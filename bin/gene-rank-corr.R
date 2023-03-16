library(tidyverse)
library(data.table)
library(survcomp)
library(broom)
library(ggpubr)
library(ggrepel)
library(plotly)
library(PerformanceAnalytics)

hum.gwa.ocd = read_tsv("~/Dropbox (UMass Medical School)/Papers/cross-species/git/cross-species-gwas/dat/gwa/hum/bat/OCD-TS.NORDiC.region-test.genes.gene.assoc.mbat")

# Set the directory path
dir_path <- "~/Dropbox (UMass Medical School)/Papers/cross-species/git/cross-species-gwas/dat/gwa/dog/bat/fastBAT/"

# Get a list of files in the directory
file_list <- list.files(path = dir_path, full.names = TRUE)

# Read in all files in the list using map_dfr()
data_list <- map_dfr(file_list, read_tsv, .id = .x)

# Pivot wider
comp = data_list %>%
  select(dog.phe= PHE,
         gene= Gene,
         dog.nsnp= No.SNPs,
         dog.chisq= `Chisq(Obs)`,
         dog.p= Pvalue) %>%
  merge((hum.gwa.ocd %>%
           select(gene = Gene,
                  hum.nsnp = No.SNPs,
                  hum.chisq = Chisq_fastBAT,
                  hum.p = P_fastBAT)),
        by = "gene",
        all.x = T) %>%
  arrange(dog.phe,gene)

# Plot
p = comp %>%
  filter(dog.phe=="fa.13-filled-by-mean") %>%
  ggplot(aes(x = -log10(dog.p),
             y = -log10(hum.p),
             label = gene)) +
  #facet_wrap(.~dog.phe) +
  geom_point(alpha=0.25) +
  theme_pubr()

ggplotly(p)
  
# Kendall's tau rank correlation: concordance vs discordance
test.p.kendall.corr = comp %>%
  group_by(dog.phe) %>%
  do(cor_test = cor.test(.$dog.p, 
                         .$hum.p,
                         method = "kendall") %>% tidy()) %>%
  unnest(cor_test)

# Fisher's combined probability
test.fish = comp %>%
  group_by(dog.phe) %>%
  rowwise() %>%
  mutate(comb.p = combine.test(p = c(dog.p,hum.p),
                               weight = c(dog.nsnp,hum.nsnp),
                               na.rm = T)) %>%
  mutate(comb.p.adj = p.adjust(comb.p, method = "BH")) %>%
  ungroup() %>%
  merge((test.p.kendall.corr %>%
           select(dog.phe,
                  cor.r=estimate,
                  cor.p=p.value)), 
        by = "dog.phe", 
        all.x = T) %>%
  mutate(all.p.adj = p.adjust(comb.p, method = "BH")) %>%
  mutate(all.cor.p = p.adjust(cor.p, method = "BH"))


test.fish %>%
  group_by(dog.phe) %>%
  arrange(comb.p.adj) %>% 
  slice_head(n=10) %>%
  ggplot(aes(x = gene,
             y = -log10(comb.p),
             fill = comb.p.adj < 0.05)) +
  facet_grid(dog.phe ~ .,
             scales = "free_y") +
  coord_flip() +
  geom_col() +
  scale_fill_discrete(type = greenfocus[1:2]) +
  theme_classic() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.position = "none")

for (i in 1:length(unique(test.fish$dog.phe))){
p = test.fish %>%
  filter(dog.phe == unique(test.fish$dog.phe)[i]) %>%
  arrange(comb.p) %>% 
  slice_head(n=25) %>%
  ggplot(aes(x = reorder(gene, comb.p.adj),
             y = -log10(comb.p),
             fill = (dog.p < 0.05 & hum.p < 0.05))) +
  coord_flip() +
  geom_col() +
  scale_fill_discrete(type = rev(greenfocus[1:2])) +
  labs(x = "region", title = unique(test.fish$dog.phe)[i]) +
  theme_classic() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))
ggsave(filename = paste("~/Dropbox (UMass Medical School)/Papers/cross-species/fig/pngs/",
                        unique(test.fish$dog.phe)[i],
                        ".fastBAT.png"),
       plot = p,
       width = 5,
       height = 5,
       device = "png")
}

p = test.p.kendall.corr %>% 
  mutate(p.adj = p.adjust(p = p.value, method = "BH")) %>%
  filter(abs(estimate)>0.01) %>%
  ggplot(aes(x = reorder(dog.phe,estimate), 
             y = estimate, 
             fill = p.adj < 0.05,
             label = scientific(p.adj))) + 
  coord_flip() +
  geom_col() +
  geom_text() +
  labs(y = "Kendall's tau", x = "dog phenotype") +
  scale_fill_discrete(type = rev(greenfocus[1:2])) +
  theme_classic()

ggsave(plot=p,filename = "~/Dropbox (UMass Medical School)/Thesis/test_kendall_genecor.png",width=6,height=6)


ggsave(plot=p,
       filename = "~/Dropbox (UMass Medical School)/Thesis/test_kendall_genecor_NORDiC.pdf",
       width=5,
       height=6,
       device='pdf',
       units='in')

test.fish = comp %>%
  filter(dog.phe %in% (test.p.kendall.corr %>%
                         mutate(p.adj = p.adjust(p = p.value, method = "BH")) %>%
                         filter(p.adj < 0.05) %>%
                         pull(dog.phe))) %>%
  group_by(dog.phe) %>%
  rowwise() %>%
  mutate(comb.p = combine.test(p = c(dog.p,hum.p),
                               weight = c(dog.nsnp,hum.nsnp),
                               na.rm = T)) %>%
  mutate(comb.p.adj = p.adjust(comb.p, method = "BH")) %>%
  ungroup() %>%
  merge((test.p.kendall.corr %>%
           select(dog.phe,
                  cor.r=estimate,
                  cor.p=p.value)), 
        by = "dog.phe", 
        all.x = T) %>%
  mutate(all.p.adj = p.adjust(comb.p, method = "BH")) %>%
  mutate(all.cor.p = p.adjust(cor.p, method = "BH"))

write_tsv(test.fish %>%
            arrange(comb.p.adj),
          "~/Dropbox (UMass Medical School)/Thesis/cross-species.kendall-corr.fisher-p.tsv")

p = test.fish %>%
  group_by(dog.phe) %>%
  arrange(comb.p.adj) %>% 
  slice_head(n=10) %>%
  ggplot(aes(x = gene,
             y = -log10(comb.p),
             fill = comb.p.adj < 0.05)) +
  facet_grid(dog.phe ~ .,
             scales = "free_y") +
  coord_flip() +
  geom_col() +
  scale_fill_discrete(type = greenfocus[1:2]) +
  theme_classic() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.position = "none")
