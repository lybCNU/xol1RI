library(stringr)
library(Seurat)
library(ggplot2)
library(reshape2)
library(patchwork)
library(SeuratObject)
library(dplyr)
library(tibble)
library(tidyr)


Seurat_gene<-readRDS(file="RNAseq.RDS")


valid_chromosomes <- c("I", "II", "III", "IV", "V", "X")

select_gene<-Seurat_gene[["RNA"]][[]] %>% 
  filter(!is.na(one2oneGeneID)) %>%
  # remove chromosome and one2oneChr prefix
  mutate(chromosome_clean = gsub("^(PX439|PX534)", "", chromosome),
         one2oneChr_clean = gsub("^(PX439|PX534)", "", one2oneChr)) %>%
  # rilter inter_pairs
  filter(chromosome_clean == one2oneChr_clean & chromosome_clean %in% valid_chromosomes) %>%
  # clean clean column
  select(-chromosome_clean, -one2oneChr_clean) %>%
  # Only keep ortholog pairs
  filter(gene_id %in% one2oneGeneID) %>%
  filter(one2oneGeneID %in% gene_id) %>%
  rownames()

Seurat_gene_selected<-Seurat_gene[select_gene,]

Seurat_gene_selected <- NormalizeData(Seurat_gene_selected,
                                      normalization.method	="RC",
                                      scale.factor = 1e6)


chr.read.fraction<-Seurat_gene_selected@assays$RNA@data %>%
  as.data.frame() %>% 
  group_by(chromosome = Seurat_gene_selected[["RNA"]]@meta.features$chromosome) %>%
  summarise(across(everything(), sum))  %>%
  column_to_rownames(var = "chromosome") %>%
  t() %>%
  as.data.frame()

chr.read.fraction_with_meta<-chr.read.fraction %>%
  rownames_to_column("sample") %>% 
  left_join(Seurat_gene_selected@meta.data, by = "sample")

chr.read.fraction_with_meta_long <- chr.read.fraction_with_meta %>%
  pivot_longer(cols = starts_with("PX439") | starts_with("PX534"), 
               names_to = "chromosome", values_to = "value") %>%
  select(code,chromosome,value) %>%
  mutate(chr_copy_number = case_when(
    code == "creXcre-M" & chromosome %in% c("PX439I", "PX439II", "PX439III", "PX439IV", "PX439V") ~ 2,
    code == "creXcre-F" & chromosome %in% c("PX439I", "PX439II", "PX439III", "PX439IV", "PX439V","PX439X") ~ 2,
    code == "creXcre-M" & chromosome == "PX439X" ~ 1,
    code == "claXcla-M" & chromosome %in% c("PX534I", "PX534II", "PX534III", "PX534IV", "PX534V") ~ 2,
    code == "claXcla-F" & chromosome %in% c("PX534I", "PX534II", "PX534III", "PX534IV", "PX534V","PX534X") ~ 2,
    code == "claXcla-M" & chromosome == "PX534X" ~ 1,
    code == "claXcre-M" & chromosome %in% c("PX439I", "PX439II", "PX439III", "PX439IV", "PX439V",
                                            "PX534I", "PX534II", "PX534III", "PX534IV", "PX534V",
                                            "PX439X") ~ 1,
    code == "claXcre-F" & chromosome %in% c("PX439I", "PX439II", "PX439III", "PX439IV", "PX439V",
                                            "PX534I", "PX534II", "PX534III", "PX534IV", "PX534V",
                                            "PX534X", "cbrX") ~ 1,
    #code == "claXcre-M" & chromosome == "PX534X" ~ 0,
    code == "creXcla-M" & chromosome %in% c("PX534I", "PX534II", "PX534III", "PX534IV", "PX534V",
                                            "PX439I", "PX439II", "PX439III", "PX439IV", "PX439V",
                                            "PX534X") ~ 1,
    code == "creXcla-F" & chromosome %in% c("PX439I", "PX439II", "PX439III", "PX439IV", "PX439V",
                                            "PX534I", "PX534II", "PX534III", "PX534IV", "PX534V",
                                            "PX534X", "PX439X") ~ 1,
    #code == "creXcla-M" & chromosome == "cbrX" ~ 0,
    TRUE ~ NA_real_  # 其他情况不计算
  )) %>%
  filter(!is.na(chr_copy_number)) %>%
  mutate(exp_per_copy=value/chr_copy_number)

# 定义安全的 t.test 函数
safe_t_test <- function(x) {
  result <- tryCatch({
    t.test(x)$conf.int
  }, error = function(e) {
    c(NA, NA) # 如果 t.test 失败，返回 NA
  })
  return(result)
}

# 对数据进行 log2 变换后计算统计量
chr_summary <- chr.read.fraction_with_meta_long %>%
  mutate(log2_value = log2(exp_per_copy)) %>%
  group_by(code, chromosome) %>%
  summarise(
    mean_log2 = mean(log2_value, na.rm = TRUE),
    conf.int1_log2 = safe_t_test(log2_value)[1],
    conf.int2_log2 = safe_t_test(log2_value)[2],
    sd = sd(exp_per_copy, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # 反变换回原始比例
  mutate(
    mean = 2^mean_log2,
    conf.int1 = 2^conf.int1_log2,
    conf.int2 = 2^conf.int2_log2,
    fill_col=case_when(
      startsWith(chromosome,"PX439") ~ "#666666",
      startsWith(chromosome,"PX534") ~ "#CCCCCC",
    ),
    chr=str_replace(chromosome,"PX439|PX534",""),
    Species=case_when(
      startsWith(chromosome,"PX439") ~ "C. remanei",
      startsWith(chromosome,"PX534") ~ "C. latens",
    )
  ) %>%
  select(code,Species,chr,mean,conf.int1,conf.int2,sd,fill_col) 


chr_summary$Species<-factor(chr_summary$Species,levels = c("C. remanei","C. latens"))
# 查看最终结果
print(chr_summary)

chr_summary$fill_col<-factor(chr_summary$fill_col,
                             levels = c("#CCCCCC","#666666"))


# 需要比较的组别
comparisons <- list(
  c("creXcla-M", "creXcre-M"),
  c("creXcla-M", "claXcla-M"),
  c("claXcre-M", "creXcre-M"),
  c("claXcre-M", "claXcla-M"),
  c("creXcla-F", "creXcre-F"),
  c("creXcla-F", "claXcla-F"),
  c("claXcre-F", "creXcre-F"),
  c("claXcre-F", "claXcla-F"),
  c("claXcla-F","claXcla-M"),
  c("creXcre-F","creXcre-M")
)

# 计算统计量的函数
calculate_statistics <- function(group1, group2, data) {
  result <- tryCatch({
    t_test_result <- t.test(
      log2(data$exp_per_copy[data$code == group1]),
      log2(data$exp_per_copy[data$code == group2]),
      na.rm = TRUE
    )
    conf_int_mean <- mean(t_test_result$conf.int)  # 计算置信区间的均值
    list(
      mean = mean(t_test_result$estimate),
      conf_int_mean = conf_int_mean,
      conf_int_min = t_test_result$conf.int[1],
      conf_int_max = t_test_result$conf.int[2]
    )
  }, error = function(e) {
    list(
      mean = NA,
      conf_int_mean = NA,
      conf_int_min = NA,
      conf_int_max = NA
    )
  })
  return(result)
}

# 按染色体分组计算所有比较的统计量
chr_ratio_transcription_results <- chr.read.fraction_with_meta_long %>%
  group_by(chromosome) %>%
  do({
    data <- .
    comparison_results <- lapply(comparisons, function(pair) {
      stats <- calculate_statistics(pair[1], pair[2], data)
      tibble(
        comparison = paste(pair[1], "vs", pair[2]),
        mean = 2^stats$mean,
        conf_int_mean = 2^stats$conf_int_mean,
        conf_int_min = 2^stats$conf_int_min,
        conf_int_max = 2^stats$conf_int_max
      )
    })
    bind_rows(comparison_results)
  }) %>%
  ungroup() %>%
  filter(!is.na(conf_int_mean)) %>%
  mutate(
    fill_col=case_when(
      startsWith(chromosome,"PX439") ~ "#666666",
      startsWith(chromosome,"PX534") ~ "#CCCCCC",
    ),
    chr=str_replace(chromosome,"PX439|PX534",""),
    Species=case_when(
      startsWith(chromosome,"PX439") ~ "C. remanei",
      startsWith(chromosome,"PX534") ~ "C. latens",
    )
  )

# 打印结果
print(chr_ratio_transcription_results)
chr_ratio_transcription_results$fill_col<-factor(chr_ratio_transcription_results$fill_col,
                                                 levels = c("#CCCCCC","#666666"))

k_format <- function(x){
  paste0(x/1000, "k")
}

 A<-  ggplot(chr_summary %>% filter(code=="creXcre-M")) +
  geom_bar(aes(x=chr, y=mean, fill = fill_col), stat = "identity" ) + 
  scale_fill_identity() +
  geom_errorbar(aes(x = chr, ymin = conf.int1, ymax = conf.int2), width=0.4) +
  ylim(0, 225000) +
  theme_bw(base_size = 7)+
  scale_y_continuous(labels = k_format,limits = c(0,230000)) +
  labs(title = "creXcre-M")+
  scale_x_discrete(labels = c("I","II","III","IV","V","X"))+
  theme_bw(base_size = 7)+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5),
        #axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank())+
  labs(y="Transcription level (CPM/chr.copy)")


B<-ggplot(chr_summary %>% filter(code=="claXcla-M")) +
  geom_bar(aes(x=chr, y=mean, fill = fill_col), stat = "identity" ) + 
  scale_fill_identity() +
  geom_errorbar(aes(x = chr, ymin = conf.int1, ymax = conf.int2), width=0.4) +
  ylim(0, 225000) +
  theme_bw(base_size = 7)+
  scale_y_continuous(labels = k_format,limits = c(0,230000)) +
  labs(title = "creXcre-M")+
  scale_x_discrete(labels = c("I","II","III","IV","V","X"))+
  theme_bw(base_size = 7)+
  theme(legend.position = "none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        #axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5),
        #axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank())+
  labs(y="Transcription level (CPM/chr.copy)")

C<-ggplot(chr_summary %>% filter(code=="creXcla-M"), 
       aes(chr, y = mean, fill = fill_col)) +
  geom_bar(position="dodge", stat = "identity", width = 0.9) + 
  scale_fill_identity() +
  geom_errorbar(aes(x = chr, ymin = conf.int1, ymax = conf.int2), width=0.4, position = position_dodge(width = 0.9)) +
  theme_bw(base_size = 7)+
  labs(title = "creXcla-M")+
  scale_y_continuous(labels = k_format,limits = c(0,230000))+
  theme(legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5))

D<-ggplot(chr_summary %>% filter(code=="claXcre-M"), 
       aes(chr, y = mean, fill = fill_col)) +
  geom_bar(position="dodge", stat = "identity", width = 0.9) + 
  scale_fill_identity() +
  geom_errorbar(aes(x = chr, ymin = conf.int1, ymax = conf.int2), width=0.4, position = position_dodge(width = 0.9)) +
  theme_bw(base_size = 7)+
  labs(title = "claXcre-M")+
  scale_y_continuous(labels = k_format,limits = c(0,230000))+
  theme(legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5))

E<-ggplot(chr_ratio_transcription_results %>% 
            filter(str_detect(comparison,"creXcla-M")) %>% 
            filter(str_detect(chromosome,"PX534")),
       aes(x = chr, y = conf_int_mean, fill = fill_col))+
  scale_fill_identity() +
  scale_y_continuous(breaks = c(0.0, 0.5, 1.0))+
  geom_bar(position="dodge", stat = "identity", width = 0.9) + 
  geom_errorbar(aes(x = chr, ymin = conf_int_min, ymax = conf_int_max), width=0.4, position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")+
  theme_bw(base_size = 7)+
  theme(legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5))+
  labs(y="Ratio of Transciption level")

F<-ggplot(chr_ratio_transcription_results %>% 
            filter(str_detect(comparison,"claXcre-M")) %>% 
            filter(str_detect(chromosome,"PX439")),
       aes(x = chr, y = conf_int_mean, fill = fill_col))+
  scale_fill_identity() +
  scale_y_continuous(breaks = c(0.0, 0.5, 1.0))+
  geom_bar(position="dodge", stat = "identity", width = 0.9) + 
  geom_errorbar(aes(x = chr, ymin = conf_int_min, ymax = conf_int_max), width=0.4, position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red")+
  theme_bw(base_size = 7)+
  theme(legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5))+
  labs(y="Ratio of Transciption level")

(B+C+E)/(A+D+F)
