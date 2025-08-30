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

#select 1vs1_only, no_inter_pairs and I II III IV V X
clean.gene <- read.delim("Data_S1.txt")
select_gene<-Seurat_gene[["RNA"]][[]] %>% 
  filter(!is.na(one2oneGeneID)) %>%
  # remove chromosome and one2oneChr cbr/cni prefix
  mutate(chromosome_clean = gsub("^(cbr|cni)", "", chromosome),
         one2oneChr_clean = gsub("^(cbr|cni)", "", one2oneChr)) %>%
  # 过滤 inter_pairs
  filter(chromosome_clean == one2oneChr_clean & chromosome_clean %in% valid_chromosomes) %>%
  # clean
  select(-chromosome_clean, -one2oneChr_clean) %>%
  # Only keep ortholog pairs whose both members have <1% inter-specific mapping reads
  filter(rownames(.) %in% rownames(clean.gene)) %>%
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
  left_join(Seurat_gene_selected@meta.data %>% rownames_to_column("sample"), by = "sample")

chr.read.fraction_with_meta_long <- chr.read.fraction_with_meta %>%
  pivot_longer(cols = starts_with("cbr") | starts_with("cni"), 
               names_to = "chromosome", values_to = "value") %>%
  select(code,chromosome,value) %>%
  mutate(chr_copy_number = case_when(
    code == "cbrXcbr-M" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV") ~ 2,
    code == "cbrXcbr-F" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV","cbrX") ~ 2,
    code == "cbrXcbr-M" & chromosome == "cbrX" ~ 1,
    code == "cniXcni-M" & chromosome %in% c("cniI", "cniII", "cniIII", "cniIV", "cniV") ~ 2,
    code == "cniXcni-F" & chromosome %in% c("cniI", "cniII", "cniIII", "cniIV", "cniV","cniX") ~ 2,
    code == "cniXcni-M" & chromosome == "cniX" ~ 1,
    code == "cniXcbr-M" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                            "cniI", "cniII", "cniIII", "cniIV", "cniV",
                                            "cbrX") ~ 1,
    code == "cniXcbr-F" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                            "cniI", "cniII", "cniIII", "cniIV", "cniV",
                                            "cniX", "cbrX") ~ 1,
    #code == "cniXcbr-M" & chromosome == "cniX" ~ 0,
    code == "cbrXcni-M" & chromosome %in% c("cniI", "cniII", "cniIII", "cniIV", "cniV",
                                            "cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                            "cniX") ~ 1,
    code == "cbrXcni-F" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                            "cniI", "cniII", "cniIII", "cniIV", "cniV",
                                            "cniX", "cbrX") ~ 1,
    code == "cniXIL3960XcbrXIL31006-F" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                                           "cniI", "cniII", "cniIII", "cniIV", "cniV",
                                                           "cniX", "cbrX") ~ 1,
    code == "cniXIL3960XcbrXIL31006-M" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                                           "cniI", "cniII", "cniIII", "cniIV", "cniV",
                                                           "cbrX") ~ 1,
    #code == "cbrXcni-M" & chromosome == "cbrX" ~ 0,
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
      startsWith(chromosome,"cbr") ~ "#666666",
      startsWith(chromosome,"cni") ~ "#CCCCCC",
    ),
    chr=str_replace(chromosome,"cbr|cni",""),
    Species=case_when(
      startsWith(chromosome,"cbr") ~ "C. briggsae",
      startsWith(chromosome,"cni") ~ "C. nigoni",
    )
  ) %>%
  select(code,Species,chr,mean,conf.int1,conf.int2,sd,fill_col) 

chr_summary$Species<-factor(chr_summary$Species,levels = c("C. nigoni","C. briggsae"))
# 查看最终结果
print(chr_summary)

chr_summary$fill_col<-factor(chr_summary$fill_col,
                             levels = c("#CCCCCC","#666666"))

# 需要比较的组别
comparisons <- list(
  c("cbrXcni-M", "cbrXcbr-M"),
  c("cbrXcni-M", "cniXcni-M"),
  c("cniXcbr-M", "cbrXcbr-M"),
  c("cniXcbr-M", "cniXcni-M"),
  c("cbrXcni-F", "cbrXcbr-F"),
  c("cbrXcni-F", "cniXcni-F"),
  c("cniXcbr-F", "cbrXcbr-F"),
  c("cniXcbr-F", "cniXcni-F"),
  c("cniXIL3960XcbrXIL31006-M","cbrXcbr-M"),
  c("cniXIL3960XcbrXIL31006-M", "cniXcni-M"),
  c("cniXcni-F","cniXcni-M"),
  c("cbrXcbr-F","cbrXcbr-M")
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
      startsWith(chromosome,"cbr") ~ "#666666",
      startsWith(chromosome,"cni") ~ "#CCCCCC",
    ),
    chr=str_replace(chromosome,"cbr|cni",""),
    Species=case_when(
      startsWith(chromosome,"cbr") ~ "C. briggsae",
      startsWith(chromosome,"cni") ~ "C. nigoni",
    )
  )

# 打印结果
print(chr_ratio_transcription_results)
chr_ratio_transcription_results$fill_col<-factor(chr_ratio_transcription_results$fill_col,
                                                 levels = c("#CCCCCC","#666666"))

k_format <- function(x){
  paste0(x/1000, "k")
}

fig1C<-ggplot(chr_summary %>% filter(code=="cniXcni-M")) +
  geom_bar(aes(x=chr, y=mean, fill = fill_col), stat = "identity" ) + 
  scale_fill_identity() +
  geom_errorbar(aes(x = chr, ymin = conf.int1, ymax = conf.int2), width=0.4) +
  #geom_errorbar(aes(x = chr, ymin = mean-sd, ymax = mean+sd), width=0.4) +
  ylim(0, 225000) +
  theme_bw(base_size = 7)+
  scale_y_continuous(labels = k_format,limits = c(0,230000)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5)
  )+
  labs(title = "cniXcni-M",y="Transcription level (CPM/chr.copy)")+
  scale_x_discrete(labels = c("I","II","III","IV","V","X"))+
  theme_bw(base_size = 7)+
  theme(legend.position = "none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        
        #axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank()
        #title = element_blank()
  )


fig1D<-ggplot(chr_summary %>% filter(code=="cbrXcbr-M")) +
  geom_bar(aes(x=chr, y=mean, fill = fill_col), stat = "identity" ) + 
  scale_fill_identity() +
  geom_errorbar(aes(x = chr, ymin = conf.int1, ymax = conf.int2), width=0.4) +
  ylim(0, 225000) +
  theme_bw(base_size = 7)+
  scale_y_continuous(labels = k_format,limits = c(0,230000)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5)
  )+
  labs(title = "cbrXcbr-M",y="Transcription level (CPM/chr.copy)")+
  scale_x_discrete(labels = c("I","II","III","IV","V","X"))+
  theme_bw(base_size = 7)+
  theme(legend.position = "none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        
        #axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank()
        #title = element_blank()
  )

fig1E<-ggplot(chr_summary %>% filter(code=="cbrXcni-M"), 
              aes(chr, y = mean, fill = fill_col)) +
  geom_bar(position="dodge", stat = "identity", width = 0.9) + 
  scale_fill_identity() +
  geom_errorbar(aes(x = chr, ymin = conf.int1, ymax = conf.int2), width=0.4, position = position_dodge(width = 0.9)) +
  theme_bw(base_size = 7)+
  labs(title = "cbrXcni-M")+
  scale_y_continuous(labels = k_format,limits = c(0,230000))+
  theme(legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank()
        #axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5)
  )

fig1F<-ggplot(chr_summary %>% filter(code=="cniXcbr-M") , 
              aes(chr, y = mean, fill = fill_col)) +
  geom_bar(position="dodge", stat = "identity", width = 0.9) + 
  #geom_hline(yintercept = 1, linetype = "dashed", color = "red")+
  scale_fill_identity() +
  geom_errorbar(aes(x = chr, ymin = conf.int1, ymax = conf.int2), width=0.4, position = position_dodge(width = 0.9)) +
  theme_bw(base_size = 7)+
  labs(title = "cniXcbr-M")+
  scale_y_continuous(labels = k_format,limits = c(0,230000))+
  theme(legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank()
        #axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5)
  )


fig1H<-ggplot(chr_ratio_transcription_results %>% filter(str_detect(comparison,"cniXcbr-M")) %>% filter(str_detect(Species,"C. briggsae")),
              aes(x = chr, y = conf_int_mean, fill = fill_col))+
  scale_fill_identity() +
  scale_y_continuous(breaks = c(0.0, 0.5, 1.0))+
  geom_bar(position="dodge", stat = "identity", width = 0.9) + 
  geom_errorbar(aes(x = chr, ymin = conf_int_min, ymax = conf_int_max), width=0.4, position = position_dodge(width = 0.9)) +
  theme_bw(base_size = 7)+
  theme(legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5))+
  labs(y="Ratio of Transciption level")+
  geom_hline(yintercept = 0.5, linetype = "dashed")

fig1G<-ggplot(chr_ratio_transcription_results %>% filter(str_detect(comparison,"cbrXcni-M")) %>% filter(str_detect(Species,"C. nigoni")),
              aes(x = chr, y = conf_int_mean, fill = fill_col))+
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")+
  scale_fill_identity() +
  geom_bar(position="dodge", stat = "identity", width = 0.9) + 
  scale_y_continuous(breaks = c(0.0, 0.5, 1.0))+
  geom_errorbar(aes(x = chr, ymin = conf_int_min, ymax = conf_int_max), width=0.4, position = position_dodge(width = 0.9)) +
  theme_bw(base_size = 7)+
  theme(legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5))+
  labs(y="Ratio of Transciption level")+
  geom_hline(yintercept = 1, linetype = "dashed")

(fig1C+fig1E+fig1G)/(fig1D+fig1F+fig1H)



figS2A<-ggplot(chr_summary %>% filter(code=="cniXcni-F")) +
  geom_bar(aes(x=chr, y=mean, fill = fill_col), stat = "identity" ) + 
  scale_fill_identity() +
  geom_errorbar(aes(x = chr, ymin = conf.int1, ymax = conf.int2), width=0.4) +
  #geom_errorbar(aes(x = chr, ymin = mean-sd, ymax = mean+sd), width=0.4) +
  ylim(0, 225000) +
  theme_bw(base_size = 7)+
  scale_y_continuous(labels = k_format,limits = c(0,130000)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5)
  )+
  labs(title = "cniXcni-F",y="Transcription level (CPM/chr.copy)")+
  scale_x_discrete(labels = c("I","II","III","IV","V","X"))+
  theme_bw(base_size = 7)+
  theme(legend.position = "none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        
        #axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank()
        #title = element_blank()
  )


figS2B<-ggplot(chr_summary %>% filter(code=="cbrXcbr-F")) +
  geom_bar(aes(x=chr, y=mean, fill = fill_col), stat = "identity" ) + 
  scale_fill_identity() +
  geom_errorbar(aes(x = chr, ymin = conf.int1, ymax = conf.int2), width=0.4) +
  ylim(0, 225000) +
  theme_bw(base_size = 7)+
  scale_y_continuous(labels = k_format,limits = c(0,130000)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.title.y = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5)
  )+
  labs(title = "cbrXcbr-F",y="Transcription level (CPM/chr.copy)")+
  scale_x_discrete(labels = c("I","II","III","IV","V","X"))+
  theme_bw(base_size = 7)+
  theme(legend.position = "none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        
        #axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank()
        #title = element_blank()
  )

figS2C<-ggplot(chr_summary %>% filter(code=="cbrXcni-F"), 
               aes(chr, y = mean, fill = fill_col)) +
  geom_bar(position="dodge", stat = "identity", width = 0.9) + 
  scale_fill_identity() +
  geom_errorbar(aes(x = chr, ymin = conf.int1, ymax = conf.int2), width=0.4, position = position_dodge(width = 0.9)) +
  theme_bw(base_size = 7)+
  labs(title = "cbrXcni-F")+
  scale_y_continuous(labels = k_format,limits = c(0,130000))+
  theme(legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank()
        #axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5)
  )

figS2D<-ggplot(chr_summary %>% filter(code=="cniXcbr-F") , 
               aes(chr, y = mean, fill = fill_col)) +
  geom_bar(position="dodge", stat = "identity", width = 0.9) + 
  #geom_hline(yintercept = 1, linetype = "dashed", color = "red")+
  scale_fill_identity() +
  geom_errorbar(aes(x = chr, ymin = conf.int1, ymax = conf.int2), width=0.4, position = position_dodge(width = 0.9)) +
  theme_bw(base_size = 7)+
  labs(title = "cniXcbr-F")+
  scale_y_continuous(labels = k_format,limits = c(0,130000))+
  theme(legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank()
        #axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5)
  )


figS2F<-ggplot(chr_ratio_transcription_results %>% filter(str_detect(comparison,"cniXcbr-F")) %>% filter(str_detect(Species,"C. briggsae")),
               aes(x = chr, y = conf_int_mean, fill = fill_col))+
  scale_fill_identity() +
  scale_y_continuous(breaks = c(0.0, 0.5, 1.0))+
  geom_bar(position="dodge", stat = "identity", width = 0.9) + 
  geom_errorbar(aes(x = chr, ymin = conf_int_min, ymax = conf_int_max), width=0.4, position = position_dodge(width = 0.9)) +
  theme_bw(base_size = 7)+
  theme(legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5))+
  labs(y="Ratio of Transciption level")+
  geom_hline(yintercept = 0.5, linetype = "dashed")

figS2E<-ggplot(chr_ratio_transcription_results %>% filter(str_detect(comparison,"cbrXcni-F")) %>% filter(str_detect(Species,"C. nigoni")),
               aes(x = chr, y = conf_int_mean, fill = fill_col))+
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = "red")+
  scale_fill_identity() +
  geom_bar(position="dodge", stat = "identity", width = 0.9) + 
  scale_y_continuous(breaks = c(0.0, 0.5, 1.0))+
  geom_errorbar(aes(x = chr, ymin = conf_int_min, ymax = conf_int_max), width=0.4, position = position_dodge(width = 0.9)) +
  theme_bw(base_size = 7)+
  theme(legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5))+
  labs(y="Ratio of Transciption level")+
  geom_hline(yintercept = 1, linetype = "dashed")

(figS2A+figS2C+figS2E)/(figS2B+figS2D+figS2F)
