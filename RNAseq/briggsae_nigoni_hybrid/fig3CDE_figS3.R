library(stringr)
library(Seurat)
library(ggplot2)
library(reshape2)
library(patchwork)
library(SeuratObject)
library(dplyr)
library(tibble)
library(tidyr)

Seurat_gene<-readRDS(file = "RNAseq.RDS")

valid_chromosomes <- c("I", "II", "III", "IV", "V", "X")

#select 1vs1_only, no_inter_pairs and I II III IV V X
clean.gene <- read.delim("Data_S1.txt")
select_gene<-Seurat_gene[["RNA"]][[]] %>% 
  filter(!is.na(one2oneGeneID)) %>%
  # 去除chromosome和one2oneChr的前缀cbr或cni
  mutate(chromosome_clean = gsub("^(cbr|cni)", "", chromosome),
         one2oneChr_clean = gsub("^(cbr|cni)", "", one2oneChr)) %>%
  # 过滤 inter_pairs
  filter(chromosome_clean == one2oneChr_clean & chromosome_clean %in% valid_chromosomes) %>%
  # 删除临时生成的clean列
  select(-chromosome_clean, -one2oneChr_clean) %>%
  # Only keep ortholog pairs whose both members have <1% inter-specific mapping reads
  filter(rownames(.) %in% rownames(clean.gene)) %>%
  # Only keep ortholog pairs
  filter(gene_id %in% one2oneGeneID) %>%
  filter(one2oneGeneID %in% gene_id) %>%
  # 输出行名
  rownames()


Seurat_gene_selected<-Seurat_gene[select_gene,]

Seurat_gene_selected <- NormalizeData(Seurat_gene_selected,
                                      normalization.method	="RC",
                                      scale.factor = 1e6)

chr.read.fraction<-Seurat_gene_selected@assays$RNA@data %>%
  as.data.frame() %>% 
  group_by(chromosome = Seurat_gene_selected[["RNA"]]@meta.features$chromosome_IR_class) %>%
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
    code == "cbrXcbr-F" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV","cbrIR10330","cbrIR10337","cbrnon-IR") ~ 2,
    code == "cbrXcbr-M" & chromosome %in% c("cbrIR10330","cbrIR10337","cbrnon-IR") ~ 1,
    code == "cniXcni-M" & chromosome %in% c("cniI", "cniII", "cniIII", "cniIV", "cniV") ~ 2,
    code == "cniXcni-F" & chromosome %in% c("cniI", "cniII", "cniIII", "cniIV", "cniV","cniIR10330","cniIR10337","cninon-IR") ~ 2,
    code == "cniXcni-M" & chromosome %in% c("cniIR10330","cniIR10337","cninon-IR") ~ 1,
    code == "cniXcbr-M" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                            "cniI", "cniII", "cniIII", "cniIV", "cniV",
                                            "cbrIR10330","cbrIR10337","cbrnon-IR") ~ 1,
    code == "cniXcbr-F" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                            "cniI", "cniII", "cniIII", "cniIV", "cniV",
                                            "cniIR10330","cniIR10337","cninon-IR", "cbrIR10330","cbrIR10337","cbrnon-IR") ~ 1,
    #code == "cniXcbr-M" & chromosome == "cniX" ~ 0,
    code == "cbrXcni-M" & chromosome %in% c("cniI", "cniII", "cniIII", "cniIV", "cniV",
                                            "cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                            "cniIR10330","cniIR10337","cninon-IR") ~ 1,
    code == "cbrXcni-F" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                            "cniI", "cniII", "cniIII", "cniIV", "cniV",
                                            "cniIR10330","cniIR10337","cninon-IR", "cbrIR10330","cbrIR10337","cbrnon-IR") ~ 1,
    code == "cniXIL3960XcbrXIL31006-F" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                                           "cniI", "cniII", "cniIII", "cniIV", "cniV",
                                                           "cniIR10330","cniIR10337","cninon-IR", "cbrIR10330","cbrIR10337","cbrnon-IR") ~ 1,
    code == "cniXIL3960XcbrXIL31006-M" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                                           "cniI", "cniII", "cniIII", "cniIV", "cniV",
                                                           "cbrIR10330","cbrIR10337","cbrnon-IR") ~ 1,
    code == "cbrXzzy10330-zzy10330F" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                                         "cniI", "cniII", "cniIII", "cniIV", "cniV",
                                                         "cniIR10337","cninon-IR","cbrIR10337","cbrnon-IR") ~ 1,
    code == "cbrXzzy10330-zzy10330F" & chromosome %in% c("cbrIR10330") ~ 2,
    code == "cbrXzzy10330-zzy10330M" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                                         "cniI", "cniII", "cniIII", "cniIV", "cniV",
                                                         "cniIR10337","cninon-IR", "cbrIR10330") ~ 1,
    code == "cbrXzzy10337-zzy1033F" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                                        "cniI", "cniII", "cniIII", "cniIV", "cniV",
                                                        "cniIR10330","cninon-IR", "cbrIR10330","cbrnon-IR") ~ 1,
    code == "cbrXzzy10330-zzy10337F" & chromosome %in% c("cbrIR10337") ~ 2,
    code == "cbrXzzy10337-zzy10337M" & chromosome %in% c("cbrI", "cbrII", "cbrIII", "cbrIV", "cbrV",
                                                         "cniI", "cniII", "cniIII", "cniIV", "cniV",
                                                         "cniIR10330","cninon-IR","cbrIR10337") ~ 1,    
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
  c("cniXcbr-M", "cbrXcbr-M"),
  c("cniXcbr-M", "cniXcni-M"),
  c("cbrXzzy10337-zzy10337M","cbrXcbr-M"),
  c("cbrXzzy10337-zzy10337M","cniXcni-M"),
  c("cbrXzzy10330-zzy10330M","cbrXcbr-M"),
  c("cbrXzzy10330-zzy10330M","cniXcni-M")
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
chr_ratio_transcription_results$chr<-factor(chr_ratio_transcription_results$chr,levels = c("I","II","III","IV","V","non-IR","IR10337","IR10330"))

figure2C<-ggplot(chr_ratio_transcription_results %>% 
                   filter(str_detect(comparison,"cbrXzzy10337-zzy10337M")),
                 aes(x = chr, y = conf_int_mean, fill = fill_col))+
  scale_fill_identity() +
  scale_y_continuous(breaks = c(0.0,0.25, 0.5,0.75, 1.0,1.25),limits = c(0,1.3))+
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
  geom_hline(yintercept = 1, linetype = "dashed")

figure2D<-ggplot(chr_ratio_transcription_results %>% 
                   filter(str_detect(comparison,"cbrXzzy10330-zzy10330M")),
                 aes(x = chr, y = conf_int_mean, fill = fill_col))+
  scale_fill_identity() +
  scale_y_continuous(breaks = c(0.0,0.25, 0.5,0.75, 1.0,1.25),limits = c(0,1.3))+
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
  geom_hline(yintercept = 1, linetype = "dashed")

figure2E<-ggplot(chr_ratio_transcription_results %>% 
                   filter(str_detect(comparison,"cniXcbr-M")),
                 aes(x = chr, y = conf_int_mean, fill = fill_col))+
  scale_fill_identity() +
  scale_y_continuous(breaks = c(0.0,0.25, 0.5,0.75, 1.0,1.25),limits = c(0,1.3))+
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

k_format <- function(x){
  paste0(x/1000, "k")
}

chr_summary$chr<-factor(chr_summary$chr,levels = c("I","II","III","IV","V","non-IR","IR10337","IR10330"))
figS3A<-ggplot(chr_summary %>% filter(code=="cbrXzzy10337-zzy10337M"), 
               aes(chr, y = mean, fill = fill_col)) +
  geom_bar(position="dodge", stat = "identity", width = 0.9) + 
  scale_fill_identity() +
  geom_errorbar(aes(x = chr, ymin = conf.int1, ymax = conf.int2), width=0.4, position = position_dodge(width = 0.9)) +
  theme_bw(base_size = 7)+
  labs(title = "cbrXzzy10337-zzy10337M",y="Transcription level (CPM/chr.copy)")+
  scale_y_continuous(labels = k_format,limits = c(0,230000))+
  theme(legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.title.y = element_blank()
        #axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5)
  )


figS3B<-ggplot(chr_summary %>% filter(code=="cbrXzzy10330-zzy10330M"), 
               aes(chr, y = mean, fill = fill_col)) +
  geom_bar(position="dodge", stat = "identity", width = 0.9) + 
  scale_fill_identity() +
  geom_errorbar(aes(x = chr, ymin = conf.int1, ymax = conf.int2), width=0.4, position = position_dodge(width = 0.9)) +
  theme_bw(base_size = 7)+
  labs(title = "cbrXzzy10330-zzy10330M",y="Transcription level (CPM/chr.copy)")+
  scale_y_continuous(labels = k_format,limits = c(0,230000))+
  theme(legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        legend.text=element_text(size=6),
        legend.key.size = unit(0.2, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.title.y = element_blank()
        #axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5)
  )

figS3C<-ggplot(chr_summary %>% filter(code=="cbrXcbr-M")) +
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
  labs(title = "cbrXcbr-M",y="Transcription level (CPM/chr.copy)")+
  theme_bw(base_size = 7)+
  theme(legend.position = "none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        
        #axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank()
        #title = element_blank()
  )

figS3D<-ggplot(chr_summary %>% filter(code=="cniXcni-M")) +
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
  theme_bw(base_size = 7)+
  theme(legend.position = "none",panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        
        #axis.text.y = element_blank(),
        # axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_blank()
        #title = element_blank()
  )

design <-"ABCD"
figS3<-figS3A+figS3B+figS3C+figS3D +
  plot_layout(design = design)

figure2CDE<-figure2C+figure2D+figure2E

figure2CDE
figS3
