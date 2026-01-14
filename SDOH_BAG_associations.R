
# 加载所有需要的R包
library(dplyr)         # 数据处理
library(lme4)          # 线性混合效应模型
library(lmerTest)      # 增强lme4，提供p值计算
library(ggplot2)       # 可视化
library(reshape2)      # 数据转换
library(stats)         # 基本统计分析
library(multcomp)      # 多重比较校正
library(data.table)    # 高效数据操作
library(mice)  # 插补缺失值的包
library(broom.mixed)  # 用于整理模型结果
library(brms)

# 读取数据
final_data <- read.csv(".../ewas_integrated4.csv", stringsAsFactors = FALSE)
dim(final_data)
colSums(is.na(final_data))  # 统计每一列的缺失值数量

# 获取所有环境变量（排除 ID 变量）
id_columns <- c("src_subject_id", "site_id_l", "rel_family_id","mri_info_manufacturer")
numeric_columns <- setdiff(names(final_data), id_columns)

# 强制转换所有非 ID 列为数值型
final_data <- final_data %>%
  mutate_at(vars(all_of(numeric_columns)), ~ as.numeric(as.character(.)))

# 检查转换后的数据类型
str(final_data)


environmental_variables <- names(final_data)[88:94]


control_variables <- c("mri_info_manufacturer",        "demo_prnt_ed_v2" ,   "demo_prtnr_ed_v2",     "asr_scr_totprob_t",      "demo_comb_income_v2")

final_data$Zcorrected_bag <- scale(final_data$corrected_bag)
final_data[environmental_variables] <- scale(final_data[environmental_variables])

# 设定结果存储路径
output_file <- ".../integrated_ewas_results_with_p_incremental3.csv"

# 如果文件存在，删除旧文件（保证新的运行不会拼接旧数据）
if (file.exists(output_file)) {
  file.remove(output_file)
}

# 遍历每个环境变量（一个一个放进去）
for (env_var in environmental_variables) {
  
  priors <- prior(normal(0, 0.1), class = "b")
  
  # 动态调整公式模板（每次只包含一个环境变量）
  formula <- bf(as.formula(
    paste(
      "Zcorrected_bag ~", env_var,  
      "+", paste(control_variables, collapse = " + "),
      "+ (1 | site_id_l/rel_family_id)"
    )
  ))
  
  # 使用 brm() 拟合模型
  model <- brm(
    formula = formula,
    data = final_data,
    family = gaussian(),  # 正态分布假设
    cores = 4,            # 降低 CPU 负担
    chains = 4,           # 运行 4 条链
    prior = priors,
    sample_prior = "yes",  # **必须开启，否则无法计算 BF**
    seed = 123,           # 设置随机种子
    iter = 4000,          # 迭代次数
    warmup = 2000,         # 预热次数
  )
  
  # 提取固定效应的结果
  fixed_effects <- as.data.frame(print(summary(model), digits = 4)$fixed)
  
  # 计算 p 值（Wald 近似）
  fixed_effects$p_value <- round(2 * pnorm(abs(fixed_effects$Estimate / fixed_effects$Est.Error), lower.tail = FALSE), 4)
  
  # 只保留环境变量对应的行
  fixed_effects_df <- fixed_effects[rownames(fixed_effects) == env_var, ]
  
  # 添加环境变量和因变量名称到结果中
  fixed_effects_df$Env_Var <- env_var
  
  # 计算贝叶斯因子（BF）
  hypothesis_test <- hypothesis(model, paste0(env_var, " = 0"))
  
  # 提取 BF 并保留 4 位小数
  fixed_effects_df$BF <- format(round(hypothesis_test$hypothesis$Evid.Ratio, 4), nsmall = 4)
  
  # **确保 `BF` 列在 `write.table()` 里**
  write.table(fixed_effects_df, file = output_file, append = TRUE, sep = ",",
              col.names = !file.exists(output_file), row.names = FALSE, quote = FALSE)
  
  # **打印当前进度**
  cat("已完成:", env_var, "\n")
}