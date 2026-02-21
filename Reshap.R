test <- read.csv("train.csv", header = TRUE, row.names = 1, check.names = FALSE)

# 转置
test_t <- t(test)

# 保存转置后的结果
write.csv(test_t, "train_transposed.csv")
