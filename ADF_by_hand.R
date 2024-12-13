#Augmented Dickey-Fuller test by hand.
k <- 6 ##設定6:(會有五個落後期做為差分後的解釋變數)

####失業率####
adf.test(x3)
x3 <- x[, 3]
x3
adf <- function(x3, k=NULL) {   #處理進行ADF的序列資料和滯後期數K
  
  if(is.null(k)) {
    k <- trunc((length(x3) - 1)^(1/3))
    k <- k + 1
  }
  
  if(k < 0) {
    k <- trunc((length(x3) - 1)^(1/3))
    k <- k + 1
  }}
adf_x3 <- diff(x3)
n <- length(adf_x3)
z <- embed(adf_x3, k)  #轉換序列資料 by設定的lag order
diff_unemt <- z[, 1]  ##被解釋變數(差分後)
unem_lag <- x3[k:n] #解釋變數:被解釋變數差分前k:n期的資料
tt <- k:n #時間
if(k > 1){
  d_unem_lag <- z[, 2:k] #解釋變數:差分過後的落後期資料(5個落後期)
  adf_x <- cbind( unem_lag, "inter"=1, tt, d_unem_lag )
  adf_be <- solve(t(adf_x) %*% adf_x) %*% t(adf_x) %*% diff_unemt
} else {adf_x <- cbind(unem_lag, "inter"=1, tt)
adf_be <- solve(t(adf_x) %*% adf_x) %*% t(adf_x) %*% diff_unemt
}
#殘差
adf_eh <- diff_unemt - adf_x %*% adf_be #
residual_variance <- sum(adf_eh^2) / (length(tt) - ncol(adf_x)) #殘差變異數
XTX_inverse <- solve(t(adf_x) %*% adf_x)
#unem_lag標準誤
element_11 <- XTX_inverse[1, 1] 
unemlag_se <- sqrt(residual_variance * element_11)
cat("unem落後期係數的標準誤差:", unemlag_se, "\n")
STAT <- adf_be[1, 1]/unemlag_se
TBL <- cbind(
  "1%" = c(-4.38, -4.15, -4.04, -3.99, -3.98, -3.96),
  "5%" = c(-3.60, -3.50, -3.45, -3.43, -3.42, -3.41))
rownames(TBL) <- paste0("T = ", c(25, 50, 100, 250, 500, Inf))
cat(strrep("#", 22))
cat("\ncritical value table\n")
print(TBL)
cat(strrep("#", 22))
cat("\n")
cat("T = ", n, "\n\n")
#Which:返回向量
num <- which(n < c(25, 50, 100, 250, 500, Inf))[1] #看樣本數滿足哪個條件
num2 <- which(STAT < TBL[num, ]) #看是否小於相應表格中的值
if(length(num2) == 0) {  
  cat("cannot reject H0: time series unem is non-stationary\n")
  cat("STAT: ", STAT, "\n")
} else {
  num3 <- num2[1]
  cat("REJECT H0: time series is stationary\n")
  cat("with significant level", colnames(TBL)[num3], "\n")
  cat("STAT: ", STAT, "\n")
}

##有套件的部分
library(tseries)
if (k > 1) {
  yt1 <- z[, 2:k] ##差分後落後期資料
  res <- lm(diff_unemt  ~ unem_lag + 1 + tt + d_unem_lag, drop.unused.levels = TRUE)
} else res <- lm(diff_unemt ~ unem_lag + 1 + tt)
res_summary <- summary(res)
res_summary


####犯罪率####
ly <- log(y)
adf <- function(ly, k=NULL) {   #處理進行ADF的序列資料和滯後期數K
  
  if(is.null(k)) {
    k <- trunc((length(ly) - 1)^(1/3))
    k <- k + 1
  }
  
  if(k < 0) {
    k <- trunc((length(ly) - 1)^(1/3))
    k <- k + 1
  }}
adf_ly <- diff(ly)
n <- length(adf_ly)
z1 <- embed(adf_ly, k)  #轉換序列資料w/落後期數
diff_cri <- z1[, 1] #ADF的被解釋變數:差分後
cri_lag <- ly[k:n] #差分前k期到n期的資料
tt <- k:n #時間
if(k > 1){
  d_cri_lag <- z1[, 2:k]
  adf_x1 <- cbind( cri_lag, "inter"=1, tt, d_cri_lag )
  adf_be1 <- solve(t(adf_x1) %*% adf_x1) %*% t(adf_x1) %*% diff_cri
} else {adf_x1 <- cbind(cri_lag, "inter"=1, tt)
adf_be1 <- solve(t(adf_x1) %*% adf_x1) %*% t(adf_x1) %*% diff_cri
}##係數結果相同
#殘差
adf_ehcri <- diff_cri - adf_x1 %*% adf_be1 #same #殘差
ecri_variance <- sum(adf_ehcri^2) / (length(tt) - ncol(adf_x1))
XTX_inverse <- solve(t(adf_x1) %*% adf_x1)
#cri標準誤
element_11 <- XTX_inverse[1, 1]
cri_lag_se <- sqrt(ecri_variance * element_11)
# 印出結果
cat("cri落後期係數的標準誤差:", cri_lag_se, "\n")
STAT2 <- adf_be1[1, 1]/cri_lag_se
STAT2
TBL <- cbind(
  "1%" = c(-4.38, -4.15, -4.04, -3.99, -3.98, -3.96),
  "5%" = c(-3.60, -3.50, -3.45, -3.43, -3.42, -3.41))
rownames(TBL) <- paste0("T = ", c(25, 50, 100, 250, 500, Inf))
cat(strrep("#", 22))
cat("\ncritical value table\n")
print(TBL)
cat(strrep("#", 22))
cat("\n")
cat("T = ", n, "\n\n")
num <- which(n < c(25, 50, 100, 250, 500, Inf))[1]
num2 <- which(STAT2 < TBL[num, ])
if(length(num2) == 0) {
  cat("cannot reject H0: time series: cri is non-stationary\n")
  cat("STAT2: ", STAT2, "\n")
} else {
  num3 <- num2[1]
  cat("REJECT H0: time series is stationary\n")
  cat("with significant level", colnames(TBL)[num3], "\n")
  cat("STAT2: ", STAT2, "\n")
}


#犯罪率套件結果
if (k > 1) {
  diff_cri <- z1[, 2:k] ##差分後落後期資料
  res <- lm(diff_cri  ~ cri_lag + 1 + tt + d_cri_lag, drop.unused.levels = TRUE)
} else res <- lm(diff_cri ~ cri_lag + 1 + tt)
res_summary <- summary(res)
res_summary


####破獲率####
x2 <- x[, 2]
adf <- function(x2, k=NULL) {   #處理進行ADF的序列資料和滯後期數K
  
  if(is.null(k)) {
    k <- trunc((length(x2) - 1)^(1/3))
    k <- k + 1
  }
  
  if(k < 0) {
    k <- trunc((length(x2) - 1)^(1/3))
    k <- k + 1
  }}
adf_x2 <- diff(x2)
n <- length(adf_x2)
z2 <- embed(adf_x2, k)  #轉換序列資料w/落後期數
diff_dec <- z2[, 1] ##被解釋變數
dec_lag <- x2[k:n] #被解釋變數差分前k:n 的資料
tt <- k:n #時間
if(k > 1){
  d_dec_lag <- z2[, 2:k] ##差分後落後期
  adf_x2 <- cbind( dec_lag, "inter"=1, tt, d_dec_lag )
  adf_be2 <- solve(t(adf_x2) %*% adf_x2) %*% t(adf_x2) %*% diff_dec
} else {adf_x2 <- cbind(dec_lag, "inter"=1, tt)
adf_be2 <- solve(t(adf_x2) %*% adf_x2) %*% t(adf_x2) %*% diff_dec
}
#殘差
adf_ehdec <- diff_dec - adf_x2 %*% adf_be2 #same
edec_variance <- sum(adf_ehdec^2) / (length(tt) - ncol(adf_x2))
XTX_inverse <- solve(t(adf_x2) %*% adf_x2)
#xt1標準誤
element_11 <- XTX_inverse[1, 1]
dec_lag_se <- sqrt(edec_variance * element_11)
dec_lag_se
# 印出結果
cat("第", 1, "個係數的標準誤差:", dec_lag_se, "\n")
STAT3 <- adf_be2[1, 1]/dec_lag_se
TBL <- cbind(
  "1%" = c(-4.38, -4.15, -4.04, -3.99, -3.98, -3.96),
  "5%" = c(-3.60, -3.50, -3.45, -3.43, -3.42, -3.41))
rownames(TBL) <- paste0("T = ", c(25, 50, 100, 250, 500, Inf))
cat(strrep("#", 22))
cat("\ncritical value table\n")
print(TBL)
cat(strrep("#", 22))
cat("\n")
cat("T = ", n, "\n\n")
num <- which(n < c(25, 50, 100, 250, 500, Inf))[1]
num2 <- which(STAT3 < TBL[num, ])
if(length(num2) == 0) {
  cat("cannot reject H0: time series dec is non-stationary\n")
  cat("STAT3: ", STAT3, "\n")
} else {
  num3 <- num2[1]
  cat("REJECT H0: time series is stationary\n")
  cat("with significant level", colnames(TBL)[num3], "\n")
  cat("STAT3: ", STAT3, "\n")
}

##破獲率套件
if (k > 1) {
  d_dec_lag <- z2[, 2:k] ##差分後落後期資料
  res <- lm(diff_dec  ~ dec_lag + 1 + tt + d_dec_lag, drop.unused.levels = TRUE)
} else res <- lm(diff_dec ~ dec_lag + 1 + tt)
res_summary <- summary(res)
res_summary

##ADF統計量表格 #xt11:犯罪率/xt12:破獲率/xt1:失業率
res_df <- cbind("ADF stat"=c(STAT2, STAT3, STAT))
res_df 
#與ADFtest套件相同


adf_data <- data.frame(
  coef=c(adf_be1[1,1], adf_be2[1,1], adf_be[1,1]),
  標準誤=c(cri_lag_se, dec_lag_se, unemlag_se),
  ADF統計=c(STAT2, STAT3, STAT),
  ADF套件=c(adf.test(ly)$statistic, adf.test(x2)$statistic, adf.test(x3)$statistic)
)
adf_data
