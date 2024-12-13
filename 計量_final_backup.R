T <- read.csv("month.csv", header=TRUE, check.names=FALSE) #2010.1~2023.8月資料
options(max.print = 10000)
T
y <- ("cri"= T[, 2])
y <- as.matrix(y)
y
x <- T[, 3:4]
x <- as.matrix(x)

#敘述統計
 ##平均
mean_y <- sum("cri"=log(y))/length(log(y))
mean_x2 <- sum("dec"=x[, 2])/148
mean_x3 <- sum("unem"=x[, 3])/148
 ##標準差
sd(log(y))
sd(x[,2])
sd(x[, 3])
summ <-cbind("nobs"= c(length(y),length(x[,2]),length(x[,3])),
            "mean"= c(mean_y, mean_x2, mean_x3), 
            "sd"= c(sd(log(y)), sd(x[, 2]), sd(x[, 3])),
            "min"= c(min(log(y)), min(x[, 2]), min(x[, 3])),
            "max"= c(max(log(y)), max(x[, 2]),  max(x[, 3])))

#放入截距
x <- cbind(rep(1), x)
x
beta <- solve(t(x) %*% x) %*% t(x) %*% y
beta
lm_f <- lm(y~x)
summary(lm_f) 


##看時間趨勢圖
plot(T[, 1], log(y), type="l", xlab="Time", ylab="犯罪率(十萬/件")
plot(T[, 1], T[, 3], type="l", xlab="Time", ylab="破獲率")
plot(T[, 1], T[, 4], type="l", xlab="Time", ylab="失業率")

##serial correlation
#手刻ACF
#犯罪率#犯罪率以取對數來分析
res <- acf(log(y))
n <- length(y)
y0 <- log(y)-mean_y
nlags <- 10
res2 <- sapply(1:nlags, \(i) {  #對每個滯後期數迭代
  a <- y0[seq_len(n-i)]  #1~(n-i)的資料(xi-mean)的值
  b <- y0[-seq_len(i)]  #剩下的值
  sum(a * b) / sum(y0 * y0) #算自相關係數
})
RES <- cbind("by_hand" = c(1, res2), 
             "acf" = res$acf[,,1][1:(nlags+1)])

rownames(RES) <- paste0("lags.", 0:nlags)
RES

acf_values1 <- acf(log(y))$acf
acf_values1

#破獲率
res_dec2 <- acf(x[, 2])
n <- length(y)
x2 <- x[, 2]-mean_x2
nlags <- 21
res_dec <- sapply(1:nlags, \(i) { 
  a2 <- x2[seq_len(n-i)]
  b2 <- x2[-seq_len(i)]
  sum(a2 * b2) / sum(x2 * x2) 
})
RES <- cbind("by_hand" = c(1, res_dec), 
             "acf" = res_dec2$acf[,,1][1:(nlags+1)])

rownames(RES) <- paste0("lags.", 0:nlags)
RES

#失業率
res_unem2 <- acf(x[, 3])
n <- length(y)
x3 <- x[, 3]-mean_x3
nlags <- 17
res_unem <- sapply(1:nlags, \(i) {
  a3 <- x3[seq_len(n-i)]
  b3 <- x3[-seq_len(i)]
  sum(a3 * b3) / sum(x3 * x3) 
})
RES <- cbind("by_hand" = c(1, res_unem), 
             "acf" = res_unem2$acf[,,1][1:(nlags+1)])


rownames(RES) <- paste0("lags.", 0:nlags)
RES
##

#套件對照
acf(log(y)) 
acf(x[, 3]) 
acf(x[, 2])


##平穩性 
library(tseries)
adf.test(y)
adf.test(x[,2])
adf.test(x[, 3])


##差分
de_y <- diff(log(y))                  #犯罪率差分
time_axis_diff <- seq_along(de_y)
plot(de_y, type="l", xlab="Time", ylab="crime")
de_x <- diff(x[, 2:3])                #失業率和破獲率差分
de_x<- cbind("Inte" = 1, de_x)
beh <- solve(t(de_x) %*% de_x) %*% t(de_x) %*% de_y
lm_d <- lm(de_y ~ 0 + de_x)  ##這邊的R squared要用centered

lm_m <- lm(de_y ~ de_x[,2]+de_x[, 3])
##差分後檢定異質性:BP test
bptest(lm_d) 
#手刻殘差
ehat <- de_y - de_x %*% beh
eh2 <- ehat^2 #被解釋變數
be3 <- solve(t(de_x) %*% de_x) %*% t(de_x) %*% eh2
be3
#做R squared
yhat <- de_x %*% be3
R2 <- drop(crossprod(yhat)/crossprod(eh2))
R2
bp1 <- nrow(de_x)*R2
bp1
p_value <- 1 - pchisq(bp1, df = 2)
p_value #有異質性
#套件對照
aux_model2 <- lm(eh_d^2 ~ 0+de_x)
summary(aux_model2)
bp_statistic2 <- nobs(lm_d) * summary(aux_model2)$r.squared ##LM
bp_statistic2
p_value2 <- 1 - pchisq(bp_statistic2, df = 2)
p_value2

##跟套件一樣(無截距)
R22<- t(yhat-mean(yhat))%*%(yhat-mean(yhat))/t(eh2-mean(eh2))%*%(eh2-mean(eh2)) #無截距
bp <- nrow(de_x)*R22
p_value22 <- 1 - pchisq(bp, df = 2)
p_value22
a_model <- lm(eh_d^2 ~ de_x[, 2]+de_x[, 3])
summary(a_model)
bp_ <- nobs(lm_d) * summary(a_model)$r.squared
bp_
bptest(lm_d) ##無截距
v3 <- diag(vcov(lm_d))
v3
my_data <- data.frame(
       Item = c("By hand", "By package"),
       BPstat = c(bp1, bp_statistic2),
       Pvalue = c(p_value, p_value2))
##

#under異質性#用white estimator
eh2 <- ehat^2
box <- 0
i <- 1
for(i in 1:nrow(de_x)){
  box　<- box+eh2[i]*de_x[i, ] %*% t(de_x[i, ])#white估計中的sum(殘差平方*X*X')
}
XtX <- t(de_x) %*% de_x
HC <- solve(XtX) %*% box %*% solve(XtX)
v <- diag(solve(XtX) %*% box %*% solve(XtX)) #white估計後的beta變異數
print(v)
#套件對照
library("lmtest")
library(sandwich)
v2 <- diag(vcovHC(lm_d, type="HC0")) ##vcov:white估計用來計算共變異數矩陣
v2
cbind("by_hand"= v, "package"= v2)
coeftest(lm_d, vcov = vcovHC(lm_d, type="HC0"))
##

#
RES <- cbind("Estimate"=c(beh),
             "std.Error"=sqrt(v),
             "t value"=c(beh)/sqrt(v),
             "Pr(>|t|)"= 2-2*pt(abs(c(beh)/sqrt(v)), df=nrow(de_x)))
round(RES, 7) ##模型結論


#wald test
R <- diag(ncol(de_x))
R <- R[-1, ]
r <- cbind(rep(0, ncol(de_x)-1))
q <- nrow(R)
N <- nrow(de_x)
Wald.test.HC <- t(R %*% beh -r) %*% solve(R %*% HC %*% t(R)) %*% (R %*% beh -r)/q
Wald.test.HC
1-pf(Wald.test.HC, q, N-q)
