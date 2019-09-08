# 必要なパッケージの読み込み
library(ks)
library(tidyverse)

# eval_component: 
# コンポーネントの中心の座標について，
# ヌルモデルとの比較によるmeaningful clusterか否か判定する関数
# 
# 入力：
#    x, # 因子スコア (回答者×因子)
#    c.center, 評価点　(コンポーネントの中心の座標; コンポーネント×因子のmatrix)
#    n.shuffle = 1000 シャッフルする回数
#
# 出力: 
#    $d.original # 評価点での元のデータの密度
#    $d.null     # 評価点でのヌルモデルでのデータの密度
#    $p.value    # p値
#    $enrichment # = d.original/d.null
eval_component <- function(x, c.center, n.shuffle = 1000) {
  
  xdim <- ncol(x) # 因子の数 (次元数)
  N <- nrow(x) # 回答者の数
  n.component <- nrow(c.center) # コンポ――ネントの数
  
  # returnする変数
  p.value <- numeric(n.component)
  enrichment <- numeric(n.component)
  
  # 元データのカーネル密度推定
  cat("Bandwidth selection...\n")
  # データ数が多いときに全データを用いると時間がかかるので，
  # 上限を1000までにして用いる。
  Hpi <- ks::Hpi.diag(x = x[1:min(N,1000),], 
                      nstage = 2)
  cat("kernel dinsity estimation for original data...\n")
  k <- ks::kde(x = x, H = Hpi,
               eval.points = c.center,
               binned = FALSE) 
  # 4次元以上の場合はbinned estimationができないということで
  # binned = FALSEとしないとerrorになる。
  density.original <- k$estimate
  
  # シャッフルデータの密度格納用行列
  density.shuffle <- matrix(0, n.shuffle, nrow(c.center))
  
  cat("kernel dinsity estimation for shuffled data...\n")
  for (idxs in 1:n.shuffle) {
    if (idxs %% 10 == 1) 
      cat("=")
    
    # シャッフルする
    sc_shuffeled <- apply(x, MARGIN = 2, sample)
    
    # シャッフルデータの密度推定
    k <- ks::kde(x = sc_shuffeled, H = Hpi, 
                 eval.points = c.center,
                 binned = FALSE)
    
    density.shuffle[idxs,] <- k$estimate
  }
  cat("\n")
  
  # ヌルモデルの密度平均
  density.null <- colMeans(density.shuffle)
  
  for (idx in 1:n.component) {
    
    # p値
    p.value[idx] <- sum(density.original[idx] < density.shuffle[,idx]) / n.shuffle
    
    # enrichment
    enrichment[idx] <- density.original[idx] /
      density.null[idx]
    
    d.null.mean <- mean(density.shuffle[,idx])
  }
  
  list(d.original = density.original,
       d.null = density.null, 
       p.value = p.value,
       enrichment = enrichment )
}


# plot_meaningful_cluster: meaningful clusterをプロットする関数
# 
# 入力：
#    res_ec: 関数eval_componentの出力
#    c.center: 評価点　(コンポーネントの中心の座標; コンポーネント×因子のmatrix)
#    p.threshold = 0.01, # 有意水準
#    enrichment.threshold = 1.25 # enrichmentの閾値
# 出力: 
#    $d.original # 評価点での元のデータの密度
#    $d.null     # 評価点でのヌルモデルでのデータの密度
#    $p.value    # p値
#    $enrichment # = d.original/d.null
plot_meaningful_cluster <- function(res_ec, c.center, 
                                    p.threshold = 0.01, 
                                    enrichment.threshold = 1.25) {
  
  n.component <- nrow(c.center) # コンポ―ネントの数
  
  df_cluster <- data.frame(c.idx = 1:n.component, #コンポーネントのID
                           res_ec) 
  
  df_c <- df_cluster %>% arrange(desc(enrichment)) %>%
    filter(enrichment > enrichment.threshold, 
           p.value < p.threshold) 
  
  ymax <- max(component.centers) + abs(max(component.centers)) * 0.2
  ymin <- min(component.centers) - abs(min(component.centers)) * 0.2
  
  if (nrow(df_c) > 0){
    for (idx in 1:nrow(df_c)) {
      barplot(height = as.matrix(component.centers[df_c$c.idx[idx],]), 
              ylim = c(ymin,ymax),
              main = sprintf("component %d,\nenrichment = %.2f\n p = %.2f", 
                             df_c$c.idx[idx], 
                             df_c$enrichment[idx],
                             df_c$p.value[idx]
              ), 
              xlab = "") 
      abline(h = 0, lty = 2)
    }
  } else {
    print("There was no meaningful cluster.\n")
  }
}
