# がんのGWAS
## GWAS による乳がん感受性多型の検出
「がんゲノムデータ解析（ 
清水　厚志, 坊農　秀雅編；メディカル・サイエンス・インターナショナル刊）」の第６章「がんのGWAS」に掲載されたコマンド、スクリプトのみ掲載しています。テキストは本文を御覧ください。

## 6.3 解析環境構築
### インストールするもの
- (M1 Mac の場合)Rosetta2
- plink 1.90 beta 6.24 64‒bit ●plink 2.00 alpha v2.00a2.3 64‒bit ●R version 4.1.0
- R studio version 1.4.1717
- wget 1.21.1
- md5sum 0.9.5

### インストール手順
#### 作業ディレクトリの準備
```bash
# bash
# 1
cd ~/Desktop
mkdir WD

# 2
cd WD
mkdir BIN
mkdir QC
mkdir GWAS
mkdir DATASET

# 3
ls
```

#### (M1 Mac の場合)Rosetta2
```bash
# bash
# 1
softwareupdate --install-rosetta
# 2-4
# (本文参照)
```

#### PLINK 1.90 beta
(本文参照)


#### PLINK 2.00 alpha
(本文参照)


#### R studio
(本文参照)

#### wget
```bash
# bash
# 1
brew install wget

# 2
wget --version
```

#### md5sum
```bash
# bash
# 1
brew install md5sha1sum

# 2
md5sum --version
```

## 6.4 データの入手
### 使用するデータセット
```bash
# bash
# 1
cd ~/Desktop/WD/DATASET
wget -N -r --no-parent -nd \
  http://www.medicalgenome.info/data/CancDAT/chapter6_GWAS/
# 1行目;DATASETフォルダへの移動
# 2行目;データセットのダウンロード

# 2
cd ~/Desktop/WD/DATASET
md5sum -c md5sum.txt
# 1行目;DATASETフォルダへの移動
# 2行目;ハッシュの計算と比較

# 3 (option)
cd ~/Desktop/WD/DATASET
wget -nd http://www.medicalgenome.info/data/CancDAT/chapter6_GWAS/(破損ファイル名)

# 1行目;DATASETフォルダへの移動 
# 2行目;データセットのダウンロード

# 4
# (本文参照)
```

## 6.5 サンプル QC
### データセットの確認
```bash
# bash
# 1
cd ~/Desktop/WD/DATASET
head 1000GP_Phase3_chr1.sample
gunzip -c 1000GP_Phase3_chr1.gen.gz | head | cut -d " " -f 1-10
# 1行目;DATASETディレクトリへの移動
# 2行目;.sampleファイルの確認。headコマンドで，上位10行のみを書き出す。
# 3行目;gen.gzファイルの確認。gz形式で圧縮されているので，gunzipでそれをいったん展開している。続いて，head で上位 10 行だけを取り出す。横にも長いファイルなので，cut コマンドで，1~10 列目までを切り出している。cut -dの後のダブルクオーテーション(")の間には半角スペースを入れる。.gen ファイルが半角スペースで区切られていることを cut コマンドに教えている。

# 2
cd ~/Desktop/WD/DATASET
head 1000GP_Phase3_chr1.sample

# 3
gunzip -c 1000GP_Phase3_chr1.gen.gz | head | cut -d " " -f 1-10
```

### ファイル形式の変換 1
```bash
# bash
# 1 
cd ~/Desktop/WD/QC
../BIN/plink \
  --gen ../DATASET/1000GP_Phase3_chr1.gen.gz \
  --sample ../DATASET/1000GP_Phase3_chr1.sample \
  --oxford-single-chr 1 \
  --make-bed \
  --keep-allele-order \
  --out 1000GP_Phase3_chr1
# 1行目;QCディレクトリへの移動
# 2行目;BINディレクトリのplinkを指定する。../は現在地より1つ上の階層のディレクトリを意味する。
# 3行目;--gen; .gen.gzファイルを指定する。
# 4行目;--sample; .gen.gzに対応するsampleファイルを指定する。
# 5行目;--oxford-single-chr; このファイルの染色体番号を指定する。 
# 6行目;--make-bed; バイナリplink形式で出力させる。

# 2
# (本文参照)

# 3
head 1000GP_Phase3_chr1.fam
head 1000GP_Phase3_chr1.bim

# 4
cd ~/Desktop/WD/QC
for chr in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22;do
../BIN/plink \
  --gen ../DATASET/1000GP_Phase3_chr$chr.gen.gz \
  --sample ../DATASET/1000GP_Phase3_chr$chr.sample \
  --oxford-single-chr $chr \
  --make-bed \
  --keep-allele-order \
  --out 1000GP_Phase3_chr$chr
done
# 2行目;forコマンド。末尾のdoから9行目のdoneまでの間に記載されたコマンドを繰り返す。このとき，指定した変数(ここでは chr)に，in 以降にリストアップされた文字や数字(ここでは 2~22 の数字) を順番に代入して繰り返す。2 から 22 の数字を記入するのが面倒な場合は seq コマンドを使って $(seq 2 22)と指定しても同じ結果が得られる。

# 5
# wc -lでファイルの行数を表示する。-lはハイフンと小文字のL。
wc -l 1000GP_Phase3_chr{1..22}.fam
```

### 低call rateによる除外
```bash
# bash
# 1
cd ~/Desktop/WD/QC
for chr in $(seq 1 22);do
../BIN/plink \
  --bfile 1000GP_Phase3_chr$chr \
  --mind 0.05 \
  --make-bed \
  --keep-allele-order \
  --out SQC1_chr$chr
done
# 2行目;forコマンド。変数chrの内容を1~22に変えながら繰り返す。
# 4行目;--bfile; 先ほど変換して作成された.bed/.bim/.famのファイル名(拡張子を除く部分)を指定する。
# 5行目;--mind;missing rateを指定する。

# 2
wc -l SQC1_chr{1..22}.fam
```

### pruning(刈りこみ)
```bash
# bash
# 1
cd ~/Desktop/WD/QC
for chr in $(seq 1 22);do 
../BIN/plink \
  --bfile SQC1_chr$chr \
  --geno 0.01 \
  --hwe 0.05 \
  --maf 0.05 \
  --make-bed \
  --keep-allele-order \
  --out SQC2_chr$chr
done
# 5行目;--geno; missing rateが指定した値より高いバリアントを除外する。
# 6行目;--hwe; ハーディ・ワインベルグ平衡からはずれたバリアントを除外する。 
# 7行目;--maf; 指定したマイナーアレル頻度より低いバリアントを除外する。

# 2
wc -l SQC2_chr{1..22}.bim

# 3
for chr in $(seq 1 22);do
../BIN/plink \
  --bfile SQC2_chr$chr \
  --indep-pairwise 1500 150 0.03 \
  --out SQC2_chr$chr
done
# 3行目;--indep-pairwise; pruning条件を指定する。1,500 kbの範囲のウィンドウを，150バリアントごとにずらして LD を計測し，その指標(r2)が>0.03 であるものを除去する。

# 4
for chr in $(seq 1 22);do
../BIN/plink \
  --bfile SQC2_chr$chr \
  --extract SQC2_chr$chr.prune.in \
  --make-bed \
  --keep-allele-order \
  --out SQC3_chr$chr
done
# 3行目;--extract; 指定されたファイル中のバリアントのみ残す。

# 5
wc -l SQC3_chr{1..22}.bim

# 6
for chr in $(seq 1 22);do
  echo SQC3_chr$chr
done > merge-list.txt
cat merge-list.txt

../BIN/plink \
  --merge-list merge-list.txt \
  --keep-allele-order \
  --out SQC3_chrALL
# 1~3行目;マージするファイルのリストを作る。
# 4行目;リストを表示して確認する。
# 7行目;--merge-list; マージするファイル名のリストを指定する。

# 7
# (本文参照)
```

### 血縁者の除去
```bash
# bash
# 1 
cd ~/Desktop/WD/QC
../BIN/plink \
  --bfile SQC3_chrALL \
  --genome \
  --min 0.1875 \
  --out SQC3_chrALL
# 4行目;--genome; PI_HATやIBDを算出する。
# 5行目;--min; 出力するPI_HATの最小値を設定する。

# 2
head SQC3_chrALL.genome

# 3
awk '{print $1, $2}' SQC3_chrALL.genome > pi-hat.txt

../BIN/plink \
  --bfile SQC3_chrALL \
  --remove pi-hat.txt \
  --make-bed \
  --keep-allele-order \
  --out SQC4_chrALL
# 1行目;血縁のあるペアのうち，片側のIDのみ抜き出す。 
# 5行目;--remove; ファイル中のIDを除く。

# 4
wc -l SQC4_chrALL.fam
```

### 集団構造の検出
```bash
# bash
# 1
cd ~/Desktop/WD/QC
awk '{ print $2 }' SQC4_chrALL.bim > extractSNP.txt

../BIN/plink \
  --bfile ../DATASET/1KG.ALL.EAS \
  --extract extractSNP.txt \
  --make-bed \
  --keep-allele-order \
  --out 1KG.ALL.EAS.ex
# 2行目;bimファイルから，pruningしたデータのバリアント名の部分のみ(bimファイルの2列目)を抜き出して，extractSNP.txt ファイルに保存する。
# 6行目;--extract; 指定したファイル中のバリアント名に合致するバリアントのみをデータセットから抜き出す。

# 2
../BIN/plink \
  --bfile SQC4_chrALL \
  --bmerge 1KG.ALL.EAS.ex \
  --keep-allele-order \
  --out SQC4_chrALL_ref
# 3行目;--bmerge; 指定したデータセット(.bed/.bim/.famの形式である必要がある)と結合する。

# 3
cd ~/Desktop/WD/QC
../BIN/plink \
  --bfile SQC4_chrALL_ref \
  --geno 0.01 \
  --make-bed \
  --out SQC4_chrALL_ref2

# 4
../BIN/plink \
  --bfile SQC4_chrALL_ref2 \
  --pca \
  --out SQC4_chrALL_ref2
# 3行目;--pca; PCAを行う。

# 5
cd ~/Desktop/WD/DATASET
wget \
  http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
```

```R
# R
# 6
install.packages("ggplot2")


# 7
# ggplot2 を読み込む
library(ggplot2)

# 作業ディレクトリを設定する
setwd("~/Desktop/WD")

# 先ほどのPCA結果を読み込む
d <- read.table("./QC/SQC4_chrALL_ref2.eigenvec")

# 列名を変更
names(d) <- c("FID","IID", paste0("PC",seq(1,20)))

# 1000 genomeデータのどの検体がどの民族に属しているかを示したデータを読み込む
i <- read.table("./DATASET/integrated_call_samples_v3.20130502.ALL.panel", header=T)

# PCA結果とマージするため，一列目の列名をIIDにして，PCA結果と列名をそろえる。
names(i)[1] <- c("IID")

# PCA結果と民族についての情報をマージする。IIDの列でマージされる。
d2 <- merge(d, i, all=T)

# 今回のデータには民族名がないため，仮に"This_study"と表記させる。
d2[is.na(d2$pop),]$pop <- "This_study"

# PC1にデータがない検体を除外
d3 <- subset(d2, !is.na(PC1))

# ggplotで作図
ggplot(data=d3, aes(x=PC1, y=PC2, color=pop)) +
 geom_point()
ggsave("pca.png")


# 8
# PC1，PC2で範囲指定し，サブセットを作成する。
d4 <- subset(d3, d3$PC1 > -0.05 & d3$PC1 < 0.01 & d3$PC2 > -0.05 & d3$PC2 < 0.05)

# もう1度作図して確認する(図6.8)。
ggplot(data=d4, aes(x=PC1, y=PC2, color=pop)) + 
  geom_point()
ggsave("pca2.png")


# 9
# IDを含む1,2列目のみ取り出す。
idlist <- d4[,2:1]

# リストを書き出す。row.names=Fを指定することで，行番号が記載されるのを防ぐ。また，quote=Fと指定して，各セルのデータがダブルクォーテーションで囲われるのを防ぐ。
write.table(idlist, "./QC/sampleQC.txt", row.names=F, quote=F)


# 10
# (本文参照)
```

### サンプル QC 結果の反映
```bash
# bash
# 1
cd ~/Desktop/WD/QC
for chr in $(seq 1 22);do 
  ../BIN/plink2 \
    --gen ../DATASET/1000GP_Phase3_chr$chr.gen.gz ref-first \
    --sample ../DATASET/1000GP_Phase3_chr$chr.sample \
    --missing-code -9,NA \
    --oxford-single-chr $chr \
    --make-pgen \
    --out 1000GP_Phase3_chr$chr
done
# 3行目;PLINK 2を使用する。
# 6行目;--missing-code; ファイル中の欠損値の記述方法を指定する。 
# 8行目;--make-pgen; PLINK 2のファイル形式で保存する。

# 2
# 1, 2行目;wc -lでファイルの行数を表示する。-lはハイフンと小文字のL
wc -l 1000GP_Phase3_chr{1..22}.psam
wc -l 1000GP_Phase3_chr{1..22}.pvar

# 3
cd ~/Desktop/WD/QC
for chr in $(seq 1 22);do
  ../BIN/plink2 \
    --pfile 1000GP_Phase3_chr$chr \
    --keep sampleQC.txt \
    --make-pgen \
    --out 1000GP_Phase3_sQC_chr$chr
done
# 4行目;--pfile; PLINK 2形式のファイルを指定する。
# 5行目;--keep; 指定したファイル内のIDリスト中にある検体のみ保持させる。

# 4
wc -l 1000GP_Phase3_sQC_chr{1..22}.psam
# 1行目;wc -lでファイルの行数を表示する。
```

## 6.6 形質の値および共変量
```R
# R
# 1
# 作業ディレクトリを設定する。
setwd("~/Desktop/WD")

# フェノタイプデータを読み込む。
d <- read.table("./DATASET/1000GP_Phase3_n2000_tags.pheno")

# 必要な部分だけ取り出す。
d <- d[,1:4]

# 列名をつける。
names(d) <- c("FID", "IID", "missing", "PHENO")

# 自分で計算したPCA結果を読み込む。
pc <- read.table("./QC/SQC4_chrALL_ref2.eigenvec")

# 列名をつける。
names(pc) <- c("FID","IID", paste0("PC",seq(1,20)))

# マージする(共通の列名でマージされる)。
d2 <- merge(d, pc)

# 保存する。
write.table(d2, "~/Desktop/WD/DATASET/Cancer.pheno", row.names = F, quote = F)

# 2
# (本文参照)
```

```bash
# bash
# 3
cd ~/Desktop/WD/DATASET

# データの総セル数
wc -w ~/Desktop/WD/DATASET/Cancer.pheno

# データの行数
wc -l ~/Desktop/WD/DATASET/Cancer.pheno
```

## 6.7 GWAS の実施
```bash
# bash
# 1
cd ~/Desktop/WD/GWAS
for chr in $(seq 1 22);do
../BIN/plink2 \
  --ci 0.95 \
  --pfile ../QC/1000GP_Phase3_sQC_chr$chr \
  --glm omit-ref hide-covar cols=chrom,pos,ref,alt,a1freq,test,nobs,beta,se,ci,p \
  --pheno ../DATASET/Cancer.pheno \
  --pheno-name PHENO \
  --covar ../DATASET/Cancer.pheno \
  --covar-name PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 \
  --covar-variance-standardize PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 \
  --maf 0.01 --hwe 1e-6 midp --geno 0.01 \
  --vif 1000 \
  --memory 4000 \
  --threads 4 \
  --out Cancer.gwas.chr$chr
done
# 4行目;--ci; 指定した信頼区間を出力する
# 6行目;--glm; 線形またはロジスティック回帰分析を行う。omit-refにより，ALTを解析するように指定し，hide-covar で調整変数の結果を出さないようにしている。cols 以下は出力のフォーマット指定。
# 7行目;--pheno;形質データのファイル名を指定。
# 8行目;--pheno-name;形質データ中の形質データの列名を指定する。
# 9行目;--covar;調整変数データのファイル名を指定する。今回は形質データとまとめられているので，同一である。
# 10行目;--covar-name; 調整変数データ中の使用する調整変数の列名を指定する。
# 11行目;--covar-variance-standardize; 調整変数を標準化して解析する。
# 13行目;--vif; 多重共線性チェックを行う際，結果を返さない閾値を指定する。今回は上げているが，通常は 50 に設定されている。
# 14行目;--memory; 使用可能なメモリサイズ。一般に大きいほど速くなる。単位はMb。
# 15行目;--threads; 使用可能なCPUコア数。一般に大きいほど速くなる。


# 2
wc -l Cancer.gwas.chr{1..22}.PHENO.glm.logistic

# 3
cd ~/Desktop/WD/GWAS
head -n 1 Cancer.gwas.chr1.PHENO.glm.logistic > header.txt
cat Cancer.gwas.chr{1..22}.PHENO.glm.logistic | sed -e '/#/d' > body.txt
cat header.txt body.txt > Cancer.gwas.chrALL.PHENO.glm.logistic
# 2行目;実行結果のヘッダー行のみ別に保存しておく。
# 3行目;各染色体の結果を結合し，#で始まるヘッダー行を除く。 
# 4行目;2，3行目で処理した結果を結合する。


# 4
head Cancer.gwas.chrALL.PHENO.glm.logistic

# 5
wc -l Cancer.gwas.chrALL.PHENO.glm.logistic
```

## 6.8 作図
### qq プロット/マンハッタンプロット
```R
# R
# 1
install.packages("qqman")
install.packages("data.table")

# 2
# ライブラリーの読み込み 
library(qqman)
library(data.table)

# 作業ディレクトリを設定する
setwd("~/Desktop/WD")

# GWAS結果の読み込み
d <- fread("./GWAS/Cancer.gwas.chrALL.PHENO.glm.logistic")

# 作図に必要な部分(染色体番号，バリアントの位置，バリアント名，リファレンスアレル，オルタナ ティブアレル，頻度，β，標準誤差，P 値)のみ抽出
d2 <- d[,c(1,2,3,4,5,7,10,11,14)]

# qqmanが各列を認識しやすいように，列名を変更する
names(d2) <- c("CHR","BP","SNP","REF","ALT","FREQ","BETA","SE","P")

# 3
manhattan(d2)

# 4
qq(d2$P)

# 5
# 観測されたP値の中央値のχ2(右側)と，本来期待されるのP値の中央値(0.5)のχ2(右側)の比をとり，λとする
lambda <- qchisq(median(d2$P), df=1, lower.tail=F)/qchisq(0.5,df=1, lower.tail=F)
lambda
```

### LocusZoom
```R
# R
# 1
# min(d2$P)でP値の最小値を算出し，その値をもつバリアントのBP，CHRの値を得る。 
# 6番染色体のピークのBP(位置)の情報を得る。
peakCHR <- 6
d3 <- subset(d2, d2$CHR == peakCHR)
peak <- subset(d3, d3$P == min(d3$P))
peakBP <- peak[,"BP"]

# 最小P値のBPから +/- 1 M(1×106)bp離れた場所のBPを計算する。as.integer()は計算結果を整数に変換するための処置。
upperBP <- as.integer(peakBP + 1e6)
lowerBP <- as.integer(peakBP - 1e6)

# +/- 1 Mbpの範囲に入るエリアをGWAS結果から抽出する。
d4 <- subset(d3, d3$CHR == peakCHR & d3$BP < upperBP & d3$BP > lowerBP)

# データを保存する。
write.table(d4, "./GWAS/hit.6.txt", row.names = F, quote = F, sep = "\t")

# 続いて16番染色体のピーク。
peakCHR <- 16
d3 <- subset(d2, d2$CHR == peakCHR)
peak <- subset(d3, d3$P == min(d3$P)) 
peakBP <- peak[,"BP"]

# 最小P値のBPから +/- 1 Mbp離れた場所のBPを計算する。
upperBP <- as.integer(peakBP + 1e6)
lowerBP <- as.integer(peakBP - 1e6)

# +/- 1M bpの範囲に入るエリアをGWAS結果から抽出する。
d4 <- subset(d3, d3$CHR == peakCHR & d3$BP < upperBP & d3$BP > lowerBP)

# データを保存する。
write.table(d4, "./GWAS/hit.16.txt", row.names = F, quote = F, sep = "\t")

# 2-8
# (本文参照)
```

## 6.9 下流解析
- GWAS catalog(https://www.ebi.ac.uk/gwas/)
- GTEx(https://gtexportal.org/home/)
- gnomAD(https://gnomad.broadinstitute.org/)
- dbSNP(https://www.ncbi.nlm.nih.gov/snp/)
- jMorp(https://jmorp.megabank.tohoku.ac.jp/202001/)
