# Scalar Wave Propagation Library
SWPLはスカラー波動光学による光波伝搬を計算するC++クラスライブラリです。

### 実装済みの伝搬手法一覧
- Rayleigh-Sommerfeld Integral (Eigenで実装 :高速)
- Rayleigh-Sommerfeld Integral (4重のforで実装 :低速)

### 依存ソフトウェア
- OpenMP
- Eigen

### ビルド
- VisualStudio
```
cl /EHsc /arch:AVX2 /O2 /Oi /openmp /I<Eigenへのインクルードパス> /source-charset:utf-8 <ソースファイル>
```

- g++
```
g++ -std=c++11 -O2 -mavx2 -fopenmp -I<Eigenへのインクルードパス>
```

- Intel Compiler
```
icpc -std=c++11 -O2 -xcore-avx2 -openmp -I<Eigenへのインクルードパス> <ソースファイル>
```

### テスト
test ディレクトリの中にあるスクリプト test.bat はテストを実行します。

```
test.bat <テスト名>
```

スクリプトは以下のテスト名を一つ引数に取り、それぞれのテストを実行します。

|引数名|テストの内容|
|---|---|
|performance|ライブラリに含まれる各伝搬計算手法ごとの実行時間を計測します。|
