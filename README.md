# 楕円フーリエ解析解析メモ

# 解析

## 手順

### 輪郭座標の抽出

\"**contour_all.py**\"で葉のスキャン画像から、葉を1枚ずつ抽出し、輪郭のxy座標をcsvファイルに出力する。他にも処理の際に生じる画像ファイルや葉のデータが以下のように出力される。これはフォルダ単位での処理となる。

| 出力ファイル                                                             | 場所とファイル名                                                   |
|--------------------------------------------------------------------------|--------------------------------------------------------------------|
| 輪郭のxy座標                                                             | ./contour/[画像フォルダ名]/contour\_[画像ファイル名]\_[葉のID].csv |
| 元画像の二値化画像に、葉ごとにIDと重心、長軸と短軸の線を引いたラベル画像 | ./label/[画像フォルダ名]/[画像ファイル名]\_label.jpg               |
| それぞれの葉の面積、周囲長、真円度、外接長方形の幅と長さのデータシート   | ./label/[画像フォルダ名]/[画像ファイル名]\_leaf_info.csv           |
| 葉ごとのトリミング画像 (縦横に任意の余白あり)                            | ./label/[画像フォルダ名]/[画像ファイル名]\_[葉のID].jpg            |

### 正規化楕円フーリエ係数の算出

"elliptic_fourier_analysis.R"での解析に入る。まずは正規化楕円フーリエ係数 (Normalized Fourier coefficients)を算出する。

### Harmonic Fourier powerの算出

計算式を書く

### 主成分分析 (principal component analysis; PCA)

# 問題点

横長の葉っぱや円に近い葉では、長軸が葉脈の方向と一致しないため、他の葉と単純に比較することができない。長軸がきれいに葉脈と平行になっているときは、楕円フーリエ解析の計算段階で90度回転する処理を入れれば対応可能と思われる。

二値化をした際に、端の部分が白くなることがあり、これが検知されてしまうことがある。対処法としては、端の部分をトリミングすれば良いと考えられる