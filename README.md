# コード
- Julialangによる2D compressible Navier-Stokes方程式のPISO法による解析コード
- コロケート格子に対応するため，Rhie-Chow補間を適用
- 解法の詳細はdoc内のPDFに記載
- 例題として衝撃波管の問題を解いてます。OpenFOAMの結果とほぼ同じ結果になるは確認しました。

# 実行
データを出力フォルダ(outfileで指定する部分）を作成してから，

```
include("Unsteady-compr-PISO_sod-colocation-relaxation-ud.jl")
para = Param()
@time main(para)
```
