### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 5be67619-2fa6-48e7-b8e4-492360b6c9c0
# Example というパッケージをロードする
using Example

# ╔═╡ bdb66039-1c23-4172-afef-ffab8452ddf6
begin
	# テストコードを記述するために必要なパッケージ
	using Test
end

# ╔═╡ d7f1cb15-6ab1-4fd6-8639-e20ded4f29ea
md"""
# Hello world
"""

# ╔═╡ 408ea404-a897-4574-97e0-5cab52b9d154
versioninfo()

# ╔═╡ 04c275d6-ddd3-11ee-0102-4ff69f77eaac
println("Hello")

# ╔═╡ aeed025d-ce2d-46ca-be64-3aea35e51932
md"""
`using Example` によって `hello` という関数が利用できる
"""

# ╔═╡ ee932197-1a89-4ada-b5eb-6e0514c2f47c
hello("YourName") # ここを変更して遊んでみよう

# ╔═╡ 4c482f7b-0f75-4b77-ad79-62ce700d5834
md"""
`Example.hello` のようにして `hello` が属するパッケージ名を明示しても良い
"""

# ╔═╡ f853d44f-cb09-4fd2-8afa-ea89821c829b
md"""
`using Example: hello` のようにすると `hello` 関数をどこから持ってきたのかを他の人にもわかるように実装できる．

Python ユーザのように名前空間を大事にしたい, 慎重に管理したい場合は下記のようにしてもよかろう:

```julia
using Example: Example # Example から Example のみをロードする
Example.hello # このような方法でしか `hello` を使用できない
```

この場合 `hello("Name")` のように直接使うことはできない． また `using Example: Example` は `Example` を二回用いてるので記述がやや冗長になる． この冗長さを回避するために Julia は `import Example` という書き方も許している:

```julia
import Example
Example.hello("Name")
```
"""

# ╔═╡ b5c5f82b-5cd6-4b71-8b9f-d4b9715d21b7
begin
	x = 2
	y = 3
end

# ╔═╡ a0f3c104-da2d-4540-ae62-84a9f2dee111
z = 5 # 上のセルで定義された x, y の変数を変更すると z も更新されることに注意

# ╔═╡ ff5b2387-b22b-423d-8b4c-e99a085f5293
"""
	Point

二次元空間の点を表すオブジェクト
"""
struct Point2D
	x::Float64
	y::Float64
end

# ╔═╡ 44641fc0-7c8c-46f6-be51-438be0ed742d
"""
	distance(p::Point2D, q::Poin2D)

点 `p`, `q` の間の距離を計算する
"""
distance(p::Point2D, q::Point2D) = √((p.x - q.x) ^ 2 + (p.y - q.y) ^ 2)

# ╔═╡ e2b487c7-7d56-4cb1-b50d-a7d58d5e5d5c
md"""
下記のようにテストを書くことができる．
"""

# ╔═╡ e8a37180-4777-4087-9155-4a3c37530d54
@testset "Point2D" begin
	o = Point2D(0, 0)
	@test o.x == 0.
	@test o.y == 0.
	p = Point2D(1, 2)
	@test p.x == 1
	@test p.y == 2
	d1 = distance(p, o)
	d2 = distance(o, p)
	@test d1 == d2 ≈ √5
end

# ╔═╡ e10f9104-9905-484a-9144-18d9467ec743
md"""
## 練習問題

- `distance` 関数があるセルを編集し敢えて間違った実装に改変してみよ． 
- このとき，上記のテストが正常に動かなくなることを確認せよ
- Point2D の実装を元に3次元版を実装してみよ
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Example = "7876af07-990d-54b4-ab0e-23690620f79a"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
Example = "~0.5.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "49929da40b6ade982f27b790233147f997c9e104"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Example]]
git-tree-sha1 = "46e44e869b4d90b96bd8ed1fdcf32244fddfb6cc"
uuid = "7876af07-990d-54b4-ab0e-23690620f79a"
version = "0.5.3"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
"""

# ╔═╡ Cell order:
# ╟─d7f1cb15-6ab1-4fd6-8639-e20ded4f29ea
# ╠═5be67619-2fa6-48e7-b8e4-492360b6c9c0
# ╠═bdb66039-1c23-4172-afef-ffab8452ddf6
# ╠═408ea404-a897-4574-97e0-5cab52b9d154
# ╠═04c275d6-ddd3-11ee-0102-4ff69f77eaac
# ╟─aeed025d-ce2d-46ca-be64-3aea35e51932
# ╠═ee932197-1a89-4ada-b5eb-6e0514c2f47c
# ╟─4c482f7b-0f75-4b77-ad79-62ce700d5834
# ╟─f853d44f-cb09-4fd2-8afa-ea89821c829b
# ╠═b5c5f82b-5cd6-4b71-8b9f-d4b9715d21b7
# ╠═a0f3c104-da2d-4540-ae62-84a9f2dee111
# ╠═ff5b2387-b22b-423d-8b4c-e99a085f5293
# ╠═44641fc0-7c8c-46f6-be51-438be0ed742d
# ╟─e2b487c7-7d56-4cb1-b50d-a7d58d5e5d5c
# ╠═e8a37180-4777-4087-9155-4a3c37530d54
# ╟─e10f9104-9905-484a-9144-18d9467ec743
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002