# 微動アレイ観測解析手法 - 技術ノート

## 1. 微動アレイ観測の基礎理論

### 1.1 微動とは
微動（microtremor）は、地盤に常時存在する微小振動であり、主に以下の要因により発生する：
- 海洋波浪（脈動：0.1-1 Hz）
- 交通振動（1-20 Hz）
- 風による構造物の振動
- 人間活動による振動

### 1.2 地震波の種類（実体波と表面波）

地震波は大きく実体波（body waves）と表面波（surface waves）に分類される。

#### 実体波
実体波は地盤内部を伝播する波動で、以下の2種類がある：

1. **P波（Primary wave、縦波）**
   - 粒子の振動方向と波の進行方向が平行
   - 速度：$v_p = \sqrt{\frac{\lambda + 2\mu}{\rho}}$
   - ここで、$\lambda$, $\mu$：ラメ定数、$\rho$：密度

2. **S波（Secondary wave or Shear Wave、横波）**
   - 粒子の振動方向と波の進行方向が垂直
   - 速度：$v_s = \sqrt{\frac{\mu}{\rho}}$
   - 液体中では伝播しない（$\mu = 0$）

#### 表面波
表面波は地表面に沿って伝播する波動で、振幅が深さとともに指数関数的に減衰する：

1. **レイリー波（Rayleigh wave）**
   - 鉛直面内での楕円運動
   - 微動の主成分（約70%）
   - 速度：約$0.92v_s$（ポアソン比0.25の場合）

2. **ラブ波（Love wave）**
   - 水平面内での振動
   - 速度構造に不連続がある場合に発生
   - SH波の重ね合わせとして表現

### 1.3 位相速度と群速度

#### 位相速度（Phase velocity）
位相速度$c$は、単一周波数の波の位相が伝播する速度：

$$c = \frac{\omega}{k} = f\lambda$$

ここで：
- $\omega = 2\pi f$：角周波数
- $k = 2\pi/\lambda$：波数
- $f$：周波数
- $\lambda$：波長

#### 群速度（Group velocity）
群速度$U$は、波束（エネルギー）が伝播する速度：

$$U = \frac{d\omega}{dk} = c - \lambda\frac{dc}{d\lambda}$$

#### 2つの正弦波の重ね合わせによる理解

位相速度と群速度の違いを理解するため、わずかに異なる周波数を持つ2つの正弦波の重ね合わせを考える：

$$u(x,t) = \sin(k_1x - \omega_1t) + \sin(k_2x - \omega_2t)$$

ここで、$k_1 = k - \Delta k$、$k_2 = k + \Delta k$、$\omega_1 = \omega - \Delta\omega$、$\omega_2 = \omega + \Delta\omega$とすると、

三角関数の和積公式を用いて：

$$u(x,t) = 2\cos(\Delta k \cdot x - \Delta\omega \cdot t) \sin(kx - \omega t)$$

この式は以下のように解釈できる：

1. **搬送波（carrier wave）**：$\sin(kx - \omega t)$
   - 位相速度で伝播：$c = \frac{\omega}{k}$
   - 個々の波の山や谷の移動速度

2. **包絡線（envelope）**：$2\cos(\Delta k \cdot x - \Delta\omega \cdot t)$
   - 群速度で伝播：$U = \frac{\Delta\omega}{\Delta k} \rightarrow \frac{d\omega}{dk}$（$\Delta k \rightarrow 0$の極限）
   - 波束（wave packet）全体の移動速度
   - エネルギーの伝播速度

#### 物理的意味

- **位相速度**：単一周波数成分の位相が進む速度。実際には観測できない概念的な速度
- **群速度**：波のエネルギーや情報が伝わる速度。実際に観測される波束の移動速度

分散性媒質では$c \neq U$となり：
- 正常分散（$\frac{dc}{d\omega} > 0$）：$U < c$
- 異常分散（$\frac{dc}{d\omega} < 0$）：$U > c$

#### 分散性
媒質が不均質な場合、位相速度は周波数に依存し（分散性）、これにより地下構造の推定が可能となる。正常分散の場合：
- 高周波数（短波長）→ 浅部の情報
- 低周波数（長波長）→ 深部の情報

### 1.4 表面波の分散性
微動の主成分は表面波（レイリー波）であり、その位相速度は周波数（または波長）に依存する。この分散性により、地下構造を推定することが可能となる。

位相速度 $c(f)$ と波長 $\lambda$ の関係：

$$c(f) = f \times \lambda$$


探査深度の目安：

$$\text{探査深度} \approx \frac{\lambda}{2} \sim \frac{\lambda}{3}$$

## 2. SPAC法（空間自己相関法）

### 2.1 理論的背景
SPAC法（Spatial Autocorrelation Method）は、Aki (1957) により提案された手法で、微動の空間的な相関を利用して位相速度を推定する。

#### なぜ空間相関から位相速度が推定できるのか

##### 直感的な理解
2つの観測点で同じ波を観測する場合を考える：
1. **観測点間距離が波長の整数倍**：波形が完全に一致（相関係数 = 1）
2. **観測点間距離が波長の半整数倍**：波形が逆位相（相関係数 = -1）
3. **その他の距離**：中間的な相関値

この相関の周期性から波長を推定でき、周波数との関係から位相速度が求まる。

##### 数学的な導出：2点間相関からベッセル関数へ

座標原点に観測点1、そこから距離$r$離れた位置に観測点2を配置する。SPAC法では各周波数成分を独立に解析するため、特定の角周波数$\omega$（周波数$f = \omega/2\pi$）の成分に着目する。方位角$\theta$から到来する平面波を考えると、各観測点での波動場は：

観測点1（原点）：
$$u_1(t, \omega) = \int_0^{2\pi} A(\theta, \omega) \exp[-i\omega t] d\theta \tag{2-1}$$

観測点2（距離$r$）：
$$u_2(t, \omega) = \int_0^{2\pi} A(\theta, \omega) \exp[ikr\cos(\theta - \alpha) - i\omega t] d\theta \tag{2-2}$$

ここで：
- $A(\theta, \omega)$：方位角$\theta$方向からの周波数$\omega$の波の複素振幅
- $k = 2\pi f/c = \omega/c$：波数（周波数依存）
- $\alpha$：観測点2の方位角
- $r\cos(\theta - \alpha)$：波の進行方向への投影距離

**ステップ1：相関関数の計算**

時間平均を用いた相関関数（特定周波数$\omega$について）：
$$\langle u_1(t, \omega) u_2^*(t, \omega) \rangle = \lim_{T \to \infty} \frac{1}{T} \int_0^T u_1(t, \omega) u_2^*(t, \omega) dt \tag{2-3}$$

式(2-1)と(2-2)を代入：
$$\langle u_1(t, \omega) u_2^*(t, \omega) \rangle = \lim_{T \to \infty} \frac{1}{T} \int_0^T \left[\int_0^{2\pi} A(\theta_1, \omega) e^{-i\omega t} d\theta_1\right] \left[\int_0^{2\pi} A^*(\theta_2, \omega) e^{-ikr\cos(\theta_2 - \alpha) + i\omega t} d\theta_2\right] dt \tag{2-4}$$

**ステップ2：時間積分の実行**

式(2-4)の時間積分を先に実行する。積分の順序を交換して：
$$\langle u_1(t, \omega) u_2^*(t, \omega) \rangle = \int_0^{2\pi} \int_0^{2\pi} A(\theta_1, \omega) A^*(\theta_2, \omega) e^{-ikr\cos(\theta_2 - \alpha)} \left[\lim_{T \to \infty} \frac{1}{T} \int_0^T e^{-i\omega t} e^{i\omega t} dt\right] d\theta_1 d\theta_2 \tag{2-4a}$$

ここで重要な点は、微動が定常確率過程であることから、異なる方向からの波は統計的に独立（非相関）であると仮定できる。すなわち、特定周波数$\omega$において：
$$A(\theta_1, \omega) A^*(\theta_2, \omega) = |A(\theta, \omega)|^2 \delta(\theta_1 - \theta_2) \tag{2-4b}$$

この仮定により、時間積分は：
$$\lim_{T \to \infty} \frac{1}{T} \int_0^T dt = 1 \quad \text{（$\theta_1 = \theta_2$のとき）} \tag{2-4c}$$

したがって：
$$\langle u_1(t, \omega) u_2^*(t, \omega) \rangle = \int_0^{2\pi} \int_0^{2\pi} |A(\theta, \omega)|^2 e^{-ikr\cos(\theta - \alpha)} \delta(\theta_1 - \theta_2) d\theta_1 d\theta_2 \tag{2-5}$$

デルタ関数により$\theta_1$についての積分が実行され：
$$\langle u_1(t, \omega) u_2^*(t, \omega) \rangle = \int_0^{2\pi} |A(\theta, \omega)|^2 e^{-ikr\cos(\theta - \alpha)} d\theta \tag{2-5a}$$

等方的な波動場の仮定（$|A(\theta, \omega)|^2 = \text{const} = A_0^2(\omega)$）より：
$$\langle u_1(t, \omega) u_2^*(t, \omega) \rangle = A_0^2(\omega) \int_0^{2\pi} e^{-ikr\cos(\theta - \alpha)} d\theta \tag{2-6}$$

**ステップ3：パワーの計算**

同様に各点でのパワー（特定周波数$\omega$について）：
$$\langle |u_1(t, \omega)|^2 \rangle = \langle |u_2(t, \omega)|^2 \rangle = 2\pi A_0^2(\omega) \tag{2-7}$$

**ステップ4：正規化された相関係数**

$$\rho(r, \omega) = \frac{\langle u_1(t, \omega) u_2^*(t, \omega) \rangle}{\sqrt{\langle |u_1(t, \omega)|^2 \rangle \langle |u_2(t, \omega)|^2 \rangle}} = \frac{A_0^2(\omega) \int_0^{2\pi} e^{-ikr\cos(\theta - \alpha)} d\theta}{2\pi A_0^2(\omega)} \tag{2-8}$$

観測点2の方位$\alpha$によらないことを示すため、変数変換$\psi = \theta - \alpha$：
$$\rho(r, f) = \frac{1}{2\pi} \int_0^{2\pi} e^{-ikr\cos\psi} d\psi \tag{2-9}$$

**ステップ5：ベッセル関数の認識**

オイラーの公式より$e^{-ikr\cos\psi} = \cos(kr\cos\psi) - i\sin(kr\cos\psi)$。
$\sin(kr\cos\psi)$の積分は奇関数のためゼロ：

$$\rho(r, f) = \frac{1}{2\pi} \int_0^{2\pi} \cos(kr\cos\psi) d\psi \tag{2-10}$$

これは第1種0次ベッセル関数の積分表現：
$$J_0(x) = \frac{1}{2\pi} \int_0^{2\pi} \cos(x\cos\psi) d\psi \tag{2-11}$$

したがって：
$$\rho(r, f) = J_0(kr) = J_0\left(\frac{2\pi rf}{c(f)}\right) \tag{2-12}$$

#### 物理的イメージ
- **波長が長い（低周波）**：観測点間の位相差が小さく、高い相関
- **波長が短い（高周波）**：観測点間の位相差が大きく、相関が振動
- **ベッセル関数の第1ゼロ点**：$2\pi rf/c \approx 2.4$のとき$\rho = 0$

これにより、相関係数がゼロになる周波数から直接位相速度を推定可能：
$$c = \frac{2\pi rf}{2.4} \approx 2.6 \cdot rf$$

### 2.2 基本原理
円形アレイ（中心点＋周辺観測点）での観測データから、空間自己相関係数を計算：

$$\rho(r, f) = \frac{2 \cdot \text{Re}[S_{12}(f)]}{\sqrt{S_{11}(f) \cdot S_{22}(f)}}$$

ここで：
- $S_{12}(f)$: クロススペクトル
- $S_{11}(f), S_{22}(f)$: オートスペクトル
- $r$: 観測点間距離

### 2.3 位相速度の推定
空間自己相関係数は第1種0次ベッセル関数で表現される：

$$\rho(r, f) = J_0\left(\frac{2\pi rf}{c(f)}\right)$$

これを逆解析することで位相速度 $c(f)$ を求める。

### 2.4 実装上の注意点
1. **時間窓の設定**：10-30分程度の連続記録を使用
2. **セグメント分割**：20-40秒程度に分割し、平均化
3. **窓関数**：Hanning窓などを適用してスペクトルリーケージを低減
4. **ノイズ除去**：過渡的なノイズを含むセグメントは除外

## 3. FK法（周波数-波数法）

### 3.1 理論的背景
FK法（Frequency-Wavenumber Method）は、アレイ観測データを周波数-波数領域で解析し、表面波の伝播方向と位相速度を推定する。

### 3.2 基本原理
観測データ $u(x, t)$ のフーリエ変換：

$$U(k, \omega) = \int\int u(x, t) \exp(-i(kx - \omega t)) \, dx \, dt$$

パワースペクトル密度：

$$P(k, \omega) = |U(k, \omega)|^2$$

### 3.3 ビームフォーミング
最大尤度法やビームフォーミングにより、各周波数での卓越波数を推定：

$$c(f) = \frac{2\pi f}{|k_{\max}|}$$

## 4. 分散曲線の算出と品質管理

### 4.1 分散曲線の結合
- 異なるアレイサイズからの結果を結合
- 周波数帯域の重複部分で連続性を確認

### 4.2 品質評価指標
1. **信号対雑音比（SNR）**
2. **空間自己相関係数の標準偏差**
3. **FK解析でのピークの明瞭度**

### 4.3 異常値の除去
- 物理的に不合理な速度値（例：表層より深部で低速）
- 隣接周波数との不連続性が大きい点

## 5. 逆解析手法

### 5.1 地下構造モデルの構築
初期モデルの設定：
- 層数、層厚、S波速度、P波速度、密度
- 既存の地質情報や物理探査結果を参考

### 5.2 順解析
レイリー波の理論分散曲線の計算：
- Thomson-Haskell法
- 反射透過係数法
- 有限要素法

### 5.3 最適化手法
1. **遺伝的アルゴリズム（GA）**
   - グローバルサーチに適している
   - 計算時間が長い

2. **シミュレーテッドアニーリング（SA）**
   - 局所解に陥りにくい
   - パラメータ調整が重要

3. **最小二乗法**
   - 収束が速い
   - 初期値依存性が高い

### 5.4 誤差評価
観測分散曲線と理論分散曲線の残差：

$$\text{misfit} = \sqrt{\frac{\sum\left[\frac{(c_{\text{obs}} - c_{\text{cal}})^2}{\sigma^2}\right]}{N}}$$

## 6. 実装上のトラブルシューティング

### 6.1 一般的な問題と対策

1. **低周波数でのSPAC係数の不安定性**
   - 原因：アレイサイズが波長に対して小さい
   - 対策：より大きなアレイを使用、または解析下限周波数を上げる

2. **高周波数でのコヒーレンスの低下**
   - 原因：局所的な不均質性、非平面波の影響
   - 対策：観測点間隔を狭める、平均化時間を長くする

3. **分散曲線の不連続**
   - 原因：モード判定の誤り、ノイズの影響
   - 対策：隣接周波数との連続性をチェック、品質の低いデータを除外

### 6.2 データ処理のベストプラクティス

1. **前処理**
   - 機器応答の補正
   - トレンド除去（線形・多項式）
   - バンドパスフィルタの適用

2. **品質管理**
   - 時系列データの目視確認
   - スペクトルの安定性チェック
   - 複数手法での結果比較

3. **解釈上の注意**
   - 高次モードの影響
   - 横方向不均質性の影響
   - 解の非唯一性

## 7. 参考文献と推奨リソース

### 基本文献
- Aki, K. (1957): Space and time spectra of stationary stochastic waves
- Okada, H. (2003): The Microtremor Survey Method, SEG

### 解析ソフトウェア
- Geopsy (オープンソース)
- MASW (商用)
- 自作コード（Python/MATLAB）

### データフォーマット
- SEG-Y, SAC, MiniSEED
- 時刻同期の重要性（GPS時刻推奨）

## 付録A: 探査深度の理論的背景

### A.1 なぜ探査深度が λ/2 ～ λ/3 なのか

表面波（特にレイリー波）の探査深度が波長の1/2から1/3程度になる理由は、波動の振幅が深さとともに指数関数的に減衰することに起因する。

### A.2 レイリー波の変位の深さ分布

レイリー波の変位は、深さ$z$に対して以下のように表現される：

水平成分：
$$u_x(z) = A \left[ e^{-k\alpha z} - \frac{2\alpha\beta}{1+\beta^2} e^{-k\beta z} \right]$$

鉛直成分：
$$u_z(z) = A \frac{i\alpha}{k} \left[ e^{-k\alpha z} - \frac{2}{1+\beta^2} e^{-k\beta z} \right]$$

ここで：
- $k = 2\pi/\lambda$：波数
- $\alpha = \sqrt{1 - (c/v_p)^2}$
- $\beta = \sqrt{1 - (c/v_s)^2}$
- $c$：レイリー波の位相速度
- $v_p, v_s$：P波、S波速度

### A.3 感度カーネル（Sensitivity Kernel）

レイリー波の位相速度は、各深さの弾性定数の重み付き平均として表現できる。この重み関数（感度カーネル）$K(z)$は：

$$\frac{\partial c}{\partial v_s(z)} = K(z) \propto \left| \frac{\partial u}{\partial z} \right|^2$$

感度カーネルは通常、深さ$z = \lambda/3$付近で最大値を取る。これは：

1. **浅部（$z < \lambda/4$）**：変位は大きいが、ひずみ（変位の空間微分）は小さい
2. **深部（$z > \lambda/2$）**：振幅が急速に減衰し、寄与が小さい
3. **中間深度（$z \approx \lambda/3$）**：変位とひずみのバランスが最適

### A.4 実用的な意味

探査深度の目安として$\lambda/2 \sim \lambda/3$を用いる理由：

1. **エネルギー分布**：レイリー波のエネルギーの約70%がこの深度範囲に集中
2. **感度の観点**：この深度の速度構造変化が位相速度に最も影響
3. **経験則**：多数の実測例から導かれた実用的な関係

### A.5 周波数と探査深度の関係

具体的な例（$v_s = 200$ m/sの場合）：

| 周波数 (Hz) | 波長 (m) | 探査深度 (m) |
|------------|----------|--------------|
| 1          | 200      | 65-100       |
| 5          | 40       | 13-20        |
| 10         | 20       | 6.5-10       |
| 20         | 10       | 3.3-5        |

この関係により、異なる周波数の分散曲線から、異なる深さの地下構造情報を得ることができる。