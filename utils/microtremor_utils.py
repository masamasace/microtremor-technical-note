"""
微動アレイ観測解析のためのユーティリティ関数

このモジュールには、微動データの前処理、SPAC解析、FK解析、
分散曲線の計算など、実践的な解析に必要な関数が含まれています。
"""

import numpy as np
from scipy import signal
from scipy.special import j0, j1
from scipy.optimize import curve_fit, minimize
import matplotlib.pyplot as plt


class MicrotremorArray:
    """微動アレイ観測データを扱うクラス"""
    
    def __init__(self, data, fs, coordinates):
        """
        Parameters:
        ----------
        data : array_like
            観測データ (n_stations × n_samples)
        fs : float
            サンプリング周波数 (Hz)
        coordinates : array_like
            観測点座標 [[x0,y0], [x1,y1], ...] (m)
        """
        self.data = np.array(data)
        self.fs = fs
        self.coordinates = np.array(coordinates)
        self.n_stations = data.shape[0]
        self.n_samples = data.shape[1]
        
    def preprocess(self, detrend='linear', taper=0.05):
        """データの前処理（トレンド除去、テーパー処理）"""
        processed_data = np.copy(self.data)
        
        for i in range(self.n_stations):
            # トレンド除去
            processed_data[i, :] = signal.detrend(processed_data[i, :], type=detrend)
            
            # テーパー処理
            if taper > 0:
                window = signal.tukey(self.n_samples, alpha=taper)
                processed_data[i, :] *= window
                
        return processed_data
    
    def bandpass_filter(self, freqmin, freqmax, corners=4):
        """バンドパスフィルタ"""
        nyquist = 0.5 * self.fs
        low = freqmin / nyquist
        high = freqmax / nyquist
        
        b, a = signal.butter(corners, [low, high], btype='band')
        filtered_data = np.zeros_like(self.data)
        
        for i in range(self.n_stations):
            filtered_data[i, :] = signal.filtfilt(b, a, self.data[i, :])
            
        return filtered_data
    
    def compute_distances(self):
        """全観測点ペア間の距離を計算"""
        distances = {}
        for i in range(self.n_stations):
            for j in range(i+1, self.n_stations):
                dx = self.coordinates[i, 0] - self.coordinates[j, 0]
                dy = self.coordinates[i, 1] - self.coordinates[j, 1]
                distances[(i, j)] = np.sqrt(dx**2 + dy**2)
        return distances


def spac_coefficient(data1, data2, fs, nperseg=None, noverlap=None):
    """
    2つの観測点間のSPAC係数を計算
    
    Parameters:
    ----------
    data1, data2 : array_like
        2つの観測点の時系列データ
    fs : float
        サンプリング周波数
    nperseg : int
        セグメント長（サンプル数）
    noverlap : int
        オーバーラップ長（サンプル数）
    
    Returns:
    -------
    frequencies : array
        周波数配列
    spac : array
        SPAC係数
    coherence : array
        コヒーレンス
    """
    if nperseg is None:
        nperseg = min(len(data1) // 8, 2048)
    if noverlap is None:
        noverlap = nperseg // 2
        
    # クロススペクトル密度
    f, Pxy = signal.csd(data1, data2, fs=fs, nperseg=nperseg, noverlap=noverlap)
    
    # オートスペクトル密度
    _, Pxx = signal.welch(data1, fs=fs, nperseg=nperseg, noverlap=noverlap)
    _, Pyy = signal.welch(data2, fs=fs, nperseg=nperseg, noverlap=noverlap)
    
    # SPAC係数（実部のみ）
    spac = np.real(Pxy) / np.sqrt(Pxx * Pyy)
    
    # コヒーレンス
    coherence = np.abs(Pxy)**2 / (Pxx * Pyy)
    
    return f, spac, coherence


def estimate_phase_velocity_spac(frequencies, spac_values, distance):
    """
    SPAC係数から位相速度を推定
    
    Parameters:
    ----------
    frequencies : array_like
        周波数配列
    spac_values : array_like
        SPAC係数配列
    distance : float
        観測点間距離
    
    Returns:
    -------
    phase_velocities : array
        推定された位相速度
    """
    phase_velocities = []
    
    for f, spac in zip(frequencies, spac_values):
        if f == 0 or abs(spac) > 1:
            phase_velocities.append(np.nan)
            continue
            
        # ベッセル関数のゼロ点を探索
        def bessel_func(c):
            return j0(2 * np.pi * f * distance / c) - spac
        
        # 探索範囲
        c_min = 50
        c_max = 2000
        
        try:
            from scipy.optimize import brentq
            c = brentq(bessel_func, c_min, c_max)
            phase_velocities.append(c)
        except:
            phase_velocities.append(np.nan)
            
    return np.array(phase_velocities)


def fk_analysis(data, coordinates, fs, freq_range, slowness_range, n_slowness=50):
    """
    FK法（周波数-波数法）による解析
    
    Parameters:
    ----------
    data : array_like
        観測データ (n_stations × n_samples)
    coordinates : array_like
        観測点座標 [[x0,y0], [x1,y1], ...]
    fs : float
        サンプリング周波数
    freq_range : tuple
        解析周波数範囲 (fmin, fmax)
    slowness_range : tuple
        スローネス範囲 (smin, smax) [s/m]
    n_slowness : int
        スローネスのグリッド数
    
    Returns:
    -------
    frequencies : array
        周波数配列
    phase_velocities : array
        推定された位相速度
    """
    n_stations = data.shape[0]
    
    # FFT
    nfft = 2**int(np.ceil(np.log2(data.shape[1])))
    fft_data = np.fft.fft(data, n=nfft, axis=1)
    frequencies = np.fft.fftfreq(nfft, 1/fs)
    
    # 正の周波数のみ
    freq_mask = (frequencies >= freq_range[0]) & (frequencies <= freq_range[1])
    frequencies = frequencies[freq_mask]
    fft_data = fft_data[:, freq_mask]
    
    # スローネスグリッド
    sx = np.linspace(-slowness_range[1], slowness_range[1], n_slowness)
    sy = np.linspace(-slowness_range[1], slowness_range[1], n_slowness)
    SX, SY = np.meshgrid(sx, sy)
    
    phase_velocities = []
    
    for freq_idx, f in enumerate(frequencies):
        if f == 0:
            phase_velocities.append(np.nan)
            continue
            
        # ビームパワー計算
        beam_power = np.zeros((n_slowness, n_slowness))
        
        for i in range(n_slowness):
            for j in range(n_slowness):
                # ステアリングベクトル
                steering = np.exp(-1j * 2 * np.pi * f * 
                                 (SX[i,j] * coordinates[:, 0] + 
                                  SY[i,j] * coordinates[:, 1]))
                
                # ビームフォーミング
                beam = np.sum(fft_data[:, freq_idx] * steering)
                beam_power[i, j] = np.abs(beam)**2
        
        # 最大パワーの位置
        max_idx = np.unravel_index(np.argmax(beam_power), beam_power.shape)
        sx_max = SX[max_idx]
        sy_max = SY[max_idx]
        
        # スローネスから位相速度へ
        slowness = np.sqrt(sx_max**2 + sy_max**2)
        if slowness > 0:
            phase_velocities.append(1.0 / slowness)
        else:
            phase_velocities.append(np.nan)
    
    return frequencies, np.array(phase_velocities)


def combine_dispersion_curves(curves_list, freq_ranges=None):
    """
    複数のアレイサイズからの分散曲線を結合
    
    Parameters:
    ----------
    curves_list : list of tuples
        [(freq1, vel1), (freq2, vel2), ...] の形式
    freq_ranges : list of tuples
        各曲線の有効周波数範囲 [(fmin1, fmax1), ...]
    
    Returns:
    -------
    combined_freq : array
        結合された周波数配列
    combined_vel : array
        結合された位相速度配列
    """
    all_freq = []
    all_vel = []
    
    for idx, (freq, vel) in enumerate(curves_list):
        # 有効範囲でフィルタリング
        if freq_ranges is not None:
            fmin, fmax = freq_ranges[idx]
            mask = (freq >= fmin) & (freq <= fmax) & ~np.isnan(vel)
        else:
            mask = ~np.isnan(vel)
            
        all_freq.extend(freq[mask])
        all_vel.extend(vel[mask])
    
    # 周波数でソート
    sorted_idx = np.argsort(all_freq)
    combined_freq = np.array(all_freq)[sorted_idx]
    combined_vel = np.array(all_vel)[sorted_idx]
    
    return combined_freq, combined_vel


def theoretical_dispersion_curve(frequencies, layer_params):
    """
    簡易的な理論分散曲線の計算（多層モデル）
    
    Parameters:
    ----------
    frequencies : array_like
        周波数配列
    layer_params : list of dicts
        各層のパラメータ [{'h': 層厚, 'vs': S波速度}, ...]
        最下層は半無限とする（'h'は不要）
    
    Returns:
    -------
    phase_velocities : array
        理論位相速度
    """
    # ここでは簡易的な近似式を使用
    # 実際にはThomson-Haskell法などを実装する必要がある
    
    phase_velocities = []
    
    for f in frequencies:
        if f == 0:
            phase_velocities.append(np.nan)
            continue
            
        # 各層の寄与を考慮
        weighted_vel = 0
        total_weight = 0
        
        depth = 0
        for i, layer in enumerate(layer_params):
            vs = layer['vs']
            
            if 'h' in layer:  # 最下層以外
                h = layer['h']
                # 波長に対する層の寄与
                wavelength = vs / f
                weight = np.exp(-4 * np.pi * depth / wavelength) * h
                weighted_vel += vs * weight
                total_weight += weight
                depth += h
            else:  # 最下層（半無限）
                wavelength = vs / f
                weight = np.exp(-4 * np.pi * depth / wavelength)
                weighted_vel += vs * weight
                total_weight += weight
        
        if total_weight > 0:
            phase_velocities.append(weighted_vel / total_weight)
        else:
            phase_velocities.append(np.nan)
    
    return np.array(phase_velocities)


def plot_dispersion_curve(frequencies, velocities, title='Dispersion Curve', 
                         save_path=None, **kwargs):
    """
    分散曲線のプロット
    
    Parameters:
    ----------
    frequencies : array_like or list
        周波数配列（複数可）
    velocities : array_like or list
        位相速度配列（複数可）
    title : str
        グラフタイトル
    save_path : str
        保存パス（Noneの場合は保存しない）
    **kwargs : dict
        追加のプロットオプション
    """
    plt.figure(figsize=(10, 6))
    
    # 単一曲線か複数曲線かを判定
    if isinstance(frequencies, list):
        for i, (f, v) in enumerate(zip(frequencies, velocities)):
            label = kwargs.get('labels', [f'Curve {i+1}' for _ in frequencies])[i]
            plt.scatter(f, v, s=30, alpha=0.7, label=label)
    else:
        plt.scatter(frequencies, velocities, s=30, alpha=0.7)
    
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Phase Velocity (m/s)')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    
    if 'xlim' in kwargs:
        plt.xlim(kwargs['xlim'])
    if 'ylim' in kwargs:
        plt.ylim(kwargs['ylim'])
    
    if isinstance(frequencies, list) or 'labels' in kwargs:
        plt.legend()
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    
    plt.show()


def quality_check(frequencies, spac_coefficients, coherence, min_coherence=0.8):
    """
    データ品質チェック
    
    Parameters:
    ----------
    frequencies : array_like
        周波数配列
    spac_coefficients : array_like
        SPAC係数
    coherence : array_like
        コヒーレンス
    min_coherence : float
        最小コヒーレンス閾値
    
    Returns:
    -------
    quality_mask : array
        品質基準を満たす要素がTrue
    quality_score : float
        全体的な品質スコア（0-1）
    """
    quality_mask = np.ones_like(frequencies, dtype=bool)
    
    # コヒーレンスチェック
    quality_mask &= (coherence >= min_coherence)
    
    # SPAC係数の範囲チェック
    quality_mask &= (np.abs(spac_coefficients) <= 1.0)
    
    # 非物理的な値を除外
    quality_mask &= ~np.isnan(spac_coefficients)
    quality_mask &= ~np.isinf(spac_coefficients)
    
    # 品質スコア
    quality_score = np.sum(quality_mask) / len(quality_mask)
    
    return quality_mask, quality_score