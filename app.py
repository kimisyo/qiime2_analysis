import matplotlib as mpl
import streamlit as st
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from functools import cmp_to_key

# df_barのサンプルIDを、カテゴリ、菌種の存在量、サンプルID順にソート


def index_compare(data1, data2):
    data1s = data1.split("|")
    data2s = data2.split("|")
    value1 = df_bar.loc[data1s[0], sort_class]
    value2 = df_bar.loc[data2s[0], sort_class]
    if data1s[1] < data2s[1]:
        return 1
    elif data1s[1] > data2s[1]:
        return -1
    else:
        # category
        if value1 < value2:
            return 1
        elif value1 > value2:
            return -1
        else:
            # Sample ID
            if data1s[0] < data2s[0]:
                return -1
            elif data1s[0] > data2s[0]:
                return 1
            else:
                return 0


# 色の設定
colors = [(56, 108, 176), (255, 255, 153), (253, 192, 134), (190, 174, 212),
          (127, 201, 127), (102, 102, 102), (191, 91, 23), (240, 2, 127)]
colors_new = []
for color in colors:
    colors_new.append(tuple(np.array(color)/255))

# グローバル変数
metadata_file = None
taxonomy_file = None
heatmap_file = None
df_metadata = None
df_bar = None
df_heatmap = None

st.sidebar.title("Qiime2向け腸内細菌解析アプリ")

metadata_file = st.sidebar.file_uploader(
    'メタデータファイルを選択してください', type=['tsv'])

if metadata_file:
    df_metadata = pd.read_table(metadata_file, index_col=0)

    taxonomy_file = st.sidebar.file_uploader(
        'taxonomyファイルを指定してください', type=['csv'])

    if taxonomy_file:
        df_bar = pd.read_csv(taxonomy_file, index_col=0)\
            .reindex(index=df_metadata.index)\
            .iloc[:, :-(len(df_metadata.columns))]

        heatmap_file = st.sidebar.file_uploader(
            'heatmapファイルを指定してください', type=['tsv', 'txt'])
        if heatmap_file:
            df_heatmap = pd.read_csv(heatmap_file, delimiter="\t", index_col=0)\
                .reindex(columns=df_metadata.index)

st.markdown("### 解析結果")
if metadata_file and taxonomy_file:

    category = st.selectbox('カテゴリIDを指定してください', df_metadata.columns)

    # サンプルIDを、サンプルID-カテゴリに変換
    sample_ids = [label + "|" +
                  str(cat) for (label, cat) in zip(df_metadata.index, df_metadata[category].values)]

    # 存在量の最も大きな菌種を求める
    df_sum_by_class = df_bar.sum(axis=0)\
        .sort_values(ascending=False)
    sort_class = df_sum_by_class.index.values[0]

    # サンプルID毎の菌種の存在量を、サンプル内の割合(%)に変換
    df_sum_by_sample = df_bar.sum(axis=1)
    for column in df_bar.columns:
        df_bar[column] = (df_bar[column] / df_sum_by_sample) * 100

    # サンプルID-菌種テーブルの列を存在量の大きな菌種順に並べ替え
    df_bar = df_bar.reindex(columns=df_sum_by_class.index)

    # サンプルID-Levelテーブルの列をサンプルID順に並べ替え
    sample_ids_sort = sample_ids.copy()
    sample_ids_sort.sort(key=cmp_to_key(index_compare))
    df_bar.index = sample_ids
    df_bar = df_bar.reindex(index=sample_ids_sort)

    # グラフの描画
    df_bottom = pd.DataFrame(index=df_bar.index)

    for i, column in enumerate(df_bar.columns):
        df_bottom[column] = np.ones(len(df_bar)) * 100
        for j, column_j in enumerate(df_bar.columns):
            if j > i:
                break
            df_bottom[column] -= df_bar.iloc[:, j]

    fig_bar = plt.figure(figsize=(7, 5))
    barWidth = 0.99
    for i, column in enumerate(df_bar.columns):
        print(column)
        plt.bar(df_bar.index, df_bar[column], bottom=df_bottom[column], color=colors_new[i % len(
            colors)], edgecolor='white', width=barWidth, label=column)

    plt.ylabel("Relative Frequency (%)")
    plt.xlabel("Sample")
    plt.xticks(rotation=90)
    plt.ylim(0, 100)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
    mpl.rcParams['axes.xmargin'] = 0
    mpl.rcParams['axes.ymargin'] = 0
    st.pyplot(fig_bar)

    st.download_button('Download CSV', df_bar.to_csv(),
                       'result.csv', 'text/csv')

    if heatmap_file:
        df_heatmap.columns = sample_ids
        df_heatmap = df_heatmap.reindex(columns=sample_ids_sort)
        df_heatmap = df_heatmap.reindex(index=df_sum_by_class.index)
        fig_heatmap = plt.figure(figsize=(8, 3))
        sns.heatmap(df_heatmap,
                    cmap='coolwarm',
                    xticklabels=True,
                    yticklabels=True,
                    cbar_kws={'label': 'log10 frequency'})

        st.pyplot(fig_heatmap)
