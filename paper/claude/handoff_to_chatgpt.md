# Claude → ChatGPT 引き継ぎ報告

**日付**: 2026-10-26
**対象**: `paper/claude/manuscript.md` (399行, Markdown原稿)

---

## 1. Claudeが完了した作業

### 1.1 SubGraph修正後の数値全面検証・更新

SubGraph128状態マッピングのバグ修正（`pattern_id := .I` → 実際のSubGraph行対応）により、47/51疾患のBasin割り当てが変化した。修正後のデータ（`patient_scores.tsv`, `energy_sanity_check_per_basin.tsv`, `basin_code_enrichment.tsv`, `state_basin_map.tsv`）を全件突合し、以下を更新済み:

| 項目 | 修正前 | 修正後 |
|---|---|---|
| TB全体HR | 1.088 | 1.069 (per-basin Cox) |
| 認知症Basin 1特徴コード | 不眠(47%) | 皮膚炎(38%) |
| 貧血Basin 1/3人数 | 2,091,203 / 46,809 | 2,167,645 / 31,362 |
| 貧血Basin 3エネルギー | +0.146 | −0.100 |
| 貧血サイズ比 | 45× | 69× |
| 甲状腺Basin 3 n | 982,324 | 1,002,031 |
| 甲状腺Basin 4 n | 139,993 | 135,495 |
| 甲状腺CKDコード頻度 | 1.8%,4.0%,7.7% etc. | 1.3%,3.2%,6.1% etc. |
| 甲状腺死亡率差 | 20.4pp (全疾患最大) | 10.0pp (最大ではない) |
| 甲状腺Basin 4特徴 | 最年少/最低死亡率 | 書き直し(76歳/47.3%) |

### 1.2 用語・方法論の修正

- **"per month" → "per step"**: 全箇所変更。Methods §4.7に「stride ≈ 1 month」の注記追加
- **survminer**: 使用していないことを確認、Methods ソフトウェアリストから削除
- **Figure S9**: ICD-10コード頻度推移ヒートマップの凡例を追加（代表4疾患）

### 1.3 確認済みで変更不要だった数値

- 45/51 死亡率per-basin majority基準 ✓
- 42/51 入院per-basin majority基準 ✓
- 認知症コホート n=929,951 ✓
- 認知症Basin分布 498,636/434,649 ✓
- Basin数分布 17/26/6/2 ✓
- PCA分散 12.5%/13.3% の説明 ✓
- TB軌跡 82.8% stable-in-basin ✓

### 1.4 実行中の計算

- **tagcloud_per_basin** (51疾患): Basin別ICD-10コード対数オッズスコアリングとタグクラウド生成。4/51完了、13:00-13:30頃に全完了予定。
  - 出力: `plot/{disease}/tagcloud_basin*.png`, `plot/{disease}/basin_code_enrichment.tsv`

### 1.5 tagcloud完了後に実行予定のスクリプト

`scratchpad/run_all_remaining.sh` に5フェーズをまとめてある:
1. **trajectory_to_death** (51疾患): 死亡前の状態遷移
2. **trajectory_to_hospitalization** (51疾患): 入院前の状態遷移
3. **stepwise_transition + summary**: 遷移率の再計算
4. **code_trajectory_per_basin** (代表4疾患): Basin別コード頻度推移ヒートマップ
5. **Fig05 data regeneration**: Figure 5のサマリーデータ再生成

---

## 2. 未完了・要再検証の項目

### 2.1 trajectory依存の数値 (run_all_remaining.sh完了後に要検証)

以下の数値は修正前のtrajectoryデータに基づいており、再計算後に変わる可能性がある:

**Figure 5関連:**
- 96.2% 1-basin疾患のfirst-to-last match
- 68% 3-4 basin疾患のmatch
- 31% 異Basin移動
- 1-4% Basin外への退出
- 60-70% Basin収束率 (非Basin開始)
- 2%未満 stable non-basin

**Table 1 分類:**
- heart_failureの分類 (現在 "Symmetric" だが rate=8.66%, asymmetry=0.268 → 閾値付近)
- 全疾患の分類カウント (Absorbing/Directional/Symmetric/Single-basin)

**個別疾患の遷移率:**
- 貧血: 22.6%/0.5% (42:1比)、78,584/76,120イベント
- 甲状腺: 10.2%/2.3%
- 認知症: 5.7%/5.9%

### 2.2 Jaccard clustering

tagcloud完了後、basin_code_enrichment.tsvから再計算が必要:
- 2大グループ (7 Acute/Frailty vs 43 Chronic/Lifestyle) の内訳
- サブクラスタの記述

### 2.3 Figure S8

Figure S8の対象figureプレースホルダーが未解決（build_docx.shで空行挿入のみ）。

### 2.4 docx再ビルド

全計算完了・数値更新後に `scratchpad/build_docx.sh` で再生成が必要。

---

## 3. 原稿で守るべき制約

1. **未実施の解析を実施済みと書かない**: 比例ハザード性検定、多重検定補正、競合リスク解析、感度解析は未実施
2. **DeSC-OnlineXXX-Landscaperの経緯は書かない**: 未公表の先行研究。本文・Discussion・参考文献に再導入不可
3. **著者・DOI・PMIDを創作しない**: 参考文献は実在確認済みのもののみ
4. **コード頻度推移から個人内発症順序を主張しない**: 集団レベルの頻度変化のみ
5. **Basin定義自体は変更しない**: 状態IDとラベルの対応のみ修正対象
6. **数値の「新しい=正しい」と判断しない**: 矛盾があれば元データに立ち返る

---

## 4. ファイル構成

```
paper/claude/
├── manuscript.md              # メイン原稿 (399行)
├── manuscript.docx            # pandocビルド済みdocx (要再ビルド)
├── author_queries.md          # VERIFY解決済み著者情報
├── evidence.tsv               # 数値根拠のエビデンストラッカー
├── revision_notes.md          # 修正履歴ノート
├── validation/                # 検証用スクリプト・結果
├── backup_20261024/           # 2024-10-24時点バックアップ
├── backup_20261025_pre_subgraph_fix/  # SubGraph修正前バックアップ
└── handoff_to_chatgpt.md      # 本ファイル

scratchpad/ (セッション一時ファイル)
├── run_all_remaining.sh       # 残り全フェーズ実行スクリプト
├── run_tagcloud_all.sh        # tagcloud実行スクリプト (実行中)
├── rebuild_transition_summary.R
├── rebuild_fig05_data.R
├── build_docx.sh              # docxビルドスクリプト
└── stale_numbers_to_verify.md # 検証状況トラッカー
```

---

## 5. 推奨アクション (ChatGPTへ)

1. **tagcloud + run_all_remaining.sh 完了を待つ** (推定 2026-10-26 15:00頃)
2. **§2.1の数値を再検証**: 特にFigure 5, Table 1分類, 各疾患遷移率
3. **Jaccard clustering記述を更新**: 新しいbasin_code_enrichment.tsvで確認
4. **Figure S8を解決**: 対象figureの特定と挿入
5. **docxを再ビルド**: `bash scratchpad/build_docx.sh`
6. **最終一貫性チェック**: 本文中の数値とデータファイルの全件突合
