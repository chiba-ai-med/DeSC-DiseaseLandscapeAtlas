# 改稿ノート (Revision Notes)

Updated: 2026-10-24

## 主要な変更

### 1. 著者情報の追加
- 4名の著者名・8所属・2名の責任著者を挿入。SETTING指定の通り。
- Lead contactは未指定のためAQ-8として確認事項に残した。

### 2. 未公表の先行研究への言及を削除
- Introduction第3段落の「DeSC-OnlineXXX-Landscaper」への言及を全削除。
- 「以前は全患者一括でELAを行い限界があった→疾患群別へ」という未公表検討の経緯を削除。
- 対応する参考文献プレースホルダー（旧ref 7, 8）を削除。
- Introductionは公開済み先行研究（Hidalgo 2009, Jensen 2014, Masuda 2025等）と疾患群別Atlasの目的から再構成。

### 3. 96.2%保持率の定義修正（§4A）
- 「throughout follow-up」→「first and last observed basin assignments matched」に修正。
- 途中の離脱・復帰を許す端点比較であることを明記。
- 死亡者限定、疾患別中央値であることを本文・Summary・Figure 5 legendで統一。
- Figure 5 legendに「endpoint comparison that permits intermediate departures」を追記。

### 4. 縦断解析の対象修正（§4B）
- Methods 4.7: 「For all patients (not restricted to deceased patients)」を削除。
- stepwise_transition_rate.shの入力がtrajectories_dead.tsvであることを確認し、「deceased patients only」と明記。
- periods_before_deathが暦上の月数ではなく「ordinal position counting backward from the last retained observation」であることを明記。
- Figure legendの "months before death" に対応する注意をAQ-6として記録。

### 5. 貧血の記述修正（§4C）
- 「virtually irreversible absorption」「near-irreversible progression」を削除。
- 遷移確率比を45:1→42:1に修正（実数値 0.226456/0.005353 = 42.30）。
- 絶対遷移件数（78,584 vs 76,120）がほぼ同数であることを追記。
- 「The asymmetry in conditional probabilities reflects the large difference in basin sizes rather than a net unidirectional flow of patients」と明記。
- Discussion: 「effectively irreversible progression」→「highly asymmetric conditional transition probabilities」。

### 6. 臨床的解釈の修正（§4D）
- 甲状腺: E00-E07が全甲状腺疾患をカバーすることを明記。「hypothyroidism-CKD progression pathway」→「known clinical association between thyroid dysfunction and CKD (Iglesias and Díez, 2009)」に変更。因果関係の示唆を除去。
- 認知症: 「autonomic dysfunction and advanced frailty」→コード名の列挙のみに変更。サブタイプの断定を除去。Discussion: 「autonomic/frailty profile」→具体的なコード記述に変更。
- 歯科コード: 「healthy population with regular healthcare access」→「a population with regular healthcare contact; it does not establish that dental care is protective」に修正。

### 7. Cox・入院関連の修正（§4E）
- Figure 7: 入院パネルが「any all-cause admission as the event and total observation duration as the time variable; it does not represent time to first admission」であることをlegendに明記。
- 入院率を「proportion with at least one all-cause admission」として定義。
- 代表レコード: 「latest available observation period」→「observation window with the maximum row identifier (max old_row_id)」。暦順との対応を断定しない。

### 8. 結核の修正（§4F）
- 32% Non-basin: 時点別の構成比と全期間Non-basin（3.3%）の違いを明記。
- HR CI: 1.082–1.089 → 1.082–1.094に修正（実データ: 1.09372）。
- PC sign比較: 「all non-basin states had PC6 = +1」→「all 10 examined non-basin states had PC6 = +1」。top-10限定であることを明記。

### 9. FAIL疾患リストの確定
- 8疾患を特定・列挙: breast_cancer, heart_failure, hip_fracture, liver_cancer, pancreatic_cancer, sepsis, tuberculosis, urinary_tract_infection。
- 各疾患のbasin数とPASS数を括弧で表示。

### 10. Table 1の不整合を明示
- liver_cancerとUTIがPredictiveに分類されているがmortality majority基準FAILであることを本文に明記。[VERIFY]タグ付き。
- 境界値付近の疾患（peptic_ulcer, delirium, AF, AKI）の数値を記載。

### 11. 参考文献の充実
- 27箇所の[REF NEEDED]をすべて解消。16件の実在文献に置換。
- DOI、PMID（可能な場合）を付与。
- 追加: Masuda et al. 2025 (ELA tutorial), Schneidman et al. 2006 (Ising/max-entropy), Ezaki et al. 2017 (ELA neuroimaging)。
- Landscaperは公開論文なしのため、container URIとMethodsでの記述のみ。

### 12. Discussion: SubGraphマッピングの限界を追加
- Limitation 8として、SubGraph.tsvの行数が128未満である問題と、0ラベルの解釈の曖昧さを追記。

## Codex版の指摘への対応

### 採用した指摘
- SubGraph<128行の問題（AQ-1として記録）
- 死亡前「月数」の定義（AQ-6）
- 入院KMのイベント定義（AQ-5）
- PC寄与率の定義差異（AQ-4）
- 貧血の遷移件数がほぼ同数
- Table 1の境界値問題
- stepwise transitionが死亡者限定であること
- 代表レコードのold_row_id選択

### 不採用とした指摘
- PC1-age median 0.358（Codex版）: 実際にはsigned medianであり、absolute median 0.407が正しい。本文は|r|=0.41を維持。
- UTI PC7-sex -0.319（Codex版のcor_energy列を参照した可能性）: 実際のcor_sex列は0.4477。本文は|r|=0.45を維持。

## 検証（2026-10-25）

### 13. AQ-1 SubGraph.tsv マッピングバグの確認

StatusNetwork.tsv（128×128有向隣接行列）を全51疾患で解析し、以下を確認：

1. 47/51疾患でSubGraph.tsvは孤立局所最小をスキップして出力（SKIP形式）
2. 4疾患（dementia, delirium, diabetes, sepsis）は128行出力（FULL形式）
3. `build_patient_scores.R` L36の `pattern_id := .I` が行番号をstate IDとして使用
4. 甲状腺のpatient_scores.tsvで誤割当を直接確認（state 55→basin 4誤付与、state 56→basin 1誤付与）
5. 30の複数basin疾患で観察の1.5〜61%が誤ったbasinに分類

**Limitation 8**: 不確定だったVERIFYタグを削除し、確認済みバグとして記述を書き直し。
**Methods §4.4**: SubGraph.tsvのSKIP形式と4疾患のFULL形式を明記。
**AQ-1**: 確認済みバグとして更新。修正コード案と再実行手順を追加。

検証結果: `paper/claude/validation/` に保存
- `verify_subgraph_v2.py`: 検証スクリプト
- `subgraph_reconciliation.tsv`: 全51疾患の照合結果
- `subgraph_verification_v2.json`: 詳細JSON
- `synthetic_example.md`: 合成8状態の例
- `validation_report.md`: 検証レポート全文

### 14. AQ-4 PCA variance 12.5% vs 13.3%の解決

- 固有値ファイル（`ooc_cov/eigenvalues.csv`）: Σ(top7)/Σ(all) = 中央値13.3%
- Figure 2c: 患者代表スコアベースの分散 = 12.5%
- 差異の原因: `which.max(old_row_id)` による患者代表化
- Methods §4.3とLimitation 3に両値を明記

### 15. AQ-5 Figure 7 入院KMの検証

- `km_per_basin_hosp.R` L52-55: first_admission_ymd → time-to-first-admission を確認
- Codex版の「any-admission flag + total observation」は実際のコードと不一致
- 認知症はfull_128形式でSubGraphバグの影響なし
- 原稿修正不要
