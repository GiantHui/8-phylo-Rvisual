#!/usr/bin/env bash
# run_chengjisihan.sh
# ====================
# 成吉思汗 Y 染色体树可视化流水线主控脚本
# 用法：bash pipe/run_chengjisihan.sh
# 需在 /mnt/d/Y-成吉思汗项目/2-MLtree/8-ggtree-visual 目录下运行

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

# ── 加载配置 ──────────────────────────────────────────────────────────
CONFIG_FILE="$BASE_DIR/conf/Config.yaml"
source "$BASE_DIR/script/load_config.sh" "$CONFIG_FILE"

# 便捷变量（来自配置文件）
PYTHON_BIN="$TOOLS__PYTHON_BIN"
RSCRIPT_BIN="$TOOLS__RSCRIPT_BIN"

TREE_FILE="$INPUTS__TREE_FILE"
META_FILE="$INPUTS__META_FILE"
DESIGN_FILE="$INPUTS__DESIGN_FILE"
[[ "$DESIGN_FILE" = /* ]] || DESIGN_FILE="$BASE_DIR/$DESIGN_FILE"

OUTPUT_ANNOT="$BASE_DIR/$PATHS__OUTPUT_ANNOTATION"
OUTPUT_FULL_TREE="$BASE_DIR/$PATHS__OUTPUT_FULL_TREE"
OUTPUT_SUBTREE="$BASE_DIR/$PATHS__OUTPUT_SUBTREE"
REPORT_DIR="$BASE_DIR/$PATHS__REPORT_DIR"

PYTHON_MODULE="$BASE_DIR/$PATHS__PYTHON_DIR/prepare_annotations.py"
R_MODULE="$BASE_DIR/$PATHS__SRC_DIR/visualize_tree.R"

FIG_W="$VISUALIZATION__FIGURE_WIDTH_MM"
FIG_H="$VISUALIZATION__FIGURE_HEIGHT_MM"
DPI="$VISUALIZATION__DPI"
INTERACTIVE="${RUNTIME__INTERACTIVE:-false}"

# ── 日志函数 ──────────────────────────────────────────────────────────
log_info()  { echo "[$(date '+%H:%M:%S')] [INFO]  $*"; }
log_ok()    { echo "[$(date '+%H:%M:%S')] [OK]    $*"; }
log_error() { echo "[$(date '+%H:%M:%S')] [ERROR] $*" >&2; }

# ── 输入检查 ──────────────────────────────────────────────────────────
log_info "检查输入文件..."
[[ -f "$TREE_FILE" ]]     || { log_error "树文件不存在: $TREE_FILE"; exit 1; }
[[ -f "$META_FILE" ]]     || { log_error "元数据文件不存在: $META_FILE"; exit 1; }
[[ -f "$DESIGN_FILE" ]]   || { log_error "设计文件不存在: $DESIGN_FILE"; exit 1; }
[[ -f "$PYTHON_MODULE" ]] || { log_error "Python 模块不存在: $PYTHON_MODULE"; exit 1; }
[[ -f "$R_MODULE" ]]      || { log_error "R 模块不存在: $R_MODULE"; exit 1; }

# ── 建立输出目录 ──────────────────────────────────────────────────────
mkdir -p "$OUTPUT_ANNOT" "$OUTPUT_FULL_TREE" "$OUTPUT_SUBTREE" "$REPORT_DIR"

# ================================================================
# 步骤 1：Python 注释准备
# ================================================================
log_info "步骤 1/2：运行 Python 注释准备..."

"$PYTHON_BIN" "$PYTHON_MODULE" \
    --tree       "$TREE_FILE" \
    --meta       "$META_FILE" \
    --design     "$DESIGN_FILE" \
    --output_dir "$OUTPUT_ANNOT"

log_ok "注释文件已生成: $OUTPUT_ANNOT"

# ================================================================
# 步骤 2：R 可视化
# ================================================================
log_info "步骤 2/2：运行 R 树可视化..."

"$RSCRIPT_BIN" "$R_MODULE" \
    --tree         "$TREE_FILE" \
    --annot        "$OUTPUT_ANNOT/tip_annotations.tsv" \
    --branch_colors "$OUTPUT_ANNOT/branch_colors.tsv" \
    --pop_colors    "$OUTPUT_ANNOT/pop_group_colors.tsv" \
    --region_colors "$OUTPUT_ANNOT/region_group_colors.tsv" \
    --lang_colors   "$OUTPUT_ANNOT/language_colors.tsv" \
    --design_config "$OUTPUT_ANNOT/design_config.tsv" \
    --out_full_tree "$OUTPUT_FULL_TREE" \
    --out_subtree   "$OUTPUT_SUBTREE" \
    --width_mm      "$FIG_W" \
    --height_mm     "$FIG_H" \
    --dpi           "$DPI" \
    --interactive   "$INTERACTIVE"

log_ok "可视化完成"

# ================================================================
# 步骤 3：汇总输出清单
# ================================================================
log_info "输出文件清单:"
find "$OUTPUT_ANNOT" "$OUTPUT_FULL_TREE" "$OUTPUT_SUBTREE" \
    -type f | sort | while read -r f; do
    echo "  $f"
done

log_ok "流水线全部完成。"
