#!/usr/bin/env Rscript
# visualize_tree.R
# ================
# 成吉思汗Y染色体系统发育树 ggtree 可视化脚本
#
# 用法（由 pipe/run_chengjisihan.sh 调用）：
#   Rscript src/visualize_tree.R \
#     --tree   <treefile> \
#     --annot  <tip_annotations.tsv> \
#     --branch_colors <branch_colors.tsv> \
#     --pop_colors    <pop_group_colors.tsv> \
#     --region_colors <region_group_colors.tsv> \
#     --lang_colors   <language_colors.tsv> \
#     --out_full_tree <path> \
#     --out_subtree   <path> \
#     --width_mm  200 \
#     --height_mm 290 \
#     --dpi       300

# ============================================================
# 工具函数
# ============================================================

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a[[1]])) a else b

pt2size <- function(pt) pt / ggplot2::.pt

# 根据画幅高度和 tip 数量自动计算标签字号（pt）
calc_font_pt <- function(fig_h_mm, n_tips, scale = 0.82, lo = 4.0, hi = 8.0) {
    pt_per_row <- (fig_h_mm / n_tips) * (72 / 25.4)
    min(hi, max(lo, pt_per_row * scale))
}

title_case_label <- function(x) {
    x <- gsub("_", " ", x)
    x <- gsub("\\b([a-z])", "\\U\\1", x, perl = TRUE)
    x
}

macro_haplogroup <- function(haplogroup) {
    hg <- ifelse(is.na(haplogroup), "", haplogroup)
    sub("^([A-Z]).*$", "\\1", hg)
}

parse_cli_args <- function(argv) {
    result <- list()
    i <- 1L
    while (i <= length(argv)) {
        if (startsWith(argv[i], "--")) {
            key <- sub("^--", "", argv[i])
            if (i + 1L <= length(argv) && !startsWith(argv[i + 1L], "--")) {
                result[[key]] <- argv[i + 1L]; i <- i + 2L
            } else {
                result[[key]] <- TRUE; i <- i + 1L
            }
        } else { i <- i + 1L }
    }
    result
}

safe_mrca <- function(tree_obj, tips) {
    in_tree <- tips[tips %in% tree_obj$tip.label]
    if (length(in_tree) >= 2) getMRCA(tree_obj, in_tree)
    else if (length(in_tree) == 1) which(tree_obj$tip.label == in_tree[1])
    else NA_integer_
}

read_design_config <- function(path) {
    defaults <- list(
        target_clade = "C2a1a3",
        output_prefix = "c2a1a3",
        full_tree_clade_labels = "B,C2,C2a1a3,C2a1a1,C2a1a2,C2a1b,C2b",
        subtree_clade_labels = "C2a1a3,C2a1a3a1,C2a1a3a2,C2a1a3a6",
        subtree_fine_clade_labels = "C2a1a3a1a1,C2a1a3a1a2,C2a1a3a6b"
    )
    if (is.null(path) || is.na(path) || !file.exists(path)) return(defaults)
    df <- readr::read_tsv(path, show_col_types = FALSE, col_types = "cc")
    if (!all(c("key", "value") %in% names(df))) return(defaults)
    cfg <- as.list(stats::setNames(df$value, df$key))
    modifyList(defaults, cfg)
}

split_design_list <- function(value) {
    out <- trimws(unlist(strsplit(value %||% "", ",", fixed = TRUE)))
    out[nzchar(out)]
}

safe_prefix <- function(value) {
    gsub("[^A-Za-z0-9_\\-]+", "_", tolower(value))
}

# ============================================================
# 依赖加载
# ============================================================
suppressPackageStartupMessages({
    library(ape)
    library(ggtree)
    library(treeio)
    library(tidytree)
    library(ggplot2)
    library(dplyr)
    library(readr)
    library(tibble)
    library(patchwork)
    library(ggnewscale)
})

# ============================================================
# 参数解析
# ============================================================
ARGS <- parse_cli_args(commandArgs(trailingOnly = TRUE))

TREE_FILE       <- ARGS[["tree"]]
ANNOT_FILE      <- ARGS[["annot"]]
BRANCH_COLORS_F <- ARGS[["branch_colors"]]
POP_COLORS_F    <- ARGS[["pop_colors"]]
REGION_COLORS_F <- ARGS[["region_colors"]]
LANG_COLORS_F   <- ARGS[["lang_colors"]]
DESIGN_CONFIG_F <- ARGS[["design_config"]]
OUT_FULL_TREE   <- ARGS[["out_full_tree"]]
OUT_SUBTREE     <- ARGS[["out_subtree"]]
FIG_W_MM        <- as.numeric(ARGS[["width_mm"]]  %||% "200")
FIG_H_MM        <- as.numeric(ARGS[["height_mm"]] %||% "290")
DPI             <- as.numeric(ARGS[["dpi"]]        %||% "300")

dir.create(OUT_FULL_TREE, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_SUBTREE,   recursive = TRUE, showWarnings = FALSE)

DESIGN_CONFIG <- read_design_config(DESIGN_CONFIG_F)
TARGET_CLADE <- DESIGN_CONFIG[["target_clade"]]
OUTPUT_PREFIX <- safe_prefix(DESIGN_CONFIG[["output_prefix"]] %||% TARGET_CLADE)
FULL_TREE_CLADE_LABELS <- split_design_list(DESIGN_CONFIG[["full_tree_clade_labels"]])
SUBTREE_CLADE_LABELS <- split_design_list(DESIGN_CONFIG[["subtree_clade_labels"]])
SUBTREE_FINE_CLADE_LABELS <- split_design_list(DESIGN_CONFIG[["subtree_fine_clade_labels"]])

# ============================================================
# 数据加载
# ============================================================
message("[INFO] 读取树文件: ", TREE_FILE)
tree <- read.tree(TREE_FILE)

message("[INFO] 读取注释文件: ", ANNOT_FILE)
annot <- read_tsv(ANNOT_FILE, show_col_types = FALSE) |>
    mutate(
        is_new_sample = as.logical(is_new_sample),
        has_metadata  = as.logical(has_metadata),
        language_detailed = ifelse(
            is.na(language_detailed) | language_detailed == "NA",
            NA_character_, language_detailed
        )
    )

# 颜色映射
branch_color_df <- read_tsv(BRANCH_COLORS_F, show_col_types = FALSE)
branch_colors   <- setNames(branch_color_df$branch_color, branch_color_df$haplogroup_group)

pop_color_df <- read_tsv(POP_COLORS_F, show_col_types = FALSE)
pop_colors   <- setNames(pop_color_df$pop_color, pop_color_df$pop_group)

region_color_df <- read_tsv(REGION_COLORS_F, show_col_types = FALSE)
region_colors   <- setNames(region_color_df$region_color, region_color_df$region_group)
region_labels <- c(
    china_north = "Northern China",
    china_south = "Southern China",
    china_unspecified = "China",
    central_asia = "Central Asia",
    west_asia = "West Asia",
    south_asia = "South Asia",
    unknown_region = "Unknown Region"
)

lang_color_df <- read_tsv(LANG_COLORS_F, show_col_types = FALSE)
lang_colors   <- setNames(lang_color_df$lang_color, lang_color_df$language_detailed)
lang_colors   <- lang_colors[!is.na(names(lang_colors))]

all_colors <- c(branch_colors, pop_colors)
color_breaks <- names(all_colors)
color_labels <- setNames(title_case_label(color_breaks), color_breaks)
custom_color_labels <- c(
    C2a1a3a1 = "C2a1a3a1",
    C2a1a3a2 = "C2a1a3a2",
    C2a1a3a4 = "C2a1a3a4",
    C2a1a3a6 = "C2a1a3a6",
    C2a1a3_base = "C2a1a3 Base",
    C2a1a1 = "C2a1a1",
    C2a1a2 = "C2a1a2",
    C2a1b = "C2a1b",
    C2b = "C2b",
    non_C2 = "Non-C2",
    mongol_china = "Mongol China",
    central_asia_c2a1a3 = "Central Asia C2a1a3",
    tungusic_china_c2a1a3 = "Tungusic China C2a1a3",
    eastasia_c2a1a3 = "East Asia C2a1a3",
    other_c2a1a3 = "Other C2a1a3",
    other_c2 = "Other C2",
    non_c2 = "Non-C2",
    new_sample_c2a1a3 = "New Sample C2a1a3",
    new_sample_other = "New Sample Other"
)
color_labels[names(custom_color_labels)] <- custom_color_labels

# ============================================================
# 树预处理：以 B 单倍群为外群重定根
# ============================================================
b_tips <- grep("(Ju_hoan|Mbuti|LP6005|_B[0-9]|_B$)",
               tree$tip.label, value = TRUE, perl = TRUE)
message("[INFO] B 外群 tip 数: ", length(b_tips))
if (length(b_tips) >= 2) {
    tree <- root(tree, node = getMRCA(tree, b_tips), resolve.root = TRUE)
    message("[INFO] 已以 B MRCA 重定根")
} else if (length(b_tips) == 1) {
    tree <- root(tree, outgroup = b_tips[1], resolve.root = TRUE)
} else {
    message("[WARN] 未找到 B 外群，保持原始定根")
}

# ============================================================
# 分支颜色分组（groupOTU）
# ============================================================
temp_grp <- annot |>
    filter(label %in% tree$tip.label) |>
    group_by(haplogroup_group) |>
    summarise(labels = list(label), .groups = "drop")
group_list   <- setNames(temp_grp$labels, temp_grp$haplogroup_group)
tree_grouped <- groupOTU(tree, group_list, group_name = "haplo_group")

# ============================================================
# 全树节点坐标 & 宏单倍群 MRCA
# ============================================================
message("[INFO] 计算节点坐标...")
p_coord_tmp <- ggtree(tree_grouped, layout = "rectangular")
tree_pos_df <- p_coord_tmp$data
tree_max_x  <- max(tree_pos_df$x[tree_pos_df$isTip], na.rm = TRUE)
message("[INFO] tree_max_x = ", signif(tree_max_x, 5))

target_all_tips <- annot |>
    filter(startsWith(haplogroup, TARGET_CLADE), label %in% tree$tip.label) |>
    pull(label)

tip_y_df <- tree_pos_df |>
    filter(isTip) |>
    select(label, y)

tips_for_clade_label <- function(clade_label, annot_df, tree_obj, b_tip_vec) {
    if (clade_label == "B") return(b_tip_vec)
    annot_df |>
        filter(label %in% tree_obj$tip.label, startsWith(haplogroup, clade_label)) |>
        pull(label)
}

macro_grp <- tibble(macro = FULL_TREE_CLADE_LABELS) |>
    rowwise() |>
    mutate(
        labels = list(tips_for_clade_label(macro, annot, tree, b_tips)),
        n = length(unlist(labels))
    ) |>
    ungroup()

macro_labels_df <- macro_grp |>
    rowwise() |>
    mutate(y = median(tip_y_df$y[tip_y_df$label %in% unlist(labels)], na.rm = TRUE)) |>
    ungroup() |>
    filter(!is.na(y), n >= 2) |>
    mutate(
        clade_label = paste0(macro, " (", n, ")"),
        x = tree_max_x * 1.24
    ) |>
    arrange(y)
message("[INFO] 宏单倍群标注节点数: ", nrow(macro_labels_df))

# ============================================================
# 字号：根据画幅自动计算
# ============================================================
n_full     <- length(tree$tip.label)
font_full  <- calc_font_pt(FIG_H_MM, n_full,  scale = 0.62, lo = 2.8, hi = 3.8)
message("[INFO] 全树字号: ", round(font_full, 2), " pt (", n_full, " tips)")

# gheatmap 偏移（让条带紧贴标签右侧，减少留白）
offset_full <- tree_max_x * 0.22

# ============================================================
# 图 1：全树概览（全部 ID，新样本彩色、已知样本灰色）
# ============================================================
message("[INFO] 绘制全树概览...")

p_full <- ggtree(
    tree_grouped,
    aes(color = haplo_group),
    layout = "rectangular",
    lwd    = 0.28
) %<+% annot +
    # 已知样本：小灰点
    geom_tippoint(
        data  = function(d) filter(d, isTip & !(is_new_sample %in% TRUE)),
        color = "grey65", shape = 16, size = 0.22, na.rm = TRUE
    ) +
    # 新样本：彩色三角（醒目）
    geom_tippoint(
        data  = function(d) filter(d, isTip & (is_new_sample %in% TRUE)),
        aes(color = pop_group),
        shape = 17, size = 0.85, na.rm = TRUE
    ) +
    # 已知样本 ID：深灰色小字
    geom_tiplab(
        data     = function(d) filter(d, isTip & !(is_new_sample %in% TRUE)),
        aes(label = display_label),
        color    = "grey40",
        align    = FALSE,
        size     = pt2size(font_full),
        family = "Arial",
        hjust    = -0.06,
        na.rm    = TRUE
    ) +
    # 新样本 ID：彩色加粗，略大
    geom_tiplab(
        data     = function(d) filter(d, isTip & (is_new_sample %in% TRUE)),
        aes(label = display_label, color = pop_group),
        align    = FALSE,
        size     = pt2size(font_full + 0.2),
        family = "Arial",
        fontface = "bold",
        hjust    = -0.06,
        na.rm    = TRUE
    ) +
    # 宏单倍群标注放在树右侧空白区，避免遮挡分支与 tip。
    geom_label(
        data          = macro_labels_df,
        aes(x = x, y = y, label = clade_label),
        inherit.aes   = FALSE,
        size          = pt2size(4),
        family = "Arial",
        fontface      = "bold",
        fill          = alpha("white", 0.88),
        color         = "#222222",
        linewidth     = 0.3,
        label.padding = unit(0.12, "lines"),
        na.rm         = TRUE
    ) +
    scale_color_manual(
        values   = all_colors,
        breaks   = color_breaks,
        labels   = color_labels,
        name     = "Haplogroup / Population",
        na.value = "#999999"
    ) +
    hexpand(0.48) +
    theme_tree(plot.margin = margin(3, 4, 3, 3, "mm")) +
    theme(
        text            = element_text(family = "Arial"),
        legend.text     = element_text(family = "Arial", size = 5),
        legend.title    = element_text(family = "Arial", size = 6, face = "bold"),
        legend.position = c(0.80, 0.18),
        legend.box      = "vertical",
        legend.background = element_rect(fill = alpha("white", 0.78), color = "grey80"),
        legend.key.size = unit(2.4, "mm")
    )

# 语系条带
lang_mat <- annot |>
    filter(label %in% tree$tip.label) |>
    select(label, language_detailed) |>
    column_to_rownames("label")

suppressWarnings({
    p_full_final <- gheatmap(
        p        = p_full,
        data     = lang_mat,
        width    = 0.04,
        offset   = offset_full,
        colnames = FALSE,
        color    = NA
    ) +
        scale_fill_manual(
            values   = lang_colors,
            name     = "Language family",
            na.value = "#CCCCCC"
        ) +
        theme(
            legend.text       = element_text(family = "Arial", size = 5),
            legend.title      = element_text(family = "Arial", size = 6, face = "bold"),
            legend.position   = c(0.80, 0.18),
            legend.box        = "vertical",
            legend.background = element_rect(fill = alpha("white", 0.78), color = "grey80"),
            legend.key.size   = unit(2.4, "mm")
        )
})

message("[INFO] 保存全树图...")
ggsave(file.path(OUT_FULL_TREE, "full_tree_overview.png"),
       p_full_final, width = FIG_W_MM, height = FIG_H_MM, units = "mm", dpi = DPI)
ggsave(file.path(OUT_FULL_TREE, "full_tree_overview.pdf"),
       p_full_final, width = FIG_W_MM, height = FIG_H_MM, units = "mm", device = cairo_pdf)
message("[OK] 全树图已保存")

# ============================================================
# 图 2：目标子树（全部 ID + 详细亚支节点标注）
# ============================================================
message("[INFO] 提取 ", TARGET_CLADE, " 子树...")

mrca_target <- getMRCA(tree, target_all_tips[target_all_tips %in% tree$tip.label])
sub_tree    <- tree_subset(tree, node = mrca_target, levels_back = 0)
message("[INFO] ", TARGET_CLADE, " 子树 tip 数: ", length(sub_tree$tip.label))

sub_annot <- annot |> filter(label %in% sub_tree$tip.label)

# 子树分支分组
temp_sub     <- sub_annot |> group_by(haplogroup_group) |>
    summarise(labels = list(label), .groups = "drop")
sub_grp_list <- setNames(temp_sub$labels, temp_sub$haplogroup_group)
sub_tree_grp <- groupOTU(sub_tree, sub_grp_list, group_name = "haplo_group")

# 子树坐标
p_sub_coord  <- ggtree(sub_tree_grp, layout = "rectangular")
sub_pos_df   <- p_sub_coord$data
sub_max_x    <- max(sub_pos_df$x[sub_pos_df$isTip], na.rm = TRUE)
message("[INFO] sub_max_x = ", signif(sub_max_x, 5))

# 子树详细亚支 MRCA
sub_detail <- setNames(
    lapply(SUBTREE_CLADE_LABELS, function(clade) {
        sub_annot |> filter(startsWith(haplogroup, clade)) |> pull(label)
    }),
    SUBTREE_CLADE_LABELS
)
sub_fine <- setNames(
    lapply(SUBTREE_FINE_CLADE_LABELS, function(clade) {
        sub_annot |> filter(startsWith(haplogroup, clade)) |> pull(label)
    }),
    SUBTREE_FINE_CLADE_LABELS
)
for (nm in names(sub_fine)) {
    if (length(sub_fine[[nm]]) >= 2) sub_detail[[nm]] <- sub_fine[[nm]]
}

sub_node_ids <- vapply(sub_detail, function(t) safe_mrca(sub_tree, t), integer(1))
sub_labels_df <- tibble(node = sub_node_ids, clade_label = names(sub_node_ids)) |>
    filter(!is.na(node)) |>
    left_join(sub_pos_df |> select(node, x, y), by = "node") |>
    filter(!is.na(x), !is.na(y)) |>
    mutate(
        label_x = x - sub_max_x * 0.23,
        label_y = y
    )
message("[INFO] 子树亚支标注节点数: ", nrow(sub_labels_df))

# 字号
n_sub        <- length(sub_tree$tip.label)
font_sub     <- calc_font_pt(FIG_H_MM, n_sub, scale = 0.82, lo = 4.5, hi = 7.0)
font_sub     <- min(6.0, max(5.0, font_sub))
message("[INFO] 子树字号: ", round(font_sub, 2), " pt (", n_sub, " tips)")

offset_sub <- sub_max_x * 0.62

p_sub <- ggtree(
    sub_tree_grp,
    aes(color = haplo_group),
    layout = "rectangular",
    lwd    = 0.35
) %<+% sub_annot +
    scale_color_manual(
        values   = branch_colors,
        breaks   = names(branch_colors),
        labels   = color_labels[names(branch_colors)],
        name     = "Haplogroup Branch",
        na.value = "#999999"
    ) +
    new_scale_color() +
    geom_tippoint(
        data  = function(d) filter(d, isTip & !(is_new_sample %in% TRUE)),
        aes(color = region_group),
        shape = 16, size = 0.35, na.rm = TRUE
    ) +
    geom_tippoint(
        data  = function(d) filter(d, isTip & (is_new_sample %in% TRUE)),
        aes(color = region_group),
        shape = 17, size = 1.1, na.rm = TRUE
    ) +
    geom_tiplab(
        data     = function(d) filter(d, isTip & !(is_new_sample %in% TRUE)),
        aes(label = display_label, color = region_group),
        align    = FALSE,
        size     = pt2size(font_sub),
        family = "Arial",
        hjust    = -0.06,
        na.rm    = TRUE
    ) +
    geom_tiplab(
        data     = function(d) filter(d, isTip & (is_new_sample %in% TRUE)),
        aes(label = display_label, color = region_group),
        align    = FALSE,
        size     = pt2size(font_sub),
        family = "Arial",
        fontface = "bold",
        hjust    = -0.06,
        na.rm    = TRUE
    ) +
    geom_segment(
        data          = sub_labels_df,
        aes(x = label_x, y = label_y, xend = x, yend = y),
        inherit.aes   = FALSE,
        color         = "black",
        linewidth     = 0.18,
        arrow         = arrow(length = unit(0.9, "mm"), type = "closed"),
        na.rm         = TRUE
    ) +
    geom_label(
        data          = sub_labels_df,
        aes(x = label_x, y = label_y, label = clade_label),
        inherit.aes   = FALSE,
        size          = pt2size(4.5),
        family = "Arial",
        fontface      = "bold",
        fill          = alpha("white", 0.88),
        color         = "#222222",
        linewidth     = 0.3,
        label.padding = unit(0.12, "lines"),
        na.rm         = TRUE
    ) +
    scale_color_manual(
        values   = region_colors,
        breaks   = names(region_colors),
        labels   = region_labels[names(region_colors)],
        name     = "Region Group",
        na.value = "#6E6E6E"
    ) +
    xlim(-sub_max_x * 1.35, NA) +
    hexpand(0.36) +
    theme_tree(plot.margin = margin(3, 4, 3, 3, "mm")) +
    theme(
        text            = element_text(family = "Arial"),
        legend.text       = element_text(family = "Arial", size = 5),
        legend.title      = element_text(family = "Arial", size = 6, face = "bold"),
        legend.position   = c(0.20, 0.16),
        legend.box        = "vertical",
        legend.background = element_rect(fill = alpha("white", 0.78), color = "grey80"),
        legend.key.size   = unit(2.4, "mm")
    )

sub_strip_mat <- sub_annot |>
    mutate(
        District = region_group,
        spacer_strip = "strip_space",
        Clade = haplogroup_group
    ) |>
    select(label, District, spacer_strip, Clade) |>
    column_to_rownames("label")
colnames(sub_strip_mat) <- c("District", " ", "Clade")

sub_strip_colors <- c(region_colors, strip_space = "#FFFFFF", branch_colors)
sub_strip_labels <- c(
    region_labels[names(region_colors)],
    strip_space = "Spacer",
    color_labels[names(branch_colors)]
)

suppressWarnings({
    p_sub_final <- gheatmap(
        p        = p_sub,
        data     = sub_strip_mat,
        width    = 0.10,
        offset   = offset_sub,
        colnames = TRUE,
        colnames_angle = 90,
        colnames_position = "top",
        font.size = 2.2,
        color    = "#FFFFFF"
    ) +
        scale_fill_manual(
            values   = sub_strip_colors,
            breaks   = names(sub_strip_colors),
            labels   = sub_strip_labels,
            name     = "Annotation Strip",
            na.value = "#CCCCCC",
            guide    = "none"
        ) +
        theme(
            legend.text       = element_text(family = "Arial", size = 5),
            legend.title      = element_text(family = "Arial", size = 6, face = "bold"),
            legend.position   = c(0.20, 0.16),
            legend.box        = "vertical",
            legend.background = element_rect(fill = alpha("white", 0.78), color = "grey80"),
            legend.key.size   = unit(2.4, "mm")
        )
})

message("[INFO] 保存子树图...")
ggsave(file.path(OUT_SUBTREE, paste0(OUTPUT_PREFIX, "_subtree.png")),
       p_sub_final, width = FIG_W_MM, height = FIG_H_MM, units = "mm", dpi = DPI)
ggsave(file.path(OUT_SUBTREE, paste0(OUTPUT_PREFIX, "_subtree.pdf")),
       p_sub_final, width = FIG_W_MM, height = FIG_H_MM, units = "mm", device = cairo_pdf)
message("[OK] 子树图已保存")

# ============================================================
# 图 3：颜色图例汇总
# ============================================================
message("[INFO] 生成颜色图例...")

make_color_legend_plot <- function(color_vec, title_str) {
    df <- data.frame(
        grp   = factor(names(color_vec), levels = names(color_vec)),
        color = unname(color_vec),
        y     = seq_along(color_vec),
        stringsAsFactors = FALSE
    )
    ggplot(df, aes(x = 0, y = y, fill = grp)) +
        geom_tile(width = 0.5, height = 0.8) +
        geom_text(aes(x = 0.4, label = grp),
                  hjust = 0, size = pt2size(7), family = "Arial") +
        scale_fill_manual(values = setNames(df$color, df$grp), guide = "none") +
        xlim(-0.5, 5) +
        labs(title = title_str, x = NULL, y = NULL) +
        theme_void() +
        theme(
            plot.title  = element_text(family = "Arial", size = 8, face = "bold"),
            plot.margin = margin(2, 2, 2, 2, "mm")
        )
}

branch_legend_colors <- branch_colors[names(branch_colors) != "0"]

p_leg_branch <- make_color_legend_plot(branch_legend_colors, "Haplogroup Branch")
p_leg_pop    <- make_color_legend_plot(pop_colors,           "Population Group")
p_leg_lang   <- make_color_legend_plot(lang_colors,          "Language")

p_legends <- p_leg_branch + p_leg_pop + p_leg_lang +
    plot_layout(ncol = 3) +
    plot_annotation(title = "Color Reference")

ggsave(file.path(OUT_FULL_TREE, "color_legend.png"),
       p_legends, width = FIG_W_MM, height = 90, units = "mm", dpi = DPI)
ggsave(file.path(OUT_FULL_TREE, "color_legend.pdf"),
       p_legends, width = FIG_W_MM, height = 90, units = "mm", device = cairo_pdf)

message("[INFO] 全部可视化完成。")
message("  全树图目录: ", OUT_FULL_TREE)
message("  子树图目录: ", OUT_SUBTREE)
