#!/usr/bin/env python3
"""
prepare_annotations.py
=======================
为成吉思汗Y染色体系统发育树可视化准备注释文件。

模块用法：
    from python.prepare_annotations import prepare_all
    result = prepare_all(tree_file, meta_file, output_dir)

脚本用法：
    python3 prepare_annotations.py \\
        --tree  <treefile> \\
        --meta  <metadata.tsv> \\
        --output_dir <annotation_output_dir>
"""

import re
import sys
import argparse
from pathlib import Path

import pandas as pd


# ============================================================
# 配色常量（单倍群分支 / 人群分组 / 语系）
# ============================================================

# 主干单倍群 → 分支颜色（用于 ggtree 分支着色）
BRANCH_COLORS: dict[str, str] = {
    "C2a1a3a1":    "#C62828",  # 深红：C2a1a3a1 星状支系（Golden Horde 候选）
    "C2a1a3a2":    "#E65100",  # 深橙：C2a1a3a2 支系
    "C2a1a3a4":    "#BF360C",  # 棕橙：C2a1a3a4 支系
    "C2a1a3a6":    "#7B1FA2",  # 紫色：C2a1a3a6 支系
    "C2a1a3_base": "#0277BD",  # 蓝色：C2a1a3 基础（未分到上游亚支的节点）
    "C2a1a1":      "#2E7D32",  # 深绿：C2a1a1 支系
    "C2a1a2":      "#00695C",  # 青色：C2a1a2 支系
    "C2a1b":       "#558B2F",  # 橄榄绿：C2a1b
    "C2b":         "#827717",  # 暗黄：C2b
    "B":           "#6D4C41",  # 棕色：B 外群
    "D":           "#5D4037",  # 深棕：D
    "E":           "#8D6E63",  # 浅棕：E
    "G":           "#455A64",  # 蓝灰：G
    "H":           "#607D8B",  # 蓝灰：H
    "I":           "#546E7A",  # 蓝灰：I
    "J":           "#1565C0",  # 蓝色：J
    "L":           "#00838F",  # 青蓝：L
    "N":           "#6A1B9A",  # 紫色：N
    "O":           "#2E7D32",  # 绿色：O
    "Q":           "#AD1457",  # 洋红：Q
    "R":           "#F9A825",  # 金黄：R
    "non_C2":      "#9E9E9E",  # 灰色：其他未识别非 C2
}

# 人群分组 → tip 标签颜色（高饱和度高对比色板，色相充分分离）
POP_GROUP_COLORS: dict[str, str] = {
    "mongol_china":           "#E53935",  # 鲜红：中国蒙古族（最核心）
    "central_asia_c2a1a3":    "#8E24AA",  # 鲜紫：中亚 C2a1a3 携带者
    "tungusic_china_c2a1a3":  "#2E9948",  # 鲜绿：中国通古斯族 C2a1a3
    "eastasia_c2a1a3":        "#FB8C00",  # 鲜橙：其他东亚 C2a1a3
    "other_c2a1a3":           "#0097A7",  # 鲜青：其他 C2a1a3
    "other_c2":               "#1E88E5",  # 鲜蓝：其他 C2 支系
    "non_c2":                 "#9E9E9E",  # 中灰：非 C2 外群
    "new_sample_c2a1a3":      "#E91E63",  # 鲜粉：新样本 C2a1a3（BDC/LP前缀）
    "new_sample_other":       "#26C6DA",  # 亮青：新样本非 C2（BDC/LP6，非C2a1a3）
}

# 地区分组 → 子树 ID 颜色。颜色选择偏深，保证 5-6 pt 标签可读。
REGION_GROUP_COLORS: dict[str, str] = {
    "china_north":      "#0072B2",
    "china_south":      "#D55E00",
    "china_unspecified":"#CC79A7",
    "central_asia":     "#009E73",
    "west_asia":        "#8E24AA",
    "south_asia":       "#E69F00",
    "unknown_region":   "#6E6E6E",
}

# 语系 → 条带颜色（用于注释条带，tidyplots/Okabe-Ito 风格离散配色）
LANGUAGE_COLORS: dict[str, str] = {
    "Sinitic":       "#0072B2",
    "Turkic":        "#56B4E9",
    "Indo-European": "#009E73",
    "Tungusic":      "#F0E442",
    "Mongolic":      "#E69F00",
    "Kazakh":        "#D55E00",
    "Khoisan":       "#8ECAE6",
    "Niger-Congo":   "#219EBC",
    "NA":            "#D9D9D9",
}

# 单倍群前缀优先级（从最特异到最宽泛）
DEFAULT_HAPLOGROUP_PRIORITY: list[tuple[str, str]] = [
    ("C2a1a3a1",    "C2a1a3a1"),
    ("C2a1a3a2",    "C2a1a3a2"),
    ("C2a1a3a4",    "C2a1a3a4"),
    ("C2a1a3a6",    "C2a1a3a6"),
    ("C2a1a3",      "C2a1a3_base"),
    ("C2a1a1",      "C2a1a1"),
    ("C2a1a2",      "C2a1a2"),
    ("C2a1b",       "C2a1b"),
    ("C2b",         "C2b"),
]

# BDC 样本按当前项目约定处理为已发表；其他缺失元数据的样本仍保留未知状态。
DEFAULT_PUBLISHED_PREFIXES_WITHOUT_METADATA: tuple[str, ...] = ("BDC",)
DEFAULT_TARGET_CLADE = "C2a1a3"


# ============================================================
# 工具函数
# ============================================================

def _clean(val) -> str:
    """将 pandas NaN / 'nan' / None 统一转为空字符串。"""
    if val is None:
        return ""
    try:
        import math
        if isinstance(val, float) and math.isnan(val):
            return ""
    except TypeError:
        pass
    s = str(val).strip()
    return "" if s.lower() == "nan" else s


def parse_design_file(design_file: str | None) -> dict[str, str]:
    """读取 Markdown 中的 visual_design_config 配置块。"""
    config: dict[str, str] = {
        "target_clade": DEFAULT_TARGET_CLADE,
        "output_prefix": DEFAULT_TARGET_CLADE.lower(),
        "published_prefixes": ",".join(DEFAULT_PUBLISHED_PREFIXES_WITHOUT_METADATA),
        "haplogroup_group_rules": ";".join(
            f"{prefix}={group}" for prefix, group in DEFAULT_HAPLOGROUP_PRIORITY
        ),
    }
    if not design_file:
        return config

    text = Path(design_file).read_text(encoding="utf-8")
    in_block = False
    found_block = False
    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith("```visual_design_config"):
            in_block = True
            found_block = True
            continue
        if in_block and line.startswith("```"):
            break
        if not in_block or not line or line.startswith("#"):
            continue
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        config[key.strip()] = value.strip()

    if not found_block and design_file.lower().endswith(".tsv"):
        df = pd.read_csv(design_file, sep="\t", dtype=str).fillna("")
        if {"key", "value"}.issubset(df.columns):
            for _, row in df.iterrows():
                config[str(row["key"]).strip()] = str(row["value"]).strip()
    return config


def parse_haplogroup_rules(rule_text: str) -> list[tuple[str, str]]:
    rules: list[tuple[str, str]] = []
    for item in (rule_text or "").split(";"):
        item = item.strip()
        if not item or "=" not in item:
            continue
        prefix, group = item.split("=", 1)
        prefix = prefix.strip()
        group = group.strip()
        if prefix and group:
            rules.append((prefix, group))
    return rules or DEFAULT_HAPLOGROUP_PRIORITY


def split_design_list(value: str) -> tuple[str, ...]:
    return tuple(x.strip() for x in (value or "").split(",") if x.strip())


def extract_tip_labels(newick_str: str) -> list[str]:
    """从 Newick 字符串中提取所有 tip 标签（排除内部节点标签）。"""
    all_labels = re.findall(r"([^(),;:\s]+):", newick_str)
    # 内部节点标签为纯数字或 数字/数字（如 92/80、100）
    tips = [lbl for lbl in all_labels
            if not re.match(r"^[\d.]+(/[\d.]+)?$", lbl)]
    return tips


def assign_haplogroup_group(
    haplogroup: str,
    haplogroup_priority: list[tuple[str, str]] | None = None,
) -> str:
    """将单倍群字符串映射到主干分组（用于分支着色）。"""
    if not haplogroup:
        return "non_C2"
    for prefix, group in (haplogroup_priority or DEFAULT_HAPLOGROUP_PRIORITY):
        if haplogroup.startswith(prefix):
            return group
    macro = haplogroup[0].upper()
    if macro in BRANCH_COLORS:
        return macro
    return "non_C2"


def assign_pop_group(
    population: str,
    country: str,
    continent: str,
    haplogroup: str,
    has_metadata: bool,
    target_clade: str = DEFAULT_TARGET_CLADE,
) -> str:
    """
    依据科学问题（成吉思汗支系）对样本进行人群分组。
    分组优先级：蒙古族中国 > 中亚C2a1a3 > 通古斯族中国 > 东亚其他 > 其他C2 > 非C2。
    """
    if not has_metadata:
        is_c2a1a3 = (haplogroup or "").startswith(target_clade)
        return "new_sample_c2a1a3" if is_c2a1a3 else "new_sample_other"

    hg = haplogroup or ""
    is_c2a1a3 = hg.startswith(target_clade)

    if not is_c2a1a3:
        return "other_c2" if hg.startswith("C2") else "non_c2"

    pop  = population or ""
    cnt  = country or ""
    cont = continent or ""

    # 层级1：中国蒙古族（与成吉思汗最直接相关）
    if pop == "Mongol" and cnt == "China":
        return "mongol_china"

    # 层级2：中亚人群 C2a1a3（哈萨克斯坦/吉尔吉斯/乌兹别克等）
    if cont == "Central_Asia":
        return "central_asia_c2a1a3"

    # 层级3：中国通古斯语系民族（满族/达斡尔/锡伯/鄂伦春）
    if pop in ("Man", "Daur", "Xibe", "Oroqen") and cnt == "China":
        return "tungusic_china_c2a1a3"

    # 层级4：其他东亚 C2a1a3（汉族/回族/哈萨克族等）
    if cont == "East_Asia":
        return "eastasia_c2a1a3"

    return "other_c2a1a3"


def assign_region_group(country: str, province: str, continent: str) -> str:
    """按大区给 tip ID 分组，用于子树标签颜色。"""
    cnt = (country or "").strip()
    prov = (province or "").strip().replace(" ", "_")
    cont = (continent or "").strip()

    central_asia_countries = {
        "Kazakhstan", "Kyrgyzstan", "Uzbekistan", "Tajikistan", "Turkmenistan"
    }
    west_asia_countries = {"Afghanistan", "Iran"}
    south_asia_countries = {"Pakistan", "India", "Nepal", "Bangladesh", "Sri_Lanka"}
    china_north_provinces = {
        "Beijing", "Tianjin", "Hebei", "Shanxi", "Inner_Mongolia",
        "InnerMongolia", "Liaoning", "Jilin", "Heilongjiang", "Shandong",
        "Shaanxi", "Gansu", "Xinjiang", "Ningxia", "Qinghai"
    }
    china_south_provinces = {
        "Shanghai", "Jiangsu", "Zhejiang", "Anhui", "Fujian", "Jiangxi",
        "Henan", "Hubei", "Hunan", "Guangdong", "Guangxi", "Hainan",
        "Chongqing", "Sichuan", "Guizhou", "Yunnan"
    }

    if cnt == "China":
        if not prov:
            return "china_unspecified"
        if prov in china_south_provinces:
            return "china_south"
        return "china_north" if prov in china_north_provinces else "china_south"
    if cnt in west_asia_countries:
        return "west_asia"
    if cnt in south_asia_countries:
        return "south_asia"
    if cnt in central_asia_countries:
        return "central_asia"
    if cont == "West_Asia":
        return "west_asia"
    if cont == "South_Asia":
        return "south_asia"
    if cont == "Central_Asia":
        return "central_asia"
    return "unknown_region"


def make_display_label(
    population: str,
    province: str,
    country: str,
    sample_id: str,
    haplogroup: str = "",
) -> str:
    """
    生成 tip 显示标签。
    统一为 population_region_fullHaplogroup。中国样本 region 用省份，
    境外样本 region 用国家；缺元数据样本保留原始 ID 以避免臆造信息。
    """
    pop  = (population or "").strip()
    prov = (province or "").strip().replace(" ", "")
    cnt  = (country or "").strip()
    hg   = (haplogroup or "").strip()

    if pop:
        if cnt == "China":
            region = prov if prov else "China"
        else:
            region = cnt if cnt else prov
        return "_".join(x for x in (pop, region, hg) if x)

    # 无民族信息时退而使用地理信息
    if cnt == "China":
        region = prov if prov else "China"
        return "_".join(x for x in (region, hg) if x)
    if cnt:
        return "_".join(x for x in (cnt, hg) if x)

    return sample_id


def _extract_haplogroup_from_label(tip: str) -> str:
    """从 tip 标签中提取单倍群（无元数据时的后备方法）。"""
    parts = tip.split("_")
    for part in reversed(parts):
        # 单倍群以大写字母+数字开头，或为单字母（B、D、E等）
        if re.match(r"^[A-Z][0-9]", part):
            return part
        if part in ("B", "D", "E", "G") and len(part) == 1:
            return part
        # 带连字符的单倍群（如 C2a1a3-M504）
        if "-" in part and re.match(r"^[A-Z][0-9]", part):
            return part
    return ""


# ============================================================
# 核心逻辑
# ============================================================

def build_annotation_table(
    tip_labels: list[str],
    meta_df: pd.DataFrame,
    design_config: dict[str, str] | None = None,
) -> pd.DataFrame:
    """
    为所有树 tip 构建注释表。
    匹配策略：对每个 tip 标签，遍历元数据 ID（按长度降序），
    找到第一个满足 tip.startswith(ID + '_') 的条目。
    """
    design_config = design_config or {}
    haplogroup_priority = parse_haplogroup_rules(
        design_config.get("haplogroup_group_rules", "")
    )
    published_prefixes = split_design_list(
        design_config.get("published_prefixes", ",".join(DEFAULT_PUBLISHED_PREFIXES_WITHOUT_METADATA))
    )
    target_clade = design_config.get("target_clade", DEFAULT_TARGET_CLADE)

    # 按 ID 长度降序排列，优先匹配最长 ID（避免短 ID 误匹配）
    meta_dict: dict[str, dict] = {
        row["ID"]: row.to_dict() for _, row in meta_df.iterrows()
    }
    sorted_ids = sorted(meta_dict.keys(), key=len, reverse=True)

    records: list[dict] = []

    for tip in tip_labels:
        matched_id: str | None = None
        for mid in sorted_ids:
            if tip.startswith(mid + "_") or tip == mid:
                matched_id = mid
                break

        if matched_id:
            row = meta_dict[matched_id]
            haplogroup    = _clean(row.get("Haplogroup"))
            population    = _clean(row.get("Population"))
            province      = _clean(row.get("Province"))
            country       = _clean(row.get("Country"))
            continent     = _clean(row.get("Continent"))
            lang_detailed = _clean(row.get("Language_detailed"))
            lang_family   = _clean(row.get("Language_family"))
            data_source   = _clean(row.get("Data"))
            data_published = _clean(row.get("data_published"))
            has_meta = True
        else:
            haplogroup    = _extract_haplogroup_from_label(tip)
            population    = ""
            province      = ""
            country       = ""
            continent     = ""
            lang_detailed = "NA"
            lang_family   = "NA"
            data_source   = "unknown"
            data_published = (
                "published"
                if tip.startswith(published_prefixes)
                else "unknown"
            )
            has_meta = False

        # 未发表：无元数据 或 元数据来源为非已发表数据集
        is_new = data_published != "published"

        hg_group      = assign_haplogroup_group(haplogroup, haplogroup_priority)
        pop_group     = assign_pop_group(
            population, country, continent, haplogroup, has_meta, target_clade=target_clade
        )
        region_group  = assign_region_group(country, province, continent)
        display_label = make_display_label(
            population, province, country, tip, haplogroup=haplogroup
        )

        records.append({
            "label":            tip,            # 与树 tip 标签一致
            "sample_id":        matched_id or tip,
            "haplogroup":       haplogroup,
            "haplogroup_group": hg_group,
            "population":       population,
            "province":         province,
            "country":          country,
            "continent":        continent,
            "language_detailed":lang_detailed,
            "language_family":  lang_family,
            "data_source":      data_source,
            "data_published":   data_published,
            "pop_group":        pop_group,
            "region_group":     region_group,
            "display_label":    display_label,
            "is_new_sample":    is_new,
            "has_metadata":     has_meta,
        })

    return pd.DataFrame(records)


def prepare_all(
    tree_file: str,
    meta_file: str,
    output_dir: str,
    design_file: str | None = None,
) -> dict:
    """
    主逻辑：读取树文件和元数据，生成所有注释 TSV 文件。

    Returns
    -------
    dict : 各输出文件路径及汇总统计
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    design_config = parse_design_file(design_file)

    # 读取树文件
    with open(tree_file, encoding="utf-8") as fh:
        newick_str = fh.read()
    tip_labels = extract_tip_labels(newick_str)

    # 读取元数据
    meta_df = pd.read_csv(meta_file, sep="\t", dtype=str)

    # 构建注释表
    annot_df = build_annotation_table(tip_labels, meta_df, design_config=design_config)

    # ── 输出：主注释表 ──────────────────────────────────────────────
    main_out = out / "tip_annotations.tsv"
    annot_df.to_csv(main_out, sep="\t", index=False)

    # ── 输出：颜色映射表 ────────────────────────────────────────────
    branch_color_df = pd.DataFrame(
        [{"haplogroup_group": k, "branch_color": v} for k, v in BRANCH_COLORS.items()]
    )
    branch_color_out = out / "branch_colors.tsv"
    branch_color_df.to_csv(branch_color_out, sep="\t", index=False)

    pop_color_df = pd.DataFrame(
        [{"pop_group": k, "pop_color": v} for k, v in POP_GROUP_COLORS.items()]
    )
    pop_color_out = out / "pop_group_colors.tsv"
    pop_color_df.to_csv(pop_color_out, sep="\t", index=False)

    lang_color_df = pd.DataFrame(
        [{"language_detailed": k, "lang_color": v} for k, v in LANGUAGE_COLORS.items()]
    )
    lang_color_out = out / "language_colors.tsv"
    lang_color_df.to_csv(lang_color_out, sep="\t", index=False)

    region_color_df = pd.DataFrame(
        [{"region_group": k, "region_color": v} for k, v in REGION_GROUP_COLORS.items()]
    )
    region_color_out = out / "region_group_colors.tsv"
    region_color_df.to_csv(region_color_out, sep="\t", index=False)

    design_config_out = out / "design_config.tsv"
    pd.DataFrame(
        [{"key": k, "value": v} for k, v in sorted(design_config.items())]
    ).to_csv(design_config_out, sep="\t", index=False)

    if design_file:
        design_source_out = out / "visual_design.md"
        design_source_out.write_text(Path(design_file).read_text(encoding="utf-8"), encoding="utf-8")

    # ── 汇总统计 ────────────────────────────────────────────────────
    summary = {
        "total_tips":       len(tip_labels),
        "matched_to_meta":  int(annot_df["has_metadata"].sum()),
        "new_samples":      int(annot_df["is_new_sample"].sum()),
        "target_clade":     design_config.get("target_clade", DEFAULT_TARGET_CLADE),
        "target_count":     int(
            annot_df["haplogroup"].fillna("").str.startswith(
                design_config.get("target_clade", DEFAULT_TARGET_CLADE)
            ).sum()
        ),
    }

    return {
        "annotations":    str(main_out),
        "branch_colors":  str(branch_color_out),
        "pop_colors":     str(pop_color_out),
        "lang_colors":    str(lang_color_out),
        "region_colors":  str(region_color_out),
        "design_config":  str(design_config_out),
        "summary":        summary,
    }


# ============================================================
# 脚本入口
# ============================================================

def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="为系统发育树可视化准备注释文件（成吉思汗Y染色体项目）。"
    )
    parser.add_argument(
        "--tree",
        required=True,
        help="输入树文件路径（Newick/IQTree .treefile 格式）",
    )
    parser.add_argument(
        "--meta",
        required=True,
        help="样本元数据 TSV 文件路径（含 ID、Haplogroup、Population 等列）",
    )
    parser.add_argument(
        "--design",
        required=False,
        help="可视化设计 Markdown 文件路径（含 visual_design_config 配置块）",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="注释文件输出目录",
    )
    return parser.parse_args(argv)


def main(args: argparse.Namespace) -> int:
    result = prepare_all(
        tree_file=args.tree,
        meta_file=args.meta,
        output_dir=args.output_dir,
        design_file=args.design,
    )
    print(f"[INFO] 注释文件已写入: {args.output_dir}")
    s = result["summary"]
    print(
        f"[INFO] 统计: 总 tip={s['total_tips']}, "
        f"匹配元数据={s['matched_to_meta']}, "
        f"新样本={s['new_samples']}, "
        f"{s['target_clade']} 总计={s['target_count']}"
    )
    return 0


if __name__ == "__main__":
    args = _parse_args()
    sys.exit(main(args))
