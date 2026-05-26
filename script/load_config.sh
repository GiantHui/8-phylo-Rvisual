#!/usr/bin/env bash
# load_config.sh — 将 conf/Config.yaml 中的 section.key 解析为 shell 变量
# 使用方法：source script/load_config.sh conf/Config.yaml
# 支持：顶层 section + 二级 key 的简单 YAML（标量值，两空格缩进）

set -euo pipefail

CONFIG_FILE="${1:?需要提供配置文件路径，例如 conf/Config.yaml}"

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "[ERROR] 配置文件不存在: $CONFIG_FILE" >&2
    exit 1
fi

current_section=""

while IFS= read -r line; do
    # 跳过空行和注释行
    [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue

    # 去除行尾注释
    line="${line%%#*}"
    # 去除尾部空白
    line="${line%"${line##*[! ]}"}"

    # 检测顶层 section（不以空格开头，以冒号结尾或冒号后有空格）
    if [[ "$line" =~ ^([A-Za-z_][A-Za-z0-9_]*)[[:space:]]*:[[:space:]]*$ ]]; then
        current_section="${BASH_REMATCH[1]}"
        continue
    fi

    # 检测二级 key: value（以两个空格开头）
    if [[ -n "$current_section" && "$line" =~ ^[[:space:]]{2}([A-Za-z_][A-Za-z0-9_]*)[[:space:]]*:[[:space:]]*(.+)[[:space:]]*$ ]]; then
        key="${BASH_REMATCH[1]}"
        value="${BASH_REMATCH[2]}"
        # 去除引号
        value="${value#\"}"
        value="${value%\"}"
        value="${value#\'}"
        value="${value%\'}"
        # 导出为 SECTION__KEY 格式（全大写）
        var_name="${current_section^^}__${key^^}"
        export "$var_name"="$value"
    fi
done < "$CONFIG_FILE"
