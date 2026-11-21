"""SDF 文本距离场核心算法。

统一契约（浮点比较容差 1e-6）：
- 距离场是二维 ``list[list[float]]``，字形内部为正、外部为负；
- 每个值的绝对值不超过 ``spread``，``spread`` 必须是有限正数；
- 抗锯齿宽度 ``width`` 以距离场单位计，过渡带覆盖 ``[-width/2, width/2]``；
- 所有函数只依赖标准库 ``math``，不联网、不使用随机数。
"""

import math


class SdfError(ValueError):
    """距离场输入、尺寸或配置非法。"""


# ---------------------------------------------------------------------------
# 输入校验辅助
# ---------------------------------------------------------------------------


def _require_positive(value, name):
    """有限正数校验，否则抛 SdfError。"""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SdfError(f"{name} 必须是有限正数，实际为 {value!r}")
    if not math.isfinite(value) or value <= 0:
        raise SdfError(f"{name} 必须是有限正数，实际为 {value!r}")
    return float(value)


def _require_non_negative(value, name):
    """有限非负数校验，否则抛 SdfError。"""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise SdfError(f"{name} 必须是有限非负数，实际为 {value!r}")
    if not math.isfinite(value) or value < 0:
        raise SdfError(f"{name} 必须是有限非负数，实际为 {value!r}")
    return float(value)


def _check_rect(field):
    """校验 field 是非空矩形二维浮点列表，返回 (行数, 列数)。"""
    if not isinstance(field, list) or not field:
        raise SdfError("距离场必须是非空二维列表")
    width = None
    for row in field:
        if not isinstance(row, list) or not row:
            raise SdfError("距离场的每一行必须是非空列表")
        if width is None:
            width = len(row)
        elif len(row) != width:
            raise SdfError("距离场各行长度不一致（非矩形）")
        for value in row:
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise SdfError("距离场的值必须是数值")
            if not math.isfinite(value):
                raise SdfError("距离场的值必须是有限数")
    return len(field), width


def _check_bits(bits):
    """校验位图是非空矩形 0/1 二维列表，返回 (行数, 列数)。"""
    if not isinstance(bits, list) or not bits:
        raise SdfError("位图必须是非空二维列表")
    width = None
    for row in bits:
        if not isinstance(row, list) or not row:
            raise SdfError("位图的每一行必须是非空列表")
        if width is None:
            width = len(row)
        elif len(row) != width:
            raise SdfError("位图各行长度不一致（非矩形）")
        for value in row:
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise SdfError("位图只能包含 0/1")
            if value not in (0, 1):
                raise SdfError("位图只能包含 0/1")
    return len(bits), width


# ---------------------------------------------------------------------------
# 距离变换
# ---------------------------------------------------------------------------


def edt_1d(values):
    """一维平方距离变换（Felzenszwalb-Huttenlocher，单趟线性扫描）。

    对每个下标 ``i``，计算 ``min_q((i - q) ** 2 + values[q])``。
    ``values`` 中允许出现 ``math.inf``（表示该位置不是种子），不得修改
    输入。二维距离场的行、列各调用一趟即可得到平方欧氏距离。

    返回同长度的 ``list[float]``；结果按绝对误差 1e-6 比较。
    """
    if not isinstance(values, list):
        raise SdfError("values 必须是一维列表")
    n = len(values)
    if n == 0:
        return []

    # 下包络用的抛物线索引 v 与交点横坐标 z；只保留有限种子。
    envelope = []
    intersections = []

    for q in range(n):
        fq = values[q]
        if not isinstance(fq, (int, float)):
            raise SdfError("values 的元素必须是数值")
        if not math.isfinite(fq):
            if fq != math.inf:
                raise SdfError("values 只允许 math.inf 作为非种子标记")
            continue

        s = -math.inf
        while envelope:
            k = envelope[-1]
            # 新抛物线 fq 与包络末抛物线 fk 的交点横坐标。
            s = ((fq + q * q) - (values[k] + k * k)) / (2.0 * (q - k))
            if s <= intersections[-1]:
                envelope.pop()
                intersections.pop()
            else:
                break
        if not envelope:
            intersections.append(-math.inf)
        else:
            intersections.append(s)
        envelope.append(q)

    result = [0.0] * n
    if not envelope:
        # 没有任何种子：所有位置距离均为无穷。
        return [math.inf] * n

    k = 0
    for i in range(n):
        while k + 1 < len(envelope) and intersections[k + 1] < i:
            k += 1
        q = envelope[k]
        delta = i - q
        result[i] = float(delta * delta + values[q])
    return result


def _squared_edt(bits, inside):
    """二维平方欧氏距离变换，返回到 ``inside`` 格子的平方距离场。

    行、列各跑一趟线性一维变换（Felzenszwalb-Huttenlocher），
    整体复杂度 O(h * w)。
    """
    h = len(bits)
    w = len(bits[0])

    # 第一趟：沿每行做一维变换。
    grid = [[0.0] * w for _ in range(h)]
    for r in range(h):
        row = bits[r]
        seeds = [0.0 if value == inside else math.inf for value in row]
        grid[r] = edt_1d(seeds)

    # 第二趟：沿每列做一维变换。
    result = [[0.0] * w for _ in range(h)]
    for c in range(w):
        column = [grid[r][c] for r in range(h)]
        transformed = edt_1d(column)
        for r in range(h):
            result[r][c] = transformed[r]
    return result


def distance_transform(bits, spread):
    """把 0/1 位图转换为有符号欧氏距离场。

    ``bits`` 是矩形二维 0/1 列表，1 表示字形内部。每个格子的取值：

    - 内部（1）：到最近外部（0）格子的欧氏距离，取正值；
    - 外部（0）：到最近内部（1）格子的欧氏距离的相反数，取负值；
    - 若某一侧不存在，则内部取 ``spread``、外部取 ``-spread``。

    距离用格子中心计算、相邻格子间距为 1；绝对值超过 ``spread`` 时
    截断到 ``[-spread, spread]``。返回与 ``bits`` 同尺寸的新二维列表。

    ``spread`` 必须为有限正数，否则抛 ``SdfError``；位图为空、非矩形
    或含 0/1 之外的值也抛 ``SdfError``。浮点比较容差 1e-6。
    """
    spread = _require_positive(spread, "spread")
    h, w = _check_bits(bits)

    inside_dist = _squared_edt(bits, 1)
    outside_dist = _squared_edt(bits, 0)

    field = [[0.0] * w for _ in range(h)]
    for r in range(h):
        for c in range(w):
            if bits[r][c] == 1:
                # 内部取正：到最近外部格子的欧氏距离。
                d = math.sqrt(outside_dist[r][c])
                field[r][c] = min(d, spread)
            else:
                # 外部取负：到最近内部格子的欧氏距离取相反数。
                d = math.sqrt(inside_dist[r][c])
                field[r][c] = -min(d, spread)
    return field


# ---------------------------------------------------------------------------
# 平滑 / 阈值 / 描边 / 采样
# ---------------------------------------------------------------------------


def smooth_field(field, sigma):
    """用高斯核对二维距离场做平滑，边界采用复制（clamp）方式。

    ``sigma`` 必须为有限非负数。``sigma > 0`` 时，核半径
    ``r = ceil(3 * sigma)``，偏移 ``k`` 的权重为
    ``exp(-k*k / (2*sigma*sigma))`` 并在 ``[-r, r]`` 上归一化；越界坐标
    夹到最近的边缘格子，因此常量区域和图像边界不会被抹糊。
    ``sigma == 0`` 时返回输入副本。

    ``field`` 必须是矩形二维浮点列表，否则抛 ``SdfError``。结果按
    绝对误差 1e-6 比较。
    """
    sigma = _require_non_negative(sigma, "sigma")
    h, w = _check_rect(field)

    if sigma == 0.0:
        return [[float(value) for value in row] for row in field]

    radius = max(1, math.ceil(3.0 * sigma))
    weights = [
        math.exp(-(k * k) / (2.0 * sigma * sigma))
        for k in range(-radius, radius + 1)
    ]
    total = sum(weights)
    weights = [weight / total for weight in weights]

    # 可分离的高斯滤波：先水平后垂直，每趟都对越界坐标做 clamp。
    horizontal = [[0.0] * w for _ in range(h)]
    for r in range(h):
        row = field[r]
        for c in range(w):
            acc = 0.0
            for offset, weight in zip(range(-radius, radius + 1), weights):
                source = min(max(c + offset, 0), w - 1)
                acc += weight * row[source]
            horizontal[r][c] = acc

    result = [[0.0] * w for _ in range(h)]
    for c in range(w):
        for r in range(h):
            acc = 0.0
            for offset, weight in zip(range(-radius, radius + 1), weights):
                source = min(max(r + offset, 0), h - 1)
                acc += weight * horizontal[source][c]
            result[r][c] = acc
    return result


def threshold_field(field, width):
    """按抗锯齿宽度把有符号距离场压到 0..1。

    ``width`` 是过渡带宽度（距离场单位），必须为有限正数，否则抛
    ``SdfError``。对每个距离 ``d`` 返回
    ``clamp((d + width/2) / width, 0, 1)``：``d >= width/2`` 为 1，
    ``d <= -width/2`` 为 0，中间线性过渡。宽度越大过渡带越宽，缩放后
    不会糊成固定的一像素。

    返回同尺寸的新二维列表；``field`` 非矩形时抛 ``SdfError``。
    结果按绝对误差 1e-6 比较。
    """
    width = _require_positive(width, "width")
    h, w = _check_rect(field)
    half = width / 2.0

    result = [[0.0] * w for _ in range(h)]
    for r in range(h):
        for c in range(w):
            value = (field[r][c] + half) / width
            result[r][c] = min(max(value, 0.0), 1.0)
    return result


def outline_field(field, thickness):
    """把描边向外扩展 ``thickness``，内部值保持不变。

    对内部（``d >= 0``）原样返回；对外部（``d < 0``）返回
    ``d + thickness``，使描边带落在原边界外侧而不把字形加粗。
    ``thickness == 0`` 时返回副本。``thickness`` 必须为有限非负数，
    ``field`` 必须为矩形二维列表，否则抛 ``SdfError``。

    结果按绝对误差 1e-6 比较。
    """
    thickness = _require_non_negative(thickness, "thickness")
    h, w = _check_rect(field)

    result = [[0.0] * w for _ in range(h)]
    for r in range(h):
        for c in range(w):
            value = float(field[r][c])
            if value < 0.0:
                value += thickness  # 只向外扩，内部（d >= 0）保持不变
            result[r][c] = value
    return result


def sample_bilinear(field, u, v):
    """在 ``(u, v)`` 处做双线性采样；``u`` 是列、``v`` 是行。

    坐标先 clamp 到 ``[0, w-1]`` 和 ``[0, h-1]``，因此越界采样返回
    最近的边缘值而不会索引越界；整数坐标返回对应格值。

    ``field`` 为空、非矩形，或坐标非有限数时抛 ``SdfError``。结果按
    绝对误差 1e-6 比较。
    """
    h, w = _check_rect(field)
    if isinstance(u, bool) or not isinstance(u, (int, float)) or not math.isfinite(u):
        raise SdfError(f"u 必须是有限数，实际为 {u!r}")
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v):
        raise SdfError(f"v 必须是有限数，实际为 {v!r}")

    # 先 clamp 再取整，越界坐标落在边缘格子上，不会索引越界。
    u = min(max(float(u), 0.0), float(w - 1))
    v = min(max(float(v), 0.0), float(h - 1))

    x0 = int(math.floor(u))
    y0 = int(math.floor(v))
    x1 = min(x0 + 1, w - 1)
    y1 = min(y0 + 1, h - 1)
    fx = u - x0
    fy = v - y0

    top = field[y0][x0] * (1.0 - fx) + field[y0][x1] * fx
    bottom = field[y1][x0] * (1.0 - fx) + field[y1][x1] * fx
    return top * (1.0 - fy) + bottom * fy


def validate_field(field, spread):
    """检查距离场，返回问题描述列表；没有问题返回空列表。

    检查项包括：空场、各行长度不一致（尺寸不符）、值非有限或超出
    ``[-spread, spread]``。``spread`` 本身必须为有限正数，否则抛
    ``SdfError``。比较时使用绝对容差 1e-6。
    """
    spread = _require_positive(spread, "spread")

    if not isinstance(field, list) or not field:
        return ["距离场为空"]

    issues = []
    width = len(field[0]) if isinstance(field[0], list) else 0
    if width == 0:
        issues.append("距离场的行不能为空")

    for r, row in enumerate(field):
        if not isinstance(row, list):
            issues.append(f"第 {r} 行不是列表")
            continue
        if len(row) != width:
            issues.append(f"第 {r} 行长度与首行不一致（非矩形）")
        for c, value in enumerate(row):
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                issues.append(f"位置 ({r}, {c}) 的值不是数值: {value!r}")
                continue
            if not math.isfinite(value):
                issues.append(f"位置 ({r}, {c}) 的值非有限: {value!r}")
            elif value < -spread - 1e-6 or value > spread + 1e-6:
                issues.append(
                    f"位置 ({r}, {c}) 的值 {value} 超出 [-{spread}, {spread}]"
                )
    return issues


def sdf_batch(glyphs, config):
    """批量生成字形距离场，返回按输入顺序排列的 ``{name: field}``。

    ``glyphs`` 是 ``dict[name, {"w": int, "h": int, "bits": list[list[int]]}]``；
    ``config`` 至少含 ``"spread"``（有限正数），可选 ``"sigma"``
    （非负，默认 0，先做高斯平滑）。

    每个字形必须独立计算，禁止共享或复用可变的距离场缓冲，避免字形
    之间串值。任一字形或配置非法时抛 ``SdfError``，并且不返回任何
    部分结果（整体回滚）。全部成功时返回 ``dict`` 的键顺序与
    ``glyphs`` 一致。字形尺寸必须满足 ``h`` 行、每行 ``w`` 个 0/1。
    结果按绝对误差 1e-6 比较。
    """
    if not isinstance(config, dict) or "spread" not in config:
        raise SdfError("config 必须包含 spread")
    spread = _require_positive(config["spread"], "spread")
    sigma = _require_non_negative(config.get("sigma", 0), "sigma")

    if not isinstance(glyphs, dict):
        raise SdfError("glyphs 必须是字典")

    # 先整体校验，任一非法直接抛错，不产出部分结果。
    for name, glyph in glyphs.items():
        if not isinstance(glyph, dict) or "w" not in glyph or "h" not in glyph or "bits" not in glyph:
            raise SdfError(f"字形 {name!r} 缺少 w/h/bits 字段")
        w, h = glyph["w"], glyph["h"]
        if isinstance(w, bool) or not isinstance(w, int) or w <= 0:
            raise SdfError(f"字形 {name!r} 的 w 必须是正整数")
        if isinstance(h, bool) or not isinstance(h, int) or h <= 0:
            raise SdfError(f"字形 {name!r} 的 h 必须是正整数")
        rows, cols = _check_bits(glyph["bits"])
        if rows != h or cols != w:
            raise SdfError(
                f"字形 {name!r} 的位图尺寸 {rows}x{cols} 与声明 {h}x{w} 不符"
            )

    # 校验通过后再逐个独立计算，每个字形拿到全新缓冲，互不共享。
    result = {}
    for name, glyph in glyphs.items():
        field = distance_transform(glyph["bits"], spread)
        if sigma > 0.0:
            field = smooth_field(field, sigma)
        result[name] = field
    return result
