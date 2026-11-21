"""sdftext.field 的验收测试：初始状态 15 个测试全部失败。"""

import math

import pytest

from sdftext.field import (
    SdfError,
    distance_transform,
    edt_1d,
    outline_field,
    sample_bilinear,
    sdf_batch,
    smooth_field,
    threshold_field,
    validate_field,
)


def _flat(field):
    """把二维场拉平，便于做容差比较。"""
    return [value for row in field for value in row]


# ---------------------------------------------------------------------------
# 10 个缺陷对应的回归测试
# ---------------------------------------------------------------------------


def test_distance_transform_uses_euclidean_distance():
    """缺陷 1：不能用曼哈顿距离，斜边应为 sqrt(2) 而不是 2。"""
    bits = [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ]
    field = distance_transform(bits, 10.0)
    assert field[2][2] == pytest.approx(1.0, abs=1e-6)
    assert field[1][1] == pytest.approx(-math.sqrt(2.0), abs=1e-6)
    assert field[3][3] == pytest.approx(-math.sqrt(2.0), abs=1e-6)
    assert field[1][2] == pytest.approx(-1.0, abs=1e-6)


def test_distance_transform_inside_positive_outside_negative():
    """缺陷 2：字形内部为正，外部为负。"""
    bits = [
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 1, 1, 1, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
    ]
    field = distance_transform(bits, 4.0)
    assert field[2][2] > 0.0
    assert field[2][1] > 0.0
    assert field[0][2] < 0.0
    assert field[4][2] < 0.0
    assert validate_field(field, 4.0) == []


def test_distance_transform_clips_far_points_to_spread():
    """缺陷 3：超过 spread 的距离必须截断，不能无界增长。"""
    bits = [[1, 0, 0, 0, 0, 0, 0]]
    field = distance_transform(bits, 2.0)
    assert field[0][6] == pytest.approx(-2.0, abs=1e-6)
    assert field[0][3] == pytest.approx(-2.0, abs=1e-6)
    assert validate_field(field, 2.0) == []


def test_edt_1d_computes_squared_distances_linearly():
    """缺陷 4：一维变换应得到正确的平方距离，而不是每趟重算的错值。"""
    values = [math.inf, 0.0, math.inf, math.inf, 0.0, math.inf]
    expected = [1.0, 0.0, 1.0, 1.0, 0.0, 1.0]
    assert edt_1d(values) == pytest.approx(expected, abs=1e-6)

    n = 64
    seeds = (0, 31, 63)
    long_values = [0.0 if index in seeds else math.inf for index in range(n)]
    brute = [float(min((index - q) ** 2 for q in seeds)) for index in range(n)]
    assert edt_1d(long_values) == pytest.approx(brute, abs=1e-6)


def test_smooth_field_matches_gaussian_kernel():
    """缺陷 5：必须是高斯核，盒滤波的台阶值会在此失败。"""
    field = [[0.0, 0.0, 0.0, 1.0, 1.0, 1.0]]
    sigma = 1.0
    radius = 3
    weights = [
        math.exp(-(offset * offset) / (2.0 * sigma * sigma))
        for offset in range(-radius, radius + 1)
    ]
    total = sum(weights)
    row = field[0]
    expected = []
    for column in range(len(row)):
        acc = 0.0
        for offset, weight in zip(range(-radius, radius + 1), weights):
            source = min(max(column + offset, 0), len(row) - 1)
            acc += weight * row[source]
        expected.append(acc / total)

    smoothed = smooth_field(field, sigma)
    assert len(smoothed) == 1
    assert smoothed[0] == pytest.approx(expected, abs=1e-6)


def test_threshold_field_respects_antialias_width():
    """缺陷 6：过渡带必须随 width 缩放，不能固定为 1 像素。"""
    field = [[-2.0, -0.5, 0.0, 0.5, 2.0]]
    narrow = threshold_field(field, 1.0)
    wide = threshold_field(field, 4.0)
    assert narrow[0] == pytest.approx([0.0, 0.0, 0.5, 1.0, 1.0], abs=1e-6)
    assert wide[0] == pytest.approx([0.0, 0.375, 0.5, 0.625, 1.0], abs=1e-6)


def test_outline_field_keeps_interior_unchanged():
    """缺陷 7：描边只向外扩，内部值不得改变。"""
    field = [[-3.0, -1.0, 0.0, 0.5, 2.0]]
    result = outline_field(field, 1.5)
    assert result[0] == pytest.approx([-1.5, 0.5, 0.0, 0.5, 2.0], abs=1e-6)


def test_sample_bilinear_clamps_out_of_bounds():
    """缺陷 8：越界坐标必须 clamp，而不是索引越界或取错格子。"""
    field = [[0.0, 10.0], [20.0, 30.0]]
    assert sample_bilinear(field, -5.0, -5.0) == pytest.approx(0.0, abs=1e-6)
    assert sample_bilinear(field, 5.0, 5.0) == pytest.approx(30.0, abs=1e-6)
    assert sample_bilinear(field, 10.0, 0.0) == pytest.approx(10.0, abs=1e-6)
    assert sample_bilinear(field, 0.5, 0.5) == pytest.approx(15.0, abs=1e-6)


def test_sdf_batch_isolates_glyphs():
    """缺陷 9：批量各字形必须使用独立缓冲；非法字形要整体回滚。"""
    glyphs = {
        "vertical": {
            "w": 5,
            "h": 5,
            "bits": [
                [0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0],
                [0, 0, 1, 0, 0],
            ],
        },
        "horizontal": {
            "w": 5,
            "h": 5,
            "bits": [
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
                [1, 1, 1, 1, 1],
                [0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0],
            ],
        },
    }
    config = {"spread": 3.0}
    result = sdf_batch(glyphs, config)
    assert _flat(result["vertical"]) == pytest.approx(
        _flat(distance_transform(glyphs["vertical"]["bits"], 3.0)), abs=1e-6
    )
    assert _flat(result["horizontal"]) == pytest.approx(
        _flat(distance_transform(glyphs["horizontal"]["bits"], 3.0)), abs=1e-6
    )
    assert _flat(result["vertical"]) != pytest.approx(
        _flat(result["horizontal"]), abs=1e-6
    )

    with pytest.raises(SdfError):
        sdf_batch(
            {
                "ok": {"w": 2, "h": 2, "bits": [[1, 1], [1, 1]]},
                "bad": {"w": 3, "h": 2, "bits": [[1, 1, 1]]},
            },
            config,
        )


def test_distance_transform_rejects_zero_spread():
    """缺陷 10：spread=0 必须抛 SdfError，不能除零。"""
    bits = [[1, 0], [0, 0]]
    with pytest.raises(SdfError):
        distance_transform(bits, 0.0)


# ---------------------------------------------------------------------------
# 5 个补充/组合测试
# ---------------------------------------------------------------------------


def test_all_ones_bitmap_has_positive_field():
    """全 1 位图没有外部，距离场应全为正且被 spread 截断。"""
    bits = [[1] * 4 for _ in range(4)]
    field = distance_transform(bits, 3.0)
    assert all(value > 0.0 for value in _flat(field))
    assert _flat(field) == pytest.approx([3.0] * 16, abs=1e-6)
    assert validate_field(field, 3.0) == []


def test_all_zeros_bitmap_has_negative_field():
    """全 0 位图没有内部，距离场应全为负且被 spread 截断。"""
    bits = [[0] * 4 for _ in range(4)]
    field = distance_transform(bits, 3.0)
    assert all(value < 0.0 for value in _flat(field))
    assert _flat(field) == pytest.approx([-3.0] * 16, abs=1e-6)


def test_sdf_batch_rejects_zero_spread():
    """批量入口同样要拒绝 spread=0。"""
    glyphs = {"a": {"w": 2, "h": 2, "bits": [[1, 0], [0, 0]]}}
    with pytest.raises(SdfError):
        sdf_batch(glyphs, {"spread": 0})


def test_sample_bilinear_at_grid_point_matches_cell_value():
    """在整数格点上采样应精确等于格值（容差 1e-6）。"""
    field = [[0.0, 10.0, 20.0], [30.0, 40.0, 50.0]]
    for row in range(2):
        for column in range(3):
            sampled = sample_bilinear(field, float(column), float(row))
            assert sampled == pytest.approx(field[row][column], abs=1e-6)


def test_sdf_batch_preserves_order_on_success():
    """批量全部成功时，返回字典的键顺序与输入一致。"""
    glyphs = {
        "first": {"w": 3, "h": 3, "bits": [[1, 1, 1], [1, 0, 1], [1, 1, 1]]},
        "second": {"w": 3, "h": 3, "bits": [[0, 0, 0], [0, 1, 0], [0, 0, 0]]},
        "third": {"w": 3, "h": 3, "bits": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]},
    }
    result = sdf_batch(glyphs, {"spread": 2.0})
    assert list(result.keys()) == list(glyphs.keys())
    for name, glyph in glyphs.items():
        assert len(result[name]) == glyph["h"]
        assert all(len(row) == glyph["w"] for row in result[name])
