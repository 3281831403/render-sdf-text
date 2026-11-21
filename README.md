# render-sdf-text

## 题面

UI 换成 SDF 字体后，斜笔画发虚、描边还把字加粗了一圈：距离场用曼哈顿距离算，斜边变成菱形；内部默认取 0，正负号反了；距离不做截断，远点距离无界；一维变换每趟重算，复杂度退化到 O(n²)；平滑用盒滤波，边缘出现台阶；抗锯齿宽度固定 1 像素，缩放后糊；描边在内部也加，字形变粗；采样越界不 clamp，索引错误；批量共用距离场缓冲，字形之间串；spread=0 不校验，除零。sdftext/field.py 是核心，tests/test_field.py 15个红。要修的是：把曼哈顿距离换成欧氏距离、内部取正值、按 spread 截断、一维变换降到线性、盒滤波换成高斯核、抗锯齿宽度随缩放走、描边只外扩、采样越界 clamp、批量独立缓冲、spread=0 显式抛 SdfError。

## 运行方式

在仓库根目录执行：

```bash
python -m pytest tests -q
```

初始状态必须恰好输出 `15 failed`，且没有 collection/import/syntax error；15 个测试都会因为调用 `sdftext/field.py` 里的桩函数而失败。修复后同一命令应变为 `15 passed`。运行离线，依赖只有标准库 `math` 与 `pytest`，不使用随机数。

## 字段表

| 字段 | 取值 |
| --- | --- |
| 所属领域 | 游戏、图形与多媒体 |
| 子方向 | 图形渲染 |
| 任务类型 | Bugfix |
| 语言/框架 | Python 3, pytest |
| 环境可复现等级 | 无外部依赖 |

## 约束

`tests/test_field.py` 不可修改。只实现 `sdftext/field.py` 中 8 个 `raise NotImplementedError` 的函数，保持 docstring 契约（距离场范围、抗锯齿宽度、容差 1e-6）与 `SdfError` 行为不变。
