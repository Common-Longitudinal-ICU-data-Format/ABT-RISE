from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd


def plot_upset(failure_df: pd.DataFrame, boolean_cols: list[str], output_path: Path) -> None:
    """Draw an UpSet-style plot of eligibility failure intersections and save as PNG."""
    import itertools

    # Count each unique combination of boolean failure flags
    combos = []
    for r in range(1, len(boolean_cols) + 1):
        for cols in itertools.combinations(boolean_cols, r):
            mask = failure_df[list(cols)].all(axis=1)
            for col in boolean_cols:
                if col not in cols:
                    mask = mask & ~failure_df[col]
            count = mask.sum()
            if count > 0:
                combos.append((set(cols), count))

    combos.sort(key=lambda x: x[1], reverse=True)
    if not combos:
        return

    # Save aggregate summary CSV (no PHI) next to the PNG
    summary_rows = []
    for cols_set, count in combos:
        row = {col: col in cols_set for col in boolean_cols}
        row["count"] = count
        summary_rows.append(row)
    pd.DataFrame(summary_rows).to_csv(output_path.with_suffix(".csv"), index=False)

    n_combos = len(combos)
    n_cats = len(boolean_cols)

    fig, (ax_bar, ax_dot) = plt.subplots(
        2, 1, figsize=(max(6, n_combos * 0.8), 3 + n_cats * 0.4),
        gridspec_kw={"height_ratios": [3, n_cats], "hspace": 0.05},
        sharex=True,
    )

    x = range(n_combos)
    counts = [c[1] for c in combos]
    ax_bar.bar(x, counts, color="#4a90d9", edgecolor="white", linewidth=0.5)
    for i, v in enumerate(counts):
        ax_bar.text(i, v + max(counts) * 0.02, f"{v:,}", ha="center", va="bottom", fontsize=7)
    ax_bar.set_ylabel("Count")
    ax_bar.set_xlim(-0.5, n_combos - 0.5)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)

    # Dot matrix
    for i, (active_cols, _) in enumerate(combos):
        for j, col in enumerate(boolean_cols):
            color = "#333333" if col in active_cols else "#cccccc"
            size = 60 if col in active_cols else 20
            ax_dot.scatter(i, j, color=color, s=size, zorder=3)
        active_indices = [boolean_cols.index(c) for c in active_cols if c in boolean_cols]
        if len(active_indices) > 1:
            ax_dot.plot(
                [i] * 2,
                [min(active_indices), max(active_indices)],
                color="#333333", linewidth=1.5, zorder=2,
            )

    ax_dot.set_yticks(range(n_cats))
    ax_dot.set_yticklabels(boolean_cols, fontsize=8)
    ax_dot.set_xlim(-0.5, n_combos - 0.5)
    ax_dot.set_ylim(-0.5, n_cats - 0.5)
    ax_dot.invert_yaxis()
    ax_dot.set_xticks([])
    ax_dot.spines["top"].set_visible(False)
    ax_dot.spines["right"].set_visible(False)
    ax_dot.spines["bottom"].set_visible(False)
    ax_dot.grid(axis="y", linestyle="--", alpha=0.3)

    fig.savefig(str(output_path), dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_consort(consort: list[dict], output_path: Path) -> None:
    """Draw a vertical CONSORT flowchart and save as PNG."""
    n_steps = len(consort)
    fig_height = max(8, n_steps * 1.6)
    fig, ax = plt.subplots(figsize=(14, fig_height))
    fig.patch.set_facecolor("#000000")
    ax.set_facecolor("#000000")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # Layout constants
    main_x = 0.38
    excl_x = 0.78
    box_w = 0.32
    box_h_main = 0.055
    box_h_excl = 0.055
    y_top = 0.95
    y_bottom = 0.03
    y_spacing = (y_top - y_bottom) / (n_steps - 1) if n_steps > 1 else 0

    # Dark theme colors
    main_color = "#2d6a4f"
    main_edge = "#40916c"
    excl_color = "#9b2226"
    excl_edge = "#c1121f"
    start_color = "#1d3557"
    start_edge = "#457b9d"
    arrow_color = "#aaaaaa"
    text_color = "#ffffff"

    def draw_box(x_center, y_center, w, h, text, facecolor, edgecolor, fontsize=8):
        x0 = x_center - w / 2
        y0 = y_center - h / 2
        box = mpatches.FancyBboxPatch(
            (x0, y0),
            w,
            h,
            boxstyle="round,pad=0.008",
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=1.5,
        )
        ax.add_patch(box)
        ax.text(
            x_center,
            y_center,
            text,
            ha="center",
            va="center",
            fontsize=fontsize,
            fontweight="bold",
            color=text_color,
            wrap=True,
        )

    def draw_arrow(x0, y0, x1, y1):
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops=dict(
                arrowstyle="-|>",
                color=arrow_color,
                lw=1.5,
                mutation_scale=12,
            ),
        )

    # Step 0 — starting box
    step0 = consort[0]
    y0 = y_top
    draw_box(
        main_x,
        y0,
        box_w,
        box_h_main,
        f"{step0['description']}\n(n = {step0['remaining']:,})",
        start_color,
        start_edge,
        fontsize=9,
    )

    prev_y = y0

    # Steps 1..N
    for i, step in enumerate(consort[1:], start=1):
        y = y_top - i * y_spacing

        # Arrow from previous main box to this one
        draw_arrow(main_x, prev_y - box_h_main / 2, main_x, y + box_h_main / 2)

        # Main (remaining) box
        draw_box(
            main_x,
            y,
            box_w,
            box_h_main,
            f"{step['description']}\n(n = {step['remaining']:,})",
            main_color,
            main_edge,
        )

        # Exclusion side branch
        excl_y = (prev_y + y) / 2
        # Horizontal arrow from main flow to exclusion box
        draw_arrow(
            main_x + box_w / 2,
            excl_y,
            excl_x - box_w / 2,
            excl_y,
        )
        reason = step["reason"] or ""
        draw_box(
            excl_x,
            excl_y,
            box_w,
            box_h_excl,
            f"Excluded: {reason}\n(n = {step['excluded']:,})",
            excl_color,
            excl_edge,
        )

        prev_y = y

    fig.savefig(str(output_path), dpi=150, bbox_inches="tight", facecolor="#000000")
    plt.close(fig)
