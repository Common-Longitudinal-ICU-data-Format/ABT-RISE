from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


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
