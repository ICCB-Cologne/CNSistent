# Get the tight bounding box of the specified ax, including tick labels, in display units
#
# Example print:
# Figure size in inches: [10.         10.17996442] (Width x Height)
# Ax size in inches: [8.41291667 0.63665135] (Width x Height)
# Ax coordinates within the figure:
#    Left: 0.07935416666666673
#    Bottom: 0.46873017812050666
#    Right: 0.9206458333333333
#    Top: 0.5312698218794933

def _get_pos(fig, ax, should_print=False):
    # Get the tight bounding box of the axes, including tick labels, in display units
    tight_bbox = ax.get_tightbbox(fig.canvas.get_renderer())
    # Transform the bounding box to figure fraction coordinates
    transformed_bbox = tight_bbox.transformed(fig.transFigure.inverted())
    # Calculate the size of the bounding box in inches
    fig_size = fig.get_size_inches()
    if should_print:
        print(f"Figure size in inches: {fig_size} (Width x Height)")
    axes_size_inches = fig_size * transformed_bbox.size
    if should_print:
        print(f"Ax size in inches: {axes_size_inches} (Width x Height)")
    # Calculate the position of the bounding box in figure fraction coordinates
    left_pos = (1 - axes_size_inches[0] / fig_size[0]) / 2
    bottom_pos = (1 - axes_size_inches[1] / fig_size[1]) / 2
    right_position = 1 - left_pos
    top_position = 1 - bottom_pos
    bound_box = [left_pos, bottom_pos, right_position, top_position]
    if should_print:
        print(f"Ax coordinates within the figure:")
        for name, val in zip(["Left", "Bottom", "Right", "Top"], bound_box):
            print(f"\t{name}: {val}")
    return fig_size, axes_size_inches, bound_box