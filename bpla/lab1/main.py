from bpla import start_bpla
from save import save_array
from plot import draw_plots

if __name__ == "__main__":
    states = start_bpla()
    # save_array(states)
    draw_plots(states)
