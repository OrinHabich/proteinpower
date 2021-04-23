import sys
sys.path.append("..")
from src.protein_plotter import ProteinPlotter


def test_plot(protein):
    ProteinPlotter.plot(protein, show_plot=False, scale_axes=False)
