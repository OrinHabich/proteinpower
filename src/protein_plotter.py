import matplotlib.pyplot as plt
import sys
sys.path.append("..")
from src.algorithms import Algorithms


class ProteinPlotter():
    """Functionality to make a 3D plot of a folded protein."""

    @classmethod
    def plot(cls, protein, show_bonds=True, scale_axes=True,
             show_plot=True):
        """Plot a protein."""
        acids = protein.acids
        x, y, z, acid_colors = cls._prepare_data_for_plot(acids)

        # Create plot
        plt.figure()
        axes = plt.axes(projection='3d')
        axes.set_xlabel('x')
        axes.set_ylabel('y')
        axes.set_zlabel('z')
        axes.plot3D(x, y, z, 'green')
        if not scale_axes:
            max_len = len(acids)
            plt.xlim([-max_len, max_len])
            plt.ylim([-max_len, max_len])
            axes.set_zlim([-max_len, max_len])

        axes.scatter3D(x, y, z, facecolor=acid_colors)

        if show_bonds:
            axes = cls._add_bonds_to_plot(axes, acids)

        if show_plot:
            plt.show()

    @staticmethod
    def _prepare_data_for_plot(acids):
        """Returns relevant data about the acids in suitable format for
        plotting.
        """
        x = []
        y = []
        z = []
        acid_colors = []
        for acid in acids:
            x.append(acid['x'])
            y.append(acid['y'])
            z.append(acid['z'])
            if acid['acid_type'] == 'H':
                acid_colors.append('r')
            elif acid['acid_type'] == 'C':
                acid_colors.append('black')
            else:
                acid_colors.append('b')
        return (x, y, z, acid_colors)

    @staticmethod
    def _add_bonds_to_plot(axes, acids):
        """Add certain bonds in a folding to the plot."""
        bonds = Algorithms.find_bonds(acids)
        linestyle = 'dotted'
        for bond_type in bonds:
            for bond in bond_type:
                x_bonds = []
                y_bonds = []
                z_bonds = []
                for acid in bond:
                    x_bonds.append(acid['x'])
                    y_bonds.append(acid['y'])
                    z_bonds.append(acid['z'])
                    axes.plot3D(x_bonds,
                                y_bonds,
                                z_bonds,
                                color='black',
                                linestyle=linestyle)
            linestyle = 'dashed'
        return axes
