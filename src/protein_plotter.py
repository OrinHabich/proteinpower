import matplotlib.pyplot as plt
import sys
sys.path.append("..")
from src.algorithms import Algorithms


class ProteinPlotter():
    """Functionality to make a 3D plot of a folded protein."""

    @classmethod
    def plot(cls, protein, show_connections=True, scale_axes=True,
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

        if show_connections:
            axes = cls._add_connections_to_plot(axes, acids)

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
    def _add_connections_to_plot(axes, acids):
        """Add certain connections in a folding to the plot."""
        weak_bonds, strong_bonds = Algorithms.find_connections(acids)
        for connection in weak_bonds:
            x_connections = []
            y_connections = []
            z_connections = []
            for acid in connection:
                x_connections.append(acid['x'])
                y_connections.append(acid['y'])
                z_connections.append(acid['z'])
                axes.plot3D(x_connections,
                            y_connections,
                            z_connections,
                            color='black',
                            linestyle='dotted')

        for connection in strong_bonds:
            x_connections = []
            y_connections = []
            z_connections = []
            for acid in connection:
                x_connections.append(acid['x'])
                y_connections.append(acid['y'])
                z_connections.append(acid['z'])
                axes.plot3D(x_connections,
                            y_connections,
                            z_connections,
                            color='black',
                            linestyle='dashed')
        return axes
