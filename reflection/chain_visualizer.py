import matplotlib.pyplot as plt
from celluloid import Camera


class ChainVisualizer:

    def __init__(self) -> None:
        self.figure = plt.figure(figsize=(24, 6))
        self.camera = Camera(self.figure)

    def plot(self, chain):
        plt.plot(list(range(len(chain.atoms))), [atom.u for atom in chain.atoms.values()], color='b', marker='o')
        self.camera.snap()

    def animate(self):
        animation = self.camera.animate()
        animation.save('anim.gif')
        
    