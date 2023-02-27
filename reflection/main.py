import matplotlib.pyplot as plt
import numpy as np


from chain import Chain
from chain_visualizer import ChainVisualizer


chain = Chain()
chain.create_atoms(n=1000)
chain.add_pair_springs()

visualizer = ChainVisualizer()

for i in range(10000):
    chain.evaluate_step(0.01)
    if i%100==0:
        visualizer.plot(chain)
visualizer.animate()