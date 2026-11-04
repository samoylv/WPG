"""Common test setup: matplotlib draws off screen, so tests that call plt.show() open no windows."""
import matplotlib

matplotlib.use('Agg')
