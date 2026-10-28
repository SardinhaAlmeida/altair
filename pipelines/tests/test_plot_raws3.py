"""Regression checks for the combined plotting entry point."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.text import Text

import Plot_Raws as plots


class PlotTests(unittest.TestCase):
    def test_filtered_exports_keep_word_type_and_rows(self):
        source = pd.DataFrame([
            [0, 1, "ACGT", 50, 4],
            [0, 2, "AAAA", 0, 4],
        ], columns=["genome_idx", "position", "raw", "gc_pct", "kmer_len"])
        with tempfile.TemporaryDirectory() as directory:
            for label, prefix in (("RAWs", "raws"), ("mRAWs", "mraws")):
                source.attrs["word_label"] = label
                plots.results_csv(source, directory, "example")
                for interval in ("45_55", "48_52"):
                    exported = pd.read_csv(Path(directory) / f"{prefix}_gc_{interval}_example.csv")
                    self.assertEqual(exported.raw.tolist(), ["ACGT"])

    def test_original_dot_colors(self):
        _, colors, _ = plots.kmer_styles(pd.DataFrame({"kmer_len": [11, 12, 13, 14]}))
        self.assertEqual(colors, {11: "red", 12: "blue", 13: "darkseagreen", 14: "magenta"})

    def test_short_profile_uses_centered_clipped_window(self):
        counts, gc = plots.window_profiles(
            np.array([1., 2.]), np.array([50., 80.]), 51
        )
        np.testing.assert_array_equal(counts, [3, 3])
        np.testing.assert_allclose(gc, [130 / 3, 130 / 3])

    def test_selected_outputs_and_labels_for_both_analyses(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "input.csv"
            source.write_text("0,1,ACGT,50,4\n1,70,ACGTA,40,5\n")
            for label in ("RAWs", "mRAWs"):
                saved = []

                def check_figure(fig, filename, **kwargs):
                    if "3d_raws_distribution_topdown" in str(filename):
                        ax = fig.axes[0]
                        self.assertEqual(ax.name, "3d")
                        self.assertFalse(ax.zaxis.get_visible())
                        self.assertGreater(ax.get_ylim()[0], ax.get_ylim()[1])
                    elif "topdown" in str(filename):
                        self.assertTrue(all(ax.name == "rectilinear" for ax in fig.axes))
                        self.assertEqual(fig.axes[0].get_xlabel(), "Position in Genome")
                        self.assertEqual(fig.axes[0].get_ylabel(), "Genome Index")
                        for collection in fig.axes[0].collections:
                            offsets = collection.get_offsets()
                            self.assertTrue(all(tuple(point) in {(1, 0), (70, 1)}
                                                for point in offsets))
                    texts = [t.get_text() for t in fig.findobj(Text)]
                    self.assertTrue(any("Pathogen: Example virus | Host: human" in t
                                        for t in texts))
                    if label == "mRAWs":
                        self.assertFalse(any("RAW" in t.replace("mRAW", "")
                                             for t in texts))
                    saved.append(Path(filename).name)

                with patch.object(Figure, "savefig", check_figure):
                    plots.run_plots(source, directory, "Example_virus", "human", label)
                self.assertEqual(len(saved), 11)
                self.assertNotIn("raw_genome_tracks.png", saved)
                self.assertIn("raw_position_topdown_k4.png", saved)
                self.assertIn("raw_position_topdown_k5.png", saved)
                self.assertIn("raw_gc_profile_window51_exact_target_count_gradient.png", saved)
                self.assertIn("raw_window51_target_count_window1.png", saved)
                self.assertFalse(any("raster" in name for name in saved))

    def test_shortest_k_sizes_stay_consistent_in_filtered_views(self):
        df = pd.DataFrame([
            [0, 1, "A" * 10, 40, 10],
            [0, 2, "A" * 15, 50, 15],
            [1, 1, "A" * 16, 50, 16],
        ], columns=["genome_idx", "position", "raw", "gc_pct", "kmer_len"])
        _, colors, _ = plots.kmer_styles(df)
        for label in ("RAWs", "mRAWs"):
            df.attrs["word_label"] = label
            for rows, expected in ((df, [1, 3]), (df.iloc[1:], [1]), (df.iloc[:1], [3])):
                fig = plots.plt.figure()
                ax = fig.add_subplot(projection="3d", computed_zorder=False)
                plots.scatter_occurrences(ax, rows, colors)
                self.assertEqual([c.get_sizes()[0] for c in ax.collections], expected)
                self.assertEqual(sum(len(c._offsets3d[0]) for c in ax.collections), len(rows))
                if len(ax.collections) == 2:
                    self.assertGreater(ax.collections[1].zorder, ax.collections[0].zorder)
                    self.assertEqual(list(ax.collections[1]._offsets3d[2]), [10])
                plots.plt.close(fig)

    def test_empty_input_skips_figures(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "empty.csv"
            source.touch()
            with patch.object(Figure, "savefig") as save:
                plots.run_plots(source, directory, "virus", "host")
            save.assert_not_called()


if __name__ == "__main__":
    unittest.main()
