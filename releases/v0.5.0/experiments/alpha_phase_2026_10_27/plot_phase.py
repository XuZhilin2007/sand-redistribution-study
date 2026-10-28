"""Plot proved sets/subsets only. Requires matplotlib, outside proof runtime."""
from fractions import Fraction as Q
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
OUT = HERE.parents[1]/'docs/alpha_phase_2026-10-27'


def main():
    data = json.loads((HERE/'exact_results.json').read_text(encoding='utf-8'))
    roots = data['central']['roots']
    a, b, a6 = [roots[k]['approximation'] for k in ('a10', 'b6', 'a6')]
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 1, figsize=(11.8, 7.2), gridspec_kw={'height_ratios': [1, 1.35]})
    ax = axes[0]
    ax.axvspan(a, b, color='#2166ac', alpha=.12)
    ax.plot([a, b], [0, 0], lw=15, color='#2166ac', solid_capstyle='butt')
    ax.plot([a, b], [0, 0], 'o', color='#2166ac', ms=8)
    oldlo, oldhi = float(Q(1999, 8000)), float(Q(2001, 8000))
    ax.plot([oldlo, oldhi], [.55, .55], lw=12, color='#d49525', solid_capstyle='butt')
    ax.text(.251, .55, 'Published v0.4.0 interval', va='center', color='#8d5d0f')
    ax.annotate('K=10: U3 = U2\nquintic boundary\n0.246154829858...', xy=(a, 0), xytext=(a-.003, -.65),
                ha='center', va='top', arrowprops={'arrowstyle': '-', 'color': '#555'})
    ax.annotate('K=6: U1 = U2\ncubic boundary\n0.291645436170...', xy=(b, 0), xytext=(b+.0007, -.65),
                ha='center', va='top', arrowprops={'arrowstyle': '-', 'color': '#555'})
    ax.set(xlim=(.239, .299), ylim=(-1.6, .95), yticks=[], xlabel='Removal fraction alpha')
    ax.set_title('Universal effective set: exactly one closed interval, determined by K=6 and K=10', loc='left', fontweight='bold')
    ax.grid(axis='x', color='#ddd', lw=.6)
    ax = axes[1]
    def bar(left, right, y, color, open_left=False):
        ax.plot([left, right], [y, y], lw=10, color=color, solid_capstyle='butt')
        ax.plot(left, y, 'o', color=color, mfc='white' if open_left else color, ms=6, mew=1.5)
        ax.plot(right, y, 'o', color=color, ms=6)
    bar(a6, b, 2, '#25846c')
    bar(a, float(Q(7, 24)), 1, '#2166ac')
    bar(float(Q(21, 25)), float(Q(43, 50)), 1, '#2166ac')
    bar(float(Q(123, 500)), float(Q(7, 24)), 0, '#2166ac')
    bar(float(Q(13, 40)), float(Q(1, 3)), 0, '#2166ac', True)
    ax.plot([.5, .9], [1, 1], 'x', color='#b13c39', ms=8, mew=2)
    ax.plot([.324, .325], [0, 0], 'x', color='#b13c39', ms=7, mew=1.5)
    ax.text(.50, 1.2, 'time 1 beats time 2', ha='center', fontsize=9, color='#953331')
    ax.text(.90, .67, 'time 4 beats time 2', ha='center', fontsize=9, color='#953331')
    ax.annotate('K=13: selector jump at 13/40\ntie point fails; right side succeeds',
                xy=(.325, 0), xytext=(.44, -.5), fontsize=9,
                arrowprops={'arrowstyle': '-', 'color': '#555'}, va='center')
    ax.set(xlim=(0, 1), ylim=(-.78, 2.5), yticks=[0, 1, 2],
           yticklabels=['K=13: certified subsets', 'K=10: certified subsets', 'K=6: complete effective set'],
           xlabel='Removal fraction alpha')
    ax.set_title('Fixed sizes can have disconnected effective sets; K=13 is not closed', loc='left', fontweight='bold')
    ax.grid(axis='x', color='#ddd', lw=.6)
    ax.legend(handles=[Line2D([0], [0], color='#2166ac', lw=7, label='Proved inclusion'),
                       Line2D([0], [0], color='#b13c39', marker='x', ls='', mew=2, label='Exact counterexample')],
              loc='upper right', frameon=False, ncol=2, fontsize=9)
    fig.text(.12, .015, 'Bars for K=10 and K=13 are certified portions, not complete classifications. White space is unclassified in those rows.', fontsize=9, color='#555')
    fig.tight_layout(rect=(0, .035, 1, 1), h_pad=2.5)
    OUT.mkdir(parents=True, exist_ok=True)
    for extension in ('png', 'svg'):
        fig.savefig(OUT/f'phase_diagram.{extension}', dpi=180, bbox_inches='tight')
    print(OUT/'phase_diagram.png')


if __name__ == '__main__':
    main()
