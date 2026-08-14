# KnockKnock
![](knockknock_logo.png)
KnockKnock is a tool to aid in the design of shRNA sequences for use in RNAi mediated knock-down of target genes. KnockKnock evaluates every possible guide sequence in one or more mRNAs using position-specific nucleotide efficacy data from [Matveeva et al. \(2012\)](https://doi.org/10.3389/fgene.2012.00163), then helps you turn the best candidates into orderable cloning oligos for use with commong shRNA expression backbones such as [pINDUCER](10.1073/pnas.1019736108) or  [pLKO](https://doi.org/10.1016/j.cell.2006.01.040).

## Installation
**Requires Python ≥ 3.10.**  (see [python.org](https://www.python.org/downloads/))
*Dependencies: biopython, pandas, xlsxwriter, typer, rich, textual.*
### Recommended (automatic installation)
either with `uv`  (see [Astral](https://docs.astral.sh/uv/getting-started/installation/))
```bash
uv tool install knockknock-shrna       
```
or with `pipx`:
```bash
pipx install knockknock-shrna
```
Either one installs KnockKnock into its own environment and puts the `knockknock` command (short alias `kk`) on your PATH. Nothing else in your Python installation is touched, and you never have to activate anything.
### Manage environment yourself
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install knockknock-shrna
```
or:
```bash
conda create -n knockknock python=3.11
conda activate knockknock
pip install knockknock-shrna
```
### Bare installation
Bare `pip install knockknock-shrna` with no environment works, but only if you know what you are doing: it can conflict with the dependencies of whatever else lives in that Python installation.
### Uninstall
Uninstall, matching how you installed: `uv tool uninstall knockknock-shrna`, `pipx uninstall knockknock-shrna`, or `pip uninstall knockknock-shrna`.

The distribution is named `knockknock-shrna` (an unrelated ML training-notification package already holds plain `knockknock` on PyPI); the import package and the command are both plain `knockknock`. In an isolated install the two never meet — only a bare system-wide `pip install` of both can collide on the `knockknock` command.

Step-by-step install instructions, including uninstall, are in
[docs/USAGE.md](docs/USAGE.md#4-requirements-and-installation).

## Quick start
Start the terminal UI with `knockknock` or its alias `kk`.
This will present a file-picker which allows you to open FASTA files containing your target mRNA sequence. You can also copy & paste sequences in basic format or FASTA format.

For other usage examples and detailed documentation consult:
| Doc                                      | Covers                                                       |
|------------------------------------------|--------------------------------------------------------------|
| [docs/CLI_MANUAL.md](docs/CLI_MANUAL.md) | Every scoring, oligo-design, and `--info` option, with examples. |
| [docs/TUI_MANUAL.md](docs/TUI_MANUAL.md) | Every interactive screen, key binding, and procedure (file picker, explorer, oligo designer). |

## Scientific background
- **Position-specific efficacy coefficients**: [Matveeva et al., PLOS ONE
  7(2): e30811 (2012)](https://doi.org/10.1371/journal.pone.0030811).
- **Nearest-neighbor ΔG**: [Xia et al., Biochemistry 37(42): 14719–14735
  (1998)](https://doi.org/10.1021/bi9809425).

## License
[MIT](LICENSE)
