"""Entry point.

    python -m nanofusion [FILE.json]     start the web interface (optionally opening a saved simulation)
    python -m nanofusion run [options]   run a simulation from the command line
    python -m nanofusion improved [...]  run the improved crystal model
    python -m nanofusion implosion [...] run a hollow-shell Coulomb implosion
    python -m nanofusion ball [...]      burn a hot fuel ball held by a tamper and a magnetic field
    python -m nanofusion sim [...]       every energy route in one simulator: list, run, search
"""

from __future__ import annotations

import sys


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "run":
        from .cli import main as cli_main

        return cli_main(argv[1:])
    if argv and argv[0] == "improved":
        from .improved.cli import main as improved_main

        return improved_main(argv[1:])
    if argv and argv[0] == "implosion":
        from .improved.implosion import main as implosion_main

        return implosion_main(argv[1:])
    if argv and argv[0] == "ball":
        from .improved.confinement import main as ball_main

        return ball_main(argv[1:])
    if argv and argv[0] == "sim":
        from .improved.simulator import main as sim_main

        return sim_main(argv[1:])
    from .web.server import main as web_main

    return web_main(argv)


if __name__ == "__main__":
    sys.exit(main())
