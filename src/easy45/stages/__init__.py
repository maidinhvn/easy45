"""Pipeline stages. Each module exposes ``run(config, state) -> dict``.

A stage reads its inputs from ``state`` (named output paths produced by
earlier stages) and/or from ``config``, writes intermediate files under
``config.workdir``, and returns a dict of named output paths that later
stages consume.
"""
