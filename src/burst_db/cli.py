import logging

import click

from .build_frame_db import create
from .create_cslc_burst_catalog import make_burst_catalog
from .query_consistent_bursts import urls_for_frame
from .query_frame_db import intersect, lookup
from .query_historical_bursts import fetch_bursts, fetch_granules

logger = logging.getLogger("burst_db")


@click.group()
def cli_app():
    """Create/interact with OPERA's burst/frame databases."""
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "[%(levelname)s|%(module)s|L%(lineno)d] %(asctime)s: %(message)s"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)


cli_app.add_command(create)
cli_app.add_command(intersect)
cli_app.add_command(lookup)
cli_app.add_command(make_burst_catalog)
cli_app.add_command(urls_for_frame)


@click.group()
def historical():
    """Sub-commands for interacting with the historical burst database."""


historical.add_command(fetch_bursts)
historical.add_command(fetch_granules)
cli_app.add_command(historical)

if __name__ == "__main__":
    cli_app()
