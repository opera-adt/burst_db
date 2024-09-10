import click

from .build_frame_db import create
from .query_frame_db import intersect, lookup
from .query_historical_bursts import fetch_bursts, fetch_granules


@click.group()
def cli_app():
    """Create/interact with OPERA's burst/frame databases."""


cli_app.add_command(create)
cli_app.add_command(intersect)
cli_app.add_command(lookup)


@click.group()
def historical():
    """Sub-commands for interacting with the historical burst database."""


historical.add_command(fetch_bursts)
historical.add_command(fetch_granules)
cli_app.add_command(historical)

if __name__ == "__main__":
    cli_app()
