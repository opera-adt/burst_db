import click

from .build_frame_db import create
from .query_frame_db import intersect, lookup
from .query_historical_bursts import fetch_bursts, fetch_granules


@click.group()
def cli_app():
    """Create/interact with OPERA's burst/frame databases."""
    pass


cli_app.add_command(create)
cli_app.add_command(intersect)
cli_app.add_command(lookup)
cli_app.add_command(fetch_bursts)
cli_app.add_command(fetch_granules)

if __name__ == "__main__":
    cli_app()
