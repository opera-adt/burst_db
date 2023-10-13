import typer

from .build_frame_db import create
from .query_frame_db import intersect, lookup
from .query_historical_bursts import fetch_granules

app = typer.Typer()
app.command()(create)
app.command()(intersect)
app.command()(lookup)
app.command()(fetch_granules)