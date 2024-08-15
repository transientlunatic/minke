import click

from .injection import injection


@click.group()
def minke():
    """
    This is the main command line program for the heron package.
    """
    pass


minke.add_command(injection)
