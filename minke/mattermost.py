"""
This is an experminetal module to allow Minke production pipelines to
communicate wasily with the Mattermost chat server which is used by
LIGO.
"""

import requests

def send_mattermost(message, hook):
    # A line to explain what the Minke-bot is when it posts.
    explainer = """I'm Minke-bot. I'm here to keep an eye on the progress of MDC job production, and report it to this Mattermost channel."""

    payload={
      "channel": "minke",
      "username": "Minke-bot",
      "icon_url": "https://camo.githubusercontent.com/d2a13f9c4e6fe16dc8021887990f4bcc4f8df076/68747470733a2f2f636f64652e64616e69656c2d77696c6c69616d732e636f2e756b2f6d696e6b652f5f696d616765732f6d696e6b652e706e67",
      "text": "{} \n {}".format(message, explainer)
      }

    code = requests.post(hook, json=payload)    
