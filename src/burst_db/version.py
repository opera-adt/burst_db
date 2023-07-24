'''
release history of `burst_db`
'''

import collections


# release history
Tag = collections.namedtuple('Tag', ['version', 'date'])
release_history = (
    Tag('0.1.2', '2023-07-24'),
    Tag('0.1.1', '2022-12-14'),
    Tag('0.1.0', '2022-12-08')
)

# latest release version number and date
release_version = release_history[0].version
release_date = release_history[0].date
