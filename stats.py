import pstats
from pstats import SortKey
p = pstats.Stats('stats_old')
p.strip_dirs().sort_stats('tottime').print_stats()