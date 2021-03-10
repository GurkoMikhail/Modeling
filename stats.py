import pstats
from pstats import SortKey
p = pstats.Stats('efg3_fix 0.0 deg projection test stats.txt')
p.strip_dirs().sort_stats('tottime').print_stats()