from src.align import AlignSave
from src.lrt import LRT
import os

families = ['drav', 'ie', 'drav_ie', 'nostratic','monkhmer', 'munda', 'monkhmer_munda', 'mayan', 'mixezoque', 'mayan_mixezoque', 'afrasian_loloburmese', 'utoaztecan', 'mayan_utoaztecan', 'monkhmer1_utoaztecan',\
            'monkhmer1_mayan']

for family in families:
    if os.path.isfile(f"results/lrt_{family}.json"):
        continue
    AlignSave(family)
    LRT(family)