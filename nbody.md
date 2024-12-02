# PCG projekt 2
- autor: xlogin00

## Měření výkonu (čas / 100 kroků simulace)

### Průběžné
|   N   | CPU [s]  | Step 0 [s] | Step 1 [s] |
|:-----:|----------|------------|------------|
|  4096 | 0.492139 |            |            |
|  8192 | 1.471328 |            |            |
| 12288 | 2.478942 |            |            |
| 16384 | 3.386801 |            |            |
| 20480 | 5.059240 |            |            |
| 24576 | 7.112179 |            |            |
| 28672 | 9.892856 |            |            |
| 32768 | 12.59829 |            |            |
| 36864 | 15.54297 |            |            |
| 40960 | 19.36099 |            |            |
| 45056 | 23.48723 |            |            |
| 49152 | 27.69359 |            |            |
| 53248 | 32.63063 |            |            |
| 57344 | 37.43660 |            |            |
| 61440 | 42.85863 |            |            |
| 65536 | 49.46104 |            |            |
| 69632 | 55.14939 |            |            |
| 73728 | 62.04446 |            |            |
| 77824 | 69.26138 |            |            |
| 81920 | 76.60071 |            |            |

### Závěrečné
|    N   |  CPU [s] | GPU [s] | Zrychlení | Propustnost [GiB/s] | Výkon [GFLOPS] |
|:------:|:--------:|:-------:|:---------:|:-------------------:|:--------------:|
|   1024 |   1.0928 |         |           |                     |                |
|   2048 |   0.5958 |         |           |                     |                |
|   4096 |   0.6652 |         |           |                     |                |
|   8192 |   1.6599 |         |           |                     |                |
|  16384 |   3.3655 |         |           |                     |                |
|  32768 |  12.7233 |         |           |                     |                |
|  65536 |  48.9732 |         |           |                     |                |
| 131072 | 195.9965 |         |           |                     |                |

## Otázky

### Krok 0: Základní implementace
**Vyskytla se nějaká anomále v naměřených časech? Pokud ano, vysvětlete:**


### Krok 1: Sloučení kernelů
**Došlo ke zrychlení?**


**Popište hlavní důvody:**

### Krok 2: Výpočet těžiště
**Kolik kernelů je nutné použít k výpočtu?**

**Kolik další paměti jste museli naalokovat?**

**Jaké je zrychelní vůči sekveční verzi? Zdůvodněte.** *(Provedu to smyčkou #pragma acc parallel loop seq)*


### Krok 4: Měření výkonu
**Jakých jste dosáhli výsledků?**

**Lze v datech pozorovat nějaké anomálie?**