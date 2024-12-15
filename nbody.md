# PCG projekt 2

-   autor: xholan11

## Měření výkonu (čas / 100 kroků simulace)

### Průběžné

|   N   | CPU [s]  | Step 0 [s] | Step 1 [s] |
| :---: | -------- | ---------- | ---------- |
| 4096  | 0.492139 | 0.126816   | 0.103891   |
| 8192  | 1.471328 | 0.249650   | 0.207355   |
| 12288 | 2.478942 | 0.372334   | 0.310779   |
| 16384 | 3.386801 | 0.521003   | 0.442817   |
| 20480 | 5.059240 | 0.650448   | 0.553413   |
| 24576 | 7.112179 | 0.779496   | 0.662052   |
| 28672 | 9.892856 | 0.943379   | 0.791775   |
| 32768 | 12.59829 | 1.123603   | 0.938089   |
| 36864 | 15.54297 | 1.273069   | 1.054807   |
| 40960 | 19.36099 | 1.414321   | 1.171734   |
| 45056 | 23.48723 | 1.945680   | 1.683561   |
| 49152 | 27.69359 | 2.129030   | 1.840568   |
| 53248 | 32.63063 | 2.307882   | 1.991440   |
| 57344 | 37.43660 | 3.014291   | 2.521638   |
| 61440 | 42.85863 | 3.238308   | 2.720996   |
| 65536 | 49.46104 | 3.450597   | 2.896619   |
| 69632 | 55.14939 | 4.235064   | 3.541781   |
| 73728 | 62.04446 | 4.536944   | 3.794009   |
| 77824 | 69.26138 | 4.788916   | 4.000583   |
| 81920 | 76.60071 | 5.038197   | 4.209727   |

### Závěrečné

|   N    | CPU [s]  |  GPU [s]  | Zrychlení | Propustnost [GiB/s] | Výkon [GFLOPS] |
| :----: | :------: | :-------: | :-------: | :-----------------: | :------------: |
|  1024  |  1.0928  | 0.031077  | 35.164269 |        0.095        |     77.53      |
|  2048  |  0.5958  | 0.057352  | 10.388478 |        0.095        |     157.00     |
|  4096  |  0.6652  | 0.111493  | 5.966294  |         0.1         |     311.6      |
|  8192  |  1.6599  | 0.218902  | 7.582845  |        0.093        |     619.98     |
| 16384  |  3.3655  | 0.471271  | 7.141326  |        0.087        |    1157.79     |
| 32768  | 12.7233  | 0.970165  | 13.114573 |        0.082        |    2180.67     |
| 65536  | 48.9732  | 2.961630  | 16.535894 |        0.052        |    2805.86     |
| 131072 | 195.9965 | 10.488941 | 18.686014 |        0.029        |    3134.71     |

## Otázky

### Krok 0: Základní implementace

**Vyskytla se nějaká anomále v naměřených časech? Pokud ano, vysvětlete:**
Větší nárůsty lze pozorovat mezi dvojicemi pro N=:

-   40960 a 45056
-   53248 a 57344
-   65536 a 69632

Příčinou mohlo být vypadnutí z cache paměti. Rovněž se mohly špatně namapovat bloky na CUDA procesory.

### Krok 1: Sloučení kernelů

**Došlo ke zrychlení?**
Ano.

**Popište hlavní důvody:**
Ze tří for cyklů zůstal jeden. Také došlo k optimalizaci některých výpočetních částí.

### Krok 2: Výpočet těžiště

**Kolik kernelů je nutné použít k výpočtu?**
Jeden kernel je spotřebován pro prvotní načtení částic do comBuffer a redukce do velikost N/2.
Následně je pro každý krok v cyklu, který běží od i=N/2 po i=0 s kroky jdoucími po polovinách z daného **i**, spuštěn jeden kernel.

**Kolik další paměti jste museli naalokovat?**
Na začátku je naalokována paměť pro comBuffer o velikosti (N/2)\*sizeof(float4).

**Jaké je zrychelní vůči sekveční verzi? Zdůvodněte.** _(Provedu to smyčkou #pragma acc parallel loop seq)_
Pro N=28672 bylo dosaženo zrychlení 49.20x (u všech ostatních N při provádění runProgressBenchmark.sh bylo zrychlení pozorováno také). Zatímco pro sekvenční
verzi byly všechny iterace prováděny postupně, u té paralelní, kde nedocházelo k závislostem pro danou iteraci, mohla být každá iterace provedena více vlákny odpovídající daným pozicím v comBuffer.

### Krok 4: Měření výkonu

**Jakých jste dosáhli výsledků?**
Pro menší množství dat je zrychlení oproti CPU největší, pak dochází k nějakému srovnávání výkonnosti a následně lze opět vidět, jak
výkon GPU předstihuje CPU.

**Lze v datech pozorovat nějaké anomálie?**
Do N=32768 lze pozorovat narůstající časy GPU po cca dvojnásobcích. Potom následují dva skoky v časech pro N=65536 a N=131072. To může být namapováním více vytvořených bloků na některé CUDA procesory.
