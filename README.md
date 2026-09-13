# CHull Tucker3 en R

Implementación documentada de una simulación Monte Carlo para seleccionar la dimensionalidad Tucker3 mediante **CHull-fp** y **CHull-S**, basada en Ceulemans y Kiers (2006).

Mantenedor: **José Luis Vásquez Pérez** ([joluvasquez](https://github.com/joluvasquez)).

## Ejecutar

Requiere R base; verificado con R 4.5.1. No necesita paquetes adicionales. Desde la carpeta raíz del proyecto:

```sh
Rscript R/CHull_Ceulemans-Kiers.R
```

Se generan **225 conjuntos de datos** y se comparan **341 candidatos** por conjunto. Los CSV se escriben en `resultados/`; se reemplazan los archivos del mismo nombre. Los gráficos utilizan el dispositivo activo de R. El script limpia su entorno al comenzar: ejecútelo en una sesión nueva.

## Pruebas rápidas

```sh
Rscript tests/verificar.R
```

Comprueban reconstrucción tensorial independiente, ajuste sin ruido, identidad de energías, número de candidatos, estabilidad de un caso de empate numérico y geometría de la frontera frente a `chull()` de R. No ejecutan la simulación completa ni escriben CSV.

## Método y resultados

El ajuste se aproxima mediante SVD; no se optimiza cada modelo mediante ALS. El filtro estructural del 98 % utiliza esos ajustes aproximados. La semilla es `20260912`.

| Método | Aciertos en la revisión de v4 | Errores | Recuperación |
|---|---:|---:|---:|
| CHull-fp | 209/225 | 16 | 92,89 % |
| CHull-S | 219/225 | 6 | 97,33 % |

Estos resultados corresponden a la ejecución revisada en R 4.5.1. Los recuentos publicados de 12 y 5 errores son referencias descriptivas, no resultados obligatorios para nuevas muestras. La semilla y el entorno numérico influyen en la reproducibilidad.

Consulte [metodología, unidades, salidas y limitaciones](docs/metodologia.md).

## Estructura

- `R/CHull_Ceulemans-Kiers.R`: simulación documentada.
- `tests/verificar.R`: pruebas independientes en R base.
- `docs/metodologia.md`: documentación ampliada.
- `CITATION.cff`: metadatos para citar esta implementación.

## Referencia metodológica

Ceulemans, E., & Kiers, H. A. L. (2006). Selecting among three-mode principal component models of different types and complexities: A numerical convex hull based method. *British Journal of Mathematical and Statistical Psychology, 59*, 133–150. [DOI: 10.1348/000711005X64817](https://doi.org/10.1348/000711005X64817).

Este repositorio contiene una implementación en R; no es el software oficial de los autores del artículo.

## Licencia

Este proyecto se distribuye bajo la [licencia MIT](LICENSE).
Copyright (c) 2026 José Luis Vásquez Pérez.
