# CHull Tucker3 — versión v4 documentada

## Uso

El archivo `R/CHull_Ceulemans-Kiers.R` conserva los cálculos de v4, añade comentarios metodológicos y utiliza una ruta de salida relativa. Requiere únicamente R base.

1. Abra el script y revise `ruta_salida`. Esta edición utiliza `resultados`, relativo al directorio de trabajo.
2. Ejecútelo en una sesión nueva. `rm(list=ls())` elimina los objetos del entorno donde se ejecuta; los CSV existentes con los mismos nombres se reemplazan.
3. Desde una terminal con Rscript disponible, ejecute:

```powershell
Rscript R/CHull_Ceulemans-Kiers.R
```

Los dos gráficos usan el dispositivo gráfico activo. En RStudio aparecen en el panel de gráficos; Rscript normalmente abre un PDF predeterminado. El script no define nombres propios para guardarlos.

## Diseño y método

Se cruzan tres tamaños (200×10×10, 50×20×20 y 27×27×27), tres dimensionalidades verdaderas ((3,2,2), (3,3,3) y (4,3,2)), cinco niveles nominales de ruido (0, 15, 30, 45 y 60 %) y cinco réplicas: 225 conjuntos.

Los 341 candidatos tienen P, Q y R entre 1 y 8 y cumplen P≤QR, Q≤PR y R≤PQ. El ajuste es la energía del subnúcleo proyectado dividida entre la energía observada. Las bases se obtienen por SVD; no se optimiza cada modelo por ALS. CHull compara ese ajuste con `fp` o con `S=P+Q+R`, obtiene la frontera superior y maximiza el cociente entre pendientes consecutivas.

El filtro del 98 % se aplica a la señal sin ruido y considera candidatos con menor `fp`. Comprueba ajustes aproximados del conjunto evaluado. No constituye una certificación mediante optimización individual de cada candidato. Si la estructura no pasa, se genera otra; el límite es 500 intentos por conjunto.

## Ruido y unidades

Se genera X=M+eE con ||E||=||M|| y e=sqrt(p/(1-p)).

| Campo | Significado | Unidad |
|---|---|---|
| error_objetivo | Nivel nominal del diseño | % |
| coef_error | Multiplicador e del ruido | Coeficiente |
| error_nominal | 100e²/(1+e²) | % |
| error_energia_separada | 100 SS_ruido/(SS_senal+SS_ruido) | % |
| error_observado | 100 SS_ruido/SS_X | % |
| termino_cruzado | 2e〈M,E〉 | Suma de productos |
| fit_fp, fit_S | Ajuste aproximado seleccionado | Proporción |
| acierto_fp, acierto_S | Coincidencia exacta de P,Q,R con el rango verdadero | 0 o 1 |
| st_fp, st_S | Cociente de pendientes; puede ser Inf o NA | Cociente |
| intentos_98 | Número de generaciones hasta aceptar la estructura | Conteo |
| error_identidad_SS | Error absoluto en SS_X=SS_senal+SS_ruido+termino_cruzado | Suma de cuadrados |

El ruido observado no tiene por qué coincidir exactamente con el nominal: SS_X incluye el término cruzado. Esta diferencia no implica que la generación sea incorrecta.

## Archivos de salida

Todos los nombres llevan el sufijo `_v4.csv`.

| Prefijo | Contenido |
|---|---|
| 00_resultados_parciales | Resultados acumulados, escritos después de cada conjunto |
| 01_resultados_MC_CHull_Tucker3 | Resultados individuales; al finalizar incluye error_identidad_SS |
| 02_validacion_Ceulemans_Kiers | Comparación descriptiva con los recuentos publicados |
| 03_resultados_por_error | Porcentajes de recuperación por ruido nominal |
| 04_resultados_por_tamano | Porcentajes de recuperación por tamaño |
| 05_resultados_por_dimensionalidad | Porcentajes de recuperación por rango verdadero |
| 06_resultados_tamano_error | Recuperación por tamaño y ruido |
| 07_resultados_dimensionalidad_error | Recuperación por rango verdadero y ruido |
| 08_errores_CHull_fp | Casos de selección incorrecta por fp |
| 09_errores_CHull_S | Casos de selección incorrecta por S |
| 10_diagnostico_ruido | Media, desviación estándar, mínimo y máximo del ruido observado, y media nominal |

El CSV parcial no incluye todas las columnas agregadas al final. Permite inspeccionar el avance, pero no implementa reanudación automática. No se guardan los tensores ni todas las tablas de candidatos/fronteras.

## Reproducibilidad y límites

La semilla es 20260912. Los datos, los intentos rechazados, la prueba inicial de 341 modelos y los desempates comparten el generador aleatorio. Cambiar su orden o número puede cambiar las muestras posteriores. Las bibliotecas numéricas también pueden producir pequeñas diferencias en casos límite.

Las tolerancias son tol_fit=1e-12, tol_hull=1e-12 y tol_slope=1e-14. El umbral de pendiente es absoluto y depende de la escala de complejidad. La implementación no valida exhaustivamente entradas vacías, no finitas o dimensiones incompatibles en todas las funciones. Los parámetros del diseño actual satisfacen las condiciones utilizadas en las pruebas.

Cuando no hay st disponible, CHull devuelve el mayor ajuste como respaldo; no debe interpretarse como un codo identificado. Si cambia max_comp, tamaños o dimensionalidades, revise también las comprobaciones que exigen 341 candidatos.

Los 12 errores de CHull-fp y 5 de CHull-S son resultados publicados de una muestra concreta. Una nueva simulación no tiene que reproducir exactamente esos recuentos. La tabla denominada validacion es una comparación descriptiva, no una prueba estadística de equivalencia.

## Verificación computacional de v4

En la revisión anterior se ejecutaron las 225 simulaciones con R 4.5.1, con los siguientes resultados:

| Comprobación | Resultado |
|---|---|
| Candidatos por conjunto | 341 |
| CHull-fp | 209 aciertos, 16 errores; 92,89 % |
| CHull-S | 219 aciertos, 6 errores; 97,33 % |
| Recuperación con ruido cero | 100 % en ambos métodos |
| Caso de perturbación 10⁻¹⁵ | Selección estable |
| Comparación con chull() de R | 100 casos aleatorios satisfactorios |
| Reconstrucción tensorial independiente | Diferencia máxima 1,87×10⁻¹⁶ |
| Identidad de sumas de cuadrados | Error absoluto máximo 1,14×10⁻¹³ |

Las pruebas reproducibles están en `tests/verificar.R`.
Estas son pruebas de funcionamiento observado, no una demostración para todas las entradas. Para entregar esta copia se comprueba adicionalmente que sus expresiones ejecutables coinciden con v4 salvo la ruta de salida. Los comentarios no añaden validaciones ni cambian el algoritmo.

## Referencia

Ceulemans, E., & Kiers, H. A. L. (2006). Selecting among three-mode principal component models of different types and complexities: A numerical convex hull based method. *British Journal of Mathematical and Statistical Psychology, 59*, 133–150. [Artículo original](https://ppw.kuleuven.be/okp/_pdf/Ceulemans2006SATMP.pdf). DOI: 10.1348/000711005X64817. La comparación solo Tucker3 corresponde a la sección 6.2.

