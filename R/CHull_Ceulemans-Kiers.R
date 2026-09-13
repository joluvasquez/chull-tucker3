# =====================================================================
# DOCUMENTACION DE LA IMPLEMENTACION v4
# =====================================================================
# Proposito: estudiar la recuperacion de la dimensionalidad Tucker3 con
# CHull-fp y CHull-S. Esta edicion adapta la salida a resultados/ y conserva las
# expresiones de calculo de la version v4.
#
# Referencia: Ceulemans, E. y Kiers, H. A. L. (2006). Selecting among
# three-mode principal component models of different types and
# complexities: A numerical convex hull based method.
# British Journal of Mathematical and Statistical Psychology, 59, 133-150.
# DOI: 10.1348/000711005X64817. Comparacion solo Tucker3: seccion 6.2.
#
# Diseno: 3 tamanos x 3 rangos verdaderos x 5 niveles x 5 replicas = 225.
# Candidatos: P,Q,R entre 1 y 8, sujetos a P<=QR, Q<=PR y R<=PQ (341).
# Ajuste: aproximacion rapida por SVD de los desplegamientos, proyectando
# X sobre las primeras columnas de cada base y sumando el nucleo retenido.
# No se realiza una optimizacion ALS independiente por candidato.
#
# Reproducibilidad: semilla 20260912; requiere R base. Los sorteos de
# empates y los intentos rechazados consumen numeros aleatorios. Cambiar
# esas operaciones puede cambiar los conjuntos generados posteriormente.
# Diferencias entre plataformas/SVD tambien pueden afectar casos limite.
#
# Ejecucion: revisar ruta_salida antes de ejecutar. El script borra objetos
# del entorno con rm(list=ls()) y escribe CSV de nombres fijos, reemplazando
# los existentes. Se recomienda ejecutarlo en una sesion nueva de Rscript.
# Los graficos usan el dispositivo activo; no se les asigna archivo aqui.
#
# Interpretacion: 12 y 5 errores son referencias del articulo para CHull-fp
# y CHull-S, no criterios automaticos de aprobacion de otra muestra.
# Esta documentacion no convierte las pruebas realizadas en una garantia
# de correccion para todas las entradas. Ver README.md para sus limites.
# =====================================================================

# =====================================================================
# REPLICACION MONTE CARLO - CEULEMANS & KIERS (2006)
# SOLO MODELOS TUCKER3
#
# VERSION v4
#
# CORRECCIONES INCORPORADAS:
#
# 1. Calibracion correcta del nivel NOMINAL de error:
#
#          e^2
#    p = ---------
#        1 + e^2
#
#    e = sqrt(p / (1-p))
#
# 2. Se diferencia entre:
#       - error objetivo
#       - error nominal
#       - error observado
#
# 3. Se guarda el termino cruzado:
#
#    ||X||^2 =
#    ||M||^2 + ||eE||^2 + 2e<M,E>
#
# 4. CHull utiliza tolerancias numericas explicitas
#    y consistentes.
#
# 5. Correccion de fp cuando una via supera el producto
#    de las otras dos.
#
# 6. 341 modelos Tucker3 candidatos.
#
# 7. Filtro estructural del 98%.
#
# 8. Replica original:
#
#    3 tamanos
#    x 3 dimensionalidades
#    x 5 niveles de error
#    x 5 replicas
#    = 225 datasets
#
# =====================================================================



# =====================================================================
# 0. LIMPIAR ENTORNO
# =====================================================================

rm(list = ls())



# =====================================================================
# 1. CONFIGURACION GENERAL
# =====================================================================

set.seed(20260912)


ruta_salida <- "resultados"


if (!dir.exists(ruta_salida)) {
  
  dir.create(
    ruta_salida,
    recursive = TRUE
  )
}


cat("\n=============================================\n")
cat("REPLICACION CEULEMANS & KIERS - VERSION v4\n")
cat("=============================================\n\n")

cat("Resultados en:\n")
cat(ruta_salida, "\n\n")



# =====================================================================
# 2. TOLERANCIAS NUMERICAS
# =====================================================================

# Tolerancia para comparar ajustes f

tol_fit <- 1e-12


# Tolerancia geometrica para determinar si
# un punto esta sobre la recta del convex hull

tol_hull <- 1e-12


# Tolerancia para pendientes cercanas a cero

tol_slope <- 1e-14


cat("Tolerancias numericas:\n")
cat("tol_fit   =", tol_fit, "\n")
cat("tol_hull  =", tol_hull, "\n")
cat("tol_slope =", tol_slope, "\n\n")



# =====================================================================
# 3. DISEÑO EXPERIMENTAL
# =====================================================================

tamanos <- list(
  
  c(200, 10, 10),
  
  c(50, 20, 20),
  
  c(27, 27, 27)
  
)


nombres_tamanos <- c(
  
  "200x10x10",
  
  "50x20x20",
  
  "27x27x27"
  
)



# ---------------------------------------------------------------------
# Dimensionalidades verdaderas
# ---------------------------------------------------------------------

dimensiones_true <- list(
  
  c(3, 2, 2),
  
  c(3, 3, 3),
  
  c(4, 3, 2)
  
)



# =====================================================================
# 4. NIVELES DE ERROR
#
# p = e^2 / (1 + e^2)
#
# e = sqrt(p / (1-p))
# =====================================================================

proporcion_error <- c(
  
  0.00,
  
  0.15,
  
  0.30,
  
  0.45,
  
  0.60
  
)


coef_error <- sqrt(
  
  proporcion_error /
    
    (1 - proporcion_error)
  
)


coef_error[
  proporcion_error == 0
] <- 0



tabla_error <- data.frame(
  
  error_objetivo = proporcion_error * 100,
  
  coeficiente_e = coef_error,
  
  error_nominal =
    
    100 *
    
    coef_error^2 /
    
    (1 + coef_error^2)
  
)


cat("Niveles de error utilizados:\n\n")

print(
  tabla_error,
  row.names = FALSE
)

cat("\n")



# =====================================================================
# 5. PARAMETROS DEL MONTE CARLO
# =====================================================================

n_replicas <- 5


max_comp <- 8



# =====================================================================
# 6. FUNCION: ORTONORMALIZACION
# =====================================================================

# CONTRATO: M es una matriz numerica finita, con filas >= columnas.
# Devuelve una base ortonormal con el mismo numero de columnas.
# No verifica explicitamente dimensiones ni rango de entrada.
ortonormalizar <- function(M) {
  
  
  qrM <- qr(M)
  
  
  Q <- qr.Q(qrM)
  
  
  Q <- Q[
    
    ,
    
    seq_len(
      ncol(M)
    ),
    
    drop = FALSE
    
  ]
  
  
  return(Q)
  
}



# =====================================================================
# 7. FUNCION: PRODUCTO N-MODO
# =====================================================================

# CONTRATO: X es un array; modo es un indice valido de dim(X).
# M tiene dim(X)[modo] columnas. El resultado conserva los otros modos
# y reemplaza la dimension elegida por nrow(M). Valida compatibilidad
# matricial, pero no todas las posibles entradas mal formadas.
producto_modo <- function(
    X,
    M,
    modo) {
  
  
  dims <- dim(X)
  
  
  N <- length(dims)
  
  
  perm <- c(
    
    modo,
    
    setdiff(
      seq_len(N),
      modo
    )
    
  )
  
  
  Xperm <- aperm(
    X,
    perm
  )
  
  
  dims_perm <- dim(Xperm)
  
  
  Xmat <- matrix(
    
    Xperm,
    
    nrow = dims_perm[1]
    
  )
  
  
  if (
    ncol(M) != nrow(Xmat)
  ) {
    
    stop(
      
      paste0(
        
        "Dimensiones incompatibles en producto_modo(). ",
        
        "ncol(M) = ",
        ncol(M),
        
        "; nrow(Xmat) = ",
        nrow(Xmat)
        
      )
    )
  }
  
  
  Ymat <- M %*% Xmat
  
  
  dims_nuevas <- c(
    
    nrow(M),
    
    dims_perm[-1]
    
  )
  
  
  Yperm <- array(
    
    Ymat,
    
    dim = dims_nuevas
    
  )
  
  
  perm_inversa <- order(perm)
  
  
  Y <- aperm(
    
    Yperm,
    
    perm_inversa
    
  )
  
  
  return(Y)
  
}



# =====================================================================
# 8. FUNCIONES: DESPLEGAMIENTOS
# =====================================================================

desplegar_modo1 <- function(X) {
  
  
  matrix(
    
    X,
    
    nrow = dim(X)[1]
    
  )
  
}



desplegar_modo2 <- function(X) {
  
  
  matrix(
    
    aperm(
      X,
      c(2, 1, 3)
    ),
    
    nrow = dim(X)[2]
    
  )
  
}



desplegar_modo3 <- function(X) {
  
  
  matrix(
    
    aperm(
      X,
      c(3, 1, 2)
    ),
    
    nrow = dim(X)[3]
    
  )
  
}



# =====================================================================
# 9. FUNCION: NUMERO DE PARAMETROS LIBRES
#
# fp =
# IP + JQ + KR + PQR - P^2 - Q^2 - R^2
#
# Si una dimension excede el producto de las otras:
#
# I* = min(I, JK)
# J* = min(J, IK)
# K* = min(K, IJ)
#
# =====================================================================

# Complejidad usada por CHull-fp:
# I*P + J*Q + K*R + P*Q*R - P^2 - Q^2 - R^2,
# con dimensiones efectivas min(I,JK), min(J,IK), min(K,IJ).
# Esta correccion sigue la convencion del articulo para modos grandes.
calcular_fp <- function(
    I,
    J,
    K,
    P,
    Q,
    R) {
  
  
  I_fp <- min(
    I,
    J * K
  )
  
  
  J_fp <- min(
    J,
    I * K
  )
  
  
  K_fp <- min(
    K,
    I * J
  )
  
  
  fp <-
    
    I_fp * P +
    
    J_fp * Q +
    
    K_fp * R +
    
    P * Q * R -
    
    P^2 -
    
    Q^2 -
    
    R^2
  
  
  return(fp)
  
}



# =====================================================================
# 10. FUNCION: GENERAR DATOS TUCKER3
#
# X = M + eE
#
# M = G x1 A x2 B x3 C
#
# A ~ N(0,1)
#
# B y C ortonormales
#
# G ~ U(-0.5, 0.5)
#
# E ~ N(0,1)
#
# Se escala E para cumplir:
#
# ||M||^2 = ||E||^2
#
# =====================================================================

# SALIDA: lista con X, estructura M, factores, nucleo y diagnosticos.
# E se escala para que ||E||=||M||; ruido=e*E y X=M+ruido.
# Los tres porcentajes devueltos estan en escala 0-100:
# error_nominal = 100*e^2/(1+e^2);
# error_energia_separada = 100*SS_ruido/(SS_senal+SS_ruido);
# error_observado = 100*SS_ruido/SS_X.
# Los dos primeros coinciden salvo redondeo. El tercero incluye en su
# denominador el termino cruzado y no es una particion aditiva de SS_X.
# Se esperan dimensiones/rangos positivos y compatibles, y e finito >=0.
generar_tucker3 <- function(
    I,
    J,
    K,
    P,
    Q,
    R,
    e) {
  
  
  # -------------------------------------------------------------------
  # A
  # -------------------------------------------------------------------
  
  A <- matrix(
    
    rnorm(
      I * P
    ),
    
    nrow = I,
    
    ncol = P
    
  )
  
  
  
  # -------------------------------------------------------------------
  # B
  # -------------------------------------------------------------------
  
  B0 <- matrix(
    
    rnorm(
      J * Q
    ),
    
    nrow = J,
    
    ncol = Q
    
  )
  
  
  B <- ortonormalizar(B0)
  
  
  
  # -------------------------------------------------------------------
  # C
  # -------------------------------------------------------------------
  
  C0 <- matrix(
    
    rnorm(
      K * R
    ),
    
    nrow = K,
    
    ncol = R
    
  )
  
  
  C <- ortonormalizar(C0)
  
  
  
  # -------------------------------------------------------------------
  # G
  # -------------------------------------------------------------------
  
  G <- array(
    
    runif(
      
      P * Q * R,
      
      min = -0.5,
      
      max = 0.5
      
    ),
    
    dim = c(
      P,
      Q,
      R
    )
    
  )
  
  
  
  # -------------------------------------------------------------------
  # Parte estructural M
  # -------------------------------------------------------------------
  
  M <- producto_modo(
    G,
    A,
    1
  )
  
  
  M <- producto_modo(
    M,
    B,
    2
  )
  
  
  M <- producto_modo(
    M,
    C,
    3
  )
  
  
  
  # -------------------------------------------------------------------
  # Error E
  # -------------------------------------------------------------------
  
  E <- array(
    
    rnorm(
      I * J * K
    ),
    
    dim = c(
      I,
      J,
      K
    )
    
  )
  
  
  
  # -------------------------------------------------------------------
  # Igualar energia de señal y error base
  #
  # ||M|| = ||E||
  # -------------------------------------------------------------------
  
  norma_M <- sqrt(
    sum(M^2)
  )
  
  
  norma_E <- sqrt(
    sum(E^2)
  )
  
  
  if (
    !is.finite(norma_M) ||
    !is.finite(norma_E) ||
    norma_E <= 0
  ) {
    
    stop(
      "Problema al normalizar la matriz de error."
    )
  }
  
  
  E <- E *
    norma_M /
    norma_E
  
  
  
  # -------------------------------------------------------------------
  # Error perturbado
  # -------------------------------------------------------------------
  
  ruido <- e * E
  
  
  
  # -------------------------------------------------------------------
  # Datos observados
  # -------------------------------------------------------------------
  
  X <- M + ruido
  
  
  
  # ===================================================================
  # DIAGNOSTICO DEL NIVEL DE ERROR
  # ===================================================================
  
  SS_senal <- sum(
    M^2
  )
  
  
  SS_error_base <- sum(
    E^2
  )
  
  
  SS_ruido <- sum(
    ruido^2
  )
  
  
  SS_X <- sum(
    X^2
  )
  
  
  
  # -------------------------------------------------------------------
  # Termino cruzado
  #
  # 2e<M,E>
  # -------------------------------------------------------------------
  
  termino_cruzado <- 2 * e * sum(
    M * E
  )
  
  
  
  # -------------------------------------------------------------------
  # Error nominal / esperado
  # -------------------------------------------------------------------
  
  error_nominal <-
    
    100 *
    
    e^2 /
    
    (1 + e^2)
  
  
  
  # -------------------------------------------------------------------
  # Proporcion respecto a las energias separadas
  #
  # Esta debe coincidir numericamente con el nivel nominal
  # debido al escalamiento ||M||^2 = ||E||^2
  # -------------------------------------------------------------------
  
  error_energia_separada <-
    
    100 *
    
    SS_ruido /
    
    (
      SS_senal +
        SS_ruido
    )
  
  
  
  # -------------------------------------------------------------------
  # Error observado en X
  #
  # Este valor NO tiene por que coincidir exactamente con
  # el nivel nominal, debido al termino cruzado.
  # -------------------------------------------------------------------
  
  error_observado <-
    
    100 *
    
    SS_ruido /
    
    SS_X
  
  
  
  return(
    
    list(
      
      X = X,
      
      M = M,
      
      A = A,
      
      B = B,
      
      C = C,
      
      G = G,
      
      E = E,
      
      ruido = ruido,
      
      SS_senal =
        SS_senal,
      
      SS_error_base =
        SS_error_base,
      
      SS_ruido =
        SS_ruido,
      
      SS_X =
        SS_X,
      
      termino_cruzado =
        termino_cruzado,
      
      error_nominal =
        error_nominal,
      
      error_energia_separada =
        error_energia_separada,
      
      error_observado =
        error_observado
      
    )
    
  )
  
}



# =====================================================================
# 11. FUNCION: AJUSTES TUCKER3
#
# Aproximacion SVD / HOSVD
# =====================================================================

# SALIDA: data.frame con P,Q,R,S,fp,fit; fit es una proporcion (no %).
# Cada ajuste usa las mismas bases SVD de X y un subnucleo distinto.
# Para el diseno actual max_comp=8 produce 341 candidatos admisibles.
# La comprobacion de 341 en el programa principal presupone ese diseno:
# modificar max_comp o usar modos pequenos exige revisar esa comprobacion.
ajustes_tucker3 <- function(
    X,
    max_comp = 8) {
  
  
  I <- dim(X)[1]
  
  J <- dim(X)[2]
  
  K <- dim(X)[3]
  
  
  
  # -------------------------------------------------------------------
  # Desplegamientos
  # -------------------------------------------------------------------
  
  X1 <- desplegar_modo1(X)
  
  X2 <- desplegar_modo2(X)
  
  X3 <- desplegar_modo3(X)
  
  
  
  # -------------------------------------------------------------------
  # SVD
  # -------------------------------------------------------------------
  
  sv1 <- svd(
    
    X1,
    
    nu = min(
      I,
      max_comp
    ),
    
    nv = 0
    
  )
  
  
  sv2 <- svd(
    
    X2,
    
    nu = min(
      J,
      max_comp
    ),
    
    nv = 0
    
  )
  
  
  sv3 <- svd(
    
    X3,
    
    nu = min(
      K,
      max_comp
    ),
    
    nv = 0
    
  )
  
  
  
  U1 <- sv1$u
  
  U2 <- sv2$u
  
  U3 <- sv3$u
  
  
  
  Pmax <- min(
    max_comp,
    ncol(U1)
  )
  
  
  Qmax <- min(
    max_comp,
    ncol(U2)
  )
  
  
  Rmax <- min(
    max_comp,
    ncol(U3)
  )
  
  
  
  # -------------------------------------------------------------------
  # Núcleo máximo
  # -------------------------------------------------------------------
  
  H <- producto_modo(
    
    X,
    
    t(
      U1[
        ,
        1:Pmax,
        drop = FALSE
      ]
    ),
    
    1
    
  )
  
  
  H <- producto_modo(
    
    H,
    
    t(
      U2[
        ,
        1:Qmax,
        drop = FALSE
      ]
    ),
    
    2
    
  )
  
  
  H <- producto_modo(
    
    H,
    
    t(
      U3[
        ,
        1:Rmax,
        drop = FALSE
      ]
    ),
    
    3
    
  )
  
  
  
  SSX <- sum(
    X^2
  )
  
  
  if (
    !is.finite(SSX) ||
    SSX <= 0
  ) {
    
    stop(
      "La suma de cuadrados de X no es valida."
    )
  }
  
  
  
  resultados_modelos <- data.frame(
    
    P = integer(),
    
    Q = integer(),
    
    R = integer(),
    
    S = integer(),
    
    fp = numeric(),
    
    fit = numeric()
    
  )
  
  
  
  contador <- 1
  
  
  
  for (P in 1:Pmax) {
    
    
    for (Q in 1:Qmax) {
      
      
      for (R in 1:Rmax) {
        
        
        
        # ---------------------------------------------------------------
        # Excluir soluciones Tucker3 redundantes
        # ---------------------------------------------------------------
        
        if (
          
          P > Q * R ||
          
          Q > P * R ||
          
          R > P * Q
          
        ) {
          
          next
          
        }
        
        
        
        # ---------------------------------------------------------------
        # Subnúcleo
        # ---------------------------------------------------------------
        
        Gsub <- H[
          
          1:P,
          
          1:Q,
          
          1:R,
          
          drop = FALSE
          
        ]
        
        
        
        # ---------------------------------------------------------------
        # Fit aproximado
        # ---------------------------------------------------------------
        
        fit <-
          
          sum(Gsub^2) /
          
          SSX
        
        
        
        # ---------------------------------------------------------------
        # Suma de componentes
        # ---------------------------------------------------------------
        
        S <-
          
          P +
          Q +
          R
        
        
        
        # ---------------------------------------------------------------
        # Numero de parametros libres
        # ---------------------------------------------------------------
        
        fp <- calcular_fp(
          
          I = I,
          
          J = J,
          
          K = K,
          
          P = P,
          
          Q = Q,
          
          R = R
          
        )
        
        
        
        resultados_modelos[
          contador,
        ] <- data.frame(
          
          P = P,
          
          Q = Q,
          
          R = R,
          
          S = S,
          
          fp = fp,
          
          fit = fit
          
        )
        
        
        
        contador <-
          contador + 1
        
      }
      
    }
    
  }
  
  
  
  rownames(
    resultados_modelos
  ) <- NULL
  
  
  
  return(
    resultados_modelos
  )
  
}



# =====================================================================
# 12. FUNCION CHULL - VERSION NUMERICAMENTE ESTABLE
# =====================================================================

# CONTRATO: tabla no vacia con fit y la complejidad solicitada, numericos
# y finitos. P,Q,R se conservan para identificar el modelo seleccionado.
# No se validan exhaustivamente tablas vacias, NA, Inf ni tolerancias.
# Devuelve mejores_por_complejidad, hull (con st) y seleccionado.
#
# tol_fit: empates de ajuste, dominancia y empates de st finito.
# tol_hull: margen vertical para eliminar puntos sobre/bajo una cuerda.
# tol_slope: umbral absoluto de pendiente posterior tratada como cero.
# Son tolerancias con usos distintos; no aseguran invariancia numerica
# ante cualquier cambio de escala de la complejidad.
#
# Los empates de ajuste y de st maximo se sortean. Si una pendiente
# posterior se considera cero, st se fija a Inf. Los extremos no tienen
# st. Si ningun st esta disponible, se devuelve el mayor fit: es una
# regla de respaldo, no la identificacion de un codo mediante el scree test.
CHull <- function(
    tabla,
    variable_complejidad = "fp",
    tol_fit = 1e-12,
    tol_hull = 1e-12,
    tol_slope = 1e-14) {
  
  
  # -------------------------------------------------------------------
  # Validacion
  # -------------------------------------------------------------------
  
  if (
    !(variable_complejidad %in% names(tabla))
  ) {
    
    stop(
      
      paste(
        
        "No existe la variable:",
        
        variable_complejidad
        
      )
      
    )
    
  }
  
  
  
  x <- tabla[[variable_complejidad]]
  
  
  
  # ===================================================================
  # PASO 1
  # Retener la mejor solucion por cada valor de complejidad
  # ===================================================================
  
  valores <- sort(
    unique(x)
  )
  
  
  mejores <- NULL
  
  
  
  for (v in valores) {
    
    
    temp <- tabla[
      
      x == v,
      
      ,
      
      drop = FALSE
      
    ]
    
    
    
    max_fit <- max(
      temp$fit
    )
    
    
    
    # -----------------------------------------------------------------
    # Empates numericos
    #
    # Cualquier solucion a distancia <= tol_fit del maximo
    # se considera numericamente empatada.
    # -----------------------------------------------------------------
    
    candidatos <- which(
      
      abs(
        temp$fit -
          max_fit
      ) <= tol_fit
      
    )
    
    
    
    if (
      length(candidatos) > 1
    ) {
      
      
      indice_mejor <- sample(
        
        candidatos,
        
        size = 1
        
      )
      
      
    } else {
      
      
      indice_mejor <-
        candidatos[1]
      
    }
    
    
    
    mejor <- temp[
      
      indice_mejor,
      
      ,
      
      drop = FALSE
      
    ]
    
    
    
    mejores <- rbind(
      
      mejores,
      
      mejor
      
    )
    
  }
  
  
  
  # -------------------------------------------------------------------
  # Ordenar por complejidad
  # -------------------------------------------------------------------
  
  mejores <- mejores[
    
    order(
      mejores[[variable_complejidad]]
    ),
    
    ,
    
    drop = FALSE
    
  ]
  
  
  
  rownames(mejores) <- NULL
  
  
  
  # ===================================================================
  # PASO 2
  # Eliminar soluciones dominadas
  #
  # Se elimina i solamente cuando existe una solucion anterior
  # cuyo ajuste es superior por mas de tol_fit.
  # ===================================================================
  
  conservar <- rep(
    
    TRUE,
    
    nrow(mejores)
    
  )
  
  
  
  if (
    nrow(mejores) >= 2
  ) {
    
    
    mejor_fit_previo <-
      mejores$fit[1]
    
    
    
    for (
      i in 2:nrow(mejores)
    ) {
      
      
      if (
        
        mejores$fit[i] <
        mejor_fit_previo -
        tol_fit
        
      ) {
        
        
        conservar[i] <-
          FALSE
        
        
      } else {
        
        
        mejor_fit_previo <- max(
          
          mejor_fit_previo,
          
          mejores$fit[i]
          
        )
        
      }
      
    }
    
  }
  
  
  
  hull <- mejores[
    
    conservar,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  
  rownames(hull) <- NULL
  
  
  
  # ===================================================================
  # PASO 3
  # FRONTERA CONVEXA SUPERIOR
  # ===================================================================
  
  repetir <- TRUE
  
  
  
  while (
    
    repetir &&
    
    nrow(hull) >= 3
    
  ) {
    
    
    repetir <- FALSE
    
    
    
    xx <- hull[[variable_complejidad]]
    
    
    yy <- hull$fit
    
    
    
    eliminar_indice <-
      NA_integer_
    
    
    
    for (
      i in 2:(nrow(hull) - 1)
    ) {
      
      
      x1 <- xx[i - 1]
      
      x2 <- xx[i]
      
      x3 <- xx[i + 1]
      
      
      
      y1 <- yy[i - 1]
      
      y2 <- yy[i]
      
      y3 <- yy[i + 1]
      
      
      
      y_linea <-
        
        y1 +
        
        (y3 - y1) *
        
        (x2 - x1) /
        
        (x3 - x1)
      
      
      
      # ----------------------------------------------------------------
      # El punto se considera debajo o sobre la recta
      # si no supera la recta por mas de tol_hull.
      # ----------------------------------------------------------------
      
      if (
        
        y2 <=
        y_linea +
        tol_hull
        
      ) {
        
        
        eliminar_indice <- i
        
        
        break
        
      }
      
    }
    
    
    
    if (
      !is.na(eliminar_indice)
    ) {
      
      
      hull <- hull[
        
        -eliminar_indice,
        
        ,
        
        drop = FALSE
        
      ]
      
      
      repetir <- TRUE
      
    }
    
  }
  
  
  
  rownames(hull) <- NULL
  
  
  
  # ===================================================================
  # PASO 4
  # SCREE TEST
  #
  #              pendiente anterior
  # st_i = -------------------------------
  #              pendiente posterior
  # ===================================================================
  
  hull$st <- NA_real_
  
  
  
  if (
    nrow(hull) >= 3
  ) {
    
    
    xx <- hull[[variable_complejidad]]
    
    
    
    for (
      i in 2:(nrow(hull) - 1)
    ) {
      
      
      pendiente_antes <-
        
        (
          hull$fit[i] -
            hull$fit[i - 1]
        ) /
        
        (
          xx[i] -
            xx[i - 1]
        )
      
      
      
      pendiente_despues <-
        
        (
          hull$fit[i + 1] -
            hull$fit[i]
        ) /
        
        (
          xx[i + 1] -
            xx[i]
        )
      
      
      
      # ----------------------------------------------------------------
      # Pendiente posterior numericamente cero
      # ----------------------------------------------------------------
      
      if (
        
        is.finite(
          pendiente_antes
        ) &&
        
        abs(
          pendiente_despues
        ) <= tol_slope
        
      ) {
        
        
        hull$st[i] <- Inf
        
        
        
      } else if (
        
        is.finite(
          pendiente_antes
        ) &&
        
        is.finite(
          pendiente_despues
        ) &&
        
        pendiente_despues >
        tol_slope
        
      ) {
        
        
        hull$st[i] <-
          
          pendiente_antes /
          
          pendiente_despues
        
        
        
      } else {
        
        
        hull$st[i] <-
          NA_real_
        
      }
      
    }
    
  }
  
  
  
  # ===================================================================
  # PASO 5
  # SELECCION DEL MAYOR st
  # ===================================================================
  
  candidatos_st <- which(
    
    !is.na(
      hull$st
    )
    
  )
  
  
  
  if (
    length(candidatos_st) == 0
  ) {
    
    
    seleccionado <- hull[
      
      which.max(
        hull$fit
      ),
      
      ,
      
      drop = FALSE
      
    ]
    
    
  } else {
    
    
    valores_st <- hull$st[
      candidatos_st
    ]
    
    
    
    max_st <- max(
      valores_st
    )
    
    
    
    # -----------------------------------------------------------------
    # Si existe Inf, conservar solamente los Inf
    # -----------------------------------------------------------------
    
    if (
      is.infinite(max_st)
    ) {
      
      
      candidatos_max <- candidatos_st[
        
        is.infinite(
          hull$st[candidatos_st]
        )
        
      ]
      
      
    } else {
      
      
      candidatos_max <- candidatos_st[
        
        abs(
          
          hull$st[candidatos_st] -
            max_st
          
        ) <= tol_fit
        
      ]
      
    }
    
    
    
    # -----------------------------------------------------------------
    # Si hubiera empate numerico de st,
    # elegir al azar para mantener la logica de empates.
    # -----------------------------------------------------------------
    
    if (
      length(candidatos_max) > 1
    ) {
      
      
      mejor_indice <- sample(
        
        candidatos_max,
        
        size = 1
        
      )
      
      
    } else {
      
      
      mejor_indice <-
        candidatos_max[1]
      
    }
    
    
    
    seleccionado <- hull[
      
      mejor_indice,
      
      ,
      
      drop = FALSE
      
    ]
    
  }
  
  
  
  return(
    
    list(
      
      mejores_por_complejidad =
        mejores,
      
      hull =
        hull,
      
      seleccionado =
        seleccionado
      
    )
    
  )
  
}



# =====================================================================
# 13. FILTRO ESTRUCTURAL DEL 98%
# =====================================================================

# ALCANCE DEL FILTRO: acepta M si todos los candidatos evaluados con fp
# menor al verdadero tienen fit aproximado <=0.98+tol_fit.
# Se aplica a M (senal sin ruido), no a X. Si no hay candidatos menores,
# devuelve TRUE. Usa los ajustes SVD de ajustes_tucker3(); no certifica
# el optimo obtenido mediante ALS para cada candidato.
# El bucle principal vuelve a generar los datos si devuelve FALSE y
# detiene la simulacion tras 500 intentos fallidos para un conjunto.
estructura_valida_98 <- function(
    M,
    P0,
    Q0,
    R0,
    max_comp = 8,
    tol_fit = 1e-12) {
  
  
  I <- dim(M)[1]
  
  J <- dim(M)[2]
  
  K <- dim(M)[3]
  
  
  
  fp_true <- calcular_fp(
    
    I = I,
    
    J = J,
    
    K = K,
    
    P = P0,
    
    Q = Q0,
    
    R = R0
    
  )
  
  
  
  tabla_M <- ajustes_tucker3(
    
    M,
    
    max_comp =
      max_comp
    
  )
  
  
  
  modelos_menores <- tabla_M[
    
    tabla_M$fp <
      fp_true,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  
  if (
    nrow(modelos_menores) == 0
  ) {
    
    return(TRUE)
    
  }
  
  
  
  max_fit_menor <- max(
    
    modelos_menores$fit,
    
    na.rm = TRUE
    
  )
  
  
  
  if (
    !is.finite(max_fit_menor)
  ) {
    
    return(FALSE)
    
  }
  
  
  
  return(
    
    max_fit_menor <=
      0.98 +
      tol_fit
    
  )
  
}



# =====================================================================
# 14. VERIFICAR EL NUMERO DE MODELOS CANDIDATOS
# =====================================================================

X_test <- array(
  
  rnorm(
    27 * 27 * 27
  ),
  
  dim = c(
    27,
    27,
    27
  )
  
)



tabla_test <- ajustes_tucker3(
  
  X_test,
  
  max_comp = 8
  
)



cat("\n=============================================\n")
cat("VERIFICACION DEL ESPACIO DE MODELOS\n")
cat("=============================================\n\n")


cat(
  
  "Numero de modelos Tucker3 =",
  
  nrow(tabla_test),
  
  "\n"
  
)



if (
  nrow(tabla_test) != 341
) {
  
  
  stop(
    
    paste0(
      
      "ERROR: se esperaban 341 modelos y se obtuvieron ",
      
      nrow(tabla_test),
      
      "."
      
    )
    
  )
  
  
  
} else {
  
  
  cat(
    "CORRECTO: se obtienen 341 modelos Tucker3.\n\n"
  )
  
}



rm(
  X_test,
  tabla_test
)



# =====================================================================
# 15. DATA FRAME PRINCIPAL
# =====================================================================

resultados <- data.frame(
  
  
  simulacion =
    integer(),
  
  
  tamano =
    character(),
  
  
  I =
    integer(),
  
  
  J =
    integer(),
  
  
  K =
    integer(),
  
  
  I_fp =
    integer(),
  
  
  J_fp =
    integer(),
  
  
  K_fp =
    integer(),
  
  
  P_true =
    integer(),
  
  
  Q_true =
    integer(),
  
  
  R_true =
    integer(),
  
  
  S_true =
    integer(),
  
  
  error_objetivo =
    numeric(),
  
  
  coef_error =
    numeric(),
  
  
  error_nominal =
    numeric(),
  
  
  error_energia_separada =
    numeric(),
  
  
  error_observado =
    numeric(),
  
  
  SS_senal =
    numeric(),
  
  
  SS_ruido =
    numeric(),
  
  
  SS_X =
    numeric(),
  
  
  termino_cruzado =
    numeric(),
  
  
  replica =
    integer(),
  
  
  intentos_98 =
    integer(),
  
  
  
  # CHull-fp
  
  P_fp =
    integer(),
  
  
  Q_fp =
    integer(),
  
  
  R_fp =
    integer(),
  
  
  S_fp =
    integer(),
  
  
  fp_seleccionado =
    numeric(),
  
  
  fit_fp =
    numeric(),
  
  
  st_fp =
    numeric(),
  
  
  acierto_fp =
    integer(),
  
  
  
  # CHull-S
  
  P_S =
    integer(),
  
  
  Q_S =
    integer(),
  
  
  R_S =
    integer(),
  
  
  S_seleccionado =
    integer(),
  
  
  fp_S =
    numeric(),
  
  
  fit_S =
    numeric(),
  
  
  st_S =
    numeric(),
  
  
  acierto_S =
    integer(),
  
  
  stringsAsFactors =
    FALSE
  
)



# =====================================================================
# 16. MONTE CARLO
# =====================================================================

contador_global <- 1



total_simulaciones <-
  
  length(tamanos) *
  
  length(dimensiones_true) *
  
  length(proporcion_error) *
  
  n_replicas



for (
  tt in seq_along(tamanos)
) {
  
  
  I <- tamanos[[tt]][1]
  
  J <- tamanos[[tt]][2]
  
  K <- tamanos[[tt]][3]
  
  
  
  # -------------------------------------------------------------------
  # Dimensiones efectivas para fp
  # -------------------------------------------------------------------
  
  I_fp <- min(
    I,
    J * K
  )
  
  
  J_fp <- min(
    J,
    I * K
  )
  
  
  K_fp <- min(
    K,
    I * J
  )
  
  
  
  for (
    dd in seq_along(dimensiones_true)
  ) {
    
    
    P0 <- dimensiones_true[[dd]][1]
    
    Q0 <- dimensiones_true[[dd]][2]
    
    R0 <- dimensiones_true[[dd]][3]
    
    
    
    for (
      ee in seq_along(proporcion_error)
    ) {
      
      
      p_error <- proporcion_error[ee]
      
      
      e_coef <- coef_error[ee]
      
      
      
      for (
        rr in seq_len(n_replicas)
      ) {
        
        
        cat("\n---------------------------------------------\n")
        
        cat(
          "Simulacion:",
          contador_global,
          "de",
          total_simulaciones,
          "\n"
        )
        
        
        cat(
          "Tamano:",
          I,
          "x",
          J,
          "x",
          K,
          "\n"
        )
        
        
        cat(
          "Tucker3 verdadero: (",
          P0,
          ",",
          Q0,
          ",",
          R0,
          ")\n",
          sep = ""
        )
        
        
        cat(
          "Error objetivo:",
          p_error * 100,
          "%\n"
        )
        
        
        cat(
          "Coeficiente e:",
          round(
            e_coef,
            6
          ),
          "\n"
        )
        
        
        cat(
          "Replica:",
          rr,
          "\n"
        )
        
        
        
        # ==============================================================
        # GENERAR ESTRUCTURA QUE CUMPLA EL FILTRO DEL 98%
        # ==============================================================
        
        valido <- FALSE
        
        
        intento <- 0
        
        
        
        while (
          !valido
        ) {
          
          
          intento <-
            intento + 1
          
          
          
          sim <- generar_tucker3(
            
            I = I,
            
            J = J,
            
            K = K,
            
            P = P0,
            
            Q = Q0,
            
            R = R0,
            
            e = e_coef
            
          )
          
          
          
          valido <- estructura_valida_98(
            
            M = sim$M,
            
            P0 = P0,
            
            Q0 = Q0,
            
            R0 = R0,
            
            max_comp =
              max_comp,
            
            tol_fit =
              tol_fit
            
          )
          
          
          
          if (
            intento >= 500 &&
            !valido
          ) {
            
            
            stop(
              
              paste0(
                
                "No se obtuvo una estructura valida tras 500 intentos. ",
                
                "Tamano = ",
                I,
                "x",
                J,
                "x",
                K,
                
                "; dimensionalidad = (",
                P0,
                ",",
                Q0,
                ",",
                R0,
                ")."
                
              )
              
            )
            
          }
          
        }
        
        
        
        # ==============================================================
        # AJUSTAR LOS 341 MODELOS
        # ==============================================================
        
        tabla_modelos <- ajustes_tucker3(
          
          X =
            sim$X,
          
          max_comp =
            max_comp
          
        )
        
        
        
        if (
          nrow(tabla_modelos) != 341
        ) {
          
          
          stop(
            
            paste0(
              
              "Se esperaban 341 modelos Tucker3 y se obtuvieron ",
              
              nrow(tabla_modelos),
              
              "."
              
            )
            
          )
          
        }
        
        
        
        # ==============================================================
        # CHull-fp
        # ==============================================================
        
        resultado_fp <- CHull(
          
          tabla =
            tabla_modelos,
          
          variable_complejidad =
            "fp",
          
          tol_fit =
            tol_fit,
          
          tol_hull =
            tol_hull,
          
          tol_slope =
            tol_slope
          
        )
        
        
        
        sel_fp <-
          resultado_fp$seleccionado
        
        
        
        # ==============================================================
        # CHull-S
        # ==============================================================
        
        resultado_S <- CHull(
          
          tabla =
            tabla_modelos,
          
          variable_complejidad =
            "S",
          
          tol_fit =
            tol_fit,
          
          tol_hull =
            tol_hull,
          
          tol_slope =
            tol_slope
          
        )
        
        
        
        sel_S <-
          resultado_S$seleccionado
        
        
        
        # ==============================================================
        # ACIERTO EXACTO
        # ==============================================================
        
        acierto_fp <- as.integer(
          
          sel_fp$P == P0 &&
            
            sel_fp$Q == Q0 &&
            
            sel_fp$R == R0
          
        )
        
        
        
        acierto_S <- as.integer(
          
          sel_S$P == P0 &&
            
            sel_S$Q == Q0 &&
            
            sel_S$R == R0
          
        )
        
        
        
        # ==============================================================
        # GUARDAR RESULTADOS
        # ==============================================================
        
        resultados[
          contador_global,
        ] <- data.frame(
          
          
          simulacion =
            contador_global,
          
          
          tamano =
            nombres_tamanos[tt],
          
          
          I =
            I,
          
          
          J =
            J,
          
          
          K =
            K,
          
          
          I_fp =
            I_fp,
          
          
          J_fp =
            J_fp,
          
          
          K_fp =
            K_fp,
          
          
          P_true =
            P0,
          
          
          Q_true =
            Q0,
          
          
          R_true =
            R0,
          
          
          S_true =
            P0 + Q0 + R0,
          
          
          error_objetivo =
            p_error * 100,
          
          
          coef_error =
            e_coef,
          
          
          error_nominal =
            sim$error_nominal,
          
          
          error_energia_separada =
            sim$error_energia_separada,
          
          
          error_observado =
            sim$error_observado,
          
          
          SS_senal =
            sim$SS_senal,
          
          
          SS_ruido =
            sim$SS_ruido,
          
          
          SS_X =
            sim$SS_X,
          
          
          termino_cruzado =
            sim$termino_cruzado,
          
          
          replica =
            rr,
          
          
          intentos_98 =
            intento,
          
          
          
          # ------------------------------------------------------------
          # CHull-fp
          # ------------------------------------------------------------
          
          P_fp =
            sel_fp$P,
          
          
          Q_fp =
            sel_fp$Q,
          
          
          R_fp =
            sel_fp$R,
          
          
          S_fp =
            sel_fp$S,
          
          
          fp_seleccionado =
            sel_fp$fp,
          
          
          fit_fp =
            sel_fp$fit,
          
          
          st_fp =
            sel_fp$st,
          
          
          acierto_fp =
            acierto_fp,
          
          
          
          # ------------------------------------------------------------
          # CHull-S
          # ------------------------------------------------------------
          
          P_S =
            sel_S$P,
          
          
          Q_S =
            sel_S$Q,
          
          
          R_S =
            sel_S$R,
          
          
          S_seleccionado =
            sel_S$S,
          
          
          fp_S =
            sel_S$fp,
          
          
          fit_S =
            sel_S$fit,
          
          
          st_S =
            sel_S$st,
          
          
          acierto_S =
            acierto_S,
          
          
          stringsAsFactors =
            FALSE
          
        )
        
        
        
        cat(
          "Error nominal:",
          round(
            sim$error_nominal,
            4
          ),
          "%\n"
        )
        
        
        cat(
          "Error observado:",
          round(
            sim$error_observado,
            4
          ),
          "%\n"
        )
        
        
        cat(
          "Termino cruzado:",
          format(
            sim$termino_cruzado,
            scientific = TRUE
          ),
          "\n"
        )
        
        
        cat(
          "CHull-fp = (",
          sel_fp$P,
          ",",
          sel_fp$Q,
          ",",
          sel_fp$R,
          ") | acierto = ",
          acierto_fp,
          "\n",
          sep = ""
        )
        
        
        cat(
          "CHull-S  = (",
          sel_S$P,
          ",",
          sel_S$Q,
          ",",
          sel_S$R,
          ") | acierto = ",
          acierto_S,
          "\n",
          sep = ""
        )
        
        
        
        # ==============================================================
        # GUARDADO PARCIAL
        # ==============================================================
        
        write.csv(
          
          resultados,
          
          file =
            file.path(
              
              ruta_salida,
              
              "00_resultados_parciales_v4.csv"
              
            ),
          
          row.names =
            FALSE
          
        )
        
        
        
        contador_global <-
          contador_global + 1
        
      }
      
    }
    
  }
  
}



# =====================================================================
# 17. ETIQUETAS DE DIMENSIONALIDAD
# =====================================================================

resultados$dim_true <- paste0(
  
  "(",
  resultados$P_true,
  ",",
  resultados$Q_true,
  ",",
  resultados$R_true,
  ")"
  
)



resultados$dim_fp <- paste0(
  
  "(",
  resultados$P_fp,
  ",",
  resultados$Q_fp,
  ",",
  resultados$R_fp,
  ")"
  
)



resultados$dim_S <- paste0(
  
  "(",
  resultados$P_S,
  ",",
  resultados$Q_S,
  ",",
  resultados$R_S,
  ")"
  
)



# =====================================================================
# 18. RESULTADOS GENERALES
# =====================================================================

errores_fp <- sum(
  
  resultados$acierto_fp == 0
  
)


errores_S <- sum(
  
  resultados$acierto_S == 0
  
)



aciertos_fp <- sum(
  
  resultados$acierto_fp == 1
  
)


aciertos_S <- sum(
  
  resultados$acierto_S == 1
  
)



porcentaje_fp <- mean(
  
  resultados$acierto_fp
  
) * 100



porcentaje_S <- mean(
  
  resultados$acierto_S
  
) * 100



# =====================================================================
# 19. VALIDACION CON LOS RESULTADOS DEL ARTICULO
# =====================================================================

validacion <- data.frame(
  
  
  Metodo = c(
    
    "CHull-fp",
    
    "CHull-P+Q+R"
    
  ),
  
  
  Errores_articulo = c(
    
    12,
    
    5
    
  ),
  
  
  Aciertos_articulo = c(
    
    213,
    
    220
    
  ),
  
  
  Porcentaje_articulo = c(
    
    213 / 225 * 100,
    
    220 / 225 * 100
    
  ),
  
  
  Errores_replicacion = c(
    
    errores_fp,
    
    errores_S
    
  ),
  
  
  Aciertos_replicacion = c(
    
    aciertos_fp,
    
    aciertos_S
    
  ),
  
  
  Porcentaje_replicacion = c(
    
    porcentaje_fp,
    
    porcentaje_S
    
  )
  
)



validacion$Diferencia_pp <-
  
  validacion$Porcentaje_replicacion -
  
  validacion$Porcentaje_articulo



# =====================================================================
# 20. RESULTADOS POR NIVEL DE ERROR
# =====================================================================

por_error_fp <- aggregate(
  
  acierto_fp ~ error_objetivo,
  
  data =
    resultados,
  
  FUN =
    mean
  
)



por_error_S <- aggregate(
  
  acierto_S ~ error_objetivo,
  
  data =
    resultados,
  
  FUN =
    mean
  
)



por_error <- merge(
  
  por_error_fp,
  
  por_error_S,
  
  by =
    "error_objetivo"
  
)



por_error$CHull_fp <-
  
  por_error$acierto_fp *
  100



por_error$CHull_S <-
  
  por_error$acierto_S *
  100



por_error <- por_error[
  
  ,
  
  c(
    
    "error_objetivo",
    
    "CHull_fp",
    
    "CHull_S"
    
  )
  
]



# =====================================================================
# 21. RESULTADOS POR TAMAÑO
# =====================================================================

por_tamano_fp <- aggregate(
  
  acierto_fp ~ tamano,
  
  data =
    resultados,
  
  FUN =
    mean
  
)



por_tamano_S <- aggregate(
  
  acierto_S ~ tamano,
  
  data =
    resultados,
  
  FUN =
    mean
  
)



por_tamano <- merge(
  
  por_tamano_fp,
  
  por_tamano_S,
  
  by =
    "tamano"
  
)



por_tamano$CHull_fp <-
  
  por_tamano$acierto_fp *
  100



por_tamano$CHull_S <-
  
  por_tamano$acierto_S *
  100



por_tamano <- por_tamano[
  
  ,
  
  c(
    
    "tamano",
    
    "CHull_fp",
    
    "CHull_S"
    
  )
  
]



# =====================================================================
# 22. RESULTADOS POR DIMENSIONALIDAD
# =====================================================================

por_dim_fp <- aggregate(
  
  acierto_fp ~ dim_true,
  
  data =
    resultados,
  
  FUN =
    mean
  
)



por_dim_S <- aggregate(
  
  acierto_S ~ dim_true,
  
  data =
    resultados,
  
  FUN =
    mean
  
)



por_dim <- merge(
  
  por_dim_fp,
  
  por_dim_S,
  
  by =
    "dim_true"
  
)



por_dim$CHull_fp <-
  
  por_dim$acierto_fp *
  100



por_dim$CHull_S <-
  
  por_dim$acierto_S *
  100



por_dim <- por_dim[
  
  ,
  
  c(
    
    "dim_true",
    
    "CHull_fp",
    
    "CHull_S"
    
  )
  
]



# =====================================================================
# 23. TAMAÑO x ERROR
# =====================================================================

tamano_error_fp <- aggregate(
  
  acierto_fp ~
    
    tamano +
    
    error_objetivo,
  
  data =
    resultados,
  
  FUN =
    mean
  
)



tamano_error_S <- aggregate(
  
  acierto_S ~
    
    tamano +
    
    error_objetivo,
  
  data =
    resultados,
  
  FUN =
    mean
  
)



tamano_error <- merge(
  
  tamano_error_fp,
  
  tamano_error_S,
  
  by = c(
    
    "tamano",
    
    "error_objetivo"
    
  )
  
)



tamano_error$CHull_fp <-
  
  tamano_error$acierto_fp *
  100



tamano_error$CHull_S <-
  
  tamano_error$acierto_S *
  100



tamano_error <- tamano_error[
  
  ,
  
  c(
    
    "tamano",
    
    "error_objetivo",
    
    "CHull_fp",
    
    "CHull_S"
    
  )
  
]



# =====================================================================
# 24. DIMENSIONALIDAD x ERROR
# =====================================================================

dim_error_fp <- aggregate(
  
  acierto_fp ~
    
    dim_true +
    
    error_objetivo,
  
  data =
    resultados,
  
  FUN =
    mean
  
)



dim_error_S <- aggregate(
  
  acierto_S ~
    
    dim_true +
    
    error_objetivo,
  
  data =
    resultados,
  
  FUN =
    mean
  
)



dim_error <- merge(
  
  dim_error_fp,
  
  dim_error_S,
  
  by = c(
    
    "dim_true",
    
    "error_objetivo"
    
  )
  
)



dim_error$CHull_fp <-
  
  dim_error$acierto_fp *
  100



dim_error$CHull_S <-
  
  dim_error$acierto_S *
  100



dim_error <- dim_error[
  
  ,
  
  c(
    
    "dim_true",
    
    "error_objetivo",
    
    "CHull_fp",
    
    "CHull_S"
    
  )
  
]



# =====================================================================
# 25. DIAGNOSTICO DEL NIVEL DE RUIDO
# =====================================================================

ruido_nominal <- aggregate(
  
  error_nominal ~ error_objetivo,
  
  data =
    resultados,
  
  FUN =
    mean
  
)



ruido_observado_media <- aggregate(
  
  error_observado ~ error_objetivo,
  
  data =
    resultados,
  
  FUN =
    mean
  
)



ruido_observado_sd <- aggregate(
  
  error_observado ~ error_objetivo,
  
  data =
    resultados,
  
  FUN =
    sd
  
)



ruido_observado_min <- aggregate(
  
  error_observado ~ error_objetivo,
  
  data =
    resultados,
  
  FUN =
    min
  
)



ruido_observado_max <- aggregate(
  
  error_observado ~ error_objetivo,
  
  data =
    resultados,
  
  FUN =
    max
  
)



diagnostico_ruido <- merge(
  
  ruido_nominal,
  
  ruido_observado_media,
  
  by =
    "error_objetivo",
  
  suffixes = c(
    
    "_nominal",
    
    "_observado"
    
  )
  
)



names(
  diagnostico_ruido
)[
  names(diagnostico_ruido) ==
    "error_nominal"
] <- "error_nominal_medio"



names(
  diagnostico_ruido
)[
  names(diagnostico_ruido) ==
    "error_observado"
] <- "error_observado_medio"



diagnostico_ruido <- merge(
  
  diagnostico_ruido,
  
  ruido_observado_sd,
  
  by =
    "error_objetivo"
  
)



names(
  diagnostico_ruido
)[
  names(diagnostico_ruido) ==
    "error_observado"
] <- "sd_error_observado"



diagnostico_ruido <- merge(
  
  diagnostico_ruido,
  
  ruido_observado_min,
  
  by =
    "error_objetivo"
  
)



names(
  diagnostico_ruido
)[
  names(diagnostico_ruido) ==
    "error_observado"
] <- "min_error_observado"



diagnostico_ruido <- merge(
  
  diagnostico_ruido,
  
  ruido_observado_max,
  
  by =
    "error_objetivo"
  
)



names(
  diagnostico_ruido
)[
  names(diagnostico_ruido) ==
    "error_observado"
] <- "max_error_observado"



# =====================================================================
# 26. IDENTIFICAR ERRORES DE SELECCION
# =====================================================================

errores_CHull_fp <- resultados[
  
  resultados$acierto_fp == 0,
  
  ,
  
  drop = FALSE
  
]



errores_CHull_S <- resultados[
  
  resultados$acierto_S == 0,
  
  ,
  
  drop = FALSE
  
]



# =====================================================================
# 27. GUARDAR ARCHIVOS
# =====================================================================

write.csv(
  
  resultados,
  
  file.path(
    
    ruta_salida,
    
    "01_resultados_MC_CHull_Tucker3_v4.csv"
    
  ),
  
  row.names = FALSE
  
)



write.csv(
  
  validacion,
  
  file.path(
    
    ruta_salida,
    
    "02_validacion_Ceulemans_Kiers_v4.csv"
    
  ),
  
  row.names = FALSE
  
)



write.csv(
  
  por_error,
  
  file.path(
    
    ruta_salida,
    
    "03_resultados_por_error_v4.csv"
    
  ),
  
  row.names = FALSE
  
)



write.csv(
  
  por_tamano,
  
  file.path(
    
    ruta_salida,
    
    "04_resultados_por_tamano_v4.csv"
    
  ),
  
  row.names = FALSE
  
)



write.csv(
  
  por_dim,
  
  file.path(
    
    ruta_salida,
    
    "05_resultados_por_dimensionalidad_v4.csv"
    
  ),
  
  row.names = FALSE
  
)



write.csv(
  
  tamano_error,
  
  file.path(
    
    ruta_salida,
    
    "06_resultados_tamano_error_v4.csv"
    
  ),
  
  row.names = FALSE
  
)



write.csv(
  
  dim_error,
  
  file.path(
    
    ruta_salida,
    
    "07_resultados_dimensionalidad_error_v4.csv"
    
  ),
  
  row.names = FALSE
  
)



write.csv(
  
  errores_CHull_fp,
  
  file.path(
    
    ruta_salida,
    
    "08_errores_CHull_fp_v4.csv"
    
  ),
  
  row.names = FALSE
  
)



write.csv(
  
  errores_CHull_S,
  
  file.path(
    
    ruta_salida,
    
    "09_errores_CHull_S_v4.csv"
    
  ),
  
  row.names = FALSE
  
)



write.csv(
  
  diagnostico_ruido,
  
  file.path(
    
    ruta_salida,
    
    "10_diagnostico_ruido_v4.csv"
    
  ),
  
  row.names = FALSE
  
)



# =====================================================================
# 28. MOSTRAR RESULTADOS
# =====================================================================

cat("\n\n=============================================\n")
cat("RESULTADOS GENERALES\n")
cat("=============================================\n\n")


print(
  validacion,
  row.names = FALSE
)



cat("\n\n=============================================\n")
cat("RESULTADOS POR ERROR\n")
cat("=============================================\n\n")


print(
  por_error,
  row.names = FALSE
)



cat("\n\n=============================================\n")
cat("RESULTADOS POR TAMAÑO\n")
cat("=============================================\n\n")


print(
  por_tamano,
  row.names = FALSE
)



cat("\n\n=============================================\n")
cat("RESULTADOS POR DIMENSIONALIDAD\n")
cat("=============================================\n\n")


print(
  por_dim,
  row.names = FALSE
)



cat("\n\n=============================================\n")
cat("DIAGNOSTICO DEL ERROR\n")
cat("=============================================\n\n")


print(
  diagnostico_ruido,
  row.names = FALSE
)



# =====================================================================
# 29. GRAFICO: ACIERTO SEGUN NIVEL DE ERROR
# =====================================================================

matplot(
  
  por_error$error_objetivo,
  
  cbind(
    
    por_error$CHull_fp,
    
    por_error$CHull_S
    
  ),
  
  type = "b",
  
  pch = c(
    16,
    17
  ),
  
  lty = c(
    1,
    2
  ),
  
  ylim = c(
    0,
    100
  ),
  
  xlab =
    "Nivel nominal de error (%)",
  
  ylab =
    "Selección correcta (%)",
  
  main =
    "CHull: recuperación de la dimensionalidad Tucker3"
  
)



legend(
  
  "bottomleft",
  
  legend = c(
    
    "CHull-fp",
    
    "CHull-P+Q+R"
    
  ),
  
  pch = c(
    16,
    17
  ),
  
  lty = c(
    1,
    2
  ),
  
  bty = "n"
  
)



# =====================================================================
# 30. GRAFICO:
# ERROR NOMINAL VS ERROR OBSERVADO MEDIO
# =====================================================================

plot(
  
  diagnostico_ruido$error_objetivo,
  
  diagnostico_ruido$error_observado_medio,
  
  type = "b",
  
  pch = 16,
  
  xlab =
    "Error nominal (%)",
  
  ylab =
    "Error observado medio (%)",
  
  main =
    "Error nominal vs. error observado"
  
)


abline(
  
  a = 0,
  
  b = 1,
  
  lty = 2
  
)



# =====================================================================
# 31. COMPROBACION DE LA IDENTIDAD DE SUMAS DE CUADRADOS
# =====================================================================

resultados$error_identidad_SS <- abs(
  
  resultados$SS_X -
    
    (
      resultados$SS_senal +
        
        resultados$SS_ruido +
        
        resultados$termino_cruzado
    )
  
)



cat("\n\n=============================================\n")
cat("COMPROBACION NUMERICA DE SS\n")
cat("=============================================\n\n")


cat(
  
  "Mayor error absoluto de la identidad:\n",
  
  max(
    resultados$error_identidad_SS,
    na.rm = TRUE
  ),
  
  "\n"
  
)



# =====================================================================
# 32. GUARDAR NUEVAMENTE EL ARCHIVO PRINCIPAL
#     INCLUYENDO LA COMPROBACION DE SS
# =====================================================================

write.csv(
  
  resultados,
  
  file.path(
    
    ruta_salida,
    
    "01_resultados_MC_CHull_Tucker3_v4.csv"
    
  ),
  
  row.names = FALSE
  
)



# =====================================================================
# 33. FINAL
# =====================================================================

cat("\n\n=============================================\n")
cat("MONTE CARLO FINALIZADO CORRECTAMENTE\n")
cat("=============================================\n\n")


cat(
  
  "Numero total de simulaciones:",
  
  nrow(resultados),
  
  "\n"
  
)


cat(
  
  "CHull-fp:",
  
  aciertos_fp,
  
  "aciertos y",
  
  errores_fp,
  
  "errores =",
  
  round(
    porcentaje_fp,
    2
  ),
  
  "%\n"
  
)


cat(
  
  "CHull-S:",
  
  aciertos_S,
  
  "aciertos y",
  
  errores_S,
  
  "errores =",
  
  round(
    porcentaje_S,
    2
  ),
  
  "%\n"
  
)


cat("\nArchivos guardados en:\n")

cat(
  ruta_salida,
  "\n"
)

