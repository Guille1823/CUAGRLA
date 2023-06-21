#' Cuadrado Grecolatino
#'
#' @param Respuesta (string) nombre de la variable respuesta
#' @param Tratamiento (string) nombre de la variable tratamiento
#' @param Filas (string) nombre de la variable filas
#' @param Columnas (string) nombre de la variable columnas
#' @param LetrasGR (string) nombre de la variable letras griegas
#' @param data (\code{data.frame}) Tabla de datos en formato largo con los datos del tratamiento, filas, columnas, letras griegas y respuesta.
#'
#' @return Devuelve una tabla en formato fata.frame con los calculos orrespondientes a una tabla Anova
#' @export
#Funcion
CGL <- function(Respuesta, Tratamiento, Filas, Columnas, LetrasGR, data) {
  # y^2/k^2
  sumtotalRESPUESTA <- sum(Respuesta)
  k <- nlevels(Tratamiento)
  s <- sumtotalRESPUESTA^2 / k^2

  # Calculamos grados de libertad, suma de cuadrado, cuadrado medio, F y significancia
  # de tratamiento, filas, columnas, letras griegas, error y total

  # TRATAMIENTO
  gl_tto <- k - 1
  sum_tto <- tapply(Respuesta, INDEX = Tratamiento, FUN = sum)
  n_tto <- tapply(Respuesta, INDEX = Tratamiento, FUN = length)
  sc_tto <- sum(sum_tto^2 / n_tto) - s
  cm_tto <- sc_tto / gl_tto

  # FILAS
  gl_filas <- k - 1
  sum_filas <- tapply(Respuesta, INDEX = Filas, FUN = sum)
  n_filas <- tapply(Respuesta, INDEX = Filas, FUN = length)
  sc_filas <- sum(sum_filas^2 / n_filas) - s
  cm_filas <- sc_filas / gl_filas

  # COLUMNAS
  gl_columnas <- k - 1
  sum_columnas <- tapply(Respuesta, INDEX = Columnas, FUN = sum)
  n_columnas <- tapply(Respuesta, INDEX = Columnas, FUN = length)
  sc_columnas <- sum(sum_columnas^2 / n_columnas) - s
  cm_columnas <- sc_columnas / gl_columnas

  # LETRAS GRIEGAS
  gl_letrasgr <- k - 1
  sum_letrasgr <- tapply(Respuesta, INDEX = LetrasGR, FUN = sum)
  n_letrasgr <- tapply(Respuesta, INDEX = LetrasGR, FUN = length)
  sc_letrasgr <- sum(sum_letrasgr^2 / n_letrasgr) - s
  cm_letrasgr <- sc_letrasgr / gl_letrasgr

  # TOTAL
  gl_total <- k^2 - 1
  sc_total <- sum(Respuesta^2) - s

  # ERROR
  gl_error <- (k - 1) * (k - 3)
  sc_error <- sc_total - sc_tto - sc_filas - sc_columnas - sc_letrasgr
  cm_error <- sc_error / gl_error

  # F
  F_tto <- cm_tto / cm_error
  F_filas <- cm_filas / cm_error
  F_columnas <- cm_columnas / cm_error
  F_letrasgr <- cm_letrasgr / cm_error

  # SIGNIFICANCIA
  p_value_tto <- pf(F_tto, gl_tto, gl_error, lower.tail = FALSE)
  p_value_filas <- pf(F_filas, gl_filas, gl_error, lower.tail = FALSE)
  p_value_columnas <- pf(F_columnas, gl_columnas, gl_error, lower.tail = FALSE)
  p_value_letrasgr <- pf(F_letrasgr, gl_letrasgr, gl_error, lower.tail = FALSE)

  # Crear la tabla de resultados
  tabla <- data.frame(
    Fuente = c("Tratamiento", "Filas", "Columnas", "Letras Griegas", "Error", "Total"),
    GL = c(gl_tto, gl_filas, gl_columnas, gl_letrasgr, gl_error, gl_total),
    SC = c(sc_tto, sc_filas, sc_columnas, sc_letrasgr, sc_error, sc_total),
    CM = c(cm_tto, cm_filas, cm_columnas, cm_letrasgr, cm_error, NA),
    F = c(F_tto, F_filas, F_columnas, F_letrasgr, NA, NA),
    `p-value` = c(p_value_tto, p_value_filas, p_value_columnas, p_value_letrasgr, NA, NA)
  )

  return(tabla)
}
