
library(ggplot2)
library(viridis)  # per scale di colori migliori
library(mvtnorm)
library(gridExtra)
library(coda)
library(grid)
set.seed(359445)  


# Carica il file .RData
load("C:/Users/sierv/Desktop/homework 2/data_for_homework2.RData")

#imposto cartella di lavoro per salvataggio dati
setwd("C:/Users/sierv/Desktop/homework 2/")

# Verifica gli oggetti caricati
print("Oggetti nell'ambiente:")
print(ls())

# Esplora la struttura di data_clean
print("Struttura di data_clean:")
str(data_clean)

# Visualizza le prime righe
print("Prime righe di data_clean:")
head(data_clean)

# Verifica i nomi delle colonne
print("Nomi delle colonne:")
names(data_clean)


inizio <- Sys.time()

# Grafico dei dati con scala di colori basata su log_PM10
# Scala invertita: colori scuri = valori alti, colori chiari = valori bassi
p1 <- ggplot(data_clean, aes(x = X, y = Y, color = log_PM10)) +
  geom_point(size = 2, alpha = 0.7) +  # Punti piccoli e leggermente trasparenti
  scale_color_viridis(
    option = "viridis",  # Scala viridis classica
    name = "log(PM10)",
    direction = -1  # Inverte la scala: scuro = alto, chiaro = basso
  ) +
  labs(
    title = "Distribuzione spaziale del log(PM10) nelle città europee",
    subtitle = "Dati rilevati il 3 dicembre 2025 alle ore 12:00",
    x = "Coordinata X",
    y = "Coordinata Y"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right"
  ) +
  coord_fixed()  # mantiene le proporzioni corrette delle coordinate

# Visualizza il grafico
print(p1)

# Salva il grafico
ggsave("punto1_mappa_pm10.png", p1, width = 10, height = 8, dpi = 300)



# Punto 2

n_sub <- 150   

N <- nrow(data_clean)

#Estraggo 150 numeri casuali evitando ripetizioni in modo da estrarre 150 righe
# dal dataset
idx_sub <- sample(1:N, size = n_sub, replace = FALSE) 
data_sub <- data_clean[idx_sub, ] #estrazione righe sfruttando i numeri estratti

# Estraggo quello che useremo nel modello
y  <- data_sub$log_PM10
S  <- as.matrix(data_sub[, c("X", "Y")])   # coordinate (s1_i, s2_i)
n  <- length(y)

#Verifico che il procedimento di estrazione sia andato a buon fine
dim(data_sub)
head(data_sub)



#Creo matrice delle distanze euclidee tra i punti (s1,s2) =(x,y) per cercare 
# distanza massima per la prior 
coordinate_xy <- data_sub[, c("X","Y")] #estrae tutte le righe ma solo le colonne
# X,Y

#Salvo le coordinate in una matrice
S <- as.matrix(coordinate_xy) #matrice n_sub * 2

#calcolo le distante tra tutte le possibili coppie di righe (che sono le coppie
# di coordinate) della matrice S e le metto in una matrice
distanze <- dist(S)
matrice_dist <- as.matrix(distanze) #matrice n_sub * n_sub

#Cerco la distanza massima minima selezionando solo la parte superiore della matrice
# intanto la diagonale è nulla e la matrice è simmetrice
dist_valori <- matrice_dist[upper.tri(matrice_dist)] #upper.tri crea una matrice di True e False che funzionano
# come criterio di selezione
dist_min  <- min(dist_valori)
dist_max  <- max(dist_valori)


#MCMC
B <- 125000
B_burn <- 25000
#costruzione matrice X
X <- cbind(1, data_sub$X, data_sub$Y) #comando cbind unisce le colonne e crea
#matrice

#funzioni per full coditional

#Costruzione matrice di correlazione spaziale R_phi
costruisci_R_phi <- function(phi, matrice_dist)
{
  return(exp(-phi * matrice_dist))
}

campiona_w <- function(y, X, beta_corrente, tau2_corrente, sigma2_corrente,
                       phi_corrente, matrice_dist)
{
  n <- length(y)
  
  # Matrice di covarianza spaziale di w
  R_phi     <- costruisci_R_phi(phi_corrente, matrice_dist)
  Sigma_w   <- sigma2_corrente * R_phi  # Σ_w = σ^2 R_phi
  
  # Prima inverto Σ_w
  Sigma_w_inv <- solve(Sigma_w)
  precisione_w <- (1 / tau2_corrente) * diag(n) + Sigma_w_inv
  
  # Varianza a posteriori V_w = precisione^{-1}
  V_w <- chol2inv(chol(precisione_w))
  
  # Media a posteriori μ_w = V_w (1/τ^2) (y - Xβ)
  residuo <- y - X %*% beta_corrente
  mu_w    <- V_w %*% ((1 / tau2_corrente) * residuo)
  
  # Campionamento da N(μ_w, V_w)
  #as.numeric mi assicura che sia un vettore
  return(as.numeric(rmvnorm(1, mean = as.vector(mu_w), sigma = V_w)))
}

campiona_beta <- function(y, w_corrente, X, tau2_corrente) {
  n <- length(y)
  p <- ncol(X) # qui p = 3
  
  # y* = y - w
  y_stella <- y - w_corrente
  
  # Precisione di beta: (1/τ^2) X'X + (1/10^4) I_3
  precisione_beta <- (1 / tau2_corrente) * t(X) %*% X +
    (1 / 1e4) * diag(p)
  
  # Varianza a posteriori
  V_beta <- chol2inv(chol(precisione_beta))
  
  # Media a posteriori μ_beta = V_beta (1/τ^2) X' y*
  mu_beta <- V_beta %*% ((1 / tau2_corrente) * t(X) %*% y_stella)
  
  # Campionamento da N(μ_beta, V_beta)
  return(as.numeric(rmvnorm(1, mean = as.vector(mu_beta), sigma = V_beta)))
}

# Campionamento da Inverse-Gamma IG(a, b) (funzione di appoggio)
campiona_invgamma <- function(forma, scala) {
  1 / rgamma(1, shape = forma, rate = scala)
}


campiona_tau2 <- function(y, X, beta_corrente, w_corrente) {
  n <- length(y)
  
  # residuo = y - Xbeta - w
  residuo <- y - X %*% beta_corrente - w_corrente
  
  # Parametri IG( A_tau, B_tau )
  A_tau <- 1 + n / 2
  B_tau <- 1 + 0.5 * sum(residuo^2)
  
  campiona_invgamma(forma = A_tau, scala = B_tau)
}

campiona_sigma2 <- function(w_corrente, phi_corrente, matrice_dist) {
  n <- length(w_corrente)
  
  # Matrice di correlazione R_phi
  R_phi <- costruisci_R_phi(phi_corrente, matrice_dist)
  
  # Inversa di R_phi
  R_inv <- chol2inv(chol(R_phi))
  
  # Parametri IG( A_sigma, B_sigma )
  A_sigma <- 1 + n / 2
  B_sigma <- 1 + 0.5 * as.numeric(t(w_corrente) %*% R_inv %*% w_corrente)
  
  return(campiona_invgamma(forma = A_sigma, scala = B_sigma))
}

#per phi non creiamo la funzione di campinamento ma solo una funzione che
# calcola la denità di un punto, dato che il campionamento avviene tramite
# Metropolis

log_phi <- function(phi, w_corrente, sigma2_corrente, matrice_dist)
{
  #non viene fatto nessun controlo sul valore di phi in entrata in quanto 
  # è gia controllato prima
  
  R_phi <- costruisci_R_phi(phi, matrice_dist)
  
  chol_R <- chol(R_phi) #decomposizione Cholesky
  
  # log |R_phi|
  log_det_R <- 2 * sum(log(diag(chol_R)))
  
  #Come da teoria risolvo Lz = w tramite sostituione in avanti
  L <- t(chol_R)
  z <- forwardsolve(L,w_corrente)
  quad <- sum(z^2)
  
  return(-0.5 * log_det_R - 0.5 * quad / sigma2_corrente)
}


#inizializzo 

#beta
mod_lm <- lm(y ~ data_sub$X + data_sub$Y)
beta_0 <- coef(mod_lm)
#matrice che salva le iterazioni sulle righe
beta_catena <- matrix(NA, nrow = B, ncol = 3)
beta_catena[1,] <- beta_0

#W
w_0 <- rep(0,n_sub) #vettore di 0
w_catena <- matrix(NA,nrow = B, ncol = n_sub)
w_catena[1,] <- w_0

#tau
residui_lm <- resid(mod_lm)
tau2_0 <- var(residui_lm) * 0.5
tau2_catena <- matrix(NA,nrow = B, ncol = 1)
tau2_catena[1,] <- tau2_0

#sigma^2
sigma2_0 <- var(residui_lm) * 0.5
sigma2_catena <- matrix(NA,nrow = B, ncol = 1)
sigma2_catena[1,] <- sigma2_0

#Phi
#limiti intervallo originali
phi_min <- 3 / dist_max
phi_max <- 3 / dist_min

#Trasformo in scala logaritmica
log_phi_min <- log(phi_min)
log_phi_max <- log(phi_max)

# Inizializzo phi in scala LOG
log_phi_0 <- (log_phi_min + log_phi_max) / 2  # Punto medio in scala log
log_phi_catena <- rep(NA, B)
log_phi_catena[1] <- log_phi_0

# Varianza proposta per log(phi)
range_log_phi <- log_phi_max - log_phi_min
varianza_log <- (range_log_phi / 10)^2  # Proposta = 10% del range log

alfa_desiderato <- 0.25
alfa_catena <- rep(NA,B)

alfa_catena[1] <- 0
c <- 10 #tasso di aggiornamento

#variabili per proposta adattiva modificati rispetto a prima
A <- 2
D <- 100

#iterazione di inizio e fine adattamento
iter_inizio_adatt <- 50
iter_fine_adatt <- B_burn - 500


for (i in 2:B)
{
  beta_corrente <- beta_catena[i-1,]
  w_corrente <- w_catena[i-1,]
  tau2_corrente <- tau2_catena[i-1,]
  sigma2_corrente <- sigma2_catena[i-1,]
  log_phi_corrente <- log_phi_catena[i-1]
  phi_corrente <- exp(log_phi_corrente)  # Trasforma per usare nelle funzioni
  # di calcolo della densità e di campionamento
  
  #campiono dalle full conditional
  beta_corrente <- campiona_beta(y, w_corrente, X, tau2_corrente)
  beta_catena[i,] <- beta_corrente
  
  w_corrente <- campiona_w(y, X, beta_corrente, tau2_corrente, sigma2_corrente,
                           phi_corrente, matrice_dist)
  w_catena[i,] <- w_corrente
  
  tau2_corrente <- campiona_tau2(y, X, beta_corrente, w_corrente)
  tau2_catena[i,] <- tau2_corrente
  
  sigma2_corrente <- campiona_sigma2(w_corrente, phi_corrente, matrice_dist)
  sigma2_catena[i,] <- sigma2_corrente
  
  #Metropolis per phi

  # Proposta in scala LOG
  log_phi_proposta <- rnorm(1, log_phi_corrente, sqrt(varianza_log))
  
  
  # Controlla se log_phi proposto è nel range ammissibile
  if (log_phi_proposta < log_phi_min || log_phi_proposta > log_phi_max)
  {
    # Rifiuta automaticamente
    alfa_catena[i] <- 0
    log_phi_catena[i] <- log_phi_corrente
  }
  else
  {
    # Trasforma in scala originale per calcolare densità
    phi_proposta <- exp(log_phi_proposta)
    
    # Calcola log-densità a posteriori
    # IMPORTANTE: Aggiungiamo Jacobiano della trasformazione sommandolo in fondo
    log_post_proposta <- log_phi(phi_proposta, w_corrente, sigma2_corrente, 
                                 matrice_dist) + log_phi_proposta
    
    log_post_corrente <- log_phi(phi_corrente, w_corrente, sigma2_corrente, 
                                 matrice_dist) + log_phi_corrente
    
    # Calcola tasso di accettazione
    log_alpha <- log_post_proposta - log_post_corrente
    alfa_corrente <- min(1, exp(log_alpha))
    
    alfa_catena[i] <- alfa_corrente
    
    # Accetta/rifiuta
    p <- runif(1,0,1)
    if (p <= alfa_corrente)
    {
      log_phi_catena[i] <- log_phi_proposta
    }
    else
    {
      log_phi_catena[i] <- log_phi_corrente
    }
  }
  
  #Aggiornamento varianza 
  if (i %% c == 0 && i >= iter_inizio_adatt && i <= iter_fine_adatt)
  {
    alfa_media_recente <- mean(alfa_catena[ (i-c+1):i])
    gamma_i <- A / (D + i)
    
    log_varianza_new <- log(varianza_log) + gamma_i * (alfa_media_recente - alfa_desiderato)
    varianza_log <- exp(log_varianza_new)
    
    # Limiti di sicurezza
    varianza_log <- min(varianza_log, (range_log_phi / 2)^2)
    varianza_log <- max(varianza_log, (range_log_phi / 100)^2)
    
    #Controllo
    if (i %% 200 == 0)
    {
      cat(sprintf("Iter %d: alfa=%.3f, var_log=%.4f, sd_log=%.3f\n", 
                  i, alfa_media_recente, varianza_log, sqrt(varianza_log)))
    }
  }
  
  
  #controllo eseguzione
  if (i %% 2000 == 0 && i > B_burn)
  {
    cat(sprintf("Iter %d: ", i ))
  }
  
}


print(Sys.time() - inizio)

# Trasforma log(phi) in phi
phi_catena <- exp(log_phi_catena) #vettore non matrice


# Rimuoviamo i primi B_burn campioni 
beta_post <- beta_catena[(B_burn + 1):B, ]
w_post <- w_catena[(B_burn + 1):B, ]
tau2_post <- tau2_catena[(B_burn + 1):B, ]
sigma2_post <- sigma2_catena[(B_burn + 1):B, ]
phi_post <- phi_catena[(B_burn + 1):B]
alfa_post <- alfa_catena[(B_burn + 1):B]

n_post <- nrow(beta_post)  # numero campioni post burn-in
cat("Campioni dopo burn-in:", n_post, "\n")

# Thinning: dato che emerge alta autorcorrelazione in tutti i parametri
# procediamo a fare del thin, ovvero a prendere un campione ogni k

thin <- 50
# indici delle iterazioni che tieni dopo il burn-in
indici_post <- seq(from = B_burn + 1, to = B, by = thin)

# catene thinnate
beta_post  <- beta_catena[indici_post, ]
w_post     <- w_catena[indici_post, ]
#metto drop = false per non far convertire la matrice in vettore
tau2_post  <- tau2_catena[indici_post, , drop = FALSE]
sigma2_post <- sigma2_catena[indici_post, , drop = FALSE]
phi_post   <- phi_catena[indici_post]
alfa_post  <- alfa_catena[indici_post]

n_post <- length(indici_post)
cat("Campioni dopo burn-in e thinning:", n_post, "\n")

# ANALISI DEI RISULTATI

#Creo una funzione semplice che mi permette di avere una visione rapida 
# dei primi risultati

riassunto_posteriori <- function(catena_post)
{
  media <- mean(catena_post)
  dev_standard <- sd(catena_post)
  
  #intervallo credibilità 95%
  q_025 <- quantile(catena_post, probs = 0.025)
  q_975 <- quantile(catena_post, probs = 0.975)
  
  #restituisco riga con i risultati assegnando nome per facile identificazione
  return(c(Media = media,
           DevStd = dev_standard,
           Q2.5 = q_025,
           Q97.5 = q_975))

}
  
# beta_post è una matrice: colonne = (mu, beta1, beta2)
# creo una tabella riassuntiva riga-per-parametro
tabella_beta <- rbind(
  mu    = riassunto_posteriori(beta_post[, 1]),
  beta1 = riassunto_posteriori(beta_post[, 2]),
  beta2 = riassunto_posteriori(beta_post[, 3])
)


print(round(tabella_beta, 4))
  
tabella_varianze_phi <- rbind(
  tau2   = riassunto_posteriori(tau2_post),
  sigma2 = riassunto_posteriori(sigma2_post),
  phi    = riassunto_posteriori(phi_post)
)


print(round(tabella_varianze_phi, 4))  


# Diagonosi grafiche, uso di IA per migliorare grafici
# 3.1 TRACEPLOTS 

# Prepariamo i dati per i traceplots
# Creiamo un dataframe con tutte le iterazioni post burn-in e thinning

iterazioni <- 1:n_post

# DataFrame per beta (mu, beta1, beta2)
df_beta <- data.frame(
  Iterazione = rep(iterazioni, 3),
  Parametro = rep(c("mu", "beta1", "beta2"), each = n_post),
  Valore = c(beta_post[, 1], beta_post[, 2], beta_post[, 3])
)

# Traceplot per beta
p_trace_beta <- ggplot(df_beta, aes(x = Iterazione, y = Valore)) +
  geom_line(color = "steelblue", alpha = 0.7) +
  facet_wrap(~ Parametro, scales = "free_y", ncol = 1) +
  labs(
    title = "Traceplots dei parametri beta",
    subtitle = "Dopo burn-in e thinning",
    x = "Iterazione",
    y = "Valore del parametro"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold")
  )

print(p_trace_beta)
ggsave("punto3_traceplots_beta.png", p_trace_beta, width = 10, height = 8, dpi = 300)


# DataFrame per varianze e phi
df_var_phi <- data.frame(
  Iterazione = rep(iterazioni, 3),
  Parametro = rep(c("tau2", "sigma2", "phi"), each = n_post),
  Valore = c(tau2_post, sigma2_post, phi_post)
)

# Traceplot per tau2, sigma2, phi
p_trace_var <- ggplot(df_var_phi, aes(x = Iterazione, y = Valore)) +
  geom_line(color = "darkgreen", alpha = 0.7) +
  facet_wrap(~ Parametro, scales = "free_y", ncol = 1) +
  labs(
    title = "Traceplots dei parametri di varianza e phi",
    subtitle = "Dopo burn-in e thinning",
    x = "Iterazione",
    y = "Valore del parametro"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    strip.text = element_text(face = "bold")
  )

print(p_trace_var)
ggsave("punto3_traceplots_varianze.png", p_trace_var, width = 10, height = 8, dpi = 300)
  

# Densità per beta
p_dens_beta <- ggplot(df_beta, aes(x = Valore, fill = Parametro)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ Parametro, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Densità a posteriori dei parametri beta",
    x = "Valore",
    y = "Densità"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

print(p_dens_beta)
ggsave("punto3_densita_beta.png", p_dens_beta, width = 10, height = 4, dpi = 300)


# Densità per varianze e phi
p_dens_var <- ggplot(df_var_phi, aes(x = Valore, fill = Parametro)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ Parametro, scales = "free", ncol = 3) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    title = "Densità a posteriori di tau2, sigma2 e phi",
    x = "Valore",
    y = "Densità"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

print(p_dens_var)
ggsave("punto3_densita_varianze.png", p_dens_var, width = 10, height = 4, dpi = 300)


# grafico riassuntivo, medie e IC


# Creiamo un dataframe con tutti i 6 parametri e le loro statistiche
df_tutti_parametri <- data.frame(
  Parametro = c("mu", "beta1", "beta2", "tau2", "sigma2", "phi"),
  Media = c(
    mean(beta_post[, 1]),
    mean(beta_post[, 2]),
    mean(beta_post[, 3]),
    mean(tau2_post),
    mean(sigma2_post),
    mean(phi_post)
  ),
  Q025 = c(
    quantile(beta_post[, 1], 0.025),
    quantile(beta_post[, 2], 0.025),
    quantile(beta_post[, 3], 0.025),
    quantile(tau2_post, 0.025),
    quantile(sigma2_post, 0.025),
    quantile(phi_post, 0.025)
  ),
  Q975 = c(
    quantile(beta_post[, 1], 0.975),
    quantile(beta_post[, 2], 0.975),
    quantile(beta_post[, 3], 0.975),
    quantile(tau2_post, 0.975),
    quantile(sigma2_post, 0.975),
    quantile(phi_post, 0.975)
  )
)

# Ordiniamo i parametri dal basso verso l'alto per migliore visualizzazione
# Così il primo parametro (mu) sarà in alto
df_tutti_parametri$Parametro <- factor(
  df_tutti_parametri$Parametro,
  levels = c("phi", "sigma2", "tau2", "beta2", "beta1", "mu")
)

# Grafico con parametri sull'asse Y e valori sull'asse X
p_stime_complete <- ggplot(df_tutti_parametri, aes(x = Media, y = Parametro)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", 
             alpha = 0.5, size = 0.8) +
  geom_segment(aes(x = Q025, xend = Q975, y = Parametro, yend = Parametro),
               color = "steelblue", size = 1.5, alpha = 0.7) +
  geom_point(size = 4, color = "darkblue") +
  labs(
    title = "Stime a posteriori e intervalli di credibilità al 95%",
    subtitle = "Media (punto) e IC 95% (segmento) per ciascun parametro",
    x = "Valore del parametro",
    y = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.y = element_text(face = "bold", size = 11),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank()
  )

print(p_stime_complete)
ggsave("punto3_stime_complete.png", p_stime_complete, 
       width = 10, height = 6, dpi = 300)


#funzione di caovarianza

# Calcoliamo la funzione di covarianza per diverse distanze
# Usiamo distanze da 0 fino alla distanza massima osservata
# Fissiamo la lunghezza di h = 100
h_seq <- seq(0, dist_max, length.out = 100)

# Per ogni campione a posteriori, calcoliamo C(h; phi) per ogni h
# C(h; phi) = sigma^2 * exp(-h * phi)

# Matrice per salvare i valori: righe = campioni post, colonne = valori di h
C_h_matrix <- matrix(NA, nrow = n_post, ncol = length(h_seq))

for (i in 1:n_post) {
  sigma2_i <- sigma2_post[i]
  phi_i <- phi_post[i]
  C_h_matrix[i, ] <- sigma2_i * exp(-h_seq * phi_i)
}

# Calcoliamo media e intervalli di credibilità per ogni h, quindi
# su ogni colonna. In pratica per ogni h prendiamo gli n_post campioni
# e calcoliamo media e IC
C_h_media <- apply(C_h_matrix, 2, mean)
C_h_q025 <- apply(C_h_matrix, 2, quantile, probs = 0.025)
C_h_q975 <- apply(C_h_matrix, 2, quantile, probs = 0.975)

# DataFrame per il grafico
df_C_h <- data.frame(
  h = h_seq,
  Media = C_h_media,
  Q025 = C_h_q025,
  Q975 = C_h_q975
)

# Grafico della funzione di covarianza
p_C_h <- ggplot(df_C_h, aes(x = h)) +
  geom_ribbon(aes(ymin = Q025, ymax = Q975), fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = Media), color = "darkblue", size = 1) +
  labs(
    title = "Funzione di covarianza spaziale C(h; phi)",
    subtitle = "Media a posteriori e intervallo di credibilità al 95%",
    x = "Distanza spaziale (h)",
    y = "C(h; phi) = sigma^2 * exp(-h * phi)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_C_h)
ggsave("punto3_covarianza_spaziale.png", p_C_h, width = 10, height = 6, dpi = 300)




# Rappresentazione spaziale W

# Calcoliamo statistiche a posteriori di W per ogni punto spaziale, quindi sulle 
# colonne del campione a posteriori
# Media a posteriori
W_media <- apply(w_post, 2, mean)

# Intervalli di credibilità al 95%
W_q025 <- apply(w_post, 2, quantile, probs = 0.025)
W_q975 <- apply(w_post, 2, quantile, probs = 0.975)

# Creiamo un dataframe con le coordinate e le statistiche di W
df_W <- data.frame(
  X = data_sub$X,
  Y = data_sub$Y,
  W_media = W_media,
  W_q025 = W_q025,
  W_q975 = W_q975,
  log_PM10 = y  # Aggiungiamo anche i valori osservati
)


# Grafico della media a posteriori di W 
p_W_media <- ggplot(df_W, aes(x = X, y = Y, color = W_media)) +
  geom_point(size = 2, alpha = 0.7) +  # Punti piccoli e leggermente trasparenti
  scale_color_viridis(
    option = "viridis",  # Scala viridis classica
    name = "W (media)",
    direction = -1  # Inverte la scala: scuro = alto, chiaro = basso
  ) +
  labs(
    title = "Media a posteriori del processo spaziale W(s)",
    subtitle = "Effetto spaziale non spiegato dalle covariate",
    x = "Coordinata X",
    y = "Coordinata Y"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right"
  ) +
  coord_fixed()  # mantiene le proporzioni corrette delle coordinate

print(p_W_media)
ggsave("punto3_W_media.png", p_W_media, width = 10, height = 8, dpi = 300)


# Grafico del limite inferiore IC 95%
p_W_q025 <- ggplot(df_W, aes(x = X, y = Y, color = W_q025)) +
  geom_point(size = 2, alpha = 0.7) +  # Punti piccoli e leggermente trasparenti
  scale_color_viridis(
    option = "viridis",  # Scala viridis classica
    name = "W (Q2.5%)",
    direction = -1  # Inverte la scala: scuro = alto, chiaro = basso
  ) +
  labs(
    title = "Limite inferiore IC 95% di W(s)",
    subtitle = "Quantile 2.5% della distribuzione a posteriori",
    x = "Coordinata X",
    y = "Coordinata Y"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right"
  ) +
  coord_fixed()  # mantiene le proporzioni corrette delle coordinate

print(p_W_q025)
ggsave("punto3_W_q025.png", p_W_q025, width = 10, height = 8, dpi = 300)


# Grafico del limite superiore IC 95%
p_W_q975 <- ggplot(df_W, aes(x = X, y = Y, color = W_q975)) +
  geom_point(size = 2, alpha = 0.7) +  # Punti piccoli e leggermente trasparenti
  scale_color_viridis(
    option = "viridis",  # Scala viridis classica
    name = "W (Q97.5%)",
    direction = -1  # Inverte la scala: scuro = alto, chiaro = basso
  ) +
  labs(
    title = "Limite superiore IC 95% di W(s)",
    subtitle = "Quantile 97.5% della distribuzione a posteriori",
    x = "Coordinata X",
    y = "Coordinata Y"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    legend.position = "right"
  ) +
  coord_fixed()  # mantiene le proporzioni corrette delle coordinate

print(p_W_q975)
ggsave("punto3_W_q975.png", p_W_q975, width = 10, height = 8, dpi = 300)

# Combina e salva i tre grafici in un'unica figura
png("punto3_W_combined.png", width = 5400, height = 1800, res = 300)
grid.arrange(p_W_media, p_W_q025, p_W_q975, 
             ncol = 3,
             top = textGrob("Distribuzione spaziale del processo W(s)",
                            gp = gpar(fontsize = 14, fontface = "bold")))
dev.off()

# Creiamo un grafico affiancato per confrontare W_media e log_PM10

# Grafico di Y osservato (per confronto)
p_Y_obs <- ggplot(df_W, aes(x = X, y = Y, color = log_PM10)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_viridis(
    option = "viridis",
    name = "log(PM10)",
    direction = -1
  ) +
  labs(
    title = "Valori osservati: log(PM10)",
    x = "Coordinata X",
    y = "Coordinata Y"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  ) +
  coord_fixed()

# Confronto affiancato
p_confronto <- grid.arrange(p_Y_obs, p_W_media, ncol = 2)

ggsave("punto3_confronto_Y_W.png", p_confronto, width = 16, height = 7, dpi = 300)
  


# Funzione per calcolare ESS usando il pacchetto coda
# L'ESS misura il numero effettivo di campioni indipendenti nella catena

# Convertiamo le catene in oggetti mcmc per usare le funzioni di coda
beta_mcmc <- as.mcmc(beta_post)
tau2_mcmc <- as.mcmc(tau2_post)
sigma2_mcmc <- as.mcmc(sigma2_post)
phi_mcmc <- as.mcmc(phi_post)

# Calcoliamo ESS per ogni parametro
ess_mu <- effectiveSize(beta_mcmc[, 1])
ess_beta1 <- effectiveSize(beta_mcmc[, 2])
ess_beta2 <- effectiveSize(beta_mcmc[, 3])
ess_tau2 <- effectiveSize(tau2_mcmc)
ess_sigma2 <- effectiveSize(sigma2_mcmc)
ess_phi <- effectiveSize(phi_mcmc)

# Creiamo una tabella riassuntiva con ESS
df_ess <- data.frame(
  Parametro = c("mu", "beta1", "beta2", "tau2", "sigma2", "phi"),
  N_campioni = rep(n_post, 6),  # Numero totale di campioni post burn-in e thinning
  ESS = c(ess_mu, ess_beta1, ess_beta2, ess_tau2, ess_sigma2, ess_phi),
  Efficienza = c(ess_mu / n_post, ess_beta1 / n_post, ess_beta2 / n_post,
                 ess_tau2 / n_post, ess_sigma2 / n_post, ess_phi / n_post) * 100
)

# Arrotondiamo per migliore leggibilità
df_ess$ESS <- round(df_ess$ESS, 1)
df_ess$Efficienza <- round(df_ess$Efficienza, 2)

# Stampiamo la tabella
cat("\n=== EFFECTIVE SAMPLE SIZE (ESS) ===\n\n")
cat("L'ESS rappresenta il numero effettivo di campioni indipendenti.\n")
cat("Un ESS più alto indica minore autocorrelazione nella catena.\n")
cat("L'efficienza è il rapporto ESS/N_campioni in percentuale.\n\n")

print(df_ess, row.names = FALSE)

# Grafici ACF per i 6 parametri - divisi in 2 gruppi
# primo gruppo: mu, beta1, beta2
png("punto3_acf_beta.png", width = 3600, height = 1200, res = 300)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 2))

acf(beta_post[, 1], 
    main = "ACF - mu",
    col = "steelblue",
    lwd = 2,
    lag.max = 50)

acf(beta_post[, 2], 
    main = "ACF - beta1",
    col = "steelblue",
    lwd = 2,
    lag.max = 50)

acf(beta_post[, 3], 
    main = "ACF - beta2",
    col = "steelblue",
    lwd = 2,
    lag.max = 50)

dev.off()
cat("Primo grafico salvato: punto3_acf_beta.png\n")


# secondo gruppo: tau2, sigma2, phi
png("punto3_acf_varianze_phi.png", width = 3600, height = 1200, res = 300)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 2))

acf(tau2_post, 
    main = "ACF - tau2",
    col = "darkgreen",
    lwd = 2,
    lag.max = 50)

acf(sigma2_post, 
    main = "ACF - sigma2",
    col = "darkgreen",
    lwd = 2,
    lag.max = 50)

acf(phi_post, 
    main = "ACF - phi",
    col = "darkgreen",
    lwd = 2,
    lag.max = 50)

dev.off()
cat("Secondo grafico salvato: punto3_acf_varianze_phi.png\n")

cat("\nEntrambi i grafici ACF sono stati salvati con successo!\n")


#Punto 5
#MCMC marginalizzato
# B e B_burn già impostati

cat(" AVVIO MCMC SUL MODELLO MARGINALIZZATO\n")
inizio_marg <- Sys.time()

# funzioni per log-densità, non in forma chiusa

# Log-densità per tau2
log_dens_tau2_marg <- function(tau2, y, X, beta, sigma2, phi, matrice_dist) 
  {
  n <- length(y)
  
  # Costruisco Sigma_y = sigma2 * R_phi + tau2 * I_n
  R_phi <- costruisci_R_phi(phi, matrice_dist)
  Sigma_y <- sigma2 * R_phi + tau2 * diag(n)
  
  # Decomposizione di Cholesky per stabilità numerica
  chol_Sigma_y <- chol(Sigma_y)
  
  # Calcolo log determinante
  log_det_Sigma_y <- 2 * sum(log(diag(chol_Sigma_y)))
  
  # Calcolo forma quadratica
  residuo <- y - X %*% beta
  L <- t(chol_Sigma_y)
  z <- forwardsolve(L, residuo)
  quad <- sum(z^2)
  
  # Log-densità: -1/2 * log|Sigma_y| - 1/2 * r' Sigma_y^{-1} r - 1/tau2 - 2*log(tau2)
  log_dens <- -0.5 * log_det_Sigma_y - 0.5 * quad - 1/tau2 - 2 * log(tau2)
  
  return(log_dens)
}

# Log-densità per sigma2 
log_dens_sigma2_marg <- function(sigma2, y, X, beta, tau2, phi, matrice_dist) {
  n <- length(y)
  
  # Costruisco Sigma_y = sigma2 * R_phi + tau2 * I_n
  R_phi <- costruisci_R_phi(phi, matrice_dist)
  Sigma_y <- sigma2 * R_phi + tau2 * diag(n)
  
  # Decomposizione di Cholesky per stabilità numerica
  chol_Sigma_y <- chol(Sigma_y)
  
  # Calcolo log determinante
  log_det_Sigma_y <- 2 * sum(log(diag(chol_Sigma_y)))
  
  # Calcolo forma quadratica
  residuo <- y - X %*% beta
  L <- t(chol_Sigma_y)
  z <- forwardsolve(L, residuo)
  quad <- sum(z^2)
  
  log_dens <- -0.5 * log_det_Sigma_y - 0.5 * quad - 1/sigma2 - 2 * log(sigma2)
  
  return(log_dens)
}

# Log-densità per phi 
log_dens_phi_marg <- function(phi, y, X, beta, sigma2, tau2, matrice_dist) {
  n <- length(y)
  
  R_phi <- costruisci_R_phi(phi, matrice_dist)
  Sigma_y <- sigma2 * R_phi + tau2 * diag(n)
  
  # Decomposizione di Cholesky per stabilità numerica
  chol_Sigma_y <- chol(Sigma_y)
  
  # Calcolo log determinante
  log_det_Sigma_y <- 2 * sum(log(diag(chol_Sigma_y)))
  
  # Calcolo forma quadratica
  residuo <- y - X %*% beta
  L <- t(chol_Sigma_y)
  z <- forwardsolve(L, residuo)
  quad <- sum(z^2)
  
  log_dens <- -0.5 * log_det_Sigma_y - 0.5 * quad
  
  return(log_dens)
}

#funzione per campionare beta in forma chiusa
campiona_beta_marg <- function(y, X, sigma2, tau2, phi, matrice_dist) {
  n <- length(y)
  p <- ncol(X)
  
  R_phi <- costruisci_R_phi(phi, matrice_dist)
  Sigma_y <- sigma2 * R_phi + tau2 * diag(n)
  
  # Inversa di Sigma_y ccon Cholesky
  chol_Sigma <- chol(Sigma_y)
  Sigma_y_inv <- chol2inv(chol_Sigma)
  
  precisione_beta <- t(X) %*% Sigma_y_inv %*% X + (1/1e4) * diag(p)
  
  # Varianza a posteriori con Cholesky
  chol_prec <- chol(precisione_beta)
  V_beta <- chol2inv(chol_prec)
  
  # Media a posteriori: V_beta * X' Sigma_y^{-1} y
  mu_beta <- V_beta %*% (t(X) %*% Sigma_y_inv %*% y)
  
  # Campionamento da N(mu_beta, V_beta)
  return(as.numeric(rmvnorm(1, mean = as.vector(mu_beta), sigma = V_beta)))
}

#Inizializzazione parametri come prima

beta_marg_0 <- coef(mod_lm)
beta_marg_catena <- matrix(NA, nrow = B, ncol = 3)
beta_marg_catena[1,] <- beta_marg_0


tau2_marg_0 <- var(residui_lm) * 0.5
tau2_marg_catena <- matrix(NA, nrow = B, ncol = 1)
tau2_marg_catena[1,] <- tau2_marg_0


sigma2_marg_0 <- var(residui_lm) * 0.5
sigma2_marg_catena <- matrix(NA, nrow = B, ncol = 1)
sigma2_marg_catena[1,] <- sigma2_marg_0

#log_phi_max e log_phi_min sono gli stessi che avevo definito per il modello
# non marginalizzato
log_phi_marg_0 <- (log_phi_min + log_phi_max) / 2
log_phi_marg_catena <- rep(NA, B)
log_phi_marg_catena[1] <- log_phi_marg_0

# Varianza proposta per log(phi)
varianza_log_marg <- (range_log_phi / 10)^2
# Varianze proposta per tau2 e sigma2 
varianza_prop_tau2 <- 0.1
varianza_prop_sigma2 <- 1
#Campio parametri funzione adattamento per avere tasso accettazione desiderato
# (trovati per mezzo di vari tentavi)
A <- 10     
D <- 50      
c <- 10     

# Catene per i tassi di accettazione
alfa_tau2_marg_catena <- rep(NA, B)
alfa_sigma2_marg_catena <- rep(NA, B)
alfa_phi_marg_catena <- rep(NA, B)

alfa_tau2_marg_catena[1] <- 0
alfa_sigma2_marg_catena[1] <- 0
alfa_phi_marg_catena[1] <- 0


for (i in 2:B) {
  # Valori correnti
  beta_marg_corrente <- beta_marg_catena[i-1,]
  tau2_marg_corrente <- tau2_marg_catena[i-1,]
  sigma2_marg_corrente <- sigma2_marg_catena[i-1,]
  log_phi_marg_corrente <- log_phi_marg_catena[i-1]
  phi_marg_corrente <- exp(log_phi_marg_corrente)
  
  #campionamento beta in forma chiusa
  beta_marg_corrente <- campiona_beta_marg(y, X, sigma2_marg_corrente, 
                                           tau2_marg_corrente, phi_marg_corrente, 
                                           matrice_dist)
  beta_marg_catena[i,] <- beta_marg_corrente
  
  #Metropolis
  #tau
  # Proposta in scala logaritmica per garantire positività
  log_tau2_corrente <- log(tau2_marg_corrente)
  log_tau2_proposta <- rnorm(1, log_tau2_corrente, sqrt(varianza_prop_tau2))
  tau2_proposta <- exp(log_tau2_proposta)
  
  # Calcolo log-densità, cosiderando anche Jacobaino trasformazione densità
  log_post_tau2_proposta <- log_dens_tau2_marg(tau2_proposta, y, X, 
                                               beta_marg_corrente, sigma2_marg_corrente, 
                                               phi_marg_corrente, matrice_dist) + log_tau2_proposta
  
  log_post_tau2_corrente <- log_dens_tau2_marg(tau2_marg_corrente, y, X, 
                                               beta_marg_corrente, sigma2_marg_corrente, 
                                               phi_marg_corrente, matrice_dist) + log_tau2_corrente
  
  # Tasso di accettazione
  log_alpha_tau2 <- log_post_tau2_proposta - log_post_tau2_corrente
  alfa_tau2 <- min(1, exp(log_alpha_tau2))
  alfa_tau2_marg_catena[i] <- alfa_tau2
  
  # Accetta/rifiuta
  if (runif(1) <= alfa_tau2)
    {
    tau2_marg_catena[i,] <- tau2_proposta
    } 
  else 
    {
    tau2_marg_catena[i,] <- tau2_marg_corrente
    }
  tau2_marg_corrente <- tau2_marg_catena[i,]
  
  #sigma2
  # Proposta in scala logaritmica per garantire positività
  log_sigma2_corrente <- log(sigma2_marg_corrente)
  log_sigma2_proposta <- rnorm(1, log_sigma2_corrente, sqrt(varianza_prop_sigma2))
  sigma2_proposta <- exp(log_sigma2_proposta)
  
  # Calcolo log-densità, considerando anche Jacobiano trasformazione densità
  log_post_sigma2_proposta <- log_dens_sigma2_marg(sigma2_proposta, y, X, 
                                                   beta_marg_corrente, tau2_marg_corrente, 
                                                   phi_marg_corrente, matrice_dist) + log_sigma2_proposta
  
  log_post_sigma2_corrente <- log_dens_sigma2_marg(sigma2_marg_corrente, y, X, 
                                                   beta_marg_corrente, tau2_marg_corrente, 
                                                   phi_marg_corrente, matrice_dist) + log_sigma2_corrente
  
  # Tasso di accettazione
  log_alpha_sigma2 <- log_post_sigma2_proposta - log_post_sigma2_corrente
  alfa_sigma2 <- min(1, exp(log_alpha_sigma2))
  alfa_sigma2_marg_catena[i] <- alfa_sigma2
  
  # Accetta/rifiuta
  if (runif(1) <= alfa_sigma2) 
    {
    sigma2_marg_catena[i,] <- sigma2_proposta
    } 
  else
    {
    sigma2_marg_catena[i,] <- sigma2_marg_corrente
    }
  sigma2_marg_corrente <- sigma2_marg_catena[i,]
  
  #phi
  # Proposta in scala logaritmica
  log_phi_proposta_marg <- rnorm(1, log_phi_marg_corrente, sqrt(varianza_log_marg))
  
  # Controlla se log_phi proposto è nel range ammissibile
  if (log_phi_proposta_marg < log_phi_min || log_phi_proposta_marg > log_phi_max)
    {
    alfa_phi_marg_catena[i] <- 0
    log_phi_marg_catena[i] <- log_phi_marg_corrente
    } 
  else
    {
    phi_proposta_marg <- exp(log_phi_proposta_marg)
    
    # Calcola log-densità a posteriori con Jacobiano
    log_post_phi_proposta <- log_dens_phi_marg(phi_proposta_marg, y, X, 
                                               beta_marg_corrente, sigma2_marg_corrente, 
                                               tau2_marg_corrente, matrice_dist) + log_phi_proposta_marg
    
    log_post_phi_corrente <- log_dens_phi_marg(phi_marg_corrente, y, X, 
                                               beta_marg_corrente, sigma2_marg_corrente, 
                                               tau2_marg_corrente, matrice_dist) + log_phi_marg_corrente
    
    # Calcola tasso di accettazione
    log_alpha_phi <- log_post_phi_proposta - log_post_phi_corrente
    alfa_phi <- min(1, exp(log_alpha_phi))
    alfa_phi_marg_catena[i] <- alfa_phi
    
    # Accetta/rifiuta
    if (runif(1) <= alfa_phi) {
      log_phi_marg_catena[i] <- log_phi_proposta_marg
    } else {
      log_phi_marg_catena[i] <- log_phi_marg_corrente
    }
  }

  
  if (i %% c == 0 && i >= iter_inizio_adatt && i <= iter_fine_adatt) 
    {
    gamma_i <- A / (D + i)
    # Adattamento per tau2
    alfa_media_tau2 <- mean(alfa_tau2_marg_catena[(i-c+1):i])
    log_varianza_tau2_new <- log(varianza_prop_tau2) + gamma_i * (alfa_media_tau2 - alfa_desiderato)
    varianza_prop_tau2 <- exp(log_varianza_tau2_new)
    
    # Adattamento per sigma2
    alfa_media_sigma2 <- mean(alfa_sigma2_marg_catena[(i-c+1):i])
    log_varianza_sigma2_new <- log(varianza_prop_sigma2) + gamma_i * (alfa_media_sigma2 - alfa_desiderato)
    varianza_prop_sigma2 <- exp(log_varianza_sigma2_new)
    
    # Adattamento per phi
    alfa_media_phi <- mean(alfa_phi_marg_catena[(i-c+1):i])
    log_varianza_phi_new <- log(varianza_log_marg) + gamma_i * (alfa_media_phi - alfa_desiderato)
    varianza_log_marg <- exp(log_varianza_phi_new)
    #Limiti di sicurezza analoghi al caso precedente
    varianza_log_marg <- min(varianza_log_marg, (range_log_phi / 2)^2)
    varianza_log_marg <- max(varianza_log_marg, (range_log_phi / 100)^2)
    
    # Controllo periodico
    if (i %% 200 == 0) {
      cat(sprintf("Iter %d: alfa_tau2=%.3f, alfa_sigma2=%.3f, alfa_phi=%.3f\n", 
                  i, alfa_media_tau2, alfa_media_sigma2, alfa_media_phi))
    }
  }
  
  # Controllo esecuzione
  if (i %% 2000 == 0 && i > B_burn) {
    cat(sprintf("Modello Marg - Iter %d\n", i))
  }
}

tempo_marg <- Sys.time() - inizio_marg
cat("\nTempo esecuzione modello marginalizzato:", tempo_marg, "\n")

# Trasforma log(phi) in phi
phi_marg_catena <- exp(log_phi_marg_catena)

# Burn-in e thinning insieme: ifatti prendo direttamente gli indici_post
# presi prima dopo burn-in e thinning
beta_marg_post  <- beta_marg_catena[indici_post, ]
tau2_marg_post  <- tau2_marg_catena[indici_post, , drop = FALSE]
sigma2_marg_post <- sigma2_marg_catena[indici_post, , drop = FALSE]
phi_marg_post   <- phi_marg_catena[indici_post]

alfa_tau2_marg_post <- alfa_tau2_marg_catena[indici_post]
alfa_sigma2_marg_post <- alfa_sigma2_marg_catena[indici_post]
alfa_phi_marg_post <- alfa_phi_marg_catena[indici_post]

cat("Campioni modello marginalizzato dopo burn-in e thinning:", nrow(beta_marg_post), "\n\n")


cat("\n=== STATISTICHE POSTERIORI MODELLO MARGINALIZZATO ===\n\n")

# Tabella per beta
tabella_beta_marg <- rbind(
  mu    = riassunto_posteriori(beta_marg_post[, 1]),
  beta1 = riassunto_posteriori(beta_marg_post[, 2]),
  beta2 = riassunto_posteriori(beta_marg_post[, 3])
)

cat("PARAMETRI BETA:\n")
print(round(tabella_beta_marg, 4))
cat("\n")

# Tabella per varianze e phi
tabella_var_phi_marg <- rbind(
  tau2   = riassunto_posteriori(tau2_marg_post),
  sigma2 = riassunto_posteriori(sigma2_marg_post),
  phi    = riassunto_posteriori(phi_marg_post)
)

cat("PARAMETRI DI VARIANZA E PHI:\n")
print(round(tabella_var_phi_marg, 4))
cat("\n")

# Tassi di accettazione
cat("TASSI DI ACCETTAZIONE METROPOLIS-HASTINGS:\n")
cat(sprintf("  tau2:   %.4f (ottimale: ~0.23-0.25)\n", mean(alfa_tau2_marg_post)))
cat(sprintf("  sigma2: %.4f (ottimale: ~0.23-0.25)\n", mean(alfa_sigma2_marg_post)))
cat(sprintf("  phi:    %.4f (ottimale: ~0.23-0.25)\n", mean(alfa_phi_marg_post)))
cat("\n")

#ESERCIZIO 2
#è stato svolto in un file diverso, copia incollo il file, commentanto le prime righe inutili,
# ma lasciandolo integro in caso che l'eseguzione di questo file completo generi 
#problemi

#set.seed(359445) 

#library(ggplot2) #Uso IA per grafici
#library(gridExtra) 
#library(coda)

#librerie per campionare
library(mvtnorm)
library(MCMCpack) # Per Dirichlet e Inverse-Wishart

# Imposto cartella di lavoro per salvataggio dati
setwd("C:/Users/sierv/Desktop/homework 2/es2")

# Carichiamo il file CSV
data_dogs <- read.csv("data_dogs2.csv")

# Calcolo del tempo di permanenza in ogni punto
# Per ogni punto calcolo la distanza dal punto successivo, aggiungo NA per avere 
# stessa lunghezza colonne data.frame perche diff restituisce vettore N-1
data_dogs$distanza <- c(sqrt(diff(data_dogs$X)^2 + diff(data_dogs$Y)^2), NA)


# Definisco una soglia per considerare il cane "fermo" (distanza molto piccola)
#questo grafico che non salvo serve a vedere la distribuzione della distanza per stabilire 
# una soglia ragionevole 
plot_distribuzione <- ggplot(data_dogs, aes(x = distanza)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Distribuzione delle distanze percorse (step-by-step)",
       x = "Distanza (unità di mappa)",
       y = "Frequenza (Conteggio)")

print(plot_distribuzione)

ggsave("distribuzione_distanza.png", plot_distribuzione, width = 10, height = 8)

soglia_fermo <- 0.15


# Creo una variabile che conta quante osservazioni consecutive il cane resta fermo
# Inizializzo il contatore di permanenza
data_dogs$tempo_permanenza <- 0

#per la prima riga inizializzo a mano
data_dogs$tempo_permanenza[1] <- 1

# Ciclo per calcolare il tempo di permanenza
contatore <- 1
for(i in 2:(nrow(data_dogs)-1)) {
  # Se la distanza è piccola (il cane è fermo), incremento il contatore
  # guardo la distanza a i-1 per come funziona diff che calcola la distanza
  # tra i-1 e i, assegnandola alla riga i-1
  if(data_dogs$distanza[i-1] < soglia_fermo) {
    contatore <- contatore + 1
    data_dogs$tempo_permanenza[i] <- contatore
  }
  else 
  {
    # Se il cane si è mosso, resetto il contatore
    contatore <- 1
    data_dogs$tempo_permanenza[i] <- 1
  }
}





# Grafico spaziale con scala di colori basata sul tempo di permanenza
plot_spaziale <- ggplot(data_dogs, aes(x = X, y = Y)) +
  geom_path(alpha = 0.3, color = "gray70") +  # Traiettoria in grigio chiaro
  geom_point(aes(color = tempo_permanenza), size = 2, alpha = 0.7) +  # Punti colorati in base alla permanenza
  scale_color_gradient(low = "yellow", high = "darkblue",  # Scala da chiaro (poco tempo) a scuro (molto tempo)
                       name = "Tempo di\npermanenza") +
  theme_minimal() +
  labs(title = "Traiettoria del cane con tempo di permanenza", 
       x = "Coordinata X", 
       y = "Coordinata Y",
       subtitle = "Colori chiari = movimento, colori scuri = sosta prolungata")

# Visualizzo il grafico
print(plot_spaziale)

# Salvo il grafico
ggsave("traiettoria_cane_permanenza.png", plot_spaziale, width = 15, height = 10)


#creo un indice temporare per il prossimo grafico
data_dogs$time <- 1:nrow(data_dogs)

# Primo grafico: tempo vs coordinata X
grafico_x <- ggplot(data_dogs, aes(x = time, y = X)) +
  geom_line(color = "blue", alpha = 0.7) +
  geom_point(color = "blue", alpha = 0.5, size = 1) +
  theme_minimal() +
  labs(title = "Coordinata X nel tempo", 
       x = "Tempo (indice)", 
       y = "Coordinata X")

# Secondo grafico: tempo vs coordinata Y
grafico_y <- ggplot(data_dogs, aes(x = time, y = Y)) +
  geom_line(color = "red", alpha = 0.7) +
  geom_point(color = "red", alpha = 0.5, size = 1) +
  theme_minimal() +
  labs(title = "Coordinata Y nel tempo", 
       x = "Tempo (indice)", 
       y = "Coordinata Y")

# Combino i due grafici affiancati
grafici_combinati <- grid.arrange(grafico_x, grafico_y, ncol = 2)

# Salvo il grafico combinato
ggsave("coordinate_nel_tempo.png", grafici_combinati, width = 14, height = 6, dpi = 300)

#MCMC

#numero di iterazioni
B <- 50000
B_burn <- 20000

#Inizializzazione
L <- 2
n <- nrow(data_dogs)
catena_m1 <- matrix(NA, nrow = B, ncol = 2)
catena_m2 <- matrix(NA, nrow = B, ncol = 2)
catena_m1[1,] <- c(mean(data_dogs$X), mean(data_dogs$Y))
catena_m2[1,] <- c(mean(data_dogs$X), mean(data_dogs$Y))

catena_sigma1 <- matrix(NA, nrow = B, ncol = 4)
catena_sigma2 <- matrix(NA, nrow = B, ncol = 4)
catena_sigma1[1,] <- c(var(data_dogs$X), cov(data_dogs$X, data_dogs$Y),
                       cov(data_dogs$X, data_dogs$Y),
                       var(data_dogs$Y))
catena_sigma2[1,] <- c(var(data_dogs$X), cov(data_dogs$X, data_dogs$Y),
                       cov(data_dogs$X, data_dogs$Y),
                       var(data_dogs$Y))                       

catena_z0 <- rep(NA, B)
catena_z0[1] <- 1 #scelgo 1 in quanto equiprobabili

catena_pi1 <- matrix(NA, nrow = B, ncol =2)
catena_pi2 <- matrix(NA, nrow = B, ncol =2)
catena_pi1[1,] <- c(0.5,0.5)
catena_pi2[1,] <- c(0.5,0.5)

catena_z <- matrix(NA,nrow = B, ncol = n)
catena_z[1,] <- sample(c(1,2), size = n, replace = TRUE)
#campiona tra 1 e 2 con probabilità entrambi 0.5


for (i in 2:B)
{
  #calcolo gli n_k e la somma degli y_k associati a quei cluster
  z_corrente <- catena_z[i-1,]
  indici_1 <- which(z_corrente == 1)
  indici_2 <- which(z_corrente == 2)
  
  n1 <- length(indici_1)
  n2 <- length(indici_2)
  
  if(n1 > 0)
  {
    somma_y1 <- colSums(data_dogs[indici_1, c("X", "Y")])
  }
  else
  {
    somma_y1 <- c(0,0)
  }
  
  if(n2 > 0)
  {
    somma_y2 <- colSums(data_dogs[indici_2, c("X", "Y")])
  }
  else
  {
    somma_y2 <- c(0,0)
  }
  
  #calcolo parametri full conditional di mu_1
  sigma_1 <- matrix(catena_sigma1[i-1,], nrow = 2, ncol = 2)
  #inverto con cholesky per capacità computazionali
  sigma_1inv <- chol2inv(chol(sigma_1))
  V_inv <- n1*sigma_1inv + 0.01
  V <- chol2inv(chol(V_inv))
  media <- V %*% (sigma_1inv %*% somma_y1)
  catena_m1[i,] <- rmvnorm(1, mean = media, sigma = V)
  
  #calcolo parametri full conditional di mu_2
  sigma_2 <- matrix(catena_sigma2[i-1,], nrow = 2, ncol = 2)
  #inverto con cholesky per capacità computazionali
  sigma_2inv <- chol2inv(chol(sigma_2))
  V_inv <- n2*sigma_2inv + diag(2)/100
  V <- chol2inv(chol(V_inv))
  media <- V %*% (sigma_2inv %*% somma_y2)
  catena_m2[i,] <- rmvnorm(1, mean = media, sigma = V)
  
  #campiono sigma_1
  v_post1 <- n1 + 4
  if (n1 > 0)
  {
    #filtro le coordinate dei punti associati al claster z = 1
    coordinata_1 <- as.matrix(data_dogs[indici_1, c("X","Y")])
    #uso sweep per sottrare mu_1 vettore a ogni riga della matrice 
    matrice_appoggio <- sweep(coordinata_1,2,catena_m1[i,], "-")
    sommatoria_finale1 <- t(matrice_appoggio) %*% matrice_appoggio
  }
  else
  {
    sommatoria_finale1 <- matrix(0, nrow = 2, ncol = 2)
  }
  
  phi_post1 <- diag(2) + sommatoria_finale1
  #campiono
  catena_sigma1[i,] <- as.vector(riwish(v_post1, phi_post1))
  
  #campiono sigma_2
  v_post2 <- n2 + 4
  if (n2 > 0)
  {
    #filtro le coordinate dei punti associati al claster z = 1
    coordinata_2 <- as.matrix(data_dogs[indici_2, c("X","Y")])
    #uso sweep per sottrare mu_1 vettore a ogni riga della matrice 
    matrice_appoggio <- sweep(coordinata_2,2,catena_m2[i,], "-")
    sommatoria_finale2 <- t(matrice_appoggio) %*% matrice_appoggio
  }
  else
  {
    sommatoria_finale2 <- matrix(0, nrow = 2, ncol = 2)
  }
  phi_post2 <- diag(2) + sommatoria_finale2
  #campiono, convertendo a vettore l output che senno sarebbe matrice
  catena_sigma2[i,] <- as.vector(riwish(v_post2, phi_post2))
  
  #conto i passaggi da uno stato all altro 
  sequenza_stati <- catena_z[i-1,]
  #il vettore da contiente tutti gli stati tranne l'ultimo
  da <- sequenza_stati[-length(sequenza_stati)]
  #il vettore a contiene tutti gli stati tranne il primo
  a <- sequenza_stati[-1]
  # table() accoppia gli elementi di 'da' e 'a' e conta le occorrenze
  #factor con levels=1:2 serve a forzare table a creare una matrice 2*2 anche se 
  # alcune transizioni hanno conteggio 0
  matrice_conteggi <- table(Da = factor(da, levels = 1:2), A = factor(a, levels = 1:2))
  
  #campiono pi_1
  catena_pi1[i,] <- rdirichlet(1, alpha = c(1 + matrice_conteggi[1,1], 
                                            1 + matrice_conteggi[1,2]))
  
  #campiono pi_2
  catena_pi2[i,] <- rdirichlet(1, alpha = c(1 + matrice_conteggi[2,1], 
                                            1 + matrice_conteggi[2,2]))
  
  
  #creo la matrice di transizione per facilitare calcolo alfa_i(k) e z0 campionamento
  pi1_corr <- catena_pi1[i,]
  pi2_corr <- catena_pi2[i,]
  P <- matrix(c(pi1_corr,pi2_corr), nrow = 2, ncol = 2, byrow = TRUE)
  
  #campiono z0
  # prendo lo stato z1 dell'iterazione precedente
  z1_old <- catena_z[i-1,1]
  #estraggo le prob di transizione verso z1
  prob_z0 <- P[, z1_old]
  #campiono z0
  catena_z0[i] <- sample(1:2, 1, prob = prob_z0)
  
  
  #calcolo degli alfa_i(k) del forward-backward
  alfa_matrice <- matrix(NA, nrow = n, ncol = 2)
  #pre calcolo le densita in una matrice dati i parametri legati a 1 e 2,
  # ovvero calcolo subito le densità f(y_t | z_t=k) per tutti i t e k, t indice
  # riga, k colonna
  mu1_corr <- catena_m1[i,]
  mu2_corr <- catena_m2[i,]
  sigma1_corr <- matrix(catena_sigma1[i,], nrow = 2, ncol = 2)
  sigma2_corr <- matrix(catena_sigma2[i,], nrow = 2, ncol = 2)
  dens_calcolate <- matrix(NA, nrow = n, ncol = 2)
  dens_calcolate[,1] <- dmvnorm(data_dogs[,c("X","Y")],
                                mean = mu1_corr, sigma = sigma1_corr)
  dens_calcolate[,2] <- dmvnorm(data_dogs[,c("X","Y")],
                                mean = mu2_corr, sigma = sigma2_corr)
  
  z0_corr <- catena_z0[i]
  
  # alfa_1(k) = P(y_1 | z_1=k) * P(z_1=k | z_0)
  #faccio prodotto elemento per elemento in modo da avere vettore 2*1
  alfa_matrice[1,] <- dens_calcolate[1,] * P[z0_corr,]
  
  #Normalizzo le probabilità per evitare che diventino numeri troppo piccoli
  alfa_matrice[1,] <- alfa_matrice[1,] / sum(alfa_matrice[1,])
  
  # t è l indice temporale/di riga (negli appunti i)
  for (t in 2:n)
  {
    for (k in 1:2) 
    {
      #calcolo sommatoria_ki1 ai1(ki1) * P_ki1_k per un k fissato, ovvero
      # alfa_t(k) propostozionale a: Likelihood_t * Somma(alfa_t-1 * Transizione)
      sommatoria <- 0
      for (k_i1 in 1:2)
      {
        sommatoria <- sommatoria + P[k_i1,k] * alfa_matrice[t-1, k_i1]
      }
      alfa_matrice[t,k] <- dens_calcolate[t,k] * sommatoria
      
    }
    
    #Normalizzo
    alfa_matrice[t,] <- alfa_matrice[t,] / sum(alfa_matrice[t,])
    
  }
  
  #procedo al campionamento
  z_new <- rep(NA, n)
  #zn, campiono da una distreta che associa a 1 prob alfa_n(1) e a 2 alfa_n(2)
  z_new[n] <- sample(1:2,1, prob = alfa_matrice[n,])
  # f(zt |...) proporsionale alfa_t(zt)*P(zt,zt+1)
  
  for (t in (n-1):1) 
  {
    #guardo lo z che ho appena campionato
    z_futuro <- z_new[t+1]
    # estraggo la COLONNA della matrice P relativa allo stato futuro.
    # questa colonna contiene: [P(1->z_future), P(2->z_future)]
    prob_transizione_verso_futuro <- P[, z_futuro]
    #calcolo il peso per z_t: Alfa * Transizione
    #è una moltiplicazione elemento per elemento (1 con 1, 2 con 2)
    prob_t <- alfa_matrice[t, ] * prob_transizione_verso_futuro
    #controllo underflow
    if(sum(prob_t) == 0)
    {
      print("underflow, prob nulle")
    }
    #campiono
    z_new[t] <- sample(1:2,1, prob = prob_t)
  }
  catena_z[i,] <- z_new
  
  #print ogni 2000 iterazioni di controllo
  if(i %% 2000 == 0) 
  {
    print(paste("Iterazione", i, "di", B))
  }
  
  
}

# Rimuovo il burn-in
catena_m1_final <- catena_m1[(B_burn+1):B,]
catena_m2_final <- catena_m2[(B_burn+1):B,]
catena_sigma1_final <- catena_sigma1[(B_burn+1):B,]
catena_sigma2_final <- catena_sigma2[(B_burn+1):B,]
catena_pi1_final <- catena_pi1[(B_burn+1):B,]
catena_pi2_final <- catena_pi2[(B_burn+1):B,]
catena_z_final <- catena_z[(B_burn+1):B,]


#Plot delle catene e anilisi dei risultati
#trace plots 
trace_mu1_x <- ggplot(data.frame(iter = 1:nrow(catena_m1_final), 
                                 value = catena_m1_final[,1]), 
                      aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "steelblue") +
  theme_minimal() +
  labs(title = "Trace plot μ₁ (coordinata X)", x = "Iterazione", y = "μ₁ₓ")

trace_mu1_y <- ggplot(data.frame(iter = 1:nrow(catena_m1_final), 
                                 value = catena_m1_final[,2]), 
                      aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "coral") +
  theme_minimal() +
  labs(title = "Trace plot μ₁ (coordinata Y)", x = "Iterazione", y = "μ₁ᵧ")

trace_mu2_x <- ggplot(data.frame(iter = 1:nrow(catena_m2_final), 
                                 value = catena_m2_final[,1]), 
                      aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "steelblue") +
  theme_minimal() +
  labs(title = "Trace plot μ₂ (coordinata X)", x = "Iterazione", y = "μ₂ₓ")

trace_mu2_y <- ggplot(data.frame(iter = 1:nrow(catena_m2_final), 
                                 value = catena_m2_final[,2]), 
                      aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "coral") +
  theme_minimal() +
  labs(title = "Trace plot μ₂ (coordinata Y)", x = "Iterazione", y = "μ₂ᵧ")

# Salvataggio
grid_mu <- grid.arrange(trace_mu1_x, trace_mu1_y, trace_mu2_x, trace_mu2_y, 
                        ncol = 2)
ggsave("trace_plots_mu_L2.png", grid_mu, width = 14, height = 10, dpi = 300)

trace_sigma1_11 <- ggplot(data.frame(iter = 1:nrow(catena_sigma1_final), 
                                     value = catena_sigma1_final[,1]), 
                          aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkgreen") +
  theme_minimal() +
  labs(title = "Trace plot Σ₁[1,1] (var X)", x = "Iterazione", y = "Σ₁[1,1]")

trace_sigma1_22 <- ggplot(data.frame(iter = 1:nrow(catena_sigma1_final), 
                                     value = catena_sigma1_final[,4]), 
                          aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkgreen") +
  theme_minimal() +
  labs(title = "Trace plot Σ₁[2,2] (var Y)", x = "Iterazione", y = "Σ₁[2,2]")

trace_sigma1_12 <- ggplot(data.frame(iter = 1:nrow(catena_sigma1_final), 
                                     value = catena_sigma1_final[,2]), 
                          aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "purple") +
  theme_minimal() +
  labs(title = "Trace plot Σ₁[1,2] (cov)", x = "Iterazione", y = "Σ₁[1,2]")

trace_sigma2_11 <- ggplot(data.frame(iter = 1:nrow(catena_sigma2_final), 
                                     value = catena_sigma2_final[,1]), 
                          aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkgreen") +
  theme_minimal() +
  labs(title = "Trace plot Σ₂[1,1] (var X)", x = "Iterazione", y = "Σ₂[1,1]")

trace_sigma2_22 <- ggplot(data.frame(iter = 1:nrow(catena_sigma2_final), 
                                     value = catena_sigma2_final[,4]), 
                          aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkgreen") +
  theme_minimal() +
  labs(title = "Trace plot Σ₂[2,2] (var Y)", x = "Iterazione", y = "Σ₂[2,2]")

trace_sigma2_12 <- ggplot(data.frame(iter = 1:nrow(catena_sigma2_final), 
                                     value = catena_sigma2_final[,2]), 
                          aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "purple") +
  theme_minimal() +
  labs(title = "Trace plot Σ₂[1,2] (cov)", x = "Iterazione", y = "Σ₂[1,2]")

grid_sigma <- grid.arrange(trace_sigma1_11, trace_sigma1_22, trace_sigma1_12,
                           trace_sigma2_11, trace_sigma2_22, trace_sigma2_12,
                           ncol = 3)
ggsave("trace_plots_sigma_L2.png", grid_sigma, width = 18, height = 8, dpi = 300)

trace_pi1_1 <- ggplot(data.frame(iter = 1:nrow(catena_pi1_final), 
                                 value = catena_pi1_final[,1]), 
                      aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkorange") +
  theme_minimal() +
  labs(title = "Trace plot π₁[1] (prob 1→1)", x = "Iterazione", y = "π₁[1]")

trace_pi1_2 <- ggplot(data.frame(iter = 1:nrow(catena_pi1_final), 
                                 value = catena_pi1_final[,2]), 
                      aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkorange") +
  theme_minimal() +
  labs(title = "Trace plot π₁[2] (prob 1→2)", x = "Iterazione", y = "π₁[2]")

trace_pi2_1 <- ggplot(data.frame(iter = 1:nrow(catena_pi2_final), 
                                 value = catena_pi2_final[,1]), 
                      aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkred") +
  theme_minimal() +
  labs(title = "Trace plot π₂[1] (prob 2→1)", x = "Iterazione", y = "π₂[1]")

trace_pi2_2 <- ggplot(data.frame(iter = 1:nrow(catena_pi2_final), 
                                 value = catena_pi2_final[,2]), 
                      aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkred") +
  theme_minimal() +
  labs(title = "Trace plot π₂[2] (prob 2→2)", x = "Iterazione", y = "π₂[2]")

grid_pi <- grid.arrange(trace_pi1_1, trace_pi1_2, trace_pi2_1, trace_pi2_2, 
                        ncol = 2)
ggsave("trace_plots_pi_L2.png", grid_pi, width = 14, height = 10, dpi = 300)


# Converto in oggetti mcmc per il pacchetto coda
mcmc_mu1_x <- mcmc(catena_m1_final[,1])
mcmc_mu1_y <- mcmc(catena_m1_final[,2])
mcmc_mu2_x <- mcmc(catena_m2_final[,1])
mcmc_mu2_y <- mcmc(catena_m2_final[,2])
mcmc_sigma1_11 <- mcmc(catena_sigma1_final[,1]) 
mcmc_sigma1_12 <- mcmc(catena_sigma1_final[,2]) 
mcmc_sigma1_22 <- mcmc(catena_sigma1_final[,4])
mcmc_sigma2_11 <- mcmc(catena_sigma2_final[,1]) 
mcmc_sigma2_12 <- mcmc(catena_sigma2_final[,2]) 
mcmc_sigma2_22 <- mcmc(catena_sigma2_final[,4])
mcmc_pi1_1 <- mcmc(catena_pi1_final[,1])
mcmc_pi1_2 <- mcmc(catena_pi1_final[,2])
mcmc_pi2_1 <- mcmc(catena_pi2_final[,1])
mcmc_pi2_2 <- mcmc(catena_pi2_final[,2])

# Calcolo ESS (Effective Sample Size)
ess_mu1_x <- effectiveSize(mcmc_mu1_x)
ess_mu1_y <- effectiveSize(mcmc_mu1_y)
ess_mu2_x <- effectiveSize(mcmc_mu2_x)
ess_mu2_y <- effectiveSize(mcmc_mu2_y)
ess_sigma1_11 <- effectiveSize(mcmc_sigma1_11)
ess_sigma1_12 <- effectiveSize(mcmc_sigma1_12)
ess_sigma1_22 <- effectiveSize(mcmc_sigma1_22)
ess_sigma2_11 <- effectiveSize(mcmc_sigma2_11)
ess_sigma2_12 <- effectiveSize(mcmc_sigma2_12)
ess_sigma2_22 <- effectiveSize(mcmc_sigma2_22)
ess_pi1_1 <- effectiveSize(mcmc_pi1_1)
ess_pi1_2 <- effectiveSize(mcmc_pi1_2)
ess_pi2_1 <- effectiveSize(mcmc_pi2_1)
ess_pi2_2 <- effectiveSize(mcmc_pi2_2)


cat("\n=== DIAGNOSTICHE MCMC (L=2) ===\n")
cat("\nEffective Sample Size (ESS):\n")
cat(sprintf("  μ₁ₓ: %.0f\n", ess_mu1_x))
cat(sprintf("  μ₁ᵧ: %.0f\n", ess_mu1_y))
cat(sprintf("  μ₂ₓ: %.0f\n", ess_mu2_x))
cat(sprintf("  μ₂ᵧ: %.0f\n", ess_mu2_y))
cat("\n--- Sigma 1 (Cluster 1) ---\n")
cat(sprintf("  Var(X): %.0f\n", ess_sigma1_11))
cat(sprintf("  Cov(X,Y): %.0f\n", ess_sigma1_12))
cat(sprintf("  Var(Y): %.0f\n", ess_sigma1_22))

cat("\n--- Sigma 2 (Cluster 2) ---\n")
cat(sprintf("  Var(X): %.0f\n", ess_sigma2_11))
cat(sprintf("  Cov(X,Y): %.0f\n", ess_sigma2_12))
cat(sprintf("  Var(Y): %.0f\n", ess_sigma2_22))

cat("\n--- Matrice di Transizione (Pi) ---\n")
cat(sprintf("  pi_1[1->1]: %.0f\n", ess_pi1_1))
cat(sprintf("  pi_1[1->2]: %.0f\n", ess_pi1_2))
cat(sprintf("  pi_2[2->1]: %.0f\n", ess_pi2_1))
cat(sprintf("  pi_2[2->2]: %.0f\n", ess_pi2_2))

# Statistiche riassuntive posteriori
cat("\n=== STATISTICHE POSTERIORI ===\n")
cat("\nμ₁ (centro cluster 1):\n")
cat(sprintf("  X: %.3f [%.3f, %.3f]\n", 
            mean(catena_m1_final[,1]),
            quantile(catena_m1_final[,1], 0.025),
            quantile(catena_m1_final[,1], 0.975)))
cat(sprintf("  Y: %.3f [%.3f, %.3f]\n", 
            mean(catena_m1_final[,2]),
            quantile(catena_m1_final[,2], 0.025),
            quantile(catena_m1_final[,2], 0.975)))

cat("\nμ₂ (centro cluster 2):\n")
cat(sprintf("  X: %.3f [%.3f, %.3f]\n", 
            mean(catena_m2_final[,1]),
            quantile(catena_m2_final[,1], 0.025),
            quantile(catena_m2_final[,1], 0.975)))
cat(sprintf("  Y: %.3f [%.3f, %.3f]\n", 
            mean(catena_m2_final[,2]),
            quantile(catena_m2_final[,2], 0.025),
            quantile(catena_m2_final[,2], 0.975)))

cat("\nProbabilità di transizione π₁ (da stato 1):\n")
cat(sprintf("  1→1: %.3f [%.3f, %.3f]\n", 
            mean(catena_pi1_final[,1]),
            quantile(catena_pi1_final[,1], 0.025),
            quantile(catena_pi1_final[,1], 0.975)))
cat(sprintf("  1→2: %.3f [%.3f, %.3f]\n", 
            mean(catena_pi1_final[,2]),
            quantile(catena_pi1_final[,2], 0.025),
            quantile(catena_pi1_final[,2], 0.975)))

cat("\nProbabilità di transizione π₂ (da stato 2):\n")
cat(sprintf("  2→1: %.3f [%.3f, %.3f]\n", 
            mean(catena_pi2_final[,1]),
            quantile(catena_pi2_final[,1], 0.025),
            quantile(catena_pi2_final[,1], 0.975)))
cat(sprintf("  2→2: %.3f [%.3f, %.3f]\n", 
            mean(catena_pi2_final[,2]),
            quantile(catena_pi2_final[,2], 0.025),
            quantile(catena_pi2_final[,2], 0.975)))

cat("\nΣ₁ (matrice di covarianza cluster 1):\n")
cat(sprintf("  Σ₁[1,1] (var X): %.4f [%.4f, %.4f]\n", 
            mean(catena_sigma1_final[,1]),
            quantile(catena_sigma1_final[,1], 0.025),
            quantile(catena_sigma1_final[,1], 0.975)))
cat(sprintf("  Σ₁[2,2] (var Y): %.4f [%.4f, %.4f]\n", 
            mean(catena_sigma1_final[,4]),
            quantile(catena_sigma1_final[,4], 0.025),
            quantile(catena_sigma1_final[,4], 0.975)))
cat(sprintf("  Σ₁[1,2] (cov XY): %.4f [%.4f, %.4f]\n", 
            mean(catena_sigma1_final[,2]),
            quantile(catena_sigma1_final[,2], 0.025),
            quantile(catena_sigma1_final[,2], 0.975)))

cat("\nΣ₂ (matrice di covarianza cluster 2):\n")
cat(sprintf("  Σ₂[1,1] (var X): %.4f [%.4f, %.4f]\n", 
            mean(catena_sigma2_final[,1]),
            quantile(catena_sigma2_final[,1], 0.025),
            quantile(catena_sigma2_final[,1], 0.975)))
cat(sprintf("  Σ₂[2,2] (var Y): %.4f [%.4f, %.4f]\n", 
            mean(catena_sigma2_final[,4]),
            quantile(catena_sigma2_final[,4], 0.025),
            quantile(catena_sigma2_final[,4], 0.975)))
cat(sprintf("  Σ₂[1,2] (cov XY): %.4f [%.4f, %.4f]\n", 
            mean(catena_sigma2_final[,2]),
            quantile(catena_sigma2_final[,2], 0.025),
            quantile(catena_sigma2_final[,2], 0.975)))


#Punto 4
#Per ogni colonna della matrice delle iterazioni di z, calcolo il 
# valore più frequente e salvo il valore che compare più spesso in
# z moda
z_moda <- rep(NA, n)

for (i in 1:n)
{
  conteggio1 <- sum(catena_z_final[,i] == 1)
  conteggio2 <- sum(catena_z_final[,i] == 2)
  
  #in caso di pareggio attribuisce il valore di moda a 1
  if(conteggio1 >= conteggio2)
  {
    z_moda[i] <- 1
  }
  else
  {
    z_moda[i] <- 2
  }
}

# Aggiungo z_moda al dataframe, salvandolo come fattore per facilitare
# colorazione
data_dogs$z_moda <- as.factor(z_moda)

# Grafico spaziale con colorazione per cluster stimato
plot_spaziale_cluster <- ggplot(data_dogs, aes(x = X, y = Y)) +
  geom_path(alpha = 0.3, color = "gray70") +  # Traiettoria in grigio
  geom_point(aes(color = z_moda), size = 2.5, alpha = 0.8) +  # Punti colorati per cluster
  scale_color_manual(values = c("1" = "steelblue", "2" = "coral"),
                     name = "Cluster stimato",
                     labels = c("Stato 1", "Stato 2")) +
  theme_minimal() +
  labs(title = "Traiettoria del cane con clustering HMM (L=2)", 
       x = "Coordinata X", 
       y = "Coordinata Y",
       subtitle = "Colorazione basata sulla moda a posteriori di z_i")

print(plot_spaziale_cluster)
ggsave("traiettoria_cane_cluster_L2.png", plot_spaziale_cluster, width = 15, height = 10, dpi = 300)


# Serie temporali con colorazione per cluster
grafico_x_cluster <- ggplot(data_dogs, aes(x = time, y = X)) +
  geom_line(color = "gray70", alpha = 0.5) +
  geom_point(aes(color = z_moda), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("1" = "steelblue", "2" = "coral"),
                     name = "Cluster",
                     labels = c("Stato 1", "Stato 2")) +
  theme_minimal() +
  labs(title = "Coordinata X nel tempo - Clustering HMM", 
       x = "Tempo (indice)", 
       y = "Coordinata X")

grafico_y_cluster <- ggplot(data_dogs, aes(x = time, y = Y)) +
  geom_line(color = "gray70", alpha = 0.5) +
  geom_point(aes(color = z_moda), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("1" = "steelblue", "2" = "coral"),
                     name = "Cluster",
                     labels = c("Stato 1", "Stato 2")) +
  theme_minimal() +
  labs(title = "Coordinata Y nel tempo - Clustering HMM", 
       x = "Tempo (indice)", 
       y = "Coordinata Y")

# Combino i due grafici
grafici_temporali_cluster <- grid.arrange(grafico_x_cluster, grafico_y_cluster, ncol = 2)

ggsave("coordinate_tempo_cluster_L2.png", grafici_temporali_cluster, width = 14, height = 6, dpi = 300)







#punto 5, L = 3
#numero di iterazioni
B <- 50000
B_burn <- 20000

#Inizializzazione
L <- 3
n <- nrow(data_dogs)
catena_m1 <- matrix(NA, nrow = B, ncol = 2)
catena_m2 <- matrix(NA, nrow = B, ncol = 2)
catena_m3 <- matrix(NA, nrow = B, ncol = 2)
catena_m1[1,] <- c(mean(data_dogs$X), mean(data_dogs$Y))
catena_m2[1,] <- c(mean(data_dogs$X), mean(data_dogs$Y))
catena_m3[1,] <- c(mean(data_dogs$X), mean(data_dogs$Y))

catena_sigma1 <- matrix(NA, nrow = B, ncol = 4)
catena_sigma2 <- matrix(NA, nrow = B, ncol = 4)
catena_sigma3 <- matrix(NA, nrow = B, ncol = 4)
catena_sigma1[1,] <- c(var(data_dogs$X), cov(data_dogs$X, data_dogs$Y),
                       cov(data_dogs$X, data_dogs$Y),
                       var(data_dogs$Y))
catena_sigma2[1,] <- c(var(data_dogs$X), cov(data_dogs$X, data_dogs$Y),
                       cov(data_dogs$X, data_dogs$Y),
                       var(data_dogs$Y))   
catena_sigma3[1,] <- c(var(data_dogs$X), cov(data_dogs$X, data_dogs$Y),
                       cov(data_dogs$X, data_dogs$Y),
                       var(data_dogs$Y)) 

catena_z0 <- rep(NA, B)
catena_z0[1] <- 1 #scelgo 1 in quanto equiprobabili

catena_pi1 <- matrix(NA, nrow = B, ncol =3)
catena_pi2 <- matrix(NA, nrow = B, ncol =3)
catena_pi3 <- matrix(NA, nrow = B, ncol =3)
catena_pi1[1,] <- c(0.33,0.33,0.34)
catena_pi2[1,] <- c(0.33,0.33,0.34)
catena_pi3[1,] <- c(0.33,0.33,0.34)


catena_z <- matrix(NA,nrow = B, ncol = n)
catena_z[1,] <- sample(c(1,2,3), size = n, replace = TRUE)
#campiona tra 1 e 2 con probabilità entrambi 0.5


for (i in 2:B)
{
  #calcolo gli n_k e la somma degli y_k associati a quei cluster
  z_corrente <- catena_z[i-1,]
  indici_1 <- which(z_corrente == 1)
  indici_2 <- which(z_corrente == 2)
  indici_3 <- which(z_corrente == 3)
  
  n1 <- length(indici_1)
  n2 <- length(indici_2)
  n3 <- length(indici_3)
  
  if(n1 > 0)
  {
    somma_y1 <- colSums(data_dogs[indici_1, c("X", "Y")])
  }
  else
  {
    somma_y1 <- c(0,0)
  }
  
  if(n2 > 0)
  {
    somma_y2 <- colSums(data_dogs[indici_2, c("X", "Y")])
  }
  else
  {
    somma_y2 <- c(0,0)
  }
  
  if(n3 > 0)
  {
    somma_y3 <- colSums(data_dogs[indici_3, c("X", "Y")])
  }
  else
  {
    somma_y3 <- c(0,0)
  }
  
  #calcolo parametri full conditional di mu_1
  sigma_1 <- matrix(catena_sigma1[i-1,], nrow = 2, ncol = 2)
  #inverto con cholesky per capacità computazionali
  sigma_1inv <- chol2inv(chol(sigma_1))
  V_inv <- n1*sigma_1inv + 0.01
  V <- chol2inv(chol(V_inv))
  media <- V %*% (sigma_1inv %*% somma_y1)
  catena_m1[i,] <- rmvnorm(1, mean = media, sigma = V)
  
  #calcolo parametri full conditional di mu_2
  sigma_2 <- matrix(catena_sigma2[i-1,], nrow = 2, ncol = 2)
  #inverto con cholesky per capacità computazionali
  sigma_2inv <- chol2inv(chol(sigma_2))
  V_inv <- n2*sigma_2inv + diag(2)/100
  V <- chol2inv(chol(V_inv))
  media <- V %*% (sigma_2inv %*% somma_y2)
  catena_m2[i,] <- rmvnorm(1, mean = media, sigma = V)
  
  #calcolo parametri full conditional di mu_3
  sigma_3 <- matrix(catena_sigma3[i-1,], nrow = 2, ncol = 2)
  #inverto con cholesky per capacità computazionali
  sigma_3inv <- chol2inv(chol(sigma_3))
  V_inv <- n3*sigma_3inv + diag(2)/100
  V <- chol2inv(chol(V_inv))
  media <- V %*% (sigma_3inv %*% somma_y3)
  catena_m3[i,] <- rmvnorm(1, mean = media, sigma = V)
  
  #campiono sigma_1
  v_post1 <- n1 + 4
  if (n1 > 0)
  {
    #filtro le coordinate dei punti associati al claster z = 1
    coordinata_1 <- as.matrix(data_dogs[indici_1, c("X","Y")])
    #uso sweep per sottrare mu_1 vettore a ogni riga della matrice 
    matrice_appoggio <- sweep(coordinata_1,2,catena_m1[i,], "-")
    sommatoria_finale1 <- t(matrice_appoggio) %*% matrice_appoggio
  }
  else
  {
    sommatoria_finale1 <- matrix(0, nrow = 2, ncol = 2)
  }
  
  phi_post1 <- diag(2) + sommatoria_finale1
  #campiono
  catena_sigma1[i,] <- as.vector(riwish(v_post1, phi_post1))
  
  #campiono sigma_2
  v_post2 <- n2 + 4
  if (n2 > 0)
  {
    #filtro le coordinate dei punti associati al claster z = 1
    coordinata_2 <- as.matrix(data_dogs[indici_2, c("X","Y")])
    #uso sweep per sottrare mu_1 vettore a ogni riga della matrice 
    matrice_appoggio <- sweep(coordinata_2,2,catena_m2[i,], "-")
    sommatoria_finale2 <- t(matrice_appoggio) %*% matrice_appoggio
  }
  else
  {
    sommatoria_finale2 <- matrix(0, nrow = 2, ncol = 2)
  }
  phi_post2 <- diag(2) + sommatoria_finale2
  #campiono, convertendo a vettore l output che senno sarebbe matrice
  catena_sigma2[i,] <- as.vector(riwish(v_post2, phi_post2))
  
  #campiono sigma_3
  v_post3 <- n3 + 4
  if (n3 > 0)
  {
    #filtro le coordinate dei punti associati al claster z = 1
    coordinata_3 <- as.matrix(data_dogs[indici_3, c("X","Y")])
    #uso sweep per sottrare mu_1 vettore a ogni riga della matrice 
    matrice_appoggio <- sweep(coordinata_3,2,catena_m3[i,], "-")
    sommatoria_finale3 <- t(matrice_appoggio) %*% matrice_appoggio
  }
  else
  {
    sommatoria_finale3 <- matrix(0, nrow = 2, ncol = 2)
  }
  
  phi_post3 <- diag(2) + sommatoria_finale3
  #campiono
  catena_sigma3[i,] <- as.vector(riwish(v_post1, phi_post3))
  
  #conto i passaggi da uno stato all altro 
  sequenza_stati <- catena_z[i-1,]
  #il vettore da contiente tutti gli stati tranne l'ultimo
  da <- sequenza_stati[-length(sequenza_stati)]
  #il vettore a contiene tutti gli stati tranne il primo
  a <- sequenza_stati[-1]
  # table() accoppia gli elementi di 'da' e 'a' e conta le occorrenze
  #factor con levels=1:3 serve a forzare table a creare una matrice 3*3 anche se 
  # alcune transizioni hanno conteggio 0
  matrice_conteggi <- table(Da = factor(da, levels = 1:3), A = factor(a, levels = 1:3))
  
  #campiono pi_1
  catena_pi1[i,] <- rdirichlet(1, alpha = c(1 + matrice_conteggi[1,1], 
                                            1 + matrice_conteggi[1,2],
                                            1 + matrice_conteggi[1,3]))
  
  #campiono pi_2
  catena_pi2[i,] <- rdirichlet(1, alpha = c(1 + matrice_conteggi[2,1], 
                                            1 + matrice_conteggi[2,2],
                                            1 + matrice_conteggi[2,3]))
  
  #campiono pi_3
  catena_pi3[i,] <- rdirichlet(1, alpha = c(1 + matrice_conteggi[3,1], 
                                            1 + matrice_conteggi[3,2],
                                            1 + matrice_conteggi[3,3]))
  
  
  
  #creo la matrice di transizione per facilitare calcolo alfa_i(k) e z0 campionamento
  pi1_corr <- catena_pi1[i,]
  pi2_corr <- catena_pi2[i,]
  pi3_corr <- catena_pi3[i,]
  P <- matrix(c(pi1_corr,pi2_corr,pi3_corr), nrow = 3, ncol = 3, byrow = TRUE)
  
  #campiono z0
  # prendo lo stato z1 dell'iterazione precedente
  z1_old <- catena_z[i-1,1]
  #estraggo le prob di transizione verso z1
  prob_z0 <- P[, z1_old]
  #campiono z0
  catena_z0[i] <- sample(1:3, 1, prob = prob_z0)
  
  
  #calcolo degli alfa_i(k) del forward-backward
  alfa_matrice <- matrix(NA, nrow = n, ncol = 3)
  #pre calcolo le densita in una matrice dati i parametri legati a 1,2,3
  # Calcolo f(y_t | z_t = k) per ogni t e ogni k=1,2,3
  mu1_corr <- catena_m1[i,]
  mu2_corr <- catena_m2[i,]
  mu3_corr <- catena_m3[i,]
  sigma1_corr <- matrix(catena_sigma1[i,], nrow = 2, ncol = 2)
  sigma2_corr <- matrix(catena_sigma2[i,], nrow = 2, ncol = 2)
  sigma3_corr <- matrix(catena_sigma3[i,], nrow = 2, ncol = 2)
  dens_calcolate <- matrix(NA, nrow = n, ncol = 3)
  dens_calcolate[,1] <- dmvnorm(data_dogs[,c("X","Y")],
                                mean = mu1_corr, sigma = sigma1_corr)
  dens_calcolate[,2] <- dmvnorm(data_dogs[,c("X","Y")],
                                mean = mu2_corr, sigma = sigma2_corr)
  dens_calcolate[,3] <- dmvnorm(data_dogs[,c("X","Y")],
                                mean = mu3_corr, sigma = sigma3_corr)
  
  z0_corr <- catena_z0[i]
  
  # Passo base (t=1): P(y_1|z_1) * P(z_1|z_0)
  #faccio prodotto elemento per elemento in modo da avere vettore 3*1
  alfa_matrice[1,] <- dens_calcolate[1,] * P[z0_corr,]
  
  #Normalizzo le probabilità per evitare che diventino numeri troppo piccoli
  alfa_matrice[1,] <- alfa_matrice[1,] / sum(alfa_matrice[1,])
  
  # t è l indice temporale/di riga (negli appunti i)
  # alpha_t(k) = f(y_t|k) * sum_j(alpha_{t-1}(j) * P(j->k))
  for (t in 2:n)
  {
    for (k in 1:3) 
    {
      #calcolo sommatoria_ki1 ai1(ki1) * P_ki1_k per un k fissato
      sommatoria <- 0
      for (k_i1 in 1:3)
      {
        sommatoria <- sommatoria + P[k_i1,k] * alfa_matrice[t-1, k_i1]
      }
      alfa_matrice[t,k] <- dens_calcolate[t,k] * sommatoria
      
    }
    
    #Normalizzo
    alfa_matrice[t,] <- alfa_matrice[t,] / sum(alfa_matrice[t,])
    
  }
  
  #procedo al campionamento
  z_new <- rep(NA, n)
  #zn, campiono da una distreta che associa a 1 prob an(1) e a 2 an(2) e a 3 a_n(3)
  z_new[n] <- sample(1:3,1, prob = alfa_matrice[n,])
  # f(zt |...) proporsionale alfa_t(zt)*P(zt,zt+1)
  
  for (t in (n-1):1) 
  {
    #guardo lo z che ho appena campionato
    z_futuro <- z_new[t+1]
    # estraggo la COLONNA della matrice P relativa allo stato futuro.
    # questa colonna contiene: [P(1->z_futuro), P(2->z_futuro), P(3->z_furuto)]
    prob_transizione_verso_futuro <- P[, z_futuro]
    #calcolo il peso per z_t: Alfa * Transizione
    #è una moltiplicazione elemento per elemento (1 con 1, 2 con 2), ovvero calcolo
    # P(z_t=k | z_{t+1}, y) propozionale a: alpha_t(k) * P(z_t=k -> z_{t+1})
    prob_t <- alfa_matrice[t, ] * prob_transizione_verso_futuro
    #controllo underflow
    if(sum(prob_t) == 0)
    {
      print("underflow, prob nulle")
    }
    #campiono
    z_new[t] <- sample(1:3,1, prob = prob_t)
  }
  catena_z[i,] <- z_new
  
  #print ogni 2000 iterazioni di controllo
  if(i %% 2000 == 0) 
  {
    print(paste("Iterazione", i, "di", B))
  }
  
  
}

# Rimuovo il burn-in
catena_m1_final <- catena_m1[(B_burn+1):B,]
catena_m2_final <- catena_m2[(B_burn+1):B,]
catena_m3_final <- catena_m3[(B_burn+1):B,]
catena_sigma1_final <- catena_sigma1[(B_burn+1):B,]
catena_sigma2_final <- catena_sigma2[(B_burn+1):B,]
catena_sigma3_final <- catena_sigma3[(B_burn+1):B,]
catena_pi1_final <- catena_pi1[(B_burn+1):B,]
catena_pi2_final <- catena_pi2[(B_burn+1):B,]
catena_pi3_final <- catena_pi3[(B_burn+1):B,]
catena_z_final <- catena_z[(B_burn+1):B,]

# diagnostica MCMC

# trace plot 

# --- Trace plots per μ (medie dei cluster) ---
trace_mu1_x_L3 <- ggplot(data.frame(iter = 1:nrow(catena_m1_final), 
                                    value = catena_m1_final[,1]), 
                         aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "steelblue") +
  theme_minimal() +
  labs(title = "Trace plot μ₁ (coordinata X) - L=3", x = "Iterazione", y = "μ₁ₓ")

trace_mu1_y_L3 <- ggplot(data.frame(iter = 1:nrow(catena_m1_final), 
                                    value = catena_m1_final[,2]), 
                         aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "coral") +
  theme_minimal() +
  labs(title = "Trace plot μ₁ (coordinata Y) - L=3", x = "Iterazione", y = "μ₁ᵧ")

trace_mu2_x_L3 <- ggplot(data.frame(iter = 1:nrow(catena_m2_final), 
                                    value = catena_m2_final[,1]), 
                         aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "steelblue") +
  theme_minimal() +
  labs(title = "Trace plot μ₂ (coordinata X) - L=3", x = "Iterazione", y = "μ₂ₓ")

trace_mu2_y_L3 <- ggplot(data.frame(iter = 1:nrow(catena_m2_final), 
                                    value = catena_m2_final[,2]), 
                         aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "coral") +
  theme_minimal() +
  labs(title = "Trace plot μ₂ (coordinata Y) - L=3", x = "Iterazione", y = "μ₂ᵧ")

trace_mu3_x_L3 <- ggplot(data.frame(iter = 1:nrow(catena_m3_final), 
                                    value = catena_m3_final[,1]), 
                         aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "steelblue") +
  theme_minimal() +
  labs(title = "Trace plot μ₃ (coordinata X) - L=3", x = "Iterazione", y = "μ₃ₓ")

trace_mu3_y_L3 <- ggplot(data.frame(iter = 1:nrow(catena_m3_final), 
                                    value = catena_m3_final[,2]), 
                         aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "coral") +
  theme_minimal() +
  labs(title = "Trace plot μ₃ (coordinata Y) - L=3", x = "Iterazione", y = "μ₃ᵧ")

# Salvataggio
grid_mu_L3 <- grid.arrange(trace_mu1_x_L3, trace_mu1_y_L3, 
                           trace_mu2_x_L3, trace_mu2_y_L3,
                           trace_mu3_x_L3, trace_mu3_y_L3, 
                           ncol = 2)
ggsave("trace_plots_mu_L3.png", grid_mu_L3, width = 14, height = 15, dpi = 300)


# --- Trace plots per Σ (matrici di covarianza) ---
trace_sigma1_11_L3 <- ggplot(data.frame(iter = 1:nrow(catena_sigma1_final), 
                                        value = catena_sigma1_final[,1]), 
                             aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkgreen") +
  theme_minimal() +
  labs(title = "Trace plot Σ₁[1,1] (var X) - L=3", x = "Iterazione", y = "Σ₁[1,1]")

trace_sigma1_22_L3 <- ggplot(data.frame(iter = 1:nrow(catena_sigma1_final), 
                                        value = catena_sigma1_final[,4]), 
                             aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkgreen") +
  theme_minimal() +
  labs(title = "Trace plot Σ₁[2,2] (var Y) - L=3", x = "Iterazione", y = "Σ₁[2,2]")

trace_sigma1_12_L3 <- ggplot(data.frame(iter = 1:nrow(catena_sigma1_final), 
                                        value = catena_sigma1_final[,2]), 
                             aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "purple") +
  theme_minimal() +
  labs(title = "Trace plot Σ₁[1,2] (cov) - L=3", x = "Iterazione", y = "Σ₁[1,2]")

trace_sigma2_11_L3 <- ggplot(data.frame(iter = 1:nrow(catena_sigma2_final), 
                                        value = catena_sigma2_final[,1]), 
                             aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkgreen") +
  theme_minimal() +
  labs(title = "Trace plot Σ₂[1,1] (var X) - L=3", x = "Iterazione", y = "Σ₂[1,1]")

trace_sigma2_22_L3 <- ggplot(data.frame(iter = 1:nrow(catena_sigma2_final), 
                                        value = catena_sigma2_final[,4]), 
                             aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkgreen") +
  theme_minimal() +
  labs(title = "Trace plot Σ₂[2,2] (var Y) - L=3", x = "Iterazione", y = "Σ₂[2,2]")

trace_sigma2_12_L3 <- ggplot(data.frame(iter = 1:nrow(catena_sigma2_final), 
                                        value = catena_sigma2_final[,2]), 
                             aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "purple") +
  theme_minimal() +
  labs(title = "Trace plot Σ₂[1,2] (cov) - L=3", x = "Iterazione", y = "Σ₂[1,2]")

trace_sigma3_11_L3 <- ggplot(data.frame(iter = 1:nrow(catena_sigma3_final), 
                                        value = catena_sigma3_final[,1]), 
                             aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkgreen") +
  theme_minimal() +
  labs(title = "Trace plot Σ₃[1,1] (var X) - L=3", x = "Iterazione", y = "Σ₃[1,1]")

trace_sigma3_22_L3 <- ggplot(data.frame(iter = 1:nrow(catena_sigma3_final), 
                                        value = catena_sigma3_final[,4]), 
                             aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkgreen") +
  theme_minimal() +
  labs(title = "Trace plot Σ₃[2,2] (var Y) - L=3", x = "Iterazione", y = "Σ₃[2,2]")

trace_sigma3_12_L3 <- ggplot(data.frame(iter = 1:nrow(catena_sigma3_final), 
                                        value = catena_sigma3_final[,2]), 
                             aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "purple") +
  theme_minimal() +
  labs(title = "Trace plot Σ₃[1,2] (cov) - L=3", x = "Iterazione", y = "Σ₃[1,2]")

grid_sigma_L3 <- grid.arrange(trace_sigma1_11_L3, trace_sigma1_22_L3, trace_sigma1_12_L3,
                              trace_sigma2_11_L3, trace_sigma2_22_L3, trace_sigma2_12_L3,
                              trace_sigma3_11_L3, trace_sigma3_22_L3, trace_sigma3_12_L3,
                              ncol = 3)
ggsave("trace_plots_sigma_L3.png", grid_sigma_L3, width = 18, height = 12, dpi = 300)


# --- Trace plots per π (probabilità di transizione) ---
trace_pi1_1_L3 <- ggplot(data.frame(iter = 1:nrow(catena_pi1_final), 
                                    value = catena_pi1_final[,1]), 
                         aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkorange") +
  theme_minimal() +
  labs(title = "Trace plot π₁[1] (prob 1→1) - L=3", x = "Iterazione", y = "π₁[1]")

trace_pi1_2_L3 <- ggplot(data.frame(iter = 1:nrow(catena_pi1_final), 
                                    value = catena_pi1_final[,2]), 
                         aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkorange") +
  theme_minimal() +
  labs(title = "Trace plot π₁[2] (prob 1→2) - L=3", x = "Iterazione", y = "π₁[2]")

trace_pi1_3_L3 <- ggplot(data.frame(iter = 1:nrow(catena_pi1_final), 
                                    value = catena_pi1_final[,3]), 
                         aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkorange") +
  theme_minimal() +
  labs(title = "Trace plot π₁[3] (prob 1→3) - L=3", x = "Iterazione", y = "π₁[3]")

trace_pi2_1_L3 <- ggplot(data.frame(iter = 1:nrow(catena_pi2_final), 
                                    value = catena_pi2_final[,1]), 
                         aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkred") +
  theme_minimal() +
  labs(title = "Trace plot π₂[1] (prob 2→1) - L=3", x = "Iterazione", y = "π₂[1]")

trace_pi2_2_L3 <- ggplot(data.frame(iter = 1:nrow(catena_pi2_final), 
                                    value = catena_pi2_final[,2]), 
                         aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkred") +
  theme_minimal() +
  labs(title = "Trace plot π₂[2] (prob 2→2) - L=3", x = "Iterazione", y = "π₂[2]")

trace_pi2_3_L3 <- ggplot(data.frame(iter = 1:nrow(catena_pi2_final), 
                                    value = catena_pi2_final[,3]), 
                         aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkred") +
  theme_minimal() +
  labs(title = "Trace plot π₂[3] (prob 2→3) - L=3", x = "Iterazione", y = "π₂[3]")

trace_pi3_1_L3 <- ggplot(data.frame(iter = 1:nrow(catena_pi3_final), 
                                    value = catena_pi3_final[,1]), 
                         aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkblue") +
  theme_minimal() +
  labs(title = "Trace plot π₃[1] (prob 3→1) - L=3", x = "Iterazione", y = "π₃[1]")

trace_pi3_2_L3 <- ggplot(data.frame(iter = 1:nrow(catena_pi3_final), 
                                    value = catena_pi3_final[,2]), 
                         aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkblue") +
  theme_minimal() +
  labs(title = "Trace plot π₃[2] (prob 3→2) - L=3", x = "Iterazione", y = "π₃[2]")

trace_pi3_3_L3 <- ggplot(data.frame(iter = 1:nrow(catena_pi3_final), 
                                    value = catena_pi3_final[,3]), 
                         aes(x = iter, y = value)) +
  geom_line(alpha = 0.7, color = "darkblue") +
  theme_minimal() +
  labs(title = "Trace plot π₃[3] (prob 3→3) - L=3", x = "Iterazione", y = "π₃[3]")

grid_pi_L3 <- grid.arrange(trace_pi1_1_L3, trace_pi1_2_L3, trace_pi1_3_L3,
                           trace_pi2_1_L3, trace_pi2_2_L3, trace_pi2_3_L3,
                           trace_pi3_1_L3, trace_pi3_2_L3, trace_pi3_3_L3,
                           ncol = 3)
ggsave("trace_plots_pi_L3.png", grid_pi_L3, width = 18, height = 12, dpi = 300)

#ESS

# Converto in oggetti mcmc
mcmc_mu1_x_L3 <- mcmc(catena_m1_final[,1])
mcmc_mu1_y_L3 <- mcmc(catena_m1_final[,2])
mcmc_mu2_x_L3 <- mcmc(catena_m2_final[,1])
mcmc_mu2_y_L3 <- mcmc(catena_m2_final[,2])
mcmc_mu3_x_L3 <- mcmc(catena_m3_final[,1])
mcmc_mu3_y_L3 <- mcmc(catena_m3_final[,2])

mcmc_sigma1_11_L3 <- mcmc(catena_sigma1_final[,1])
mcmc_sigma1_12_L3 <- mcmc(catena_sigma1_final[,2])
mcmc_sigma1_22_L3 <- mcmc(catena_sigma1_final[,4])
mcmc_sigma2_11_L3 <- mcmc(catena_sigma2_final[,1])
mcmc_sigma2_12_L3 <- mcmc(catena_sigma2_final[,2])
mcmc_sigma2_22_L3 <- mcmc(catena_sigma2_final[,4])
mcmc_sigma3_11_L3 <- mcmc(catena_sigma3_final[,1])
mcmc_sigma3_12_L3 <- mcmc(catena_sigma3_final[,2])
mcmc_sigma3_22_L3 <- mcmc(catena_sigma3_final[,4])

mcmc_pi1_1_L3 <- mcmc(catena_pi1_final[,1])
mcmc_pi1_2_L3 <- mcmc(catena_pi1_final[,2])
mcmc_pi1_3_L3 <- mcmc(catena_pi1_final[,3])
mcmc_pi2_1_L3 <- mcmc(catena_pi2_final[,1])
mcmc_pi2_2_L3 <- mcmc(catena_pi2_final[,2])
mcmc_pi2_3_L3 <- mcmc(catena_pi2_final[,3])
mcmc_pi3_1_L3 <- mcmc(catena_pi3_final[,1])
mcmc_pi3_2_L3 <- mcmc(catena_pi3_final[,2])
mcmc_pi3_3_L3 <- mcmc(catena_pi3_final[,3])

# Calcolo ESS
ess_mu1_x_L3 <- effectiveSize(mcmc_mu1_x_L3)
ess_mu1_y_L3 <- effectiveSize(mcmc_mu1_y_L3)
ess_mu2_x_L3 <- effectiveSize(mcmc_mu2_x_L3)
ess_mu2_y_L3 <- effectiveSize(mcmc_mu2_y_L3)
ess_mu3_x_L3 <- effectiveSize(mcmc_mu3_x_L3)
ess_mu3_y_L3 <- effectiveSize(mcmc_mu3_y_L3)

ess_sigma1_11_L3 <- effectiveSize(mcmc_sigma1_11_L3)
ess_sigma1_12_L3 <- effectiveSize(mcmc_sigma1_12_L3)
ess_sigma1_22_L3 <- effectiveSize(mcmc_sigma1_22_L3)
ess_sigma2_11_L3 <- effectiveSize(mcmc_sigma2_11_L3)
ess_sigma2_12_L3 <- effectiveSize(mcmc_sigma2_12_L3)
ess_sigma2_22_L3 <- effectiveSize(mcmc_sigma2_22_L3)
ess_sigma3_11_L3 <- effectiveSize(mcmc_sigma3_11_L3)
ess_sigma3_12_L3 <- effectiveSize(mcmc_sigma3_12_L3)
ess_sigma3_22_L3 <- effectiveSize(mcmc_sigma3_22_L3)

ess_pi1_1_L3 <- effectiveSize(mcmc_pi1_1_L3)
ess_pi1_2_L3 <- effectiveSize(mcmc_pi1_2_L3)
ess_pi1_3_L3 <- effectiveSize(mcmc_pi1_3_L3)
ess_pi2_1_L3 <- effectiveSize(mcmc_pi2_1_L3)
ess_pi2_2_L3 <- effectiveSize(mcmc_pi2_2_L3)
ess_pi2_3_L3 <- effectiveSize(mcmc_pi2_3_L3)
ess_pi3_1_L3 <- effectiveSize(mcmc_pi3_1_L3)
ess_pi3_2_L3 <- effectiveSize(mcmc_pi3_2_L3)
ess_pi3_3_L3 <- effectiveSize(mcmc_pi3_3_L3)


#Output diagnostiche e statistiche a posteriori


cat("=== DIAGNOSTICHE MCMC (L=3) ===\n")


cat("\n--- Effective Sample Size (ESS) ---\n")
cat("\nMedie dei cluster (μ):\n")
cat(sprintf("  μ₁ₓ: %.0f\n", ess_mu1_x_L3))
cat(sprintf("  μ₁ᵧ: %.0f\n", ess_mu1_y_L3))
cat(sprintf("  μ₂ₓ: %.0f\n", ess_mu2_x_L3))
cat(sprintf("  μ₂ᵧ: %.0f\n", ess_mu2_y_L3))
cat(sprintf("  μ₃ₓ: %.0f\n", ess_mu3_x_L3))
cat(sprintf("  μ₃ᵧ: %.0f\n", ess_mu3_y_L3))

cat("\nMatrici di covarianza (Σ):\n")
cat("  Sigma 1:\n")
cat(sprintf("    Var(X): %.0f\n", ess_sigma1_11_L3))
cat(sprintf("    Cov(X,Y): %.0f\n", ess_sigma1_12_L3))
cat(sprintf("    Var(Y): %.0f\n", ess_sigma1_22_L3))
cat("  Sigma 2:\n")
cat(sprintf("    Var(X): %.0f\n", ess_sigma2_11_L3))
cat(sprintf("    Cov(X,Y): %.0f\n", ess_sigma2_12_L3))
cat(sprintf("    Var(Y): %.0f\n", ess_sigma2_22_L3))
cat("  Sigma 3:\n")
cat(sprintf("    Var(X): %.0f\n", ess_sigma3_11_L3))
cat(sprintf("    Cov(X,Y): %.0f\n", ess_sigma3_12_L3))
cat(sprintf("    Var(Y): %.0f\n", ess_sigma3_22_L3))

cat("\nMatrice di transizione (π):\n")
cat("  Da stato 1:\n")
cat(sprintf("    π₁[1→1]: %.0f\n", ess_pi1_1_L3))
cat(sprintf("    π₁[1→2]: %.0f\n", ess_pi1_2_L3))
cat(sprintf("    π₁[1→3]: %.0f\n", ess_pi1_3_L3))
cat("  Da stato 2:\n")
cat(sprintf("    π₂[2→1]: %.0f\n", ess_pi2_1_L3))
cat(sprintf("    π₂[2→2]: %.0f\n", ess_pi2_2_L3))
cat(sprintf("    π₂[2→3]: %.0f\n", ess_pi2_3_L3))
cat("  Da stato 3:\n")
cat(sprintf("    π₃[3→1]: %.0f\n", ess_pi3_1_L3))
cat(sprintf("    π₃[3→2]: %.0f\n", ess_pi3_2_L3))
cat(sprintf("    π₃[3→3]: %.0f\n", ess_pi3_3_L3))



cat("=== STATISTICHE POSTERIORI (L=3) ===\n")


cat("\n--- μ₁ (centro cluster 1) ---\n")
cat(sprintf("  X: %.3f [%.3f, %.3f]\n", 
            mean(catena_m1_final[,1]),
            quantile(catena_m1_final[,1], 0.025),
            quantile(catena_m1_final[,1], 0.975)))
cat(sprintf("  Y: %.3f [%.3f, %.3f]\n", 
            mean(catena_m1_final[,2]),
            quantile(catena_m1_final[,2], 0.025),
            quantile(catena_m1_final[,2], 0.975)))

cat("\n--- μ₂ (centro cluster 2) ---\n")
cat(sprintf("  X: %.3f [%.3f, %.3f]\n", 
            mean(catena_m2_final[,1]),
            quantile(catena_m2_final[,1], 0.025),
            quantile(catena_m2_final[,1], 0.975)))
cat(sprintf("  Y: %.3f [%.3f, %.3f]\n", 
            mean(catena_m2_final[,2]),
            quantile(catena_m2_final[,2], 0.025),
            quantile(catena_m2_final[,2], 0.975)))

cat("\n--- μ₃ (centro cluster 3) ---\n")
cat(sprintf("  X: %.3f [%.3f, %.3f]\n", 
            mean(catena_m3_final[,1]),
            quantile(catena_m3_final[,1], 0.025),
            quantile(catena_m3_final[,1], 0.975)))
cat(sprintf("  Y: %.3f [%.3f, %.3f]\n", 
            mean(catena_m3_final[,2]),
            quantile(catena_m3_final[,2], 0.025),
            quantile(catena_m3_final[,2], 0.975)))

cat("\n--- Probabilità di transizione π₁ (da stato 1) ---\n")
cat(sprintf("  1→1: %.3f [%.3f, %.3f]\n", 
            mean(catena_pi1_final[,1]),
            quantile(catena_pi1_final[,1], 0.025),
            quantile(catena_pi1_final[,1], 0.975)))
cat(sprintf("  1→2: %.3f [%.3f, %.3f]\n", 
            mean(catena_pi1_final[,2]),
            quantile(catena_pi1_final[,2], 0.025),
            quantile(catena_pi1_final[,2], 0.975)))
cat(sprintf("  1→3: %.3f [%.3f, %.3f]\n", 
            mean(catena_pi1_final[,3]),
            quantile(catena_pi1_final[,3], 0.025),
            quantile(catena_pi1_final[,3], 0.975)))

cat("\n--- Probabilità di transizione π₂ (da stato 2) ---\n")
cat(sprintf("  2→1: %.3f [%.3f, %.3f]\n", 
            mean(catena_pi2_final[,1]),
            quantile(catena_pi2_final[,1], 0.025),
            quantile(catena_pi2_final[,1], 0.975)))
cat(sprintf("  2→2: %.3f [%.3f, %.3f]\n", 
            mean(catena_pi2_final[,2]),
            quantile(catena_pi2_final[,2], 0.025),
            quantile(catena_pi2_final[,2], 0.975)))
cat(sprintf("  2→3: %.3f [%.3f, %.3f]\n", 
            mean(catena_pi2_final[,3]),
            quantile(catena_pi2_final[,3], 0.025),
            quantile(catena_pi2_final[,3], 0.975)))

cat("\n--- Probabilità di transizione π₃ (da stato 3) ---\n")
cat(sprintf("  3→1: %.3f [%.3f, %.3f]\n", 
            mean(catena_pi3_final[,1]),
            quantile(catena_pi3_final[,1], 0.025),
            quantile(catena_pi3_final[,1], 0.975)))
cat(sprintf("  3→2: %.3f [%.3f, %.3f]\n", 
            mean(catena_pi3_final[,2]),
            quantile(catena_pi3_final[,2], 0.025),
            quantile(catena_pi3_final[,2], 0.975)))
cat(sprintf("  3→3: %.3f [%.3f, %.3f]\n", 
            mean(catena_pi3_final[,3]),
            quantile(catena_pi3_final[,3], 0.025),
            quantile(catena_pi3_final[,3], 0.975)))

cat("\n--- Σ₁ (matrice di covarianza cluster 1) ---\n")
cat(sprintf("  Σ₁[1,1] (var X): %.4f [%.4f, %.4f]\n", 
            mean(catena_sigma1_final[,1]),
            quantile(catena_sigma1_final[,1], 0.025),
            quantile(catena_sigma1_final[,1], 0.975)))
cat(sprintf("  Σ₁[2,2] (var Y): %.4f [%.4f, %.4f]\n", 
            mean(catena_sigma1_final[,4]),
            quantile(catena_sigma1_final[,4], 0.025),
            quantile(catena_sigma1_final[,4], 0.975)))
cat(sprintf("  Σ₁[1,2] (cov XY): %.4f [%.4f, %.4f]\n", 
            mean(catena_sigma1_final[,2]),
            quantile(catena_sigma1_final[,2], 0.025),
            quantile(catena_sigma1_final[,2], 0.975)))

cat("\n--- Σ₂ (matrice di covarianza cluster 2) ---\n")
cat(sprintf("  Σ₂[1,1] (var X): %.4f [%.4f, %.4f]\n", 
            mean(catena_sigma2_final[,1]),
            quantile(catena_sigma2_final[,1], 0.025),
            quantile(catena_sigma2_final[,1], 0.975)))
cat(sprintf("  Σ₂[2,2] (var Y): %.4f [%.4f, %.4f]\n", 
            mean(catena_sigma2_final[,4]),
            quantile(catena_sigma2_final[,4], 0.025),
            quantile(catena_sigma2_final[,4], 0.975)))
cat(sprintf("  Σ₂[1,2] (cov XY): %.4f [%.4f, %.4f]\n", 
            mean(catena_sigma2_final[,2]),
            quantile(catena_sigma2_final[,2], 0.025),
            quantile(catena_sigma2_final[,2], 0.975)))

cat("\n--- Σ₃ (matrice di covarianza cluster 3) ---\n")
cat(sprintf("  Σ₃[1,1] (var X): %.4f [%.4f, %.4f]\n", 
            mean(catena_sigma3_final[,1]),
            quantile(catena_sigma3_final[,1], 0.025),
            quantile(catena_sigma3_final[,1], 0.975)))
cat(sprintf("  Σ₃[2,2] (var Y): %.4f [%.4f, %.4f]\n", 
            mean(catena_sigma3_final[,4]),
            quantile(catena_sigma3_final[,4], 0.025),
            quantile(catena_sigma3_final[,4], 0.975)))
cat(sprintf("  Σ₃[1,2] (cov XY): %.4f [%.4f, %.4f]\n", 
            mean(catena_sigma3_final[,2]),
            quantile(catena_sigma3_final[,2], 0.025),
            quantile(catena_sigma3_final[,2], 0.975)))


# Grafici per clustering

# Calcolo la moda a posteriori per ogni z_i
z_moda_L3 <- rep(NA, n)

for (i in 1:n) {
  conteggio1 <- sum(catena_z_final[,i] == 1)
  conteggio2 <- sum(catena_z_final[,i] == 2)
  conteggio3 <- sum(catena_z_final[,i] == 3)
  
  # Determino quale stato ha il conteggio maggiore
  conteggi <- c(conteggio1, conteggio2, conteggio3)
  #Prima avevo eseguito il conto a mano ora uso whith.max per evitare di fare
  # tutta la casistica
  z_moda_L3[i] <- which.max(conteggi)
}

# Aggiungo al dataframe
data_dogs$z_moda_L3 <- as.factor(z_moda_L3)

# Grafico spaziale con colorazione per cluster stimato (L=3)
plot_spaziale_cluster_L3 <- ggplot(data_dogs, aes(x = X, y = Y)) +
  geom_path(alpha = 0.3, color = "gray70") +
  geom_point(aes(color = z_moda_L3), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("1" = "steelblue", "2" = "coral", "3" = "darkgreen"),
                     name = "Cluster stimato",
                     labels = c("Stato 1", "Stato 2", "Stato 3")) +
  theme_minimal() +
  labs(title = "Traiettoria del cane con clustering HMM (L=3)", 
       x = "Coordinata X", 
       y = "Coordinata Y",
       subtitle = "Colorazione basata sulla moda a posteriori di z_i")

print(plot_spaziale_cluster_L3)
ggsave("traiettoria_cane_cluster_L3.png", plot_spaziale_cluster_L3, 
       width = 15, height = 10, dpi = 300)


# Serie temporali con colorazione per cluster (L=3)
grafico_x_cluster_L3 <- ggplot(data_dogs, aes(x = time, y = X)) +
  geom_line(color = "gray70", alpha = 0.5) +
  geom_point(aes(color = z_moda_L3), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("1" = "steelblue", "2" = "coral", "3" = "darkgreen"),
                     name = "Cluster",
                     labels = c("Stato 1", "Stato 2", "Stato 3")) +
  theme_minimal() +
  labs(title = "Coordinata X nel tempo - Clustering HMM (L=3)", 
       x = "Tempo (indice)", 
       y = "Coordinata X")

grafico_y_cluster_L3 <- ggplot(data_dogs, aes(x = time, y = Y)) +
  geom_line(color = "gray70", alpha = 0.5) +
  geom_point(aes(color = z_moda_L3), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("1" = "steelblue", "2" = "coral", "3" = "darkgreen"),
                     name = "Cluster",
                     labels = c("Stato 1", "Stato 2", "Stato 3")) +
  theme_minimal() +
  labs(title = "Coordinata Y nel tempo - Clustering HMM (L=3)", 
       x = "Tempo (indice)", 
       y = "Coordinata Y")

# Combino i due grafici
grafici_temporali_cluster_L3 <- grid.arrange(grafico_x_cluster_L3, grafico_y_cluster_L3, 
                                             ncol = 2)

ggsave("coordinate_tempo_cluster_L3.png", grafici_temporali_cluster_L3, 
       width = 14, height = 6, dpi = 300)







