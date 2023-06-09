### S�ries Temporais com R: Analise ARIMA do Consumo do Varejo em MA
### Elabora��o: Paulo Roberto Carneiro de Sa

library(BETS)

# Pegando as series a partir do site do Banco Central do Brasil
# �ndice de volume de vendas no varejo Total do maranh�o
# mensal a partir de jan/2000 at� fev/2020 
# 242 observa��es mensais
# pacotes necessarios:
library(BETS)
library(dygraphs)
library(ggplot2)
library(fpp2)
library(FitAR)
library(tseries)
library(forecast)

######################################################################################################
varejoma <- BETSget(1463) 
print(varejoma)
class(varejoma)
dput(varejoma)  # op��o para ter os dados como na structure abaixo

#A rotina de dados obtidos pelo BETS j� retorna a s�rie em formato ts, ou seja, s�rie temporal. 
#Farei ent�o a cria��o de uma s�rie em diferen�as para observar o comportamento da s�rie em n�vel e em diferen�as.
#Inicialmente olharei as estat�sticas descritivas da s�rie. Em seguida farei um plot b�sico da s�rie e o plot pelo pacote dygraphs,
#�til para ver os pontos de picos e momentos espec�ficos.

dvarejo <- diff(varejoma)
# estatisticas basicas
summary(varejoma)

### plot basico lembrar que em class(), ele j� indicou que era ts = serie
# temporal
plot(varejoma)

# pelo pacote dygraph d� mais op��es
library(dygraphs)
help("dygraph")
dygraph(varejoma, main = "�ndice de volume de vendas no varejo total do Maranh�o <br> (Mensal)  (2011=100) BCB 1463") %>% 
  dyAxis("x", drawGrid = TRUE) %>% dyEvent("2005-1-01", "2005", labelLoc = "bottom") %>% 
  dyEvent("2015-1-01", "2015", labelLoc = "bottom") %>% dyEvent("2018-1-01", 
                                                                "2018", labelLoc = "bottom") %>% dyEvent("2019-2-20", "2020", labelLoc = "bottom") %>% 
  dyOptions(drawPoints = TRUE, pointSize = 2)

### Fun��o de Autocorrela��o (FAC) e Autocorrela��o parcial (FACp) com defasagem 36
# S�rie em n�vel - Usarei a rotina do Hyndman e Athanasopoulos (2018).

library(ggplot2)
varejo <- varejoma
varejo %>% ggtsdisplay(main = "")

### S�rie em primeira diferen�a

#Como a s�rie apresenta varia��es sazonais importantes, 
#assim como tend�ncia importante, claramente n�o-estacion�rias, 
#vou olhar tamb�m em primeira diferen�a.

varejo %>% diff() %>% ggtsdisplay(main = "S�rie varejo do MA em primeira diferen�a")

### Primeira diferen�a sazonal e ACF e PACF
#Farei agora a diferen�a sazonal.

varejo %>% diff(lag = 12) %>% ggtsdisplay(main = "S�rie varejo do MA em primeira diferen�a sazonal")

### Ajuste sazonal e ACF e PACF
#Aplicaremos o ajuste sazonal tipo STL (Seasonal Decomposition of Time Series by Loess) aos dados.

library(fpp2)
varejoadj <- varejo %>% stl(s.window = "periodic") %>% seasadj()
autoplot(varejoadj)

########################### analise ACF e PACF ############################################################# 
#E agora os plots de ACF e PACF na s�rie ajustada sazonalmente.
varejoadj %>% diff() %>% ggtsdisplay(main = "S�rie varejo do MA em primeira diferen�a e ajuste sazonal")
#Parece indicar para algum AR de ordem 4 na PACF,
#alguma sazonalidade marcante na ACF (indicando um termo SMA), 
#e indicacao de um ARIMA (4,1,1)(0,1,1)[12].

#Uma estima��o inicial desse modelo daria algo como: 
(fit <- Arima(varejoadj, order = c(4, 1, 1), seasonal = c(0, 1, 1)))
 
#############################################################################################################

autoplot(fit)

################### analise ARIMA ################################

checkresiduals(fit)

##################################################################

library(FitAR)  # pratico para plotar graficos de Q
# LBQPlot(res, lag.max = 30) res � a s�rie que temos interesse de realizar
# os testes lag.max � o n�mero de lags que se deseja imprimir o gr�fico
LBQPlot(residuals(fit), 36)

# checar normalidade dos res�duos
library(tseries)
jarque.bera.test(residuals(fit))

########################## ARIMA #########################################
#Uma estima��o automatizada pela fun��o auto.arima {forecast} padr�o indica um modelo diferente.
library(forecast)
fit2 <- auto.arima(varejo)
summary(fit2)

# p-values dos coeficientes
pvalues <- (1 - pnorm(abs(fit2$coef)/sqrt(diag(fit2$var.coef)))) * 2
pvalues

library(FitAR)  # pratico para plotar graficos de Q
LBQPlot(residuals(fit2), 36)

# checar normalidade dos res�duos
library(tseries)
jarque.bera.test(residuals(fit2))

######################## Modelo ARIMA autom�tico mais geral ###########################################
require(fpp2)
fit3 <- auto.arima(varejoma, stepwise = FALSE, approximation = FALSE)
summary(fit3)

# p-values dos coeficientes
pvalues <- (1 - pnorm(abs(fit3$coef)/sqrt(diag(fit3$var.coef)))) * 2
pvalues

##################### Checando res�duos ##################################################
#Os modelos podem ser considerados white noise, mas n�o normais.

checkresiduals(fit3)

require(FitAR)  # pratico para plotar graficos de Q
LBQPlot(residuals(fit3), 36)

# checar normalidade dos res�duos
require(tseries)
jarque.bera.test(residuals(fit3))

#Pode-se fazer o gr�fico dos forecasts deste modelo:
autoplot(forecast(fit3, h = 24), title = "Forecasts de volume de vendas no varejo total do Maranh�o", 
         xlab = "Ano", ylab = "�ndice (2011=100)")

######################################################################################################
# Teste de raiz unit�ria
# O teste de Dickey-Fuller aumentado (ADF) � nossa primeira escolha, 
# O padr�o � que H0: s�rie n�o estacion�ria ou s�rie tem raiz unit�ria.
ADF.test <- adf.test(varejoma)
ADF.test

######################################################################################################
# Teste ADF invertendo a hip�tese nula
#Neste caso, o padr�o � que H0: a s�rie � estacion�ria ou a s�rie n�o tem raiz unit�ria, 
#Neste caso, a hip�tese alternativa � que a s�rie � "explosiva".

ADF.test12 <- adf.test(varejoma, k = 12, alternative = c("explosive"))
ADF.test12

######################################################################################################
# ARIMA - Codigo do Hyndman e Athanasopoulos 
x <- varejoma
ns <- nsdiffs(x)  #numero de diferencas sazonais
ns

if (ns > 0) {
  xstar <- diff(x, lag = frequency(x), differences = ns)
} else {
  xstar <- x
}

nd <- ndiffs(xstar)  # numero de diferencas
nd

if (nd > 0) {
  xstar <- diff(xstar, differences = nd)
}
# xstar ser� a s�rie devidamente em diferen�as sazonais e diferen�as normais
autoplot(xstar)

auto.arima(xstar) # ver que o resultado � igual, mas agora nao tem diferencas

modelo <- auto.arima(xstar, stepwise = FALSE, approximation = FALSE)
summary(modelo)

autoplot(modelo)








































































































































