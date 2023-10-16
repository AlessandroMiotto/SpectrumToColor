# Spectrum To Color

Il *colorspace XYZ* è uno spazio-colori simile a quello *RGB* con la differenza di avere valori di risposta nella *color matching function* sempre positivi. La *Y* rappresenta la luminanza (che è la luminosità emessa per unità di superficie), la *X* rappresenta uno stimolo blu e la *Z* un'opportuna mistura di colori.  

Nel *colorspace* XYZ i colori sono ottenuti, nel caso riflessivo e trasmissivo, come gli integrali nella lunghezza d'onda $\lambda\in \left[380,780\right]$:
$$X=\frac{1}{N}\int_\lambda S(\lambda )I(\lambda )\bar{x}(\lambda )d\lambda$$
$$Y=\frac{1}{N}\int_\lambda S(\lambda )I(\lambda )\bar{y}(\lambda )d\lambda$$
$$Z=\frac{1}{N}\int_\lambda S(\lambda )I(\lambda )\bar{z}(\lambda )d\lambda$$
dove $I(\lambda )$ è la distribuzione spettrale dell'**illuminante** (nel nostro caso consideremo la luce solare, l'illuminante standard D65), $S(\lambda )$ è la **riflettanza/trasmittanza** dell'oggetto e infine $\bar{x}(\lambda )$ è la **color matching function** e rappresenta la risposta dell'occhio umano allo spettro elettromagnatico. $N$ è una costante di normalizzazione ed è pari a $N=\int_\lambda I(\lambda)\bar{y}(\lambda )d\lambda$.

Nel nostro caso tratteremo la trasmittanza e allora $S(\lambda)=T(\lambda)=\Phi_e^t/\Phi_e^i$ dove $\Phi_e^t$ è il flusso radiante trasmesso e $\Phi_e^i$ è il flusso radente incidente. La legge di Beer-Lambert permette di calcolare la trasmittanza che è $T=e^{-\tau}$ dove $\tau=\sum_i\sigma_i\int_0^l n_i(z)dz$ con $\sigma_i$ la sezione d'urto d'atenuazione e $n_i$ la densità nel mezzo dell'*i*-esima specie. Assumendo un mezzo isotropo si ottiene che $\tau=\sigma n l$ (dove $l$ è lo spessore del materiale). In generale troviamo che $\sigma n=\alpha(\lambda)$ è il coefficiente di assorbimento e allora troveremo che:
$$X=\frac{1}{N}\int_{\lambda}e^{-\alpha(\lambda)l}I(\lambda)\bar{x}(\lambda)d\lambda$$

Il passaggio dal color space XYZ a quello RGB (in particolare quello sRGB) è ottenibile tramite trasformazioni lineari introdotte da opportune [matrici](http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html) che dipendono anche dall'illuminante scelto. Per lo spazio colore sRGB con l'illuminante standard D65, la trasformazione è:
```math
\begin{pmatrix} r \\ g \\ b \end{pmatrix}=\begin{pmatrix} 3.2404542 & -1.5371385 & -0.4985314  \\ -0.9692660 & 1.8760108 & 0.0415560 \\ 0.0556434 & -0.2040259 & 1.0572252 \end{pmatrix} \begin{pmatrix} X \\ Y \\ Z \end{pmatrix}
```
I valori $(r,g,b)$ sono detti canali RGB lineari e per passare ai non lineari $(R,G,B)$, è necessario applicare la cosiddetta *gamma correction*:
```math
C_\text{sRGB} = \begin{cases}
12.92C_\text{linear}, & C_\text{linear} \le 0.0031308 \\[5mu]
1.055C_\text{linear}^{1/2.4}-0.055, & C_\text{linear} > 0.0031308
\end{cases}
```
dove $C_\text{sRGB}=\{R,G,B\}$ e $C_\text{linear}=\{r,g,b\}$. Generalmente $C_{sRGB}\in[0,1]$ ma è possibile ottenere valori negativi e ciò è dovuto al fatto che la gamma di colore (*gamut*) dello spazio XYZ è più ampia rispetto a quella sRGB. Valori che eccedono lo spazio sRGB non sono accuratamente rappresentabili su uno schermo.
