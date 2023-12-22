# Projet2 BICHON Pierre & GONTIER Lucien

Ce projet a été réalisé par la collaboration de Pierre Bichon et Lucien Gontier.
$ \rho_0(x) = \left{\begin{cases} 0 & \text{si } x < 1 \ 1 & \text{si } 1 < x < 2 \ \end{cases}\right . $


Il s'inscrit dans le cadre du cours de calcul scientifique numérique du Master 2 Modélisation, Analyse numérique et Calcul Scientifique de Nantes université.

Il vise à retrouver numériquement une frontière d'un problème de Cauchy à l'aide des conditions aux bords partiellement surdéterminées, et ce en implémentant trois différentes méthodes pour résoudre une intégrale de Fredholm essentiel à la résolution du problème.

Pour en savoir plus nous vous invitons à lire notre rapport qui se trouve dans le dossier intitulé "rapport".

Ici nous ne parlerons que du programme en lui même.



## Structure du projet
Notre projet est constitué de plusiseurs dossiers et fichiers.

Une fois dans le dossier Projet2_BICHON_GONTIER vous devriez trouver les dossiers:

- headers
- prog
- dat
- images
- obj
- gnuplot
- rapport

Les fichiers :

- makefile
- README.MD

Et l'executable :

- exe

### Dossier "headers"
Le dossier headers comme son nom l'indique comporte tous les headers .h des fichier .c du dossier "prog" de notre programme, dans cette rubrique, nous détaillerons chaques fonctions simplifier la lecture et l'utilisation de ces dernières.

Comme il risque d'il y avoir beaucoup dinformations, nous allons faire des sous-sections pour chaques fichier.h :

#### Constante.h
Ce fichier contient nos appels aux librairies extérieurs et le paramétrage de nos constantes :
- pi : $\pi$ en double précision.
- H : le H définissant la hauteur selon y à laquelle la frontière $\Gamma_3$ de $D$ se trouve.
- L : le L définissant la largeur selon x à laquelle la frontière $\Gamma_2$ de $\Omega$ se trouve.
- N : le nombre de points et poids utilisé pour notre quadrature de Gauss Legendre.
- N1 : le nombre d'itération pour la somme dans k(x,t) et h(x).
- N2 : nombre d'itérations pour le calcul de $\tilde T(x,y)$
- k : permet de changer le $f_0$ et $U_{ex}$
- num : permet aussi un changement sur $f_0$ et $U_{ex}$

Exemple pour num et k, premier cas si num = 0 :

$U_{ex}(x,y)=cos(\frac{k\pi}{L}x)cosh(\frac{k\pi}{L}y)$ et $f_0(x)=cos(\frac{k\pi}{L}x)$

Et si num = 1 :

$U_{ex}(x,y)=cos(\frac{k\pi}{L}x)cosh(\frac{k\pi}{L}y)\frac{x}{k}$ et $f_0(x)=cos(\frac{k\pi}{L}x)\frac{x}{k}$

#### fonctions.h
Ce fichier contient toutes les fonctions mathématiques utilisées dans ce projet :
```c
double Uex(double x, double y);
```
Cette fonction renvoie la valeur de $U_{ex}$ aux points x et y.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
| `x`       | `double`  |x de [0,L]     |
| `y`       | `double`  | y de [O,H]    |

```c
double f0(double x);
```
Cette fonction renvoie la valeur de $f_0$ au point x.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
| `x`       | `double`  |x de [0,L]     |

```c
double q0(double x);
```
Cette fonction renvoie la valeur de $q_0$ au point x.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
| `x`       | `double`  |x de [0,L]     |


```c
 double alpha_i(double x,int i);
```
Cette fonction renvoie la valeur $\alpha_i(x)=\frac{2\pi}{L^{2}}\frac{icos(\frac{i\pi}{L}x)}{sinh(\frac{i\pi}{L}H)}$ pour $1\leq i \leq N1$, et $\alpha_{N1+1}(x) =\frac{1}{HL}$, pour la méthode des noyaux séparables.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
| `x`       | `double`  |x de [0,L] de $\alpha_i(x)$    |
| `i`       |`int`    | i entre 1 et N1+1 de  $\alpha_i$ |

```c
double beta_i(double t, int i);
```
Cette fonction renvoie la valeur $\beta_i(t)=cos(\frac{i\pi}{L}t)$ pour $1\leq i \leq N1$ et  $\beta_{N1+1}(x) =1$, pour la méthode des noyaux séparables.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
| `t`       | `double`  |x de [0,L] de $\beta_i(t)$    |
| `i`       |`int`    | i entre 1 et N1+1 de  $\beta_i$ |


```c
double h(double x);
```
Cette fonction renvoie la valeur de $h$ au point x.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
| `x`       | `double`  |x de [0,L]     |

#### Gauss_Legendre.h
```c
double integrale_GL(double (*f)(double), double a, double b,int m);
```
Cette fonction renvoie l'approximation de valeur de l'intégrale $\displaystyle\int_a^b f(x)cos(\frac{m\pi}{L}x)dx$ par la quadrature de Gauss Legendre.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
| `f`       | `fonction de double`  |la fonction que l'on souhaite intégrer contre  $cos(\frac{m\pi}{L}x)$  |
| `a`       |`double` | borne inférieur de l'intégrale |
| `b`       |`double` | borne supérieur de l'intégrale |
| `m`       |`int`    | m de $cos(\frac{m\pi}{L}x)$ |

```c
double GL2( double tab[N], double a, double b,int m);
```
Cette fonction renvoie l'approximation de valeur de l'intégrale $\displaystyle\int_a^b f(x)cos(\frac{m\pi}{L}x)dx$ par la quadrature de Gauss Legendre, avec tab les valeur de f calculés aux points $\frac{b-a}{2}X[i] +\frac{b+a}{2}$, avec X[i] les points de Gauss Legendre originaux pour N points.

Cette fonction demande un calcul au préalable de notre tab[N] aux bons points.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
| `tab`     | `tableau de N double`  |le tableau evalué aux bon points de la fonction qu'on souhaite intégrer contre  $cos(\frac{m\pi}{L}x)$  |
| `a`       |`double` | borne inférieur de l'intégrale |
|`b`        |`double` | borne supérieur de l'intégrale |
| `m`       |`int`    | m de $cos(\frac{m\pi}{L}x)$ |

#### OperationsMatrices.h
Dans ce programme se trouve les fonctions permettant les opérations sur les matrices nécessaires au fonctionnement de la méthode d'approximation succesive.
```c
void mult_mat_scal(double lambda,double A[N1+1][N1+1]);
```
Cette fonction renvoie la la valeur du produit matriciel $\lambda\times A$, dans $A$.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
|`lambda`   |`double` | scalaire multipliant la matrice A |
|`A`        |`tableau de (N+1)x(N+1) doubles` | matrice de taille (N+1)²               |

```c
void sous_matrice(double In[N1+1][N1+1],double A[N1+1][N1+1]);
```
Cette fonction renvoie la la valeur du produit matriciel $I_n - A$, dans $A$.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
|`In`       |`tableau de (N+1)x(N+1) doubles` | matrice identité de taille  (N+1)²           |
|`A`        |`tableau de (N+1)x(N+1) doubles` | matrice de taille (N+1)²            |


#### methodes.h

Dans ce programme on va retrouver les méthodes que nous avons utiliser pour approximer $f_3^\alpha$.

```c
double adomain_v1( double x,double alpha, int n);
```
Cette fonction renvoit la valeur de $f_3^\alpha$ au point x, en utilisant la méthode d'Adomain, elle correspond à la première des deux versions. (voir rapport)
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
|`x`        |`double` | x entre [0,L] de $f_3^\alpha (x)$              |
| `alpha`   |`double` | $\alpha$ de $f_3^\alpha$  |
| `n`       |`int`    | nombre d'itération pour la méthode|


```c
double Kn( double x, int s[], int n);
```
Cette fonction renvoie K(s0,s1,...,sn) définie dans le rapport aux points de Gauss Legendre definie par le vecteur s, le tout multiplie par le produit des poids w[s[i]]
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
|`x`        |`double` | x entre [0,L] de $f_3^\alpha (x)$              |
| `s`[]     |`int`    | points pour Gauss-Legendre |
| `n`       |`int`    | numéro de l'itération|


```c
double U_n( double x, int n);
```
Cette fontion renvoie la valeur de Un(x) sans le terme dépendant de alpha pour chaque n
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
|`x`        |`double` | x entre [0,L] de $f_3^\alpha (x)$              |
| `n`       |`int`    | numéro de l'itération|


```c
double adomain_v2( double x,double alpha, int n);
```
Cette fonction renvoit la valeur de $f_3^\alpha$ au point x, en utilisant la méthode d'Adomain, elle correspond à la seconde des deux versions. (voir rapport)
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
|`x`        |`double` | x entre [0,L] de $f_3^\alpha (x)$              |
| `alpha`   |`double` | $\alpha$ de $f_3^\alpha$  |
| `n`       |`int`    | nombre d'itération pour la méthode|


```c
double approximation_succesive(double (*fonction_initial)(double), double x,double alpha, int n);
```
Cette fonction renvoit la valeur de $f_3^\alpha$ au point x, en utilisant la méthode d'approximations succesives.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
|`fonction_initial`| `fonction de double` | fonction initial de la méthode d'apporximations succesives|
|`x`        |`double` | x entre [0,L] de $f_3^\alpha (x)$              |
| `alpha`   |`double` | $\alpha$ de $f_3^\alpha$  |
| `n`       |`int`    | nombre d'itération pour la méthode|

```c
double ker_sep( double x,double alpha );
```
Cette fonction renvoit la valeur de $f_3^\alpha$ au point x, en utilisant la méthode des noyaux séparés.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
|`x`        |`double` | x entre [0,L] de $f_3^\alpha (x)$  |
| `alpha`   |`double` | $\alpha$ de $f_3^\alpha$  |

#### err.h
Ce programme est celui dans lequel on retrouve les fonctions qui calcul les erreurs en norme $L^{2}$ des différentes méthodes.

```c
double L2_err_f3_adov1(double alpha, int n);
```
Cette fonction renvoie $\|\| f_3^\alpha - U_{ex}(.,H)\|\| _{L^{2}} $, avec la première version de l'approximation de la méthode d'Adomain.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
| `alpha`   |`double` | $\alpha$ de $f_3^\alpha$  |
| `n`       |`int`    | nombre d'itération pour la méthode d'Adomain| 

```c
double L2_err_f3_adov2(double alpha, int n);
```
Cette fonction renvoie $\|\| f_3^\alpha - U_{ex}(.,H)\|\| _{L^{2}} $, avec la seconde version de l'approximation de la méthode d'Adomain.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
| `alpha`   |`double` | $\alpha$ de $f_3^\alpha$  |
| `n`       |`int`    | nombre d'itération pour la méthode d'Adomain| 

```c
double L2_err_f3_approx_succesive(double (*f)(double),double alpha,int n);
```
Cette fonction renvoie $\|\| f_3^\alpha - U_{ex}(.,H)\|\| _{L^{2}} $, avec la méthode d'approximation succesive.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |
|`fonction_initial`| `fonction de double` | fonction initial de la méthode d'apporximation succesive.|          
| `alpha`   |`double` | $\alpha$ de $f_3^\alpha$  |
| `n`       |`int`    | nombre d'itération pour la méthode d'approximation succesive|

```c
double L2_err_ker_sep_f3(double alpha);
```
Cette fonction renvoie $\|\| f_3^\alpha - U_{ex}(.,H)\|\| _{L^{2}} $, avec la méthode des noyaux séparables.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |         
| `alpha`   |`double` | $\alpha$ de $f_3^\alpha$  |

#### plot.h
Dans ce fichier on retrouve toutes les fonctions nous permettant de générer nos fichier .dat , qui nous permettrons ensuite de créer les figures associées.
```c
double plot_err_alpha_adov1(double alpha_deb,double alpha_fin, double pas,int n);
```
Cette fonction renvoie le $\alpha$ pour lequel l'erreur est minime parmis les points testés par la première version de la méthode d'Adomain, et écris toutes les valeur de $\alpha$ et l'erreur $L^{2}$ associé dans le fichier "dat/err_alpha_adov1.dat".
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |         
| `alpha_deb`|`double` |  premier $\alpha$ testé |
| `alpha_fin`|`double` |  dernier $\alpha$ testé |
| `pas`      |`double` | pas sur $\alpha$ |
| `n`       |`int`    | nombre d'itération pour la méthode|

```c
double plot_err_alpha_adov2(double alpha_deb,double alpha_fin, double pas,int n);
```
Cette fonction renvoie le $\alpha$ pour lequel l'erreur est minime parmis les points testés par la première version de la méthode d'Adomain, et écris toutes les valeur de $\alpha$ et l'erreur $L^{2}$ associé dans le fichier "dat/err_alpha_adov2.dat".
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- |         
| `alpha_deb`|`double` |  premier $\alpha$ testé |
| `alpha_fin`|`double` |  dernier $\alpha$ testé |
| `pas`      |`double` | pas sur $\alpha$ |
| `n`       |`int`    | nombre d'itération pour la méthode|

```c
double plot_err_alpha_approx_succesive(double (*f)(double),double alpha_deb,double alpha_fin,double pas,int n);
```
Cette fonction renvoie le $\alpha$ pour lequel l'erreur est minime parmis les points testés par la méthode d'approximation succesive, et écris toutes les valeur de $\alpha$ et l'erreur $L^{2}$ associé dans le fichier "dat/err_alpha_approx_succesive.dat".
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- | 
|`fonction_initial`| `fonction de double` | fonction initial de la méthode d'apporximation succesive|
| `alpha_deb`|`double` |  premier $\alpha$ testé |
| `alpha_fin`|`double` |  dernier $\alpha$ testé |
| `pas`      |`double` | pas sur $\alpha$ |
| `n`       |`int`    | nombre d'itération pour la méthode|

```c
double plot_err_alpha_ker_sep(double alpha_deb,double alpha_fin,double pas);
```
Cette fonction renvoie le $\alpha$ pour lequel l'erreur est minime parmis les points testés par la méthode de noyaux séparés, et écris toutes les valeur de $\alpha$ et l'erreur $L^{2}$ associé dans le fichier "dat/err_alpha_ker_sep.dat".
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- | 
| `alpha_deb`|`double` |  premier $\alpha$ testé |
| `alpha_fin`|`double` |  dernier $\alpha$ testé |
| `pas`      |`double` | pas sur $\alpha$ |

```c
int plot_f3_adov1( int nx,int n,double alpha);
```
Cette fonction écris dans le fichier "dat/f3_adov1.dat" les valeur de x et le $f_3^\alpha (x) $ associées par la première version de la méthode d'Adomain.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- | 
|`nx`| `int` | nombre de points selon l'axe des abcisses|
|`n`| `int` | nombre d'itération pour la méthode|
| `alpha`   |`double` | $\alpha$ (on utilisera le $\alpha$ optimal trouvé par la fonction plot_err_alpha_adov1) |


```c
int plot_f3_adov2( int nx,int n,double alpha);
```
Cette fonction écris dans le fichier "dat/f3_adov2.dat" les valeur de x et le $f_3^\alpha (x) $ associées par la seconde version de la méthode d'Adomain.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- | 
|`nx`| `int` | nombre de points selon l'axe des abcisses|
|`n`| `int` | nombre d'itération pour la méthode|
| `alpha`   |`double` | $\alpha$ (on utilisera le $\alpha$ optimal trouvé par la fonction plot_err_alpha_adov2) |


```c
int plot_f3ex( int nx);
```
Cette fonction écris dans le fichier "dat/f3ex.dat" les valeur de x et le $U_{ex} (x,H)$ associées.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- | 
|`nx`| `int` | nombre de points selon l'axe des abcisses|

```c
int plot_f3_app_succ(double (*f)(double), int nx,int n,double alpha);
```
Cette fonction écris dans le fichier "dat/f3_app_succ.dat" les valeur de x et le $f_3^\alpha (x) $ associées par la méthode d'approximations succesives.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- | 
|`fonction_initial`| `fonction de double` | fonction initial de la méthode d'apporximation succesive|
|`nx`| `int` | nombre de points selon l'axe des abcisses|
|`n`| `int` | nombre d'itération pour la méthode|
| `alpha`   |`double` | $\alpha$ (on utilisera le $\alpha$ optimal trouvé par la fonction plot_err_alpha_approx_succesive) |

```c
int plot_f3_ker_sep(int nx,double alpha);
```
Cette fonction écrit dans le fichier "dat/f3_ker_sep.dat" les valeur de x et le $f_3^\alpha (x) $ associées par la méthode des noyaux séparables.
| Paramètre | Type   | Description                   |
| --------- | ------ | ----------------------------- | 
|`nx`| `int` | nombre de points selon l'axe des abcisses|
| `alpha`   |`double` | $\alpha$ (on utilisera le $\alpha$ optimal trouvé par la fonction plot_err_alpha_ker_sep) |


### Dossier "prog"
Dans ce dossier, nous y retrouverons tous les fichiers .c de notre projet, y compris le "main.c", ou se trouve les instructions de notre programme.

Les détails concernants ces fichiers et leurs fonctions se trouvent dans la rubrique concernant les headers juste au dessus.

Voici la liste des programmes .c que vous y trouverez :
- main.c
- fonctions.c
- Gauss_Legendre.c
- plot.c
- methodes.c
- err.c
- OperationsMatrice.c

### Dossier "dat"
Dans ce dossier vous trouverez tous les fichiers .dat de données qui nous servent à tracer les figures, vous y trouverez :
- err_alpha_adov1.dat              : fichier où se trouve les donnnées permettant de tracer l'erreur selon $\alpha$ pour notre première version de la méthode de Adomain.

- err_alpha_approx_succesive.dat   : fichier où se trouve les donnnées permettant de tracer l'erreur selon $\alpha$ pour la méthode d'approximation succesive.

- err_alpha_ker_sep.dat            : fichier où se trouve les donnnées permettant de tracer l'erreur selon $\alpha$ pour la méthode des noyaux séparables.

- f3_adov1.dat                     : fichier où se trouve les donnnées permettant de tracer $f_3^\alpha$, l'approximation calculée par la méthode d'Adomain.

- f3_app_succ.dat                  : fichier où se trouve les donnnées permettant de tracer $f_3^\alpha$, l'approximation calculée par la méthode des approximations successives.

- f3_ker_sep.dat                   : fichier où se trouve les donnnées permettant de tracer $f_3^\alpha$, l'approximation calculée par la méthode des noyaux séparables.

- f3ex.dat                         : fichier où se trouve les donnnées permettant de tracer $f_3$, la solution exacte.






### Dossier "images"
Dans ce dossier vous trouverez toutes les figures en .png généré par gnuplot en utilisant les fichiers .dat corresponants :
- f3_adov1.png
- f3_approx_succesive.png
- f3_ker_sep.png
- f3_comparaison.png
- f3_comparaison_ss_app_suc.png
- T_err_alpha_adov1.png
- T_err_alpha_approx_succesive.png
- T_err_alpha_ker_sep.png

### Dossier "obj"
Ce dossier contiendra tous les .o correspondant aux .c du dossier "prog".
### Dossier "gnuplot"
Dans ce dossier, vous trouverez les scripts gnuplots créer pour générer les images du dossier "images", à l'aide des données du dossiers "dat", vous y trouverez :
- f3_adov1.plt : trace $f_3^\alpha$ calculée par notre première version de la méthode de Adomain.

- f3_approx_succesive.plt : trace $f_3^\alpha$ calculée par la méthode des approximations succesives.

- f3_ker_sep.plt : trace $f_3^\alpha$ calculée par la méthode des noyaux séparables.

- f3_comparaison.plt : trace les $f_3^\alpha$ calculé par nos différentes méthodes, et aussi $U_{ex}(x,H)$ pour visualiser les différences.

- T_err_alpha_adomainv1.plt : trace l'erreur en norme $L^{2}$ de $\|f_3^\alpha -U_{ex}(x,H)\|$ selon $\alpha$, avec $f_3^\alpha$ calculé par notre première version de la méthode d'Adomain.

- T_err_alpha_approx_succesive.plt : trace l'erreur en norme $L^{2}$ de $\|f_3^\alpha -U_{ex}(x,H)/|$ selon $\alpha$, avec $f_3^\alpha$ calculé par notre première version de la méthode d'approximation succesive.

- T_err_alpha_ker_sep.plt : trace l'erreur en norme $L^{2}$ de $\|f_3^\alpha -U_{ex}(x,H)/|$ selon $\alpha$, avec $f_3^\alpha$ calculée par la méthode des noyaux séparables.


### Dossier "rapport"
Dans ce dossier vous trouverez, un dossier "images" contenant les images utilisées pour notre rapport, et deux fichiers :

- Rapport_Projet2_Bichon_Gontier.pdf : contenant le rapport de notre projet.

- Rapport_Projet2_Bichon_Gontier.tex : contenant le fichier .tex avec lequel nousa vons écris ce rapport.



## Fonctionnement
Vous vous trouver dans le dossier du projet: "Projet2_BICHON_GONTIER
En premier lieux, le programme est déja compilé, mais si vous voulez le recompiler complétement ou faire des modifications pour afficher tester tel ou tel fonction, il peut être bon de faire la commande :
```
make clean
```

Cela supprimera tous les fichiers .o et exe de la précédente compilation, ensuite vous pourrez recompiler le programe, avec les instructions qui se trouve dans le fichier "prog/main.c" (attention cependant, des lignes sont mises en commentaires pour éviter de les recompiler à chaque fois, décommenter la ligne dont vous avez besoin) avec la commande :
```
make
```
Puis vous pourrez lancer le programme par l'une des deux commande :
```
make run
```
ou 
```
./exe
```
Notez que les toutes les fonctions sont listées dans la partie structure de ce README, plus particulièrement dans la rubrique headers.



### Afficher avec script gnuplot

Si vous souhaitez afficher une courbe ou un schéma en particulier.

Nous vous invitons à regarder la fonction du fichier "headers/plot" qui vous convient, la décommenter ou ajouter une ligne l'appelant dans le fichier "prog/main.c".

Puis ensuite relancer la compilation et l'execution du programme, avec les commandes précédentes.

Maintenant, toujours dans le dossier "Projet2_BICHON_GONTIER"(important car les appels et les écriture des script gnuplot se font de ce dossier), faites la commande :
```
gnuplot
```
Cela vous ouvrira l'interface gnuplot dans le terminal.
Ensuite vous pourrez lancer la création du fichier .png correspondant à votre figure par la commande : (ici pour le script T_err_alpha_adomainv1.plt qui affiche l'erreur en norme $L^{2}$ entre $f_3^\alpha$ et la solution exacte en H, $U_{ex}(x,H)$ )
```
load "gnuplot/T_err_alpha_adomainv1.plt"
```

Ensuite vous pourrez ouvrir l'image correspondante dans le dossier "images", vous pourrez de la manière que vous voulez, par exemple depuis depuis le terminal dans le dossier principal du projet "Projet2_BICHON_GONTIER" par la commande :
```
eog images/T_err_alpha_adomainv1.png
```
si vous avez eog sur votre machine.
