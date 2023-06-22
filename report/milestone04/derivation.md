# Derivation of forces for Lennard-Jones potential
The force acting upon a single atom $k$ is given by:
$$\vec{f_k} = \sum_{j} \vec{f_{kj}} 
= \sum_{j} \frac{\partial{V(r_{kj})}}{\partial{r_{kj}}} \cdot \frac{\vec{r_{kj}}}{r_{kj}}$$
Where $r_{kj} = |\vec{r_k} - \vec{r_j}|$ is the Euclidian distance between positions of atom $k$ and atom $j$, and $\vec{r_{kj}}$ is the vector pointing from atom $k$ to atom $j$, i.e. $\vec{r_{kj}} = \vec{r_k} - \vec{r_j}$.
$$\frac{\partial{V(r_{kj})}}{\partial{r_{kj}}} 
= \frac{\partial}{\partial{r_{kj}}} 4 \epsilon \left[ \frac{\sigma^{12}}{r_{kj}^{12}} - \frac{\sigma^6}{r_{kj}^6} \right]
= 4 \epsilon \left[ \frac{\partial}{\partial{r_{kj}}} \frac{\sigma^{12}}{r_{kj}^{12}} - \frac{\partial}{\partial{r_{kj}}} \frac{\sigma^6}{r_{kj}^6} \right]$$
$$= 4 \epsilon \left[ 12 \frac{\sigma^{12}}{r_{kj}^{13}} - 6 \frac{\sigma^6}{r_{kj}^7} \right]$$
So finally the forces acting upon atom $k$ are:
$$\vec{f_k} = \sum_{j} 4 \epsilon \left[ 12 \frac{\sigma^{12}}{r_{kj}^{13}} - 6 \frac{\sigma^6}{r_{kj}^7} \right] \cdot \frac{\vec{r_{kj}}}{r_{kj}}$$

**Note: some confusion with signs and $r_{kj}$ bzw. $r_{jk}$, not sure if what I have is right b ut if I change signs in the code the tests fail**