Let  ${\Theta}= (\Theta_1, \Theta_2, \dots, \Theta_p)$ be a vector of random angles, it has a multivariate von Mises distribution; i.e. such that:

$f_{{\Theta}}(\Theta)= C^{-1}_p({\kappa},{\Lambda})\exp\{{\kappa}^Tc({\Theta},{\mu})+\frac{1}{2}s({\Theta},{\mu})^T{\Lambda}s({\Theta},{\mu})\}$

where $-\pi<\theta_j\leq\pi$  , $-\pi<\mu_j\leq\pi$ , $\kappa_j \geq 0$ ,  $-\infty<\lambda_{jl}<\infty$ , and  $C^{-1}_p({\kappa},{\Lambda})$ is a normalizing constant, 
$c({\Theta}, {\mu})^T=(\cos(\theta_1-\mu_1),\cos(\theta_2-\mu_2),\dots,\cos(\theta_p-\mu_p)),$ 

$s({\Theta},{\mu})^T=(\sin(\theta_1- \mu_1),\sin(\theta_2- \mu_2),\dots, \sin(\theta_p - \mu_p) ),$

${\mu}^T=(\mu_1, \mu_2, \dots, \mu_p)\\qquad\qquad{\kappa}^T=(\kappa_1,\kappa_2,\dots,\kappa_p),$

$\lambda_{jl}=\lambda_{lj}\,\qquad\qquad\lambda_{jj}=0$

It is possible to build a graphical model, $(G,Q)$, where $G=(V,E)$, and the set $V$ contains $p$ vertices, each of this is referred about each variables $\Theta=(\Theta_1,\Theta_2,\dots,\Theta_p)$, and the set $E$ is such that a missing edge indicates a conditional independence among variables.


The literature about this kind of graphical models is really poor and it focuses the attention on protein structure problem. Some examples are: \cite{boomsma2008generative} that use a dynamic Bayesian network, specifically, they adapt a generalization of hidden Markov model, to localize the protein structure in naive status; \cite{lennox2009density} that present a Bayesian approach to density estimation for bivariate data which combines a Dirichlet process mixture model and a bivariate von Mises.\\

Most important work for our aim is one of \cite{razavian2011mises}. They would predict the three dimensional structure of a protein given its amino acid sequence, and adopting an undirected graphical model.\\ They use a factor graph for study the independence structure, furthermore, they assume that angles of the amino acid sequence have a von Mises distribution.\\ The biggest problem of this approach is that the likelihood, and consequentially the log-likelihood, are not in closed form. To overcome this problem, they use a pseudo-likelihood to parameter estimate.

\section{Introduction to protein structure}
In structural biology there is an unsolved problem called \emph{protein folding problem}. In short, this problem is really important, in fact,since 1994, there is a biennial global competition called the Community Wide Experiment on the Critical Assessment of Techniques for Protein Structure Prediction (CASP). This competition has become the gold standard for assessing techniques, and the last year it has been won by a department of Google, Deep Mind, which has presented an Artificial Intelligence called AlphaFold. (\cite{AlphaFold})\\
This is to predict the three dimensional structure of a protein given its amino acid sequence.\\


Proteins are large and complex molecules, they are essential in sustaining life. In this paragraph, we give a brief outline of the fundamental aspects of their structure. For furthermore information can be found, for example in \cite{branden1999introduction, lesk2001introduction}.\\
\begin{figure}
	\centering
	\includegraphics[scale=0.6]{README_files/figure-html/pep}
	\caption{Peptide bond arising from the carboxyl group of the amino acid 1 condenses with anime group of amino acid 2. This drops a water molecule ($H_2O$) and form a peptide bond between the carbon atom of first amino acid and the nitrogen atom of the second amino acid. The process continues to produce a polypeptide.  }
	\label{fig:pep}
\end{figure}
Specific atomic groups, so-called \emph{amino acids}, give rise to a protein. There are twenty commonly occurring amino acid, each with the same structure. In particular, an amino acid has a central carbon acid, denoting by $C^\alpha$, which bind to a hydrogen atom, $H$, an amine group $NH_2$ , a carboxyl group ($COOH$) and a side chain. The last identify the amino acid kind, while the carboxyl group makes a link among different amino acids, and is called peptide bond. Figure \ref{fig:pep} shows an example of this bond. One or more peptide bond form a polypeptide chains, and one or more of those form a protein.

The amino acid sequence is called primary structure, the peptide chains are named secondary structure, instead, polypeptide chains are defined as tertiary structures. 

\paragraph{The protein backbone}
\cite{moss1996basic} gives the following definition of backbone :\\
\begin{center}
	"\emph{ That linear chain to which all other chains, long or short or both, may be regarded as being pendant. Note: where two or more chains could equally be considered to be the main chain, that one is selected which leads to the simplest representation of the molecule. }" 
\end{center}
For comprehensive accounts on protein structure see, for example, \cite{lesk2001introduction}.\\

 An example of protein backbone is:
$$ N_1 - C^\alpha_1 - C_1 - N_2 - C^\alpha_2 - C_2 - \cdots - N_p - C^\alpha_p -C_p$$
 \begin{figure}
	\centering
	\includegraphics[scale=0.2]{img/dihedral1}  $\qquad$   \includegraphics[scale=0.15]{img/treangoli}
	\caption{Left panel shows a dihedral angle $\theta$ defined in terms of atoms $A_{i-2}, A_{i-1}, A_i , A_{i+1}$. Note that $A_{i-1}$ and $A_{i+1}$ are coplanar, it is the same for $A_{i}$ and $A_{i+1}$. The right panel shows a dihedral angles $\psi_i , \phi_i$, and $\omega_i$ in terms of atoms. }
	\label{fig:dihedral1}
\end{figure}
We consider a dihedral angle, which is the angle between planes through two sets of three atoms, having two atoms in common. 
\begin{figure}
	\centering
	\includegraphics[scale=0.5]{img/dihedral2}
	\caption{Dihedral angles $\phi$, $\psi$, and $\omega$.}
	\label{fig:dihedral2}
\end{figure}

Figure \ref{fig:dihedral1} shows a dihedral angle $\theta$ defined in terms of atoms $A_{i-2}, A_{i-1}, A_i , A_{i+1}$. It is important to know that the plane $\Pi_1$ passing through the atoms $A_i$ and $A_{i-2}$, as the plane $\Pi_2$ passing through the atoms $A_{i-1}$ and $A_{i+1}$. The angle $\theta$ is measured between $-\pi$ and $\pi$, the zero direction is observed when the plane $\Pi_1$ and $\Pi_2$ collapses. Clearly, if $A_{i}=N_i$, then $A_{i+1}=C_{i}$, $A_{i-1}=C^\alpha_i$ and $A_{i-2}=C_i$, while $\theta$ is $\phi_i$. This situation is clarifying in figure \ref{fig:dihedral2}.\\

 However, most combinations of $\phi$ and $\psi$ angles are impossible because they would result in steric collisions between backbone and side chain atoms. Hence, it is useful to look at the Ramachandran plot. (\cite{ramachandran1968conformation}) \\
 
 This graph plots $\psi$ versus $\phi$ angles for amino acid proteins, in this way the combination observables come out, and they can be classified as alpha helices and beta strands. The pairs ($\phi$, $\psi$)  commonly tends to make an asymmetric clusters: if $\psi$ angle is positive and $\phi$ is negative, then the pairs contributes to a beta stand shape, while if both of them are negative, they contributes to a right-handed alpha helices, and if both of them are positive, we can look a left-handed alpha helices.
 
 
 Figure \ref{fig:Rama} clarifies how to read a Ramachandran plot, and highlight each zone that corresponds to the different shape of the secondary structure.\\
 
 \subsection{Protein data}
\begin{figure}
	\centering
	\includegraphics[scale=0.5]{img/prova}
	\caption{Ramachandran plot, that draws a $\psi$ angle versus $\phi$ angle.}
	\label{fig:Rama}
\end{figure}

Protein structure can be surveyed  through an atomic level by X-ray diffraction and neutron-diffraction studies of crystallized proteins, and more recently through nuclear magnetic resonance (NMR) spectroscopy of proteins in solution. In this way, we observe the $N-C_\alpha$ and $C_\alpha - C$ bonds that are free rotation, and are represented by the torsion angle $\phi$ and $\psi$, respectively.\\

Fortunately, there is a Protein Structure Databank (PDB) (\cite{berman2000protein}), that collects over ten thousand angle pairs, which provide ample data for even relatively studies.  For our purpose, we choose an existing protein codifying as 1yt6, which is unclassified protein. This protein is cited for the first time in \cite{murata2005structure}.

Nine amino acids compose 1yt6, and its structure is characterized by only a beta sheet.  
\begin{figure}
	\centering
	\includegraphics[scale=0.3]{img/scorpion}
	\caption{Three dimensional shape of 1yt6 protein.}
	\label{fig:scorpion}
\end{figure}

In Figure \ref{fig:scorpion} we can see the three dimensional shape of protein, that is given by PDB file, this also contains a matrix data, with fifty rows, each corresponding to an observation, and nine columns, one for each amino acid of the sequence.


```r
summary(cars)
```

```
##      speed           dist       
##  Min.   : 4.0   Min.   :  2.00  
##  1st Qu.:12.0   1st Qu.: 26.00  
##  Median :15.0   Median : 36.00  
##  Mean   :15.4   Mean   : 42.98  
##  3rd Qu.:19.0   3rd Qu.: 56.00  
##  Max.   :25.0   Max.   :120.00
```

## Including Plots

You can also embed plots, for example:

<img src="README_files/figure-html/pressure-1.png" width="672" />

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
