#  What to do next?

## Readings on advanced sampling methods
 - The big review: [Enhanced Sampling Methods for Molecular Dynamics Simulations [Article v1.0]](https://livecomsjournal.org/index.php/livecoms/article/view/v4i1e1583)
   - This review is part of LiveCoMS, and attemps to systemetize the differences between different methods.  It is rather extensive, with over 400 references so far!  -    - Because it is a "live" article, you can file issues on the GitHub repository associated with the article asking to have things explained better, or have and additional method we forgot to add!
- Other reviews about enhanced sampling: 
  - [Enhanced sampling without borders: on global biasing functions and how to reweight them](https://doi.org/10.1039/d1cp04809k)
  - [Enhancing Important Fluctuations: Rare Events and Metadynamics from a Conceptual Viewpoint](https://doi.org/10.1146/annurev-physchem-040215-112229)
  - [Using metadynamics to explore complex free-energy landscapes](https://doi.org/10.1038/s42254-020-0153-0)
  - [Metadynamics: A Unified Framework for Accelerating Rare Events and Sampling Thermodynamics and Kinetics](https://doi.org/10.1007/978-3-319-44677-6_49)
  - [Analyzing and Biasing Simulations with PLUMED](https://doi.org/10.1007/978-1-4939-9608-7_21) (also https://arxiv.org/abs/1812.08213)
  - [Best Practices for Alchemical Free Energy Calculations v1.0](https://livecomsjournal.org/index.php/livecoms/article/view/v2i1e18378) 

   
## Useful Software Tools

- [PLUMED](https://www.plumed.org/)
  - PLUMED is the metadynamics code that was discussed earlier.  Note that is bigger than just metadynamics - there are other related methods that are supported as well. The official manual for PLUMED version 2.9 can be found [here](https://www.plumed.org/doc-v2.9/user-doc/html/index.html). A number of tutorials can also be found in the [manual](https://www.plumed.org/doc-v2.9/user-doc/html/tutorials.html)
  - [PLUMED-Nest](https://www.plumed-nest.org/) 
    - A repository of "recipies" (input file sets) written for PLUMED; it's possible that something you want to try to do, or something SIMILAR to what you want to do, is already implemented, and you can borrow that!
   - [PLUMED Masterclass](https://www.plumed.org/masterclass): Corresponding YouTube videos with lectures and solutions found [here](https://www.youtube.com/watch?v=2eGhMSdIJEs&list=PLmdKEn2znJEld8l6Hp9PXf4EursC4-8nC)
    - [Masterclass 21.1: PLUMED syntax and analysis](https://www.plumed.org/doc-v2.9/user-doc/html/masterclass-21-1.html)
    - [Masterclass 21.2: Statistical errors in MD](https://www.plumed.org/doc-v2.9/user-doc/html/masterclass-21-2.html)
    - [Masterclass 21.3: Umbrella sampling](https://www.plumed.org/doc-v2.9/user-doc/html/masterclass-21-3.html)
    - [Masterclass 21.4: Metadynamics](https://www.plumed.org/doc-v2.9/user-doc/html/masterclass-21-4.html)
    - [Masterclass 21.5: Simulations with multiple replicas](https://www.plumed.org/doc-v2.9/user-doc/html/masterclass-21-5.html)
    - [Masterclass 21.6: Dimensionality reduction](https://www.plumed.org/doc-v2.9/user-doc/html/masterclass-21-6.html)
    - [Masterclass 21.7: Optimizing PLUMED performances](https://www.plumed.org/doc-v2.9/user-doc/html/masterclass-21-7.html)
    - [Masterclass 22.3: OPES method](https://www.plumed.org/doc-v2.9/user-doc/html/masterclass-22-03.html)
    - [Masterclass 22.11: Variationally enhanced sampling with PLUMED](https://www.plumed.org/doc-v2.9/user-doc/html/masterclass-22-11.html)
- [pymbar](https://pymbar.readthedocs.io/)
  - A set of functions for calculating free energies between states and removing biases using multistate reweighting. [You can see a video about multistate reweighting here](https://www.youtube.com/watch?v=yGyQa8opfi0). 
- [alchemlyb](https://alchemlyb.readthedocs.io/)
  - a wrapper library to make it easy to parse output from molecular simulation programs into pymbar and other free energy calculation routines, as well as interact usefully with the analysis coming out of this software.
- [physical_validation](https://physical-validation.readthedocs.io/)
  - A set of tools to determine if your simulation obeys certain physical laws, such as equipartition and consistency with the Boltzamm distribution. You can see a [video about the approaches to do this here](https://www.youtube.com/watch?v=-Zxvi7EQwE4).
- [WESTPA](https://westpa.readthedocs.io/)
  - Software for a method we didn't really talk about, _weighted ensemble_ which biases sampling by maintaining "copies" of the system equally populated along the collective variable of interest.  
  - A LiveCoMS article describing some WESTPA tutorials [linked here](https://livecomsjournal.org/index.php/livecoms/article/view/v1i2e10607). 

## Cutting edge questions

How do I find good collective variabies?

### Reviews and other useful papers about collective variables
- [Reaction coordinates in complex systems-a perspective](https://doi.org/10.1140/epjb/s10051-021-00233-5)
- [Collective Variables for Crystallization Simulations─from Early Developments to Recent Advances](https://doi.org/10.1021/acsomega.2c06310)
- [Critical comparison of general-purpose collective variables for crystal nucleation](https://doi.org/10.1103/physreve.107.l012601)

### Review papers about machine learning for finding collective variables
- [Manifold Learning in Atomistic Simulations: A Conceptual Review](https://doi.org/10.48550/arXiv.2303.08486)
- [Machine Learning Force Fields and Coarse-Grained Variables in Molecular Dynamics: Application to Materials and Biological Systems
](https://doi.org/10.1021/acs.jctc.0c00355)
- [Machine learning for collective variable discovery and enhanced sampling in biomolecular simulation](https://doi.org/10.1080/00268976.2020.1737742)
- [Chasing collective variables using temporal data-driven strategies](https://doi.org/10.1017/qrd.2022.23)
- [Machine learning approaches for analyzing and enhancing molecular dynamics simulations](https://doi.org/10.1016/j.sbi.2019.12.016)
  



