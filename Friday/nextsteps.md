#  What to do next?

## Readings on advanced sampling methods
 - The big review: [Enhanced Sampling Methods for Molecular Dynamics Simulations [Article v1.0]](https://livecomsjournal.org/index.php/livecoms/article/view/v4i1e1583)
   - This review is part of LiveCoMS, and attemps to systemetize the differences between different methods.  It is rather extensive, with over 400 references so far!  -    - Because it is a "live" article, you can file issues on the GitHub repository associated with the article asking to have things explained better, or have and additional method we forgot to add!
 - Other reviews:
   -
## Useful Software Tools

- [PLUMED](https://www.plumed.org/)
  - PLUMED is the metadynamics code that was discussed earlier.  Note that is bigger than just metadynamics - there are other related methods that are supported as well. 
- [PLUMED-Nest](https://www.plumed-nest.org/) 
  - A repository of "recipies" written for PLUMED; it's possible that something you want to try to do, or something SIMILAR to what you wanted to do, is already implemented, and you can borrow that!
- [WESTPA](https://westpa.readthedocs.io/)
  - Software for a method we didn't really talk about, _weighted ensemble_ which biases sampling by maintaining "copies" of the system equally populated along the collective variable of interest.  
  - A LiveCoMS article describing some WESTPA tutorials [linked here](https://livecomsjournal.org/index.php/livecoms/article/view/v1i2e10607). 
- [pymbar](https://pymbar.readthedocs.io/)
  - A set of functions for calculating free energies between states and removing biases using multistate reweighting. [You can see a video about multistate reweighting here](https://www.youtube.com/watch?v=yGyQa8opfi0). 
- [alchemlyb](https://alchemlyb.readthedocs.io/)
  - a wrapper library to make it easy to parse output from molecular simulation programs into pymbar and other free energy calculation routines, as well as interact usefully with the analysis coming out of this software.
- [physical_validation](https://physical-validation.readthedocs.io/)
  - A set of tools to determine if your simulation obeys certain physical laws, such as equipartition and consistency with the Boltzamm distribution. You can see a [video about the approaches to do this here](https://www.youtube.com/watch?v=-Zxvi7EQwE4).

## Cutting edge questions

- How to find better collective variables?
