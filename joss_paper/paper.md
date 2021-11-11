---
title: 'sunraster: Handling Solar Scanning Slit Spectrograph Observations in Python'
tags:
  - Python
  - astronomy
  - UV
  - solar physics
  - spectroscopy
authors:
  - name: Daniel F. Ryan^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0001-8661-3825
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Author Without ORCID^[co-first author] # note this makes a footnote saying 'co-first author'
    affiliation: 2
  - name: Author with no affiliation^[corresponding author]
    affiliation: 3
affiliations:
 - name: University of Applied Sciences Northwest Switzerland
   index: 1
 - name: Institution Name
   index: 2
 - name: Independent Researcher
   index: 3
date: 12 August 2021
bibliography: paper.bib

---

# Summary

sunraster is a free, open-source, community-developed Python package for inspecting,
manipulating, and visualizing solar observations from slit-spectrograph instruments.
It builds upon the data classes of the ndcube package to provide slicing, coordinate
transformation, and visualization APIs.
These APIs come in two flavours ("sns" and "raster") that allow scientists to treat
the data as if they were 3-D (time-space-wavelength) or 4-D (time-space-space-wavelength).
Due to the design of slit-spectrographs both these models can be simultaneously valid
representations of the data.
(See the Introduction to Slit-spectrographs section for more detail.)
Providing both these APIs on the same data objects increase the flexibility of
scientists in their analysis while reducing the complexity of their analysis workflows.
This enhances scientific output and reduces the likelihood of mistakes.
sunraster also provides data readers for commonly used solar slit-spectrograph instruments.
The current version at time of writing (v0.3) provides readers for Solar Orbiter/SPICE
and IRIS.

# Introduction to Slit-spectrographs

Slit-spectrographs disperse light passing through a narrow slit perpendicular to the
slit's long axis.
A single "image" produced by a slit-spectrograph therefore has space along the long slit
axis and wavelength along the short one.
Such an image is called a spectrogram.
When spectrograms are stacked chronologically, we can extend the data into a 3-D
time-space-wavelength spectrogram cube.
While Python tools already exist for handling N-D astronomical data, namely ndcube,
the nature of slit-spectrographs make their data unusually complex and require more
specialized data classes.
To overcome their very narrow field of view in the x-direction, slit-spectrographs
often scan their slits across a region of the Sun to build up a spectral image.
This is known as rastering.
Raster scans are typically repeated numerous times during an observing campaign,
sometimes at rapid cadence.
This results in 4-D data cubes with axes corresponding to
the number of the raster scan, the slit position within the scan,
the position along the slit, and wavelength.
These axes corresponds to physical types of time-time/space-space-wavelength
where time and space are folded together in the x-axis due to the fact that
spectrograms at each position in the scan are taken sequentially in time.
Depending on the scientist's use-case, the data may need to be represented as a
3-D time-space-wavelength cube or a 4-D time-space-space-wavelength cube.
Moreover scientists may need to switch back and forth between these paradigms at
different stages of their analysis.

# Statement of Need

The traditional solution to providing both the 3-D and 4-D representations of
slit-spectrograph observations is to create two copies of the data with a 3-D and
4-D shape.
Scientists can then use whichever version is most appropriate at the time.
However this approach is inefficient for the computer memory and unwieldly for the user,
as the two versions must be kept consistent.
It adds complexity to the analysis workflow and increases the chance of mistakes.
By contrast sunraster provides data classes with equivalent sets
of manipulation and visualization APIs for the 3-D and 4-D pardigms.
The tools underlying these APIs translate between them and manipulate the data in
a self-consistent way, thereby removing the need to duplicate the data.
This enables users to seemlessly and instantly transition between the 3- and 4-D paradigms.
This makes analysis simpler, clearer, and less prone to mistakes.

Solar slit-spectrographs, like Hinode/EIS and IRIS, have long provided valuable
spectral measurements of the solar atmosphere and helped us better understand many
phenomena including solar flares, jets, and active region heating, to name but a few.
Proposed and newly launched instruments, such as Solar Orbiter/SPICE,
continue to be designed on the slit-spectrograph principle.
Therefore the tools provided by sunraster are important for current and future solar
physics research.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
