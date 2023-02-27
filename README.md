# AGN-in-Mergers

## Background

### Merging galaxies

Galaxies are large gravitationally bound systems of stars, gas, dust, and dark matter. Contrary to first impressions, these objects are not stationary. Galaxies also tend not to exist in isolation, but in clusters of galaxies. Within these clusters, galaxies frequently interact and merge like the example below. 

Merging galaxies are not only beautiful, but they are an important part of galaxy evolution. Most massive galaxies, including our galaxy, the Milky Way, are a product of several past galaxy mergers. In fact, the Milky Way will collide with the Andromeda galaxies in about 4 billion years. 

![Figure 1](https://github.com/jlsteffen/AGN-in-Mergers/blob/main/images/merger_ex.jpg)

### Active Galactic Nuclei (AGN)

<img align="left" width="400" height="350" src="https://github.com/jlsteffen/AGN-in-Mergers/blob/main/images/AGN_NASA.jpg"> 

When large amounts of material fall onto a supermassive black hole, the material does not fall directly into black hole. Imagine water draining from a sink. The water slowly drains through a whirlpool which forms around the opening of the sink. Material falls through a black hole in a fundamentally similar way, though the extreme scales involved create some differences. The gases and dust in the whirlpool, which we call the accretion disk, have extreme velocities and the friction between the gases cause intense amounts of heat. The heated accretion disk then produces large amounts of high energy photons. The black hole, an entity so massive that light cannot escape it, ironically becomes one of the brightest objects in the universe (or at least the accretion disk become the brightest object). We refer to these objects as active galactic nuclei (AGN). 

This phenomena only occurs when a large amount of material falls onto the black hole. Where does this material come from? One possible source of this fuel could be from galaxy mergers. When two galaxies interact, gases in the the disks of the galaxies will lose angular momentum. As a consequence, these gases will fall inward to smaller orbits. If some of these gases loose enough angular momentum, they may fall into the very center of the galaxy where they can fuel the central black hole. 

## The Goal
Do merging galaxies induce AGN activity? To do this we will take a sample of galaxies and split it into a merging galaxies and isolated galaxies. Then we will identify AGN within both samples. If we find that the merging galaxies host a higher rate of AGN than the isolated sample, we can conclude that galaxy mergers can induce AGN activity. 

## The Survey

<img align="left" width="400" height="350" src="https://github.com/jlsteffen/AGN-in-Mergers/blob/main/images/manga_v3.jpg"> 

This work is based on the data from the SDSS's (Sloan Digital Sky Survey) MaNGA (Mapping Nearby Galaxies at Apache Point Observatory) survey. MaNGA is an integral field spectroscopic (IFS) survey which has observed 10,000 nearby galaxies. Spectroscopy takes the light from an object and splits it by wavelength. We can tell many things about an object by its spectra. You can tell what king of object you are looking at, the ages of stars, the chemical composition, et cetera. 

Traditionally, for large spectroscopic surveys, light is collected through fiber optic cables to prevent spectra from overlapping with other light sources. A single fiber optic cable is placed on a single source, and a single spectrum is collected. IFS takes a whole bundle of fiber optic cables and places them onto a single target. This gives astronomers several spectra for a single target, each of which covers a different part of the object. With MaNGA, IFS allows us to simultaneously study a galaxy's center and disk with a single observation. 

## Identifying Merging Galaxies

While we can easily tell that some galaxies are undergoing a merger event due to the presence of clear tidal features, like tidal tails or warped disks, many merging galaxies are hard to spot by visual inspection alone. Many galaxies may be gravitationally paired with one another but display no obvious tidal features in their images. We need to come up with a robust method for identifying paired galaxies. 

The first step is to come up with a catalog of galaxies to inspect. The MaNGA survey sets a good starting point, having observed 10,000 galaxies, but this sample can be improved. Each of the MaNGA survey's 10,000 observations target a single galaxy; however, in many observations there will be other galaxies, stars, et cetera that also fall within the observation. Some of these other galaxies may be gravitationally bound with the target galaxy. To obtain these galaxies, we need to build a list of all object within the MaNGA survey's fields-of-view. I have constructed this catalog in a separate project (https://github.com/jlsteffen/MaNGAObj). This catalog contains 15,000 objects which have been classified (star, galaxy , quasar, etc.) based on their spectra. 

We can now build a sample of paired galaxies by selecting galaxies which are physically near each other. We do this by selecting galaxies that are nearby on the sky and which have similar redshifts and line-of-sight velocities (which tells us that the objects are nearby along our line-of-sight). Finally, since we want to tell how paired galaxies are different from other galaxies, we create a control sample of galaxies which have no nearby companion. In total we have 391 galaxy pairs and 7811 control galaxies.

## Identifying Active Galactic Nuclei

The next step is to determine which galaxies are hosting an active galactic nucleus. This can be done by studying the spectrum from a galaxy's center. The high energy photons generated by the AGN will ionize nearby gas cloud which then emit distinct emission lines which we can detect in their spectrum. Below is the spectrum from the center of an ordinary galaxy (with active star formation).

image.png

Now, here is another spectrum but for a galaxy hosting an AGN.

image.png

One can see the strong emission lines at xxx nm and xxx nm which are the signature tell of an AGN. As an objective way to identify AGN, we take the emission line ratios [O III]/H$\beta$ and [N II]/H$\alpha$ and plot them on the below diagram. 

image.png

Each data point represents the emission line ratios from a single galaxy. The dashed line separates galaxies whose emission lines are purely from star formation (young massive stars; below the line) and galaxies whose emission lines are from AGN (above the line). Galaxies between the dashed line and the solid line have contributions from star formation and AGN. Galaxies above the solid line have emission lines that are dominated by AGN activity. The dashed-dot line separates Seyfert galaxies (strong AGN) and LINERS (weak AGN). For our work, we will use galaxies whose emission line ratios are above the dashed line. 

Finally, we are also aware that some old but massive stars can create AGN-like emission line ratios. We know that the equivalent width of the H$\alpha$ line is a good tracer for these objects so we decide that objects in our AGN catalog need to have a H$\alpha$ equivalent width above 6$\AA$, as is shown in the below figure.

In our pair sample, we find 
