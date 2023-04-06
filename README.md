# AGN-in-Mergers

## Background

<img align="left" width="500" height="350" src="https://github.com/jlsteffen/AGN-in-Mergers/blob/main/images/merger_ex.jpg"> 

### Merging galaxies

Galaxies are large gravitationally bound systems of stars, gas, dust, and dark matter. Contrary to first impressions, these objects are not stationary. Galaxies also tend not to exist in isolation, but in clusters of galaxies. Within these clusters, galaxies frequently interact and merge like the example on the left. 

Merging galaxies are not only beautiful, but they are an important part of galaxy evolution. Most massive galaxies, including our galaxy, the Milky Way, are a product of several past galaxy mergers. In fact, the Milky Way will collide with the Andromeda galaxies in about 4 billion years. 

#

<img align="left" width="500" height="350" src="https://github.com/jlsteffen/AGN-in-Mergers/blob/main/images/AGN_NASA.jpg"> 

### Active Galactic Nuclei (AGN)

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

The next step is to determine which galaxies are hosting an active galactic nucleus. This can be done by studying the spectrum from a galaxy's center. The high energy photons generated by the AGN will ionize nearby gas cloud which then emit distinct emission lines which we can detect in their spectrum. Below is the spectrum from a galaxy.

![Figure 1](https://github.com/jlsteffen/AGN-in-Mergers/blob/main/images/7959-12704.png)

The emission lines around 5000 angstroms and 6500 angstroms are created by ionized gases. The typical sources of ionization in a galaxy is through star formation (i.e. young massive stars) or AGN. An AGN will produce a strong power-law spectrum which extend far into the far-UV in comparison to the spectrum of OB stars. As a consequence, AGN will produce stronger low-ionization forbidden lines relative to the H alpha line in comparison to OB stars. Thus, the ionization sources for these emission lines can be revealed by studying the ratios between forbidden lines and the Balmer series. 

![Figure 2](https://github.com/jlsteffen/AGN-in-Mergers/blob/main/images/8341-12704.png)

The panel on the right shows what is called a BPT diagram. Each data point represents an emission line ratio from a points on the galaxy on the left panel. The dashed line separates galaxies whose emission lines are purely from star formation (young massive stars; below the line) and galaxies whose emission lines are from AGN (above the line). Galaxies between the dashed line and the solid line have contributions from star formation and AGN. Galaxies above the solid line have emission lines that are dominated by AGN activity. The dashed-dot line separates Seyfert galaxies (strong AGN) and LINERS (weak AGN). For our work, we will use galaxies whose emission line ratios are above the dashed line. The middle panel shows the image of the galaxy with the BPT classifications overlaid. The region ionized by an AGN is centrally concentrated and the outlying regions are dominated by star formation. 

I perform this analysis for each galaxy in our sample; however, only for the centers of the galaxies since that is where an AGN would be found. I find 105 AGN in the pair sample and 872 AGN in the control sample. 

## AGN Volume Density

Now answering our question seems simple. If the paired galaxies have a higher fraction of AGN than the control galaxies, then we can say that galaxy merging can enhance the rate of AGN activity? The issue is that AGN duty cycle in MaNGA is mass dependent and peaks around log(M/Msun) = 10.5 and the pair fraction in MaNGA prefers higher mass galaxies. Further, the MaNGA sample itself has a flat mass distribution so high mass galaxies are being oversampled and low mass galaxies are being undersampled. 

To resolve this, I will model the mass dependent AGN fraction in the control galaxies and use this model to estimate expected volume density of AGN in the pair sample assuming stochastic fueling alone. If the measured volume density of AGN is higher than the expected stochastic volume density, we can assume that the excess volume density is from merger induced fueling.

Below is the distribution of stellar masses for the whole control sample (grey) and the control sample hosting AGN (blue). The AGN fraction is simply the AGN sample divided by the whole sample, shown in the black outlined diamonds. 

![Figure 3](https://github.com/jlsteffen/AGN-in-Mergers/blob/main/images/20agn_hist.png)

I then model the mass dependent AGN fraction with a Normal distribution in the below formula. 

![Figure 4](https://github.com/jlsteffen/AGN-in-Mergers/blob/main/images/Equation7.png)

Here, M is the stellar mass of a galaxy, z is the galaxy's redshift, and f_0, sigma, and b are amplitude, variance, and mean of the Normal distribution that I fit for. I fit these values with a minimization routine and find that f_0 = 0.12, sigma = 0.44, and b = 10.56. The AGN fraction here will have no significant evolution in redshift within the range included in the MaNGA sample, but I include it here for completeness since we know that AGN fraction is redshift dependent across broader redshift ranges. 

This model represents the AGN induced through the stochastic, or random fueling of the AGN. The AGN in the paired galaxies will likely still be subject to the stochastic means of AGN fueling. If the merger process does not influence AGN activity, we should be able to predict the volume density of AGN in the pair sample using the model above. If we find that the observed volume density of AGN is higher than what is predicted by stochastic fueling, we may assume that some of the observed AGN are being induced by the galaxy merger. 

I break the AGN in pairs sample into three sub-groups; offset AGN in which one galaxy hosts an AGN, dual AGN in which both galaxies host an AGN, and the combined offset+dual AGN sample. Below I show the equations for the volume densities for the dual, offset, and offset+dual AGN samples respectively. W_j represents the volume weights from the MaNGA sample and f^t and f^c are the stochastic AGN probabilities for the MaNGA target galaxy and its companion respectively. These equations may initially look intimidating, they are similar in form to the probability of flipping two coins. The first equation is the probability of getting two heads up, the second equation is the probability of getting a single heads up in either coin, and the last equation is the probability of getting any heads up. 

![Figure 5](https://github.com/jlsteffen/AGN-in-Mergers/blob/main/images/Equation1012)

Below I plot the AGN volume density for my pair sample as a function of the pair separation. The grey scatter points are the volume weights for individual paired galaxy while the blue scatter points are the volume weights for paired galaxies hosting an AGN. The grey diamonds represent the expected volume density while the black diamonds show the measured volume density. Regions where the measured volume density is higher than what is expected is highlighted in green while the inverse is highlighted in red. It is apparent that there are excess AGN in the dual and offset+dual subsamples and that the excess is separation dependent. 

![Figure 6](https://github.com/jlsteffen/AGN-in-Mergers/blob/main/images/21vol_den_sep.png)

![Figure 7](https://github.com/jlsteffen/AGN-in-Mergers/blob/main/images/21vol_den_sep_model.png.png)

![Figure 8](https://github.com/jlsteffen/AGN-in-Mergers/blob/main/images/21vol_den_sep_fit.png.png)

This work has been published in [Steffen et al. 2023](https://ui.adsabs.harvard.edu/abs/2023ApJ...942..107S/abstract).

