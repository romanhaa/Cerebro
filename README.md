# Cerebro

![Screenshot Cerebro: overview panel](screenshots/Screenshot_2019-03-01_at_15.22.15.png?raw=true "Cerebro: overview panel")

![Video Cerebro](screenshots/Cerebro_2019-03-12.2019-03-12_14_37_20.gif?raw=true "Cerebro")

## Description

Motivation ...

Visualize and explore complex data sets from single-cell RNA-sequencing experiments.

Currently available only for Windows and macOS.
Linux users have the option to launch the application through a Docker container or hosting the application using the `app.R` file and a local R installation.

## Installation

Download release for your OS, unpack and run.

## Usage

* Export data set using the `cerebroPrepare` R package.
* Launch the browser and load the data set.

## How it works

* 

## FAQ

* 

## Package from source

### macOS

To package Cerebro you need [Git](https://git-scm.com) and [Node.js](https://nodejs.org/en/download/) (which comes with [npm](http://npmjs.com)) installed on your computer.
Then, from the command line, run:

```bash
# clone this repository
git clone https://gitlab.com/romanhaa/Cerebro.git
# install Electron packager
npm install electron-packager --global
# go into the repository
cd Cerebro
# install dependencies
npm install
# run the app
npm start
# build the app
npm run package-win
# or
npm run package-mac
```

### Windows

If you're using Linux Bash for Windows, [see this guide](https://www.howtogeek.com/261575/how-to-run-graphical-linux-desktop-applications-from-windows-10s-bash-shell/) or use `node` from the command prompt.

## Release History

* version X.X
  * ...

## Contributions

* <https://github.com/ColumbusCollaboratory/electron-quick-start>

## To Do

* [ ] Prepare Docker container and instructions on how to use it.
* [ ] Add contributions that are currently mentioned inside the browser.
* [ ] Logo inside app?
* [ ] Is there a more robust way to match gene (set) expression values to cells without assuming that columns in `expression` are in the same order as rows in `cells`?
* [ ] Include detailed sample information (path and data type) on 'Sample info' tab.
* [ ] Check where hallmark gene sets would appear in pathway enrichement.
* [ ] Add parameters of marker gene detection to info tab.
* [ ] Make panel with links to useful pages.
  * How UMAP works: <https://umap-learn.readthedocs.io/en/latest/how_umap_works.html>
  * How to interpret distances in UMAP : <https://github.com/lmcinnes/umap/issues/92>
  * How to effectively use t-SNE: <https://distill.pub/2016/misread-tsne/>
  * How to interpret distances in t-SNE: <https://stats.stackexchange.com/questions/263539/clustering-on-the-output-of-t-sne>
  * Problems of t-SNE: <https://stats.stackexchange.com/questions/270391/should-dimensionality-reduction-for-visualization-be-considered-a-closed-probl/270414>
* [ ] Check if `require()` could be a way to make startup faster.
* [ ] Gene expression tab: Update `textAreaInput` after checking which genes are present in data set? Or at least remove empty lines?
* [ ] Make `hoverinfo` background white, like in scatterD3.
* [ ] Check which packages can be removed. / Clean up R libraries.
* [ ] Common function that filters user-specified genes for the ones present in the data set?
* [ ] Common plotly layout for all box plots?
* [ ] Find out how to make separate object for alternative text obsolete.
* [ ] `per` to `by`?
* [ ] Custom hoverinfo for all plots.
* [ ] Info box for gene ID conversion.
* [ ] Display only the sample information that is available.
* [ ] Be more consistent with "dot", "cell" and "point".
* [ ] Check how well dimensional reductions scale in plotly.
  * 3D is doable up to 150,000 cells. At 300,000 cells it is barely usable on my MacBook Pro.
* [ ] Cell size by variable doesn't work anymore.
* [ ] Use loading animations.
* [ ] Plot dimensional reduction from overview panel with sample/cluster labels.
* [ ] Allow user to assign custom colors.
* [ ] Change name of box plots because now they are violin plots.
* [ ] Exporting dimensional reductions is misleading when using 3D because a plot is created but only the first 2 dimensions used. Either deactivate or find a way to plot an angled 3D with ggplot.
  * Add warning message when exporting.
* [ ] Do thorough performance check. What step could/should be improved? Gene (set) expression?
  * colMeans runs in C, super fast already.
