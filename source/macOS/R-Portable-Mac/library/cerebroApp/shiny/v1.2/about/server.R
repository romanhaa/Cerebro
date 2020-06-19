##----------------------------------------------------------------------------##
## Tab: About.
##----------------------------------------------------------------------------##

output[["about"]] <- renderText({
  paste0(
    '<b>Version of cerebroApp</b><br>',
    'v1.2.2<br>
    <br>
    <b>Author</b><br>
    Roman Hillje<br>
    Department of Experimental Oncology<br>
    IEO, European Institute of Oncology IRCCS, Milan<br>
    <br>
    <b>Links</b><br>
    <ul>
      <li><a href=https://romanhaa.github.io/cerebroApp/ title="Official cerebroApp website" target="_blank"><b>Official cerebroApp website</b></a></li>
      <li><a href=https://github.com/romanhaa/Cerebro title="Official Cerebro repository on GitHub" target="_blank"><b>Official Cerebro repository on GitHub</b></a></li>
      <li><a href=https://github.com/romanhaa/Cerebro/releases title="Cerebro releases" target="_blank"><b>Cerebro releases</b></a></li>
      <li><a href=https://github.com/romanhaa/Cerebro/tree/master/examples title="Cerebro example data sets" target="_blank"><b>Cerebro example data sets</b></a></li>
    </ul>
    <br>
    <b>Citation</b><br>
    If you used Cerebro for your research, please cite the following publication:
    <br>
    Roman Hillje, Pier Giuseppe Pelicci, Lucilla Luzi, Cerebro: Interactive visualization of scRNA-seq data, Bioinformatics, btz877, <a href=https://doi.org/10.1093/bioinformatics/btz877 title="DOI" target="_blank">https://doi.org/10.1093/bioinformatics/btz877</a><br>
    <br>
    <b>Contact</b><br>
    <a href="mailto:roman.hillje@ieo.it?subject=Cerebro">roman.hillje@ieo.it</a><br>
    <br>
    <b>License</b><br>
    Cerebro is distributed under the terms of the <a href=https://github.com/romanhaa/Cerebro/blob/master/LICENSE.md title="MIT license" target="_blank">MIT license.</a><br>
    <br>
    <b>Credit where credit is due</b><br>
    <ul>
      <li>Sample and cluster color palettes taken from <a href="https://flatuicolors.com/" title="Flat UI Colors 2" target="_blank">https://flatuicolors.com/</a></li>
    </ul>
    <br>
    <b>Preferences</b>'
  )
})

output[["preferences"]] <- renderUI({
  tagList(
    checkboxInput("webgl_checkbox", label = "Use WebGL", value = TRUE)
  )
})

observeEvent(input[["webgl_checkbox"]], {
  preferences$use_webgl <- input[["webgl_checkbox"]]
  print(paste0("WebGL status is now: ", preferences$use_webgl))
})

observeEvent(input[["browser"]], {
  browser()
})

output[["logo_Cerebro"]] <- renderImage({
  list(
    src = system.file('extdata/logo_Cerebro.png', package = 'cerebroApp'),
    contentType = 'image/png',
    width = 350,
    height = 405,
    alt = "Cerebro logo",
    align = "right"
  )},
  deleteFile = FALSE
)
