##----------------------------------------------------------------------------##
## Tab: About.
##----------------------------------------------------------------------------##

output[["about"]] <- renderText({
  paste0(
    '<b>Version of cerebroApp:</b><br>',
    packageVersion("cerebroApp"), '<br>
    <br>
    <b>Author:</b><br>
    Roman Hillje<br>
    Department of Experimental Oncology<br>
    IEO, European Institute of Oncology IRCCS, Milan<br>
    <br>
    <b>Links:</b><br>
    <ul>
      <li><a href=https://github.com/romanhaa/Cerebro title="Official Cerebro repository on GitHub" target="_blank"><b>Official Cerebro repository on GitHub</b></a></li>
      <li><a href=https://github.com/romanhaa/Cerebro/releases title="Cerebro releases" target="_blank"><b>Cerebro releases</b></a></li>
      <li><a href=https://github.com/romanhaa/Cerebro/tree/master/examples title="Cerebro example data sets" target="_blank"><b>Cerebro example data sets</b></a></li>
    </ul>
    <br>
    <b>Contact:</b><br>
    <a href="mailto:roman.hillje@ieo.it?subject=Cerebro">roman.hillje@ieo.it</a><br>
    <br>
    <b>License:</b><br>
    Cerebro is distributed under the terms of the <a href=https://github.com/romanhaa/Cerebro/blob/master/LICENSE.md title="MIT license" target="_blank">MIT license.</a><br>
    <br>
    <b>Credit where credit is due:</b><br>
    <ul>
      <li>App icon made by <a href="https://www.flaticon.com/authors/kiranshastry" title="Kiranshastry" target="_blank">Kiranshastry</a> from <a href="https://www.flaticon.com/" title="Flaticon" target="_blank">www.flaticon.com</a> is licensed by <a href="http://creativecommons.org/licenses/by/3.0/" title="Creative Commons BY 3.0" target="_blank">CC 3.0 BY</a></li>
      <li>Sample and cluster color palettes taken from <a href="https://flatuicolors.com/" title="Flat UI Colors 2" target="_blank">https://flatuicolors.com/</a></li>
    </ul>'
  )
})