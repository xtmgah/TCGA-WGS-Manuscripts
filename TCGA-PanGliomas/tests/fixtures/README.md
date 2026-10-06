# Cairo blank-space regression fixture

`cairo_type3_spaces.pdf` contains only synthetic text and a red line. It was generated using R 4.5.1 and Cairo 1.18.4, with Roboto Condensed (Apache License 2.0; see `functions/fonts/NOTICE.txt`). It contains no study data.

```r
cairo_pdf("cairo_type3_spaces.pdf", family="Roboto Condensed", width=4, height=2)
par(mar=rep(0,4))
plot.new()
text(.5,.7,"Multivariable Cox model for OS",cex=.7)
text(.5,.3,"P ≥ 0.05; 95% CI",cex=.7)
segments(.1,.5,.9,.5,col="red")
dev.off()
```

This Cairo build represents spaces in empty Type3 fonts. The regression test verifies that normalization preserves all visible pixels at 150 and 300 dpi, text content, and vector drawing count after saving and reopening the PDF. The checked-in fixture keeps this regression test effective on older Cairo builds that do not emit this representation.
