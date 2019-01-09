
all: clean bitstream-vera


BITSTREAM_VERA_VERSION= `cat bitstream-vera-VERSION`
BITSTREAM_VERA_FILE= ttf-bitstream-vera-$(BITSTREAM_VERA_VERSION).tar.bz2
BITSTREAM_VERA_URL= http://ftp.gnome.org/pub/GNOME/sources/ttf-bitstream-vera/$(BITSTREAM_VERA_VERSION)/$(BITSTREAM_VERA_FILE)

bitstream-vera:
	wget $(BITSTREAM_VERA_URL)
	mkdir bitstream-vera-fonts
	tar xjvf $(BITSTREAM_VERA_FILE) -C bitstream-vera-fonts --strip-components=1
	rm $(BITSTREAM_VERA_FILE)
	cp bitstream-vera-fonts/COPYRIGHT.TXT ../../LICENCE
	webify --no-eot --no-svg bitstream-vera-fonts/*.ttf


clean:
	find . \
		! -name "." ! -name ".." \
		! -name "README.md" \
		! -name "Makefile" \
		! -name "*-VERSION" | \
		xargs rm -rf
