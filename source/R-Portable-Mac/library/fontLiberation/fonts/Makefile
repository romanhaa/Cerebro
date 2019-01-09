
all: clean liberation

LIBERATION_VERSION= `cat liberation-VERSION`
LIBERATION_FILE= liberation-fonts-ttf-$(LIBERATION_VERSION).tar.gz
LIBERATION_URL= https://fedorahosted.org/releases/l/i/liberation-fonts/$(LIBERATION_FILE)

liberation:
	wget $(LIBERATION_URL)
	mkdir liberation-fonts
	tar xjvf $(LIBERATION_FILE) -C liberation-fonts --strip-components=1
	rm $(LIBERATION_FILE)
	cp liberation-fonts/LICENSE ../../LICENSE
	webify --no-eot --no-svg liberation-fonts/*.ttf

clean:
	find . \
		! -name "." ! -name ".." \
		! -name "README.md" \
		! -name "Makefile" \
		! -name "*-VERSION" | \
		xargs rm -rf
