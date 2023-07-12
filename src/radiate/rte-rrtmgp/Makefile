#
# Top-level Makefile
#
.PHONY: rams docs
all: rams docs

rams:
	make -C build for_rams
	make -C interface all

docs:
	@cd doc; ./build_documentation.sh

clean:
	make -C build clean
	make -C interface clean
	rm -rf public
