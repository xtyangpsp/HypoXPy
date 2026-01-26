MK	= make
SUBDIRS	= hypoDD ncsn2pha ph2dt hista2ddsta
BINDIR = ~/bin
all:
	for d in $(SUBDIRS); do (cd $$d; echo $$d; $(MK) $@; cp $$d $(BINDIR)); done
	
clean:
	for d in $(SUBDIRS); do (cd $$d; echo $$d; $(MK) $@; rm $(BINDIR)/$$d); done


