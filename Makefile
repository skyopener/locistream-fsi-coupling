include Stream.conf

all: fluidPhysics stream MMStest cavitation flamelet

fluidPhysics: FRC
	cd fluidPhysics; $(MAKE)

stream: fluidPhysics
	cd src; $(MAKE)

cavitation: FRC
	cd cavitation; $(MAKE)

flamelet: FRC
	cd flamelet; $(MAKE)

MMStest: fluidPhysics
	cd MMStest; $(MAKE)

clean: FRC
	cd fluidPhysics; $(MAKE) clean
	cd src; $(MAKE) clean
	cd cavitation; $(MAKE) clean
	cd flamelet; $(MAKE) clean
	cd MMStest; $(MAKE) clean
	rm -rf bin lib

FRC:
