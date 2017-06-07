#
# Makefile to convert DocBook CLM Users-Guide into html and/or pdf
# (rtf, txt, ps, tex, man, dvi, and texi are also valid docbook formats)
#
VPATH := ../../tools . .. ../../bld ../../tools/ncl_scripts \
	../../tools/mksurfdata_map ../../tools/mkmapdata ../../test/system ../../bld/namelist_files \
	../../bld/config_files ../../../../../scripts/ccsm_utils/Tools \
        ../../src/main

PDFUG     := clm_ug.pdf
HTMLUG    := book1.html
DOCBKUG   := clm_ug.xml
CFGLOG    := config_help
CPRLOG    := cprnc_readme
TOOLOG    := tools_readme
CO2STR    := datm.streams.txt.co2tseries
BNMLOG    := buildnml_help
RESLOG    := buildnml_resolutions
USCLOG    := buildnml_usecases
QCKLOG    := quickstart_guide
COPLOG    := filecopies
MKSLOG    := mksurfdata.pl
MKGLOG    := mknoocnmap.pl
MKMLOG    := mkmapdata.sh
USRLOG    := quickstart_usrdat
TDRLOG    := test_driver.sh
GETREG    := getregional_datasets
DATLOG    := build_date
NMLDFTBL  := namelist_definition_table
NMLDLTBL  := namelist_defaults_clm_table
HSFLDTBL  := history_fields_table
CFGDFTBL  := config_definition_table
COMPLIST  := compsets_list_ofIcases.xml
MODNLIDX  := index.html
SOURCES   := $(DOCBKUG) $(COMPLIST) $(CFGLOG).xml $(BNMLOG).xml \
             $(RESLOG).xml $(USCLOG).xml $(QCKLOG).xml $(COPLOG).xml \
             $(USRLOG).xml $(GETREG).xml preface.xml custom.xml special_cases.xml \
             tools.xml adding_files.xml single_point.xml \
             appendix.xml trouble_shooting.xml $(CO2STR).xml \
             $(MKSLOG).xml $(MKGLOG).xml $(MKMLOG).xml $(TDRLOG).xml $(DATLOG).xml $(CPRLOG).xml $(TOOLOG).xml \
	     $(NMLDFTBL).xml $(NMLDLTBL).xml $(CFGDFTBL).xml $(HSFLDTBL).xml modelnl/$(MODNLIDX)

CONVAMP    := sed 's/\&/\&amp;/g'
CONVSIGNS  := sed 's/>/\&gt;/g' | sed 's/</\&lt;/g'

all: $(HTMLUG) $(PDFUG) modelnl/$(MODNLIDX)

.SUFFIXES:
.SUFFIXES: .xml .pdf .html .tlog .user_nl .diff .xhtml .xsl

RM := /bin/rm

$(CFGLOG).xml: $(CFGLOG).tlog
$(BNMLOG).xml: $(BNMLOG).tlog
$(RESLOG).xml: $(RESLOG).tlog
$(USCLOG).xml: $(USCLOG).tlog
$(QCKLOG).xml: $(QCKLOG).tlog
$(COPLOG).xml: $(COPLOG).tlog
$(USRLOG).xml: $(USRLOG).tlog
$(GETREG).xml: $(GETREG).tlog
$(MKSLOG).xml: $(MKSLOG).tlog
$(MKGLOG).xml: $(MKGLOG).tlog
$(MKMLOG).xml: $(MKMLOG).tlog
$(CPRLOG).xml: $(CPRLOG).tlog
$(TOOLOG).xml: $(TOOLOG).tlog
$(DATLOG).xml: $(DATLOG).tlog
$(CO2STR).xml: $(CO2STR).tlog

$(NMLDFTBL).xhtml: namelist_definition.xml namelist_definition.xsl
	xsltproc ../../bld/namelist_files/namelist_definition.xsl $< > tmpFile.txt
	./addxhtmlhead.pl tmpFile.txt > $@
	$(RM) tmpFile.txt

$(NMLDLTBL).xhtml: namelist_defaults_clm.xml namelist_defaults.xsl
	xsltproc ../../bld/namelist_files/namelist_defaults.xsl $< > tmpFile.txt
	./addxhtmlhead.pl tmpFile.txt > $@
	$(RM) tmpFile.txt

$(CFGDFTBL).xhtml: config_definition.xml config_definition.xsl
	xsltproc ../../bld/config_files/config_definition.xsl $< > tmpFile.txt
	./addxhtmlhead.pl tmpFile.txt > $@
	$(RM) tmpFile.txt

history_fields.xml: findHistFields.pl
	cd ../../src/main; ./findHistFields.pl >& tmpFile.txt
	$(RM) tmpFile.txt

$(HSFLDTBL).xhtml: history_fields.xml history_fields.xsl
	xsltproc ../../bld/namelist_files/history_fields.xsl $< > tmpFile.txt
	./addxhtmlhead.pl tmpFile.txt > $@
	$(RM) tmpFile.txt

.xhtml.xml:
	xsltproc stylesheethtml2docbook.xsl $< > $@

.tlog.xml:
	$(CONVAMP) $< | $(CONVSIGNS) > tempFile.txt
	./limitLineLen.pl tempFile.txt  > $@
	$(RM) tempFile.txt

.diff.xml:
	$(CONVAMP) $< | $(CONVSIGNS) > $@

.user_nl.xml:
	$(CONVAMP) $< | $(CONVSIGNS) > $@

debug:
	@echo "SOURCES:   $(SOURCES)"
	@echo "CONVAMP:   $(CONVAMP)"
	@echo "CONVSIGNS: $(CONVSIGNS)"

$(COMPLIST): 
	./get_Icaselist.pl > $@

$(HTMLUG): $(SOURCES)
	docbook2html --dsl clm_stylesheet.dsl#html $<

$(PDFUG): $(SOURCES)
	docbook2pdf --dsl clm_stylesheet.dsl#print $<

modelnl/$(MODNLIDX): 
	cd modelnl && $(MAKE) $(MAKEFLAGS)

$(BNMLOG).tlog: build-namelist
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../bld/build-namelist -help >& $@

$(DATLOG).tlog:
	@echo "Get current build date"
	date  +%b-%d-%Y >& $@

$(RESLOG).tlog: build-namelist
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../bld/build-namelist -res list >& $@

$(USCLOG).tlog: build-namelist
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../bld/build-namelist -use_case list >& $@

$(CFGLOG).tlog: configure
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../bld/configure -help >& $@

$(MKSLOG).tlog: mksurfdata.pl
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../tools/mksurfdata_map/mksurfdata.pl -help >& $@

$(MKGLOG).tlog: mknoocnmap.pl
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../tools/mkmapdata/mknoocnmap.pl -help >& $@

$(MKMLOG).tlog: mkmapdata.sh
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../tools/mkmapdata/mkmapdata.sh --help >& $@

$(TDRLOG).tlog: test_driver.sh
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $@ is good and redo your make"
	../../test/system/test_driver.sh -help >& $@

$(QCKLOG).tlog: Quickstart.GUIDE
	cp $<  $@

$(COPLOG).tlog: README.filecopies
	cp $<  $@

$(CPRLOG).tlog: README
	cp $<  $@

$(CO2STR).tlog: datm.streams.txt.co2tseries
	cp $<  $@

$(TOOLOG).tlog: README
	cp $<  $@

$(USRLOG).tlog: Quickstart.userdatasets
	cp $<  $@

$(GETREG).tlog:  getregional_datasets.pl
	@echo "The following line will fail in the make as it calls die -- but that is expected"
	@echo "Check that the output $(GETREG) is good and redo your make"
	../../tools/ncl_scripts/getregional_datasets.pl -help >& $@

clean:
	$(RM) -f $(HTMLUG) $(PDFUG) *.tlog $(DATLOG).xml *.xhtml *.tex

realclean: clean
	$(RM) -f f*.html c*.html x*.html a*.html i*.html $(COMPLIST) $(CFGLOG).xml \
	$(BNMLOG).xml \
	$(RESLOG).xml $(USCLOG).xml $(USRLOG).xml $(GETREG).xml $(QCKLOG).xml \
	*.tlog $(MKSLOG).xml $(MKGLOG).xml $(MKMLOG).xml $(TDRLOG).xml $(DATLOG).xml \
	$(NMLDFTBL).xml $(NMLDLTBL).xml $(CFGDFTBL).xml $(CPRLOG).xml $(TOOLOG).xml \
	$(COPLOG).xml $(CO2STR).xml

