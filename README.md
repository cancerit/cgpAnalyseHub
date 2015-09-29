cgpAnalyseHub
=========

cgpAnalyseHub has been created to try and simplify bulk retrieval of data from [CGHub](https://cghub.ucsc.edu/), remapping
through BWA-mem and analysis using the CGP variant calling algoritms [CaVEMan](http://cancerit.github.io/CaVEMan/)
and [cgpPindel](http://cancerit.github.io/cgpPindel/) followed by annotation using [Vagrent](http://cancerit.github.io/VAGrENT/).

To do this it utilises the [CGHub cart](https://browser.cghub.ucsc.edu/search/) summary file and arranges data in a uniform manner.

To download data from CGHub you will need to apply for access directly with the administrators of the resource, we are not able to help you with this process.

Code is also included to simplify bulk generation of commands and helper scripts for farm submissions if this facility is available to you.
A reasonable understanding of linux and your compute infrastructure is required.

### Non-Hub data

It is not necessary for data to be sourced from CGHub, the subsequent callers and helper code can be used
with any paired WXS data (tumour/normal) provided data is arranged in an expected format.  This is detailed
in the [wiki pages](../../wiki).

---

###Dependencies/Install

Please install the following first:

* [PCAP-core](http://github.com/ICGC-TCGA-PanCancer/PCAP-core/releases)

Please see this for any child dependencies.

Once complete please run:

./setup.sh /some/install/location

Note: The environment `PERL5LIB` variable is modified to only use libraries installed under the specified path,
please ensure that all items are installed to the same area.

To use the download functionality you will need to ensure that [genetorrent](https://cghub.ucsc.edu/software/downloads.html) is installed and available in your path environment.

If you wish to use the callers and not just the download facility you also need to install the following tools:

* [cgpVcf](http://cancerit.github.io/cgpVcf/)
* [cgpPindel](http://cancerit.github.io/cgpPindel/)
* [Vagrent](http://cancerit.github.io/VAGrENT/)
* [cgpCaVEManPostProcessing](http://cancerit.github.io/cgpCaVEManPostProcessing/)
* [cgpCaVEManWrapper](http://cancerit.github.io/cgpCaVEManWrapper/)

---

## Developer info

This is a hubflow managed project.  Please use appropriate feature, release and hotfix tools.

####Cutting the release
1. Update `lib/Sanger/CGP/AnalyseHub.pm` to the correct version (adding rc/beta to end if applicable).
2. Update `Changes.md` to show major items.
3. Run `./prerelease.sh`
4. Check all tests and coverage reports are acceptable.
5. Commit the updated docs tree and updated module/version.
6. Push commits.
7. Use the GitHub tools to draft a release.

LICENCE
=======

Copyright (c) 2015 Genome Research Ltd.

Author: Cancer Genome Project <cgpit@sanger.ac.uk>

This file is part of cgpAnalyseHub.

cgpAnalyseHub is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
2009, 2011-2012’ should be interpreted as being identical to a statement that
reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012’."
