# SyntenyPlot

	usage: synteny_plot.pl global.config synteny.config annotation.config output.svg

## Descriptions:

	Produce publishable vector-format synteny views using user-defined parameter and datasets

	SVG format further improvement could be achieved using Adobe Illustrator on Windows or Inkscape on Linux. 

	>Click *block4.svg* for an example, to see what the output svg looks like.

## Requirements: 

	Perl Modules: SVG, Data::Dumper

#HOW-TO
### 1. **user.config**

user.config contains plot parameter for the plot, including length, height, and color of line and fillins

>		Better to defined the sequence length here as autodetected length by this script might not be the real length

### 2. **synteny.config**

synteny.config specifies the syntonic region between two sequences

>		*Only first 6 column are necessart
>		*END not necessarily larger than start
>		*Could use whole.genome.alignments.pipeline.sh to generate
>		*whole.genome.alignments.pipeline.sh use LAST for pairwise alinment and ChainNet to join close alignments into larger blocks
>		*Final format:	
>			#REF1	START1	END1	REF2	START2	END2	...

### 3. **anntation.config**

anntation.config specifies the gene/feature/repeat locations

>		*Feature ID of non-EXPRESSION tracks must be unique
>		*For EXPRESSION track, note the FEATURE_ID and TISSUE change in col2 and col3, and the STDEV is before the AVERAGE;
>		*check /examples/block4/annotation.config for examples
>		*Curently CDS might not supplorted very vell as can not link the individual CDS box using lines.

### TO-DO list;

+ [ ] 1. Support CDS TRACKs
+ [ ] 2. Make delicate legends.
+ [ ] 3. Avoid Tick marks out of main plot
+ [ ] 4. Regenerate examples

## Citation:

To be published.

## Author:

---------------------------------------------------------------------

>	**Fu-Hao Lu**

>	Post-Doctoral Scientist in Micheal Bevan laboratory

>	Cell and Developmental Department, John Innes Centre

>	Norwich NR4 7UH, United Kingdom

>	E-mail: <Fu-Hao.Lu@jic.ac.uk>

---------------------------------------------------------------------

## Copyright

Copyright (c) 2016-2018 Fu-Hao Lu

Under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version, Permission is hereby granted, 
free of charge, to any person obtaining a copy of this software and 
associated documentation files (the "Software"), to deal in the Software 
without restriction, including without limitation the rights to use, 
copy, modify, merge, publish, distribute, sublicense, and/or sell 
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
