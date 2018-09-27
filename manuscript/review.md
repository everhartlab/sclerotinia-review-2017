---
title: "Something in the agar does not compute: On the discriminatory power of mycelial compatibility in *Sclerotinia sclerotiorum*"
mainfont: Georgia
fontsize: 12pt
geometry:
  margin=1in
author:
  - Zhian N. Kamvar
  - Sydney E. Everhart
institute: "University of Nebraska-Lincoln, Lincoln, NE"
bibliography: references.bib
header-includes:
  - \usepackage{lineno}
  - \linenumbers
csl: springer-basic-author-date-no-et-al.csl
abstract: |
  Mycelial compatibility, the ability for fungal isolates to grow together and form one single colony, was defined for *Sclerotinia sclerotiorum* nearly 30 years ago and has since been used as a marker to describe clonal variation in population genetic studies.
  While evidence suggests an associative relationship between mycelial compatibility and vegetative compatibility, contemporary research has treated these traits as analogous.
  As molecular markers have been developed to describe genetic variation, researchers combined these with the mycelial compatibility groups to assess and define clonal lineages.
  However, several inconsistent relationships between mycelial compatibility groups, haplotypes, and even vegetative compatibility groups have been observed throughout the literature, suggesting that mycelial compatibility may not accurately reflect self-recognition.
  We argue that the *Sclerotinia* community needs to move beyond using MCG data in
  population genetic studies.
---

Key Words: vegetative compatibility, microsatellite genotyping

# A brief introduction to vegetative compatibility

Vegetative compatibility is a trait in eukaryotes that allows two individuals to fuse and grow as one; it is often seen as a self-recognition system [@Leslie1993-hj; @Glass2000-cg].
In fungi, this process of fusion between hyphae is known as anastomosis.
Incompatible reactions from non-self recognition, result in a barrage line: a band of melanized and/or dead tissue at the contact zone between isolates [@Leslie1993-hj; @Glass1992-af; @Glass2000-cg; @Strom2016-di].
The barrage line, however is one of four possible states in hyphal contact—C0, no contact; C1, hyphal contact without fusion; C2, fusion with subsequent killing reaction; and C4, fusion with a stable heterokaryon [@Carling1988-xk].
To address this ambiguity, a method using complementary auxotrophic nitrate nonutilizing mutant strains grown on minimal media was developed to accurately assess stable heterokaryon formation [@Puhalla1985-bq; @Correll1987-hh; @Kohn1990-po].

For several decades, fungal biologists have used this system to classify intra-species groups called Vegetative Compatibility Groups (VCG) as a means to assess the diversity of natural fungal populations.
It has been shown that this trait has a genetic component in which two isolates must have identical alleles across several unlinked incompatibility loci in order to be compatible and successfully fuse.
The number of loci that regulate this interaction within each species of fungi is highly variable, and since only a single mismatch at one locus is required to induce an incompatible reaction, it is possible for a species to have anywhere from 4 (e.g. 2 loci = 2^2^ combinations) to over 1000 (e.g. 10 loci = 2^10^ combinations) VCG.

Because of the genetic component, the use of VCG to delineate intraspecific variation was an important addition to the existing classifications, often based on severity, which was affected by environmental factors [@Correll1987-hh].
This was used in *Fusarium oxysporum* to collate over 100 forma specialis based on severity into 21 VCG, giving researchers an indirect method for determining pathovar that was rapid compared to the use of differential plant lines able to differentiate these pathovars [@Correll1987-hh; @Puhalla1985-bq].
One of the most powerful applications of VCG has been shown in the clonal pathogen *Aspergillus flavus*.
This pathogen showed little to no gene flow between VCG [@Grubisha2010-ld].
As it was found that VCG could delineate aflatoxin-producing strains, researchers showed that non-aflatoxin-producing strains could reliably be used as a biocontrol agent [@Atehnkeng2016-qb; @Strom2016-di].
For many fungal species, before the advent of molecular genetic markers, VCG were routinely used to identify clonal lineages within populations.

# Mycelial compatibility in *Sclerotinia sclerotiorum*

The term mycelial compatibility was coined for *Sclerotinia sclerotiorum* in 1990 when an attempt was made to develop a protocol for determining the compatibility groups within populations [@Kohn1990-po].
Mycelial compatibility was defined as the ability for two isolates to grow together without the formation of a barrage reaction.
Attempts at assessing vegetative compatibility directly through mutagenesis were made, but only one was successful with a small sub-sample of isolates [@Kohn1990-po; @Ford1995-wk].
Despite an initial definition indicating that mycelial compatibility was a subset of vegetative compatibility, several contemporary publications treat them as analogous [@Kohn1990-po; @Schafer2006-ph; @Lehner2017-ny].

Not long after the initial definition of Mycelial Compatibility Groups (MCG) for *S. sclerotiorum*, these were paired with restriction fragment length polymorphism (RFLP) fingerprint profiles, and showed near concordance between the two [@Kohn1991-wq].
Due to the diversity of MCG, subsequent studies have also found concordance, but it did not take long for uncoupling of genotype and compatibility group to be reported (Table \ref{tab:meta}).
While the genetic component of non-self recognition is still not fully understood in *S. sclerotiorum*, 44 candidate genes have recently been identified, which may explain the diversity of MCG found across studies [@Wu2017-dx].

Table: \label{tab:meta} A selection of studies that included both MCG and MLH characterization of *Sclerotinia sclerotiorum* in their analysis and results of their work.
N = number of samples, MLH = number of multilocus haplotypes, MCG = number of mycelial compatibility groups, MCG concordance = number of MCG that belonged to a single MLH.

Study                       N  MLH  MCG  MCG concordance
----------------------- ----- ---- ---- ----------------
@Kohn1991-wq               63   33   32            31/32
@Kohli1992-pe              66   39   36            33/36
@Cubeta1997-rr             84   50   41            32/41
@Carbone2001-hi            64   19   15             7/15
@Hambleton2002-an         213   41   21            16/21
@Atallah2004-es           167  145   82               NA
@Malvarez2007-jo          294  152  141          130/141
@Attanayake2012-mq         40   16   15             7/15
@Barari2012-dn             65   44   26            22/26
@Lehner2015-oj            118   70   14             4/14
@Lehner2017-mm            187   46   34            24/34
@Kamvar2017-cl            366  165   87            15/87
--------------------------------------------------------


In this review, we will address inconsistencies between haplotype and MCG assignments, discuss proposed hypotheses that explain these patterns, and suggest a path forward for future population genetic studies of *S. sclerotiorum*.
The first step is to acknowledge the fact that mycelial compatibility is not equivalent to vegetative compatibility.

# Mycelial compatibility is not vegetative compatibility

Mycelial compatibility in *S. sclerotiorum* is considered by some to be synonymous with vegetative compatibility [@Schafer2006-ph; @Kohn1990-po], yet evidence to the contrary exists in *S. sclerotiorum* [@Ford1995-wk], *S. homoeocarpa* [@Jo2008-ft; @Chang2014-rn], *Verticillium dalhiae* [@Papaioannou2014-pe], and *Neurospora crassa* [@Micali2003-li].
The crux in these cases is that barrage formation may only serve as indirect evidence for vegetative incompatibility.
In order to be absolutely certain anastomosis has taken place, the most efficient method is to pair complementary auxotrophic nitrate non-utilizing (*nit*-) mutants on minimal media [@Kohn1990-po; @Ford1995-wk; @Leslie1993-hj; @Glass1992-af].
Vegetatively compatible strains anastomose and produce aerial hyphae and regain wildtype growth patterns, whereas incompatible strains fail to grow larger than the control strains, which are spindly and fail to produce aerial hyphae.
Unfortunately, *S. sclerotiorum* has proven notoriously difficult to mutate, which is the reason this technique has not yet been applied to determine mycelial compatibility [@Kohn1990-po; @Ford1995-wk].
To date, only one study has been able to create auxotrophic mutants to test for anastomosis and heterokaryosis in *S. sclerotiorum* [@Ford1995-wk].

One explanation for the lack of congruence between MCG and VCG is that MCG represent subsets of VCG, such that one MCG may represent more than one VCG.
While there is evidence that MCG are nested within VCG in *S. homoeocarpa* [@Chang2014-rn], the relationship between mycelial compatibility and vegetative compatibility in *S. sclerotiorum* is more complex, as shown in Table 4 of Ford et al. [-@Ford1995-wk] where one MCG was represented by two VCG and, in turn, both of these VCG represented two MCG.
Microscopically, however, it was shown that the hyphae remained in a stable, heterokaryotic state, as evidenced by 4,6-diamidino-2-phenylindole (DAPI) staining and stable transfers of colonies.
While this is the only study that was able to perform this experiment, there are several important conclusions: 1) vegetatively compatible strains are not necessarily identical, 2) mycelially compatible strains do not necessarily anastomose, 3) in *S. sclerotiorum*, vegetative compatibility may serve as an avenue to expand genetic variability through the parasexual cycle [@Ford1995-wk; @Strom2016-di], and 4) neither mycelial compatibility nor vegetative compatibility alone have the power to detect clones.
In the next section, we will examine the relationship between MCG and haplotype.

# Do mycelial compatibility groups represent clonal lineages?

Both MCG and haplotypes are controlled by multiple, independent loci and are thus expected to reflect clonal lineages—groups of closely related fungal strains derived from a single ancestral strain via asexual replication.
With a sufficiently large number of loci, the probability for two unrelated strains to converge on the same MCG or haplotype is extremely low [@Kohli1992-pe].
However, when MCG and haplotypes have been combined to detect clonal lineages in the past, the results revealed a disconnect between MCG and haplotype [@Cubeta1997-rr; @Atallah2004-es; @Attanayake2012-mq; @Kamvar2017-cl].

The conflict between MCG and VCG documented by Ford et al. [-@Ford1995-wk], while perplexing, may serve to explain complex relationships that continue to be observed between MCG and haplotypes [@Schafer2006-ph; @Attanayake2012-mq; @Kamvar2017-cl].
There are three possible relationships we can observe between a set of MCG and haplotypes in a population (we will be referring to these relationships throughout the manuscript): 1) A one-to-one relationship of a single haplotype to a single MCG (congruent, Fig. \ref{fig:mcg-figure}A), which is the expected outcome; 2) Several closely-related haplotypes belong to one MCG (trans-haplotype, Fig. \ref{fig:mcg-figure}B), which would be expected in a population that has undergone a clonal expansion; and 3) A single haplotype belongs to more than one MCG (trans-MCG, Fig. \ref{fig:mcg-figure}C), which would be expected for genotypes produced vi independent recombination events [@Parks1993-nv; @Arnaud-Haond2007-zo; @Cubeta1997-rr; @Schafer2006-ph].
All of these relationships have been observed in the nearly 30 years since MCG were first elucidated for *S. sclerotiorum*.
Below, we highlight a few examples that described these relationships.

![
Schematic representation of three possible interactions between haplotype (shaded cartoon DNA strands) and MCG (as depicted by cartoon mycelia on petri plates; grey plates represent redundant pairings).
An incompatible reaction is seen as a thick, black line down the middle of the plate.
A) Congruent relationship: two haplotypes, two MCG showing 1:1 relationship between haplotype and MCG.
Unlike haplotypes produce barrage reaction.
B) trans-haplotype relationship: two haplotypes, one MCG showing closely related haplotypes sharing the same MCG given by the lack of a barrage reaction.
C) trans-MCG relationship: one haplotype, two MCG showing a barrage reaction between identical genotypes. ](../figures/mcg-figure-small.png){#fig:mcg-figure width="100%"}

The first study to compare *S. sclerotiorum* MCG with haplotypes analyzed 63 isolates via RFLP probes in both nuclear and mitochondrial
DNA, detected 32 MCG and 33 haplotypes[^1] [@Kohn1991-wq].
Nearly all relationships observed were congruent, with the exception of one MCG that exhibited a trans-haplotype relationship; the two isolates that comprised that MCG expressed different haplotypes for the mitochondrial probe.
Ultimately, it was concluded that MCG represented unique clones.
In 1992, a study of 66 isolates comprising 36 MCG and 39 haplotypes showed three MCG exhibited trans-haplotype relationships [@Kohli1992-pe].
The authors suggested that because it was unlikely that the same MCG would arise by chance in different isolates AND the probability of sharing the same haplotype was 1 in a million, clones should be determined by the combination of MCG and haplotype.
Since it was impossible to ascertain inheritance from the probes used (and thus confirm that haplotypes within an MCG are closely related), they inferred that the observed trans-haplotype relationships were caused by recombination or mutation that did not occur at the *het* loci [@Kohli1992-pe].

Inconsistencies began to appear in subsequent studies of genetic diversity using both molecular markers and MCG.
For example, Cubeta et al. [-@Cubeta1997-rr] surveyed 100 isolates from North Carolina and Louisiana, USA.
They found all three MCG to haplotype relationships in North Carolina, where—out of 41 MCG—32 were congruent relationships, 10 MCG had a trans-haplotype relationship with more than one haplotype and two haplotypes had a trans-MCG relationship with more than one MCG identified and confirmed across two laboratories.
The trans-haplotype relationships were hypothesized to be due to transposition in one of the probes and the trans-MCG relationships were thought to represent partial compatibility, driven by a match at some, but not all *het* loci [@Cubeta1997-rr].

No studies since Kohn et al. [-@Kohn1991-wq] have shown a 1:1 relationship between haplotype and MCG (Table 1).
However, the assumption that strains within an MCG should be genetically identical has never been strict [@Leslie1993-hj].
Three mechanisms have been hypothesized for trans-haplotype relationships: recombination, transposable elements, and mutation.
We know that there are a minimum of 8 and a maximum of 44 loci controlling the determination of vegetative compatibility, which would make trans-haplotype relationships unlikely by recombination alone [@Malvarez2007-jo; @Wu2017-dx; @Kohli1992-pe].
Moreover, because several marker types have been used to investigate haplotypes, it is unlikely that transposable elements are responsible for changing the haplotypes.
As was observed in Kamvar et al. [-@Kamvar2017-cl], all three relationships were present, with congruent relationships comprising only 15 of 87 MCG.
The average genetic distance between haplotypes within MCG (trans-haplotype) were not consistently smaller than between MCG, suggesting that simple mutation could not be the cause of these incongruencies.
Moreover, both trans-haplotype and trans-MCG relationships were entangled (Figure 1 in Kamvar et al. [-@Kamvar2017-cl]); there were MCG comprised of several haplotypes that were further comprised of several MCG such that, in one extreme case, one single haplotype could be connected with 78 other haplotypes through four MCG.

# Why are there inconsistencies?

The inconsistencies between MCG and haplotypes have been acknowledged throughout the literature [@Atallah2004-es; @Malvarez2007-jo; @Cubeta1997-rr; @Kamvar2017-cl].
A few hypotheses to explain these patterns have been proposed.
Here, we will highlight two of them and briefly explore their strengths and weaknesses.

Perhaps one of the most widely recognized hypotheses is that a decoupling of MCG and haplotype is a sign of recombination [@Kohli1992-pe; @Milgroom1996-we; @Phillips2002-pq; @Schafer2006-ph].
This concept has been demonstrated in the chestnut blight fungus *Cryphonectria parasitica* by comparing genetic similarity within VCG to randomly-assigned groups [@Milgroom1996-we; @Liu1996-dr].
Using this hypothesis, Schafer and Kohn [@Schafer2006-ph] specified expected relationships between MCG and haplotype under different reproductive scenarios.
Clonal populations are expected to correlate 1:1 or have closely related haplotypes within MCG due to clonal expansion.
However, in recombining populations, “...except for isolates from adjacent plants, the expectation is that each isolate sampled is either incompatible with all other isolates or is part of an intransitive[^2] MCG and each isolate either has a unique fingerprint or a fingerprint associated with more than one MCG” [@Schafer2006-ph].

![
Sampling from disparate populations can show signatures of linkage equilibrium ($\bar{r}_d$) in populations that are otherwise clonal, which is demonstrated in these 20 simulated populations.
A) Twenty simulated clonal populations (large shaded circles) with no migration between them.
Each population has 20 individuals sampled (white points) one filled circle in each population represents one sample to create the overlaid tree in B.
B) A neighbor-joining tree for all individuals in A (light grey branches).
Tip color represents population of origin.
Larger tips and dark shaded branches represent the tree created from sampling one individual from each population.
C) Distributions of the index of association calculated from 20 random samples within (shaded densities) and among (black density) populations with 99 replicates.
](../results/figures/iatree.png){#fig:iatree width="75%"}

Decoupling due to recombination is possible and expected in fungal species with few genes controlling vegetative compatibility, or when comparing mtDNA haplotypes, since transfer has been documented between fungal strains of different compatibility groups [@Gordon1992-fs; @Leslie1993-hj].
However, in *S. sclerotiorum*, we are more likely to observe entirely unique genotypes and MCG.
Consider that there are as few as eight and as many as 44 genes controlling mycelial compatibility; the odds of expressing the exact same MCG or genotype through independent recombination events are low [@Kohli1992-pe; @Wu2017-dx].
If this is the case, why do we consider recombination to be the cause of these incongruent MCG/haplotype relationships?
These decoupling patterns via recombination may be plausible in small populations with fixed alleles at several loci controlling mycelial compatibility.
For instance, a recent study found a total of 15 MCG, but only seven of them were congruent, which may suggest four loci controlling mycelial compatibility segregating in the population [@Attanayake2012-mq].
However, considering that identical MCG from both disparate regions and different haplotypes have been found, this explanation does not scale well---even considering the potential for migration via agricultural mechanisms [@Kohli1992-pe; @Ford1995-wk; @Kamvar2017-cl].

A recent hypothesis suggests that MCG constitute distinct populations.
In their recent review, Lehner and Mizubuti [@Lehner2017-ny]—while addressing the issue of apparent recombination within tropical regions—proposed that MCG would impose a barrier for genetic exchange and thus, MCG represent groups of recombining individuals.
For *S. sclerotiorum*, which has a homothallic mating system, recombination between identical genotypes would appear to be no different than clonal reproduction.
If populations were analyzed for linkage without taking MCG into account, the results could give signatures of random mating due to sampling from several independent clonal populations [@Prugnolle2010-yb; @Lehner2017-ny].
However, it has been shown in other ascomycetes that vegetative incompatibility is not a barrier to recombination and, in *S. sclerotiorum*, mycelial incompatibility is not a barrier to vegetative compatibility [@Perkins1988-mt; @Ford1995-wk; @Leslie1993-hj; @Milgroom1996-we].
Thus, if members of a single MCG were sampled from several disparate populations and treated as one, assessment of linkage within MCG may suffer from the same error (Fig. \ref{fig:iatree}).
This hypothesis, however can only address haplotype discordance within MCG and does not explain the identical haplotypes across MCG.

![
Diagram representing processes (A) governing the expression of observed patterns (B) of haplotypes, VCG, and MCG.
Population genetic processes like mutation, migration, and recombination directly affect haplotype while a combination of haplotype and unknown environmental processes (i.e. viruses) indirectly affect MCG and VCG expression.
Black, solid arrows represent direct causal relationships and grey, dashed arrows represent indirect causal relationships. ](../figures/process-pattern.png){#fig:process width="40%"}

One possibility for identical haplotypes across MCG could be explained by exogenous forces capable of altering mycelial compatibility reactions.
A virus capable of switching on mycelial compatibility was recently discovered [@Wu2017-dx].
*Sclerotinia sclerotiorum* mycoreovirus 4 (SsMYRV4) has been demonstrated to force two incompatible strains to become compatible by modifying expression of *het* genes and preventing apoptosis at the junction between the two strains.
It’s hypothesized that this facilitates viral transmission between hosts.
This study identified up to 44 putative genes involved in incompatibility reactions, but showed that the reaction can be modified through expression.
Thus, while mycelial compatibility may utilize the machinery of vegetative compatibility, there is evidence that it may still be altered by exogenous forces (Fig. \ref{fig:process}).
The extent to which this virus or other exogenous forces may play a role in mycelial compatibility interactions remains uncharacterized.

# The way forward in genetic studies

In order to ascertain the evolutionary processes that lead to observed population structure, one must use multiple independent, quantitative, selectively neutral, and heritable markers.
RFLP, microsatellite, and SNP data all meet these criteria; mycelial compatibility, on the other hand, does not.
It is impossible to infer anything about migration, mutation, or recombination with MCG alone [@Milgroom2003-mu; @McDonald1997-ob; @Goss2015-ue].
Mycelial compatibility groups in the *S. sclerotiorum* literature have roughly reflected genetic diversity, but there are a few key areas (unclear patterns of inheritance, immeasurable uncertainty, and no standard panel of strains for testing) that lead the authors of this paper to suggest the discontinued use of MCG as a marker for population genetic inference.

## Are haplotype and MCG better together?

Because genetic markers are heritable, we are able to assess the relationships between our samples, which can tell us something about the underlying genetic structure of the populations.
As *S. sclerotiorum* has the capacity for clonal reproduction, taking into account the haplotypic structure is also very important in assessing the clonal dynamics [@Arnaud-Haond2007-zo; @Grunwald2017-wd].
This is usually performed by borrowing tools from ecology to analyze the density and distribution of duplicated haplotypes in a population.
Both haplotypes and MCG were previously considered to be independent and, thus could be combined in order to elucidate the true distribution of clones [@Kohn1991-wq].
In practice, depending on the type of MLG/haplotype interactions observed, this technique renders the exact same distribution as haplotypes (in the case of trans-haplotype relationships, Fig. \ref{fig:mcg-figure}B) or a more diverse distribution (in the case of trans-MCG relationships, Fig. \ref{fig:mcg-figure}C).
However, due to the potentially non-independent nature of MCG, we strongly advise against joining haplotypes and MCG to determine clones, even if it means sacrificing perceived specificity [@Leslie1993-hj].

Both MCG and genetic markers can be used to identify clones in a population.
While the main use of MCG has been for self/nonself determination, there is clear documented evidence of non-self isolates producing a compatible reaction and vice-versa [@Ford1995-wk; @Leslie1996-di; @Wu2017-dx].
While this problem exists with molecular markers, the implications differ.
Identical haplotypes are generally considered to be the same individual, however, haplotypes can be identical in state without being identical by descent or mutation rendering two otherwise identical genotypes different [@Parks1993-nv; @Arnaud-Haond2007-zo].

While both of these markers suffer from uncertainty, haplotypes have a clear advantage over MCG as there are ways to quantify the uncertainty.
For example, to handle the presence of highly variable loci or missing data, one can use a genetic distance threshold to group similar haplotypes that may be part of a clonal lineage [@Arnaud-Haond2007-zo; @Kamvar2015-ff; @Bailleul2016-lw].
The metric *P~sex~* can be used to quantify the probability of obtaining identical genotypes by chance, allowing the researcher to subset identical haplotypes into different groups [@Parks1993-nv; @Arnaud-Haond2007-zo; @Bailleul2016-lw].
Unfortunately, these operations cannot be performed when MCG are included.
Thus, MCG should not be combined with genetic markers to determine self from non-self due to the qualitative nature of MCG.

## Reproducibility across studies

Even if we assume that MCG and haplotype measure the same phenomena and error could be calculated with equal confidence, there is still a difference in access and reproducibility across studies.
With the advent of RFLP fingerprinting, it was possible to develop a reference set of fingerprints for comparison with future studies [@Kohn1991-wq; @Hambleton2002-an; @Kohli1998-hh].
However, unlike other economically important plant pathogens such as *Aspergillus flavus*, no such database exists for mycelial compatibility groups [@Grubisha2010-ld].
Each lab thus performs their own analysis of mycelial compatibility, producing their own bespoke group assignments.
While these experiments are replicated using standardized media, there are rare cases of replication efforts across laboratories [@Cubeta1997-rr].
Of course, the same can be said of microsatellite genotypes, which came to replace the RFLP probes; no reference panel yet exists for microsatellite loci [@Sirjusingh2001-sq].

While both measures may be equally deficient in this respect, there is a measurable difference in the amount of effort necessary to constructing a reference panel of either marker.
Microsatellite loci are known to differ from study to study, but these differences are based on the machine used, primer modifications (ie. direct or indirect fluorescent labeling and poly-A tails), and choice in binning ambiguous alleles.
This means that 1) the difference should be relative across all allele calls between studies and 2) if representative isolates—the minimal number of isolates that would capture all alleles observed in a study, which could represent a small fraction of the number of unique haplotypes—from several studies were re-sequenced together, it is possible to adjust historical and future data to match the allele calls of the reference panel.
Thus, while there is no current reference panel of microsatellite haplotypes, it is possible to achieve this with a large enough selection of isolates representing known variation at each locus and requiring only a minimal amount of additional cost and effort to include as standards in any genotyping study.
In contrast, a MCG representative panel would require all represented groups from all studies to be compared together.
Moreover, MCG tester strains would need to be shared as living biological material, limiting access due to hazardous material import policies; genotype panels, on the other hand, can be shared as purified genomic DNA.

## Recommended Practices

For population genetic studies, it is recommended to follow best practices of population genetic analysis [@Grunwald2017-wd].
Researchers should genotype populations over several unlinked, selectively neutral genetic loci [@Sirjusingh2001-sq].
When assessing multilocus genotypes, researchers should calculate *P~sex~* to account for the probability of encountering a genotype more than once via independent meiotic events [@Arnaud-Haond2007-zo; @Parks1993-nv].
However, not all genetic loci are created equal.
The future of *S. sclerotiorum* population genetic research needs to use the high quality genome available to generate new and better marker systems [@Derbyshire2017-mx].
To enable reuse, haplotype data including all relevant metadata (population assignments) and analysis code should be deposited in public repositories such as the Open Science Framework ([*https://osf.io*](https://osf.io)).

## Are there still uses for Mycelial Compatibility?

Readers may be forgiven for coming to the conclusion that the authors wish to completely do away with using mycelial compatibility altogether.
While we strongly oppose their use in inference of population genetic processes, there are still instances that merit use of MCG testing.
While we have highlighted inconsistencies in the ability of MCG to determine clones, most of the assignments are largely congruent with haplotypes.
Thus, MCG are well suited as a rough measure of standing genotypic diversity for research groups that do not have the resources available for genotyping.
Because of the more rapid nature of MCG, these are ideal for situations in which a grower needs information about management decisions or pilot studies to determine appropriate sampling design for a population genetic study.

While mycelial compatibility may not have a direct relationship with genotype, it is still a phenotype that is related to self-recognition and mating.
Thus it may influence gene flow within and between *S. sclerotiorum* populations through either sexual or parasexual recombination [@Strom2016-di].
The environmental cues or processes that influence barrage formation in *S. sclerotiorum* have yet to be elucidated and may serve as fertile ground for future research.

# Conclusions

All genetic markers are imprecise, but imprecision is not the same as uselessness.
Researchers can use genetic markers for applications with different degrees of tolerance for imprecision.
For pilot studies to assess standing clonal diversity in a given population, sequencing full genomes of all isolates would be a waste of time and money when MCG could be rapidly employed with relatively little effort and cost.
In this review, we argue that if the goal of a particular research endeavor is to assess the evolutionary processes that drive diversity and population subdivision, MCG are too coarse a marker for such inference.
We have shown that MCG can be used to roughly detect genetically identical isolates, but cannot be used for more in-depth analysis due to the immeasurable type-I and type-II errors inherent to the system.

The benefits of DNA-based markers for population genetic inference far outweigh the cost burden: relatedness and differentiation can be estimated, error can be quantified, and—more importantly—data can easily be shared and harmonized in meta-population genetic studies.
Thus, given these advantages, combined with the falling cost of molecular genetic reagents and analyses and availability of previously described microsatellites and a high-quality reference genome, we recommend that investigators should employ neutral DNA-based markers for population genetic inference [@Sirjusingh2001-sq; @Attanayake2014-uy; @Derbyshire2017-mx].
For the past 20 years, we have been using mycelial compatibility as a tool when instead, we should view it as a biological pattern that may allow us to understand the natural history of this pathogen.

# Acknowledgements

We would like to thank Gerard Adams for stimulating discussions on the nature of vegetative compatibility in filamentous fungi.
We additionally would like to thank two anonymous reviewers for their comments that improved the clarity of the final manuscript.

# Funding Sources

Funding was provided for salaries and previous research on this topic that was also reviewed here.
This includes partial support from the Nebraska Agricultural Experiment Station with funding from the Hatch Act (Accession Number 1007272) through the USDA National Institute of Food and Agriculture, grant \#58-5442-2-209 from the USDA-ARS National Sclerotinia Initiative to SEE, and start-up funds from the University of Nebraska-Lincoln (UNL) to SEE.
Funders had no role in study design, data collection and analysis, decision to publish, or preparation of the manuscript.

 [^1]: We use haplotype here to refer to the combination of the nuclear and mitochondrial probes.
    Kohn et al. describes each of these as “fingerprints” or (confusingly) “phenotypes”.

 [^2]: The idea of transitivity here refers to the mathematical term.
    If MCG A=B; B=C; then A=C.
    An intransitive relationship would be A=B; B=C; A$\neq$C.
    However, intransitive MCG may not be the result of sexual reproduction as the mechanism controlling this would be that of an imprecise recognition system in which only a percentage of loci would need to match [@Neigel1983-gt].
    In our analysis of the literature, all *S. sclerotiorum* MCG have been transitive.

# References
