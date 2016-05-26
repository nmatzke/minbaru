# minbaru
Replicating the "baraminic distance correlation" calculations that creationists use on cladistic morphology datasets. Should be entertaining.

"Baraminology" is a field invented by Young-Earth Creationists (YECs) to attempt to discern the "created kinds" in biology.  Their logic, such as it is, is that since the book of Genesis says that God created organisms and commanded them to reproduce "after their kind", these "kinds" must be the originally created groups.  Creationists have added the idea that evolution can occur within "kinds", including speciation, substantial morphological change, etc. -- things that would be called "macroevolution" in standard evolutionary biology.  For example, many YECs will agree that all cats (lions, tigers, house cats, etc.) might be the same kind, or all canids (dogs, wolves, foxes, etc.)  

This solves several problems for YECs -- it reduces the number of species they have to imagine putting on Noah's Ark, and it allows them to dismiss much of the evidence for evolution with "that's just microevolution within the kind."  There are several problems with this whole strategy -- it goes far beyond a literal reading of Genesis, for example -- but the biggest one is that creationists have never been able to come up with a rigorous way to determine what the Really Really Distinct And Not Evolutionarily Connected kinds are.  The only real rule is "whatever the Bible says", but the Bible doesn't discuss most species.  The Bible does say humans were specially created, but even this "fact" doesn't tell creationists which hominin transitional fossils are humans, and which are apes -- and creationists disagree about which categories to put which fossils in!

The rational strategy would have been to agree "it looks like evolution happened, doesn't it?", but creationists invented a different strategy.  They took an old creationst term "baramin" and constructed "baraminology". Baraminology attempts to quantify the differences between species, and, using clustering, determine what the created kinds are.  This basically amounts to phenetics with clustering, but it is interesting [1] that creationists would even attempt to come up with a quantitative, replicable method based on their views.

The program the baraminologists use is called bdist.  It originally existed as a Perl script, bdist.pl, and later was made into a web service, bdistmds, that added bootstrapping and graphical MDS (multidimensional scaling) outputs. (See references/links below.)

The R code here replicates the basic bdist calculations and derivatives. My basic goal was to see exactly what the baraminologists are doing, since I suspect that cladistic datasets contain objective tree structure that is poorly represented in a phenetic analysis and MDS plots. I also suspect that typical baraminology practices, such as excluding data below a cutoff, and looking for clusters only within a narrowly-defined cladistic dataset, are objectively bad, even given their assumptions and goals.

Eventually I may get this up to a full R package on CRAN, and/or a publication -- although it is way down on my list for the moment. I may also use it just to practice R package construction and unit testing on a simple (pure R) code base.

Here is the story behind the name, "minbaru": The term "baramin" for created (bara) kind (min) is, I gather, ungrammatical Hebrew.  Various websites suggest that "minbaru" would be the correct way to say what the creationists are trying to say. I have no idea, personally, but I think the reversal of the term captures what my goals are here.

Notes

[1] Interesting, that is, at least to the community of "creationism watchers" that observe creationist shenanigans with a mixture of horror and fascination.


References

The main person supporting and promoting BDIST is Todd Wood at Bryan College (yes, that's William Jennings Bryan College, located in Dayton, TN, home of the 1925 Scopes Monkey Trial). Wood is one of the very few creationists with a nonzero understanding of modern evolutionary biology and its quantitative methods, and he has at various points admitted the incompetence of much creationist work, and the empirical difficulties that creationism has.  He also seems to take the most honest position that a creationist can take, which is that he is not a creationist because of the empirical data, he is a creationist because of his literalist reading of the Bible. (Kurt Wise is another one of these rarities.)

Anyway, we have Todd Wood to thank for putting BDIST online.

The original bdist.pl can be found at:
http://web.archive.org/web/20060616052853/http://www.bryancore.org/bsg/bdist.html
http://web.archive.org/web/20050905145521/http://www.bryancore.org/bsg/bdist.pl

The online version of bdistmds (no source code, but the key calculations seem to be identical with bdist.pl) can be found at:
http://www.coresci.org/bdist.html

And limited documentation at:
http://www.coresci.org/bdistinfo.html

The original creationist paper on baraminic distances was: 

Robinson, D. Ashley; Cavanaugh, David P. (1998). A Quantitative Approach to Baraminology With Examples from the Catarrhine Primates. Creation Research Society Quarterly, 34(4),  196-208.
