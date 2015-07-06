#! /bin/sh
#run combiner.C to make the plots for all variables

root -l -b << EOF
.L combiner.C

combiner("mgg");
combiner("t1pfmet");

combiner("r91");
combiner("r92");
combiner("hoe1");
combiner("hoe2");
combiner("sieie1");
combiner("sieie2");
combiner("chiso1");
combiner("chiso2");
combiner("neuiso1");
combiner("neuiso2");
combiner("phoiso1");
combiner("phoiso2");

.q

EOF

