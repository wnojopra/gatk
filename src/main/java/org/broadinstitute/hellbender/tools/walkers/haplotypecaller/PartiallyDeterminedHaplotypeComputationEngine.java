package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.Tuple;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.ArrayUtils;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.utils.haplotype.EventMap;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.haplotype.PartiallyDeterminedHaplotype;
import org.broadinstitute.hellbender.utils.read.CigarBuilder;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Class that manages the complicated steps involved in generating artifical haplotypes for the PDHMM:
 *
 * The primoary method ton this class is {@link #generatePDHaplotypes} that gets called with an existing assembly result set
 * as well as the list of alleles that were found from the pileupcaller code and it attempts to replace the existing haplotypes
 * from the provided result set object with new haplotyeps that were artificially generated from merging the assembled alleles
 * with those found by the pileup caller. The method might hit one of a few heuristic barriers and chose to fallback on only
 * the assembled haplotypes if the processing would become too complicated.
 */
public class PartiallyDeterminedHaplotypeComputationEngine {
    final static int MAX_PD_HAPS_TO_GENERATE = 256*2;; //(2048 is illuminas #) (without optimizing the hmm to some degree this is probably unattainable)
    final static int MAX_BRANCH_PD_HAPS = 128; //(128 is illuminas #)
    final static int MAX_VAR_IN_EVENT_GROUP = 17; // (20 is illuminas #)

    //To make this somewhat cleaner of a port from Illumina, we have two base spaces. R and U space. R space is vcf coordinate space,
    //U is a 0->N (N = region size) based space where Insertions are +0.5 and deletions are + 1 from their original position

    // We use this comparator for haplotype construcion to make it slightly easier to build compound haplotypes (i.e. snp+insertion/deletion at the same anchor base)
    public static final Comparator<VariantContext> HAPLOTYPE_SNP_FIRST_COMPARATOR = Comparator.comparingInt(VariantContext::getStart)
            // Decide arbitrarily so as not to accidentally throw away overlapping variants
            .thenComparingInt(vc -> vc.getReference().length())
            .thenComparingInt(vc -> vc.getAlternateAllele(0).length())
            .thenComparing(vc -> vc.getAlternateAllele(0));


    /**
     * The workhorse method for the PDHMM branching code. This method is responsible for the lionshare of the work in taking a list of
     * assembly haplotypes/alleles and the alleles found in the ColumnwiseDetection code and converting them into the artificial haplotypes
     * that are necessary for variant calling. This code might fail due to having too complex a site, if that happens it will return
     * the input sourceSet object without anything being changed.
     *
     * Overall pattern is as follows:
     *  1. Generate the combined list of variants from assembly and pileup-caller that pass the heuristic filers
     *  2. Generate lists of variants that overlap in groups
     *  3. Use SmithWaterman to realign pairs of 2/3 events that contain an indel, recording pairs that realign in a better
     *     set of sub-events or result in the reference.
     *  4. Merge event groups that have equivalent events.
     *  5. Iterate over every start position and attempt to build determined variants:
     *      5a. Select an allele at a position (ref or one of the alts)
     *      5b. Ask each event group for branches of mutually exclusive events (overlapping or from the equivalent event code) for the determined variants.
     *      5c. Merge all event groups into distinct branches (combinatorially)
     *      5d. For each branch, build a PDHaplotype with the determined allele and the allowed variants for that branch.
     *
     * @param sourceSet AssemblyResultSet to be modified with the new haplotypes
     * @param callingSpan Calling span to subset determined events to (to handle padding)
     * @param referenceHaplotype Reference haplotype against which to build artifical haps
     * @param assemblyVariants Assembly variants.
     * @param pileupAllelesFoundShouldFilter Pileup alleles that should be filtered if they are part of the assembly
     * @param pileupAllelesPassingFilters Pileup alleles that pass the heuristics to be included in genotyping
     * @param snpAdjacentToIndelLimit If pileup allele snps are too close to assembled indels we thorw them out.
     * @param aligner SmithWatermanAligner to use for filtering out equivalent event sets
     * @param swParameters Parameters for hap-hap comparisons
     * @param makeDeterminedHapsInstead If true, generate regular haplotypes to be handed to the PairHMM (note this will result in a lot more failures due to combinitorial expansion)
     * @param debugSite If true, print out all the details of how this code is functioning for debugging.
     * @return unchanged assembly result set if failed, updated haplotyeps otherwise
     */
    public static AssemblyResultSet generatePDHaplotypes(final AssemblyResultSet sourceSet,
                                                         final Locatable callingSpan,
                                                         final Haplotype referenceHaplotype,
                                                         final SortedSet<VariantContext> assemblyVariants,
                                                         final List<VariantContext> pileupAllelesFoundShouldFilter,
                                                         final List<VariantContext> pileupAllelesPassingFilters,
                                                         final int snpAdjacentToIndelLimit,
                                                         final SmithWatermanAligner aligner,
                                                         final SWParameters swParameters,
                                                         final boolean makeDeterminedHapsInstead,
                                                         final boolean debugSite) {

        final TreeSet<VariantContext> variantsInOrder = new TreeSet<>(
                HAPLOTYPE_SNP_FIRST_COMPARATOR);

        // First we filter the assembly variants based on badness from the graph
        for (VariantContext delVariant : pileupAllelesFoundShouldFilter) {

            List<VariantContext> variantsToRemove = assemblyVariants.stream().filter(
                    v -> v.getStart() == delVariant.getStart() &&
                            delVariant.getReference().equals(v.getReference()) &&
                            delVariant.getAlternateAllele(0).equals(v.getAlternateAllele(0))).collect(Collectors.toList());

            if (!variantsToRemove.isEmpty()) {
                if (debugSite) System.out.println("Removing assembly variants due to columnwise heurisits: " + variantsToRemove);
                variantsToRemove.forEach(assemblyVariants::remove);
            }
        }

        // Ignore any snps from pileups that were close to indels
        final List<VariantContext> givenAllelesFiltered = pileupAllelesPassingFilters.stream()
                .filter(vc -> vc.isIndel() ||
                        assemblyVariants.stream().filter(VariantContext::isIndel).noneMatch(indel -> vc.withinDistanceOf(indel, snpAdjacentToIndelLimit))).collect(Collectors.toList());

        // TODO make sure this TREE-SET is properly comparing the VCs
        variantsInOrder.addAll(assemblyVariants);
        variantsInOrder.addAll(givenAllelesFiltered);

        if (debugSite) System.out.println("Variants to PDHapDetermination:\n"+
                variantsInOrder.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("\n")));

        // TODO this is where we filter out if indels > 32
        // TODO this is where we do eventTestSet work
        List<VariantContext> vcsAsList = new ArrayList<>(variantsInOrder);

        // NOTE: we iterate over this several times and expect it to be sorted.
        Map<Double, List<VariantContext>> eventsByDRAGENCoordinates = new LinkedHashMap<>();
        SortedMap<Integer, List<VariantContext>> variantsByStartPos = new TreeMap<>();
        List<EventGroup> eventGroups = new ArrayList<>();
        int lastEventEnd = -1;
        for (VariantContext vc : variantsInOrder) {
            // Break everything into independent groups (don't worry about transitivitiy right now)
            Double eventKey = vc.getStart() + (vc.isSimpleInsertion()? 0.5:0) + (vc.isSimpleDeletion()? 1 : 0) - referenceHaplotype.getStartPosition();
            eventsByDRAGENCoordinates.putIfAbsent(eventKey, new ArrayList<>());
            eventsByDRAGENCoordinates.get(eventKey).add(vc);
            variantsByStartPos.putIfAbsent(vc.getStart(), new ArrayList<>());
            variantsByStartPos.get(vc.getStart()).add(vc);
            if (debugSite) System.out.println("testing:"+vc);
            if (debugSite) System.out.println("EventKey:"+eventKey);
            if (eventKey <= lastEventEnd + 0.5) {
                eventGroups.get(eventGroups.size()-1).addEvent(vc);
            } else {
                eventGroups.add(new EventGroup(vc));
            }
            int newEnd = (int) (vc.getEnd() - referenceHaplotype.getStartPosition());
            if (debugSite) System.out.println("LastEventEnd:"+lastEventEnd+"    newEventEnd:"+newEnd);
            lastEventEnd = Math.max(newEnd, lastEventEnd);
        }
        //Print the event groups
        if (debugSite) eventsByDRAGENCoordinates.entrySet().stream().map(e -> {
            return String.format("%.1f", e.getKey()) + " -> " + e.getValue().stream()
                    .map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition()))
                    .collect(Collectors.joining(","));
        }).forEach(System.out::println);

        // Iterate over all events starting with all indels

        List<List<VariantContext>> dissalowedPairs = smithWatermanRealignPairsOfVariantsForEquivalentEvents(referenceHaplotype, aligner, swParameters, debugSite, variantsInOrder, vcsAsList);
        if (debugSite) {
            System.out.println("Dissalowed Variant pairs:");
            dissalowedPairs.stream().map(l -> l.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->"))).forEach(System.out::println);
        }

        if (debugSite) {
            System.out.println("Event groups before merging:\n"+eventGroups.stream().map(eg -> eg.toDisplayString((int)referenceHaplotype.getStartPosition())).collect(Collectors.joining("\n")));
        }
        //Now that we have the dissalowed groups, lets merge any of them from seperate groups:
        //TODO this is not an efficient way of doing this
        for (List<VariantContext> pair : dissalowedPairs) {
            EventGroup eventGrpLeft = null;
            for (VariantContext event : pair) {
                EventGroup grpForEvent = eventGroups.stream().filter(grp -> grp.contains(event)).findFirst().get();
                // If the event isn't in the same event group as its predicessor, merge this group with that one and
                if (eventGrpLeft != grpForEvent) {
                    if (eventGrpLeft == null) {
                        eventGrpLeft = grpForEvent;
                    } else {
                        eventGrpLeft.mergeEvent(grpForEvent);
                        eventGroups.remove(grpForEvent);
                    }
                }
            }
        }
        if (debugSite) {
            System.out.println("Event groups after merging:\n"+eventGroups.stream().map(eg -> eg.toDisplayString((int)referenceHaplotype.getStartPosition())).collect(Collectors.joining("\n")));
        }

        //Now we have finished with the work of merging event groups transitively by position and mutually exclusiveness. Now every group should be entirely independant of one another:
        if (eventGroups.stream().map(eg -> eg.populateBitset(dissalowedPairs)).anyMatch(b->!b)) {
            // if any of our event groups is too large, abort.
            if (debugSite ) System.out.println("Found event group with too many variants! Aborting haplotype building");
            return sourceSet;
        };

        Set<Haplotype> outputHaplotypes = new LinkedHashSet<>();
        if (makeDeterminedHapsInstead) {
            // only add the reference haplotype if we are producing regular haplotype objects (not PD haplotypes for the haplotype alleles)
            outputHaplotypes.add(referenceHaplotype);
        }

        //Iterate over very VCF start position in R space
        List<Map.Entry<Integer, List<VariantContext>>> entriesRInOrder = new ArrayList<>(variantsByStartPos.entrySet());
        /**
         * Overall Loop:
         * Iterate over every cluster of variants with the same start position.
         */
        for (int indexOfDeterminedInR = 0; indexOfDeterminedInR < entriesRInOrder.size(); indexOfDeterminedInR++) {
            Map.Entry<Integer, List<VariantContext>> variantSiteGroup = entriesRInOrder.get(indexOfDeterminedInR);
            if (debugSite) System.out.println("working with variants of the group: " + variantSiteGroup);
            // Skip
            if (entriesRInOrder.get(indexOfDeterminedInR).getKey() < callingSpan.getStart() || entriesRInOrder.get(indexOfDeterminedInR).getKey() > callingSpan.getEnd()) {
                if (debugSite) System.out.println("Skipping determined hap construction! otside of span: "+callingSpan);
                continue;
            }

            final List<VariantContext> determinedVariants = variantSiteGroup.getValue();

            /**
             * Determined Event Loop:
             * We iterate over every ref position and select single alleles (including ref) from that reference position (in R space) to be "determined"
             *
             * NOTE: we skip the reference allele in the event that we are making determined haplotypes instead of undetermined haplotypes
             */
            for (int IndexOfAllele = (makeDeterminedHapsInstead?0:-1); IndexOfAllele < determinedVariants.size(); IndexOfAllele++) { //note -1 for I here corresponds to the reference allele at this site
                if (debugSite) System.out.println("Working with allele at site: "+(IndexOfAllele ==-1? "[ref:"+(variantSiteGroup.getKey()-referenceHaplotype.getStartPosition())+"]" : PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int)referenceHaplotype.getStartPosition()).apply(determinedVariants.get(IndexOfAllele))));
                // This corresponds to the DRAGEN code for
                // 0 0
                // 0 1
                // 1 0
                final boolean isRef = IndexOfAllele == -1;
                final VariantContext determinedEventToTest = determinedVariants.get(isRef ? 0 : IndexOfAllele);

                /*
                 * Here we handle any of the necessary work to deal with the event groups and maybe forming compund branches out of the groups
                 */
                //Set Determined pairs:
                List<Tuple<VariantContext, Boolean>> determinedPairs = new ArrayList<>();
                for(int j = 0; j < determinedVariants.size(); j++) {
                    determinedPairs.add(new Tuple<>(determinedVariants.get(j), IndexOfAllele == j));
                }

                // Loop over eventGroups, have each of them return a list of VariantContexts
                List<Set<VariantContext>> branchExcludeAlleles = new ArrayList<>();
                branchExcludeAlleles.add(new HashSet<>()); // Add the null branch (assuming no exclusions)

                /* Note for future posterity:
                 * An assembly region could potentially have any number of (within some limitations) of event groups. When we are constructing
                 * haplotypes out of the assembled variants we want to take the dot product of the branches for each set of event groups that
                 * we find. I.E. if an event group with mutex variants (B,C) requires two branches for Variant A and Variant A also leads to two branches in
                 * another event group with mutex variants (D,E). Then we would want to ultimately generate the branches A,B,D -> A,B,E -> A,C,D -> A,C,E.
                 * This is why we iterate over branchExcludeAlleles internally here.
                 */
                for(EventGroup group : eventGroups ) {
                    if (group.causesBranching()) {
                        List<List<Tuple<VariantContext, Boolean>>> groupVCs = group.getVariantGroupsForEvent(determinedPairs, true);
                        // Combinatorially expand the branches as necessary
                        List<Set<VariantContext>> newBranchesToAdd = new ArrayList<>();
                        for (Set<VariantContext> excludedVars : branchExcludeAlleles) {
                            //For every exclude group, fork it by each subset we have:
                            for (int i = 1; i < groupVCs.size(); i++) { //NOTE: iterate starting at 1 here because we special case that branch at the end
                                Set<VariantContext> newSet = new HashSet<>(excludedVars);
                                groupVCs.get(i).stream().filter(t -> !t.b).forEach(t -> newSet.add(t.a));
                                newBranchesToAdd.add(newSet);
                            }
                            // Be careful since this event group might have returned nothing
                            if (!groupVCs.isEmpty()) {
                                groupVCs.get(0).stream().filter(t -> !t.b).forEach(t -> excludedVars.add(t.a));
                            }
                        }
                        branchExcludeAlleles.addAll(newBranchesToAdd);

                        if (branchExcludeAlleles.size() > MAX_BRANCH_PD_HAPS) {
                            if (debugSite ) System.out.println("Found too many branches for variants at: "+determinedEventToTest.getStart()+" aborting and falling back to Assembly Varinats!");
                            return sourceSet;
                        }
                    }
                }

                if (debugSite) {
                    System.out.println("Branches:");
                    for (int i = 0; i < branchExcludeAlleles.size(); i++) {
                        final int ifinal = i;
                        System.out.println("Branch "+i+" VCs:");
                        System.out.println("exclude:"+branchExcludeAlleles.get(i).stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int)referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                        //to match dragen debug output for personal sanity
                        System.out.println("include:"+variantsInOrder.stream().filter(variantContext -> !branchExcludeAlleles.get(ifinal).contains(variantContext)).map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int)referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                    }
                }


                /**
                 * Now handle each branch independently of the others. (the logic is the same in every case except we must ensure that none of the excluded alleles get included when constructing haps.
                 */
                for (Set<VariantContext> excludeEvents : branchExcludeAlleles) {

                    List<Haplotype> branchHaps = new ArrayList<>();
                    // iterate over all of the
                    List<VariantContext> newBranch = new ArrayList<>();

                    // Handle the simple case of making PD haplotypes
                    if (!makeDeterminedHapsInstead) {
                        for (int secondRIndex = 0; secondRIndex < entriesRInOrder.size(); secondRIndex++) {
                            if (secondRIndex != indexOfDeterminedInR) {
                                // We know here that nothing illegally overlaps because there are no groups.
                                // Also exclude any events that overlap the determined allele since we cant construct them (also this stops compound alleles from being formed)
                                // NOTE: it is important that we allow reference alleles to overlap undetermined variants as it leads to mismatches against DRAGEN otherwise.
                                entriesRInOrder.get(secondRIndex).getValue().stream()
                                        .filter(vc -> !excludeEvents.contains(vc))
                                        .forEach(newBranch::add);
                            } else {
                                newBranch.add(determinedEventToTest);
                            }
                        }
                        newBranch.sort(HAPLOTYPE_SNP_FIRST_COMPARATOR);
                        PartiallyDeterminedHaplotype newPDHaplotypeFromEvents = createNewPDHaplotypeFromEvents(referenceHaplotype, determinedEventToTest, isRef, newBranch);
                        newPDHaplotypeFromEvents.setAllDeterminedEventsAtThisSite(determinedVariants); // accounting for determined variants for later in case we are in optimization mode
                        branchHaps.add(newPDHaplotypeFromEvents);

                    } else {
                        // TODO currently this approach doesn't properly handle a bunch of duplicate events...
                        // If we are producing determined bases, then we want to enforce that every new event at least has THIS event as a variant.
                        List<List<VariantContext>> variantGroupsCombinatorialExpansion = new ArrayList<>();
                        variantGroupsCombinatorialExpansion.add(new ArrayList<>());
                        // We can drastically cut down on combinatorial expansion here by saying each allele is the FIRST variant in the list, thus preventing double counting.
                        for (int secondRIndex = indexOfDeterminedInR; secondRIndex < entriesRInOrder.size(); secondRIndex++) {
                            if (variantGroupsCombinatorialExpansion.size() > MAX_BRANCH_PD_HAPS) {
                                if(debugSite ) System.out.println("Too many branch haplotypes ["+variantGroupsCombinatorialExpansion.size()+"] generated from site, falling back on assebmly variants!");
                                return sourceSet;
                            }
                            // Iterate through the growing combinatorial expansion of haps, split it into either having or not having the variant.
                            if (secondRIndex == indexOfDeterminedInR) {
                                for (List<VariantContext> hclist : variantGroupsCombinatorialExpansion) {
                                    hclist.add(determinedEventToTest);
                                }
                            // Othewise make sure to include the combinatorial expansion of events at the other site
                            } else {
                                List<List<VariantContext>> hapsPerVCsAtRSite = new ArrayList<>();
                                for (VariantContext vc : entriesRInOrder.get(secondRIndex).getValue()) {
                                    for (List<VariantContext> hclist : variantGroupsCombinatorialExpansion) {
                                        if (!excludeEvents.contains(vc)) {
                                            List<VariantContext> newList = new ArrayList<>(hclist);
                                            newList.add(vc);
                                            hapsPerVCsAtRSite.add(newList);
                                        }
                                    }
                                }
                                //Add them after to prevent accidentally adding duplicates of events at a site
                                variantGroupsCombinatorialExpansion.addAll(hapsPerVCsAtRSite);
                            }
                        }

                        for (List<VariantContext> subset : variantGroupsCombinatorialExpansion) {
                            subset.sort(HAPLOTYPE_SNP_FIRST_COMPARATOR);
                            if (debugSite) System.out.println("Construcing Hap From Events:"+ subset.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                            branchHaps.add(constructHaplotypeFromVariants(referenceHaplotype, subset, true));
                        }
                    }
                    // Add the branch haps to the results:
                    if (debugSite) {
                        System.out.println("Constructed Haps for Branch"+excludeEvents.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining(",")) + ":");
                        System.out.println(branchHaps.stream().map(h -> h.getCigar() + " " + h.toString()).collect(Collectors.joining("\n")));
                    }

                    outputHaplotypes.addAll(branchHaps);
                    if (outputHaplotypes.size() > MAX_PD_HAPS_TO_GENERATE) {
                        if (debugSite) System.out.println("Too many Haps ["+outputHaplotypes.size()+"] generated at this site! Aborting!");
                        return sourceSet;
                    }
                }
            }
        }

        if (outputHaplotypes.size() > MAX_PD_HAPS_TO_GENERATE) {
            if (debugSite ) System.out.println("Too many branch haplotypes found, aborting ["+outputHaplotypes.size()+"]");
            return sourceSet;
        }
        sourceSet.storeAssemblyHaplotypes();
        Set<Haplotype> result = new LinkedHashSet<>();
        // TODO this is an entirely unnecessary step that can be done away with but i leave in because it makes debugging against previous versions much easier.
        outputHaplotypes.stream().sorted(new Comparator<Haplotype>() {
            @Override
            public int compare(Haplotype o1, Haplotype o2) {
                return new String(o1.getBases()).compareTo(new String(o2.getBases()));
            }
        }).forEach(h -> result.add(h));
        sourceSet.replaceAllHaplotypes(result);
        if (debugSite) System.out.println("Constructed Haps for Branch"+sourceSet.getHaplotypeList().stream().map(Haplotype::toString).collect(Collectors.joining("\n")));
        if (!makeDeterminedHapsInstead) {
            // Setting a boolean on the source-set to indicate to downstream code that we have PD haplotypes
            sourceSet.setPartiallyDeterminedMode();
        }
        if (debugSite ) System.out.println("Returning "+outputHaplotypes.size()+" to the HMM");
        return sourceSet;
    }

    /**
     * Helper method that handles one of the Heuristics baked by DRAGEN into this artifical haplotype genration code.
     *
     * To help mitigate the risk of generating combinatorial haplotypes with SNPs/Indels that that might or might not add
     * up to equivalent events, DRAGEN enforces that events are NOT allowed to be present in the same haplotype if they
     * (when run through smith waterman) add up to other events that were found at this assembly region.
     *
     * To cut down on the complexity of the task; we (and DRAGEN) follow this procedure:
     * 1. look at all sets of 2 variants where at least one is an indel and none overlap.
     *    a) for each set construct an artifical haplotype with only those two variants on it
     *    b) SmithWtarman align it against the reference to generate the cheapest cigar string representation
     *    c) Construct the event map for the new artificial haplotype, if any events in the new event map are in our list of variants
     *       but are NOT the constituent events that were used to construct the haplotype then dissalow the pair
     * 2. Look at all sets of 3 variants that do not contain disallowed pairs found in step 1.
     *    a-b-c) repeat steps 1a,1b,and 1c on the 3 evetn sets
     *
     * @return A list of lists of variant contexts that correspond to dissalowed groups. This list may be empty if none are found.
     */
    private static List<List<VariantContext>> smithWatermanRealignPairsOfVariantsForEquivalentEvents(Haplotype referenceHaplotype, SmithWatermanAligner aligner, SWParameters swParameters, boolean debugSite, TreeSet<VariantContext> variantsInOrder, List<VariantContext> vcsAsList) {
        List<List<VariantContext>> dissalowedPairs = new ArrayList<>();

        //Iterate over all 2 element permutations in which one element is an indel and test for alignments
        for (int i = 0; i < vcsAsList.size(); i++) {
            final VariantContext firstEvent = vcsAsList.get(i);
            if (firstEvent.isIndel()) {
                // For every indel, make every 2-3 element subset (without overlapping) of variants to test for equivalency
                for (int j = 0; j < vcsAsList.size(); j++) {
                    final VariantContext secondEvent = vcsAsList.get(j);
                    // Don't compare myslef, any overlappers to myself, or indels i've already examined me (to prevent double counting)
                    if (j != i && !eventsOverlapForPDHapsCode(firstEvent, secondEvent, true) && ((!secondEvent.isIndel()) || j > i)) {
                        final List<VariantContext> events = new ArrayList<>(Arrays.asList(firstEvent, secondEvent));
                        events.sort(HAPLOTYPE_SNP_FIRST_COMPARATOR);
                        if (debugSite) System.out.println("Testing events: "+ events.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                        if (constructArtificialHaplotypeAndTestEquivalentEvents(referenceHaplotype, aligner, swParameters, variantsInOrder, events, debugSite)) {
                            dissalowedPairs.add(events);
                        }
                    }
                }
            }
        }

        //TODO NOTE: there are some discrepancies with the iteration over 3x variants in some complicated cases involving
        //TODO       lots of transitively dissalowed pairs. Hopefully this is a minor effect.
        //Now iterate over all 3 element pairs and make sure none of the
        for (int i = 0; i < vcsAsList.size(); i++) {
            final VariantContext firstEvent = vcsAsList.get(i);
            if (firstEvent.isIndel()) {
                // For every indel, make every 2-3 element subset (without overlapping) of variants to test for equivalency
                for (int j = 0; j < vcsAsList.size(); j++) {
                    final VariantContext secondEvent = vcsAsList.get(j);
                    // Don't compare myslef, any overlappers to myself, or indels i've already examined me (to prevent double counting)
                    if (j != i && !eventsOverlapForPDHapsCode(firstEvent, secondEvent, true) && ((!secondEvent.isIndel()) || j > i)) {
                        // if i and j area lready disalowed keep going
                        if (dissalowedPairs.stream().anyMatch(p -> p.contains(firstEvent) && p.contains(secondEvent))) {
                            continue;
                        }
                        final List<VariantContext> events = new ArrayList<>(Arrays.asList(firstEvent, secondEvent));
                        // If our 2 element arrays weren't inequivalent, test subsets of 3 including this:
                        for (int k = j+1; k < vcsAsList.size(); k++) {
                            final VariantContext thirdEvent = vcsAsList.get(k);
                            if (k != i && !eventsOverlapForPDHapsCode(thirdEvent, firstEvent, true) && !eventsOverlapForPDHapsCode(thirdEvent, secondEvent, true)) {
                                // if k and j or k and i are dissalowed, keep looking
                                if (dissalowedPairs.stream().anyMatch(p -> (p.contains(firstEvent) && p.contains(thirdEvent)) || (p.contains(secondEvent) && p.contains(thirdEvent)))) {
                                    continue;
                                }
                                List<VariantContext> subList = new ArrayList<>(events);
                                subList.add(thirdEvent);
                                subList.sort(HAPLOTYPE_SNP_FIRST_COMPARATOR);
                                if (debugSite) System.out.println("Testing events: " + subList.stream().map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition())).collect(Collectors.joining("->")));
                                if (constructArtificialHaplotypeAndTestEquivalentEvents(referenceHaplotype, aligner, swParameters, variantsInOrder, subList, debugSite)) {
                                    dissalowedPairs.add(subList);
                                }
                            }
                        }
                    }
                }
            }
        }

        return dissalowedPairs;
    }

    /**
     * Overlaps method to handle indels and snps correctly. Specifically for this branching codes purposes, SNPS never overlap other SNPS,
     * indels don't overlap on their anchor bases and insertions don't overlap anything except deletions spanning them or other insertions
     * at the same base.
     *
     * @param snpsOverlap if false, don't ever evaluate snps as overlapping other snps (we do this because sometimes we need to construct artifical haps where we don't allow overlapping)
     */
    static boolean eventsOverlapForPDHapsCode(VariantContext vc1, VariantContext vc2, boolean snpsOverlap){
        if (!snpsOverlap && vc2.isSNP() && vc1.isSNP()) {
            return false;
        }
        if (!vc1.getContig().equals(vc2.getContig())) {
            return false;
        }
        double vc1start = vc1.isIndel() ? (vc1.isSimpleDeletion() ? vc1.getStart() + 1 : vc1.getStart() + 0.5) : vc1.getStart();
        double vc1end = vc1.isSimpleInsertion() ? vc1.getEnd() + 0.5 : vc1.getEnd();
        double vc2start = vc2.isIndel() ? (vc2.isSimpleDeletion() ? vc2.getStart() + 1 : vc2.getStart() + 0.5) : vc2.getStart();
        double vc2end = vc2.isSimpleInsertion() ? vc2.getEnd() + 0.5 : vc2.getEnd();

        //Pulled directly from CoordMath.java (not using here because of doubles)
        return (vc2start >= vc1start && vc2start <= vc1end) || (vc2end >=vc1start && vc2end <= vc1end) || vc1start >= vc2start && vc1end <= vc2end;
    }


    /**
     * This method is the helper that manages the actual SW alignment and testing of a group of variants vs the reference haplotype.
     *
     * The method is as follows, construct and artificial haplotype of the provided events, then realign it vs the reference and test
     * if any of the resulting variants are present in the inputs (but doesn't match)
     *
     * @return
     */
    @VisibleForTesting
    private static boolean constructArtificialHaplotypeAndTestEquivalentEvents(Haplotype referenceHaplotype, SmithWatermanAligner aligner, SWParameters swParameters, TreeSet<VariantContext> vcs, List<VariantContext> eventsToTest, boolean debugSite) {
        final Haplotype realignHap = constructHaplotypeFromVariants(referenceHaplotype, eventsToTest, false);
        //Special case to capture events that equal the reference (and thus have empty event maps).
        if (Arrays.equals(realignHap.getBases(), referenceHaplotype.getBases())) {
            if (debugSite) System.out.println("Events add up to the reference! Dissalowing pair");
            return true;
        }
        //ALIGN!
        realignHap.setCigar(CigarUtils.calculateCigar(referenceHaplotype.getBases(), realignHap.getBases(), aligner, swParameters, SWOverhangStrategy.INDEL));
        EventMap.buildEventMapsForHaplotypes(Collections.singletonList(realignHap), referenceHaplotype.getBases(), referenceHaplotype.getGenomeLocation(), false,0);
        //TODO this isn't exactly how it works in DRAGEN, or at least we've lost the score here...
        final boolean wasEquivalentEvent = realignHap.getEventMap().getVariantContexts().stream().filter(eMapVC ->
                // Are there any variants NOT in our initial list
                eventsToTest.stream().noneMatch(v -> {
                    return doVariantsMatch(eMapVC, v);
                }))
                // Do any of those variants appear in our overall list of alleles
                .anyMatch(eMapVc -> vcs.stream().anyMatch(v -> {
                    return doVariantsMatch(eMapVc, v);
                }));
        if (debugSite) System.out.println(
                realignHap.getEventMap().getVariantContexts().stream()
                .map(PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString((int) referenceHaplotype.getStartPosition()))
                .collect(Collectors.joining("\n")));
        if (wasEquivalentEvent) {
            if (debugSite) System.out.println("Events mismatched!");
        }

        return wasEquivalentEvent;
    }

    // A helper method to assert that variant contexts in the event map match those outside of it.
    private static boolean doVariantsMatch(VariantContext eMapVC, VariantContext v) {
        return eMapVC.getStart() == v.getStart() &&
                eMapVC.getReference().equals(v.getReference()) &&
                eMapVC.getAlternateAllele(0).equals(v.getAlternateAllele(0)) &&
                eMapVC.getAlternateAlleles().size() == v.getAlternateAlleles().size();
    }


    /**
     *
     *
     * NOTE: this accepts multiple alleles stacked up at the same base (assuming the order is SNP -> INDEL)
     * @param refHap
     * @param variantContexts
     * @return
     */
    @VisibleForTesting
    public static Haplotype constructHaplotypeFromVariants(final Haplotype refHap, final List<VariantContext> variantContexts, final boolean setEventMap) {
        //ASSERT that the base is ref and cool
        if (!refHap.isReference() || refHap.getCigar().numCigarElements() > 1) {
            throw new RuntimeException("This is not a valid base haplotype for construction");
        }
        final long genomicStartPosition = refHap.getStartPosition();
        long refOffsetOfNextBaseToAdd = genomicStartPosition;

        byte[] refbases = refHap.getBases();
        CigarBuilder runningCigar = new CigarBuilder();
        byte[] newRefBases = {};

        //ASSUME sorted for now
        // use the reverse list to save myself figuring out cigars for right now
        for (VariantContext vc : variantContexts) {
            if (vc.getAlternateAlleles().size() > 1) {
                throw new RuntimeException("too may alt alleles");
            }
            Allele refAllele = vc.getReference();
            Allele altAllele = vc.getAlternateAllele(0);

            int intermediateRefStartPosition = (int) (refOffsetOfNextBaseToAdd - genomicStartPosition);
            int intermediateRefEndPosition = Math.toIntExact(vc.getStart() - genomicStartPosition);

            if ((vc.isIndel() && intermediateRefEndPosition - intermediateRefStartPosition < -1) || (!vc.isIndel() && intermediateRefEndPosition - intermediateRefStartPosition < 0)) {//todo clean this up
                throw new RuntimeException("Variant "+vc+" is out of order in the PD event list: "+variantContexts);
            }
            if (intermediateRefEndPosition - intermediateRefStartPosition > 0) { // Append the cigar element for the anchor base if necessary.
                runningCigar.add(new CigarElement(intermediateRefEndPosition - intermediateRefStartPosition, CigarOperator.M));
            }
            // Include the ref base for indel if the base immediately proceeding this event is not already tracked
            boolean includeRefBaseForIndel = vc.isIndel() && (intermediateRefStartPosition <= intermediateRefEndPosition);

            CigarElement newCigarElement;
            if (refAllele.length() == altAllele.length()) {
                newCigarElement = new CigarElement(refAllele.length(), CigarOperator.X);
            } else {
                // If the last base was filled by another event, don't attempt to fill in the indel ref base.
                if (includeRefBaseForIndel) {
                    runningCigar.add(new CigarElement(1, CigarOperator.M)); //When we add an indel we end up inserting a matching base
                }
                newCigarElement = new CigarElement(Math.abs(altAllele.length() - refAllele.length()),
                        refAllele.length() > altAllele.length() ?
                                CigarOperator.D : CigarOperator.I);
            }
            runningCigar.add(newCigarElement);

            if (intermediateRefEndPosition - intermediateRefStartPosition > 0) {
                newRefBases = ArrayUtils.addAll(newRefBases, ArrayUtils.subarray(refbases, intermediateRefStartPosition, (int) (vc.getStart() - genomicStartPosition))); // bases before the variant
            }
            // Handle the ref base for indels that exlcude their ref bases
            if (refAllele.length() != altAllele.length() && !includeRefBaseForIndel) {
                newRefBases = ArrayUtils.addAll(newRefBases, Arrays.copyOfRange(altAllele.getBases(),1, altAllele.length()));
            // else add the snp
            } else {
                newRefBases = ArrayUtils.addAll(newRefBases, altAllele.getBases()); // refbases added
            }
            refOffsetOfNextBaseToAdd = vc.getEnd() + 1; //TODO this is probably not set for future reference
        }

        // Finish off the haplotype with the final bases
        int refStartIndex = (int) (refOffsetOfNextBaseToAdd - genomicStartPosition);
        newRefBases = ArrayUtils.addAll(newRefBases, ArrayUtils.subarray(refbases, refStartIndex, refbases.length));
        runningCigar.add(new CigarElement(refbases.length - refStartIndex, CigarOperator.M));

        final Haplotype outHaplotype = new Haplotype(newRefBases, false, refHap.getGenomeLocation(), runningCigar.make());
        if (setEventMap) {
            EventMap.buildEventMapsForHaplotypes(Collections.singletonList(outHaplotype), refHap.getBases(), refHap.getGenomeLocation(), false,0);
            // NOTE: we set this AFTER generating the event maps because hte event map code above is being generated from the ref hap so this offset will cause out of bounds errors
            outHaplotype.setAlignmentStartHapwrtRef(refHap.getAlignmentStartHapwrtRef()); //TODO better logic here
        }
        return outHaplotype;
    }

    /**
     * Construct a PD haplotype from scratch
     *
     *
     * Generally we are constructing a new haplotype with all of the reference bases for SNP events and with the longest possible allele for INDEL events.
     * For deleteions, we extend the haplotype by the ref lenght
     *
     * NOTE: we assume each provided VC
     *
     * //TODO these both need guardrails against exceptions
     *
     * @param base
     * @param eventWithVariant
     * @param useRef
     * @param constituentEvents
     */
    @VisibleForTesting
    static PartiallyDeterminedHaplotype createNewPDHaplotypeFromEvents(final Haplotype base, final VariantContext eventWithVariant, final boolean useRef, final List<VariantContext> constituentEvents) {

        //ASSERT that the base is ref and cool
        if (!base.isReference() || base.getCigar().numCigarElements() > 1) {
            throw new RuntimeException("This is not a valid base haplotype for construction");
        }
        long genomicStartPosition = base.getStartPosition();
        long refOffsetOfNextBaseToAdd = genomicStartPosition;

        byte[] refBasesToAddTo = base.getBases();
        CigarBuilder runningCigar = new CigarBuilder(false); // NOTE: in some incredibly rare edge cases involving the legacy assembly region trimmer a deletion can hang past the edge of an active window.
        byte[] newHaplotypeBasees = {};
        byte[] pdBytes = {};

        //ASSUME sorted for now
        // use the reverse list to save myself figuring out cigars for right now
        for (VariantContext vc : constituentEvents) {
            int intermediateRefStartPosition = (int) (refOffsetOfNextBaseToAdd - genomicStartPosition);
            int intermediateRefEndPosition = Math.toIntExact(vc.getStart() - genomicStartPosition);

            // An extra special case if we are a SNP following a SNP
            if (vc.isSNP() && intermediateRefEndPosition - intermediateRefStartPosition == -1 && ((pdBytes[pdBytes.length-1] & PartiallyDeterminedHaplotype.SNP) != 0)  ) {
                byte[] array = PartiallyDeterminedHaplotype.getPDBytesForHaplotypes(vc.getReference(), vc.getAlternateAllele(0));
                pdBytes[pdBytes.length-1] = (byte) (pdBytes[pdBytes.length-1] | array[0]); // adding any partial bases if necessary
                continue;
            }

            // Ref alleles (even if they overlap undetermined events) should be skipped
            if (vc.getStart()==eventWithVariant.getStart() && useRef) {
                continue;
            }

            //Check that we are allowed to add this event (and if not we are)
            if ((vc.isIndel() && intermediateRefEndPosition - intermediateRefStartPosition < -1) || (!vc.isIndel() && intermediateRefEndPosition - intermediateRefStartPosition < 0)) {//todo clean this up
                throw new RuntimeException("Variant "+vc+" is out of order in the PD event list: "+constituentEvents);
            }

            // Add the cigar for bases we skip over
            if (intermediateRefEndPosition - intermediateRefStartPosition > 0) {
                runningCigar.add(new CigarElement(intermediateRefEndPosition - intermediateRefStartPosition, CigarOperator.M));
            }

            // Include the ref base for indel if the base immediately proceeding this event is not already tracked
            boolean includeRefBaseForIndel = vc.isIndel() && (intermediateRefStartPosition <= intermediateRefEndPosition);

            // Determine the alleles to add
            Allele refAllele = vc.getReference();
            Allele altAllele = vc.getAlternateAllele(0);
            boolean isInsertion = altAllele.length() > refAllele.length(); // If its an insertion we flip to "ADD" the bases to the ref.
            boolean isEvent = false;
            byte[] basesToAdd;
            // If this is the blessed variant, add
            if (vc.getStart()==eventWithVariant.getStart()) {
                if (!useRef && eventWithVariant.getAlternateAlleles().size() > 1) {
                    throw new RuntimeException("the Blessed variant must be monoallelic");
                }
                isEvent = true;
                basesToAdd = useRef? refAllele.getBases() : altAllele.getBases();
            // Otherwise make sure we are adding the longest allele (for indels) or the ref allele for snps.
            } else {
                basesToAdd = refAllele.length() >= altAllele.length() ? refAllele.getBases() : altAllele.getBases();
            }

            // Remove anchor base if necessary
            if (vc.isIndel() && !includeRefBaseForIndel) {
                basesToAdd = Arrays.copyOfRange(basesToAdd, 1, basesToAdd.length);
            }

            // Figure out the cigar to add:
            // - If we are in the ref, simply add the cigar corresponding to the allele we are using
            // -
            CigarElement newCigarElement;
            // if this is the event special case
            if (isEvent) {
                if (vc.isSNP()) {
                    newCigarElement = new CigarElement(refAllele.length(), useRef? CigarOperator.M : CigarOperator.X);
                } else {
                    if (vc.isIndel() && includeRefBaseForIndel) {
                        runningCigar.add(new CigarElement( 1, CigarOperator.M));
                    }
                    // For Insertions: mark the Cigar as I if we aren't in ref
                    if (isInsertion) {
                        newCigarElement = new CigarElement(useRef ? 0 : Math.max(refAllele.length(), altAllele.length()) - 1, CigarOperator.I);
                    // For Deletions: Always include the bases, however mark them as M or D accordingly
                    } else {
                        newCigarElement = new CigarElement(Math.max(refAllele.length(), altAllele.length()) - 1, useRef ? CigarOperator.M : CigarOperator.D);
                    }
                }
           // If we aren't in the blessed variant, add a match and make sure the array is set accordingly
            } else {
                if (!vc.isIndel()) {
                    newCigarElement = new CigarElement(refAllele.length() , CigarOperator.M);
                } else {
                    // Maybe add the cigar for the anchor base
                    if (includeRefBaseForIndel) {
                        runningCigar.add(new CigarElement(1, CigarOperator.M));
                    }
                    // Add the cigar for the indel allele bases
                    if (isInsertion) {
                        // Insertions stay in the cigar since they are added relative to the reference
                        newCigarElement = new CigarElement(Math.abs(altAllele.length() - refAllele.length()), CigarOperator.I);
                    } else {
                        // Deletions become matches because they still exist as bases on the reference
                        newCigarElement = new CigarElement(Math.abs(altAllele.length() - refAllele.length()), CigarOperator.M);
                    }
                }
            }
            runningCigar.add(newCigarElement);

            // Add ref basses up to this if necessary
            if (intermediateRefEndPosition - intermediateRefStartPosition > 0) {
                newHaplotypeBasees = ArrayUtils.addAll(newHaplotypeBasees, ArrayUtils.subarray(refBasesToAddTo, intermediateRefStartPosition, (int) (vc.getStart() - genomicStartPosition))); // bases before the variant
                pdBytes = ArrayUtils.addAll(pdBytes, new byte[vc.getStart() - (int)refOffsetOfNextBaseToAdd]); // bases before the variant
            }
            newHaplotypeBasees = ArrayUtils.addAll(newHaplotypeBasees, basesToAdd); // refbases added
            if (includeRefBaseForIndel) {
                pdBytes = ArrayUtils.add(pdBytes, (byte)0);
            }
            pdBytes = ArrayUtils.addAll(pdBytes, isEvent?
                    new byte[basesToAdd.length - (includeRefBaseForIndel?1:0)] :
                    PartiallyDeterminedHaplotype.getPDBytesForHaplotypes(isInsertion?
                            altAllele :
                            refAllele,
                            isInsertion? refAllele :
                                    altAllele)); // refbases added
            refOffsetOfNextBaseToAdd = vc.getEnd() + 1; //TODO this is probably not set for future reference
        }

        // Finish off the haplotype with the final bases
        int refStartIndex = (int) (refOffsetOfNextBaseToAdd - genomicStartPosition);
        newHaplotypeBasees = ArrayUtils.addAll(newHaplotypeBasees, ArrayUtils.subarray(refBasesToAddTo, refStartIndex, refBasesToAddTo.length));
        pdBytes = ArrayUtils.addAll(pdBytes, new byte[refBasesToAddTo.length - refStartIndex]);
        runningCigar.add(new CigarElement(refBasesToAddTo.length - refStartIndex, CigarOperator.M));

        return new PartiallyDeterminedHaplotype(
                new Haplotype(newHaplotypeBasees, false, base.getGenomeLocation(), runningCigar.make()),
                useRef,
                pdBytes,
                constituentEvents,
                eventWithVariant,
                runningCigar.make(),
                eventWithVariant.getStart(),
                base.getAlignmentStartHapwrtRef());
    }

    // A helper class for managing mutually exclusive event clusters and the logic arround forming valid events vs eachother.
    private static class EventGroup {
        List<VariantContext> variantsInBitmapOrder;
        HashSet<VariantContext> variantContextSet;
        //From Illumina (there is a LOT of math that will eventually go into these)/
        BitSet allowedEvents = null;

        // Optimizaiton to save ourselves recomputing the subsets at every point its necessary to do so.
        List<List<Tuple<VariantContext,Boolean>>> cachedEventLitsts = null;

        public EventGroup(final VariantContext variantContext) {
            variantsInBitmapOrder = new ArrayList<>();
            variantContextSet = new HashSet<>();
            variantsInBitmapOrder.add(variantContext);
            variantContextSet.add(variantContext);
        }
        public EventGroup() {
            variantsInBitmapOrder = new ArrayList<>();
            variantContextSet = new HashSet<>();
        }

        /**
         * This is the primary method for handling mutually exclusive events in this subgroup. This code amd methods comes directly from DRAGEN:
         *
         * Create a #Variants bitset to store valid pairings:
         *      - The index of each element corresponds to an enumerated subset of alleles in this group
         *      - Each bit in the index corresponds to the presence or absence of a variant from the vcf list.
         *          - For example with variants [A,B,C] the number 5 corresponds to subset [A,C]
         *      - A false in the bitset corresponds to a dissalowed pair.
         *      - NOTE: we can use 32bit ints for these bitshift operations by virtue of the fact that we limit ourselves to at most 22 variants per group.
         * Iterate through pairs of Variants that overlap and mark off any pairings including this.
         * Iterate through the mutex variants and ensure pairs containing all mutex variant groups are marked as true
         *
         * @param dissalowedEvents Pairs of events dissalowed
         * @return false if the event group is too large to process
         */
        public boolean populateBitset(List<List<VariantContext>> dissalowedEvents) {
            if (variantsInBitmapOrder.size() > MAX_VAR_IN_EVENT_GROUP) {
                return false;
            }
            if (variantsInBitmapOrder.size() < 2) {
                return true;
            }

            allowedEvents = new BitSet(variantsInBitmapOrder.size());
            allowedEvents.flip(1, 1 << variantsInBitmapOrder.size());
            // initialize all events as being allowed and then dissalow them in turn .

            // Ensure the list is in positional order before commencing.
            variantsInBitmapOrder.sort(HAPLOTYPE_SNP_FIRST_COMPARATOR);
            List<Integer> bitmasks = new ArrayList<>();
            // Mark as disallowed all events that overlap eachother
            for (int i = 0; i < variantsInBitmapOrder.size(); i++) {
                for (int j = i+1; j < variantsInBitmapOrder.size(); j++) {
                    if (eventsOverlapForPDHapsCode(variantsInBitmapOrder.get(i), variantsInBitmapOrder.get(j), false)) {
                        bitmasks.add(1 << i | 1 << j);
                    }
                }
            }
            // mark as dissalowed any sets of variants from the bitmask.
            for (List<VariantContext> disallowed : dissalowedEvents) {
                //
                if (disallowed.stream().anyMatch(v -> variantContextSet.contains(v))){
                    int bitmask = 0;
                    for (VariantContext v : disallowed) {
                        int indexOfV = variantsInBitmapOrder.indexOf(v);
                        if (indexOfV < 0) {
                            throw new RuntimeException("Something went wrong in event group merging, variant "+v+" is missing from the event group despite being in a mutex pair: "+disallowed+"\n"+this);
                        }
                        bitmask += 1 << variantsInBitmapOrder.indexOf(v);
                    }
                    bitmasks.add(bitmask);
                }
            }

            // Now iterate through the list and dissalow all events with every bitmask
            if (!bitmasks.isEmpty()) {
                events:
                for (int i = 1; i < allowedEvents.length(); i++) {
                    for (final int mask : bitmasks) {
                        if ((i & mask) == mask) { // are the bits form the mask true?
                            allowedEvents.set(i, false);
                            continue events;
                            // Once i is set we don't need to keep checking bitmasks
                        }
                    }
                }
            }

            return true;
        }

        /**
         * This method handles the logic involved in getting all of the allowed subsets of alleles for this event group.
         *
         * @param dissalowSubsets
         * @return
         */
        public List<List<Tuple<VariantContext,Boolean>>> getVariantGroupsForEvent(final List<Tuple<VariantContext, Boolean>> eventsForMask, final boolean dissalowSubsets) {
            // If we are dealing with an external to this list event
            int eventMask = 0;
            int maskValues = 0;
            for(Tuple<VariantContext, Boolean> event : eventsForMask) {
                if (variantContextSet.contains(event.a)) {
                    int index = variantsInBitmapOrder.indexOf(event.a);
                    eventMask = eventMask | (1 << index);
                    maskValues = maskValues | ((event.b ? 1 : 0) << index);
                }
            }
            // Special case (if we are determining bases outside of this mutex cluster we can reuse the work from previous iterations)
            if (eventMask == 0 && cachedEventLitsts != null) {
                return cachedEventLitsts;
            }

            List<Integer> ints = new ArrayList<>();
            // Iterate from the BACK of the list (i.e. ~supersets -> subsets)
            // NOTE: we skip over 0 here since that corresponds to ref-only events, handle those externally to this code
            outerLoop:
            for (int i = allowedEvents.length(); i > 0; i--) {
                // If the event is allowed AND if we are looking for a particular event to be present or absent.
                if (allowedEvents.get(i) && (eventMask == 0 || ((i & eventMask) == maskValues))) {
                    // Only check for subsets if we need to
                    if (dissalowSubsets) {
                        for (Integer group : ints) {
                            // if the current i is a subset of an existing group
                            if ((i | group) == group) {
                                continue outerLoop;
                            }
                        }
                    }
                    ints.add(i);
                }
            }

            // Now that we have all the mutex groups, unpack them into lists of variants
            List<List<Tuple<VariantContext,Boolean>>> output = new ArrayList<>();
            for (Integer grp : ints) {
                List<Tuple<VariantContext,Boolean>> newGrp = new ArrayList<>();
                for (int i = 0; i < variantsInBitmapOrder.size(); i++) {
                    // if the corresponding bit is 1, set it as such, otherwise set it as 0.
                    newGrp.add(new Tuple<>(variantsInBitmapOrder.get(i), ((1<<i) & grp) != 0));
                }
                output.add(newGrp);
            }
            // Cache the result
            if(eventMask==0) {
                cachedEventLitsts = Collections.unmodifiableList(output);
            }
            return output;
        }

        public boolean causesBranching() {
            return variantsInBitmapOrder.size() > 1;
        }

        //Print The event group in Illumina indexed ordering:
        public String toDisplayString(int startPos) {
            return "EventGroup: " + variantsInBitmapOrder.stream().map(vc -> PartiallyDeterminedHaplotype.getDRAGENDebugVariantContextString(startPos).apply(vc)).collect(Collectors.joining("->"));
        }

        public boolean contains(final VariantContext event) {
            return variantContextSet.contains(event);
        }

        public void addEvent(final VariantContext event) {
            variantsInBitmapOrder.add(event);
            variantContextSet.add(event);
            allowedEvents = null;
        }

        public EventGroup mergeEvent(final EventGroup other) {
            variantsInBitmapOrder.addAll(other.variantsInBitmapOrder);
            variantContextSet.addAll(other.variantsInBitmapOrder);
            allowedEvents = null;
            return this;
        }
    }
}
