#!/usr/bin/env python3

import os
import sys
import glob
import math
import ROOT

# helper function to read branches
# note: several branches are nested RVec<RVec<...>> (per jet, per vertex), so a single
#       [0] does not reach the scalars. convert the whole thing to plain python lists:
#       comparing two RVecs with != is element-wise and returns an RVec, which is always
#       truthy when non-empty, i.e. every such branch would be reported as different.

def get_value(obj):
    if hasattr(obj, "__len__") and not isinstance(obj, (str, bytes)):
        return [get_value(x) for x in obj]
    return obj


# helper function to compare two values read by get_value.
# floats are compared with a relative tolerance: the two codes do the same arithmetic in a
# different order, so last-bit differences are expected and are not physics.

def values_equal(val1, val2, rtol=1e-4):
    if isinstance(val1, list) != isinstance(val2, list):
        return False
    if isinstance(val1, list):
        if len(val1) != len(val2):
            return False
        return all(values_equal(a, b, rtol) for a, b in zip(val1, val2))
    if isinstance(val1, float) or isinstance(val2, float):
        if math.isnan(val1) and math.isnan(val2):
            return True
        return abs(val1 - val2) <= rtol * max(1.0, abs(val1), abs(val2))
    return val1 == val2


# helper to get a tree from either a single .root file or a directory of chunk_*.root
# (the batch productions are chunked, e.g. <tag>/Zbb/QQB/chunk_*.root, so a plain TFile
# is not enough). Returns the tree; the caller must keep the returned owner alive, else
# ROOT garbage-collects the file/chain underneath the tree.
def open_tree(path, tree_name):
    if os.path.isdir(path):
        chain = ROOT.TChain(tree_name)
        chunks = sorted(glob.glob(os.path.join(path, "**", "chunk_*.root"), recursive=True))
        if not chunks:
            raise RuntimeError(f"No chunk_*.root found under directory: {path}")
        for c in chunks:
            chain.Add(c)
        print(f"  chained {len(chunks)} chunks from {path}")
        return chain, chain

    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open file: {path}")
    tree = f.Get(tree_name)
    if not tree:
        raise RuntimeError(f"Cannot find tree '{tree_name}' in {path}")
    return f, tree


# helper function which returns the unique run & event_number combinations in a file as python set
def load_events(filename, tree_name, run_branch, event_branch, do_flavour_filter = False, flavour_val = 5):
    owner, tree = open_tree(filename, tree_name)   # keep `owner` referenced until we are done

    events = set()

    for entry in tree:
        run_vec = getattr(entry, run_branch)
        evt_vec = getattr(entry, event_branch)

        # Safety checks
        if len(run_vec) == 0 or len(evt_vec) == 0:
            continue

        if len(run_vec) != 1 or len(evt_vec) != 1:
            raise RuntimeError(
                f"Expected exactly 1 element per event, "
                f"got {len(run_vec)} (run) and {len(evt_vec)} (event)"
            )
        
        # luka's files aren't filtered by flavour, ours are, add the filter here:
        if do_flavour_filter:
            evt_type = getattr(entry, "genEventType")

            if evt_type != flavour_val:
                continue

        run = int(run_vec[0])
        evt = int(evt_vec[0])

        events.add((int(run), int(evt)))

    del owner   # release the TFile / TChain
    return events

# function that compares the common events branch by branch
# translation of branch names needs to be passed 

def compare_events(filepath1, filepath2,
                   common_events,
                   branch_map,
                   max_print=10, tree_name = "events"):

    #debug:
    print("file1:", filepath1)
    print("file2:", filepath2)

    # open again - either a single file or a directory of chunks
    # (owner1/owner2 must stay in scope, they keep the TFile/TChain alive)
    owner1, tree1 = open_tree(filepath1, tree_name)
    owner2, tree2 = open_tree(filepath2, tree_name)

    # Build ROOT internal indices
    print("Building indices...")
    tree1.BuildIndex(branch_map["run"][0], branch_map["event"][0])
    tree2.BuildIndex(branch_map["run"][1], branch_map["event"][1])

    n_diff_events = 0
    # count in how many events each branch differs, so the summary shows where the
    # disagreement actually sits instead of only which branches ever differed
    n_diff_per_branch = {name: 0 for name in branch_map}

    for run_number, event_number in common_events:
        entry_number_file1 = tree1.GetEntryNumberWithIndex(run_number, event_number)
        tree1.GetEntry(entry_number_file1)

        entry_number_file2 = tree2.GetEntryNumberWithIndex(run_number, event_number)
        tree2.GetEntry(entry_number_file2)

        event_has_diff = False

        for branch_name, (branch_1, branch_2) in branch_map.items():

            val1 = get_value(getattr(tree1, branch_1))
            val2 = get_value(getattr(tree2, branch_2))

            if values_equal(val1, val2):
                continue

            n_diff_per_branch[branch_name] += 1

            if not event_has_diff:
                event_has_diff = True
                n_diff_events += 1
                if n_diff_events <= max_print:
                    print(f"\nDifference in event {run_number} {event_number}")

            if n_diff_events <= max_print:
                print(f"  Branch mismatch: {branch_1} vs {branch_2} in event {event_number}")
                print(f"    file1 (Luka): {val1}")
                print(f"    file2 (ours): {val2}")

    n_common = len(common_events)
    print(f"\nTotal events with differences: {n_diff_events} / {n_common}")
    print("\nEvents differing per branch:")
    for branch_name, n_diff in sorted(n_diff_per_branch.items(), key=lambda kv: -kv[1]):
        frac = 100. * n_diff / n_common if n_common else 0.
        print(f"  {branch_name:24s} {n_diff:6d}  ({frac:5.1f}%)")


def main():

    tree_name = "events"

    # First we get the sets of run and event number from each file and find the overlap (=common events)
    file1 = "/eos/user/l/llambrec/aleph-data/ntuples-withks/eventlevel/mc/output_qqb_1.root" # for training? 
    # file1 = "/eos/user/l/llambrec/aleph-data/ntuples-withksloose/eventlevel/mc/output_qqb_1.root" #for V0 plots, no pointing angle cut!
    run_branch_1 = "runNumber"
    event_branch_1 = "eventNumber"

    print("Checking events in file {}".format(file1))
    # Note: Luka's file is not filtered by flavour, ours is, so apply the flavour filter here
    events1 = load_events(file1, tree_name, run_branch_1, event_branch_1, do_flavour_filter = True, flavour_val = 5.)
    print(f"File1: {len(events1)} events")

    # file2 = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/1994/zqq/stage1/v09_ntuple_valid/ntuple_valid_tester_5.root" 
    # file2 = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/1994/zqq/stage1/v15_SV_test/Zbb.root" 
    # either a single .root file, or a directory of chunks (batch productions are chunked,
    # the directory is chained automatically - see open_tree)
    file2 = "/eos/experiment/fcc/ee/analyses/case-studies/aleph/processedMC/1994/zqq/stage1/v17_SVsFix_w_data/Zbb"
    run_branch_2 = "run_number"
    event_branch_2 = "event_number"

    print("Checking events in file {}".format(file2))
    events2 = load_events(file2, tree_name, run_branch_2, event_branch_2, do_flavour_filter = False, flavour_val = 5.)
    print(f"File2: {len(events2)} events")

    only_in_1 = events1 - events2
    only_in_2 = events2 - events1
    common = events1 & events2

    print("\nComparison results:")
    print(f"Common events: {len(common)}")
    print(f"Only in file1: {len(only_in_1)}") # This should be 0 !!!
    print(f"Only in file2: {len(only_in_2)}") # Note: out file processed larger fraction of input, so will be  > 0


    # print(common)

    # now check the common events branch by branch

    # names are not unified, so need translation map:
    branch_translation = {
        #"name":"(name_file1, name_file2)",
        "run":("runNumber", "run_number"),
        "event":("eventNumber", "event_number"),
        # input track collections (what goes INTO the primary vertex fit).
        # These should agree in every event - if they don't, the disagreement is upstream of
        # the vertex fit and nothing below is meaningful.
        "n_tracks_all":("Event_nTracks", "n_tracks_all"),
        "n_selected_tracks":("Event_nSelectedTracks", "n_tracks_sel"),
        # note: our n_tracks_sel_vertexfit / n_trackstates_sel (tracks passing the |D0|,|Z0|
        # preselection, i.e. the actual fit input) have no equivalent in Luka's ntuples -
        # he applies that cut inside getPrimaryTracks and never stores the intermediate count.
        # output of primary vertex fit
        "n_primary_tracks":("Event_nPrimaryTracks", "n_primary_tracks"),
        "n_secondary_tracks":("Event_nSecondaryTracks", "n_secondary_tracks"),
        # # vertex position
        "vertex_x":("PV_x", "Vertex_refit_x"),
        "vertex_y":("PV_y", "Vertex_refit_y"),
        "vertex_z":("PV_z", "Vertex_refit_z"),
        # #truth vertex
        # "gen_vertex_x":("GenPV_x", "gen_vertex_x"),
        # "gen_vertex_y":("GenPV_y", "gen_vertex_y"),
        # "gen_vertex_z":("GenPV_z", "gen_vertex_z"),
        #secondary vertices:
        "n_sv_event":("Event_nSV", "n_sv_event"),
        "n_sv_jets":("Jets_nSV", "n_sv_jets"),
        "sv_ntracks":("SecondaryVertices_nTracks", "sv_ntracks"),
        "n_v0_event":("Event_nV0Candidates", "n_v0_event"),
        "v0_invM":("V0Candidates_mass", "v0_invM"),
        # 



    }

    compare_events(file1, file2, common, branch_translation)
    print(f"Reminder: Common events: {len(common)}")



if __name__ == "__main__":
    main()
