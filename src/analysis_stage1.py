from argparse import ArgumentParser
from datetime import datetime

class Analysis():
    def __init__(self, cmdline_args):
        parser = ArgumentParser(
        description='Additional analysis arguments',
        usage='Provide additional arguments after analysis script path')

        self.ana_args, _ = parser.parse_known_args(cmdline_args['remaining'])

        # Define config and sample list here if needed
        self.process_list = {
            "QQB": {"fraction": 1, "output": "Zqq", "chunks": 70}
            # "QQB": {"fraction": 1/1000., "output": "Zqq"}
        }
        self.output_dir = "../test/Zqq_vertex/"
        self.output_dir_eos = "/eos/user/m/mdefranc/aleph_vertex/Zqq/1994_stage1/checkSV"
        self.input_dir = "/eos/experiment/aleph/EDM4HEP/MC/1994/"
        self.test_file = "/eos/experiment/aleph/EDM4HEP/MC/1994/QQB/ZM4212_39_AL.root"
        self.comp_group = 'group_u_FCC.local_gen'

        # self.batch_queue = 'workday'
        self.batch_queue = 'tomorrow'
        #self.batch_queue = 'testmatch'
        
        self.n_threads = 4
        self.include_paths = ["analyzer.h"]

    def analyzers(self, df):
        coll = {
            "GenParticles": "MCParticles",
            "PFParticles": "RecoParticles",
            "PFTracks": "EFlowTrack",
            "PFPhotons": "EFlowPhoton",
            "PFNeutralHadrons": "EFlowNeutralHadron",
            "TrackState": "_Tracks_trackStates",
            "TrackerHits": "TrackerHits",
            "CalorimeterHits": "CalorimeterHits",
            "PathLength": "EFlowTrack_L",
            "Bz": "magFieldBz",
        }

        # MC EVENT FILTERING
        df = df.Define("event_type", "AlephSelection::get_EventType({})".format(coll["GenParticles"]))
        # df = df.Filter("event_type == 2") # example filter

        df = df.Define("TrackStateFlipped", f"AlephSelection::flipD0_copy( {coll['TrackState']} )")
        df = df.Filter("AlephSelection::sel_class_filter(16)(ClassBitset) ")

        #used_collection = "TrackStateFlipped"
        used_collection = coll["TrackState"]

        # Vertex
        df = df.Filter("Vertices.size() > 0")  # to remove eventually
        df = df.Define("pv", "TLorentzVector(Vertices[0].position.x, Vertices[0].position.y, Vertices[0].position.z, 0.0)",)
        df = df.Define("pv_true", "AlephSelection::get_EventPrimaryVertexP4()({})".format(coll["GenParticles"]))

        df = df.Define("ntracks", "ReconstructedParticle2Track::getTK_n({})".format(used_collection))

        # values from Matteo based on sim-only
        res_x = 150. # in um
        res_y = 5. # in um
        res_z = 1. # in cm

        df = df.Define("RecoedPrimaryTracks", "VertexFitterSimple::get_PrimaryTracks({}, true, {},{},{},0.,0.,0.)".format(used_collection, res_x/10., res_y/10., res_z*1E03)) # 10um as unit (x,y), 1cm as unit (z)
        df = df.Define("VertexObject", "VertexFitterSimple::VertexFitter_Tk(1, RecoedPrimaryTracks, true, {},{},{},0.,0.,0.)".format(res_x/10., res_y/10., res_z*1E03)) # 10um as unit (x,y), 1cm as unit (z)
        df = df.Define("Vertex", "VertexingUtils::get_VertexData(VertexObject)")
        df = df.Define("SecondaryTracks", "VertexFitterSimple::get_NonPrimaryTracks({}, RecoedPrimaryTracks)".format(used_collection))

        df = df.Define("n_RecoedPrimaryTracks", "ReconstructedParticle2Track::getTK_n(RecoedPrimaryTracks)")
        df = df.Define("n_SecondaryTracks", "ReconstructedParticle2Track::getTK_n(SecondaryTracks)")

        # values from Luka taking into account performance on data
        res_x_loose = 200. # in um
        res_y_loose = 100. # in um
        res_z_loose = 2. # in cm

        df = df.Define("RecoedPrimaryTracks_looseBS", "VertexFitterSimple::get_PrimaryTracks({}, true, {},{},{},0.,0.,0.)".format(used_collection, res_x_loose/10., res_y_loose/10., res_z_loose*1E03)) # 10um as unit (x,y), 1cm as unit (z)
        df = df.Define("n_RecoedPrimaryTracks_looseBS", "ReconstructedParticle2Track::getTK_n(RecoedPrimaryTracks_looseBS)")
        df = df.Define("SecondaryTracks_looseBS", "VertexFitterSimple::get_NonPrimaryTracks({}, RecoedPrimaryTracks_looseBS)".format(used_collection))
        df = df.Define("n_SecondaryTracks_looseBS", "ReconstructedParticle2Track::getTK_n(SecondaryTracks_looseBS)")



        #df = df.Define("PV", "ROOT::VecOps::RVec<edm4hep::VertexData> v; v.push_back(Vertex); return v;")

        # Return the dataframe for further chaining or final output
        return df

    def output(self):
        # List of branches to save in output
        return [
                "event_type", "ntracks",
                "n_RecoedPrimaryTracks", "n_SecondaryTracks",
                "n_RecoedPrimaryTracks_looseBS", "n_SecondaryTracks_looseBS",
        ]
