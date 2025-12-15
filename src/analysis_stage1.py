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
        }
        self.output_dir = "../test/Zqq_vertex/"
        self.output_dir_eos = "/eos/user/m/mdefranc/aleph_vertex/Zqq/1994_stage1/{}_xy_flip".format(datetime.now().strftime("%Y_%m_%d"))
        self.input_dir = "/eos/experiment/aleph/EDM4HEP/MC/1994/"
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

        #df = df.Define("TrackStateFlipped", f"AlephSelection::flipD0_copy( {coll['TrackState']} )")
        df = df.Filter("AlephSelection::sel_class_filter(16)(ClassBitset) ")

        #used_collection = "TrackStateFlipped"
        used_collection = coll["TrackState"]

        # Vertex
        df = df.Filter("Vertices.size() == 1")  # to remove eventually
        df = df.Define("pv", "TLorentzVector(Vertices[0].position.x, Vertices[0].position.y, Vertices[0].position.z, 0.0)",)
        df = df.Define("pv_true", "AlephSelection::get_EventPrimaryVertexP4()({})".format(coll["GenParticles"]))

        df = df.Define("pv_true_x", "pv_true.X()")
        df = df.Define("pv_true_y", "pv_true.Y()")
        df = df.Define("pv_true_z", "pv_true.Z()")

        df = df.Define("V0_x", "Vertices[0].position.x")
        df = df.Define("V0_y", "Vertices[0].position.y")
        df = df.Define("V0_z", "Vertices[0].position.z")

        df = df.Define("V0_res_x", "V0_x - pv_true_x")
        df = df.Define("V0_res_y", "V0_y - pv_true_y")
        df = df.Define("V0_res_z", "V0_z - pv_true_z")

        df = df.Define("V0_ratio_x", "V0_x / pv_true_x")
        df = df.Define("V0_ratio_y", "V0_y / pv_true_y")
        df = df.Define("V0_ratio_z", "V0_z / pv_true_z")

        df = df.Define("ntracks", "ReconstructedParticle2Track::getTK_n({})".format(used_collection))

        #df = df.Define("VertexObject_allTracks_noBS", "VertexFitterSimple::VertexFitter_Tk(1, {}, false)".format(coll["TrackState"]))
        #df = df.Define("Vertex_allTracks_noBS", "VertexingUtils::get_VertexData(VertexObject_allTracks_noBS)")

        #df = df.Define("VertexObject_allTracks", "VertexFitterSimple::VertexFitter_Tk(1, {}, true, 12.5,.5,720.,0.,0.,0.)".format(used_collection)) # 10um as unit
        #df = df.Define("Vertex_allTracks", "VertexingUtils::get_VertexData(VertexObject_allTracks)")

        #df = df.Define("PV_allTracks", "ROOT::VecOps::RVec<edm4hep::VertexData> v; v.push_back(Vertex_allTracks); return v;")


        df = df.Define("RecoedPrimaryTracks", "VertexFitterSimple::get_PrimaryTracks({}, true, 12.5,.5,720.,0.,0.,0.)".format(used_collection)) # 10um as unit
        df = df.Define("VertexObject", "VertexFitterSimple::VertexFitter_Tk(1, RecoedPrimaryTracks, true, 12.5,.5,720.,0.,0.,0.)") # 10um as unit
        df = df.Define("Vertex", "VertexingUtils::get_VertexData(VertexObject)")
        df = df.Define("SecondaryTracks", "VertexFitterSimple::get_NonPrimaryTracks({}, RecoedPrimaryTracks)".format(used_collection))

        df = df.Define("n_RecoedPrimaryTracks", "ReconstructedParticle2Track::getTK_n(RecoedPrimaryTracks)")
        df = df.Define("n_SecondaryTracks", "ReconstructedParticle2Track::getTK_n(SecondaryTracks)")

        df = df.Define("PV_x", "Vertex.position.x")
        df = df.Define("PV_y", "Vertex.position.y")
        df = df.Define("PV_z", "Vertex.position.z")

        df = df.Define("PV_res_true_x", "PV_x - pv_true_x")
        df = df.Define("PV_res_true_y", "PV_y - pv_true_y")
        df = df.Define("PV_res_true_z", "PV_z - pv_true_z")

        df = df.Define("PV_ratio_true_x", "PV_x / pv_true_x")
        df = df.Define("PV_ratio_true_y", "PV_y / pv_true_y")
        df = df.Define("PV_ratio_true_z", "PV_z / pv_true_z")

        df = df.Define("PV_res_V0_x", "PV_x - V0_x")
        df = df.Define("PV_res_V0_y", "PV_y - V0_y")
        df = df.Define("PV_res_V0_z", "PV_z - V0_z")

        df = df.Define("PV_ratio_V0_x", "PV_x / V0_x")
        df = df.Define("PV_ratio_V0_y", "PV_y / V0_y")
        df = df.Define("PV_ratio_V0_z", "PV_z / V0_z")

        #df = df.Define("PV", "ROOT::VecOps::RVec<edm4hep::VertexData> v; v.push_back(Vertex); return v;")

        # Return the dataframe for further chaining or final output
        return df

    def output(self):
        # List of branches to save in output
        return [
            "event_type", "ntracks",
            "V0_x", "V0_y", "V0_z",
            "pv_true_x", "pv_true_y", "pv_true_z",
            "V0_res_x", "V0_res_y", "V0_res_z",
            "V0_ratio_x", "V0_ratio_y", "V0_ratio_z",
            "n_RecoedPrimaryTracks", "n_SecondaryTracks",
            "PV_x", "PV_y", "PV_z",
            "PV_res_true_x", "PV_res_true_y", "PV_res_true_z",
            "PV_ratio_true_x", "PV_ratio_true_y", "PV_ratio_true_z",
            "PV_res_V0_x", "PV_res_V0_y", "PV_res_V0_z",
            "PV_ratio_V0_x", "PV_ratio_V0_y", "PV_ratio_V0_z",
        ]
