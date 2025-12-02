

processList = {
"QQB" : {"fraction" : 1, "output": "Zqq"},
}



outputDir = "/eos/user/m/mdefranc/aleph_vertex/Zqq/1994_stage1/"
inputDir = "/eos/experiment/aleph/EDM4HEP/MC/1994/"
nCPUS = 4


includePaths = ["analyzer.h"]


class RDFanalysis:
    def analysers(df):
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
        #df = df.Filter("event_type == 2") # d-quark: 1, u-quark:2, s-quark:3, c-quark:4, b-quark: 5
        ########################################################################################################################


        df = df.Define("TrackStateFlipped",f"AlephSelection::flipD0_copy( {coll['TrackState']} )")
        df = df.Filter("AlephSelection::sel_class_filter(16)(ClassBitset) ")


        # Define RP kinematics
        ####################################################################################################
        df = df.Define("RP_px", "ReconstructedParticle::get_px(RecoParticles)")
        df = df.Define("RP_py", "ReconstructedParticle::get_py(RecoParticles)")
        df = df.Define("RP_pz", "ReconstructedParticle::get_pz(RecoParticles)")
        df = df.Define("RP_e", "ReconstructedParticle::get_e(RecoParticles)")
        df = df.Define("RP_m", "ReconstructedParticle::get_mass(RecoParticles)")

        # Define pseudo-jets
        ####################################################################################################
        df = df.Define("pjetc", "JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")
        # Anti-kt clustering and jet constituents
        ####################################################################################################
        df = df.Define("_jet", "JetClustering::clustering_ee_kt(2, 2, 1, 0)(pjetc)")
        df = df.Define("jets","JetClusteringUtils::get_pseudoJets(_jet)" )
        df = df.Define("_jetc", "JetClusteringUtils::get_constituents(_jet)") 
        df = df.Define("jetc", "JetConstituentsUtils::build_constituents_cluster(RecoParticles, _jetc)")
        df = df.Define("jetConstitutentsTypes", f"AlephSelection::build_constituents_Types()(ParticleID, _jetc)")


        ############################################# Event Level Variables #######################################################
        df = df.Define("jet_p4", "JetConstituentsUtils::compute_tlv_jets(jets)" )
        df = df.Define("event_invariant_mass", "JetConstituentsUtils::InvariantMass(jet_p4[0], jet_p4[1])")
        ################################################### MC Filtering ############################################################

        # ===== VERTEX
        df = df.Filter("Vertices.size() == 1") # to remove eventually
        df = df.Define(
            "pv", "TLorentzVector(Vertices[0].position.x, Vertices[0].position.y, Vertices[0].position.z, 0.0)",
        )
        df = df.Define("pv_true", "AlephSelection::get_EventPrimaryVertexP4()({})".format(coll["GenParticles"]))

        df = df.Define("pv_x", "Vertices[0].position.x")
        df = df.Define("pv_y", "Vertices[0].position.y")
        df = df.Define("pv_z", "Vertices[0].position.z")



        df = df.Define("pv_true_x", "pv_true.X()")
        df = df.Define("pv_true_y", "pv_true.Y()")
        df = df.Define("pv_true_z", "pv_true.Z()")
        df = df.Define("pv_res_x", "pv_x - pv_true_x")
        df = df.Define("pv_res_y", "pv_y - pv_true_y")
        df = df.Define("pv_res_z", "pv_z - pv_true_z")
        


        df = df.Define("ntracks", "ReconstructedParticle2Track::getTK_n({})".format(coll["TrackState"]))

        # fit vertex with all tracks (first step)
        df = df.Define("VertexObject_allTracks_noBS", "VertexFitterSimple::VertexFitter_Tk(1, {}, false)".format(coll["TrackState"]))
        df = df.Define("Vertex_allTracks_noBS", "VertexingUtils::get_VertexData(VertexObject_allTracks_noBS)")


        df = df.Define("VertexObject_allTracks", "VertexFitterSimple::VertexFitter_Tk(1, {}, true, 150e-3,50e-3,10.,0.,0.,0.)".format(coll["TrackState"]))
        df = df.Define("Vertex_allTracks", "VertexingUtils::get_VertexData(VertexObject_allTracks)")
        

        df = df.Define("PV_allTracks_noBS_x", "Vertex_allTracks_noBS.position.x")
        df = df.Define("PV_allTracks_noBS_y", "Vertex_allTracks_noBS.position.y")
        df = df.Define("PV_allTracks_noBS_z", "Vertex_allTracks_noBS.position.z")

        df = df.Define("PV_res_allTracks_noBS_true_x", "PV_allTracks_noBS_x - pv_true_x")
        df = df.Define("PV_res_allTracks_noBS_true_y", "PV_allTracks_noBS_y - pv_true_y")
        df = df.Define("PV_res_allTracks_noBS_true_z", "PV_allTracks_noBS_z - pv_true_z")

        df = df.Define("PV_ratio_allTracks_noBS_true_x", "PV_allTracks_noBS_x / pv_true_x")
        df = df.Define("PV_ratio_allTracks_noBS_true_y", "PV_allTracks_noBS_y / pv_true_y")
        df = df.Define("PV_ratio_allTracks_noBS_true_z", "PV_allTracks_noBS_z / pv_true_z")

        df = df.Define("PV_res_allTracks_noBS_V0_x", "PV_allTracks_noBS_x - pv_x")
        df = df.Define("PV_res_allTracks_noBS_V0_y", "PV_allTracks_noBS_y - pv_y")
        df = df.Define("PV_res_allTracks_noBS_V0_z", "PV_allTracks_noBS_z - pv_z")

        df = df.Define("PV_ratio_allTracks_noBS_V0_x", "PV_allTracks_noBS_x / pv_x")
        df = df.Define("PV_ratio_allTracks_noBS_V0_y", "PV_allTracks_noBS_y / pv_y")
        df = df.Define("PV_ratio_allTracks_noBS_V0_z", "PV_allTracks_noBS_z / pv_z")

        df = df.Define("PV_allTracks_noBS", "ROOT::VecOps::RVec<edm4hep::VertexData> v; v.push_back(Vertex_allTracks_noBS); return v;")


        df = df.Define("PV_allTracks_x", "Vertex_allTracks.position.x")
        df = df.Define("PV_allTracks_y", "Vertex_allTracks.position.y")
        df = df.Define("PV_allTracks_z", "Vertex_allTracks.position.z")

        df = df.Define("PV_res_allTracks_true_x", "PV_allTracks_x - pv_true_x")
        df = df.Define("PV_res_allTracks_true_y", "PV_allTracks_y - pv_true_y")
        df = df.Define("PV_res_allTracks_true_z", "PV_allTracks_z - pv_true_z")

        df = df.Define("PV_ratio_allTracks_true_x", "PV_allTracks_x / pv_true_x")
        df = df.Define("PV_ratio_allTracks_true_y", "PV_allTracks_y / pv_true_y")
        df = df.Define("PV_ratio_allTracks_true_z", "PV_allTracks_z / pv_true_z")

        df = df.Define("PV_res_allTracks_V0_x", "PV_allTracks_x - pv_x")
        df = df.Define("PV_res_allTracks_V0_y", "PV_allTracks_y - pv_y")
        df = df.Define("PV_res_allTracks_V0_z", "PV_allTracks_z - pv_z")

        df = df.Define("PV_ratio_allTracks_V0_x", "PV_allTracks_x / pv_x")
        df = df.Define("PV_ratio_allTracks_V0_y", "PV_allTracks_y / pv_y")
        df = df.Define("PV_ratio_allTracks_V0_z", "PV_allTracks_z / pv_z")

        df = df.Define("PV_allTracks", "ROOT::VecOps::RVec<edm4hep::VertexData> v; v.push_back(Vertex_allTracks); return v;")   # remove SV tracks from PV fit


        df = df.Define("PV_res_allTracks_noBS_BS_x", "PV_allTracks_noBS_x - PV_allTracks_x")
        df = df.Define("PV_res_allTracks_noBS_BS_y", "PV_allTracks_noBS_y - PV_allTracks_y")
        df = df.Define("PV_res_allTracks_noBS_BS_z", "PV_allTracks_noBS_z - PV_allTracks_z")



        return df


        ############################################# Particle Flow Level Variables #######################################################
        df = df.Define("pfcand_isMu",     "AlephSelection::get_isType(jetConstitutentsTypes,2)")
        df = df.Define("pfcand_isEl",     "AlephSelection::get_isType(jetConstitutentsTypes,1)")
        df = df.Define("pfcand_isGamma",  "AlephSelection::get_isType(jetConstitutentsTypes,4)")
        df = df.Define("pfcand_isChargedHad", "AlephSelection::get_isType(jetConstitutentsTypes,0)")
        df = df.Define("pfcand_isNeutralHad", "AlephSelection::get_isType(jetConstitutentsTypes,5)")
        ############################################# Kinematics and PID #######################################################
        df = df.Define("pfcand_e",        "JetConstituentsUtils::get_e(jetc)") 
        df = df.Define("pfcand_p",        "JetConstituentsUtils::get_p(jetc)") 
        df = df.Define("pfcand_theta",    "JetConstituentsUtils::get_theta(jetc)") 
        df = df.Define("pfcand_phi",      "JetConstituentsUtils::get_phi(jetc)") 
        df = df.Define("pfcand_charge",   "JetConstituentsUtils::get_charge(jetc)") 
        df = df.Define("pfcand_type",     "JetConstituentsUtils::get_type(jetc)") 
        df = df.Define("pfcand_erel",     "JetConstituentsUtils::get_erel_cluster(jets, jetc)")
        df = df.Define("pfcand_erel_log", "JetConstituentsUtils::get_erel_log_cluster(jets, jetc)")
        df = df.Define("pfcand_thetarel", "JetConstituentsUtils::get_thetarel_cluster(jets, jetc)")
        df = df.Define("pfcand_phirel",   "JetConstituentsUtils::get_phirel_cluster(jets, jetc)")
        df = df.Define("Bz", '1.5')
############################################# Track Parameters and Covariance #######################################################

        df = df.Define("pfcand_dxy",        f'JetConstituentsUtils::XPtoPar_dxy(jetc, TrackStateFlipped, pv, Bz)') 
        df = df.Define("pfcand_dz",         f'JetConstituentsUtils::XPtoPar_dz(jetc, TrackStateFlipped, pv, Bz)') 
        df = df.Define("pfcand_phi0",       f'JetConstituentsUtils::XPtoPar_phi(jetc, TrackStateFlipped, pv, Bz)') 
        df = df.Define("pfcand_C",          f'JetConstituentsUtils::XPtoPar_C(jetc, TrackStateFlipped, Bz)') 
        df = df.Define("pfcand_ct",         f'JetConstituentsUtils::XPtoPar_ct(jetc, TrackStateFlipped, Bz)') 
        df = df.Define("pfcand_dptdpt",     f'JetConstituentsUtils::get_omega_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dxydxy",     f'JetConstituentsUtils::get_d0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dzdz",       f'JetConstituentsUtils::get_z0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dphidphi",   f'JetConstituentsUtils::get_phi0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_detadeta",   f'JetConstituentsUtils::get_tanlambda_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dxydz",      f'JetConstituentsUtils::get_d0_z0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dphidxy",    f'JetConstituentsUtils::get_phi0_d0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_phidz",      f'JetConstituentsUtils::get_phi0_z0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_phictgtheta",f'JetConstituentsUtils::get_tanlambda_phi0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dxyctgtheta",f'JetConstituentsUtils::get_tanlambda_d0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dlambdadz",  f'JetConstituentsUtils::get_tanlambda_z0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_cctgtheta",  f'JetConstituentsUtils::get_omega_tanlambda_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_phic",       f'JetConstituentsUtils::get_omega_phi0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_dxyc",       f'JetConstituentsUtils::get_omega_d0_cov(jetc, TrackStateFlipped)') 
        df = df.Define("pfcand_cdz",        f'JetConstituentsUtils::get_omega_z0_cov(jetc, TrackStateFlipped)')
############################################# Btag Variables #######################################################
        df = df.Define("pfcand_btagSip2dVal",   "JetConstituentsUtils::get_Sip2dVal_clusterV(jets, pfcand_dxy, pfcand_phi0, Bz)") 
        df = df.Define("pfcand_btagSip2dSig",   "JetConstituentsUtils::get_Sip2dSig(pfcand_btagSip2dVal, pfcand_dxydxy)") 
        df = df.Define("pfcand_btagSip3dVal",   "JetConstituentsUtils::get_Sip3dVal_clusterV(jets, pfcand_dxy, pfcand_dz, pfcand_phi0, Bz)") 
        df = df.Define("pfcand_btagSip3dSig",   "JetConstituentsUtils::get_Sip3dSig(pfcand_btagSip3dVal, pfcand_dxydxy, pfcand_dzdz)") 
        df = df.Define("pfcand_btagJetDistVal","JetConstituentsUtils::get_JetDistVal_clusterV(jets, jetc, pfcand_dxy, pfcand_dz, pfcand_phi0, Bz)") 
        df = df.Define("pfcand_btagJetDistSig","JetConstituentsUtils::get_JetDistSig(pfcand_btagJetDistVal, pfcand_dxydxy, pfcand_dzdz)")
        ############################################# Jet Level Variables #######################################################
        df=df.Define("event_njet",   "JetConstituentsUtils::count_jets(jetc)")
        df.Filter("event_njet > 1")
        ##############################################################################################################
        df = df.Define("jet_p", "JetClusteringUtils::get_p(jets)")
        df = df.Define("jet_pT", "JetClusteringUtils::get_pt(jets)")
        df = df.Define("jet_e", "JetClusteringUtils::get_e(jets)")
        df = df.Define("jet_mass", "JetClusteringUtils::get_m(jets)")
        df = df.Define("jet_phi", "JetClusteringUtils::get_phi(jets)")
        df = df.Define("jet_theta", "JetClusteringUtils::get_theta(jets)")
        df = df.Define("jet_eta", "JetClusteringUtils::get_eta(jets)")
        ##############################################################################################################
        df = df.Define("jet_nconst", "JetConstituentsUtils::count_consts(jetc)") 
        df = df.Define(f"jet_nmu",    f"JetConstituentsUtils::count_type(pfcand_isMu)") 
        df = df.Define(f"jet_nel",    f"JetConstituentsUtils::count_type(pfcand_isEl)") 
        df = df.Define(f"jet_nchad",  f"JetConstituentsUtils::count_type(pfcand_isChargedHad)") 
        df = df.Define(f"jet_ngamma", f"JetConstituentsUtils::count_type(pfcand_isGamma)") 
        df = df.Define(f"jet_nnhad",  f"JetConstituentsUtils::count_type(pfcand_isNeutralHad)")
        ##############################################################################################################
        df = df.Define("jet_constituents_dEdx_pads_objs", "AlephSelection::build_constituents_dEdx()(RecoParticles, _RecoParticles_tracks.index, dEdxPads, _dEdxPads_track.index, _jetc)" )
        df = df.Define("pfcand_dEdx_pads_type", "AlephSelection::get_dEdx_type(jet_constituents_dEdx_pads_objs)")
        df = df.Define("pfcand_dEdx_pads_value", "AlephSelection::get_dEdx_value(jet_constituents_dEdx_pads_objs)")
        df = df.Define("pfcand_dEdx_pads_error", "AlephSelection::get_dEdx_error(jet_constituents_dEdx_pads_objs)")
        df = df.Define("jet_constituents_dEdx_wires_objs", "AlephSelection::build_constituents_dEdx()(RecoParticles, _RecoParticles_tracks.index, dEdxWires, _dEdxWires_track.index, _jetc)" )
        df = df.Define("pfcand_dEdx_wires_type", "AlephSelection::get_dEdx_type(jet_constituents_dEdx_wires_objs)")
        df = df.Define("pfcand_dEdx_wires_value", "AlephSelection::get_dEdx_value(jet_constituents_dEdx_wires_objs)")
        df = df.Define("pfcand_dEdx_wires_error", "AlephSelection::get_dEdx_error(jet_constituents_dEdx_wires_objs)")
        ##############################################################################################################


        return df

    def output():

        return [
            "event_type",
            #"event_invariant_mass","event_njet", 
            # Jet Level Variables  
            #"jet_mass","jet_p","jet_e", "jet_phi", "jet_theta", "jet_pT","jet_eta",
            #                "jet_nnhad","jet_ngamma","jet_nchad","jet_nel", "jet_nmu", "jet_nconst",


            #the dEdX values associated to the jet constituents:
            #"pfcand_dEdx_pads_type", "pfcand_dEdx_pads_value", "pfcand_dEdx_pads_error",
            #    "pfcand_dEdx_wires_type", "pfcand_dEdx_wires_value", "pfcand_dEdx_wires_error",
            #    "pfcand_isMu", "pfcand_isEl", "pfcand_isChargedHad", "pfcand_isGamma", "pfcand_isNeutralHad",
            #    "pfcand_e", "pfcand_p", "pfcand_theta", "pfcand_phi", "pfcand_charge", "pfcand_type",
            #    "pfcand_erel", "pfcand_erel_log", "pfcand_thetarel", "pfcand_phirel", 
            #    "pfcand_dxy", "pfcand_dz", "pfcand_phi0", "pfcand_C", "pfcand_ct",
            #    "pfcand_dptdpt", "pfcand_dxydxy", "pfcand_dzdz", "pfcand_dphidphi", "pfcand_detadeta",
            #    "pfcand_dxydz", "pfcand_dphidxy", "pfcand_phidz", "pfcand_phictgtheta", "pfcand_dxyctgtheta",
            #    "pfcand_dlambdadz", "pfcand_cctgtheta", "pfcand_phic", "pfcand_dxyc", "pfcand_cdz",
            #    "pfcand_btagSip2dVal", "pfcand_btagSip2dSig", "pfcand_btagSip3dVal", "pfcand_btagSip3dSig", 
            #    "pfcand_btagJetDistVal", "pfcand_btagJetDistSig",
                "pv_x", "pv_y", "pv_z", "pv_true_x", "pv_true_y", "pv_true_z", "pv_res_x", "pv_res_y", "pv_res_z",
                "ntracks",
                "PV_allTracks_noBS_x", "PV_allTracks_noBS_y", "PV_allTracks_noBS_z", 
                "PV_res_allTracks_noBS_true_x", "PV_res_allTracks_noBS_true_y", "PV_res_allTracks_noBS_true_z",
                "PV_ratio_allTracks_noBS_true_x", "PV_ratio_allTracks_noBS_true_y", "PV_ratio_allTracks_noBS_true_z",
                "PV_res_allTracks_noBS_V0_x", "PV_res_allTracks_noBS_V0_y", "PV_res_allTracks_noBS_V0_z",
                "PV_ratio_allTracks_noBS_V0_x", "PV_ratio_allTracks_noBS_V0_y", "PV_ratio_allTracks_noBS_V0_z",
                "PV_allTracks_x", "PV_allTracks_y", "PV_allTracks_z", 
                "PV_res_allTracks_true_x", "PV_res_allTracks_true_y", "PV_res_allTracks_true_z",
                "PV_ratio_allTracks_true_x", "PV_ratio_allTracks_true_y", "PV_ratio_allTracks_true_z",
                "PV_res_allTracks_V0_x", "PV_res_allTracks_V0_y", "PV_res_allTracks_V0_z",
                "PV_ratio_allTracks_V0_x", "PV_ratio_allTracks_V0_y", "PV_ratio_allTracks_V0_z",
                "PV_res_allTracks_noBS_BS_x", "PV_res_allTracks_noBS_BS_y", "PV_res_allTracks_noBS_BS_z",
        ]
