
include_directories ( $ENV{G4INCLUDE} )

art_make(  
          MODULE_LIBRARIES  larsim_Simulation 
			    lardataobj_Simulation
			    lardata_RecoObjects
			    lardataobj_RecoBase
			    larcorealg_Geometry
			    larcore_Geometry_Geometry_service
			    lardataobj_AnalysisBase
			    lardata_Utilities
			    nutools_MagneticField_MagneticField_service
			    nusimdata_SimulationBase
			    ${ART_FRAMEWORK_CORE}
			    ${ART_FRAMEWORK_PRINCIPAL}
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${ART_FRAMEWORK_SERVICES_BASIC}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
			    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			    art_Persistency_Common
			    art_Utilities
			    canvas
			    ${MF_MESSAGELOGGER}
			    ${MF_UTILITIES}
			    ${FHICLCPP}
			    cetlib cetlib_except
			    ${CLHEP}
			    ${ROOT_BASIC_LIB_LIST}
			    ${G4_LIB_LIST}
			    larreco_RecoAlg
			    larreco_RecoAlg_ClusterRecoUtil
			    larreco_ClusterFinder
			    ${ROOT_GENVECTOR}
			    larsim_MCCheater_BackTrackerService_service
         )

install_headers()
install_fhicl()
install_source()


