
include_directories ( $ENV{G4INCLUDE} )
include_directories( $ENV{PANDORA_INC} )
include_directories( $ENV{LARPANDORACONTENT_INC} )

art_make(  
          MODULE_LIBRARIES  larsim_Simulation 
			    lardataobj_Simulation
			    lardata_RecoObjects
			    lardataobj_RecoBase
			    larcorealg_Geometry
                            larpandora_LArPandoraInterface
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
			    art_Persistency_Provenance
			    art_Utilities
			    art_Framework_Services_Registry
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
                            ${PANDORASDK}
         )

install_headers()
install_fhicl()
install_source()


