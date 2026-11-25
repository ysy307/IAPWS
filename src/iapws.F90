module iapws
    use :: module_iapws, only:type_iapws_property, IAPWS_BATCH_OK, IAPWS_BATCH_INVALID
    use :: module_iapws95, only:type_iapws95
    use :: module_iapws97, only:type_iapws97
    use :: module_iapws06_iceIh, only:type_iapws06
    implicit none

    public :: IAPWS_BATCH_OK, IAPWS_BATCH_INVALID
    public :: type_iapws_property
    public :: type_iapws95
    public :: type_iapws97
    public :: type_iapws06

end module iapws
