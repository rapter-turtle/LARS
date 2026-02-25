

/*
WARNING: THIS FILE IS AUTO-GENERATED. DO NOT MODIFY.

This file was generated from pils.idl
using RTI Code Generator (rtiddsgen) version 4.3.0.
The rtiddsgen tool is part of the RTI Connext DDS distribution.
For more information, type 'rtiddsgen -help' at a command shell
or consult the Code Generator User's Manual.
*/

#ifndef pilsPlugin_475505405_h
#define pilsPlugin_475505405_h

#include "pils.hpp"

struct RTICdrStream;

#ifndef pres_typePlugin_h
#include "pres/pres_typePlugin.h"
#endif

#if (defined(RTI_WIN32) || defined (RTI_WINCE) || defined(RTI_INTIME)) && defined(NDDS_USER_DLL_EXPORT)
/* If the code is building on Windows, start exporting symbols.
*/
#undef NDDSUSERDllExport
#define NDDSUSERDllExport __declspec(dllexport)
#endif

namespace ULARS {
    namespace EO {

        #define USVControlTypePlugin_get_sample PRESTypePluginDefaultEndpointData_getSample

        #define USVControlTypePlugin_get_buffer PRESTypePluginDefaultEndpointData_getBuffer 
        #define USVControlTypePlugin_return_buffer PRESTypePluginDefaultEndpointData_returnBuffer

        #define USVControlTypePlugin_create_sample PRESTypePluginDefaultEndpointData_createSample 
        #define USVControlTypePlugin_destroy_sample PRESTypePluginDefaultEndpointData_deleteSample 

        /* --------------------------------------------------------------------------------------
        Support functions:
        * -------------------------------------------------------------------------------------- */

        NDDSUSERDllExport extern USVControlType*
        USVControlTypePluginSupport_create_data_w_params(
            const struct DDS_TypeAllocationParams_t * alloc_params);

        NDDSUSERDllExport extern USVControlType*
        USVControlTypePluginSupport_create_data_ex(RTIBool allocate_pointers);

        NDDSUSERDllExport extern USVControlType*
        USVControlTypePluginSupport_create_data(void);

        NDDSUSERDllExport extern RTIBool 
        USVControlTypePluginSupport_copy_data(
            USVControlType *out,
            const USVControlType *in);

        NDDSUSERDllExport extern void 
        USVControlTypePluginSupport_destroy_data_w_params(
            USVControlType *sample,
            const struct DDS_TypeDeallocationParams_t * dealloc_params);

        NDDSUSERDllExport extern void 
        USVControlTypePluginSupport_destroy_data_ex(
            USVControlType *sample,RTIBool deallocate_pointers);

        NDDSUSERDllExport extern void 
        USVControlTypePluginSupport_destroy_data(
            USVControlType *sample);

        NDDSUSERDllExport extern void 
        USVControlTypePluginSupport_print_data(
            const USVControlType *sample,
            const char *desc,
            unsigned int indent);

        /* ----------------------------------------------------------------------------
        Callback functions:
        * ---------------------------------------------------------------------------- */

        NDDSUSERDllExport extern PRESTypePluginParticipantData 
        USVControlTypePlugin_on_participant_attached(
            void *registration_data, 
            const struct PRESTypePluginParticipantInfo *participant_info,
            RTIBool top_level_registration, 
            void *container_plugin_context,
            RTICdrTypeCode *typeCode);

        NDDSUSERDllExport extern void 
        USVControlTypePlugin_on_participant_detached(
            PRESTypePluginParticipantData participant_data);

        NDDSUSERDllExport extern PRESTypePluginEndpointData 
        USVControlTypePlugin_on_endpoint_attached(
            PRESTypePluginParticipantData participant_data,
            const struct PRESTypePluginEndpointInfo *endpoint_info,
            RTIBool top_level_registration, 
            void *container_plugin_context);

        NDDSUSERDllExport extern void 
        USVControlTypePlugin_on_endpoint_detached(
            PRESTypePluginEndpointData endpoint_data);

        NDDSUSERDllExport extern void    
        USVControlTypePlugin_return_sample(
            PRESTypePluginEndpointData endpoint_data,
            USVControlType *sample,
            void *handle);    

        NDDSUSERDllExport extern RTIBool 
        USVControlTypePlugin_copy_sample(
            PRESTypePluginEndpointData endpoint_data,
            USVControlType *out,
            const USVControlType *in);

        /* ----------------------------------------------------------------------------
        (De)Serialize functions:
        * ------------------------------------------------------------------------- */

        NDDSUSERDllExport extern RTIBool
        USVControlTypePlugin_serialize_to_cdr_buffer(
            char * buffer,
            unsigned int * length,
            const USVControlType *sample,
            ::dds::core::policy::DataRepresentationId representation
            = ::dds::core::policy::DataRepresentation::xcdr()); 

        NDDSUSERDllExport extern RTIBool 
        USVControlTypePlugin_deserialize(
            PRESTypePluginEndpointData endpoint_data,
            USVControlType **sample, 
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_sample, 
            void *endpoint_plugin_qos);

        NDDSUSERDllExport extern RTIBool
        USVControlTypePlugin_deserialize_from_cdr_buffer(
            USVControlType *sample,
            const char * buffer,
            unsigned int length);    

        NDDSUSERDllExport extern unsigned int 
        USVControlTypePlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        /* --------------------------------------------------------------------------------------
        Key Management functions:
        * -------------------------------------------------------------------------------------- */
        NDDSUSERDllExport extern PRESTypePluginKeyKind 
        USVControlTypePlugin_get_key_kind(void);

        NDDSUSERDllExport extern unsigned int 
        USVControlTypePlugin_get_serialized_key_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        NDDSUSERDllExport extern unsigned int 
        USVControlTypePlugin_get_serialized_key_max_size_for_keyhash(
            PRESTypePluginEndpointData endpoint_data,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        NDDSUSERDllExport extern RTIBool 
        USVControlTypePlugin_deserialize_key(
            PRESTypePluginEndpointData endpoint_data,
            USVControlType ** sample,
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_key,
            void *endpoint_plugin_qos);

        /* Plugin Functions */
        NDDSUSERDllExport extern struct PRESTypePlugin*
        USVControlTypePlugin_new(void);

        NDDSUSERDllExport extern void
        USVControlTypePlugin_delete(struct PRESTypePlugin *);

    } /* namespace EO  */
    namespace C2 {

        #define WaypointPlugin_get_sample PRESTypePluginDefaultEndpointData_getSample

        #define WaypointPlugin_get_buffer PRESTypePluginDefaultEndpointData_getBuffer 
        #define WaypointPlugin_return_buffer PRESTypePluginDefaultEndpointData_returnBuffer

        #define WaypointPlugin_create_sample PRESTypePluginDefaultEndpointData_createSample 
        #define WaypointPlugin_destroy_sample PRESTypePluginDefaultEndpointData_deleteSample 

        /* --------------------------------------------------------------------------------------
        Support functions:
        * -------------------------------------------------------------------------------------- */

        NDDSUSERDllExport extern Waypoint*
        WaypointPluginSupport_create_data_w_params(
            const struct DDS_TypeAllocationParams_t * alloc_params);

        NDDSUSERDllExport extern Waypoint*
        WaypointPluginSupport_create_data_ex(RTIBool allocate_pointers);

        NDDSUSERDllExport extern Waypoint*
        WaypointPluginSupport_create_data(void);

        NDDSUSERDllExport extern RTIBool 
        WaypointPluginSupport_copy_data(
            Waypoint *out,
            const Waypoint *in);

        NDDSUSERDllExport extern void 
        WaypointPluginSupport_destroy_data_w_params(
            Waypoint *sample,
            const struct DDS_TypeDeallocationParams_t * dealloc_params);

        NDDSUSERDllExport extern void 
        WaypointPluginSupport_destroy_data_ex(
            Waypoint *sample,RTIBool deallocate_pointers);

        NDDSUSERDllExport extern void 
        WaypointPluginSupport_destroy_data(
            Waypoint *sample);

        NDDSUSERDllExport extern void 
        WaypointPluginSupport_print_data(
            const Waypoint *sample,
            const char *desc,
            unsigned int indent);

        /* ----------------------------------------------------------------------------
        Callback functions:
        * ---------------------------------------------------------------------------- */

        NDDSUSERDllExport extern PRESTypePluginParticipantData 
        WaypointPlugin_on_participant_attached(
            void *registration_data, 
            const struct PRESTypePluginParticipantInfo *participant_info,
            RTIBool top_level_registration, 
            void *container_plugin_context,
            RTICdrTypeCode *typeCode);

        NDDSUSERDllExport extern void 
        WaypointPlugin_on_participant_detached(
            PRESTypePluginParticipantData participant_data);

        NDDSUSERDllExport extern PRESTypePluginEndpointData 
        WaypointPlugin_on_endpoint_attached(
            PRESTypePluginParticipantData participant_data,
            const struct PRESTypePluginEndpointInfo *endpoint_info,
            RTIBool top_level_registration, 
            void *container_plugin_context);

        NDDSUSERDllExport extern void 
        WaypointPlugin_on_endpoint_detached(
            PRESTypePluginEndpointData endpoint_data);

        NDDSUSERDllExport extern void    
        WaypointPlugin_return_sample(
            PRESTypePluginEndpointData endpoint_data,
            Waypoint *sample,
            void *handle);    

        NDDSUSERDllExport extern RTIBool 
        WaypointPlugin_copy_sample(
            PRESTypePluginEndpointData endpoint_data,
            Waypoint *out,
            const Waypoint *in);

        /* ----------------------------------------------------------------------------
        (De)Serialize functions:
        * ------------------------------------------------------------------------- */

        NDDSUSERDllExport extern RTIBool
        WaypointPlugin_serialize_to_cdr_buffer(
            char * buffer,
            unsigned int * length,
            const Waypoint *sample,
            ::dds::core::policy::DataRepresentationId representation
            = ::dds::core::policy::DataRepresentation::xcdr()); 

        NDDSUSERDllExport extern RTIBool 
        WaypointPlugin_deserialize(
            PRESTypePluginEndpointData endpoint_data,
            Waypoint **sample, 
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_sample, 
            void *endpoint_plugin_qos);

        NDDSUSERDllExport extern RTIBool
        WaypointPlugin_deserialize_from_cdr_buffer(
            Waypoint *sample,
            const char * buffer,
            unsigned int length);    

        NDDSUSERDllExport extern unsigned int 
        WaypointPlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        /* --------------------------------------------------------------------------------------
        Key Management functions:
        * -------------------------------------------------------------------------------------- */
        NDDSUSERDllExport extern PRESTypePluginKeyKind 
        WaypointPlugin_get_key_kind(void);

        NDDSUSERDllExport extern unsigned int 
        WaypointPlugin_get_serialized_key_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        NDDSUSERDllExport extern unsigned int 
        WaypointPlugin_get_serialized_key_max_size_for_keyhash(
            PRESTypePluginEndpointData endpoint_data,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        NDDSUSERDllExport extern RTIBool 
        WaypointPlugin_deserialize_key(
            PRESTypePluginEndpointData endpoint_data,
            Waypoint ** sample,
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_key,
            void *endpoint_plugin_qos);

        /* Plugin Functions */
        NDDSUSERDllExport extern struct PRESTypePlugin*
        WaypointPlugin_new(void);

        NDDSUSERDllExport extern void
        WaypointPlugin_delete(struct PRESTypePlugin *);

        #define GlobalWaypointTypePlugin_get_sample PRESTypePluginDefaultEndpointData_getSample

        #define GlobalWaypointTypePlugin_get_buffer PRESTypePluginDefaultEndpointData_getBuffer 
        #define GlobalWaypointTypePlugin_return_buffer PRESTypePluginDefaultEndpointData_returnBuffer

        #define GlobalWaypointTypePlugin_create_sample PRESTypePluginDefaultEndpointData_createSample 
        #define GlobalWaypointTypePlugin_destroy_sample PRESTypePluginDefaultEndpointData_deleteSample 

        /* --------------------------------------------------------------------------------------
        Support functions:
        * -------------------------------------------------------------------------------------- */

        NDDSUSERDllExport extern GlobalWaypointType*
        GlobalWaypointTypePluginSupport_create_data_w_params(
            const struct DDS_TypeAllocationParams_t * alloc_params);

        NDDSUSERDllExport extern GlobalWaypointType*
        GlobalWaypointTypePluginSupport_create_data_ex(RTIBool allocate_pointers);

        NDDSUSERDllExport extern GlobalWaypointType*
        GlobalWaypointTypePluginSupport_create_data(void);

        NDDSUSERDllExport extern RTIBool 
        GlobalWaypointTypePluginSupport_copy_data(
            GlobalWaypointType *out,
            const GlobalWaypointType *in);

        NDDSUSERDllExport extern void 
        GlobalWaypointTypePluginSupport_destroy_data_w_params(
            GlobalWaypointType *sample,
            const struct DDS_TypeDeallocationParams_t * dealloc_params);

        NDDSUSERDllExport extern void 
        GlobalWaypointTypePluginSupport_destroy_data_ex(
            GlobalWaypointType *sample,RTIBool deallocate_pointers);

        NDDSUSERDllExport extern void 
        GlobalWaypointTypePluginSupport_destroy_data(
            GlobalWaypointType *sample);

        NDDSUSERDllExport extern void 
        GlobalWaypointTypePluginSupport_print_data(
            const GlobalWaypointType *sample,
            const char *desc,
            unsigned int indent);

        /* ----------------------------------------------------------------------------
        Callback functions:
        * ---------------------------------------------------------------------------- */

        NDDSUSERDllExport extern PRESTypePluginParticipantData 
        GlobalWaypointTypePlugin_on_participant_attached(
            void *registration_data, 
            const struct PRESTypePluginParticipantInfo *participant_info,
            RTIBool top_level_registration, 
            void *container_plugin_context,
            RTICdrTypeCode *typeCode);

        NDDSUSERDllExport extern void 
        GlobalWaypointTypePlugin_on_participant_detached(
            PRESTypePluginParticipantData participant_data);

        NDDSUSERDllExport extern PRESTypePluginEndpointData 
        GlobalWaypointTypePlugin_on_endpoint_attached(
            PRESTypePluginParticipantData participant_data,
            const struct PRESTypePluginEndpointInfo *endpoint_info,
            RTIBool top_level_registration, 
            void *container_plugin_context);

        NDDSUSERDllExport extern void 
        GlobalWaypointTypePlugin_on_endpoint_detached(
            PRESTypePluginEndpointData endpoint_data);

        NDDSUSERDllExport extern void    
        GlobalWaypointTypePlugin_return_sample(
            PRESTypePluginEndpointData endpoint_data,
            GlobalWaypointType *sample,
            void *handle);    

        NDDSUSERDllExport extern RTIBool 
        GlobalWaypointTypePlugin_copy_sample(
            PRESTypePluginEndpointData endpoint_data,
            GlobalWaypointType *out,
            const GlobalWaypointType *in);

        /* ----------------------------------------------------------------------------
        (De)Serialize functions:
        * ------------------------------------------------------------------------- */

        NDDSUSERDllExport extern RTIBool
        GlobalWaypointTypePlugin_serialize_to_cdr_buffer(
            char * buffer,
            unsigned int * length,
            const GlobalWaypointType *sample,
            ::dds::core::policy::DataRepresentationId representation
            = ::dds::core::policy::DataRepresentation::xcdr()); 

        NDDSUSERDllExport extern RTIBool 
        GlobalWaypointTypePlugin_deserialize(
            PRESTypePluginEndpointData endpoint_data,
            GlobalWaypointType **sample, 
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_sample, 
            void *endpoint_plugin_qos);

        NDDSUSERDllExport extern RTIBool
        GlobalWaypointTypePlugin_deserialize_from_cdr_buffer(
            GlobalWaypointType *sample,
            const char * buffer,
            unsigned int length);    

        NDDSUSERDllExport extern unsigned int 
        GlobalWaypointTypePlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        /* --------------------------------------------------------------------------------------
        Key Management functions:
        * -------------------------------------------------------------------------------------- */
        NDDSUSERDllExport extern PRESTypePluginKeyKind 
        GlobalWaypointTypePlugin_get_key_kind(void);

        NDDSUSERDllExport extern unsigned int 
        GlobalWaypointTypePlugin_get_serialized_key_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        NDDSUSERDllExport extern unsigned int 
        GlobalWaypointTypePlugin_get_serialized_key_max_size_for_keyhash(
            PRESTypePluginEndpointData endpoint_data,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        NDDSUSERDllExport extern RTIBool 
        GlobalWaypointTypePlugin_deserialize_key(
            PRESTypePluginEndpointData endpoint_data,
            GlobalWaypointType ** sample,
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_key,
            void *endpoint_plugin_qos);

        /* Plugin Functions */
        NDDSUSERDllExport extern struct PRESTypePlugin*
        GlobalWaypointTypePlugin_new(void);

        NDDSUSERDllExport extern void
        GlobalWaypointTypePlugin_delete(struct PRESTypePlugin *);

        #define StatusTypePlugin_get_sample PRESTypePluginDefaultEndpointData_getSample

        #define StatusTypePlugin_get_buffer PRESTypePluginDefaultEndpointData_getBuffer 
        #define StatusTypePlugin_return_buffer PRESTypePluginDefaultEndpointData_returnBuffer

        #define StatusTypePlugin_create_sample PRESTypePluginDefaultEndpointData_createSample 
        #define StatusTypePlugin_destroy_sample PRESTypePluginDefaultEndpointData_deleteSample 

        /* --------------------------------------------------------------------------------------
        Support functions:
        * -------------------------------------------------------------------------------------- */

        NDDSUSERDllExport extern StatusType*
        StatusTypePluginSupport_create_data_w_params(
            const struct DDS_TypeAllocationParams_t * alloc_params);

        NDDSUSERDllExport extern StatusType*
        StatusTypePluginSupport_create_data_ex(RTIBool allocate_pointers);

        NDDSUSERDllExport extern StatusType*
        StatusTypePluginSupport_create_data(void);

        NDDSUSERDllExport extern RTIBool 
        StatusTypePluginSupport_copy_data(
            StatusType *out,
            const StatusType *in);

        NDDSUSERDllExport extern void 
        StatusTypePluginSupport_destroy_data_w_params(
            StatusType *sample,
            const struct DDS_TypeDeallocationParams_t * dealloc_params);

        NDDSUSERDllExport extern void 
        StatusTypePluginSupport_destroy_data_ex(
            StatusType *sample,RTIBool deallocate_pointers);

        NDDSUSERDllExport extern void 
        StatusTypePluginSupport_destroy_data(
            StatusType *sample);

        NDDSUSERDllExport extern void 
        StatusTypePluginSupport_print_data(
            const StatusType *sample,
            const char *desc,
            unsigned int indent);

        /* ----------------------------------------------------------------------------
        Callback functions:
        * ---------------------------------------------------------------------------- */

        NDDSUSERDllExport extern PRESTypePluginParticipantData 
        StatusTypePlugin_on_participant_attached(
            void *registration_data, 
            const struct PRESTypePluginParticipantInfo *participant_info,
            RTIBool top_level_registration, 
            void *container_plugin_context,
            RTICdrTypeCode *typeCode);

        NDDSUSERDllExport extern void 
        StatusTypePlugin_on_participant_detached(
            PRESTypePluginParticipantData participant_data);

        NDDSUSERDllExport extern PRESTypePluginEndpointData 
        StatusTypePlugin_on_endpoint_attached(
            PRESTypePluginParticipantData participant_data,
            const struct PRESTypePluginEndpointInfo *endpoint_info,
            RTIBool top_level_registration, 
            void *container_plugin_context);

        NDDSUSERDllExport extern void 
        StatusTypePlugin_on_endpoint_detached(
            PRESTypePluginEndpointData endpoint_data);

        NDDSUSERDllExport extern void    
        StatusTypePlugin_return_sample(
            PRESTypePluginEndpointData endpoint_data,
            StatusType *sample,
            void *handle);    

        NDDSUSERDllExport extern RTIBool 
        StatusTypePlugin_copy_sample(
            PRESTypePluginEndpointData endpoint_data,
            StatusType *out,
            const StatusType *in);

        /* ----------------------------------------------------------------------------
        (De)Serialize functions:
        * ------------------------------------------------------------------------- */

        NDDSUSERDllExport extern RTIBool
        StatusTypePlugin_serialize_to_cdr_buffer(
            char * buffer,
            unsigned int * length,
            const StatusType *sample,
            ::dds::core::policy::DataRepresentationId representation
            = ::dds::core::policy::DataRepresentation::xcdr()); 

        NDDSUSERDllExport extern RTIBool 
        StatusTypePlugin_deserialize(
            PRESTypePluginEndpointData endpoint_data,
            StatusType **sample, 
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_sample, 
            void *endpoint_plugin_qos);

        NDDSUSERDllExport extern RTIBool
        StatusTypePlugin_deserialize_from_cdr_buffer(
            StatusType *sample,
            const char * buffer,
            unsigned int length);    

        NDDSUSERDllExport extern unsigned int 
        StatusTypePlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        /* --------------------------------------------------------------------------------------
        Key Management functions:
        * -------------------------------------------------------------------------------------- */
        NDDSUSERDllExport extern PRESTypePluginKeyKind 
        StatusTypePlugin_get_key_kind(void);

        NDDSUSERDllExport extern unsigned int 
        StatusTypePlugin_get_serialized_key_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        NDDSUSERDllExport extern unsigned int 
        StatusTypePlugin_get_serialized_key_max_size_for_keyhash(
            PRESTypePluginEndpointData endpoint_data,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        NDDSUSERDllExport extern RTIBool 
        StatusTypePlugin_deserialize_key(
            PRESTypePluginEndpointData endpoint_data,
            StatusType ** sample,
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_key,
            void *endpoint_plugin_qos);

        /* Plugin Functions */
        NDDSUSERDllExport extern struct PRESTypePlugin*
        StatusTypePlugin_new(void);

        NDDSUSERDllExport extern void
        StatusTypePlugin_delete(struct PRESTypePlugin *);

        #define ObstacleTypePlugin_get_sample PRESTypePluginDefaultEndpointData_getSample

        #define ObstacleTypePlugin_get_buffer PRESTypePluginDefaultEndpointData_getBuffer 
        #define ObstacleTypePlugin_return_buffer PRESTypePluginDefaultEndpointData_returnBuffer

        #define ObstacleTypePlugin_create_sample PRESTypePluginDefaultEndpointData_createSample 
        #define ObstacleTypePlugin_destroy_sample PRESTypePluginDefaultEndpointData_deleteSample 

        /* --------------------------------------------------------------------------------------
        Support functions:
        * -------------------------------------------------------------------------------------- */

        NDDSUSERDllExport extern ObstacleType*
        ObstacleTypePluginSupport_create_data_w_params(
            const struct DDS_TypeAllocationParams_t * alloc_params);

        NDDSUSERDllExport extern ObstacleType*
        ObstacleTypePluginSupport_create_data_ex(RTIBool allocate_pointers);

        NDDSUSERDllExport extern ObstacleType*
        ObstacleTypePluginSupport_create_data(void);

        NDDSUSERDllExport extern RTIBool 
        ObstacleTypePluginSupport_copy_data(
            ObstacleType *out,
            const ObstacleType *in);

        NDDSUSERDllExport extern void 
        ObstacleTypePluginSupport_destroy_data_w_params(
            ObstacleType *sample,
            const struct DDS_TypeDeallocationParams_t * dealloc_params);

        NDDSUSERDllExport extern void 
        ObstacleTypePluginSupport_destroy_data_ex(
            ObstacleType *sample,RTIBool deallocate_pointers);

        NDDSUSERDllExport extern void 
        ObstacleTypePluginSupport_destroy_data(
            ObstacleType *sample);

        NDDSUSERDllExport extern void 
        ObstacleTypePluginSupport_print_data(
            const ObstacleType *sample,
            const char *desc,
            unsigned int indent);

        /* ----------------------------------------------------------------------------
        Callback functions:
        * ---------------------------------------------------------------------------- */

        NDDSUSERDllExport extern PRESTypePluginParticipantData 
        ObstacleTypePlugin_on_participant_attached(
            void *registration_data, 
            const struct PRESTypePluginParticipantInfo *participant_info,
            RTIBool top_level_registration, 
            void *container_plugin_context,
            RTICdrTypeCode *typeCode);

        NDDSUSERDllExport extern void 
        ObstacleTypePlugin_on_participant_detached(
            PRESTypePluginParticipantData participant_data);

        NDDSUSERDllExport extern PRESTypePluginEndpointData 
        ObstacleTypePlugin_on_endpoint_attached(
            PRESTypePluginParticipantData participant_data,
            const struct PRESTypePluginEndpointInfo *endpoint_info,
            RTIBool top_level_registration, 
            void *container_plugin_context);

        NDDSUSERDllExport extern void 
        ObstacleTypePlugin_on_endpoint_detached(
            PRESTypePluginEndpointData endpoint_data);

        NDDSUSERDllExport extern void    
        ObstacleTypePlugin_return_sample(
            PRESTypePluginEndpointData endpoint_data,
            ObstacleType *sample,
            void *handle);    

        NDDSUSERDllExport extern RTIBool 
        ObstacleTypePlugin_copy_sample(
            PRESTypePluginEndpointData endpoint_data,
            ObstacleType *out,
            const ObstacleType *in);

        /* ----------------------------------------------------------------------------
        (De)Serialize functions:
        * ------------------------------------------------------------------------- */

        NDDSUSERDllExport extern RTIBool
        ObstacleTypePlugin_serialize_to_cdr_buffer(
            char * buffer,
            unsigned int * length,
            const ObstacleType *sample,
            ::dds::core::policy::DataRepresentationId representation
            = ::dds::core::policy::DataRepresentation::xcdr()); 

        NDDSUSERDllExport extern RTIBool 
        ObstacleTypePlugin_deserialize(
            PRESTypePluginEndpointData endpoint_data,
            ObstacleType **sample, 
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_sample, 
            void *endpoint_plugin_qos);

        NDDSUSERDllExport extern RTIBool
        ObstacleTypePlugin_deserialize_from_cdr_buffer(
            ObstacleType *sample,
            const char * buffer,
            unsigned int length);    

        NDDSUSERDllExport extern unsigned int 
        ObstacleTypePlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        /* --------------------------------------------------------------------------------------
        Key Management functions:
        * -------------------------------------------------------------------------------------- */
        NDDSUSERDllExport extern PRESTypePluginKeyKind 
        ObstacleTypePlugin_get_key_kind(void);

        NDDSUSERDllExport extern unsigned int 
        ObstacleTypePlugin_get_serialized_key_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        NDDSUSERDllExport extern unsigned int 
        ObstacleTypePlugin_get_serialized_key_max_size_for_keyhash(
            PRESTypePluginEndpointData endpoint_data,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        NDDSUSERDllExport extern RTIBool 
        ObstacleTypePlugin_deserialize_key(
            PRESTypePluginEndpointData endpoint_data,
            ObstacleType ** sample,
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_key,
            void *endpoint_plugin_qos);

        /* Plugin Functions */
        NDDSUSERDllExport extern struct PRESTypePlugin*
        ObstacleTypePlugin_new(void);

        NDDSUSERDllExport extern void
        ObstacleTypePlugin_delete(struct PRESTypePlugin *);

        #define CommandInitialTypePlugin_get_sample PRESTypePluginDefaultEndpointData_getSample

        #define CommandInitialTypePlugin_get_buffer PRESTypePluginDefaultEndpointData_getBuffer 
        #define CommandInitialTypePlugin_return_buffer PRESTypePluginDefaultEndpointData_returnBuffer

        #define CommandInitialTypePlugin_create_sample PRESTypePluginDefaultEndpointData_createSample 
        #define CommandInitialTypePlugin_destroy_sample PRESTypePluginDefaultEndpointData_deleteSample 

        /* --------------------------------------------------------------------------------------
        Support functions:
        * -------------------------------------------------------------------------------------- */

        NDDSUSERDllExport extern CommandInitialType*
        CommandInitialTypePluginSupport_create_data_w_params(
            const struct DDS_TypeAllocationParams_t * alloc_params);

        NDDSUSERDllExport extern CommandInitialType*
        CommandInitialTypePluginSupport_create_data_ex(RTIBool allocate_pointers);

        NDDSUSERDllExport extern CommandInitialType*
        CommandInitialTypePluginSupport_create_data(void);

        NDDSUSERDllExport extern RTIBool 
        CommandInitialTypePluginSupport_copy_data(
            CommandInitialType *out,
            const CommandInitialType *in);

        NDDSUSERDllExport extern void 
        CommandInitialTypePluginSupport_destroy_data_w_params(
            CommandInitialType *sample,
            const struct DDS_TypeDeallocationParams_t * dealloc_params);

        NDDSUSERDllExport extern void 
        CommandInitialTypePluginSupport_destroy_data_ex(
            CommandInitialType *sample,RTIBool deallocate_pointers);

        NDDSUSERDllExport extern void 
        CommandInitialTypePluginSupport_destroy_data(
            CommandInitialType *sample);

        NDDSUSERDllExport extern void 
        CommandInitialTypePluginSupport_print_data(
            const CommandInitialType *sample,
            const char *desc,
            unsigned int indent);

        /* ----------------------------------------------------------------------------
        Callback functions:
        * ---------------------------------------------------------------------------- */

        NDDSUSERDllExport extern PRESTypePluginParticipantData 
        CommandInitialTypePlugin_on_participant_attached(
            void *registration_data, 
            const struct PRESTypePluginParticipantInfo *participant_info,
            RTIBool top_level_registration, 
            void *container_plugin_context,
            RTICdrTypeCode *typeCode);

        NDDSUSERDllExport extern void 
        CommandInitialTypePlugin_on_participant_detached(
            PRESTypePluginParticipantData participant_data);

        NDDSUSERDllExport extern PRESTypePluginEndpointData 
        CommandInitialTypePlugin_on_endpoint_attached(
            PRESTypePluginParticipantData participant_data,
            const struct PRESTypePluginEndpointInfo *endpoint_info,
            RTIBool top_level_registration, 
            void *container_plugin_context);

        NDDSUSERDllExport extern void 
        CommandInitialTypePlugin_on_endpoint_detached(
            PRESTypePluginEndpointData endpoint_data);

        NDDSUSERDllExport extern void    
        CommandInitialTypePlugin_return_sample(
            PRESTypePluginEndpointData endpoint_data,
            CommandInitialType *sample,
            void *handle);    

        NDDSUSERDllExport extern RTIBool 
        CommandInitialTypePlugin_copy_sample(
            PRESTypePluginEndpointData endpoint_data,
            CommandInitialType *out,
            const CommandInitialType *in);

        /* ----------------------------------------------------------------------------
        (De)Serialize functions:
        * ------------------------------------------------------------------------- */

        NDDSUSERDllExport extern RTIBool
        CommandInitialTypePlugin_serialize_to_cdr_buffer(
            char * buffer,
            unsigned int * length,
            const CommandInitialType *sample,
            ::dds::core::policy::DataRepresentationId representation
            = ::dds::core::policy::DataRepresentation::xcdr()); 

        NDDSUSERDllExport extern RTIBool 
        CommandInitialTypePlugin_deserialize(
            PRESTypePluginEndpointData endpoint_data,
            CommandInitialType **sample, 
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_sample, 
            void *endpoint_plugin_qos);

        NDDSUSERDllExport extern RTIBool
        CommandInitialTypePlugin_deserialize_from_cdr_buffer(
            CommandInitialType *sample,
            const char * buffer,
            unsigned int length);    

        NDDSUSERDllExport extern unsigned int 
        CommandInitialTypePlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        /* --------------------------------------------------------------------------------------
        Key Management functions:
        * -------------------------------------------------------------------------------------- */
        NDDSUSERDllExport extern PRESTypePluginKeyKind 
        CommandInitialTypePlugin_get_key_kind(void);

        NDDSUSERDllExport extern unsigned int 
        CommandInitialTypePlugin_get_serialized_key_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        NDDSUSERDllExport extern unsigned int 
        CommandInitialTypePlugin_get_serialized_key_max_size_for_keyhash(
            PRESTypePluginEndpointData endpoint_data,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        NDDSUSERDllExport extern RTIBool 
        CommandInitialTypePlugin_deserialize_key(
            PRESTypePluginEndpointData endpoint_data,
            CommandInitialType ** sample,
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_key,
            void *endpoint_plugin_qos);

        /* Plugin Functions */
        NDDSUSERDllExport extern struct PRESTypePlugin*
        CommandInitialTypePlugin_new(void);

        NDDSUSERDllExport extern void
        CommandInitialTypePlugin_delete(struct PRESTypePlugin *);

    } /* namespace C2  */
} /* namespace ULARS  */

#if (defined(RTI_WIN32) || defined (RTI_WINCE) || defined(RTI_INTIME)) && defined(NDDS_USER_DLL_EXPORT)
/* If the code is building on Windows, stop exporting symbols.
*/
#undef NDDSUSERDllExport
#define NDDSUSERDllExport
#endif

#endif /* pilsPlugin_475505405_h */

