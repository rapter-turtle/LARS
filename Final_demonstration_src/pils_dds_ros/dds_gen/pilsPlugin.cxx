

/*
WARNING: THIS FILE IS AUTO-GENERATED. DO NOT MODIFY.

This file was generated from pils.idl
using RTI Code Generator (rtiddsgen) version 4.3.0.
The rtiddsgen tool is part of the RTI Connext DDS distribution.
For more information, type 'rtiddsgen -help' at a command shell
or consult the Code Generator User's Manual.
*/

#include <string.h>

#ifndef ndds_c_h
#include "ndds/ndds_c.h"
#endif

#ifndef osapi_type_h
#include "osapi/osapi_type.h"
#endif
#ifndef osapi_heap_h
#include "osapi/osapi_heap.h"
#endif

#ifndef osapi_utility_h
#include "osapi/osapi_utility.h"
#endif

#ifndef cdr_log_h
#include "cdr/cdr_log.h"
#endif

#ifndef cdr_type_h
#include "cdr/cdr_type.h"
#endif

#ifndef cdr_type_object_h
#include "cdr/cdr_typeObject.h"
#endif

#ifndef cdr_encapsulation_h
#include "cdr/cdr_encapsulation.h"
#endif

#ifndef cdr_stream_h
#include "cdr/cdr_stream.h"
#endif

#ifndef cdr_log_h
#include "cdr/cdr_log.h"
#endif

#ifndef pres_typePlugin_h
#include "pres/pres_typePlugin.h"
#endif

#include "dds_c/dds_c_typecode_impl.h"

#include "rti/topic/cdr/Serialization.hpp"

#define RTI_CDR_CURRENT_SUBMODULE RTI_CDR_SUBMODULE_MASK_STREAM

#include "pilsPlugin.hpp"

namespace ULARS {
    namespace EO {

        /* ----------------------------------------------------------------------------
        *  Type USVControlType
        * -------------------------------------------------------------------------- */

        /* -----------------------------------------------------------------------------
        Support functions:
        * -------------------------------------------------------------------------- */

        USVControlType *
        USVControlTypePluginSupport_create_data(void)
        {
            try {
                USVControlType *sample = new USVControlType();
                ::rti::topic::allocate_sample(*sample);
                return sample;
            } catch (...) {
                return NULL;
            }
        }

        void 
        USVControlTypePluginSupport_destroy_data(
            USVControlType *sample) 
        {
            delete sample;
        }

        RTIBool 
        USVControlTypePluginSupport_copy_data(
            USVControlType *dst,
            const USVControlType *src)
        {
            try {
                *dst = *src;
            } catch (...) {
                return RTI_FALSE;
            }

            return RTI_TRUE;
        }

        /* ----------------------------------------------------------------------------
        Callback functions:
        * ---------------------------------------------------------------------------- */

        PRESTypePluginParticipantData 
        USVControlTypePlugin_on_participant_attached(
            void *registration_data,
            const struct PRESTypePluginParticipantInfo *participant_info,
            RTIBool top_level_registration,
            void *container_plugin_context,
            RTICdrTypeCode *type_code)
        {
            try {
                struct RTIXCdrInterpreterPrograms *programs = NULL;
                struct PRESTypePluginDefaultParticipantData *pd = NULL;
                struct RTIXCdrInterpreterProgramsGenProperty programProperty =
                RTIXCdrInterpreterProgramsGenProperty_INITIALIZER;
                if (registration_data) {} /* To avoid warnings */
                if (participant_info) {} /* To avoid warnings */
                if (top_level_registration) {} /* To avoid warnings */
                if (container_plugin_context) {} /* To avoid warnings */
                if (type_code) {} /* To avoid warnings */
                pd = (struct PRESTypePluginDefaultParticipantData *)
                PRESTypePluginDefaultParticipantData_new(participant_info);

                programProperty.generateV1Encapsulation = RTI_XCDR_TRUE;
                programProperty.generateV2Encapsulation = RTI_XCDR_TRUE;
                programProperty.resolveAlias = RTI_XCDR_TRUE;
                programProperty.inlineStruct = RTI_XCDR_TRUE;
                programProperty.optimizeEnum = RTI_XCDR_TRUE;
                programProperty.unboundedSize = RTIXCdrLong_MAX;

                programProperty.externalReferenceSize = 
                (RTIXCdrUnsignedShort) sizeof(::dds::core::external<char>);
                programProperty.getExternalRefPointerFcn = 
                ::rti::topic::interpreter::get_external_value_pointer;

                programs = DDS_TypeCodeFactory_assert_programs_in_global_list(
                    DDS_TypeCodeFactory_get_instance(),
                    (DDS_TypeCode *) (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< USVControlType >::get().native()
                    ,
                    &programProperty,
                    RTI_XCDR_PROGRAM_MASK_TYPEPLUGIN);

                if (programs == NULL) {
                    PRESTypePluginDefaultParticipantData_delete(
                        (PRESTypePluginParticipantData)pd);
                    return NULL;
                }

                pd->programs = programs;
                return (PRESTypePluginParticipantData)pd;
            } catch (...) {
                return NULL;
            }
        }

        void 
        USVControlTypePlugin_on_participant_detached(
            PRESTypePluginParticipantData participant_data)
        {
            if (participant_data != NULL) {
                struct PRESTypePluginDefaultParticipantData *pd = 
                (struct PRESTypePluginDefaultParticipantData *)participant_data;

                if (pd->programs != NULL) {
                    DDS_TypeCodeFactory_remove_programs_from_global_list(
                        DDS_TypeCodeFactory_get_instance(),
                        pd->programs);
                    pd->programs = NULL;
                }
                PRESTypePluginDefaultParticipantData_delete(participant_data);
            }
        }

        PRESTypePluginEndpointData
        USVControlTypePlugin_on_endpoint_attached(
            PRESTypePluginParticipantData participant_data,
            const struct PRESTypePluginEndpointInfo *endpoint_info,
            RTIBool top_level_registration, 
            void *containerPluginContext)
        {
            try {
                PRESTypePluginEndpointData epd = NULL;
                unsigned int serializedSampleMaxSize = 0;

                if (top_level_registration) {} /* To avoid warnings */
                if (containerPluginContext) {} /* To avoid warnings */

                if (participant_data == NULL) {
                    return NULL;
                } 

                epd = PRESTypePluginDefaultEndpointData_new(
                    participant_data,
                    endpoint_info,
                    (PRESTypePluginDefaultEndpointDataCreateSampleFunction)
                    USVControlTypePluginSupport_create_data,
                    (PRESTypePluginDefaultEndpointDataDestroySampleFunction)
                    USVControlTypePluginSupport_destroy_data,
                    NULL , NULL );

                if (epd == NULL) {
                    return NULL;
                } 

                if (endpoint_info->endpointKind == PRES_TYPEPLUGIN_ENDPOINT_WRITER) {
                    serializedSampleMaxSize = ::ULARS::EO::USVControlTypePlugin_get_serialized_sample_max_size(
                        epd,RTI_FALSE,RTI_CDR_ENCAPSULATION_ID_CDR_BE,0);
                    PRESTypePluginDefaultEndpointData_setMaxSizeSerializedSample(epd, serializedSampleMaxSize);

                    if (PRESTypePluginDefaultEndpointData_createWriterPool(
                        epd,
                        endpoint_info,
                        (PRESTypePluginGetSerializedSampleMaxSizeFunction)
                        ::ULARS::EO::USVControlTypePlugin_get_serialized_sample_max_size, epd,
                        (PRESTypePluginGetSerializedSampleSizeFunction)
                        PRESTypePlugin_interpretedGetSerializedSampleSize,
                        epd) == RTI_FALSE) {
                        PRESTypePluginDefaultEndpointData_delete(epd);
                        return NULL;
                    }
                }

                return epd;
            } catch (...) {
                return NULL;
            }
        }

        void 
        USVControlTypePlugin_on_endpoint_detached(
            PRESTypePluginEndpointData endpoint_data)
        {  
            PRESTypePluginDefaultEndpointData_delete(endpoint_data);
        }

        void    
        USVControlTypePlugin_return_sample(
            PRESTypePluginEndpointData endpoint_data,
            USVControlType *sample,
            void *handle)
        {
            try {
                ::rti::topic::reset_sample(*sample);
            } catch(const std::exception& ex) {
                RTICdrLog_logWithFunctionName(
                    RTI_LOG_BIT_EXCEPTION,
                    "USVControlTypePlugin_return_sample",
                    &RTI_LOG_ANY_FAILURE_ss,
                    "exception: ",
                    ex.what());
            }

            PRESTypePluginDefaultEndpointData_returnSample(
                endpoint_data, sample, handle);
        }

        RTIBool 
        USVControlTypePlugin_copy_sample(
            PRESTypePluginEndpointData,
            USVControlType *dst,
            const USVControlType *src)
        {
            return ::ULARS::EO::USVControlTypePluginSupport_copy_data(dst,src);
        }

        /* ----------------------------------------------------------------------------
        (De)Serialize functions:
        * ------------------------------------------------------------------------- */
        unsigned int 
        USVControlTypePlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        RTIBool
        USVControlTypePlugin_serialize_to_cdr_buffer(
            char * buffer,
            unsigned int * length,
            const USVControlType *sample,
            ::dds::core::policy::DataRepresentationId representation)
        {
            using namespace ::dds::core::policy;

            try{
                RTIEncapsulationId encapsulationId = RTI_CDR_ENCAPSULATION_ID_INVALID;
                struct RTICdrStream cdrStream;
                struct PRESTypePluginDefaultEndpointData epd;
                RTIBool result;
                struct PRESTypePluginDefaultParticipantData pd;
                struct RTIXCdrTypePluginProgramContext defaultProgramContext =
                RTIXCdrTypePluginProgramContext_INTIALIZER;
                struct PRESTypePlugin plugin = PRES_TYPEPLUGIN_DEFAULT;

                if (length == NULL) {
                    return RTI_FALSE;
                }

                RTIOsapiMemory_zero(&epd, sizeof(struct PRESTypePluginDefaultEndpointData));
                epd.programContext = defaultProgramContext;
                epd._participantData = &pd;
                epd.typePlugin = &plugin;
                epd.programContext.endpointPluginData = &epd;
                try {
                    plugin.typeCode = (struct RTICdrTypeCode *)
                    (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< USVControlType >::get().native()
                    ;
                } catch (...) {
                    return RTI_FALSE;
                }
                pd.programs = ::rti::topic::interpreter::get_cdr_serialization_programs<
                USVControlType, 
                true, true, true>();

                encapsulationId = DDS_TypeCode_get_native_encapsulation(
                    (DDS_TypeCode *) plugin.typeCode,
                    representation);

                if (encapsulationId == RTI_CDR_ENCAPSULATION_ID_INVALID) {
                    return RTI_FALSE;
                }

                epd._maxSizeSerializedSample =
                USVControlTypePlugin_get_serialized_sample_max_size(
                    (PRESTypePluginEndpointData)&epd, 
                    RTI_TRUE, 
                    encapsulationId,
                    0);

                if (buffer == NULL) {
                    *length = 
                    PRESTypePlugin_interpretedGetSerializedSampleSize(
                        (PRESTypePluginEndpointData)&epd,
                        RTI_TRUE,
                        encapsulationId,
                        0,
                        sample);

                    if (*length == 0) {
                        return RTI_FALSE;
                    }

                    return RTI_TRUE;
                }    

                RTICdrStream_init(&cdrStream);
                RTICdrStream_set(&cdrStream, buffer, *length);

                result = PRESTypePlugin_interpretedSerialize(
                    (PRESTypePluginEndpointData)&epd, 
                    sample, 
                    &cdrStream, 
                    RTI_TRUE, 
                    encapsulationId,
                    RTI_TRUE, 
                    NULL);  

                *length = (unsigned int) RTICdrStream_getCurrentPositionOffset(&cdrStream);
                return result;
            } catch (...) {
                return RTI_FALSE;
            }
        }

        RTIBool
        USVControlTypePlugin_deserialize_from_cdr_buffer(
            USVControlType *sample,
            const char * buffer,
            unsigned int length)
        {
            struct RTICdrStream cdrStream;
            struct PRESTypePluginDefaultParticipantData pd;
            struct RTIXCdrTypePluginProgramContext defaultProgramContext =
            RTIXCdrTypePluginProgramContext_INTIALIZER;
            struct PRESTypePlugin plugin;
            struct PRESTypePluginDefaultEndpointData epd;

            RTICdrStream_init(&cdrStream);
            /*
            * The buffer is constant because this is a deserialization function
            * (the buffer is an input parameter, not an output parameter).
            * However, the buffer member in the stream is a (char *) so coverity
            * complains in case something else modifies the buffer's contents later.
            *
            * We don't currently have a stream type with a constant buffer.
            * Therefore, we suppress the warning after making sure that this function
            * doesn't modify the contents of the stream's buffer.
            */
            /* coverity[cert_exp40_c_violation : FALSE] */
            RTICdrStream_set(&cdrStream, (char *) buffer, length);

            epd.programContext = defaultProgramContext;
            epd._participantData = &pd;
            epd.typePlugin = &plugin;
            epd.programContext.endpointPluginData = &epd;
            try {
                plugin.typeCode = (struct RTICdrTypeCode *)
                (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< USVControlType >::get().native()
                ;
            } catch (...) {
                return RTI_FALSE;
            }
            pd.programs = ::rti::topic::interpreter::get_cdr_serialization_programs<
            USVControlType, 
            true, true, true>();

            epd._assignabilityProperty.acceptUnknownEnumValue = RTI_XCDR_TRUE;
            epd._assignabilityProperty.acceptUnknownUnionDiscriminator = 
            RTI_XCDR_ACCEPT_UNKNOWN_DISCRIMINATOR_AND_SELECT_DEFAULT;

            ::rti::topic::reset_sample(*sample);
            return PRESTypePlugin_interpretedDeserialize( 
                (PRESTypePluginEndpointData)&epd,
                sample,
                &cdrStream,
                RTI_TRUE, 
                RTI_TRUE, 
                NULL);
        }

        unsigned int 
        USVControlTypePlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            try {
                unsigned int size;
                RTIBool overflow = RTI_FALSE;

                size = PRESTypePlugin_interpretedGetSerializedSampleMaxSize(
                    endpoint_data,&overflow,include_encapsulation,encapsulation_id,current_alignment);

                if (overflow) {
                    size = RTI_CDR_MAX_SERIALIZED_SIZE;
                }

                return size;
            } catch (...) {
                return 0;
            }
        }

        /* --------------------------------------------------------------------------------------
        Key Management functions:
        * -------------------------------------------------------------------------------------- */

        PRESTypePluginKeyKind 
        USVControlTypePlugin_get_key_kind(void)
        {
            return PRES_TYPEPLUGIN_NO_KEY;
        }

        RTIBool USVControlTypePlugin_deserialize_key(
            PRESTypePluginEndpointData endpoint_data,
            USVControlType **sample, 
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_key,
            void *endpoint_plugin_qos)
        {
            try {
                RTIBool result;
                if (drop_sample) {} /* To avoid warnings */
                cdrStream->_xTypesState.unassignable = RTI_FALSE;
                result= PRESTypePlugin_interpretedDeserializeKey( 
                    endpoint_data, (sample != NULL)?*sample:NULL, cdrStream,
                    deserialize_encapsulation, deserialize_key, endpoint_plugin_qos);
                if (result) {
                    if (cdrStream->_xTypesState.unassignable) {
                        result = RTI_FALSE;
                    }
                }
                return result;    
            } catch (...) {
                return RTI_FALSE;
            }     
        }

        unsigned int
        USVControlTypePlugin_get_serialized_key_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            try {
                unsigned int size;
                RTIBool overflow = RTI_FALSE;

                size = PRESTypePlugin_interpretedGetSerializedKeyMaxSize(
                    endpoint_data,&overflow,include_encapsulation,encapsulation_id,current_alignment);
                if (overflow) {
                    size = RTI_CDR_MAX_SERIALIZED_SIZE;
                }

                return size;
            } catch (...) {
                return RTI_FALSE;
            }    
        }

        unsigned int
        USVControlTypePlugin_get_serialized_key_max_size_for_keyhash(
            PRESTypePluginEndpointData endpoint_data,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            unsigned int size;
            RTIBool overflow = RTI_FALSE;

            size = PRESTypePlugin_interpretedGetSerializedKeyMaxSizeForKeyhash(
                endpoint_data,
                &overflow,
                encapsulation_id,
                current_alignment);
            if (overflow) {
                size = RTI_CDR_MAX_SERIALIZED_SIZE;
            }

            return size;
        }

        /* ------------------------------------------------------------------------
        * Plug-in Installation Methods
        * ------------------------------------------------------------------------ */
        struct PRESTypePlugin *USVControlTypePlugin_new(void) 
        { 
            struct PRESTypePlugin *plugin = NULL;
            const struct PRESTypePluginVersion PLUGIN_VERSION = 
            PRES_TYPE_PLUGIN_VERSION_2_0;

            RTIOsapiHeap_allocateStructure(
                &plugin, struct PRESTypePlugin);
            if (plugin == NULL) {
                return NULL;
            }

            plugin->version = PLUGIN_VERSION;

            /* set up parent's function pointers */
            plugin->onParticipantAttached =
            (PRESTypePluginOnParticipantAttachedCallback)
            ::ULARS::EO::USVControlTypePlugin_on_participant_attached;
            plugin->onParticipantDetached =
            (PRESTypePluginOnParticipantDetachedCallback)
            ::ULARS::EO::USVControlTypePlugin_on_participant_detached;
            plugin->onEndpointAttached =
            (PRESTypePluginOnEndpointAttachedCallback)
            ::ULARS::EO::USVControlTypePlugin_on_endpoint_attached;
            plugin->onEndpointDetached =
            (PRESTypePluginOnEndpointDetachedCallback)
            ::ULARS::EO::USVControlTypePlugin_on_endpoint_detached;

            plugin->copySampleFnc =
            (PRESTypePluginCopySampleFunction)
            ::ULARS::EO::USVControlTypePlugin_copy_sample;
            plugin->createSampleFnc =
            (PRESTypePluginCreateSampleFunction)
            USVControlTypePlugin_create_sample;
            plugin->destroySampleFnc =
            (PRESTypePluginDestroySampleFunction)
            USVControlTypePlugin_destroy_sample;

            plugin->serializeFnc = 
            (PRESTypePluginSerializeFunction) PRESTypePlugin_interpretedSerialize;
            plugin->deserializeFnc =
            (PRESTypePluginDeserializeFunction) PRESTypePlugin_interpretedDeserializeWithAlloc;
            plugin->getSerializedSampleMaxSizeFnc =
            (PRESTypePluginGetSerializedSampleMaxSizeFunction)
            ::ULARS::EO::USVControlTypePlugin_get_serialized_sample_max_size;
            plugin->getSerializedSampleMinSizeFnc =
            (PRESTypePluginGetSerializedSampleMinSizeFunction)
            PRESTypePlugin_interpretedGetSerializedSampleMinSize;
            plugin->getDeserializedSampleMaxSizeFnc = NULL; 
            plugin->getSampleFnc =
            (PRESTypePluginGetSampleFunction)
            USVControlTypePlugin_get_sample;
            plugin->returnSampleFnc =
            (PRESTypePluginReturnSampleFunction)
            USVControlTypePlugin_return_sample;
            plugin->getKeyKindFnc =
            (PRESTypePluginGetKeyKindFunction)
            ::ULARS::EO::USVControlTypePlugin_get_key_kind;

            /* These functions are only used for keyed types. As this is not a keyed
            type they are all set to NULL
            */
            plugin->serializeKeyFnc = NULL ;    
            plugin->deserializeKeyFnc = NULL;  
            plugin->getKeyFnc = NULL;
            plugin->returnKeyFnc = NULL;
            plugin->instanceToKeyFnc = NULL;
            plugin->keyToInstanceFnc = NULL;
            plugin->getSerializedKeyMaxSizeFnc = NULL;
            plugin->instanceToKeyHashFnc = NULL;
            plugin->serializedSampleToKeyHashFnc = NULL;
            plugin->serializedKeyToKeyHashFnc = NULL;    
            #ifdef NDDS_STANDALONE_TYPE
            plugin->typeCode = NULL; 
            #else
            try {
                plugin->typeCode = (struct RTICdrTypeCode *)
                &::rti::topic::dynamic_type< ::ULARS::EO::USVControlType >::get().native();
            } catch (...) {
                ::ULARS::EO::USVControlTypePlugin_delete(plugin);
                return NULL;
            }
            #endif
            plugin->languageKind = PRES_TYPEPLUGIN_CPPSTL_LANG;

            /* Serialized buffer */
            plugin->getBuffer = 
            (PRESTypePluginGetBufferFunction)
            USVControlTypePlugin_get_buffer;
            plugin->returnBuffer = 
            (PRESTypePluginReturnBufferFunction)
            USVControlTypePlugin_return_buffer;
            plugin->getBufferWithParams = NULL;
            plugin->returnBufferWithParams = NULL;
            plugin->getSerializedSampleSizeFnc =
            (PRESTypePluginGetSerializedSampleSizeFunction)
            PRESTypePlugin_interpretedGetSerializedSampleSize;

            plugin->getWriterLoanedSampleFnc = NULL; 
            plugin->returnWriterLoanedSampleFnc = NULL;
            plugin->returnWriterLoanedSampleFromCookieFnc = NULL;
            plugin->validateWriterLoanedSampleFnc = NULL;
            plugin->setWriterLoanedSampleSerializedStateFnc = NULL;

            static const char * TYPE_NAME = "ULARS::EO::USVControlType";
            plugin->endpointTypeName = TYPE_NAME;
            plugin->isMetpType = RTI_FALSE;
            return plugin;
        }

        void
        USVControlTypePlugin_delete(struct PRESTypePlugin *plugin)
        {
            RTIOsapiHeap_freeStructure(plugin);
        } 
    } /* namespace EO  */
    namespace C2 {

        /* ----------------------------------------------------------------------------
        *  Type Waypoint
        * -------------------------------------------------------------------------- */

        /* -----------------------------------------------------------------------------
        Support functions:
        * -------------------------------------------------------------------------- */

        Waypoint *
        WaypointPluginSupport_create_data(void)
        {
            try {
                Waypoint *sample = new Waypoint();
                ::rti::topic::allocate_sample(*sample);
                return sample;
            } catch (...) {
                return NULL;
            }
        }

        void 
        WaypointPluginSupport_destroy_data(
            Waypoint *sample) 
        {
            delete sample;
        }

        RTIBool 
        WaypointPluginSupport_copy_data(
            Waypoint *dst,
            const Waypoint *src)
        {
            try {
                *dst = *src;
            } catch (...) {
                return RTI_FALSE;
            }

            return RTI_TRUE;
        }

        /* ----------------------------------------------------------------------------
        Callback functions:
        * ---------------------------------------------------------------------------- */

        PRESTypePluginParticipantData 
        WaypointPlugin_on_participant_attached(
            void *registration_data,
            const struct PRESTypePluginParticipantInfo *participant_info,
            RTIBool top_level_registration,
            void *container_plugin_context,
            RTICdrTypeCode *type_code)
        {
            try {
                struct RTIXCdrInterpreterPrograms *programs = NULL;
                struct PRESTypePluginDefaultParticipantData *pd = NULL;
                struct RTIXCdrInterpreterProgramsGenProperty programProperty =
                RTIXCdrInterpreterProgramsGenProperty_INITIALIZER;
                if (registration_data) {} /* To avoid warnings */
                if (participant_info) {} /* To avoid warnings */
                if (top_level_registration) {} /* To avoid warnings */
                if (container_plugin_context) {} /* To avoid warnings */
                if (type_code) {} /* To avoid warnings */
                pd = (struct PRESTypePluginDefaultParticipantData *)
                PRESTypePluginDefaultParticipantData_new(participant_info);

                programProperty.generateV1Encapsulation = RTI_XCDR_TRUE;
                programProperty.generateV2Encapsulation = RTI_XCDR_TRUE;
                programProperty.resolveAlias = RTI_XCDR_TRUE;
                programProperty.inlineStruct = RTI_XCDR_TRUE;
                programProperty.optimizeEnum = RTI_XCDR_TRUE;
                programProperty.unboundedSize = RTIXCdrLong_MAX;

                programProperty.externalReferenceSize = 
                (RTIXCdrUnsignedShort) sizeof(::dds::core::external<char>);
                programProperty.getExternalRefPointerFcn = 
                ::rti::topic::interpreter::get_external_value_pointer;

                programs = DDS_TypeCodeFactory_assert_programs_in_global_list(
                    DDS_TypeCodeFactory_get_instance(),
                    (DDS_TypeCode *) (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< Waypoint >::get().native()
                    ,
                    &programProperty,
                    RTI_XCDR_PROGRAM_MASK_TYPEPLUGIN);

                if (programs == NULL) {
                    PRESTypePluginDefaultParticipantData_delete(
                        (PRESTypePluginParticipantData)pd);
                    return NULL;
                }

                pd->programs = programs;
                return (PRESTypePluginParticipantData)pd;
            } catch (...) {
                return NULL;
            }
        }

        void 
        WaypointPlugin_on_participant_detached(
            PRESTypePluginParticipantData participant_data)
        {
            if (participant_data != NULL) {
                struct PRESTypePluginDefaultParticipantData *pd = 
                (struct PRESTypePluginDefaultParticipantData *)participant_data;

                if (pd->programs != NULL) {
                    DDS_TypeCodeFactory_remove_programs_from_global_list(
                        DDS_TypeCodeFactory_get_instance(),
                        pd->programs);
                    pd->programs = NULL;
                }
                PRESTypePluginDefaultParticipantData_delete(participant_data);
            }
        }

        PRESTypePluginEndpointData
        WaypointPlugin_on_endpoint_attached(
            PRESTypePluginParticipantData participant_data,
            const struct PRESTypePluginEndpointInfo *endpoint_info,
            RTIBool top_level_registration, 
            void *containerPluginContext)
        {
            try {
                PRESTypePluginEndpointData epd = NULL;
                unsigned int serializedSampleMaxSize = 0;

                if (top_level_registration) {} /* To avoid warnings */
                if (containerPluginContext) {} /* To avoid warnings */

                if (participant_data == NULL) {
                    return NULL;
                } 

                epd = PRESTypePluginDefaultEndpointData_new(
                    participant_data,
                    endpoint_info,
                    (PRESTypePluginDefaultEndpointDataCreateSampleFunction)
                    WaypointPluginSupport_create_data,
                    (PRESTypePluginDefaultEndpointDataDestroySampleFunction)
                    WaypointPluginSupport_destroy_data,
                    NULL , NULL );

                if (epd == NULL) {
                    return NULL;
                } 

                if (endpoint_info->endpointKind == PRES_TYPEPLUGIN_ENDPOINT_WRITER) {
                    serializedSampleMaxSize = ::ULARS::C2::WaypointPlugin_get_serialized_sample_max_size(
                        epd,RTI_FALSE,RTI_CDR_ENCAPSULATION_ID_CDR_BE,0);
                    PRESTypePluginDefaultEndpointData_setMaxSizeSerializedSample(epd, serializedSampleMaxSize);

                    if (PRESTypePluginDefaultEndpointData_createWriterPool(
                        epd,
                        endpoint_info,
                        (PRESTypePluginGetSerializedSampleMaxSizeFunction)
                        ::ULARS::C2::WaypointPlugin_get_serialized_sample_max_size, epd,
                        (PRESTypePluginGetSerializedSampleSizeFunction)
                        PRESTypePlugin_interpretedGetSerializedSampleSize,
                        epd) == RTI_FALSE) {
                        PRESTypePluginDefaultEndpointData_delete(epd);
                        return NULL;
                    }
                }

                return epd;
            } catch (...) {
                return NULL;
            }
        }

        void 
        WaypointPlugin_on_endpoint_detached(
            PRESTypePluginEndpointData endpoint_data)
        {  
            PRESTypePluginDefaultEndpointData_delete(endpoint_data);
        }

        void    
        WaypointPlugin_return_sample(
            PRESTypePluginEndpointData endpoint_data,
            Waypoint *sample,
            void *handle)
        {
            try {
                ::rti::topic::reset_sample(*sample);
            } catch(const std::exception& ex) {
                RTICdrLog_logWithFunctionName(
                    RTI_LOG_BIT_EXCEPTION,
                    "WaypointPlugin_return_sample",
                    &RTI_LOG_ANY_FAILURE_ss,
                    "exception: ",
                    ex.what());
            }

            PRESTypePluginDefaultEndpointData_returnSample(
                endpoint_data, sample, handle);
        }

        RTIBool 
        WaypointPlugin_copy_sample(
            PRESTypePluginEndpointData,
            Waypoint *dst,
            const Waypoint *src)
        {
            return ::ULARS::C2::WaypointPluginSupport_copy_data(dst,src);
        }

        /* ----------------------------------------------------------------------------
        (De)Serialize functions:
        * ------------------------------------------------------------------------- */
        unsigned int 
        WaypointPlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        RTIBool
        WaypointPlugin_serialize_to_cdr_buffer(
            char * buffer,
            unsigned int * length,
            const Waypoint *sample,
            ::dds::core::policy::DataRepresentationId representation)
        {
            using namespace ::dds::core::policy;

            try{
                RTIEncapsulationId encapsulationId = RTI_CDR_ENCAPSULATION_ID_INVALID;
                struct RTICdrStream cdrStream;
                struct PRESTypePluginDefaultEndpointData epd;
                RTIBool result;
                struct PRESTypePluginDefaultParticipantData pd;
                struct RTIXCdrTypePluginProgramContext defaultProgramContext =
                RTIXCdrTypePluginProgramContext_INTIALIZER;
                struct PRESTypePlugin plugin = PRES_TYPEPLUGIN_DEFAULT;

                if (length == NULL) {
                    return RTI_FALSE;
                }

                RTIOsapiMemory_zero(&epd, sizeof(struct PRESTypePluginDefaultEndpointData));
                epd.programContext = defaultProgramContext;
                epd._participantData = &pd;
                epd.typePlugin = &plugin;
                epd.programContext.endpointPluginData = &epd;
                try {
                    plugin.typeCode = (struct RTICdrTypeCode *)
                    (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< Waypoint >::get().native()
                    ;
                } catch (...) {
                    return RTI_FALSE;
                }
                pd.programs = ::rti::topic::interpreter::get_cdr_serialization_programs<
                Waypoint, 
                true, true, true>();

                encapsulationId = DDS_TypeCode_get_native_encapsulation(
                    (DDS_TypeCode *) plugin.typeCode,
                    representation);

                if (encapsulationId == RTI_CDR_ENCAPSULATION_ID_INVALID) {
                    return RTI_FALSE;
                }

                epd._maxSizeSerializedSample =
                WaypointPlugin_get_serialized_sample_max_size(
                    (PRESTypePluginEndpointData)&epd, 
                    RTI_TRUE, 
                    encapsulationId,
                    0);

                if (buffer == NULL) {
                    *length = 
                    PRESTypePlugin_interpretedGetSerializedSampleSize(
                        (PRESTypePluginEndpointData)&epd,
                        RTI_TRUE,
                        encapsulationId,
                        0,
                        sample);

                    if (*length == 0) {
                        return RTI_FALSE;
                    }

                    return RTI_TRUE;
                }    

                RTICdrStream_init(&cdrStream);
                RTICdrStream_set(&cdrStream, buffer, *length);

                result = PRESTypePlugin_interpretedSerialize(
                    (PRESTypePluginEndpointData)&epd, 
                    sample, 
                    &cdrStream, 
                    RTI_TRUE, 
                    encapsulationId,
                    RTI_TRUE, 
                    NULL);  

                *length = (unsigned int) RTICdrStream_getCurrentPositionOffset(&cdrStream);
                return result;
            } catch (...) {
                return RTI_FALSE;
            }
        }

        RTIBool
        WaypointPlugin_deserialize_from_cdr_buffer(
            Waypoint *sample,
            const char * buffer,
            unsigned int length)
        {
            struct RTICdrStream cdrStream;
            struct PRESTypePluginDefaultParticipantData pd;
            struct RTIXCdrTypePluginProgramContext defaultProgramContext =
            RTIXCdrTypePluginProgramContext_INTIALIZER;
            struct PRESTypePlugin plugin;
            struct PRESTypePluginDefaultEndpointData epd;

            RTICdrStream_init(&cdrStream);
            /*
            * The buffer is constant because this is a deserialization function
            * (the buffer is an input parameter, not an output parameter).
            * However, the buffer member in the stream is a (char *) so coverity
            * complains in case something else modifies the buffer's contents later.
            *
            * We don't currently have a stream type with a constant buffer.
            * Therefore, we suppress the warning after making sure that this function
            * doesn't modify the contents of the stream's buffer.
            */
            /* coverity[cert_exp40_c_violation : FALSE] */
            RTICdrStream_set(&cdrStream, (char *) buffer, length);

            epd.programContext = defaultProgramContext;
            epd._participantData = &pd;
            epd.typePlugin = &plugin;
            epd.programContext.endpointPluginData = &epd;
            try {
                plugin.typeCode = (struct RTICdrTypeCode *)
                (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< Waypoint >::get().native()
                ;
            } catch (...) {
                return RTI_FALSE;
            }
            pd.programs = ::rti::topic::interpreter::get_cdr_serialization_programs<
            Waypoint, 
            true, true, true>();

            epd._assignabilityProperty.acceptUnknownEnumValue = RTI_XCDR_TRUE;
            epd._assignabilityProperty.acceptUnknownUnionDiscriminator = 
            RTI_XCDR_ACCEPT_UNKNOWN_DISCRIMINATOR_AND_SELECT_DEFAULT;

            ::rti::topic::reset_sample(*sample);
            return PRESTypePlugin_interpretedDeserialize( 
                (PRESTypePluginEndpointData)&epd,
                sample,
                &cdrStream,
                RTI_TRUE, 
                RTI_TRUE, 
                NULL);
        }

        unsigned int 
        WaypointPlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            try {
                unsigned int size;
                RTIBool overflow = RTI_FALSE;

                size = PRESTypePlugin_interpretedGetSerializedSampleMaxSize(
                    endpoint_data,&overflow,include_encapsulation,encapsulation_id,current_alignment);

                if (overflow) {
                    size = RTI_CDR_MAX_SERIALIZED_SIZE;
                }

                return size;
            } catch (...) {
                return 0;
            }
        }

        /* --------------------------------------------------------------------------------------
        Key Management functions:
        * -------------------------------------------------------------------------------------- */

        PRESTypePluginKeyKind 
        WaypointPlugin_get_key_kind(void)
        {
            return PRES_TYPEPLUGIN_NO_KEY;
        }

        RTIBool WaypointPlugin_deserialize_key(
            PRESTypePluginEndpointData endpoint_data,
            Waypoint **sample, 
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_key,
            void *endpoint_plugin_qos)
        {
            try {
                RTIBool result;
                if (drop_sample) {} /* To avoid warnings */
                cdrStream->_xTypesState.unassignable = RTI_FALSE;
                result= PRESTypePlugin_interpretedDeserializeKey( 
                    endpoint_data, (sample != NULL)?*sample:NULL, cdrStream,
                    deserialize_encapsulation, deserialize_key, endpoint_plugin_qos);
                if (result) {
                    if (cdrStream->_xTypesState.unassignable) {
                        result = RTI_FALSE;
                    }
                }
                return result;    
            } catch (...) {
                return RTI_FALSE;
            }     
        }

        unsigned int
        WaypointPlugin_get_serialized_key_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            try {
                unsigned int size;
                RTIBool overflow = RTI_FALSE;

                size = PRESTypePlugin_interpretedGetSerializedKeyMaxSize(
                    endpoint_data,&overflow,include_encapsulation,encapsulation_id,current_alignment);
                if (overflow) {
                    size = RTI_CDR_MAX_SERIALIZED_SIZE;
                }

                return size;
            } catch (...) {
                return RTI_FALSE;
            }    
        }

        unsigned int
        WaypointPlugin_get_serialized_key_max_size_for_keyhash(
            PRESTypePluginEndpointData endpoint_data,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            unsigned int size;
            RTIBool overflow = RTI_FALSE;

            size = PRESTypePlugin_interpretedGetSerializedKeyMaxSizeForKeyhash(
                endpoint_data,
                &overflow,
                encapsulation_id,
                current_alignment);
            if (overflow) {
                size = RTI_CDR_MAX_SERIALIZED_SIZE;
            }

            return size;
        }

        /* ------------------------------------------------------------------------
        * Plug-in Installation Methods
        * ------------------------------------------------------------------------ */
        struct PRESTypePlugin *WaypointPlugin_new(void) 
        { 
            struct PRESTypePlugin *plugin = NULL;
            const struct PRESTypePluginVersion PLUGIN_VERSION = 
            PRES_TYPE_PLUGIN_VERSION_2_0;

            RTIOsapiHeap_allocateStructure(
                &plugin, struct PRESTypePlugin);
            if (plugin == NULL) {
                return NULL;
            }

            plugin->version = PLUGIN_VERSION;

            /* set up parent's function pointers */
            plugin->onParticipantAttached =
            (PRESTypePluginOnParticipantAttachedCallback)
            ::ULARS::C2::WaypointPlugin_on_participant_attached;
            plugin->onParticipantDetached =
            (PRESTypePluginOnParticipantDetachedCallback)
            ::ULARS::C2::WaypointPlugin_on_participant_detached;
            plugin->onEndpointAttached =
            (PRESTypePluginOnEndpointAttachedCallback)
            ::ULARS::C2::WaypointPlugin_on_endpoint_attached;
            plugin->onEndpointDetached =
            (PRESTypePluginOnEndpointDetachedCallback)
            ::ULARS::C2::WaypointPlugin_on_endpoint_detached;

            plugin->copySampleFnc =
            (PRESTypePluginCopySampleFunction)
            ::ULARS::C2::WaypointPlugin_copy_sample;
            plugin->createSampleFnc =
            (PRESTypePluginCreateSampleFunction)
            WaypointPlugin_create_sample;
            plugin->destroySampleFnc =
            (PRESTypePluginDestroySampleFunction)
            WaypointPlugin_destroy_sample;

            plugin->serializeFnc = 
            (PRESTypePluginSerializeFunction) PRESTypePlugin_interpretedSerialize;
            plugin->deserializeFnc =
            (PRESTypePluginDeserializeFunction) PRESTypePlugin_interpretedDeserializeWithAlloc;
            plugin->getSerializedSampleMaxSizeFnc =
            (PRESTypePluginGetSerializedSampleMaxSizeFunction)
            ::ULARS::C2::WaypointPlugin_get_serialized_sample_max_size;
            plugin->getSerializedSampleMinSizeFnc =
            (PRESTypePluginGetSerializedSampleMinSizeFunction)
            PRESTypePlugin_interpretedGetSerializedSampleMinSize;
            plugin->getDeserializedSampleMaxSizeFnc = NULL; 
            plugin->getSampleFnc =
            (PRESTypePluginGetSampleFunction)
            WaypointPlugin_get_sample;
            plugin->returnSampleFnc =
            (PRESTypePluginReturnSampleFunction)
            WaypointPlugin_return_sample;
            plugin->getKeyKindFnc =
            (PRESTypePluginGetKeyKindFunction)
            ::ULARS::C2::WaypointPlugin_get_key_kind;

            /* These functions are only used for keyed types. As this is not a keyed
            type they are all set to NULL
            */
            plugin->serializeKeyFnc = NULL ;    
            plugin->deserializeKeyFnc = NULL;  
            plugin->getKeyFnc = NULL;
            plugin->returnKeyFnc = NULL;
            plugin->instanceToKeyFnc = NULL;
            plugin->keyToInstanceFnc = NULL;
            plugin->getSerializedKeyMaxSizeFnc = NULL;
            plugin->instanceToKeyHashFnc = NULL;
            plugin->serializedSampleToKeyHashFnc = NULL;
            plugin->serializedKeyToKeyHashFnc = NULL;    
            #ifdef NDDS_STANDALONE_TYPE
            plugin->typeCode = NULL; 
            #else
            try {
                plugin->typeCode = (struct RTICdrTypeCode *)
                &::rti::topic::dynamic_type< ::ULARS::C2::Waypoint >::get().native();
            } catch (...) {
                ::ULARS::C2::WaypointPlugin_delete(plugin);
                return NULL;
            }
            #endif
            plugin->languageKind = PRES_TYPEPLUGIN_CPPSTL_LANG;

            /* Serialized buffer */
            plugin->getBuffer = 
            (PRESTypePluginGetBufferFunction)
            WaypointPlugin_get_buffer;
            plugin->returnBuffer = 
            (PRESTypePluginReturnBufferFunction)
            WaypointPlugin_return_buffer;
            plugin->getBufferWithParams = NULL;
            plugin->returnBufferWithParams = NULL;
            plugin->getSerializedSampleSizeFnc =
            (PRESTypePluginGetSerializedSampleSizeFunction)
            PRESTypePlugin_interpretedGetSerializedSampleSize;

            plugin->getWriterLoanedSampleFnc = NULL; 
            plugin->returnWriterLoanedSampleFnc = NULL;
            plugin->returnWriterLoanedSampleFromCookieFnc = NULL;
            plugin->validateWriterLoanedSampleFnc = NULL;
            plugin->setWriterLoanedSampleSerializedStateFnc = NULL;

            static const char * TYPE_NAME = "ULARS::C2::Waypoint";
            plugin->endpointTypeName = TYPE_NAME;
            plugin->isMetpType = RTI_FALSE;
            return plugin;
        }

        void
        WaypointPlugin_delete(struct PRESTypePlugin *plugin)
        {
            RTIOsapiHeap_freeStructure(plugin);
        } 

        /* ----------------------------------------------------------------------------
        *  Type GlobalWaypointType
        * -------------------------------------------------------------------------- */

        /* -----------------------------------------------------------------------------
        Support functions:
        * -------------------------------------------------------------------------- */

        GlobalWaypointType *
        GlobalWaypointTypePluginSupport_create_data(void)
        {
            try {
                GlobalWaypointType *sample = new GlobalWaypointType();
                ::rti::topic::allocate_sample(*sample);
                return sample;
            } catch (...) {
                return NULL;
            }
        }

        void 
        GlobalWaypointTypePluginSupport_destroy_data(
            GlobalWaypointType *sample) 
        {
            delete sample;
        }

        RTIBool 
        GlobalWaypointTypePluginSupport_copy_data(
            GlobalWaypointType *dst,
            const GlobalWaypointType *src)
        {
            try {
                *dst = *src;
            } catch (...) {
                return RTI_FALSE;
            }

            return RTI_TRUE;
        }

        /* ----------------------------------------------------------------------------
        Callback functions:
        * ---------------------------------------------------------------------------- */

        PRESTypePluginParticipantData 
        GlobalWaypointTypePlugin_on_participant_attached(
            void *registration_data,
            const struct PRESTypePluginParticipantInfo *participant_info,
            RTIBool top_level_registration,
            void *container_plugin_context,
            RTICdrTypeCode *type_code)
        {
            try {
                struct RTIXCdrInterpreterPrograms *programs = NULL;
                struct PRESTypePluginDefaultParticipantData *pd = NULL;
                struct RTIXCdrInterpreterProgramsGenProperty programProperty =
                RTIXCdrInterpreterProgramsGenProperty_INITIALIZER;
                if (registration_data) {} /* To avoid warnings */
                if (participant_info) {} /* To avoid warnings */
                if (top_level_registration) {} /* To avoid warnings */
                if (container_plugin_context) {} /* To avoid warnings */
                if (type_code) {} /* To avoid warnings */
                pd = (struct PRESTypePluginDefaultParticipantData *)
                PRESTypePluginDefaultParticipantData_new(participant_info);

                programProperty.generateV1Encapsulation = RTI_XCDR_TRUE;
                programProperty.generateV2Encapsulation = RTI_XCDR_TRUE;
                programProperty.resolveAlias = RTI_XCDR_TRUE;
                programProperty.inlineStruct = RTI_XCDR_TRUE;
                programProperty.optimizeEnum = RTI_XCDR_TRUE;
                programProperty.unboundedSize = RTIXCdrLong_MAX;

                programProperty.externalReferenceSize = 
                (RTIXCdrUnsignedShort) sizeof(::dds::core::external<char>);
                programProperty.getExternalRefPointerFcn = 
                ::rti::topic::interpreter::get_external_value_pointer;

                programs = DDS_TypeCodeFactory_assert_programs_in_global_list(
                    DDS_TypeCodeFactory_get_instance(),
                    (DDS_TypeCode *) (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< GlobalWaypointType >::get().native()
                    ,
                    &programProperty,
                    RTI_XCDR_PROGRAM_MASK_TYPEPLUGIN);

                if (programs == NULL) {
                    PRESTypePluginDefaultParticipantData_delete(
                        (PRESTypePluginParticipantData)pd);
                    return NULL;
                }

                pd->programs = programs;
                return (PRESTypePluginParticipantData)pd;
            } catch (...) {
                return NULL;
            }
        }

        void 
        GlobalWaypointTypePlugin_on_participant_detached(
            PRESTypePluginParticipantData participant_data)
        {
            if (participant_data != NULL) {
                struct PRESTypePluginDefaultParticipantData *pd = 
                (struct PRESTypePluginDefaultParticipantData *)participant_data;

                if (pd->programs != NULL) {
                    DDS_TypeCodeFactory_remove_programs_from_global_list(
                        DDS_TypeCodeFactory_get_instance(),
                        pd->programs);
                    pd->programs = NULL;
                }
                PRESTypePluginDefaultParticipantData_delete(participant_data);
            }
        }

        PRESTypePluginEndpointData
        GlobalWaypointTypePlugin_on_endpoint_attached(
            PRESTypePluginParticipantData participant_data,
            const struct PRESTypePluginEndpointInfo *endpoint_info,
            RTIBool top_level_registration, 
            void *containerPluginContext)
        {
            try {
                PRESTypePluginEndpointData epd = NULL;
                unsigned int serializedSampleMaxSize = 0;

                if (top_level_registration) {} /* To avoid warnings */
                if (containerPluginContext) {} /* To avoid warnings */

                if (participant_data == NULL) {
                    return NULL;
                } 

                epd = PRESTypePluginDefaultEndpointData_new(
                    participant_data,
                    endpoint_info,
                    (PRESTypePluginDefaultEndpointDataCreateSampleFunction)
                    GlobalWaypointTypePluginSupport_create_data,
                    (PRESTypePluginDefaultEndpointDataDestroySampleFunction)
                    GlobalWaypointTypePluginSupport_destroy_data,
                    NULL , NULL );

                if (epd == NULL) {
                    return NULL;
                } 

                if (endpoint_info->endpointKind == PRES_TYPEPLUGIN_ENDPOINT_WRITER) {
                    serializedSampleMaxSize = ::ULARS::C2::GlobalWaypointTypePlugin_get_serialized_sample_max_size(
                        epd,RTI_FALSE,RTI_CDR_ENCAPSULATION_ID_CDR_BE,0);
                    PRESTypePluginDefaultEndpointData_setMaxSizeSerializedSample(epd, serializedSampleMaxSize);

                    if (PRESTypePluginDefaultEndpointData_createWriterPool(
                        epd,
                        endpoint_info,
                        (PRESTypePluginGetSerializedSampleMaxSizeFunction)
                        ::ULARS::C2::GlobalWaypointTypePlugin_get_serialized_sample_max_size, epd,
                        (PRESTypePluginGetSerializedSampleSizeFunction)
                        PRESTypePlugin_interpretedGetSerializedSampleSize,
                        epd) == RTI_FALSE) {
                        PRESTypePluginDefaultEndpointData_delete(epd);
                        return NULL;
                    }
                }

                return epd;
            } catch (...) {
                return NULL;
            }
        }

        void 
        GlobalWaypointTypePlugin_on_endpoint_detached(
            PRESTypePluginEndpointData endpoint_data)
        {  
            PRESTypePluginDefaultEndpointData_delete(endpoint_data);
        }

        void    
        GlobalWaypointTypePlugin_return_sample(
            PRESTypePluginEndpointData endpoint_data,
            GlobalWaypointType *sample,
            void *handle)
        {
            try {
                ::rti::topic::reset_sample(*sample);
            } catch(const std::exception& ex) {
                RTICdrLog_logWithFunctionName(
                    RTI_LOG_BIT_EXCEPTION,
                    "GlobalWaypointTypePlugin_return_sample",
                    &RTI_LOG_ANY_FAILURE_ss,
                    "exception: ",
                    ex.what());
            }

            PRESTypePluginDefaultEndpointData_returnSample(
                endpoint_data, sample, handle);
        }

        RTIBool 
        GlobalWaypointTypePlugin_copy_sample(
            PRESTypePluginEndpointData,
            GlobalWaypointType *dst,
            const GlobalWaypointType *src)
        {
            return ::ULARS::C2::GlobalWaypointTypePluginSupport_copy_data(dst,src);
        }

        /* ----------------------------------------------------------------------------
        (De)Serialize functions:
        * ------------------------------------------------------------------------- */
        unsigned int 
        GlobalWaypointTypePlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        RTIBool
        GlobalWaypointTypePlugin_serialize_to_cdr_buffer(
            char * buffer,
            unsigned int * length,
            const GlobalWaypointType *sample,
            ::dds::core::policy::DataRepresentationId representation)
        {
            using namespace ::dds::core::policy;

            try{
                RTIEncapsulationId encapsulationId = RTI_CDR_ENCAPSULATION_ID_INVALID;
                struct RTICdrStream cdrStream;
                struct PRESTypePluginDefaultEndpointData epd;
                RTIBool result;
                struct PRESTypePluginDefaultParticipantData pd;
                struct RTIXCdrTypePluginProgramContext defaultProgramContext =
                RTIXCdrTypePluginProgramContext_INTIALIZER;
                struct PRESTypePlugin plugin = PRES_TYPEPLUGIN_DEFAULT;

                if (length == NULL) {
                    return RTI_FALSE;
                }

                RTIOsapiMemory_zero(&epd, sizeof(struct PRESTypePluginDefaultEndpointData));
                epd.programContext = defaultProgramContext;
                epd._participantData = &pd;
                epd.typePlugin = &plugin;
                epd.programContext.endpointPluginData = &epd;
                try {
                    plugin.typeCode = (struct RTICdrTypeCode *)
                    (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< GlobalWaypointType >::get().native()
                    ;
                } catch (...) {
                    return RTI_FALSE;
                }
                pd.programs = ::rti::topic::interpreter::get_cdr_serialization_programs<
                GlobalWaypointType, 
                true, true, true>();

                encapsulationId = DDS_TypeCode_get_native_encapsulation(
                    (DDS_TypeCode *) plugin.typeCode,
                    representation);

                if (encapsulationId == RTI_CDR_ENCAPSULATION_ID_INVALID) {
                    return RTI_FALSE;
                }

                epd._maxSizeSerializedSample =
                GlobalWaypointTypePlugin_get_serialized_sample_max_size(
                    (PRESTypePluginEndpointData)&epd, 
                    RTI_TRUE, 
                    encapsulationId,
                    0);

                if (buffer == NULL) {
                    *length = 
                    PRESTypePlugin_interpretedGetSerializedSampleSize(
                        (PRESTypePluginEndpointData)&epd,
                        RTI_TRUE,
                        encapsulationId,
                        0,
                        sample);

                    if (*length == 0) {
                        return RTI_FALSE;
                    }

                    return RTI_TRUE;
                }    

                RTICdrStream_init(&cdrStream);
                RTICdrStream_set(&cdrStream, buffer, *length);

                result = PRESTypePlugin_interpretedSerialize(
                    (PRESTypePluginEndpointData)&epd, 
                    sample, 
                    &cdrStream, 
                    RTI_TRUE, 
                    encapsulationId,
                    RTI_TRUE, 
                    NULL);  

                *length = (unsigned int) RTICdrStream_getCurrentPositionOffset(&cdrStream);
                return result;
            } catch (...) {
                return RTI_FALSE;
            }
        }

        RTIBool
        GlobalWaypointTypePlugin_deserialize_from_cdr_buffer(
            GlobalWaypointType *sample,
            const char * buffer,
            unsigned int length)
        {
            struct RTICdrStream cdrStream;
            struct PRESTypePluginDefaultParticipantData pd;
            struct RTIXCdrTypePluginProgramContext defaultProgramContext =
            RTIXCdrTypePluginProgramContext_INTIALIZER;
            struct PRESTypePlugin plugin;
            struct PRESTypePluginDefaultEndpointData epd;

            RTICdrStream_init(&cdrStream);
            /*
            * The buffer is constant because this is a deserialization function
            * (the buffer is an input parameter, not an output parameter).
            * However, the buffer member in the stream is a (char *) so coverity
            * complains in case something else modifies the buffer's contents later.
            *
            * We don't currently have a stream type with a constant buffer.
            * Therefore, we suppress the warning after making sure that this function
            * doesn't modify the contents of the stream's buffer.
            */
            /* coverity[cert_exp40_c_violation : FALSE] */
            RTICdrStream_set(&cdrStream, (char *) buffer, length);

            epd.programContext = defaultProgramContext;
            epd._participantData = &pd;
            epd.typePlugin = &plugin;
            epd.programContext.endpointPluginData = &epd;
            try {
                plugin.typeCode = (struct RTICdrTypeCode *)
                (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< GlobalWaypointType >::get().native()
                ;
            } catch (...) {
                return RTI_FALSE;
            }
            pd.programs = ::rti::topic::interpreter::get_cdr_serialization_programs<
            GlobalWaypointType, 
            true, true, true>();

            epd._assignabilityProperty.acceptUnknownEnumValue = RTI_XCDR_TRUE;
            epd._assignabilityProperty.acceptUnknownUnionDiscriminator = 
            RTI_XCDR_ACCEPT_UNKNOWN_DISCRIMINATOR_AND_SELECT_DEFAULT;

            ::rti::topic::reset_sample(*sample);
            return PRESTypePlugin_interpretedDeserialize( 
                (PRESTypePluginEndpointData)&epd,
                sample,
                &cdrStream,
                RTI_TRUE, 
                RTI_TRUE, 
                NULL);
        }

        unsigned int 
        GlobalWaypointTypePlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            try {
                unsigned int size;
                RTIBool overflow = RTI_FALSE;

                size = PRESTypePlugin_interpretedGetSerializedSampleMaxSize(
                    endpoint_data,&overflow,include_encapsulation,encapsulation_id,current_alignment);

                if (overflow) {
                    size = RTI_CDR_MAX_SERIALIZED_SIZE;
                }

                return size;
            } catch (...) {
                return 0;
            }
        }

        /* --------------------------------------------------------------------------------------
        Key Management functions:
        * -------------------------------------------------------------------------------------- */

        PRESTypePluginKeyKind 
        GlobalWaypointTypePlugin_get_key_kind(void)
        {
            return PRES_TYPEPLUGIN_NO_KEY;
        }

        RTIBool GlobalWaypointTypePlugin_deserialize_key(
            PRESTypePluginEndpointData endpoint_data,
            GlobalWaypointType **sample, 
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_key,
            void *endpoint_plugin_qos)
        {
            try {
                RTIBool result;
                if (drop_sample) {} /* To avoid warnings */
                cdrStream->_xTypesState.unassignable = RTI_FALSE;
                result= PRESTypePlugin_interpretedDeserializeKey( 
                    endpoint_data, (sample != NULL)?*sample:NULL, cdrStream,
                    deserialize_encapsulation, deserialize_key, endpoint_plugin_qos);
                if (result) {
                    if (cdrStream->_xTypesState.unassignable) {
                        result = RTI_FALSE;
                    }
                }
                return result;    
            } catch (...) {
                return RTI_FALSE;
            }     
        }

        unsigned int
        GlobalWaypointTypePlugin_get_serialized_key_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            try {
                unsigned int size;
                RTIBool overflow = RTI_FALSE;

                size = PRESTypePlugin_interpretedGetSerializedKeyMaxSize(
                    endpoint_data,&overflow,include_encapsulation,encapsulation_id,current_alignment);
                if (overflow) {
                    size = RTI_CDR_MAX_SERIALIZED_SIZE;
                }

                return size;
            } catch (...) {
                return RTI_FALSE;
            }    
        }

        unsigned int
        GlobalWaypointTypePlugin_get_serialized_key_max_size_for_keyhash(
            PRESTypePluginEndpointData endpoint_data,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            unsigned int size;
            RTIBool overflow = RTI_FALSE;

            size = PRESTypePlugin_interpretedGetSerializedKeyMaxSizeForKeyhash(
                endpoint_data,
                &overflow,
                encapsulation_id,
                current_alignment);
            if (overflow) {
                size = RTI_CDR_MAX_SERIALIZED_SIZE;
            }

            return size;
        }

        /* ------------------------------------------------------------------------
        * Plug-in Installation Methods
        * ------------------------------------------------------------------------ */
        struct PRESTypePlugin *GlobalWaypointTypePlugin_new(void) 
        { 
            struct PRESTypePlugin *plugin = NULL;
            const struct PRESTypePluginVersion PLUGIN_VERSION = 
            PRES_TYPE_PLUGIN_VERSION_2_0;

            RTIOsapiHeap_allocateStructure(
                &plugin, struct PRESTypePlugin);
            if (plugin == NULL) {
                return NULL;
            }

            plugin->version = PLUGIN_VERSION;

            /* set up parent's function pointers */
            plugin->onParticipantAttached =
            (PRESTypePluginOnParticipantAttachedCallback)
            ::ULARS::C2::GlobalWaypointTypePlugin_on_participant_attached;
            plugin->onParticipantDetached =
            (PRESTypePluginOnParticipantDetachedCallback)
            ::ULARS::C2::GlobalWaypointTypePlugin_on_participant_detached;
            plugin->onEndpointAttached =
            (PRESTypePluginOnEndpointAttachedCallback)
            ::ULARS::C2::GlobalWaypointTypePlugin_on_endpoint_attached;
            plugin->onEndpointDetached =
            (PRESTypePluginOnEndpointDetachedCallback)
            ::ULARS::C2::GlobalWaypointTypePlugin_on_endpoint_detached;

            plugin->copySampleFnc =
            (PRESTypePluginCopySampleFunction)
            ::ULARS::C2::GlobalWaypointTypePlugin_copy_sample;
            plugin->createSampleFnc =
            (PRESTypePluginCreateSampleFunction)
            GlobalWaypointTypePlugin_create_sample;
            plugin->destroySampleFnc =
            (PRESTypePluginDestroySampleFunction)
            GlobalWaypointTypePlugin_destroy_sample;

            plugin->serializeFnc = 
            (PRESTypePluginSerializeFunction) PRESTypePlugin_interpretedSerialize;
            plugin->deserializeFnc =
            (PRESTypePluginDeserializeFunction) PRESTypePlugin_interpretedDeserializeWithAlloc;
            plugin->getSerializedSampleMaxSizeFnc =
            (PRESTypePluginGetSerializedSampleMaxSizeFunction)
            ::ULARS::C2::GlobalWaypointTypePlugin_get_serialized_sample_max_size;
            plugin->getSerializedSampleMinSizeFnc =
            (PRESTypePluginGetSerializedSampleMinSizeFunction)
            PRESTypePlugin_interpretedGetSerializedSampleMinSize;
            plugin->getDeserializedSampleMaxSizeFnc = NULL; 
            plugin->getSampleFnc =
            (PRESTypePluginGetSampleFunction)
            GlobalWaypointTypePlugin_get_sample;
            plugin->returnSampleFnc =
            (PRESTypePluginReturnSampleFunction)
            GlobalWaypointTypePlugin_return_sample;
            plugin->getKeyKindFnc =
            (PRESTypePluginGetKeyKindFunction)
            ::ULARS::C2::GlobalWaypointTypePlugin_get_key_kind;

            /* These functions are only used for keyed types. As this is not a keyed
            type they are all set to NULL
            */
            plugin->serializeKeyFnc = NULL ;    
            plugin->deserializeKeyFnc = NULL;  
            plugin->getKeyFnc = NULL;
            plugin->returnKeyFnc = NULL;
            plugin->instanceToKeyFnc = NULL;
            plugin->keyToInstanceFnc = NULL;
            plugin->getSerializedKeyMaxSizeFnc = NULL;
            plugin->instanceToKeyHashFnc = NULL;
            plugin->serializedSampleToKeyHashFnc = NULL;
            plugin->serializedKeyToKeyHashFnc = NULL;    
            #ifdef NDDS_STANDALONE_TYPE
            plugin->typeCode = NULL; 
            #else
            try {
                plugin->typeCode = (struct RTICdrTypeCode *)
                &::rti::topic::dynamic_type< ::ULARS::C2::GlobalWaypointType >::get().native();
            } catch (...) {
                ::ULARS::C2::GlobalWaypointTypePlugin_delete(plugin);
                return NULL;
            }
            #endif
            plugin->languageKind = PRES_TYPEPLUGIN_CPPSTL_LANG;

            /* Serialized buffer */
            plugin->getBuffer = 
            (PRESTypePluginGetBufferFunction)
            GlobalWaypointTypePlugin_get_buffer;
            plugin->returnBuffer = 
            (PRESTypePluginReturnBufferFunction)
            GlobalWaypointTypePlugin_return_buffer;
            plugin->getBufferWithParams = NULL;
            plugin->returnBufferWithParams = NULL;
            plugin->getSerializedSampleSizeFnc =
            (PRESTypePluginGetSerializedSampleSizeFunction)
            PRESTypePlugin_interpretedGetSerializedSampleSize;

            plugin->getWriterLoanedSampleFnc = NULL; 
            plugin->returnWriterLoanedSampleFnc = NULL;
            plugin->returnWriterLoanedSampleFromCookieFnc = NULL;
            plugin->validateWriterLoanedSampleFnc = NULL;
            plugin->setWriterLoanedSampleSerializedStateFnc = NULL;

            static const char * TYPE_NAME = "ULARS::C2::GlobalWaypointType";
            plugin->endpointTypeName = TYPE_NAME;
            plugin->isMetpType = RTI_FALSE;
            return plugin;
        }

        void
        GlobalWaypointTypePlugin_delete(struct PRESTypePlugin *plugin)
        {
            RTIOsapiHeap_freeStructure(plugin);
        } 

        /* ----------------------------------------------------------------------------
        *  Type StatusType
        * -------------------------------------------------------------------------- */

        /* -----------------------------------------------------------------------------
        Support functions:
        * -------------------------------------------------------------------------- */

        StatusType *
        StatusTypePluginSupport_create_data(void)
        {
            try {
                StatusType *sample = new StatusType();
                ::rti::topic::allocate_sample(*sample);
                return sample;
            } catch (...) {
                return NULL;
            }
        }

        void 
        StatusTypePluginSupport_destroy_data(
            StatusType *sample) 
        {
            delete sample;
        }

        RTIBool 
        StatusTypePluginSupport_copy_data(
            StatusType *dst,
            const StatusType *src)
        {
            try {
                *dst = *src;
            } catch (...) {
                return RTI_FALSE;
            }

            return RTI_TRUE;
        }

        /* ----------------------------------------------------------------------------
        Callback functions:
        * ---------------------------------------------------------------------------- */

        PRESTypePluginParticipantData 
        StatusTypePlugin_on_participant_attached(
            void *registration_data,
            const struct PRESTypePluginParticipantInfo *participant_info,
            RTIBool top_level_registration,
            void *container_plugin_context,
            RTICdrTypeCode *type_code)
        {
            try {
                struct RTIXCdrInterpreterPrograms *programs = NULL;
                struct PRESTypePluginDefaultParticipantData *pd = NULL;
                struct RTIXCdrInterpreterProgramsGenProperty programProperty =
                RTIXCdrInterpreterProgramsGenProperty_INITIALIZER;
                if (registration_data) {} /* To avoid warnings */
                if (participant_info) {} /* To avoid warnings */
                if (top_level_registration) {} /* To avoid warnings */
                if (container_plugin_context) {} /* To avoid warnings */
                if (type_code) {} /* To avoid warnings */
                pd = (struct PRESTypePluginDefaultParticipantData *)
                PRESTypePluginDefaultParticipantData_new(participant_info);

                programProperty.generateV1Encapsulation = RTI_XCDR_TRUE;
                programProperty.generateV2Encapsulation = RTI_XCDR_TRUE;
                programProperty.resolveAlias = RTI_XCDR_TRUE;
                programProperty.inlineStruct = RTI_XCDR_TRUE;
                programProperty.optimizeEnum = RTI_XCDR_TRUE;
                programProperty.unboundedSize = RTIXCdrLong_MAX;

                programProperty.externalReferenceSize = 
                (RTIXCdrUnsignedShort) sizeof(::dds::core::external<char>);
                programProperty.getExternalRefPointerFcn = 
                ::rti::topic::interpreter::get_external_value_pointer;

                programs = DDS_TypeCodeFactory_assert_programs_in_global_list(
                    DDS_TypeCodeFactory_get_instance(),
                    (DDS_TypeCode *) (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< StatusType >::get().native()
                    ,
                    &programProperty,
                    RTI_XCDR_PROGRAM_MASK_TYPEPLUGIN);

                if (programs == NULL) {
                    PRESTypePluginDefaultParticipantData_delete(
                        (PRESTypePluginParticipantData)pd);
                    return NULL;
                }

                pd->programs = programs;
                return (PRESTypePluginParticipantData)pd;
            } catch (...) {
                return NULL;
            }
        }

        void 
        StatusTypePlugin_on_participant_detached(
            PRESTypePluginParticipantData participant_data)
        {
            if (participant_data != NULL) {
                struct PRESTypePluginDefaultParticipantData *pd = 
                (struct PRESTypePluginDefaultParticipantData *)participant_data;

                if (pd->programs != NULL) {
                    DDS_TypeCodeFactory_remove_programs_from_global_list(
                        DDS_TypeCodeFactory_get_instance(),
                        pd->programs);
                    pd->programs = NULL;
                }
                PRESTypePluginDefaultParticipantData_delete(participant_data);
            }
        }

        PRESTypePluginEndpointData
        StatusTypePlugin_on_endpoint_attached(
            PRESTypePluginParticipantData participant_data,
            const struct PRESTypePluginEndpointInfo *endpoint_info,
            RTIBool top_level_registration, 
            void *containerPluginContext)
        {
            try {
                PRESTypePluginEndpointData epd = NULL;
                unsigned int serializedSampleMaxSize = 0;

                if (top_level_registration) {} /* To avoid warnings */
                if (containerPluginContext) {} /* To avoid warnings */

                if (participant_data == NULL) {
                    return NULL;
                } 

                epd = PRESTypePluginDefaultEndpointData_new(
                    participant_data,
                    endpoint_info,
                    (PRESTypePluginDefaultEndpointDataCreateSampleFunction)
                    StatusTypePluginSupport_create_data,
                    (PRESTypePluginDefaultEndpointDataDestroySampleFunction)
                    StatusTypePluginSupport_destroy_data,
                    NULL , NULL );

                if (epd == NULL) {
                    return NULL;
                } 

                if (endpoint_info->endpointKind == PRES_TYPEPLUGIN_ENDPOINT_WRITER) {
                    serializedSampleMaxSize = ::ULARS::C2::StatusTypePlugin_get_serialized_sample_max_size(
                        epd,RTI_FALSE,RTI_CDR_ENCAPSULATION_ID_CDR_BE,0);
                    PRESTypePluginDefaultEndpointData_setMaxSizeSerializedSample(epd, serializedSampleMaxSize);

                    if (PRESTypePluginDefaultEndpointData_createWriterPool(
                        epd,
                        endpoint_info,
                        (PRESTypePluginGetSerializedSampleMaxSizeFunction)
                        ::ULARS::C2::StatusTypePlugin_get_serialized_sample_max_size, epd,
                        (PRESTypePluginGetSerializedSampleSizeFunction)
                        PRESTypePlugin_interpretedGetSerializedSampleSize,
                        epd) == RTI_FALSE) {
                        PRESTypePluginDefaultEndpointData_delete(epd);
                        return NULL;
                    }
                }

                return epd;
            } catch (...) {
                return NULL;
            }
        }

        void 
        StatusTypePlugin_on_endpoint_detached(
            PRESTypePluginEndpointData endpoint_data)
        {  
            PRESTypePluginDefaultEndpointData_delete(endpoint_data);
        }

        void    
        StatusTypePlugin_return_sample(
            PRESTypePluginEndpointData endpoint_data,
            StatusType *sample,
            void *handle)
        {
            try {
                ::rti::topic::reset_sample(*sample);
            } catch(const std::exception& ex) {
                RTICdrLog_logWithFunctionName(
                    RTI_LOG_BIT_EXCEPTION,
                    "StatusTypePlugin_return_sample",
                    &RTI_LOG_ANY_FAILURE_ss,
                    "exception: ",
                    ex.what());
            }

            PRESTypePluginDefaultEndpointData_returnSample(
                endpoint_data, sample, handle);
        }

        RTIBool 
        StatusTypePlugin_copy_sample(
            PRESTypePluginEndpointData,
            StatusType *dst,
            const StatusType *src)
        {
            return ::ULARS::C2::StatusTypePluginSupport_copy_data(dst,src);
        }

        /* ----------------------------------------------------------------------------
        (De)Serialize functions:
        * ------------------------------------------------------------------------- */
        unsigned int 
        StatusTypePlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        RTIBool
        StatusTypePlugin_serialize_to_cdr_buffer(
            char * buffer,
            unsigned int * length,
            const StatusType *sample,
            ::dds::core::policy::DataRepresentationId representation)
        {
            using namespace ::dds::core::policy;

            try{
                RTIEncapsulationId encapsulationId = RTI_CDR_ENCAPSULATION_ID_INVALID;
                struct RTICdrStream cdrStream;
                struct PRESTypePluginDefaultEndpointData epd;
                RTIBool result;
                struct PRESTypePluginDefaultParticipantData pd;
                struct RTIXCdrTypePluginProgramContext defaultProgramContext =
                RTIXCdrTypePluginProgramContext_INTIALIZER;
                struct PRESTypePlugin plugin = PRES_TYPEPLUGIN_DEFAULT;

                if (length == NULL) {
                    return RTI_FALSE;
                }

                RTIOsapiMemory_zero(&epd, sizeof(struct PRESTypePluginDefaultEndpointData));
                epd.programContext = defaultProgramContext;
                epd._participantData = &pd;
                epd.typePlugin = &plugin;
                epd.programContext.endpointPluginData = &epd;
                try {
                    plugin.typeCode = (struct RTICdrTypeCode *)
                    (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< StatusType >::get().native()
                    ;
                } catch (...) {
                    return RTI_FALSE;
                }
                pd.programs = ::rti::topic::interpreter::get_cdr_serialization_programs<
                StatusType, 
                true, true, true>();

                encapsulationId = DDS_TypeCode_get_native_encapsulation(
                    (DDS_TypeCode *) plugin.typeCode,
                    representation);

                if (encapsulationId == RTI_CDR_ENCAPSULATION_ID_INVALID) {
                    return RTI_FALSE;
                }

                epd._maxSizeSerializedSample =
                StatusTypePlugin_get_serialized_sample_max_size(
                    (PRESTypePluginEndpointData)&epd, 
                    RTI_TRUE, 
                    encapsulationId,
                    0);

                if (buffer == NULL) {
                    *length = 
                    PRESTypePlugin_interpretedGetSerializedSampleSize(
                        (PRESTypePluginEndpointData)&epd,
                        RTI_TRUE,
                        encapsulationId,
                        0,
                        sample);

                    if (*length == 0) {
                        return RTI_FALSE;
                    }

                    return RTI_TRUE;
                }    

                RTICdrStream_init(&cdrStream);
                RTICdrStream_set(&cdrStream, buffer, *length);

                result = PRESTypePlugin_interpretedSerialize(
                    (PRESTypePluginEndpointData)&epd, 
                    sample, 
                    &cdrStream, 
                    RTI_TRUE, 
                    encapsulationId,
                    RTI_TRUE, 
                    NULL);  

                *length = (unsigned int) RTICdrStream_getCurrentPositionOffset(&cdrStream);
                return result;
            } catch (...) {
                return RTI_FALSE;
            }
        }

        RTIBool
        StatusTypePlugin_deserialize_from_cdr_buffer(
            StatusType *sample,
            const char * buffer,
            unsigned int length)
        {
            struct RTICdrStream cdrStream;
            struct PRESTypePluginDefaultParticipantData pd;
            struct RTIXCdrTypePluginProgramContext defaultProgramContext =
            RTIXCdrTypePluginProgramContext_INTIALIZER;
            struct PRESTypePlugin plugin;
            struct PRESTypePluginDefaultEndpointData epd;

            RTICdrStream_init(&cdrStream);
            /*
            * The buffer is constant because this is a deserialization function
            * (the buffer is an input parameter, not an output parameter).
            * However, the buffer member in the stream is a (char *) so coverity
            * complains in case something else modifies the buffer's contents later.
            *
            * We don't currently have a stream type with a constant buffer.
            * Therefore, we suppress the warning after making sure that this function
            * doesn't modify the contents of the stream's buffer.
            */
            /* coverity[cert_exp40_c_violation : FALSE] */
            RTICdrStream_set(&cdrStream, (char *) buffer, length);

            epd.programContext = defaultProgramContext;
            epd._participantData = &pd;
            epd.typePlugin = &plugin;
            epd.programContext.endpointPluginData = &epd;
            try {
                plugin.typeCode = (struct RTICdrTypeCode *)
                (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< StatusType >::get().native()
                ;
            } catch (...) {
                return RTI_FALSE;
            }
            pd.programs = ::rti::topic::interpreter::get_cdr_serialization_programs<
            StatusType, 
            true, true, true>();

            epd._assignabilityProperty.acceptUnknownEnumValue = RTI_XCDR_TRUE;
            epd._assignabilityProperty.acceptUnknownUnionDiscriminator = 
            RTI_XCDR_ACCEPT_UNKNOWN_DISCRIMINATOR_AND_SELECT_DEFAULT;

            ::rti::topic::reset_sample(*sample);
            return PRESTypePlugin_interpretedDeserialize( 
                (PRESTypePluginEndpointData)&epd,
                sample,
                &cdrStream,
                RTI_TRUE, 
                RTI_TRUE, 
                NULL);
        }

        unsigned int 
        StatusTypePlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            try {
                unsigned int size;
                RTIBool overflow = RTI_FALSE;

                size = PRESTypePlugin_interpretedGetSerializedSampleMaxSize(
                    endpoint_data,&overflow,include_encapsulation,encapsulation_id,current_alignment);

                if (overflow) {
                    size = RTI_CDR_MAX_SERIALIZED_SIZE;
                }

                return size;
            } catch (...) {
                return 0;
            }
        }

        /* --------------------------------------------------------------------------------------
        Key Management functions:
        * -------------------------------------------------------------------------------------- */

        PRESTypePluginKeyKind 
        StatusTypePlugin_get_key_kind(void)
        {
            return PRES_TYPEPLUGIN_NO_KEY;
        }

        RTIBool StatusTypePlugin_deserialize_key(
            PRESTypePluginEndpointData endpoint_data,
            StatusType **sample, 
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_key,
            void *endpoint_plugin_qos)
        {
            try {
                RTIBool result;
                if (drop_sample) {} /* To avoid warnings */
                cdrStream->_xTypesState.unassignable = RTI_FALSE;
                result= PRESTypePlugin_interpretedDeserializeKey( 
                    endpoint_data, (sample != NULL)?*sample:NULL, cdrStream,
                    deserialize_encapsulation, deserialize_key, endpoint_plugin_qos);
                if (result) {
                    if (cdrStream->_xTypesState.unassignable) {
                        result = RTI_FALSE;
                    }
                }
                return result;    
            } catch (...) {
                return RTI_FALSE;
            }     
        }

        unsigned int
        StatusTypePlugin_get_serialized_key_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            try {
                unsigned int size;
                RTIBool overflow = RTI_FALSE;

                size = PRESTypePlugin_interpretedGetSerializedKeyMaxSize(
                    endpoint_data,&overflow,include_encapsulation,encapsulation_id,current_alignment);
                if (overflow) {
                    size = RTI_CDR_MAX_SERIALIZED_SIZE;
                }

                return size;
            } catch (...) {
                return RTI_FALSE;
            }    
        }

        unsigned int
        StatusTypePlugin_get_serialized_key_max_size_for_keyhash(
            PRESTypePluginEndpointData endpoint_data,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            unsigned int size;
            RTIBool overflow = RTI_FALSE;

            size = PRESTypePlugin_interpretedGetSerializedKeyMaxSizeForKeyhash(
                endpoint_data,
                &overflow,
                encapsulation_id,
                current_alignment);
            if (overflow) {
                size = RTI_CDR_MAX_SERIALIZED_SIZE;
            }

            return size;
        }

        /* ------------------------------------------------------------------------
        * Plug-in Installation Methods
        * ------------------------------------------------------------------------ */
        struct PRESTypePlugin *StatusTypePlugin_new(void) 
        { 
            struct PRESTypePlugin *plugin = NULL;
            const struct PRESTypePluginVersion PLUGIN_VERSION = 
            PRES_TYPE_PLUGIN_VERSION_2_0;

            RTIOsapiHeap_allocateStructure(
                &plugin, struct PRESTypePlugin);
            if (plugin == NULL) {
                return NULL;
            }

            plugin->version = PLUGIN_VERSION;

            /* set up parent's function pointers */
            plugin->onParticipantAttached =
            (PRESTypePluginOnParticipantAttachedCallback)
            ::ULARS::C2::StatusTypePlugin_on_participant_attached;
            plugin->onParticipantDetached =
            (PRESTypePluginOnParticipantDetachedCallback)
            ::ULARS::C2::StatusTypePlugin_on_participant_detached;
            plugin->onEndpointAttached =
            (PRESTypePluginOnEndpointAttachedCallback)
            ::ULARS::C2::StatusTypePlugin_on_endpoint_attached;
            plugin->onEndpointDetached =
            (PRESTypePluginOnEndpointDetachedCallback)
            ::ULARS::C2::StatusTypePlugin_on_endpoint_detached;

            plugin->copySampleFnc =
            (PRESTypePluginCopySampleFunction)
            ::ULARS::C2::StatusTypePlugin_copy_sample;
            plugin->createSampleFnc =
            (PRESTypePluginCreateSampleFunction)
            StatusTypePlugin_create_sample;
            plugin->destroySampleFnc =
            (PRESTypePluginDestroySampleFunction)
            StatusTypePlugin_destroy_sample;

            plugin->serializeFnc = 
            (PRESTypePluginSerializeFunction) PRESTypePlugin_interpretedSerialize;
            plugin->deserializeFnc =
            (PRESTypePluginDeserializeFunction) PRESTypePlugin_interpretedDeserializeWithAlloc;
            plugin->getSerializedSampleMaxSizeFnc =
            (PRESTypePluginGetSerializedSampleMaxSizeFunction)
            ::ULARS::C2::StatusTypePlugin_get_serialized_sample_max_size;
            plugin->getSerializedSampleMinSizeFnc =
            (PRESTypePluginGetSerializedSampleMinSizeFunction)
            PRESTypePlugin_interpretedGetSerializedSampleMinSize;
            plugin->getDeserializedSampleMaxSizeFnc = NULL; 
            plugin->getSampleFnc =
            (PRESTypePluginGetSampleFunction)
            StatusTypePlugin_get_sample;
            plugin->returnSampleFnc =
            (PRESTypePluginReturnSampleFunction)
            StatusTypePlugin_return_sample;
            plugin->getKeyKindFnc =
            (PRESTypePluginGetKeyKindFunction)
            ::ULARS::C2::StatusTypePlugin_get_key_kind;

            /* These functions are only used for keyed types. As this is not a keyed
            type they are all set to NULL
            */
            plugin->serializeKeyFnc = NULL ;    
            plugin->deserializeKeyFnc = NULL;  
            plugin->getKeyFnc = NULL;
            plugin->returnKeyFnc = NULL;
            plugin->instanceToKeyFnc = NULL;
            plugin->keyToInstanceFnc = NULL;
            plugin->getSerializedKeyMaxSizeFnc = NULL;
            plugin->instanceToKeyHashFnc = NULL;
            plugin->serializedSampleToKeyHashFnc = NULL;
            plugin->serializedKeyToKeyHashFnc = NULL;    
            #ifdef NDDS_STANDALONE_TYPE
            plugin->typeCode = NULL; 
            #else
            try {
                plugin->typeCode = (struct RTICdrTypeCode *)
                &::rti::topic::dynamic_type< ::ULARS::C2::StatusType >::get().native();
            } catch (...) {
                ::ULARS::C2::StatusTypePlugin_delete(plugin);
                return NULL;
            }
            #endif
            plugin->languageKind = PRES_TYPEPLUGIN_CPPSTL_LANG;

            /* Serialized buffer */
            plugin->getBuffer = 
            (PRESTypePluginGetBufferFunction)
            StatusTypePlugin_get_buffer;
            plugin->returnBuffer = 
            (PRESTypePluginReturnBufferFunction)
            StatusTypePlugin_return_buffer;
            plugin->getBufferWithParams = NULL;
            plugin->returnBufferWithParams = NULL;
            plugin->getSerializedSampleSizeFnc =
            (PRESTypePluginGetSerializedSampleSizeFunction)
            PRESTypePlugin_interpretedGetSerializedSampleSize;

            plugin->getWriterLoanedSampleFnc = NULL; 
            plugin->returnWriterLoanedSampleFnc = NULL;
            plugin->returnWriterLoanedSampleFromCookieFnc = NULL;
            plugin->validateWriterLoanedSampleFnc = NULL;
            plugin->setWriterLoanedSampleSerializedStateFnc = NULL;

            static const char * TYPE_NAME = "ULARS::C2::StatusType";
            plugin->endpointTypeName = TYPE_NAME;
            plugin->isMetpType = RTI_FALSE;
            return plugin;
        }

        void
        StatusTypePlugin_delete(struct PRESTypePlugin *plugin)
        {
            RTIOsapiHeap_freeStructure(plugin);
        } 

        /* ----------------------------------------------------------------------------
        *  Type ObstacleType
        * -------------------------------------------------------------------------- */

        /* -----------------------------------------------------------------------------
        Support functions:
        * -------------------------------------------------------------------------- */

        ObstacleType *
        ObstacleTypePluginSupport_create_data(void)
        {
            try {
                ObstacleType *sample = new ObstacleType();
                ::rti::topic::allocate_sample(*sample);
                return sample;
            } catch (...) {
                return NULL;
            }
        }

        void 
        ObstacleTypePluginSupport_destroy_data(
            ObstacleType *sample) 
        {
            delete sample;
        }

        RTIBool 
        ObstacleTypePluginSupport_copy_data(
            ObstacleType *dst,
            const ObstacleType *src)
        {
            try {
                *dst = *src;
            } catch (...) {
                return RTI_FALSE;
            }

            return RTI_TRUE;
        }

        /* ----------------------------------------------------------------------------
        Callback functions:
        * ---------------------------------------------------------------------------- */

        PRESTypePluginParticipantData 
        ObstacleTypePlugin_on_participant_attached(
            void *registration_data,
            const struct PRESTypePluginParticipantInfo *participant_info,
            RTIBool top_level_registration,
            void *container_plugin_context,
            RTICdrTypeCode *type_code)
        {
            try {
                struct RTIXCdrInterpreterPrograms *programs = NULL;
                struct PRESTypePluginDefaultParticipantData *pd = NULL;
                struct RTIXCdrInterpreterProgramsGenProperty programProperty =
                RTIXCdrInterpreterProgramsGenProperty_INITIALIZER;
                if (registration_data) {} /* To avoid warnings */
                if (participant_info) {} /* To avoid warnings */
                if (top_level_registration) {} /* To avoid warnings */
                if (container_plugin_context) {} /* To avoid warnings */
                if (type_code) {} /* To avoid warnings */
                pd = (struct PRESTypePluginDefaultParticipantData *)
                PRESTypePluginDefaultParticipantData_new(participant_info);

                programProperty.generateV1Encapsulation = RTI_XCDR_TRUE;
                programProperty.generateV2Encapsulation = RTI_XCDR_TRUE;
                programProperty.resolveAlias = RTI_XCDR_TRUE;
                programProperty.inlineStruct = RTI_XCDR_TRUE;
                programProperty.optimizeEnum = RTI_XCDR_TRUE;
                programProperty.unboundedSize = RTIXCdrLong_MAX;

                programProperty.externalReferenceSize = 
                (RTIXCdrUnsignedShort) sizeof(::dds::core::external<char>);
                programProperty.getExternalRefPointerFcn = 
                ::rti::topic::interpreter::get_external_value_pointer;

                programs = DDS_TypeCodeFactory_assert_programs_in_global_list(
                    DDS_TypeCodeFactory_get_instance(),
                    (DDS_TypeCode *) (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< ObstacleType >::get().native()
                    ,
                    &programProperty,
                    RTI_XCDR_PROGRAM_MASK_TYPEPLUGIN);

                if (programs == NULL) {
                    PRESTypePluginDefaultParticipantData_delete(
                        (PRESTypePluginParticipantData)pd);
                    return NULL;
                }

                pd->programs = programs;
                return (PRESTypePluginParticipantData)pd;
            } catch (...) {
                return NULL;
            }
        }

        void 
        ObstacleTypePlugin_on_participant_detached(
            PRESTypePluginParticipantData participant_data)
        {
            if (participant_data != NULL) {
                struct PRESTypePluginDefaultParticipantData *pd = 
                (struct PRESTypePluginDefaultParticipantData *)participant_data;

                if (pd->programs != NULL) {
                    DDS_TypeCodeFactory_remove_programs_from_global_list(
                        DDS_TypeCodeFactory_get_instance(),
                        pd->programs);
                    pd->programs = NULL;
                }
                PRESTypePluginDefaultParticipantData_delete(participant_data);
            }
        }

        PRESTypePluginEndpointData
        ObstacleTypePlugin_on_endpoint_attached(
            PRESTypePluginParticipantData participant_data,
            const struct PRESTypePluginEndpointInfo *endpoint_info,
            RTIBool top_level_registration, 
            void *containerPluginContext)
        {
            try {
                PRESTypePluginEndpointData epd = NULL;
                unsigned int serializedSampleMaxSize = 0;

                if (top_level_registration) {} /* To avoid warnings */
                if (containerPluginContext) {} /* To avoid warnings */

                if (participant_data == NULL) {
                    return NULL;
                } 

                epd = PRESTypePluginDefaultEndpointData_new(
                    participant_data,
                    endpoint_info,
                    (PRESTypePluginDefaultEndpointDataCreateSampleFunction)
                    ObstacleTypePluginSupport_create_data,
                    (PRESTypePluginDefaultEndpointDataDestroySampleFunction)
                    ObstacleTypePluginSupport_destroy_data,
                    NULL , NULL );

                if (epd == NULL) {
                    return NULL;
                } 

                if (endpoint_info->endpointKind == PRES_TYPEPLUGIN_ENDPOINT_WRITER) {
                    serializedSampleMaxSize = ::ULARS::C2::ObstacleTypePlugin_get_serialized_sample_max_size(
                        epd,RTI_FALSE,RTI_CDR_ENCAPSULATION_ID_CDR_BE,0);
                    PRESTypePluginDefaultEndpointData_setMaxSizeSerializedSample(epd, serializedSampleMaxSize);

                    if (PRESTypePluginDefaultEndpointData_createWriterPool(
                        epd,
                        endpoint_info,
                        (PRESTypePluginGetSerializedSampleMaxSizeFunction)
                        ::ULARS::C2::ObstacleTypePlugin_get_serialized_sample_max_size, epd,
                        (PRESTypePluginGetSerializedSampleSizeFunction)
                        PRESTypePlugin_interpretedGetSerializedSampleSize,
                        epd) == RTI_FALSE) {
                        PRESTypePluginDefaultEndpointData_delete(epd);
                        return NULL;
                    }
                }

                return epd;
            } catch (...) {
                return NULL;
            }
        }

        void 
        ObstacleTypePlugin_on_endpoint_detached(
            PRESTypePluginEndpointData endpoint_data)
        {  
            PRESTypePluginDefaultEndpointData_delete(endpoint_data);
        }

        void    
        ObstacleTypePlugin_return_sample(
            PRESTypePluginEndpointData endpoint_data,
            ObstacleType *sample,
            void *handle)
        {
            try {
                ::rti::topic::reset_sample(*sample);
            } catch(const std::exception& ex) {
                RTICdrLog_logWithFunctionName(
                    RTI_LOG_BIT_EXCEPTION,
                    "ObstacleTypePlugin_return_sample",
                    &RTI_LOG_ANY_FAILURE_ss,
                    "exception: ",
                    ex.what());
            }

            PRESTypePluginDefaultEndpointData_returnSample(
                endpoint_data, sample, handle);
        }

        RTIBool 
        ObstacleTypePlugin_copy_sample(
            PRESTypePluginEndpointData,
            ObstacleType *dst,
            const ObstacleType *src)
        {
            return ::ULARS::C2::ObstacleTypePluginSupport_copy_data(dst,src);
        }

        /* ----------------------------------------------------------------------------
        (De)Serialize functions:
        * ------------------------------------------------------------------------- */
        unsigned int 
        ObstacleTypePlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        RTIBool
        ObstacleTypePlugin_serialize_to_cdr_buffer(
            char * buffer,
            unsigned int * length,
            const ObstacleType *sample,
            ::dds::core::policy::DataRepresentationId representation)
        {
            using namespace ::dds::core::policy;

            try{
                RTIEncapsulationId encapsulationId = RTI_CDR_ENCAPSULATION_ID_INVALID;
                struct RTICdrStream cdrStream;
                struct PRESTypePluginDefaultEndpointData epd;
                RTIBool result;
                struct PRESTypePluginDefaultParticipantData pd;
                struct RTIXCdrTypePluginProgramContext defaultProgramContext =
                RTIXCdrTypePluginProgramContext_INTIALIZER;
                struct PRESTypePlugin plugin = PRES_TYPEPLUGIN_DEFAULT;

                if (length == NULL) {
                    return RTI_FALSE;
                }

                RTIOsapiMemory_zero(&epd, sizeof(struct PRESTypePluginDefaultEndpointData));
                epd.programContext = defaultProgramContext;
                epd._participantData = &pd;
                epd.typePlugin = &plugin;
                epd.programContext.endpointPluginData = &epd;
                try {
                    plugin.typeCode = (struct RTICdrTypeCode *)
                    (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< ObstacleType >::get().native()
                    ;
                } catch (...) {
                    return RTI_FALSE;
                }
                pd.programs = ::rti::topic::interpreter::get_cdr_serialization_programs<
                ObstacleType, 
                true, true, true>();

                encapsulationId = DDS_TypeCode_get_native_encapsulation(
                    (DDS_TypeCode *) plugin.typeCode,
                    representation);

                if (encapsulationId == RTI_CDR_ENCAPSULATION_ID_INVALID) {
                    return RTI_FALSE;
                }

                epd._maxSizeSerializedSample =
                ObstacleTypePlugin_get_serialized_sample_max_size(
                    (PRESTypePluginEndpointData)&epd, 
                    RTI_TRUE, 
                    encapsulationId,
                    0);

                if (buffer == NULL) {
                    *length = 
                    PRESTypePlugin_interpretedGetSerializedSampleSize(
                        (PRESTypePluginEndpointData)&epd,
                        RTI_TRUE,
                        encapsulationId,
                        0,
                        sample);

                    if (*length == 0) {
                        return RTI_FALSE;
                    }

                    return RTI_TRUE;
                }    

                RTICdrStream_init(&cdrStream);
                RTICdrStream_set(&cdrStream, buffer, *length);

                result = PRESTypePlugin_interpretedSerialize(
                    (PRESTypePluginEndpointData)&epd, 
                    sample, 
                    &cdrStream, 
                    RTI_TRUE, 
                    encapsulationId,
                    RTI_TRUE, 
                    NULL);  

                *length = (unsigned int) RTICdrStream_getCurrentPositionOffset(&cdrStream);
                return result;
            } catch (...) {
                return RTI_FALSE;
            }
        }

        RTIBool
        ObstacleTypePlugin_deserialize_from_cdr_buffer(
            ObstacleType *sample,
            const char * buffer,
            unsigned int length)
        {
            struct RTICdrStream cdrStream;
            struct PRESTypePluginDefaultParticipantData pd;
            struct RTIXCdrTypePluginProgramContext defaultProgramContext =
            RTIXCdrTypePluginProgramContext_INTIALIZER;
            struct PRESTypePlugin plugin;
            struct PRESTypePluginDefaultEndpointData epd;

            RTICdrStream_init(&cdrStream);
            /*
            * The buffer is constant because this is a deserialization function
            * (the buffer is an input parameter, not an output parameter).
            * However, the buffer member in the stream is a (char *) so coverity
            * complains in case something else modifies the buffer's contents later.
            *
            * We don't currently have a stream type with a constant buffer.
            * Therefore, we suppress the warning after making sure that this function
            * doesn't modify the contents of the stream's buffer.
            */
            /* coverity[cert_exp40_c_violation : FALSE] */
            RTICdrStream_set(&cdrStream, (char *) buffer, length);

            epd.programContext = defaultProgramContext;
            epd._participantData = &pd;
            epd.typePlugin = &plugin;
            epd.programContext.endpointPluginData = &epd;
            try {
                plugin.typeCode = (struct RTICdrTypeCode *)
                (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< ObstacleType >::get().native()
                ;
            } catch (...) {
                return RTI_FALSE;
            }
            pd.programs = ::rti::topic::interpreter::get_cdr_serialization_programs<
            ObstacleType, 
            true, true, true>();

            epd._assignabilityProperty.acceptUnknownEnumValue = RTI_XCDR_TRUE;
            epd._assignabilityProperty.acceptUnknownUnionDiscriminator = 
            RTI_XCDR_ACCEPT_UNKNOWN_DISCRIMINATOR_AND_SELECT_DEFAULT;

            ::rti::topic::reset_sample(*sample);
            return PRESTypePlugin_interpretedDeserialize( 
                (PRESTypePluginEndpointData)&epd,
                sample,
                &cdrStream,
                RTI_TRUE, 
                RTI_TRUE, 
                NULL);
        }

        unsigned int 
        ObstacleTypePlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            try {
                unsigned int size;
                RTIBool overflow = RTI_FALSE;

                size = PRESTypePlugin_interpretedGetSerializedSampleMaxSize(
                    endpoint_data,&overflow,include_encapsulation,encapsulation_id,current_alignment);

                if (overflow) {
                    size = RTI_CDR_MAX_SERIALIZED_SIZE;
                }

                return size;
            } catch (...) {
                return 0;
            }
        }

        /* --------------------------------------------------------------------------------------
        Key Management functions:
        * -------------------------------------------------------------------------------------- */

        PRESTypePluginKeyKind 
        ObstacleTypePlugin_get_key_kind(void)
        {
            return PRES_TYPEPLUGIN_NO_KEY;
        }

        RTIBool ObstacleTypePlugin_deserialize_key(
            PRESTypePluginEndpointData endpoint_data,
            ObstacleType **sample, 
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_key,
            void *endpoint_plugin_qos)
        {
            try {
                RTIBool result;
                if (drop_sample) {} /* To avoid warnings */
                cdrStream->_xTypesState.unassignable = RTI_FALSE;
                result= PRESTypePlugin_interpretedDeserializeKey( 
                    endpoint_data, (sample != NULL)?*sample:NULL, cdrStream,
                    deserialize_encapsulation, deserialize_key, endpoint_plugin_qos);
                if (result) {
                    if (cdrStream->_xTypesState.unassignable) {
                        result = RTI_FALSE;
                    }
                }
                return result;    
            } catch (...) {
                return RTI_FALSE;
            }     
        }

        unsigned int
        ObstacleTypePlugin_get_serialized_key_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            try {
                unsigned int size;
                RTIBool overflow = RTI_FALSE;

                size = PRESTypePlugin_interpretedGetSerializedKeyMaxSize(
                    endpoint_data,&overflow,include_encapsulation,encapsulation_id,current_alignment);
                if (overflow) {
                    size = RTI_CDR_MAX_SERIALIZED_SIZE;
                }

                return size;
            } catch (...) {
                return RTI_FALSE;
            }    
        }

        unsigned int
        ObstacleTypePlugin_get_serialized_key_max_size_for_keyhash(
            PRESTypePluginEndpointData endpoint_data,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            unsigned int size;
            RTIBool overflow = RTI_FALSE;

            size = PRESTypePlugin_interpretedGetSerializedKeyMaxSizeForKeyhash(
                endpoint_data,
                &overflow,
                encapsulation_id,
                current_alignment);
            if (overflow) {
                size = RTI_CDR_MAX_SERIALIZED_SIZE;
            }

            return size;
        }

        /* ------------------------------------------------------------------------
        * Plug-in Installation Methods
        * ------------------------------------------------------------------------ */
        struct PRESTypePlugin *ObstacleTypePlugin_new(void) 
        { 
            struct PRESTypePlugin *plugin = NULL;
            const struct PRESTypePluginVersion PLUGIN_VERSION = 
            PRES_TYPE_PLUGIN_VERSION_2_0;

            RTIOsapiHeap_allocateStructure(
                &plugin, struct PRESTypePlugin);
            if (plugin == NULL) {
                return NULL;
            }

            plugin->version = PLUGIN_VERSION;

            /* set up parent's function pointers */
            plugin->onParticipantAttached =
            (PRESTypePluginOnParticipantAttachedCallback)
            ::ULARS::C2::ObstacleTypePlugin_on_participant_attached;
            plugin->onParticipantDetached =
            (PRESTypePluginOnParticipantDetachedCallback)
            ::ULARS::C2::ObstacleTypePlugin_on_participant_detached;
            plugin->onEndpointAttached =
            (PRESTypePluginOnEndpointAttachedCallback)
            ::ULARS::C2::ObstacleTypePlugin_on_endpoint_attached;
            plugin->onEndpointDetached =
            (PRESTypePluginOnEndpointDetachedCallback)
            ::ULARS::C2::ObstacleTypePlugin_on_endpoint_detached;

            plugin->copySampleFnc =
            (PRESTypePluginCopySampleFunction)
            ::ULARS::C2::ObstacleTypePlugin_copy_sample;
            plugin->createSampleFnc =
            (PRESTypePluginCreateSampleFunction)
            ObstacleTypePlugin_create_sample;
            plugin->destroySampleFnc =
            (PRESTypePluginDestroySampleFunction)
            ObstacleTypePlugin_destroy_sample;

            plugin->serializeFnc = 
            (PRESTypePluginSerializeFunction) PRESTypePlugin_interpretedSerialize;
            plugin->deserializeFnc =
            (PRESTypePluginDeserializeFunction) PRESTypePlugin_interpretedDeserializeWithAlloc;
            plugin->getSerializedSampleMaxSizeFnc =
            (PRESTypePluginGetSerializedSampleMaxSizeFunction)
            ::ULARS::C2::ObstacleTypePlugin_get_serialized_sample_max_size;
            plugin->getSerializedSampleMinSizeFnc =
            (PRESTypePluginGetSerializedSampleMinSizeFunction)
            PRESTypePlugin_interpretedGetSerializedSampleMinSize;
            plugin->getDeserializedSampleMaxSizeFnc = NULL; 
            plugin->getSampleFnc =
            (PRESTypePluginGetSampleFunction)
            ObstacleTypePlugin_get_sample;
            plugin->returnSampleFnc =
            (PRESTypePluginReturnSampleFunction)
            ObstacleTypePlugin_return_sample;
            plugin->getKeyKindFnc =
            (PRESTypePluginGetKeyKindFunction)
            ::ULARS::C2::ObstacleTypePlugin_get_key_kind;

            /* These functions are only used for keyed types. As this is not a keyed
            type they are all set to NULL
            */
            plugin->serializeKeyFnc = NULL ;    
            plugin->deserializeKeyFnc = NULL;  
            plugin->getKeyFnc = NULL;
            plugin->returnKeyFnc = NULL;
            plugin->instanceToKeyFnc = NULL;
            plugin->keyToInstanceFnc = NULL;
            plugin->getSerializedKeyMaxSizeFnc = NULL;
            plugin->instanceToKeyHashFnc = NULL;
            plugin->serializedSampleToKeyHashFnc = NULL;
            plugin->serializedKeyToKeyHashFnc = NULL;    
            #ifdef NDDS_STANDALONE_TYPE
            plugin->typeCode = NULL; 
            #else
            try {
                plugin->typeCode = (struct RTICdrTypeCode *)
                &::rti::topic::dynamic_type< ::ULARS::C2::ObstacleType >::get().native();
            } catch (...) {
                ::ULARS::C2::ObstacleTypePlugin_delete(plugin);
                return NULL;
            }
            #endif
            plugin->languageKind = PRES_TYPEPLUGIN_CPPSTL_LANG;

            /* Serialized buffer */
            plugin->getBuffer = 
            (PRESTypePluginGetBufferFunction)
            ObstacleTypePlugin_get_buffer;
            plugin->returnBuffer = 
            (PRESTypePluginReturnBufferFunction)
            ObstacleTypePlugin_return_buffer;
            plugin->getBufferWithParams = NULL;
            plugin->returnBufferWithParams = NULL;
            plugin->getSerializedSampleSizeFnc =
            (PRESTypePluginGetSerializedSampleSizeFunction)
            PRESTypePlugin_interpretedGetSerializedSampleSize;

            plugin->getWriterLoanedSampleFnc = NULL; 
            plugin->returnWriterLoanedSampleFnc = NULL;
            plugin->returnWriterLoanedSampleFromCookieFnc = NULL;
            plugin->validateWriterLoanedSampleFnc = NULL;
            plugin->setWriterLoanedSampleSerializedStateFnc = NULL;

            static const char * TYPE_NAME = "ULARS::C2::ObstacleType";
            plugin->endpointTypeName = TYPE_NAME;
            plugin->isMetpType = RTI_FALSE;
            return plugin;
        }

        void
        ObstacleTypePlugin_delete(struct PRESTypePlugin *plugin)
        {
            RTIOsapiHeap_freeStructure(plugin);
        } 

        /* ----------------------------------------------------------------------------
        *  Type CommandInitialType
        * -------------------------------------------------------------------------- */

        /* -----------------------------------------------------------------------------
        Support functions:
        * -------------------------------------------------------------------------- */

        CommandInitialType *
        CommandInitialTypePluginSupport_create_data(void)
        {
            try {
                CommandInitialType *sample = new CommandInitialType();
                ::rti::topic::allocate_sample(*sample);
                return sample;
            } catch (...) {
                return NULL;
            }
        }

        void 
        CommandInitialTypePluginSupport_destroy_data(
            CommandInitialType *sample) 
        {
            delete sample;
        }

        RTIBool 
        CommandInitialTypePluginSupport_copy_data(
            CommandInitialType *dst,
            const CommandInitialType *src)
        {
            try {
                *dst = *src;
            } catch (...) {
                return RTI_FALSE;
            }

            return RTI_TRUE;
        }

        /* ----------------------------------------------------------------------------
        Callback functions:
        * ---------------------------------------------------------------------------- */

        PRESTypePluginParticipantData 
        CommandInitialTypePlugin_on_participant_attached(
            void *registration_data,
            const struct PRESTypePluginParticipantInfo *participant_info,
            RTIBool top_level_registration,
            void *container_plugin_context,
            RTICdrTypeCode *type_code)
        {
            try {
                struct RTIXCdrInterpreterPrograms *programs = NULL;
                struct PRESTypePluginDefaultParticipantData *pd = NULL;
                struct RTIXCdrInterpreterProgramsGenProperty programProperty =
                RTIXCdrInterpreterProgramsGenProperty_INITIALIZER;
                if (registration_data) {} /* To avoid warnings */
                if (participant_info) {} /* To avoid warnings */
                if (top_level_registration) {} /* To avoid warnings */
                if (container_plugin_context) {} /* To avoid warnings */
                if (type_code) {} /* To avoid warnings */
                pd = (struct PRESTypePluginDefaultParticipantData *)
                PRESTypePluginDefaultParticipantData_new(participant_info);

                programProperty.generateV1Encapsulation = RTI_XCDR_TRUE;
                programProperty.generateV2Encapsulation = RTI_XCDR_TRUE;
                programProperty.resolveAlias = RTI_XCDR_TRUE;
                programProperty.inlineStruct = RTI_XCDR_TRUE;
                programProperty.optimizeEnum = RTI_XCDR_TRUE;
                programProperty.unboundedSize = RTIXCdrLong_MAX;

                programProperty.externalReferenceSize = 
                (RTIXCdrUnsignedShort) sizeof(::dds::core::external<char>);
                programProperty.getExternalRefPointerFcn = 
                ::rti::topic::interpreter::get_external_value_pointer;

                programs = DDS_TypeCodeFactory_assert_programs_in_global_list(
                    DDS_TypeCodeFactory_get_instance(),
                    (DDS_TypeCode *) (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< CommandInitialType >::get().native()
                    ,
                    &programProperty,
                    RTI_XCDR_PROGRAM_MASK_TYPEPLUGIN);

                if (programs == NULL) {
                    PRESTypePluginDefaultParticipantData_delete(
                        (PRESTypePluginParticipantData)pd);
                    return NULL;
                }

                pd->programs = programs;
                return (PRESTypePluginParticipantData)pd;
            } catch (...) {
                return NULL;
            }
        }

        void 
        CommandInitialTypePlugin_on_participant_detached(
            PRESTypePluginParticipantData participant_data)
        {
            if (participant_data != NULL) {
                struct PRESTypePluginDefaultParticipantData *pd = 
                (struct PRESTypePluginDefaultParticipantData *)participant_data;

                if (pd->programs != NULL) {
                    DDS_TypeCodeFactory_remove_programs_from_global_list(
                        DDS_TypeCodeFactory_get_instance(),
                        pd->programs);
                    pd->programs = NULL;
                }
                PRESTypePluginDefaultParticipantData_delete(participant_data);
            }
        }

        PRESTypePluginEndpointData
        CommandInitialTypePlugin_on_endpoint_attached(
            PRESTypePluginParticipantData participant_data,
            const struct PRESTypePluginEndpointInfo *endpoint_info,
            RTIBool top_level_registration, 
            void *containerPluginContext)
        {
            try {
                PRESTypePluginEndpointData epd = NULL;
                unsigned int serializedSampleMaxSize = 0;

                if (top_level_registration) {} /* To avoid warnings */
                if (containerPluginContext) {} /* To avoid warnings */

                if (participant_data == NULL) {
                    return NULL;
                } 

                epd = PRESTypePluginDefaultEndpointData_new(
                    participant_data,
                    endpoint_info,
                    (PRESTypePluginDefaultEndpointDataCreateSampleFunction)
                    CommandInitialTypePluginSupport_create_data,
                    (PRESTypePluginDefaultEndpointDataDestroySampleFunction)
                    CommandInitialTypePluginSupport_destroy_data,
                    NULL , NULL );

                if (epd == NULL) {
                    return NULL;
                } 

                if (endpoint_info->endpointKind == PRES_TYPEPLUGIN_ENDPOINT_WRITER) {
                    serializedSampleMaxSize = ::ULARS::C2::CommandInitialTypePlugin_get_serialized_sample_max_size(
                        epd,RTI_FALSE,RTI_CDR_ENCAPSULATION_ID_CDR_BE,0);
                    PRESTypePluginDefaultEndpointData_setMaxSizeSerializedSample(epd, serializedSampleMaxSize);

                    if (PRESTypePluginDefaultEndpointData_createWriterPool(
                        epd,
                        endpoint_info,
                        (PRESTypePluginGetSerializedSampleMaxSizeFunction)
                        ::ULARS::C2::CommandInitialTypePlugin_get_serialized_sample_max_size, epd,
                        (PRESTypePluginGetSerializedSampleSizeFunction)
                        PRESTypePlugin_interpretedGetSerializedSampleSize,
                        epd) == RTI_FALSE) {
                        PRESTypePluginDefaultEndpointData_delete(epd);
                        return NULL;
                    }
                }

                return epd;
            } catch (...) {
                return NULL;
            }
        }

        void 
        CommandInitialTypePlugin_on_endpoint_detached(
            PRESTypePluginEndpointData endpoint_data)
        {  
            PRESTypePluginDefaultEndpointData_delete(endpoint_data);
        }

        void    
        CommandInitialTypePlugin_return_sample(
            PRESTypePluginEndpointData endpoint_data,
            CommandInitialType *sample,
            void *handle)
        {
            try {
                ::rti::topic::reset_sample(*sample);
            } catch(const std::exception& ex) {
                RTICdrLog_logWithFunctionName(
                    RTI_LOG_BIT_EXCEPTION,
                    "CommandInitialTypePlugin_return_sample",
                    &RTI_LOG_ANY_FAILURE_ss,
                    "exception: ",
                    ex.what());
            }

            PRESTypePluginDefaultEndpointData_returnSample(
                endpoint_data, sample, handle);
        }

        RTIBool 
        CommandInitialTypePlugin_copy_sample(
            PRESTypePluginEndpointData,
            CommandInitialType *dst,
            const CommandInitialType *src)
        {
            return ::ULARS::C2::CommandInitialTypePluginSupport_copy_data(dst,src);
        }

        /* ----------------------------------------------------------------------------
        (De)Serialize functions:
        * ------------------------------------------------------------------------- */
        unsigned int 
        CommandInitialTypePlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment);

        RTIBool
        CommandInitialTypePlugin_serialize_to_cdr_buffer(
            char * buffer,
            unsigned int * length,
            const CommandInitialType *sample,
            ::dds::core::policy::DataRepresentationId representation)
        {
            using namespace ::dds::core::policy;

            try{
                RTIEncapsulationId encapsulationId = RTI_CDR_ENCAPSULATION_ID_INVALID;
                struct RTICdrStream cdrStream;
                struct PRESTypePluginDefaultEndpointData epd;
                RTIBool result;
                struct PRESTypePluginDefaultParticipantData pd;
                struct RTIXCdrTypePluginProgramContext defaultProgramContext =
                RTIXCdrTypePluginProgramContext_INTIALIZER;
                struct PRESTypePlugin plugin = PRES_TYPEPLUGIN_DEFAULT;

                if (length == NULL) {
                    return RTI_FALSE;
                }

                RTIOsapiMemory_zero(&epd, sizeof(struct PRESTypePluginDefaultEndpointData));
                epd.programContext = defaultProgramContext;
                epd._participantData = &pd;
                epd.typePlugin = &plugin;
                epd.programContext.endpointPluginData = &epd;
                try {
                    plugin.typeCode = (struct RTICdrTypeCode *)
                    (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< CommandInitialType >::get().native()
                    ;
                } catch (...) {
                    return RTI_FALSE;
                }
                pd.programs = ::rti::topic::interpreter::get_cdr_serialization_programs<
                CommandInitialType, 
                true, true, true>();

                encapsulationId = DDS_TypeCode_get_native_encapsulation(
                    (DDS_TypeCode *) plugin.typeCode,
                    representation);

                if (encapsulationId == RTI_CDR_ENCAPSULATION_ID_INVALID) {
                    return RTI_FALSE;
                }

                epd._maxSizeSerializedSample =
                CommandInitialTypePlugin_get_serialized_sample_max_size(
                    (PRESTypePluginEndpointData)&epd, 
                    RTI_TRUE, 
                    encapsulationId,
                    0);

                if (buffer == NULL) {
                    *length = 
                    PRESTypePlugin_interpretedGetSerializedSampleSize(
                        (PRESTypePluginEndpointData)&epd,
                        RTI_TRUE,
                        encapsulationId,
                        0,
                        sample);

                    if (*length == 0) {
                        return RTI_FALSE;
                    }

                    return RTI_TRUE;
                }    

                RTICdrStream_init(&cdrStream);
                RTICdrStream_set(&cdrStream, buffer, *length);

                result = PRESTypePlugin_interpretedSerialize(
                    (PRESTypePluginEndpointData)&epd, 
                    sample, 
                    &cdrStream, 
                    RTI_TRUE, 
                    encapsulationId,
                    RTI_TRUE, 
                    NULL);  

                *length = (unsigned int) RTICdrStream_getCurrentPositionOffset(&cdrStream);
                return result;
            } catch (...) {
                return RTI_FALSE;
            }
        }

        RTIBool
        CommandInitialTypePlugin_deserialize_from_cdr_buffer(
            CommandInitialType *sample,
            const char * buffer,
            unsigned int length)
        {
            struct RTICdrStream cdrStream;
            struct PRESTypePluginDefaultParticipantData pd;
            struct RTIXCdrTypePluginProgramContext defaultProgramContext =
            RTIXCdrTypePluginProgramContext_INTIALIZER;
            struct PRESTypePlugin plugin;
            struct PRESTypePluginDefaultEndpointData epd;

            RTICdrStream_init(&cdrStream);
            /*
            * The buffer is constant because this is a deserialization function
            * (the buffer is an input parameter, not an output parameter).
            * However, the buffer member in the stream is a (char *) so coverity
            * complains in case something else modifies the buffer's contents later.
            *
            * We don't currently have a stream type with a constant buffer.
            * Therefore, we suppress the warning after making sure that this function
            * doesn't modify the contents of the stream's buffer.
            */
            /* coverity[cert_exp40_c_violation : FALSE] */
            RTICdrStream_set(&cdrStream, (char *) buffer, length);

            epd.programContext = defaultProgramContext;
            epd._participantData = &pd;
            epd.typePlugin = &plugin;
            epd.programContext.endpointPluginData = &epd;
            try {
                plugin.typeCode = (struct RTICdrTypeCode *)
                (RTIXCdrTypeCode *)&::rti::topic::dynamic_type< CommandInitialType >::get().native()
                ;
            } catch (...) {
                return RTI_FALSE;
            }
            pd.programs = ::rti::topic::interpreter::get_cdr_serialization_programs<
            CommandInitialType, 
            true, true, true>();

            epd._assignabilityProperty.acceptUnknownEnumValue = RTI_XCDR_TRUE;
            epd._assignabilityProperty.acceptUnknownUnionDiscriminator = 
            RTI_XCDR_ACCEPT_UNKNOWN_DISCRIMINATOR_AND_SELECT_DEFAULT;

            ::rti::topic::reset_sample(*sample);
            return PRESTypePlugin_interpretedDeserialize( 
                (PRESTypePluginEndpointData)&epd,
                sample,
                &cdrStream,
                RTI_TRUE, 
                RTI_TRUE, 
                NULL);
        }

        unsigned int 
        CommandInitialTypePlugin_get_serialized_sample_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            try {
                unsigned int size;
                RTIBool overflow = RTI_FALSE;

                size = PRESTypePlugin_interpretedGetSerializedSampleMaxSize(
                    endpoint_data,&overflow,include_encapsulation,encapsulation_id,current_alignment);

                if (overflow) {
                    size = RTI_CDR_MAX_SERIALIZED_SIZE;
                }

                return size;
            } catch (...) {
                return 0;
            }
        }

        /* --------------------------------------------------------------------------------------
        Key Management functions:
        * -------------------------------------------------------------------------------------- */

        PRESTypePluginKeyKind 
        CommandInitialTypePlugin_get_key_kind(void)
        {
            return PRES_TYPEPLUGIN_NO_KEY;
        }

        RTIBool CommandInitialTypePlugin_deserialize_key(
            PRESTypePluginEndpointData endpoint_data,
            CommandInitialType **sample, 
            RTIBool * drop_sample,
            struct RTICdrStream *cdrStream,
            RTIBool deserialize_encapsulation,
            RTIBool deserialize_key,
            void *endpoint_plugin_qos)
        {
            try {
                RTIBool result;
                if (drop_sample) {} /* To avoid warnings */
                cdrStream->_xTypesState.unassignable = RTI_FALSE;
                result= PRESTypePlugin_interpretedDeserializeKey( 
                    endpoint_data, (sample != NULL)?*sample:NULL, cdrStream,
                    deserialize_encapsulation, deserialize_key, endpoint_plugin_qos);
                if (result) {
                    if (cdrStream->_xTypesState.unassignable) {
                        result = RTI_FALSE;
                    }
                }
                return result;    
            } catch (...) {
                return RTI_FALSE;
            }     
        }

        unsigned int
        CommandInitialTypePlugin_get_serialized_key_max_size(
            PRESTypePluginEndpointData endpoint_data,
            RTIBool include_encapsulation,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            try {
                unsigned int size;
                RTIBool overflow = RTI_FALSE;

                size = PRESTypePlugin_interpretedGetSerializedKeyMaxSize(
                    endpoint_data,&overflow,include_encapsulation,encapsulation_id,current_alignment);
                if (overflow) {
                    size = RTI_CDR_MAX_SERIALIZED_SIZE;
                }

                return size;
            } catch (...) {
                return RTI_FALSE;
            }    
        }

        unsigned int
        CommandInitialTypePlugin_get_serialized_key_max_size_for_keyhash(
            PRESTypePluginEndpointData endpoint_data,
            RTIEncapsulationId encapsulation_id,
            unsigned int current_alignment)
        {
            unsigned int size;
            RTIBool overflow = RTI_FALSE;

            size = PRESTypePlugin_interpretedGetSerializedKeyMaxSizeForKeyhash(
                endpoint_data,
                &overflow,
                encapsulation_id,
                current_alignment);
            if (overflow) {
                size = RTI_CDR_MAX_SERIALIZED_SIZE;
            }

            return size;
        }

        /* ------------------------------------------------------------------------
        * Plug-in Installation Methods
        * ------------------------------------------------------------------------ */
        struct PRESTypePlugin *CommandInitialTypePlugin_new(void) 
        { 
            struct PRESTypePlugin *plugin = NULL;
            const struct PRESTypePluginVersion PLUGIN_VERSION = 
            PRES_TYPE_PLUGIN_VERSION_2_0;

            RTIOsapiHeap_allocateStructure(
                &plugin, struct PRESTypePlugin);
            if (plugin == NULL) {
                return NULL;
            }

            plugin->version = PLUGIN_VERSION;

            /* set up parent's function pointers */
            plugin->onParticipantAttached =
            (PRESTypePluginOnParticipantAttachedCallback)
            ::ULARS::C2::CommandInitialTypePlugin_on_participant_attached;
            plugin->onParticipantDetached =
            (PRESTypePluginOnParticipantDetachedCallback)
            ::ULARS::C2::CommandInitialTypePlugin_on_participant_detached;
            plugin->onEndpointAttached =
            (PRESTypePluginOnEndpointAttachedCallback)
            ::ULARS::C2::CommandInitialTypePlugin_on_endpoint_attached;
            plugin->onEndpointDetached =
            (PRESTypePluginOnEndpointDetachedCallback)
            ::ULARS::C2::CommandInitialTypePlugin_on_endpoint_detached;

            plugin->copySampleFnc =
            (PRESTypePluginCopySampleFunction)
            ::ULARS::C2::CommandInitialTypePlugin_copy_sample;
            plugin->createSampleFnc =
            (PRESTypePluginCreateSampleFunction)
            CommandInitialTypePlugin_create_sample;
            plugin->destroySampleFnc =
            (PRESTypePluginDestroySampleFunction)
            CommandInitialTypePlugin_destroy_sample;

            plugin->serializeFnc = 
            (PRESTypePluginSerializeFunction) PRESTypePlugin_interpretedSerialize;
            plugin->deserializeFnc =
            (PRESTypePluginDeserializeFunction) PRESTypePlugin_interpretedDeserializeWithAlloc;
            plugin->getSerializedSampleMaxSizeFnc =
            (PRESTypePluginGetSerializedSampleMaxSizeFunction)
            ::ULARS::C2::CommandInitialTypePlugin_get_serialized_sample_max_size;
            plugin->getSerializedSampleMinSizeFnc =
            (PRESTypePluginGetSerializedSampleMinSizeFunction)
            PRESTypePlugin_interpretedGetSerializedSampleMinSize;
            plugin->getDeserializedSampleMaxSizeFnc = NULL; 
            plugin->getSampleFnc =
            (PRESTypePluginGetSampleFunction)
            CommandInitialTypePlugin_get_sample;
            plugin->returnSampleFnc =
            (PRESTypePluginReturnSampleFunction)
            CommandInitialTypePlugin_return_sample;
            plugin->getKeyKindFnc =
            (PRESTypePluginGetKeyKindFunction)
            ::ULARS::C2::CommandInitialTypePlugin_get_key_kind;

            /* These functions are only used for keyed types. As this is not a keyed
            type they are all set to NULL
            */
            plugin->serializeKeyFnc = NULL ;    
            plugin->deserializeKeyFnc = NULL;  
            plugin->getKeyFnc = NULL;
            plugin->returnKeyFnc = NULL;
            plugin->instanceToKeyFnc = NULL;
            plugin->keyToInstanceFnc = NULL;
            plugin->getSerializedKeyMaxSizeFnc = NULL;
            plugin->instanceToKeyHashFnc = NULL;
            plugin->serializedSampleToKeyHashFnc = NULL;
            plugin->serializedKeyToKeyHashFnc = NULL;    
            #ifdef NDDS_STANDALONE_TYPE
            plugin->typeCode = NULL; 
            #else
            try {
                plugin->typeCode = (struct RTICdrTypeCode *)
                &::rti::topic::dynamic_type< ::ULARS::C2::CommandInitialType >::get().native();
            } catch (...) {
                ::ULARS::C2::CommandInitialTypePlugin_delete(plugin);
                return NULL;
            }
            #endif
            plugin->languageKind = PRES_TYPEPLUGIN_CPPSTL_LANG;

            /* Serialized buffer */
            plugin->getBuffer = 
            (PRESTypePluginGetBufferFunction)
            CommandInitialTypePlugin_get_buffer;
            plugin->returnBuffer = 
            (PRESTypePluginReturnBufferFunction)
            CommandInitialTypePlugin_return_buffer;
            plugin->getBufferWithParams = NULL;
            plugin->returnBufferWithParams = NULL;
            plugin->getSerializedSampleSizeFnc =
            (PRESTypePluginGetSerializedSampleSizeFunction)
            PRESTypePlugin_interpretedGetSerializedSampleSize;

            plugin->getWriterLoanedSampleFnc = NULL; 
            plugin->returnWriterLoanedSampleFnc = NULL;
            plugin->returnWriterLoanedSampleFromCookieFnc = NULL;
            plugin->validateWriterLoanedSampleFnc = NULL;
            plugin->setWriterLoanedSampleSerializedStateFnc = NULL;

            static const char * TYPE_NAME = "ULARS::C2::CommandInitialType";
            plugin->endpointTypeName = TYPE_NAME;
            plugin->isMetpType = RTI_FALSE;
            return plugin;
        }

        void
        CommandInitialTypePlugin_delete(struct PRESTypePlugin *plugin)
        {
            RTIOsapiHeap_freeStructure(plugin);
        } 
    } /* namespace C2  */
} /* namespace ULARS  */
#undef RTI_CDR_CURRENT_SUBMODULE 
