

/*
WARNING: THIS FILE IS AUTO-GENERATED. DO NOT MODIFY.

This file was generated from pils.idl
using RTI Code Generator (rtiddsgen) version 4.3.0.
The rtiddsgen tool is part of the RTI Connext DDS distribution.
For more information, type 'rtiddsgen -help' at a command shell
or consult the Code Generator User's Manual.
*/

#include <iosfwd>
#include <iomanip>
#include <cmath>
#include <limits>

#ifndef NDDS_STANDALONE_TYPE
#include "rti/topic/cdr/Serialization.hpp"
#include "pilsPlugin.hpp"
#else
#include "rti/topic/cdr/SerializationHelpers.hpp"
#endif

#include "pils.hpp"

#include <rti/util/ostream_operators.hpp>

namespace ULARS {

    namespace EO {

        // ---- USVControlType: 

        USVControlType::USVControlType() :
            m_integratedTargetRPM_ (0) ,
            m_integratedTargetSteering_ (0) ,
            m_integratedBowThrust_ (0)  {

        }   

        USVControlType::USVControlType (int16_t integratedTargetRPM_,int16_t integratedTargetSteering_,int16_t integratedBowThrust_):
            m_integratedTargetRPM_(integratedTargetRPM_), 
            m_integratedTargetSteering_(integratedTargetSteering_), 
            m_integratedBowThrust_(integratedBowThrust_) {
        }

        void USVControlType::swap(USVControlType& other_)  noexcept 
        {
            using std::swap;
            swap(m_integratedTargetRPM_, other_.m_integratedTargetRPM_);
            swap(m_integratedTargetSteering_, other_.m_integratedTargetSteering_);
            swap(m_integratedBowThrust_, other_.m_integratedBowThrust_);
        }  

        bool USVControlType::operator == (const USVControlType& other_) const {
            if (m_integratedTargetRPM_ != other_.m_integratedTargetRPM_) {
                return false;
            }
            if (m_integratedTargetSteering_ != other_.m_integratedTargetSteering_) {
                return false;
            }
            if (m_integratedBowThrust_ != other_.m_integratedBowThrust_) {
                return false;
            }
            return true;
        }

        bool USVControlType::operator != (const USVControlType& other_) const {
            return !this->operator ==(other_);
        }

        std::ostream& operator << (std::ostream& o,const USVControlType& sample)
        {
            ::rti::util::StreamFlagSaver flag_saver (o);
            o <<"[";
            o << "integratedTargetRPM: " << sample.integratedTargetRPM ()<<", ";
            o << "integratedTargetSteering: " << sample.integratedTargetSteering ()<<", ";
            o << "integratedBowThrust: " << sample.integratedBowThrust ();
            o <<"]";
            return o;
        }

    } // namespace EO  

    namespace C2 {

        // ---- Waypoint: 

        Waypoint::Waypoint() :
            m_latitude_ (0.0) ,
            m_longitude_ (0.0) ,
            m_lateralErrorTolerance_ (0u) ,
            m_goalErrorTolerance_ (0u) ,
            m_maxVelocity_ (0u)  {

        }   

        Waypoint::Waypoint (double latitude_,double longitude_,uint32_t lateralErrorTolerance_,uint32_t goalErrorTolerance_,uint32_t maxVelocity_):
            m_latitude_(latitude_), 
            m_longitude_(longitude_), 
            m_lateralErrorTolerance_(lateralErrorTolerance_), 
            m_goalErrorTolerance_(goalErrorTolerance_), 
            m_maxVelocity_(maxVelocity_) {
        }

        void Waypoint::swap(Waypoint& other_)  noexcept 
        {
            using std::swap;
            swap(m_latitude_, other_.m_latitude_);
            swap(m_longitude_, other_.m_longitude_);
            swap(m_lateralErrorTolerance_, other_.m_lateralErrorTolerance_);
            swap(m_goalErrorTolerance_, other_.m_goalErrorTolerance_);
            swap(m_maxVelocity_, other_.m_maxVelocity_);
        }  

        bool Waypoint::operator == (const Waypoint& other_) const {
            if (std::fabs(m_latitude_ - other_.m_latitude_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_latitude_ - other_.m_latitude_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_longitude_ - other_.m_longitude_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_longitude_ - other_.m_longitude_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (m_lateralErrorTolerance_ != other_.m_lateralErrorTolerance_) {
                return false;
            }
            if (m_goalErrorTolerance_ != other_.m_goalErrorTolerance_) {
                return false;
            }
            if (m_maxVelocity_ != other_.m_maxVelocity_) {
                return false;
            }
            return true;
        }

        bool Waypoint::operator != (const Waypoint& other_) const {
            return !this->operator ==(other_);
        }

        std::ostream& operator << (std::ostream& o,const Waypoint& sample)
        {
            ::rti::util::StreamFlagSaver flag_saver (o);
            o <<"[";
            o << "latitude: " << std::setprecision(15) << sample.latitude ()<<", ";
            o << "longitude: " << std::setprecision(15) << sample.longitude ()<<", ";
            o << "lateralErrorTolerance: " << sample.lateralErrorTolerance ()<<", ";
            o << "goalErrorTolerance: " << sample.goalErrorTolerance ()<<", ";
            o << "maxVelocity: " << sample.maxVelocity ();
            o <<"]";
            return o;
        }

        // ---- GlobalWaypointType: 

        GlobalWaypointType::GlobalWaypointType() :
            m_commandID_ (0)  {

        }   

        GlobalWaypointType::GlobalWaypointType (int16_t commandID_,const ::rti::core::bounded_sequence< ::ULARS::C2::Waypoint, 100L >& waypointData_):
            m_commandID_(commandID_), 
            m_waypointData_(waypointData_) {
        }

        void GlobalWaypointType::swap(GlobalWaypointType& other_)  noexcept 
        {
            using std::swap;
            swap(m_commandID_, other_.m_commandID_);
            swap(m_waypointData_, other_.m_waypointData_);
        }  

        bool GlobalWaypointType::operator == (const GlobalWaypointType& other_) const {
            if (m_commandID_ != other_.m_commandID_) {
                return false;
            }
            if (m_waypointData_ != other_.m_waypointData_) {
                return false;
            }
            return true;
        }

        bool GlobalWaypointType::operator != (const GlobalWaypointType& other_) const {
            return !this->operator ==(other_);
        }

        std::ostream& operator << (std::ostream& o,const GlobalWaypointType& sample)
        {
            ::rti::util::StreamFlagSaver flag_saver (o);
            o <<"[";
            o << "commandID: " << sample.commandID ()<<", ";
            o << "waypointData: " << sample.waypointData ();
            o <<"]";
            return o;
        }

        // ---- StatusType: 

        StatusType::StatusType() :
            m_USV_x_ (0.0) ,
            m_USV_y_ (0.0) ,
            m_USV_h_ (0.0) ,
            m_MSHIP_x_ (0.0) ,
            m_MSHIP_y_ (0.0) ,
            m_MSHIP_h_ (0.0) ,
            m_Craddle_x1_ (0.0) ,
            m_Craddle_y1_ (0.0) ,
            m_Craddle_x2_ (0.0) ,
            m_Craddle_y2_ (0.0) ,
            m_Craddle_x3_ (0.0) ,
            m_Craddle_y3_ (0.0) ,
            m_Craddle_x4_ (0.0) ,
            m_Craddle_y4_ (0.0) ,
            m_Craddle_h_ (0.0)  {

        }   

        StatusType::StatusType (double USV_x_,double USV_y_,double USV_h_,double MSHIP_x_,double MSHIP_y_,double MSHIP_h_,double Craddle_x1_,double Craddle_y1_,double Craddle_x2_,double Craddle_y2_,double Craddle_x3_,double Craddle_y3_,double Craddle_x4_,double Craddle_y4_,double Craddle_h_):
            m_USV_x_(USV_x_), 
            m_USV_y_(USV_y_), 
            m_USV_h_(USV_h_), 
            m_MSHIP_x_(MSHIP_x_), 
            m_MSHIP_y_(MSHIP_y_), 
            m_MSHIP_h_(MSHIP_h_), 
            m_Craddle_x1_(Craddle_x1_), 
            m_Craddle_y1_(Craddle_y1_), 
            m_Craddle_x2_(Craddle_x2_), 
            m_Craddle_y2_(Craddle_y2_), 
            m_Craddle_x3_(Craddle_x3_), 
            m_Craddle_y3_(Craddle_y3_), 
            m_Craddle_x4_(Craddle_x4_), 
            m_Craddle_y4_(Craddle_y4_), 
            m_Craddle_h_(Craddle_h_) {
        }

        void StatusType::swap(StatusType& other_)  noexcept 
        {
            using std::swap;
            swap(m_USV_x_, other_.m_USV_x_);
            swap(m_USV_y_, other_.m_USV_y_);
            swap(m_USV_h_, other_.m_USV_h_);
            swap(m_MSHIP_x_, other_.m_MSHIP_x_);
            swap(m_MSHIP_y_, other_.m_MSHIP_y_);
            swap(m_MSHIP_h_, other_.m_MSHIP_h_);
            swap(m_Craddle_x1_, other_.m_Craddle_x1_);
            swap(m_Craddle_y1_, other_.m_Craddle_y1_);
            swap(m_Craddle_x2_, other_.m_Craddle_x2_);
            swap(m_Craddle_y2_, other_.m_Craddle_y2_);
            swap(m_Craddle_x3_, other_.m_Craddle_x3_);
            swap(m_Craddle_y3_, other_.m_Craddle_y3_);
            swap(m_Craddle_x4_, other_.m_Craddle_x4_);
            swap(m_Craddle_y4_, other_.m_Craddle_y4_);
            swap(m_Craddle_h_, other_.m_Craddle_h_);
        }  

        bool StatusType::operator == (const StatusType& other_) const {
            if (std::fabs(m_USV_x_ - other_.m_USV_x_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_USV_x_ - other_.m_USV_x_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_USV_y_ - other_.m_USV_y_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_USV_y_ - other_.m_USV_y_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_USV_h_ - other_.m_USV_h_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_USV_h_ - other_.m_USV_h_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_MSHIP_x_ - other_.m_MSHIP_x_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_MSHIP_x_ - other_.m_MSHIP_x_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_MSHIP_y_ - other_.m_MSHIP_y_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_MSHIP_y_ - other_.m_MSHIP_y_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_MSHIP_h_ - other_.m_MSHIP_h_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_MSHIP_h_ - other_.m_MSHIP_h_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_Craddle_x1_ - other_.m_Craddle_x1_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_Craddle_x1_ - other_.m_Craddle_x1_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_Craddle_y1_ - other_.m_Craddle_y1_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_Craddle_y1_ - other_.m_Craddle_y1_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_Craddle_x2_ - other_.m_Craddle_x2_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_Craddle_x2_ - other_.m_Craddle_x2_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_Craddle_y2_ - other_.m_Craddle_y2_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_Craddle_y2_ - other_.m_Craddle_y2_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_Craddle_x3_ - other_.m_Craddle_x3_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_Craddle_x3_ - other_.m_Craddle_x3_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_Craddle_y3_ - other_.m_Craddle_y3_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_Craddle_y3_ - other_.m_Craddle_y3_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_Craddle_x4_ - other_.m_Craddle_x4_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_Craddle_x4_ - other_.m_Craddle_x4_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_Craddle_y4_ - other_.m_Craddle_y4_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_Craddle_y4_ - other_.m_Craddle_y4_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_Craddle_h_ - other_.m_Craddle_h_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_Craddle_h_ - other_.m_Craddle_h_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            return true;
        }

        bool StatusType::operator != (const StatusType& other_) const {
            return !this->operator ==(other_);
        }

        std::ostream& operator << (std::ostream& o,const StatusType& sample)
        {
            ::rti::util::StreamFlagSaver flag_saver (o);
            o <<"[";
            o << "USV_x: " << std::setprecision(15) << sample.USV_x ()<<", ";
            o << "USV_y: " << std::setprecision(15) << sample.USV_y ()<<", ";
            o << "USV_h: " << std::setprecision(15) << sample.USV_h ()<<", ";
            o << "MSHIP_x: " << std::setprecision(15) << sample.MSHIP_x ()<<", ";
            o << "MSHIP_y: " << std::setprecision(15) << sample.MSHIP_y ()<<", ";
            o << "MSHIP_h: " << std::setprecision(15) << sample.MSHIP_h ()<<", ";
            o << "Craddle_x1: " << std::setprecision(15) << sample.Craddle_x1 ()<<", ";
            o << "Craddle_y1: " << std::setprecision(15) << sample.Craddle_y1 ()<<", ";
            o << "Craddle_x2: " << std::setprecision(15) << sample.Craddle_x2 ()<<", ";
            o << "Craddle_y2: " << std::setprecision(15) << sample.Craddle_y2 ()<<", ";
            o << "Craddle_x3: " << std::setprecision(15) << sample.Craddle_x3 ()<<", ";
            o << "Craddle_y3: " << std::setprecision(15) << sample.Craddle_y3 ()<<", ";
            o << "Craddle_x4: " << std::setprecision(15) << sample.Craddle_x4 ()<<", ";
            o << "Craddle_y4: " << std::setprecision(15) << sample.Craddle_y4 ()<<", ";
            o << "Craddle_h: " << std::setprecision(15) << sample.Craddle_h ();
            o <<"]";
            return o;
        }

        // ---- ObstacleType: 

        ObstacleType::ObstacleType() :
            m_Obs_x1_ (0.0) ,
            m_Obs_y1_ (0.0) ,
            m_Obs_x2_ (0.0) ,
            m_Obs_y2_ (0.0) ,
            m_commandID_ (0)  {

        }   

        ObstacleType::ObstacleType (double Obs_x1_,double Obs_y1_,double Obs_x2_,double Obs_y2_,int16_t commandID_):
            m_Obs_x1_(Obs_x1_), 
            m_Obs_y1_(Obs_y1_), 
            m_Obs_x2_(Obs_x2_), 
            m_Obs_y2_(Obs_y2_), 
            m_commandID_(commandID_) {
        }

        void ObstacleType::swap(ObstacleType& other_)  noexcept 
        {
            using std::swap;
            swap(m_Obs_x1_, other_.m_Obs_x1_);
            swap(m_Obs_y1_, other_.m_Obs_y1_);
            swap(m_Obs_x2_, other_.m_Obs_x2_);
            swap(m_Obs_y2_, other_.m_Obs_y2_);
            swap(m_commandID_, other_.m_commandID_);
        }  

        bool ObstacleType::operator == (const ObstacleType& other_) const {
            if (std::fabs(m_Obs_x1_ - other_.m_Obs_x1_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_Obs_x1_ - other_.m_Obs_x1_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_Obs_y1_ - other_.m_Obs_y1_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_Obs_y1_ - other_.m_Obs_y1_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_Obs_x2_ - other_.m_Obs_x2_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_Obs_x2_ - other_.m_Obs_x2_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (std::fabs(m_Obs_y2_ - other_.m_Obs_y2_) > std::numeric_limits< double>::epsilon()
            && !(std::fabs(m_Obs_y2_ - other_.m_Obs_y2_) < (std::numeric_limits< double>::min)())) {
                return false;
            }
            if (m_commandID_ != other_.m_commandID_) {
                return false;
            }
            return true;
        }

        bool ObstacleType::operator != (const ObstacleType& other_) const {
            return !this->operator ==(other_);
        }

        std::ostream& operator << (std::ostream& o,const ObstacleType& sample)
        {
            ::rti::util::StreamFlagSaver flag_saver (o);
            o <<"[";
            o << "Obs_x1: " << std::setprecision(15) << sample.Obs_x1 ()<<", ";
            o << "Obs_y1: " << std::setprecision(15) << sample.Obs_y1 ()<<", ";
            o << "Obs_x2: " << std::setprecision(15) << sample.Obs_x2 ()<<", ";
            o << "Obs_y2: " << std::setprecision(15) << sample.Obs_y2 ()<<", ";
            o << "commandID: " << sample.commandID ();
            o <<"]";
            return o;
        }

        // ---- CommandInitialType: 

        CommandInitialType::CommandInitialType() :
            m_commandID_ (0)  {

        }   

        CommandInitialType::CommandInitialType (int16_t commandID_):
            m_commandID_(commandID_) {
        }

        void CommandInitialType::swap(CommandInitialType& other_)  noexcept 
        {
            using std::swap;
            swap(m_commandID_, other_.m_commandID_);
        }  

        bool CommandInitialType::operator == (const CommandInitialType& other_) const {
            if (m_commandID_ != other_.m_commandID_) {
                return false;
            }
            return true;
        }

        bool CommandInitialType::operator != (const CommandInitialType& other_) const {
            return !this->operator ==(other_);
        }

        std::ostream& operator << (std::ostream& o,const CommandInitialType& sample)
        {
            ::rti::util::StreamFlagSaver flag_saver (o);
            o <<"[";
            o << "commandID: " << sample.commandID ();
            o <<"]";
            return o;
        }

    } // namespace C2  

} // namespace ULARS  

#ifdef NDDS_STANDALONE_TYPE
namespace rti {
    namespace topic {
    }
}

#else
// --- Type traits: -------------------------------------------------

namespace rti { 
    namespace topic {

        template<>
        struct native_type_code< ::ULARS::EO::USVControlType > {
            static DDS_TypeCode * get()
            {
                using namespace ::rti::topic::interpreter;

                static RTIBool is_initialized = RTI_FALSE;

                static DDS_TypeCode_Member USVControlType_g_tc_members[3]=
                {

                    {
                        (char *)"integratedTargetRPM",/* Member name */
                        {
                            0,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"integratedTargetSteering",/* Member name */
                        {
                            1,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"integratedBowThrust",/* Member name */
                        {
                            2,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }
                };

                static DDS_TypeCode USVControlType_g_tc =
                {{
                        DDS_TK_STRUCT, /* Kind */
                        DDS_BOOLEAN_FALSE, /* Ignored */
                        -1, /*Ignored*/
                        (char *)"ULARS::EO::USVControlType", /* Name */
                        NULL, /* Ignored */      
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        3, /* Number of members */
                        USVControlType_g_tc_members, /* Members */
                        DDS_VM_NONE, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER,
                        DDS_BOOLEAN_TRUE, /* _isCopyable */
                        NULL, /* _sampleAccessInfo: assigned later */
                        NULL /* _typePlugin: assigned later */
                    }}; /* Type code for USVControlType*/

                if (is_initialized) {
                    return &USVControlType_g_tc;
                }

                is_initialized = RTI_TRUE;

                USVControlType_g_tc._data._annotations._allowedDataRepresentationMask = 5;

                USVControlType_g_tc_members[0]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_short;
                USVControlType_g_tc_members[1]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_short;
                USVControlType_g_tc_members[2]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_short;

                /* Initialize the values for member annotations. */
                USVControlType_g_tc_members[0]._annotations._defaultValue._d = RTI_XCDR_TK_SHORT;
                USVControlType_g_tc_members[0]._annotations._defaultValue._u.short_value = 0;
                USVControlType_g_tc_members[0]._annotations._minValue._d = RTI_XCDR_TK_SHORT;
                USVControlType_g_tc_members[0]._annotations._minValue._u.short_value = RTIXCdrShort_MIN;
                USVControlType_g_tc_members[0]._annotations._maxValue._d = RTI_XCDR_TK_SHORT;
                USVControlType_g_tc_members[0]._annotations._maxValue._u.short_value = RTIXCdrShort_MAX;
                USVControlType_g_tc_members[1]._annotations._defaultValue._d = RTI_XCDR_TK_SHORT;
                USVControlType_g_tc_members[1]._annotations._defaultValue._u.short_value = 0;
                USVControlType_g_tc_members[1]._annotations._minValue._d = RTI_XCDR_TK_SHORT;
                USVControlType_g_tc_members[1]._annotations._minValue._u.short_value = RTIXCdrShort_MIN;
                USVControlType_g_tc_members[1]._annotations._maxValue._d = RTI_XCDR_TK_SHORT;
                USVControlType_g_tc_members[1]._annotations._maxValue._u.short_value = RTIXCdrShort_MAX;
                USVControlType_g_tc_members[2]._annotations._defaultValue._d = RTI_XCDR_TK_SHORT;
                USVControlType_g_tc_members[2]._annotations._defaultValue._u.short_value = 0;
                USVControlType_g_tc_members[2]._annotations._minValue._d = RTI_XCDR_TK_SHORT;
                USVControlType_g_tc_members[2]._annotations._minValue._u.short_value = RTIXCdrShort_MIN;
                USVControlType_g_tc_members[2]._annotations._maxValue._d = RTI_XCDR_TK_SHORT;
                USVControlType_g_tc_members[2]._annotations._maxValue._u.short_value = RTIXCdrShort_MAX;

                USVControlType_g_tc._data._sampleAccessInfo = sample_access_info();
                USVControlType_g_tc._data._typePlugin = type_plugin_info();    

                return &USVControlType_g_tc;
            }

            static RTIXCdrSampleAccessInfo * sample_access_info()
            {
                static RTIBool is_initialized = RTI_FALSE;

                ::ULARS::EO::USVControlType *sample;

                static RTIXCdrMemberAccessInfo USVControlType_g_memberAccessInfos[3] =
                {RTIXCdrMemberAccessInfo_INITIALIZER};

                static RTIXCdrSampleAccessInfo USVControlType_g_sampleAccessInfo = 
                RTIXCdrSampleAccessInfo_INITIALIZER;

                if (is_initialized) {
                    return (RTIXCdrSampleAccessInfo*) &USVControlType_g_sampleAccessInfo;
                }

                RTIXCdrHeap_allocateStruct(
                    &sample, 
                    ::ULARS::EO::USVControlType);
                if (sample == NULL) {
                    return NULL;
                }

                USVControlType_g_memberAccessInfos[0].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->integratedTargetRPM() - (char *)sample);

                USVControlType_g_memberAccessInfos[1].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->integratedTargetSteering() - (char *)sample);

                USVControlType_g_memberAccessInfos[2].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->integratedBowThrust() - (char *)sample);

                USVControlType_g_sampleAccessInfo.memberAccessInfos = 
                USVControlType_g_memberAccessInfos;

                {
                    size_t candidateTypeSize = sizeof(::ULARS::EO::USVControlType);

                    if (candidateTypeSize > RTIXCdrLong_MAX) {
                        USVControlType_g_sampleAccessInfo.typeSize[0] =
                        RTIXCdrLong_MAX;
                    } else {
                        USVControlType_g_sampleAccessInfo.typeSize[0] =
                        (RTIXCdrUnsignedLong) candidateTypeSize;
                    }
                }

                USVControlType_g_sampleAccessInfo.useGetMemberValueOnlyWithRef =
                RTI_XCDR_TRUE;

                USVControlType_g_sampleAccessInfo.getMemberValuePointerFcn = 
                interpreter::get_aggregation_value_pointer< ::ULARS::EO::USVControlType >;

                USVControlType_g_sampleAccessInfo.languageBinding = 
                RTI_XCDR_TYPE_BINDING_CPP_11_STL ;

                RTIXCdrHeap_freeStruct(sample);
                is_initialized = RTI_TRUE;
                return (RTIXCdrSampleAccessInfo*) &USVControlType_g_sampleAccessInfo;
            }
            static RTIXCdrTypePlugin * type_plugin_info()
            {
                static RTIXCdrTypePlugin USVControlType_g_typePlugin = 
                {
                    NULL, /* serialize */
                    NULL, /* serialize_key */
                    NULL, /* deserialize_sample */
                    NULL, /* deserialize_key_sample */
                    NULL, /* skip */
                    NULL, /* get_serialized_sample_size */
                    NULL, /* get_serialized_sample_max_size_ex */
                    NULL, /* get_serialized_key_max_size_ex */
                    NULL, /* get_serialized_sample_min_size */
                    NULL, /* serialized_sample_to_key */
                    NULL,
                    NULL,
                    NULL,
                    NULL,
                    NULL
                };

                return &USVControlType_g_typePlugin;
            }
        }; // native_type_code

        const ::dds::core::xtypes::StructType& dynamic_type< ::ULARS::EO::USVControlType >::get()
        {
            return static_cast<const ::dds::core::xtypes::StructType&>(
                ::rti::core::native_conversions::cast_from_native< ::dds::core::xtypes::DynamicType >(
                    *(native_type_code< ::ULARS::EO::USVControlType >::get())));
        }

        template<>
        struct native_type_code< ::ULARS::C2::Waypoint > {
            static DDS_TypeCode * get()
            {
                using namespace ::rti::topic::interpreter;

                static RTIBool is_initialized = RTI_FALSE;

                static DDS_TypeCode_Member Waypoint_g_tc_members[5]=
                {

                    {
                        (char *)"latitude",/* Member name */
                        {
                            0,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"longitude",/* Member name */
                        {
                            1,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"lateralErrorTolerance",/* Member name */
                        {
                            2,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"goalErrorTolerance",/* Member name */
                        {
                            3,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"maxVelocity",/* Member name */
                        {
                            4,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }
                };

                static DDS_TypeCode Waypoint_g_tc =
                {{
                        DDS_TK_STRUCT, /* Kind */
                        DDS_BOOLEAN_FALSE, /* Ignored */
                        -1, /*Ignored*/
                        (char *)"ULARS::C2::Waypoint", /* Name */
                        NULL, /* Ignored */      
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        5, /* Number of members */
                        Waypoint_g_tc_members, /* Members */
                        DDS_VM_NONE, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER,
                        DDS_BOOLEAN_TRUE, /* _isCopyable */
                        NULL, /* _sampleAccessInfo: assigned later */
                        NULL /* _typePlugin: assigned later */
                    }}; /* Type code for Waypoint*/

                if (is_initialized) {
                    return &Waypoint_g_tc;
                }

                is_initialized = RTI_TRUE;

                Waypoint_g_tc._data._annotations._allowedDataRepresentationMask = 5;

                Waypoint_g_tc_members[0]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                Waypoint_g_tc_members[1]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                Waypoint_g_tc_members[2]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_ulong;
                Waypoint_g_tc_members[3]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_ulong;
                Waypoint_g_tc_members[4]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_ulong;

                /* Initialize the values for member annotations. */
                Waypoint_g_tc_members[0]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                Waypoint_g_tc_members[0]._annotations._defaultValue._u.double_value = 0.0;
                Waypoint_g_tc_members[0]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                Waypoint_g_tc_members[0]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                Waypoint_g_tc_members[0]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                Waypoint_g_tc_members[0]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                Waypoint_g_tc_members[1]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                Waypoint_g_tc_members[1]._annotations._defaultValue._u.double_value = 0.0;
                Waypoint_g_tc_members[1]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                Waypoint_g_tc_members[1]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                Waypoint_g_tc_members[1]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                Waypoint_g_tc_members[1]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                Waypoint_g_tc_members[2]._annotations._defaultValue._d = RTI_XCDR_TK_ULONG;
                Waypoint_g_tc_members[2]._annotations._defaultValue._u.ulong_value = 0u;
                Waypoint_g_tc_members[2]._annotations._minValue._d = RTI_XCDR_TK_ULONG;
                Waypoint_g_tc_members[2]._annotations._minValue._u.ulong_value = RTIXCdrUnsignedLong_MIN;
                Waypoint_g_tc_members[2]._annotations._maxValue._d = RTI_XCDR_TK_ULONG;
                Waypoint_g_tc_members[2]._annotations._maxValue._u.ulong_value = RTIXCdrUnsignedLong_MAX;
                Waypoint_g_tc_members[3]._annotations._defaultValue._d = RTI_XCDR_TK_ULONG;
                Waypoint_g_tc_members[3]._annotations._defaultValue._u.ulong_value = 0u;
                Waypoint_g_tc_members[3]._annotations._minValue._d = RTI_XCDR_TK_ULONG;
                Waypoint_g_tc_members[3]._annotations._minValue._u.ulong_value = RTIXCdrUnsignedLong_MIN;
                Waypoint_g_tc_members[3]._annotations._maxValue._d = RTI_XCDR_TK_ULONG;
                Waypoint_g_tc_members[3]._annotations._maxValue._u.ulong_value = RTIXCdrUnsignedLong_MAX;
                Waypoint_g_tc_members[4]._annotations._defaultValue._d = RTI_XCDR_TK_ULONG;
                Waypoint_g_tc_members[4]._annotations._defaultValue._u.ulong_value = 0u;
                Waypoint_g_tc_members[4]._annotations._minValue._d = RTI_XCDR_TK_ULONG;
                Waypoint_g_tc_members[4]._annotations._minValue._u.ulong_value = RTIXCdrUnsignedLong_MIN;
                Waypoint_g_tc_members[4]._annotations._maxValue._d = RTI_XCDR_TK_ULONG;
                Waypoint_g_tc_members[4]._annotations._maxValue._u.ulong_value = RTIXCdrUnsignedLong_MAX;

                Waypoint_g_tc._data._sampleAccessInfo = sample_access_info();
                Waypoint_g_tc._data._typePlugin = type_plugin_info();    

                return &Waypoint_g_tc;
            }

            static RTIXCdrSampleAccessInfo * sample_access_info()
            {
                static RTIBool is_initialized = RTI_FALSE;

                ::ULARS::C2::Waypoint *sample;

                static RTIXCdrMemberAccessInfo Waypoint_g_memberAccessInfos[5] =
                {RTIXCdrMemberAccessInfo_INITIALIZER};

                static RTIXCdrSampleAccessInfo Waypoint_g_sampleAccessInfo = 
                RTIXCdrSampleAccessInfo_INITIALIZER;

                if (is_initialized) {
                    return (RTIXCdrSampleAccessInfo*) &Waypoint_g_sampleAccessInfo;
                }

                RTIXCdrHeap_allocateStruct(
                    &sample, 
                    ::ULARS::C2::Waypoint);
                if (sample == NULL) {
                    return NULL;
                }

                Waypoint_g_memberAccessInfos[0].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->latitude() - (char *)sample);

                Waypoint_g_memberAccessInfos[1].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->longitude() - (char *)sample);

                Waypoint_g_memberAccessInfos[2].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->lateralErrorTolerance() - (char *)sample);

                Waypoint_g_memberAccessInfos[3].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->goalErrorTolerance() - (char *)sample);

                Waypoint_g_memberAccessInfos[4].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->maxVelocity() - (char *)sample);

                Waypoint_g_sampleAccessInfo.memberAccessInfos = 
                Waypoint_g_memberAccessInfos;

                {
                    size_t candidateTypeSize = sizeof(::ULARS::C2::Waypoint);

                    if (candidateTypeSize > RTIXCdrLong_MAX) {
                        Waypoint_g_sampleAccessInfo.typeSize[0] =
                        RTIXCdrLong_MAX;
                    } else {
                        Waypoint_g_sampleAccessInfo.typeSize[0] =
                        (RTIXCdrUnsignedLong) candidateTypeSize;
                    }
                }

                Waypoint_g_sampleAccessInfo.useGetMemberValueOnlyWithRef =
                RTI_XCDR_TRUE;

                Waypoint_g_sampleAccessInfo.getMemberValuePointerFcn = 
                interpreter::get_aggregation_value_pointer< ::ULARS::C2::Waypoint >;

                Waypoint_g_sampleAccessInfo.languageBinding = 
                RTI_XCDR_TYPE_BINDING_CPP_11_STL ;

                RTIXCdrHeap_freeStruct(sample);
                is_initialized = RTI_TRUE;
                return (RTIXCdrSampleAccessInfo*) &Waypoint_g_sampleAccessInfo;
            }
            static RTIXCdrTypePlugin * type_plugin_info()
            {
                static RTIXCdrTypePlugin Waypoint_g_typePlugin = 
                {
                    NULL, /* serialize */
                    NULL, /* serialize_key */
                    NULL, /* deserialize_sample */
                    NULL, /* deserialize_key_sample */
                    NULL, /* skip */
                    NULL, /* get_serialized_sample_size */
                    NULL, /* get_serialized_sample_max_size_ex */
                    NULL, /* get_serialized_key_max_size_ex */
                    NULL, /* get_serialized_sample_min_size */
                    NULL, /* serialized_sample_to_key */
                    NULL,
                    NULL,
                    NULL,
                    NULL,
                    NULL
                };

                return &Waypoint_g_typePlugin;
            }
        }; // native_type_code

        const ::dds::core::xtypes::StructType& dynamic_type< ::ULARS::C2::Waypoint >::get()
        {
            return static_cast<const ::dds::core::xtypes::StructType&>(
                ::rti::core::native_conversions::cast_from_native< ::dds::core::xtypes::DynamicType >(
                    *(native_type_code< ::ULARS::C2::Waypoint >::get())));
        }

        template<>
        struct native_type_code< ::ULARS::C2::GlobalWaypointType > {
            static DDS_TypeCode * get()
            {
                using namespace ::rti::topic::interpreter;

                static RTIBool is_initialized = RTI_FALSE;

                static DDS_TypeCode GlobalWaypointType_g_tc_waypointData_sequence;

                static DDS_TypeCode_Member GlobalWaypointType_g_tc_members[2]=
                {

                    {
                        (char *)"commandID",/* Member name */
                        {
                            0,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"waypointData",/* Member name */
                        {
                            1,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }
                };

                static DDS_TypeCode GlobalWaypointType_g_tc =
                {{
                        DDS_TK_STRUCT, /* Kind */
                        DDS_BOOLEAN_FALSE, /* Ignored */
                        -1, /*Ignored*/
                        (char *)"ULARS::C2::GlobalWaypointType", /* Name */
                        NULL, /* Ignored */      
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        2, /* Number of members */
                        GlobalWaypointType_g_tc_members, /* Members */
                        DDS_VM_NONE, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER,
                        DDS_BOOLEAN_TRUE, /* _isCopyable */
                        NULL, /* _sampleAccessInfo: assigned later */
                        NULL /* _typePlugin: assigned later */
                    }}; /* Type code for GlobalWaypointType*/

                if (is_initialized) {
                    return &GlobalWaypointType_g_tc;
                }

                is_initialized = RTI_TRUE;

                GlobalWaypointType_g_tc_waypointData_sequence = initialize_sequence_typecode< ::rti::core::bounded_sequence< ::ULARS::C2::Waypoint, 100L > >((100L));

                GlobalWaypointType_g_tc._data._annotations._allowedDataRepresentationMask = 5;

                GlobalWaypointType_g_tc_waypointData_sequence._data._typeCode = (RTICdrTypeCode *)&::rti::topic::dynamic_type< ::ULARS::C2::Waypoint>::get().native();
                GlobalWaypointType_g_tc_members[0]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_short;
                GlobalWaypointType_g_tc_members[1]._representation._typeCode = (RTICdrTypeCode *)& GlobalWaypointType_g_tc_waypointData_sequence;

                /* Initialize the values for member annotations. */
                GlobalWaypointType_g_tc_members[0]._annotations._defaultValue._d = RTI_XCDR_TK_SHORT;
                GlobalWaypointType_g_tc_members[0]._annotations._defaultValue._u.short_value = 0;
                GlobalWaypointType_g_tc_members[0]._annotations._minValue._d = RTI_XCDR_TK_SHORT;
                GlobalWaypointType_g_tc_members[0]._annotations._minValue._u.short_value = RTIXCdrShort_MIN;
                GlobalWaypointType_g_tc_members[0]._annotations._maxValue._d = RTI_XCDR_TK_SHORT;
                GlobalWaypointType_g_tc_members[0]._annotations._maxValue._u.short_value = RTIXCdrShort_MAX;

                GlobalWaypointType_g_tc._data._sampleAccessInfo = sample_access_info();
                GlobalWaypointType_g_tc._data._typePlugin = type_plugin_info();    

                return &GlobalWaypointType_g_tc;
            }

            static RTIXCdrSampleAccessInfo * sample_access_info()
            {
                static RTIBool is_initialized = RTI_FALSE;

                ::ULARS::C2::GlobalWaypointType *sample;

                static RTIXCdrMemberAccessInfo GlobalWaypointType_g_memberAccessInfos[2] =
                {RTIXCdrMemberAccessInfo_INITIALIZER};

                static RTIXCdrSampleAccessInfo GlobalWaypointType_g_sampleAccessInfo = 
                RTIXCdrSampleAccessInfo_INITIALIZER;

                if (is_initialized) {
                    return (RTIXCdrSampleAccessInfo*) &GlobalWaypointType_g_sampleAccessInfo;
                }

                RTIXCdrHeap_allocateStruct(
                    &sample, 
                    ::ULARS::C2::GlobalWaypointType);
                if (sample == NULL) {
                    return NULL;
                }

                GlobalWaypointType_g_memberAccessInfos[0].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->commandID() - (char *)sample);

                GlobalWaypointType_g_memberAccessInfos[1].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->waypointData() - (char *)sample);

                GlobalWaypointType_g_sampleAccessInfo.memberAccessInfos = 
                GlobalWaypointType_g_memberAccessInfos;

                {
                    size_t candidateTypeSize = sizeof(::ULARS::C2::GlobalWaypointType);

                    if (candidateTypeSize > RTIXCdrLong_MAX) {
                        GlobalWaypointType_g_sampleAccessInfo.typeSize[0] =
                        RTIXCdrLong_MAX;
                    } else {
                        GlobalWaypointType_g_sampleAccessInfo.typeSize[0] =
                        (RTIXCdrUnsignedLong) candidateTypeSize;
                    }
                }

                GlobalWaypointType_g_sampleAccessInfo.useGetMemberValueOnlyWithRef =
                RTI_XCDR_TRUE;

                GlobalWaypointType_g_sampleAccessInfo.getMemberValuePointerFcn = 
                interpreter::get_aggregation_value_pointer< ::ULARS::C2::GlobalWaypointType >;

                GlobalWaypointType_g_sampleAccessInfo.languageBinding = 
                RTI_XCDR_TYPE_BINDING_CPP_11_STL ;

                RTIXCdrHeap_freeStruct(sample);
                is_initialized = RTI_TRUE;
                return (RTIXCdrSampleAccessInfo*) &GlobalWaypointType_g_sampleAccessInfo;
            }
            static RTIXCdrTypePlugin * type_plugin_info()
            {
                static RTIXCdrTypePlugin GlobalWaypointType_g_typePlugin = 
                {
                    NULL, /* serialize */
                    NULL, /* serialize_key */
                    NULL, /* deserialize_sample */
                    NULL, /* deserialize_key_sample */
                    NULL, /* skip */
                    NULL, /* get_serialized_sample_size */
                    NULL, /* get_serialized_sample_max_size_ex */
                    NULL, /* get_serialized_key_max_size_ex */
                    NULL, /* get_serialized_sample_min_size */
                    NULL, /* serialized_sample_to_key */
                    NULL,
                    NULL,
                    NULL,
                    NULL,
                    NULL
                };

                return &GlobalWaypointType_g_typePlugin;
            }
        }; // native_type_code

        const ::dds::core::xtypes::StructType& dynamic_type< ::ULARS::C2::GlobalWaypointType >::get()
        {
            return static_cast<const ::dds::core::xtypes::StructType&>(
                ::rti::core::native_conversions::cast_from_native< ::dds::core::xtypes::DynamicType >(
                    *(native_type_code< ::ULARS::C2::GlobalWaypointType >::get())));
        }

        template<>
        struct native_type_code< ::ULARS::C2::StatusType > {
            static DDS_TypeCode * get()
            {
                using namespace ::rti::topic::interpreter;

                static RTIBool is_initialized = RTI_FALSE;

                static DDS_TypeCode_Member StatusType_g_tc_members[15]=
                {

                    {
                        (char *)"USV_x",/* Member name */
                        {
                            0,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"USV_y",/* Member name */
                        {
                            1,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"USV_h",/* Member name */
                        {
                            2,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"MSHIP_x",/* Member name */
                        {
                            3,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"MSHIP_y",/* Member name */
                        {
                            4,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"MSHIP_h",/* Member name */
                        {
                            5,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"Craddle_x1",/* Member name */
                        {
                            6,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"Craddle_y1",/* Member name */
                        {
                            7,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"Craddle_x2",/* Member name */
                        {
                            8,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"Craddle_y2",/* Member name */
                        {
                            9,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"Craddle_x3",/* Member name */
                        {
                            10,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"Craddle_y3",/* Member name */
                        {
                            11,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"Craddle_x4",/* Member name */
                        {
                            12,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"Craddle_y4",/* Member name */
                        {
                            13,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"Craddle_h",/* Member name */
                        {
                            14,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }
                };

                static DDS_TypeCode StatusType_g_tc =
                {{
                        DDS_TK_STRUCT, /* Kind */
                        DDS_BOOLEAN_FALSE, /* Ignored */
                        -1, /*Ignored*/
                        (char *)"ULARS::C2::StatusType", /* Name */
                        NULL, /* Ignored */      
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        15, /* Number of members */
                        StatusType_g_tc_members, /* Members */
                        DDS_VM_NONE, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER,
                        DDS_BOOLEAN_TRUE, /* _isCopyable */
                        NULL, /* _sampleAccessInfo: assigned later */
                        NULL /* _typePlugin: assigned later */
                    }}; /* Type code for StatusType*/

                if (is_initialized) {
                    return &StatusType_g_tc;
                }

                is_initialized = RTI_TRUE;

                StatusType_g_tc._data._annotations._allowedDataRepresentationMask = 5;

                StatusType_g_tc_members[0]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                StatusType_g_tc_members[1]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                StatusType_g_tc_members[2]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                StatusType_g_tc_members[3]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                StatusType_g_tc_members[4]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                StatusType_g_tc_members[5]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                StatusType_g_tc_members[6]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                StatusType_g_tc_members[7]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                StatusType_g_tc_members[8]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                StatusType_g_tc_members[9]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                StatusType_g_tc_members[10]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                StatusType_g_tc_members[11]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                StatusType_g_tc_members[12]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                StatusType_g_tc_members[13]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                StatusType_g_tc_members[14]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;

                /* Initialize the values for member annotations. */
                StatusType_g_tc_members[0]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[0]._annotations._defaultValue._u.double_value = 0.0;
                StatusType_g_tc_members[0]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[0]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                StatusType_g_tc_members[0]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[0]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                StatusType_g_tc_members[1]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[1]._annotations._defaultValue._u.double_value = 0.0;
                StatusType_g_tc_members[1]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[1]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                StatusType_g_tc_members[1]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[1]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                StatusType_g_tc_members[2]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[2]._annotations._defaultValue._u.double_value = 0.0;
                StatusType_g_tc_members[2]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[2]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                StatusType_g_tc_members[2]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[2]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                StatusType_g_tc_members[3]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[3]._annotations._defaultValue._u.double_value = 0.0;
                StatusType_g_tc_members[3]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[3]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                StatusType_g_tc_members[3]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[3]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                StatusType_g_tc_members[4]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[4]._annotations._defaultValue._u.double_value = 0.0;
                StatusType_g_tc_members[4]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[4]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                StatusType_g_tc_members[4]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[4]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                StatusType_g_tc_members[5]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[5]._annotations._defaultValue._u.double_value = 0.0;
                StatusType_g_tc_members[5]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[5]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                StatusType_g_tc_members[5]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[5]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                StatusType_g_tc_members[6]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[6]._annotations._defaultValue._u.double_value = 0.0;
                StatusType_g_tc_members[6]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[6]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                StatusType_g_tc_members[6]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[6]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                StatusType_g_tc_members[7]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[7]._annotations._defaultValue._u.double_value = 0.0;
                StatusType_g_tc_members[7]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[7]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                StatusType_g_tc_members[7]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[7]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                StatusType_g_tc_members[8]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[8]._annotations._defaultValue._u.double_value = 0.0;
                StatusType_g_tc_members[8]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[8]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                StatusType_g_tc_members[8]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[8]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                StatusType_g_tc_members[9]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[9]._annotations._defaultValue._u.double_value = 0.0;
                StatusType_g_tc_members[9]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[9]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                StatusType_g_tc_members[9]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[9]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                StatusType_g_tc_members[10]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[10]._annotations._defaultValue._u.double_value = 0.0;
                StatusType_g_tc_members[10]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[10]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                StatusType_g_tc_members[10]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[10]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                StatusType_g_tc_members[11]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[11]._annotations._defaultValue._u.double_value = 0.0;
                StatusType_g_tc_members[11]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[11]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                StatusType_g_tc_members[11]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[11]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                StatusType_g_tc_members[12]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[12]._annotations._defaultValue._u.double_value = 0.0;
                StatusType_g_tc_members[12]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[12]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                StatusType_g_tc_members[12]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[12]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                StatusType_g_tc_members[13]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[13]._annotations._defaultValue._u.double_value = 0.0;
                StatusType_g_tc_members[13]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[13]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                StatusType_g_tc_members[13]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[13]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                StatusType_g_tc_members[14]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[14]._annotations._defaultValue._u.double_value = 0.0;
                StatusType_g_tc_members[14]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[14]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                StatusType_g_tc_members[14]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                StatusType_g_tc_members[14]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;

                StatusType_g_tc._data._sampleAccessInfo = sample_access_info();
                StatusType_g_tc._data._typePlugin = type_plugin_info();    

                return &StatusType_g_tc;
            }

            static RTIXCdrSampleAccessInfo * sample_access_info()
            {
                static RTIBool is_initialized = RTI_FALSE;

                ::ULARS::C2::StatusType *sample;

                static RTIXCdrMemberAccessInfo StatusType_g_memberAccessInfos[15] =
                {RTIXCdrMemberAccessInfo_INITIALIZER};

                static RTIXCdrSampleAccessInfo StatusType_g_sampleAccessInfo = 
                RTIXCdrSampleAccessInfo_INITIALIZER;

                if (is_initialized) {
                    return (RTIXCdrSampleAccessInfo*) &StatusType_g_sampleAccessInfo;
                }

                RTIXCdrHeap_allocateStruct(
                    &sample, 
                    ::ULARS::C2::StatusType);
                if (sample == NULL) {
                    return NULL;
                }

                StatusType_g_memberAccessInfos[0].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->USV_x() - (char *)sample);

                StatusType_g_memberAccessInfos[1].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->USV_y() - (char *)sample);

                StatusType_g_memberAccessInfos[2].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->USV_h() - (char *)sample);

                StatusType_g_memberAccessInfos[3].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->MSHIP_x() - (char *)sample);

                StatusType_g_memberAccessInfos[4].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->MSHIP_y() - (char *)sample);

                StatusType_g_memberAccessInfos[5].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->MSHIP_h() - (char *)sample);

                StatusType_g_memberAccessInfos[6].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->Craddle_x1() - (char *)sample);

                StatusType_g_memberAccessInfos[7].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->Craddle_y1() - (char *)sample);

                StatusType_g_memberAccessInfos[8].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->Craddle_x2() - (char *)sample);

                StatusType_g_memberAccessInfos[9].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->Craddle_y2() - (char *)sample);

                StatusType_g_memberAccessInfos[10].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->Craddle_x3() - (char *)sample);

                StatusType_g_memberAccessInfos[11].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->Craddle_y3() - (char *)sample);

                StatusType_g_memberAccessInfos[12].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->Craddle_x4() - (char *)sample);

                StatusType_g_memberAccessInfos[13].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->Craddle_y4() - (char *)sample);

                StatusType_g_memberAccessInfos[14].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->Craddle_h() - (char *)sample);

                StatusType_g_sampleAccessInfo.memberAccessInfos = 
                StatusType_g_memberAccessInfos;

                {
                    size_t candidateTypeSize = sizeof(::ULARS::C2::StatusType);

                    if (candidateTypeSize > RTIXCdrLong_MAX) {
                        StatusType_g_sampleAccessInfo.typeSize[0] =
                        RTIXCdrLong_MAX;
                    } else {
                        StatusType_g_sampleAccessInfo.typeSize[0] =
                        (RTIXCdrUnsignedLong) candidateTypeSize;
                    }
                }

                StatusType_g_sampleAccessInfo.useGetMemberValueOnlyWithRef =
                RTI_XCDR_TRUE;

                StatusType_g_sampleAccessInfo.getMemberValuePointerFcn = 
                interpreter::get_aggregation_value_pointer< ::ULARS::C2::StatusType >;

                StatusType_g_sampleAccessInfo.languageBinding = 
                RTI_XCDR_TYPE_BINDING_CPP_11_STL ;

                RTIXCdrHeap_freeStruct(sample);
                is_initialized = RTI_TRUE;
                return (RTIXCdrSampleAccessInfo*) &StatusType_g_sampleAccessInfo;
            }
            static RTIXCdrTypePlugin * type_plugin_info()
            {
                static RTIXCdrTypePlugin StatusType_g_typePlugin = 
                {
                    NULL, /* serialize */
                    NULL, /* serialize_key */
                    NULL, /* deserialize_sample */
                    NULL, /* deserialize_key_sample */
                    NULL, /* skip */
                    NULL, /* get_serialized_sample_size */
                    NULL, /* get_serialized_sample_max_size_ex */
                    NULL, /* get_serialized_key_max_size_ex */
                    NULL, /* get_serialized_sample_min_size */
                    NULL, /* serialized_sample_to_key */
                    NULL,
                    NULL,
                    NULL,
                    NULL,
                    NULL
                };

                return &StatusType_g_typePlugin;
            }
        }; // native_type_code

        const ::dds::core::xtypes::StructType& dynamic_type< ::ULARS::C2::StatusType >::get()
        {
            return static_cast<const ::dds::core::xtypes::StructType&>(
                ::rti::core::native_conversions::cast_from_native< ::dds::core::xtypes::DynamicType >(
                    *(native_type_code< ::ULARS::C2::StatusType >::get())));
        }

        template<>
        struct native_type_code< ::ULARS::C2::ObstacleType > {
            static DDS_TypeCode * get()
            {
                using namespace ::rti::topic::interpreter;

                static RTIBool is_initialized = RTI_FALSE;

                static DDS_TypeCode_Member ObstacleType_g_tc_members[5]=
                {

                    {
                        (char *)"Obs_x1",/* Member name */
                        {
                            0,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"Obs_y1",/* Member name */
                        {
                            1,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"Obs_x2",/* Member name */
                        {
                            2,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"Obs_y2",/* Member name */
                        {
                            3,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }, 
                    {
                        (char *)"commandID",/* Member name */
                        {
                            4,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }
                };

                static DDS_TypeCode ObstacleType_g_tc =
                {{
                        DDS_TK_STRUCT, /* Kind */
                        DDS_BOOLEAN_FALSE, /* Ignored */
                        -1, /*Ignored*/
                        (char *)"ULARS::C2::ObstacleType", /* Name */
                        NULL, /* Ignored */      
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        5, /* Number of members */
                        ObstacleType_g_tc_members, /* Members */
                        DDS_VM_NONE, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER,
                        DDS_BOOLEAN_TRUE, /* _isCopyable */
                        NULL, /* _sampleAccessInfo: assigned later */
                        NULL /* _typePlugin: assigned later */
                    }}; /* Type code for ObstacleType*/

                if (is_initialized) {
                    return &ObstacleType_g_tc;
                }

                is_initialized = RTI_TRUE;

                ObstacleType_g_tc._data._annotations._allowedDataRepresentationMask = 5;

                ObstacleType_g_tc_members[0]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                ObstacleType_g_tc_members[1]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                ObstacleType_g_tc_members[2]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                ObstacleType_g_tc_members[3]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_double;
                ObstacleType_g_tc_members[4]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_short;

                /* Initialize the values for member annotations. */
                ObstacleType_g_tc_members[0]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                ObstacleType_g_tc_members[0]._annotations._defaultValue._u.double_value = 0.0;
                ObstacleType_g_tc_members[0]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                ObstacleType_g_tc_members[0]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                ObstacleType_g_tc_members[0]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                ObstacleType_g_tc_members[0]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                ObstacleType_g_tc_members[1]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                ObstacleType_g_tc_members[1]._annotations._defaultValue._u.double_value = 0.0;
                ObstacleType_g_tc_members[1]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                ObstacleType_g_tc_members[1]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                ObstacleType_g_tc_members[1]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                ObstacleType_g_tc_members[1]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                ObstacleType_g_tc_members[2]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                ObstacleType_g_tc_members[2]._annotations._defaultValue._u.double_value = 0.0;
                ObstacleType_g_tc_members[2]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                ObstacleType_g_tc_members[2]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                ObstacleType_g_tc_members[2]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                ObstacleType_g_tc_members[2]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                ObstacleType_g_tc_members[3]._annotations._defaultValue._d = RTI_XCDR_TK_DOUBLE;
                ObstacleType_g_tc_members[3]._annotations._defaultValue._u.double_value = 0.0;
                ObstacleType_g_tc_members[3]._annotations._minValue._d = RTI_XCDR_TK_DOUBLE;
                ObstacleType_g_tc_members[3]._annotations._minValue._u.double_value = RTIXCdrDouble_MIN;
                ObstacleType_g_tc_members[3]._annotations._maxValue._d = RTI_XCDR_TK_DOUBLE;
                ObstacleType_g_tc_members[3]._annotations._maxValue._u.double_value = RTIXCdrDouble_MAX;
                ObstacleType_g_tc_members[4]._annotations._defaultValue._d = RTI_XCDR_TK_SHORT;
                ObstacleType_g_tc_members[4]._annotations._defaultValue._u.short_value = 0;
                ObstacleType_g_tc_members[4]._annotations._minValue._d = RTI_XCDR_TK_SHORT;
                ObstacleType_g_tc_members[4]._annotations._minValue._u.short_value = RTIXCdrShort_MIN;
                ObstacleType_g_tc_members[4]._annotations._maxValue._d = RTI_XCDR_TK_SHORT;
                ObstacleType_g_tc_members[4]._annotations._maxValue._u.short_value = RTIXCdrShort_MAX;

                ObstacleType_g_tc._data._sampleAccessInfo = sample_access_info();
                ObstacleType_g_tc._data._typePlugin = type_plugin_info();    

                return &ObstacleType_g_tc;
            }

            static RTIXCdrSampleAccessInfo * sample_access_info()
            {
                static RTIBool is_initialized = RTI_FALSE;

                ::ULARS::C2::ObstacleType *sample;

                static RTIXCdrMemberAccessInfo ObstacleType_g_memberAccessInfos[5] =
                {RTIXCdrMemberAccessInfo_INITIALIZER};

                static RTIXCdrSampleAccessInfo ObstacleType_g_sampleAccessInfo = 
                RTIXCdrSampleAccessInfo_INITIALIZER;

                if (is_initialized) {
                    return (RTIXCdrSampleAccessInfo*) &ObstacleType_g_sampleAccessInfo;
                }

                RTIXCdrHeap_allocateStruct(
                    &sample, 
                    ::ULARS::C2::ObstacleType);
                if (sample == NULL) {
                    return NULL;
                }

                ObstacleType_g_memberAccessInfos[0].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->Obs_x1() - (char *)sample);

                ObstacleType_g_memberAccessInfos[1].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->Obs_y1() - (char *)sample);

                ObstacleType_g_memberAccessInfos[2].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->Obs_x2() - (char *)sample);

                ObstacleType_g_memberAccessInfos[3].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->Obs_y2() - (char *)sample);

                ObstacleType_g_memberAccessInfos[4].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->commandID() - (char *)sample);

                ObstacleType_g_sampleAccessInfo.memberAccessInfos = 
                ObstacleType_g_memberAccessInfos;

                {
                    size_t candidateTypeSize = sizeof(::ULARS::C2::ObstacleType);

                    if (candidateTypeSize > RTIXCdrLong_MAX) {
                        ObstacleType_g_sampleAccessInfo.typeSize[0] =
                        RTIXCdrLong_MAX;
                    } else {
                        ObstacleType_g_sampleAccessInfo.typeSize[0] =
                        (RTIXCdrUnsignedLong) candidateTypeSize;
                    }
                }

                ObstacleType_g_sampleAccessInfo.useGetMemberValueOnlyWithRef =
                RTI_XCDR_TRUE;

                ObstacleType_g_sampleAccessInfo.getMemberValuePointerFcn = 
                interpreter::get_aggregation_value_pointer< ::ULARS::C2::ObstacleType >;

                ObstacleType_g_sampleAccessInfo.languageBinding = 
                RTI_XCDR_TYPE_BINDING_CPP_11_STL ;

                RTIXCdrHeap_freeStruct(sample);
                is_initialized = RTI_TRUE;
                return (RTIXCdrSampleAccessInfo*) &ObstacleType_g_sampleAccessInfo;
            }
            static RTIXCdrTypePlugin * type_plugin_info()
            {
                static RTIXCdrTypePlugin ObstacleType_g_typePlugin = 
                {
                    NULL, /* serialize */
                    NULL, /* serialize_key */
                    NULL, /* deserialize_sample */
                    NULL, /* deserialize_key_sample */
                    NULL, /* skip */
                    NULL, /* get_serialized_sample_size */
                    NULL, /* get_serialized_sample_max_size_ex */
                    NULL, /* get_serialized_key_max_size_ex */
                    NULL, /* get_serialized_sample_min_size */
                    NULL, /* serialized_sample_to_key */
                    NULL,
                    NULL,
                    NULL,
                    NULL,
                    NULL
                };

                return &ObstacleType_g_typePlugin;
            }
        }; // native_type_code

        const ::dds::core::xtypes::StructType& dynamic_type< ::ULARS::C2::ObstacleType >::get()
        {
            return static_cast<const ::dds::core::xtypes::StructType&>(
                ::rti::core::native_conversions::cast_from_native< ::dds::core::xtypes::DynamicType >(
                    *(native_type_code< ::ULARS::C2::ObstacleType >::get())));
        }

        template<>
        struct native_type_code< ::ULARS::C2::CommandInitialType > {
            static DDS_TypeCode * get()
            {
                using namespace ::rti::topic::interpreter;

                static RTIBool is_initialized = RTI_FALSE;

                static DDS_TypeCode_Member CommandInitialType_g_tc_members[1]=
                {

                    {
                        (char *)"commandID",/* Member name */
                        {
                            0,/* Representation ID */
                            DDS_BOOLEAN_FALSE,/* Is a pointer? */
                            -1, /* Bitfield bits */
                            NULL/* Member type code is assigned later */
                        },
                        0, /* Ignored */
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        RTI_CDR_REQUIRED_MEMBER, /* Is a key? */
                        DDS_PUBLIC_MEMBER,/* Member visibility */
                        1,
                        NULL, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER
                    }
                };

                static DDS_TypeCode CommandInitialType_g_tc =
                {{
                        DDS_TK_STRUCT, /* Kind */
                        DDS_BOOLEAN_FALSE, /* Ignored */
                        -1, /*Ignored*/
                        (char *)"ULARS::C2::CommandInitialType", /* Name */
                        NULL, /* Ignored */      
                        0, /* Ignored */
                        0, /* Ignored */
                        NULL, /* Ignored */
                        1, /* Number of members */
                        CommandInitialType_g_tc_members, /* Members */
                        DDS_VM_NONE, /* Ignored */
                        RTICdrTypeCodeAnnotations_INITIALIZER,
                        DDS_BOOLEAN_TRUE, /* _isCopyable */
                        NULL, /* _sampleAccessInfo: assigned later */
                        NULL /* _typePlugin: assigned later */
                    }}; /* Type code for CommandInitialType*/

                if (is_initialized) {
                    return &CommandInitialType_g_tc;
                }

                is_initialized = RTI_TRUE;

                CommandInitialType_g_tc._data._annotations._allowedDataRepresentationMask = 5;

                CommandInitialType_g_tc_members[0]._representation._typeCode = (RTICdrTypeCode *)&DDS_g_tc_short;

                /* Initialize the values for member annotations. */
                CommandInitialType_g_tc_members[0]._annotations._defaultValue._d = RTI_XCDR_TK_SHORT;
                CommandInitialType_g_tc_members[0]._annotations._defaultValue._u.short_value = 0;
                CommandInitialType_g_tc_members[0]._annotations._minValue._d = RTI_XCDR_TK_SHORT;
                CommandInitialType_g_tc_members[0]._annotations._minValue._u.short_value = RTIXCdrShort_MIN;
                CommandInitialType_g_tc_members[0]._annotations._maxValue._d = RTI_XCDR_TK_SHORT;
                CommandInitialType_g_tc_members[0]._annotations._maxValue._u.short_value = RTIXCdrShort_MAX;

                CommandInitialType_g_tc._data._sampleAccessInfo = sample_access_info();
                CommandInitialType_g_tc._data._typePlugin = type_plugin_info();    

                return &CommandInitialType_g_tc;
            }

            static RTIXCdrSampleAccessInfo * sample_access_info()
            {
                static RTIBool is_initialized = RTI_FALSE;

                ::ULARS::C2::CommandInitialType *sample;

                static RTIXCdrMemberAccessInfo CommandInitialType_g_memberAccessInfos[1] =
                {RTIXCdrMemberAccessInfo_INITIALIZER};

                static RTIXCdrSampleAccessInfo CommandInitialType_g_sampleAccessInfo = 
                RTIXCdrSampleAccessInfo_INITIALIZER;

                if (is_initialized) {
                    return (RTIXCdrSampleAccessInfo*) &CommandInitialType_g_sampleAccessInfo;
                }

                RTIXCdrHeap_allocateStruct(
                    &sample, 
                    ::ULARS::C2::CommandInitialType);
                if (sample == NULL) {
                    return NULL;
                }

                CommandInitialType_g_memberAccessInfos[0].bindingMemberValueOffset[0] = 
                (RTIXCdrUnsignedLong) ((char *)&sample->commandID() - (char *)sample);

                CommandInitialType_g_sampleAccessInfo.memberAccessInfos = 
                CommandInitialType_g_memberAccessInfos;

                {
                    size_t candidateTypeSize = sizeof(::ULARS::C2::CommandInitialType);

                    if (candidateTypeSize > RTIXCdrLong_MAX) {
                        CommandInitialType_g_sampleAccessInfo.typeSize[0] =
                        RTIXCdrLong_MAX;
                    } else {
                        CommandInitialType_g_sampleAccessInfo.typeSize[0] =
                        (RTIXCdrUnsignedLong) candidateTypeSize;
                    }
                }

                CommandInitialType_g_sampleAccessInfo.useGetMemberValueOnlyWithRef =
                RTI_XCDR_TRUE;

                CommandInitialType_g_sampleAccessInfo.getMemberValuePointerFcn = 
                interpreter::get_aggregation_value_pointer< ::ULARS::C2::CommandInitialType >;

                CommandInitialType_g_sampleAccessInfo.languageBinding = 
                RTI_XCDR_TYPE_BINDING_CPP_11_STL ;

                RTIXCdrHeap_freeStruct(sample);
                is_initialized = RTI_TRUE;
                return (RTIXCdrSampleAccessInfo*) &CommandInitialType_g_sampleAccessInfo;
            }
            static RTIXCdrTypePlugin * type_plugin_info()
            {
                static RTIXCdrTypePlugin CommandInitialType_g_typePlugin = 
                {
                    NULL, /* serialize */
                    NULL, /* serialize_key */
                    NULL, /* deserialize_sample */
                    NULL, /* deserialize_key_sample */
                    NULL, /* skip */
                    NULL, /* get_serialized_sample_size */
                    NULL, /* get_serialized_sample_max_size_ex */
                    NULL, /* get_serialized_key_max_size_ex */
                    NULL, /* get_serialized_sample_min_size */
                    NULL, /* serialized_sample_to_key */
                    NULL,
                    NULL,
                    NULL,
                    NULL,
                    NULL
                };

                return &CommandInitialType_g_typePlugin;
            }
        }; // native_type_code

        const ::dds::core::xtypes::StructType& dynamic_type< ::ULARS::C2::CommandInitialType >::get()
        {
            return static_cast<const ::dds::core::xtypes::StructType&>(
                ::rti::core::native_conversions::cast_from_native< ::dds::core::xtypes::DynamicType >(
                    *(native_type_code< ::ULARS::C2::CommandInitialType >::get())));
        }
    }
}

namespace dds { 
    namespace topic {
        void topic_type_support< ::ULARS::EO::USVControlType >:: register_type(
            ::dds::domain::DomainParticipant& participant,
            const std::string& type_name) 
        {

            ::rti::domain::register_type_plugin(
                participant,
                type_name,
                ::ULARS::EO::USVControlTypePlugin_new,
                ::ULARS::EO::USVControlTypePlugin_delete);
        }

        std::vector<char>& topic_type_support< ::ULARS::EO::USVControlType >::to_cdr_buffer(
            std::vector<char>& buffer, 
            const ::ULARS::EO::USVControlType& sample,
            ::dds::core::policy::DataRepresentationId representation)
        {
            // First get the length of the buffer
            unsigned int length = 0;
            RTIBool ok = USVControlTypePlugin_serialize_to_cdr_buffer(
                NULL, 
                &length,
                &sample,
                representation);
            ::rti::core::check_return_code(
                ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
                "Failed to calculate cdr buffer size");

            // Create a vector with that size and copy the cdr buffer into it
            buffer.resize(length);
            ok = USVControlTypePlugin_serialize_to_cdr_buffer(
                &buffer[0], 
                &length, 
                &sample,
                representation);
            ::rti::core::check_return_code(
                ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
                "Failed to copy cdr buffer");

            return buffer;
        }

        void topic_type_support< ::ULARS::EO::USVControlType >::from_cdr_buffer(::ULARS::EO::USVControlType& sample, 
        const std::vector<char>& buffer)
        {

            RTIBool ok  = USVControlTypePlugin_deserialize_from_cdr_buffer(
                &sample, 
                &buffer[0], 
                static_cast<unsigned int>(buffer.size()));
            ::rti::core::check_return_code(ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
            "Failed to create ::ULARS::EO::USVControlType from cdr buffer");
        }

        void topic_type_support< ::ULARS::EO::USVControlType >::reset_sample(::ULARS::EO::USVControlType& sample) 
        {
            sample.integratedTargetRPM(0);
            sample.integratedTargetSteering(0);
            sample.integratedBowThrust(0);
        }

        void topic_type_support< ::ULARS::EO::USVControlType >::allocate_sample(::ULARS::EO::USVControlType& sample, int, int) 
        {
            RTIOsapiUtility_unusedParameter(sample);
        }
        void topic_type_support< ::ULARS::C2::Waypoint >:: register_type(
            ::dds::domain::DomainParticipant& participant,
            const std::string& type_name) 
        {

            ::rti::domain::register_type_plugin(
                participant,
                type_name,
                ::ULARS::C2::WaypointPlugin_new,
                ::ULARS::C2::WaypointPlugin_delete);
        }

        std::vector<char>& topic_type_support< ::ULARS::C2::Waypoint >::to_cdr_buffer(
            std::vector<char>& buffer, 
            const ::ULARS::C2::Waypoint& sample,
            ::dds::core::policy::DataRepresentationId representation)
        {
            // First get the length of the buffer
            unsigned int length = 0;
            RTIBool ok = WaypointPlugin_serialize_to_cdr_buffer(
                NULL, 
                &length,
                &sample,
                representation);
            ::rti::core::check_return_code(
                ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
                "Failed to calculate cdr buffer size");

            // Create a vector with that size and copy the cdr buffer into it
            buffer.resize(length);
            ok = WaypointPlugin_serialize_to_cdr_buffer(
                &buffer[0], 
                &length, 
                &sample,
                representation);
            ::rti::core::check_return_code(
                ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
                "Failed to copy cdr buffer");

            return buffer;
        }

        void topic_type_support< ::ULARS::C2::Waypoint >::from_cdr_buffer(::ULARS::C2::Waypoint& sample, 
        const std::vector<char>& buffer)
        {

            RTIBool ok  = WaypointPlugin_deserialize_from_cdr_buffer(
                &sample, 
                &buffer[0], 
                static_cast<unsigned int>(buffer.size()));
            ::rti::core::check_return_code(ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
            "Failed to create ::ULARS::C2::Waypoint from cdr buffer");
        }

        void topic_type_support< ::ULARS::C2::Waypoint >::reset_sample(::ULARS::C2::Waypoint& sample) 
        {
            sample.latitude(0.0);
            sample.longitude(0.0);
            sample.lateralErrorTolerance(0u);
            sample.goalErrorTolerance(0u);
            sample.maxVelocity(0u);
        }

        void topic_type_support< ::ULARS::C2::Waypoint >::allocate_sample(::ULARS::C2::Waypoint& sample, int, int) 
        {
            RTIOsapiUtility_unusedParameter(sample);
        }
        void topic_type_support< ::ULARS::C2::GlobalWaypointType >:: register_type(
            ::dds::domain::DomainParticipant& participant,
            const std::string& type_name) 
        {

            ::rti::domain::register_type_plugin(
                participant,
                type_name,
                ::ULARS::C2::GlobalWaypointTypePlugin_new,
                ::ULARS::C2::GlobalWaypointTypePlugin_delete);
        }

        std::vector<char>& topic_type_support< ::ULARS::C2::GlobalWaypointType >::to_cdr_buffer(
            std::vector<char>& buffer, 
            const ::ULARS::C2::GlobalWaypointType& sample,
            ::dds::core::policy::DataRepresentationId representation)
        {
            // First get the length of the buffer
            unsigned int length = 0;
            RTIBool ok = GlobalWaypointTypePlugin_serialize_to_cdr_buffer(
                NULL, 
                &length,
                &sample,
                representation);
            ::rti::core::check_return_code(
                ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
                "Failed to calculate cdr buffer size");

            // Create a vector with that size and copy the cdr buffer into it
            buffer.resize(length);
            ok = GlobalWaypointTypePlugin_serialize_to_cdr_buffer(
                &buffer[0], 
                &length, 
                &sample,
                representation);
            ::rti::core::check_return_code(
                ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
                "Failed to copy cdr buffer");

            return buffer;
        }

        void topic_type_support< ::ULARS::C2::GlobalWaypointType >::from_cdr_buffer(::ULARS::C2::GlobalWaypointType& sample, 
        const std::vector<char>& buffer)
        {

            RTIBool ok  = GlobalWaypointTypePlugin_deserialize_from_cdr_buffer(
                &sample, 
                &buffer[0], 
                static_cast<unsigned int>(buffer.size()));
            ::rti::core::check_return_code(ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
            "Failed to create ::ULARS::C2::GlobalWaypointType from cdr buffer");
        }

        void topic_type_support< ::ULARS::C2::GlobalWaypointType >::reset_sample(::ULARS::C2::GlobalWaypointType& sample) 
        {
            sample.commandID(0);
            ::rti::topic::reset_sample(sample.waypointData());
        }

        void topic_type_support< ::ULARS::C2::GlobalWaypointType >::allocate_sample(::ULARS::C2::GlobalWaypointType& sample, int, int) 
        {
            ::rti::topic::allocate_sample(sample.waypointData(),  100L, -1);
        }
        void topic_type_support< ::ULARS::C2::StatusType >:: register_type(
            ::dds::domain::DomainParticipant& participant,
            const std::string& type_name) 
        {

            ::rti::domain::register_type_plugin(
                participant,
                type_name,
                ::ULARS::C2::StatusTypePlugin_new,
                ::ULARS::C2::StatusTypePlugin_delete);
        }

        std::vector<char>& topic_type_support< ::ULARS::C2::StatusType >::to_cdr_buffer(
            std::vector<char>& buffer, 
            const ::ULARS::C2::StatusType& sample,
            ::dds::core::policy::DataRepresentationId representation)
        {
            // First get the length of the buffer
            unsigned int length = 0;
            RTIBool ok = StatusTypePlugin_serialize_to_cdr_buffer(
                NULL, 
                &length,
                &sample,
                representation);
            ::rti::core::check_return_code(
                ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
                "Failed to calculate cdr buffer size");

            // Create a vector with that size and copy the cdr buffer into it
            buffer.resize(length);
            ok = StatusTypePlugin_serialize_to_cdr_buffer(
                &buffer[0], 
                &length, 
                &sample,
                representation);
            ::rti::core::check_return_code(
                ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
                "Failed to copy cdr buffer");

            return buffer;
        }

        void topic_type_support< ::ULARS::C2::StatusType >::from_cdr_buffer(::ULARS::C2::StatusType& sample, 
        const std::vector<char>& buffer)
        {

            RTIBool ok  = StatusTypePlugin_deserialize_from_cdr_buffer(
                &sample, 
                &buffer[0], 
                static_cast<unsigned int>(buffer.size()));
            ::rti::core::check_return_code(ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
            "Failed to create ::ULARS::C2::StatusType from cdr buffer");
        }

        void topic_type_support< ::ULARS::C2::StatusType >::reset_sample(::ULARS::C2::StatusType& sample) 
        {
            sample.USV_x(0.0);
            sample.USV_y(0.0);
            sample.USV_h(0.0);
            sample.MSHIP_x(0.0);
            sample.MSHIP_y(0.0);
            sample.MSHIP_h(0.0);
            sample.Craddle_x1(0.0);
            sample.Craddle_y1(0.0);
            sample.Craddle_x2(0.0);
            sample.Craddle_y2(0.0);
            sample.Craddle_x3(0.0);
            sample.Craddle_y3(0.0);
            sample.Craddle_x4(0.0);
            sample.Craddle_y4(0.0);
            sample.Craddle_h(0.0);
        }

        void topic_type_support< ::ULARS::C2::StatusType >::allocate_sample(::ULARS::C2::StatusType& sample, int, int) 
        {
            RTIOsapiUtility_unusedParameter(sample);
        }
        void topic_type_support< ::ULARS::C2::ObstacleType >:: register_type(
            ::dds::domain::DomainParticipant& participant,
            const std::string& type_name) 
        {

            ::rti::domain::register_type_plugin(
                participant,
                type_name,
                ::ULARS::C2::ObstacleTypePlugin_new,
                ::ULARS::C2::ObstacleTypePlugin_delete);
        }

        std::vector<char>& topic_type_support< ::ULARS::C2::ObstacleType >::to_cdr_buffer(
            std::vector<char>& buffer, 
            const ::ULARS::C2::ObstacleType& sample,
            ::dds::core::policy::DataRepresentationId representation)
        {
            // First get the length of the buffer
            unsigned int length = 0;
            RTIBool ok = ObstacleTypePlugin_serialize_to_cdr_buffer(
                NULL, 
                &length,
                &sample,
                representation);
            ::rti::core::check_return_code(
                ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
                "Failed to calculate cdr buffer size");

            // Create a vector with that size and copy the cdr buffer into it
            buffer.resize(length);
            ok = ObstacleTypePlugin_serialize_to_cdr_buffer(
                &buffer[0], 
                &length, 
                &sample,
                representation);
            ::rti::core::check_return_code(
                ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
                "Failed to copy cdr buffer");

            return buffer;
        }

        void topic_type_support< ::ULARS::C2::ObstacleType >::from_cdr_buffer(::ULARS::C2::ObstacleType& sample, 
        const std::vector<char>& buffer)
        {

            RTIBool ok  = ObstacleTypePlugin_deserialize_from_cdr_buffer(
                &sample, 
                &buffer[0], 
                static_cast<unsigned int>(buffer.size()));
            ::rti::core::check_return_code(ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
            "Failed to create ::ULARS::C2::ObstacleType from cdr buffer");
        }

        void topic_type_support< ::ULARS::C2::ObstacleType >::reset_sample(::ULARS::C2::ObstacleType& sample) 
        {
            sample.Obs_x1(0.0);
            sample.Obs_y1(0.0);
            sample.Obs_x2(0.0);
            sample.Obs_y2(0.0);
            sample.commandID(0);
        }

        void topic_type_support< ::ULARS::C2::ObstacleType >::allocate_sample(::ULARS::C2::ObstacleType& sample, int, int) 
        {
            RTIOsapiUtility_unusedParameter(sample);
        }
        void topic_type_support< ::ULARS::C2::CommandInitialType >:: register_type(
            ::dds::domain::DomainParticipant& participant,
            const std::string& type_name) 
        {

            ::rti::domain::register_type_plugin(
                participant,
                type_name,
                ::ULARS::C2::CommandInitialTypePlugin_new,
                ::ULARS::C2::CommandInitialTypePlugin_delete);
        }

        std::vector<char>& topic_type_support< ::ULARS::C2::CommandInitialType >::to_cdr_buffer(
            std::vector<char>& buffer, 
            const ::ULARS::C2::CommandInitialType& sample,
            ::dds::core::policy::DataRepresentationId representation)
        {
            // First get the length of the buffer
            unsigned int length = 0;
            RTIBool ok = CommandInitialTypePlugin_serialize_to_cdr_buffer(
                NULL, 
                &length,
                &sample,
                representation);
            ::rti::core::check_return_code(
                ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
                "Failed to calculate cdr buffer size");

            // Create a vector with that size and copy the cdr buffer into it
            buffer.resize(length);
            ok = CommandInitialTypePlugin_serialize_to_cdr_buffer(
                &buffer[0], 
                &length, 
                &sample,
                representation);
            ::rti::core::check_return_code(
                ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
                "Failed to copy cdr buffer");

            return buffer;
        }

        void topic_type_support< ::ULARS::C2::CommandInitialType >::from_cdr_buffer(::ULARS::C2::CommandInitialType& sample, 
        const std::vector<char>& buffer)
        {

            RTIBool ok  = CommandInitialTypePlugin_deserialize_from_cdr_buffer(
                &sample, 
                &buffer[0], 
                static_cast<unsigned int>(buffer.size()));
            ::rti::core::check_return_code(ok ? DDS_RETCODE_OK : DDS_RETCODE_ERROR,
            "Failed to create ::ULARS::C2::CommandInitialType from cdr buffer");
        }

        void topic_type_support< ::ULARS::C2::CommandInitialType >::reset_sample(::ULARS::C2::CommandInitialType& sample) 
        {
            sample.commandID(0);
        }

        void topic_type_support< ::ULARS::C2::CommandInitialType >::allocate_sample(::ULARS::C2::CommandInitialType& sample, int, int) 
        {
            RTIOsapiUtility_unusedParameter(sample);
        }
    }
}  

#endif // NDDS_STANDALONE_TYPE
