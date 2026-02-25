

/*
WARNING: THIS FILE IS AUTO-GENERATED. DO NOT MODIFY.

This file was generated from pils.idl
using RTI Code Generator (rtiddsgen) version 4.3.0.
The rtiddsgen tool is part of the RTI Connext DDS distribution.
For more information, type 'rtiddsgen -help' at a command shell
or consult the Code Generator User's Manual.
*/

#ifndef pils_475505405_hpp
#define pils_475505405_hpp

#include <iosfwd>

#if (defined(RTI_WIN32) || defined (RTI_WINCE) || defined(RTI_INTIME)) && defined(NDDS_USER_DLL_EXPORT)
/* If the code is building on Windows, start exporting symbols.
*/
#undef RTIUSERDllExport
#define RTIUSERDllExport __declspec(dllexport)
#endif

#include "dds/core/SafeEnumeration.hpp"
#include "dds/core/String.hpp"
#include "dds/core/array.hpp"
#include "dds/core/vector.hpp"
#include "dds/core/External.hpp"
#include "rti/core/LongDouble.hpp"
#include "rti/core/Pointer.hpp"
#include "rti/core/array.hpp"
#include "rti/topic/TopicTraits.hpp"

#include "omg/types/string_view.hpp"

#include "rti/core/BoundedSequence.hpp"
#include "dds/core/Optional.hpp"

#ifndef NDDS_STANDALONE_TYPE
#include "dds/domain/DomainParticipant.hpp"
#include "dds/topic/TopicTraits.hpp"
#include "dds/core/xtypes/DynamicType.hpp"
#include "dds/core/xtypes/StructType.hpp"
#include "dds/core/xtypes/UnionType.hpp"
#include "dds/core/xtypes/EnumType.hpp"
#include "dds/core/xtypes/AliasType.hpp"
#include "rti/util/StreamFlagSaver.hpp"
#include "rti/domain/PluginSupport.hpp"
#endif

#if (defined(RTI_WIN32) || defined (RTI_WINCE) || defined(RTI_INTIME)) && defined(NDDS_USER_DLL_EXPORT)
/* If the code is building on Windows, stop exporting symbols.
*/
#undef RTIUSERDllExport
#define RTIUSERDllExport
#endif

#if (defined(RTI_WIN32) || defined (RTI_WINCE) || defined(RTI_INTIME)) && defined(NDDS_USER_DLL_EXPORT)
/* If the code is building on Windows, start exporting symbols.
*/
#undef NDDSUSERDllExport
#define NDDSUSERDllExport __declspec(dllexport)
#endif

namespace ULARS {
    namespace EO {

        class NDDSUSERDllExport USVControlType {
          public:

            USVControlType();

            USVControlType(int16_t integratedTargetRPM_,int16_t integratedTargetSteering_,int16_t integratedBowThrust_);

            int16_t& integratedTargetRPM() noexcept {
                return m_integratedTargetRPM_;
            }

            const int16_t& integratedTargetRPM() const noexcept {
                return m_integratedTargetRPM_;
            }

            void integratedTargetRPM(int16_t value) {

                m_integratedTargetRPM_ = value;
            }

            int16_t& integratedTargetSteering() noexcept {
                return m_integratedTargetSteering_;
            }

            const int16_t& integratedTargetSteering() const noexcept {
                return m_integratedTargetSteering_;
            }

            void integratedTargetSteering(int16_t value) {

                m_integratedTargetSteering_ = value;
            }

            int16_t& integratedBowThrust() noexcept {
                return m_integratedBowThrust_;
            }

            const int16_t& integratedBowThrust() const noexcept {
                return m_integratedBowThrust_;
            }

            void integratedBowThrust(int16_t value) {

                m_integratedBowThrust_ = value;
            }

            bool operator == (const USVControlType& other_) const;
            bool operator != (const USVControlType& other_) const;

            void swap(USVControlType& other_) noexcept ;

          private:

            int16_t m_integratedTargetRPM_;
            int16_t m_integratedTargetSteering_;
            int16_t m_integratedBowThrust_;

        };

        inline void swap(USVControlType& a, USVControlType& b)  noexcept 
        {
            a.swap(b);
        }

        NDDSUSERDllExport std::ostream& operator<<(std::ostream& o, const USVControlType& sample);

    } // namespace EO  
    namespace C2 {

        class NDDSUSERDllExport Waypoint {
          public:

            Waypoint();

            Waypoint(double latitude_,double longitude_,uint32_t lateralErrorTolerance_,uint32_t goalErrorTolerance_,uint32_t maxVelocity_);

            double& latitude() noexcept {
                return m_latitude_;
            }

            const double& latitude() const noexcept {
                return m_latitude_;
            }

            void latitude(double value) {

                m_latitude_ = value;
            }

            double& longitude() noexcept {
                return m_longitude_;
            }

            const double& longitude() const noexcept {
                return m_longitude_;
            }

            void longitude(double value) {

                m_longitude_ = value;
            }

            uint32_t& lateralErrorTolerance() noexcept {
                return m_lateralErrorTolerance_;
            }

            const uint32_t& lateralErrorTolerance() const noexcept {
                return m_lateralErrorTolerance_;
            }

            void lateralErrorTolerance(uint32_t value) {

                m_lateralErrorTolerance_ = value;
            }

            uint32_t& goalErrorTolerance() noexcept {
                return m_goalErrorTolerance_;
            }

            const uint32_t& goalErrorTolerance() const noexcept {
                return m_goalErrorTolerance_;
            }

            void goalErrorTolerance(uint32_t value) {

                m_goalErrorTolerance_ = value;
            }

            uint32_t& maxVelocity() noexcept {
                return m_maxVelocity_;
            }

            const uint32_t& maxVelocity() const noexcept {
                return m_maxVelocity_;
            }

            void maxVelocity(uint32_t value) {

                m_maxVelocity_ = value;
            }

            bool operator == (const Waypoint& other_) const;
            bool operator != (const Waypoint& other_) const;

            void swap(Waypoint& other_) noexcept ;

          private:

            double m_latitude_;
            double m_longitude_;
            uint32_t m_lateralErrorTolerance_;
            uint32_t m_goalErrorTolerance_;
            uint32_t m_maxVelocity_;

        };

        inline void swap(Waypoint& a, Waypoint& b)  noexcept 
        {
            a.swap(b);
        }

        NDDSUSERDllExport std::ostream& operator<<(std::ostream& o, const Waypoint& sample);

        #if (defined(RTI_WIN32) || defined (RTI_WINCE)) && defined(NDDS_USER_DLL_EXPORT)
        // On Windows, dll-export template instantiations of standard types used by
        // other dll-exported types
        template class NDDSUSERDllExport std::allocator< ::ULARS::C2::Waypoint >;
        template class NDDSUSERDllExport std::vector< ::ULARS::C2::Waypoint >;
        #endif
        class NDDSUSERDllExport GlobalWaypointType {
          public:

            GlobalWaypointType();

            GlobalWaypointType(int16_t commandID_,const ::rti::core::bounded_sequence< ::ULARS::C2::Waypoint, 100L >& waypointData_);

            int16_t& commandID() noexcept {
                return m_commandID_;
            }

            const int16_t& commandID() const noexcept {
                return m_commandID_;
            }

            void commandID(int16_t value) {

                m_commandID_ = value;
            }

            ::rti::core::bounded_sequence< ::ULARS::C2::Waypoint, 100L >& waypointData() noexcept {
                return m_waypointData_;
            }

            const ::rti::core::bounded_sequence< ::ULARS::C2::Waypoint, 100L >& waypointData() const noexcept {
                return m_waypointData_;
            }

            void waypointData(const ::rti::core::bounded_sequence< ::ULARS::C2::Waypoint, 100L >& value) {

                m_waypointData_ = value;
            }

            void waypointData(::rti::core::bounded_sequence< ::ULARS::C2::Waypoint, 100L >&& value) {
                m_waypointData_ = std::move(value);
            }
            bool operator == (const GlobalWaypointType& other_) const;
            bool operator != (const GlobalWaypointType& other_) const;

            void swap(GlobalWaypointType& other_) noexcept ;

          private:

            int16_t m_commandID_;
            ::rti::core::bounded_sequence< ::ULARS::C2::Waypoint, 100L > m_waypointData_;

        };

        inline void swap(GlobalWaypointType& a, GlobalWaypointType& b)  noexcept 
        {
            a.swap(b);
        }

        NDDSUSERDllExport std::ostream& operator<<(std::ostream& o, const GlobalWaypointType& sample);

        class NDDSUSERDllExport StatusType {
          public:

            StatusType();

            StatusType(double USV_x_,double USV_y_,double USV_h_,double MSHIP_x_,double MSHIP_y_,double MSHIP_h_,double Craddle_x1_,double Craddle_y1_,double Craddle_x2_,double Craddle_y2_,double Craddle_x3_,double Craddle_y3_,double Craddle_x4_,double Craddle_y4_,double Craddle_h_);

            double& USV_x() noexcept {
                return m_USV_x_;
            }

            const double& USV_x() const noexcept {
                return m_USV_x_;
            }

            void USV_x(double value) {

                m_USV_x_ = value;
            }

            double& USV_y() noexcept {
                return m_USV_y_;
            }

            const double& USV_y() const noexcept {
                return m_USV_y_;
            }

            void USV_y(double value) {

                m_USV_y_ = value;
            }

            double& USV_h() noexcept {
                return m_USV_h_;
            }

            const double& USV_h() const noexcept {
                return m_USV_h_;
            }

            void USV_h(double value) {

                m_USV_h_ = value;
            }

            double& MSHIP_x() noexcept {
                return m_MSHIP_x_;
            }

            const double& MSHIP_x() const noexcept {
                return m_MSHIP_x_;
            }

            void MSHIP_x(double value) {

                m_MSHIP_x_ = value;
            }

            double& MSHIP_y() noexcept {
                return m_MSHIP_y_;
            }

            const double& MSHIP_y() const noexcept {
                return m_MSHIP_y_;
            }

            void MSHIP_y(double value) {

                m_MSHIP_y_ = value;
            }

            double& MSHIP_h() noexcept {
                return m_MSHIP_h_;
            }

            const double& MSHIP_h() const noexcept {
                return m_MSHIP_h_;
            }

            void MSHIP_h(double value) {

                m_MSHIP_h_ = value;
            }

            double& Craddle_x1() noexcept {
                return m_Craddle_x1_;
            }

            const double& Craddle_x1() const noexcept {
                return m_Craddle_x1_;
            }

            void Craddle_x1(double value) {

                m_Craddle_x1_ = value;
            }

            double& Craddle_y1() noexcept {
                return m_Craddle_y1_;
            }

            const double& Craddle_y1() const noexcept {
                return m_Craddle_y1_;
            }

            void Craddle_y1(double value) {

                m_Craddle_y1_ = value;
            }

            double& Craddle_x2() noexcept {
                return m_Craddle_x2_;
            }

            const double& Craddle_x2() const noexcept {
                return m_Craddle_x2_;
            }

            void Craddle_x2(double value) {

                m_Craddle_x2_ = value;
            }

            double& Craddle_y2() noexcept {
                return m_Craddle_y2_;
            }

            const double& Craddle_y2() const noexcept {
                return m_Craddle_y2_;
            }

            void Craddle_y2(double value) {

                m_Craddle_y2_ = value;
            }

            double& Craddle_x3() noexcept {
                return m_Craddle_x3_;
            }

            const double& Craddle_x3() const noexcept {
                return m_Craddle_x3_;
            }

            void Craddle_x3(double value) {

                m_Craddle_x3_ = value;
            }

            double& Craddle_y3() noexcept {
                return m_Craddle_y3_;
            }

            const double& Craddle_y3() const noexcept {
                return m_Craddle_y3_;
            }

            void Craddle_y3(double value) {

                m_Craddle_y3_ = value;
            }

            double& Craddle_x4() noexcept {
                return m_Craddle_x4_;
            }

            const double& Craddle_x4() const noexcept {
                return m_Craddle_x4_;
            }

            void Craddle_x4(double value) {

                m_Craddle_x4_ = value;
            }

            double& Craddle_y4() noexcept {
                return m_Craddle_y4_;
            }

            const double& Craddle_y4() const noexcept {
                return m_Craddle_y4_;
            }

            void Craddle_y4(double value) {

                m_Craddle_y4_ = value;
            }

            double& Craddle_h() noexcept {
                return m_Craddle_h_;
            }

            const double& Craddle_h() const noexcept {
                return m_Craddle_h_;
            }

            void Craddle_h(double value) {

                m_Craddle_h_ = value;
            }

            bool operator == (const StatusType& other_) const;
            bool operator != (const StatusType& other_) const;

            void swap(StatusType& other_) noexcept ;

          private:

            double m_USV_x_;
            double m_USV_y_;
            double m_USV_h_;
            double m_MSHIP_x_;
            double m_MSHIP_y_;
            double m_MSHIP_h_;
            double m_Craddle_x1_;
            double m_Craddle_y1_;
            double m_Craddle_x2_;
            double m_Craddle_y2_;
            double m_Craddle_x3_;
            double m_Craddle_y3_;
            double m_Craddle_x4_;
            double m_Craddle_y4_;
            double m_Craddle_h_;

        };

        inline void swap(StatusType& a, StatusType& b)  noexcept 
        {
            a.swap(b);
        }

        NDDSUSERDllExport std::ostream& operator<<(std::ostream& o, const StatusType& sample);

        class NDDSUSERDllExport ObstacleType {
          public:

            ObstacleType();

            ObstacleType(double Obs_x1_,double Obs_y1_,double Obs_x2_,double Obs_y2_,int16_t commandID_);

            double& Obs_x1() noexcept {
                return m_Obs_x1_;
            }

            const double& Obs_x1() const noexcept {
                return m_Obs_x1_;
            }

            void Obs_x1(double value) {

                m_Obs_x1_ = value;
            }

            double& Obs_y1() noexcept {
                return m_Obs_y1_;
            }

            const double& Obs_y1() const noexcept {
                return m_Obs_y1_;
            }

            void Obs_y1(double value) {

                m_Obs_y1_ = value;
            }

            double& Obs_x2() noexcept {
                return m_Obs_x2_;
            }

            const double& Obs_x2() const noexcept {
                return m_Obs_x2_;
            }

            void Obs_x2(double value) {

                m_Obs_x2_ = value;
            }

            double& Obs_y2() noexcept {
                return m_Obs_y2_;
            }

            const double& Obs_y2() const noexcept {
                return m_Obs_y2_;
            }

            void Obs_y2(double value) {

                m_Obs_y2_ = value;
            }

            int16_t& commandID() noexcept {
                return m_commandID_;
            }

            const int16_t& commandID() const noexcept {
                return m_commandID_;
            }

            void commandID(int16_t value) {

                m_commandID_ = value;
            }

            bool operator == (const ObstacleType& other_) const;
            bool operator != (const ObstacleType& other_) const;

            void swap(ObstacleType& other_) noexcept ;

          private:

            double m_Obs_x1_;
            double m_Obs_y1_;
            double m_Obs_x2_;
            double m_Obs_y2_;
            int16_t m_commandID_;

        };

        inline void swap(ObstacleType& a, ObstacleType& b)  noexcept 
        {
            a.swap(b);
        }

        NDDSUSERDllExport std::ostream& operator<<(std::ostream& o, const ObstacleType& sample);

        class NDDSUSERDllExport CommandInitialType {
          public:

            CommandInitialType();

            explicit CommandInitialType(int16_t commandID_);

            int16_t& commandID() noexcept {
                return m_commandID_;
            }

            const int16_t& commandID() const noexcept {
                return m_commandID_;
            }

            void commandID(int16_t value) {

                m_commandID_ = value;
            }

            bool operator == (const CommandInitialType& other_) const;
            bool operator != (const CommandInitialType& other_) const;

            void swap(CommandInitialType& other_) noexcept ;

          private:

            int16_t m_commandID_;

        };

        inline void swap(CommandInitialType& a, CommandInitialType& b)  noexcept 
        {
            a.swap(b);
        }

        NDDSUSERDllExport std::ostream& operator<<(std::ostream& o, const CommandInitialType& sample);

    } // namespace C2  
} // namespace ULARS  

#ifdef NDDS_STANDALONE_TYPE
namespace rti { 
    namespace topic {
    }
}
#else

namespace rti {
    namespace flat {
        namespace topic {
        }
    }
}
namespace dds {
    namespace topic {

        template<>
        struct topic_type_name< ::ULARS::EO::USVControlType > {
            NDDSUSERDllExport static std::string value() {
                return "ULARS::EO::USVControlType";
            }
        };

        template<>
        struct is_topic_type< ::ULARS::EO::USVControlType > : public ::dds::core::true_type {};

        template<>
        struct topic_type_support< ::ULARS::EO::USVControlType > {
            NDDSUSERDllExport 
            static void register_type(
                ::dds::domain::DomainParticipant& participant,
                const std::string & type_name);

            NDDSUSERDllExport 
            static std::vector<char>& to_cdr_buffer(
                std::vector<char>& buffer, 
                const ::ULARS::EO::USVControlType& sample,
                ::dds::core::policy::DataRepresentationId representation 
                = ::dds::core::policy::DataRepresentation::auto_id());

            NDDSUSERDllExport 
            static void from_cdr_buffer(::ULARS::EO::USVControlType& sample, const std::vector<char>& buffer);
            NDDSUSERDllExport 
            static void reset_sample(::ULARS::EO::USVControlType& sample);

            NDDSUSERDllExport 
            static void allocate_sample(::ULARS::EO::USVControlType& sample, int, int);

            static const ::rti::topic::TypePluginKind::type type_plugin_kind = 
            ::rti::topic::TypePluginKind::STL;
        };
        template<>
        struct topic_type_name< ::ULARS::C2::Waypoint > {
            NDDSUSERDllExport static std::string value() {
                return "ULARS::C2::Waypoint";
            }
        };

        template<>
        struct is_topic_type< ::ULARS::C2::Waypoint > : public ::dds::core::true_type {};

        template<>
        struct topic_type_support< ::ULARS::C2::Waypoint > {
            NDDSUSERDllExport 
            static void register_type(
                ::dds::domain::DomainParticipant& participant,
                const std::string & type_name);

            NDDSUSERDllExport 
            static std::vector<char>& to_cdr_buffer(
                std::vector<char>& buffer, 
                const ::ULARS::C2::Waypoint& sample,
                ::dds::core::policy::DataRepresentationId representation 
                = ::dds::core::policy::DataRepresentation::auto_id());

            NDDSUSERDllExport 
            static void from_cdr_buffer(::ULARS::C2::Waypoint& sample, const std::vector<char>& buffer);
            NDDSUSERDllExport 
            static void reset_sample(::ULARS::C2::Waypoint& sample);

            NDDSUSERDllExport 
            static void allocate_sample(::ULARS::C2::Waypoint& sample, int, int);

            static const ::rti::topic::TypePluginKind::type type_plugin_kind = 
            ::rti::topic::TypePluginKind::STL;
        };
        template<>
        struct topic_type_name< ::ULARS::C2::GlobalWaypointType > {
            NDDSUSERDllExport static std::string value() {
                return "ULARS::C2::GlobalWaypointType";
            }
        };

        template<>
        struct is_topic_type< ::ULARS::C2::GlobalWaypointType > : public ::dds::core::true_type {};

        template<>
        struct topic_type_support< ::ULARS::C2::GlobalWaypointType > {
            NDDSUSERDllExport 
            static void register_type(
                ::dds::domain::DomainParticipant& participant,
                const std::string & type_name);

            NDDSUSERDllExport 
            static std::vector<char>& to_cdr_buffer(
                std::vector<char>& buffer, 
                const ::ULARS::C2::GlobalWaypointType& sample,
                ::dds::core::policy::DataRepresentationId representation 
                = ::dds::core::policy::DataRepresentation::auto_id());

            NDDSUSERDllExport 
            static void from_cdr_buffer(::ULARS::C2::GlobalWaypointType& sample, const std::vector<char>& buffer);
            NDDSUSERDllExport 
            static void reset_sample(::ULARS::C2::GlobalWaypointType& sample);

            NDDSUSERDllExport 
            static void allocate_sample(::ULARS::C2::GlobalWaypointType& sample, int, int);

            static const ::rti::topic::TypePluginKind::type type_plugin_kind = 
            ::rti::topic::TypePluginKind::STL;
        };
        template<>
        struct topic_type_name< ::ULARS::C2::StatusType > {
            NDDSUSERDllExport static std::string value() {
                return "ULARS::C2::StatusType";
            }
        };

        template<>
        struct is_topic_type< ::ULARS::C2::StatusType > : public ::dds::core::true_type {};

        template<>
        struct topic_type_support< ::ULARS::C2::StatusType > {
            NDDSUSERDllExport 
            static void register_type(
                ::dds::domain::DomainParticipant& participant,
                const std::string & type_name);

            NDDSUSERDllExport 
            static std::vector<char>& to_cdr_buffer(
                std::vector<char>& buffer, 
                const ::ULARS::C2::StatusType& sample,
                ::dds::core::policy::DataRepresentationId representation 
                = ::dds::core::policy::DataRepresentation::auto_id());

            NDDSUSERDllExport 
            static void from_cdr_buffer(::ULARS::C2::StatusType& sample, const std::vector<char>& buffer);
            NDDSUSERDllExport 
            static void reset_sample(::ULARS::C2::StatusType& sample);

            NDDSUSERDllExport 
            static void allocate_sample(::ULARS::C2::StatusType& sample, int, int);

            static const ::rti::topic::TypePluginKind::type type_plugin_kind = 
            ::rti::topic::TypePluginKind::STL;
        };
        template<>
        struct topic_type_name< ::ULARS::C2::ObstacleType > {
            NDDSUSERDllExport static std::string value() {
                return "ULARS::C2::ObstacleType";
            }
        };

        template<>
        struct is_topic_type< ::ULARS::C2::ObstacleType > : public ::dds::core::true_type {};

        template<>
        struct topic_type_support< ::ULARS::C2::ObstacleType > {
            NDDSUSERDllExport 
            static void register_type(
                ::dds::domain::DomainParticipant& participant,
                const std::string & type_name);

            NDDSUSERDllExport 
            static std::vector<char>& to_cdr_buffer(
                std::vector<char>& buffer, 
                const ::ULARS::C2::ObstacleType& sample,
                ::dds::core::policy::DataRepresentationId representation 
                = ::dds::core::policy::DataRepresentation::auto_id());

            NDDSUSERDllExport 
            static void from_cdr_buffer(::ULARS::C2::ObstacleType& sample, const std::vector<char>& buffer);
            NDDSUSERDllExport 
            static void reset_sample(::ULARS::C2::ObstacleType& sample);

            NDDSUSERDllExport 
            static void allocate_sample(::ULARS::C2::ObstacleType& sample, int, int);

            static const ::rti::topic::TypePluginKind::type type_plugin_kind = 
            ::rti::topic::TypePluginKind::STL;
        };
        template<>
        struct topic_type_name< ::ULARS::C2::CommandInitialType > {
            NDDSUSERDllExport static std::string value() {
                return "ULARS::C2::CommandInitialType";
            }
        };

        template<>
        struct is_topic_type< ::ULARS::C2::CommandInitialType > : public ::dds::core::true_type {};

        template<>
        struct topic_type_support< ::ULARS::C2::CommandInitialType > {
            NDDSUSERDllExport 
            static void register_type(
                ::dds::domain::DomainParticipant& participant,
                const std::string & type_name);

            NDDSUSERDllExport 
            static std::vector<char>& to_cdr_buffer(
                std::vector<char>& buffer, 
                const ::ULARS::C2::CommandInitialType& sample,
                ::dds::core::policy::DataRepresentationId representation 
                = ::dds::core::policy::DataRepresentation::auto_id());

            NDDSUSERDllExport 
            static void from_cdr_buffer(::ULARS::C2::CommandInitialType& sample, const std::vector<char>& buffer);
            NDDSUSERDllExport 
            static void reset_sample(::ULARS::C2::CommandInitialType& sample);

            NDDSUSERDllExport 
            static void allocate_sample(::ULARS::C2::CommandInitialType& sample, int, int);

            static const ::rti::topic::TypePluginKind::type type_plugin_kind = 
            ::rti::topic::TypePluginKind::STL;
        };
    }
}

namespace rti { 
    namespace topic {

        template<>
        struct dynamic_type< ::ULARS::EO::USVControlType > {
            typedef ::dds::core::xtypes::StructType type;
            NDDSUSERDllExport static const ::dds::core::xtypes::StructType& get();
        };

        template <>
        struct extensibility< ::ULARS::EO::USVControlType > {
            static const ::dds::core::xtypes::ExtensibilityKind::type kind =
            ::dds::core::xtypes::ExtensibilityKind::EXTENSIBLE;    };

        template<>
        struct dynamic_type< ::ULARS::C2::Waypoint > {
            typedef ::dds::core::xtypes::StructType type;
            NDDSUSERDllExport static const ::dds::core::xtypes::StructType& get();
        };

        template <>
        struct extensibility< ::ULARS::C2::Waypoint > {
            static const ::dds::core::xtypes::ExtensibilityKind::type kind =
            ::dds::core::xtypes::ExtensibilityKind::EXTENSIBLE;    };

        template<>
        struct dynamic_type< ::ULARS::C2::GlobalWaypointType > {
            typedef ::dds::core::xtypes::StructType type;
            NDDSUSERDllExport static const ::dds::core::xtypes::StructType& get();
        };

        template <>
        struct extensibility< ::ULARS::C2::GlobalWaypointType > {
            static const ::dds::core::xtypes::ExtensibilityKind::type kind =
            ::dds::core::xtypes::ExtensibilityKind::EXTENSIBLE;    };

        template<>
        struct dynamic_type< ::ULARS::C2::StatusType > {
            typedef ::dds::core::xtypes::StructType type;
            NDDSUSERDllExport static const ::dds::core::xtypes::StructType& get();
        };

        template <>
        struct extensibility< ::ULARS::C2::StatusType > {
            static const ::dds::core::xtypes::ExtensibilityKind::type kind =
            ::dds::core::xtypes::ExtensibilityKind::EXTENSIBLE;    };

        template<>
        struct dynamic_type< ::ULARS::C2::ObstacleType > {
            typedef ::dds::core::xtypes::StructType type;
            NDDSUSERDllExport static const ::dds::core::xtypes::StructType& get();
        };

        template <>
        struct extensibility< ::ULARS::C2::ObstacleType > {
            static const ::dds::core::xtypes::ExtensibilityKind::type kind =
            ::dds::core::xtypes::ExtensibilityKind::EXTENSIBLE;    };

        template<>
        struct dynamic_type< ::ULARS::C2::CommandInitialType > {
            typedef ::dds::core::xtypes::StructType type;
            NDDSUSERDllExport static const ::dds::core::xtypes::StructType& get();
        };

        template <>
        struct extensibility< ::ULARS::C2::CommandInitialType > {
            static const ::dds::core::xtypes::ExtensibilityKind::type kind =
            ::dds::core::xtypes::ExtensibilityKind::EXTENSIBLE;    };

    }
}

#endif // NDDS_STANDALONE_TYPE
#if (defined(RTI_WIN32) || defined (RTI_WINCE) || defined(RTI_INTIME)) && defined(NDDS_USER_DLL_EXPORT)
/* If the code is building on Windows, stop exporting symbols.
*/
#undef NDDSUSERDllExport
#define NDDSUSERDllExport
#endif

#endif // pils_475505405_hpp

