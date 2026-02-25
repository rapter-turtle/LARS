#ifndef DDS_COMMUNICATOR_HPP
#define DDS_COMMUNICATOR_HPP

#include "pils.hpp" // Generated from pils.idl
#include <dds/domain/DomainParticipant.hpp>
#include <dds/topic/Topic.hpp>
#include <dds/pub/Publisher.hpp>
#include <dds/sub/Subscriber.hpp>
#include <dds/pub/DataWriter.hpp>
#include <dds/sub/DataReader.hpp>
#include <dds/sub/Sample.hpp>

class DDSCommunicator
{
public:
    DDSCommunicator(int domain_id = 0);

    // No special destructor needed, C++11 objects manage their own lifecycle.
    ~DDSCommunicator() = default;

    bool init();

    // Publish methods
    bool publish_waypoints(const ULARS::C2::GlobalWaypointType& waypoints);
    bool publish_control(const ULARS::EO::USVControlType& control);

    // Subscriber access
    // Reads available samples from the DDS queue and returns them.
    std::vector<ULARS::C2::StatusType> read_status();

private:
    int domain_id_;

    // DDS C++11 API Entities
    dds::domain::DomainParticipant participant_{dds::core::null};
    dds::pub::Publisher publisher_{dds::core::null};
    dds::sub::Subscriber subscriber_{dds::core::null};

    // Topics
    dds::topic::Topic<ULARS::C2::GlobalWaypointType> waypoints_topic_{dds::core::null};
    dds::topic::Topic<ULARS::EO::USVControlType> control_topic_{dds::core::null};
    dds::topic::Topic<ULARS::C2::StatusType> status_topic_{dds::core::null};

    // DataWriters
    dds::pub::DataWriter<ULARS::C2::GlobalWaypointType> waypoints_writer_{dds::core::null};
    dds::pub::DataWriter<ULARS::EO::USVControlType> control_writer_{dds::core::null};

    // DataReader
    dds::sub::DataReader<ULARS::C2::StatusType> status_reader_{dds::core::null};
};

#endif // DDS_COMMUNICATOR_HPP
