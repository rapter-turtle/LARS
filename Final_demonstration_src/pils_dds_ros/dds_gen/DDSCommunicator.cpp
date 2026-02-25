#include "DDSCommunicator.hpp"
#include <iostream>

DDSCommunicator::DDSCommunicator(int domain_id) : domain_id_(domain_id) {}

bool DDSCommunicator::init()
{
    try {
        // --- Create Participant, Publisher, and Subscriber ---
        participant_ = dds::domain::DomainParticipant(domain_id_);
        publisher_ = dds::pub::Publisher(participant_);
        subscriber_ = dds::sub::Subscriber(participant_);

        // --- Create Topics ---
        // The modern C++ API finds the type support automatically.
        waypoints_topic_ = dds::topic::Topic<ULARS::C2::GlobalWaypointType>(participant_, "GlobalWaypointType");
        control_topic_ = dds::topic::Topic<ULARS::EO::USVControlType>(participant_, "USVControlType");
        status_topic_ = dds::topic::Topic<ULARS::C2::StatusType>(participant_, "StatusType");

        // --- Create DataWriters ---
        waypoints_writer_ = dds::pub::DataWriter<ULARS::C2::GlobalWaypointType>(publisher_, waypoints_topic_);
        control_writer_ = dds::pub::DataWriter<ULARS::EO::USVControlType>(publisher_, control_topic_);

        // --- Create DataReader ---
        status_reader_ = dds::sub::DataReader<ULARS::C2::StatusType>(subscriber_, status_topic_);

    } catch (const dds::core::Exception& e) {
        std::cerr << "DDSCommunicator initialization failed: " << e.what() << std::endl;
        return false;
    }

    std::cout << "DDSCommunicator initialized successfully." << std::endl;
    return true;
}

bool DDSCommunicator::publish_waypoints(const ULARS::C2::GlobalWaypointType& waypoints)
{
    try {
        waypoints_writer_.write(waypoints);
        return true;
    } catch (const dds::core::Exception& e) {
        std::cerr << "Failed to publish waypoints: " << e.what() << std::endl;
        return false;
    }
}

bool DDSCommunicator::publish_control(const ULARS::EO::USVControlType& control)
{
    try {
        control_writer_.write(control);
        return true;
    } catch (const dds::core::Exception& e) {
        std::cerr << "Failed to publish control: " << e.what() << std::endl;
        return false;
    }
}

std::vector<ULARS::C2::StatusType> DDSCommunicator::read_status()
{
    std::vector<ULARS::C2::StatusType> received_data;
    
    // Take all available samples from the reader
    dds::sub::LoanedSamples<ULARS::C2::StatusType> samples = status_reader_.take();

    for (const auto& sample : samples) {
        if (sample.info().valid()) {
            received_data.push_back(sample.data());
        }
    }

    // LoanedSamples automatically returns the loan upon destruction
    return received_data;
}
