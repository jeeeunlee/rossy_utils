#include <Eigen/Dense>

#include <osg/Geode>
#include <osg/Group>
#include <osg/ShapeDrawable>
#include <osg/Node>
#include <osg/PositionAttitudeTransform>

//osg_viewer
#include <osgViewer/Viewer>
#include <osgViewer/Renderer>
#include <osgViewer/ViewerEventHandlers>


class osgViewerWrapper {
    public:
        osgViewerWrapper();
        ~osgViewerWrapper();
        
        void addBox(const Eigen::VectorXf &center,
                    const Eigen::VectorXf &half_length,
                    const Eigen::Quaternionf &ori = Eigen::Quaternionf::Identity());
        void addCylinder(const Eigen::VectorXf &center,
                    const double &radius, 
                    const double &height, 
                    const Eigen::Quaternionf &ori = Eigen::Quaternionf::Identity());
        void addSphere(const Eigen::VectorXf &center,
                    const double &radius);
        void addSphereRed(const Eigen::VectorXf & center, 
                    const double & radius);

        int show();

    private:
        osg::ref_ptr<osgViewer::Viewer> viewer_;
        osg::ref_ptr<osg::Group> root_;
};

// Example code
// // create a simple shape
// osg::ref_ptr<osg::Geode> myGeode = new osg::Geode;
// myGeode->addDrawable(new osg::ShapeDrawable(new osg::Sphere()));

// // load obj
// osg::ref_ptr<osg::Node> myNode
//     = osgDB::readNodeFile("MyModel.obj");

// // attach nodes
// someNode->addChild(otherNode.get()); // if otherNode is a ref_ptr
// someNode->addChild(otherNode); // if otherNode is a raw pointer

// // Create transformation node
// osg::ref_ptr<osg::PositionAttitudeTransform> myPat
//     = new osg::PositionAttitudeTransform;
// myPat->setPosition(osg::Vec3(1.0, ydist, zdist));
// myPat->setScale(osg::Vec3(xscale, yscale, zscale));
// myPat->setAttitude(osg::Quat(angle_in_rads, axis_of_rot_vec));
// // Attach node to be transformed
// myPat->addChild(myNode.get());