#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/Node>
#include <osg/PositionAttitudeTransform>

#include <osgViewer/Viewer>
#include <osgViewer/ViewerEventHandlers>
#include "osgViewerWrapper.hpp"

//osgGA
#include <osgGA/GUIEventAdapter>
#include <osgGA/TrackballManipulator>
#include <osgGA/StateSetManipulator>

osgViewerWrapper::osgViewerWrapper()
{
    // rootNode_ = osg::ref_ptr<osg::Node>
    root_ = new osg::Group();
    viewer_ = new osgViewer::Viewer();
}

osgViewerWrapper::~osgViewerWrapper()
{
    for(unsigned int j=root_->getNumChildren(); j>0; j--){
        osg::ref_ptr<osg::Node> children = root_->getChild(j-1);
        root_->removeChild(children);
    }
    delete viewer_;
}


void osgViewerWrapper::addBox(const Eigen::VectorXf & center, 
                    const Eigen::VectorXf & hl, 
                    const Eigen::Quaternionf & ori)
{
    // osg::Vec4 red(0.90f,0.30f,0.30f,1.0f);
    osg::Vec4 lightblue(0.30f,0.6f,0.90f,0.3f);
    // osg::Vec4 blue(0.10f,0.30f,0.40f,1.0f);

    osg::Geode* geode = new osg::Geode;
    osg::ref_ptr<osg::TessellationHints> hints = new osg::TessellationHints;
    hints->setDetailRatio(2.0f);

    osg::Box* box = new osg::Box();
    box->setCenter(osg::Vec3(center(0),center(1),center(2)));
    box->setHalfLengths(osg::Vec3(hl(0),hl(1),hl(2)));
    box->setRotation(osg::Quat(ori.x(),ori.y(),ori.z(),ori.w()));
    osg::ShapeDrawable* sd = new osg::ShapeDrawable(box, hints);    
    sd->setColor(lightblue);    
    sd->getOrCreateStateSet()->setMode( GL_BLEND, osg::StateAttribute::ON );
    sd->getOrCreateStateSet()->setRenderingHint( osg::StateSet::TRANSPARENT_BIN );
    
    geode->addDrawable(sd);    
    root_->addChild(geode);
}

void osgViewerWrapper::addCylinder(const Eigen::VectorXf & center, 
                                const double & radius, 
                                const double & height, 
                                const Eigen::Quaternionf & ori)
{
    osg::Geode* geode = new osg::Geode;

    osg::Cylinder* cyl = new osg::Cylinder();
    cyl->setCenter(osg::Vec3(center(0),center(1),center(2)));
    cyl->setRadius(radius);
    cyl->setHeight(height);
    cyl->setRotation(osg::Quat(ori.x(),ori.y(),ori.z(),ori.w()));
    osg::ShapeDrawable* sd = new osg::ShapeDrawable(cyl);
    geode->addDrawable(sd);
    root_->addChild(geode);
}

void osgViewerWrapper::addSphere(const Eigen::VectorXf & center, 
                            const double & radius)
{
    osg::Geode* geode = new osg::Geode;
    osg::Sphere* sph = new osg::Sphere();
    sph->setCenter(osg::Vec3(center(0),center(1),center(2)));
    sph->setRadius(radius);
    osg::ShapeDrawable* sd = new osg::ShapeDrawable(sph);
    geode->addDrawable(sd);
    root_->addChild(geode);
}

void osgViewerWrapper::addSphereRed(const Eigen::VectorXf & center, 
                            const double & radius)
{
    osg::Vec4 red(0.90f,0.30f,0.30f,1.0f);
    // osg::Vec4 lightblue(0.30f,0.6f,0.90f,1.0f);
    // osg::Vec4 blue(0.10f,0.30f,0.40f,1.0f);

    osg::Geode* geode = new osg::Geode;
    osg::ref_ptr<osg::TessellationHints> hints = new osg::TessellationHints;
    hints->setDetailRatio(2.0f);

    osg::Sphere* sph = new osg::Sphere();
    sph->setCenter(osg::Vec3(center(0),center(1),center(2)));
    sph->setRadius(radius);
    osg::ShapeDrawable* sd = new osg::ShapeDrawable(sph, hints);    
    sd->setColor(red);    
    geode->addDrawable(sd);
    root_->addChild(geode);

}

int osgViewerWrapper::show()
{    
    viewer_->setSceneData(root_.get());
    viewer_->addEventHandler(new osgViewer::StatsHandler);
    viewer_->addEventHandler(new osgViewer::WindowSizeHandler);
    viewer_->addEventHandler(
    new osgGA::StateSetManipulator(
    viewer_->getCamera()->getOrCreateStateSet()));
    viewer_->setCameraManipulator(new osgGA::TrackballManipulator());
    viewer_->realize();

    return (viewer_->run());
}