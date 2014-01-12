#ifndef RIGIDPATCHSIMILARITY_H
#define RIGIDPATCHSIMILARITY_H

#include "DenseVector.h"
#include "DenseMatrix.h"

#include "EuclideanMetric.h"

#include "SVD.h"
#include "RandomSVD.h"

#include <algorithm>
#include <map>

#include "GWT.h"
#include "Wasserstein.h"

#include "itkCastImageFilter.h"
#include "itkEllipseSpatialObject.h"
#include "itkImage.h"
#include "itkImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkCenteredRigid2DTransform.h"
#include <itkNormalizedCorrelationImageToImageMetric.h>
#include <itkConjugateGradientOptimizer.h>
#include <itkExhaustiveOptimizer.h>



class CommandIterationUpdateExhaust : public itk::Command 
{
  public:
    typedef  CommandIterationUpdateExhaust   Self;
    typedef  itk::Command             Superclass;
    typedef itk::SmartPointer<Self>  Pointer;
    typedef itk::Rigid2DTransform< double > TransformType;

    itkNewMacro( Self );

    double bestValue;
    TransformType::Pointer transform;

  protected:
    CommandIterationUpdateExhaust() 
    {
      bestValue = std::numeric_limits<double>::max();
    };

  public:

    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    };

    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      typedef itk::ExhaustiveOptimizer     OptimizerType;
      typedef const OptimizerType * OptimizerPointer;
      OptimizerPointer optimizer = 
        dynamic_cast< OptimizerPointer >( object );
     // if( itk::IterationEvent().CheckEvent( &event ) ){
        if(bestValue > optimizer->GetCurrentValue() ){
          bestValue = optimizer->GetCurrentValue();
          //std::cout << transform->GetParameters() << std::endl;
        }
      //}
    };

    void SetTransform(TransformType::Pointer t){
      transform = t;
    };

};


template <typename TPrecision>
class RigidPatchSimilarity : public SimilarityVisitor<TPrecision>{
private:

  std::map<int, GWTInfo<TPrecision> * > nodes;
  int maxID;
  EuclideanMetric<TPrecision> l2;


  typedef typename itk::Image<TPrecision, 2> ImageType;
  typedef typename ImageType::Pointer ImagePointer;

public:
  
   RigidPatchSimilarity<TPrecision>(){
     maxID = 0;
   };

   void visit(GMRANode<TPrecision> *node){
    GWTInfo<TPrecision> *gwt = (GWTInfo<TPrecision>*)node->getNodeInfo();
    
    int ID = gwt->ID;
    if(ID > maxID){
      maxID = ID;
    }
    nodes[ID] = gwt;  
  };

  std::map<int, GWTInfo<TPrecision> * > &getNodes(){
    return nodes;
  };

  FortranLinalg::DenseMatrix<TPrecision> distances(){
    using namespace FortranLinalg;

    DenseMatrix<TPrecision> D(maxID+1, maxID+1);
    for(int i=0; i< D.M(); i++){
       D(i, i) = 0;
       GWTInfo<TPrecision> *v1 = nodes[i];
       for(int j=i+1; j<D.N(); j++){
         GWTInfo<TPrecision> *v2 = nodes[j];
         D(i, j) = 1+distance(v1, v2);
         D(j, i) = D(i, j);
      }
       std::cout << i << std::endl;
    }
    return D;
  };


  TPrecision distance(GWTInfo<TPrecision> *n1, GWTInfo<TPrecision> *n2){
    using namespace FortranLinalg;
    ImagePointer m1 = getImage(n1->center);
    ImagePointer m2 = getImage(n2->center);

    //RigidRegistration
    //  The transform that will map the fixed image into the moving image.
    //typedef itk::TranslationTransform< double, 2 > TransformType;
    //typedef itk::CenteredRigid2DTransform< double > TransformType;
    typedef itk::Rigid2DTransform< double > TransformType;

    //  An optimizer is required to explore the parameter space of the transform
    //  in search of optimal values of the metric.
    //typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;
    typedef itk::ExhaustiveOptimizer     OptimizerType;
    //typedef itk::ConjugateGradientOptimizer OptimizerType;


    //  The metric will compare how well the two images match each other. Metric
    //  types are usually parameterized by the image types as it can be seen in
    //  the following type declaration.
   //typedef itk::MeanSquaresImageToImageMetric< ImageType, ImageType >    MetricType;
    typedef itk::NormalizedCorrelationImageToImageMetric< ImageType, ImageType >    MetricType;

    //  Finally, the type of the interpolator is declared. The interpolator will
    //  evaluate the intensities of the moving image at non-grid positions.
    typedef itk:: LinearInterpolateImageFunction<
      ImageType,
      double          >    InterpolatorType;

    //  The registration method type is instantiated using the types of the
    //  fixed and moving images. This class is responsible for interconnecting
    //  all the components that we have described so far.
    typedef itk::ImageRegistrationMethod<
      ImageType,
      ImageType >    RegistrationType;

    // Create components
    typename MetricType::Pointer         metric        = MetricType::New();
    typename TransformType::Pointer      transform     = TransformType::New();
    typename OptimizerType::Pointer      optimizer     = OptimizerType::New();
    typename InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
    typename RegistrationType::Pointer   registration  = RegistrationType::New();

    // Each component is now connected to the instance of the registration method.
    registration->SetMetric(        metric        );
    registration->SetOptimizer(     optimizer     );
    registration->SetTransform(     transform     );
    registration->SetInterpolator(  interpolator  ); 

    typename TransformType::CenterType center;
    center.Fill(sqrt(n1->center.N())/2.0);
    transform->SetCenter(center);

    //  Initialize the transform
    typedef typename RegistrationType::ParametersType ParametersType;
    ParametersType initialParameters( transform->GetNumberOfParameters() );

    initialParameters.Fill(0.0);  

    registration->SetInitialTransformParameters( initialParameters );


    /*
    optimizer->SetMaximumStepLength( 1.0 );
    optimizer->SetMinimumStepLength( 0.1 );
    // Set a stopping criterion
    optimizer->SetNumberOfIterations( 50 );
    */

    typename OptimizerType::StepsType steps( transform->GetNumberOfParameters() );
    steps[0] = 8;
    steps[1] = 4;
    steps[2] = 4;
    optimizer->SetNumberOfSteps( steps );
    optimizer->SetStepLength( 1 );
    typename OptimizerType::ScalesType scale( transform->GetNumberOfParameters() );
    scale[0] = 3.1416/steps[0];
    scale[1] = 0.5;
    scale[2] = 0.5;
    optimizer->SetScales( scale );

    CommandIterationUpdateExhaust::Pointer observer = 
      CommandIterationUpdateExhaust::New();
    observer->SetTransform(transform);
    optimizer->AddObserver( itk::IterationEvent(), observer );


    // Set the registration inputs
    registration->SetFixedImage(m1);
    registration->SetMovingImage(m2);

    registration->SetFixedImageRegion(
        m1->GetLargestPossibleRegion() );

    
    //std::cout << optimizer->GetValue() << " / ";
    try
    {
      registration->Update();
    }
    catch( itk::ExceptionObject & err )
    {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      //std::cout << optimizer->GetValue() << " *" << std::endl;
      //return optimizer->GetValue();
      std::cout << observer->bestValue << " *" << std::endl;
      return observer->bestValue;
      //return EXIT_FAILURE;
    }
    //std::cout << optimizer->GetValue() << std::endl;
    //return optimizer->GetValue();
    return observer->bestValue;
  };




private:


   static ImagePointer getImage(FortranLinalg::DenseVector<TPrecision> &v){
     int n = sqrt(v.N());
     // Create a black image with 2 white regions
     ImagePointer image = ImageType::New();
     typename ImageType::IndexType start;
     start.Fill(0);

     typename ImageType::SizeType size;
     size.Fill(n);

     typename ImageType::RegionType region(start,size);
     image->SetRegions(region);
     image->Allocate();
     image->FillBuffer(0);

     // Make a square
     int index =0;
     for(unsigned int r = 0; r < n; r++)
     {
       for(unsigned int c = 0; c < n; c++)
       {
         typename ImageType::IndexType pixelIndex;
         pixelIndex[0] = r;
         pixelIndex[1] = c;

         image->SetPixel(pixelIndex, v(index) );
         index++;

       }
     }
     return image;
   };



};




#endif
