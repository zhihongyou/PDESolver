#ifndef FIELDFUNCTION_CU
#define FIELDFUNCTION_CU

#include "finiteDifferenceCentralO2Isotropic2D.h"
#include "finiteDifferenceCentralO4Isotropic2D.h"
#include "fieldFunctionOnsite.h"
// #include "userDefinedFunction.h"


// ---------------------------------------------------------------------
void initFFuncMapAll() {
    f_func_map_all=
        {
            // Central difference O2 Isotropic in 2D            
            {{"d1x","CentralDifferenceO2Iso2D"},d1xCO2I2D},
            {{"d1y","CentralDifferenceO2Iso2D"},d1yCO2I2D},
            {{"d2x","CentralDifferenceO2Iso2D"},d2xCO2I2D},
            {{"d2y","CentralDifferenceO2Iso2D"},d2yCO2I2D},
            {{"d1x1y","CentralDifferenceO2Iso2D"},d1x1yCO2I2D},
            {{"d2x2y","CentralDifferenceO2Iso2D"},d2x2yCO2I2D},
            {{"laplace","CentralDifferenceO2Iso2D"},laplaceCO2I2D},
            {{"biLaplace","CentralDifferenceO2Iso2D"},biLaplaceCO2I2D},
            
            // Central difference O4 Isotropic in 2D            
            {{"d1x","CentralDifferenceO4Iso2D"},d1xCO4I2D},
            {{"d1y","CentralDifferenceO4Iso2D"},d1yCO4I2D},
            {{"d2x","CentralDifferenceO4Iso2D"},d2xCO4I2D},
            {{"d2y","CentralDifferenceO4Iso2D"},d2yCO4I2D},
            {{"d1x1y","CentralDifferenceO4Iso2D"},d1x1yCO4I2D},
            {{"d2x2y","CentralDifferenceO4Iso2D"},d2x2yCO4I2D},
            {{"laplace","CentralDifferenceO4Iso2D"},laplaceCO4I2D},
            {{"biLaplace","CentralDifferenceO4Iso2D"},biLaplaceCO4I2D},
            
            // Onsite functions            
            {{"1",""},oneF},
            {{"1/f",""},oneOverF},
            {{"sin",""},sinF},
            {{"cos",""},cosF},
            {{"exp",""},expF},
            {{"abs",""},absF},
            {{"sign",""},signF},
            {{"pow",""},powF},
        };
    
    // One can also use to following to add new functions
    //  to the list:
    // f_func_map_all[{"d1x","CentralDifferenceO2Iso2D"}]=d1xCO2I2D;
};


// ----------------------------------------------------------------------
void setFFuncMapDev () {
    f_func_map_all_dev[{"d1x","CentralDifferenceO2Iso2D"}]=getFFuncDevPtr(&d1xCO2I2D_dev);
    f_func_map_all_dev[{"d1y","CentralDifferenceO2Iso2D"}]=getFFuncDevPtr(&d1yCO2I2D_dev);
    f_func_map_all_dev[{"d2x","CentralDifferenceO2Iso2D"}]=getFFuncDevPtr(&d2xCO2I2D_dev);
    f_func_map_all_dev[{"d2y","CentralDifferenceO2Iso2D"}]=getFFuncDevPtr(&d2yCO2I2D_dev);
    f_func_map_all_dev[{"d1x1y","CentralDifferenceO2Iso2D"}]=getFFuncDevPtr(&d1x1yCO2I2D_dev);
    f_func_map_all_dev[{"d2x2y","CentralDifferenceO2Iso2D"}]=getFFuncDevPtr(&d2x2yCO2I2D_dev);
    f_func_map_all_dev[{"laplace","CentralDifferenceO2Iso2D"}]=getFFuncDevPtr(&laplaceCO2I2D_dev);
    f_func_map_all_dev[{"biLaplace","CentralDifferenceO2Iso2D"}]=getFFuncDevPtr(&biLaplaceCO2I2D_dev);

    f_func_map_all_dev[{"d1x","CentralDifferenceO4Iso2D"}]=getFFuncDevPtr(&d1xCO4I2D_dev);
    f_func_map_all_dev[{"d1y","CentralDifferenceO4Iso2D"}]=getFFuncDevPtr(&d1yCO4I2D_dev);
    f_func_map_all_dev[{"d2x","CentralDifferenceO4Iso2D"}]=getFFuncDevPtr(&d2xCO4I2D_dev);
    f_func_map_all_dev[{"d2y","CentralDifferenceO4Iso2D"}]=getFFuncDevPtr(&d2yCO4I2D_dev);
    f_func_map_all_dev[{"d1x1y","CentralDifferenceO4Iso2D"}]=getFFuncDevPtr(&d1x1yCO4I2D_dev);
    f_func_map_all_dev[{"d2x2y","CentralDifferenceO4Iso2D"}]=getFFuncDevPtr(&d2x2yCO4I2D_dev);
    f_func_map_all_dev[{"laplace","CentralDifferenceO4Iso2D"}]=getFFuncDevPtr(&laplaceCO4I2D_dev);
    f_func_map_all_dev[{"biLaplace","CentralDifferenceO4Iso2D"}]=getFFuncDevPtr(&biLaplaceCO4I2D_dev);

    f_func_map_all_dev[{"1",""}]=getFFuncDevPtr(&oneF_dev);
    f_func_map_all_dev[{"1/f",""}]=getFFuncDevPtr(&oneOverF_dev);
    f_func_map_all_dev[{"sin",""}]=getFFuncDevPtr(&sinF_dev);
    f_func_map_all_dev[{"cos",""}]=getFFuncDevPtr(&cosF_dev);
    f_func_map_all_dev[{"exp",""}]=getFFuncDevPtr(&expF_dev);
    f_func_map_all_dev[{"abs",""}]=getFFuncDevPtr(&absF_dev);
    f_func_map_all_dev[{"sign",""}]=getFFuncDevPtr(&signF_dev);
    f_func_map_all_dev[{"pow",""}]=getFFuncDevPtr(&powF_dev);
};


// ----------------------------------------------------------------------
void setFFuncMap () {
    initFFuncMapAll();
    f_func_map_all_dev=f_func_map_all;
    setFFuncMapDev();
    addUserDefinedFuncs();
};


#endif
