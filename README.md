# GENIE Reweight

The GENIE Reweight product is a selection of tools to propagate model uncertainties and to support generator-related 
analysis tasks. Users should note the inherent limitations of the reweighting procedure and be aware that this product
does not include weight calculators for several important systematics, and that important modelling aspects are not
reweightable in principle. The GENIE Reweight product does *not* provide the full systematic error for any GENIE 
comprehensive model or tune and, indeed, the **GENIE tuning procedure makes no use of the ReWeight product** 
but it relies on weighting functions from brute-force parameter scans made with the aid of the Professor tool [https://professor.hepforge.org/].

The GENIE Collaboration has medium-term plans to release Professor/YODA Generator response functions for all important 
modelling uncertainties, as well as to release all covariance matrices from the GENIE global fits to scattering data.
In the mean time, please be aware of the numerous caveats with the use of the Reweight product.

For more information, visit http://www.genie-mc.org

<pre>
                                                   .oooooo.    oooooooooooo ooooo      ooo ooooo oooooooooooo  
                                                  d8P'  `Y8b   `888'     `8 `888b.     `8' `888' `888'     `8  
                                                 888            888          8 `88b.    8   888   888           
                        Ndyooym          dN      888            888oooo8     8   `88b.  8   888   888oooo8     
                     Nds//+sdmoy       d+m       888     ooooo  888    "     8     `88b.8   888   888    "      
                   Nh+//ohN  m+s      N//syyN    `88.    .88'   888       o  8       `888   888   888       o  
                 Ny+//od   Nh+oN       o///+      `Y8bood8P'   o888ooooood8 o8o        `8  o888o o888ooooood8  
               Nh+//om   Nh+/yN       o///s                                                                    
              d+//+d   my+/smmyhN    m///h                                                           REWEIGHT    
            Ns///yN NdyoshNNs///d   h////yN                                                                    
           mo//om        ms///+m   d///////oyhmN                                                               
          N+//yN       ms////+N    h////////////oym                                                         
          s//h       ho/+///sN     N///////////////od                                                          
         N+/h     my++yh+//y       s/////////////////oN                                                        
         Ns/m Nmy++ymNh//+d        s//////////////////+m                                                       
          NhssoshN  Ny//sN         m///////////////////+                                                       
                  Nmo/+ohdmN        mo//////////////////h            NmmN                                      
              ms/+s//o/-----:+sdN     mhso+ooys/////////o      mhs+/------/ohN                                 
    Nhd     N+-/o+/o+------------/shN     Ndy+//////////+ mhs+:---------------om                               
  mo/sN    m:/oo/oo:-----------------:://////////////////-----------------------y                              
  y//d    Noo++o+:----------:------------:://///////////:-------:/+osyso/:-------h                             
  Nyo+syysooo+/--------------::--------------:://///////:----::////+ossyso/------:                             
     NNNNNs--------------------:/::-------------:://///:-:://////////oosyo+:------                             
          y---------------------:////:::-----------:////////////////+o++++/:-----/                             
          m-----------------------://////////////////////////////+so/--:---------y                             
           +------------------------://////////////////////////+yy:---+oh/--:y+-/N                             
           N/-------------------------://////////////////////ohy/-------yy--os-/m                              
            No--------------------------://///////////////+shs/---------/d++/-oN                               
              mo--------------------------:////////////+shyo:------------sy/sm                                 
                Ny+:------------------------::://///oyhyo:--------------/sd                                    
                    mhso+/:------------------:/+oyhyo/----------:/+syhm                                        
                           NmddhyyyssssyyyhdmmNNNmhhhyyyyhhddmN                                                
</pre>

## Authors

<pre>
Luis Alvarez-Ruso [8], Costas Andreopoulos (*) [4,6], Christopher Barry [4], Francis Bench [4], Steve Dennis [4], 
 Steve Dytman [5], Hugh Gallagher [7], Steven Gardiner[3], Walter Giele [3], Robert Hatcher [3], Libo Jiang [5], 
      Rhiannon Jones [4], Igor Kakorin [2], Konstantin Kuzmin [2], Anselmo Meregaglia [1], Donna Naples [5], 
       Vadim Naumov [2], Gabriel Perdue [3], Marco Roda [4], Vladyslav Syrotenko [7], Julia Tena Vidal [4], 
                                    Jeremy Wolcott [7], and Julia Yarba [3]

                                           (The GENIE Collaboration)

(1) CENBG, Université de Bordeaux, CNRS/IN2P3, 33175 Gradignan, France
(2) Joint Institute for Nuclear Research (JINR), Dubna, Moscow region, 141980, Russia
(3) Fermi National Accelerator Laboratory, Batavia, Illinois 60510, USA
(4) University of Liverpool, Dept. of Physics, Liverpool L69 7ZE, UK 
(5) University of Pittsburgh, Dept. of Physics and Astronomy, Pittsburgh PA 15260, USA
(6) STFC Rutherford Appleton Laboratory, Particle Physics Dept., Oxfordshire OX11 0QX, UK
(7) Tufts University, Dept. of Physics and Astronomy, Medford MA 02155, USA
(8) University of Valencia, Valencia, Spain

--------------------
(*) Corresponding Author:

 Dr. Costas Andreopoulos < constantinos.andreopoulos \at cern.ch >
    
 University of Liverpool          |  U.K. Research & Innovation (UKRI)
 Faculty of Science & Engineering |  Science & Technology Facilities Council (STFC)
 School of Physical Sciences      |  Rutherford Appleton Laboratory 
 Department of Physics            |  Particle Physics Department
 Oliver Lodge Lab 316             |  Harwell Oxford Campus, R1 2.89
 Liverpool L69 7ZE, UK            |  Oxfordshire OX11 0QX, UK          
 tel: +44-(0)1517-943201          |  tel: +44-(0)1235-445091 
</pre>
 

## Copyright

Copyright (c) 2003-2018, The GENIE Collaboration. For information, visit http://copyright.genie-mc.org 


## Physics & User manual

For installation and usage information, as well as information on the GENIE framework, event generator modules and tuning, 
see the GENIE Physics & User Manual in the public section of the GENIE Document Database:
https://genie-docdb.pp.rl.ac.uk/cgi-bin/ShowDocument?docid=2


## Public releases and physics tunes

For a list of public releases and a summary information, see http://releases.genie-mc.org

For a list of model configurations and tunes supported in each release, see http://tunes.genie-mc.org

The naming conventions for releases, model configurations and tunes are outlined in the GENIE web page and in the Physics & User manual.


## Contribution guidelines

GENIE welcomes community contributions through its Incubator. An Incubator Project is the unique route for any physics or software development into any of the GENIE suite products (Generator, Comparisons, Tuning). Details on the Incubator Project Phases (Launch, Research and Development, Graduation, Integration and Deployment) can be found in the GENIE Bylaws in the public section of the GENIE Document Database: https://genie-docdb.pp.rl.ac.uk/cgi-bin/ShowDocument?docid=1


## Citing GENIE

If you use GENIE, please **always** cite the following reference: 

<pre>
@article{Andreopoulos:2009rq,
      author         = "Andreopoulos, C. and others",
      title          = "{The GENIE Neutrino Monte Carlo Generator}",
      journal        = "Nucl. Instrum. Meth.",
      volume         = "A614",
      year           = "2010",
      pages          = "87-104",
      doi            = "10.1016/j.nima.2009.12.009",
      eprint         = "0905.2517",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "FERMILAB-PUB-09-418-CD",
      SLACcitation   = "%%CITATION = ARXIV:0905.2517;%%"
}
</pre>

If you used any of the new model configurations and tunes provided in the GENIE v3* series, please **add the following reference**:
<pre>
</pre>

Finally, if you used any of the standard GENIE applications, built-in flux and geometry drivers, or if you used any of its event reweightng and error propagation tools, please **add the following reference**:
<pre>
@article{Andreopoulos:2015wxa,
      author         = "Andreopoulos, Costas and Barry, Christopher and Dytman,
                        Steve and Gallagher, Hugh and Golan, Tomasz and Hatcher,
                        Robert and Perdue, Gabriel and Yarba, Julia",
      title          = "{The GENIE Neutrino Monte Carlo Generator: Physics and
                        User Manual}",
      year           = "2015",
      eprint         = "1510.05494",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "FERMILAB-FN-1004-CD",
      SLACcitation   = "%%CITATION = ARXIV:1510.05494;%%"
}
</pre>

Please notice that the GENIE authors endorse the **MCNET guidelines for fair academic use** which can be found in http://www.montecarlonet.org/GUIDELINES. We invite users to consider which GENIE components are important for a particular analysis and cite them, in addition to the main references. A list of such references in maintained in the official GENIE web page.
