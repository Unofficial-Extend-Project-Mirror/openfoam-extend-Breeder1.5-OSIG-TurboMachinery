/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008 Franz Blaim All rights reserved
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    PrimitiveRadialPatch class 

Description
    Implements all the Sorting functions which are needed within the rotational patch

TODO   Maybe some of the functions must be made global in the following days 
\*---------------------------------------------------------------------------*/



#include "primitiveRotationalPatch.H"
namespace Foam

{

List< primitiveRotationalPatch::Corner > primitiveRotationalPatch::mergesort(const List< primitiveRotationalPatch::Corner >& a)
{    
    if (a.size() == 1)
        return a;
    
    List< Corner> l1(a);    
    l1.clear(); 
    
    label n = a.size()-a.size()/2;
    l1.setSize(n);    
    for (label i = 0; i < n; i++)
    {
        l1[i] = a[i];      
    }
    
    List< Corner> l2(a);
    l2.clear();
        
    l2.setSize(a.size()/2);    
    for (label i = n; i < a.size(); i++)
    {        
        l2[i-n] = a[i];
    }  
    l1 = mergesort(l1);
    l2 = mergesort(l2);
  
    return merge(l1, l2);
}

List< primitiveRotationalPatch::Corner> primitiveRotationalPatch::merge(List< Corner>& a, List< Corner>& b)
{
    List< Corner> c(a);
    c.clear();
    c.setSize(a.size()+b.size());
  
    int i = 0;
    while (a.size() > 0 && b.size() > 0)
    {
        if (a[0].Angle > b[0].Angle)
        {
            c[i++] = b[0];
            removeFirst(b);           
        } 
        else
        {
            c[i++] = a[0];
            removeFirst(a);
            
        }        
    }   
  
    while (a.size() != 0)
    {     
        c[i++] = a[0];
        removeFirst(a);
    }
    
    while (b.size() != 0)
    {
        c[i++] = b[0];
        removeFirst(b);
    }
    return c;
}

void primitiveRotationalPatch::removeFirst(List< Corner>& rT)
{
     List< Corner> Temp(rT);
     Temp.clear();
     Temp.setSize(rT.size()-1);
         
     forAll(rT,iT)
     {         
         if (iT != 0)
             Temp[iT-1]= rT[iT];    
     }
     rT = Temp;     
}         
         
List< scalar > primitiveRotationalPatch::mergesort(const List< scalar >& a)
{    
    if (a.size() == 1)
        return a;
    
    List< scalar > l1(a);    
    l1.clear(); 
    
    label n = a.size()-a.size()/2;
    l1.setSize(n);    
    for (label i = 0; i < n; i++)
    {
        l1[i] = a[i];      
    }
    
    List< scalar > l2(a);
    l2.clear();
        
    l2.setSize(a.size()/2);    
    for (label i = n; i < a.size(); i++)
    {        
        l2[i-n] = a[i];
    }  
    l1 = mergesort(l1);
    l2 = mergesort(l2);
  
    return merge(l1, l2);
}

List< scalar > primitiveRotationalPatch::merge(List< scalar >& a, List< scalar >& b)
{
    List< scalar > c(a);
    c.clear();
    c.setSize(a.size() + b.size());
  
    int i = 0;
    while (a.size() > 0 && b.size() > 0)
    {
        if (a[0] > b[0])
        {
            c[i++] = b[0];
            removeFirst(b);           
        } 
        else
        {
            c[i++] = a[0];
            removeFirst(a);
            
        }        
    }   
  
    while (a.size() != 0)
    {     
        c[i++] = a[0];
        removeFirst(a);
    }
    
    while (b.size() != 0)
    {
        c[i++] = b[0];
        removeFirst(b);
    }
    return c;
}

void primitiveRotationalPatch::removeFirst(List< scalar>& rT)
{
     List<scalar> Temp(rT);
     Temp.clear();
     Temp.setSize(rT.size()-1);
         
     forAll(rT,iT)
     {         
         if (iT != 0)
             Temp[iT-1]= rT[iT];    
     }
     rT = Temp;     
} 



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
         
} // End namespace Foam
        
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
