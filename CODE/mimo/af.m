## Copyright (C) 2021 pudja
## 
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see
## <https://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {} {@var{retval} =} array_factor (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: pudja <pudja@almanshurin>
## Created: 2021-01-07

function retval = af (data, range, point)
step = range/point;

for ph0=0:step:(range)-step
  i = round(ph0/step)+1;
    
  tmp = 0;
  for n=1:length(data)
    amplitude = abs(data(n));
    phase = angle(data(n));

    if n==1
      phase_normalized = 0;
    else 
      phase_normalized = phase - ref_phase;
    endif
    ref_phase = phase
    
##    amplitude = 1;
##    phase_normalized = 150*(pi/180);

    a = -pi*cos(phase_normalized);
    ps0 = pi*cos(ph0) + a;
    
    tmp += amplitude*exp(j*ps0*n);
  endfor
  AF(i) = abs(tmp);
endfor

retval = AF;

endfunction
