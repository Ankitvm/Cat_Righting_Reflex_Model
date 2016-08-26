function msfcn_newtonmethod(block)
% Level-2 M-file S-function to find beta from differential equation using
% Newton-Raphson method

% Author: Ankit Manerikar, Tanvi Anandpara

% This S-function receives as input the current values of alpha and theta
% and solves the differential equation for the Kane-Scher model of
% cat-righting reflex using numerical techniques (Newton-Raphson method) to
% compute the value of beta for the model.

setup(block);

%endfunction


function setup(block)

% Register number of ports
block.NumInputPorts  = 2;
block.NumOutputPorts = 1;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
  block.InputPort(1).Complexity   = 'Real'; 
  block.InputPort(1).DataTypeId   = 0;
  block.InputPort(1).SamplingMode = 'Sample';
  block.InputPort(1).Dimensions   = 1;
  
  block.InputPort(2).Complexity   = 'Real'; 
  block.InputPort(2).DataTypeId   = 0;
  block.InputPort(2).SamplingMode = 'Sample';
  block.InputPort(2).Dimensions   = 1;
  
  block.OutputPort(1).Complexity   = 'Real';
  block.OutputPort(1).DataTypeId   = 0;
  block.OutputPort(1).SamplingMode = 'Sample';
  block.OutputPort(1).Dimensions   = 1;

% Register sample times
block.SampleTimes = [-1 0];

% Register methods
block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);
block.RegBlockMethod('Update', @Update);
block.RegBlockMethod('Terminate', @Terminate);

%endfunction


function DoPostPropSetup(block)

% Initialize the Dwork vector
block.NumDworks = 1;

% Dwork(1) stores the initial value of beta for the numerical method
block.Dwork(1).Name            = 'b0';
block.Dwork(1).Dimensions      = 1;
block.Dwork(1).DatatypeID      = 0;      % double
block.Dwork(1).Complexity      = 'Real'; % real
block.Dwork(1).UsedAsDiscState = true;

%endfunction

function Start(block)

% Populate the Dwork vector
block.Dwork(1).Data = 0.785;

%endfunction

function InitializeConditions(block)

% Set the initial pulse width value
% set_param(block.Dwork(2).Data, 'PulseWidth', num2str(50));

%endfunction

function Outputs(block)

%finding beta using Newton-Raphson method

%initializing work parameter for finding beta 
b0 = block.Dwork(1).Data; 

%taking current input values for alpha and theta  
a = block.InputPort(1).Data;
d = block.InputPort(2).Data;

bn = b0;
b_diff = b0; 
bn0 = b0;
    
while (abs(b_diff)>= 0.000001) % putting low errro margin
        
        bn0 = bn;
        
        %values for S and T in the Kane-Scher equation
        t = cos(a)*cos(bn) - sin(a)*sin(bn)*cos(d);
        s = (-sqrt(2))*(sin(bn)*(cos(a)*sin(bn) + sin(a)*cos(bn)*cos(d)));

        %derivatives of S and T
        td = -(cos(a).*sin(bn)+ sin(a).*cos(bn).*cos(d));
        sd = -sqrt(2)*(cos(a).*sin(2*bn) + sin(a).*cos(2*bn).*cos(d));

        %simplified expression for Kane-Scher Equation by taking
        %d(phi)/d(theta) = 0.5
        fcheck = (((1+ t).^2).*((5*t+3).^2).*(1-t)) - 4*(s.^2); 
        %derivative for expression
        fcheckd = (-2.*(1-t).*(1+t).*((5*t+3).^2) + 10*((1-t).^2).*(1+t).*(5*t+3) + ((1-t).^2).*(5*t+3)).*(td) - (8.*s.*(sd));

        %formula for Newton-Raphson method
        bnew = bn -(fcheck./fcheckd);
        
        b_diff = bnew - bn0;    % check for error difference
        bn = bnew;              % update value for next iteration
        
end  

 block.OutputPort(1).Data = (3.14157 -(bnew));   %put final value as output
 block.Dwork(1).Data = (bnew);
%endfunction


function Update(block)

% Store the input value in the Dwork(1)
%block.Dwork(1).Data = block.InputPort(1).Data;

%endfunction


function Terminate(block)

%endfunction
