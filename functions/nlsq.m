function [x,CostFunction,JAC,EXITFLAG,OUTPUT,msg] = nlsq(funfcn,x,verbosity,options,defaultopt,CostFunction,JAC,XDATA,YDATA,caller,varargin)
%NLSQ Helper function that solves non-linear least squares problems.
%   NLSQ is the core code for solving problems of the form:
%   min sum {FUN(X).^2} where FUN and X may be vectors or matrices.   
%             x
%

%   Copyright 1990-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2007/05/23 18:59:12 $

% ------------Initialization----------------
XOUT = x(:);
% numberOfVariables must be the name of this variable
numberOfVariables = length(XOUT);
msg = [];
how = [];
OUTPUT = [];
iter = 0;
EXITFLAG = 1;  %assume convergence
currstepsize = 0;

% Handle the output
if isfield(options,'OutputFcn')
    outputfcn = optimget(options,'OutputFcn',defaultopt,'fast');
else
    outputfcn = defaultopt.OutputFcn;
end
if isempty(outputfcn)
    haveoutputfcn = false;
else
    haveoutputfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse OutputFcn which is needed to support cell array syntax for OutputFcn.
    outputfcn = createCellArrayOfFunctions(outputfcn,'OutputFcn');
end
stop = false;

% Handle the output
if isfield(options,'PlotFcns')
    plotfcns = optimget(options,'PlotFcns',defaultopt,'fast');
else
    plotfcns = defaultopt.PlotFcns;
end
if isempty(plotfcns)
    haveplotfcn = false;
else
    haveplotfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse PlotFcns which is needed to support cell array syntax for PlotFcns.
    plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
end

formatstrFirstIter = ' %5.0f       %5.0f   %13.6g';
formatstr = ' %5.0f       %5.0f   %13.6g %12.3g    %12.3g';

% If lsqcurvefit, then add XDATA to objective's input list.
% xargin{1} will be updated at each iteration
if ~isempty(XDATA)
    xargin = {x,XDATA};
else
    xargin = {x};
end

% options
gradflag =  strcmp(optimget(options,'Jacobian',defaultopt,'fast'),'on');
tolX = optimget(options,'TolX',defaultopt,'fast');
lineSearchType = strcmp(optimget(options,'LineSearchType',defaultopt,'fast'),'cubicpoly');
% If caller is fsolve, LevenbergMarquardt is not a member of the defaultopt structure,
% but options.LevenbergMarquardt will exist and will have been set to either 'on' or 'off'
% in fsolve. Thus, the call to optimget below will work as expected. 
levMarq = strcmp(optimget(options,'LevenbergMarquardt',defaultopt,'fast'),'on');
tolFun = optimget(options,'TolFun',defaultopt,'fast');
DiffMinChange = optimget(options,'DiffMinChange',defaultopt,'fast');
DiffMaxChange = optimget(options,'DiffMaxChange',defaultopt,'fast');
DerivativeCheck = strcmp(optimget(options,'DerivativeCheck',defaultopt,'fast'),'on');
typicalx = optimget(options,'TypicalX',defaultopt,'fast') ;
if ischar(typicalx)
   if isequal(lower(typicalx),'ones(numberofvariables,1)')
      typicalx = ones(numberOfVariables,1);
   else
      error('optim:nlsq:InvalidTypicalX','Option ''TypicalX'' must be a numerical value if not the default.')
   end
end
maxFunEvals = optimget(options,'MaxFunEvals',defaultopt,'fast');
maxIter = optimget(options,'MaxIter',defaultopt,'fast');
if ischar(maxFunEvals)
   if isequal(lower(maxFunEvals),'100*numberofvariables')
      maxFunEvals = 100*numberOfVariables;
   else
      error('optim:nlsq:MaxFunEvals','Option ''MaxFunEvals'' must be a numeric value if not the default.')
   end
end

nfun=length(CostFunction);
numFunEvals = 1;
numGradEvals = 0;
MATX=zeros(3,1);
MATL=[CostFunction'*CostFunction;0;0];
FIRSTF=CostFunction'*CostFunction;
PCNT = 0;
EstSum=0.5;
% system of equations or overdetermined
if nfun >= numberOfVariables
   GradFactor = 0;  
else % underdetermined: singularity in JAC'*JAC or GRAD*GRAD' 
   GradFactor = 1;
end
done = false;
NewF = CostFunction'*CostFunction;


while ~done  
   % Work Out Gradients
   if ~(gradflag) || DerivativeCheck
      JACFD = zeros(nfun, numberOfVariables);  % set to correct size
      JACFD = finitedifferences(XOUT,x,funfcn,[],[],[],CostFunction,[],XDATA,YDATA,...
                              DiffMinChange,DiffMaxChange,typicalx,[],'all',[],[],[],[],[],[],...
                              false,varargin{:});       
      numFunEvals=numFunEvals+numberOfVariables;
      % In the particular case when there is only one function in more
      % than one variable, finitedifferences() returns JACFD in a column
      % vector, whereas nlsq() expects a single-row Jacobian:
      if (nfun == 1) && (numberOfVariables > 1)
         JACFD = JACFD';
      end
      % Gradient check
      if DerivativeCheck && gradflag         
         if isa(funfcn{3},'inline') 
            % if using inlines, the gradient is in funfcn{4}
            graderr(JACFD, JAC, formula(funfcn{4})); %
         else 
            % otherwise fun/grad in funfcn{3}
            graderr(JACFD, JAC,  funfcn{3});
         end
         DerivativeCheck = 0;
      else
         JAC = JACFD;
      end
   else
      x(:) = XOUT;
   end
   GradF = 2*(CostFunction'*JAC)'; %2*GRAD*CostFunction;
   
   %---------------Initialization of Search Direction------------------
   JacTJac = JAC'*JAC;
   if iter == 0
       OLDX = XOUT;
       OLDF = CostFunction;
       displayHeaders(verbosity,levMarq);
       
       % Initialize the output function.
       if haveoutputfcn || haveplotfcn
           [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,caller,XOUT,xOutputfcn,'init',iter,numFunEvals, ...
               CostFunction,NewF,[],[],GradF,[],[],varargin{:});
           if stop
               [x,CostFunction,JAC,EXITFLAG,OUTPUT,msg] = cleanUpInterrupt(xOutputfcn,optimValues,levMarq,caller);
               return;
           end
       end
       
       if verbosity > 1
           disp(sprintf(formatstrFirstIter,iter,numFunEvals,NewF)); 
       end
       
       % 0th iteration
       if haveoutputfcn || haveplotfcn
           [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,caller,XOUT,xOutputfcn,'iter',iter,numFunEvals, ...
               CostFunction,NewF,[],[],GradF,[],[],varargin{:});
           if stop
               [x,CostFunction,JAC,EXITFLAG,OUTPUT,msg] = cleanUpInterrupt(xOutputfcn,optimValues,levMarq,caller);
               return;
           end
       end

      % Disable the warnings about conditioning for singular and
      % nearly singular matrices
      warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
      warningstate2 = warning('off', 'MATLAB:singularMatrix');
      
      if condest(JacTJac)>1e16 
         SD=-(JacTJac+(norm(full(JAC))+1)*(eye(numberOfVariables,numberOfVariables)))\(CostFunction'*JAC)';
         if levMarq 
            GradFactor=norm(full(JAC))+1; 
         end
         how='COND';
      else
         SD=-(JacTJac+GradFactor*(eye(numberOfVariables,numberOfVariables)))\(CostFunction'*JAC)';
      end
      FIRSTF=NewF;
      OLDJ = JAC;
      GDOLD=GradF'*SD;
      % currstepsize controls the initial starting step-size.
      % If currstepsize has been set externally then it will
      % be non-zero, otherwise set to 1. Right now it's not
      % possible to set it externally.
      if currstepsize == 0, 
         currstepsize=1; 
      end
      LMorSwitchFromGNtoLM = false;
      if levMarq
         newf=JAC*SD+CostFunction;
         GradFactor=newf'*newf;
         SD=-(JacTJac+GradFactor*(eye(numberOfVariables,numberOfVariables)))\(CostFunction'*JAC)'; 
         LMorSwitchFromGNtoLM = true;
      end

      % Restore the warning states to their original settings
      warning(warningstate1)
      warning(warningstate2)

      newf=JAC*SD+CostFunction;
      EstSum=newf'*newf;
      status=0;
      if lineSearchType==0; 
         PCNT=1; 
      end
     
  else % iter >= 1
      gdnew=GradF'*SD;
      if verbosity > 1 
          num=sprintf(formatstr,iter,numFunEvals,NewF,currstepsize,gdnew);
          if LMorSwitchFromGNtoLM
              num=[num,sprintf('   %12.6g ',GradFactor)]; 
          end
          if isinf(verbosity)
              disp([num,'       ',how])
          else
              disp(num)
          end
      end
      
      if haveoutputfcn || haveplotfcn
          [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,caller,XOUT,xOutputfcn,'iter',iter,numFunEvals, ...
              CostFunction,NewF,currstepsize,gdnew,GradF,SD,GradFactor,varargin{:});
          if stop
              [x,CostFunction,JAC,EXITFLAG,OUTPUT,msg] = cleanUpInterrupt(xOutputfcn,optimValues,levMarq,caller);
              return;
          end
      end
      
      %-------------Direction Update------------------


      if gdnew>0 && NewF>FIRSTF
         % Case 1: New function is bigger than last and gradient w.r.t. SD -ve
         % ... interpolate. 
         how='inter';
         [stepsize]=cubici1(NewF,FIRSTF,gdnew,GDOLD,currstepsize);
         currstepsize=0.9*stepsize;
      elseif NewF<FIRSTF
         %  New function less than old fun. and OK for updating 
         %         .... update and calculate new direction. 
         [newstep,fbest] =cubici3(NewF,FIRSTF,gdnew,GDOLD,currstepsize);
         if fbest>NewF,
            fbest=0.9*NewF; 
         end 
         if gdnew<0
            how='incstep';
            if newstep<currstepsize,  
               newstep=(2*currstepsize+1e-4); how=[how,'IF']; 
            end
            currstepsize=abs(newstep);
         else 
            if currstepsize>0.9
               how='int_step';
               currstepsize=min([1,abs(newstep)]);
            end
         end
         % SET DIRECTION.      
         % Gauss-Newton Method    
         LMorSwitchFromGNtoLM = true ;

         % Disable the warnings about conditioning for singular and
         % nearly singular matrices
         warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
         warningstate2 = warning('off', 'MATLAB:singularMatrix');

         if ~levMarq
            if currstepsize>1e-8 && condest(JacTJac)<1e16
               SD=JAC\(JAC*XOUT-CostFunction)-XOUT;
               if SD'*GradF>eps,
                  how='ERROR- GN not descent direction';
               end
               LMorSwitchFromGNtoLM = false;
            else
               if verbosity > 1
                   if currstepsize <= 1e-8
                       fprintf('Step size too small - Switching to LM method.\n');
                   else % condest(JacTJac) >= 1e16
                       fprintf('Iteration matrix ill-conditioned - Switching to LM method.\n');
                   end
               end
               how='CHG2LM';
               levMarq = true;
               currstepsize=abs(currstepsize);               
            end
         end
         
         if (LMorSwitchFromGNtoLM)      
            % Levenberg_marquardt Method N.B. EstSum is the estimated sum of squares.
            %                                 GradFactor is the value of lambda.
            % Estimated Residual:
            if EstSum>fbest
               GradFactor=GradFactor/(1+currstepsize);
            else
               GradFactor=GradFactor+(fbest-EstSum)/(currstepsize+eps);
            end
            SD=-(JacTJac+GradFactor*(eye(numberOfVariables,numberOfVariables)))\(CostFunction'*JAC)'; 
            currstepsize=1; 
            estf=JAC*SD+CostFunction;
            EstSum=estf'*estf;
         end
         gdnew=GradF'*SD;

         % Restore the warning states to their original settings
         warning(warningstate1)
         warning(warningstate2)
         
         OLDX=XOUT;
         % Save Variables
         FIRSTF=NewF;
         OLDG=GradF;
         GDOLD=gdnew;    
         
         % If quadratic interpolation set PCNT
         if lineSearchType==0, 
            PCNT=1; MATX=zeros(3,1); MATL(1)=NewF; 
         end
      else 
         % Halve Step-length
         how='Red_Step';
         if NewF==FIRSTF,
            msg = sprintf('No improvement in search direction: Terminating.');
            done = true;
            EXITFLAG = -4;
         else
            currstepsize=currstepsize/8;
            if currstepsize<1e-8
               currstepsize=-currstepsize;
            end
         end
      end
   end % if iter==0 
   
   %----------End of Direction Update-------------------
   iter = iter + 1;

   if lineSearchType==0, 
      PCNT=1; MATX=zeros(3,1);  MATL(1)=NewF; 
   end
   % Check Termination 
   if (GradF'*SD) < tolFun && ...
           max(abs(GradF)) < 10*(tolFun+tolX)
       msg = sprintf(['Optimization terminated: directional derivative along\n' ... 
                          ' search direction less than TolFun and infinity-norm of\n' ...
                          ' gradient less than 10*(TolFun+TolX).']);
       done = true; EXITFLAG = 1;
   elseif max(abs(SD))< tolX 
       msg = sprintf('Optimization terminated: magnitude of search direction less than TolX.');     
       done = true; EXITFLAG = 4;

   elseif numFunEvals > maxFunEvals
      msg = sprintf(['Maximum number of function evaluations exceeded.',...
                 ' Increase OPTIONS.MaxFunEvals.']);
      done = true;
      EXITFLAG = 0;
   elseif iter > maxIter
      msg = sprintf(['Maximum number of iterations exceeded.', ...
                        ' Increase OPTIONS.MaxIter.']);
      done = true;
      EXITFLAG = 0;
   else % continue
      XOUT=OLDX+currstepsize*SD;
      % Line search using mixed polynomial interpolation and extrapolation.
      if PCNT~=0
         % initialize OX and OLDF 
         OX = XOUT; OLDF = CostFunction;
         while PCNT > 0 && numFunEvals <= maxFunEvals
            x(:) = XOUT; 
            xargin{1} = x; % update x in list of input arguments to objective
            CostFunction = feval(funfcn{3},xargin{:},varargin{:});
            if ~isempty(YDATA)
               CostFunction = CostFunction - YDATA;
            end
            CostFunction = CostFunction(:);
            numFunEvals=numFunEvals+1;
            NewF = CostFunction'*CostFunction;
            % <= used in case when no improvement found.
            if NewF <= OLDF'*OLDF, 
               OX = XOUT; 
               OLDF=CostFunction; 
            end
            [PCNT,MATL,MATX,steplen,NewF,how]=searchq(PCNT,NewF,OLDX,MATL,MATX,SD,GDOLD,currstepsize,how);
            currstepsize=steplen;
            XOUT=OLDX+steplen*SD;
            if NewF==FIRSTF,  
               PCNT=0; 
            end
            if haveoutputfcn % Call output functions via callAllOptimOutputFcns wrapper
                stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,'interrupt',varargin{:});
                if stop
                    [x,CostFunction,JAC,EXITFLAG,OUTPUT,msg] = cleanUpInterrupt(xOutputfcn,optimValues,levMarq,caller);
                    return;
                end
            end
         end % end while
         XOUT = OX;
         CostFunction=OLDF;
         if numFunEvals>maxFunEvals 
            msg = sprintf(['Maximum number of function evaluations exceeded.',...
                            ' Increase OPTIONS.MaxFunEvals.']);
            done = true; 
            EXITFLAG = 0;
         end
      end % if PCNT~=0
   
      % Evaluate objective
      x(:)=XOUT;
      xargin{1} = x; % update x in list of input arguments to objective
      switch funfcn{1}
          case 'fun'
              CostFunction = feval(funfcn{3},xargin{:},varargin{:});
              if ~isempty(YDATA)
                  CostFunction = CostFunction - YDATA;
              end
              CostFunction = CostFunction(:);
              nfun=length(CostFunction);
              % JAC will be updated when it is finite-differenced
          case 'fungrad'
              [CostFunction,JAC] = feval(funfcn{3},xargin{:},varargin{:});
              if ~isempty(YDATA)
                  CostFunction = CostFunction - YDATA;
              end
              CostFunction = CostFunction(:);
              numGradEvals=numGradEvals+1;
          case 'fun_then_grad'
              CostFunction = feval(funfcn{3},xargin{:},varargin{:});
              if ~isempty(YDATA)
                  CostFunction = CostFunction - YDATA;
              end
              CostFunction = CostFunction(:);
              JAC = feval(funfcn{4},xargin{:},varargin{:});
              numGradEvals=numGradEvals+1;
          otherwise
              error('optim:nlsq:InvalidCalltype','Undefined calltype in LSQNONLIN.')
      end
      numFunEvals=numFunEvals+1;
      NewF = CostFunction'*CostFunction;

   end % Convergence testing
end  % while

NewF = CostFunction'*CostFunction;
gdnew=GradF'*SD;

if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,caller,XOUT,xOutputfcn,'iter',iter,numFunEvals, ...
        CostFunction,NewF,currstepsize,GDOLD,GradF,SD,GradFactor,varargin{:});
    if stop
        [x,CostFunction,JAC,EXITFLAG,OUTPUT,msg] = cleanUpInterrupt(xOutputfcn,optimValues,levMarq,caller);
        return;
    end
end

x(:)=XOUT;

OUTPUT.iterations = iter;
OUTPUT.funcCount = numFunEvals;
OUTPUT.stepsize=currstepsize;
OUTPUT.cgiterations = [];
OUTPUT.firstorderopt = [];

if levMarq
   OUTPUT.algorithm='medium-scale: Levenberg-Marquardt, line-search';
else
   OUTPUT.algorithm='medium-scale: Gauss-Newton, line-search';
end

if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues] = callOutputAndPlotFcns(outputfcn,plotfcns,caller,XOUT,xOutputfcn,'done',iter,numFunEvals, ...
        CostFunction,NewF,currstepsize,GDOLD,GradF,SD,GradFactor,varargin{:});
    % Optimization done, so ignore "stop"
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function displayHeaders(verbosity,levMarq)

if verbosity>1
   if ~levMarq
      if isinf(verbosity)
         header = sprintf(['\n                                                     Directional \n',...
                             ' Iteration  Func-count    Residual     Step-size      derivative   Line-search']);
      else
         header = sprintf(['\n                                                     Directional \n',...
                             ' Iteration  Func-count    Residual     Step-size      derivative ']);
      end
   else
      if isinf(verbosity)
         header = sprintf(['\n                                                     Directional \n',...
                             ' Iteration  Func-count    Residual     Step-size      derivative   Lambda       Line-search']);
      else
         header = sprintf(['\n                                                     Directional \n',...
                             ' Iteration  Func-count    Residual     Step-size      derivative    Lambda']);
      end
   end
   disp(header)
end
%--------------------------------------------------------------------------
function [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,caller,x,xOutputfcn,state,iter,numFunEvals, ...
    CostFunction,NewF,currstepsize,gdnew,GradF,SD,GradFactor,varargin)
% CALLOUTPUTANDPLOTFCNS assigns values to the struct OptimValues and then calls the
% outputfcn/plotfcns.  
%
% state - can have the values 'init','iter', or 'done'. 

% For the 'done' state we do not check the value of 'stop' because the
% optimization is already done.
optimValues.iteration = iter;
optimValues.funccount = numFunEvals;
optimValues.stepsize = currstepsize;
optimValues.directionalderivative = gdnew;
optimValues.gradient = GradF;
optimValues.firstorderopt = norm(GradF,Inf);
optimValues.searchdirection = SD;
optimValues.lambda = GradFactor;
if isequal(caller,'fsolve') 
   optimValues.fval = CostFunction; 
else % lsqnonlin, lsqcurvefit 
   optimValues.residual = CostFunction; 
   optimValues.resnorm = NewF; 
end 
xOutputfcn(:) = x;  % Set x to have user expected size
stop = false;
% Call output function
if ~isempty(outputfcn)
    switch state
        case {'iter','init'}
            stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:});
        case 'interrupt'
            % No 'interrupt' case in nlsq
        otherwise
            error('optim:nlsq:UnknownStateInCALLOUTPUTANDPLOTFCNS','Unknown state in CALLOUTPUTANDPLOTFCNS.')
    end
end
% Call plot functions
if ~isempty(plotfcns)
    switch state
        case {'iter','init'}
            stop = callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error('optim:nlsq:UnknownStateInCALLOUTPUTANDPLOTFCNS','Unknown state in CALLOUTPUTANDPLOTFCNS.')
    end
end
%--------------------------------------------------------------------------
function [x,CostFunction,JAC,EXITFLAG,OUTPUT,msg] = cleanUpInterrupt(xOutputfcn,optimValues,levMarq,caller)

x = xOutputfcn;
% CostFunction can be either 'fval' (fsolve) or 'residual'
if isequal(caller,'fsolve') 
    CostFunction = optimValues.fval;
else 
    CostFunction = optimValues.residual;
end
EXITFLAG = -1; 
OUTPUT.iterations = optimValues.iteration;
OUTPUT.funcCount = optimValues.funccount;
OUTPUT.stepsize = optimValues.stepsize;
OUTPUT.cgiterations = [];
OUTPUT.firstorderopt = [];

if levMarq
   OUTPUT.algorithm='medium-scale: Levenberg-Marquardt, line-search';
else
   OUTPUT.algorithm='medium-scale: Gauss-Newton, line-search';
end

JAC = []; % May be in an inconsistent state
msg = 'Optimization terminated prematurely by user.';
