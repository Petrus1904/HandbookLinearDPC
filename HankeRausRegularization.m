function [x, y] = HankeRausRegularization(A, b, range, markValue)
    %HankeRausRegularization(A, b, range, markValue)
    % plots the regularization score for various values
    % of lambda regarding the Tikhonov regularization problem:
    %
    % min_g ||Ag-b||^2_2 + lambda ||g||^2_2
    %
    % arguments
    %   A, b -- problem matrix and vector
    %   range -- range of lambdas to evaluate and number of points. 
    %            Point are logaritmically spaced, but use the desired values, not
    %            the logaritmic exponents. - Default: [1e-3, 1e6, 1000];
    %   markValue -- if you want to highlight the chosen value (optional)
    %
    % Returns
    %   plot -- if nargout == 0. Plots the Hanke and Raus heuristic for the
    %           regularization problem.
    %   x,y --  the x axis values and y axis values of the plot (if you like
    %           to make a better one).

    %% Compute hanke and raus heuristic
    if(nargin < 3 || isempty(range))
        range = [1e-3, 1e6, 1000];
    end
    
    lamLog = logspace(log10(range(1)), log10(range(2)), range(3));
    fh=zeros(1, length(lamLog));
    %Hanke and Rause heuristic     
    for i = 1:length(lamLog)
        lam = lamLog(i);
        aux=inv(A*A'+lam*eye(size(A, 1)));
        flam=sqrt(abs(b'*lam^2*aux^3*b));
        fh(:,i)=flam;
    end
    
    if(nargin == 4 && ~isempty(markValue))
        aux=inv(A*A'+markValue*eye(size(A, 1)));
        markY=sqrt(abs(b'*markValue^2*aux^3*b));
    end
    
    %% plot or return

    if(nargout == 0)
        figure();
        semilogx(lamLog,fh,'r', 'LineWidth', 2, 'DisplayName', 'Hanke and Raus');
        hold on;
        if(nargin == 4 && ~isempty(markValue))
            semilogx(markValue, markY, 'w*', 'LineWidth', 2, 'HandleVisibility', 'off');
            semilogx(markValue, markY, 'ro', 'LineWidth', 2, 'HandleVisibility', 'off');
        end
%         semilogx(lamLog,sqrt(b'*b)*lamLog.^(-0.5),'b--', 'LineWidth', 2, 'DisplayName', '$lim_{inf}$');
%         if(rank(A) < size(A,1))
%             semilogx(lamLog,sqrt(b'*(eye(size(A,1))-A*pinv(A))*b)*lamLog.^(-0.5),'b--', 'LineWidth', 2, 'DisplayName', '$lim_0$');
%         else
%             semilogx(lamLog,sqrt(b'*inv(A*A')^3*b)*lamLog,'b--', 'LineWidth', 2, 'DisplayName', '$lim_0$');
%         end
%         ylim([0, 2*max(fh)]);
%         semilogx(lamLog,fhT,'r', 'LineWidth', 2, 'DisplayName', 'Tikhonov');
        title('Logarithmic regularization score plot');
        ylabel('Regularization score');
        xlabel('lambda');
        grid on;
        legend;
    else
        x = lamLog;
        y = fh;
    end

end

