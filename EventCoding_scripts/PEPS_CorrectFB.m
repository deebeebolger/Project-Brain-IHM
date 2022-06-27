function [fbcurr, incorr] = PEPS_CorrectFB(fborig, fbcurr, currentdir, fname)
% Programmed by: D. Bolger .                     Date: 5-02-2020
% Function to correct the feedbacks based on a comparison of ITI between
% the original feedback list and the feedback list of the current subject.
% It takes in the table resuming the feedbacks information of the current subject and
% the a table resuming the defined feedback properties.
%***************************************************************
% Create a text file in which to write the feedback information for the
% current subject.

fID = fopen(fullfile(currentdir,[fname,'.txt']),'a+');


% Begin by comparing the ITIs of the original and the current subject.
% If the ITI difference < 1second, then check if both are assigned the same
% trial number and feedback type. If this is not the case, assign the
% correct feedback category and trial number to the current feedback.
% Here we are presuming that the length of the current subject feedback
% list is same or greater than that of the original feedback list.
incorr=0;
cntr_orig=1;
cntr_curr =1;
while cntr_orig<=length(fborig.ITI)
    
    if fbcurr.ITI(cntr_curr)<0  % Sometimes weird negative ITIs happen...
        fbcurr.ITI(cntr_curr)=0;
    end
    
    diffcurr = abs(fborig.ITI(cntr_orig)-fbcurr.ITI(cntr_curr));   % Find the difference between the original ITI and the current ITI
    if diffcurr<1         % If the ITIs match
        
        fbcurr.Feedbacks{cntr_curr}=fborig.Feedbacks{cntr_orig};            % Assign the current feedback type.
        fbcurr.Trial_numbers(cntr_curr) = fborig.Trial_numbers(cntr_orig);  % Assign the current feedback trial number.
        cntr_orig = cntr_orig+1;                                            % Increment the original counter
        cntr_curr = cntr_curr+1;                                            % Increment the current subject feedback counter
    else
        fbcurr.Feedbacks{cntr_curr}='incorrect-fb';         % If the ITIs do not match for the current index assign "incorrect" title
        incorr = incorr+1;                                  % Increment the number of incorrect feedbacks.
        fbcurr.incorrectindx(incorr) = cntr_curr;           % Add the index to the current subject feedback table.
        cntr_orig = cntr_orig+1;                            % Increment only the original feedback information table. 
    end
    
end



 %Write the missing feedbacks to table and save.


end