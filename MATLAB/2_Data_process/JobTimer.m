function t = JobTimer(jobname)
t = timer;
t.StartFcn = @(~,~)disp('A timer is started!');
t.TimerFcn = @(~,~)system(['cd ~/work/job/ && sbatch ' jobname]);% resubmit the job if not done
t.stopFcn = @(~,~)disp('A timer is stopped!');
t.StartDelay = 60 * 60 * 2 - 90;
t.ExecutionMode = 'singleShot';
start(t);
end

