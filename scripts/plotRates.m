clear;
clc;
close all;

% dir = '/home/mnowicki/Projects/Tardos/KITTI/orbslam2_klt/klt0.85_patchSize9/inliers/sequence_10/';
% dir = '/home/mnowicki/Projects/Tardos/KITTI/orbslam2_klt/klt0.9_run1/inliers/sequence_10/';
% dir = '/home/mnowicki/Projects/Tardos/KITTI/orbslam2_klt/klt0.95_run1/inliers/sequence_01/';
%dir = '/home/mnowicki/Projects/Tardos/KITTI/orbslam2_klt/orig_run1/inliers/sequence_01/';
dir = '/home/mnowicki/Projects/Tardos/ORB_SLAM2/logs/';

%
% candidates (1), matches (2), tracks (3), extra tracks over matches (4), inliers (5),
% inliersTracking (6), zncc (7)
%
voInliers = importdata(strcat(dir, 'voInlierCount.txt'));

%
% matches, tracks, candidates, inliers, total
%
mapInliers = importdata(strcat(dir, 'mapInlierCount.txt'));

figure;
plot(voInliers(:,1), 'r', 'LineWidth',2);
hold on;
plot(voInliers(:,2), 'g', 'LineWidth',2);
plot(voInliers(:,3), 'b', 'LineWidth',2);
plot(voInliers(:,2) + voInliers(:,4), 'c', 'LineWidth',2);
legend('Last inliers (candidates)', 'VO Matches', 'VO Trackings', 'VO Matches + extra tracks');
title('VO - frame2frame');

figure;
plot(voInliers(:,5), 'r', 'LineWidth',2);
hold on;
plot(voInliers(:,6), 'g', 'LineWidth',2);
plot(voInliers(:,5) - voInliers(:,6), 'b', 'LineWidth',2);
legend('Total inliers', 'Inliers matches', 'Inliers trackings');
title('VO - inlier analysis');

figure;
plot(voInliers(:,6), 'r', 'LineWidth',2);
hold on;
yyaxis right;
plot(voInliers(:,7), 'g', 'LineWidth',2);
legend('Inliers tracking', 'ZNCC');
title('Tracking inliers vs avg zncc');

clc
fprintf(['VO\n\tCand: %.2f \n\tMatches: %.2f \n\tTracks: %.2f \n\tMatches + tracks: %.2f (%.2f) ',...
    '\n\tInliers: %.2f (t:%.2f, m:%.2f)\n\tInlier perc. - t:%.2f m:%.2f\n\tAvg znnc of inliners: %.2f\n'], ...
    mean(voInliers(:,1)), mean(voInliers(:,2)), mean(voInliers(:,3)), mean(voInliers(:,2) + voInliers(:,4)), ...p[
    mean(voInliers(:,4)), mean(voInliers(:,5)), mean(voInliers(:,6)), mean(voInliers(:,5) - voInliers(:,6)), ...
    mean(voInliers(:,6)) / mean(voInliers(:,4)), ...
    (mean(voInliers(:,5)) - mean(voInliers(:,6))) / mean(voInliers(:,2)), ...
    mean(voInliers(:,7)));

figure;
plot(mapInliers(:,1), 'r', 'LineWidth',2);
hold on;
plot(mapInliers(:,2), 'g', 'LineWidth',2);
plot(mapInliers(:,3), 'b', 'LineWidth',2);
plot(mapInliers(:,4), 'c', 'LineWidth',2);
plot(mapInliers(:,5), 'm', 'LineWidth',2);
legend('Map matches', 'Map tracks', 'Map add cand', 'Map inliers', 'Map total cand');
title('Map - frame2localFrames');

fprintf('Map\n\tAdd cand: %f \n\tAdd matches: %f \n\tAdd tracks: %f \n\tInliers: %f \n\tTotal: %f\n\n', ...
    mean(mapInliers(:,3)), mean(mapInliers(:,1)), mean(mapInliers(:,2)), ...
    mean(mapInliers(:,4)), mean(mapInliers(:,5)))

% figure;
% plot(voInliers(:,3) ./voInliers(:,4), 'r', 'LineWidth',2);
% hold on;
% plot(mapInliers(:,4) ./ mapInliers(:,5), 'g', 'LineWidth',2);
% legend('VO Inlier Rate', 'Map Inlier Rate');
% title('Inlier rate');

% figure;
% plot(voInliers(:,3), 'r', 'LineWidth',2);
% hold on;
% plot(mapInliers(:,4), 'g', 'LineWidth',2);
% legend('VO Inlier Count', 'Map Inlier Count');
% title('Inlier count');
