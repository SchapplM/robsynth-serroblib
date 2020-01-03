% Calculate vector of inverse dynamics joint torques for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPRRR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:04:18
% EndTime: 2019-12-31 19:04:25
% DurationCPUTime: 4.25s
% Computational Cost: add. (2233->389), mult. (4669->532), div. (0->0), fcn. (3202->14), ass. (0->172)
t383 = cos(qJ(3));
t375 = sin(pkin(9));
t356 = pkin(1) * t375 + pkin(6);
t346 = t356 * qJD(1);
t379 = sin(qJ(3));
t439 = qJD(3) * t383;
t342 = t356 * qJDD(1);
t481 = -qJD(2) * qJD(3) - t342;
t426 = -t346 * t439 + t481 * t379;
t266 = -qJDD(3) * pkin(3) - qJDD(2) * t383 - t426;
t444 = qJD(1) * t383;
t355 = -qJD(4) + t444;
t371 = qJ(1) + pkin(9);
t362 = sin(t371);
t363 = cos(t371);
t409 = g(1) * t363 + g(2) * t362;
t391 = -g(3) * t383 + t379 * t409;
t484 = qJD(4) * pkin(7) * t355 - t266 + t391;
t378 = sin(qJ(4));
t382 = cos(qJ(4));
t433 = t382 * qJD(3);
t445 = qJD(1) * t379;
t330 = t378 * t445 - t433;
t381 = cos(qJ(5));
t441 = qJD(3) * t378;
t332 = t382 * t445 + t441;
t377 = sin(qJ(5));
t464 = t332 * t377;
t281 = t381 * t330 + t464;
t353 = -qJD(5) + t355;
t483 = t281 * t353;
t404 = t330 * t377 - t381 * t332;
t482 = t353 * t404;
t306 = qJD(2) * t383 - t379 * t346;
t480 = qJD(3) * t306;
t437 = qJD(4) * t379;
t479 = -qJD(1) * t437 + qJDD(3);
t307 = t379 * qJD(2) + t383 * t346;
t297 = qJD(3) * pkin(7) + t307;
t376 = cos(pkin(9));
t357 = -pkin(1) * t376 - pkin(2);
t323 = -pkin(3) * t383 - pkin(7) * t379 + t357;
t298 = t323 * qJD(1);
t262 = t297 * t382 + t298 * t378;
t257 = -pkin(8) * t330 + t262;
t435 = qJD(5) * t377;
t255 = t257 * t435;
t296 = -qJD(3) * pkin(3) - t306;
t277 = pkin(4) * t330 + t296;
t374 = qJ(4) + qJ(5);
t369 = sin(t374);
t370 = cos(t374);
t460 = t370 * t383;
t292 = -t362 * t460 + t363 * t369;
t294 = t362 * t369 + t363 * t460;
t471 = g(3) * t379;
t478 = g(1) * t294 - g(2) * t292 + t277 * t281 + t370 * t471 + t255;
t461 = t369 * t383;
t291 = t362 * t461 + t363 * t370;
t293 = t362 * t370 - t363 * t461;
t265 = qJDD(3) * pkin(7) + qJDD(2) * t379 + t342 * t383 + t480;
t410 = pkin(3) * t379 - pkin(7) * t383;
t339 = t410 * qJD(3);
t278 = qJD(1) * t339 + qJDD(1) * t323;
t272 = t382 * t278;
t432 = qJD(1) * qJD(3);
t419 = t383 * t432;
t430 = qJDD(1) * t379;
t275 = qJD(4) * t433 + (t419 + t430) * t382 + t479 * t378;
t366 = t383 * qJDD(1);
t326 = t379 * t432 + qJDD(4) - t366;
t244 = pkin(4) * t326 - pkin(8) * t275 - qJD(4) * t262 - t265 * t378 + t272;
t276 = (qJD(3) * (qJD(4) + t444) + t430) * t378 - t479 * t382;
t436 = qJD(4) * t382;
t428 = t382 * t265 + t378 * t278 + t298 * t436;
t438 = qJD(4) * t378;
t397 = -t297 * t438 + t428;
t245 = -pkin(8) * t276 + t397;
t417 = t381 * t244 - t377 * t245;
t477 = -g(1) * t293 + g(2) * t291 + t277 * t404 + t369 * t471 + t417;
t322 = qJDD(5) + t326;
t476 = t322 * MDP(23) + (-t281 ^ 2 + t404 ^ 2) * MDP(20) - t281 * MDP(19) * t404;
t334 = t377 * t382 + t378 * t381;
t308 = t334 * t379;
t395 = -t378 * t437 + t383 * t433;
t474 = qJD(4) + qJD(5);
t416 = t275 * t377 + t381 * t276;
t249 = -qJD(5) * t404 + t416;
t472 = pkin(7) + pkin(8);
t261 = -t297 * t378 + t382 * t298;
t256 = -pkin(8) * t332 + t261;
t253 = -pkin(4) * t355 + t256;
t469 = t253 * t381;
t468 = t257 * t381;
t467 = t275 * t378;
t466 = t330 * t355;
t465 = t332 * t355;
t463 = t355 * t382;
t462 = t356 * t378;
t459 = t378 * t379;
t458 = t378 * t383;
t457 = t379 * t382;
t456 = t382 * t383;
t455 = qJDD(2) - g(3);
t420 = t378 * t439;
t259 = -t435 * t459 + (t474 * t457 + t420) * t381 + t395 * t377;
t454 = t259 * t353 - t308 * t322;
t333 = t377 * t378 - t381 * t382;
t398 = t333 * t383;
t453 = qJD(1) * t398 - t474 * t333;
t452 = (-t444 + t474) * t334;
t336 = t410 * qJD(1);
t451 = t382 * t306 + t378 * t336;
t450 = t323 * t436 + t378 * t339;
t440 = qJD(3) * t379;
t449 = t382 * t339 + t440 * t462;
t335 = t356 * t456;
t448 = t378 * t323 + t335;
t372 = t379 ^ 2;
t447 = -t383 ^ 2 + t372;
t347 = qJD(1) * t357;
t442 = qJD(3) * t330;
t434 = qJD(5) * t381;
t427 = t381 * t275 - t377 * t276 - t330 * t434;
t425 = qJD(4) * t472;
t424 = t378 * t444;
t423 = t355 * t441;
t415 = t355 * t356 + t297;
t414 = -qJD(4) * t298 - t265;
t413 = qJD(5) * t253 + t245;
t411 = -t307 + (-t424 + t438) * pkin(4);
t380 = sin(qJ(1));
t384 = cos(qJ(1));
t408 = g(1) * t380 - g(2) * t384;
t320 = t382 * t336;
t349 = t472 * t382;
t403 = pkin(4) * t379 - pkin(8) * t456;
t407 = qJD(1) * t403 + qJD(5) * t349 - t306 * t378 + t382 * t425 + t320;
t348 = t472 * t378;
t406 = pkin(8) * t424 - qJD(5) * t348 - t378 * t425 - t451;
t247 = t253 * t377 + t468;
t258 = -qJD(3) * t398 - t474 * t308;
t309 = t333 * t379;
t405 = t258 * t353 + t309 * t322;
t401 = -t326 * t378 + t355 * t436;
t400 = -t326 * t382 - t355 * t438;
t248 = -t332 * t435 + t427;
t396 = -qJD(1) * t347 + t409;
t394 = t379 * t436 + t420;
t393 = -pkin(7) * t326 - t296 * t355;
t392 = 0.2e1 * qJD(3) * t347 - qJDD(3) * t356;
t385 = qJD(3) ^ 2;
t388 = g(1) * t362 - g(2) * t363 - 0.2e1 * qJDD(1) * t357 - t356 * t385;
t361 = -pkin(4) * t382 - pkin(3);
t341 = qJDD(3) * t383 - t379 * t385;
t340 = qJDD(3) * t379 + t383 * t385;
t314 = (pkin(4) * t378 + t356) * t379;
t312 = t332 * t440;
t311 = t382 * t323;
t304 = t362 * t378 + t363 * t456;
t303 = t362 * t382 - t363 * t458;
t302 = -t362 * t456 + t363 * t378;
t301 = t362 * t458 + t363 * t382;
t287 = pkin(4) * t394 + t356 * t439;
t274 = -pkin(8) * t459 + t448;
t273 = t404 * t440;
t267 = -pkin(8) * t457 + t311 + (-pkin(4) - t462) * t383;
t252 = (-t379 * t433 - t383 * t438) * t356 - t394 * pkin(8) + t450;
t251 = pkin(4) * t276 + t266;
t250 = t403 * qJD(3) + (-t335 + (pkin(8) * t379 - t323) * t378) * qJD(4) + t449;
t246 = -t257 * t377 + t469;
t1 = [(-t275 * t383 + t326 * t457 - t355 * t395 + t312) * MDP(14) + ((t276 + t423) * t383 + (t401 - t442) * t379) * MDP(15) + (-t326 * t383 - t355 * t440) * MDP(16) + (-(-t323 * t438 + t449) * t355 + t311 * t326 - g(1) * t302 - g(2) * t304 + (t356 * t442 - t272 + t415 * t436 + (qJD(3) * t296 - t326 * t356 - t414) * t378) * t383 + (qJD(3) * t261 + t266 * t378 + t276 * t356 + t296 * t436) * t379) * MDP(17) + (t450 * t355 - t448 * t326 - g(1) * t301 - g(2) * t303 + (-t415 * t438 + (t296 * t382 + t332 * t356) * qJD(3) + t428) * t383 + (-t296 * t438 + t266 * t382 + t356 * t275 + (-t356 * t463 - t262) * qJD(3)) * t379) * MDP(18) + (-t248 * t309 - t258 * t404) * MDP(19) + (-t248 * t308 + t249 * t309 - t258 * t281 + t259 * t404) * MDP(20) + (-t248 * t383 - t273 - t405) * MDP(21) + (t249 * t383 - t281 * t440 + t454) * MDP(22) + (-t322 * t383 - t353 * t440) * MDP(23) + (-(t250 * t381 - t252 * t377) * t353 + (t267 * t381 - t274 * t377) * t322 - t417 * t383 + t246 * t440 + t287 * t281 + t314 * t249 + t251 * t308 + t277 * t259 - g(1) * t292 - g(2) * t294 + (-(-t267 * t377 - t274 * t381) * t353 + t247 * t383) * qJD(5)) * MDP(24) + (-t247 * t440 - g(1) * t291 - g(2) * t293 + t314 * t248 - t251 * t309 - t255 * t383 + t277 * t258 - t287 * t404 + ((-qJD(5) * t274 + t250) * t353 - t267 * t322 + t244 * t383) * t377 + ((qJD(5) * t267 + t252) * t353 - t274 * t322 + t413 * t383) * t381) * MDP(25) + qJDD(1) * MDP(1) + (t408 + (t375 ^ 2 + t376 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t372 + 0.2e1 * t379 * t419) * MDP(5) + 0.2e1 * (t366 * t379 - t432 * t447) * MDP(6) + t340 * MDP(7) + t341 * MDP(8) + (t379 * t392 + t383 * t388) * MDP(10) + (-t379 * t388 + t383 * t392) * MDP(11) + (t275 * t457 + t332 * t395) * MDP(12) + ((-t330 * t382 - t332 * t378) * t439 + (-t467 - t276 * t382 + (t330 * t378 - t332 * t382) * qJD(4)) * t379) * MDP(13) + t408 * MDP(2) + (g(1) * t384 + g(2) * t380) * MDP(3); t455 * MDP(4) + t341 * MDP(10) - t340 * MDP(11) + t312 * MDP(18) + t454 * MDP(24) + (-t273 + t405) * MDP(25) + ((-t276 + t423) * MDP(17) + (t355 * t433 - t275) * MDP(18) - t249 * MDP(24) - t248 * MDP(25)) * t383 + ((t401 + t442) * MDP(17) + t400 * MDP(18) + qJD(3) * t281 * MDP(24)) * t379; MDP(7) * t430 + MDP(8) * t366 + qJDD(3) * MDP(9) + (qJD(3) * t307 + t379 * t396 + t383 * t455 + t426) * MDP(10) + (t480 + (qJD(3) * t346 - t455) * t379 + (t396 + t481) * t383) * MDP(11) + (-t332 * t463 + t467) * MDP(12) + ((t275 + t466) * t382 + (-t276 + t465) * t378) * MDP(13) + ((-t332 * t379 + t355 * t456) * qJD(1) - t401) * MDP(14) + ((t330 * t379 - t355 * t458) * qJD(1) - t400) * MDP(15) + (-pkin(3) * t276 - t307 * t330 + t320 * t355 + (-t306 * t355 + t393) * t378 + t484 * t382) * MDP(17) + (-pkin(3) * t275 - t307 * t332 - t451 * t355 - t484 * t378 + t393 * t382) * MDP(18) + (t248 * t334 - t404 * t453) * MDP(19) + (-t248 * t333 - t249 * t334 - t281 * t453 + t404 * t452) * MDP(20) + (t322 * t334 - t353 * t453) * MDP(21) + (-t322 * t333 + t353 * t452) * MDP(22) + ((-t348 * t381 - t349 * t377) * t322 + t361 * t249 + t251 * t333 + (t377 * t406 + t381 * t407) * t353 + t411 * t281 + t452 * t277 + t391 * t370) * MDP(24) + (-(-t348 * t377 + t349 * t381) * t322 + t361 * t248 + t251 * t334 + (-t377 * t407 + t381 * t406) * t353 - t411 * t404 + t453 * t277 - t391 * t369) * MDP(25) + (MDP(16) * t355 - MDP(17) * t261 + MDP(18) * t262 + MDP(21) * t404 + MDP(22) * t281 + MDP(23) * t353 - MDP(24) * t246 + MDP(25) * t247) * t445 + (-MDP(5) * t379 * t383 + MDP(6) * t447) * qJD(1) ^ 2; t332 * t330 * MDP(12) + (-t330 ^ 2 + t332 ^ 2) * MDP(13) + (t275 - t466) * MDP(14) + (-t276 - t465) * MDP(15) + t326 * MDP(16) + (-t297 * t436 - g(1) * t303 + g(2) * t301 - t262 * t355 - t296 * t332 + t272 + (t414 + t471) * t378) * MDP(17) + (g(1) * t304 - g(2) * t302 + g(3) * t457 - t261 * t355 + t296 * t330 - t397) * MDP(18) + (t248 - t483) * MDP(21) + (-t249 + t482) * MDP(22) + ((-t256 * t377 - t468) * t353 - t247 * qJD(5) + (-t281 * t332 + t322 * t381 + t353 * t435) * pkin(4) + t477) * MDP(24) + ((t257 * t353 - t244) * t377 + (-t256 * t353 - t413) * t381 + (-t322 * t377 + t332 * t404 + t353 * t434) * pkin(4) + t478) * MDP(25) + t476; (t427 - t483) * MDP(21) + (-t416 + t482) * MDP(22) + (-t247 * t353 + t477) * MDP(24) + (-t377 * t244 - t381 * t245 - t246 * t353 + t478) * MDP(25) + (-MDP(21) * t464 + MDP(22) * t404 - MDP(24) * t247 - MDP(25) * t469) * qJD(5) + t476;];
tau = t1;
