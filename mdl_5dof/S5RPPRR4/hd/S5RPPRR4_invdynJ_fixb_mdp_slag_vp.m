% Calculate vector of inverse dynamics joint torques for
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RPPRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:31:38
% EndTime: 2020-01-03 11:31:49
% DurationCPUTime: 4.77s
% Computational Cost: add. (2093->339), mult. (5244->460), div. (0->0), fcn. (4071->14), ass. (0->172)
t466 = qJDD(1) * pkin(1);
t362 = qJDD(2) - t466;
t383 = sin(qJ(1));
t386 = cos(qJ(1));
t447 = -g(2) * t386 - g(3) * t383;
t481 = t447 - t362;
t378 = sin(pkin(8));
t377 = sin(pkin(9));
t379 = cos(pkin(9));
t382 = sin(qJ(4));
t385 = cos(qJ(4));
t334 = t377 * t385 + t379 * t382;
t395 = qJD(1) * t334;
t306 = t378 * t395;
t384 = cos(qJ(5));
t443 = qJD(1) * t378;
t426 = t377 * t443;
t411 = t382 * t426;
t455 = t379 * t385;
t429 = t378 * t455;
t309 = qJD(1) * t429 - t411;
t381 = sin(qJ(5));
t459 = t309 * t381;
t263 = t384 * t306 + t459;
t380 = cos(pkin(8));
t442 = qJD(1) * t380;
t479 = qJD(4) - t442;
t351 = -qJD(5) - t479;
t462 = t263 * t351;
t439 = qJD(3) * t378;
t440 = qJD(2) * t380;
t326 = -t377 * t440 - t379 * t439;
t327 = -t377 * t439 + t379 * t440;
t338 = -pkin(2) * t380 - qJ(3) * t378 - pkin(1);
t332 = t379 * t338;
t393 = -pkin(6) * t378 * t379 + (-qJ(2) * t377 - pkin(3)) * t380;
t284 = t332 + t393;
t300 = t379 * t380 * qJ(2) + t377 * t338;
t458 = t377 * t378;
t289 = -pkin(6) * t458 + t300;
t451 = t382 * t284 + t385 * t289;
t480 = qJD(4) * t451 - t385 * t326 + t327 * t382;
t324 = qJD(1) * t338 + qJD(2);
t314 = t379 * t324;
t276 = qJD(1) * t393 + t314;
t427 = qJ(2) * t442;
t288 = t377 * t324 + t379 * t427;
t281 = -pkin(6) * t426 + t288;
t404 = -t276 * t382 - t281 * t385;
t251 = -pkin(7) * t306 - t404;
t436 = qJD(5) * t381;
t249 = t251 * t436;
t352 = qJ(2) * t443 + qJD(3);
t325 = pkin(3) * t426 + t352;
t279 = pkin(4) * t306 + t325;
t376 = pkin(9) + qJ(4);
t366 = qJ(5) + t376;
t361 = cos(t366);
t360 = sin(t366);
t452 = t386 * t360;
t454 = t380 * t383;
t302 = t361 * t454 - t452;
t453 = t380 * t386;
t304 = t360 * t383 + t361 * t453;
t469 = g(1) * t378;
t478 = g(2) * t302 - g(3) * t304 + t279 * t263 + t361 * t469 + t249;
t336 = qJD(4) * t411;
t437 = qJD(4) * t385;
t425 = t379 * t437;
t275 = (qJD(1) * t425 + qJDD(1) * t334) * t378 - t336;
t402 = -t306 * t381 + t384 * t309;
t477 = t263 * MDP(19) * t402 + (-t263 ^ 2 + t402 ^ 2) * MDP(20);
t329 = t334 * qJD(4);
t339 = qJDD(1) * t429;
t456 = t377 * t382;
t274 = t339 + (-qJD(1) * t329 - qJDD(1) * t456) * t378;
t430 = qJDD(1) * t380;
t354 = -qJDD(4) + t430;
t297 = -qJD(1) * t439 + qJDD(1) * t338 + qJDD(2);
t292 = t379 * t297;
t432 = qJD(1) * qJD(2);
t424 = t380 * t432;
t259 = qJDD(1) * t393 - t377 * t424 + t292;
t422 = t379 * t430;
t273 = qJ(2) * t422 + t377 * t297 + t379 * t424;
t431 = qJDD(1) * t378;
t423 = t377 * t431;
t261 = -pkin(6) * t423 + t273;
t390 = qJD(4) * t404 + t385 * t259 - t382 * t261;
t242 = -pkin(4) * t354 - pkin(7) * t274 + t390;
t475 = t381 * t242;
t463 = t402 * t351;
t348 = -qJDD(5) + t354;
t337 = t348 * MDP(23);
t474 = t354 * MDP(16) + t337;
t417 = t385 * t276 - t281 * t382;
t250 = -pkin(7) * t309 + t417;
t248 = pkin(4) * t479 + t250;
t464 = t251 * t384;
t405 = -t248 * t381 - t464;
t472 = t405 * qJD(5);
t433 = qJ(2) * qJDD(1);
t301 = -t360 * t454 - t361 * t386;
t303 = -t361 * t383 + t380 * t452;
t438 = qJD(4) * t382;
t397 = -t382 * t259 - t385 * t261 - t276 * t437 + t281 * t438;
t243 = -pkin(7) * t275 - t397;
t420 = t384 * t242 - t381 * t243;
t471 = -g(2) * t301 - g(3) * t303 - t279 * t402 + t360 * t469 + t420;
t418 = t274 * t381 + t384 * t275;
t245 = qJD(5) * t402 + t418;
t468 = g(2) * t383;
t467 = g(3) * t386;
t465 = t248 * t384;
t461 = t306 * t479;
t460 = t309 * t479;
t457 = t377 * t380;
t450 = -t380 * t395 + t329;
t333 = t455 - t456;
t449 = t479 * t333;
t335 = pkin(3) * t458 + t378 * qJ(2);
t448 = t386 * pkin(1) + t383 * qJ(2);
t373 = t378 ^ 2;
t446 = t380 ^ 2 + t373;
t445 = MDP(17) * t385;
t441 = qJD(2) * t378;
t435 = qJD(5) * t384;
t428 = t384 * t274 - t381 * t275 - t306 * t435;
t330 = qJ(2) * t431 + t378 * t432 + qJDD(3);
t387 = qJD(1) ^ 2;
t421 = t446 * t387;
t416 = t385 * t284 - t289 * t382;
t414 = MDP(10) * (-t377 ^ 2 - t379 ^ 2);
t413 = qJD(5) * t248 + t243;
t412 = 0.2e1 * t446;
t298 = pkin(3) * t423 + t330;
t408 = -t467 + t468;
t407 = qJD(5) * t334 + t450;
t406 = qJD(5) * t333 + t449;
t322 = t334 * t378;
t323 = t333 * t378;
t277 = t384 * t322 + t323 * t381;
t278 = -t322 * t381 + t323 * t384;
t401 = t432 + t433;
t272 = -t401 * t457 + t292;
t299 = -qJ(2) * t457 + t332;
t400 = -t326 * qJD(1) - t299 * qJDD(1) - t272;
t399 = t327 * qJD(1) + t300 * qJDD(1) + t273;
t396 = t284 * t437 - t289 * t438 + t382 * t326 + t385 * t327;
t244 = -t309 * t436 + t428;
t394 = t466 + t481;
t391 = t412 * t432 + t467;
t368 = t383 * pkin(1);
t364 = cos(t376);
t363 = sin(t376);
t318 = t363 * t383 + t364 * t453;
t317 = t363 * t453 - t364 * t383;
t316 = -t363 * t386 + t364 * t454;
t315 = -t363 * t454 - t364 * t386;
t312 = t378 * t425 - t438 * t458;
t311 = t378 * t329;
t290 = pkin(4) * t312 + t441;
t287 = -t377 * t427 + t314;
t286 = pkin(4) * t322 + t335;
t256 = pkin(4) * t275 + t298;
t255 = -pkin(7) * t322 + t451;
t254 = -pkin(4) * t380 - pkin(7) * t323 + t416;
t253 = qJD(5) * t278 - t311 * t381 + t384 * t312;
t252 = -qJD(5) * t277 - t311 * t384 - t312 * t381;
t247 = pkin(7) * t311 - t480;
t246 = -pkin(7) * t312 + t396;
t1 = [(t273 * t300 + t288 * t327 + t272 * t299 + t287 * t326 - g(2) * (pkin(2) * t453 + t448) - g(3) * (pkin(2) * t454 - qJ(2) * t386 + t368)) * MDP(11) + t447 * MDP(2) + (-t362 * pkin(1) - g(2) * t448 - g(3) * t368 + (t446 * t433 + t391) * qJ(2)) * MDP(7) + (g(2) * t317 - g(3) * t315 + t335 * t274 + t298 * t323 + t309 * t441 - t325 * t311 + t354 * t451 - t396 * t479) * MDP(18) + (-g(2) * t318 - g(3) * t316 + t335 * t275 + t298 * t322 + t306 * t441 + t325 * t312 - t416 * t354 - t479 * t480) * MDP(17) + (t412 * t433 + t391 - t468) * MDP(6) + (-(t254 * t384 - t255 * t381) * t348 + t290 * t263 + t286 * t245 + t256 * t277 + t279 * t253 - g(2) * t304 - g(3) * t302 + (t246 * t381 - t247 * t384 - (-t254 * t381 - t255 * t384) * qJD(5)) * t351) * MDP(24) + (g(2) * t303 - g(3) * t301 + t286 * t244 + t279 * t252 + t256 * t278 + t290 * t402 + ((-qJD(5) * t255 + t247) * t351 + t254 * t348) * t381 + ((qJD(5) * t254 + t246) * t351 + t255 * t348) * t384) * MDP(25) + qJDD(1) * MDP(1) + t408 * MDP(3) + (-t312 * t479 + t322 * t354) * MDP(15) + (-t311 * t479 - t323 * t354) * MDP(14) + (t253 * t351 + t277 * t348) * MDP(22) + (-t252 * t351 - t278 * t348) * MDP(21) + (-t274 * t322 - t275 * t323 + t306 * t311 - t309 * t312) * MDP(13) + (t244 * t278 + t252 * t402) * MDP(19) + (-t244 * t277 - t245 * t278 - t252 * t263 - t253 * t402) * MDP(20) + (t274 * t323 - t309 * t311) * MDP(12) + (t377 * MDP(8) + MDP(9) * t379) * (t330 * t378 + t373 * t401 - t408) + ((t330 * qJ(2) + qJ(3) * t447 + t352 * qJD(2)) * MDP(11) - t394 * MDP(5) + (-t377 * t399 + t379 * t400 + t447) * MDP(10)) * t378 + (-t397 * MDP(18) - t390 * MDP(17) + t394 * MDP(4) + (-t420 - t472) * MDP(24) + (t413 * t384 - t249 + t475) * MDP(25) + (-t377 * t447 + t399) * MDP(9) + (t379 * t447 + t400) * MDP(8) + t275 * MDP(15) - t274 * MDP(14) + t245 * MDP(22) - t244 * MDP(21) + t474) * t380; -MDP(4) * t430 - MDP(6) * t421 + (-qJ(2) * t421 - t481) * MDP(7) + (-t377 * t421 - t422) * MDP(8) + (t377 * t430 - t379 * t421) * MDP(9) + (t272 * t379 + t273 * t377 + (-t352 * t378 + (t287 * t377 - t288 * t379) * t380) * qJD(1) - t447) * MDP(11) + (-t306 * t443 - t333 * t354 - t450 * t479) * MDP(17) + (-t309 * t443 + t334 * t354 - t449 * t479) * MDP(18) + (-(t333 * t384 - t334 * t381) * t348 - t263 * t443 + (t381 * t406 + t384 * t407) * t351) * MDP(24) + ((t333 * t381 + t334 * t384) * t348 - t402 * t443 + (-t381 * t407 + t384 * t406) * t351) * MDP(25) + (MDP(5) + t414) * t431; (g(1) * t380 + t330) * MDP(11) + (-t336 + t460) * MDP(17) + (t339 - t461) * MDP(18) + (t245 - t463) * MDP(24) + (t244 + t462) * MDP(25) + t387 * t373 * t414 + (-t408 * MDP(11) + (-MDP(8) * t379 + MDP(9) * t377) * t387 * t380 + ((MDP(17) * t382 + MDP(9)) * t379 + (-MDP(18) * t382 + MDP(8) + t445) * t377) * qJDD(1) + ((t287 * t379 + t288 * t377) * MDP(11) + (-MDP(18) * t334 + t379 * t445) * qJD(4)) * qJD(1)) * t378; t309 * t306 * MDP(12) + (-t306 ^ 2 + t309 ^ 2) * MDP(13) + (t274 + t461) * MDP(14) + (-t275 + t460) * MDP(15) + (-g(2) * t315 - g(3) * t317 - t325 * t309 + t363 * t469 - t404 * t479 + t390) * MDP(17) + (g(2) * t316 - g(3) * t318 + t325 * t306 + t364 * t469 + t417 * t479 + t397) * MDP(18) + (t244 - t462) * MDP(21) + (-t245 - t463) * MDP(22) + ((-t250 * t381 - t464) * t351 + t472 + (-t263 * t309 - t348 * t384 + t351 * t436) * pkin(4) + t471) * MDP(24) + ((t251 * t351 - t242) * t381 + (-t250 * t351 - t413) * t384 + (-t309 * t402 + t348 * t381 + t351 * t435) * pkin(4) + t478) * MDP(25) - t474 + t477; (t428 - t462) * MDP(21) + (-t418 - t463) * MDP(22) - t337 + (t351 * t405 + t471) * MDP(24) + (-t384 * t243 - t475 - (-t251 * t381 + t465) * t351 + t478) * MDP(25) + (-MDP(21) * t459 - MDP(22) * t402 + t405 * MDP(24) - MDP(25) * t465) * qJD(5) + t477;];
tau = t1;
