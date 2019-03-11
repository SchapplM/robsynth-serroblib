% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:42:45
% EndTime: 2019-03-09 01:42:51
% DurationCPUTime: 4.10s
% Computational Cost: add. (2465->382), mult. (5125->455), div. (0->0), fcn. (3742->14), ass. (0->180)
t367 = sin(pkin(10));
t373 = sin(qJ(4));
t441 = qJD(4) * t373;
t422 = t367 * t441;
t369 = cos(pkin(10));
t472 = cos(qJ(4));
t424 = t472 * t369;
t414 = qJD(1) * t424;
t420 = qJDD(1) * t472;
t429 = qJDD(1) * t373;
t425 = qJD(4) * t414 + t367 * t420 + t369 * t429;
t287 = qJD(1) * t422 - t425;
t286 = -qJDD(6) + t287;
t375 = cos(qJ(6));
t279 = t375 * t286;
t336 = t367 * t472 + t373 * t369;
t482 = t336 * qJD(1);
t485 = qJD(6) + t482;
t372 = sin(qJ(6));
t489 = t485 * t372;
t491 = -t485 * t489 - t279;
t446 = t367 * t373;
t423 = qJD(1) * t446;
t325 = -t414 + t423;
t306 = qJD(4) * t372 - t375 * t325;
t490 = t306 * t485;
t418 = t375 * t485;
t308 = qJD(4) * t375 + t325 * t372;
t417 = t485 * t308;
t370 = cos(pkin(9));
t471 = pkin(1) * t370;
t352 = -pkin(2) - t471;
t431 = qJDD(1) * t352;
t337 = qJDD(3) + t431;
t366 = qJ(1) + pkin(9);
t358 = sin(t366);
t360 = cos(t366);
t483 = -g(1) * t358 + g(2) * t360;
t488 = -t337 - t483;
t487 = MDP(14) - MDP(17);
t365 = pkin(10) + qJ(4);
t357 = sin(t365);
t359 = cos(t365);
t413 = g(1) * t360 + g(2) * t358;
t385 = -g(3) * t357 - t359 * t413;
t368 = sin(pkin(9));
t347 = pkin(1) * t368 + qJ(3);
t333 = qJD(1) * qJD(3) + qJDD(1) * t347;
t354 = t369 * qJDD(2);
t298 = t354 + (-pkin(7) * qJDD(1) - t333) * t367;
t310 = t367 * qJDD(2) + t369 * t333;
t430 = qJDD(1) * t369;
t299 = pkin(7) * t430 + t310;
t339 = t347 * qJD(1);
t356 = t369 * qJD(2);
t462 = pkin(7) * qJD(1);
t304 = t356 + (-t339 - t462) * t367;
t318 = t367 * qJD(2) + t369 * t339;
t305 = t369 * t462 + t318;
t421 = qJD(4) * t472;
t415 = -t373 * t298 - t472 * t299 - t304 * t421 + t305 * t441;
t486 = t385 - t415;
t444 = t472 * t304 - t373 * t305;
t481 = qJD(5) - t444;
t267 = t373 * t304 + t472 * t305;
t265 = -qJD(4) * qJ(5) - t267;
t469 = pkin(5) * t325;
t255 = -t265 - t469;
t473 = pkin(4) + pkin(8);
t480 = t473 * t286 + (t255 - t267 + t469) * t485;
t330 = t336 * qJD(4);
t335 = -t424 + t446;
t438 = qJD(6) * t375;
t394 = t330 * t372 + t335 * t438;
t452 = t335 * t372;
t479 = t286 * t452 - t394 * t485;
t348 = g(3) * t359;
t478 = -t357 * t413 + t348;
t463 = pkin(7) + t347;
t331 = t463 * t367;
t332 = t463 * t369;
t268 = (qJD(3) * t367 + qJD(4) * t332) * t373 - qJD(3) * t424 + t331 * t421;
t281 = -t373 * t331 + t332 * t472;
t477 = qJD(4) * t268 - qJDD(4) * t281 + t357 * t483;
t269 = qJD(3) * t336 + qJD(4) * t281;
t280 = t331 * t472 + t373 * t332;
t476 = -qJD(4) * t269 - qJDD(4) * t280 - t359 * t483;
t475 = t325 ^ 2;
t474 = t482 ^ 2;
t408 = t367 * t429 - t369 * t420;
t288 = qJD(1) * t330 + t408;
t470 = pkin(4) * t288;
t374 = sin(qJ(1));
t466 = g(1) * t374;
t461 = qJ(5) * t325;
t460 = qJDD(4) * pkin(4);
t426 = t375 * qJDD(4) + t372 * t288 + t325 * t438;
t432 = qJD(4) * qJD(6);
t260 = -t372 * t432 + t426;
t459 = t260 * t375;
t351 = pkin(3) * t369 + pkin(2);
t338 = -t351 - t471;
t389 = -qJ(5) * t336 + t338;
t270 = t335 * t473 + t389;
t458 = t270 * t286;
t457 = t286 * t372;
t456 = t306 * t325;
t455 = t308 * t325;
t453 = t325 * t482;
t451 = t358 * t372;
t450 = t358 * t375;
t449 = t360 * t372;
t448 = t360 * t375;
t447 = t367 * MDP(6);
t329 = -t369 * t421 + t422;
t445 = t260 * t336 - t308 * t329;
t442 = t367 ^ 2 + t369 ^ 2;
t324 = qJD(1) * t338 + qJD(3);
t384 = -qJ(5) * t482 + t324;
t262 = t325 * t473 + t384;
t440 = qJD(6) * t262;
t439 = qJD(6) * t335;
t437 = t267 * qJD(4);
t435 = pkin(5) * t482 + t481;
t428 = qJDD(4) * qJ(5);
t321 = qJDD(1) * t338 + qJDD(3);
t383 = qJ(5) * t287 + t321;
t382 = -qJD(5) * t482 + t383;
t249 = t288 * t473 + t382;
t254 = -qJD(4) * t473 + t435;
t416 = qJD(6) * t254 + t249;
t376 = cos(qJ(1));
t411 = -g(2) * t376 + t466;
t410 = pkin(4) * t359 + qJ(5) * t357;
t407 = -t440 + t348;
t248 = t254 * t372 + t262 * t375;
t283 = t375 * t288;
t261 = qJD(6) * t308 + qJDD(4) * t372 - t283;
t405 = -t261 * t336 + t306 * t329;
t404 = -t287 * t335 + t330 * t482;
t403 = -t288 * t336 + t325 * t329;
t309 = -t333 * t367 + t354;
t402 = -t309 * t367 + t310 * t369;
t401 = (-t339 * t367 + t356) * t367 - t318 * t369;
t400 = qJ(5) * t329 - qJD(5) * t336;
t289 = -qJD(4) * t329 + qJDD(4) * t336;
t396 = qJD(4) * t330 + qJDD(4) * t335;
t395 = t351 + t410;
t393 = -t298 * t472 + t373 * t299 + t304 * t441 + t305 * t421;
t392 = -t335 * t279 + t330 * t418 - t439 * t489;
t250 = -qJD(4) * qJD(5) + t415 - t428;
t246 = -pkin(5) * t288 - t250;
t272 = t336 * pkin(5) + t280;
t391 = t246 * t335 + t255 * t330 + t272 * t286;
t390 = -t418 * t485 + t457;
t387 = qJDD(5) + t393;
t381 = t246 + (t473 * t485 + t461) * t485 + t385;
t275 = pkin(4) * t325 + t384;
t380 = t275 * t482 + t387 + t478;
t371 = -pkin(7) - qJ(3);
t362 = t376 * pkin(1);
t320 = qJD(4) * t325;
t314 = -t357 * t451 + t448;
t313 = t357 * t450 + t449;
t312 = t357 * t449 + t450;
t311 = t357 * t448 - t451;
t285 = pkin(4) * t482 + t461;
t278 = pkin(4) * t335 + t389;
t276 = pkin(4) * t330 + t400;
t273 = -t335 * pkin(5) + t281;
t264 = -qJD(4) * pkin(4) + t481;
t263 = t330 * t473 + t400;
t259 = -t329 * pkin(5) + t269;
t258 = -pkin(5) * t330 - t268;
t252 = t382 + t470;
t251 = t387 - t460;
t247 = t254 * t375 - t262 * t372;
t245 = -t287 * pkin(5) - qJDD(4) * t473 + t387;
t244 = t375 * t245;
t1 = [(t250 * t335 + t251 * t336 - t264 * t329 + t265 * t330 + t268 * t325 + t269 * t482 - t280 * t287 - t281 * t288 - t413) * MDP(16) + qJDD(1) * MDP(1) + t289 * MDP(11) - t396 * MDP(12) + (-t252 * t335 - t275 * t330 - t276 * t325 - t278 * t288 - t476) * MDP(17) + (t288 * t338 + t321 * t335 + t324 * t330 + t476) * MDP(14) + (-t252 * t336 + t275 * t329 - t276 * t482 + t278 * t287 - t477) * MDP(18) + (-t287 * t338 + t321 * t336 - t324 * t329 + t477) * MDP(15) + (t445 - t479) * MDP(22) + (t411 + (t368 ^ 2 + t370 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t333 * t442 + t402 - t413) * MDP(7) + (-t286 * t336 - t329 * t485) * MDP(24) + (-t287 * t336 - t329 * t482) * MDP(9) + (t337 * t352 - g(1) * (-pkin(1) * t374 - pkin(2) * t358 + qJ(3) * t360) - g(2) * (pkin(2) * t360 + qJ(3) * t358 + t362) + t402 * t347 - t401 * qJD(3)) * MDP(8) + (t403 - t404) * MDP(10) + (t392 + t405) * MDP(23) + (g(1) * t376 + g(2) * t374) * MDP(3) + t411 * MDP(2) + (t260 * t452 + t308 * t394) * MDP(20) + (g(1) * t313 - g(2) * t311 + t248 * t329 + t258 * t308 + t273 * t260 + (-(qJD(6) * t272 + t263) * t485 + t458 - t416 * t336 + t255 * t439) * t375 + (-(-qJD(6) * t270 + t259) * t485 - (t245 - t440) * t336 + t391) * t372) * MDP(26) + (-g(1) * t314 - g(2) * t312 + t244 * t336 - t247 * t329 + t258 * t306 + t273 * t261 + (-t249 * t336 - t263 * t485 + t458) * t372 + (t259 * t485 - t391) * t375 + ((-t270 * t375 - t272 * t372) * t485 - t248 * t336 + t255 * t452) * qJD(6)) * MDP(25) + ((-t306 * t372 + t308 * t375) * t330 + (t459 - t261 * t372 + (-t306 * t375 - t308 * t372) * qJD(6)) * t335) * MDP(21) + (pkin(1) * t466 - g(2) * t362 - t250 * t281 + t251 * t280 + t252 * t278 + t264 * t269 + t265 * t268 + t275 * t276 + (g(1) * t371 - g(2) * t395) * t360 + (g(1) * t395 + g(2) * t371) * t358) * MDP(19) + (t369 * MDP(5) - t447) * (-t431 + t488); (qJDD(2) - g(3)) * MDP(4) + (t309 * t369 + t310 * t367 - g(3)) * MDP(8) + (t403 + t404) * MDP(16) + (-t250 * t336 + t251 * t335 + t264 * t330 + t265 * t329 - g(3)) * MDP(19) + (t392 - t405) * MDP(25) + (t445 + t479) * MDP(26) + (-MDP(15) + MDP(18)) * t289 - t487 * t396; -MDP(5) * t430 + qJDD(1) * t447 - t442 * MDP(7) * qJD(1) ^ 2 + (qJD(1) * t401 - t488) * MDP(8) + ((-t325 - t423) * qJD(4) + t425) * MDP(15) + (-t474 - t475) * MDP(16) + (t287 + t320) * MDP(18) + (t470 - t265 * t325 + (-qJD(5) - t264) * t482 + t383 + t483) * MDP(19) + (t390 + t456) * MDP(25) + (t455 - t491) * MDP(26) + t487 * (0.2e1 * t482 * qJD(4) + t408); MDP(9) * t453 + (t474 - t475) * MDP(10) + ((t325 - t423) * qJD(4) + t425) * MDP(11) - t408 * MDP(12) + qJDD(4) * MDP(13) + (-t324 * t482 - t393 + t437 - t478) * MDP(14) + (qJD(4) * t444 + t324 * t325 - t486) * MDP(15) + (pkin(4) * t287 - qJ(5) * t288 + (-t265 - t267) * t482 + (t264 - t481) * t325) * MDP(16) + (t285 * t325 + t380 - t437 - 0.2e1 * t460) * MDP(17) + (0.2e1 * t428 - t275 * t325 + t285 * t482 + (0.2e1 * qJD(5) - t444) * qJD(4) + t486) * MDP(18) + (-t251 * pkin(4) - g(3) * t410 - t250 * qJ(5) - t264 * t267 - t481 * t265 - t275 * t285 + t413 * (pkin(4) * t357 - qJ(5) * t359)) * MDP(19) + (-t372 * t417 + t459) * MDP(20) + ((-t261 - t417) * t375 + (-t260 + t490) * t372) * MDP(21) + (t455 + t491) * MDP(22) + (t390 - t456) * MDP(23) + t485 * t325 * MDP(24) + (qJ(5) * t261 + t247 * t325 + t435 * t306 + t381 * t372 + t375 * t480) * MDP(25) + (qJ(5) * t260 - t248 * t325 + t435 * t308 - t372 * t480 + t381 * t375) * MDP(26); (-t287 + t320) * MDP(16) + (qJDD(4) - t453) * MDP(17) + (-qJD(4) ^ 2 - t474) * MDP(18) + (t265 * qJD(4) + t380 - t460) * MDP(19) + (-qJD(4) * t306 - t279) * MDP(25) + (-qJD(4) * t308 + t457) * MDP(26) + (-MDP(25) * t489 - MDP(26) * t418) * t485; t308 * t306 * MDP(20) + (-t306 ^ 2 + t308 ^ 2) * MDP(21) + (t426 + t490) * MDP(22) + (t283 + t417) * MDP(23) - t286 * MDP(24) + (-g(1) * t311 - g(2) * t313 + t248 * t485 - t255 * t308 + t244) * MDP(25) + (g(1) * t312 - g(2) * t314 + t247 * t485 + t255 * t306) * MDP(26) + (-MDP(23) * t432 + MDP(25) * t407 - MDP(26) * t416) * t375 + (-MDP(22) * t432 + (-qJD(6) * t325 - qJDD(4)) * MDP(23) - t416 * MDP(25) + (-t245 - t407) * MDP(26)) * t372;];
tau  = t1;
