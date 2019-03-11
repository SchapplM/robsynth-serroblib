% Calculate vector of inverse dynamics joint torques for
% S6RPPRRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRRP7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRRP7_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:14:03
% EndTime: 2019-03-09 02:14:09
% DurationCPUTime: 4.59s
% Computational Cost: add. (3564->406), mult. (7094->500), div. (0->0), fcn. (4961->10), ass. (0->174)
t367 = sin(pkin(9));
t373 = sin(qJ(4));
t431 = qJD(1) * t373;
t415 = t367 * t431;
t368 = cos(pkin(9));
t376 = cos(qJ(4));
t430 = qJD(1) * t376;
t418 = t368 * t430;
t323 = -t415 + t418;
t375 = cos(qJ(5));
t372 = sin(qJ(5));
t429 = qJD(4) * t372;
t302 = t323 * t375 + t429;
t329 = t367 * t376 + t368 * t373;
t474 = t329 * qJD(1);
t484 = qJD(5) + t474;
t489 = t302 * t484;
t371 = -pkin(1) - qJ(3);
t471 = -qJD(1) * qJD(3) + qJDD(1) * t371;
t330 = qJDD(2) + t471;
t408 = -pkin(7) * qJDD(1) + t330;
t306 = t408 * t367;
t307 = t408 * t368;
t338 = qJD(1) * t371 + qJD(2);
t410 = -pkin(7) * qJD(1) + t338;
t311 = t410 * t367;
t312 = t410 * t368;
t427 = qJD(4) * t376;
t428 = qJD(4) * t373;
t393 = -t373 * t306 + t307 * t376 - t311 * t427 - t312 * t428;
t257 = -qJDD(4) * pkin(4) - t393;
t363 = pkin(9) + qJ(4);
t350 = sin(t363);
t351 = cos(t363);
t374 = sin(qJ(1));
t377 = cos(qJ(1));
t477 = g(1) * t374 - g(2) * t377;
t383 = -g(3) * t350 + t351 * t477;
t426 = qJD(5) * t484;
t488 = pkin(8) * t426 + t257 + t383;
t328 = t367 * t373 - t376 * t368;
t386 = t328 * qJDD(1);
t424 = qJD(5) * t375;
t267 = -t375 * qJDD(4) + t323 * t424 - t372 * t386 + t429 * (qJD(5) - t474);
t399 = t306 * t376 + t307 * t373;
t479 = -t373 * t311 + t312 * t376;
t256 = qJDD(4) * pkin(8) + qJD(4) * t479 + t399;
t278 = t376 * t311 + t373 * t312;
t275 = qJD(4) * pkin(8) + t278;
t349 = qJD(1) * qJ(2) + qJD(3);
t355 = t367 * pkin(3);
t333 = qJD(1) * t355 + t349;
t276 = pkin(4) * t474 - pkin(8) * t323 + t333;
t259 = t275 * t375 + t276 * t372;
t417 = t368 * t427;
t337 = qJD(4) * t415;
t472 = -t329 * qJDD(1) + t337;
t293 = qJD(1) * t417 - t472;
t364 = qJDD(1) * qJ(2);
t365 = qJD(1) * qJD(2);
t476 = t364 + t365;
t336 = qJDD(3) + t476;
t327 = qJDD(1) * t355 + t336;
t380 = -qJD(4) * t474 - t386;
t263 = t293 * pkin(4) - pkin(8) * t380 + t327;
t262 = t375 * t263;
t423 = t375 * qJD(4);
t425 = qJD(5) * t372;
t266 = -qJD(5) * t423 - t372 * qJDD(4) + t323 * t425 - t375 * t380;
t290 = qJDD(5) + t293;
t244 = pkin(5) * t290 + qJ(6) * t266 - qJD(5) * t259 - qJD(6) * t302 - t256 * t372 + t262;
t300 = t323 * t372 - t423;
t251 = -qJ(6) * t300 + t259;
t487 = t251 * t484 + t244;
t388 = t375 * t256 + t372 * t263 - t275 * t425 + t276 * t424;
t245 = -qJ(6) * t267 - qJD(6) * t300 + t388;
t258 = -t275 * t372 + t375 * t276;
t250 = -qJ(6) * t302 + t258;
t249 = pkin(5) * t484 + t250;
t486 = -t249 * t484 + t245;
t402 = t367 * MDP(7) + t368 * MDP(8);
t405 = g(1) * t377 + g(2) * t374;
t481 = t405 * t351;
t478 = t336 - t405;
t432 = t367 ^ 2 + t368 ^ 2;
t443 = t375 * t377;
t446 = t372 * t374;
t315 = -t350 * t446 + t443;
t444 = t374 * t375;
t445 = t372 * t377;
t317 = t350 * t445 + t444;
t475 = -g(1) * t315 - g(2) * t317;
t464 = g(3) * t351;
t473 = t350 * t477 + t464;
t470 = t302 ^ 2;
t469 = 0.2e1 * t365;
t468 = pkin(5) * t372;
t463 = -pkin(7) + t371;
t369 = -qJ(6) - pkin(8);
t462 = pkin(1) * qJDD(1);
t461 = t266 * t372;
t460 = t267 * t375;
t459 = t300 * t474;
t458 = t300 * t323;
t457 = t300 * t372;
t456 = t300 * t375;
t455 = t302 * t323;
t454 = t302 * t372;
t453 = t302 * t375;
t451 = t474 * t372;
t450 = t328 * t372;
t449 = t328 * t375;
t283 = t375 * t290;
t331 = t463 * t367;
t332 = t463 * t368;
t296 = t331 * t376 + t332 * t373;
t294 = t375 * t296;
t343 = qJ(2) + t355;
t442 = -t250 + t249;
t441 = -t372 * t267 - t300 * t424;
t292 = pkin(4) * t323 + pkin(8) * t474;
t440 = t372 * t292 + t375 * t479;
t439 = -t451 * t484 + t283;
t291 = pkin(4) * t329 + pkin(8) * t328 + t343;
t438 = t372 * t291 + t294;
t324 = -t367 * t427 - t368 * t428;
t437 = t324 * qJD(4) - t328 * qJDD(4);
t411 = qJD(5) * t369;
t436 = -qJ(6) * t451 + qJD(6) * t375 + t372 * t411 - t440;
t286 = t375 * t292;
t435 = -pkin(5) * t323 - t286 + (-qJ(6) * t474 + t411) * t375 + (-qJD(6) + t479) * t372;
t434 = t377 * pkin(1) + t374 * qJ(2);
t295 = t331 * t373 - t332 * t376;
t271 = -qJD(3) * t329 - qJD(4) * t295;
t325 = -t367 * t428 + t417;
t288 = pkin(4) * t325 - pkin(8) * t324 + qJD(2);
t419 = t375 * t271 + t372 * t288 + t291 * t424;
t416 = t328 * t424;
t414 = g(2) * t434;
t413 = qJDD(2) - t477;
t412 = pkin(7) + qJ(3) + t468;
t409 = t432 * MDP(9);
t407 = t375 * t484;
t406 = -qJD(5) * t276 - t256;
t403 = -t275 * t424 + t262;
t401 = -t249 * t375 - t251 * t372;
t400 = -t249 * t372 + t251 * t375;
t348 = pkin(5) * t375 + pkin(4);
t396 = -t348 * t350 - t351 * t369;
t395 = -qJ(6) * t324 + qJD(6) * t328;
t394 = -t372 * MDP(23) - t375 * MDP(24);
t274 = -qJD(4) * pkin(4) - t479;
t392 = -t324 * t372 + t416;
t391 = t324 * t375 + t328 * t425;
t390 = -pkin(8) * t290 + t274 * t484;
t385 = -t396 + t355;
t248 = pkin(5) * t267 + qJDD(6) + t257;
t272 = -qJD(3) * t328 + qJD(4) * t296;
t379 = (t453 + t457) * MDP(25) + t401 * MDP(26) + (-t375 * MDP(23) + t372 * MDP(24)) * t484;
t378 = qJD(1) ^ 2;
t357 = t377 * qJ(2);
t335 = t369 * t375;
t334 = t369 * t372;
t318 = t350 * t443 - t446;
t316 = t350 * t444 + t445;
t299 = t300 ^ 2;
t284 = t375 * t291;
t281 = t375 * t288;
t264 = pkin(5) * t300 + qJD(6) + t274;
t260 = qJ(6) * t450 + t438;
t254 = pkin(5) * t329 + qJ(6) * t449 - t296 * t372 + t284;
t247 = qJ(6) * t416 + (-qJD(5) * t296 + t395) * t372 + t419;
t246 = pkin(5) * t325 - t271 * t372 + t281 + t395 * t375 + (-t294 + (-qJ(6) * t328 - t291) * t372) * qJD(5);
t1 = [(qJD(2) * t323 - t296 * qJDD(4) + t333 * t324 - t327 * t328 - t481 - t343 * t386 + (-t343 * t474 - t271) * qJD(4)) * MDP(17) + (-t246 * t302 - t247 * t300 + t254 * t266 - t260 * t267 + t481 + t401 * t324 + (qJD(5) * t400 + t244 * t375 + t245 * t372) * t328) * MDP(25) + (t477 + t432 * (-t330 - t471)) * MDP(9) + (t336 * qJ(2) + t349 * qJD(2) - g(1) * (t371 * t374 + t357) - g(2) * (qJ(3) * t377 + t434) + (-qJD(3) * t338 + t330 * t371) * t432) * MDP(10) + (-t266 * t329 - t283 * t328 + t302 * t325 + t391 * t484) * MDP(20) + ((-t454 - t456) * t324 + (-t461 + t460 + (t453 - t457) * qJD(5)) * t328) * MDP(19) + t477 * MDP(2) + t437 * MDP(13) + ((-t296 * t424 + t281) * t484 + t284 * t290 + t403 * t329 + t258 * t325 + t272 * t300 + t295 * t267 - t274 * t416 - g(1) * t318 - g(2) * t316 + ((-qJD(5) * t291 - t271) * t484 - t296 * t290 + t406 * t329 - t257 * t328 + t274 * t324) * t372) * MDP(23) + (t413 - 0.2e1 * t462) * MDP(4) + (t328 * t293 - t323 * t325 - t324 * t474 - t329 * t380) * MDP(12) + (qJD(2) * t474 - qJD(4) * t272 - qJDD(4) * t295 + t293 * t343 + t325 * t333 + t327 * t329 - t350 * t405) * MDP(16) + (-(qJDD(2) - t462) * pkin(1) - g(1) * (-pkin(1) * t374 + t357) - t414 + (t364 + t469) * qJ(2)) * MDP(6) + (0.2e1 * t364 + t469 - t405) * MDP(5) + (t323 * t324 - t328 * t380) * MDP(11) + (-(-t296 * t425 + t419) * t484 - t438 * t290 - t388 * t329 - t259 * t325 + t272 * t302 - t295 * t266 - t257 * t449 + g(1) * t317 - g(2) * t315 + t391 * t274) * MDP(24) + (t266 * t449 + t302 * t391) * MDP(18) + (t245 * t260 + t251 * t247 + t244 * t254 + t249 * t246 + t248 * (-pkin(5) * t450 + t295) + t264 * (-pkin(5) * t392 + t272) - g(1) * t357 - t414 + (-g(1) * t385 - g(2) * t412) * t377 + (-g(1) * (-pkin(1) - t412) - g(2) * t385) * t374) * MDP(26) + (-t267 * t329 + t290 * t450 - t300 * t325 + t392 * t484) * MDP(21) + qJDD(1) * MDP(1) + t405 * MDP(3) + (t290 * t329 + t325 * t484) * MDP(22) + (-qJD(4) * t325 - qJDD(4) * t329) * MDP(14) + t402 * (t476 + t478); t413 * MDP(6) + (t330 * t432 - t477) * MDP(10) + t437 * MDP(16) + (t267 * t328 - t300 * t324) * MDP(23) + (-t266 * t328 - t302 * t324) * MDP(24) + (t248 * t328 - t264 * t324 - t477) * MDP(26) + (-qJD(4) * MDP(17) + (t454 - t456) * MDP(25) + t400 * MDP(26) + t394 * t484) * t325 + (-t349 * MDP(10) - MDP(16) * t474 - t323 * MDP(17) + t379) * qJD(1) + (-MDP(6) * qJ(2) - MDP(5) - t402) * t378 + (-pkin(1) * MDP(6) + MDP(4) - t409) * qJDD(1) + (-qJDD(4) * MDP(17) + (-t460 - t461) * MDP(25) + (-t244 * t372 + t245 * t375) * MDP(26) + t394 * t290 + t379 * qJD(5)) * t329; (qJD(1) * t338 * t432 + t478) * MDP(10) - t337 * MDP(16) + (t439 - t458) * MDP(23) - MDP(24) * t455 + t441 * MDP(25) + (-t264 * t323 - t405) * MDP(26) - t378 * t409 + (-MDP(23) * t426 - t290 * MDP(24) + MDP(25) * t489 + MDP(26) * t486) * t372 + ((t266 - t459) * MDP(25) + t487 * MDP(26) - t484 ^ 2 * MDP(24)) * t375 + ((t323 + t418) * MDP(16) + (-t367 * t430 - t368 * t431 - t474) * MDP(17)) * qJD(4) + (MDP(16) * t329 - MDP(17) * t328 + t402) * qJDD(1); t323 * t474 * MDP(11) + (t323 ^ 2 - t474 ^ 2) * MDP(12) - t386 * MDP(13) + ((t323 - t418) * qJD(4) + t472) * MDP(14) + qJDD(4) * MDP(15) + (qJD(4) * t278 - t323 * t333 - t383 + t393) * MDP(16) + (t474 * t333 - t399 + t473) * MDP(17) + (t302 * t407 - t461) * MDP(18) + ((-t266 - t459) * t375 - t484 * t454 + t441) * MDP(19) + (t290 * t372 + t407 * t484 - t455) * MDP(20) + (-t425 * t484 + t439 + t458) * MDP(21) - t484 * t323 * MDP(22) + (-pkin(4) * t267 - t258 * t323 - t278 * t300 - t286 * t484 + (t479 * t484 + t390) * t372 - t488 * t375) * MDP(23) + (pkin(4) * t266 + t259 * t323 - t278 * t302 + t372 * t488 + t390 * t375 + t440 * t484) * MDP(24) + (t266 * t334 + t267 * t335 - t436 * t300 - t435 * t302 - t372 * t487 + t486 * t375 - t473) * MDP(25) + (-t245 * t335 + t244 * t334 - t248 * t348 - g(3) * t396 + (t468 * t484 - t278) * t264 + t436 * t251 + t435 * t249 - t477 * (t348 * t351 - t350 * t369)) * MDP(26); t302 * t300 * MDP(18) + (-t299 + t470) * MDP(19) + (t300 * t484 - t266) * MDP(20) + (-t267 + t489) * MDP(21) + t290 * MDP(22) + (t259 * t484 - t274 * t302 + (t406 + t464) * t372 + t403 + t475) * MDP(23) + (g(1) * t316 - g(2) * t318 + t258 * t484 + t274 * t300 + t375 * t464 - t388) * MDP(24) + (pkin(5) * t266 - t300 * t442) * MDP(25) + (t442 * t251 + (-t264 * t302 + t372 * t464 + t244 + t475) * pkin(5)) * MDP(26); (-t299 - t470) * MDP(25) + (t249 * t302 + t251 * t300 + t248 + t383) * MDP(26);];
tau  = t1;
