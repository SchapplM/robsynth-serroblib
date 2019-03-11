% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S6RPRPRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:03:12
% EndTime: 2019-03-09 03:03:17
% DurationCPUTime: 3.91s
% Computational Cost: add. (3963->409), mult. (8491->524), div. (0->0), fcn. (5856->14), ass. (0->191)
t398 = sin(pkin(9));
t376 = pkin(1) * t398 + pkin(7);
t475 = qJ(4) + t376;
t407 = cos(qJ(3));
t388 = t407 * qJDD(2);
t404 = sin(qJ(3));
t365 = t376 * qJDD(1);
t416 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t365;
t440 = t475 * qJD(1);
t427 = t440 * qJD(3);
t293 = qJDD(3) * pkin(3) - t404 * t416 - t407 * t427 + t388;
t297 = (qJDD(2) - t427) * t404 + t416 * t407;
t397 = sin(pkin(10));
t399 = cos(pkin(10));
t273 = t397 * t293 + t399 * t297;
t271 = qJDD(3) * pkin(8) + t273;
t338 = t407 * qJD(2) - t404 * t440;
t496 = qJD(3) * pkin(3);
t331 = t338 + t496;
t339 = qJD(2) * t404 + t407 * t440;
t479 = t399 * t339;
t296 = t397 * t331 + t479;
t290 = qJD(3) * pkin(8) + t296;
t400 = cos(pkin(9));
t378 = -t400 * pkin(1) - pkin(2);
t391 = t407 * pkin(3);
t513 = t378 - t391;
t350 = qJD(1) * t513 + qJD(4);
t461 = qJD(1) * t404;
t478 = t399 * t407;
t351 = qJD(1) * t478 - t397 * t461;
t361 = t397 * t407 + t399 * t404;
t353 = t361 * qJD(1);
t305 = -pkin(4) * t351 - pkin(8) * t353 + t350;
t403 = sin(qJ(5));
t406 = cos(qJ(5));
t276 = t290 * t406 + t305 * t403;
t454 = qJDD(1) * t407;
t455 = qJDD(1) * t404;
t433 = -t397 * t455 + t399 * t454;
t325 = -qJD(3) * t353 + t433;
t456 = qJD(1) * qJD(3);
t447 = t404 * t456;
t369 = t397 * t447;
t446 = t407 * t456;
t426 = t399 * t446 - t369;
t453 = pkin(3) * t447 + qJDD(4);
t516 = -pkin(8) * t361 + t513;
t281 = -t325 * pkin(4) - t426 * pkin(8) + qJDD(1) * t516 + t453;
t280 = t406 * t281;
t457 = t406 * qJD(3);
t459 = qJD(5) * t403;
t509 = qJDD(1) * t361 + t426;
t286 = -qJD(5) * t457 - t403 * qJDD(3) + t353 * t459 - t406 * t509;
t352 = t361 * qJD(3);
t324 = qJD(1) * t352 + qJDD(5) - t433;
t336 = qJD(3) * t403 + t353 * t406;
t260 = pkin(5) * t324 + qJ(6) * t286 - qJD(5) * t276 - qJD(6) * t336 - t271 * t403 + t280;
t334 = t353 * t403 - t457;
t267 = -qJ(6) * t334 + t276;
t349 = qJD(5) - t351;
t519 = t267 * t349 + t260;
t287 = qJD(5) * t336 - t406 * qJDD(3) + t403 * t509;
t458 = qJD(5) * t406;
t423 = t406 * t271 + t403 * t281 - t290 * t459 + t305 * t458;
t261 = -qJ(6) * t287 - qJD(6) * t334 + t423;
t275 = -t290 * t403 + t406 * t305;
t266 = -qJ(6) * t336 + t275;
t265 = pkin(5) * t349 + t266;
t518 = -t265 * t349 + t261;
t375 = pkin(3) * t397 + pkin(8);
t474 = qJ(6) + t375;
t517 = qJ(6) * t351 - qJD(5) * t474;
t394 = qJ(1) + pkin(9);
t383 = sin(t394);
t385 = cos(t394);
t444 = -g(1) * t383 + g(2) * t385;
t512 = qJDD(1) * t513 + t453;
t437 = g(1) * t385 + g(2) * t383;
t393 = qJ(3) + pkin(10);
t384 = cos(t393);
t480 = t385 * t406;
t483 = t383 * t403;
t344 = t384 * t483 + t480;
t481 = t385 * t403;
t482 = t383 * t406;
t346 = -t384 * t481 + t482;
t511 = -g(1) * t346 + g(2) * t344;
t272 = t293 * t399 - t397 * t297;
t270 = -qJDD(3) * pkin(4) - t272;
t382 = sin(t393);
t417 = g(3) * t384 - t382 * t437;
t460 = qJD(5) * t349;
t510 = t375 * t460 + t270 + t417;
t508 = t336 ^ 2;
t507 = pkin(3) * t399;
t405 = sin(qJ(1));
t502 = g(1) * t405;
t499 = g(3) * t382;
t497 = g(3) * t407;
t494 = t286 * t403;
t493 = t334 * t351;
t492 = t334 * t353;
t491 = t334 * t403;
t490 = t336 * t349;
t489 = t336 * t353;
t488 = t336 * t406;
t360 = t397 * t404 - t478;
t355 = t360 * qJD(3);
t487 = t355 * t403;
t486 = t355 * t406;
t485 = t361 * t403;
t484 = t361 * t406;
t328 = t397 * t339;
t477 = t403 * t324;
t476 = t403 * t349;
t357 = t475 * t404;
t359 = t475 * t407;
t322 = -t357 * t397 + t359 * t399;
t317 = t406 * t322;
t319 = t406 * t324;
t473 = qJDD(2) - g(3);
t472 = -t266 + t265;
t471 = -t287 * t484 + t334 * t486;
t470 = -t403 * t287 - t334 * t458;
t299 = t338 * t399 - t328;
t314 = pkin(3) * t461 + pkin(4) * t353 - pkin(8) * t351;
t469 = t406 * t299 + t403 * t314;
t468 = -t286 * t360 + t336 * t352;
t318 = pkin(4) * t360 + t516;
t467 = t403 * t318 + t317;
t466 = t351 * t476 + t319;
t465 = qJD(6) * t406 + t403 * t517 - t469;
t308 = t406 * t314;
t464 = -pkin(5) * t353 - t308 + t517 * t406 + (-qJD(6) + t299) * t403;
t381 = t391 + pkin(2);
t408 = cos(qJ(1));
t463 = t408 * pkin(1) + t385 * t381;
t395 = t404 ^ 2;
t462 = -t407 ^ 2 + t395;
t368 = qJD(1) * t378;
t452 = t404 * t496;
t451 = t336 * t487;
t443 = qJD(3) * t475;
t340 = qJD(4) * t407 - t404 * t443;
t341 = -qJD(4) * t404 - t407 * t443;
t304 = t340 * t399 + t341 * t397;
t315 = pkin(4) * t352 + pkin(8) * t355 + t452;
t450 = t406 * t304 + t403 * t315 + t318 * t458;
t380 = pkin(5) * t406 + pkin(4);
t449 = t361 * t458;
t402 = -qJ(4) - pkin(7);
t445 = pkin(5) * t403 - t402;
t295 = t331 * t399 - t328;
t298 = t397 * t338 + t479;
t303 = t340 * t397 - t399 * t341;
t321 = t399 * t357 + t359 * t397;
t441 = t349 * t406;
t439 = -qJD(5) * t305 - t271;
t435 = -g(2) * t408 + t502;
t434 = -t290 * t458 + t280;
t432 = -t265 * t406 - t267 * t403;
t431 = -t287 * t360 - t334 * t352;
t401 = -qJ(6) - pkin(8);
t430 = t380 * t384 - t382 * t401;
t429 = qJ(6) * t355 - qJD(6) * t361;
t289 = -qJD(3) * pkin(4) - t295;
t425 = t449 - t487;
t424 = -t361 * t459 - t486;
t422 = t289 * t349 - t375 * t324;
t419 = -qJD(1) * t368 - t365 + t437;
t418 = 0.2e1 * qJD(3) * t368 - qJDD(3) * t376;
t264 = pkin(5) * t287 + qJDD(6) + t270;
t409 = qJD(3) ^ 2;
t415 = -0.2e1 * qJDD(1) * t378 - t376 * t409 - t444;
t377 = -pkin(4) - t507;
t364 = qJDD(3) * t407 - t404 * t409;
t363 = qJDD(3) * t404 + t407 * t409;
t358 = t474 * t406;
t356 = t474 * t403;
t347 = t384 * t480 + t483;
t345 = -t384 * t482 + t481;
t333 = t334 ^ 2;
t313 = t406 * t318;
t309 = t406 * t315;
t282 = pkin(5) * t334 + qJD(6) + t289;
t278 = -qJ(6) * t485 + t467;
t277 = pkin(5) * t360 - qJ(6) * t484 - t322 * t403 + t313;
t263 = -qJ(6) * t449 + (-qJD(5) * t322 + t429) * t403 + t450;
t262 = pkin(5) * t352 - t304 * t403 + t309 + t429 * t406 + (-t317 + (qJ(6) * t361 - t318) * t403) * qJD(5);
t1 = [qJDD(1) * MDP(1) + t435 * MDP(2) + (g(1) * t408 + g(2) * t405) * MDP(3) + (t435 + (t398 ^ 2 + t400 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t395 + 0.2e1 * t404 * t446) * MDP(5) + 0.2e1 * (t404 * t454 - t456 * t462) * MDP(6) + t363 * MDP(7) + t364 * MDP(8) + (t404 * t418 + t407 * t415) * MDP(10) + (-t404 * t415 + t407 * t418) * MDP(11) + (-t272 * t361 - t273 * t360 + t295 * t355 - t296 * t352 + t303 * t353 + t304 * t351 + t321 * t509 + t322 * t325 - t437) * MDP(12) + (t273 * t322 + t296 * t304 - t272 * t321 - t295 * t303 + t512 * t513 + t350 * t452 - g(1) * (-pkin(1) * t405 - t381 * t383 - t385 * t402) - g(2) * (-t383 * t402 + t463)) * MDP(13) + (-t286 * t484 + t336 * t424) * MDP(14) + (t451 + (t494 + (-t488 + t491) * qJD(5)) * t361 + t471) * MDP(15) + (t319 * t361 + t349 * t424 + t468) * MDP(16) + (-t349 * t425 - t361 * t477 + t431) * MDP(17) + (t324 * t360 + t349 * t352) * MDP(18) + ((-t322 * t458 + t309) * t349 + t313 * t324 + t434 * t360 + t275 * t352 + t303 * t334 + t321 * t287 + t289 * t449 - g(1) * t345 - g(2) * t347 + ((-qJD(5) * t318 - t304) * t349 - t322 * t324 + t439 * t360 + t270 * t361 - t289 * t355) * t403) * MDP(19) + (-(-t322 * t459 + t450) * t349 - t467 * t324 - t423 * t360 - t276 * t352 + t303 * t336 - t321 * t286 + t270 * t484 - g(1) * t344 - g(2) * t346 + t424 * t289) * MDP(20) + (-t262 * t336 - t263 * t334 + t277 * t286 - t278 * t287 - t444 * t382 - t432 * t355 + (-t260 * t406 - t261 * t403 + (t265 * t403 - t267 * t406) * qJD(5)) * t361) * MDP(21) + (t261 * t278 + t267 * t263 + t260 * t277 + t265 * t262 + t264 * (pkin(5) * t485 + t321) + t282 * (pkin(5) * t425 + t303) + pkin(1) * t502 - g(2) * t463 + (-g(1) * t445 - g(2) * t430) * t385 + (-g(1) * (-t381 - t430) - g(2) * t445) * t383) * MDP(22); t473 * MDP(4) + t364 * MDP(10) - t363 * MDP(11) + (-t355 * t351 + t352 * t353 + t360 * t509) * MDP(12) + (-t272 * t360 - t295 * t352 - t296 * t355 - g(3)) * MDP(13) + (t355 * t476 - t431) * MDP(19) + (t349 * t486 + t468) * MDP(20) + (-t451 + t471) * MDP(21) + (t264 * t360 + t265 * t487 - t267 * t486 + t282 * t352 - g(3)) * MDP(22) + (t325 * MDP(12) + t273 * MDP(13) - MDP(21) * t494 + (-t260 * t403 + t261 * t406) * MDP(22) + (-t403 * MDP(19) - t406 * MDP(20)) * t324 + ((t488 + t491) * MDP(21) + t432 * MDP(22) + (-t406 * MDP(19) + t403 * MDP(20)) * t349) * qJD(5)) * t361; MDP(7) * t455 + MDP(8) * t454 + qJDD(3) * MDP(9) + (t404 * t419 + t388 - t497) * MDP(10) + (-t404 * t473 + t407 * t419) * MDP(11) + ((t296 - t298) * t353 + (-t299 + t295) * t351 + (t397 * t325 + (-t397 * t454 + t369 + (-t446 - t455) * t399) * t399) * pkin(3)) * MDP(12) + (t295 * t298 - t296 * t299 + (-t497 + t272 * t399 + t273 * t397 + (-qJD(1) * t350 + t437) * t404) * pkin(3)) * MDP(13) + (t336 * t441 - t494) * MDP(14) + ((-t286 + t493) * t406 - t336 * t476 + t470) * MDP(15) + (t349 * t441 + t477 - t489) * MDP(16) + (-t349 * t459 + t466 + t492) * MDP(17) - t349 * t353 * MDP(18) + (-t275 * t353 + t377 * t287 - t298 * t334 - t308 * t349 + (t299 * t349 + t422) * t403 - t510 * t406) * MDP(19) + (t276 * t353 - t377 * t286 - t298 * t336 + t469 * t349 + t403 * t510 + t422 * t406) * MDP(20) + (-t286 * t356 - t287 * t358 - t465 * t334 - t464 * t336 - t437 * t384 - t403 * t519 + t406 * t518 - t499) * MDP(21) + (t261 * t358 - t260 * t356 + t264 * (-t380 - t507) - g(3) * (t391 + t430) + (pkin(5) * t476 - t298) * t282 + t465 * t267 + t464 * t265 + t437 * (pkin(3) * t404 + t380 * t382 + t384 * t401)) * MDP(22) + (-MDP(5) * t404 * t407 + MDP(6) * t462) * qJD(1) ^ 2; (-t351 ^ 2 - t353 ^ 2) * MDP(12) + (t466 - t492) * MDP(19) - MDP(20) * t489 + t470 * MDP(21) + (-t282 * t353 + t444) * MDP(22) + (-MDP(19) * t460 - t324 * MDP(20) + MDP(21) * t490 + MDP(22) * t518) * t403 + ((t286 + t493) * MDP(21) + t519 * MDP(22) - t349 ^ 2 * MDP(20)) * t406 + (t295 * t353 - t296 * t351 + t444 + t512) * MDP(13); t336 * t334 * MDP(14) + (-t333 + t508) * MDP(15) + (t334 * t349 - t286) * MDP(16) + (-t287 + t490) * MDP(17) + t324 * MDP(18) + (t276 * t349 - t289 * t336 + (t439 + t499) * t403 + t434 + t511) * MDP(19) + (g(1) * t347 - g(2) * t345 + t275 * t349 + t289 * t334 + t406 * t499 - t423) * MDP(20) + (pkin(5) * t286 - t334 * t472) * MDP(21) + (t472 * t267 + (-t282 * t336 + t403 * t499 + t260 + t511) * pkin(5)) * MDP(22); (-t333 - t508) * MDP(21) + (t265 * t336 + t267 * t334 + t264 + t417) * MDP(22);];
tau  = t1;
