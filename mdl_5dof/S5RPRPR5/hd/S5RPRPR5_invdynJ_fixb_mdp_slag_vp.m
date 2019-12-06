% Calculate vector of inverse dynamics joint torques for
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:57:25
% EndTime: 2019-12-05 17:57:32
% DurationCPUTime: 3.57s
% Computational Cost: add. (2387->351), mult. (5936->473), div. (0->0), fcn. (4375->12), ass. (0->195)
t400 = sin(pkin(8));
t395 = t400 ^ 2;
t511 = 0.2e1 * t395;
t408 = cos(qJ(3));
t467 = qJD(3) * t408;
t519 = qJ(2) * t467;
t399 = sin(pkin(9));
t401 = cos(pkin(9));
t405 = sin(qJ(3));
t361 = t399 * t408 + t401 * t405;
t423 = qJD(1) * t361;
t335 = t400 * t423;
t407 = cos(qJ(5));
t327 = t407 * t335;
t470 = qJD(1) * t408;
t449 = t400 * t470;
t472 = qJD(1) * t400;
t450 = t405 * t472;
t338 = -t399 * t450 + t401 * t449;
t404 = sin(qJ(5));
t495 = t338 * t404;
t295 = t327 + t495;
t402 = cos(pkin(8));
t471 = qJD(1) * t402;
t381 = -qJD(3) + t471;
t375 = -qJD(5) + t381;
t496 = t295 * t375;
t364 = -pkin(2) * t402 - pkin(6) * t400 - pkin(1);
t348 = t364 * qJDD(1) + qJDD(2);
t344 = t408 * t348;
t457 = qJDD(1) * t402;
t379 = -qJDD(3) + t457;
t499 = qJ(2) * t408;
t380 = t402 * t499;
t466 = qJD(4) * t400;
t469 = qJD(2) * t402;
t418 = -t405 * t469 - t408 * t466;
t491 = t400 * t408;
t452 = qJ(4) * t491;
t500 = qJ(2) * t405;
t454 = t402 * t500;
t420 = -t452 - t454;
t492 = t400 * t405;
t453 = qJ(4) * t492;
t349 = t364 * qJD(1) + qJD(2);
t468 = qJD(3) * t349;
t273 = -t405 * t468 - pkin(3) * t379 + t344 + t420 * qJDD(1) + ((-t380 + t453) * qJD(3) + t418) * qJD(1);
t412 = t420 * qJD(3) - t405 * t466;
t461 = qJD(1) * qJD(2);
t448 = t408 * t461;
t440 = qJDD(1) * t380 + t405 * t348 + t349 * t467 + t402 * t448;
t456 = qJDD(1) * t405;
t446 = t400 * t456;
t279 = -qJ(4) * t446 + t412 * qJD(1) + t440;
t262 = t401 * t273 - t279 * t399;
t422 = t361 * qJD(3);
t455 = qJDD(1) * t408;
t308 = (qJD(1) * t422 + t399 * t456 - t401 * t455) * t400;
t260 = -pkin(4) * t379 + pkin(7) * t308 + t262;
t345 = t408 * t349;
t316 = t420 * qJD(1) + t345;
t306 = -pkin(3) * t381 + t316;
t317 = -qJ(4) * t450 + qJD(1) * t380 + t405 * t349;
t490 = t401 * t317;
t284 = t399 * t306 + t490;
t508 = pkin(7) * t335;
t270 = t284 - t508;
t465 = qJD(5) * t404;
t518 = t404 * t260 - t270 * t465;
t263 = t399 * t273 + t401 * t279;
t430 = t399 * t405 - t401 * t408;
t421 = t430 * qJD(3);
t307 = (qJD(1) * t421 - t361 * qJDD(1)) * t400;
t261 = pkin(7) * t307 + t263;
t350 = pkin(3) * t450 + qJ(2) * t472 + qJD(4);
t315 = pkin(4) * t335 + t350;
t393 = qJ(3) + pkin(9) + qJ(5);
t386 = sin(t393);
t409 = cos(qJ(1));
t387 = cos(t393);
t406 = sin(qJ(1));
t486 = t406 * t387;
t331 = -t386 * t409 + t402 * t486;
t489 = t402 * t409;
t494 = t386 * t406;
t333 = -t387 * t489 - t494;
t506 = g(1) * t400;
t517 = -g(2) * t331 - g(3) * t333 - t407 * t261 + t315 * t295 + t387 * t506 - t518;
t370 = -qJDD(5) + t379;
t362 = t370 * MDP(21);
t432 = -t335 * t404 + t407 * t338;
t516 = t295 * t432 * MDP(17) + (-t295 ^ 2 + t432 ^ 2) * MDP(18) - t362;
t501 = g(3) * t409;
t459 = qJDD(1) * qJ(2);
t427 = t459 + t461;
t514 = t427 * t402;
t497 = t432 * t375;
t483 = t408 * t409;
t488 = t405 * t406;
t353 = t402 * t488 + t483;
t485 = t406 * t408;
t487 = t405 * t409;
t355 = t402 * t487 - t485;
t513 = -g(2) * t353 + g(3) * t355;
t330 = t387 * t409 + t402 * t494;
t332 = t386 * t489 - t486;
t444 = t407 * t260 - t404 * t261;
t512 = -g(2) * t330 + g(3) * t332 - t315 * t432 + t386 * t506 + t444;
t443 = t407 * t307 + t308 * t404;
t265 = t432 * qJD(5) - t443;
t396 = t402 ^ 2;
t510 = pkin(3) * t399;
t509 = pkin(3) * t405;
t507 = pkin(7) * t338;
t504 = g(2) * t406;
t502 = qJ(2) * t501;
t498 = qJDD(1) * pkin(1);
t410 = qJD(1) ^ 2;
t493 = t395 * t410;
t311 = t399 * t317;
t283 = t401 * t306 - t311;
t268 = -pkin(4) * t381 + t283 - t507;
t484 = t407 * t268;
t480 = t364 * t467 + t408 * t469;
t304 = t412 + t480;
t305 = (-t380 + (qJ(4) * t400 - t364) * t405) * qJD(3) + t418;
t277 = t401 * t304 + t399 * t305;
t287 = t401 * t316 - t311;
t359 = t408 * t364;
t321 = -t452 + t359 + (-pkin(3) - t500) * t402;
t479 = t405 * t364 + t380;
t325 = -t453 + t479;
t290 = t399 * t321 + t401 * t325;
t482 = -t402 * t423 + t422;
t481 = t430 * t471 - t421;
t377 = t400 * pkin(3) * t467;
t478 = t400 * qJD(2) + t377;
t477 = pkin(3) * t492 + t400 * qJ(2);
t476 = t395 + t396;
t398 = t408 ^ 2;
t475 = t405 ^ 2 - t398;
t474 = MDP(10) * t400;
t473 = MDP(11) * t400;
t463 = t379 * MDP(12);
t462 = qJD(3) + t381;
t460 = qJD(1) * qJD(3);
t458 = qJDD(1) * t400;
t451 = -qJD(5) * t327 + t404 * t307 - t407 * t308;
t447 = t405 * t460;
t445 = t476 * t410;
t276 = -t304 * t399 + t401 * t305;
t286 = -t316 * t399 - t490;
t289 = t401 * t321 - t325 * t399;
t442 = qJD(1) * t462;
t441 = t379 + t457;
t439 = 0.2e1 * t476;
t438 = qJD(3) * t454;
t437 = g(2) * t409 + g(3) * t406;
t436 = -t501 + t504;
t435 = qJD(5) * t361 + t482;
t434 = -qJD(5) * t430 + t481;
t433 = -t404 * t268 - t407 * t270;
t346 = t361 * t400;
t347 = t430 * t400;
t431 = -t407 * t346 + t347 * t404;
t310 = -t346 * t404 - t347 * t407;
t429 = pkin(3) * t446 + qJ(2) * t458 + qJD(1) * t377 + t400 * t461 + qJDD(4);
t428 = qJD(3) * (t381 + t471);
t388 = pkin(3) * t401 + pkin(4);
t426 = t388 * t404 + t407 * t510;
t425 = t388 * t407 - t404 * t510;
t424 = (pkin(3) * t408 + pkin(2)) * t402 - t400 * (-qJ(4) - pkin(6)) + pkin(1);
t264 = -t338 * t465 + t451;
t419 = -t437 - t498;
t390 = qJDD(2) - t498;
t417 = -t390 - t419;
t416 = -t381 ^ 2 - t493;
t414 = t439 * t461 + t504;
t356 = -t402 * t483 - t488;
t354 = t402 * t485 - t487;
t341 = t400 * t422;
t337 = t400 * t421;
t324 = pkin(3) * t449 + pkin(4) * t338;
t322 = pkin(4) * t346 + t477;
t318 = -pkin(4) * t337 + t478;
t288 = -pkin(4) * t307 + t429;
t285 = -pkin(7) * t346 + t290;
t282 = -pkin(4) * t402 + pkin(7) * t347 + t289;
t281 = t310 * qJD(5) - t407 * t337 - t341 * t404;
t280 = t431 * qJD(5) + t337 * t404 - t341 * t407;
t275 = t287 - t507;
t274 = t286 + t508;
t267 = pkin(7) * t337 + t277;
t266 = pkin(7) * t341 + t276;
t1 = [qJDD(1) * MDP(1) + t437 * MDP(2) - t436 * MDP(3) - t417 * t400 * MDP(5) + (t439 * t459 + t414 - t501) * MDP(6) + (-t502 + (-t390 + t437) * pkin(1) + (t476 * t459 + t414) * qJ(2)) * MDP(7) + (qJDD(1) * t398 - 0.2e1 * t408 * t447) * t395 * MDP(8) + (-t405 * t455 + t475 * t460) * MDP(9) * t511 + (t405 * t428 - t441 * t408) * t474 + (t441 * t405 + t408 * t428) * t473 + (-g(2) * t356 + g(3) * t354 - t359 * t379 + (t511 + t396) * qJD(1) * t519 + (qJD(3) * t364 * t381 + t427 * t511) * t405) * MDP(13) + ((-t438 + t480) * t381 + t479 * t379 - g(2) * t355 - g(3) * t353 + (t448 + (-t447 + t455) * qJ(2)) * t511) * MDP(14) + (t262 * t347 - t263 * t346 - t276 * t338 - t277 * t335 + t283 * t341 + t284 * t337 + t289 * t308 + t290 * t307 + t437 * t400) * MDP(15) + (t263 * t290 + t284 * t277 + t262 * t289 + t283 * t276 + t429 * t477 + t350 * t478 - t502 + (g(2) * t424 - g(3) * t509) * t409 + (-g(2) * (-qJ(2) - t509) + g(3) * t424) * t406) * MDP(16) + (t264 * t310 + t280 * t432) * MDP(17) + (t264 * t431 - t265 * t310 - t280 * t295 - t281 * t432) * MDP(18) + (-t280 * t375 - t310 * t370) * MDP(19) + (t281 * t375 - t370 * t431) * MDP(20) + (-(t282 * t407 - t285 * t404) * t370 + t318 * t295 + t322 * t265 - t288 * t431 + t315 * t281 - g(2) * t333 + g(3) * t331 + (-t266 * t407 + t267 * t404 - (-t282 * t404 - t285 * t407) * qJD(5)) * t375) * MDP(22) + (-g(2) * t332 - g(3) * t330 + t322 * t264 + t315 * t280 + t288 * t310 + t318 * t432 + ((-qJD(5) * t285 + t266) * t375 + t282 * t370) * t404 + ((qJD(5) * t282 + t267) * t375 + t285 * t370) * t407) * MDP(23) + (t417 * MDP(4) + t463 + (-t344 + t381 * t519 + (qJ(2) * t379 + qJD(2) * t381 + t468 + t514) * t405) * MDP(13) + (-qJD(1) * t438 + t440) * MDP(14) - t264 * MDP(19) + t265 * MDP(20) + t362 + (-t433 * qJD(5) - t444) * MDP(22) + ((qJD(5) * t268 + t261) * t407 + t518) * MDP(23)) * t402; -MDP(4) * t457 + MDP(5) * t458 - MDP(6) * t445 + (-qJ(2) * t445 + qJDD(2) + t419) * MDP(7) + (-t379 * t408 + t416 * t405) * MDP(13) + (t405 * t379 + t416 * t408) * MDP(14) + (t307 * t361 - t308 * t430 - t481 * t335 + t482 * t338) * MDP(15) + (-t262 * t430 + t263 * t361 - t482 * t283 + t481 * t284 - t350 * t472 - t437) * MDP(16) + (-(-t361 * t404 - t407 * t430) * t370 - t295 * t472 + (t434 * t404 + t435 * t407) * t375) * MDP(22) + ((t361 * t407 - t404 * t430) * t370 - t432 * t472 + (-t435 * t404 + t434 * t407) * t375) * MDP(23); t408 * t405 * MDP(8) * t493 - t475 * MDP(9) * t493 + (-t405 * t442 + t455) * t474 + (-t408 * t442 - t456) * t473 - t463 + (t344 + (-t402 * t442 - t493) * t499 + (-t462 * t349 + t506 - t514) * t405 + t513) * MDP(13) + (g(1) * t491 - g(2) * t354 - g(3) * t356 - t345 * t381 + (t462 * t471 + t493) * t500 - t440) * MDP(14) + ((t284 + t286) * t338 - (t283 - t287) * t335 + (t307 * t399 + t308 * t401) * pkin(3)) * MDP(15) + (-t283 * t286 - t284 * t287 + (t263 * t399 + t262 * t401 + (g(1) * t405 - t350 * t470) * t400 + t513) * pkin(3)) * MDP(16) + (t264 - t496) * MDP(19) + (-t265 - t497) * MDP(20) + (-t425 * t370 + (t274 * t407 - t275 * t404) * t375 - t324 * t295 + (t426 * t375 + t433) * qJD(5) + t512) * MDP(22) + (t426 * t370 - (t274 * t404 + t275 * t407) * t375 - t324 * t432 + (t425 * t375 - t484) * qJD(5) + t517) * MDP(23) + t516; (-t335 ^ 2 - t338 ^ 2) * MDP(15) + (g(1) * t402 + t283 * t338 + t284 * t335 + t436 * t400 + t429) * MDP(16) + (t265 - t497) * MDP(22) + (t264 + t496) * MDP(23); (t451 - t496) * MDP(19) + (t443 - t497) * MDP(20) + (t433 * t375 + t512) * MDP(22) + (-(-t270 * t404 + t484) * t375 + t517) * MDP(23) + (-MDP(19) * t495 - t432 * MDP(20) + t433 * MDP(22) - MDP(23) * t484) * qJD(5) + t516;];
tau = t1;
