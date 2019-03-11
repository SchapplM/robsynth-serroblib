% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:40:12
% EndTime: 2019-03-09 01:40:21
% DurationCPUTime: 6.28s
% Computational Cost: add. (4134->430), mult. (8966->557), div. (0->0), fcn. (6835->18), ass. (0->189)
t432 = cos(pkin(10));
t518 = cos(qJ(4));
t479 = t518 * t432;
t402 = qJD(1) * t479;
t429 = sin(pkin(10));
t436 = sin(qJ(4));
t495 = t436 * t429;
t478 = qJD(1) * t495;
t378 = -t402 + t478;
t374 = qJD(6) + t378;
t392 = t429 * t518 + t436 * t432;
t380 = t392 * qJD(1);
t428 = sin(pkin(11));
t431 = cos(pkin(11));
t362 = -t431 * qJD(4) + t380 * t428;
t438 = cos(qJ(6));
t364 = qJD(4) * t428 + t380 * t431;
t435 = sin(qJ(6));
t503 = t364 * t435;
t524 = -t438 * t362 - t503;
t526 = t374 * t524;
t430 = sin(pkin(9));
t403 = pkin(1) * t430 + qJ(3);
t387 = qJD(1) * qJD(3) + qJDD(1) * t403;
t412 = t432 * qJDD(2);
t352 = t412 + (-pkin(7) * qJDD(1) - t387) * t429;
t369 = t429 * qJDD(2) + t432 * t387;
t484 = qJDD(1) * t432;
t353 = pkin(7) * t484 + t369;
t397 = t403 * qJD(1);
t415 = t432 * qJD(2);
t509 = pkin(7) * qJD(1);
t366 = t415 + (-t397 - t509) * t429;
t372 = t429 * qJD(2) + t432 * t397;
t367 = t432 * t509 + t372;
t313 = t436 * t366 + t518 * t367;
t523 = t313 * qJD(4);
t449 = -t352 * t518 + t436 * t353 + t523;
t281 = -qJDD(4) * pkin(4) + qJDD(5) + t449;
t426 = pkin(10) + qJ(4);
t417 = sin(t426);
t420 = cos(t426);
t427 = qJ(1) + pkin(9);
t418 = sin(t427);
t421 = cos(t427);
t473 = g(1) * t421 + g(2) * t418;
t447 = -g(3) * t420 + t417 * t473;
t442 = t281 - t447;
t433 = cos(pkin(9));
t517 = pkin(1) * t433;
t408 = -pkin(2) - t517;
t485 = qJDD(1) * t408;
t393 = qJDD(3) + t485;
t477 = -g(1) * t418 + g(2) * t421;
t525 = -t393 - t477;
t460 = t362 * t435 - t364 * t438;
t522 = t374 * t460;
t391 = t428 * t438 + t431 * t435;
t382 = t391 * qJD(6);
t491 = t391 * t378 + t382;
t520 = t477 * t417;
t456 = t477 * t420;
t454 = -t518 * t366 + t436 * t367;
t384 = t392 * qJD(4);
t476 = qJDD(1) * t518;
t483 = qJDD(1) * t436;
t469 = t429 * t483 - t432 * t476;
t344 = qJD(1) * t384 + t469;
t342 = qJDD(6) + t344;
t389 = t428 * t435 - t438 * t431;
t492 = t374 * t389;
t519 = -t342 * t391 + t374 * t492;
t513 = g(3) * t417;
t445 = -t473 * t420 - t513;
t377 = t378 ^ 2;
t516 = pkin(8) * t431;
t437 = sin(qJ(1));
t514 = g(1) * t437;
t511 = pkin(7) + t403;
t510 = pkin(8) + qJ(5);
t508 = t524 * t380;
t507 = t460 * t380;
t505 = t344 * t428;
t504 = t344 * t431;
t502 = t378 * t428;
t452 = t479 - t495;
t383 = t452 * qJD(4);
t501 = t383 * t428;
t500 = t392 * t428;
t499 = t392 * t431;
t498 = t418 * t420;
t497 = t420 * t421;
t496 = t429 * MDP(6);
t481 = qJD(4) * t402 + t429 * t476 + t432 * t483;
t343 = -qJD(4) * t478 + t481;
t331 = -t431 * qJDD(4) + t343 * t428;
t332 = qJDD(4) * t428 + t343 * t431;
t488 = qJD(6) * t438;
t482 = -t435 * t331 + t438 * t332 - t362 * t488;
t489 = qJD(6) * t435;
t273 = -t364 * t489 + t482;
t494 = -t273 * t452 - t384 * t460;
t455 = t436 * t352 + t353 * t518;
t280 = qJDD(4) * qJ(5) + (qJD(5) - t454) * qJD(4) + t455;
t407 = pkin(3) * t432 + pkin(2);
t394 = -t407 - t517;
t373 = qJDD(1) * t394 + qJDD(3);
t290 = pkin(4) * t344 - qJ(5) * t343 - qJD(5) * t380 + t373;
t269 = t431 * t280 + t428 * t290;
t297 = t383 * t391 + t488 * t499 - t489 * t500;
t337 = t391 * t392;
t493 = -t297 * t374 - t337 * t342;
t307 = qJD(4) * qJ(5) + t313;
t376 = qJD(1) * t394 + qJD(3);
t321 = pkin(4) * t378 - qJ(5) * t380 + t376;
t284 = t431 * t307 + t428 * t321;
t341 = pkin(4) * t380 + qJ(5) * t378;
t294 = t428 * t341 - t431 * t454;
t385 = t511 * t429;
t386 = t511 * t432;
t453 = -t385 * t518 - t436 * t386;
t314 = qJD(3) * t452 + qJD(4) * t453;
t322 = pkin(4) * t384 - qJ(5) * t383 - qJD(5) * t392;
t287 = t431 * t314 + t428 * t322;
t334 = -pkin(4) * t452 - qJ(5) * t392 + t394;
t340 = -t436 * t385 + t386 * t518;
t299 = t428 * t334 + t431 * t340;
t490 = t429 ^ 2 + t432 ^ 2;
t305 = -qJD(4) * pkin(4) + qJD(5) + t454;
t487 = -qJD(5) + t305;
t268 = -t280 * t428 + t431 * t290;
t264 = pkin(5) * t344 - pkin(8) * t332 + t268;
t267 = -pkin(8) * t331 + t269;
t475 = t438 * t264 - t267 * t435;
t283 = -t307 * t428 + t431 * t321;
t293 = t431 * t341 + t428 * t454;
t286 = -t314 * t428 + t431 * t322;
t474 = t438 * t331 + t435 * t332;
t298 = t431 * t334 - t340 * t428;
t439 = cos(qJ(1));
t471 = -g(2) * t439 + t514;
t470 = -t389 * t342 - t374 * t491;
t468 = t264 * t435 + t267 * t438;
t467 = -t268 * t431 - t269 * t428;
t272 = pkin(5) * t378 - pkin(8) * t364 + t283;
t276 = -pkin(8) * t362 + t284;
t265 = t272 * t438 - t276 * t435;
t266 = t272 * t435 + t276 * t438;
t274 = -qJD(6) * t460 + t474;
t466 = t274 * t452 + t384 * t524;
t465 = -t283 * t428 + t284 * t431;
t285 = -pkin(5) * t452 - pkin(8) * t499 + t298;
t292 = -pkin(8) * t500 + t299;
t464 = t285 * t438 - t292 * t435;
t463 = t285 * t435 + t292 * t438;
t296 = -t382 * t392 - t383 * t389;
t338 = t389 * t392;
t462 = -t296 * t374 + t338 * t342;
t461 = -t344 * t392 - t378 * t383;
t368 = -t387 * t429 + t412;
t459 = -t368 * t429 + t369 * t432;
t458 = (-t397 * t429 + t415) * t429 - t372 * t432;
t457 = pkin(4) * t420 + qJ(5) * t417 + t407;
t396 = t510 * t431;
t451 = pkin(5) * t380 + qJD(5) * t428 + qJD(6) * t396 + t378 * t516 + t293;
t395 = t510 * t428;
t450 = pkin(8) * t502 - qJD(5) * t431 + qJD(6) * t395 + t294;
t444 = t281 * t392 + t305 * t383 - t473;
t315 = qJD(3) * t392 + qJD(4) * t340;
t434 = -pkin(7) - qJ(3);
t425 = pkin(11) + qJ(6);
t422 = t439 * pkin(1);
t419 = cos(t425);
t416 = sin(t425);
t409 = -pkin(5) * t431 - pkin(4);
t361 = t416 * t418 + t419 * t497;
t360 = -t416 * t497 + t418 * t419;
t359 = t416 * t421 - t419 * t498;
t358 = t416 * t498 + t419 * t421;
t346 = -qJD(4) * t384 + qJDD(4) * t452;
t345 = qJD(4) * t383 + qJDD(4) * t392;
t318 = pkin(5) * t500 - t453;
t301 = pkin(5) * t501 + t315;
t300 = -pkin(5) * t502 + t313;
t295 = t362 * pkin(5) + t305;
t279 = -pkin(8) * t501 + t287;
t275 = pkin(5) * t384 - t383 * t516 + t286;
t270 = t331 * pkin(5) + t281;
t1 = [(t387 * t490 + t459 - t473) * MDP(7) + (-t273 * t337 + t274 * t338 + t296 * t524 + t297 * t460) * MDP(21) + (pkin(1) * t514 - g(2) * t422 + t268 * t298 + t269 * t299 - t281 * t453 + t283 * t286 + t284 * t287 + t305 * t315 + (g(1) * t434 - g(2) * t457) * t421 + (g(1) * t457 + g(2) * t434) * t418) * MDP(19) + qJDD(1) * MDP(1) + (t269 * t452 - t284 * t384 - t287 * t378 - t299 * t344 + t315 * t364 - t332 * t453 + t428 * t456 + t431 * t444) * MDP(17) + (-qJD(4) * t314 - qJDD(4) * t340 + t343 * t394 + t373 * t392 + t376 * t383 + t520) * MDP(15) + (-t286 * t364 - t287 * t362 - t298 * t332 - t299 * t331 - t520 + t467 * t392 + (-t283 * t431 - t284 * t428) * t383) * MDP(18) + (-t268 * t452 + t283 * t384 + t286 * t378 + t298 * t344 + t315 * t362 - t331 * t453 + t428 * t444 - t431 * t456) * MDP(16) + (-qJD(4) * t315 + qJDD(4) * t453 + t344 * t394 - t373 * t452 + t376 * t384 - t456) * MDP(14) + (-(t275 * t435 + t279 * t438) * t374 - t463 * t342 + t468 * t452 - t266 * t384 - t301 * t460 + t318 * t273 - t270 * t338 + t295 * t296 - g(1) * t358 - g(2) * t360 + (t265 * t452 - t374 * t464) * qJD(6)) * MDP(26) + ((t275 * t438 - t279 * t435) * t374 + t464 * t342 - t475 * t452 + t265 * t384 - t301 * t524 + t318 * t274 + t270 * t337 + t295 * t297 - g(1) * t359 - g(2) * t361 + (t266 * t452 - t374 * t463) * qJD(6)) * MDP(25) + (t343 * t452 - t380 * t384 + t461) * MDP(10) + (-t342 * t452 + t374 * t384) * MDP(24) + (-t273 * t338 - t296 * t460) * MDP(20) + (t466 + t493) * MDP(23) + (-t462 + t494) * MDP(22) + (t471 + (t430 ^ 2 + t433 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + t471 * MDP(2) + (t393 * t408 - g(1) * (-pkin(1) * t437 - pkin(2) * t418 + qJ(3) * t421) - g(2) * (pkin(2) * t421 + qJ(3) * t418 + t422) + t459 * t403 - t458 * qJD(3)) * MDP(8) + t345 * MDP(11) + t346 * MDP(12) + (t343 * t392 + t380 * t383) * MDP(9) + (g(1) * t439 + g(2) * t437) * MDP(3) + (t432 * MDP(5) - t496) * (-t485 + t525); (qJDD(2) - g(3)) * MDP(4) + (t368 * t432 + t369 * t429 - g(3)) * MDP(8) + t346 * MDP(14) - t345 * MDP(15) + (-t331 * t452 + t362 * t384) * MDP(16) + (-t332 * t452 + t364 * t384) * MDP(17) + (-t281 * t452 + t305 * t384 - g(3)) * MDP(19) + (-t466 + t493) * MDP(25) + (t462 + t494) * MDP(26) + (t461 * MDP(17) + (-t331 * t392 - t362 * t383) * MDP(18) + (t269 * t392 + t284 * t383) * MDP(19)) * t431 + (t461 * MDP(16) + (t332 * t392 + t364 * t383) * MDP(18) + (-t268 * t392 - t283 * t383) * MDP(19)) * t428; -MDP(5) * t484 + qJDD(1) * t496 - t490 * MDP(7) * qJD(1) ^ 2 + (qJD(1) * t458 - t525) * MDP(8) + (0.2e1 * qJD(4) * t380 + t469) * MDP(14) + ((-t378 - t478) * qJD(4) + t481) * MDP(15) + (-t362 * t380 - t377 * t428 + t504) * MDP(16) + (-t364 * t380 - t377 * t431 - t505) * MDP(17) + (-t331 * t428 - t332 * t431 + (-t362 * t431 + t364 * t428) * t378) * MDP(18) + (-t305 * t380 + t378 * t465 - t467 + t477) * MDP(19) + (t470 + t508) * MDP(25) + (t507 + t519) * MDP(26); -t377 * MDP(10) + ((t378 - t478) * qJD(4) + t481) * MDP(11) - t469 * MDP(12) + qJDD(4) * MDP(13) + (t447 - t449 + t523) * MDP(14) + (t376 * t378 - t445 - t455) * MDP(15) + (-qJ(5) * t505 - pkin(4) * t331 - t313 * t362 + (t428 * t487 - t293) * t378 - t442 * t431) * MDP(16) + (-qJ(5) * t504 - pkin(4) * t332 - t313 * t364 + (t431 * t487 + t294) * t378 + t442 * t428) * MDP(17) + (t293 * t364 + t294 * t362 + (-qJ(5) * t331 - qJD(5) * t362 - t283 * t378 + t269) * t431 + (qJ(5) * t332 + qJD(5) * t364 - t284 * t378 - t268) * t428 + t445) * MDP(18) + (-t283 * t293 - t284 * t294 - t305 * t313 + t465 * qJD(5) - t442 * pkin(4) + (-t268 * t428 + t269 * t431 + t445) * qJ(5)) * MDP(19) + (t273 * t391 + t460 * t492) * MDP(20) + (-t273 * t389 - t274 * t391 + t460 * t491 - t492 * t524) * MDP(21) + (t507 - t519) * MDP(22) + (t470 - t508) * MDP(23) + ((-t395 * t438 - t396 * t435) * t342 + t409 * t274 + t270 * t389 + t300 * t524 + (t435 * t450 - t438 * t451) * t374 + t491 * t295 + t447 * t419) * MDP(25) + (-(-t395 * t435 + t396 * t438) * t342 + t409 * t273 + t270 * t391 + t300 * t460 + (t435 * t451 + t438 * t450) * t374 - t492 * t295 - t447 * t416) * MDP(26) + (MDP(10) * t380 - MDP(14) * t376 - MDP(16) * t283 + MDP(17) * t284 - MDP(24) * t374 - MDP(25) * t265 + MDP(26) * t266 + MDP(9) * t378) * t380; (t364 * t378 + t331) * MDP(16) + (-t362 * t378 + t332) * MDP(17) + (-t362 ^ 2 - t364 ^ 2) * MDP(18) + (t283 * t364 + t284 * t362 + t442) * MDP(19) + (t274 - t522) * MDP(25) + (t273 + t526) * MDP(26); t460 * t524 * MDP(20) + (t460 ^ 2 - t524 ^ 2) * MDP(21) + (t482 - t526) * MDP(22) + (-t474 - t522) * MDP(23) + t342 * MDP(24) + (-g(1) * t360 + g(2) * t358 + t266 * t374 + t295 * t460 + t416 * t513 + t475) * MDP(25) + (g(1) * t361 - g(2) * t359 + t265 * t374 - t295 * t524 + t419 * t513 - t468) * MDP(26) + (-MDP(22) * t503 + MDP(23) * t460 - MDP(25) * t266 - MDP(26) * t265) * qJD(6);];
tau  = t1;
