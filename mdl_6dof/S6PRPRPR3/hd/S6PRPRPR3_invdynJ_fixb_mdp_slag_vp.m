% Calculate vector of inverse dynamics joint torques for
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:37:48
% EndTime: 2019-03-08 19:37:54
% DurationCPUTime: 4.20s
% Computational Cost: add. (1961->402), mult. (4538->535), div. (0->0), fcn. (3774->12), ass. (0->193)
t405 = cos(qJ(2));
t396 = sin(pkin(6));
t485 = qJD(1) * t396;
t459 = t405 * t485;
t359 = qJD(2) * pkin(2) + t459;
t397 = cos(pkin(11));
t402 = sin(qJ(2));
t460 = t402 * t485;
t364 = t397 * t460;
t394 = sin(pkin(11));
t321 = t394 * t359 + t364;
t319 = qJD(2) * pkin(8) + t321;
t399 = cos(pkin(6));
t380 = qJD(1) * t399 + qJD(3);
t401 = sin(qJ(4));
t404 = cos(qJ(4));
t488 = -t401 * t319 + t404 * t380;
t530 = -qJD(5) + t488;
t484 = qJD(2) * t396;
t455 = qJD(1) * t484;
t471 = qJDD(1) * t396;
t533 = t402 * t471 + t405 * t455;
t294 = -qJD(4) * pkin(4) - t530;
t472 = qJD(2) * qJD(4);
t454 = t401 * t472;
t469 = qJDD(2) * t404;
t532 = t454 - t469;
t482 = qJD(2) * t404;
t453 = t404 * t472;
t470 = qJDD(2) * t401;
t429 = t453 + t470;
t348 = qJDD(6) + t429;
t403 = cos(qJ(6));
t341 = t403 * t348;
t400 = sin(qJ(6));
t351 = qJD(4) * t400 + t403 * t482;
t480 = qJD(4) * t351;
t531 = (t480 - t341) * MDP(22);
t466 = MDP(11) - MDP(14);
t465 = MDP(12) - MDP(15);
t327 = t394 * t459 + t364;
t479 = qJD(4) * t401;
t529 = -pkin(4) * t479 + t327;
t395 = sin(pkin(10));
t398 = cos(pkin(10));
t496 = t399 * t405;
t528 = -t395 * t496 - t398 * t402;
t301 = t404 * t319 + t401 * t380;
t297 = -qJD(4) * qJ(5) - t301;
t483 = qJD(2) * t401;
t473 = pkin(5) * t483 - t530;
t376 = t405 * t471;
t331 = qJDD(2) * pkin(2) - t402 * t455 + t376;
t296 = t394 * t331 + t397 * t533;
t291 = qJDD(2) * pkin(8) + t296;
t477 = qJD(4) * t404;
t526 = t401 * t291 + t319 * t477 + t380 * t479;
t494 = t405 * t397;
t349 = t394 * t402 - t494;
t430 = t349 * t399;
t443 = t394 * t405 + t397 * t402;
t304 = -t395 * t443 - t398 * t430;
t307 = t395 * t430 - t398 * t443;
t501 = t396 * t402;
t332 = t394 * t501 - t396 * t494;
t525 = -g(1) * t307 - g(2) * t304 + g(3) * t332;
t334 = t443 * t399;
t305 = t334 * t398 - t349 * t395;
t500 = t396 * t404;
t285 = t305 * t401 + t398 * t500;
t306 = t395 * t334 + t349 * t398;
t287 = -t306 * t401 - t395 * t500;
t333 = t443 * t396;
t316 = t333 * t401 - t399 * t404;
t524 = -g(1) * t287 - g(2) * t285 - g(3) * t316;
t387 = pkin(5) * t482;
t289 = -t297 + t387;
t382 = qJD(6) + t483;
t406 = -pkin(4) - pkin(9);
t523 = t406 * t348 + (t289 - t387 - t301) * t382;
t476 = qJD(5) * t401;
t432 = -qJ(5) * t477 - t476;
t295 = t331 * t397 - t394 * t533;
t437 = pkin(4) * t454 - t295;
t451 = -qJ(5) * t401 - pkin(3);
t438 = pkin(4) * t404 - t451;
t277 = qJD(2) * t432 - qJDD(2) * t438 + t437;
t521 = pkin(2) * t397;
t343 = -t438 - t521;
t385 = pkin(2) * t394 + pkin(8);
t407 = qJD(4) ^ 2;
t420 = t385 * t407 - t525;
t490 = t432 - t529;
t522 = qJD(2) * t490 + qJDD(2) * t343 + t277 + t420;
t514 = pkin(5) + t385;
t513 = qJ(5) * t404;
t512 = qJDD(4) * pkin(4);
t436 = t403 * qJDD(4) + t400 * t532;
t310 = -qJD(6) * t351 + t436;
t511 = t310 * t403;
t456 = t400 * t482;
t445 = -qJD(6) * t456 + qJDD(4) * t400 - t403 * t454;
t311 = (qJD(4) * qJD(6) + t469) * t403 + t445;
t510 = t311 * t401;
t509 = t348 * t400;
t508 = t351 * t382;
t478 = qJD(4) * t403;
t353 = -t456 + t478;
t507 = t353 * t382;
t377 = t399 * qJDD(1) + qJDD(3);
t506 = t377 * t404;
t505 = t382 * t400;
t504 = t382 * t403;
t503 = t395 * t402;
t502 = t396 * t401;
t499 = t396 * t405;
t497 = t399 * t402;
t495 = t400 * t404;
t493 = qJDD(1) - g(3);
t492 = t310 * t401 + t353 * t477;
t446 = pkin(9) * t401 - t513;
t418 = qJD(4) * t446 - t476;
t491 = -t418 + t529;
t475 = qJD(6) * t400;
t457 = t382 * t475;
t464 = t401 * t504;
t489 = qJD(4) * t464 + t404 * t457;
t392 = t401 ^ 2;
t393 = t404 ^ 2;
t487 = t392 - t393;
t474 = qJD(6) * t404;
t468 = qJDD(4) * qJ(5);
t467 = qJDD(4) * t385;
t462 = t398 * t496;
t408 = qJD(2) ^ 2;
t461 = t401 * t404 * t408;
t458 = t400 * t479;
t345 = t514 * t404;
t363 = t394 * t460;
t320 = t359 * t397 - t363;
t424 = qJDD(5) - t506 + t526;
t273 = t424 - t512;
t449 = -qJD(4) * t297 - t273;
t448 = -t404 * t291 + t319 * t479 - t401 * t377 - t380 * t477;
t280 = qJD(4) * t406 + t473;
t428 = t404 * t406 + t451;
t298 = qJD(2) * t428 - t320;
t274 = t280 * t403 - t298 * t400;
t275 = t280 * t400 + t298 * t403;
t281 = t316 * t403 - t332 * t400;
t282 = t316 * t400 + t332 * t403;
t317 = t333 * t404 + t399 * t401;
t433 = qJD(6) * t504 + t509;
t431 = -g(3) * t399 + (-g(1) * t395 + g(2) * t398) * t396;
t286 = t305 * t404 - t398 * t502;
t288 = -t306 * t404 + t395 * t502;
t427 = -g(1) * t288 - g(2) * t286 - g(3) * t317;
t426 = g(1) * t306 - g(2) * t305 - g(3) * t333;
t423 = -t403 * t474 + t458;
t276 = qJD(2) * t418 + qJDD(2) * t428 + t437;
t422 = -t276 + t525;
t344 = t514 * t401;
t421 = t344 * t348 + t426;
t419 = qJD(4) * qJD(5) - t448 + t468;
t318 = -qJD(2) * pkin(3) - t320;
t330 = t397 * t459 - t363;
t386 = -pkin(3) - t521;
t417 = -t467 + (qJD(2) * t386 + t318 + t330) * qJD(4);
t302 = -qJD(2) * t438 - t320;
t416 = t467 + (-qJD(2) * t343 - t302 - t330) * qJD(4);
t415 = t427 - t448;
t414 = -g(1) * t528 - g(3) * t499;
t413 = t419 * t404 + t273 * t401 + (t294 * t404 + t297 * t401) * qJD(4);
t412 = -qJD(4) * t301 + t524 + t526;
t271 = -pkin(5) * t532 + t419;
t389 = pkin(4) * t483;
t410 = t271 + (qJD(2) * t446 - qJD(6) * t406 + t389) * t382 + t427;
t409 = -qJD(2) * t327 - t295 + t420 + (-pkin(3) + t386) * qJDD(2);
t369 = pkin(2) * t462;
t368 = qJDD(4) * t404 - t401 * t407;
t367 = qJDD(4) * t401 + t404 * t407;
t354 = -qJ(5) * t482 + t389;
t339 = qJD(4) * t345;
t338 = t514 * t479;
t335 = t428 - t521;
t329 = t349 * t484;
t328 = qJD(2) * t333;
t299 = t302 * t483;
t279 = -t333 * t479 + (qJD(4) * t399 - t329) * t404;
t278 = qJD(4) * t317 - t329 * t401;
t270 = t429 * pkin(5) + qJDD(4) * t406 + t424;
t269 = t403 * t270;
t1 = [t493 * MDP(1) + (-t295 * t332 + t296 * t333 - t320 * t328 - t321 * t329 + t377 * t399 - g(3)) * MDP(5) + (t273 * t316 + t277 * t332 + t278 * t294 - t279 * t297 + t302 * t328 + t317 * t419 - g(3)) * MDP(16) + ((-qJD(6) * t282 + t278 * t403 - t328 * t400) * t382 + t281 * t348 + t279 * t351 + t317 * t311) * MDP(22) + (-(qJD(6) * t281 + t278 * t400 + t328 * t403) * t382 - t282 * t348 + t279 * t353 + t317 * t310) * MDP(23) + ((t316 * t401 + t317 * t404) * MDP(13) + (t401 * t465 - t404 * t466) * t332) * qJDD(2) + ((qJDD(2) * t405 - t402 * t408) * MDP(3) + (-qJDD(2) * t402 - t405 * t408) * MDP(4)) * t396 + ((t278 * t401 + t279 * t404 + t316 * t477 - t317 * t479) * MDP(13) + t465 * (t328 * t401 + t332 * t477) + t466 * (-t328 * t404 + t332 * t479)) * qJD(2) - t466 * (qJD(4) * t278 + qJDD(4) * t316) - t465 * (qJD(4) * t279 + qJDD(4) * t317); qJDD(2) * MDP(2) + (t376 - g(2) * (t462 - t503) + t414) * MDP(3) + (-g(1) * (t395 * t497 - t398 * t405) - g(2) * (-t395 * t405 - t398 * t497) - t493 * t501) * MDP(4) + (-g(2) * t369 + t320 * t327 - t321 * t330 + (g(2) * t503 + t295 * t397 + t296 * t394 + t414) * pkin(2)) * MDP(5) + (qJDD(2) * t392 + 0.2e1 * t401 * t453) * MDP(6) + 0.2e1 * (t401 * t469 - t472 * t487) * MDP(7) + t367 * MDP(8) + t368 * MDP(9) + (t401 * t417 - t404 * t409) * MDP(11) + (t401 * t409 + t404 * t417) * MDP(12) + (t413 + t426 + (-qJD(2) * t330 + qJDD(2) * t385) * (t392 + t393)) * MDP(13) + (t416 * t401 + t404 * t522) * MDP(14) + (-t401 * t522 + t416 * t404) * MDP(15) + (t277 * t343 - g(1) * (pkin(2) * t528 - pkin(8) * t306) - g(2) * (-pkin(2) * t503 + pkin(8) * t305 + t369) - g(3) * (pkin(2) * t499 + pkin(8) * t333) + (-t294 * t401 + t297 * t404) * t330 + t490 * t302 + t413 * t385 + t525 * t438) * MDP(16) + (-t310 * t495 + t353 * t423) * MDP(17) + ((-t351 * t400 + t353 * t403) * t479 + (-t511 + t311 * t400 + (t351 * t403 + t353 * t400) * qJD(6)) * t404) * MDP(18) + (-t348 * t495 + t382 * t423 + t492) * MDP(19) + (-t510 + (-t480 - t341) * t404 + t489) * MDP(20) + (t348 * t401 + t382 * t477) * MDP(21) + (-t335 * t509 + t345 * t311 - t338 * t351 + t421 * t403 + (-t289 * t478 + t400 * t422 + t269) * t401 + ((-t330 * t401 + t339) * t403 + t491 * t400) * t382 + ((-t335 * t403 - t344 * t400) * t382 - t275 * t401) * qJD(6) + (qJD(4) * t274 + t271 * t403 - t289 * t475 - t330 * t351) * t404) * MDP(22) + (-t275 * t477 + t345 * t310 + (-t330 * t404 - t338) * t353 + (-t289 * t474 - t335 * t348 + (-qJD(6) * t344 + t491) * t382 + (-qJD(6) * t280 + t422) * t401) * t403 + (-(-qJD(6) * t335 + t339) * t382 - t271 * t404 + (qJD(4) * t289 + qJD(6) * t298 + t330 * t382 - t270) * t401 - t421) * t400) * MDP(23); (t431 + t377) * MDP(5) + (t294 * t479 + t401 * t419 + t431) * MDP(16) + (t489 + t510) * MDP(22) + (-t382 * t458 + t492) * MDP(23) + t466 * t368 - t465 * t367 + (t449 * MDP(16) + t433 * MDP(23) + t531) * t404; -MDP(6) * t461 + t487 * t408 * MDP(7) + MDP(8) * t470 + MDP(9) * t469 + qJDD(4) * MDP(10) + (-t318 * t483 - t412 + t506) * MDP(11) + (qJD(4) * t488 - t318 * t482 - t415) * MDP(12) + (-pkin(4) * t401 + t513) * qJDD(2) * MDP(13) + (-0.2e1 * t512 + qJDD(5) + t299 + (-qJD(2) * t354 - t377) * t404 + t412) * MDP(14) + (0.2e1 * t468 + (0.2e1 * qJD(5) - t488) * qJD(4) + (t302 * t404 + t354 * t401) * qJD(2) + t415) * MDP(15) + (t419 * qJ(5) - t273 * pkin(4) - t302 * t354 - t294 * t301 - g(1) * (-pkin(4) * t287 + qJ(5) * t288) - g(2) * (-pkin(4) * t285 + qJ(5) * t286) - g(3) * (-pkin(4) * t316 + qJ(5) * t317) + t530 * t297) * MDP(16) + (-t353 * t505 + t511) * MDP(17) + ((-t311 - t507) * t403 + (-t310 + t508) * t400) * MDP(18) + (-t457 + t341 + (-t353 * t404 - t401 * t505) * qJD(2)) * MDP(19) + ((t351 * t404 - t464) * qJD(2) - t433) * MDP(20) - t382 * MDP(21) * t482 + (qJ(5) * t311 - t274 * t482 + t473 * t351 + t410 * t400 + t403 * t523) * MDP(22) + (qJ(5) * t310 + t275 * t482 + t473 * t353 - t400 * t523 + t410 * t403) * MDP(23); MDP(13) * t470 + (qJDD(4) + t461) * MDP(14) + (-t392 * t408 - t407) * MDP(15) + (t299 - t449 + t524) * MDP(16) - t531 + (-qJD(4) * t353 - t509) * MDP(23) + (-MDP(22) * t505 - MDP(23) * t504) * t382; t353 * t351 * MDP(17) + (-t351 ^ 2 + t353 ^ 2) * MDP(18) + (t436 + t508) * MDP(19) + (-t403 * t469 - t445 + t507) * MDP(20) + t348 * MDP(21) + (-t400 * t276 + t269 + t275 * t382 - t289 * t353 - g(1) * (t287 * t403 + t307 * t400) - g(2) * (t285 * t403 + t304 * t400) - g(3) * t281) * MDP(22) + (-t403 * t276 - t400 * t270 + t274 * t382 + t289 * t351 - g(1) * (-t287 * t400 + t307 * t403) - g(2) * (-t285 * t400 + t304 * t403) + g(3) * t282) * MDP(23) + (-MDP(19) * t351 - MDP(20) * t478 - MDP(22) * t275 - MDP(23) * t274) * qJD(6);];
tau  = t1;
