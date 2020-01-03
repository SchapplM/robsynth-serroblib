% Calculate vector of inverse dynamics joint torques for
% S5RRPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:06:03
% EndTime: 2020-01-03 12:06:08
% DurationCPUTime: 3.42s
% Computational Cost: add. (2533->341), mult. (3859->465), div. (0->0), fcn. (2535->14), ass. (0->193)
t388 = qJDD(1) + qJDD(2);
t391 = qJD(1) + qJD(2);
t396 = sin(pkin(9));
t402 = cos(qJ(5));
t403 = cos(qJ(4));
t468 = qJD(4) + qJD(5);
t440 = t468 * t403;
t398 = sin(qJ(5));
t495 = t398 * t403;
t399 = sin(qJ(4));
t507 = t388 * t399;
t474 = qJD(4) * t399;
t450 = t396 * t474;
t496 = t398 * t399;
t458 = t396 * t496;
t526 = -qJD(5) * t458 - t398 * t450;
t274 = t526 * t391 + (t388 * t495 + (t391 * t440 + t507) * t402) * t396;
t400 = sin(qJ(2));
t477 = qJD(2) * t400;
t463 = pkin(1) * t477;
t404 = cos(qJ(2));
t519 = pkin(1) * t404;
t485 = -qJD(1) * t463 + qJDD(1) * t519;
t447 = qJDD(3) - t485;
t518 = pkin(2) * t388;
t321 = t447 - t518;
t395 = qJ(1) + qJ(2);
t384 = sin(t395);
t386 = cos(t395);
t434 = g(2) * t386 + g(3) * t384;
t529 = t321 + t434;
t514 = pkin(1) * qJD(1);
t465 = t400 * t514;
t349 = qJ(3) * t391 + t465;
t389 = t396 ^ 2;
t397 = cos(pkin(9));
t503 = t391 * t397;
t359 = -qJD(4) + t503;
t470 = qJD(4) + t359;
t527 = t349 * (t389 * t391 + t470 * t397);
t478 = qJD(1) * t404;
t435 = -pkin(1) * t478 + qJD(3);
t504 = t391 * t396;
t352 = -pkin(3) * t397 - pkin(7) * t396 - pkin(2);
t473 = qJD(4) * t403;
t475 = qJD(3) * t397;
t497 = t397 * t404;
t523 = -(t399 * t400 + t403 * t497) * t514 + t352 * t473 + t403 * t475;
t336 = t352 - t519;
t375 = pkin(1) * t400 + qJ(3);
t498 = t397 * t403;
t489 = t399 * t336 + t375 * t498;
t522 = qJ(3) * t498 + t399 * t352;
t520 = t522 * qJD(4) + (-t399 * t497 + t400 * t403) * t514 + t399 * t475;
t517 = g(1) * t396;
t378 = g(3) * t386;
t513 = MDP(7) * t397;
t306 = t352 * t391 + t435;
t461 = t349 * t498;
t422 = -t306 * t399 - t461;
t501 = t396 * t399;
t467 = pkin(8) * t501;
t287 = -t391 * t467 - t422;
t512 = t287 * t402;
t351 = -qJD(5) + t359;
t511 = t351 * t397;
t510 = t384 * t397;
t509 = t386 * t397;
t508 = t388 * t397;
t506 = t388 * t403;
t469 = qJDD(1) * t400;
t476 = qJD(2) * t404;
t314 = qJ(3) * t388 + qJD(3) * t391 + (qJD(1) * t476 + t469) * pkin(1);
t307 = t389 * t314;
t502 = t391 * t399;
t500 = t396 * t403;
t499 = t397 * t399;
t356 = -qJDD(4) + t508;
t494 = t399 * t356;
t493 = t399 * t403;
t492 = t402 * t403;
t487 = t386 * pkin(2) + t384 * qJ(3);
t486 = g(2) * t384 - t378;
t390 = t397 ^ 2;
t484 = t389 + t390;
t393 = t403 ^ 2;
t483 = t399 ^ 2 - t393;
t482 = MDP(11) * t389;
t481 = MDP(12) * t389;
t480 = MDP(13) * t396;
t479 = MDP(14) * t396;
t472 = qJD(5) * t398;
t350 = -qJDD(5) + t356;
t343 = t350 * MDP(22);
t471 = t356 * MDP(15);
t466 = pkin(8) * t500;
t462 = qJ(3) * t499;
t459 = t391 * t500;
t457 = t396 * t492;
t298 = t352 * t388 + t447;
t456 = t399 * t298 + t306 * t473 + t314 * t498;
t455 = t529 * t396;
t369 = pkin(1) * t476 + qJD(3);
t454 = t336 * t473 + t369 * t498 + t399 * t463;
t453 = t391 * t477;
t451 = t391 * t473;
t449 = t397 * t474;
t448 = t349 * t474;
t446 = t369 * t484;
t445 = t484 * t314;
t444 = t484 * t388;
t443 = t356 + t508;
t442 = t391 * t470;
t415 = -t397 * t448 + t456;
t421 = t451 + t507;
t268 = -t421 * t396 * pkin(8) + t415;
t301 = t403 * t306;
t286 = -pkin(8) * t459 - t349 * t499 + t301;
t276 = -pkin(4) * t359 + t286;
t441 = qJD(5) * t276 + t268;
t439 = t391 * t465;
t438 = t390 * t314 + t307 - t486;
t295 = t403 * t298;
t433 = -t314 * t499 + t295;
t340 = t403 * t352;
t299 = -t466 + t340 + (-qJ(3) * t399 - pkin(4)) * t397;
t432 = qJD(5) * t299 + (-t462 - t466) * qJD(4) + t523;
t305 = -t467 + t522;
t364 = pkin(8) * t450;
t431 = qJD(5) * t305 - t364 + t520;
t430 = -t276 * t398 - t512;
t333 = t403 * t336;
t289 = -t466 + t333 + (-t375 * t399 - pkin(4)) * t397;
t297 = -t467 + t489;
t429 = t289 * t402 - t297 * t398;
t428 = t289 * t398 + t297 * t402;
t427 = t399 * t402 + t495;
t426 = -t492 + t496;
t425 = qJD(4) * (t359 + t503);
t365 = t396 * pkin(4) * t473;
t424 = -t435 * t396 - t365;
t325 = -t384 * t499 - t386 * t403;
t327 = -t384 * t403 + t386 * t499;
t420 = g(2) * t327 - g(3) * t325 + t403 * t307 + t415 * t397;
t326 = t384 * t498 - t386 * t399;
t328 = t384 * t399 + t386 * t498;
t419 = t389 * t349 * t473 - g(2) * t328 - g(3) * t326 + t399 * t307;
t418 = t434 - t485;
t380 = -pkin(2) - t519;
t417 = pkin(1) * t453 + t380 * t388;
t267 = -t388 * t466 - pkin(4) * t356 + (-t461 + (pkin(8) * t504 - t306) * t399) * qJD(4) + t433;
t281 = t287 * t472;
t288 = (t421 * pkin(4) + t314) * t396;
t409 = t468 * t427;
t290 = t409 * t396;
t309 = (pkin(4) * t502 + t349) * t396;
t394 = qJ(4) + qJ(5);
t383 = sin(t394);
t385 = cos(t394);
t316 = -t383 * t510 - t385 * t386;
t318 = t383 * t509 - t384 * t385;
t331 = t426 * t396;
t416 = g(2) * t318 - g(3) * t316 + (t398 * t267 + t441 * t402 - t281) * t397 - t288 * t331 - t309 * t290;
t273 = t388 * t457 + (-t388 * t496 - t409 * t391) * t396;
t311 = (-t457 + t458) * t391;
t312 = t427 * t504;
t414 = -t311 * t312 * MDP(18) + (-t312 * t351 + t273) * MDP(20) + (t311 * t351 - t274) * MDP(21) + (t311 ^ 2 - t312 ^ 2) * MDP(19) - t343;
t387 = t391 ^ 2;
t413 = -t359 ^ 2 - t387 * t389;
t412 = -t434 + t439;
t261 = t430 * qJD(5) + t402 * t267 - t398 * t268;
t291 = t402 * t396 * t440 + t526;
t317 = -t386 * t383 + t385 * t510;
t319 = t383 * t384 + t385 * t509;
t330 = t427 * t396;
t411 = -g(2) * t319 - g(3) * t317 - t261 * t397 + t288 * t330 + t309 * t291;
t410 = t435 * t484;
t408 = t281 + t385 * t517 + g(2) * t317 + (t287 * t351 - t267) * t398 - g(3) * t319 + t309 * t312;
t407 = (-t273 * t330 + t274 * t331 + t290 * t312 + t291 * t311) * MDP(19) + (-t273 * t397 + t290 * t351 + t331 * t350) * MDP(20) + (t274 * t397 + t291 * t351 + t330 * t350) * MDP(21) + (-t273 * t331 + t290 * t311) * MDP(18) + (t399 * t425 - t443 * t403) * t480 + (t443 * t399 + t403 * t425) * t479 + 0.2e1 * (t483 * t391 * qJD(4) - t388 * t493) * t481 + (t388 * t393 - 0.2e1 * t399 * t451) * t482 + t388 * MDP(4) + (t343 + t471) * t397;
t406 = -g(2) * t316 - g(3) * t318 + t309 * t311 + t383 * t517 + t261;
t405 = cos(qJ(1));
t401 = sin(qJ(1));
t376 = t384 * pkin(2);
t372 = pkin(4) * t501;
t367 = t403 * t463;
t345 = qJ(3) * t396 + t372;
t344 = -pkin(2) * t391 + t435;
t334 = t375 * t396 + t372;
t320 = t369 * t396 + t365;
t283 = -qJD(4) * t489 - t369 * t499 + t364 + t367;
t282 = (-t375 * t499 - t466) * qJD(4) + t454;
t270 = t422 * qJD(4) + t433;
t1 = [(-g(2) * t405 - g(3) * t401) * MDP(2) + (g(2) * t401 - g(3) * t405) * MDP(3) + qJDD(1) * MDP(1) + ((t388 * t404 - t453) * pkin(1) - t418) * MDP(5) + t407 + (-(-t336 * t474 + t367) * t359 - t333 * t356 + (-(-t369 * t399 - t375 * t473) * t359 + t375 * t494 - t270) * t397 + (t369 * t502 + t421 * t375) * t389 + t419) * MDP(16) + ((t429 * qJD(5) + t282 * t402 + t283 * t398) * t351 + t428 * t350 - t320 * t311 + t334 * t273 + t416) * MDP(24) + (t417 * t396 + t455) * MDP(8) + (t375 * t444 + t391 * t446 + t438) * MDP(9) + (-(-t428 * qJD(5) - t282 * t398 + t283 * t402) * t351 - t429 * t350 + t320 * t312 + t334 * t274 + t411) * MDP(23) + (((-qJDD(1) - t388) * t400 + (-qJD(1) - t391) * t476) * pkin(1) + t486) * MDP(6) + (t321 * t380 + t344 * t463 - g(2) * (pkin(1) * t405 + t487) - g(3) * (pkin(1) * t401 - qJ(3) * t386 + t376) + t349 * t446 + t375 * t445) * MDP(10) + ((-t375 * t449 + t454) * t359 + t489 * t356 + ((t369 * t391 + t375 * t388) * t403 + (-t375 * t391 - t349) * t474) * t389 + t420) * MDP(17) + (-t417 - t529) * t513; (t412 + t485) * MDP(5) + ((-t469 + (-qJD(2) + t391) * t478) * pkin(1) + t486) * MDP(6) + t407 + (-(t340 - t462) * t356 - t270 * t397 + t520 * t359 + (qJ(3) * t507 + (qJ(3) * t473 + t435 * t399) * t391) * t389 + t419) * MDP(16) + (qJ(3) * t444 + t410 * t391 + t438) * MDP(9) + (-(t299 * t402 - t305 * t398) * t350 + t345 * t274 + (t432 * t398 + t431 * t402) * t351 - t424 * t312 + t411) * MDP(23) + (t522 * t356 + (-qJ(3) * t449 + t523) * t359 + (qJ(3) * t506 - t448 + (-qJ(3) * t474 + t435 * t403) * t391) * t389 + t420) * MDP(17) + (-t321 * pkin(2) - t344 * t465 - g(2) * t487 - g(3) * t376 + (t445 + t378) * qJ(3) + t410 * t349) * MDP(10) + ((-t439 - t518) * t396 + t455) * MDP(8) + ((t299 * t398 + t305 * t402) * t350 + t345 * t273 + (-t431 * t398 + t432 * t402) * t351 + t424 * t311 + t416) * MDP(24) + (-t321 + t412 + t518) * t513; -MDP(7) * t508 + t396 * t388 * MDP(8) - t484 * t387 * MDP(9) + (-t484 * t391 * t349 + qJDD(3) + t418 - t518) * MDP(10) + (-t356 * t403 + t413 * t399) * MDP(16) + (t413 * t403 + t494) * MDP(17) + (t409 * t351 + t426 * t350 + (-t396 * t312 - t427 * t511) * t391) * MDP(23) + (t311 * t504 + t427 * t350 + (-t351 * t468 + t511 * t391) * t426) * MDP(24); (-t399 * t442 + t506) * t480 + (-t403 * t442 - t507) * t479 - t471 + (-g(2) * t325 - g(3) * t327 + t295 - t403 * t527 + (-t470 * t306 - t397 * t314 + t517) * t399) * MDP(16) + (g(1) * t500 + g(2) * t326 - g(3) * t328 - t301 * t359 + t399 * t527 - t456) * MDP(17) + ((-t286 * t398 - t512) * t351 + (-t312 * t459 - t402 * t350 + t351 * t472) * pkin(4) + t406) * MDP(23) + ((-t286 * t351 - t441) * t402 + (qJD(5) * t351 * t402 + t311 * t459 + t398 * t350) * pkin(4) + t408) * MDP(24) + t414 + (-t483 * t481 + t482 * t493) * t387; (t430 * t351 + t406) * MDP(23) + ((-t268 + (-qJD(5) - t351) * t276) * t402 + t408) * MDP(24) + t414;];
tau = t1;
