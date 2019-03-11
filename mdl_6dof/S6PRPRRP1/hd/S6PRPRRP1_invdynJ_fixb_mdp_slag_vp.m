% Calculate vector of inverse dynamics joint torques for
% S6PRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S6PRPRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:58:39
% EndTime: 2019-03-08 19:58:45
% DurationCPUTime: 4.50s
% Computational Cost: add. (2817->416), mult. (6517->570), div. (0->0), fcn. (5276->12), ass. (0->195)
t412 = cos(pkin(11));
t418 = sin(qJ(2));
t411 = sin(pkin(6));
t485 = qJD(1) * t411;
t461 = t418 * t485;
t379 = t412 * t461;
t409 = sin(pkin(11));
t421 = cos(qJ(2));
t460 = t421 * t485;
t339 = t409 * t460 + t379;
t417 = sin(qJ(4));
t420 = cos(qJ(4));
t449 = pkin(4) * t417 - pkin(9) * t420;
t373 = t449 * qJD(4);
t539 = t339 - t373;
t484 = qJD(2) * t411;
t458 = qJD(1) * t484;
t466 = t411 * qJDD(1);
t538 = t418 * t466 + t421 * t458;
t474 = qJD(5) * t417;
t537 = qJD(2) * t474 - qJDD(4);
t416 = sin(qJ(5));
t419 = cos(qJ(5));
t467 = qJDD(2) * t417;
t482 = qJD(2) * t420;
t326 = (qJD(4) * (qJD(5) + t482) + t467) * t416 + t537 * t419;
t401 = pkin(5) * t419 + pkin(4);
t536 = t401 * t420 + pkin(3);
t523 = pkin(5) * t416;
t535 = pkin(8) + t523;
t378 = t409 * t461;
t342 = t412 * t460 - t378;
t478 = qJD(4) * t417;
t502 = t416 * t420;
t398 = pkin(2) * t409 + pkin(8);
t512 = t398 * t416;
t534 = -t342 * t502 + t539 * t419 - t478 * t512;
t444 = -pkin(4) * t420 - pkin(9) * t417 - pkin(3);
t525 = pkin(2) * t412;
t359 = t444 - t525;
t473 = qJD(5) * t419;
t500 = t419 * t420;
t533 = -t342 * t500 + t359 * t473 - t539 * t416;
t375 = qJD(2) * pkin(2) + t460;
t335 = t409 * t375 + t379;
t333 = qJD(2) * pkin(8) + t335;
t414 = cos(pkin(6));
t394 = qJD(1) * t414 + qJD(3);
t532 = -t417 * t333 + t394 * t420;
t445 = t409 * t421 + t412 * t418;
t346 = t445 * t411;
t331 = t346 * t420 + t414 * t417;
t499 = t421 * t412;
t508 = t411 * t418;
t345 = t409 * t508 - t411 * t499;
t295 = -t331 * t416 + t345 * t419;
t347 = t445 * t414;
t363 = t409 * t418 - t499;
t410 = sin(pkin(10));
t413 = cos(pkin(10));
t320 = t347 * t413 - t363 * t410;
t509 = t411 * t417;
t301 = t320 * t420 - t413 * t509;
t321 = t410 * t347 + t363 * t413;
t303 = -t321 * t420 + t410 * t509;
t438 = t363 * t414;
t319 = -t410 * t445 - t413 * t438;
t322 = t410 * t438 - t413 * t445;
t531 = -g(1) * (-t303 * t416 - t322 * t419) - g(2) * (-t301 * t416 - t319 * t419) - g(3) * t295;
t503 = t414 * t421;
t530 = -t410 * t503 - t413 * t418;
t529 = pkin(5) * t326 + qJDD(6);
t313 = t420 * t333 + t417 * t394;
t310 = qJD(4) * pkin(9) + t313;
t395 = -qJD(5) + t482;
t433 = g(1) * t322 + g(2) * t319 - g(3) * t345;
t527 = (t395 * t398 + t310) * qJD(5) - t433;
t479 = qJD(4) * t416;
t483 = qJD(2) * t417;
t367 = t419 * t483 + t479;
t526 = t367 ^ 2;
t415 = -qJ(6) - pkin(9);
t521 = qJ(6) * t417;
t520 = qJDD(4) * pkin(4);
t465 = t420 * qJDD(2);
t468 = qJD(2) * qJD(4);
t362 = t417 * t468 + qJDD(5) - t465;
t519 = t362 * t416;
t469 = t419 * qJD(4);
t365 = t416 * t483 - t469;
t518 = t365 * t395;
t517 = t367 * t395;
t516 = t367 * t417;
t391 = t414 * qJDD(1) + qJDD(3);
t515 = t391 * t417;
t513 = t395 * t419;
t510 = t410 * t418;
t507 = t411 * t420;
t506 = t411 * t421;
t504 = t414 * t418;
t501 = t417 * t419;
t498 = qJDD(1) - g(3);
t334 = t375 * t412 - t378;
t317 = qJD(2) * t444 - t334;
t283 = -t310 * t416 + t419 * t317;
t281 = -qJ(6) * t367 + t283;
t278 = -pkin(5) * t395 + t281;
t497 = -t281 + t278;
t371 = t398 * t500;
t442 = pkin(5) * t417 - qJ(6) * t500;
t472 = qJD(6) * t419;
t496 = -t417 * t472 + t442 * qJD(4) + (-t371 + (-t359 + t521) * t416) * qJD(5) - t534;
t480 = qJD(4) * t398;
t495 = (-qJ(6) * qJD(5) - t480) * t501 + (-qJD(6) * t417 + (-qJ(6) * qJD(4) - qJD(5) * t398) * t420) * t416 + t533;
t372 = t449 * qJD(2);
t494 = t416 * t372 + t419 * t532;
t459 = t420 * t469;
t493 = -t326 * t501 - t365 * t459;
t453 = qJD(5) * t415;
t491 = t472 - t494 + (qJ(6) * t482 + t453) * t416;
t356 = t419 * t372;
t490 = -qJD(2) * t442 + t419 * t453 - t356 + (-qJD(6) + t532) * t416;
t488 = t416 * t359 + t371;
t407 = t417 ^ 2;
t487 = -t420 ^ 2 + t407;
t486 = MDP(20) * t416;
t481 = qJD(4) * t365;
t477 = qJD(4) * t420;
t476 = qJD(5) * t395;
t475 = qJD(5) * t416;
t309 = -qJD(4) * pkin(4) - t532;
t292 = pkin(5) * t365 + qJD(6) + t309;
t471 = t292 * qJD(4);
t470 = t419 * MDP(19);
t463 = t413 * t503;
t390 = t421 * t466;
t343 = qJDD(2) * pkin(2) - t418 * t458 + t390;
t307 = t409 * t343 + t538 * t412;
t305 = qJDD(2) * pkin(8) + t307;
t279 = qJDD(4) * pkin(9) + qJD(4) * t532 + t305 * t420 + t515;
t306 = t343 * t412 - t538 * t409;
t290 = qJD(2) * t373 + qJDD(2) * t444 - t306;
t462 = t419 * t279 + t416 * t290 + t317 * t473;
t457 = t420 * t468;
t454 = t398 + t523;
t284 = t310 * t419 + t317 * t416;
t282 = -qJ(6) * t365 + t284;
t448 = -t278 * t419 - t282 * t416;
t447 = t278 * t416 - t282 * t419;
t296 = t331 * t419 + t345 * t416;
t330 = t346 * t417 - t414 * t420;
t443 = t417 * t305 + t333 * t477 - t391 * t420 + t394 * t478;
t440 = -t395 * t473 + t519;
t439 = -g(3) * t414 + (-g(1) * t410 + g(2) * t413) * t411;
t437 = t310 * t475 - t462;
t300 = t320 * t417 + t413 * t507;
t302 = -t321 * t417 - t410 * t507;
t436 = g(1) * t302 + g(2) * t300 + g(3) * t330;
t435 = g(1) * t303 + g(2) * t301 + g(3) * t331;
t434 = g(1) * t321 - g(2) * t320 - g(3) * t346;
t432 = -t416 * t474 + t459;
t280 = t443 - t520;
t431 = -pkin(9) * t362 - t309 * t395;
t325 = -qJD(5) * t469 + (-t457 - t467) * t419 + t537 * t416;
t332 = -qJD(2) * pkin(3) - t334;
t399 = -pkin(3) - t525;
t430 = -qJDD(4) * t398 + (qJD(2) * t399 + t332 + t342) * qJD(4);
t288 = t419 * t290;
t429 = -t284 * qJD(5) - t416 * t279 + t288;
t428 = -g(1) * t530 - g(3) * t506;
t427 = pkin(9) * t476 - t280 + t436;
t426 = t436 - t443;
t422 = qJD(4) ^ 2;
t424 = -qJD(2) * t339 + t398 * t422 - t306 + t433 + (-pkin(3) + t399) * qJDD(2);
t423 = qJD(2) ^ 2;
t385 = t415 * t419;
t384 = t415 * t416;
t383 = pkin(2) * t463;
t382 = qJDD(4) * t420 - t417 * t422;
t381 = qJDD(4) * t417 + t420 * t422;
t361 = t365 ^ 2;
t351 = t367 * t478;
t350 = t419 * t359;
t341 = t363 * t484;
t340 = qJD(2) * t346;
t324 = -t416 * t521 + t488;
t316 = -qJ(6) * t501 + t350 + (-pkin(5) - t512) * t420;
t294 = -qJD(4) * t330 - t341 * t420;
t293 = qJD(4) * t331 - t341 * t417;
t276 = t295 * qJD(5) + t294 * t419 + t340 * t416;
t275 = -t296 * qJD(5) - t294 * t416 + t340 * t419;
t274 = t280 + t529;
t273 = -qJ(6) * t326 - qJD(6) * t365 - t437;
t272 = pkin(5) * t362 + qJ(6) * t325 - qJD(6) * t367 + t429;
t1 = [t498 * MDP(1) + (-t306 * t345 + t307 * t346 - t334 * t340 - t335 * t341 + t391 * t414 - g(3)) * MDP(5) + (-qJD(4) * t293 - qJDD(4) * t330 - t345 * t465) * MDP(11) + (-qJD(4) * t294 - qJDD(4) * t331 + t345 * t467) * MDP(12) + (-t275 * t395 + t293 * t365 + t295 * t362 + t326 * t330) * MDP(18) + (t276 * t395 + t293 * t367 - t296 * t362 - t325 * t330) * MDP(19) + (-t275 * t367 - t276 * t365 + t295 * t325 - t296 * t326) * MDP(20) + (t272 * t295 + t273 * t296 + t274 * t330 + t275 * t278 + t276 * t282 + t292 * t293 - g(3)) * MDP(21) + ((-t340 * t420 + t345 * t478) * MDP(11) + (t340 * t417 + t345 * t477) * MDP(12)) * qJD(2) + ((qJDD(2) * t421 - t418 * t423) * MDP(3) + (-qJDD(2) * t418 - t421 * t423) * MDP(4)) * t411; qJDD(2) * MDP(2) + (t390 - g(2) * (t463 - t510) + t428) * MDP(3) + (-g(1) * (t410 * t504 - t413 * t421) - g(2) * (-t410 * t421 - t413 * t504) - t498 * t508) * MDP(4) + (-g(2) * t383 + t334 * t339 - t335 * t342 + (g(2) * t510 + t306 * t412 + t307 * t409 + t428) * pkin(2)) * MDP(5) + (qJDD(2) * t407 + 0.2e1 * t417 * t457) * MDP(6) + 0.2e1 * (t417 * t465 - t487 * t468) * MDP(7) + t381 * MDP(8) + t382 * MDP(9) + (t417 * t430 - t420 * t424) * MDP(11) + (t417 * t424 + t420 * t430) * MDP(12) + (-t325 * t501 + t367 * t432) * MDP(13) + (-t473 * t516 + (-t367 * t477 + (qJD(5) * t365 + t325) * t417) * t416 + t493) * MDP(14) + (t325 * t420 + t362 * t501 - t395 * t432 + t351) * MDP(15) + ((t395 * t479 + t326) * t420 + (-t440 - t481) * t417) * MDP(16) + (-t362 * t420 - t395 * t478) * MDP(17) + (t350 * t362 + t534 * t395 + (t359 * t476 + t434) * t416 + (t365 * t480 - t288 + (qJD(4) * t309 + qJD(5) * t317 - t362 * t398 + t279) * t416 + t527 * t419) * t420 + (qJD(4) * t283 + t280 * t416 + t309 * t473 + t326 * t398 - t342 * t365) * t417) * MDP(18) + (-t488 * t362 + t533 * t395 + t434 * t419 + ((t309 * t419 + t367 * t398) * qJD(4) - t527 * t416 + t462) * t420 + (-t309 * t475 + t280 * t419 - t398 * t325 - t342 * t367 + (-t398 * t513 - t284) * qJD(4)) * t417) * MDP(19) + (t316 * t325 - t324 * t326 - t496 * t367 - t495 * t365 + t448 * t477 + (t447 * qJD(5) - t272 * t419 - t273 * t416 - t433) * t417) * MDP(20) + (t273 * t324 + t272 * t316 - g(1) * (t530 * pkin(2) - t535 * t321 + t536 * t322) - g(2) * (-pkin(2) * t510 + t536 * t319 + t320 * t535 + t383) - g(3) * (pkin(2) * t506 - t536 * t345 + t535 * t346) + t495 * t282 + t496 * t278 + t454 * t420 * t471 + (t274 * t454 + (pkin(5) * t473 - t342) * t292 + t433 * t415) * t417) * MDP(21); (t439 + t391) * MDP(5) + t382 * MDP(11) - t381 * MDP(12) + t351 * MDP(19) + t493 * MDP(20) + t439 * MDP(21) + (-t326 * MDP(18) + t325 * MDP(19) - t274 * MDP(21) + (t367 * t486 - t447 * MDP(21) + (t416 * MDP(18) + t470) * t395) * qJD(4)) * t420 + ((t481 - t519) * MDP(18) - t362 * t470 - t325 * t486 + (-t272 * t416 + t273 * t419 + t471) * MDP(21) + ((t365 * t416 + t367 * t419) * MDP(20) + t448 * MDP(21) + (t419 * MDP(18) - t416 * MDP(19)) * t395) * qJD(5)) * t417; MDP(8) * t467 + MDP(9) * t465 + qJDD(4) * MDP(10) + (qJD(4) * t313 - t332 * t483 + t426) * MDP(11) + (-t515 + (-qJD(2) * t332 - t305) * t420 + t435) * MDP(12) + (-t325 * t416 - t367 * t513) * MDP(13) + ((-t325 + t518) * t419 + (-t326 + t517) * t416) * MDP(14) + ((t395 * t500 - t516) * qJD(2) + t440) * MDP(15) + (t395 * t475 + t362 * t419 + (t365 * t417 - t395 * t502) * qJD(2)) * MDP(16) + t395 * MDP(17) * t483 + (-t283 * t483 - pkin(4) * t326 - t313 * t365 + t356 * t395 + (-t395 * t532 + t431) * t416 + t427 * t419) * MDP(18) + (pkin(4) * t325 + t284 * t483 - t313 * t367 - t494 * t395 - t427 * t416 + t431 * t419) * MDP(19) + (t325 * t384 + t326 * t385 - t490 * t367 - t491 * t365 + (t278 * t395 + t273) * t419 + (t282 * t395 - t272) * t416 - t435) * MDP(20) + (-t273 * t385 + t272 * t384 - t274 * t401 - g(1) * (-t302 * t401 - t303 * t415) - g(2) * (-t300 * t401 - t301 * t415) - g(3) * (-t330 * t401 - t331 * t415) + (-t395 * t523 - t313) * t292 + t491 * t282 + t490 * t278) * MDP(21) + (-t417 * t420 * MDP(6) + t487 * MDP(7)) * t423; t367 * t365 * MDP(13) + (-t361 + t526) * MDP(14) + (-t325 - t518) * MDP(15) + (-t326 - t517) * MDP(16) + t362 * MDP(17) + (-t284 * t395 - t309 * t367 + t429 + t531) * MDP(18) + (-t283 * t395 + t309 * t365 - g(1) * (-t303 * t419 + t322 * t416) - g(2) * (-t301 * t419 + t319 * t416) + g(3) * t296 + t437) * MDP(19) + (pkin(5) * t325 - t497 * t365) * MDP(20) + (t497 * t282 + (-t292 * t367 + t272 + t531) * pkin(5)) * MDP(21); (-t361 - t526) * MDP(20) + (t278 * t367 + t282 * t365 - t426 - t520 + t529) * MDP(21);];
tau  = t1;
