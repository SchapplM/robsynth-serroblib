% Calculate vector of inverse dynamics joint torques for
% S5RRPPR10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:52
% EndTime: 2019-12-31 19:44:59
% DurationCPUTime: 6.00s
% Computational Cost: add. (2118->474), mult. (4848->603), div. (0->0), fcn. (3272->8), ass. (0->206)
t423 = cos(qJ(2));
t414 = g(3) * t423;
t420 = sin(qJ(2));
t488 = qJDD(1) * t420;
t400 = pkin(6) * t488;
t489 = qJD(1) * qJD(2);
t476 = t423 * t489;
t447 = qJDD(2) * pkin(2) - pkin(6) * t476 - qJDD(3) - t400;
t473 = t447 - t414;
t421 = sin(qJ(1));
t424 = cos(qJ(1));
t462 = g(1) * t424 + g(2) * t421;
t542 = t462 * t420;
t427 = -t542 - t473;
t504 = t423 * pkin(2) + t420 * qJ(3);
t543 = -pkin(1) - t504;
t498 = qJD(2) * t420;
t540 = qJ(4) * t498 - qJD(4) * t423;
t407 = t423 * qJDD(1);
t477 = t420 * t489;
t437 = -t477 + t407;
t501 = qJD(1) * t420;
t417 = sin(pkin(8));
t482 = t417 * t501;
t418 = cos(pkin(8));
t490 = t418 * qJD(2);
t362 = t482 - t490;
t480 = t418 * t501;
t499 = qJD(2) * t417;
t364 = t480 + t499;
t419 = sin(qJ(5));
t422 = cos(qJ(5));
t310 = t362 * t419 + t364 * t422;
t405 = t418 * qJDD(2);
t438 = t476 + t488;
t332 = t417 * t438 - t405;
t333 = qJDD(2) * t417 + t418 * t438;
t471 = -t332 * t422 + t419 * t333;
t281 = qJD(5) * t310 + t471;
t495 = qJD(4) * t364;
t527 = qJ(4) * t333;
t429 = t447 + t495 + t527;
t536 = pkin(3) + pkin(4);
t279 = -t332 * t536 + t429;
t538 = t279 + t542;
t537 = -0.2e1 * pkin(1);
t360 = t364 ^ 2;
t535 = pkin(3) * t332;
t534 = pkin(6) * t364;
t533 = g(1) * t421;
t530 = g(2) * t424;
t529 = g(3) * t420;
t528 = -pkin(7) + qJ(3);
t526 = qJ(4) * t418;
t525 = qJ(4) * t423;
t460 = pkin(2) * t420 - qJ(3) * t423;
t348 = qJD(2) * t460 - qJD(3) * t420;
t306 = qJD(1) * t348 + qJDD(1) * t543;
t339 = pkin(6) * t437 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t286 = t306 * t417 + t339 * t418;
t282 = qJ(4) * t477 + (-qJ(4) * qJDD(1) - qJD(1) * qJD(4)) * t423 + t286;
t524 = t282 * t418;
t523 = t286 * t418;
t519 = t364 * t419;
t308 = -t362 * t422 + t519;
t500 = qJD(1) * t423;
t393 = qJD(5) + t500;
t522 = t308 * t393;
t521 = t310 * t393;
t520 = t348 * t418;
t371 = t460 * qJD(1);
t518 = t371 * t418;
t517 = t417 * t420;
t516 = t417 * t422;
t515 = t417 * t423;
t514 = t418 * t420;
t513 = t418 * t423;
t512 = t420 * t421;
t511 = t420 * t424;
t510 = t421 * t418;
t509 = t421 * t423;
t508 = t423 * t424;
t507 = t424 * t417;
t367 = t417 * t419 + t418 * t422;
t439 = t367 * t423;
t506 = -qJD(1) * t439 - qJD(5) * t367;
t481 = t417 * t500;
t486 = t419 * t513;
t491 = qJD(5) * t422;
t492 = qJD(5) * t419;
t505 = -qJD(1) * t486 + t417 * t491 - t418 * t492 + t422 * t481;
t358 = t543 * qJD(1);
t403 = pkin(6) * t500;
t381 = qJD(2) * qJ(3) + t403;
t313 = t358 * t417 + t381 * t418;
t335 = pkin(6) * t513 + t417 * t543;
t415 = t420 ^ 2;
t416 = t423 ^ 2;
t503 = t415 - t416;
t502 = qJD(1) * t418;
t497 = qJD(2) * t423;
t496 = qJD(3) * t418;
t494 = qJD(4) * t417;
t487 = pkin(6) * t498;
t485 = t332 * t419 + t333 * t422 + t362 * t491;
t484 = -pkin(6) * t417 - pkin(3);
t483 = qJ(3) * t498;
t479 = qJD(4) * t514;
t478 = qJ(3) * t407;
t475 = qJ(4) * t417 + pkin(2);
t283 = -t429 + t535;
t474 = -t283 - t414;
t285 = t306 * t418 - t339 * t417;
t446 = pkin(3) * t407 + qJDD(4) - t285;
t284 = -pkin(3) * t477 + t446;
t277 = pkin(4) * t437 - pkin(7) * t333 + t284;
t278 = pkin(7) * t332 + t282;
t472 = t277 * t422 - t419 * t278;
t312 = t358 * t418 - t381 * t417;
t387 = pkin(6) * t515;
t334 = t418 * t543 - t387;
t444 = -t417 * t536 + t526;
t470 = -t444 * t500 + t403 + t494;
t469 = g(1) * t418 * t511 + g(2) * t420 * t510 + qJD(3) * t481 + t417 * t478;
t468 = pkin(1) * t424 + pkin(2) * t508 + pkin(6) * t421 + qJ(3) * t511;
t467 = t418 * t478;
t466 = -qJ(3) * t332 * t418 - t362 * t496 - t529;
t465 = t484 * t420;
t352 = t417 * t509 + t418 * t424;
t354 = t423 * t507 - t510;
t464 = -g(1) * t352 + g(2) * t354;
t353 = t418 * t509 - t507;
t355 = t417 * t421 + t418 * t508;
t463 = g(1) * t353 - g(2) * t355;
t461 = qJD(2) * pkin(2) - pkin(6) * t501 - qJD(3);
t322 = t335 - t525;
t459 = pkin(3) * t417 - t526;
t458 = pkin(6) * t362 - t417 * t461;
t456 = t419 * t277 + t422 * t278;
t299 = pkin(3) * t500 + qJD(4) - t312;
t287 = pkin(4) * t500 - pkin(7) * t364 + t299;
t304 = -qJ(4) * t500 + t313;
t289 = pkin(7) * t362 + t304;
t274 = t287 * t422 - t289 * t419;
t275 = t287 * t419 + t289 * t422;
t410 = t423 * pkin(3);
t305 = pkin(4) * t423 + t387 + t410 + (-pkin(7) * t420 - t543) * t418;
t311 = pkin(7) * t517 + t322;
t455 = t305 * t422 - t311 * t419;
t454 = t305 * t419 + t311 * t422;
t453 = t352 * t422 - t353 * t419;
t452 = t352 * t419 + t353 * t422;
t368 = -t418 * t419 + t516;
t451 = qJ(3) * t333 + qJD(3) * t364;
t356 = t417 * t371;
t328 = -pkin(6) * t480 + t356;
t338 = t417 * t348;
t319 = -t418 * t487 + t338;
t449 = pkin(3) * t418 + t475;
t448 = pkin(6) + t459;
t443 = -pkin(6) * qJDD(2) + t489 * t537;
t378 = t528 * t418;
t430 = -pkin(7) * t513 + (-pkin(4) + t484) * t420;
t442 = -qJD(1) * t430 + qJD(3) * t417 - qJD(5) * t378 + t518;
t377 = t528 * t417;
t398 = qJ(4) * t501;
t440 = -pkin(6) * t514 + pkin(7) * t515;
t441 = qJD(1) * t440 - qJD(5) * t377 + t356 + t398 - t496;
t343 = t367 * t420;
t436 = t364 * t492 - t485;
t435 = -pkin(6) + t444;
t426 = qJD(1) ^ 2;
t434 = pkin(1) * t426 + t462;
t425 = qJD(2) ^ 2;
t433 = pkin(6) * t425 + qJDD(1) * t537 + t530;
t432 = qJ(4) * t364 + t461;
t431 = t543 * t533;
t412 = t424 * pkin(6);
t395 = g(1) * t512;
t391 = qJ(3) * t508;
t388 = qJ(3) * t509;
t366 = -qJDD(5) - t437;
t357 = t418 * t536 + t475;
t346 = t362 * t500;
t342 = t419 * t514 - t420 * t516;
t340 = t448 * t420;
t331 = t459 * t500 + t403;
t327 = pkin(6) * t482 + t518;
t326 = -t334 + t410;
t321 = t435 * t420;
t318 = t417 * t487 + t520;
t317 = qJD(1) * t465 - t518;
t316 = t328 + t398;
t315 = t448 * t497 - t479;
t307 = qJD(2) * t465 - t520;
t303 = t354 * t419 + t355 * t422;
t302 = t354 * t422 - t355 * t419;
t298 = t435 * t497 + t479;
t297 = t319 + t540;
t295 = pkin(3) * t362 - t432;
t294 = qJD(5) * t368 * t420 + qJD(2) * t439;
t293 = qJD(2) * t486 + qJD(5) * t343 - t497 * t516;
t291 = qJD(2) * t440 + t338 + t540;
t290 = qJD(2) * t430 - t520;
t288 = -t362 * t536 + t432;
t1 = [0.2e1 * (t407 * t420 - t489 * t503) * MDP(5) + (qJDD(2) * t420 + t423 * t425) * MDP(6) + (qJDD(2) * t423 - t420 * t425) * MDP(7) + (t443 * t420 + (-t433 + t533) * t423) * MDP(9) + (t420 * t433 + t423 * t443 - t395) * MDP(10) + ((pkin(6) * t332 - t447 * t417 + (qJD(1) * t334 + t312) * qJD(2)) * t420 + (-qJD(1) * t318 + qJD(2) * t458 - qJDD(1) * t334 - t285) * t423 + t463) * MDP(11) + ((pkin(6) * t333 - t447 * t418 + (-qJD(1) * t335 - t313) * qJD(2)) * t420 + (qJD(1) * t319 + qJDD(1) * t335 + t286 + (-t418 * t461 + t534) * qJD(2)) * t423 + t464) * MDP(12) + (-t318 * t364 - t319 * t362 - t332 * t335 - t333 * t334 + t395 + (-t312 * t418 - t313 * t417) * t497 + (-t285 * t418 - t286 * t417 - t530) * t420) * MDP(13) + (t286 * t335 + t313 * t319 + t285 * t334 + t312 * t318 - g(1) * t412 - g(2) * t468 - t431 + (-t420 * t447 - t461 * t497) * pkin(6)) * MDP(14) + (t315 * t362 + t332 * t340 + (t283 * t417 + (-qJD(1) * t326 - t299) * qJD(2)) * t420 + (qJD(1) * t307 + qJDD(1) * t326 + t295 * t499 + t284) * t423 + t463) * MDP(15) + (-t297 * t362 + t307 * t364 - t322 * t332 + t326 * t333 + t395 + (t299 * t418 - t304 * t417) * t497 + (-t282 * t417 + t284 * t418 - t530) * t420) * MDP(16) + (-t315 * t364 - t333 * t340 + (-t283 * t418 + (qJD(1) * t322 + t304) * qJD(2)) * t420 + (-qJD(1) * t297 - qJDD(1) * t322 - t295 * t490 - t282) * t423 - t464) * MDP(17) + (t282 * t322 + t304 * t297 + t283 * t340 + t295 * t315 + t284 * t326 + t299 * t307 - g(1) * (-pkin(3) * t353 - qJ(4) * t352 + t412) - g(2) * (pkin(3) * t355 + qJ(4) * t354 + t468) - t431) * MDP(18) + (t294 * t310 - t343 * t436) * MDP(19) + (-t281 * t343 - t293 * t310 - t294 * t308 + t342 * t436) * MDP(20) + (t294 * t393 - t310 * t498 - t343 * t366 - t423 * t436) * MDP(21) + (-t281 * t423 - t293 * t393 + t308 * t498 + t342 * t366) * MDP(22) + (-t366 * t423 - t393 * t498) * MDP(23) + ((t290 * t422 - t291 * t419) * t393 - t455 * t366 + t472 * t423 - t274 * t498 + t298 * t308 + t321 * t281 + t279 * t342 + t288 * t293 + g(1) * t452 - g(2) * t303 + (-t275 * t423 - t393 * t454) * qJD(5)) * MDP(24) + (-(t290 * t419 + t291 * t422) * t393 + t454 * t366 - t456 * t423 + t275 * t498 + t298 * t310 - t321 * t436 + t279 * t343 + t288 * t294 + g(1) * t453 - g(2) * t302 + (-t274 * t423 - t393 * t455) * qJD(5)) * MDP(25) + qJDD(1) * MDP(1) + (qJDD(1) * t415 + 0.2e1 * t420 * t476) * MDP(4) + (-t530 + t533) * MDP(2) + t462 * MDP(3); MDP(6) * t488 + MDP(7) * t407 + qJDD(2) * MDP(8) + (t420 * t434 - t400 - t414) * MDP(9) + (t529 + (-pkin(6) * qJDD(1) + t434) * t423) * MDP(10) + (-pkin(2) * t332 + t473 * t418 + ((-qJ(3) * t499 - t312) * t420 + (t327 - t458) * t423) * qJD(1) + t469) * MDP(11) + (t467 - pkin(2) * t333 + t427 * t417 + ((-qJ(3) * t490 + t313) * t420 + (-t534 - t328 + (qJD(3) + t461) * t418) * t423) * qJD(1)) * MDP(12) + (t523 + t327 * t364 + t328 * t362 + (t312 * t502 - t462) * t423 + (t313 * t500 - t285 + t451) * t417 + t466) * MDP(13) + (t447 * pkin(2) - t313 * t328 - t312 * t327 + t461 * t403 - g(1) * (-pkin(2) * t511 + t391) - g(2) * (-pkin(2) * t512 + t388) - g(3) * t504 + (-t312 * t417 + t313 * t418) * qJD(3) + (-t285 * t417 + t523) * qJ(3)) * MDP(14) + (-t332 * t449 + t474 * t418 + (-t331 - t494) * t362 + (t299 * t420 - t317 * t423 + (-t295 * t423 - t483) * t417) * qJD(1) + t469) * MDP(15) + (t524 + t316 * t362 - t317 * t364 + (-t299 * t502 - t462) * t423 + (t304 * t500 + t284 + t451) * t417 + t466) * MDP(16) + (-t467 + t331 * t364 + t333 * t449 + (t474 + t495 + t542) * t417 + (-t304 * t420 + t316 * t423 + (t483 + (-qJD(3) + t295) * t423) * t418) * qJD(1)) * MDP(17) + (qJ(3) * t524 - t295 * t331 - t299 * t317 - g(1) * t391 - g(2) * t388 - g(3) * (pkin(3) * t513 + t504) + (-t316 + t496) * t304 + (-g(3) * t525 + qJ(3) * t284 + qJD(3) * t299 - qJD(4) * t295) * t417 + (-t283 + t542) * t449) * MDP(18) + (t310 * t506 - t368 * t436) * MDP(19) + (-t281 * t368 - t308 * t506 - t310 * t505 + t367 * t436) * MDP(20) + (t310 * t501 - t366 * t368 + t393 * t506) * MDP(21) + (-t308 * t501 + t366 * t367 - t393 * t505) * MDP(22) + t393 * MDP(23) * t501 + (-(t377 * t422 - t378 * t419) * t366 + t357 * t281 - g(3) * t439 + (t419 * t441 + t422 * t442) * t393 + t470 * t308 + t505 * t288 + t274 * t501 + t538 * t367) * MDP(24) + ((t377 * t419 + t378 * t422) * t366 - t357 * t436 + (-t419 * t442 + t422 * t441) * t393 + t470 * t310 + t506 * t288 - t275 * t501 + (-t414 + t538) * t368) * MDP(25) + (-MDP(4) * t420 * t423 + MDP(5) * t503) * t426; (t312 * t364 + t313 * t362 + t427) * MDP(14) + (t535 - t527 + t304 * t362 + (-qJD(4) - t299) * t364 + t427) * MDP(18) + (-t281 - t521) * MDP(24) + (t436 + t522) * MDP(25) + (MDP(11) + MDP(15)) * (t417 * t488 - t405 + (-t364 + t499) * t500) + (MDP(12) - MDP(17)) * (t346 + t333) + (MDP(13) + MDP(16)) * (-t362 ^ 2 - t360); (t362 * t364 + t437) * MDP(15) + (-t346 + t333) * MDP(16) + (-t416 * t426 - t360) * MDP(17) + (-g(3) * t517 - g(1) * t354 - g(2) * t352 + t295 * t364 + (-pkin(3) * t498 + t304 * t423) * qJD(1) + t446) * MDP(18) + (-t308 * t364 - t366 * t422) * MDP(24) + (-t310 * t364 + t366 * t419) * MDP(25) + (-MDP(24) * t419 - MDP(25) * t422) * t393 ^ 2; t310 * t308 * MDP(19) + (-t308 ^ 2 + t310 ^ 2) * MDP(20) + (t485 + t522) * MDP(21) + (-t471 + t521) * MDP(22) - t366 * MDP(23) + (-g(1) * t302 - g(2) * t453 + g(3) * t342 + t275 * t393 - t288 * t310 + t472) * MDP(24) + (g(1) * t303 + g(2) * t452 + g(3) * t343 + t274 * t393 + t288 * t308 - t456) * MDP(25) + (-MDP(21) * t519 - MDP(22) * t310 - MDP(24) * t275 - MDP(25) * t274) * qJD(5);];
tau = t1;
