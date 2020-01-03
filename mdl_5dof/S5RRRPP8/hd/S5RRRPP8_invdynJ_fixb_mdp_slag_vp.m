% Calculate vector of inverse dynamics joint torques for
% S5RRRPP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP8_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:35
% EndTime: 2019-12-31 21:09:43
% DurationCPUTime: 5.56s
% Computational Cost: add. (2843->515), mult. (6102->604), div. (0->0), fcn. (3767->6), ass. (0->206)
t429 = cos(qJ(2));
t413 = t429 * qJDD(1);
t426 = sin(qJ(2));
t503 = qJD(1) * qJD(2);
t568 = -t426 * t503 + t413;
t358 = qJDD(3) - t568;
t550 = pkin(3) + qJ(5);
t490 = t550 * t358;
t430 = cos(qJ(1));
t533 = t426 * t430;
t427 = sin(qJ(1));
t535 = t426 * t427;
t574 = g(1) * t533 + g(2) * t535;
t415 = t426 * pkin(7);
t419 = t429 * pkin(2);
t498 = -pkin(1) - t419;
t461 = t498 - t415;
t352 = t461 * qJD(1);
t515 = qJD(1) * t429;
t408 = pkin(6) * t515;
t376 = qJD(2) * pkin(7) + t408;
t425 = sin(qJ(3));
t428 = cos(qJ(3));
t319 = -t428 * t352 + t376 * t425;
t513 = qJD(2) * t425;
t516 = qJD(1) * t426;
t362 = t428 * t516 + t513;
t456 = pkin(4) * t362 + t319;
t505 = qJD(4) + t456;
t509 = qJD(3) * t426;
t573 = qJD(1) * t509 - qJDD(2);
t502 = qJDD(1) * t426;
t313 = (qJD(2) * (qJD(3) + t515) + t502) * t425 + t573 * t428;
t350 = t358 * qJ(4);
t391 = -qJD(3) + t515;
t379 = qJD(4) * t391;
t571 = t379 - t350;
t540 = t362 * t391;
t570 = -t313 - t540;
t536 = t425 * t429;
t501 = pkin(3) * t536;
t510 = qJD(3) * t425;
t569 = pkin(3) * t510 - qJD(1) * t501 - qJD(4) * t425 - t408;
t470 = g(1) * t430 + g(2) * t427;
t567 = -pkin(4) * t313 + qJDD(5);
t565 = -0.2e1 * pkin(1);
t507 = t428 * qJD(2);
t360 = t425 * t516 - t507;
t564 = t360 ^ 2;
t357 = t362 ^ 2;
t387 = t391 ^ 2;
t563 = 0.2e1 * t350;
t562 = pkin(4) + pkin(7);
t561 = pkin(3) * t358;
t488 = t429 * t503;
t312 = -qJD(3) * t507 + (-t488 - t502) * t428 + t573 * t425;
t560 = pkin(4) * t312;
t558 = pkin(4) * t360;
t557 = pkin(7) * t358;
t556 = g(1) * t427;
t553 = g(2) * t430;
t552 = g(3) * t426;
t551 = g(3) * t429;
t549 = pkin(7) * qJD(3);
t548 = qJ(4) * t313;
t547 = qJ(4) * t360;
t546 = qJ(4) * t428;
t320 = t425 * t352 + t428 * t376;
t309 = qJ(4) * t391 - t320;
t545 = t309 * t391;
t544 = t312 * t425;
t543 = t320 * t391;
t542 = t360 * t362;
t541 = t360 * t391;
t539 = t362 * t428;
t473 = pkin(2) * t426 - pkin(7) * t429;
t365 = t473 * qJD(1);
t538 = t365 * t428;
t537 = t425 * t426;
t534 = t426 * t428;
t532 = t427 * t429;
t531 = t428 * t429;
t530 = t429 * t430;
t529 = t430 * t425;
t465 = qJ(5) * t425 - t546;
t451 = t465 * t429;
t528 = -qJD(1) * t451 + qJD(3) * t465 - qJD(5) * t428 + t569;
t508 = qJD(3) * t428;
t527 = -qJ(4) * t508 + t515 * t546 + t569;
t368 = t473 * qJD(2);
t518 = -t415 - t419;
t370 = -pkin(1) + t518;
t526 = t425 * t368 + t370 * t508;
t348 = t425 * t365;
t485 = pkin(6) * t428 - qJ(4);
t443 = -pkin(4) * t536 - t426 * t485;
t525 = -qJD(1) * t443 - t562 * t510 - t348;
t378 = t562 * t428;
t497 = -pkin(6) * t425 - pkin(3);
t441 = pkin(4) * t531 + (-qJ(5) + t497) * t426;
t524 = -qJD(1) * t441 + qJD(3) * t378 + t538;
t523 = t574 * t425;
t522 = t574 * t428;
t520 = pkin(3) * t537 - qJ(4) * t534;
t396 = pkin(6) * t531;
t519 = t425 * t370 + t396;
t422 = t426 ^ 2;
t517 = -t429 ^ 2 + t422;
t514 = qJD(2) * t362;
t512 = qJD(2) * t426;
t511 = qJD(2) * t429;
t506 = -qJD(4) - t319;
t305 = t320 - t558;
t504 = -qJD(5) - t305;
t394 = pkin(6) * t536;
t496 = t391 * t513;
t495 = t391 * t507;
t375 = -qJD(2) * pkin(2) + pkin(6) * t516;
t455 = -qJ(4) * t362 + t375;
t301 = t360 * t550 + t455;
t494 = t301 * t510;
t493 = t301 * t508;
t492 = t391 * t510;
t491 = t425 * t509;
t486 = -t425 * qJ(4) - pkin(2);
t406 = pkin(6) * t502;
t343 = -qJDD(2) * pkin(2) + pkin(6) * t488 + t406;
t439 = qJ(4) * t312 - qJD(4) * t362 + t343;
t289 = qJD(5) * t360 + t313 * t550 + t439;
t484 = -t289 - t551;
t344 = t425 * t532 + t428 * t430;
t345 = t427 * t531 - t529;
t483 = -t344 * pkin(3) + qJ(4) * t345;
t346 = -t427 * t428 + t429 * t529;
t347 = t425 * t427 + t428 * t530;
t482 = -t346 * pkin(3) + qJ(4) * t347;
t481 = t370 * t428 - t394;
t324 = qJD(1) * t368 + qJDD(1) * t461;
t342 = t568 * pkin(6) + qJDD(2) * pkin(7);
t478 = -t428 * t324 + t425 * t342 + t352 * t510 + t376 * t508;
t477 = -t425 * t324 - t428 * t342 - t352 * t508 + t376 * t510;
t476 = t426 * pkin(3) * t508 + pkin(6) * t511 + qJ(4) * t491 + qJD(2) * t501;
t475 = pkin(3) * t531 + qJ(4) * t536 - t518;
t474 = t497 * t426;
t472 = g(1) * t344 - g(2) * t346;
t471 = g(1) * t345 - g(2) * t347;
t469 = -t345 * pkin(3) + t430 * pkin(6) - qJ(4) * t344;
t327 = qJ(4) * t429 - t519;
t468 = -qJD(3) * t396 + t368 * t428 - t370 * t510;
t467 = -qJD(4) * t429 + t526;
t466 = qJD(3) * t375 - t557;
t297 = t391 * t550 + t505;
t300 = qJD(5) - t309 - t558;
t464 = t297 * t428 - t300 * t425;
t308 = pkin(3) * t391 - t506;
t463 = t308 * t428 + t309 * t425;
t462 = -qJDD(4) - t478;
t460 = pkin(3) * t428 - t486;
t459 = -qJ(5) * t537 - t520;
t458 = t391 * t549 - t551;
t291 = t477 + t571;
t454 = -pkin(6) * qJDD(2) + t503 * t565;
t453 = t358 * t425 - t391 * t508;
t452 = t358 * t428 + t492;
t292 = pkin(3) * t313 + t439;
t450 = -t292 + t458;
t432 = qJD(1) ^ 2;
t449 = pkin(1) * t432 + t470;
t448 = t428 * t550 - t486;
t431 = qJD(2) ^ 2;
t447 = pkin(6) * t431 + qJDD(1) * t565 + t553;
t446 = t430 * pkin(1) + pkin(2) * t530 + t347 * pkin(3) + t427 * pkin(6) + pkin(7) * t533 + qJ(4) * t346;
t311 = pkin(3) * t360 + t455;
t445 = t311 * t391 + t557;
t440 = g(1) * t346 + g(2) * t344 + g(3) * t537 - t478;
t437 = -qJDD(4) + t440;
t299 = -t312 - t541;
t436 = g(1) * t347 + g(2) * t345 + g(3) * t534 + t477;
t435 = t311 * t362 - t437;
t434 = t301 * t362 - t437 - t560;
t433 = -t301 * t360 - t436 + t567;
t418 = t429 * pkin(3);
t416 = t426 * pkin(6);
t402 = g(1) * t535;
t399 = pkin(7) * t530;
t395 = pkin(7) * t532;
t377 = t562 * t425;
t334 = t416 + t520;
t328 = t418 - t481;
t326 = t416 - t459;
t323 = pkin(3) * t362 + t547;
t322 = qJD(1) * t474 - t538;
t321 = t485 * t516 - t348;
t317 = -pkin(4) * t537 - t327;
t314 = qJ(5) * t429 + t394 + t418 + (pkin(4) * t426 - t370) * t428;
t306 = t362 * t550 + t547;
t303 = (-qJ(4) * t511 - qJD(4) * t426) * t428 + t476;
t302 = qJD(2) * t474 - t468;
t298 = -qJ(4) * t512 + (t426 * t507 + t429 * t510) * pkin(6) - t467;
t296 = qJD(2) * t451 + (qJD(5) * t425 + (qJ(5) * qJD(3) - qJD(4)) * t428) * t426 + t476;
t295 = (-pkin(4) * t534 - t394) * qJD(3) + t443 * qJD(2) + t467;
t294 = -pkin(4) * t491 + qJD(2) * t441 + qJD(5) * t429 - t468;
t293 = -t462 - t561;
t290 = -t291 + t567;
t288 = qJD(5) * t391 - t462 - t490 - t560;
t1 = [((t312 - t495) * t429 + (t452 + t514) * t426) * MDP(13) + ((t313 + t496) * t429 + (-qJD(2) * t360 - t453) * t426) * MDP(14) + (-t358 * t429 - t391 * t512) * MDP(15) + (-t468 * t391 + t481 * t358 + ((pkin(6) * t360 + t375 * t425) * qJD(2) + t478) * t429 + (t375 * t508 - t319 * qJD(2) + t343 * t425 + (t313 - t496) * pkin(6)) * t426 + t471) * MDP(16) + (t526 * t391 - t519 * t358 + (t375 * t507 + (-t492 + t514) * pkin(6) - t477) * t429 + (-t375 * t510 - t320 * qJD(2) + t343 * t428 + (-t312 - t495) * pkin(6)) * t426 - t472) * MDP(17) + (t298 * t360 + t302 * t362 - t312 * t328 + t313 * t327 + t402 + t463 * t511 + (-t553 + t291 * t425 + t293 * t428 + (-t308 * t425 + t309 * t428) * qJD(3)) * t426) * MDP(18) + (-t302 * t391 - t303 * t360 - t313 * t334 + t328 * t358 + (-t311 * t513 - t293) * t429 + (qJD(2) * t308 - t292 * t425 - t311 * t508) * t426 - t471) * MDP(19) + (t298 * t391 - t303 * t362 + t312 * t334 - t327 * t358 + (-t311 * t507 + t291) * t429 + (-qJD(2) * t309 - t292 * t428 + t311 * t510) * t426 + t472) * MDP(20) + (-g(1) * t469 - g(2) * t446 + t291 * t327 + t292 * t334 + t293 * t328 + t309 * t298 + t308 * t302 + t311 * t303 - t461 * t556) * MDP(21) + (t294 * t362 - t295 * t360 - t312 * t314 - t313 * t317 + t402 + t464 * t511 + (-t553 + t288 * t428 - t290 * t425 + (-t297 * t425 - t300 * t428) * qJD(3)) * t426) * MDP(22) + (-t295 * t391 - t296 * t362 + t312 * t326 + t317 * t358 + (-t301 * t507 - t290) * t429 + (qJD(2) * t300 - t289 * t428 + t494) * t426 + t472) * MDP(23) + (t294 * t391 + t296 * t360 + t313 * t326 - t314 * t358 + (t301 * t513 + t288) * t429 + (-qJD(2) * t297 + t289 * t425 + t493) * t426 + t471) * MDP(24) + (t289 * t326 + t301 * t296 + t288 * t314 + t297 * t294 + t290 * t317 + t300 * t295 - g(1) * (-qJ(5) * t345 + t469) - g(2) * (pkin(4) * t533 + qJ(5) * t347 + t446) - (-t426 * t562 + t498) * t556) * MDP(25) + qJDD(1) * MDP(1) + (qJDD(1) * t422 + 0.2e1 * t426 * t488) * MDP(4) + 0.2e1 * (t413 * t426 - t503 * t517) * MDP(5) + (qJDD(2) * t426 + t429 * t431) * MDP(6) + (qJDD(2) * t429 - t426 * t431) * MDP(7) + (t454 * t426 + (-t447 + t556) * t429) * MDP(9) + (t426 * t447 + t429 * t454 - t402) * MDP(10) + (-t312 * t534 + (t429 * t507 - t491) * t362) * MDP(11) + ((-t360 * t428 - t362 * t425) * t511 + (t544 - t313 * t428 + (t360 * t425 - t539) * qJD(3)) * t426) * MDP(12) + (-t553 + t556) * MDP(2) + t470 * MDP(3); MDP(6) * t502 + MDP(7) * t413 + qJDD(2) * MDP(8) + (t426 * t449 - t406 - t551) * MDP(9) + (t552 + (-pkin(6) * qJDD(1) + t449) * t429) * MDP(10) + (-t391 * t539 - t544) * MDP(11) + ((-t312 + t541) * t428 + (-t313 + t540) * t425) * MDP(12) + ((-t362 * t426 + t391 * t531) * qJD(1) + t453) * MDP(13) + ((t360 * t426 - t391 * t536) * qJD(1) + t452) * MDP(14) + t391 * MDP(15) * t516 + (-pkin(2) * t313 + t466 * t425 + (-t551 - t343 + (t365 + t549) * t391) * t428 + (-t375 * t536 + t319 * t426 + (-t360 * t429 + t391 * t537) * pkin(6)) * qJD(1) + t522) * MDP(16) + (pkin(2) * t312 - t348 * t391 + t466 * t428 + (t343 - t458) * t425 + (-t375 * t531 + t320 * t426 + (-t362 * t429 + t391 * t534) * pkin(6)) * qJD(1) - t523) * MDP(17) + (-t552 - t321 * t360 - t322 * t362 - t470 * t429 + (-t291 - t391 * t308 + (qJD(3) * t362 - t313) * pkin(7)) * t428 + (t293 - t545 + (qJD(3) * t360 - t312) * pkin(7)) * t425) * MDP(18) + (-t308 * t516 + t313 * t460 + t322 * t391 - t360 * t527 + t425 * t445 - t428 * t450 - t522) * MDP(19) + (t309 * t516 - t312 * t460 - t321 * t391 - t362 * t527 + t425 * t450 + t428 * t445 + t523) * MDP(20) + (-t309 * t321 - t308 * t322 - g(1) * t399 - g(2) * t395 - g(3) * t475 + t527 * t311 + (qJD(3) * t463 - t291 * t428 + t293 * t425) * pkin(7) + (t470 * t426 - t292) * t460) * MDP(21) + (-t552 + t288 * t425 + t290 * t428 - t312 * t377 - t313 * t378 + t524 * t362 - t525 * t360 + t464 * qJD(3) + (-qJD(1) * t464 - t470) * t429) * MDP(22) + (-t493 - t312 * t448 + t358 * t378 + t484 * t425 - t525 * t391 - t528 * t362 + (-t300 * t426 + t301 * t531) * qJD(1) + t523) * MDP(23) + (t494 - t313 * t448 - t358 * t377 + t484 * t428 + t524 * t391 + t528 * t360 + (t297 * t426 - t301 * t536) * qJD(1) + t522) * MDP(24) + (-t289 * t448 + t288 * t377 + t290 * t378 - g(1) * (pkin(4) * t530 + t399) - g(2) * (pkin(4) * t532 + t395) - g(3) * (qJ(5) * t531 + t475) + t528 * t301 + t525 * t300 + t524 * t297 + (-g(3) * pkin(4) + t470 * t448) * t426) * MDP(25) + (-t426 * t429 * MDP(4) + t517 * MDP(5)) * t432; MDP(11) * t542 + (t357 - t564) * MDP(12) + t299 * MDP(13) + t570 * MDP(14) + t358 * MDP(15) + (-t362 * t375 + t440 - t543) * MDP(16) + (t319 * t391 + t360 * t375 + t436) * MDP(17) + (pkin(3) * t312 - t548 + (-t309 - t320) * t362 + (t308 + t506) * t360) * MDP(18) + (t323 * t360 + t435 + t543 - 0.2e1 * t561) * MDP(19) + (-t311 * t360 + t323 * t362 + t391 * t506 - t379 - t436 + t563) * MDP(20) + (-t293 * pkin(3) - g(1) * t482 - g(2) * t483 + g(3) * t520 - t291 * qJ(4) - t308 * t320 + t309 * t506 - t311 * t323) * MDP(21) + (-t548 + t312 * t550 + (t300 + t504) * t362 + (t297 - t505) * t360) * MDP(22) + (t306 * t362 - t391 * t456 - 0.2e1 * t379 + t433 + t563) * MDP(23) + (-t306 * t360 + (-0.2e1 * qJD(5) - t305) * t391 + 0.2e1 * t490 - t434) * MDP(24) + (-t288 * t550 + t290 * qJ(4) - t301 * t306 - g(1) * (-qJ(5) * t346 + t482) - g(2) * (-qJ(5) * t344 + t483) - g(3) * t459 + t505 * t300 + t504 * t297) * MDP(25); (t435 - t545 - t561) * MDP(21) + ((qJD(5) + t300) * t391 - t490 + t434) * MDP(25) + (MDP(19) - MDP(24)) * (t358 - t542) + (MDP(18) + MDP(22)) * t299 + (MDP(20) + MDP(23)) * (-t357 - t387); (t358 + t542) * MDP(23) + (-t387 - t564) * MDP(24) + (-t297 * t391 + t433 - t571) * MDP(25) + t570 * MDP(22);];
tau = t1;
