% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:10:27
% EndTime: 2019-03-09 05:10:41
% DurationCPUTime: 8.19s
% Computational Cost: add. (8669->435), mult. (23006->588), div. (0->0), fcn. (18570->10), ass. (0->192)
t500 = cos(pkin(10));
t505 = cos(qJ(3));
t564 = t500 * t505;
t498 = sin(pkin(10));
t503 = sin(qJ(3));
t565 = t498 * t503;
t524 = -t564 + t565;
t461 = t524 * qJD(1);
t473 = t498 * t505 + t500 * t503;
t462 = t473 * qJD(1);
t502 = sin(qJ(4));
t588 = cos(qJ(4));
t443 = t588 * t461 + t462 * t502;
t441 = qJD(6) + t443;
t499 = cos(pkin(11));
t504 = cos(qJ(6));
t563 = t504 * t499;
t497 = sin(pkin(11));
t501 = sin(qJ(6));
t566 = t497 * t501;
t470 = -t563 + t566;
t622 = t441 * t470;
t545 = qJD(1) * qJD(3);
t539 = t505 * t545;
t484 = t500 * t539;
t540 = t503 * t545;
t458 = -t498 * t540 + t484;
t521 = -t502 * t461 + t462 * t588;
t556 = t498 * t539 + t500 * t540;
t408 = qJD(4) * t521 + t502 * t458 + t588 * t556;
t472 = t497 * t504 + t499 * t501;
t624 = t472 * t408 - t622 * t441;
t623 = t441 * t472;
t496 = qJD(3) + qJD(4);
t433 = t496 * t497 + t499 * t521;
t541 = qJD(4) * t588;
t552 = qJD(4) * t502;
t407 = t588 * t458 - t461 * t541 - t462 * t552 - t502 * t556;
t581 = t407 * t497;
t523 = -qJD(6) * t433 - t581;
t431 = -t499 * t496 + t497 * t521;
t548 = qJD(6) * t504;
t557 = t407 * t563 - t431 * t548;
t334 = t523 * t501 + t557;
t604 = t431 * t501 - t433 * t504;
t335 = -qJD(6) * t604 + t407 * t472;
t608 = t504 * t431;
t384 = t433 * t501 + t608;
t531 = -t470 * t408 - t623 * t441;
t575 = t443 * t496;
t578 = t521 * t496;
t582 = t604 * t521;
t583 = t384 * t521;
t621 = (t531 + t583) * MDP(29) + (-t408 + t578) * MDP(18) - t443 ^ 2 * MDP(16) + (t443 * MDP(15) + MDP(16) * t521 - MDP(30) * t441) * t521 + (t407 + t575) * MDP(17) + t334 * t472 * MDP(26) + (t582 + t624) * MDP(28) + (-t334 * t470 - t472 * t335 + t622 * t384) * MDP(27) + (t622 * MDP(26) + t623 * MDP(27)) * t604;
t612 = t443 * t497;
t620 = pkin(5) * t612;
t619 = pkin(9) * t612;
t611 = t443 * t499;
t617 = pkin(5) * t521 + pkin(9) * t611;
t586 = pkin(7) + qJ(2);
t479 = t586 * t498;
t474 = qJD(1) * t479;
t481 = t586 * t500;
t475 = qJD(1) * t481;
t596 = -t505 * t474 - t475 * t503;
t416 = -t556 * pkin(8) - qJD(2) * t461 + qJD(3) * t596;
t515 = t473 * qJD(2);
t513 = qJD(1) * t515;
t526 = t474 * t503 - t475 * t505;
t417 = -pkin(8) * t458 + qJD(3) * t526 - t513;
t435 = -pkin(8) * t462 + t596;
t430 = qJD(3) * pkin(3) + t435;
t436 = -pkin(8) * t461 - t526;
t508 = t416 * t588 + t502 * t417 + t430 * t541 - t436 * t552;
t339 = t496 * qJD(5) + t508;
t349 = t556 * pkin(3) + t408 * pkin(4) - t407 * qJ(5) - qJD(5) * t521;
t322 = t499 * t339 + t497 * t349;
t319 = t322 * t499;
t429 = t588 * t436;
t382 = t502 * t430 + t429;
t376 = qJ(5) * t496 + t382;
t491 = -pkin(2) * t500 - pkin(1);
t476 = qJD(1) * t491 + qJD(2);
t449 = pkin(3) * t461 + t476;
t383 = pkin(4) * t443 - qJ(5) * t521 + t449;
t350 = -t376 * t497 + t499 * t383;
t616 = -t350 * t611 + t319;
t328 = pkin(5) * t443 - pkin(9) * t433 + t350;
t351 = t499 * t376 + t497 * t383;
t336 = -pkin(9) * t431 + t351;
t316 = t328 * t501 + t336 * t504;
t341 = t502 * t416 - t588 * t417 + t430 * t552 + t436 * t541;
t327 = pkin(5) * t581 + t341;
t428 = t502 * t436;
t381 = t430 * t588 - t428;
t375 = -t496 * pkin(4) + qJD(5) - t381;
t367 = t431 * pkin(5) + t375;
t615 = t316 * t521 + t327 * t472 - t622 * t367;
t315 = t328 * t504 - t336 * t501;
t614 = -t315 * t521 + t327 * t470 + t623 * t367;
t409 = pkin(4) * t521 + qJ(5) * t443;
t321 = -t339 * t497 + t499 * t349;
t606 = -t351 * t443 - t321;
t605 = t449 * t443 - t508;
t388 = t435 * t588 - t428;
t391 = pkin(3) * t462 + t409;
t356 = -t388 * t497 + t499 * t391;
t485 = pkin(3) * t541 + qJD(5);
t598 = -t485 * t497 - t356;
t357 = t499 * t388 + t497 * t391;
t537 = t485 * t499 - t357;
t362 = t499 * t381 + t497 * t409;
t597 = -qJD(5) * t499 + t362;
t387 = t502 * t435 + t429;
t532 = pkin(3) * t552 - t387;
t569 = t479 * t505;
t438 = -pkin(8) * t473 - t481 * t503 - t569;
t525 = t479 * t503 - t481 * t505;
t439 = -pkin(8) * t524 - t525;
t595 = t588 * t438 - t502 * t439;
t594 = (t498 ^ 2 + t500 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t591 = t341 * t497 + t351 * t521;
t590 = -t521 * t449 - t341;
t589 = -t341 * t499 - t350 * t521;
t587 = t499 * pkin(5);
t493 = t499 * pkin(9);
t580 = t407 * t499;
t465 = t524 * qJD(3);
t466 = t473 * qJD(3);
t520 = -t502 * t473 - t524 * t588;
t414 = t520 * qJD(4) - t465 * t588 - t502 * t466;
t579 = t414 * t497;
t448 = t473 * t588 - t502 * t524;
t572 = t448 * t497;
t571 = t448 * t499;
t511 = -qJD(3) * t569 + qJD(2) * t564 + (-qJD(2) * t498 - qJD(3) * t481) * t503;
t418 = -pkin(8) * t466 + t511;
t507 = qJD(3) * t525 - t515;
t419 = pkin(8) * t465 + t507;
t354 = qJD(4) * t595 + t588 * t418 + t502 * t419;
t415 = t448 * qJD(4) - t502 * t465 + t466 * t588;
t360 = pkin(3) * t466 + pkin(4) * t415 - qJ(5) * t414 - qJD(5) * t448;
t324 = t499 * t354 + t497 * t360;
t453 = pkin(3) * t524 + t491;
t398 = -pkin(4) * t520 - qJ(5) * t448 + t453;
t404 = t502 * t438 + t439 * t588;
t364 = t497 * t398 + t499 * t404;
t553 = qJD(3) * t462;
t550 = qJD(6) * t336;
t549 = qJD(6) * t448;
t546 = qJD(1) * qJD(2);
t542 = qJD(1) * t565;
t323 = -t354 * t497 + t499 * t360;
t361 = -t381 * t497 + t499 * t409;
t363 = t499 * t398 - t404 * t497;
t317 = -pkin(9) * t581 + t322;
t535 = -qJD(6) * t328 - t317;
t492 = -pkin(3) * t588 - pkin(4);
t533 = t532 + t620;
t530 = -t321 * t497 + t319;
t342 = -pkin(5) * t520 - pkin(9) * t571 + t363;
t353 = -pkin(9) * t572 + t364;
t529 = t342 * t504 - t353 * t501;
t528 = t342 * t501 + t353 * t504;
t527 = t350 * t497 - t351 * t499;
t489 = pkin(3) * t502 + qJ(5);
t468 = t489 * t499 + t493;
t519 = qJD(6) * t468 - t598 + t617;
t467 = (-pkin(9) - t489) * t497;
t518 = -qJD(6) * t467 - t537 + t619;
t480 = qJ(5) * t499 + t493;
t517 = qJD(5) * t497 + qJD(6) * t480 + t361 + t617;
t478 = (-pkin(9) - qJ(5)) * t497;
t516 = -qJD(6) * t478 + t597 + t619;
t514 = t341 * t448 + t375 * t414 - t407 * t595;
t512 = -pkin(4) * t407 - qJ(5) * t408 + (-qJD(5) + t375) * t443;
t510 = t407 * t492 - t408 * t489 + (t375 - t485) * t443;
t355 = t404 * qJD(4) + t502 * t418 - t419 * t588;
t490 = -pkin(4) - t587;
t477 = t492 - t587;
t411 = t470 * t448;
t410 = t472 * t448;
t370 = pkin(5) * t572 - t595;
t368 = t382 - t620;
t345 = t414 * t472 + t548 * t571 - t549 * t566;
t344 = -t414 * t470 - t472 * t549;
t332 = pkin(5) * t579 + t355;
t320 = -pkin(9) * t579 + t324;
t318 = pkin(5) * t415 - t414 * t493 + t323;
t314 = pkin(5) * t408 - pkin(9) * t580 + t321;
t313 = t504 * t314;
t1 = [(-t323 * t433 - t324 * t431 + (-t321 * t448 - t350 * t414 - t363 * t407) * t499 + (-t322 * t448 - t351 * t414 - t364 * t407) * t497) * MDP(24) + ((t318 * t504 - t320 * t501) * t441 + t529 * t408 - (-t317 * t501 + t313) * t520 + t315 * t415 + t332 * t384 + t370 * t335 + t327 * t410 + t367 * t345 + (t316 * t520 - t441 * t528) * qJD(6)) * MDP(31) + (-(t318 * t501 + t320 * t504) * t441 - t528 * t408 + (t314 * t501 + t317 * t504) * t520 - t316 * t415 - t332 * t604 + t370 * t334 - t327 * t411 + t367 * t344 + (t315 * t520 - t441 * t529) * qJD(6)) * MDP(32) + (t476 * t466 + t491 * t556) * MDP(13) + (t453 * t407 + t449 * t414 + (t448 * t556 + t466 * t521) * pkin(3)) * MDP(21) + (t453 * t408 + t449 * t415 + (t466 * t443 - t520 * t556) * pkin(3)) * MDP(20) + (-t458 * t524 + t465 * t461 - t462 * t466 - t473 * t556) * MDP(9) + (t322 * t520 - t324 * t443 - t351 * t415 + t355 * t433 - t364 * t408 + t499 * t514) * MDP(23) + (-t321 * t520 + t323 * t443 + t350 * t415 + t355 * t431 + t363 * t408 + t497 * t514) * MDP(22) + (t458 * t473 - t462 * t465) * MDP(8) + (t335 * t520 - t345 * t441 - t384 * t415 - t408 * t410) * MDP(29) + (-t334 * t520 + t344 * t441 - t408 * t411 - t415 * t604) * MDP(28) + (-t408 * t520 + t415 * t441) * MDP(30) + (t407 * t520 - t408 * t448 - t414 * t443 - t415 * t521) * MDP(16) + (t407 * t448 + t414 * t521) * MDP(15) + (-t334 * t410 + t335 * t411 - t344 * t384 + t345 * t604) * MDP(27) + (-t334 * t411 - t344 * t604) * MDP(26) + (t321 * t363 + t322 * t364 + t323 * t350 + t324 * t351 - t341 * t595 + t355 * t375) * MDP(25) + (t491 * t458 - t476 * t465) * MDP(14) + 0.2e1 * t546 * t594 + (t414 * MDP(17) - t415 * MDP(18) - t355 * MDP(20) - t354 * MDP(21)) * t496 + (-t465 * MDP(10) - t466 * MDP(11) + MDP(13) * t507 - MDP(14) * t511) * qJD(3); (t553 + t556) * MDP(13) + (t484 + (-t461 - t542) * qJD(3)) * MDP(14) + (t408 + t578) * MDP(20) + (t407 - t575) * MDP(21) + (t408 * t499 - t431 * t521 - t443 * t612) * MDP(22) + (-t408 * t497 - t433 * t521 - t443 * t611) * MDP(23) + (-(t431 * t499 - t433 * t497) * t443 + (-t497 ^ 2 - t499 ^ 2) * t407) * MDP(24) + (t321 * t499 + t322 * t497 - t375 * t521 - t443 * t527) * MDP(25) + (t531 - t583) * MDP(31) + (t582 - t624) * MDP(32) - qJD(1) ^ 2 * t594; ((t467 * t504 - t468 * t501) * t408 + t477 * t335 + (t501 * t518 - t504 * t519) * t441 + t533 * t384 + t614) * MDP(31) + (t356 * t433 - t537 * t431 + (t433 * t485 + t606) * t497 + t616) * MDP(24) + (-(t467 * t501 + t468 * t504) * t408 + t477 * t334 + (t501 * t519 + t504 * t518) * t441 - t533 * t604 + t615) * MDP(32) + (t484 + (t461 - t542) * qJD(3)) * MDP(10) + (t476 * t461 + t524 * t546) * MDP(14) + (t388 * t496 + (-t462 * t521 - t496 * t541) * pkin(3) + t605) * MDP(21) + (t341 * t492 + t350 * t598 + t537 * t351 + t532 * t375 + t530 * t489) * MDP(25) + (-t476 * t462 - t513) * MDP(13) + (-t356 * t443 + t431 * t532 + t497 * t510 + t589) * MDP(22) + (t387 * t496 + (-t443 * t462 - t496 * t552) * pkin(3) + t590) * MDP(20) + (t357 * t443 + t433 * t532 + t499 * t510 + t591) * MDP(23) + (t553 - t556) * MDP(11) + (-t461 ^ 2 + t462 ^ 2) * MDP(9) + t462 * t461 * MDP(8) + t621; (t382 * t496 + t590) * MDP(20) + (t381 * t496 + t605) * MDP(21) + (-t361 * t443 - t382 * t431 + t497 * t512 + t589) * MDP(22) + (t362 * t443 - t382 * t433 + t499 * t512 + t591) * MDP(23) + (t361 * t433 + t597 * t431 + (qJD(5) * t433 + t606) * t497 + t616) * MDP(24) + (-pkin(4) * t341 + qJ(5) * t530 - qJD(5) * t527 - t350 * t361 - t351 * t362 - t375 * t382) * MDP(25) + ((t478 * t504 - t480 * t501) * t408 + t490 * t335 - t368 * t384 + (t501 * t516 - t504 * t517) * t441 + t614) * MDP(31) + (-(t478 * t501 + t480 * t504) * t408 + t490 * t334 + t368 * t604 + (t501 * t517 + t504 * t516) * t441 + t615) * MDP(32) + t621; (t433 * t443 + t581) * MDP(22) + (-t431 * t443 + t580) * MDP(23) + (-t431 ^ 2 - t433 ^ 2) * MDP(24) + (t350 * t433 + t351 * t431 + t341) * MDP(25) + (-t441 * t604 + t335) * MDP(31) + (-t441 * t608 + (-t433 * t441 + t523) * t501 + t557) * MDP(32); -t384 ^ 2 * MDP(27) + (t384 * t441 + t557) * MDP(28) + t408 * MDP(30) + (t316 * t441 + t313) * MDP(31) + (t315 * t441 + t367 * t384) * MDP(32) - (MDP(26) * t384 - MDP(27) * t604 + t441 * MDP(29) - t367 * MDP(31)) * t604 + (MDP(29) * t523 - MDP(31) * t550 + MDP(32) * t535) * t504 + (t523 * MDP(28) + (qJD(6) * t431 - t580) * MDP(29) + t535 * MDP(31) + (-t314 + t550) * MDP(32)) * t501;];
tauc  = t1;
