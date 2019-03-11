% Calculate vector of inverse dynamics joint torques for
% S6PRPRPR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRPRPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:41:30
% EndTime: 2019-03-08 19:41:42
% DurationCPUTime: 9.82s
% Computational Cost: add. (4555->503), mult. (10815->678), div. (0->0), fcn. (9090->18), ass. (0->222)
t498 = cos(pkin(11));
t615 = cos(qJ(4));
t561 = t615 * t498;
t478 = qJD(2) * t561;
t501 = sin(qJ(4));
t495 = sin(pkin(11));
t575 = qJD(2) * t495;
t558 = t501 * t575;
t446 = -t478 + t558;
t438 = qJD(6) + t446;
t457 = t495 * t615 + t501 * t498;
t448 = t457 * qJD(2);
t494 = sin(pkin(12));
t497 = cos(pkin(12));
t426 = -t497 * qJD(4) + t448 * t494;
t503 = cos(qJ(6));
t428 = qJD(4) * t494 + t448 * t497;
t500 = sin(qJ(6));
t598 = t428 * t500;
t627 = -t503 * t426 - t598;
t630 = t627 * t438;
t521 = -t501 * t495 + t561;
t451 = t521 * qJD(4);
t452 = t457 * qJD(4);
t502 = sin(qJ(2));
t496 = sin(pkin(6));
t576 = qJD(1) * t496;
t560 = t502 * t576;
t629 = pkin(4) * t452 - qJ(5) * t451 - qJD(5) * t457 - t560;
t504 = cos(qJ(2));
t588 = t496 * t504;
t511 = t521 * t588;
t610 = pkin(8) + qJ(3);
t463 = t610 * t495;
t465 = t610 * t498;
t522 = -t463 * t615 - t501 * t465;
t628 = -qJD(1) * t511 + qJD(3) * t521 + qJD(4) * t522;
t459 = qJD(2) * qJ(3) + t560;
t608 = cos(pkin(6));
t551 = qJD(1) * t608;
t474 = t498 * t551;
t417 = t474 + (-pkin(8) * qJD(2) - t459) * t495;
t432 = t498 * t459 + t495 * t551;
t574 = qJD(2) * t498;
t418 = pkin(8) * t574 + t432;
t360 = t501 * t417 + t615 * t418;
t626 = t360 * qJD(4);
t533 = t426 * t500 - t428 * t503;
t625 = t438 * t533;
t456 = t494 * t503 + t497 * t500;
t450 = t456 * qJD(6);
t578 = t456 * t446 + t450;
t590 = t495 * MDP(6);
t623 = t498 * MDP(5) - t590;
t577 = t495 ^ 2 + t498 ^ 2;
t622 = MDP(7) * t577;
t493 = pkin(11) + qJ(4);
t487 = sin(t493);
t607 = cos(pkin(10));
t544 = t608 * t607;
t606 = sin(pkin(10));
t442 = t502 * t606 - t504 * t544;
t543 = t608 * t606;
t444 = t502 * t607 + t504 * t543;
t547 = g(1) * t444 + g(2) * t442;
t620 = g(3) * t588 - t547;
t621 = t620 * t487;
t583 = -t628 * t494 + t629 * t497;
t582 = t629 * t494 + t628 * t497;
t423 = -t501 * t463 + t465 * t615;
t512 = t457 * t588;
t580 = -qJD(1) * t512 + qJD(3) * t457 + qJD(4) * t423;
t524 = -t615 * t417 + t501 * t418;
t559 = t504 * t576;
t542 = qJD(3) - t559;
t573 = qJD(2) * t502;
t556 = qJD(1) * t573;
t619 = t496 * t556 + qJDD(3);
t618 = -qJDD(4) * pkin(4) + qJDD(5);
t554 = qJDD(2) * t615;
t566 = qJDD(2) * t501;
t541 = t495 * t566 - t498 * t554;
t406 = qJD(2) * t452 + t541;
t401 = qJDD(6) + t406;
t454 = t494 * t500 - t503 * t497;
t579 = t438 * t454;
t617 = -t401 * t456 + t438 * t579;
t439 = t446 ^ 2;
t614 = pkin(9) * t497;
t611 = g(3) * t496;
t609 = pkin(9) + qJ(5);
t605 = qJDD(2) * pkin(2);
t603 = t627 * t448;
t602 = t533 * t448;
t600 = t406 * t494;
t599 = t406 * t497;
t597 = t446 * t494;
t596 = t451 * t494;
t595 = t457 * t494;
t594 = t457 * t497;
t492 = pkin(12) + qJ(6);
t486 = sin(t492);
t489 = cos(t493);
t593 = t486 * t489;
t488 = cos(t492);
t592 = t488 * t489;
t591 = t489 * t504;
t589 = t496 * t502;
t586 = qJDD(1) - g(3);
t568 = qJDD(2) * qJ(3);
t569 = qJDD(1) * t496;
t430 = t502 * t569 + t568 + (qJD(3) + t559) * qJD(2);
t548 = qJDD(1) * t608;
t472 = t498 * t548;
t386 = t472 + (-pkin(8) * qJDD(2) - t430) * t495;
t399 = t498 * t430 + t495 * t548;
t567 = qJDD(2) * t498;
t387 = pkin(8) * t567 + t399;
t525 = t501 * t386 + t387 * t615;
t325 = qJDD(4) * qJ(5) + (qJD(5) - t524) * qJD(4) + t525;
t563 = qJD(4) * t478 + t495 * t554 + t498 * t566;
t405 = -qJD(4) * t558 + t563;
t482 = pkin(3) * t498 + pkin(2);
t555 = t504 * t569;
t526 = -t555 + t619;
t420 = -qJDD(2) * t482 + t526;
t338 = pkin(4) * t406 - qJ(5) * t405 - qJD(5) * t448 + t420;
t321 = t497 * t325 + t494 * t338;
t585 = pkin(5) * t452 - t451 * t614 + t583;
t584 = pkin(9) * t596 - t582;
t356 = qJD(4) * qJ(5) + t360;
t437 = -qJD(2) * t482 + t542;
t369 = pkin(4) * t446 - qJ(5) * t448 + t437;
t333 = t497 * t356 + t494 * t369;
t400 = pkin(4) * t448 + qJ(5) * t446;
t340 = t494 * t400 - t497 * t524;
t581 = pkin(5) * t596 + t580;
t402 = -pkin(4) * t521 - qJ(5) * t457 - t482;
t354 = t494 * t402 + t497 * t423;
t572 = qJD(6) * t500;
t571 = qJD(6) * t503;
t352 = -qJD(4) * pkin(4) + qJD(5) + t524;
t570 = -qJD(5) + t352;
t565 = g(3) * t589;
t384 = -t497 * qJDD(4) + t405 * t494;
t385 = qJDD(4) * t494 + t405 * t497;
t564 = -t500 * t384 + t503 * t385 - t426 * t571;
t557 = t496 * t573;
t553 = t496 * t607;
t552 = t496 * t606;
t320 = -t325 * t494 + t497 * t338;
t316 = pkin(5) * t406 - pkin(9) * t385 + t320;
t319 = -pkin(9) * t384 + t321;
t550 = t503 * t316 - t500 * t319;
t332 = -t356 * t494 + t497 * t369;
t339 = t497 * t400 + t494 * t524;
t549 = t503 * t384 + t500 * t385;
t353 = t497 * t402 - t423 * t494;
t443 = t502 * t544 + t504 * t606;
t445 = -t502 * t543 + t504 * t607;
t546 = g(1) * t445 + g(2) * t443;
t545 = -t454 * t401 - t578 * t438;
t540 = t500 * t316 + t503 * t319;
t539 = -t320 * t497 - t321 * t494;
t323 = pkin(5) * t446 - pkin(9) * t428 + t332;
t327 = -pkin(9) * t426 + t333;
t317 = t323 * t503 - t327 * t500;
t318 = t323 * t500 + t327 * t503;
t538 = -t332 * t494 + t333 * t497;
t343 = -pkin(5) * t521 - pkin(9) * t594 + t353;
t345 = -pkin(9) * t595 + t354;
t537 = t343 * t503 - t345 * t500;
t536 = t343 * t500 + t345 * t503;
t440 = -t495 * t589 + t498 * t608;
t441 = t495 * t608 + t498 * t589;
t391 = t501 * t440 + t441 * t615;
t370 = -t391 * t494 - t497 * t588;
t371 = t391 * t497 - t494 * t588;
t535 = t370 * t503 - t371 * t500;
t534 = t370 * t500 + t371 * t503;
t431 = -t459 * t495 + t474;
t532 = t431 * t495 - t432 * t498;
t433 = t526 - t605;
t531 = -t433 + t547;
t527 = MDP(3) + t623;
t523 = t440 * t615 - t501 * t441;
t464 = t609 * t497;
t520 = pkin(5) * t448 + qJD(5) * t494 + qJD(6) * t464 + t446 * t614 + t339;
t462 = t609 * t494;
t519 = pkin(9) * t597 - qJD(5) * t497 + qJD(6) * t462 + t340;
t518 = -t386 * t615 + t501 * t387 + t626;
t328 = -t428 * t572 + t564;
t517 = g(1) * (t445 * t487 - t489 * t552) + g(2) * (t443 * t487 + t489 * t553) + g(3) * (t487 * t589 - t489 * t608);
t409 = t443 * t489 - t487 * t553;
t411 = t445 * t489 + t487 * t552;
t436 = t487 * t608 + t489 * t589;
t516 = g(1) * t411 + g(2) * t409 + g(3) * t436;
t326 = t518 + t618;
t514 = -t326 + t517;
t510 = t620 * t489;
t509 = -t620 + t555;
t508 = t326 * t457 + t352 * t451 - t546;
t398 = -t430 * t495 + t472;
t507 = -t398 * t495 + t399 * t498 - t546;
t329 = -qJD(6) * t533 + t549;
t506 = t517 - t518;
t505 = qJD(2) ^ 2;
t483 = -pkin(5) * t497 - pkin(4);
t458 = -qJD(2) * pkin(2) + t542;
t395 = t454 * t457;
t394 = t456 * t457;
t381 = pkin(5) * t595 - t522;
t358 = qJD(2) * t512 + qJD(4) * t391;
t357 = qJD(2) * t511 + qJD(4) * t523;
t350 = t357 * t497 + t494 * t557;
t349 = -t357 * t494 + t497 * t557;
t348 = t451 * t456 + t571 * t594 - t572 * t595;
t347 = -t450 * t457 - t451 * t454;
t346 = -pkin(5) * t597 + t360;
t344 = t426 * pkin(5) + t352;
t322 = t384 * pkin(5) + t326;
t1 = [t586 * MDP(1) + (t398 * t440 + t399 * t441 - g(3)) * MDP(8) + (-qJD(4) * t358 + qJDD(4) * t523) * MDP(14) + (-qJD(4) * t357 - qJDD(4) * t391) * MDP(15) + (t349 * t446 + t358 * t426 + t370 * t406 - t384 * t523) * MDP(16) + (-t350 * t446 + t358 * t428 - t371 * t406 - t385 * t523) * MDP(17) + (-t349 * t428 - t350 * t426 - t370 * t385 - t371 * t384) * MDP(18) + (t320 * t370 + t321 * t371 - t326 * t523 + t332 * t349 + t333 * t350 + t352 * t358 - g(3)) * MDP(19) + ((-qJD(6) * t534 + t349 * t503 - t350 * t500) * t438 + t535 * t401 - t358 * t627 - t523 * t329) * MDP(25) + (-(qJD(6) * t535 + t349 * t500 + t350 * t503) * t438 - t534 * t401 - t358 * t533 - t523 * t328) * MDP(26) + (-t440 * t495 + t441 * t498) * MDP(7) * qJDD(2) + ((-qJDD(2) * MDP(4) - t527 * t505 + (t446 * MDP(14) + t448 * MDP(15) + t458 * MDP(8)) * qJD(2)) * t502 + ((-t431 * t575 + t432 * t574 - t433) * MDP(8) - t406 * MDP(14) - t405 * MDP(15) + (-MDP(4) + t622) * t505 + t527 * qJDD(2)) * t504) * t496; qJDD(2) * MDP(2) + t509 * MDP(3) + (-t586 * t589 + t546) * MDP(4) + (-t565 + t507 + (qJD(2) * t542 + t568) * t577) * MDP(7) + (-t532 * qJD(3) + t531 * pkin(2) + t507 * qJ(3) + (-g(3) * (pkin(2) * t504 + qJ(3) * t502) + (-t458 * t502 + t504 * t532) * qJD(1)) * t496) * MDP(8) + (t405 * t457 + t448 * t451) * MDP(9) + (t405 * t521 - t406 * t457 - t446 * t451 - t448 * t452) * MDP(10) + (qJD(4) * t451 + qJDD(4) * t457) * MDP(11) + (-qJD(4) * t452 + qJDD(4) * t521) * MDP(12) + (-qJD(4) * t580 + qJDD(4) * t522 - t406 * t482 - t420 * t521 + t437 * t452 - t446 * t560 - t510) * MDP(14) + (-t628 * qJD(4) - qJDD(4) * t423 - t405 * t482 + t420 * t457 + t437 * t451 - t448 * t560 + t621) * MDP(15) + (-t320 * t521 + t332 * t452 + t353 * t406 - t522 * t384 - t497 * t510 + (t508 - t565) * t494 + t583 * t446 + t580 * t426) * MDP(16) + (t321 * t521 - t333 * t452 - t354 * t406 - t522 * t385 - t547 * t494 * t489 + t508 * t497 - (-t494 * t591 + t497 * t502) * t611 - t582 * t446 + t580 * t428) * MDP(17) + (-t353 * t385 - t354 * t384 + t539 * t457 + (-t332 * t497 - t333 * t494) * t451 - t583 * t428 - t582 * t426 - t621) * MDP(18) + (t320 * t353 + t321 * t354 - t326 * t522 + t583 * t332 + t582 * t333 + t580 * t352 - (t502 * t611 + t546) * t610 + (-t504 * t611 + t547) * (pkin(4) * t489 + qJ(5) * t487 + t482)) * MDP(19) + (-t328 * t395 - t347 * t533) * MDP(20) + (-t328 * t394 + t329 * t395 + t347 * t627 + t348 * t533) * MDP(21) + (-t328 * t521 + t347 * t438 - t395 * t401 - t452 * t533) * MDP(22) + (t329 * t521 - t348 * t438 - t394 * t401 + t452 * t627) * MDP(23) + (-t401 * t521 + t438 * t452) * MDP(24) + (t537 * t401 - t550 * t521 + t317 * t452 + t381 * t329 + t322 * t394 + t344 * t348 - g(1) * (-t444 * t592 + t445 * t486) - g(2) * (-t442 * t592 + t443 * t486) - (t486 * t502 + t488 * t591) * t611 + (t500 * t584 + t503 * t585) * t438 - t581 * t627 + (t318 * t521 - t438 * t536) * qJD(6)) * MDP(25) + (-t536 * t401 + t540 * t521 - t318 * t452 + t381 * t328 - t322 * t395 + t344 * t347 - g(1) * (t444 * t593 + t445 * t488) - g(2) * (t442 * t593 + t443 * t488) - (-t486 * t591 + t488 * t502) * t611 + (-t500 * t585 + t503 * t584) * t438 - t581 * t533 + (t317 * t521 - t438 * t537) * qJD(6)) * MDP(26) + t623 * (t496 * (-g(3) * t504 + t556) + t531 + t605); -MDP(5) * t567 + qJDD(2) * t590 - t505 * t622 + (qJD(2) * t532 - t509 - t605 + t619) * MDP(8) + (0.2e1 * qJD(4) * t448 + t541) * MDP(14) + ((-t446 - t558) * qJD(4) + t563) * MDP(15) + (-t426 * t448 - t439 * t494 + t599) * MDP(16) + (-t428 * t448 - t439 * t497 - t600) * MDP(17) + (-t384 * t494 - t385 * t497 + (-t426 * t497 + t428 * t494) * t446) * MDP(18) + (-t352 * t448 + t446 * t538 - t539 + t620) * MDP(19) + (t545 + t603) * MDP(25) + (t602 + t617) * MDP(26); -t439 * MDP(10) + ((t446 - t558) * qJD(4) + t563) * MDP(11) - t541 * MDP(12) + qJDD(4) * MDP(13) + (t506 + t626) * MDP(14) + (t437 * t446 + t516 - t525) * MDP(15) + (-qJ(5) * t600 - pkin(4) * t384 - t360 * t426 + (t494 * t570 - t339) * t446 + t514 * t497) * MDP(16) + (-qJ(5) * t599 - pkin(4) * t385 - t360 * t428 + (t497 * t570 + t340) * t446 - t514 * t494) * MDP(17) + (t339 * t428 + t340 * t426 + (-qJ(5) * t384 - qJD(5) * t426 - t332 * t446 + t321) * t497 + (qJ(5) * t385 + qJD(5) * t428 - t333 * t446 - t320) * t494 - t516) * MDP(18) + (-t332 * t339 - t333 * t340 - t352 * t360 + t538 * qJD(5) + t514 * pkin(4) + (-t320 * t494 + t321 * t497 - t516) * qJ(5)) * MDP(19) + (t328 * t456 + t533 * t579) * MDP(20) + (-t328 * t454 - t329 * t456 + t533 * t578 - t579 * t627) * MDP(21) + (t602 - t617) * MDP(22) + (t545 - t603) * MDP(23) + ((-t462 * t503 - t464 * t500) * t401 + t483 * t329 + t322 * t454 + t346 * t627 + (t500 * t519 - t503 * t520) * t438 + t578 * t344 + t517 * t488) * MDP(25) + (-(-t462 * t500 + t464 * t503) * t401 + t483 * t328 + t322 * t456 + t346 * t533 + (t500 * t520 + t503 * t519) * t438 - t579 * t344 - t517 * t486) * MDP(26) + (t448 * MDP(10) - t437 * MDP(14) - t332 * MDP(16) + t333 * MDP(17) - t438 * MDP(24) - t317 * MDP(25) + t318 * MDP(26) + MDP(9) * t446) * t448; (t428 * t446 + t384) * MDP(16) + (-t426 * t446 + t385) * MDP(17) + (-t426 ^ 2 - t428 ^ 2) * MDP(18) + (t332 * t428 + t333 * t426 - t506 + t618) * MDP(19) + (t329 - t625) * MDP(25) + (t328 + t630) * MDP(26); t533 * t627 * MDP(20) + (t533 ^ 2 - t627 ^ 2) * MDP(21) + (t564 - t630) * MDP(22) + (-t549 - t625) * MDP(23) + t401 * MDP(24) + (t318 * t438 + t344 * t533 - g(1) * (-t411 * t486 + t444 * t488) - g(2) * (-t409 * t486 + t442 * t488) - g(3) * (-t436 * t486 - t488 * t588) + t550) * MDP(25) + (t317 * t438 - t344 * t627 - g(1) * (-t411 * t488 - t444 * t486) - g(2) * (-t409 * t488 - t442 * t486) - g(3) * (-t436 * t488 + t486 * t588) - t540) * MDP(26) + (-MDP(22) * t598 + MDP(23) * t533 - MDP(25) * t318 - MDP(26) * t317) * qJD(6);];
tau  = t1;
