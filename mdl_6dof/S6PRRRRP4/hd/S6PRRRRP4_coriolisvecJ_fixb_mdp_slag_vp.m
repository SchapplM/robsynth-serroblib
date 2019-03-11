% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:18:06
% EndTime: 2019-03-09 00:18:16
% DurationCPUTime: 6.66s
% Computational Cost: add. (5755->512), mult. (14003->663), div. (0->0), fcn. (10221->10), ass. (0->219)
t498 = sin(qJ(3));
t501 = cos(qJ(3));
t529 = pkin(3) * t498 - pkin(9) * t501;
t454 = t529 * qJD(3);
t460 = -pkin(3) * t501 - pkin(9) * t498 - pkin(2);
t497 = sin(qJ(4));
t499 = sin(qJ(2));
t500 = cos(qJ(4));
t573 = qJD(4) * t500;
t494 = sin(pkin(6));
t586 = qJD(1) * t494;
t502 = cos(qJ(2));
t604 = t501 * t502;
t654 = -(t497 * t499 + t500 * t604) * t586 + t497 * t454 + t460 * t573;
t578 = qJD(3) * t498;
t626 = pkin(8) * t497;
t653 = t500 * t454 + t578 * t626 - (-t497 * t604 + t499 * t500) * t586;
t583 = qJD(2) * t498;
t552 = t497 * t583;
t577 = qJD(3) * t500;
t447 = -t552 + t577;
t579 = qJD(3) * t497;
t448 = t500 * t583 + t579;
t496 = sin(qJ(5));
t627 = cos(qJ(5));
t394 = -t627 * t447 + t448 * t496;
t520 = t496 * t447 + t448 * t627;
t652 = t394 * t520;
t606 = t500 * t501;
t484 = pkin(8) * t606;
t527 = pkin(4) * t498 - pkin(10) * t606;
t512 = t527 * qJD(3);
t651 = t512 + (-t484 + (pkin(10) * t498 - t460) * t497) * qJD(4) + t653;
t575 = qJD(4) * t497;
t553 = t498 * t573;
t576 = qJD(3) * t501;
t556 = t497 * t576;
t644 = t553 + t556;
t650 = -t644 * pkin(10) + (-t498 * t577 - t501 * t575) * pkin(8) + t654;
t562 = t499 * t586;
t457 = qJD(2) * pkin(8) + t562;
t495 = cos(pkin(6));
t585 = qJD(1) * t501;
t424 = -t498 * t457 + t495 * t585;
t451 = t529 * qJD(2);
t543 = -t424 * t497 + t500 * t451;
t628 = -pkin(10) - pkin(9);
t564 = qJD(4) * t628;
t649 = qJD(2) * t527 - t500 * t564 + t543;
t582 = qJD(2) * t501;
t557 = t497 * t582;
t592 = t500 * t424 + t497 * t451;
t648 = pkin(10) * t557 + t497 * t564 - t592;
t574 = qJD(4) * t498;
t554 = t497 * t574;
t511 = t500 * t576 - t554;
t568 = qJD(3) * qJD(4);
t545 = t500 * t568;
t507 = qJD(2) * t511 + t545;
t550 = t627 * qJD(5);
t547 = qJD(2) * t574;
t569 = qJD(2) * qJD(3);
t548 = t501 * t569;
t565 = t500 * t547 + (t548 + t568) * t497;
t572 = qJD(5) * t496;
t348 = -t447 * t550 + t448 * t572 + t496 * t565 - t627 * t507;
t482 = -qJD(4) + t582;
t469 = -qJD(5) + t482;
t343 = -t394 * t469 - t348;
t349 = qJD(5) * t520 + t496 * t507 + t627 * t565;
t549 = t498 * t569;
t629 = t520 ^ 2;
t647 = MDP(23) * t549 + (-t469 * t520 - t349) * MDP(22) + MDP(19) * t652 + (-t394 ^ 2 + t629) * MDP(20) + t343 * MDP(21);
t622 = qJD(3) * pkin(3);
t416 = -t424 - t622;
t387 = -pkin(4) * t447 + t416;
t346 = pkin(5) * t394 - qJ(6) * t520 + t387;
t646 = t346 * t394;
t645 = t387 * t394;
t611 = t496 * t497;
t518 = t500 * t627 - t611;
t632 = qJD(4) + qJD(5);
t633 = t627 * qJD(4) + t550;
t594 = -t500 * t633 + t518 * t582 + t611 * t632;
t450 = t496 * t500 + t497 * t627;
t403 = t632 * t450;
t593 = -t450 * t582 + t403;
t365 = pkin(5) * t520 + qJ(6) * t394;
t640 = MDP(5) * t498;
t492 = t498 ^ 2;
t639 = MDP(6) * (-t501 ^ 2 + t492);
t621 = t346 * t520;
t638 = t387 * t520;
t464 = t628 * t497;
t465 = t628 * t500;
t519 = t464 * t627 + t496 * t465;
t637 = qJD(5) * t519 - t496 * t649 + t627 * t648;
t420 = t496 * t464 - t465 * t627;
t636 = qJD(5) * t420 + t496 * t648 + t627 * t649;
t613 = t495 * t498;
t476 = qJD(1) * t613;
t425 = t501 * t457 + t476;
t417 = qJD(3) * pkin(9) + t425;
t561 = t502 * t586;
t427 = qJD(2) * t460 - t561;
t610 = t497 * t427;
t380 = t500 * t417 + t610;
t364 = pkin(10) * t447 + t380;
t446 = t500 * t460;
t608 = t498 * t500;
t399 = -pkin(10) * t608 + t446 + (-pkin(4) - t626) * t501;
t589 = t497 * t460 + t484;
t609 = t497 * t498;
t408 = -pkin(10) * t609 + t589;
t635 = t399 * t550 - t408 * t572 + t496 * t651 + t627 * t650;
t531 = -t425 + (-t557 + t575) * pkin(4);
t634 = t496 * t399 + t627 * t408;
t631 = qJD(5) * t634 + t496 * t650 - t627 * t651;
t625 = pkin(9) * t482;
t623 = qJD(2) * pkin(2);
t584 = qJD(2) * t494;
t558 = t502 * t584;
t536 = t498 * t558;
t389 = qJD(1) * t536 + qJD(3) * t476 + t457 * t576;
t620 = t389 * t497;
t619 = t389 * t500;
t618 = t416 * t497;
t617 = t448 * t482;
t616 = t494 * t499;
t615 = t494 * t502;
t504 = qJD(2) ^ 2;
t614 = t494 * t504;
t612 = t496 * t364;
t503 = qJD(3) ^ 2;
t607 = t498 * t503;
t605 = t501 * t482;
t603 = t501 * t503;
t602 = qJ(6) * t578 - qJD(6) * t501 + t635;
t601 = -pkin(5) * t578 + t631;
t379 = -t417 * t497 + t500 * t427;
t363 = -pkin(10) * t448 + t379;
t342 = t363 * t627 - t612;
t600 = -pkin(4) * t550 - qJD(6) + t342;
t599 = qJ(6) * t583 - t637;
t598 = pkin(5) * t583 + t636;
t597 = pkin(5) * t593 + qJ(6) * t594 - qJD(6) * t450 + t531;
t455 = pkin(4) * t609 + t498 * pkin(8);
t581 = qJD(3) * t519;
t580 = qJD(3) * t420;
t571 = t492 * qJD(2);
t354 = -pkin(4) * t482 + t363;
t336 = t354 * t627 - t612;
t570 = qJD(6) - t336;
t566 = t499 * t614;
t426 = pkin(4) * t644 + pkin(8) * t576;
t489 = -pkin(4) * t500 - pkin(3);
t563 = t627 * t364;
t559 = t499 * t584;
t555 = t482 * t575;
t546 = MDP(16) * t578;
t388 = -t457 * t578 + (qJD(3) * t495 + t558) * t585;
t423 = (t454 + t562) * qJD(2);
t544 = t497 * t388 - t500 * t423;
t542 = t482 + t582;
t541 = -t447 + t577;
t540 = qJD(4) + t582;
t539 = pkin(5) * t549;
t335 = qJD(2) * t512 - qJD(4) * t364 - t544;
t514 = t500 * t388 - t417 * t575 + t497 * t423 + t427 * t573;
t340 = -pkin(10) * t565 + t514;
t538 = -t496 * t335 - t627 * t340 - t354 * t550 + t364 * t572;
t537 = -t627 * t335 + t496 * t340 + t354 * t572 + t364 * t550;
t535 = t501 * t558;
t534 = t394 * t561;
t533 = t520 * t561;
t532 = t627 * t576;
t341 = t496 * t363 + t563;
t530 = pkin(4) * t572 - t341;
t458 = -t561 - t623;
t528 = -t458 - t561;
t459 = t469 * qJD(6);
t478 = qJ(6) * t549;
t323 = t478 - t459 - t538;
t437 = t501 * t616 + t613;
t406 = -t437 * t497 - t500 * t615;
t526 = -t437 * t500 + t497 * t615;
t436 = -t495 * t501 + t498 * t616;
t524 = t399 * t627 - t496 * t408;
t337 = t496 * t354 + t563;
t522 = t406 * t627 + t496 * t526;
t369 = t496 * t406 - t526 * t627;
t517 = t540 * t579;
t516 = -t336 * t469 + t538;
t515 = -t337 * t469 - t537;
t371 = pkin(4) * t565 + t389;
t324 = t537 - t539;
t510 = qJD(3) * (-t528 - t623);
t488 = -pkin(4) * t627 - pkin(5);
t485 = pkin(4) * t496 + qJ(6);
t431 = t518 * t498;
t430 = t450 * t498;
t405 = qJD(3) * t437 + t536;
t404 = -qJD(3) * t436 + t535;
t392 = -pkin(5) * t518 - qJ(6) * t450 + t489;
t384 = pkin(5) * t430 - qJ(6) * t431 + t455;
t374 = t497 * t532 - t496 * t554 - t572 * t609 + (t496 * t576 + t498 * t633) * t500;
t373 = t403 * t498 + t496 * t556 - t500 * t532;
t362 = t501 * pkin(5) - t524;
t361 = -qJ(6) * t501 + t634;
t359 = qJD(4) * t406 + t404 * t500 + t497 * t559;
t358 = qJD(4) * t526 - t404 * t497 + t500 * t559;
t351 = pkin(4) * t448 + t365;
t332 = -t469 * qJ(6) + t337;
t331 = pkin(5) * t374 + qJ(6) * t373 - qJD(6) * t431 + t426;
t330 = t469 * pkin(5) + t570;
t327 = qJD(5) * t369 - t358 * t627 + t496 * t359;
t326 = qJD(5) * t522 + t496 * t358 + t359 * t627;
t325 = t349 * pkin(5) + t348 * qJ(6) - qJD(6) * t520 + t371;
t1 = [-MDP(3) * t566 - t502 * MDP(4) * t614 + (-t501 * t566 + (-t405 - t536) * qJD(3)) * MDP(10) + (t498 * t566 + (-t404 - t535) * qJD(3)) * MDP(11) + (-t358 * t482 - t405 * t447 + t406 * t549 + t436 * t565) * MDP(17) + (t436 * t545 + t359 * t482 + t405 * t448 + (t436 * t511 + t526 * t578) * qJD(2)) * MDP(18) + (-t326 * t394 + t327 * t520 + t348 * t522 - t349 * t369) * MDP(27) + (t323 * t369 - t324 * t522 + t325 * t436 + t326 * t332 + t327 * t330 + t346 * t405) * MDP(29) + (MDP(24) + MDP(26)) * (t327 * t469 + t436 * t349 + t405 * t394 + t522 * t549) + (-MDP(25) + MDP(28)) * (-t326 * t469 + t348 * t436 + t369 * t549 - t405 * t520); -0.2e1 * t569 * t639 + (-pkin(8) * t603 + t498 * t510) * MDP(10) + (pkin(8) * t607 + t501 * t510) * MDP(11) + (t448 * t511 + t507 * t608) * MDP(12) + ((t447 * t500 - t448 * t497) * t576 + ((-t447 + t552) * t575 + (-t448 * qJD(4) - t517 - t565) * t500) * t498) * MDP(13) + (t542 * t554 + (t448 * t498 + (t571 + (-t482 - t540) * t501) * t500) * qJD(3)) * MDP(14) + (t482 * t553 + t565 * t501 + (t447 * t498 + (-t571 + t605) * t497) * qJD(3)) * MDP(15) - t542 * t546 + ((t460 * t575 - t653) * t482 + ((-pkin(8) * t447 + t618) * qJD(3) + (t610 + (pkin(8) * t482 + t417) * t500) * qJD(4) + t544) * t501 + (pkin(8) * t565 + t620 + t416 * t573 + t447 * t561 + ((-t501 * t626 + t446) * qJD(2) + t379) * qJD(3)) * t498) * MDP(17) + (t654 * t482 + (t416 * t577 + (qJD(3) * t448 - t555) * pkin(8) + t514) * t501 + (-t448 * t561 + t619 + (-pkin(8) * t583 - t416) * t575 + (-t589 * qJD(2) - t380 + (-t482 + t540) * pkin(8) * t500) * qJD(3)) * t498) * MDP(18) + (-t348 * t431 - t373 * t520) * MDP(19) + (t348 * t430 - t349 * t431 + t373 * t394 - t374 * t520) * MDP(20) + (t348 * t501 + t373 * t469 + (qJD(2) * t431 + t520) * t578) * MDP(21) + (t349 * t501 + t374 * t469 + (-qJD(2) * t430 - t394) * t578) * MDP(22) + (-t469 - t582) * MDP(23) * t578 + (t336 * t578 + t455 * t349 + t371 * t430 + t387 * t374 + t426 * t394 + t469 * t631 - t498 * t534 + t537 * t501 + t524 * t549) * MDP(24) + (-t538 * t501 + t426 * t520 - t455 * t348 + t371 * t431 - t387 * t373 + t635 * t469 + (-t533 + (-qJD(2) * t634 - t337) * qJD(3)) * t498) * MDP(25) + (t324 * t501 + t325 * t430 + t331 * t394 + t346 * t374 + t349 * t384 + t601 * t469 + (-t534 + (-qJD(2) * t362 - t330) * qJD(3)) * t498) * MDP(26) + (-t323 * t430 + t324 * t431 - t330 * t373 - t332 * t374 - t348 * t362 - t349 * t361 - t394 * t602 + t520 * t601) * MDP(27) + (-t323 * t501 - t325 * t431 - t331 * t520 + t346 * t373 + t348 * t384 - t602 * t469 + (t533 + (qJD(2) * t361 + t332) * qJD(3)) * t498) * MDP(28) + (t323 * t361 + t324 * t362 + t325 * t384 + (-t498 * t561 + t331) * t346 + t602 * t332 + t601 * t330) * MDP(29) + MDP(7) * t603 - MDP(8) * t607 + 0.2e1 * t548 * t640; (qJD(3) * t425 - t389) * MDP(10) + t528 * t582 * MDP(11) + (-t497 ^ 2 * t547 + (t517 - t617) * t500) * MDP(12) + ((-t565 + t617) * t497 + ((t447 + t577) * qJD(4) + (t501 * t541 - t554) * qJD(2)) * t500) * MDP(13) + (-t482 * t573 + (t500 * t605 + (-t448 + t579) * t498) * qJD(2)) * MDP(14) + (t555 + (-t497 * t605 + t498 * t541) * qJD(2)) * MDP(15) + (-pkin(3) * t565 - t619 + t543 * t482 + t425 * t447 + (t500 * t625 + t618) * qJD(4) + (-t379 * t498 + (-pkin(9) * t578 - t416 * t501) * t497) * qJD(2)) * MDP(17) + (t620 - t592 * t482 - t425 * t448 + (-t497 * t625 + (t416 - t622) * t500) * qJD(4) + ((-t416 - t622) * t606 + (pkin(3) * t575 - pkin(9) * t577 + t380) * t498) * qJD(2)) * MDP(18) + (-t348 * t450 - t520 * t594) * MDP(19) + (-t348 * t518 - t349 * t450 + t394 * t594 - t520 * t593) * MDP(20) + (t489 * t349 - t371 * t518 + t593 * t387 + t531 * t394) * MDP(24) + (-t489 * t348 + t371 * t450 - t594 * t387 + t520 * t531) * MDP(25) + (-t325 * t518 + t593 * t346 + t349 * t392 + t597 * t394) * MDP(26) + (t323 * t518 + t324 * t450 - t330 * t594 - t332 * t593 + t348 * t519 - t349 * t420 + t394 * t599 + t520 * t598) * MDP(27) + (-t325 * t450 + t594 * t346 + t348 * t392 - t520 * t597) * MDP(28) + (t323 * t420 - t324 * t519 + t325 * t392 + t330 * t598 - t332 * t599 + t346 * t597) * MDP(29) + (t594 * MDP(21) + t593 * MDP(22) + MDP(24) * t636 + MDP(25) * t637 + t598 * MDP(26) + t599 * MDP(28)) * t469 + (-t458 * MDP(10) + t482 * MDP(16) + (qJD(3) * t450 - t520) * MDP(21) + (qJD(3) * t518 + t394) * MDP(22) + t469 * MDP(23) + (-t336 + t581) * MDP(24) + (t337 - t580) * MDP(25) + (t330 + t581) * MDP(26) + (-t332 + t580) * MDP(28)) * t583 + (-t501 * t640 + t639) * t504; -t448 * t447 * MDP(12) + (-t447 ^ 2 + t448 ^ 2) * MDP(13) + (t447 * t482 + t507) * MDP(14) + (-t565 - t617) * MDP(15) + qJD(2) * t546 + (-t416 * t448 - t544 + (-qJD(4) - t482) * t380) * MDP(17) + (-t379 * t482 - t416 * t447 - t514) * MDP(18) + (-t341 * t469 - t638 + (-t394 * t448 + t469 * t572 + t549 * t627) * pkin(4) - t537) * MDP(24) + (-t342 * t469 + t645 + (-t448 * t520 + t469 * t550 - t496 * t549) * pkin(4) + t538) * MDP(25) + (-t621 - t351 * t394 + t530 * t469 + (pkin(5) - t488) * t549 - t537) * MDP(26) + (-t348 * t488 - t349 * t485 + (t332 + t530) * t520 + (t330 + t600) * t394) * MDP(27) + (t351 * t520 + t469 * t600 + t485 * t549 + t323 - t646) * MDP(28) + (t323 * t485 + t324 * t488 + t330 * t530 - t332 * t600 - t346 * t351) * MDP(29) + t647; (t515 - t638) * MDP(24) + (t516 + t645) * MDP(25) + (-t365 * t394 + t515 + 0.2e1 * t539 - t621) * MDP(26) + (pkin(5) * t348 - qJ(6) * t349 + (t332 - t337) * t520 + (t330 - t570) * t394) * MDP(27) + (t365 * t520 - 0.2e1 * t459 + 0.2e1 * t478 - t516 - t646) * MDP(28) + (-pkin(5) * t324 + qJ(6) * t323 - t330 * t337 + t332 * t570 - t346 * t365) * MDP(29) + t647; (-t549 + t652) * MDP(26) + t343 * MDP(27) + (-t469 ^ 2 - t629) * MDP(28) + (t332 * t469 + t324 + t621) * MDP(29);];
tauc  = t1;
