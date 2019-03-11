% Calculate vector of inverse dynamics joint torques for
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:14:14
% EndTime: 2019-03-08 23:14:23
% DurationCPUTime: 6.93s
% Computational Cost: add. (4332->510), mult. (9685->651), div. (0->0), fcn. (7602->14), ass. (0->231)
t506 = cos(qJ(3));
t650 = cos(qJ(4));
t592 = t650 * t506;
t562 = qJD(2) * t592;
t502 = sin(qJ(4));
t503 = sin(qJ(3));
t610 = qJD(2) * t503;
t586 = t502 * t610;
t444 = -t562 + t586;
t505 = cos(qJ(6));
t493 = qJD(3) + qJD(4);
t501 = sin(qJ(6));
t629 = t493 * t501;
t422 = -t505 * t444 + t629;
t453 = t502 * t506 + t503 * t650;
t446 = t453 * qJD(2);
t664 = qJD(6) + t446;
t667 = t422 * t664;
t424 = t444 * t501 + t493 * t505;
t568 = t664 * t424;
t666 = pkin(9) + pkin(8);
t462 = t666 * t503;
t463 = t666 * t506;
t426 = -t502 * t462 + t463 * t650;
t492 = qJDD(3) + qJDD(4);
t497 = qJ(3) + qJ(4);
t490 = sin(t497);
t500 = cos(pkin(6));
t507 = cos(qJ(2));
t641 = cos(pkin(11));
t579 = t641 * t507;
t498 = sin(pkin(11));
t504 = sin(qJ(2));
t627 = t498 * t504;
t438 = -t500 * t579 + t627;
t580 = t641 * t504;
t626 = t498 * t507;
t440 = t500 * t626 + t580;
t559 = g(1) * t440 + g(2) * t438;
t499 = sin(pkin(6));
t623 = t499 * t507;
t531 = g(3) * t623 - t559;
t529 = t531 * t490;
t593 = qJD(3) * t666;
t454 = t503 * t593;
t455 = t506 * t593;
t585 = qJD(4) * t650;
t612 = qJD(1) * t499;
t590 = t507 * t612;
t608 = qJD(4) * t502;
t621 = t502 * t503;
t660 = t592 - t621;
t616 = t650 * t454 + t502 * t455 + t462 * t585 + t463 * t608 + t660 * t590;
t665 = t426 * t492 - t493 * t616 - t529;
t625 = t499 * t504;
t511 = qJD(2) ^ 2;
t534 = (qJDD(2) * t507 - t504 * t511) * t499;
t591 = t504 * t612;
t578 = t666 * qJD(2) + t591;
t611 = qJD(1) * t500;
t419 = t503 * t611 + t506 * t578;
t414 = t650 * t419;
t418 = -t578 * t503 + t506 * t611;
t376 = t502 * t418 + t414;
t560 = pkin(3) * t608 - t376;
t413 = t502 * t419;
t377 = t418 * t650 - t413;
t663 = -pkin(3) * t585 - qJD(5) + t377;
t615 = qJD(4) * t426 - t453 * t590 - t502 * t454 + t455 * t650;
t624 = t499 * t506;
t443 = t500 * t503 + t504 * t624;
t540 = -t500 * t506 + t503 * t625;
t394 = t502 * t443 + t540 * t650;
t622 = t501 * t507;
t465 = t499 * t622;
t662 = t394 * t505 + t465;
t486 = t492 * qJ(5);
t661 = -t493 * qJD(5) - t486;
t659 = qJD(3) * t540;
t642 = qJD(3) * pkin(3);
t415 = t418 + t642;
t374 = -t650 * t415 + t413;
t604 = qJD(5) + t374;
t555 = t493 * t621;
t582 = qJDD(2) * t650;
t599 = qJDD(2) * t506;
t565 = -t493 * t562 - t502 * t599 - t503 * t582;
t387 = qJD(2) * t555 + t565;
t384 = -qJDD(6) + t387;
t484 = -pkin(3) * t650 - pkin(4);
t480 = -pkin(10) + t484;
t645 = t444 * pkin(5);
t658 = t664 * (t645 + t560) - t480 * t384;
t510 = qJD(3) ^ 2;
t602 = qJD(1) * qJD(2);
t584 = t504 * t602;
t554 = -qJDD(1) * t623 + t499 * t584;
t646 = g(3) * t507;
t656 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t510 + t499 * (t584 - t646) - t554 + t559;
t425 = t462 * t650 + t502 * t463;
t491 = cos(t497);
t653 = t425 * t492 + t491 * t531 + t493 * t615;
t652 = t446 ^ 2;
t651 = pkin(4) + pkin(10);
t649 = pkin(4) * t492;
t644 = t446 * pkin(5);
t643 = qJD(2) * pkin(2);
t417 = t493 * t453;
t600 = qJDD(2) * t503;
t553 = t502 * t600 - t506 * t582;
t388 = qJD(2) * t417 + t553;
t606 = qJD(6) * t505;
t594 = t501 * t388 + t444 * t606 + t505 * t492;
t607 = qJD(6) * t501;
t355 = -t493 * t607 + t594;
t639 = t355 * t505;
t375 = t502 * t415 + t414;
t638 = t375 * t493;
t637 = t384 * t501;
t485 = pkin(3) * t506 + pkin(2);
t546 = -qJ(5) * t453 - t485;
t390 = -t651 * t660 + t546;
t636 = t390 * t384;
t632 = t444 * t446;
t631 = t660 * t501;
t628 = t498 * t499;
t381 = t505 * t384;
t620 = qJDD(1) - g(3);
t619 = -pkin(5) * t417 - t616;
t416 = -qJD(3) * t592 - t506 * t585 + t555;
t618 = -t416 * pkin(5) + t615;
t614 = t644 - t663;
t495 = t503 ^ 2;
t613 = -t506 ^ 2 + t495;
t609 = qJD(2) * t504;
t605 = t644 + t604;
t601 = qJD(2) * qJD(3);
t598 = t500 * qJDD(1);
t597 = g(3) * t625;
t488 = t503 * t642;
t595 = t505 * t623;
t588 = t499 * t609;
t587 = qJD(2) * t623;
t583 = t503 * t601;
t581 = t499 * t641;
t429 = qJDD(2) * pkin(8) + (qJDD(1) * t504 + t507 * t602) * t499;
t577 = pkin(9) * qJDD(2) + t429;
t439 = t500 * t580 + t626;
t405 = t439 * t490 + t491 * t581;
t406 = t439 * t491 - t490 * t581;
t576 = -t405 * pkin(4) + t406 * qJ(5);
t441 = -t500 * t627 + t579;
t407 = t441 * t490 - t491 * t628;
t408 = t441 * t491 + t490 * t628;
t575 = -t407 * pkin(4) + qJ(5) * t408;
t432 = t490 * t625 - t500 * t491;
t433 = t490 * t500 + t491 * t625;
t574 = -t432 * pkin(4) + qJ(5) * t433;
t469 = t506 * t598;
t368 = qJDD(3) * pkin(3) - qJD(3) * t419 - t577 * t503 + t469;
t372 = qJD(3) * t418 + t503 * t598 + t577 * t506;
t566 = -t502 * t368 - t650 * t372 - t415 * t585 + t419 * t608;
t338 = t566 + t661;
t337 = -pkin(5) * t388 - t338;
t351 = -t493 * t651 + t605;
t435 = -qJD(2) * t485 - t590;
t525 = -qJ(5) * t446 + t435;
t373 = t444 * t651 + t525;
t343 = t351 * t501 + t373 * t505;
t573 = t337 * t505 - t343 * t444;
t572 = -t505 * t388 + t492 * t501;
t571 = t501 * t664;
t570 = t505 * t664;
t369 = -qJ(5) * t493 - t375;
t357 = -t369 - t645;
t569 = t664 * t357;
t567 = -t650 * t368 + t502 * t372 + t415 * t608 + t419 * t585;
t564 = t503 * t587;
t563 = t506 * t587;
t558 = g(1) * t441 + g(2) * t439;
t536 = qJ(5) * t416 - qJD(5) * t453 + t488;
t557 = -t417 * t651 - t536 + t591;
t358 = pkin(4) * t417 + t536;
t556 = -t358 + t591;
t398 = pkin(4) * t446 + qJ(5) * t444;
t342 = t351 * t505 - t373 * t501;
t549 = -(t375 - t645) * t664 + t651 * t384;
t548 = qJDD(5) + t567;
t545 = t337 * t501 + t342 * t444 + (t446 * t505 + t606) * t357;
t391 = pkin(3) * t610 + t398;
t543 = -g(1) * t498 + g(2) * t641;
t541 = -t394 * t501 + t595;
t395 = t443 * t650 - t502 * t540;
t537 = t417 * t501 - t606 * t660;
t535 = -g(1) * t408 - g(2) * t406 - g(3) * t433;
t533 = -t591 + t488;
t532 = t443 * qJD(3);
t457 = -t590 - t643;
t530 = -qJD(2) * t457 - t429 + t558;
t527 = -t535 + t566;
t526 = g(1) * t407 + g(2) * t405 + g(3) * t432 - t567;
t404 = pkin(3) * t583 - qJDD(2) * t485 + t554;
t524 = -pkin(8) * qJDD(3) + (t457 + t590 - t643) * qJD(3);
t392 = t453 * pkin(5) + t425;
t523 = -t337 * t660 + t357 * t417 + t392 * t384 + t558;
t436 = t446 * pkin(10);
t522 = (-qJD(6) * t480 + t391 + t436) * t664 + t535;
t521 = (qJD(6) * t651 + t398 + t436) * t664 + t535;
t520 = t532 + t564;
t519 = t563 - t659;
t518 = -t565 + (t444 - t586) * t493;
t517 = t435 * t444 + t527;
t516 = -t435 * t446 + t526;
t386 = pkin(4) * t444 + t525;
t515 = t386 * t446 + qJDD(5) - t526;
t356 = qJD(6) * t424 + t572;
t514 = t664 * t444 * MDP(27) + MDP(12) * t632 + ((-t356 - t568) * t505 + (-t355 + t667) * t501) * MDP(24) + (-t501 * t568 + t639) * MDP(23) + (t424 * t444 - t571 * t664 - t381) * MDP(25) + (-t422 * t444 - t570 * t664 + t637) * MDP(26) + t518 * MDP(14) - t553 * MDP(15) + (-t444 ^ 2 + t652) * MDP(13) + t492 * MDP(16);
t513 = -t386 * t444 - t527 - t661;
t512 = qJ(5) * t387 - qJD(5) * t446 + t404;
t481 = pkin(3) * t502 + qJ(5);
t403 = -pkin(4) * t660 + t546;
t393 = pkin(5) * t660 + t426;
t365 = -pkin(4) * t493 + t604;
t350 = qJD(4) * t395 + t502 * t519 + t520 * t650;
t349 = t443 * t608 + t502 * t520 - t519 * t650 + t540 * t585;
t346 = pkin(4) * t388 + t512;
t340 = t388 * t651 + t512;
t339 = t548 - t649;
t336 = -pkin(5) * t387 - t492 * t651 + t548;
t333 = t505 * t336;
t1 = [t620 * MDP(1) + MDP(3) * t534 + (-qJDD(2) * t504 - t507 * t511) * t499 * MDP(4) + (-t540 * qJDD(3) + t506 * t534 + (-t532 - 0.2e1 * t564) * qJD(3)) * MDP(10) + (-t443 * qJDD(3) - t503 * t534 + (-0.2e1 * t563 + t659) * qJD(3)) * MDP(11) + (t349 * t444 + t350 * t446 - t387 * t394 - t388 * t395) * MDP(19) + (-t338 * t395 + t339 * t394 + t349 * t369 + t350 * t365 - g(3) + (-t346 * t507 + t386 * t609) * t499) * MDP(22) + ((qJD(6) * t541 + t350 * t505 - t501 * t588) * t664 - t662 * t384 - t349 * t422 + t395 * t356) * MDP(28) + (-(qJD(6) * t662 + t350 * t501 + t505 * t588) * t664 - t541 * t384 - t349 * t424 + t395 * t355) * MDP(29) + (MDP(17) - MDP(20)) * (t499 * (-t388 * t507 + t444 * t609) - t350 * t493 - t394 * t492) + (-MDP(18) + MDP(21)) * (t499 * (-t387 * t507 - t446 * t609) - t349 * t493 + t395 * t492); (t333 * t453 - t342 * t416 + t393 * t356 + (-t340 * t453 + t490 * t559 + t636) * t501 - t523 * t505 - g(3) * (t490 * t622 + t504 * t505) * t499 + (t557 * t501 + t618 * t505) * t664 + t619 * t422 + ((-t390 * t505 - t392 * t501) * t664 - t343 * t453 - t357 * t631) * qJD(6)) * MDP(28) + (t355 * t453 + t384 * t631 - t416 * t424 + t537 * t664) * MDP(25) + (-t384 * t453 - t416 * t664) * MDP(27) + (t660 * t381 - t356 * t453 + t416 * t422 + (t417 * t505 + t607 * t660) * t664) * MDP(26) + (t343 * t416 + t393 * t355 + t619 * t424 + (t636 - (qJD(6) * t351 + t340) * t453 - t357 * qJD(6) * t660 + (-qJD(6) * t392 + t557) * t664 - t529) * t505 + (-(-qJD(6) * t373 + t336) * t453 + t597 + (qJD(6) * t390 - t618) * t664 + t523) * t501) * MDP(29) + (-t387 * t660 - t388 * t453 + t416 * t444 - t417 * t446) * MDP(13) + ((-t422 * t501 + t424 * t505) * t417 - (t639 - t356 * t501 + (-t422 * t505 - t424 * t501) * qJD(6)) * t660) * MDP(24) + (-t338 * t660 + t339 * t453 - t365 * t416 + t369 * t417 - t387 * t425 - t388 * t426 + t444 * t616 + t446 * t615 - t558 - t597) * MDP(19) + (-t388 * t485 - t404 * t660 + t417 * t435 + t444 * t533 - t653) * MDP(17) + (t346 * t660 - t386 * t417 - t388 * t403 + t444 * t556 + t653) * MDP(20) + (-t417 * t493 + t492 * t660) * MDP(15) + (-t338 * t426 + t339 * t425 + t346 * t403 + t386 * t358 - t558 * t666 + t616 * t369 + t615 * t365 + (-g(3) * t666 - qJD(1) * t386) * t625 + (-t499 * t646 + t559) * (pkin(4) * t491 + qJ(5) * t490 + t485)) * MDP(22) + (-t620 * t625 + t558) * MDP(4) + (t620 * t623 + t559) * MDP(3) + 0.2e1 * (t503 * t599 - t601 * t613) * MDP(6) + (-t355 * t631 + t424 * t537) * MDP(23) + (qJDD(2) * t495 + 0.2e1 * t506 * t583) * MDP(5) + (t387 * t485 + t404 * t453 - t416 * t435 + t446 * t533 - t665) * MDP(18) + (-t346 * t453 + t386 * t416 + t387 * t403 + t446 * t556 + t665) * MDP(21) + qJDD(2) * MDP(2) + (t524 * t503 + t506 * t656) * MDP(10) + (-t503 * t656 + t524 * t506) * MDP(11) + (qJDD(3) * t506 - t503 * t510) * MDP(8) + (qJDD(3) * t503 + t506 * t510) * MDP(7) + (-t416 * t493 + t453 * t492) * MDP(14) + (-t387 * t453 - t416 * t446) * MDP(12); (-t387 * t484 - t388 * t481 + (-t369 + t560) * t446 + (t365 + t663) * t444) * MDP(19) + (t391 * t444 + t560 * t493 + (-pkin(4) + t484) * t492 + t515) * MDP(20) + (t481 * t355 + t614 * t424 + t522 * t505 + (-t569 - t658) * t501 + t573) * MDP(29) + (t377 * t493 + (-t446 * t610 - t492 * t502 - t493 * t585) * pkin(3) + t517) * MDP(18) + (t376 * t493 + (-t444 * t610 + t492 * t650 - t493 * t608) * pkin(3) + t516) * MDP(17) + (-t338 * t481 + t339 * t484 - t386 * t391 - g(1) * ((-t441 * t503 + t498 * t624) * pkin(3) + t575) - g(2) * ((-t439 * t503 - t506 * t581) * pkin(3) + t576) - g(3) * (-pkin(3) * t540 + t574) + t663 * t369 + t560 * t365) * MDP(22) + (t481 * t356 + t614 * t422 + t522 * t501 + t505 * t658 + t545) * MDP(28) + (t391 * t446 + t481 * t492 - t493 * t663 + t513) * MDP(21) + MDP(8) * t599 + (g(3) * t443 + (-t499 * t543 - t598) * t503 + t530 * t506) * MDP(11) + (g(3) * t540 + t503 * t530 + t543 * t624 + t469) * MDP(10) + qJDD(3) * MDP(9) + t514 + MDP(7) * t600 + (-MDP(5) * t503 * t506 + MDP(6) * t613) * t511; (qJ(5) * t355 + t605 * t424 + (-t569 - t549) * t501 + t521 * t505 + t573) * MDP(29) + (-t339 * pkin(4) - g(1) * t575 - g(2) * t576 - g(3) * t574 - t338 * qJ(5) - t365 * t375 - t369 * t604 - t386 * t398) * MDP(22) + (t516 + t638) * MDP(17) + (pkin(4) * t387 - qJ(5) * t388 + (-t369 - t375) * t446 + (t365 - t604) * t444) * MDP(19) + (t398 * t444 + t515 - t638 - 0.2e1 * t649) * MDP(20) + (qJ(5) * t356 + t422 * t605 + t501 * t521 + t505 * t549 + t545) * MDP(28) + (t398 * t446 + t493 * t604 + t486 + t513) * MDP(21) + t514 + (-t374 * t493 + t517) * MDP(18); t518 * MDP(19) + (t492 - t632) * MDP(20) + (-t493 ^ 2 - t652) * MDP(21) + (t369 * t493 + t515 - t649) * MDP(22) + (-t422 * t493 - t381) * MDP(28) + (-t424 * t493 + t637) * MDP(29) + (-MDP(28) * t571 - MDP(29) * t570) * t664; t424 * t422 * MDP(23) + (-t422 ^ 2 + t424 ^ 2) * MDP(24) + (t594 + t667) * MDP(25) + (-t572 + t568) * MDP(26) - t384 * MDP(27) + (-t501 * t340 + t333 + t343 * t664 - t357 * t424 - g(1) * (t407 * t505 - t440 * t501) - g(2) * (t405 * t505 - t438 * t501) - g(3) * (t432 * t505 + t465)) * MDP(28) + (-t505 * t340 - t501 * t336 + t342 * t664 + t357 * t422 - g(1) * (-t407 * t501 - t440 * t505) - g(2) * (-t405 * t501 - t438 * t505) - g(3) * (-t432 * t501 + t595)) * MDP(29) + (-MDP(25) * t629 - MDP(26) * t424 - MDP(28) * t343 - MDP(29) * t342) * qJD(6);];
tau  = t1;
