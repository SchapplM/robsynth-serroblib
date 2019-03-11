% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRPR2
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
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:02:43
% EndTime: 2019-03-09 05:02:56
% DurationCPUTime: 6.67s
% Computational Cost: add. (4734->432), mult. (11300->597), div. (0->0), fcn. (7898->10), ass. (0->206)
t519 = cos(qJ(6));
t517 = sin(qJ(4));
t520 = cos(qJ(4));
t572 = t520 * qJD(3);
t518 = sin(qJ(3));
t585 = qJD(1) * t518;
t482 = t517 * t585 - t572;
t582 = qJD(3) * t517;
t484 = t520 * t585 + t582;
t512 = sin(pkin(11));
t514 = cos(pkin(11));
t537 = -t482 * t514 - t484 * t512;
t600 = t519 * t537;
t425 = t482 * t512 - t484 * t514;
t516 = sin(qJ(6));
t618 = t425 * t516;
t373 = t600 + t618;
t521 = cos(qJ(3));
t584 = qJD(1) * t521;
t501 = -qJD(4) + t584;
t496 = -qJD(6) + t501;
t619 = t373 * t496;
t571 = qJD(1) * qJD(3);
t557 = t521 * t571;
t578 = qJD(4) * t518;
t560 = t517 * t578;
t570 = qJD(3) * qJD(4);
t434 = -qJD(1) * t560 + (t557 + t570) * t520;
t503 = sin(pkin(10)) * pkin(1) + pkin(7);
t490 = t503 * qJD(1);
t508 = t518 * qJD(2);
t456 = t521 * t490 + t508;
t441 = qJD(3) * pkin(8) + t456;
t505 = -cos(pkin(10)) * pkin(1) - pkin(2);
t472 = -pkin(3) * t521 - pkin(8) * t518 + t505;
t446 = t472 * qJD(1);
t606 = t517 * t446;
t401 = t441 * t520 + t606;
t626 = qJD(2) * t521 - t518 * t490;
t444 = t626 * qJD(3);
t543 = pkin(3) * t518 - pkin(8) * t521;
t487 = t543 * qJD(3);
t471 = qJD(1) * t487;
t547 = t517 * t444 - t520 * t471;
t524 = -qJD(4) * t401 - t547;
t558 = t518 * t571;
t344 = pkin(4) * t558 - qJ(5) * t434 - qJD(5) * t484 + t524;
t577 = qJD(4) * t520;
t559 = t518 * t577;
t580 = qJD(3) * t521;
t562 = t517 * t580;
t632 = t559 + t562;
t435 = qJD(1) * t632 + t517 * t570;
t567 = t520 * t444 + t446 * t577 + t517 * t471;
t579 = qJD(4) * t517;
t527 = -t441 * t579 + t567;
t350 = -qJ(5) * t435 - qJD(5) * t482 + t527;
t326 = t514 * t344 - t350 * t512;
t390 = t434 * t514 - t435 * t512;
t324 = pkin(5) * t558 - pkin(9) * t390 + t326;
t327 = t512 * t344 + t514 * t350;
t389 = -t434 * t512 - t435 * t514;
t325 = pkin(9) * t389 + t327;
t400 = -t441 * t517 + t520 * t446;
t382 = -qJ(5) * t484 + t400;
t375 = -pkin(4) * t501 + t382;
t383 = -qJ(5) * t482 + t401;
t607 = t514 * t383;
t343 = t512 * t375 + t607;
t630 = pkin(9) * t537;
t336 = t343 + t630;
t575 = qJD(6) * t516;
t335 = t336 * t575;
t440 = -qJD(3) * pkin(3) - t626;
t419 = pkin(4) * t482 + qJD(5) + t440;
t376 = -pkin(5) * t537 + t419;
t641 = -t516 * t324 - t519 * t325 - t376 * t373 + t335;
t581 = qJD(3) * t518;
t555 = MDP(25) * t581;
t631 = -t519 * t425 + t516 * t537;
t640 = qJD(1) * t555 + (-t373 ^ 2 + t631 ^ 2) * MDP(22) - t373 * t631 * MDP(21);
t620 = t631 * t496;
t599 = t520 * t521;
t532 = pkin(4) * t518 - qJ(5) * t599;
t486 = t543 * qJD(1);
t546 = t520 * t486 - t517 * t626;
t623 = -qJ(5) - pkin(8);
t553 = qJD(4) * t623;
t638 = -qJD(1) * t532 - qJD(5) * t517 + t520 * t553 - t546;
t563 = t517 * t584;
t576 = qJD(5) * t520;
t591 = t517 * t486 + t520 * t626;
t637 = -qJ(5) * t563 - t517 * t553 - t576 + t591;
t552 = t519 * t324 - t516 * t325;
t636 = -t376 * t631 + t552;
t635 = pkin(9) * t425;
t476 = t512 * t520 + t514 * t517;
t528 = t476 * t521;
t634 = qJD(1) * t528 - t476 * qJD(4);
t536 = t512 * t517 - t514 * t520;
t633 = t501 * t536;
t629 = MDP(5) * t518;
t510 = t518 ^ 2;
t628 = MDP(6) * (-t521 ^ 2 + t510);
t594 = t512 * t637 + t514 * t638;
t593 = t512 * t638 - t514 * t637;
t625 = -t456 + (-t563 + t579) * pkin(4);
t550 = -t519 * t389 + t390 * t516;
t331 = qJD(6) * t631 + t550;
t624 = pkin(4) * t512;
t622 = t331 * t521;
t409 = -qJD(3) * t528 + t536 * t578;
t561 = t521 * t572;
t410 = t476 * t578 + t512 * t562 - t514 * t561;
t453 = t476 * t518;
t454 = t536 * t518;
t539 = -t519 * t453 + t454 * t516;
t346 = qJD(6) * t539 + t409 * t516 - t410 * t519;
t621 = t346 * t496;
t617 = t434 * t517;
t616 = t435 * t521;
t615 = t440 * t517;
t614 = t440 * t520;
t445 = qJD(3) * t508 + t490 * t580;
t613 = t445 * t517;
t612 = t445 * t520;
t611 = t482 * t501;
t610 = t484 * t501;
t609 = t501 * t520;
t608 = t503 * t517;
t378 = t512 * t383;
t605 = t517 * t518;
t604 = t517 * t521;
t603 = t518 * t520;
t522 = qJD(3) ^ 2;
t602 = t518 * t522;
t342 = t514 * t375 - t378;
t332 = -pkin(5) * t501 + t342 + t635;
t601 = t519 * t332;
t598 = t521 * t522;
t408 = -t453 * t516 - t454 * t519;
t347 = qJD(6) * t408 - t519 * t409 - t410 * t516;
t597 = t347 * t496 + t539 * t558;
t485 = t503 * t599;
t589 = t520 * t487 + t581 * t608;
t360 = -t518 * t576 + t532 * qJD(3) + (-t485 + (qJ(5) * t518 - t472) * t517) * qJD(4) + t589;
t590 = t472 * t577 + t517 * t487;
t364 = (-qJ(5) * qJD(4) - qJD(3) * t503) * t603 + (-qJD(5) * t518 + (-qJ(5) * qJD(3) - qJD(4) * t503) * t521) * t517 + t590;
t334 = t512 * t360 + t514 * t364;
t538 = -t476 * t516 - t519 * t536;
t596 = qJD(6) * t538 + t516 * t634 + t519 * t633;
t424 = t476 * t519 - t516 * t536;
t595 = qJD(6) * t424 + t516 * t633 - t519 * t634;
t349 = t514 * t382 - t378;
t460 = t520 * t472;
t411 = -qJ(5) * t603 + t460 + (-pkin(4) - t608) * t521;
t588 = t517 * t472 + t485;
t418 = -qJ(5) * t605 + t588;
t363 = t512 * t411 + t514 * t418;
t592 = -pkin(5) * t634 + t625;
t492 = t623 * t517;
t493 = t623 * t520;
t430 = t512 * t492 - t514 * t493;
t587 = pkin(4) * t605 + t518 * t503;
t491 = qJD(1) * t505;
t573 = t440 * qJD(4);
t568 = qJD(6) * t600 + t516 * t389 + t519 * t390;
t565 = pkin(4) * t632 + t503 * t580;
t564 = -pkin(4) * t520 - pkin(3);
t554 = MDP(16) * t581;
t330 = t425 * t575 + t568;
t551 = -t330 * t521 + t581 * t631;
t333 = t514 * t360 - t364 * t512;
t348 = -t382 * t512 - t607;
t362 = t514 * t411 - t418 * t512;
t549 = -t434 * t521 + t484 * t581;
t548 = t501 * t503 + t441;
t429 = t514 * t492 + t493 * t512;
t545 = t501 * t560;
t544 = t501 * t559;
t402 = pkin(4) * t435 + t445;
t415 = -pkin(9) * t536 + t430;
t542 = pkin(5) * t585 + pkin(9) * t633 + qJD(6) * t415 - t594;
t414 = -pkin(9) * t476 + t429;
t541 = pkin(9) * t634 + qJD(6) * t414 + t593;
t322 = t516 * t332 + t519 * t336;
t353 = -pkin(5) * t521 + pkin(9) * t454 + t362;
t354 = -pkin(9) * t453 + t363;
t540 = t353 * t516 + t354 * t519;
t535 = qJD(1) * t510 - t501 * t521;
t533 = 0.2e1 * qJD(3) * t491;
t504 = pkin(4) * t514 + pkin(5);
t531 = t504 * t516 + t519 * t624;
t530 = t504 * t519 - t516 * t624;
t526 = t535 * t517;
t447 = pkin(5) * t536 + t564;
t420 = pkin(5) * t453 + t587;
t395 = pkin(4) * t484 - pkin(5) * t425;
t377 = -pkin(5) * t409 + t565;
t357 = -pkin(5) * t389 + t402;
t338 = t349 + t635;
t337 = t348 - t630;
t329 = pkin(9) * t409 + t334;
t328 = pkin(5) * t581 + pkin(9) * t410 + t333;
t321 = -t336 * t516 + t601;
t1 = [0.2e1 * t557 * t629 - 0.2e1 * t571 * t628 + MDP(7) * t598 - MDP(8) * t602 + (-t503 * t598 + t518 * t533) * MDP(10) + (t503 * t602 + t521 * t533) * MDP(11) + (t434 * t603 + (-t560 + t561) * t484) * MDP(12) + ((-t482 * t520 - t484 * t517) * t580 + (-t617 - t435 * t520 + (t482 * t517 - t484 * t520) * qJD(4)) * t518) * MDP(13) + (t535 * t572 + t545 + t549) * MDP(14) + (t544 + t616 + (-t482 * t518 - t526) * qJD(3)) * MDP(15) + (-t501 - t584) * t554 + (-(-t472 * t579 + t589) * t501 + ((t482 * t503 + t615) * qJD(3) + (t520 * t548 + t606) * qJD(4) + t547) * t521 + (t520 * t573 + t503 * t435 + t613 + ((-t503 * t604 + t460) * qJD(1) + t400) * qJD(3)) * t518) * MDP(17) + (t590 * t501 + (-t548 * t579 + (t484 * t503 + t614) * qJD(3) + t567) * t521 + (-t517 * t573 + t503 * t434 + t612 + (-qJD(1) * t588 - t503 * t609 - t401) * qJD(3)) * t518) * MDP(18) + (t326 * t454 - t327 * t453 + t333 * t425 + t334 * t537 + t342 * t410 + t343 * t409 - t362 * t390 + t363 * t389) * MDP(19) + (t326 * t362 + t327 * t363 + t342 * t333 + t343 * t334 + t402 * t587 + t419 * t565) * MDP(20) + (t330 * t408 + t346 * t631) * MDP(21) + (t330 * t539 - t331 * t408 + t346 * t373 - t347 * t631) * MDP(22) + (t408 * t558 + t551 - t621) * MDP(23) + (t373 * t581 + t597 + t622) * MDP(24) + (-t496 - t584) * t555 + (-(t328 * t519 - t329 * t516) * t496 - t552 * t521 - t377 * t373 + t420 * t331 - t357 * t539 + t376 * t347 + (t322 * t521 + t496 * t540) * qJD(6) + ((t353 * t519 - t354 * t516) * qJD(1) + t321) * t581) * MDP(26) + (t420 * t330 - t335 * t521 + t376 * t346 + t357 * t408 + t377 * t631 + ((-qJD(6) * t354 + t328) * t496 + t324 * t521) * t516 + ((qJD(6) * t353 + t329) * t496 + (qJD(6) * t332 + t325) * t521) * t519 + (-qJD(1) * t540 - t322) * t581) * MDP(27); (t544 - t616) * MDP(17) + (-t545 + t549) * MDP(18) + (-t389 * t454 + t390 * t453 + t409 * t425 - t410 * t537) * MDP(19) + (-t326 * t453 - t327 * t454 + t342 * t409 - t343 * t410 - t402 * t521) * MDP(20) + (t597 - t622) * MDP(26) + (t551 + t621) * MDP(27) + (-MDP(10) * t518 - MDP(11) * t521) * t522 + (-t535 * MDP(18) * t520 - MDP(17) * t526 + (-MDP(27) * qJD(1) * t408 + t482 * MDP(17) + MDP(20) * t419 - MDP(26) * t373) * t518) * qJD(3); (qJD(3) * t456 - t445) * MDP(10) - t491 * t584 * MDP(11) + (-t484 * t609 + t617) * MDP(12) + ((t434 + t611) * t520 + (-t435 + t610) * t517) * MDP(13) + (-t501 * t577 + (t501 * t599 + (-t484 + t582) * t518) * qJD(1)) * MDP(14) + (t501 * t579 + (-t501 * t604 + (t482 + t572) * t518) * qJD(1)) * MDP(15) + (-pkin(3) * t435 - t612 + t546 * t501 - t456 * t482 + (pkin(8) * t609 + t615) * qJD(4) + (-t400 * t518 + (-pkin(8) * t581 - t440 * t521) * t517) * qJD(1)) * MDP(17) + (-pkin(3) * t434 + t613 - t591 * t501 - t456 * t484 + (-pkin(8) * t501 * t517 + t614) * qJD(4) + (-t440 * t599 + (-pkin(8) * t572 + t401) * t518) * qJD(1)) * MDP(18) + (-t326 * t476 - t327 * t536 - t342 * t633 + t343 * t634 + t389 * t430 - t390 * t429 + t594 * t425 + t593 * t537) * MDP(19) + (t326 * t429 + t327 * t430 + t594 * t342 + t593 * t343 + t402 * t564 + t419 * t625) * MDP(20) + (t330 * t424 + t596 * t631) * MDP(21) + (t330 * t538 - t331 * t424 + t373 * t596 - t595 * t631) * MDP(22) + (t447 * t331 - t357 * t538 - t373 * t592 + t595 * t376) * MDP(26) + (t447 * t330 + t357 * t424 + t596 * t376 + t592 * t631) * MDP(27) + (-t596 * MDP(23) + t595 * MDP(24) + (t516 * t541 + t519 * t542) * MDP(26) + (-t516 * t542 + t519 * t541) * MDP(27)) * t496 + (-t491 * MDP(10) + t501 * MDP(16) + (qJD(3) * t424 - t631) * MDP(23) + (qJD(3) * t538 - t373) * MDP(24) + t496 * MDP(25) + ((t414 * t519 - t415 * t516) * qJD(3) - t321) * MDP(26) + (-(t414 * t516 + t415 * t519) * qJD(3) + t322) * MDP(27)) * t585 + (-t521 * t629 + t628) * qJD(1) ^ 2; t484 * t482 * MDP(12) + (-t482 ^ 2 + t484 ^ 2) * MDP(13) + (t434 - t611) * MDP(14) + (-t435 - t610) * MDP(15) + qJD(1) * t554 + (-t401 * t501 - t440 * t484 + t524) * MDP(17) + (-t400 * t501 + t440 * t482 - t527) * MDP(18) + ((t389 * t512 - t390 * t514) * pkin(4) + (t342 - t349) * t537 + (-t343 - t348) * t425) * MDP(19) + (-t342 * t348 - t343 * t349 + (t326 * t514 + t327 * t512 - t419 * t484) * pkin(4)) * MDP(20) + (t330 + t619) * MDP(23) + (-t331 - t620) * MDP(24) + (t530 * t558 + (t337 * t519 - t338 * t516) * t496 + t395 * t373 + (t496 * t531 - t322) * qJD(6) + t636) * MDP(26) + (-t531 * t558 - (t337 * t516 + t338 * t519) * t496 - t395 * t631 + (t496 * t530 - t601) * qJD(6) + t641) * MDP(27) + t640; (-t425 ^ 2 - t537 ^ 2) * MDP(19) + (-t342 * t425 - t343 * t537 + t402) * MDP(20) + (t331 - t620) * MDP(26) + (t330 - t619) * MDP(27); (t568 + t619) * MDP(23) + (-t550 - t620) * MDP(24) + (-t322 * t496 + t636) * MDP(26) + (-t321 * t496 + t641) * MDP(27) + (MDP(23) * t618 - MDP(24) * t631 - MDP(26) * t322 - MDP(27) * t601) * qJD(6) + t640;];
tauc  = t1;
