% Calculate vector of inverse dynamics joint torques for
% S5RRRRP10
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRRP10_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:12:41
% EndTime: 2019-12-31 22:12:56
% DurationCPUTime: 10.28s
% Computational Cost: add. (5502->562), mult. (13720->756), div. (0->0), fcn. (10669->10), ass. (0->228)
t532 = cos(qJ(2));
t523 = sin(pkin(5));
t621 = qJD(1) * t523;
t511 = t532 * t621;
t561 = t511 - qJD(3);
t528 = sin(qJ(2));
t524 = cos(pkin(5));
t620 = qJD(1) * t524;
t603 = pkin(1) * t620;
t467 = pkin(7) * t511 + t528 * t603;
t527 = sin(qJ(3));
t531 = cos(qJ(3));
t684 = t467 + t561 * (pkin(3) * t527 - pkin(9) * t531);
t610 = qJDD(1) * t524;
t512 = qJDD(2) + t610;
t611 = qJD(1) * qJD(2);
t587 = t532 * t611;
t568 = t523 * t587;
t574 = qJD(2) * t603;
t609 = qJDD(1) * t528;
t586 = t523 * t609;
t601 = pkin(1) * t610;
t572 = t528 * t574 - t532 * t601 + (t568 + t586) * pkin(7);
t409 = -t512 * pkin(2) + t572;
t668 = cos(qJ(1));
t594 = t668 * t532;
t529 = sin(qJ(1));
t640 = t528 * t529;
t474 = -t524 * t594 + t640;
t595 = t668 * t528;
t639 = t529 * t532;
t476 = t524 * t639 + t595;
t645 = t523 * t532;
t544 = g(1) * t476 + g(2) * t474 - g(3) * t645;
t683 = t409 - t544;
t475 = t524 * t595 + t639;
t477 = -t524 * t640 + t594;
t564 = g(1) * t477 + g(2) * t475;
t648 = t523 * t528;
t543 = -g(3) * t648 - t564;
t577 = qJD(2) + t620;
t593 = t528 * t621;
t452 = t527 * t593 - t531 * t577;
t525 = -qJ(5) - pkin(9);
t682 = -qJ(5) * t452 + qJD(4) * t525;
t464 = -pkin(7) * t593 + t532 * t603;
t455 = t527 * t464;
t681 = pkin(3) * t593 - t455;
t596 = t523 * t668;
t429 = t475 * t531 - t527 * t596;
t526 = sin(qJ(4));
t530 = cos(qJ(4));
t680 = t429 * t526 - t474 * t530;
t679 = t429 * t530 + t474 * t526;
t559 = qJD(3) * t577;
t618 = qJD(2) * t532;
t591 = t527 * t618;
t616 = qJD(3) * t531;
t392 = -t531 * t512 + t523 * (qJD(1) * (t528 * t616 + t591) + t527 * t609) + t527 * t559;
t588 = t528 * t611;
t569 = t523 * t588;
t607 = qJDD(1) * t532;
t510 = t523 * t607;
t606 = qJDD(3) - t510;
t546 = t569 + t606;
t677 = pkin(8) * t546;
t513 = pkin(7) * t648;
t461 = t513 + (-pkin(1) * t532 - pkin(2)) * t524;
t472 = -t524 * t531 + t527 * t648;
t646 = t523 * t531;
t473 = t524 * t527 + t528 * t646;
t395 = pkin(3) * t472 - pkin(9) * t473 + t461;
t667 = pkin(1) * t528;
t624 = pkin(7) * t645 + t524 * t667;
t462 = pkin(8) * t524 + t624;
t560 = -pkin(2) * t532 - pkin(8) * t528 - pkin(1);
t463 = t560 * t523;
t630 = t531 * t462 + t527 * t463;
t397 = -pkin(9) * t645 + t630;
t632 = t526 * t395 + t530 * t397;
t617 = qJD(3) * t527;
t602 = pkin(8) * t617;
t675 = t526 * t602 - t684 * t530;
t590 = t531 * t618;
t608 = qJDD(1) * t531;
t642 = t527 * t512;
t535 = (t528 * t608 + (-t528 * t617 + t590) * qJD(1)) * t523 + t531 * t559 + t642;
t570 = t531 * t511;
t674 = t570 - t616;
t567 = pkin(2) * t528 - pkin(8) * t532;
t465 = t567 * t621;
t629 = t531 * t464 + t527 * t465;
t399 = pkin(9) * t593 + t629;
t491 = -pkin(3) * t531 - pkin(9) * t527 - pkin(2);
t613 = qJD(4) * t530;
t673 = t530 * t399 - t491 * t613 + t684 * t526;
t647 = t523 * t529;
t433 = t477 * t531 + t527 * t647;
t400 = -t433 * t526 + t476 * t530;
t637 = t530 * t532;
t426 = t473 * t526 + t523 * t637;
t671 = -g(1) * t400 + g(2) * t680 + g(3) * t426;
t454 = t527 * t577 + t531 * t593;
t418 = t530 * t454 - t526 * t561;
t367 = qJD(4) * t418 + t526 * t535 - t530 * t546;
t669 = t418 ^ 2;
t533 = qJD(1) ^ 2;
t666 = pkin(4) * t526;
t659 = pkin(8) * qJD(3);
t657 = qJ(5) * t527;
t581 = t530 * t561;
t614 = qJD(4) * t526;
t366 = qJD(4) * t581 + t454 * t614 - t526 * t546 - t530 * t535;
t656 = t366 * t526;
t387 = qJDD(4) + t392;
t655 = t387 * t530;
t416 = t454 * t526 + t581;
t446 = qJD(4) + t452;
t654 = t416 * t446;
t653 = t418 * t446;
t579 = t446 * t530;
t650 = t512 * MDP(8);
t520 = t523 ^ 2;
t649 = t520 * t533;
t644 = t524 * t532;
t643 = t526 * t387;
t641 = t527 * t530;
t638 = t530 * t531;
t438 = -pkin(2) * t577 - t464;
t380 = t452 * pkin(3) - t454 * pkin(9) + t438;
t439 = pkin(8) * t577 + t467;
t445 = qJD(1) * t463;
t391 = t531 * t439 + t527 * t445;
t383 = -pkin(9) * t561 + t391;
t360 = t530 * t380 - t383 * t526;
t350 = -qJ(5) * t418 + t360;
t349 = pkin(4) * t446 + t350;
t636 = -t350 + t349;
t444 = (t526 * t528 + t531 * t637) * t621;
t517 = pkin(8) * t638;
t571 = t527 * t511;
t612 = qJD(5) * t530;
t635 = -pkin(4) * t571 + qJ(5) * t444 + t399 * t526 - t527 * t612 + (pkin(4) * t527 - qJ(5) * t638) * qJD(3) + (-t517 + (-t491 + t657) * t526) * qJD(4) + t675;
t443 = t526 * t570 - t530 * t593;
t634 = qJ(5) * t443 + (-qJ(5) * qJD(4) - t659) * t641 + (-qJD(5) * t527 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t531) * t526 - t673;
t390 = -t527 * t439 + t531 * t445;
t410 = pkin(3) * t454 + pkin(9) * t452;
t633 = t530 * t390 + t526 * t410;
t628 = t526 * t682 + t612 - t633;
t405 = t530 * t410;
t627 = -pkin(4) * t454 - t405 + t682 * t530 + (-qJD(5) + t390) * t526;
t623 = t526 * t491 + t517;
t521 = t528 ^ 2;
t622 = -t532 ^ 2 + t521;
t619 = qJD(2) * t528;
t615 = qJD(4) * t446;
t605 = 0.2e1 * t520;
t600 = t532 * t649;
t599 = t526 * t645;
t598 = -pkin(7) * t510 - t528 * t601 - t532 * t574;
t597 = pkin(8) + t666;
t592 = t523 * t619;
t589 = t523 * t524 * t533;
t583 = t530 * t395 - t397 * t526;
t582 = -t527 * t462 + t463 * t531;
t580 = t532 * t561;
t578 = qJD(3) * t561;
t576 = qJD(2) + 0.2e1 * t620;
t575 = t512 + t610;
t541 = -pkin(7) * t569 - t598;
t408 = pkin(8) * t512 + t541;
t556 = t567 * qJD(2);
t411 = (qJD(1) * t556 + qJDD(1) * t560) * t523;
t573 = t527 * t408 - t531 * t411 + t439 * t616 + t445 * t617;
t428 = t475 * t527 + t531 * t596;
t432 = t477 * t527 - t529 * t646;
t565 = -g(1) * t428 + g(2) * t432;
t396 = pkin(3) * t645 - t582;
t563 = -t526 * t616 + t443;
t562 = t530 * t616 - t444;
t361 = t380 * t526 + t383 * t530;
t466 = t523 * t556;
t468 = (t644 * pkin(1) - t513) * qJD(2);
t558 = -t462 * t616 - t463 * t617 + t466 * t531 - t527 * t468;
t382 = pkin(3) * t561 - t390;
t554 = -pkin(9) * t387 + t382 * t446;
t553 = t587 + t609;
t552 = -t531 * t408 - t527 * t411 + t439 * t617 - t445 * t616;
t551 = -t462 * t617 + t463 * t616 + t527 * t466 + t531 * t468;
t353 = pkin(9) * t546 - t552;
t357 = t392 * pkin(3) - pkin(9) * t535 + t409;
t550 = -t530 * t353 - t526 * t357 - t380 * t613 + t383 * t614;
t370 = pkin(9) * t592 + t551;
t424 = qJD(3) * t473 + t523 * t591;
t425 = -qJD(3) * t472 + t523 * t590;
t469 = t624 * qJD(2);
t374 = pkin(3) * t424 - pkin(9) * t425 + t469;
t549 = t530 * t370 + t526 * t374 + t395 * t613 - t397 * t614;
t548 = g(1) * t432 + g(2) * t428 + g(3) * t472;
t547 = -g(1) * t433 - g(2) * t429 - g(3) * t473;
t371 = -pkin(3) * t592 - t558;
t539 = -qJD(4) * t632 - t370 * t526 + t530 * t374;
t345 = -qJD(4) * t361 - t526 * t353 + t530 * t357;
t538 = pkin(8) * t615 - t544;
t354 = -pkin(3) * t546 + t573;
t537 = pkin(9) * t615 + t354 - t548;
t348 = t367 * pkin(4) + qJDD(5) + t354;
t519 = pkin(4) * t530 + pkin(3);
t496 = t525 * t530;
t495 = t525 * t526;
t483 = t530 * t491;
t436 = -t526 * t657 + t623;
t427 = t473 * t530 - t599;
t420 = -qJ(5) * t641 + t483 + (-pkin(8) * t526 - pkin(4)) * t531;
t415 = t416 ^ 2;
t401 = t433 * t530 + t476 * t526;
t398 = -t465 * t531 - t681;
t379 = -qJD(4) * t426 + t425 * t530 + t526 * t592;
t378 = -qJD(4) * t599 + t425 * t526 + t473 * t613 - t530 * t592;
t368 = t416 * pkin(4) + qJD(5) + t382;
t362 = -qJ(5) * t426 + t632;
t358 = pkin(4) * t472 - qJ(5) * t427 + t583;
t351 = -qJ(5) * t416 + t361;
t347 = -qJ(5) * t378 - qJD(5) * t426 + t549;
t346 = pkin(4) * t424 - qJ(5) * t379 - qJD(5) * t427 + t539;
t343 = -qJ(5) * t367 - qJD(5) * t416 - t550;
t342 = pkin(4) * t387 + qJ(5) * t366 - qJD(5) * t418 + t345;
t1 = [(t387 * t472 + t424 * t446) * MDP(22) + (-t367 * t472 - t378 * t446 - t387 * t426 - t416 * t424) * MDP(21) + (-t366 * t472 + t379 * t446 + t387 * t427 + t418 * t424) * MDP(20) + (t366 * t426 - t367 * t427 - t378 * t418 - t379 * t416) * MDP(19) + (-t366 * t427 + t379 * t418) * MDP(18) + t524 * t650 + (-t473 * t392 - t454 * t424 - t425 * t452 - t472 * t535) * MDP(12) + (t454 * t425 + t473 * t535) * MDP(11) + (g(1) * t529 - g(2) * t668) * MDP(2) + (g(1) * t668 + g(2) * t529) * MDP(3) + (t343 * t362 + t351 * t347 + t342 * t358 + t349 * t346 + t348 * (pkin(4) * t426 + t396) + t368 * (pkin(4) * t378 + t371) - g(1) * (-t529 * pkin(1) - t475 * pkin(2) + pkin(7) * t596 + t428 * t525 - t429 * t519 - t474 * t597) - g(2) * (pkin(1) * t668 + t477 * pkin(2) + pkin(7) * t647 - t432 * t525 + t433 * t519 + t476 * t597)) * MDP(26) + (-g(1) * t680 - g(2) * t400 + t354 * t427 - t361 * t424 - t396 * t366 + t371 * t418 + t382 * t379 - t632 * t387 - t549 * t446 + t550 * t472) * MDP(24) + (g(1) * t679 - g(2) * t401 + t345 * t472 + t354 * t426 + t360 * t424 + t396 * t367 + t371 * t416 + t382 * t378 + t583 * t387 + t539 * t446) * MDP(23) + (-t469 * t577 - t513 * t512 - t572 * t524 + g(1) * t475 - g(2) * t477 + (t512 * t644 + (-t588 + t607) * t605) * pkin(1)) * MDP(9) + (-t425 * t561 + t473 * t606) * MDP(13) + (-t391 * t592 + t409 * t473 + t438 * t425 + t469 * t454 + t461 * t535 - t546 * t630 + t551 * t561 - t552 * t645 + t565) * MDP(17) + (g(1) * t429 - g(2) * t433 + t390 * t592 + t461 * t392 + t409 * t472 + t438 * t424 + t469 * t452 + t546 * t582 - t558 * t561 + t573 * t645) * MDP(16) + (-pkin(1) * t553 * t605 - g(1) * t474 + g(2) * t476 - t468 * t577 - t512 * t624 - t524 * t541) * MDP(10) + (t424 * t561 - t472 * t606) * MDP(14) + qJDD(1) * MDP(1) + (-t342 * t427 - t343 * t426 - t346 * t418 - t347 * t416 - t349 * t379 - t351 * t378 + t358 * t366 - t362 * t367 - t565) * MDP(25) + ((qJDD(1) * t521 + 0.2e1 * t528 * t587) * MDP(4) + 0.2e1 * (t528 * t607 - t611 * t622) * MDP(5)) * t520 + (((-t642 + (-t559 - t568) * t531) * t532 + (-(-qJD(1) * t617 + t608) * t645 + (qJD(1) * t473 + t454) * qJD(2)) * t528) * MDP(13) + (t392 * t532 + (-qJD(1) * t472 - t452) * t619) * MDP(14) + (t528 * t575 + t576 * t618) * MDP(6) + (-t606 * t532 + (-t511 - t561) * t619) * MDP(15) + (t532 * t575 - t576 * t619) * MDP(7)) * t523; -t528 * MDP(4) * t600 + t622 * MDP(5) * t649 + (-t532 * t589 + t586) * MDP(6) + (t528 * t589 + t510) * MDP(7) + t650 + (t467 * t577 + t649 * t667 + t544 - t572) * MDP(9) + (pkin(1) * t600 + t464 * t577 + (pkin(7) * t611 + g(3)) * t648 + t564 + t598) * MDP(10) + ((-qJD(3) * t593 + t512) * t527 ^ 2 + ((t523 * t553 + t559) * t527 - t561 * t454) * t531) * MDP(11) + (-t527 * t392 + t531 * t535 + (t571 - t617) * t454 + t674 * t452) * MDP(12) + (-t531 * t578 + t527 * t606 + (t531 * t580 + (qJD(2) * t527 - t454) * t528) * t621) * MDP(13) + (t527 * t578 + t531 * t606 + (-t527 * t580 + (qJD(2) * t531 + t452) * t528) * t621) * MDP(14) + t561 * MDP(15) * t593 + (-pkin(2) * t392 - t455 * t561 - t390 * t593 - t467 * t452 + (-t438 * t561 - t677) * t527 + (pkin(8) * t578 + t465 * t561 - t683) * t531) * MDP(16) + (-t531 * t677 - pkin(2) * t535 + t391 * t593 - t467 * t454 + (-t602 - t629) * t561 - t674 * t438 + t683 * t527) * MDP(17) + (-t366 * t641 + (-t527 * t614 + t562) * t418) * MDP(18) + (t416 * t444 + t418 * t443 + (-t416 * t530 - t418 * t526) * t616 + (t656 - t367 * t530 + (t416 * t526 - t418 * t530) * qJD(4)) * t527) * MDP(19) + (t366 * t531 + t562 * t446 + (-t418 * t561 - t446 * t614 + t655) * t527) * MDP(20) + (t367 * t531 + t563 * t446 + (t416 * t561 - t446 * t613 - t643) * t527) * MDP(21) + (-t446 * t527 * t561 - t387 * t531) * MDP(22) + (-t382 * t443 + t483 * t387 - t398 * t416 + t675 * t446 + ((-qJD(4) * t491 + t399) * t446 + t543) * t526 + (t382 * t526 * qJD(3) - t345 + (qJD(3) * t416 - t643) * pkin(8) - t538 * t530) * t531 + (pkin(8) * t367 + t354 * t526 - t360 * t561 + t382 * t613) * t527) * MDP(23) + (-t623 * t387 - t398 * t418 - t382 * t444 + t673 * t446 + t543 * t530 + (-t550 + (pkin(8) * t418 + t382 * t530) * qJD(3) + t538 * t526) * t531 + (-t382 * t614 + t354 * t530 + t561 * t361 + (qJD(3) * t579 - t366) * pkin(8)) * t527) * MDP(24) + (t349 * t444 + t351 * t443 + t366 * t420 - t367 * t436 - t635 * t418 - t634 * t416 + (-t349 * t530 - t351 * t526) * t616 + (-t342 * t530 - t343 * t526 + (t349 * t526 - t351 * t530) * qJD(4) + t544) * t527) * MDP(25) + (t342 * t420 + t343 * t436 + t635 * t349 + t634 * t351 + ((t465 + t659) * t531 + (t527 * t613 - t563) * pkin(4) + t681) * t368 + (t348 * t527 + t543) * t597 + t544 * (t519 * t531 - t525 * t527 + pkin(2))) * MDP(26); -t452 ^ 2 * MDP(12) + (-t452 * t561 + t535) * MDP(13) - t392 * MDP(14) + t546 * MDP(15) + (-t391 * t561 + t548 - t573) * MDP(16) + (-t390 * t561 + t438 * t452 - t547 + t552) * MDP(17) + (t418 * t579 - t656) * MDP(18) + ((-t366 - t654) * t530 + (-t367 - t653) * t526) * MDP(19) + (t446 * t579 + t643) * MDP(20) + (-t446 ^ 2 * t526 + t655) * MDP(21) + (-pkin(3) * t367 - t391 * t416 - t405 * t446 + (t390 * t446 + t554) * t526 - t537 * t530) * MDP(23) + (pkin(3) * t366 - t391 * t418 + t446 * t633 + t526 * t537 + t530 * t554) * MDP(24) + (t366 * t495 + t367 * t496 - t627 * t418 - t628 * t416 + (-t349 * t446 + t343) * t530 + (-t351 * t446 - t342) * t526 + t547) * MDP(25) + (-t343 * t496 + t342 * t495 - t348 * t519 - g(1) * (-t432 * t519 - t433 * t525) - g(2) * (-t428 * t519 - t429 * t525) - g(3) * (-t472 * t519 - t473 * t525) + (t446 * t666 - t391) * t368 + t628 * t351 + t627 * t349) * MDP(26) + (MDP(11) * t452 + MDP(12) * t454 - t561 * MDP(14) - t438 * MDP(16) - t418 * MDP(20) + t416 * MDP(21) - t446 * MDP(22) - t360 * MDP(23) + t361 * MDP(24)) * t454; t418 * t416 * MDP(18) + (-t415 + t669) * MDP(19) + (-t366 + t654) * MDP(20) + (-t367 + t653) * MDP(21) + t387 * MDP(22) + (t361 * t446 - t382 * t418 + t345 + t671) * MDP(23) + (g(1) * t401 + g(2) * t679 + g(3) * t427 + t360 * t446 + t382 * t416 + t550) * MDP(24) + (pkin(4) * t366 - t416 * t636) * MDP(25) + (t636 * t351 + (-t368 * t418 + t342 + t671) * pkin(4)) * MDP(26); (-t415 - t669) * MDP(25) + (t349 * t418 + t351 * t416 + t348 - t548) * MDP(26);];
tau = t1;
