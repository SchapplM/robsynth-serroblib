% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S5RRRRR11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:44:01
% EndTime: 2019-12-31 22:44:20
% DurationCPUTime: 10.10s
% Computational Cost: add. (6527->550), mult. (17092->762), div. (0->0), fcn. (13472->10), ass. (0->222)
t540 = sin(pkin(5));
t549 = cos(qJ(2));
t618 = qJD(1) * t549;
t601 = t540 * t618;
t524 = -qJD(3) + t601;
t545 = sin(qJ(2));
t541 = cos(pkin(5));
t619 = qJD(1) * t541;
t606 = pkin(1) * t619;
t493 = pkin(7) * t601 + t545 * t606;
t544 = sin(qJ(3));
t548 = cos(qJ(3));
t678 = t493 + t524 * (pkin(3) * t544 - pkin(9) * t548);
t543 = sin(qJ(4));
t547 = cos(qJ(4));
t620 = qJD(1) * t540;
t602 = t545 * t620;
t463 = t543 * t548 * t601 - t547 * t602;
t613 = qJD(3) * t548;
t671 = -t543 * t613 + t463;
t636 = t548 * t549;
t464 = (t543 * t545 + t547 * t636) * t620;
t571 = t547 * t613 - t464;
t612 = qJD(4) * t543;
t677 = -t544 * t612 + t571;
t582 = qJD(2) + t619;
t458 = pkin(8) * t582 + t493;
t489 = (-pkin(2) * t549 - pkin(8) * t545 - pkin(1)) * t540;
t469 = qJD(1) * t489;
t411 = -t544 * t458 + t469 * t548;
t400 = pkin(3) * t524 - t411;
t477 = t544 * t582 + t548 * t602;
t435 = t477 * t543 + t547 * t524;
t370 = pkin(4) * t435 + t400;
t546 = cos(qJ(5));
t437 = t477 * t547 - t524 * t543;
t542 = sin(qJ(5));
t650 = t437 * t542;
t383 = t546 * t435 + t650;
t676 = t370 * t383;
t564 = t435 * t542 - t546 * t437;
t675 = t370 * t564;
t669 = -t544 * t602 + t548 * t582;
t470 = qJD(4) - t669;
t459 = qJD(5) + t470;
t674 = t383 * t459;
t673 = t459 * t564;
t611 = qJD(4) * t547;
t672 = t544 * t611 - t671;
t663 = qJD(4) + qJD(5);
t670 = -t669 + t663;
t608 = qJD(1) * qJD(2);
t592 = t540 * t608;
t576 = t549 * t592;
t563 = t544 * t576;
t440 = qJD(3) * t477 + t563;
t668 = t440 * MDP(29) + (-t383 ^ 2 + t564 ^ 2) * MDP(26) - t383 * MDP(25) * t564;
t537 = t540 ^ 2;
t667 = -0.2e1 * t537 * t608;
t666 = MDP(5) * (t545 ^ 2 - t549 ^ 2);
t508 = t542 * t547 + t543 * t546;
t496 = t508 * t544;
t641 = t540 * t545;
t531 = pkin(7) * t641;
t660 = pkin(1) * t549;
t487 = t531 + (-pkin(2) - t660) * t541;
t498 = -t541 * t548 + t544 * t641;
t499 = t541 * t544 + t548 * t641;
t419 = pkin(3) * t498 - pkin(9) * t499 + t487;
t640 = t540 * t549;
t607 = pkin(7) * t640;
t661 = pkin(1) * t545;
t488 = t607 + (pkin(8) + t661) * t541;
t628 = t548 * t488 + t544 * t489;
t421 = -pkin(9) * t640 + t628;
t632 = t543 * t419 + t547 * t421;
t615 = qJD(3) * t544;
t658 = pkin(8) * t543;
t665 = -t547 * t678 + t615 * t658;
t490 = -pkin(7) * t602 + t549 * t606;
t561 = t540 * (pkin(2) * t545 - pkin(8) * t549);
t491 = qJD(1) * t561;
t626 = t548 * t490 + t544 * t491;
t424 = pkin(9) * t602 + t626;
t521 = -pkin(3) * t548 - pkin(9) * t544 - pkin(2);
t664 = t547 * t424 - t521 * t611 + t543 * t678;
t439 = qJD(3) * t669 + t548 * t576;
t577 = t545 * t592;
t381 = qJD(4) * t437 + t439 * t543 - t547 * t577;
t380 = t547 * t439 - t477 * t612 - t524 * t611 + t543 * t577;
t589 = t380 * t542 + t546 * t381;
t350 = -qJD(5) * t564 + t589;
t550 = qJD(1) ^ 2;
t662 = pkin(9) + pkin(10);
t659 = pkin(4) * t544;
t457 = -pkin(2) * t582 - t490;
t398 = -pkin(3) * t669 - t477 * pkin(9) + t457;
t412 = t548 * t458 + t544 * t469;
t401 = -pkin(9) * t524 + t412;
t361 = t547 * t398 - t401 * t543;
t356 = -pkin(10) * t437 + t361;
t353 = pkin(4) * t470 + t356;
t657 = t353 * t546;
t362 = t398 * t543 + t401 * t547;
t357 = -pkin(10) * t435 + t362;
t656 = t357 * t546;
t492 = qJD(2) * t561;
t484 = qJD(1) * t492;
t494 = (t541 * t660 - t531) * qJD(2);
t485 = qJD(1) * t494;
t579 = t458 * t613 + t469 * t615 - t548 * t484 + t544 * t485;
t369 = -pkin(3) * t577 + t579;
t655 = t369 * t543;
t654 = t369 * t547;
t653 = t380 * t543;
t652 = t435 * t470;
t651 = t437 * t470;
t649 = t440 * t543;
t648 = t440 * t547;
t647 = t440 * t548;
t646 = t669 * t524;
t645 = t669 * t543;
t644 = t477 * t524;
t558 = t524 * t544;
t643 = t524 * t548;
t642 = t537 * t550;
t639 = t543 * t544;
t638 = t544 * t547;
t637 = t547 * t548;
t427 = pkin(3) * t477 - pkin(9) * t669;
t635 = t547 * t411 + t543 * t427;
t507 = t542 * t543 - t546 * t547;
t634 = t463 * t542 - t464 * t546 - t496 * t663 - t507 * t613;
t610 = qJD(5) * t542;
t633 = -t610 * t639 + (t638 * t663 - t671) * t546 + t677 * t542;
t630 = t670 * t507;
t629 = t670 * t508;
t585 = -t544 * t490 + t491 * t548;
t423 = -pkin(3) * t602 - t585;
t627 = pkin(4) * t672 + pkin(8) * t613 - t423;
t533 = pkin(8) * t637;
t623 = t543 * t521 + t533;
t617 = qJD(2) * t545;
t616 = qJD(2) * t548;
t614 = qJD(3) * t547;
t609 = qJD(5) * t546;
t605 = t543 * t640;
t604 = t546 * t380 - t542 * t381 - t435 * t609;
t603 = qJD(4) * t662;
t600 = t540 * t617;
t599 = qJD(2) * t640;
t596 = t470 * t612;
t556 = -t458 * t615 + t469 * t613 + t544 * t484 + t548 * t485;
t368 = pkin(9) * t577 + t556;
t495 = (t541 * t661 + t607) * qJD(2);
t486 = qJD(1) * t495;
t376 = pkin(3) * t440 - pkin(9) * t439 + t486;
t348 = -qJD(4) * t362 - t368 * t543 + t547 * t376;
t341 = pkin(4) * t440 - pkin(10) * t380 + t348;
t347 = t547 * t368 + t543 * t376 + t398 * t611 - t401 * t612;
t342 = -pkin(10) * t381 + t347;
t591 = t546 * t341 - t542 * t342;
t355 = t357 * t610;
t590 = t542 * t341 - t355;
t588 = t547 * t419 - t421 * t543;
t586 = -t544 * t488 + t489 * t548;
t584 = t470 * t547;
t583 = qJD(5) * t353 + t342;
t580 = t537 * t545 * t549 * MDP(4);
t575 = MDP(15) * t602;
t574 = -t412 + (t612 - t645) * pkin(4);
t573 = pkin(1) * t667;
t420 = pkin(3) * t640 - t586;
t426 = t547 * t427;
t526 = t662 * t547;
t570 = pkin(4) * t477 + qJD(5) * t526 - t411 * t543 + t426 + (-pkin(10) * t669 + t603) * t547;
t525 = t662 * t543;
t569 = -pkin(10) * t645 + qJD(5) * t525 + t543 * t603 + t635;
t454 = -pkin(10) * t639 + t623;
t568 = -pkin(10) * t464 + qJD(5) * t454 - t424 * t543 + t601 * t659 - (-pkin(10) * t637 + t659) * qJD(3) - (-t533 + (pkin(10) * t544 - t521) * t543) * qJD(4) - t665;
t505 = t547 * t521;
t442 = -pkin(10) * t638 + t505 + (-pkin(4) - t658) * t548;
t567 = -qJD(5) * t442 - (-t544 * t614 - t548 * t612) * pkin(8) + t664 + t672 * pkin(10);
t345 = t353 * t542 + t656;
t451 = t499 * t547 - t605;
t359 = pkin(4) * t498 - pkin(10) * t451 + t588;
t450 = t499 * t543 + t547 * t640;
t363 = -pkin(10) * t450 + t632;
t566 = t359 * t546 - t363 * t542;
t565 = t359 * t542 + t363 * t546;
t402 = t546 * t450 + t451 * t542;
t403 = -t450 * t542 + t451 * t546;
t562 = -t488 * t613 - t489 * t615 + t492 * t548 - t544 * t494;
t560 = -t470 * t611 - t649;
t557 = -pkin(9) * t440 + t400 * t470;
t555 = -t488 * t615 + t489 * t613 + t544 * t492 + t548 * t494;
t372 = pkin(9) * t600 + t555;
t448 = qJD(3) * t499 + t544 * t599;
t449 = -qJD(3) * t498 + t548 * t599;
t391 = pkin(3) * t448 - pkin(9) * t449 + t495;
t554 = t547 * t372 + t543 * t391 + t419 * t611 - t421 * t612;
t349 = -t437 * t610 + t604;
t552 = pkin(1) * (-t541 * t608 + t642);
t373 = -pkin(3) * t600 - t562;
t551 = -qJD(4) * t632 - t372 * t543 + t547 * t391;
t339 = -qJD(5) * t345 + t591;
t536 = -pkin(4) * t547 - pkin(3);
t513 = (pkin(4) * t543 + pkin(8)) * t544;
t497 = t507 * t544;
t422 = t440 * t498;
t397 = -qJD(4) * t450 + t449 * t547 + t543 * t600;
t396 = -qJD(4) * t605 + t449 * t543 + t499 * t611 - t547 * t600;
t387 = pkin(4) * t450 + t420;
t358 = pkin(4) * t396 + t373;
t354 = pkin(4) * t381 + t369;
t352 = qJD(5) * t403 + t546 * t396 + t397 * t542;
t351 = -qJD(5) * t402 - t396 * t542 + t397 * t546;
t346 = -pkin(10) * t396 + t554;
t344 = -t357 * t542 + t657;
t343 = pkin(4) * t448 - pkin(10) * t397 + t551;
t338 = t583 * t546 + t590;
t1 = [(t448 * t524 + (t440 * t549 + (-qJD(1) * t498 + t669) * t617) * t540) * MDP(14) + (-t562 * t524 - t495 * t669 + t487 * t440 + t486 * t498 + t457 * t448 + (t579 * t549 + (qJD(1) * t586 + t411) * t617) * t540) * MDP(16) + (-t439 * t498 - t440 * t499 - t448 * t477 + t449 * t669) * MDP(12) + (t348 * t498 + t361 * t448 + t369 * t450 + t373 * t435 + t420 * t381 + t400 * t396 + t440 * t588 + t470 * t551) * MDP(23) + (-t485 * t541 - t494 * t582 + t549 * t573) * MDP(10) + (-t486 * t541 - t495 * t582 + t545 * t573) * MDP(9) + (t358 * t383 + t387 * t350 + t354 * t402 + t370 * t352 + (-qJD(5) * t565 + t343 * t546 - t346 * t542) * t459 + t566 * t440 + t339 * t498 + t344 * t448) * MDP(30) + t666 * t667 + (-t358 * t564 + t387 * t349 + t354 * t403 + t370 * t351 - (qJD(5) * t566 + t343 * t542 + t346 * t546) * t459 - t565 * t440 - t338 * t498 - t345 * t448) * MDP(31) + (t349 * t403 - t351 * t564) * MDP(25) + (-t349 * t402 - t350 * t403 - t351 * t383 + t352 * t564) * MDP(26) + (t349 * t498 + t351 * t459 + t403 * t440 - t448 * t564) * MDP(27) + 0.2e1 * t580 * t608 + (-t347 * t498 - t362 * t448 + t369 * t451 + t373 * t437 + t420 * t380 + t400 * t397 - t440 * t632 - t470 * t554) * MDP(24) + (MDP(6) * t599 - MDP(7) * t600) * (qJD(2) + 0.2e1 * t619) + (-t449 * t524 + (-t439 * t549 + (qJD(1) * t499 + t477) * t617) * t540) * MDP(13) + (-t524 * t540 - t537 * t618) * MDP(15) * t617 + (t555 * t524 + t495 * t477 + t487 * t439 + t486 * t499 + t457 * t449 + (t556 * t549 + (-qJD(1) * t628 - t412) * t617) * t540) * MDP(17) + (t380 * t451 + t397 * t437) * MDP(18) + (-t380 * t450 - t381 * t451 - t396 * t437 - t397 * t435) * MDP(19) + (t448 * t459 + t422) * MDP(29) + (t448 * t470 + t422) * MDP(22) + (-t350 * t498 - t352 * t459 - t383 * t448 - t402 * t440) * MDP(28) + (-t381 * t498 - t396 * t470 - t435 * t448 - t440 * t450) * MDP(21) + (t380 * t498 + t397 * t470 + t437 * t448 + t440 * t451) * MDP(20) + (t439 * t499 + t449 * t477) * MDP(11); (-pkin(2) * t440 - t486 * t548 + t585 * t524 + t493 * t669 + (pkin(8) * t643 + t457 * t544) * qJD(3) + (-t411 * t545 + (-pkin(8) * t617 - t457 * t549) * t544) * t620) * MDP(16) + (t524 * t615 + (-t549 * t558 + (-t669 + t616) * t545) * t620) * MDP(14) + (t381 * t548 + t671 * t470 + (t435 * t524 + t560) * t544) * MDP(21) + ((t439 - t646) * t548 + (-t440 + t644) * t544) * MDP(12) + (-t400 * t463 - t423 * t435 + t505 * t440 + ((-qJD(4) * t521 + t424) * t543 + t665) * t470 + (t400 * t543 * qJD(3) - t348 + (qJD(3) * t435 + t560) * pkin(8)) * t548 + (pkin(8) * t381 - t361 * t524 + t400 * t611 + t655) * t544) * MDP(23) + (-t623 * t440 - t423 * t437 - t400 * t464 + t664 * t470 + (t400 * t614 + t347 + (qJD(3) * t437 + t596) * pkin(8)) * t548 + (-t400 * t612 + t654 + t524 * t362 + (t470 * t614 + t380) * pkin(8)) * t544) * MDP(24) + (-pkin(7) * t576 + t493 * t582 + t545 * t552) * MDP(9) + (pkin(7) * t577 + t490 * t582 + t549 * t552) * MDP(10) + (t380 * t638 + t437 * t677) * MDP(18) + t642 * t666 + (-t349 * t497 - t564 * t634) * MDP(25) + (t513 * t349 - t354 * t497 - (t442 * t542 + t454 * t546) * t440 + t338 * t548 + (t542 * t568 + t546 * t567) * t459 - t627 * t564 + t634 * t370 + t345 * t558) * MDP(31) + (-t349 * t496 + t350 * t497 - t383 * t634 + t564 * t633) * MDP(26) + (-t349 * t548 - t440 * t497 + t459 * t634 + t558 * t564) * MDP(27) + (t513 * t350 + t354 * t496 + (t442 * t546 - t454 * t542) * t440 - t339 * t548 + (t542 * t567 - t546 * t568) * t459 + t627 * t383 + t633 * t370 - t344 * t558) * MDP(30) + (-t470 * t558 - t647) * MDP(22) + (-t459 * t558 - t647) * MDP(29) + (-t524 * t613 + (t524 * t636 + (qJD(2) * t544 - t477) * t545) * t620) * MDP(13) + (t350 * t548 + t383 * t558 - t440 * t496 - t459 * t633) * MDP(28) + (-t380 * t548 + t571 * t470 + (-t437 * t524 - t596 + t648) * t544) * MDP(20) + t524 * t575 + (t439 * t544 - t477 * t643) * MDP(11) + (-pkin(2) * t439 + t486 * t544 - t626 * t524 - t493 * t477 + (-pkin(8) * t558 + t457 * t548) * qJD(3) + (-t457 * t636 + (-pkin(8) * t616 + t412) * t545) * t620) * MDP(17) + (-t580 + (-MDP(6) * t549 + MDP(7) * t545) * t540 * t541) * t550 + (t435 * t464 + t437 * t463 + (-t435 * t547 - t437 * t543) * t613 + (-t653 - t381 * t547 + (t435 * t543 - t437 * t547) * qJD(4)) * t544) * MDP(19); -t669 ^ 2 * MDP(12) + (t439 + t646) * MDP(13) + (-t563 - t644) * MDP(14) + qJD(2) * t575 + (-t412 * t524 - t579) * MDP(16) + (-t411 * t524 - t457 * t669 - t556) * MDP(17) + (t437 * t584 + t653) * MDP(18) + ((t380 - t652) * t547 + (-t381 - t651) * t543) * MDP(19) + (t470 * t584 + t649) * MDP(20) + (-t470 ^ 2 * t543 + t648) * MDP(21) + (-pkin(3) * t381 - t654 - t412 * t435 + (-pkin(9) * t611 - t426) * t470 + (t411 * t470 + t557) * t543) * MDP(23) + (-pkin(3) * t380 + t655 - t412 * t437 + (pkin(9) * t612 + t635) * t470 + t557 * t547) * MDP(24) + (t349 * t508 + t564 * t630) * MDP(25) + (-t349 * t507 - t350 * t508 + t383 * t630 + t564 * t629) * MDP(26) + (t440 * t508 - t459 * t630) * MDP(27) + (-t440 * t507 - t459 * t629) * MDP(28) + (t536 * t350 + t354 * t507 + (-t525 * t546 - t526 * t542) * t440 + (t542 * t569 - t546 * t570) * t459 + t574 * t383 + t629 * t370) * MDP(30) + (t536 * t349 + t354 * t508 - (-t525 * t542 + t526 * t546) * t440 + (t542 * t570 + t546 * t569) * t459 - t574 * t564 - t630 * t370) * MDP(31) + (-MDP(11) * t669 + t477 * MDP(12) - MDP(14) * qJD(3) - t457 * MDP(16) - t437 * MDP(20) + t435 * MDP(21) - t470 * MDP(22) - t361 * MDP(23) + t362 * MDP(24) + MDP(27) * t564 + t383 * MDP(28) - t459 * MDP(29) - t344 * MDP(30) + t345 * MDP(31)) * t477; t437 * t435 * MDP(18) + (-t435 ^ 2 + t437 ^ 2) * MDP(19) + (t380 + t652) * MDP(20) + (-t381 + t651) * MDP(21) + t440 * MDP(22) + (t362 * t470 - t400 * t437 + t348) * MDP(23) + (t361 * t470 + t400 * t435 - t347) * MDP(24) + (t349 + t674) * MDP(27) + (-t350 - t673) * MDP(28) + (t675 - (-t356 * t542 - t656) * t459 + (-t383 * t437 + t546 * t440 - t459 * t610) * pkin(4) + t339) * MDP(30) + (t676 + t355 + (-t357 * t459 - t341) * t542 + (t356 * t459 - t583) * t546 + (t437 * t564 - t542 * t440 - t459 * t609) * pkin(4)) * MDP(31) + t668; (t604 + t674) * MDP(27) + (-t589 - t673) * MDP(28) + (t345 * t459 + t591 + t675) * MDP(30) + (-t546 * t342 + t344 * t459 - t590 + t676) * MDP(31) + (-MDP(27) * t650 + MDP(28) * t564 - MDP(30) * t345 - MDP(31) * t657) * qJD(5) + t668;];
tauc = t1;
