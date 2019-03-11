% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:45:36
% EndTime: 2019-03-09 00:45:47
% DurationCPUTime: 6.78s
% Computational Cost: add. (5695->463), mult. (13789->644), div. (0->0), fcn. (10832->12), ass. (0->221)
t667 = qJD(5) + qJD(6);
t543 = sin(qJ(4));
t548 = cos(qJ(3));
t664 = cos(qJ(4));
t599 = qJD(2) * t664;
t544 = sin(qJ(3));
t623 = qJD(2) * t544;
t684 = -t543 * t623 + t548 * t599;
t687 = -t684 + t667;
t545 = sin(qJ(2));
t539 = sin(pkin(6));
t626 = qJD(1) * t539;
t606 = t545 * t626;
t660 = qJD(3) * pkin(3);
t558 = t544 * t660 - t606;
t536 = qJD(3) + qJD(4);
t565 = -t543 * t544 + t548 * t664;
t469 = t536 * t565;
t641 = t543 * t548;
t505 = t544 * t664 + t641;
t470 = t536 * t505;
t686 = pkin(4) * t470 - pkin(10) * t469 + t558;
t666 = pkin(8) + pkin(9);
t609 = qJD(3) * t666;
t507 = t544 * t609;
t508 = t548 * t609;
t516 = t666 * t544;
t518 = t666 * t548;
t669 = -t664 * t516 - t543 * t518;
t421 = qJD(4) * t669 - t664 * t507 - t543 * t508;
t549 = cos(qJ(2));
t605 = t549 * t626;
t479 = t565 * t605;
t685 = -t421 + t479;
t594 = qJD(2) * t666 + t606;
t540 = cos(pkin(6));
t625 = qJD(1) * t540;
t477 = t544 * t625 + t548 * t594;
t463 = t543 * t477;
t476 = -t594 * t544 + t548 * t625;
t414 = t476 * t664 - t463;
t598 = qJD(4) * t664;
t678 = -pkin(3) * t598 + t414;
t499 = -qJD(2) * t641 - t544 * t599;
t459 = -pkin(4) * t499 - pkin(10) * t684;
t443 = pkin(3) * t623 + t459;
t542 = sin(qJ(5));
t547 = cos(qJ(5));
t683 = -t547 * t443 + t542 * t678;
t465 = t476 + t660;
t411 = t465 * t664 - t463;
t401 = -t536 * pkin(4) - t411;
t480 = -t499 * t542 - t547 * t536;
t389 = t480 * pkin(5) + t401;
t546 = cos(qJ(6));
t541 = sin(qJ(6));
t569 = t499 * t547 - t536 * t542;
t650 = t569 * t541;
t417 = t546 * t480 - t650;
t682 = t389 * t417;
t493 = qJD(5) - t684;
t490 = qJD(6) + t493;
t681 = t417 * t490;
t570 = t480 * t541 + t546 * t569;
t680 = t490 * t570;
t643 = t541 * t547;
t504 = t542 * t546 + t643;
t679 = t687 * t504;
t502 = t541 * t542 - t546 * t547;
t631 = t687 * t502;
t456 = t684 * t536;
t620 = qJD(5) * t547;
t621 = qJD(5) * t542;
t408 = t547 * t456 + t499 * t621 + t536 * t620;
t457 = t470 * qJD(2);
t624 = qJD(2) * t539;
t597 = qJD(1) * t624;
t586 = t549 * t597;
t434 = qJD(3) * t476 + t548 * t586;
t435 = -qJD(3) * t477 - t544 * t586;
t622 = qJD(4) * t543;
t368 = t434 * t664 + t543 * t435 + t465 * t598 - t477 * t622;
t464 = t664 * t477;
t412 = t543 * t465 + t464;
t402 = pkin(10) * t536 + t412;
t534 = -pkin(3) * t548 - pkin(2);
t491 = qJD(2) * t534 - t605;
t429 = -pkin(4) * t684 + pkin(10) * t499 + t491;
t384 = t402 * t547 + t429 * t542;
t615 = qJD(2) * qJD(3);
t596 = t544 * t615;
t494 = pkin(3) * t596 + t545 * t597;
t393 = pkin(4) * t457 - pkin(10) * t456 + t494;
t392 = t547 * t393;
t554 = -qJD(5) * t384 - t368 * t542 + t392;
t348 = pkin(5) * t457 - pkin(11) * t408 + t554;
t409 = -qJD(5) * t569 + t456 * t542;
t560 = t547 * t368 + t542 * t393 - t402 * t621 + t429 * t620;
t349 = -pkin(11) * t409 + t560;
t593 = t546 * t348 - t541 * t349;
t677 = t389 * t570 + t593;
t676 = t457 * MDP(30) + (-t417 ^ 2 + t570 ^ 2) * MDP(27) - t417 * MDP(26) * t570;
t675 = MDP(5) * t548;
t674 = MDP(6) * (t544 ^ 2 - t548 ^ 2);
t450 = t504 * t505;
t673 = t479 * t542 + t547 * t686;
t413 = t543 * t476 + t464;
t584 = pkin(3) * t622 - t413;
t484 = -t543 * t516 + t518 * t664;
t633 = qJD(4) * t484 - t505 * t605 - t543 * t507 + t508 * t664;
t460 = -pkin(4) * t565 - pkin(10) * t505 + t534;
t672 = -t460 * t620 + t484 * t621 - t542 * t686 + t685 * t547;
t649 = t684 * t542;
t671 = (t621 - t649) * pkin(5);
t670 = t542 * t443 + t547 * t678;
t668 = t544 * MDP(10) + t548 * MDP(11);
t591 = t408 * t541 + t546 * t409;
t362 = -qJD(6) * t570 + t591;
t665 = -pkin(10) - pkin(11);
t663 = t547 * pkin(5);
t531 = pkin(3) * t543 + pkin(10);
t662 = -pkin(11) - t531;
t661 = qJD(2) * pkin(2);
t383 = -t402 * t542 + t547 * t429;
t371 = pkin(11) * t569 + t383;
t363 = pkin(5) * t493 + t371;
t659 = t363 * t546;
t372 = -pkin(11) * t480 + t384;
t658 = t372 * t546;
t657 = t408 * t542;
t655 = t457 * t547;
t654 = t469 * t542;
t653 = t469 * t547;
t652 = t480 * t493;
t651 = t569 * t493;
t648 = t684 * t547;
t647 = t505 * t542;
t646 = t505 * t547;
t645 = t539 * t545;
t644 = t539 * t549;
t642 = t542 * t457;
t550 = qJD(3) ^ 2;
t640 = t544 * t550;
t475 = t547 * t484;
t639 = t548 * t550;
t638 = -pkin(11) * t653 + pkin(5) * t470 - t421 * t542 + (-t475 + (pkin(11) * t505 - t460) * t542) * qJD(5) + t673;
t564 = t505 * t620 + t654;
t637 = pkin(11) * t564 + t672;
t636 = pkin(5) * t564 + t633;
t635 = t547 * t411 + t542 * t459;
t630 = t542 * t460 + t475;
t629 = t671 + t584;
t619 = qJD(6) * t541;
t618 = qJD(6) * t546;
t614 = pkin(11) * t649;
t610 = t546 * t408 - t541 * t409 - t480 * t618;
t608 = qJD(5) * t665;
t603 = t545 * t624;
t602 = t549 * t624;
t600 = t505 * t621;
t398 = t401 * t620;
t595 = qJD(5) * t662;
t370 = t372 * t619;
t592 = t541 * t348 - t370;
t590 = -t411 * t542 + t547 * t459;
t589 = t493 * t547;
t588 = qJD(6) * t363 + t349;
t369 = t543 * t434 - t664 * t435 + t465 * t622 + t477 * t598;
t532 = -pkin(3) * t664 - pkin(4);
t585 = -t499 * pkin(5) - pkin(11) * t648;
t583 = -t412 + t671;
t582 = t369 * t542 - t384 * t499 + t398;
t535 = t547 * pkin(11);
t501 = t531 * t547 + t535;
t580 = qJD(6) * t501 - t547 * t595 + t585 - t683;
t517 = pkin(10) * t547 + t535;
t579 = qJD(6) * t517 - t547 * t608 + t585 + t590;
t500 = t662 * t542;
t578 = -qJD(6) * t500 - t542 * t595 - t614 + t670;
t515 = t665 * t542;
t577 = -qJD(6) * t515 - t542 * t608 - t614 + t635;
t352 = t363 * t541 + t658;
t454 = t547 * t460;
t390 = -pkin(5) * t565 - pkin(11) * t646 - t484 * t542 + t454;
t396 = -pkin(11) * t647 + t630;
t574 = t390 * t541 + t396 * t546;
t573 = -t401 * t684 - t457 * t531;
t495 = t540 * t548 - t544 * t645;
t496 = t540 * t544 + t548 * t645;
t445 = t543 * t495 + t496 * t664;
t427 = -t445 * t542 - t547 * t644;
t568 = -t445 * t547 + t542 * t644;
t572 = t427 * t546 + t541 * t568;
t571 = t427 * t541 - t546 * t568;
t567 = -t369 * t547 + t383 * t499 + t401 * t621;
t566 = t495 * t664 - t543 * t496;
t563 = -t600 + t653;
t561 = t491 * t499 - t369;
t361 = t569 * t619 + t610;
t351 = -t372 * t541 + t659;
t360 = pkin(5) * t409 + t369;
t557 = t351 * t499 + t360 * t502 + t389 * t679;
t556 = -t352 * t499 + t360 * t504 - t389 * t631;
t555 = -0.2e1 * qJD(3) * t661;
t553 = -t491 * t684 - t368;
t552 = (-t361 * t502 - t362 * t504 + t417 * t631 + t570 * t679) * MDP(27) + (t361 * t504 + t570 * t631) * MDP(26) + ((t408 - t652) * t547 + (-t409 + t651) * t542) * MDP(20) + (t457 * t504 - t490 * t631 - t499 * t570) * MDP(28) + (-t417 * t499 - t457 * t502 - t490 * t679) * MDP(29) + (-t569 * t589 + t657) * MDP(19) + (-t493 ^ 2 * t542 - t480 * t499 + t655) * MDP(22) + (t493 * t589 - t499 * t569 + t642) * MDP(21) + t456 * MDP(14) + (t499 ^ 2 - t684 ^ 2) * MDP(13) + (MDP(12) * t684 + MDP(23) * t493 + MDP(30) * t490) * t499 + (-t684 * MDP(14) + (-qJD(2) * t505 - t499) * MDP(15)) * t536;
t551 = qJD(2) ^ 2;
t533 = -pkin(4) - t663;
t514 = t532 - t663;
t474 = -qJD(3) * t496 - t544 * t602;
t473 = qJD(3) * t495 + t548 * t602;
t451 = t502 * t505;
t440 = pkin(5) * t647 - t669;
t431 = t457 * t565;
t388 = qJD(4) * t445 + t543 * t473 - t474 * t664;
t387 = qJD(4) * t566 + t473 * t664 + t543 * t474;
t374 = t469 * t643 - t541 * t600 - t619 * t647 + (t646 * t667 + t654) * t546;
t373 = -t450 * t667 - t502 * t469;
t365 = qJD(5) * t568 - t387 * t542 + t547 * t603;
t364 = qJD(5) * t427 + t387 * t547 + t542 * t603;
t1 = [(t365 * t493 + t388 * t480 - t409 * t566 + t427 * t457) * MDP(24) + (-t364 * t493 - t388 * t569 - t408 * t566 + t457 * t568) * MDP(25) + ((-qJD(6) * t571 - t364 * t541 + t365 * t546) * t490 + t572 * t457 + t388 * t417 - t566 * t362) * MDP(31) + (-(qJD(6) * t572 + t364 * t546 + t365 * t541) * t490 - t571 * t457 - t388 * t570 - t566 * t361) * MDP(32) + (-MDP(17) * t388 - MDP(18) * t387) * t536 + (MDP(10) * t474 - MDP(11) * t473) * qJD(3) + ((-MDP(17) * t457 - MDP(18) * t456) * t549 + ((-MDP(17) * t684 - t499 * MDP(18)) * t545 - t668 * t549 * qJD(3)) * qJD(2) + (-t549 * MDP(4) + (-MDP(10) * t548 + MDP(11) * t544 - MDP(3)) * t545) * t551) * t539; (t369 * t646 - t384 * t470 + t563 * t401 - t408 * t669 - t630 * t457 + t493 * t672 + t560 * t565 - t569 * t633) * MDP(25) + (t454 * t457 - (-t402 * t620 + t392) * t565 + t383 * t470 - t669 * t409 + t505 * t398 + (-t484 * t620 + t673) * t493 + t633 * t480 + ((-qJD(5) * t460 - t421) * t493 - t484 * t457 - (-qJD(5) * t429 - t368) * t565 + t369 * t505 + t401 * t469) * t542) * MDP(24) + ((-t480 * t547 + t542 * t569) * t469 + (-t657 - t409 * t547 + (t480 * t542 + t547 * t569) * qJD(5)) * t505) * MDP(20) + (t408 * t646 - t563 * t569) * MDP(19) + (-t408 * t565 + t457 * t646 - t470 * t569 + t493 * t563) * MDP(21) + (-t361 * t565 + t373 * t490 - t451 * t457 - t470 * t570) * MDP(28) + (-t574 * t457 + (t546 * t588 + t592) * t565 - t352 * t470 + t440 * t361 - t360 * t451 + t389 * t373 + ((-qJD(6) * t390 + t637) * t546 + (qJD(6) * t396 - t638) * t541) * t490 - t636 * t570) * MDP(32) + (t362 * t565 - t374 * t490 - t417 * t470 - t450 * t457) * MDP(29) + (t409 * t565 - t470 * t480 - t493 * t564 - t505 * t642) * MDP(22) + ((t390 * t546 - t396 * t541) * t457 - t593 * t565 + t351 * t470 + t440 * t362 + t360 * t450 + t389 * t374 + (t541 * t637 + t546 * t638) * t490 + t636 * t417 + (t352 * t565 - t490 * t574) * qJD(6)) * MDP(31) + (t470 * t493 - t431) * MDP(23) + (t470 * t490 - t431) * MDP(30) + (-t361 * t450 + t362 * t451 - t373 * t417 + t374 * t570) * MDP(27) + (-t361 * t451 - t373 * t570) * MDP(26) + (t456 * t534 + t469 * t491 + t494 * t505 - t558 * t499) * MDP(18) + (t469 * MDP(14) - t470 * MDP(15) - t633 * MDP(17) + MDP(18) * t685) * t536 - 0.2e1 * t615 * t674 + 0.2e1 * t596 * t675 + (t456 * t505 - t469 * t499) * MDP(12) + (-pkin(8) * t639 + t544 * t555) * MDP(10) - MDP(8) * t640 + (pkin(8) * t640 + t548 * t555) * MDP(11) + (t457 * t534 + t470 * t491 - t494 * t565 - t558 * t684) * MDP(17) + (t456 * t565 - t457 * t505 + t469 * t684 + t470 * t499) * MDP(13) + MDP(7) * t639; t552 + (t414 * t536 + (t499 * t623 - t536 * t598) * pkin(3) + t553) * MDP(18) + (t413 * t536 + (-t536 * t622 + t623 * t684) * pkin(3) + t561) * MDP(17) + (t532 * t408 + t573 * t547 - t584 * t569 + (t531 * t621 + t670) * t493 + t582) * MDP(25) + (t532 * t409 + t573 * t542 + t584 * t480 + (-t531 * t620 + t683) * t493 + t567) * MDP(24) + ((t500 * t546 - t501 * t541) * t457 + t514 * t362 + (t541 * t578 - t546 * t580) * t490 + t629 * t417 + t557) * MDP(31) + (-(t500 * t541 + t501 * t546) * t457 + t514 * t361 + (t541 * t580 + t546 * t578) * t490 - t629 * t570 + t556) * MDP(32) + t668 * qJD(2) * t661 + (-t544 * t675 + t674) * t551; t552 + (-pkin(4) * t409 - t590 * t493 - t412 * t480 - t401 * t649 + (-t493 * t620 - t642) * pkin(10) + t567) * MDP(24) + (-pkin(4) * t408 + t635 * t493 + t412 * t569 - t401 * t648 + (t493 * t621 - t655) * pkin(10) + t582) * MDP(25) + ((t515 * t546 - t517 * t541) * t457 + t533 * t362 + (t541 * t577 - t546 * t579) * t490 + t583 * t417 + t557) * MDP(31) + (-(t515 * t541 + t517 * t546) * t457 + t533 * t361 + (t541 * t579 + t546 * t577) * t490 - t583 * t570 + t556) * MDP(32) + (t411 * t536 + t553) * MDP(18) + (t412 * t536 + t561) * MDP(17); -t569 * t480 * MDP(19) + (-t480 ^ 2 + t569 ^ 2) * MDP(20) + (t408 + t652) * MDP(21) + (-t409 - t651) * MDP(22) + t457 * MDP(23) + (t384 * t493 + t401 * t569 + t554) * MDP(24) + (t383 * t493 + t401 * t480 - t560) * MDP(25) + (t361 + t681) * MDP(28) + (-t362 - t680) * MDP(29) + (-(-t371 * t541 - t658) * t490 - t352 * qJD(6) + (t417 * t569 + t457 * t546 - t490 * t619) * pkin(5) + t677) * MDP(31) + (t682 + t370 + (-t372 * t490 - t348) * t541 + (t371 * t490 - t588) * t546 + (-t457 * t541 - t490 * t618 - t569 * t570) * pkin(5)) * MDP(32) + t676; (t610 + t681) * MDP(28) + (-t591 - t680) * MDP(29) + (t352 * t490 + t677) * MDP(31) + (-t546 * t349 + t351 * t490 - t592 + t682) * MDP(32) + (MDP(28) * t650 + MDP(29) * t570 - MDP(31) * t352 - MDP(32) * t659) * qJD(6) + t676;];
tauc  = t1;
