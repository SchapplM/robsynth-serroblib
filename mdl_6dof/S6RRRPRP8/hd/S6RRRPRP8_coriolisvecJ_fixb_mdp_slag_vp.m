% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:19:52
% EndTime: 2019-03-09 17:20:06
% DurationCPUTime: 7.97s
% Computational Cost: add. (4733->537), mult. (10968->695), div. (0->0), fcn. (7050->6), ass. (0->221)
t540 = cos(qJ(3));
t618 = qJD(3) * t540;
t541 = cos(qJ(2));
t623 = qJD(1) * t541;
t691 = -t540 * t623 + t618;
t537 = sin(qJ(3));
t597 = t537 * t623;
t619 = qJD(3) * t537;
t690 = t597 - t619;
t538 = sin(qJ(2));
t491 = -pkin(2) * t541 - pkin(8) * t538 - pkin(1);
t460 = t491 * qJD(1);
t526 = pkin(7) * t623;
t498 = qJD(2) * pkin(8) + t526;
t419 = t540 * t460 - t537 * t498;
t610 = qJD(4) - t419;
t611 = t540 * qJD(2);
t624 = qJD(1) * t538;
t473 = t537 * t624 - t611;
t596 = t540 * t624;
t622 = qJD(2) * t537;
t475 = t596 + t622;
t536 = sin(qJ(5));
t539 = cos(qJ(5));
t560 = -t539 * t473 + t475 * t536;
t411 = t560 ^ 2;
t612 = t538 * MDP(26);
t583 = qJD(2) * t612;
t416 = t473 * t536 + t475 * t539;
t668 = t416 ^ 2;
t689 = -(t411 - t668) * MDP(23) - qJD(1) * t583;
t666 = pkin(8) - pkin(9);
t500 = t666 * t540;
t664 = pkin(7) * t537;
t598 = -pkin(3) - t664;
t647 = t540 * t541;
t547 = -pkin(9) * t647 + (-pkin(4) + t598) * t538;
t564 = pkin(2) * t538 - pkin(8) * t541;
t485 = t564 * qJD(1);
t656 = t485 * t540;
t688 = qJD(1) * t547 - qJD(3) * t500 - t656;
t456 = t537 * t485;
t632 = qJ(4) * t624 + t456;
t649 = t538 * t540;
t650 = t537 * t541;
t682 = (-t649 * pkin(7) + pkin(9) * t650) * qJD(1) + t632 + t666 * t619;
t686 = pkin(9) * t475 - t610;
t685 = qJ(6) * t416;
t497 = -qJD(2) * pkin(2) + pkin(7) * t624;
t408 = t473 * pkin(3) - t475 * qJ(4) + t497;
t387 = -pkin(4) * t473 - t408;
t684 = t387 * t416;
t516 = -qJD(3) + t623;
t670 = qJD(5) + t516;
t683 = t416 * t670;
t477 = t536 * t537 + t539 * t540;
t422 = t477 * qJD(5) - t536 * t619 - t539 * t618;
t555 = t541 * t477;
t638 = qJD(1) * t555 + t422;
t615 = qJD(5) * t539;
t616 = qJD(5) * t536;
t637 = t691 * t536 + t537 * t615 + t690 * t539 - t540 * t616;
t681 = -t416 * MDP(22) - MDP(24) * t670 - t387 * MDP(28);
t608 = qJD(1) * qJD(2);
t680 = -0.2e1 * t608;
t679 = MDP(4) * t538;
t533 = t538 ^ 2;
t678 = MDP(5) * (-t541 ^ 2 + t533);
t677 = qJ(6) * t560;
t488 = t564 * qJD(2);
t462 = qJD(1) * t488;
t587 = t538 * t608;
t569 = pkin(7) * t587;
t568 = -t460 * t619 + t540 * t462 - t498 * t618 + t537 * t569;
t570 = pkin(3) * t587;
t377 = -t568 - t570;
t420 = t537 * t460 + t540 * t498;
t503 = t516 * qJ(4);
t402 = -t503 + t420;
t676 = t402 * t516 + t377;
t518 = pkin(7) * t650;
t532 = t541 * pkin(3);
t409 = pkin(4) * t541 + t518 + t532 + (-pkin(9) * t538 - t491) * t540;
t519 = pkin(7) * t647;
t628 = t537 * t491 + t519;
t433 = -qJ(4) * t541 + t628;
t652 = t537 * t538;
t418 = pkin(9) * t652 + t433;
t639 = t536 * t409 + t539 * t418;
t395 = pkin(9) * t473 + t420;
t667 = pkin(3) + pkin(4);
t581 = -qJ(4) * t536 - t539 * t667;
t675 = -qJD(5) * t581 + t536 * t395 + t539 * t686;
t489 = qJ(4) * t539 - t536 * t667;
t674 = -qJD(5) * t489 - t539 * t395 + t536 * t686;
t673 = t688 * t539;
t499 = t666 * t537;
t630 = t536 * t499 + t539 * t500;
t672 = -t499 * t615 + t500 * t616 + t536 * t688 + t682 * t539;
t671 = t691 * qJ(4) + t537 * qJD(4) + t526;
t515 = qJ(4) * t649;
t602 = t667 * t537;
t567 = -pkin(7) - t602;
t432 = t567 * t538 + t515;
t590 = t538 * t618;
t620 = qJD(2) * t541;
t594 = t537 * t620;
t549 = t590 + t594;
t607 = qJD(2) * qJD(3);
t436 = qJD(1) * t549 + t537 * t607;
t380 = t516 * t667 - t686;
t385 = t395 - t503;
t354 = t380 * t536 + t385 * t539;
t501 = t516 * qJD(4);
t509 = qJ(4) * t587;
t557 = t460 * t618 + t537 * t462 - t498 * t619;
t546 = t540 * t569 - t557;
t370 = -t501 + t509 - t546;
t362 = pkin(9) * t436 + t370;
t586 = t541 * t608;
t592 = t538 * t619;
t435 = qJD(1) * t592 + (-t586 - t607) * t540;
t363 = pkin(9) * t435 - t587 * t667 - t568;
t580 = t536 * t362 - t539 * t363;
t669 = t354 * qJD(5) + t580;
t529 = t537 * qJ(4);
t593 = t541 * t611;
t629 = qJ(4) * t593 + qJD(4) * t649;
t381 = qJD(3) * (-t540 * t667 - t529) * t538 + t567 * t620 + t629;
t665 = pkin(3) * t516;
t375 = t436 * pkin(3) + pkin(7) * t586 + t435 * qJ(4) - t475 * qJD(4);
t662 = t375 * t537;
t661 = t375 * t540;
t659 = t435 * t537;
t658 = t473 * t516;
t657 = t475 * t516;
t655 = t497 * t537;
t654 = t497 * t540;
t653 = t516 * t540;
t651 = t537 * t539;
t543 = qJD(2) ^ 2;
t648 = t538 * t543;
t646 = t541 * t543;
t544 = qJD(1) ^ 2;
t645 = t541 * t544;
t353 = t539 * t380 - t385 * t536;
t349 = t353 - t685;
t348 = pkin(5) * t670 + t349;
t644 = t348 - t349;
t643 = -t637 * qJ(6) - qJD(6) * t477 - t672;
t559 = t536 * t540 - t651;
t642 = pkin(5) * t624 + t638 * qJ(6) - qJD(5) * t630 + qJD(6) * t559 + t682 * t536 - t673;
t636 = t690 * pkin(3) + t671;
t635 = -qJD(3) * t602 + t597 * t667 + t671;
t634 = -t675 - t685;
t633 = t674 + t677;
t631 = t537 * t488 + t491 * t618;
t426 = t475 * pkin(3) + t473 * qJ(4);
t621 = qJD(2) * t538;
t617 = qJD(5) * t416;
t613 = t538 * MDP(15);
t606 = pkin(8) * t516 * t537;
t605 = pkin(8) * t653;
t490 = -t540 * pkin(3) - pkin(2) - t529;
t604 = pkin(8) * t621;
t603 = pkin(8) * t611;
t601 = t539 * t362 + t536 * t363 + t380 * t615;
t600 = -t539 * t435 + t536 * t436 + t473 * t615;
t599 = qJ(4) * t621 + t631;
t591 = t541 * t619;
t584 = qJD(2) * t613;
t582 = pkin(1) * t680;
t563 = -qJD(3) * t519 + t488 * t540 - t491 * t619;
t374 = pkin(9) * t592 + qJD(2) * t547 - t563;
t376 = (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t649 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t537) * t541 + t599;
t579 = t539 * t374 - t376 * t536;
t577 = t539 * t409 - t418 * t536;
t576 = -t435 * t536 - t539 * t436;
t575 = t539 * t499 - t500 * t536;
t574 = t491 * t540 - t518;
t573 = t670 ^ 2;
t467 = t540 * pkin(4) - t490;
t572 = t473 + t611;
t571 = -t475 + t622;
t397 = -pkin(4) * t475 - t426;
t565 = t598 * t538;
t400 = t665 + t610;
t562 = t400 * t540 - t402 * t537;
t558 = qJD(1) * t533 - t516 * t541;
t556 = -t420 * t516 + t568;
t554 = -t385 * t616 + t601;
t553 = t536 * t374 + t539 * t376 + t409 * t615 - t418 * t616;
t367 = t475 * t616 - t600;
t364 = -pkin(4) * t436 - t375;
t368 = t576 + t617;
t347 = pkin(5) * t368 + t364;
t545 = -t419 * t516 + t546;
t484 = -pkin(5) + t581;
t451 = t477 * t538;
t450 = t536 * t649 - t538 * t651;
t443 = -t515 + (pkin(3) * t537 + pkin(7)) * t538;
t434 = t532 - t574;
t425 = qJD(1) * t565 - t656;
t424 = -pkin(7) * t596 + t632;
t404 = -qJ(6) * t477 + t630;
t403 = qJ(6) * t559 + t575;
t396 = -t435 - t658;
t391 = pkin(3) * t549 + pkin(7) * t620 + qJ(4) * t592 - t629;
t386 = qJD(2) * t565 - t563;
t384 = -qJD(4) * t541 + (-t538 * t611 - t591) * pkin(7) + t599;
t383 = qJD(2) * t555 + (qJD(3) - qJD(5)) * t538 * t559;
t382 = t422 * t538 + t536 * t593 - t539 * t594;
t369 = pkin(5) * t560 + qJD(6) + t387;
t366 = -qJ(6) * t450 + t639;
t365 = pkin(5) * t541 - qJ(6) * t451 + t577;
t350 = t354 - t677;
t346 = -qJ(6) * t382 - qJD(6) * t450 + t553;
t345 = -pkin(5) * t621 - qJ(6) * t383 - qJD(5) * t639 - qJD(6) * t451 + t579;
t344 = -qJ(6) * t368 - qJD(6) * t560 + t554;
t343 = -pkin(5) * t587 + qJ(6) * t367 - qJD(6) * t416 - t669;
t1 = [((-t473 * t540 - t475 * t537) * t620 + (t659 - t436 * t540 + (t473 * t537 - t475 * t540) * qJD(3)) * t538) * MDP(12) + ((-pkin(7) * t591 + t631) * t516 + t557 * t541 + (-pkin(7) * t435 - t497 * t619) * t538 + ((pkin(7) * t475 + t654) * t541 + (-pkin(7) * t653 - qJD(1) * t628 - t420) * t538) * qJD(2)) * MDP(17) + (-t435 * t649 + (-t592 + t593) * t475) * MDP(11) + (-t367 * t451 + t383 * t416) * MDP(22) + (t370 * t433 + t375 * t443 + t377 * t434 + t384 * t402 + t386 * t400 + t391 * t408) * MDP(21) + (t579 * t670 - t580 * t541 + t381 * t560 + t432 * t368 + t364 * t450 + t387 * t382 + (-t354 * t541 - t639 * t670) * qJD(5) + (-qJD(1) * t577 - t353) * t621) * MDP(27) + (-t670 - t623) * t583 + (-t368 * t541 - t382 * t670 + (qJD(1) * t450 + t560) * t621) * MDP(25) + (-t553 * t670 - t554 * t541 + t381 * t416 - t432 * t367 + t364 * t451 + t387 * t383 + (qJD(1) * t639 + t354) * t621) * MDP(28) + (-t367 * t541 + t383 * t670 + (-qJD(1) * t451 - t416) * t621) * MDP(24) + (-pkin(7) * t646 + t538 * t582) * MDP(9) + (-t563 * t516 - t568 * t541 + (pkin(7) * t436 + t497 * t618) * t538 + ((pkin(7) * t473 + t655) * t541 + (t574 * qJD(1) + t419 + (-t516 + t623) * t664) * t538) * qJD(2)) * MDP(16) - MDP(7) * t648 + (pkin(7) * t648 + t541 * t582) * MDP(10) + (-t516 - t623) * t584 + (-t384 * t516 - t391 * t475 + t435 * t443 + (-t408 * t611 - t370) * t541 + (t408 * t619 - t661 + (qJD(1) * t433 + t402) * qJD(2)) * t538) * MDP(20) + (t386 * t516 + t391 * t473 + t436 * t443 + (t408 * t622 + t377) * t541 + (t408 * t618 + t662 + (-qJD(1) * t434 - t400) * qJD(2)) * t538) * MDP(18) + (t367 * t450 - t368 * t451 - t382 * t416 - t383 * t560) * MDP(23) + (-t343 * t451 - t344 * t450 - t345 * t416 - t346 * t560 - t348 * t383 - t350 * t382 + t365 * t367 - t366 * t368) * MDP(29) + 0.2e1 * t586 * t679 + (t343 * t365 + t344 * t366 + t348 * t345 + t350 * t346 + (pkin(5) * t450 + t432) * t347 + (pkin(5) * t382 + t381) * t369) * MDP(30) + (t516 * t590 + t436 * t541 + (-t473 * t538 - t537 * t558) * qJD(2)) * MDP(14) + (t516 * t592 + t435 * t541 + (t475 * t538 + t540 * t558) * qJD(2)) * MDP(13) + (-t384 * t473 + t386 * t475 - t433 * t436 - t434 * t435 + t562 * t620 + (-t370 * t537 + t377 * t540 + (-t400 * t537 - t402 * t540) * qJD(3)) * t538) * MDP(19) + MDP(6) * t646 + t678 * t680; -t645 * t679 + t544 * t678 + (-t475 * t653 - t659) * MDP(11) + ((-t435 + t658) * t540 + (-t436 + t657) * t537) * MDP(12) - t516 * t618 * MDP(13) + t516 * t619 * MDP(14) + (t485 * t653 - pkin(2) * t436 + (t605 + t655) * qJD(3)) * MDP(16) + (pkin(2) * t435 - t456 * t516 + (-t606 + t654) * qJD(3)) * MDP(17) + (-t661 - t425 * t516 + t436 * t490 - t636 * t473 + (t408 * t537 + t605) * qJD(3)) * MDP(18) + (t424 * t473 - t425 * t475 + (t370 - t516 * t400 + (qJD(3) * t475 - t436) * pkin(8)) * t540 + ((qJD(3) * t473 - t435) * pkin(8) + t676) * t537) * MDP(19) + (-t662 + t424 * t516 + t435 * t490 + t636 * t475 + (-t408 * t540 + t606) * qJD(3)) * MDP(20) + (t375 * t490 - t400 * t425 - t402 * t424 - t636 * t408 + (qJD(3) * t562 + t370 * t540 + t377 * t537) * pkin(8)) * MDP(21) + (t367 * t559 - t416 * t638) * MDP(22) + (t367 * t477 + t368 * t559 - t416 * t637 + t560 * t638) * MDP(23) + (-t638 * t670 + (qJD(2) * t559 + t416) * t624) * MDP(24) + (-t637 * t670 + (qJD(2) * t477 - t560) * t624) * MDP(25) + (t364 * t477 + t467 * t368 + (-t500 * t615 + (-qJD(5) * t499 + t682) * t536 - t673) * t670 + t635 * t560 + t637 * t387 + (-qJD(2) * t575 + t353) * t624) * MDP(27) + (-t364 * t559 - t467 * t367 + t672 * t670 + t635 * t416 - t638 * t387 + (qJD(2) * t630 - t354) * t624) * MDP(28) + (t343 * t559 - t344 * t477 + t348 * t638 - t350 * t637 + t367 * t403 - t368 * t404 - t416 * t642 - t560 * t643) * MDP(29) + (t344 * t404 + t343 * t403 + t347 * (pkin(5) * t477 + t467) + (t637 * pkin(5) + (pkin(4) * t516 + t665) * t537 + t671) * t369 + t643 * t350 + t642 * t348) * MDP(30) + (t516 * t613 + t670 * t612 + (t516 * t647 + t538 * t571) * MDP(13) + (-t516 * t650 + t538 * t572) * MDP(14) + (-t419 * t538 + (-t497 * t541 - t604) * t537 + (t516 * t652 - t541 * t572) * pkin(7)) * MDP(16) + (-t497 * t647 + (t420 - t603) * t538 + (t516 * t649 + t541 * t571) * pkin(7)) * MDP(17) + (t400 * t538 + (-t408 * t541 - t604) * t537) * MDP(18) + (t408 * t647 + (-t402 + t603) * t538) * MDP(20)) * qJD(1) + (MDP(9) * t538 * t544 + MDP(10) * t645) * pkin(1); t475 * t473 * MDP(11) + (-t473 ^ 2 + t475 ^ 2) * MDP(12) + t396 * MDP(13) + (-t436 - t657) * MDP(14) + qJD(1) * t584 + (-t475 * t497 + t556) * MDP(16) + (t473 * t497 + t545) * MDP(17) + (-t408 * t475 - t426 * t473 + t556 + 0.2e1 * t570) * MDP(18) + (pkin(3) * t435 - qJ(4) * t436 + (t402 - t420) * t475 + (t400 - t610) * t473) * MDP(19) + (-t408 * t473 + t426 * t475 - 0.2e1 * t501 + 0.2e1 * t509 - t545) * MDP(20) + (-pkin(3) * t377 + qJ(4) * t370 - t400 * t420 + t402 * t610 - t408 * t426) * MDP(21) + t367 * MDP(24) + (t368 - t683) * MDP(25) + (-t581 * t587 + t670 * t674 + t669 + t684) * MDP(27) + (-t397 * t416 + t489 * t587 + t670 * t675 + t554) * MDP(28) + (t367 * t484 - t368 * t489 + (-t350 - t633) * t416) * MDP(29) + (t344 * t489 + t343 * t484 - t369 * (-pkin(5) * t416 + t397) + t634 * t350 + t633 * t348) * MDP(30) + (-t397 * MDP(27) + (t348 - t634) * MDP(29) + t681) * t560 - t689; -MDP(18) * t587 + t396 * MDP(19) - t516 ^ 2 * MDP(20) + t676 * MDP(21) + (MDP(18) * t473 - MDP(20) * t475 + MDP(21) * t408 - MDP(27) * t560 - MDP(28) * t416 - MDP(30) * t369) * t475 + (-MDP(27) * t587 + (-t560 * t670 + t367) * MDP(29) + (t350 * t670 + t343) * MDP(30) - MDP(28) * t573) * t539 + (MDP(28) * t587 + (t416 * t516 - t368 + t617) * MDP(29) + (-t348 * t670 + t344) * MDP(30) - MDP(27) * t573) * t536; t600 * MDP(24) + (-t576 + t683) * MDP(25) + (t354 * t670 - t580 - t684) * MDP(27) + (t353 * t670 - t601) * MDP(28) + t644 * MDP(30) * t350 + (t367 * MDP(29) + (-t369 * t416 + t343) * MDP(30)) * pkin(5) + (-MDP(29) * t644 - t681) * t560 + ((-MDP(25) * t475 - MDP(27) * t385) * t539 + (-MDP(24) * t475 - MDP(25) * t473 - MDP(27) * t380 + MDP(28) * t385) * t536) * qJD(5) + t689; (-t411 - t668) * MDP(29) + (t348 * t416 + t350 * t560 + t347) * MDP(30);];
tauc  = t1;
