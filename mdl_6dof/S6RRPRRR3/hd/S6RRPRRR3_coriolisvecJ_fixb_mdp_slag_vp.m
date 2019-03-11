% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:25:35
% EndTime: 2019-03-09 13:25:52
% DurationCPUTime: 11.11s
% Computational Cost: add. (9620->538), mult. (24354->712), div. (0->0), fcn. (19121->10), ass. (0->239)
t590 = sin(pkin(11));
t591 = cos(pkin(11));
t595 = sin(qJ(2));
t599 = cos(qJ(2));
t568 = t590 * t599 + t591 * t595;
t557 = t568 * qJD(1);
t594 = sin(qJ(4));
t598 = cos(qJ(4));
t647 = t598 * qJD(2);
t520 = t557 * t594 - t647;
t522 = qJD(2) * t594 + t557 * t598;
t593 = sin(qJ(5));
t597 = cos(qJ(5));
t455 = t597 * t520 + t522 * t593;
t596 = cos(qJ(6));
t592 = sin(qJ(6));
t614 = t520 * t593 - t597 * t522;
t695 = t614 * t592;
t421 = -t596 * t455 + t695;
t556 = t568 * qJD(2);
t542 = qJD(1) * t556;
t616 = t455 * t592 + t596 * t614;
t732 = t542 * MDP(31) + (-t421 ^ 2 + t616 ^ 2) * MDP(28) + t421 * MDP(27) * t616;
t654 = qJD(1) * t595;
t679 = t591 * t599;
t555 = qJD(1) * t679 - t590 * t654;
t706 = qJD(4) + qJD(5);
t731 = t555 - t706;
t646 = qJD(1) * qJD(2);
t638 = t599 * t646;
t639 = t595 * t646;
t543 = -t590 * t639 + t591 * t638;
t653 = qJD(4) * t594;
t472 = qJD(4) * t647 + t598 * t543 - t557 * t653;
t473 = qJD(4) * t522 + t543 * t594;
t650 = qJD(5) * t597;
t651 = qJD(5) * t593;
t400 = t597 * t472 - t593 * t473 - t520 * t650 - t522 * t651;
t602 = qJD(5) * t614 - t472 * t593 - t597 * t473;
t648 = qJD(6) * t596;
t643 = t596 * t400 - t455 * t648 + t592 * t602;
t649 = qJD(6) * t592;
t369 = t614 * t649 + t643;
t549 = qJD(4) - t555;
t541 = qJD(5) + t549;
t631 = t400 * t592 - t596 * t602;
t603 = qJD(6) * t616 - t631;
t537 = qJD(6) + t541;
t721 = t537 * t616;
t722 = t421 * t537;
t730 = t542 * MDP(24) + (-t455 ^ 2 + t614 ^ 2) * MDP(21) + (t455 * t541 + t400) * MDP(22) + (-t541 * t614 + t602) * MDP(23) - t455 * MDP(20) * t614 + (t603 - t721) * MDP(30) + (t369 - t722) * MDP(29) + t732;
t642 = -pkin(2) * t599 - pkin(1);
t620 = t642 * qJD(1);
t574 = qJD(3) + t620;
t478 = -pkin(3) * t555 - pkin(8) * t557 + t574;
t702 = -qJ(3) - pkin(7);
t576 = t702 * t595;
t572 = qJD(1) * t576;
t701 = qJD(2) * pkin(2);
t564 = t572 + t701;
t577 = t702 * t599;
t573 = qJD(1) * t577;
t680 = t591 * t573;
t504 = t590 * t564 - t680;
t498 = qJD(2) * pkin(8) + t504;
t440 = t598 * t478 - t498 * t594;
t424 = -pkin(9) * t522 + t440;
t404 = pkin(4) * t549 + t424;
t441 = t478 * t594 + t498 * t598;
t425 = -pkin(9) * t520 + t441;
t415 = t597 * t425;
t386 = t404 * t593 + t415;
t724 = pkin(10) * t455;
t376 = t386 - t724;
t374 = t376 * t649;
t560 = t590 * t573;
t503 = t564 * t591 + t560;
t497 = -qJD(2) * pkin(3) - t503;
t449 = pkin(4) * t520 + t497;
t408 = pkin(5) * t455 + t449;
t717 = -t408 * t421 + t374;
t570 = t593 * t594 - t597 * t598;
t661 = t731 * t570;
t676 = t593 * t598;
t571 = t594 * t597 + t676;
t660 = t731 * t571;
t580 = pkin(2) * t639;
t476 = pkin(3) * t542 - pkin(8) * t543 + t580;
t466 = t598 * t476;
t637 = qJD(2) * t702;
t553 = qJD(3) * t599 + t595 * t637;
t527 = t553 * qJD(1);
t554 = -qJD(3) * t595 + t599 * t637;
t528 = t554 * qJD(1);
t469 = t527 * t591 + t528 * t590;
t605 = -qJD(4) * t441 - t469 * t594 + t466;
t383 = pkin(4) * t542 - pkin(9) * t472 + t605;
t652 = qJD(4) * t598;
t610 = t598 * t469 + t594 * t476 + t478 * t652 - t498 * t653;
t389 = -pkin(9) * t473 + t610;
t633 = t597 * t383 - t593 * t389;
t604 = -t386 * qJD(5) + t633;
t363 = pkin(5) * t542 - pkin(10) * t400 + t604;
t624 = -t593 * t383 - t597 * t389 - t404 * t650 + t425 * t651;
t364 = pkin(10) * t602 - t624;
t634 = t596 * t363 - t592 * t364;
t716 = t408 * t616 + t634;
t490 = pkin(2) * t654 + pkin(3) * t557 - pkin(8) * t555;
t482 = t598 * t490;
t509 = t572 * t591 + t560;
t582 = pkin(2) * t590 + pkin(8);
t703 = pkin(9) + t582;
t636 = qJD(4) * t703;
t727 = pkin(4) * t557 - t509 * t594 + t482 + (-pkin(9) * t555 + t636) * t598;
t663 = t594 * t490 + t598 * t509;
t685 = t555 * t594;
t726 = -pkin(9) * t685 + t594 * t636 + t663;
t725 = t653 - t685;
t723 = pkin(10) * t614;
t510 = -t570 * t592 + t571 * t596;
t666 = qJD(6) * t510 + t661 * t592 - t660 * t596;
t640 = t568 * t652;
t567 = t590 * t595 - t679;
t559 = t567 * qJD(2);
t674 = t594 * t559;
t718 = t640 - t674;
t715 = t449 * t455 + t624;
t714 = t449 * t614 + t604;
t711 = -0.2e1 * t646;
t710 = MDP(4) * t595;
t709 = MDP(5) * (t595 ^ 2 - t599 ^ 2);
t493 = t571 * t568;
t502 = pkin(3) * t567 - pkin(8) * t568 + t642;
t496 = t598 * t502;
t519 = t576 * t590 - t577 * t591;
t682 = t568 * t598;
t434 = pkin(4) * t567 - pkin(9) * t682 - t519 * t594 + t496;
t512 = t598 * t519;
t662 = t594 * t502 + t512;
t683 = t568 * t594;
t443 = -pkin(9) * t683 + t662;
t665 = t593 * t434 + t597 * t443;
t507 = t572 * t590 - t680;
t623 = pkin(4) * t725 - t507;
t708 = t727 * t597;
t565 = t703 * t594;
t566 = t703 * t598;
t658 = -t593 * t565 + t597 * t566;
t707 = t565 * t650 + t566 * t651 + t593 * t727 + t726 * t597;
t508 = t596 * t570 + t571 * t592;
t667 = -qJD(6) * t508 + t592 * t660 + t596 * t661;
t705 = -t510 * t542 - t537 * t667;
t704 = -t541 * t661 - t542 * t571;
t413 = t593 * t425;
t385 = t597 * t404 - t413;
t375 = t385 + t723;
t373 = pkin(5) * t541 + t375;
t700 = t373 * t596;
t699 = t421 * t557;
t698 = t616 * t557;
t697 = t455 * t557;
t696 = t614 * t557;
t694 = t472 * t594;
t692 = t520 * t549;
t691 = t520 * t557;
t690 = t522 * t549;
t689 = t522 * t557;
t687 = t542 * t593;
t684 = t559 * t598;
t586 = pkin(4) * t597 + pkin(5);
t681 = t586 * t542;
t678 = t592 * t542;
t677 = t593 * t596;
t675 = t594 * t542;
t600 = qJD(2) ^ 2;
t673 = t595 * t600;
t672 = t596 * t376;
t525 = t598 * t542;
t671 = t599 * t600;
t601 = qJD(1) ^ 2;
t670 = t599 * t601;
t669 = t597 * t424 - t413;
t664 = -pkin(5) * t660 + t623;
t645 = t595 * t701;
t583 = -pkin(2) * t591 - pkin(3);
t641 = t568 * t653;
t635 = pkin(1) * t711;
t491 = pkin(3) * t556 + pkin(8) * t559 + t645;
t483 = t598 * t491;
t489 = t553 * t591 + t554 * t590;
t393 = pkin(9) * t684 + pkin(4) * t556 - t489 * t594 + t483 + (-t512 + (pkin(9) * t568 - t502) * t594) * qJD(4);
t609 = t598 * t489 + t594 * t491 + t502 * t652 - t519 * t653;
t396 = -pkin(9) * t718 + t609;
t632 = t597 * t393 - t396 * t593;
t630 = -t424 * t593 - t415;
t629 = t597 * t434 - t443 * t593;
t468 = t527 * t590 - t591 * t528;
t488 = t553 * t590 - t591 * t554;
t627 = -t597 * t565 - t566 * t593;
t518 = -t591 * t576 - t577 * t590;
t626 = t549 * t598;
t625 = qJD(6) * t373 + t364;
t622 = -t508 * t542 - t537 * t666;
t621 = t541 * t660 - t570 * t542;
t486 = pkin(4) * t683 + t518;
t471 = -pkin(10) * t570 + t658;
t619 = pkin(5) * t557 + pkin(10) * t661 + qJD(5) * t658 + qJD(6) * t471 - t593 * t726 + t708;
t470 = -pkin(10) * t571 + t627;
t618 = -pkin(10) * t660 - qJD(6) * t470 + t707;
t368 = t592 * t373 + t672;
t615 = t468 * t568 - t519 * t542;
t494 = t570 * t568;
t444 = t596 * t493 - t494 * t592;
t445 = -t493 * t592 - t494 * t596;
t575 = -pkin(4) * t598 + t583;
t448 = pkin(4) * t718 + t488;
t613 = -t549 * t725 + t525;
t435 = pkin(4) * t473 + t468;
t611 = -t641 - t684;
t608 = t593 * t393 + t597 * t396 + t434 * t650 - t443 * t651;
t606 = t497 * t549 - t582 * t542;
t524 = pkin(5) * t570 + t575;
t506 = t542 * t567;
t447 = pkin(5) * t493 + t486;
t446 = pkin(4) * t522 - pkin(5) * t614;
t407 = -t559 * t676 - t593 * t641 - t651 * t683 + (t682 * t706 - t674) * t597;
t406 = -t493 * t706 + t570 * t559;
t394 = pkin(5) * t407 + t448;
t390 = -pkin(10) * t493 + t665;
t388 = pkin(5) * t567 + pkin(10) * t494 + t629;
t384 = -pkin(5) * t602 + t435;
t378 = t669 + t723;
t377 = t630 + t724;
t372 = qJD(6) * t445 + t406 * t592 + t596 * t407;
t371 = -qJD(6) * t444 + t406 * t596 - t407 * t592;
t367 = -t376 * t592 + t700;
t366 = -pkin(10) * t407 + t608;
t365 = pkin(5) * t556 - pkin(10) * t406 - qJD(5) * t665 + t632;
t1 = [(-t469 * t567 + t488 * t557 + t489 * t555 + t503 * t559 - t504 * t556 + t518 * t543 + t615) * MDP(11) + ((-t519 * t652 + t483) * t549 + t496 * t542 + (-t498 * t652 + t466) * t567 + t440 * t556 + t488 * t520 + t518 * t473 + t497 * t640 + ((-qJD(4) * t502 - t489) * t549 + (-qJD(4) * t478 - t469) * t567 - t497 * t559 + t615) * t594) * MDP(18) + (-(-t520 * t598 - t522 * t594) * t559 + (-t694 - t473 * t598 + (t520 * t594 - t522 * t598) * qJD(4)) * t568) * MDP(14) + (t400 * t567 + t406 * t541 - t494 * t542 - t556 * t614) * MDP(22) + (-t400 * t494 - t406 * t614) * MDP(20) + (-t386 * t556 + t486 * t400 + t449 * t406 - t435 * t494 - t448 * t614 - t541 * t608 - t542 * t665 + t567 * t624) * MDP(26) + (-t368 * t556 + t447 * t369 + t408 * t371 + t374 * t567 + t384 * t445 - t394 * t616 + (-(-qJD(6) * t390 + t365) * t537 - t388 * t542 - t363 * t567) * t592 + (-(qJD(6) * t388 + t366) * t537 - t390 * t542 - t625 * t567) * t596) * MDP(33) + (t369 * t567 + t371 * t537 + t445 * t542 - t556 * t616) * MDP(29) + (t369 * t445 - t371 * t616) * MDP(27) + 0.2e1 * t638 * t710 + (-t400 * t493 - t406 * t455 + t407 * t614 - t494 * t602) * MDP(21) + (-t407 * t541 - t455 * t556 - t493 * t542 + t567 * t602) * MDP(23) + (t448 * t455 - t486 * t602 + t435 * t493 + t449 * t407 + t632 * t541 + t629 * t542 + t633 * t567 + t385 * t556 + (-t386 * t567 - t541 * t665) * qJD(5)) * MDP(25) + (-t369 * t444 + t371 * t421 + t372 * t616 + t445 * t603) * MDP(28) + ((t365 * t596 - t366 * t592) * t537 + (t388 * t596 - t390 * t592) * t542 + t634 * t567 + t367 * t556 - t394 * t421 - t447 * t603 + t384 * t444 + t408 * t372 + ((-t388 * t592 - t390 * t596) * t537 - t368 * t567) * qJD(6)) * MDP(32) + (-t372 * t537 + t421 * t556 - t444 * t542 + t567 * t603) * MDP(30) + (-t473 * t567 - t520 * t556 - t549 * t718 - t568 * t675) * MDP(16) + MDP(6) * t671 + (t549 * t556 + t506) * MDP(17) + (t541 * t556 + t506) * MDP(24) + (t537 * t556 + t506) * MDP(31) + (t468 * t518 + t469 * t519 - t503 * t488 + t504 * t489 + (t574 + t620) * t645) * MDP(12) + t709 * t711 + (-pkin(7) * t671 + t595 * t635) * MDP(9) + (t472 * t567 + t522 * t556 + t525 * t568 + t549 * t611) * MDP(15) - MDP(7) * t673 + (pkin(7) * t673 + t599 * t635) * MDP(10) + (t472 * t682 + t522 * t611) * MDP(13) + (-t441 * t556 + t468 * t682 + t518 * t472 + t488 * t522 + t497 * t611 - t542 * t662 - t549 * t609 - t567 * t610) * MDP(19); ((t503 - t509) * t555 + (-t542 * t590 - t543 * t591) * pkin(2)) * MDP(11) + (-t468 * t598 + t583 * t473 - t507 * t520 + (-t582 * t652 - t482) * t549 + (t509 * t549 + t606) * t594) * MDP(18) + (-t575 * t602 + t435 * t570 + t627 * t542 + (-t566 * t650 + (qJD(5) * t565 + t726) * t593 - t708) * t541 + t623 * t455 - t660 * t449) * MDP(25) + (t468 * t594 + t583 * t472 - t507 * t522 + (t582 * t653 + t663) * t549 + t606 * t598) * MDP(19) + ((t470 * t596 - t471 * t592) * t542 - t524 * t603 + t384 * t508 + (t592 * t618 - t596 * t619) * t537 - t664 * t421 + t666 * t408) * MDP(32) + (-(t470 * t592 + t471 * t596) * t542 + t524 * t369 + t384 * t510 + (t592 * t619 + t596 * t618) * t537 - t664 * t616 + t667 * t408) * MDP(33) + (t575 * t400 + t435 * t571 + t661 * t449 + t541 * t707 - t658 * t542 - t614 * t623) * MDP(26) + t601 * t709 - t670 * t710 + (t503 * t507 - t504 * t509 + (-t468 * t591 + t469 * t590 - t574 * t654) * pkin(2)) * MDP(12) + (t400 * t571 - t614 * t661) * MDP(20) + (-t400 * t570 - t455 * t661 + t571 * t602 - t614 * t660) * MDP(21) + (t369 * t510 - t616 * t667) * MDP(27) + (-t369 * t508 + t421 * t667 + t510 * t603 + t616 * t666) * MDP(28) + (t549 * t626 + t675 - t689) * MDP(15) + (t613 + t691) * MDP(16) + ((t472 - t692) * t598 + (-t473 - t690) * t594) * MDP(14) + (t522 * t626 + t694) * MDP(13) + (t696 - t704) * MDP(22) + (t621 + t697) * MDP(23) + (t698 - t705) * MDP(29) + (t622 - t699) * MDP(30) + (MDP(9) * t595 * t601 + MDP(10) * t670) * pkin(1) - ((-t504 + t507) * MDP(11) + t440 * MDP(18) + t385 * MDP(25) - t441 * MDP(19) + t367 * MDP(32) - t368 * MDP(33) - t386 * MDP(26) + t537 * MDP(31) + t541 * MDP(24) + t549 * MDP(17)) * t557; (-t555 ^ 2 - t557 ^ 2) * MDP(11) + (t503 * t557 - t504 * t555 + t580) * MDP(12) + (t613 - t691) * MDP(18) + (-t549 ^ 2 * t598 - t675 - t689) * MDP(19) + (t621 - t697) * MDP(25) + (t696 + t704) * MDP(26) + (t622 + t699) * MDP(32) + (t698 + t705) * MDP(33); (t669 * t541 + (t522 * t614 - t541 * t650 - t687) * pkin(4) + t715) * MDP(26) + (t446 * t616 + (-t681 - t363 + (t377 - (-qJD(5) - qJD(6)) * t593 * pkin(4)) * t537) * t592 + (-pkin(4) * t687 + (-pkin(4) * t650 - qJD(6) * t586 + t378) * t537 - t625) * t596 + t717) * MDP(33) + (-t630 * t541 + (-t455 * t522 - t541 * t651 + t542 * t597) * pkin(4) + t714) * MDP(25) + (-t473 + t690) * MDP(16) + (t441 * t549 - t497 * t522 + t605) * MDP(18) + t542 * MDP(17) + (-t520 ^ 2 + t522 ^ 2) * MDP(14) + (t440 * t549 + t497 * t520 - t610) * MDP(19) + (t596 * t681 - (t377 * t596 - t378 * t592) * t537 + t446 * t421 + (-t593 * t678 + (-t592 * t597 - t677) * t537 * qJD(5)) * pkin(4) + ((-pkin(4) * t677 - t586 * t592) * t537 - t368) * qJD(6) + t716) * MDP(32) + (t472 + t692) * MDP(15) + t522 * t520 * MDP(13) + t730; (t386 * t541 + t714) * MDP(25) + (t385 * t541 + t715) * MDP(26) + (-(-t375 * t592 - t672) * t537 - t368 * qJD(6) + (-t421 * t614 - t537 * t649 + t542 * t596) * pkin(5) + t716) * MDP(32) + ((-t376 * t537 - t363) * t592 + (t375 * t537 - t625) * t596 + (-t537 * t648 - t614 * t616 - t678) * pkin(5) + t717) * MDP(33) + t730; (t643 - t722) * MDP(29) + (-t631 - t721) * MDP(30) + (t368 * t537 + t716) * MDP(32) + (-t592 * t363 - t596 * t364 + t367 * t537 + t717) * MDP(33) + (MDP(29) * t695 + MDP(30) * t616 - MDP(32) * t368 - MDP(33) * t700) * qJD(6) + t732;];
tauc  = t1;
