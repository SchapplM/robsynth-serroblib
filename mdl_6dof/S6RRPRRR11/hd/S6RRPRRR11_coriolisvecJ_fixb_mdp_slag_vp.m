% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:33:06
% EndTime: 2019-03-09 14:33:24
% DurationCPUTime: 11.13s
% Computational Cost: add. (6084->543), mult. (13850->731), div. (0->0), fcn. (9493->8), ass. (0->245)
t592 = cos(qJ(2));
t593 = -pkin(2) - pkin(8);
t588 = sin(qJ(2));
t645 = -qJ(3) * t588 - pkin(1);
t525 = t592 * t593 + t645;
t492 = t525 * qJD(1);
t674 = qJD(1) * t588;
t572 = pkin(7) * t674;
t727 = qJD(3) + t572;
t660 = pkin(3) * t674 + t727;
t498 = qJD(2) * t593 + t660;
t587 = sin(qJ(4));
t591 = cos(qJ(4));
t441 = -t492 * t587 + t591 * t498;
t673 = qJD(1) * t592;
t648 = t587 * t673;
t670 = qJD(2) * t591;
t534 = -t648 + t670;
t430 = -pkin(9) * t534 + t441;
t568 = qJD(4) + t674;
t423 = pkin(4) * t568 + t430;
t442 = t492 * t591 + t498 * t587;
t672 = qJD(2) * t587;
t532 = t591 * t673 + t672;
t431 = -pkin(9) * t532 + t442;
t590 = cos(qJ(5));
t427 = t590 * t431;
t586 = sin(qJ(5));
t388 = t423 * t586 + t427;
t466 = t590 * t532 + t534 * t586;
t734 = pkin(10) * t466;
t382 = t388 - t734;
t585 = sin(qJ(6));
t662 = qJD(6) * t585;
t380 = t382 * t662;
t589 = cos(qJ(6));
t618 = t532 * t586 - t590 * t534;
t711 = t618 * t585;
t421 = -t589 * t466 + t711;
t573 = pkin(7) * t673;
t542 = pkin(3) * t673 + t573;
t582 = qJD(2) * qJ(3);
t518 = t582 + t542;
t475 = pkin(4) * t532 + t518;
t428 = pkin(5) * t466 + t475;
t742 = -t428 * t421 + t380;
t620 = t466 * t585 + t589 * t618;
t659 = qJD(1) * qJD(2);
t646 = t592 * t659;
t741 = MDP(33) * t646 + (-t421 ^ 2 + t620 ^ 2) * MDP(30) + t421 * MDP(29) * t620;
t647 = t588 * t659;
t482 = -qJD(4) * t532 + t587 * t647;
t666 = qJD(4) * t591;
t483 = qJD(2) * t666 - qJD(4) * t648 - t591 * t647;
t663 = qJD(5) * t590;
t664 = qJD(5) * t586;
t411 = t590 * t482 - t586 * t483 - t532 * t663 - t534 * t664;
t567 = pkin(2) * t647;
t623 = pkin(8) * t588 - qJ(3) * t592;
t668 = qJD(3) * t588;
t601 = qJD(2) * t623 - t668;
t477 = qJD(1) * t601 + t567;
t566 = pkin(7) * t646;
t524 = pkin(3) * t646 + t566;
t634 = -t477 * t587 + t591 * t524;
t600 = -qJD(4) * t442 + t634;
t392 = pkin(4) * t646 - pkin(9) * t482 + t600;
t654 = -t591 * t477 - t498 * t666 - t587 * t524;
t667 = qJD(4) * t587;
t605 = -t492 * t667 - t654;
t396 = -pkin(9) * t483 + t605;
t639 = t590 * t392 - t586 * t396;
t599 = -qJD(5) * t388 + t639;
t371 = pkin(5) * t646 - pkin(10) * t411 + t599;
t597 = qJD(5) * t618 - t482 * t586 - t590 * t483;
t628 = -t586 * t392 - t590 * t396 - t423 * t663 + t431 * t664;
t372 = pkin(10) * t597 - t628;
t740 = -t585 * t371 - t589 * t372 + t742;
t661 = qJD(6) * t589;
t655 = t589 * t411 - t466 * t661 + t585 * t597;
t377 = t618 * t662 + t655;
t559 = qJD(5) + t568;
t637 = t411 * t585 - t589 * t597;
t598 = qJD(6) * t620 - t637;
t550 = qJD(6) + t559;
t730 = t550 * t620;
t731 = t421 * t550;
t739 = MDP(26) * t646 + (-t466 ^ 2 + t618 ^ 2) * MDP(23) + (t466 * t559 + t411) * MDP(24) + (-t559 * t618 + t597) * MDP(25) - t466 * MDP(22) * t618 + (t598 - t730) * MDP(32) + (t377 - t731) * MDP(31) + t741;
t640 = t589 * t371 - t585 * t372;
t726 = t428 * t620 + t640;
t698 = t587 * t588;
t613 = pkin(4) * t592 - pkin(9) * t698;
t577 = pkin(2) * t674;
t503 = qJD(1) * t623 + t577;
t631 = -t503 * t587 + t591 * t542;
t713 = pkin(9) - t593;
t736 = qJD(1) * t613 - t713 * t667 + t631;
t547 = t713 * t591;
t652 = t591 * t674;
t682 = t591 * t503 + t587 * t542;
t735 = pkin(9) * t652 + qJD(4) * t547 + t682;
t535 = t586 * t591 + t587 * t590;
t606 = t535 * t588;
t716 = qJD(4) + qJD(5);
t684 = -qJD(1) * t606 - t716 * t535;
t714 = pkin(3) + pkin(7);
t733 = pkin(10) * t618;
t694 = t590 * t591;
t700 = t586 * t587;
t683 = -t586 * t667 - t587 * t664 + t590 * t652 - t674 * t700 + t694 * t716;
t725 = t466 * t475 + t628;
t724 = t475 * t618 + t599;
t721 = -0.2e1 * t659;
t583 = t588 ^ 2;
t584 = t592 ^ 2;
t720 = MDP(5) * (t583 - t584);
t552 = t714 * t588;
t538 = t591 * t552;
t644 = pkin(9) * t592 - t525;
t453 = pkin(4) * t588 + t587 * t644 + t538;
t537 = t587 * t552;
t681 = t591 * t525 + t537;
t693 = t591 * t592;
t457 = -pkin(9) * t693 + t681;
t686 = t586 * t453 + t590 * t457;
t719 = t684 * t589;
t718 = t736 * t590;
t617 = -t694 + t700;
t471 = -t535 * t585 - t589 * t617;
t546 = t713 * t587;
t680 = -t590 * t546 - t586 * t547;
t653 = -pkin(4) * t591 - pkin(3);
t679 = pkin(4) * t666 - t653 * t674 + t727;
t717 = -t546 * t664 + t547 * t663 + t586 * t736 + t735 * t590;
t715 = t592 * t716;
t712 = qJD(2) * pkin(2);
t710 = t482 * t591;
t671 = qJD(2) * t588;
t541 = t714 * t671;
t581 = qJD(2) * qJD(3);
t500 = -qJD(1) * t541 + t581;
t709 = t500 * t587;
t708 = t500 * t591;
t707 = t532 * t568;
t706 = t534 * t568;
t705 = t534 * t592;
t703 = t568 * t588;
t702 = t568 * t593;
t701 = t585 * t586;
t425 = t586 * t431;
t699 = t586 * t589;
t594 = qJD(2) ^ 2;
t697 = t588 * t594;
t387 = t590 * t423 - t425;
t381 = t387 + t733;
t379 = pkin(5) * t559 + t381;
t696 = t589 * t379;
t695 = t589 * t382;
t692 = t592 * t594;
t595 = qJD(1) ^ 2;
t691 = t592 * t595;
t570 = t587 * pkin(4) + qJ(3);
t470 = t589 * t535 - t585 * t617;
t690 = -qJD(6) * t470 - t585 * t683 + t719;
t689 = qJD(6) * t471 + t585 * t684 + t589 * t683;
t688 = t590 * t430 - t425;
t685 = pkin(5) * t683 + t679;
t553 = t714 * t592;
t548 = -pkin(2) * t592 + t645;
t519 = qJD(1) * t548;
t669 = qJD(2) * t592;
t665 = qJD(4) * t592;
t658 = pkin(4) * qJD(5) * t550;
t657 = t591 * t703;
t656 = t588 * t691;
t514 = pkin(4) * t693 + t553;
t651 = t568 * t666;
t650 = t587 * t665;
t649 = t591 * t665;
t643 = t683 * t559;
t642 = pkin(1) * t721;
t641 = qJD(3) - t712;
t576 = pkin(2) * t671;
t487 = t576 + t601;
t543 = t714 * t669;
t632 = -t487 * t587 + t591 * t543;
t403 = t613 * qJD(2) + (t591 * t644 - t537) * qJD(4) + t632;
t604 = t591 * t487 - t525 * t667 + t587 * t543 + t552 * t666;
t410 = (t588 * t670 + t650) * pkin(9) + t604;
t638 = t590 * t403 - t410 * t586;
t636 = -t430 * t586 - t427;
t635 = t590 * t453 - t457 * t586;
t630 = t546 * t586 - t590 * t547;
t629 = qJD(6) * t379 + t372;
t555 = t588 * t646;
t627 = t559 * t684 - t617 * t646;
t451 = -pkin(10) * t535 + t680;
t626 = pkin(5) * t673 + pkin(10) * t684 + qJD(5) * t680 + qJD(6) * t451 - t586 * t735 + t718;
t450 = pkin(10) * t617 + t630;
t625 = pkin(10) * t683 - qJD(6) * t450 + t717;
t624 = qJD(6) * t617 - t683;
t374 = t585 * t379 + t695;
t505 = t535 * t592;
t397 = pkin(5) * t588 + pkin(10) * t505 + t635;
t504 = t617 * t592;
t398 = pkin(10) * t504 + t686;
t622 = t397 * t585 + t398 * t589;
t619 = t589 * t504 + t505 * t585;
t449 = t504 * t585 - t505 * t589;
t616 = -qJD(1) * t584 + t703;
t615 = -0.2e1 * qJD(2) * t519;
t614 = t568 * t587;
t571 = pkin(4) * t590 + pkin(5);
t612 = pkin(4) * t699 + t571 * t585;
t611 = -pkin(4) * t701 + t571 * t589;
t607 = -qJ(3) * t669 - t668;
t489 = qJD(1) * t607 + t567;
t508 = t576 + t607;
t610 = pkin(7) * t594 + qJD(1) * t508 + t489;
t608 = t518 * t588 + t593 * t669;
t603 = t586 * t403 + t590 * t410 + t453 * t663 - t457 * t664;
t447 = pkin(4) * t483 + t500;
t476 = -pkin(4) * t650 + (-pkin(7) + t653) * t671;
t544 = pkin(7) * t647 - t581;
t545 = t572 + t641;
t551 = -t573 - t582;
t596 = -t544 * t592 + (t545 * t592 + (t551 + t573) * t588) * qJD(2);
t557 = t591 * t646;
t539 = -qJ(3) * t673 + t577;
t497 = t519 * t674;
t495 = pkin(5) * t535 + t570;
t462 = -pkin(5) * t504 + t514;
t437 = pkin(4) * t534 - pkin(5) * t618;
t433 = t535 * t715 - t617 * t671;
t432 = qJD(2) * t606 + t617 * t715;
t413 = -pkin(5) * t433 + t476;
t389 = -pkin(5) * t597 + t447;
t386 = t688 + t733;
t385 = t636 + t734;
t384 = qJD(6) * t449 + t432 * t585 - t589 * t433;
t383 = qJD(6) * t619 + t432 * t589 + t433 * t585;
t376 = pkin(10) * t433 + t603;
t375 = pkin(5) * t669 - pkin(10) * t432 - qJD(5) * t686 + t638;
t373 = -t382 * t585 + t696;
t1 = [(t597 * t588 + t433 * t559 + (qJD(1) * t504 - t466) * t669) * MDP(25) + (t411 * t504 - t432 * t466 - t433 * t618 - t505 * t597) * MDP(23) + (t638 * t559 + t639 * t588 + t476 * t466 - t514 * t597 - t447 * t504 - t475 * t433 + (-t388 * t588 - t559 * t686) * qJD(5) + (qJD(1) * t635 + t387) * t669) * MDP(27) + (t377 * t619 + t383 * t421 + t384 * t620 + t449 * t598) * MDP(30) + ((t375 * t589 - t376 * t585) * t550 + t640 * t588 - t413 * t421 - t462 * t598 - t389 * t619 + t428 * t384 + (-t374 * t588 - t550 * t622) * qJD(6) + ((t397 * t589 - t398 * t585) * qJD(1) + t373) * t669) * MDP(34) + (t598 * t588 - t384 * t550 + (qJD(1) * t619 + t421) * t669) * MDP(32) + (-t588 * t610 + t592 * t615) * MDP(13) + (t588 * t615 + t592 * t610) * MDP(12) - MDP(7) * t697 + (pkin(7) * t697 + t592 * t642) * MDP(10) + (t568 * t650 - t483 * t588 + (-t532 * t592 + t591 * t616) * qJD(2)) * MDP(18) + (t568 * t669 + t555) * MDP(19) + (t559 * t669 + t555) * MDP(26) + (t550 * t669 + t555) * MDP(33) + (-t482 * t587 * t592 + (t587 * t671 - t649) * t534) * MDP(15) + (-t568 * t649 + t482 * t588 + (t587 * t616 + t705) * qJD(2)) * MDP(17) + ((-t532 * t587 + t534 * t591) * t671 + (-t710 + t483 * t587 + (t532 * t591 + t534 * t587) * qJD(4)) * t592) * MDP(16) + t720 * t721 + MDP(6) * t692 + (t632 * t568 - t541 * t532 + t553 * t483 + (-t518 * t670 + t634) * t588 + (-t442 * t588 - t568 * t681) * qJD(4) + (-t518 * t667 + t708 + ((-t525 * t587 + t538) * qJD(1) + t441) * qJD(2)) * t592) * MDP(20) + 0.2e1 * MDP(4) * t555 + (-pkin(7) * t692 + t588 * t642) * MDP(9) + (t377 * t588 + t383 * t550 + (qJD(1) * t449 - t620) * t669) * MDP(31) + (t462 * t377 + t380 * t588 + t428 * t383 + t389 * t449 - t413 * t620 + (-(-qJD(6) * t398 + t375) * t550 - t371 * t588) * t585 + (-(qJD(6) * t397 + t376) * t550 - t629 * t588) * t589 + (-qJD(1) * t622 - t374) * t669) * MDP(35) + (t377 * t449 - t383 * t620) * MDP(29) + (-t411 * t505 - t432 * t618) * MDP(22) + (t411 * t588 + t432 * t559 + (-qJD(1) * t505 - t618) * t669) * MDP(24) + (-t603 * t559 + t628 * t588 - t476 * t618 + t514 * t411 - t447 * t505 + t475 * t432 + (-qJD(1) * t686 - t388) * t669) * MDP(28) + (pkin(7) * t596 + t489 * t548 + t508 * t519) * MDP(14) + t596 * MDP(11) + (-t604 * t568 - t541 * t534 + t553 * t482 + ((qJD(2) * t518 + qJD(4) * t492) * t587 + t654) * t588 + (-t518 * t666 - t709 + (-qJD(1) * t681 - t442) * qJD(2)) * t592) * MDP(21); (-t411 * t535 - t466 * t684 - t597 * t617 + t618 * t683) * MDP(23) + ((t585 * t626 + t589 * t625) * MDP(35) + (t585 * t625 - t589 * t626) * MDP(34) - t689 * MDP(32) + t690 * MDP(31)) * t550 + (t389 * t470 - t421 * t685 + t689 * t428 - t495 * t598) * MDP(34) + (-t377 * t470 + t421 * t690 + t471 * t598 + t620 * t689) * MDP(30) + ((-qJD(2) * t470 - t421) * MDP(32) + (-(t450 * t585 + t451 * t589) * qJD(2) + t374) * MDP(35) + (qJD(2) * t471 + t620) * MDP(31) + (-qJD(2) * t680 + t388) * MDP(28) + ((t450 * t589 - t451 * t585) * qJD(2) - t373) * MDP(34) + (-qJD(2) * t535 + t466) * MDP(25) - t539 * MDP(12) + t618 * MDP(24) + (qJD(2) * t630 - t387) * MDP(27) - t568 * MDP(19) - t559 * MDP(26) - t550 * MDP(33)) * t673 + (-t570 * t597 + t447 * t535 + (t546 * t663 + (qJD(5) * t547 + t735) * t586 - t718) * t559 + t683 * t475 + t679 * t466) * MDP(27) + t497 * MDP(12) + (-t568 * t667 + t557 + (-t568 * t698 - t705) * qJD(1)) * MDP(17) + ((-t483 - t706) * t591 + (-t482 + t707) * t587) * MDP(16) + (qJ(3) * t482 + t708 + t682 * t568 + t660 * t534 + (-t518 * t587 - t591 * t702) * qJD(4) + (t442 * t592 - t587 * t608) * qJD(1)) * MDP(21) + (qJ(3) * t483 + t709 - t631 * t568 + t660 * t532 + (t518 * t591 - t587 * t702) * qJD(4) + (-t441 * t592 + t591 * t608) * qJD(1)) * MDP(20) + t595 * t720 + t627 * MDP(24) + (0.2e1 * t581 + (t519 * t592 + t539 * t588) * qJD(1)) * MDP(13) + (t377 * t471 - t620 * t690) * MDP(29) + (-t411 * t617 - t618 * t684) * MDP(22) + (t570 * t411 - t447 * t617 + t684 * t475 + t717 * t559 - t618 * t679) * MDP(28) - MDP(4) * t656 + (-qJ(3) * t544 - qJD(3) * t551 - t519 * t539 + (-t551 * t588 + (-t545 - t712) * t592) * qJD(1) * pkin(7)) * MDP(14) + (t495 * t377 + t389 * t471 + t690 * t428 - t620 * t685) * MDP(35) + (MDP(9) * t588 * t595 + MDP(10) * t691) * pkin(1) + (-t651 + (-t657 + (t532 - t672) * t592) * qJD(1)) * MDP(18) + ((-t551 - t582) * t588 + (-t545 + t641) * t592) * qJD(1) * MDP(11) - t643 * MDP(25) + (-t534 * t614 + t710) * MDP(15); MDP(12) * t656 + (-t583 * t595 - t594) * MDP(13) + (qJD(2) * t551 + t497 + t566) * MDP(14) + (-qJD(2) * t532 - t568 * t614 + t557) * MDP(20) + (-t651 - qJD(2) * t534 + (-t587 * t669 - t657) * qJD(1)) * MDP(21) + (-qJD(2) * t466 + t627) * MDP(27) + (-t643 + (-t535 * t673 + t618) * qJD(2)) * MDP(28) + ((-t535 * t661 + t585 * t624 + t719) * t550 + (t471 * t673 + t421) * qJD(2)) * MDP(34) + ((t624 * t589 + (qJD(6) * t535 - t684) * t585) * t550 + (-t470 * t673 + t620) * qJD(2)) * MDP(35); MDP(19) * t646 + (-t532 ^ 2 + t534 ^ 2) * MDP(16) + (t611 * t646 - (t385 * t589 - t386 * t585) * t550 + t437 * t421 + (-t585 * t590 - t699) * t658 + (-t550 * t612 - t374) * qJD(6) + t726) * MDP(34) + (-t483 + t706) * MDP(18) + (t482 + t707) * MDP(17) + (-t636 * t559 + (-t466 * t534 - t559 * t664 + t590 * t646) * pkin(4) + t724) * MDP(27) + (t688 * t559 + (t534 * t618 - t559 * t663 - t586 * t646) * pkin(4) + t725) * MDP(28) + t534 * t532 * MDP(15) + (-t612 * t646 + (t385 * t585 + t386 * t589) * t550 + t437 * t620 - (t589 * t590 - t701) * t658 + (-t550 * t611 - t696) * qJD(6) + t740) * MDP(35) + (t441 * t568 + t518 * t532 - t605) * MDP(21) + (t442 * t568 - t518 * t534 + t600) * MDP(20) + t739; (t388 * t559 + t724) * MDP(27) + (t387 * t559 + t725) * MDP(28) + (-(-t381 * t585 - t695) * t550 - t374 * qJD(6) + (-t421 * t618 - t550 * t662 + t589 * t646) * pkin(5) + t726) * MDP(34) + ((-t382 * t550 - t371) * t585 + (t381 * t550 - t629) * t589 + (-t550 * t661 - t585 * t646 - t618 * t620) * pkin(5) + t742) * MDP(35) + t739; (t655 - t731) * MDP(31) + (-t637 - t730) * MDP(32) + (t374 * t550 + t726) * MDP(34) + (t373 * t550 + t740) * MDP(35) + (MDP(31) * t711 + MDP(32) * t620 - MDP(34) * t374 - MDP(35) * t696) * qJD(6) + t741;];
tauc  = t1;
