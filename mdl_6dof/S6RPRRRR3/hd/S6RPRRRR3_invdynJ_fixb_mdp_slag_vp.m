% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:02:49
% EndTime: 2019-03-09 07:03:07
% DurationCPUTime: 12.84s
% Computational Cost: add. (7259->594), mult. (15439->788), div. (0->0), fcn. (11177->18), ass. (0->251)
t643 = sin(qJ(4));
t648 = cos(qJ(4));
t717 = t648 * qJD(3);
t644 = sin(qJ(3));
t731 = qJD(1) * t644;
t578 = t643 * t731 - t717;
t727 = qJD(3) * t643;
t580 = t648 * t731 + t727;
t642 = sin(qJ(5));
t647 = cos(qJ(5));
t513 = t647 * t578 + t580 * t642;
t646 = cos(qJ(6));
t641 = sin(qJ(6));
t676 = t578 * t642 - t647 * t580;
t768 = t676 * t641;
t455 = -t646 * t513 + t768;
t649 = cos(qJ(3));
t629 = t649 * qJDD(1);
t716 = qJD(1) * qJD(3);
t574 = t644 * t716 + qJDD(4) - t629;
t569 = qJDD(5) + t574;
t558 = qJDD(6) + t569;
t677 = t513 * t641 + t646 * t676;
t801 = t558 * MDP(30) + (-t455 ^ 2 + t677 ^ 2) * MDP(27) + t455 * MDP(26) * t677;
t701 = t649 * t716;
t714 = qJDD(1) * t644;
t723 = qJD(4) * t644;
t786 = -qJD(1) * t723 + qJDD(3);
t502 = qJD(4) * t717 + (t701 + t714) * t648 + t786 * t643;
t730 = qJD(1) * t649;
t503 = t643 * (qJD(3) * (qJD(4) + t730) + t714) - t786 * t648;
t720 = qJD(5) * t647;
t721 = qJD(5) * t642;
t432 = t647 * t502 - t642 * t503 - t578 * t720 - t580 * t721;
t656 = qJD(5) * t676 - t502 * t642 - t647 * t503;
t718 = qJD(6) * t646;
t711 = t646 * t432 - t513 * t718 + t641 * t656;
t719 = qJD(6) * t641;
t408 = t676 * t719 + t711;
t614 = -qJD(4) + t730;
t609 = -qJD(5) + t614;
t697 = t432 * t641 - t646 * t656;
t657 = qJD(6) * t677 - t697;
t600 = -qJD(6) + t609;
t792 = t600 * t677;
t793 = t455 * t600;
t800 = t569 * MDP(23) + (-t513 ^ 2 + t676 ^ 2) * MDP(20) + (-t513 * t609 + t432) * MDP(21) + (t609 * t676 + t656) * MDP(22) - t513 * MDP(19) * t676 + (t657 + t792) * MDP(29) + (t408 + t793) * MDP(28) + t801;
t639 = sin(pkin(11));
t616 = pkin(1) * t639 + pkin(7);
t601 = t616 * qJD(1);
t549 = t644 * qJD(2) + t649 * t601;
t536 = qJD(3) * pkin(8) + t549;
t640 = cos(pkin(11));
t617 = -pkin(1) * t640 - pkin(2);
t570 = -pkin(3) * t649 - pkin(8) * t644 + t617;
t538 = t570 * qJD(1);
t473 = -t536 * t643 + t648 * t538;
t464 = -pkin(9) * t580 + t473;
t457 = -pkin(4) * t614 + t464;
t474 = t536 * t648 + t538 * t643;
t465 = -pkin(9) * t578 + t474;
t463 = t647 * t465;
t426 = t457 * t642 + t463;
t795 = pkin(10) * t513;
t419 = t426 - t795;
t416 = t419 * t719;
t548 = qJD(2) * t649 - t644 * t601;
t535 = -qJD(3) * pkin(3) - t548;
t508 = pkin(4) * t578 + t535;
t458 = pkin(5) * t513 + t508;
t638 = qJ(4) + qJ(5);
t634 = qJ(6) + t638;
t621 = sin(t634);
t622 = cos(t634);
t635 = qJ(1) + pkin(11);
t626 = cos(t635);
t625 = sin(t635);
t763 = t625 * t649;
t524 = t621 * t626 - t622 * t763;
t762 = t626 * t649;
t526 = t621 * t625 + t622 * t762;
t772 = g(3) * t644;
t785 = g(1) * t526 - g(2) * t524 - t458 * t455 + t622 * t772 + t416;
t523 = t621 * t763 + t622 * t626;
t525 = -t621 * t762 + t622 * t625;
t594 = t616 * qJDD(1);
t787 = qJD(3) * t548;
t487 = qJDD(3) * pkin(8) + qJDD(2) * t644 + t594 * t649 + t787;
t685 = pkin(3) * t644 - pkin(8) * t649;
t590 = t685 * qJD(3);
t509 = qJD(1) * t590 + qJDD(1) * t570;
t497 = t648 * t509;
t413 = pkin(4) * t574 - pkin(9) * t502 - qJD(4) * t474 - t487 * t643 + t497;
t722 = qJD(4) * t648;
t710 = t648 * t487 + t643 * t509 + t538 * t722;
t724 = qJD(4) * t643;
t669 = -t536 * t724 + t710;
t420 = -pkin(9) * t503 + t669;
t698 = t647 * t413 - t642 * t420;
t658 = -t426 * qJD(5) + t698;
t402 = pkin(5) * t569 - pkin(10) * t432 + t658;
t687 = -t642 * t413 - t647 * t420 - t457 * t720 + t465 * t721;
t403 = pkin(10) * t656 - t687;
t699 = t646 * t402 - t641 * t403;
t784 = -g(1) * t525 + g(2) * t523 + t458 * t677 + t621 * t772 + t699;
t587 = t685 * qJD(1);
t567 = t648 * t587;
t752 = t648 * t649;
t675 = pkin(4) * t644 - pkin(9) * t752;
t774 = pkin(8) + pkin(9);
t708 = qJD(4) * t774;
t797 = qJD(1) * t675 - t548 * t643 + t648 * t708 + t567;
t707 = t643 * t730;
t740 = t648 * t548 + t643 * t587;
t790 = pkin(9) * t707 - t643 * t708 - t740;
t725 = qJD(3) * t649;
t789 = -qJD(2) * qJD(3) - t594;
t709 = -t601 * t725 + t644 * t789;
t488 = -qJDD(3) * pkin(3) - qJDD(2) * t649 - t709;
t684 = g(1) * t626 + g(2) * t625;
t662 = -g(3) * t649 + t644 * t684;
t796 = qJD(4) * pkin(8) * t614 - t488 + t662;
t794 = pkin(10) * t676;
t581 = t642 * t643 - t647 * t648;
t670 = t581 * t649;
t776 = qJD(4) + qJD(5);
t742 = qJD(1) * t670 - t776 * t581;
t582 = t642 * t648 + t643 * t647;
t741 = (-t730 + t776) * t582;
t705 = t643 * t725;
t788 = t644 * t722 + t705;
t632 = sin(t638);
t633 = cos(t638);
t760 = t633 * t649;
t531 = -t625 * t760 + t626 * t632;
t533 = t625 * t632 + t626 * t760;
t783 = g(1) * t533 - g(2) * t531 + t508 * t513 + t633 * t772 + t687;
t761 = t632 * t649;
t530 = t625 * t761 + t626 * t633;
t532 = t625 * t633 - t626 * t761;
t782 = -g(1) * t532 + g(2) * t530 + t508 * t676 + t632 * t772 + t658;
t550 = t582 * t644;
t779 = t797 * t647;
t553 = t648 * t570;
t755 = t644 * t648;
t764 = t616 * t643;
t489 = -pkin(9) * t755 + t553 + (-pkin(4) - t764) * t649;
t586 = t616 * t752;
t735 = t643 * t570 + t586;
t757 = t643 * t644;
t501 = -pkin(9) * t757 + t735;
t744 = t642 * t489 + t647 * t501;
t686 = -t549 + (-t707 + t724) * pkin(4);
t604 = t774 * t643;
t605 = t774 * t648;
t736 = -t642 * t604 + t647 * t605;
t777 = -t604 * t720 - t605 * t721 - t797 * t642 + t790 * t647;
t665 = -t643 * t723 + t649 * t717;
t773 = pkin(4) * t642;
t461 = t642 * t465;
t425 = t647 * t457 - t461;
t418 = t425 + t794;
t414 = -pkin(5) * t609 + t418;
t770 = t414 * t646;
t769 = t502 * t643;
t767 = t578 * t614;
t766 = t580 * t614;
t765 = t614 * t648;
t759 = t641 * t558;
t758 = t642 * t646;
t756 = t643 * t649;
t754 = t646 * t419;
t753 = t646 * t558;
t751 = qJDD(2) - g(3);
t467 = -qJD(3) * t670 - t550 * t776;
t468 = -t721 * t757 + (t755 * t776 + t705) * t647 + t665 * t642;
t551 = t581 * t644;
t486 = -t550 * t641 - t551 * t646;
t422 = qJD(6) * t486 + t467 * t641 + t646 * t468;
t485 = t646 * t550 - t551 * t641;
t750 = t422 * t600 - t485 * t558;
t517 = t646 * t581 + t582 * t641;
t749 = -qJD(6) * t517 - t641 * t741 + t646 * t742;
t518 = -t581 * t641 + t582 * t646;
t748 = qJD(6) * t518 + t641 * t742 + t646 * t741;
t747 = t468 * t609 - t550 * t569;
t746 = t647 * t464 - t461;
t743 = pkin(5) * t741 + t686;
t739 = t570 * t722 + t643 * t590;
t726 = qJD(3) * t644;
t737 = t648 * t590 + t726 * t764;
t556 = pkin(4) * t757 + t644 * t616;
t636 = t644 ^ 2;
t734 = -t649 ^ 2 + t636;
t602 = qJD(1) * t617;
t728 = qJD(3) * t578;
t521 = pkin(4) * t788 + t616 * t725;
t624 = -pkin(4) * t648 - pkin(3);
t706 = t614 * t727;
t443 = t675 * qJD(3) + (-t586 + (pkin(9) * t644 - t570) * t643) * qJD(4) + t737;
t446 = (-t644 * t717 - t724 * t649) * t616 - t788 * pkin(9) + t739;
t696 = t647 * t443 - t446 * t642;
t695 = -t464 * t642 - t463;
t694 = t647 * t489 - t501 * t642;
t692 = t614 * t616 + t536;
t691 = -t647 * t604 - t605 * t642;
t690 = -qJD(4) * t538 - t487;
t689 = qJD(6) * t414 + t403;
t645 = sin(qJ(1));
t650 = cos(qJ(1));
t683 = g(1) * t645 - g(2) * t650;
t495 = -pkin(10) * t581 + t736;
t682 = pkin(5) * t731 + pkin(10) * t742 + qJD(5) * t736 + qJD(6) * t495 + t642 * t790 + t779;
t494 = -pkin(10) * t582 + t691;
t681 = -pkin(10) * t741 + qJD(6) * t494 + t777;
t405 = t641 * t414 + t754;
t421 = -qJD(6) * t485 + t467 * t646 - t468 * t641;
t680 = t421 * t600 - t486 * t558;
t679 = t467 * t609 + t551 * t569;
t673 = -t574 * t643 + t614 * t722;
t672 = -t574 * t648 - t614 * t724;
t668 = t642 * t443 + t647 * t446 + t489 * t720 - t501 * t721;
t666 = -qJD(1) * t602 + t684;
t664 = -pkin(8) * t574 - t535 * t614;
t663 = 0.2e1 * qJD(3) * t602 - qJDD(3) * t616;
t445 = pkin(4) * t503 + t488;
t651 = qJD(3) ^ 2;
t655 = g(1) * t625 - g(2) * t626 - 0.2e1 * qJDD(1) * t617 - t616 * t651;
t623 = pkin(4) * t647 + pkin(5);
t593 = qJDD(3) * t649 - t644 * t651;
t592 = qJDD(3) * t644 + t649 * t651;
t554 = t580 * t726;
t545 = t625 * t643 + t626 * t752;
t544 = t625 * t648 - t626 * t756;
t543 = -t625 * t752 + t626 * t643;
t542 = t625 * t756 + t626 * t648;
t541 = pkin(5) * t581 + t624;
t504 = pkin(5) * t550 + t556;
t499 = t676 * t726;
t475 = pkin(4) * t580 - pkin(5) * t676;
t447 = t677 * t726;
t440 = pkin(5) * t468 + t521;
t435 = -pkin(10) * t550 + t744;
t434 = -pkin(5) * t649 + pkin(10) * t551 + t694;
t424 = t746 + t794;
t423 = t695 + t795;
t410 = -pkin(5) * t656 + t445;
t407 = -pkin(10) * t468 + t668;
t406 = pkin(5) * t726 - pkin(10) * t467 - qJD(5) * t744 + t696;
t404 = -t419 * t641 + t770;
t1 = [0.2e1 * (t629 * t644 - t716 * t734) * MDP(6) + (-(-t570 * t724 + t737) * t614 + t553 * t574 - g(1) * t543 - g(2) * t545 + (t616 * t728 - t497 + t692 * t722 + (qJD(3) * t535 - t574 * t616 - t690) * t643) * t649 + (qJD(3) * t473 + t488 * t643 + t503 * t616 + t535 * t722) * t644) * MDP(17) + (-t405 * t726 - g(1) * t523 - g(2) * t525 + t504 * t408 + t410 * t486 - t416 * t649 + t458 * t421 - t440 * t677 + ((-qJD(6) * t435 + t406) * t600 - t434 * t558 + t402 * t649) * t641 + ((qJD(6) * t434 + t407) * t600 - t435 * t558 + t689 * t649) * t646) * MDP(32) + (t408 * t486 - t421 * t677) * MDP(26) + (-g(1) * t530 - g(2) * t532 - t426 * t726 + t556 * t432 - t445 * t551 + t508 * t467 - t521 * t676 - t569 * t744 + t609 * t668 - t649 * t687) * MDP(25) + (-t432 * t551 - t467 * t676) * MDP(19) + (-t513 * t726 - t649 * t656 + t747) * MDP(22) + (-t432 * t550 - t467 * t513 + t468 * t676 - t551 * t656) * MDP(20) + (-t696 * t609 + t694 * t569 - t698 * t649 + t425 * t726 + t521 * t513 - t556 * t656 + t445 * t550 + t508 * t468 - g(1) * t531 - g(2) * t533 + (t426 * t649 + t609 * t744) * qJD(5)) * MDP(24) + (-t502 * t649 + t574 * t755 - t614 * t665 + t554) * MDP(14) + (t502 * t755 + t580 * t665) * MDP(12) + (-t408 * t649 - t447 - t680) * MDP(28) + (-t432 * t649 - t499 - t679) * MDP(21) + (-t574 * t649 - t614 * t726) * MDP(16) + (-t569 * t649 - t609 * t726) * MDP(23) + (-t558 * t649 - t600 * t726) * MDP(30) + ((t503 + t706) * t649 + (t673 - t728) * t644) * MDP(15) + ((-t578 * t648 - t580 * t643) * t725 + (-t769 - t503 * t648 + (t578 * t643 - t580 * t648) * qJD(4)) * t644) * MDP(13) + (g(1) * t650 + g(2) * t645) * MDP(3) + (t455 * t726 - t649 * t657 + t750) * MDP(29) + (-t408 * t485 + t421 * t455 + t422 * t677 + t486 * t657) * MDP(27) + (-(t406 * t646 - t407 * t641) * t600 + (t434 * t646 - t435 * t641) * t558 - t699 * t649 + t404 * t726 - t440 * t455 - t504 * t657 + t410 * t485 + t458 * t422 - g(1) * t524 - g(2) * t526 + (-(-t434 * t641 - t435 * t646) * t600 + t405 * t649) * qJD(6)) * MDP(31) + (t644 * t663 + t649 * t655) * MDP(10) + (-t644 * t655 + t649 * t663) * MDP(11) + t683 * MDP(2) + (t739 * t614 - t735 * t574 - g(1) * t542 - g(2) * t544 + (-t692 * t724 + (t535 * t648 + t580 * t616) * qJD(3) + t710) * t649 + (-t535 * t724 + t488 * t648 + t616 * t502 + (-t616 * t765 - t474) * qJD(3)) * t644) * MDP(18) + t592 * MDP(7) + t593 * MDP(8) + (qJDD(1) * t636 + 0.2e1 * t644 * t701) * MDP(5) + (t683 + (t639 ^ 2 + t640 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + qJDD(1) * MDP(1); t751 * MDP(4) + t593 * MDP(10) - t592 * MDP(11) + t554 * MDP(18) + t747 * MDP(24) + (-t499 + t679) * MDP(25) + t750 * MDP(31) + (-t447 + t680) * MDP(32) + ((-t503 + t706) * MDP(17) + (t614 * t717 - t502) * MDP(18) + t656 * MDP(24) - t432 * MDP(25) + t657 * MDP(31) - t408 * MDP(32)) * t649 + (t673 * MDP(17) + t672 * MDP(18) + (t578 * MDP(17) + MDP(24) * t513 - MDP(31) * t455) * qJD(3)) * t644; (qJD(3) * t549 + t644 * t666 + t649 * t751 + t709) * MDP(10) + ((-t580 * t644 + t614 * t752) * qJD(1) - t673) * MDP(14) + (t408 * t518 - t677 * t749) * MDP(26) + (t432 * t582 - t676 * t742) * MDP(19) + (t691 * t569 - t624 * t656 + t445 * t581 + (t605 * t720 + (-qJD(5) * t604 + t790) * t642 + t779) * t609 + t686 * t513 + t741 * t508 + t662 * t633) * MDP(24) + (-(t494 * t641 + t495 * t646) * t558 + t541 * t408 + t410 * t518 + (-t641 * t682 + t646 * t681) * t600 + t749 * t458 - t743 * t677 - t662 * t621) * MDP(32) + (t624 * t432 + t445 * t582 + t742 * t508 - t736 * t569 + t609 * t777 - t632 * t662 - t676 * t686) * MDP(25) + (-MDP(5) * t644 * t649 + MDP(6) * t734) * qJD(1) ^ 2 + ((t578 * t644 - t614 * t756) * qJD(1) - t672) * MDP(15) + ((t502 + t767) * t648 + (-t503 + t766) * t643) * MDP(13) + (-t580 * t765 + t769) * MDP(12) + (t614 * MDP(16) - t473 * MDP(17) + MDP(18) * t474 + MDP(21) * t676 + t513 * MDP(22) + t609 * MDP(23) - t425 * MDP(24) + t426 * MDP(25) + MDP(28) * t677 - MDP(29) * t455 + t600 * MDP(30) - t404 * MDP(31) + t405 * MDP(32)) * t731 + ((t494 * t646 - t495 * t641) * t558 - t541 * t657 + t410 * t517 + (t641 * t681 + t646 * t682) * t600 + t748 * t458 - t743 * t455 + t662 * t622) * MDP(31) + (-t408 * t517 + t455 * t749 + t518 * t657 + t677 * t748) * MDP(27) + (-t432 * t581 - t513 * t742 + t582 * t656 + t676 * t741) * MDP(20) + (-t517 * t558 + t600 * t748) * MDP(29) + (t569 * t582 - t609 * t742) * MDP(21) + (-t569 * t581 + t609 * t741) * MDP(22) + (t787 + (qJD(3) * t601 - t751) * t644 + (t666 + t789) * t649) * MDP(11) + (t518 * t558 - t600 * t749) * MDP(28) + MDP(8) * t629 + MDP(7) * t714 + qJDD(3) * MDP(9) + (-pkin(3) * t503 - t549 * t578 + t567 * t614 + (-t548 * t614 + t664) * t643 + t796 * t648) * MDP(17) + (-pkin(3) * t502 - t549 * t580 - t740 * t614 - t796 * t643 + t664 * t648) * MDP(18); (t475 * t677 + (-t623 * t558 - t402 + (-t423 + (-qJD(5) - qJD(6)) * t773) * t600) * t641 + (-t558 * t773 + (pkin(4) * t720 + qJD(6) * t623 - t424) * t600 - t689) * t646 + t785) * MDP(32) + (t695 * t609 + (-t513 * t580 + t569 * t647 + t609 * t721) * pkin(4) + t782) * MDP(24) + (-t746 * t609 + (-t569 * t642 + t580 * t676 + t609 * t720) * pkin(4) + t783) * MDP(25) + (t623 * t753 + (t423 * t646 - t424 * t641) * t600 + t475 * t455 + (-t642 * t759 - (-t641 * t647 - t758) * t600 * qJD(5)) * pkin(4) + (-(-pkin(4) * t758 - t623 * t641) * t600 - t405) * qJD(6) + t784) * MDP(31) + (-t536 * t722 - g(1) * t544 + g(2) * t542 - t474 * t614 - t535 * t580 + t497 + (t690 + t772) * t643) * MDP(17) + (g(1) * t545 - g(2) * t543 + g(3) * t755 - t473 * t614 + t535 * t578 - t669) * MDP(18) + (-t503 - t766) * MDP(15) + (t502 - t767) * MDP(14) + (-t578 ^ 2 + t580 ^ 2) * MDP(13) + t574 * MDP(16) + t580 * t578 * MDP(12) + t800; (-t426 * t609 + t782) * MDP(24) + (-t425 * t609 + t783) * MDP(25) + ((-t418 * t641 - t754) * t600 - t405 * qJD(6) + (-t455 * t676 + t600 * t719 + t753) * pkin(5) + t784) * MDP(31) + ((t419 * t600 - t402) * t641 + (-t418 * t600 - t689) * t646 + (t600 * t718 - t676 * t677 - t759) * pkin(5) + t785) * MDP(32) + t800; (t711 + t793) * MDP(28) + (-t697 + t792) * MDP(29) + (-t405 * t600 + t784) * MDP(31) + (-t641 * t402 - t646 * t403 - t404 * t600 + t785) * MDP(32) + (MDP(28) * t768 + MDP(29) * t677 - MDP(31) * t405 - MDP(32) * t770) * qJD(6) + t801;];
tau  = t1;
