% Calculate vector of inverse dynamics joint torques for
% S6PRRRPR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:09:17
% EndTime: 2019-03-08 23:09:32
% DurationCPUTime: 11.74s
% Computational Cost: add. (7050->580), mult. (15977->788), div. (0->0), fcn. (12765->18), ass. (0->256)
t613 = sin(qJ(4));
t617 = cos(qJ(3));
t748 = cos(qJ(4));
t682 = qJD(2) * t748;
t614 = sin(qJ(3));
t702 = qJD(2) * t614;
t762 = -t613 * t702 + t617 * t682;
t544 = qJD(6) - t762;
t716 = t613 * t617;
t555 = -qJD(2) * t716 - t614 * t682;
t605 = qJD(3) + qJD(4);
t609 = sin(pkin(12));
t611 = cos(pkin(12));
t528 = t555 * t609 + t611 * t605;
t616 = cos(qJ(6));
t526 = t555 * t611 - t605 * t609;
t612 = sin(qJ(6));
t731 = t526 * t612;
t761 = t528 * t616 + t731;
t766 = t544 * t761;
t765 = t605 * t762;
t615 = sin(qJ(2));
t610 = sin(pkin(6));
t704 = qJD(1) * t610;
t687 = t615 * t704;
t738 = qJD(3) * pkin(3);
t638 = t614 * t738 - t687;
t646 = -t613 * t614 + t617 * t748;
t518 = t605 * t646;
t561 = t614 * t748 + t716;
t519 = t605 * t561;
t764 = pkin(4) * t519 - qJ(5) * t518 - qJD(5) * t561 + t638;
t749 = pkin(9) + pkin(8);
t688 = qJD(3) * t749;
t562 = t614 * t688;
t563 = t617 * t688;
t573 = t749 * t614;
t574 = t749 * t617;
t647 = -t573 * t748 - t613 * t574;
t618 = cos(qJ(2));
t686 = t618 * t704;
t763 = qJD(4) * t647 - t562 * t748 - t613 * t563 - t646 * t686;
t651 = t526 * t616 - t528 * t612;
t759 = t544 * t651;
t737 = cos(pkin(6));
t673 = qJD(1) * t737;
t674 = qJD(2) * t749 + t687;
t524 = t614 * t673 + t617 * t674;
t669 = t737 * qJDD(1);
t579 = t617 * t669;
t696 = qJD(1) * qJD(2);
t537 = qJDD(2) * pkin(8) + (qJDD(1) * t615 + t618 * t696) * t610;
t672 = pkin(9) * qJDD(2) + t537;
t447 = qJDD(3) * pkin(3) - qJD(3) * t524 - t672 * t614 + t579;
t523 = -t674 * t614 + t617 * t673;
t452 = qJD(3) * t523 + t614 * t669 + t672 * t617;
t515 = t523 + t738;
t681 = qJD(4) * t748;
t701 = qJD(4) * t613;
t668 = -t748 * t447 + t613 * t452 + t515 * t701 + t524 * t681;
t603 = qJDD(3) + qJDD(4);
t751 = -pkin(4) * t603 + qJDD(5);
t402 = t668 + t751;
t736 = cos(pkin(11));
t661 = t737 * t736;
t735 = sin(pkin(11));
t546 = t615 * t661 + t618 * t735;
t608 = qJ(3) + qJ(4);
t600 = sin(t608);
t601 = cos(t608);
t676 = t610 * t736;
t507 = t546 * t600 + t601 * t676;
t660 = t737 * t735;
t548 = -t615 * t660 + t618 * t736;
t675 = t610 * t735;
t509 = t548 * t600 - t601 * t675;
t719 = t610 * t615;
t540 = t600 * t719 - t601 * t737;
t640 = g(1) * t509 + g(2) * t507 + g(3) * t540;
t758 = -t402 + t640;
t559 = t609 * t616 + t611 * t612;
t550 = t559 * qJD(6);
t757 = -t559 * t762 + t550;
t558 = t609 * t612 - t616 * t611;
t706 = t544 * t558;
t545 = t615 * t735 - t618 * t661;
t547 = t615 * t736 + t618 * t660;
t663 = g(1) * t547 + g(2) * t545;
t717 = t610 * t618;
t636 = g(3) * t717 - t663;
t756 = t636 * t600;
t627 = t613 * t447 + t452 * t748 + t515 * t681 - t524 * t701;
t400 = t603 * qJ(5) + t605 * qJD(5) + t627;
t677 = qJDD(2) * t748;
t693 = qJDD(2) * t617;
t477 = t613 * t693 + t614 * t677 + t765;
t694 = qJDD(2) * t614;
t658 = t613 * t694 - t617 * t677;
t478 = qJD(2) * t519 + t658;
t597 = pkin(3) * t617 + pkin(2);
t680 = t615 * t696;
t659 = -qJDD(1) * t717 + t610 * t680;
t695 = qJD(2) * qJD(3);
t679 = t614 * t695;
t506 = pkin(3) * t679 - qJDD(2) * t597 + t659;
t411 = pkin(4) * t478 - qJ(5) * t477 + qJD(5) * t555 + t506;
t390 = t611 * t400 + t609 * t411;
t508 = t546 * t601 - t600 * t676;
t510 = t548 * t601 + t600 * t675;
t541 = t600 * t737 + t601 * t719;
t639 = g(1) * t510 + g(2) * t508 + g(3) * t541;
t755 = t390 * t611 - t639;
t712 = -t763 * t609 + t611 * t764;
t711 = t609 * t764 + t763 * t611;
t513 = t613 * t524;
t456 = t515 * t748 - t513;
t501 = -pkin(4) * t555 - qJ(5) * t762;
t431 = t611 * t456 + t609 * t501;
t754 = -qJD(5) * t611 + t431;
t514 = t748 * t524;
t461 = t613 * t523 + t514;
t665 = pkin(3) * t701 - t461;
t533 = -t613 * t573 + t574 * t748;
t708 = qJD(4) * t533 - t561 * t686 - t613 * t562 + t563 * t748;
t462 = t523 * t748 - t513;
t487 = pkin(3) * t702 + t501;
t429 = t611 * t462 + t609 * t487;
t580 = pkin(3) * t681 + qJD(5);
t753 = t580 * t611 - t429;
t428 = -t462 * t609 + t611 * t487;
t752 = -t580 * t609 - t428;
t620 = qJD(3) ^ 2;
t750 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t620 + t610 * (-g(3) * t618 + t680) - t659 + t663;
t741 = g(3) * t610;
t740 = t611 * pkin(5);
t602 = t611 * pkin(10);
t739 = qJD(2) * pkin(2);
t734 = qJ(5) * t611;
t732 = t518 * t609;
t730 = t762 * t609;
t729 = t762 * t611;
t728 = t561 * t609;
t727 = t561 * t611;
t592 = pkin(3) * t613 + qJ(5);
t724 = t592 * t609;
t723 = t592 * t611;
t604 = pkin(12) + qJ(6);
t598 = sin(t604);
t722 = t598 * t601;
t599 = cos(t604);
t721 = t599 * t601;
t720 = t601 * t618;
t718 = t610 * t617;
t715 = qJDD(1) - g(3);
t714 = pkin(5) * t519 - t518 * t602 + t712;
t713 = pkin(10) * t732 - t711;
t710 = pkin(5) * t732 + t708;
t457 = t613 * t515 + t514;
t448 = qJ(5) * t605 + t457;
t542 = -qJD(2) * t597 - t686;
t476 = -pkin(4) * t762 + qJ(5) * t555 + t542;
t420 = t611 * t448 + t609 * t476;
t445 = -t605 * pkin(4) + qJD(5) - t456;
t709 = t445 - t580;
t505 = -pkin(4) * t646 - qJ(5) * t561 - t597;
t455 = t609 * t505 + t611 * t533;
t606 = t614 ^ 2;
t705 = -t617 ^ 2 + t606;
t703 = qJD(2) * t610;
t699 = qJD(6) * t612;
t698 = qJD(6) * t616;
t697 = -qJD(5) + t445;
t692 = pkin(10) * t730;
t463 = t477 * t609 - t611 * t603;
t464 = t477 * t611 + t603 * t609;
t689 = -t612 * t463 + t616 * t464 + t528 * t698;
t685 = t615 * t703;
t684 = t618 * t703;
t678 = t617 * t695;
t389 = -t400 * t609 + t611 * t411;
t386 = pkin(5) * t478 - pkin(10) * t464 + t389;
t387 = -pkin(10) * t463 + t390;
t671 = t616 * t386 - t612 * t387;
t419 = -t448 * t609 + t611 * t476;
t430 = -t456 * t609 + t611 * t501;
t670 = t616 * t463 + t612 * t464;
t454 = t611 * t505 - t533 * t609;
t596 = -pkin(3) * t748 - pkin(4);
t667 = -t555 * pkin(5) - pkin(10) * t729;
t539 = pkin(5) * t730;
t666 = -t539 + t665;
t662 = g(1) * t548 + g(2) * t546;
t656 = t612 * t386 + t616 * t387;
t407 = -pkin(5) * t762 + pkin(10) * t526 + t419;
t408 = pkin(10) * t528 + t420;
t392 = t407 * t616 - t408 * t612;
t393 = t407 * t612 + t408 * t616;
t435 = -pkin(5) * t646 - pkin(10) * t727 + t454;
t440 = -pkin(10) * t728 + t455;
t655 = t435 * t616 - t440 * t612;
t654 = t435 * t612 + t440 * t616;
t551 = -t614 * t719 + t617 * t737;
t552 = t614 * t737 + t615 * t718;
t491 = t613 * t551 + t552 * t748;
t473 = -t491 * t609 - t611 * t717;
t474 = t491 * t611 - t609 * t717;
t653 = t473 * t616 - t474 * t612;
t652 = t473 * t612 + t474 * t616;
t648 = t551 * t748 - t613 * t552;
t557 = t602 + t723;
t645 = qJD(6) * t557 + t667 - t752;
t556 = (-pkin(10) - t592) * t609;
t644 = -qJD(6) * t556 - t692 - t753;
t570 = t602 + t734;
t643 = qJD(5) * t609 + qJD(6) * t570 + t430 + t667;
t569 = (-pkin(10) - qJ(5)) * t609;
t642 = -qJD(6) * t569 - t692 + t754;
t641 = -t420 * t555 - t609 * t758;
t404 = t526 * t699 + t689;
t637 = g(1) * t735 - g(2) * t736;
t565 = -t686 - t739;
t635 = -t565 * qJD(2) - t537 + t662;
t634 = t636 * t601;
t633 = t402 * t561 + t445 * t518 - t662;
t632 = t419 * t729 + t420 * t730 + t755;
t396 = pkin(5) * t463 + t402;
t434 = -pkin(5) * t528 + t445;
t631 = -t393 * t555 + t396 * t559 - t434 * t706 - t598 * t640;
t629 = t640 - t668;
t405 = -qJD(6) * t651 + t670;
t628 = -pkin(8) * qJDD(3) + (t565 + t686 - t739) * qJD(3);
t626 = t419 * t555 + t611 * t758;
t625 = t542 * t555 + t629;
t475 = qJDD(6) + t478;
t624 = (-t404 * t558 - t405 * t559 + t651 * t757 - t706 * t761) * MDP(24) + (t404 * t559 + t651 * t706) * MDP(23) + (t475 * t559 - t544 * t706 - t555 * t651) * MDP(25) + (-t475 * t558 - t544 * t757 + t555 * t761) * MDP(26) + (t477 - t765) * MDP(14) + (-t658 + (-qJD(2) * t561 - t555) * t605) * MDP(15) + (t555 ^ 2 - t762 ^ 2) * MDP(13) + t603 * MDP(16) + (MDP(12) * t762 + MDP(27) * t544) * t555;
t623 = t392 * t555 + t396 * t558 + t434 * t757 + t599 * t640;
t622 = -t542 * t762 - t627 + t639;
t621 = qJD(2) ^ 2;
t593 = -pkin(4) - t740;
t568 = t596 - t740;
t535 = t540 * pkin(4);
t522 = -qJD(3) * t552 - t614 * t684;
t521 = qJD(3) * t551 + t617 * t684;
t504 = t509 * pkin(4);
t503 = t507 * pkin(4);
t495 = t558 * t561;
t494 = t559 * t561;
t484 = pkin(5) * t728 - t647;
t439 = t539 + t457;
t433 = qJD(4) * t491 + t613 * t521 - t522 * t748;
t432 = qJD(4) * t648 + t521 * t748 + t613 * t522;
t427 = t518 * t559 + t698 * t727 - t699 * t728;
t426 = -t518 * t558 - t550 * t561;
t423 = t432 * t611 + t609 * t685;
t422 = -t432 * t609 + t611 * t685;
t1 = [t715 * MDP(1) + (qJD(3) * t522 + qJDD(3) * t551) * MDP(10) + (-qJD(3) * t521 - qJDD(3) * t552) * MDP(11) + (-t433 * t605 + t603 * t648) * MDP(17) + (-t432 * t605 - t491 * t603) * MDP(18) + (-t422 * t762 - t433 * t528 - t463 * t648 + t473 * t478) * MDP(19) + (t423 * t762 - t433 * t526 - t464 * t648 - t474 * t478) * MDP(20) + (t422 * t526 + t423 * t528 - t463 * t474 - t464 * t473) * MDP(21) + (t389 * t473 + t390 * t474 - t402 * t648 + t419 * t422 + t420 * t423 + t433 * t445 - g(3)) * MDP(22) + ((-qJD(6) * t652 + t422 * t616 - t423 * t612) * t544 + t653 * t475 - t433 * t761 - t648 * t405) * MDP(28) + (-(qJD(6) * t653 + t422 * t612 + t423 * t616) * t544 - t652 * t475 - t433 * t651 - t648 * t404) * MDP(29) + ((-qJDD(2) * MDP(4) + (-MDP(17) * t762 - t555 * MDP(18)) * qJD(2) + (-MDP(10) * t617 + MDP(11) * t614 - MDP(3)) * t621) * t615 + (qJDD(2) * MDP(3) - t621 * MDP(4) + (-t679 + t693) * MDP(10) + (-t678 - t694) * MDP(11) - t478 * MDP(17) - t477 * MDP(18)) * t618) * t610; (t477 * t561 - t518 * t555) * MDP(12) + (t655 * t475 - t671 * t646 + t392 * t519 + t484 * t405 + t396 * t494 + t434 * t427 - g(1) * (-t547 * t721 + t548 * t598) - g(2) * (-t545 * t721 + t546 * t598) - (t598 * t615 + t599 * t720) * t741 + (t612 * t713 + t616 * t714) * t544 - t710 * t761 + (t393 * t646 - t544 * t654) * qJD(6)) * MDP(28) + (t405 * t646 - t427 * t544 - t475 * t494 + t519 * t761) * MDP(26) + (-t404 * t494 + t405 * t495 + t426 * t761 + t427 * t651) * MDP(24) + (t477 * t646 - t478 * t561 + t518 * t762 + t519 * t555) * MDP(13) + (-t478 * t597 - t506 * t646 + t519 * t542 + t603 * t647 - t605 * t708 - t638 * t762 - t634) * MDP(17) + (-t389 * t646 + t419 * t519 + t454 * t478 - t647 * t463 - t611 * t634 + (-g(3) * t719 + t633) * t609 - t712 * t762 - t708 * t528) * MDP(19) + (t390 * t646 - t420 * t519 - t455 * t478 - t647 * t464 - t663 * t609 * t601 + t633 * t611 - (-t609 * t720 + t611 * t615) * t741 + t711 * t762 - t708 * t526) * MDP(20) + (t518 * t605 + t561 * t603) * MDP(14) + (-t454 * t464 - t455 * t463 + (-t389 * t611 - t390 * t609) * t561 + t712 * t526 + t711 * t528 + (-t419 * t611 - t420 * t609) * t518 - t756) * MDP(21) + (t389 * t454 + t390 * t455 - t402 * t647 + t712 * t419 + t711 * t420 + t708 * t445 - (t615 * t741 + t662) * t749 + (-t618 * t741 + t663) * (pkin(4) * t601 + qJ(5) * t600 + t597)) * MDP(22) + (-t519 * t605 + t603 * t646) * MDP(15) + (-t475 * t646 + t519 * t544) * MDP(27) + (-t404 * t646 + t426 * t544 - t475 * t495 - t519 * t651) * MDP(25) + (-t654 * t475 + t656 * t646 - t393 * t519 + t484 * t404 - t396 * t495 + t434 * t426 - g(1) * (t547 * t722 + t548 * t599) - g(2) * (t545 * t722 + t546 * t599) - (-t598 * t720 + t599 * t615) * t741 + (-t612 * t714 + t616 * t713) * t544 - t710 * t651 + (t392 * t646 - t544 * t655) * qJD(6)) * MDP(29) + qJDD(2) * MDP(2) + (-t404 * t495 - t426 * t651) * MDP(23) + (qJDD(3) * t614 + t617 * t620) * MDP(7) + (qJDD(3) * t617 - t614 * t620) * MDP(8) + (-t477 * t597 + t506 * t561 + t518 * t542 - t533 * t603 - t638 * t555 - t605 * t763 + t756) * MDP(18) + (-t715 * t719 + t662) * MDP(4) + (t715 * t717 + t663) * MDP(3) + 0.2e1 * (t614 * t693 - t695 * t705) * MDP(6) + (qJDD(2) * t606 + 0.2e1 * t614 * t678) * MDP(5) + (t628 * t614 + t617 * t750) * MDP(10) + (-t614 * t750 + t628 * t617) * MDP(11); (-(t556 * t612 + t557 * t616) * t475 + t568 * t404 + (t612 * t645 + t616 * t644) * t544 - t666 * t651 + t631) * MDP(29) + (-t478 * t723 + t464 * t596 - t665 * t526 - (t611 * t709 + t429) * t762 + t641) * MDP(20) + (t390 * t723 - t389 * t724 + t402 * t596 - g(1) * (t510 * qJ(5) - t504 + (-t548 * t614 + t617 * t675) * pkin(3)) - g(2) * (t508 * qJ(5) - t503 + (-t546 * t614 - t617 * t676) * pkin(3)) - g(3) * (pkin(3) * t551 + t541 * qJ(5) - t535) + t665 * t445 + t753 * t420 + t752 * t419) * MDP(22) + (-t463 * t723 - t428 * t526 + t753 * t528 + (t464 * t592 - t526 * t580 - t389) * t609 + t632) * MDP(21) + MDP(7) * t694 + (-t478 * t724 + t463 * t596 - t665 * t528 - (t609 * t709 - t428) * t762 + t626) * MDP(19) + (t462 * t605 + (t555 * t702 - t603 * t613 - t605 * t681) * pkin(3) + t622) * MDP(18) + qJDD(3) * MDP(9) + MDP(8) * t693 + (t461 * t605 + (t603 * t748 - t605 * t701 + t702 * t762) * pkin(3) + t625) * MDP(17) + (-g(3) * t551 + t614 * t635 - t637 * t718 + t579) * MDP(10) + ((t556 * t616 - t557 * t612) * t475 + t568 * t405 + (t612 * t644 - t616 * t645) * t544 - t666 * t761 + t623) * MDP(28) + t624 + (g(3) * t552 + (t610 * t637 - t669) * t614 + t635 * t617) * MDP(11) + (-MDP(5) * t614 * t617 + MDP(6) * t705) * t621; (-(t569 * t612 + t570 * t616) * t475 + t593 * t404 + t439 * t651 + (t612 * t643 + t616 * t642) * t544 + t631) * MDP(29) + (-t478 * t734 - pkin(4) * t464 + t457 * t526 - (t611 * t697 + t431) * t762 + t641) * MDP(20) + (-t402 * pkin(4) + g(1) * t504 + g(2) * t503 + g(3) * t535 - t419 * t430 - t420 * t431 - t445 * t457 + (-t419 * t609 + t420 * t611) * qJD(5) + (-t389 * t609 + t755) * qJ(5)) * MDP(22) + (-t463 * t734 - t430 * t526 - t754 * t528 + (qJ(5) * t464 - qJD(5) * t526 - t389) * t609 + t632) * MDP(21) + (t457 * t605 + t625) * MDP(17) + (t456 * t605 + t622) * MDP(18) + (-qJ(5) * t478 * t609 - pkin(4) * t463 + t457 * t528 - (t609 * t697 - t430) * t762 + t626) * MDP(19) + ((t569 * t616 - t570 * t612) * t475 + t593 * t405 + t439 * t761 + (t612 * t642 - t616 * t643) * t544 + t623) * MDP(28) + t624; (t526 * t762 + t463) * MDP(19) + (-t528 * t762 + t464) * MDP(20) + (-t526 ^ 2 - t528 ^ 2) * MDP(21) + (-t419 * t526 - t420 * t528 - t629 + t751) * MDP(22) + (t405 - t759) * MDP(28) + (t404 + t766) * MDP(29); t651 * t761 * MDP(23) + (t651 ^ 2 - t761 ^ 2) * MDP(24) + (t689 - t766) * MDP(25) + (-t670 - t759) * MDP(26) + t475 * MDP(27) + (t393 * t544 + t434 * t651 - g(1) * (-t510 * t598 + t547 * t599) - g(2) * (-t508 * t598 + t545 * t599) - g(3) * (-t541 * t598 - t599 * t717) + t671) * MDP(28) + (t392 * t544 - t434 * t761 - g(1) * (-t510 * t599 - t547 * t598) - g(2) * (-t508 * t599 - t545 * t598) - g(3) * (-t541 * t599 + t598 * t717) - t656) * MDP(29) + (MDP(25) * t731 + MDP(26) * t651 - MDP(28) * t393 - MDP(29) * t392) * qJD(6);];
tau  = t1;
