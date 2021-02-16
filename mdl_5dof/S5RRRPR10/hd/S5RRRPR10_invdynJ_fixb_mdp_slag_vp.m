% Calculate vector of inverse dynamics joint torques for
% S5RRRPR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:43
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRPR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:41:43
% EndTime: 2021-01-15 23:42:16
% DurationCPUTime: 14.75s
% Computational Cost: add. (7531->610), mult. (19002->836), div. (0->0), fcn. (15063->14), ass. (0->268)
t619 = cos(qJ(2));
t700 = qJD(1) * qJD(2);
t673 = t619 * t700;
t615 = sin(qJ(2));
t698 = qJDD(1) * t615;
t781 = t673 + t698;
t609 = sin(pkin(5));
t709 = qJD(1) * t619;
t681 = t609 * t709;
t780 = qJD(3) - t681;
t613 = sin(qJ(5));
t617 = cos(qJ(5));
t611 = cos(pkin(5));
t710 = qJD(1) * t611;
t591 = qJD(2) + t710;
t618 = cos(qJ(3));
t614 = sin(qJ(3));
t711 = qJD(1) * t609;
t682 = t615 * t711;
t658 = t614 * t682;
t520 = -t618 * t591 + t658;
t522 = t591 * t614 + t618 * t682;
t608 = sin(pkin(10));
t610 = cos(pkin(10));
t647 = -t520 * t608 + t610 * t522;
t454 = t613 * t647 - t617 * t780;
t475 = t610 * t520 + t522 * t608;
t771 = qJD(5) + t475;
t779 = t454 * t771;
t663 = t771 * t617;
t693 = pkin(1) * t710;
t538 = -pkin(7) * t682 + t619 * t693;
t653 = pkin(2) * t615 - pkin(8) * t619;
t539 = t653 * t711;
t664 = -t538 * t614 + t618 * t539;
t612 = qJ(4) + pkin(8);
t670 = qJD(3) * t612;
t722 = t618 * t619;
t778 = -(pkin(3) * t615 - qJ(4) * t722) * t711 - t664 - qJD(4) * t614 - t618 * t670;
t657 = t614 * t681;
t716 = t618 * t538 + t614 * t539;
t777 = -qJ(4) * t657 - qJD(4) * t618 + t614 * t670 + t716;
t706 = qJD(3) * t614;
t776 = -t657 + t706;
t558 = t608 * t614 - t610 * t618;
t775 = t780 * t558;
t559 = t608 * t618 + t610 * t614;
t714 = t780 * t559;
t774 = t781 * t609;
t697 = qJDD(1) * t619;
t589 = t609 * t697;
t674 = t615 * t700;
t656 = t609 * t674;
t536 = qJDD(3) - t589 + t656;
t705 = qJD(3) * t618;
t773 = -t614 * t536 - t705 * t780;
t699 = qJDD(1) * t611;
t590 = qJDD(2) + t699;
t467 = -qJD(3) * t658 + t614 * t590 + t591 * t705 + t618 * t774;
t707 = qJD(2) * t619;
t678 = t614 * t707;
t468 = -t618 * t590 + t591 * t706 + t609 * (qJD(1) * (t615 * t705 + t678) + t614 * t698);
t440 = t467 * t608 + t610 * t468;
t439 = qJDD(5) + t440;
t731 = t613 * t439;
t772 = -t663 * t771 - t731;
t762 = cos(qJ(1));
t684 = t762 * t615;
t616 = sin(qJ(1));
t725 = t616 * t619;
t634 = -t611 * t725 - t684;
t734 = t609 * t619;
t644 = g(1) * t634 + g(3) * t734;
t683 = t762 * t619;
t727 = t615 * t616;
t552 = t611 * t727 - t683;
t605 = qJ(3) + pkin(10);
t602 = sin(t605);
t603 = cos(t605);
t735 = t609 * t616;
t495 = t552 * t602 + t603 * t735;
t553 = -t611 * t684 - t725;
t685 = t609 * t762;
t632 = t553 * t602 - t603 * t685;
t736 = t609 * t615;
t770 = g(2) * t632 + g(3) * (-t602 * t736 + t603 * t611) + g(1) * t495;
t721 = t608 * t777 + t610 * t778;
t719 = t608 * t778 - t610 * t777;
t541 = pkin(7) * t681 + t615 * t693;
t654 = pkin(3) * t776 - t541;
t761 = pkin(1) * t615;
t713 = pkin(7) * t734 + t611 * t761;
t533 = pkin(8) * t611 + t713;
t645 = -pkin(2) * t619 - pkin(8) * t615 - pkin(1);
t534 = t645 * t609;
t717 = t618 * t533 + t614 * t534;
t435 = t617 * t439;
t704 = qJD(5) * t613;
t768 = t704 * t771 - t435;
t766 = pkin(3) * t468 + qJDD(4);
t565 = t602 * t685;
t741 = t553 * t603;
t633 = t565 + t741;
t554 = -t611 * t683 + t727;
t724 = t617 * t554;
t765 = t613 * t633 + t724;
t763 = (pkin(3) * t522 + pkin(4) * t647 + pkin(9) * t475) * t771 + t770;
t760 = pkin(2) * t590;
t757 = g(2) * t552;
t756 = g(2) * t634;
t755 = g(2) * t616;
t754 = g(2) * t618;
t752 = g(3) * t609;
t751 = MDP(6) * t609;
t509 = pkin(8) * t591 + t541;
t513 = qJD(1) * t534;
t466 = t509 * t618 + t513 * t614;
t692 = pkin(1) * qJD(2) * t611;
t660 = qJD(1) * t692;
t690 = pkin(1) * t699;
t686 = -pkin(7) * t589 - t615 * t690 - t619 * t660;
t626 = -pkin(7) * t656 - t686;
t482 = pkin(8) * t590 + t626;
t637 = t653 * qJD(2);
t484 = (qJD(1) * t637 + qJDD(1) * t645) * t609;
t622 = -qJD(3) * t466 - t614 * t482 + t618 * t484;
t409 = pkin(3) * t536 - qJ(4) * t467 - qJD(4) * t522 + t622;
t636 = -t618 * t482 - t614 * t484 + t509 * t706 - t513 * t705;
t411 = -qJ(4) * t468 - qJD(4) * t520 - t636;
t397 = t610 * t409 - t411 * t608;
t395 = -pkin(4) * t536 - t397;
t750 = t395 * t613;
t441 = t467 * t610 - t468 * t608;
t703 = qJD(5) * t617;
t412 = t617 * t441 + t613 * t536 - t647 * t704 + t703 * t780;
t749 = t412 * t613;
t451 = -qJ(4) * t520 + t466;
t748 = t451 * t608;
t747 = t454 * t647;
t456 = t613 * t780 + t617 * t647;
t746 = t456 * t647;
t745 = t475 * t613;
t744 = t520 * t780;
t743 = t522 * t780;
t742 = t552 * t603;
t740 = t634 * t617;
t739 = t559 * t617;
t738 = t590 * MDP(8);
t604 = t609 ^ 2;
t737 = t604 * qJD(1) ^ 2;
t448 = t610 * t451;
t733 = t611 * t618;
t732 = t611 * t619;
t537 = t613 * t554;
t730 = t613 * t619;
t728 = t614 * t615;
t726 = t615 * t618;
t723 = t617 * t619;
t398 = t608 * t409 + t610 * t411;
t549 = t609 * t728 - t733;
t679 = t609 * t707;
t499 = -qJD(3) * t549 + t618 * t679;
t550 = t609 * t726 + t611 * t614;
t540 = t609 * t637;
t592 = pkin(7) * t736;
t542 = (t732 * pkin(1) - t592) * qJD(2);
t623 = -qJD(3) * t717 + t618 * t540 - t542 * t614;
t708 = qJD(2) * t615;
t680 = t609 * t708;
t425 = pkin(3) * t680 - qJ(4) * t499 - qJD(4) * t550 + t623;
t498 = qJD(3) * t550 + t609 * t678;
t635 = -t533 * t706 + t534 * t705 + t614 * t540 + t618 * t542;
t429 = -qJ(4) * t498 - qJD(4) * t549 + t635;
t406 = t608 * t425 + t610 * t429;
t465 = -t509 * t614 + t618 * t513;
t450 = -qJ(4) * t522 + t465;
t446 = pkin(3) * t780 + t450;
t418 = t608 * t446 + t448;
t665 = -t533 * t614 + t618 * t534;
t453 = -pkin(3) * t734 - qJ(4) * t550 + t665;
t461 = -qJ(4) * t549 + t717;
t428 = t608 * t453 + t610 * t461;
t720 = pkin(4) * t682 - t721;
t715 = -t617 * t565 + t537;
t543 = pkin(7) * t679 + t615 * t692;
t606 = t615 ^ 2;
t712 = -t619 ^ 2 + t606;
t702 = qJD(2) - t591;
t696 = 0.2e1 * t604;
t695 = g(3) * t736;
t689 = t619 * t737;
t688 = t609 * t730;
t687 = t609 * t723;
t601 = pkin(3) * t618 + pkin(2);
t675 = t612 * t614;
t396 = pkin(9) * t536 + t398;
t659 = pkin(7) * t774 + t615 * t660 - t619 * t690;
t483 = t659 - t760;
t443 = t483 + t766;
t402 = pkin(4) * t440 - pkin(9) * t441 + t443;
t669 = -t613 * t396 + t617 * t402;
t668 = t441 * t613 - t617 * t536;
t667 = t613 * t775 - t617 * t682;
t666 = t613 * t682 + t617 * t775;
t662 = t591 + t710;
t661 = t590 + t699;
t479 = pkin(3) * t498 + t543;
t652 = -g(1) * t552 - g(2) * t553;
t651 = -g(1) * t554 - t756;
t493 = pkin(4) * t558 - pkin(9) * t559 - t601;
t649 = pkin(9) * t682 - qJD(5) * t493 - t719;
t577 = t612 * t618;
t507 = t610 * t577 - t608 * t675;
t648 = -t714 * pkin(4) - pkin(9) * t775 + qJD(5) * t507 - t654;
t415 = pkin(9) * t780 + t418;
t508 = -pkin(2) * t591 - t538;
t473 = pkin(3) * t520 + qJD(4) + t508;
t421 = pkin(4) * t475 - pkin(9) * t647 + t473;
t401 = t415 * t617 + t421 * t613;
t400 = -t415 * t613 + t421 * t617;
t405 = t425 * t610 - t429 * t608;
t417 = t446 * t610 - t748;
t427 = t453 * t610 - t461 * t608;
t646 = t601 * t619 + t612 * t615;
t551 = -t601 * t615 + t612 * t619;
t643 = g(1) * t762 + t755;
t642 = -g(1) * t616 + g(2) * t762;
t641 = -t745 * t771 - t768;
t532 = t592 + (-pkin(1) * t619 - pkin(2)) * t611;
t491 = -t549 * t608 + t550 * t610;
t470 = t491 * t613 + t687;
t640 = t602 * t735 - t742;
t639 = -t552 * t618 + t614 * t735;
t638 = t609 * t618 + t611 * t728;
t631 = t559 * t703 - t667;
t630 = -t559 * t704 - t666;
t627 = -g(2) * t554 + t644;
t492 = pkin(3) * t549 + t532;
t531 = t602 * t611 + t603 * t736;
t625 = g(1) * t640 + g(3) * t531;
t624 = -t613 * t640 - t740;
t621 = -t627 - t659;
t598 = -pkin(3) * t610 - pkin(4);
t597 = pkin(3) * t608 + pkin(9);
t546 = pkin(1) + t646;
t519 = t646 * t611;
t506 = t577 * t608 + t610 * t675;
t501 = -t553 * t614 + t618 * t685;
t497 = t609 * (pkin(3) * t614 + pkin(7)) + t551 * t611;
t490 = t610 * t549 + t550 * t608;
t471 = t491 * t617 - t688;
t460 = -t498 * t608 + t499 * t610;
t459 = t610 * t498 + t499 * t608;
t442 = pkin(4) * t490 - pkin(9) * t491 + t492;
t431 = -qJD(5) * t688 + t460 * t613 + t491 * t703 - t617 * t680;
t430 = -qJD(5) * t470 + t460 * t617 + t613 * t680;
t424 = -pkin(9) * t734 + t428;
t423 = pkin(4) * t734 - t427;
t420 = t450 * t610 - t748;
t419 = t450 * t608 + t448;
t416 = pkin(4) * t459 - pkin(9) * t460 + t479;
t414 = -pkin(4) * t780 - t417;
t413 = qJD(5) * t456 + t668;
t404 = pkin(9) * t680 + t406;
t403 = -pkin(4) * t680 - t405;
t394 = -qJD(5) * t401 + t669;
t393 = qJD(5) * t400 + t617 * t396 + t613 * t402;
t1 = [(g(1) * t632 - g(2) * t495 - t406 * t780 - t428 * t536 + t441 * t492 + t443 * t491 + t460 * t473 + t647 * t479 + (t398 * t619 - t418 * t708) * t609) * MDP(19) + (t405 * t780 + t427 * t536 + t479 * t475 + t492 * t440 + t443 * t490 + t473 * t459 - g(1) * t633 + g(2) * t742 + (-t397 * t619 + t417 * t708 - t602 * t755) * t609) * MDP(18) + (-t635 * t780 - t717 * t536 + t543 * t522 + t532 * t467 + t483 * t550 + t508 * t499 - g(1) * t501 - t614 * t757 + (-t466 * t708 - t616 * t754 - t619 * t636) * t609) * MDP(17) + (-t498 * t780 - t536 * t549 + (t468 * t619 - t520 * t708) * t609) * MDP(14) + (t499 * t780 + t536 * t550 + (-t467 * t619 + t522 * t708) * t609) * MDP(13) + (t623 * t780 + t665 * t536 - t622 * t734 + t465 * t680 + t543 * t520 + t532 * t468 + t483 * t549 + t508 * t498 - g(1) * ((t609 * t614 - t611 * t726) * t762 - t616 * t722) - g(2) * t639) * MDP(16) - t642 * MDP(2) + t643 * MDP(3) + (-t781 * pkin(1) * t696 - t542 * t591 - t713 * t590 - t626 * t611 + t651) * MDP(10) + (-t397 * t491 - t398 * t490 - t405 * t647 - t406 * t475 - t417 * t460 - t418 * t459 - t427 * t441 - t428 * t440 - t651) * MDP(20) + (t398 * t428 + t418 * t406 + t397 * t427 + t417 * t405 + t443 * t492 + t473 * t479 - g(1) * (t497 * t762 - t546 * t616) - g(2) * (t497 * t616 + t546 * t762)) * MDP(21) + (t619 * t661 - t662 * t708) * t609 * MDP(7) + (-t536 * t619 + t708 * t780) * t609 * MDP(15) + (-(t404 * t617 + t416 * t613 + (-t424 * t613 + t442 * t617) * qJD(5)) * t771 - (t424 * t617 + t442 * t613) * t439 - t393 * t490 - t401 * t459 + t403 * t456 + t423 * t412 + t395 * t471 + t414 * t430 + g(1) * t765 - g(2) * t624) * MDP(28) + (t394 * t490 + t400 * t459 + t403 * t454 + t423 * t413 + t395 * t470 + t414 * t431 + g(1) * t715 + ((-qJD(5) * t442 - t404) * t771 - t424 * t439 + t756) * t613 + ((-qJD(5) * t424 + t416) * t771 + t442 * t439 - g(1) * t741 - g(2) * t640) * t617) * MDP(27) + (t439 * t490 + t459 * t771) * MDP(26) + (-t413 * t490 - t431 * t771 - t439 * t470 - t454 * t459) * MDP(25) + (t412 * t490 + t430 * t771 + t439 * t471 + t456 * t459) * MDP(24) + qJDD(1) * MDP(1) + (0.2e1 * (t615 * t697 - t700 * t712) * MDP(5) + (qJDD(1) * t606 + 0.2e1 * t615 * t673) * MDP(4)) * t604 + (t615 * t661 + t662 * t707) * t751 + t611 * t738 + (-t543 * t591 - t592 * t590 - t659 * t611 - g(1) * t553 + t757 + (t590 * t732 + (-t674 + t697) * t696) * pkin(1)) * MDP(9) + (t412 * t471 + t430 * t456) * MDP(22) + (-t412 * t470 - t413 * t471 - t430 * t454 - t431 * t456) * MDP(23) + (t467 * t550 + t499 * t522) * MDP(11) + (-t467 * t549 - t468 * t550 - t498 * t522 - t499 * t520) * MDP(12); (-t397 * t559 - t398 * t558 + t417 * t775 - t418 * t714 - t440 * t507 + t441 * t506 - t475 * t719 - t647 * t721 - t652 - t695) * MDP(20) + (pkin(1) * t689 + t538 * t591 + (pkin(7) * t700 + g(3)) * t736 + t652 + t686) * MDP(10) + ((-pkin(8) * t536 + t508 * t780) * MDP(17) + t536 * MDP(14) + t743 * MDP(11) + (t467 - t744) * MDP(12) + (-t483 - t644) * MDP(16)) * t618 + (t418 * t682 - t441 * t601 + t443 * t559 - t473 * t775 - t507 * t536 + t602 * t627 + t647 * t654 - t719 * t780) * MDP(19) + (-t780 * t706 + (t614 * t619 * t780 + t520 * t615) * t711) * MDP(14) + (-t417 * t682 - t440 * t601 + t443 * t558 + t473 * t714 + t475 * t654 - t506 * t536 - t603 * t627 + t721 * t780) * MDP(18) + (-pkin(2) * t468 + t773 * pkin(8) - t465 * t682 + t508 * t776 - t541 * t520 + t554 * t754 - t664 * t780) * MDP(16) + ((-t522 * t615 - t722 * t780) * t711 - t773) * MDP(13) + (-pkin(2) * t467 + t716 * t780 + t466 * t682 - t541 * t522 + (pkin(8) * qJD(3) * t780 + t483 + t627) * t614) * MDP(17) - t780 * MDP(15) * t682 + (-t468 - t743) * t614 * MDP(12) + t467 * t614 * MDP(11) + t712 * MDP(5) * t737 + (-t682 * t702 + t589) * MDP(7) + t738 + (t541 * t591 + t737 * t761 + t621) * MDP(9) + (t398 * t507 - t397 * t506 - t443 * t601 - g(1) * (-t519 * t616 + t551 * t762) - g(2) * (t519 * t762 + t551 * t616) - t646 * t752 + t654 * t473 + t719 * t418 + t721 * t417) * MDP(21) + (t412 * t739 + t456 * t630) * MDP(22) + (t667 * t456 + t666 * t454 + (-t749 - t413 * t617 + (t454 * t613 - t456 * t617) * qJD(5)) * t559) * MDP(23) + (t412 * t558 + t435 * t559 + t456 * t714 + t630 * t771) * MDP(24) + (t439 * t558 + t714 * t771) * MDP(26) + (-(t493 * t613 + t507 * t617) * t439 - t393 * t558 + t506 * t412 + t395 * t739 - g(1) * (-t603 * t613 * t634 - t552 * t617) - g(2) * (t537 * t603 - t553 * t617) - (-t603 * t730 + t615 * t617) * t752 + (t613 * t648 + t617 * t649) * t771 + t720 * t456 - t714 * t401 + t630 * t414) * MDP(28) + ((t493 * t617 - t507 * t613) * t439 + t394 * t558 + t506 * t413 + t559 * t750 - g(1) * (-t552 * t613 + t603 * t740) - g(2) * (-t553 * t613 - t603 * t724) - (t603 * t723 + t613 * t615) * t752 + (t613 * t649 - t617 * t648) * t771 + t720 * t454 + t714 * t400 + t631 * t414) * MDP(27) + (-t413 * t558 - t454 * t714 - t559 * t731 - t631 * t771) * MDP(25) + (t702 * t709 + t698) * t751 - t615 * MDP(4) * t689; t522 * t520 * MDP(11) + (-t520 ^ 2 + t522 ^ 2) * MDP(12) + (t467 + t744) * MDP(13) + (-t468 + t743) * MDP(14) + t536 * MDP(15) + (t466 * t780 - t508 * t522 - g(1) * (-t614 * t683 + t616 * t638) + g(2) * t501 + g(3) * t549 + t622) * MDP(16) + (t465 * t780 + t508 * t520 + g(1) * t639 - g(2) * (t553 * t618 + t614 * t685) + g(3) * t550 + t636) * MDP(17) + (t419 * t780 - t473 * t647 + (-t475 * t522 + t536 * t610) * pkin(3) + t397 - t770) * MDP(18) + (t420 * t780 + t473 * t475 - g(2) * t633 + (-t522 * t647 - t536 * t608) * pkin(3) + t625 - t398) * MDP(19) + ((-t440 * t608 - t441 * t610) * pkin(3) + (t418 - t419) * t647 + (-t417 + t420) * t475) * MDP(20) + (t417 * t419 - t418 * t420 + (t642 * t638 - g(3) * t733 + t397 * t610 + t398 * t608 - t473 * t522 + (t619 * t643 + t695) * t614) * pkin(3)) * MDP(21) + (t456 * t663 + t749) * MDP(22) + ((t412 - t779) * t617 + (-t456 * t771 - t413) * t613) * MDP(23) + (-t746 - t772) * MDP(24) + (t641 + t747) * MDP(25) - t771 * t647 * MDP(26) + (-t400 * t647 + t598 * t413 - t419 * t454 - t597 * t731 + (t420 * t613 - t597 * t703) * t771 + (t704 + t745) * t414 + (-t395 - t763) * t617) * MDP(27) + (t420 * t663 + t401 * t647 + t598 * t412 - t419 * t456 + t750 + t768 * t597 + (t475 * t617 + t703) * t414 + t763 * t613) * MDP(28); (t647 * t780 + t440) * MDP(18) + (-t475 * t780 + t441) * MDP(19) + (-t475 ^ 2 - t647 ^ 2) * MDP(20) + (t417 * t647 + t418 * t475 - t621 - t760 + t766) * MDP(21) + (t641 - t747) * MDP(27) + (-t746 + t772) * MDP(28); -t454 ^ 2 * MDP(23) + (t412 + t779) * MDP(24) - t668 * MDP(25) + t439 * MDP(26) + (-t415 * t703 - t421 * t704 + t401 * t771 - g(1) * t624 - g(2) * t765 - g(3) * (-t531 * t613 - t687) + t669) * MDP(27) + (t400 * t771 + t414 * t454 + g(2) * t715 + (qJD(5) * t415 - t402 - t644) * t613 + (-g(2) * t741 - qJD(5) * t421 - t396 + t625) * t617) * MDP(28) + (t454 * MDP(22) + (-qJD(5) + t771) * MDP(25) - t414 * MDP(27) + MDP(23) * t456) * t456;];
tau = t1;
