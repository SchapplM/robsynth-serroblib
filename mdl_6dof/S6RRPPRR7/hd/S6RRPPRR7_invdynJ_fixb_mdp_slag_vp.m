% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:21:19
% EndTime: 2019-03-09 09:21:33
% DurationCPUTime: 11.37s
% Computational Cost: add. (4739->667), mult. (11245->855), div. (0->0), fcn. (8359->10), ass. (0->279)
t586 = sin(pkin(6));
t595 = cos(qJ(2));
t707 = qJD(1) * qJD(2);
t684 = t595 * t707;
t591 = sin(qJ(2));
t704 = qJDD(1) * t591;
t623 = t684 + t704;
t780 = t586 * t623;
t723 = qJD(1) * t586;
t691 = t591 * t723;
t587 = cos(pkin(6));
t713 = t587 * qJD(1);
t700 = pkin(1) * t713;
t730 = pkin(8) * t691 - t595 * t700;
t778 = qJD(3) + t730;
t783 = -qJ(4) * t691 + t778;
t532 = qJD(5) + t691;
t593 = cos(qJ(6));
t563 = qJD(2) + t713;
t590 = sin(qJ(5));
t594 = cos(qJ(5));
t722 = qJD(1) * t595;
t690 = t586 * t722;
t491 = t590 * t563 + t594 * t690;
t589 = sin(qJ(6));
t757 = t491 * t589;
t453 = -t593 * t532 - t757;
t665 = t590 * t690;
t490 = -t563 * t594 + t665;
t708 = qJD(6) - t490;
t782 = t453 * t708;
t641 = t491 * t593 - t532 * t589;
t781 = t641 * t708;
t685 = t591 * t707;
t663 = t586 * t685;
t703 = qJDD(1) * t595;
t682 = t586 * t703;
t779 = t663 - t682;
t596 = cos(qJ(1));
t735 = t595 * t596;
t592 = sin(qJ(1));
t742 = t591 * t592;
t507 = -t587 * t735 + t742;
t746 = t586 * t596;
t471 = t507 * t594 + t590 * t746;
t739 = t592 * t595;
t740 = t591 * t596;
t508 = t587 * t740 + t739;
t777 = t471 * t589 - t508 * t593;
t776 = t471 * t593 + t508 * t589;
t775 = t641 * t532;
t705 = qJDD(1) * t587;
t561 = qJDD(2) + t705;
t597 = -pkin(2) - pkin(3);
t698 = t561 * t597;
t579 = t586 * pkin(1);
t747 = t586 * t595;
t749 = t586 * t591;
t729 = pkin(2) * t747 + qJ(3) * t749;
t493 = -t579 - t729;
t567 = pkin(3) * t747;
t477 = t567 - t493;
t661 = pkin(4) * t591 + pkin(9) * t595;
t634 = t661 * t586;
t450 = t634 + t477;
t565 = pkin(8) * t749;
t653 = -qJ(4) * t749 + t565;
t692 = -pkin(1) * t595 - pkin(2);
t672 = -pkin(3) + t692;
t457 = (-pkin(9) + t672) * t587 + t653;
t774 = t590 * t450 + t594 * t457;
t509 = t587 * t739 + t740;
t773 = g(1) * t509 + g(2) * t507 - g(3) * t747;
t495 = qJDD(5) + t780;
t583 = -pkin(9) + t597;
t435 = t563 * t583 + t783;
t483 = -pkin(1) * t723 - pkin(2) * t690 - qJ(3) * t691;
t461 = pkin(3) * t690 + qJD(4) - t483;
t440 = qJD(1) * t634 + t461;
t411 = t435 * t594 + t440 * t590;
t626 = pkin(4) * t595 + t583 * t591;
t612 = qJD(2) * t626;
t559 = qJD(3) * t749;
t650 = pkin(2) * t682 + qJ(3) * t780 + qJD(1) * t559 + qJDD(1) * t579;
t615 = pkin(3) * t682 + qJDD(4) + t650;
t412 = (qJD(1) * t612 + qJDD(1) * t661) * t586 + t615;
t699 = pkin(1) * qJD(2) * t587;
t669 = qJD(1) * t699;
t697 = pkin(1) * t705;
t666 = pkin(8) * t780 + t591 * t669 - t595 * t697;
t639 = -qJDD(3) - t666;
t706 = qJD(1) * qJD(4);
t599 = (-qJ(4) * t623 - t591 * t706) * t586 - t639;
t415 = t561 * t583 + t599;
t604 = -qJD(5) * t411 + t594 * t412 - t415 * t590;
t397 = -pkin(5) * t495 - t604;
t748 = t586 * t592;
t473 = -t509 * t590 - t594 * t748;
t505 = -t587 * t594 + t590 * t747;
t736 = t594 * t596;
t756 = t507 * t590;
t619 = g(1) * t473 + g(2) * (t586 * t736 - t756) + g(3) * t505;
t772 = t708 * (-pkin(5) * t491 + pkin(10) * t708) + t397 + t619;
t716 = qJD(5) * t594;
t717 = qJD(5) * t590;
t431 = -t561 * t594 + t563 * t717 - t590 * t663 + (qJD(1) * t716 + qJDD(1) * t590) * t747;
t429 = -qJDD(6) + t431;
t588 = qJ(3) + pkin(4);
t520 = pkin(5) * t594 + pkin(10) * t590 + t588;
t771 = t708 * (-t783 + t532 * (pkin(5) * t590 - pkin(10) * t594)) + t520 * t429;
t558 = t595 * t699;
t576 = t587 * qJD(3);
t721 = qJD(2) * t591;
t689 = t586 * t721;
t451 = t586 * (pkin(8) * t721 + qJD(4) * t595) - qJ(4) * t689 - t558 - t576;
t533 = t561 * qJ(3);
t534 = t563 * qJD(3);
t667 = pkin(8) * t779 - t591 * t697 - t595 * t669;
t432 = t533 + t534 - t667;
t523 = qJ(4) * t663;
t423 = (qJ(4) * qJDD(1) + t706) * t747 - t432 - t523;
t770 = (qJD(2) - t563) * t722 + t704;
t720 = qJD(2) * t595;
t688 = t586 * t720;
t731 = qJ(3) * t688 + t559;
t437 = t586 * t612 + t731;
t767 = pkin(1) * t591;
t571 = t587 * t767;
t460 = -qJD(4) * t749 + (t571 + (pkin(8) - qJ(4)) * t747) * qJD(2);
t769 = -qJD(5) * t774 + t437 * t594 - t460 * t590;
t768 = 0.2e1 * t533;
t766 = pkin(2) * t561;
t765 = g(3) * t591;
t764 = qJ(3) * t595;
t430 = qJD(5) * t665 - t590 * t561 - t563 * t716 + t594 * t779;
t714 = qJD(6) * t593;
t694 = t593 * t430 + t589 * t495 + t532 * t714;
t715 = qJD(6) * t589;
t404 = t491 * t715 + t694;
t763 = t404 * t589;
t762 = t641 * t593;
t761 = t708 * t589;
t760 = t708 * t593;
t759 = t490 * t532;
t758 = t491 * t532;
t627 = t532 * t590;
t752 = t532 * t594;
t751 = t561 * t587;
t582 = t586 ^ 2;
t750 = t582 * qJD(1) ^ 2;
t745 = t589 * t429;
t744 = t590 * t404;
t743 = t590 * t495;
t741 = t591 * t594;
t738 = t593 * t429;
t737 = t594 * t495;
t540 = qJ(3) * t690;
t449 = t626 * t723 + t540;
t499 = pkin(8) * t690 + t591 * t700;
t481 = -qJ(4) * t690 + t499;
t732 = t590 * t449 + t594 * t481;
t728 = pkin(8) * t747 + t571;
t584 = t591 ^ 2;
t585 = t595 ^ 2;
t727 = t584 - t585;
t726 = MDP(12) * t586;
t725 = MDP(17) * t586;
t724 = qJ(3) * qJD(2);
t719 = qJD(5) * t532;
t718 = qJD(5) * t583;
t709 = -qJD(4) - t461;
t702 = 0.2e1 * t582;
t701 = g(3) * t749;
t696 = t595 * t750;
t695 = t594 * t747;
t492 = t587 * qJ(3) + t728;
t686 = qJ(4) * t704;
t680 = -t507 * pkin(2) + qJ(3) * t508;
t510 = -t587 * t742 + t735;
t679 = -t509 * pkin(2) + qJ(3) * t510;
t621 = t590 * t412 + t594 * t415 - t435 * t717 + t440 * t716;
t396 = pkin(10) * t495 + t621;
t416 = pkin(4) * t561 - t423;
t401 = -pkin(5) * t431 - pkin(10) * t430 + t416;
t678 = -t589 * t396 + t593 * t401;
t676 = t430 * t589 - t593 * t495;
t675 = t709 * t591;
t674 = qJD(1) * t493 + t483;
t673 = t563 + t713;
t671 = t561 + t705;
t670 = t597 * t749;
t668 = t591 * t696;
t658 = g(1) * t507 - g(2) * t509;
t657 = -g(1) * t510 - g(2) * t508;
t656 = -g(1) * t508 + g(2) * t510;
t655 = g(1) * t596 + g(2) * t592;
t654 = g(1) * t592 - g(2) * t596;
t484 = t589 * t594 * t691 - t593 * t690;
t652 = -t589 * t716 - t484;
t485 = (t589 * t595 + t593 * t741) * t723;
t651 = t593 * t716 + t485;
t649 = qJD(2) * t670;
t648 = t593 * t396 + t589 * t401;
t407 = pkin(10) * t532 + t411;
t537 = t563 * qJ(3);
t458 = -t481 - t537;
t448 = pkin(4) * t563 - t458;
t419 = -pkin(5) * t490 + pkin(10) * t491 + t448;
t399 = t407 * t593 + t419 * t589;
t647 = t407 * t589 - t419 * t593;
t418 = pkin(10) * t749 + t774;
t476 = qJ(4) * t747 - t492;
t463 = t587 * pkin(4) - t476;
t506 = t587 * t590 + t695;
t427 = -pkin(5) * t505 + pkin(10) * t506 + t463;
t646 = t418 * t593 + t427 * t589;
t645 = -t418 * t589 + t427 * t593;
t410 = -t435 * t590 + t440 * t594;
t643 = t450 * t594 - t457 * t590;
t638 = qJD(2) * (qJD(1) * t477 + t461);
t637 = -pkin(8) * t689 + t558;
t636 = t596 * pkin(1) + t510 * pkin(2) + pkin(8) * t748 + qJ(3) * t509;
t632 = t506 * t589 + t593 * t749;
t469 = -t506 * t593 + t589 * t749;
t631 = t708 * t714 - t745;
t630 = t708 * t715 + t738;
t628 = t532 * t453;
t424 = qJD(1) * t649 + t615;
t459 = t649 + t731;
t625 = qJD(1) * t459 + qJDD(1) * t477 + t424;
t436 = pkin(2) * t663 - t650;
t482 = pkin(2) * t689 - t731;
t624 = -qJD(1) * t482 - qJDD(1) * t493 - t436;
t622 = -t685 + t703;
t620 = t590 * t437 + t450 * t716 - t457 * t717 + t594 * t460;
t618 = -pkin(1) * t592 - t508 * pkin(2) + pkin(8) * t746 - qJ(3) * t507;
t500 = t728 * qJD(2);
t616 = -t500 * t563 - t656;
t614 = -t657 + t701;
t611 = -t590 * t715 + t651;
t610 = t657 - t667;
t609 = t416 - t614;
t406 = -pkin(5) * t532 - t410;
t608 = pkin(10) * t429 + (t406 + t410) * t708;
t607 = -t448 * t532 - t583 * t495;
t606 = -t666 + t773;
t605 = -qJDD(3) + t606;
t603 = t563 * t730 + t610;
t602 = qJD(6) * t583 * t708 + t614;
t601 = t499 * t563 + t606;
t600 = -(-pkin(10) * t690 + qJD(6) * t520 - t732) * t708 + t773;
t511 = t563 * t691;
t496 = pkin(2) * t691 - t540;
t494 = t587 * t692 + t565;
t487 = t576 + t637;
t480 = qJD(1) * t670 + t540;
t478 = t537 + t499;
t475 = -pkin(2) * t563 + t778;
t474 = t509 * t594 - t590 * t748;
t467 = qJD(5) * t505 + t594 * t689;
t466 = -qJD(5) * t695 - t587 * t717 + t590 * t689;
t464 = t587 * t672 + t653;
t462 = t770 * t586;
t446 = t563 * t597 + t783;
t442 = -t639 - t766;
t439 = t474 * t593 + t510 * t589;
t438 = -t474 * t589 + t510 * t593;
t426 = qJD(6) * t632 + t467 * t593 + t589 * t688;
t425 = qJD(6) * t469 + t467 * t589 - t593 * t688;
t422 = t698 + t599;
t420 = -pkin(5) * t690 - t449 * t594 + t481 * t590;
t417 = -pkin(5) * t749 - t643;
t413 = pkin(5) * t466 - pkin(10) * t467 - t451;
t405 = -qJD(6) * t641 + t676;
t403 = -pkin(5) * t688 - t769;
t402 = pkin(10) * t688 + t620;
t395 = -qJD(6) * t399 + t678;
t394 = -t647 * qJD(6) + t648;
t1 = [(-t565 * t561 - t666 * t587 + (t595 * t751 + t622 * t702) * pkin(1) + t616) * MDP(9) + ((-qJD(2) * t446 + qJDD(1) * t476 + t423 + (-qJD(2) * t464 + t451) * qJD(1)) * t595 + (-qJD(2) * t458 - qJDD(1) * t464 - t422 + (-qJD(2) * t476 - t460) * qJD(1)) * t591 + t655) * t725 + ((qJD(2) * t475 + qJDD(1) * t492 + t432 + (qJD(2) * t494 + t487) * qJD(1)) * t595 + (-qJD(2) * t478 + qJDD(1) * t494 + t442 + (-qJD(2) * t492 + t500) * qJD(1)) * t591 - t655) * t726 + (-t466 * t532 + t495 * t505) * MDP(22) + (t467 * t532 - t495 * t506) * MDP(21) + (-pkin(1) * t623 * t702 - t561 * t728 - t563 * t637 + t587 * t667 - t658) * MDP(10) + (0.2e1 * (t591 * t703 - t707 * t727) * MDP(5) + (qJDD(1) * t584 + 0.2e1 * t591 * t684) * MDP(4)) * t582 + (-g(1) * t756 - g(2) * t473 - t416 * t506 + t463 * t430 + t448 * t467 + t451 * t491 - t495 * t774 - t620 * t532) * MDP(25) + (t432 * t587 + t487 * t563 + t492 * t561 + t658) * MDP(13) + (-t442 * t587 - t494 * t561 + t616) * MDP(11) + (t422 * t587 + t460 * t563 + t464 * t561 + t656) * MDP(16) + (-t423 * t587 - t451 * t563 - t476 * t561 + t658) * MDP(15) + ((g(1) * t736 - t411 * t720 - t591 * t621) * MDP(25) + (t591 * t624 - t674 * t720) * MDP(13) + (t430 * t591 - t491 * t720) * MDP(21) + (t410 * t720 + t591 * t604) * MDP(24) + (t431 * t591 + t490 * t720) * MDP(22) + (t595 * t624 + t674 * t721) * MDP(11) + (t591 * t638 - t595 * t625) * MDP(16) + (t591 * t625 + t595 * t638) * MDP(15) + (t591 * t671 + t673 * t720) * MDP(6) + (t495 * t591 + t532 * t720) * MDP(23) + (t595 * t671 - t673 * t721) * MDP(7)) * t586 + (t422 * t464 + t446 * t460 + t423 * t476 + t458 * t451 + t424 * t477 + t461 * t459 - g(1) * (-pkin(3) * t508 - qJ(4) * t746 + t618) - g(2) * (pkin(3) * t510 - qJ(4) * t748 + t636)) * MDP(18) + qJDD(1) * MDP(1) + (-g(1) * t618 - g(2) * t636 + t432 * t492 + t436 * t493 + t442 * t494 + t475 * t500 + t478 * t487 + t483 * t482) * MDP(14) + t654 * MDP(2) + t655 * MDP(3) + ((-qJD(6) * t646 - t402 * t589 + t413 * t593) * t708 - t645 * t429 - t395 * t505 - t647 * t466 + t403 * t453 + t417 * t405 - t397 * t632 + t406 * t425 + g(1) * t776 - g(2) * t439) * MDP(31) + (-(qJD(6) * t645 + t402 * t593 + t413 * t589) * t708 + t646 * t429 + t394 * t505 - t399 * t466 - t403 * t641 + t417 * t404 + t397 * t469 + t406 * t426 - g(1) * t777 - g(2) * t438) * MDP(32) + (-t430 * t506 - t467 * t491) * MDP(19) + (t430 * t505 - t431 * t506 + t466 * t491 + t467 * t490) * MDP(20) + MDP(8) * t751 + (t404 * t469 - t426 * t641) * MDP(26) + (t404 * t632 - t405 * t469 + t425 * t641 - t426 * t453) * MDP(27) + (g(1) * t471 - g(2) * t474 - t416 * t505 - t463 * t431 + t448 * t466 + t451 * t490 + t643 * t495 + t769 * t532) * MDP(24) + (t429 * t505 + t466 * t708) * MDP(30) + (-t404 * t505 + t426 * t708 - t429 * t469 - t466 * t641) * MDP(28) + (t405 * t505 - t425 * t708 - t429 * t632 - t453 * t466) * MDP(29); (t422 * t597 - t423 * qJ(3) - t446 * t481 - t461 * t480 - g(1) * (-pkin(3) * t509 + t679) - g(2) * (-pkin(3) * t507 + t680) - g(3) * (t567 + t729) - t783 * t458) * MDP(18) + (t588 * t430 + t732 * t532 + t411 * t690 - t783 * t491 + t607 * t594 + (t532 * t718 - t609) * t590) * MDP(25) + (-t410 * t690 - t588 * t431 - t783 * t490 + (t481 * t532 + t607) * t590 + ((-t449 - t718) * t532 + t609) * t594) * MDP(24) + ((-t591 * t597 - t764) * qJDD(1) + ((t458 + t481 + t724) * t591 + (-qJD(2) * t597 + t446 - t783) * t595) * qJD(1)) * t725 + (t768 + t523 + t534 + t783 * t563 + (-qJ(4) * t703 - t765 + (-t480 * t591 + t595 * t709) * qJD(1)) * t586 + t610) * MDP(15) + (-t430 * t590 + t491 * t752) * MDP(19) + (t532 * t717 - t737 + (-t490 * t595 + t591 * t627) * t723) * MDP(22) + (-t442 * pkin(2) - g(1) * t679 - g(2) * t680 - g(3) * t729 + t432 * qJ(3) - t475 * t499 + t478 * t778 - t483 * t496) * MDP(14) + ((-pkin(2) * t591 + t764) * qJDD(1) + ((t478 - t499 - t724) * t591 + (-pkin(2) * qJD(2) - t475 + t778) * t595) * qJD(1)) * t726 - MDP(4) * t668 + ((-t430 - t759) * t594 + (-t431 - t758) * t590) * MDP(20) + (pkin(1) * t696 - t603 + t701) * MDP(10) + (t586 * t622 + t511) * MDP(7) + (-t532 * t716 - t743 + (t491 * t595 - t532 * t741) * t723) * MDP(21) - t532 * MDP(23) * t690 + (t768 + 0.2e1 * t534 + (-t765 + (t483 * t595 + t496 * t591) * qJD(1)) * t586 + t603) * MDP(13) + t462 * MDP(6) + (0.2e1 * t766 - qJDD(3) + (-t483 * t591 + t496 * t595) * t723 + t601) * MDP(11) + (t750 * t767 + t601) * MDP(9) + (-t481 * t563 + 0.2e1 * t698 + (-t686 + ((-qJ(4) * qJD(2) + t480) * t595 + t675) * qJD(1)) * t586 - t605) * MDP(16) + t727 * MDP(5) * t750 + t561 * MDP(8) + (-t406 * t484 - t420 * t453 - t771 * t593 + t600 * t589 + (t583 * t745 + t395 + (-t406 * t589 + t453 * t583) * qJD(5) - t602 * t593) * t594 + (t647 * t691 - t406 * t714 - t397 * t589 + t583 * t405 + (t583 * t761 + t647) * qJD(5)) * t590) * MDP(31) + (-t593 * t744 + t611 * t641) * MDP(26) + (t453 * t485 - t641 * t484 + (t453 * t593 - t589 * t641) * t716 + (t763 + t405 * t593 + (-t453 * t589 - t762) * qJD(6)) * t590) * MDP(27) + (-t406 * t485 + t420 * t641 + t771 * t589 + t600 * t593 + (t583 * t738 - t394 + (-t406 * t593 - t583 * t641) * qJD(5) + t602 * t589) * t594 + (t399 * t691 + t406 * t715 - t397 * t593 + t583 * t404 + (t583 * t760 + t399) * qJD(5)) * t590) * MDP(32) + (-t405 * t594 - t652 * t708 + (t628 + t631) * t590) * MDP(29) + (-t429 * t594 - t627 * t708) * MDP(30) + (t404 * t594 - t651 * t708 + (t630 + t775) * t590) * MDP(28); t462 * MDP(12) + (-t478 * t563 + t483 * t691 - t605 - t766) * MDP(14) - t770 * t725 + (t458 * t563 + t698 + (-t686 + (-qJ(4) * t720 + t675) * qJD(1)) * t586 - t605) * MDP(18) + (-t532 ^ 2 * t594 + t490 * t563 - t743) * MDP(24) + (t491 * t563 + t532 * t627 - t737) * MDP(25) + (t590 * t405 - (t593 * t563 - t589 * t627) * t708 + (t628 - t631) * t594) * MDP(31) + (t744 - (-t589 * t563 - t593 * t627) * t708 + (t630 - t775) * t594) * MDP(32) + (MDP(13) + MDP(15)) * (-t563 ^ 2 - t584 * t750) + (-MDP(11) + MDP(16)) * (t561 + t668); t511 * MDP(16) + (g(3) * t587 + t615) * MDP(18) + (MDP(24) * t495 - MDP(25) * t719 - MDP(31) * t405 - MDP(32) * t404) * t594 + (-MDP(24) * t719 - t495 * MDP(25) + (qJD(5) * t453 + t745) * MDP(31) + (-qJD(5) * t641 + t738) * MDP(32)) * t590 + (-t584 - t585) * MDP(17) * t750 - ((t590 * t714 - t652) * MDP(31) + t611 * MDP(32)) * t708 + (t654 * MDP(18) + (t591 * MDP(15) - MDP(16) * t595) * qJDD(1) + (((qJD(2) + t563) * MDP(15) - t458 * MDP(18) - t490 * MDP(24) - t491 * MDP(25)) * t595 + (-MDP(25) * t752 + t446 * MDP(18) + (MDP(18) * t597 + MDP(16)) * qJD(2) + (-MDP(24) * t532 + MDP(31) * t453 - MDP(32) * t641) * t590) * t591) * qJD(1)) * t586; -t490 ^ 2 * MDP(20) + (t430 - t759) * MDP(21) + (t431 - t758) * MDP(22) + t495 * MDP(23) + (t411 * t532 + t604 - t619) * MDP(24) + (g(1) * t474 + g(2) * t471 - g(3) * t506 + t410 * t532 - t448 * t490 - t621) * MDP(25) + (-t708 * t762 + t763) * MDP(26) + ((t404 - t782) * t593 + (-t405 + t781) * t589) * MDP(27) + (t708 * t760 - t745) * MDP(28) + (-t708 * t761 - t738) * MDP(29) + (-pkin(5) * t405 - t411 * t453 + t608 * t589 - t593 * t772) * MDP(31) + (-pkin(5) * t404 + t411 * t641 + t589 * t772 + t608 * t593) * MDP(32) + (MDP(19) * t490 + t491 * MDP(20) + t448 * MDP(24) - MDP(28) * t641 - t453 * MDP(29) + MDP(30) * t708 - MDP(31) * t647 - MDP(32) * t399) * t491; -t641 * t453 * MDP(26) + (-t453 ^ 2 + t641 ^ 2) * MDP(27) + (t694 + t782) * MDP(28) + (-t676 - t781) * MDP(29) - t429 * MDP(30) + (-g(1) * t438 + g(2) * t777 - g(3) * t632 + t399 * t708 + t406 * t641 + t678) * MDP(31) + (g(1) * t439 + g(2) * t776 + g(3) * t469 + t406 * t453 - t647 * t708 - t648) * MDP(32) + (MDP(28) * t757 + MDP(29) * t641 - MDP(31) * t399 + MDP(32) * t647) * qJD(6);];
tau  = t1;
