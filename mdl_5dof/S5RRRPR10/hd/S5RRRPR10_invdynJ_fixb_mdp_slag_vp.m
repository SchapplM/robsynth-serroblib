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
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(26,1)}
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
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:30:09
% EndTime: 2019-12-31 21:30:25
% DurationCPUTime: 10.54s
% Computational Cost: add. (6611->537), mult. (16703->753), div. (0->0), fcn. (13288->14), ass. (0->254)
t590 = cos(qJ(2));
t668 = qJD(1) * qJD(2);
t645 = t590 * t668;
t586 = sin(qJ(2));
t666 = qJDD(1) * t586;
t746 = t645 + t666;
t580 = sin(pkin(5));
t677 = qJD(1) * t590;
t651 = t580 * t677;
t745 = qJD(3) - t651;
t623 = pkin(2) * t586 - pkin(8) * t590;
t679 = qJD(1) * t580;
t516 = t623 * t679;
t589 = cos(qJ(3));
t504 = t589 * t516;
t652 = t586 * t679;
t582 = cos(pkin(5));
t678 = qJD(1) * t582;
t663 = pkin(1) * t678;
t515 = -pkin(7) * t652 + t590 * t663;
t585 = sin(qJ(3));
t583 = -qJ(4) - pkin(8);
t642 = qJD(3) * t583;
t690 = t589 * t590;
t744 = -t504 - (pkin(3) * t586 - qJ(4) * t690) * t679 + t589 * t642 + (t515 - qJD(4)) * t585;
t627 = t585 * t651;
t683 = t589 * t515 + t585 * t516;
t743 = -qJ(4) * t627 - qJD(4) * t589 - t585 * t642 + t683;
t588 = cos(qJ(5));
t584 = sin(qJ(5));
t564 = qJD(2) + t678;
t628 = t585 * t652;
t498 = -t589 * t564 + t628;
t500 = t564 * t585 + t589 * t652;
t579 = sin(pkin(10));
t581 = cos(pkin(10));
t613 = -t498 * t579 + t581 * t500;
t713 = t613 * t584;
t431 = -t588 * t745 + t713;
t637 = -t581 * t498 - t500 * t579;
t732 = qJD(5) - t637;
t742 = t431 * t732;
t433 = t584 * t745 + t588 * t613;
t741 = t433 * t732;
t531 = t579 * t585 - t581 * t589;
t740 = t745 * t531;
t532 = t579 * t589 + t581 * t585;
t682 = t745 * t532;
t739 = t746 * t580;
t667 = qJDD(1) * t582;
t563 = qJDD(2) + t667;
t673 = qJD(3) * t589;
t444 = -qJD(3) * t628 + t585 * t563 + t564 * t673 + t739 * t589;
t675 = qJD(2) * t590;
t648 = t585 * t675;
t674 = qJD(3) * t585;
t445 = -t589 * t563 + t564 * t674 + t580 * (qJD(1) * (t586 * t673 + t648) + t585 * t666);
t417 = -t444 * t579 - t581 * t445;
t416 = qJDD(5) - t417;
t634 = t732 * t588;
t738 = -t416 * t584 - t732 * t634;
t587 = sin(qJ(1));
t692 = t587 * t590;
t591 = cos(qJ(1));
t693 = t586 * t591;
t527 = t582 * t693 + t692;
t576 = qJ(3) + pkin(10);
t573 = sin(t576);
t574 = cos(t576);
t697 = t580 * t591;
t473 = -t527 * t574 + t573 * t697;
t688 = t590 * t591;
t694 = t586 * t587;
t526 = -t582 * t688 + t694;
t737 = t473 * t584 + t526 * t588;
t736 = t473 * t588 - t526 * t584;
t575 = t580 ^ 2;
t664 = 0.2e1 * t575;
t687 = t743 * t579 + t744 * t581;
t685 = t744 * t579 - t743 * t581;
t698 = t580 * t590;
t727 = pkin(1) * t586;
t681 = pkin(7) * t698 + t582 * t727;
t511 = pkin(8) * t582 + t681;
t612 = -pkin(2) * t590 - pkin(8) * t586 - pkin(1);
t512 = t612 * t580;
t684 = t589 * t511 + t585 * t512;
t518 = pkin(7) * t651 + t586 * t663;
t734 = -t518 + (-t627 + t674) * pkin(3);
t733 = g(1) * t591 + g(2) * t587;
t731 = pkin(3) * t445 + qJDD(4);
t529 = -t582 * t694 + t688;
t699 = t580 * t589;
t479 = -t529 * t585 + t587 * t699;
t701 = t580 * t586;
t524 = -t582 * t589 + t585 * t701;
t689 = t589 * t591;
t708 = t527 * t585;
t730 = g(3) * t524 - g(2) * (-t580 * t689 - t708) - g(1) * t479;
t665 = qJDD(1) * t590;
t562 = t580 * t665;
t646 = t586 * t668;
t626 = t580 * t646;
t514 = qJDD(3) - t562 + t626;
t488 = pkin(8) * t564 + t518;
t492 = qJD(1) * t512;
t443 = t488 * t589 + t492 * t585;
t662 = pkin(1) * qJD(2) * t582;
t630 = qJD(1) * t662;
t660 = pkin(1) * t667;
t654 = -pkin(7) * t562 - t586 * t660 - t590 * t630;
t599 = -pkin(7) * t626 - t654;
t460 = pkin(8) * t563 + t599;
t610 = t623 * qJD(2);
t462 = (qJD(1) * t610 + qJDD(1) * t612) * t580;
t595 = -qJD(3) * t443 - t585 * t460 + t589 * t462;
t386 = pkin(3) * t514 - qJ(4) * t444 - qJD(4) * t500 + t595;
t609 = -t589 * t460 - t585 * t462 + t488 * t674 - t492 * t673;
t388 = -qJ(4) * t445 - qJD(4) * t498 - t609;
t375 = t386 * t581 - t388 * t579;
t373 = -pkin(4) * t514 - t375;
t570 = pkin(3) * t579 + pkin(9);
t700 = t580 * t587;
t729 = t732 * (pkin(3) * t500 + pkin(4) * t613 - pkin(9) * t637 + qJD(5) * t570) + g(1) * (-t529 * t573 + t574 * t700) + g(2) * (-t527 * t573 - t574 * t697) + g(3) * (-t573 * t701 + t574 * t582) + t373;
t726 = pkin(2) * t563;
t720 = g(3) * t580;
t719 = MDP(6) * t580;
t418 = t444 * t581 - t445 * t579;
t671 = qJD(5) * t588;
t656 = t588 * t418 + t584 * t514 + t671 * t745;
t672 = qJD(5) * t584;
t389 = -t613 * t672 + t656;
t718 = t389 * t584;
t428 = -qJ(4) * t498 + t443;
t716 = t428 * t579;
t715 = t431 * t613;
t714 = t433 * t613;
t712 = t498 * t745;
t711 = t500 * t745;
t707 = t532 * t584;
t706 = t532 * t588;
t705 = t563 * MDP(8);
t704 = t574 * t584;
t703 = t574 * t588;
t702 = t575 * qJD(1) ^ 2;
t425 = t581 * t428;
t696 = t582 * t590;
t695 = t584 * t590;
t412 = t588 * t416;
t691 = t588 * t590;
t376 = t579 * t386 + t581 * t388;
t649 = t580 * t675;
t478 = -qJD(3) * t524 + t589 * t649;
t525 = t582 * t585 + t586 * t699;
t517 = t580 * t610;
t565 = pkin(7) * t701;
t519 = (t696 * pkin(1) - t565) * qJD(2);
t596 = -qJD(3) * t684 + t589 * t517 - t519 * t585;
t676 = qJD(2) * t586;
t650 = t580 * t676;
t402 = pkin(3) * t650 - qJ(4) * t478 - qJD(4) * t525 + t596;
t477 = qJD(3) * t525 + t580 * t648;
t608 = -t511 * t674 + t512 * t673 + t585 * t517 + t589 * t519;
t406 = -qJ(4) * t477 - qJD(4) * t524 + t608;
t384 = t579 * t402 + t581 * t406;
t442 = -t488 * t585 + t589 * t492;
t427 = -qJ(4) * t500 + t442;
t423 = pkin(3) * t745 + t427;
t395 = t579 * t423 + t425;
t636 = -t511 * t585 + t589 * t512;
t430 = -pkin(3) * t698 - qJ(4) * t525 + t636;
t438 = -qJ(4) * t524 + t684;
t405 = t579 * t430 + t581 * t438;
t686 = pkin(4) * t652 - t687;
t520 = pkin(7) * t649 + t586 * t662;
t577 = t586 ^ 2;
t680 = -t590 ^ 2 + t577;
t670 = qJD(2) - t564;
t659 = t590 * t702;
t658 = t580 * t695;
t657 = t580 * t691;
t572 = pkin(3) * t589 + pkin(2);
t647 = t583 * t585;
t374 = pkin(9) * t514 + t376;
t629 = t739 * pkin(7) + t586 * t630 - t590 * t660;
t461 = t629 - t726;
t420 = t461 + t731;
t380 = -pkin(4) * t417 - pkin(9) * t418 + t420;
t641 = -t584 * t374 + t588 * t380;
t640 = t418 * t584 - t588 * t514;
t639 = t740 * t584 - t588 * t652;
t638 = t584 * t652 + t740 * t588;
t635 = t527 * t589 - t585 * t697;
t633 = t564 + t678;
t631 = t563 + t667;
t624 = pkin(3) * t477 + t520;
t528 = t582 * t692 + t693;
t622 = -g(1) * t526 + g(2) * t528;
t621 = g(1) * t529 + g(2) * t527;
t470 = pkin(4) * t531 - pkin(9) * t532 - t572;
t619 = pkin(9) * t652 - qJD(5) * t470 - t685;
t549 = t583 * t589;
t486 = -t581 * t549 + t579 * t647;
t618 = -t682 * pkin(4) - t740 * pkin(9) + qJD(5) * t486 - t734;
t617 = t588 * t374 + t584 * t380;
t392 = pkin(9) * t745 + t395;
t487 = -pkin(2) * t564 - t515;
t452 = pkin(3) * t498 + qJD(4) + t487;
t398 = -pkin(4) * t637 - pkin(9) * t613 + t452;
t379 = t392 * t588 + t398 * t584;
t616 = t392 * t584 - t398 * t588;
t401 = -pkin(9) * t698 + t405;
t468 = t581 * t524 + t525 * t579;
t469 = -t524 * t579 + t525 * t581;
t510 = t565 + (-pkin(1) * t590 - pkin(2)) * t582;
t600 = pkin(3) * t524 + t510;
t419 = pkin(4) * t468 - pkin(9) * t469 + t600;
t615 = t401 * t588 + t419 * t584;
t614 = -t401 * t584 + t419 * t588;
t383 = t402 * t581 - t406 * t579;
t394 = t423 * t581 - t716;
t404 = t430 * t581 - t438 * t579;
t611 = t412 + (t584 * t637 - t672) * t732;
t449 = t469 * t584 + t657;
t607 = t532 * t671 - t639;
t606 = -t532 * t672 - t638;
t602 = -g(1) * t528 - g(2) * t526 + g(3) * t698;
t598 = -pkin(8) * t514 + t487 * t745;
t391 = -pkin(4) * t745 - t394;
t397 = t427 * t581 - t716;
t597 = -t570 * t416 + (t391 + t397) * t732;
t594 = -t602 - t629;
t593 = -pkin(8) * qJD(3) * t745 - t461 - t602;
t571 = -pkin(3) * t581 - pkin(4);
t509 = t573 * t582 + t574 * t701;
t485 = -t549 * t579 - t581 * t647;
t480 = t529 * t589 + t585 * t700;
t475 = t529 * t574 + t573 * t700;
t450 = t469 * t588 - t658;
t447 = t475 * t588 + t528 * t584;
t446 = -t475 * t584 + t528 * t588;
t437 = -t477 * t579 + t478 * t581;
t436 = t581 * t477 + t478 * t579;
t408 = -qJD(5) * t658 + t437 * t584 + t469 * t671 - t588 * t650;
t407 = -qJD(5) * t449 + t437 * t588 + t584 * t650;
t400 = pkin(4) * t698 - t404;
t396 = t427 * t579 + t425;
t393 = pkin(4) * t436 - pkin(9) * t437 + t624;
t390 = qJD(5) * t433 + t640;
t382 = pkin(9) * t650 + t384;
t381 = -pkin(4) * t650 - t383;
t372 = -t379 * qJD(5) + t641;
t371 = -qJD(5) * t616 + t617;
t1 = [(t416 * t468 + t436 * t732) * MDP(24) + (-t390 * t468 - t408 * t732 - t416 * t449 - t431 * t436) * MDP(23) + (t389 * t468 + t407 * t732 + t416 * t450 + t433 * t436) * MDP(22) + (-t520 * t564 - t565 * t563 - t629 * t582 + g(1) * t527 - g(2) * t529 + (t563 * t696 + (-t646 + t665) * t664) * pkin(1)) * MDP(9) + (t376 * t405 + t395 * t384 + t375 * t404 + t394 * t383 + t420 * t600 + t452 * t624 - g(1) * (-pkin(1) * t587 + t526 * t583 - t527 * t572) - g(2) * (pkin(1) * t591 - t528 * t583 + t529 * t572)) * MDP(19) + (-g(1) * t708 - g(2) * t479 + t510 * t444 + t461 * t525 + t487 * t478 + t520 * t500 - t684 * t514 - t608 * t745) * MDP(17) + (g(1) * t635 - g(2) * t480 + t510 * t445 + t461 * t524 + t487 * t477 + t520 * t498 + t636 * t514 + t596 * t745) * MDP(16) + (t478 * t745 + t514 * t525) * MDP(13) + (g(1) * t587 - g(2) * t591) * MDP(2) + (qJDD(1) * t577 + 0.2e1 * t586 * t645) * t575 * MDP(4) + (-t444 * t524 - t445 * t525 - t477 * t500 - t478 * t498) * MDP(12) + (t444 * t525 + t478 * t500) * MDP(11) + (t586 * t665 - t668 * t680) * MDP(5) * t664 + t582 * t705 + (t586 * t631 + t633 * t675) * t719 + t733 * MDP(3) + (-t477 * t745 - t514 * t524) * MDP(14) + ((-g(1) * t689 - t443 * t676 - t590 * t609) * MDP(17) + (t445 * t590 - t498 * t676) * MDP(14) + (t442 * t676 - t590 * t595) * MDP(16) + (-t444 * t590 + t500 * t676) * MDP(13) + (t590 * t631 - t633 * t676) * MDP(7) + (-t514 * t590 + t676 * t745) * MDP(15) - t733 * (pkin(3) * t585 + pkin(7)) * MDP(19)) * t580 + (-t746 * pkin(1) * t664 - t519 * t564 - t681 * t563 - t599 * t582 + t622) * MDP(10) + qJDD(1) * MDP(1) + (-t375 * t469 - t376 * t468 - t383 * t613 + t384 * t637 - t394 * t437 - t395 * t436 - t404 * t418 + t405 * t417 - t622) * MDP(18) + (t389 * t450 + t407 * t433) * MDP(20) + (-t389 * t449 - t390 * t450 - t407 * t431 - t408 * t433) * MDP(21) + ((-qJD(5) * t615 - t382 * t584 + t393 * t588) * t732 + t614 * t416 + t372 * t468 - t616 * t436 + t381 * t431 + t400 * t390 + t373 * t449 + t391 * t408 - g(1) * t736 - g(2) * t447) * MDP(25) + (-(qJD(5) * t614 + t382 * t588 + t393 * t584) * t732 - t615 * t416 - t371 * t468 - t379 * t436 + t381 * t433 + t400 * t389 + t373 * t450 + t391 * t407 + g(1) * t737 - g(2) * t446) * MDP(26); -t586 * MDP(4) * t659 + t680 * MDP(5) * t702 + (t670 * t677 + t666) * t719 + (-t652 * t670 + t562) * MDP(7) + t705 + (t518 * t564 + t702 * t727 + t594) * MDP(9) + (pkin(1) * t659 + t515 * t564 + (pkin(7) * t668 + g(3)) * t701 + t621 + t654) * MDP(10) + (t444 * t585 + t589 * t711) * MDP(11) + ((t444 - t712) * t589 + (-t445 - t711) * t585) * MDP(12) + (t745 * t673 + t514 * t585 + (-t500 * t586 - t690 * t745) * t679) * MDP(13) + (-t745 * t674 + t514 * t589 + (t585 * t590 * t745 + t498 * t586) * t679) * MDP(14) - t745 * MDP(15) * t652 + (-t442 * t652 - pkin(2) * t445 - t518 * t498 - t504 * t745 + (t515 * t745 + t598) * t585 + t593 * t589) * MDP(16) + (-pkin(2) * t444 + t443 * t652 - t518 * t500 - t585 * t593 + t589 * t598 + t683 * t745) * MDP(17) + (-g(3) * t701 - t375 * t532 - t376 * t531 + t740 * t394 - t682 * t395 + t417 * t486 + t418 * t485 - t687 * t613 + t685 * t637 - t621) * MDP(18) + (t376 * t486 - t375 * t485 - t420 * t572 - g(1) * (-t528 * t572 - t529 * t583) - g(2) * (-t526 * t572 - t527 * t583) - (t572 * t590 - t583 * t586) * t720 + t734 * t452 + t685 * t395 + t687 * t394) * MDP(19) + (t389 * t706 + t433 * t606) * MDP(20) + (t639 * t433 + t638 * t431 + (-t718 - t390 * t588 + (t431 * t584 - t433 * t588) * qJD(5)) * t532) * MDP(21) + (t389 * t531 + t412 * t532 + t433 * t682 + t606 * t732) * MDP(22) + (-t390 * t531 - t416 * t707 - t431 * t682 - t607 * t732) * MDP(23) + (t416 * t531 + t682 * t732) * MDP(24) + ((t470 * t588 - t486 * t584) * t416 + t372 * t531 + t485 * t390 + t373 * t707 - g(1) * (-t528 * t703 + t529 * t584) - g(2) * (-t526 * t703 + t527 * t584) - (t574 * t691 + t584 * t586) * t720 + (t584 * t619 - t588 * t618) * t732 + t686 * t431 - t682 * t616 + t607 * t391) * MDP(25) + (-(t470 * t584 + t486 * t588) * t416 - t371 * t531 + t485 * t389 + t373 * t706 - g(1) * (t528 * t704 + t529 * t588) - g(2) * (t526 * t704 + t527 * t588) - (-t574 * t695 + t586 * t588) * t720 + (t584 * t618 + t588 * t619) * t732 + t686 * t433 - t682 * t379 + t606 * t391) * MDP(26); t500 * t498 * MDP(11) + (-t498 ^ 2 + t500 ^ 2) * MDP(12) + (t444 + t712) * MDP(13) + (-t445 + t711) * MDP(14) + t514 * MDP(15) + (t443 * t745 - t487 * t500 + t595 + t730) * MDP(16) + (g(1) * t480 + g(2) * t635 + g(3) * t525 + t442 * t745 + t487 * t498 + t609) * MDP(17) + ((t417 * t579 - t418 * t581) * pkin(3) + (t394 - t397) * t637 + (t395 - t396) * t613) * MDP(18) + (t394 * t396 - t395 * t397 + (t375 * t581 + t376 * t579 - t452 * t500 + t730) * pkin(3)) * MDP(19) + (t433 * t634 + t718) * MDP(20) + ((t389 - t742) * t588 + (-t390 - t741) * t584) * MDP(21) + (-t714 - t738) * MDP(22) + (t611 + t715) * MDP(23) - t732 * t613 * MDP(24) + (t571 * t390 - t396 * t431 + t597 * t584 - t588 * t729 + t613 * t616) * MDP(25) + (t379 * t613 + t571 * t389 - t396 * t433 + t584 * t729 + t597 * t588) * MDP(26); (-t613 ^ 2 - t637 ^ 2) * MDP(18) + (t394 * t613 - t395 * t637 - t594 - t726 + t731) * MDP(19) + (t611 - t715) * MDP(25) + (-t714 + t738) * MDP(26); t433 * t431 * MDP(20) + (-t431 ^ 2 + t433 ^ 2) * MDP(21) + (t656 + t742) * MDP(22) + (-t640 + t741) * MDP(23) + t416 * MDP(24) + (t379 * t732 - t391 * t433 - g(1) * t446 - g(2) * t737 - g(3) * (-t509 * t584 - t657) + t641) * MDP(25) + (-t616 * t732 + t391 * t431 + g(1) * t447 - g(2) * t736 - g(3) * (-t509 * t588 + t658) - t617) * MDP(26) + (-MDP(22) * t713 - MDP(23) * t433 - MDP(25) * t379 + MDP(26) * t616) * qJD(5);];
tau = t1;
