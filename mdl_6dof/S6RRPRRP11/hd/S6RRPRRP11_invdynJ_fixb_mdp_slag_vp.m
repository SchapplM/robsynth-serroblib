% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP11
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP11_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP11_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:48:53
% EndTime: 2019-03-09 12:49:04
% DurationCPUTime: 8.46s
% Computational Cost: add. (5697->592), mult. (11693->734), div. (0->0), fcn. (7493->10), ass. (0->269)
t589 = sin(qJ(2));
t690 = qJD(1) * t589;
t565 = pkin(2) * t690;
t593 = cos(qJ(2));
t726 = qJ(3) * t593;
t634 = pkin(8) * t589 - t726;
t477 = qJD(1) * t634 + t565;
t689 = qJD(1) * t593;
t561 = pkin(7) * t689;
t521 = pkin(3) * t689 + t561;
t592 = cos(qJ(4));
t500 = t592 * t521;
t588 = sin(qJ(4));
t717 = t588 * t589;
t627 = pkin(4) * t593 - pkin(9) * t717;
t681 = qJD(4) * t588;
t595 = -pkin(2) - pkin(8);
t728 = pkin(9) - t595;
t759 = qJD(1) * t627 - t477 * t588 - t728 * t681 + t500;
t527 = t728 * t592;
t663 = t592 * t690;
t697 = t592 * t477 + t588 * t521;
t756 = pkin(9) * t663 + qJD(4) * t527 + t697;
t674 = t588 * qJD(2);
t511 = -t592 * t689 - t674;
t672 = qJD(1) * qJD(2);
t659 = t589 * t672;
t670 = qJDD(1) * t593;
t444 = t511 * qJD(4) + t592 * qJDD(2) + (t659 - t670) * t588;
t660 = t588 * t689;
t685 = qJD(2) * t592;
t512 = -t660 + t685;
t587 = sin(qJ(5));
t591 = cos(qJ(5));
t679 = qJD(4) * t593;
t662 = t588 * t679;
t612 = t685 * t589 + t662;
t680 = qJD(4) * t592;
t667 = qJD(2) * t680 + t588 * qJDD(2) + t592 * t670;
t600 = qJD(1) * t612 - t667;
t676 = qJD(5) * t591;
t677 = qJD(5) * t587;
t397 = -t591 * t444 - t511 * t676 + t512 * t677 - t587 * t600;
t632 = t511 * t587 + t591 * t512;
t398 = qJD(5) * t632 + t444 * t587 - t591 * t600;
t450 = -t591 * t511 + t512 * t587;
t447 = t450 ^ 2;
t658 = t593 * t672;
t671 = qJDD(1) * t589;
t618 = t658 + t671;
t510 = qJDD(4) + t618;
t502 = qJDD(5) + t510;
t547 = qJD(4) + t690;
t535 = qJD(5) + t547;
t740 = t632 ^ 2;
t758 = t502 * MDP(26) + (t535 * t632 - t398) * MDP(25) + t450 * MDP(22) * t632 + (t450 * t535 - t397) * MDP(24) + (-t447 + t740) * MDP(23);
t739 = pkin(3) + pkin(7);
t757 = qJ(6) * t450;
t573 = t592 * pkin(4);
t664 = -pkin(3) - t573;
t686 = qJD(2) * t589;
t458 = (-pkin(7) + t664) * t686 - pkin(4) * t662;
t558 = pkin(7) * t670;
t581 = qJDD(2) * qJ(3);
t582 = qJD(2) * qJD(3);
t666 = t558 + t581 + t582;
t639 = pkin(3) * t670 + t666;
t410 = pkin(4) * t667 + t458 * qJD(1) + t639;
t385 = t398 * pkin(5) + qJDD(6) + t410;
t513 = t587 * t592 + t588 * t591;
t747 = qJD(4) + qJD(5);
t455 = t747 * t513;
t620 = t513 * t589;
t699 = qJD(1) * t620 + t455;
t710 = t591 * t592;
t719 = t587 * t588;
t698 = -t587 * t681 - t588 * t677 + t591 * t663 - t690 * t719 + t710 * t747;
t560 = pkin(7) * t690;
t755 = qJD(3) + t560;
t583 = qJD(2) * qJ(3);
t492 = t583 + t521;
t457 = -pkin(4) * t511 + t492;
t586 = qJ(4) + qJ(5);
t569 = sin(t586);
t570 = cos(t586);
t590 = sin(qJ(1));
t594 = cos(qJ(1));
t714 = t589 * t594;
t473 = t569 * t714 + t570 * t590;
t716 = t589 * t590;
t475 = -t569 * t716 + t570 * t594;
t579 = g(3) * t593;
t571 = t589 * qJ(3);
t656 = -pkin(1) - t571;
t614 = t593 * t595 + t656;
t465 = t614 * qJD(1);
t673 = pkin(3) * t690 + t755;
t471 = qJD(2) * t595 + t673;
t428 = t465 * t592 + t471 * t588;
t546 = pkin(2) * t659;
t683 = qJD(3) * t589;
t607 = qJD(2) * t634 - t683;
t425 = qJD(1) * t607 + qJDD(1) * t614 + t546;
t545 = pkin(7) * t658;
t557 = pkin(7) * t671;
t657 = qJDD(3) + t545 + t557;
t446 = pkin(3) * t618 + qJDD(2) * t595 + t657;
t650 = -t425 * t588 + t592 * t446;
t604 = -qJD(4) * t428 + t650;
t384 = pkin(4) * t510 - pkin(9) * t444 + t604;
t668 = t592 * t425 + t588 * t446 + t471 * t680;
t388 = pkin(9) * t600 - t465 * t681 + t668;
t427 = -t465 * t588 + t592 * t471;
t417 = -pkin(9) * t512 + t427;
t411 = pkin(4) * t547 + t417;
t418 = pkin(9) * t511 + t428;
t642 = -t587 * t384 - t591 * t388 - t411 * t676 + t418 * t677;
t754 = g(1) * t473 - g(2) * t475 + t450 * t457 - t569 * t579 + t642;
t752 = qJ(6) * t632;
t413 = pkin(5) * t450 + qJD(6) + t457;
t751 = t413 * t632;
t629 = -t710 + t719;
t478 = t629 * t593;
t531 = t739 * t589;
t516 = t592 * t531;
t576 = t593 * pkin(2);
t693 = t576 + t571;
t528 = -pkin(1) - t693;
t503 = -pkin(8) * t593 + t528;
t655 = pkin(9) * t593 - t503;
t435 = pkin(4) * t589 + t588 * t655 + t516;
t515 = t588 * t531;
t696 = t592 * t503 + t515;
t709 = t592 * t593;
t441 = -pkin(9) * t709 + t696;
t700 = t587 * t435 + t591 * t441;
t750 = t759 * t591;
t526 = t728 * t588;
t695 = -t591 * t526 - t587 * t527;
t749 = pkin(4) * t680 + t755;
t748 = -t526 * t677 + t527 * t676 + t759 * t587 + t756 * t591;
t637 = g(1) * t594 + g(2) * t590;
t472 = -t569 * t590 + t570 * t714;
t474 = t569 * t594 + t570 * t716;
t745 = -g(1) * t472 - g(2) * t474 + t570 * t579;
t416 = t591 * t418;
t394 = t411 * t587 + t416;
t653 = t591 * t384 - t587 * t388;
t606 = -t394 * qJD(5) + t653;
t744 = -t457 * t632 + t606 + t745;
t743 = t397 * t629 - t632 * t699;
t378 = pkin(5) * t502 + qJ(6) * t397 - qJD(6) * t632 + t606;
t379 = -qJ(6) * t398 - qJD(6) * t450 - t642;
t414 = t587 * t418;
t393 = t591 * t411 - t414;
t389 = t393 - t752;
t386 = pkin(5) * t535 + t389;
t390 = t394 - t757;
t665 = -g(1) * t714 - g(2) * t716 + t579;
t742 = -t378 * t629 + t379 * t513 - t386 * t699 + t390 * t698 + t665;
t738 = pkin(4) * t588;
t736 = g(1) * t590;
t732 = g(2) * t594;
t731 = g(3) * t589;
t727 = pkin(7) * qJDD(2);
t725 = qJDD(2) * pkin(2);
t723 = t444 * t588;
t722 = t444 * t592;
t721 = t510 * t588;
t580 = -qJ(6) - pkin(9) - pkin(8);
t720 = t580 * t593;
t718 = t588 * t511;
t715 = t589 * t592;
t597 = qJD(1) ^ 2;
t713 = t589 * t597;
t712 = t590 * t592;
t711 = t590 * t593;
t487 = t592 * t510;
t708 = t592 * t594;
t707 = t593 * t594;
t706 = t595 * t510;
t705 = -t389 + t386;
t704 = -qJ(6) * t698 - qJD(6) * t513 - t748;
t703 = -pkin(5) * t689 + qJ(6) * t699 - t695 * qJD(5) + qJD(6) * t629 + t587 * t756 - t750;
t702 = t591 * t417 - t414;
t635 = qJD(1) * t664;
t694 = -t589 * t635 + t749;
t524 = pkin(5) * t570 + t573;
t532 = t739 * t593;
t584 = t589 ^ 2;
t585 = t593 ^ 2;
t692 = t584 - t585;
t688 = qJD(2) * t511;
t687 = qJD(2) * t512;
t684 = qJD(2) * t593;
t682 = qJD(4) * t465;
t678 = qJD(4) * t595;
t675 = t492 * qJD(4);
t669 = t593 * t713;
t490 = pkin(4) * t709 + t532;
t661 = t592 * t679;
t553 = qJ(3) + t738;
t654 = -qJD(2) * pkin(2) + qJD(3);
t564 = pkin(2) * t686;
t462 = t564 + t607;
t522 = t739 * t684;
t647 = -t462 * t588 + t592 * t522;
t405 = t627 * qJD(2) + (t592 * t655 - t515) * qJD(4) + t647;
t617 = t592 * t462 - t503 * t681 + t588 * t522 + t531 * t680;
t407 = pkin(9) * t612 + t617;
t652 = t591 * t405 - t407 * t587;
t651 = -t417 * t587 - t416;
t649 = t591 * t435 - t441 * t587;
t646 = t526 * t587 - t591 * t527;
t645 = -qJD(1) * t532 - t492;
t644 = t547 + t690;
t641 = t594 * pkin(1) + pkin(2) * t707 + t590 * pkin(7) + qJ(3) * t714;
t640 = -t557 - t665;
t520 = t739 * t686;
t596 = qJD(2) ^ 2;
t638 = pkin(7) * t596 + t732;
t636 = -t502 * t629 - t535 * t699;
t631 = t512 * t592 + t718;
t525 = t560 + t654;
t530 = -t561 - t583;
t630 = t525 * t593 + t530 * t589;
t628 = t547 * t588;
t626 = t656 - t576;
t625 = -0.2e1 * pkin(1) * t672 - t727;
t493 = t626 * qJD(1);
t624 = t493 * t690 + qJDD(3) - t640;
t623 = -t547 * t680 - t721;
t621 = -qJ(3) * t684 - t683;
t619 = -t502 * t513 - t535 * t698;
t616 = t587 * t405 + t591 * t407 + t435 * t676 - t441 * t677;
t613 = 0.2e1 * qJDD(1) * pkin(1) - t638;
t610 = t727 + (-qJD(1) * t528 - t493) * qJD(2);
t608 = -t637 * t593 - t731;
t445 = qJD(1) * t621 + qJDD(1) * t626 + t546;
t482 = t564 + t621;
t603 = qJD(1) * t482 + qJDD(1) * t528 + t445 + t638;
t468 = pkin(7) * t659 - t666;
t476 = t657 - t725;
t601 = qJD(2) * t630 - t468 * t593 + t476 * t589;
t577 = t594 * pkin(7);
t556 = pkin(4) * t591 + pkin(5);
t551 = g(1) * t711;
t544 = qJ(3) * t707;
t542 = qJ(3) * t711;
t523 = pkin(5) * t569 + t738;
t518 = pkin(3) + t524;
t517 = -qJ(3) * t689 + t565;
t497 = -t588 * t716 + t708;
t496 = t588 * t594 + t589 * t712;
t495 = t588 * t714 + t712;
t494 = -t588 * t590 + t589 * t708;
t479 = t513 * t593;
t448 = -qJD(1) * t520 + t639;
t434 = -qJ(6) * t513 + t695;
t433 = qJ(6) * t629 + t646;
t420 = t455 * t593 - t629 * t686;
t419 = qJD(2) * t620 + t478 * t747;
t402 = qJ(6) * t478 + t700;
t401 = pkin(5) * t589 + qJ(6) * t479 + t649;
t392 = t702 - t752;
t391 = t651 + t757;
t381 = qJ(6) * t420 + qJD(6) * t478 + t616;
t380 = pkin(5) * t684 - qJ(6) * t419 - qJD(5) * t700 + qJD(6) * t479 + t652;
t1 = [(t589 * t610 + t593 * t603 - t551) * MDP(12) + (-t397 * t478 + t398 * t479 - t419 * t450 + t420 * t632) * MDP(23) + (t397 * t479 + t419 * t632) * MDP(22) + (-g(2) * t707 + t378 * t479 + t379 * t478 - t380 * t632 - t381 * t450 - t386 * t419 + t390 * t420 + t397 * t401 - t398 * t402 + t551) * MDP(29) + (g(1) * t474 - g(2) * t472 - t394 * t684 - t490 * t397 - t410 * t479 + t457 * t419 + t458 * t632 - t502 * t700 - t535 * t616 + t589 * t642) * MDP(28) + (t379 * t402 + t390 * t381 + t378 * t401 + t386 * t380 + t385 * (-pkin(5) * t478 + t490) - g(1) * (t518 * t594 + t577) - g(2) * (t523 * t714 - t580 * t707 + t641) + (-g(1) * (-t523 * t589 + t626 + t720) - g(2) * t518) * t590 + (-pkin(5) * t420 + t458) * t413) * MDP(30) + (qJDD(2) * t589 + t593 * t596) * MDP(6) + (qJDD(2) * t593 - t589 * t596) * MDP(7) + qJDD(1) * MDP(1) + (t652 * t535 + t649 * t502 + t653 * t589 + t393 * t684 + t458 * t450 + t490 * t398 - t410 * t478 - t457 * t420 - g(1) * t475 - g(2) * t473 + (-t394 * t589 - t535 * t700) * qJD(5)) * MDP(27) + (pkin(7) * t601 - g(1) * t577 - g(2) * t641 + t445 * t528 + t493 * t482 - t626 * t736) * MDP(14) + (t625 * t593 + (-t613 - t736) * t589) * MDP(10) + (t610 * t593 + (-t603 + t736) * t589) * MDP(13) + (-t732 + t736) * MDP(2) + (t631 * t686 + (-t588 * (t592 * t659 - t667) - t722 + (-t592 * t511 + (t512 - t660) * t588) * qJD(4)) * t593) * MDP(16) + (-t593 * t723 + (t589 * t674 - t661) * t512) * MDP(15) + ((t644 * t685 - t667) * t589 + (t644 * t681 - t487 + t688) * t593) * MDP(18) + 0.2e1 * (t589 * t670 - t672 * t692) * MDP(5) + (-t398 * t589 + t420 * t535 - t450 * t684 + t478 * t502) * MDP(25) + (t510 * t589 + t547 * t684) * MDP(19) + (t502 * t589 + t535 * t684) * MDP(26) + ((t547 * t674 + t444) * t589 + (t623 + t687) * t593) * MDP(17) + (qJDD(1) * t584 + 0.2e1 * t589 * t658) * MDP(4) + (t589 * t625 + t593 * t613 + t551) * MDP(9) + ((t584 + t585) * qJDD(1) * pkin(7) + t601 - t637) * MDP(11) + t637 * MDP(3) + (-t397 * t589 + t419 * t535 - t479 * t502 + t632 * t684) * MDP(24) + (t647 * t547 + (-t503 * t588 + t516) * t510 + t650 * t589 + t520 * t511 + t532 * t667 + t448 * t709 - g(1) * t497 - g(2) * t495 + (t427 * t593 + t645 * t715) * qJD(2) + ((-t465 * t589 - t503 * t547) * t592 + (-t471 * t589 - t531 * t547 + t593 * t645) * t588) * qJD(4)) * MDP(20) + (-t617 * t547 - t696 * t510 - t520 * t512 + t532 * t444 + g(1) * t496 - g(2) * t494 + ((qJD(2) * t492 + t682) * t588 - t668) * t589 + (-qJD(2) * t428 - t448 * t588 - t592 * t675) * t593) * MDP(21); (-t553 * t397 - t410 * t629 - t699 * t457 - t695 * t502 + t535 * t748 + t608 * t570 + t694 * t632) * MDP(28) + (t397 * t513 + t398 * t629 + t450 * t699 - t632 * t698) * MDP(23) + t636 * MDP(24) + t743 * MDP(22) + (t379 * t434 + t378 * t433 + t385 * (pkin(5) * t513 + t553) - g(1) * (t523 * t707 + t544) - g(2) * (t523 * t711 + t542) - g(3) * (t693 - t720) + (pkin(5) * t698 + t749) * t413 + t704 * t390 + t703 * t386 + (-g(3) * t523 - t413 * t635 + t637 * (pkin(2) - t580)) * t589) * MDP(30) + (t646 * t502 + t553 * t398 + t410 * t513 + (t526 * t676 + (qJD(5) * t527 + t756) * t587 - t750) * t535 + t698 * t457 + t694 * t450 + t608 * t569) * MDP(27) + (-t517 * MDP(12) - t547 * MDP(19) - t427 * MDP(20) + t428 * MDP(21) - MDP(24) * t632 + t450 * MDP(25) - t535 * MDP(26) - t393 * MDP(27) + t394 * MDP(28)) * t689 + qJDD(2) * MDP(8) + (t397 * t433 - t398 * t434 - t450 * t704 - t632 * t703 - t742) * MDP(29) + t619 * MDP(25) + (t624 - 0.2e1 * t725) * MDP(12) + (t731 - t558 + (pkin(1) * t597 + t637) * t593) * MDP(10) + ((-pkin(2) * t589 + t726) * qJDD(1) + ((-t530 - t583) * t589 + (-t525 + t654) * t593) * qJD(1)) * MDP(11) + (-t592 * t667 - t723 - t631 * qJD(4) + (t588 * t661 + (-t718 + (-t512 + t685) * t592) * t589) * qJD(1)) * MDP(16) + (-t547 * t681 + t487 + (-t512 * t593 - t547 * t717) * qJD(1)) * MDP(17) + ((-t511 * t593 - t547 * t715) * qJD(1) + t623) * MDP(18) + (-t468 * qJ(3) - t530 * qJD(3) - t476 * pkin(2) - t493 * t517 - g(1) * (-pkin(2) * t714 + t544) - g(2) * (-pkin(2) * t716 + t542) - g(3) * t693 - t630 * qJD(1) * pkin(7)) * MDP(14) + (pkin(1) * t713 + t640) * MDP(9) - MDP(4) * t669 + t692 * MDP(5) * t597 + MDP(7) * t670 + MDP(6) * t671 + (t558 + 0.2e1 * t581 + 0.2e1 * t582 + (qJD(1) * t517 - g(3)) * t589 + (qJD(1) * t493 - t637) * t593) * MDP(13) + (-t512 * t628 + t722) * MDP(15) + (qJ(3) * t667 - t500 * t547 - t673 * t511 + (t675 + t706 + (t492 - t583) * t690) * t592 + (-t731 + t448 + (t477 - t678) * t547 + (-qJ(3) * qJD(1) * qJD(4) - t637) * t593) * t588) * MDP(20) + (qJ(3) * t444 + t697 * t547 + t673 * t512 + (-t492 * t547 - t706) * t588 + (-t547 * t678 + t448 + t608) * t592) * MDP(21); MDP(11) * t671 + (qJDD(2) + t669) * MDP(12) + (-t584 * t597 - t596) * MDP(13) + (qJD(2) * t530 + t545 + t624 - t725) * MDP(14) + (t487 + t688) * MDP(20) + (-t687 - t721) * MDP(21) + (-qJD(2) * t450 + t636) * MDP(27) + (-qJD(2) * t632 + t619) * MDP(28) + (-t398 * t513 - t450 * t698 - t743) * MDP(29) + (-t413 * qJD(2) + t742) * MDP(30) + (-MDP(21) * t547 * t592 - MDP(20) * t628) * t547; -t512 * t511 * MDP(15) + (-t511 ^ 2 + t512 ^ 2) * MDP(16) + (-t511 * t547 + t444) * MDP(17) + (t512 * t547 + t600) * MDP(18) + t510 * MDP(19) + (-g(1) * t494 - g(2) * t496 + g(3) * t709 + t428 * t547 - t492 * t512 + t604) * MDP(20) + (g(1) * t495 - g(2) * t497 + t427 * t547 - t492 * t511 + (t682 - t579) * t588 - t668) * MDP(21) + (-t651 * t535 + (-t450 * t512 + t502 * t591 - t535 * t677) * pkin(4) + t744) * MDP(27) + (t702 * t535 + (-t502 * t587 - t512 * t632 - t535 * t676) * pkin(4) + t754) * MDP(28) + (-t386 * t450 + t390 * t632 + t391 * t632 + t392 * t450 + t397 * t556 + (-t398 * t587 + (-t450 * t591 + t587 * t632) * qJD(5)) * pkin(4)) * MDP(29) + (t378 * t556 - t390 * t392 - t386 * t391 - pkin(5) * t751 - g(1) * (-t523 * t590 + t524 * t714) - g(2) * (t523 * t594 + t524 * t716) + t524 * t579 + (t379 * t587 - t413 * t512 + (-t386 * t587 + t390 * t591) * qJD(5)) * pkin(4)) * MDP(30) + t758; (t394 * t535 + t744) * MDP(27) + (t393 * t535 + t754) * MDP(28) + (pkin(5) * t397 - t450 * t705) * MDP(29) + (t705 * t390 + (t378 + t745 - t751) * pkin(5)) * MDP(30) + t758; (-t447 - t740) * MDP(29) + (-g(1) * t707 - g(2) * t711 + t386 * t632 + t390 * t450 + t385 - t731) * MDP(30);];
tau  = t1;
