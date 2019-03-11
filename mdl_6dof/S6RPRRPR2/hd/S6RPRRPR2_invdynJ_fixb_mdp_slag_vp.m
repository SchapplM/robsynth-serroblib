% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:02:46
% EndTime: 2019-03-09 05:03:00
% DurationCPUTime: 8.91s
% Computational Cost: add. (6031->530), mult. (12979->705), div. (0->0), fcn. (9163->16), ass. (0->239)
t603 = cos(qJ(6));
t600 = sin(qJ(4));
t604 = cos(qJ(4));
t669 = t604 * qJD(3);
t601 = sin(qJ(3));
t681 = qJD(1) * t601;
t548 = t600 * t681 - t669;
t677 = qJD(3) * t600;
t550 = t604 * t681 + t677;
t594 = sin(pkin(11));
t596 = cos(pkin(11));
t635 = -t548 * t596 - t550 * t594;
t697 = t603 * t635;
t487 = t548 * t594 - t550 * t596;
t599 = sin(qJ(6));
t709 = t487 * t599;
t433 = t697 + t709;
t605 = cos(qJ(3));
t680 = qJD(1) * t605;
t573 = -qJD(4) + t680;
t568 = -qJD(6) + t573;
t711 = t433 * t568;
t595 = sin(pkin(10));
t577 = pkin(1) * t595 + pkin(7);
t562 = t577 * qJD(1);
t521 = t601 * qJD(2) + t605 * t562;
t504 = qJD(3) * pkin(8) + t521;
t597 = cos(pkin(10));
t579 = -pkin(1) * t597 - pkin(2);
t536 = -pkin(3) * t605 - pkin(8) * t601 + t579;
t507 = t536 * qJD(1);
t456 = t504 * t604 + t507 * t600;
t560 = t577 * qJDD(1);
t520 = qJD(2) * t605 - t601 * t562;
t731 = qJD(3) * t520;
t464 = qJDD(3) * pkin(8) + qJDD(2) * t601 + t560 * t605 + t731;
t644 = pkin(3) * t601 - pkin(8) * t605;
t553 = t644 * qJD(3);
t484 = qJD(1) * t553 + qJDD(1) * t536;
t475 = t604 * t484;
t668 = qJD(1) * qJD(3);
t652 = t605 * t668;
t666 = qJDD(1) * t601;
t673 = qJD(4) * t601;
t730 = -qJD(1) * t673 + qJDD(3);
t482 = qJD(4) * t669 + (t652 + t666) * t604 + t730 * t600;
t587 = t605 * qJDD(1);
t539 = t601 * t668 + qJDD(4) - t587;
t388 = pkin(4) * t539 - qJ(5) * t482 - qJD(4) * t456 - qJD(5) * t550 - t464 * t600 + t475;
t483 = t600 * (qJD(3) * (qJD(4) + t680) + t666) - t730 * t604;
t672 = qJD(4) * t604;
t661 = t604 * t464 + t600 * t484 + t507 * t672;
t674 = qJD(4) * t600;
t620 = -t504 * t674 + t661;
t390 = -qJ(5) * t483 - qJD(5) * t548 + t620;
t378 = t596 * t388 - t390 * t594;
t424 = t482 * t596 - t483 * t594;
t376 = pkin(5) * t539 - pkin(9) * t424 + t378;
t379 = t594 * t388 + t596 * t390;
t423 = -t482 * t594 - t483 * t596;
t377 = pkin(9) * t423 + t379;
t455 = -t504 * t600 + t604 * t507;
t442 = -qJ(5) * t550 + t455;
t435 = -pkin(4) * t573 + t442;
t443 = -qJ(5) * t548 + t456;
t702 = t596 * t443;
t401 = t594 * t435 + t702;
t728 = pkin(9) * t635;
t395 = t401 + t728;
t670 = qJD(6) * t599;
t394 = t395 * t670;
t503 = -qJD(3) * pkin(3) - t520;
t478 = pkin(4) * t548 + qJD(5) + t503;
t436 = -pkin(5) * t635 + t478;
t588 = qJ(4) + pkin(11) + qJ(6);
t575 = sin(t588);
t576 = cos(t588);
t591 = qJ(1) + pkin(10);
t584 = cos(t591);
t583 = sin(t591);
t704 = t583 * t605;
t494 = t575 * t584 - t576 * t704;
t703 = t584 * t605;
t496 = t575 * t583 + t576 * t703;
t715 = g(3) * t601;
t742 = g(1) * t496 - g(2) * t494 - t599 * t376 - t603 * t377 - t436 * t433 + t576 * t715 + t394;
t535 = qJDD(6) + t539;
t729 = -t603 * t487 + t599 * t635;
t741 = t535 * MDP(25) + (-t433 ^ 2 + t729 ^ 2) * MDP(22) - t433 * MDP(21) * t729;
t675 = qJD(3) * t605;
t732 = -qJD(2) * qJD(3) - t560;
t660 = -t562 * t675 + t601 * t732;
t631 = -qJDD(3) * pkin(3) - t660;
t665 = qJDD(2) * t605;
t465 = t631 - t665;
t643 = g(1) * t584 + g(2) * t583;
t624 = t643 * t601;
t615 = -g(3) * t605 + t624;
t740 = qJD(4) * pkin(8) * t573 - t465 + t615;
t712 = t729 * t568;
t552 = t644 * qJD(1);
t533 = t604 * t552;
t696 = t604 * t605;
t720 = pkin(4) * t601;
t629 = -qJ(5) * t696 + t720;
t713 = qJ(5) + pkin(8);
t650 = qJD(4) * t713;
t738 = -qJD(1) * t629 - t604 * t650 - t533 + (-qJD(5) + t520) * t600;
t656 = t600 * t680;
t671 = qJD(5) * t604;
t688 = t604 * t520 + t600 * t552;
t737 = -qJ(5) * t656 + t600 * t650 - t671 + t688;
t493 = t575 * t704 + t576 * t584;
t495 = -t575 * t703 + t576 * t583;
t649 = t603 * t376 - t599 * t377;
t736 = -g(1) * t495 + g(2) * t493 - t436 * t729 + t575 * t715 + t649;
t735 = pkin(9) * t487;
t541 = t594 * t604 + t596 * t600;
t621 = t541 * t605;
t734 = qJD(1) * t621 - t541 * qJD(4);
t632 = t594 * t600 - t596 * t604;
t733 = t573 * t632;
t695 = qJDD(2) - g(3);
t726 = t605 * t695;
t691 = t594 * t737 + t596 * t738;
t690 = t594 * t738 - t596 * t737;
t724 = -t521 + (-t656 + t674) * pkin(4);
t700 = t600 * t605;
t512 = t583 * t700 + t584 * t604;
t514 = t583 * t604 - t584 * t700;
t723 = -g(1) * t514 + g(2) * t512;
t648 = -t603 * t423 + t424 * t599;
t383 = qJD(6) * t729 + t648;
t721 = pkin(4) * t594;
t710 = t482 * t600;
t708 = t548 * t573;
t707 = t550 * t573;
t706 = t573 * t604;
t705 = t577 * t600;
t438 = t594 * t443;
t701 = t600 * t601;
t699 = t601 * t604;
t400 = t596 * t435 - t438;
t391 = -pkin(5) * t573 + t400 + t735;
t698 = t603 * t391;
t518 = t541 * t601;
t519 = t632 * t601;
t463 = -t518 * t599 - t519 * t603;
t466 = -qJD(3) * t621 + t632 * t673;
t653 = t605 * t669;
t654 = t600 * t675;
t467 = t541 * t673 + t594 * t654 - t596 * t653;
t403 = qJD(6) * t463 - t603 * t466 - t467 * t599;
t637 = -t603 * t518 + t519 * t599;
t694 = t403 * t568 + t535 * t637;
t551 = t577 * t696;
t676 = qJD(3) * t601;
t686 = t604 * t553 + t676 * t705;
t414 = -t601 * t671 + t629 * qJD(3) + (-t551 + (qJ(5) * t601 - t536) * t600) * qJD(4) + t686;
t678 = qJD(3) * t577;
t687 = t536 * t672 + t600 * t553;
t418 = (-qJ(5) * qJD(4) - t678) * t699 + (-qJD(5) * t601 + (-qJ(5) * qJD(3) - qJD(4) * t577) * t605) * t600 + t687;
t393 = t594 * t414 + t596 * t418;
t636 = -t541 * t599 - t603 * t632;
t693 = qJD(6) * t636 + t599 * t734 + t603 * t733;
t486 = t541 * t603 - t599 * t632;
t692 = qJD(6) * t486 + t599 * t733 - t603 * t734;
t405 = t596 * t442 - t438;
t523 = t604 * t536;
t468 = -qJ(5) * t699 + t523 + (-pkin(4) - t705) * t605;
t685 = t600 * t536 + t551;
t477 = -qJ(5) * t701 + t685;
t417 = t594 * t468 + t596 * t477;
t689 = -pkin(5) * t734 + t724;
t564 = t713 * t600;
t565 = t713 * t604;
t492 = -t594 * t564 + t596 * t565;
t684 = pkin(4) * t701 + t601 * t577;
t592 = t601 ^ 2;
t683 = -t605 ^ 2 + t592;
t563 = qJD(1) * t579;
t662 = qJD(6) * t697 + t599 * t423 + t603 * t424;
t658 = pkin(4) * t654 + t577 * t675 + t672 * t720;
t582 = pkin(4) * t604 + pkin(3);
t657 = pkin(4) * t600 + pkin(7);
t655 = t573 * t677;
t392 = t596 * t414 - t418 * t594;
t404 = -t442 * t594 - t702;
t416 = t596 * t468 - t477 * t594;
t647 = t573 * t577 + t504;
t491 = -t596 * t564 - t565 * t594;
t646 = -qJD(4) * t507 - t464;
t642 = g(1) * t583 - g(2) * t584;
t602 = sin(qJ(1));
t606 = cos(qJ(1));
t641 = g(1) * t602 - g(2) * t606;
t472 = -pkin(9) * t632 + t492;
t640 = pkin(5) * t681 + pkin(9) * t733 + qJD(6) * t472 - t691;
t471 = -pkin(9) * t541 + t491;
t639 = pkin(9) * t734 + qJD(6) * t471 + t690;
t381 = t599 * t391 + t603 * t395;
t402 = qJD(6) * t637 + t466 * t599 - t467 * t603;
t638 = t402 * t568 - t463 * t535;
t634 = t582 * t605 + t601 * t713;
t630 = t641 * pkin(1);
t578 = pkin(4) * t596 + pkin(5);
t628 = t578 * t599 + t603 * t721;
t627 = t578 * t603 - t599 * t721;
t625 = pkin(2) + t634;
t623 = -t539 * t600 + t573 * t672;
t622 = -t539 * t604 - t573 * t674;
t382 = t487 * t670 + t662;
t619 = -qJD(1) * t563 + t643;
t618 = -t600 * t673 + t653;
t617 = -pkin(8) * t539 - t503 * t573;
t616 = 0.2e1 * qJD(3) * t563 - qJDD(3) * t577;
t612 = pkin(4) * t483 + qJDD(5) + t631;
t607 = qJD(3) ^ 2;
t611 = -0.2e1 * qJDD(1) * t579 - t577 * t607 + t642;
t422 = t612 - t665;
t556 = qJDD(3) * t605 - t601 * t607;
t555 = qJDD(3) * t601 + t605 * t607;
t524 = t550 * t676;
t515 = t583 * t600 + t584 * t696;
t513 = -t583 * t696 + t584 * t600;
t508 = pkin(5) * t632 - t582;
t479 = pkin(5) * t518 + t684;
t450 = pkin(4) * t550 - pkin(5) * t487;
t437 = -pkin(5) * t466 + t658;
t425 = t729 * t676;
t409 = -pkin(9) * t518 + t417;
t408 = -pkin(5) * t605 + pkin(9) * t519 + t416;
t398 = -pkin(5) * t423 + t422;
t397 = t405 + t735;
t396 = t404 - t728;
t385 = pkin(9) * t466 + t393;
t384 = pkin(5) * t676 + pkin(9) * t467 + t392;
t380 = -t395 * t599 + t698;
t1 = [(t378 * t519 - t379 * t518 + t392 * t487 + t393 * t635 + t400 * t467 + t401 * t466 - t416 * t424 + t417 * t423 + t601 * t642) * MDP(19) + (-(t384 * t603 - t385 * t599) * t568 + (t408 * t603 - t409 * t599) * t535 - t649 * t605 + t380 * t676 - t437 * t433 + t479 * t383 - t398 * t637 + t436 * t403 - g(1) * t494 - g(2) * t496 + (-(-t408 * t599 - t409 * t603) * t568 + t381 * t605) * qJD(6)) * MDP(26) + (t382 * t637 - t383 * t463 + t402 * t433 - t403 * t729) * MDP(22) + (t383 * t605 + t433 * t676 + t694) * MDP(24) + (-t381 * t676 - g(1) * t493 - g(2) * t495 + t479 * t382 - t394 * t605 + t398 * t463 + t436 * t402 + t437 * t729 + ((-qJD(6) * t409 + t384) * t568 - t408 * t535 + t376 * t605) * t599 + ((qJD(6) * t408 + t385) * t568 - t409 * t535 + (qJD(6) * t391 + t377) * t605) * t603) * MDP(27) + (t382 * t463 + t402 * t729) * MDP(21) + ((-t548 * t604 - t550 * t600) * t675 + (-t710 - t483 * t604 + (t548 * t600 - t550 * t604) * qJD(4)) * t601) * MDP(13) + (t687 * t573 - t685 * t539 - g(1) * t512 - g(2) * t514 + (-t647 * t674 + (t503 * t604 + t550 * t577) * qJD(3) + t661) * t605 + (-t503 * t674 + t465 * t604 + t577 * t482 + (-t577 * t706 - t456) * qJD(3)) * t601) * MDP(18) + (-t482 * t605 + t539 * t699 - t573 * t618 + t524) * MDP(14) + (t482 * t699 + t550 * t618) * MDP(12) + (t379 * t417 + t401 * t393 + t378 * t416 + t400 * t392 + t422 * t684 + t478 * t658 + t630 + (-g(1) * t657 - g(2) * t625) * t584 + (g(1) * t625 - g(2) * t657) * t583) * MDP(20) + (-(-t536 * t674 + t686) * t573 + t523 * t539 - g(1) * t513 - g(2) * t515 + (t548 * t678 - t475 + t647 * t672 + (qJD(3) * t503 - t539 * t577 - t646) * t600) * t605 + (qJD(3) * t455 + t465 * t600 + t483 * t577 + t503 * t672) * t601) * MDP(17) + 0.2e1 * (t587 * t601 - t668 * t683) * MDP(6) + (-t539 * t605 - t573 * t676) * MDP(16) + (-t535 * t605 - t568 * t676) * MDP(25) + ((t483 + t655) * t605 + (-qJD(3) * t548 + t623) * t601) * MDP(15) + (qJDD(1) * t592 + 0.2e1 * t601 * t652) * MDP(5) + (g(1) * t606 + g(2) * t602) * MDP(3) + qJDD(1) * MDP(1) + t641 * MDP(2) + (-t382 * t605 + t425 - t638) * MDP(23) + (t601 * t616 + t605 * t611) * MDP(10) + (-t601 * t611 + t605 * t616) * MDP(11) + t555 * MDP(7) + t556 * MDP(8) + ((t595 ^ 2 + t597 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t630) * MDP(4); t695 * MDP(4) + t556 * MDP(10) - t555 * MDP(11) + t524 * MDP(18) + (-t423 * t519 + t424 * t518 + t466 * t487 - t467 * t635) * MDP(19) + (-t378 * t518 - t379 * t519 + t400 * t466 - t401 * t467 - g(3)) * MDP(20) + t694 * MDP(26) + (t425 + t638) * MDP(27) + ((-t483 + t655) * MDP(17) + (t573 * t669 - t482) * MDP(18) - t422 * MDP(20) - t383 * MDP(26) - t382 * MDP(27)) * t605 + (t623 * MDP(17) + t622 * MDP(18) + (MDP(17) * t548 + MDP(20) * t478 - MDP(26) * t433) * qJD(3)) * t601; MDP(7) * t666 + MDP(8) * t587 + qJDD(3) * MDP(9) + (qJD(3) * t521 + t601 * t619 + t660 + t726) * MDP(10) + (t731 + (qJD(3) * t562 - t695) * t601 + (t619 + t732) * t605) * MDP(11) + (-t550 * t706 + t710) * MDP(12) + ((t482 + t708) * t604 + (-t483 + t707) * t600) * MDP(13) + ((-t550 * t601 + t573 * t696) * qJD(1) - t623) * MDP(14) + ((t548 * t601 - t573 * t700) * qJD(1) - t622) * MDP(15) + (-pkin(3) * t483 - t521 * t548 + t533 * t573 + (-t520 * t573 + t617) * t600 + t740 * t604) * MDP(17) + (-pkin(3) * t482 - t521 * t550 - t688 * t573 - t600 * t740 + t617 * t604) * MDP(18) + (-t378 * t541 - t379 * t632 - t400 * t733 + t401 * t734 + t423 * t492 - t424 * t491 + t691 * t487 - t643 * t605 + t690 * t635 - t715) * MDP(19) + (t379 * t492 + t378 * t491 - t422 * t582 - g(3) * t634 + t724 * t478 + t690 * t401 + t691 * t400 + t643 * (t582 * t601 - t605 * t713)) * MDP(20) + (t382 * t486 + t693 * t729) * MDP(21) + (t382 * t636 - t383 * t486 + t433 * t693 - t692 * t729) * MDP(22) + (t486 * t535 - t568 * t693) * MDP(23) + (t535 * t636 + t568 * t692) * MDP(24) + ((t471 * t603 - t472 * t599) * t535 + t508 * t383 - t398 * t636 + (t599 * t639 + t603 * t640) * t568 + t692 * t436 - t689 * t433 + t615 * t576) * MDP(26) + (-(t471 * t599 + t472 * t603) * t535 + t508 * t382 + t398 * t486 + (-t599 * t640 + t603 * t639) * t568 + t693 * t436 + t689 * t729 - t615 * t575) * MDP(27) + (t573 * MDP(16) - t455 * MDP(17) + MDP(18) * t456 - MDP(23) * t729 - MDP(24) * t433 + t568 * MDP(25) - t380 * MDP(26) + t381 * MDP(27)) * t681 + (-MDP(5) * t601 * t605 + MDP(6) * t683) * qJD(1) ^ 2; t550 * t548 * MDP(12) + (-t548 ^ 2 + t550 ^ 2) * MDP(13) + (t482 - t708) * MDP(14) + (-t483 - t707) * MDP(15) + t539 * MDP(16) + (-t504 * t672 - t456 * t573 - t503 * t550 + t475 + (t646 + t715) * t600 + t723) * MDP(17) + (g(1) * t515 - g(2) * t513 + g(3) * t699 - t455 * t573 + t503 * t548 - t620) * MDP(18) + ((t423 * t594 - t424 * t596) * pkin(4) + (t400 - t405) * t635 + (-t401 - t404) * t487) * MDP(19) + (-t400 * t404 - t401 * t405 + (g(3) * t701 + t378 * t596 + t379 * t594 - t478 * t550 + t723) * pkin(4)) * MDP(20) + (t382 + t711) * MDP(23) + (-t383 - t712) * MDP(24) + (t627 * t535 + (t396 * t603 - t397 * t599) * t568 + t450 * t433 + (t568 * t628 - t381) * qJD(6) + t736) * MDP(26) + (-t628 * t535 - (t396 * t599 + t397 * t603) * t568 - t450 * t729 + (t568 * t627 - t698) * qJD(6) + t742) * MDP(27) + t741; (-t487 ^ 2 - t635 ^ 2) * MDP(19) + (-t400 * t487 - t401 * t635 + t612 - t624 - t726) * MDP(20) + (t383 - t712) * MDP(26) + (t382 - t711) * MDP(27); (t662 + t711) * MDP(23) + (-t648 - t712) * MDP(24) + (-t381 * t568 + t736) * MDP(26) + (-t380 * t568 + t742) * MDP(27) + (MDP(23) * t709 - MDP(24) * t729 - MDP(26) * t381 - MDP(27) * t698) * qJD(6) + t741;];
tau  = t1;
