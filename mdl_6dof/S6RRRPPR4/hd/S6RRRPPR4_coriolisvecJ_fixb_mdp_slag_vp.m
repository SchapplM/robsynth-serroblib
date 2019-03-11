% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:36:15
% EndTime: 2019-03-09 15:36:28
% DurationCPUTime: 8.98s
% Computational Cost: add. (5961->536), mult. (14764->705), div. (0->0), fcn. (10274->8), ass. (0->233)
t590 = cos(qJ(3));
t591 = cos(qJ(2));
t650 = qJD(1) * qJD(2);
t634 = t591 * t650;
t588 = sin(qJ(2));
t587 = sin(qJ(3));
t659 = qJD(3) * t587;
t638 = t588 * t659;
t649 = qJD(2) * qJD(3);
t500 = -qJD(1) * t638 + (t634 + t649) * t590;
t658 = qJD(3) * t590;
t637 = t588 * t658;
t660 = qJD(2) * t591;
t640 = t587 * t660;
t600 = t637 + t640;
t501 = qJD(1) * t600 + t587 * t649;
t584 = sin(pkin(10));
t585 = cos(pkin(10));
t447 = t500 * t584 + t585 * t501;
t448 = t500 * t585 - t501 * t584;
t664 = qJD(1) * t588;
t642 = t587 * t664;
t654 = t590 * qJD(2);
t536 = t642 - t654;
t662 = qJD(2) * t587;
t538 = t590 * t664 + t662;
t484 = t585 * t536 + t538 * t584;
t586 = sin(qJ(6));
t589 = cos(qJ(6));
t608 = -t536 * t584 + t585 * t538;
t655 = qJD(6) * t589;
t656 = qJD(6) * t586;
t378 = t586 * t447 + t589 * t448 + t484 * t655 - t608 * t656;
t432 = t484 * t586 + t589 * t608;
t661 = qJD(2) * t588;
t631 = MDP(28) * t661;
t663 = qJD(1) * t591;
t566 = -qJD(3) + t663;
t651 = -qJD(6) - t566;
t717 = -t589 * t484 + t586 * t608;
t702 = t717 * t651;
t727 = MDP(24) * t432 * t717 + (t432 ^ 2 - t717 ^ 2) * MDP(25) + (t378 - t702) * MDP(26) - qJD(1) * t631;
t703 = t432 * t651;
t550 = -qJD(2) * pkin(2) + pkin(7) * t664;
t601 = -pkin(3) * t536 - qJD(4) - t550;
t595 = qJ(5) * t608 + t601;
t707 = pkin(4) + pkin(5);
t393 = -t484 * t707 + t595;
t553 = t566 * qJD(5);
t577 = pkin(7) * t663;
t551 = qJD(2) * pkin(8) + t577;
t546 = -pkin(2) * t591 - pkin(8) * t588 - pkin(1);
t526 = t546 * qJD(1);
t691 = t587 * t526;
t493 = t551 * t590 + t691;
t621 = pkin(2) * t588 - pkin(8) * t591;
t540 = t621 * qJD(2);
t527 = qJD(1) * t540;
t635 = t588 * t650;
t624 = pkin(7) * t635;
t673 = -t590 * t527 - t587 * t624;
t596 = -qJD(3) * t493 - t673;
t401 = pkin(3) * t635 - qJ(4) * t500 - qJD(4) * t538 + t596;
t603 = t526 * t658 + t587 * t527 - t551 * t659;
t594 = -t590 * t624 + t603;
t407 = -qJ(4) * t501 - qJD(4) * t536 + t594;
t375 = t584 * t401 + t585 * t407;
t645 = qJ(5) * t635 + t375;
t372 = -t553 + t645;
t369 = pkin(9) * t447 + t372;
t647 = t707 * t588;
t623 = qJD(2) * t647;
t683 = -t585 * t401 + t584 * t407;
t370 = -pkin(9) * t448 - qJD(1) * t623 + t683;
t628 = t586 * t369 - t589 * t370;
t724 = t393 * t432 + t628;
t636 = t584 * t659;
t641 = t587 * t663;
t693 = t585 * t590;
t674 = -t584 * t641 - t585 * t658 + t663 * t693 + t636;
t492 = t590 * t526 - t551 * t587;
t459 = -qJ(4) * t538 + t492;
t449 = -pkin(3) * t566 + t459;
t460 = -qJ(4) * t536 + t493;
t695 = t584 * t460;
t402 = t449 * t585 - t695;
t617 = qJD(5) - t402;
t714 = pkin(9) * t608;
t385 = t566 * t707 + t617 - t714;
t646 = t589 * t369 + t586 * t370 + t385 * t655;
t722 = -t393 * t717 + t646;
t720 = pkin(9) * t484;
t719 = t484 * t566;
t530 = t584 * t590 + t585 * t587;
t520 = t530 * qJD(3);
t675 = t530 * t663 - t520;
t718 = -t577 + (-t641 + t659) * pkin(3);
t716 = t608 ^ 2;
t715 = -0.2e1 * t650;
t713 = MDP(4) * t588;
t582 = t588 ^ 2;
t712 = MDP(5) * (-t591 ^ 2 + t582);
t411 = pkin(4) * t484 - t595;
t711 = t411 * t608;
t686 = t590 * t591;
t602 = pkin(3) * t588 - qJ(4) * t686;
t539 = t621 * qJD(1);
t670 = pkin(7) * t642 + t590 * t539;
t471 = qJD(1) * t602 + t670;
t522 = t587 * t539;
t688 = t588 * t590;
t689 = t587 * t591;
t488 = t522 + (-pkin(7) * t688 - qJ(4) * t689) * qJD(1);
t420 = t584 * t471 + t585 * t488;
t414 = qJ(5) * t664 + t420;
t705 = -qJ(4) - pkin(8);
t630 = qJD(3) * t705;
t657 = qJD(4) * t590;
t515 = t587 * t630 + t657;
t516 = -qJD(4) * t587 + t590 * t630;
t467 = t585 * t515 + t584 * t516;
t682 = t414 - t467;
t676 = (-t471 + t516) * t585 + (t488 - t515) * t584;
t708 = -qJ(5) * t674 + qJD(5) * t530 - t718;
t409 = t585 * t459 - t695;
t652 = qJD(5) - t409;
t627 = -t589 * t447 + t448 * t586;
t379 = qJD(6) * t432 + t627;
t706 = pkin(7) * t587;
t694 = t585 * t460;
t408 = t459 * t584 + t694;
t704 = t408 * t608;
t701 = t500 * t587;
t700 = t536 * t566;
t699 = t538 * t566;
t698 = t550 * t587;
t697 = t550 * t590;
t696 = t566 * t590;
t403 = t584 * t449 + t694;
t397 = -t566 * qJ(5) + t403;
t388 = t397 + t720;
t692 = t586 * t388;
t690 = t587 * t588;
t592 = qJD(2) ^ 2;
t687 = t588 * t592;
t685 = t591 * t592;
t593 = qJD(1) ^ 2;
t684 = t591 * t593;
t681 = pkin(4) * t664 - t676;
t529 = t584 * t587 - t693;
t609 = t589 * t529 - t530 * t586;
t680 = qJD(6) * t609 - t586 * t675 - t589 * t674;
t481 = t529 * t586 + t530 * t589;
t679 = qJD(6) * t481 - t586 * t674 + t589 * t675;
t569 = pkin(7) * t686;
t671 = t590 * t540 + t661 * t706;
t421 = -t588 * t657 + t602 * qJD(2) + (-t569 + (qJ(4) * t588 - t546) * t587) * qJD(3) + t671;
t672 = t587 * t540 + t546 * t658;
t436 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t688 + (-qJD(4) * t588 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t591) * t587 + t672;
t387 = t584 * t421 + t585 * t436;
t678 = t675 * t707 + t708;
t677 = pkin(4) * t675 + t708;
t532 = t590 * t546;
t489 = -qJ(4) * t688 + t532 + (-pkin(3) - t706) * t591;
t668 = t587 * t546 + t569;
t494 = -qJ(4) * t690 + t668;
t438 = t584 * t489 + t585 * t494;
t548 = t705 * t587;
t549 = t705 * t590;
t497 = t584 * t548 - t585 * t549;
t667 = pkin(3) * t690 + t588 * pkin(7);
t653 = -t652 + t714;
t644 = pkin(3) * t600 + pkin(7) * t660;
t643 = -pkin(3) * t590 - pkin(2);
t572 = -pkin(3) * t585 - pkin(4);
t639 = t591 * t654;
t632 = MDP(15) * t661;
t482 = pkin(3) * t501 + pkin(7) * t634;
t629 = pkin(1) * t715;
t386 = t421 * t585 - t584 * t436;
t437 = t489 * t585 - t584 * t494;
t496 = -t585 * t548 - t549 * t584;
t626 = t536 + t654;
t625 = -t538 + t662;
t433 = -qJ(5) * t591 + t438;
t511 = -t584 * t690 + t585 * t688;
t619 = qJ(5) * t511 - t667;
t618 = -t497 * t447 + t448 * t496 - t467 * t484;
t435 = t591 * pkin(4) - t437;
t470 = pkin(9) * t529 + t497;
t616 = -pkin(9) * t674 - qJD(1) * t647 + qJD(6) * t470 + t676;
t469 = -pkin(9) * t530 + t496;
t615 = pkin(9) * t675 - qJD(6) * t469 + t682;
t614 = -pkin(3) * t538 - qJ(5) * t484;
t365 = t586 * t385 + t589 * t388;
t410 = pkin(5) * t591 - pkin(9) * t511 + t435;
t510 = t530 * t588;
t412 = pkin(9) * t510 + t433;
t613 = t410 * t589 - t412 * t586;
t612 = t410 * t586 + t412 * t589;
t610 = t589 * t510 - t511 * t586;
t463 = t510 * t586 + t511 * t589;
t564 = -pkin(5) + t572;
t570 = pkin(3) * t584 + qJ(5);
t607 = t564 * t589 - t570 * t586;
t606 = t564 * t586 + t570 * t589;
t605 = qJD(1) * t582 - t566 * t591;
t604 = qJ(5) * t530 - t643;
t382 = qJ(5) * t661 - qJD(5) * t591 + t387;
t373 = -pkin(4) * t635 + t683;
t465 = t520 * t588 + t584 * t640 - t585 * t639;
t599 = -qJ(5) * t465 + qJD(5) * t511 - t644;
t598 = qJ(5) * t448 + qJD(5) * t608 - t482;
t380 = pkin(4) * t447 - t598;
t478 = pkin(4) * t529 - t604;
t464 = -t584 * t639 - t585 * t600 + t588 * t636;
t451 = -t529 * t707 + t604;
t450 = pkin(4) * t510 - t619;
t434 = -t510 * t707 + t619;
t413 = pkin(4) * t608 - t614;
t396 = pkin(4) * t566 + t617;
t395 = -t608 * t707 + t614;
t394 = -pkin(4) * t464 - t599;
t391 = t408 + t720;
t390 = qJD(6) * t463 + t589 * t464 - t465 * t586;
t389 = qJD(6) * t610 - t464 * t586 - t465 * t589;
t384 = -pkin(4) * t661 - t386;
t381 = t464 * t707 + t599;
t377 = -pkin(9) * t464 + t382;
t376 = pkin(9) * t465 - t386 - t623;
t371 = -t447 * t707 + t598;
t364 = t385 * t589 - t692;
t1 = [(-t372 * t591 - t380 * t511 - t382 * t566 - t394 * t608 + t411 * t465 - t448 * t450 + (qJD(1) * t433 + t397) * t661) * MDP(22) + (-t372 * t510 + t373 * t511 - t382 * t484 + t384 * t608 - t396 * t465 + t397 * t464 - t433 * t447 + t435 * t448) * MDP(21) + ((-pkin(7) * t591 * t659 + t672) * t566 + t603 * t591 + (pkin(7) * t500 - t550 * t659) * t588 + ((pkin(7) * t538 + t697) * t591 + (-pkin(7) * t696 - qJD(1) * t668 - t493) * t588) * qJD(2)) * MDP(17) + (-(-t546 * t659 + t671) * t566 + (t550 * t658 + pkin(7) * t501 + (qJD(1) * t532 + t492) * qJD(2)) * t588 + ((pkin(7) * t536 + t698) * qJD(2) + (t691 + (pkin(7) * t566 + t551) * t590) * qJD(3) + t673) * t591) * MDP(16) - MDP(7) * t687 + (pkin(7) * t687 + t591 * t629) * MDP(10) + (t500 * t688 + (-t638 + t639) * t538) * MDP(11) + (-pkin(7) * t685 + t588 * t629) * MDP(9) + (-t566 - t663) * t632 + (t373 * t591 + t380 * t510 + t384 * t566 + t394 * t484 - t411 * t464 + t447 * t450 + (-qJD(1) * t435 - t396) * t661) * MDP(20) + (t566 * t638 - t500 * t591 + (t538 * t588 + t590 * t605) * qJD(2)) * MDP(13) + (t566 * t637 + t501 * t591 + (-t536 * t588 - t587 * t605) * qJD(2)) * MDP(14) + ((-t536 * t590 - t538 * t587) * t660 + (-t701 - t501 * t590 + (t536 * t587 - t538 * t590) * qJD(3)) * t588) * MDP(12) + (-t375 * t510 - t386 * t608 - t387 * t484 + t402 * t465 + t403 * t464 - t437 * t448 - t438 * t447 + t511 * t683) * MDP(18) + (t375 * t438 + t402 * t386 + t403 * t387 - t437 * t683 + t482 * t667 - t601 * t644) * MDP(19) + (t651 - t663) * t631 + (t378 * t591 - t389 * t651 + (-qJD(1) * t463 - t432) * t661) * MDP(26) + ((qJD(6) * t613 + t376 * t586 + t377 * t589) * t651 - (-t388 * t656 + t646) * t591 + t381 * t432 + t434 * t378 + t371 * t463 + t393 * t389 + (qJD(1) * t612 + t365) * t661) * MDP(30) + t712 * t715 + (t378 * t463 + t389 * t432) * MDP(24) + (t372 * t433 + t373 * t435 + t380 * t450 + t382 * t397 + t384 * t396 + t394 * t411) * MDP(23) + (t378 * t610 - t379 * t463 - t389 * t717 - t390 * t432) * MDP(25) + (-t379 * t591 + t390 * t651 + (-qJD(1) * t610 + t717) * t661) * MDP(27) + (-(t376 * t589 - t377 * t586) * t651 - t628 * t591 + t381 * t717 + t434 * t379 - t371 * t610 + t393 * t390 + (-t365 * t591 + t612 * t651) * qJD(6) + (-qJD(1) * t613 - t364) * t661) * MDP(29) + 0.2e1 * t634 * t713 + MDP(6) * t685; -t684 * t713 + t593 * t712 + (-t538 * t696 + t701) * MDP(11) + ((t500 + t700) * t590 + (-t501 + t699) * t587) * MDP(12) + (-t566 * t658 + (t566 * t686 + t588 * t625) * qJD(1)) * MDP(13) + (t566 * t659 + (-t566 * t689 + t588 * t626) * qJD(1)) * MDP(14) + (-pkin(2) * t501 + t670 * t566 + (pkin(8) * t696 + t698) * qJD(3) + ((-pkin(8) * t662 - t492) * t588 + (-pkin(7) * t626 - t698) * t591) * qJD(1)) * MDP(16) + (-pkin(2) * t500 - t522 * t566 + (-pkin(8) * t566 * t587 + t697) * qJD(3) + (-t550 * t686 + (-pkin(8) * t654 + t493) * t588 + (t566 * t688 + t591 * t625) * pkin(7)) * qJD(1)) * MDP(17) + (-t375 * t529 + t402 * t674 + t403 * t675 + t420 * t484 + t530 * t683 - t608 * t676 + t618) * MDP(18) + (t375 * t497 + t683 * t496 + t482 * t643 - t718 * t601 + (t467 - t420) * t403 + t676 * t402) * MDP(19) + (t380 * t529 - t675 * t411 + t447 * t478 - t677 * t484 + t681 * t566) * MDP(20) + (-t372 * t529 + t373 * t530 - t396 * t674 + t397 * t675 + t414 * t484 + t608 * t681 + t618) * MDP(21) + (-t380 * t530 + t674 * t411 - t448 * t478 + t682 * t566 + t608 * t677) * MDP(22) + (t372 * t497 + t373 * t496 + t380 * t478 + t396 * t681 - t397 * t682 - t411 * t677) * MDP(23) + (t378 * t481 + t432 * t680) * MDP(24) + (t378 * t609 - t379 * t481 - t432 * t679 - t680 * t717) * MDP(25) + (-t371 * t609 + t451 * t379 + t679 * t393 + t678 * t717) * MDP(29) + (t371 * t481 + t451 * t378 + t680 * t393 + t678 * t432) * MDP(30) - (t680 * MDP(26) - t679 * MDP(27) + (t586 * t615 - t589 * t616) * MDP(29) + (t586 * t616 + t589 * t615) * MDP(30)) * t651 + (t566 * MDP(15) - t651 * MDP(28) + (-qJD(2) * t496 + t396) * MDP(20) + (qJD(2) * t497 - t397) * MDP(22) + (-qJD(2) * t481 + t432) * MDP(26) + (-qJD(2) * t609 - t717) * MDP(27) + (-(t469 * t589 - t470 * t586) * qJD(2) + t364) * MDP(29) + ((t469 * t586 + t470 * t589) * qJD(2) - t365) * MDP(30)) * t664 + (MDP(9) * t588 * t593 + MDP(10) * t684) * pkin(1); t538 * t536 * MDP(11) + (-t536 ^ 2 + t538 ^ 2) * MDP(12) + (t500 - t700) * MDP(13) + (-t501 - t699) * MDP(14) + qJD(1) * t632 + (-t493 * t566 - t538 * t550 + t596) * MDP(16) + (-t492 * t566 + t536 * t550 - t594) * MDP(17) + (t403 * t608 - t704 + (-t447 * t584 - t448 * t585) * pkin(3) + (-t402 + t409) * t484) * MDP(18) + (t402 * t408 - t403 * t409 + (t375 * t584 + t538 * t601 - t585 * t683) * pkin(3)) * MDP(19) + (-t408 * t566 - t711 - t413 * t484 + (pkin(4) - t572) * t635 - t683) * MDP(20) + (t397 * t608 - t447 * t570 + t448 * t572 - t704 + (t396 - t652) * t484) * MDP(21) + (t409 * t566 - t411 * t484 + t413 * t608 + t570 * t635 - 0.2e1 * t553 + t645) * MDP(22) + (t372 * t570 + t373 * t572 - t396 * t408 + t397 * t652 - t411 * t413) * MDP(23) + (t379 + t703) * MDP(27) + (-t607 * t635 - t395 * t717 - (-t589 * t391 + t586 * t653) * t651 + (t606 * t651 + t365) * qJD(6) + t724) * MDP(29) + (t606 * t635 - t395 * t432 - (t586 * t391 + t589 * t653) * t651 + (t607 * t651 - t692) * qJD(6) + t722) * MDP(30) - t727; (t402 * t608 + t403 * t484 + t482) * MDP(19) + (-t566 * t608 + t447) * MDP(20) + (-t448 - t719) * MDP(22) + (-t396 * t608 + t397 * t484 + t380) * MDP(23) + (-t379 + t703) * MDP(29) + (-t378 - t702) * MDP(30) + (MDP(18) + MDP(21)) * (-t484 ^ 2 - t716); (t484 * t608 - t635) * MDP(20) + (t448 - t719) * MDP(21) + (-t566 ^ 2 - t716) * MDP(22) + (t397 * t566 + t373 + t711) * MDP(23) + (-t589 * t635 - t608 * t717) * MDP(29) + (-t432 * t608 + t586 * t635) * MDP(30) - (MDP(29) * t586 + MDP(30) * t589) * t651 ^ 2; (-t627 - t703) * MDP(27) + (-t365 * t651 - t724) * MDP(29) + (-t364 * t651 - t722) * MDP(30) + (-MDP(27) * t432 - MDP(29) * t365 + MDP(30) * t692) * qJD(6) + t727;];
tauc  = t1;
