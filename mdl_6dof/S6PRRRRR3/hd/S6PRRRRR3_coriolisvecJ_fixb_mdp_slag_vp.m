% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRRR3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:51:49
% EndTime: 2019-03-09 00:52:07
% DurationCPUTime: 11.56s
% Computational Cost: add. (6024->546), mult. (14957->756), div. (0->0), fcn. (11509->12), ass. (0->238)
t581 = sin(qJ(4));
t586 = cos(qJ(4));
t654 = t586 * qJD(3);
t582 = sin(qJ(3));
t665 = qJD(2) * t582;
t533 = t581 * t665 - t654;
t663 = qJD(3) * t581;
t535 = t586 * t665 + t663;
t580 = sin(qJ(5));
t585 = cos(qJ(5));
t470 = t585 * t533 + t535 * t580;
t584 = cos(qJ(6));
t579 = sin(qJ(6));
t607 = t533 * t580 - t585 * t535;
t708 = t607 * t579;
t424 = -t584 * t470 + t708;
t609 = t470 * t579 + t584 * t607;
t653 = qJD(2) * qJD(3);
t636 = t582 * t653;
t743 = MDP(30) * t636 + (-t424 ^ 2 + t609 ^ 2) * MDP(27) + t424 * MDP(26) * t609;
t587 = cos(qJ(3));
t635 = t587 * t653;
t660 = qJD(4) * t581;
t638 = t582 * t660;
t652 = qJD(3) * qJD(4);
t497 = -qJD(2) * t638 + (t635 + t652) * t586;
t659 = qJD(4) * t586;
t637 = t582 * t659;
t661 = qJD(3) * t587;
t641 = t581 * t661;
t596 = t637 + t641;
t498 = qJD(2) * t596 + t581 * t652;
t591 = qJD(5) * t607 - t497 * t580 - t585 * t498;
t583 = sin(qJ(2));
t577 = sin(pkin(6));
t668 = qJD(1) * t577;
t647 = t583 * t668;
t547 = qJD(2) * pkin(8) + t647;
t578 = cos(pkin(6));
t699 = t578 * t582;
t562 = qJD(1) * t699;
t504 = t587 * t547 + t562;
t495 = qJD(3) * pkin(9) + t504;
t549 = -pkin(3) * t587 - pkin(9) * t582 - pkin(2);
t588 = cos(qJ(2));
t646 = t588 * t668;
t506 = qJD(2) * t549 - t646;
t705 = t506 * t581;
t443 = t495 * t586 + t705;
t666 = qJD(2) * t577;
t643 = t588 * t666;
t662 = qJD(3) * t582;
t667 = qJD(1) * t587;
t459 = -t547 * t662 + (qJD(3) * t578 + t643) * t667;
t616 = pkin(3) * t582 - pkin(9) * t587;
t544 = t616 * qJD(3);
t502 = (t544 + t647) * qJD(2);
t625 = t459 * t581 - t586 * t502;
t594 = -qJD(4) * t443 - t625;
t389 = pkin(4) * t636 - pkin(10) * t497 + t594;
t599 = t586 * t459 - t495 * t660 + t581 * t502 + t506 * t659;
t394 = -pkin(10) * t498 + t599;
t442 = -t495 * t581 + t586 * t506;
t428 = -pkin(10) * t535 + t442;
t664 = qJD(2) * t587;
t565 = -qJD(4) + t664;
t411 = -pkin(4) * t565 + t428;
t429 = -pkin(10) * t533 + t443;
t657 = qJD(5) * t585;
t658 = qJD(5) * t580;
t619 = -t580 * t389 - t585 * t394 - t411 * t657 + t429 * t658;
t370 = pkin(11) * t591 - t619;
t402 = t585 * t497 - t580 * t498 - t533 * t657 - t535 * t658;
t418 = t585 * t429;
t391 = t411 * t580 + t418;
t629 = t585 * t389 - t580 * t394;
t593 = -t391 * qJD(5) + t629;
t369 = pkin(5) * t636 - pkin(11) * t402 + t593;
t733 = pkin(11) * t470;
t382 = t391 - t733;
t656 = qJD(6) * t579;
t380 = t382 * t656;
t630 = t579 * t369 - t380;
t503 = -t582 * t547 + t578 * t667;
t494 = -qJD(3) * pkin(3) - t503;
t455 = pkin(4) * t533 + t494;
t413 = pkin(5) * t470 + t455;
t731 = t413 * t424;
t742 = -t584 * t370 - t630 - t731;
t655 = qJD(6) * t584;
t649 = t584 * t402 - t470 * t655 + t579 * t591;
t375 = t607 * t656 + t649;
t557 = -qJD(5) + t565;
t628 = t402 * t579 - t584 * t591;
t592 = qJD(6) * t609 - t628;
t632 = MDP(23) * t662;
t552 = -qJD(6) + t557;
t729 = t552 * t609;
t730 = t424 * t552;
t741 = qJD(2) * t632 + (-t470 ^ 2 + t607 ^ 2) * MDP(20) + (-t470 * t557 + t402) * MDP(21) + (t557 * t607 + t591) * MDP(22) - t470 * MDP(19) * t607 + (t592 + t729) * MDP(29) + (t375 + t730) * MDP(28) + t743;
t689 = t587 * t588;
t738 = -(t581 * t583 + t586 * t689) * t668 + t581 * t544 + t549 * t659;
t712 = pkin(8) * t581;
t737 = t586 * t544 + t662 * t712 - (-t581 * t689 + t583 * t586) * t668;
t631 = t584 * t369 - t579 * t370;
t726 = t413 * t609 + t631;
t690 = t586 * t587;
t567 = pkin(8) * t690;
t605 = pkin(4) * t582 - pkin(10) * t690;
t736 = t605 * qJD(3) + (-t567 + (pkin(10) * t582 - t549) * t581) * qJD(4) + t737;
t735 = -t596 * pkin(10) + (-t654 * t582 - t587 * t660) * pkin(8) + t738;
t541 = t616 * qJD(2);
t622 = -t503 * t581 + t586 * t541;
t713 = pkin(9) + pkin(10);
t648 = qJD(4) * t713;
t734 = qJD(2) * t605 + t586 * t648 + t622;
t642 = t581 * t664;
t677 = t586 * t503 + t581 * t541;
t727 = pkin(10) * t642 - t581 * t648 - t677;
t732 = pkin(11) * t607;
t536 = t580 * t581 - t585 * t586;
t600 = t536 * t587;
t714 = qJD(4) + qJD(5);
t679 = qJD(2) * t600 - t714 * t536;
t537 = t580 * t586 + t581 * t585;
t678 = (-t664 + t714) * t537;
t725 = t455 * t470 + t619;
t724 = t455 * t607 + t593;
t721 = MDP(5) * t582;
t575 = t582 ^ 2;
t720 = MDP(6) * (-t587 ^ 2 + t575);
t511 = t537 * t582;
t719 = t734 * t585;
t718 = t735 * t580 - t585 * t736;
t532 = t586 * t549;
t694 = t582 * t586;
t476 = -pkin(10) * t694 + t532 + (-pkin(4) - t712) * t587;
t672 = t581 * t549 + t567;
t696 = t581 * t582;
t486 = -pkin(10) * t696 + t672;
t717 = t476 * t657 - t486 * t658 + t580 * t736 + t735 * t585;
t617 = -t504 + (-t642 + t660) * pkin(4);
t680 = t580 * t476 + t585 * t486;
t553 = t713 * t581;
t554 = t713 * t586;
t674 = -t580 * t553 + t585 * t554;
t716 = -t553 * t657 - t554 * t658 - t580 * t734 + t585 * t727;
t715 = t587 * t654 - t638;
t711 = qJD(2) * pkin(2);
t618 = t582 * t643;
t460 = qJD(1) * t618 + qJD(3) * t562 + t547 * t661;
t710 = t460 * t581;
t709 = t460 * t586;
t707 = t494 * t581;
t706 = t497 * t581;
t704 = t533 * t565;
t703 = t535 * t565;
t702 = t565 * t586;
t701 = t577 * t583;
t700 = t577 * t588;
t698 = t579 * t580;
t416 = t580 * t429;
t697 = t580 * t584;
t695 = t581 * t587;
t589 = qJD(3) ^ 2;
t693 = t582 * t589;
t390 = t585 * t411 - t416;
t381 = t390 + t732;
t379 = -pkin(5) * t557 + t381;
t692 = t584 * t379;
t691 = t584 * t382;
t688 = t587 * t589;
t437 = -qJD(3) * t600 - t511 * t714;
t687 = -pkin(5) * t662 + pkin(11) * t437 + qJD(5) * t680 + t718;
t438 = -t658 * t696 + (t694 * t714 + t641) * t585 + t715 * t580;
t686 = -pkin(11) * t438 + t717;
t474 = t584 * t536 + t537 * t579;
t685 = -qJD(6) * t474 - t579 * t678 + t584 * t679;
t475 = -t536 * t579 + t537 * t584;
t684 = qJD(6) * t475 + t579 * t679 + t584 * t678;
t683 = t585 * t428 - t416;
t681 = pkin(5) * t678 + t617;
t545 = pkin(4) * t696 + t582 * pkin(8);
t650 = pkin(4) * qJD(5) * t552;
t505 = pkin(4) * t596 + pkin(8) * t661;
t571 = -pkin(4) * t586 - pkin(3);
t644 = t583 * t666;
t639 = t565 * t660;
t633 = MDP(16) * t662;
t627 = -t428 * t580 - t418;
t624 = t585 * t476 - t486 * t580;
t621 = -t585 * t553 - t554 * t580;
t620 = qJD(6) * t379 + t370;
t548 = -t646 - t711;
t615 = -t548 - t646;
t457 = -pkin(11) * t536 + t674;
t614 = pkin(5) * t665 + pkin(11) * t679 + qJD(5) * t674 + qJD(6) * t457 + t580 * t727 + t719;
t456 = -pkin(11) * t537 + t621;
t613 = -pkin(11) * t678 + qJD(6) * t456 + t716;
t372 = t579 * t379 + t691;
t512 = t536 * t582;
t406 = -pkin(5) * t587 + pkin(11) * t512 + t624;
t407 = -pkin(11) * t511 + t680;
t612 = t406 * t579 + t407 * t584;
t518 = t587 * t701 + t699;
t484 = -t518 * t581 - t586 * t700;
t602 = -t518 * t586 + t581 * t700;
t432 = t484 * t585 + t580 * t602;
t433 = t484 * t580 - t585 * t602;
t611 = t432 * t584 - t433 * t579;
t610 = t432 * t579 + t433 * t584;
t452 = t584 * t511 - t512 * t579;
t453 = -t511 * t579 - t512 * t584;
t606 = qJD(2) * t575 - t565 * t587;
t436 = pkin(4) * t498 + t460;
t570 = pkin(4) * t585 + pkin(5);
t604 = pkin(4) * t697 + t570 * t579;
t603 = -pkin(4) * t698 + t570 * t584;
t517 = -t578 * t587 + t582 * t701;
t595 = qJD(3) * (-t615 - t711);
t590 = qJD(2) ^ 2;
t510 = pkin(5) * t536 + t571;
t483 = qJD(3) * t518 + t618;
t482 = -qJD(3) * t517 + t587 * t643;
t478 = pkin(5) * t511 + t545;
t449 = pkin(4) * t535 - pkin(5) * t607;
t415 = qJD(4) * t484 + t482 * t586 + t581 * t644;
t414 = qJD(4) * t602 - t482 * t581 + t586 * t644;
t408 = pkin(5) * t438 + t505;
t392 = -pkin(5) * t591 + t436;
t386 = qJD(6) * t453 + t437 * t579 + t584 * t438;
t385 = -qJD(6) * t452 + t437 * t584 - t438 * t579;
t384 = t683 + t732;
t383 = t627 + t733;
t378 = -qJD(5) * t433 + t414 * t585 - t415 * t580;
t377 = qJD(5) * t432 + t414 * t580 + t415 * t585;
t371 = -t382 * t579 + t692;
t1 = [(-t414 * t565 + t483 * t533 + t498 * t517) * MDP(17) + (t415 * t565 + t483 * t535 + t497 * t517) * MDP(18) + (-t378 * t557 + t470 * t483 - t517 * t591) * MDP(24) + (t377 * t557 + t402 * t517 - t483 * t607) * MDP(25) + (-(-qJD(6) * t610 - t377 * t579 + t378 * t584) * t552 - t483 * t424 - t517 * t592) * MDP(31) + ((qJD(6) * t611 + t377 * t584 + t378 * t579) * t552 - t483 * t609 + t517 * t375) * MDP(32) + (-t483 * MDP(10) - t482 * MDP(11) + (t484 * MDP(17) + MDP(18) * t602 + t432 * MDP(24) - t433 * MDP(25) + MDP(31) * t611 - MDP(32) * t610) * t665) * qJD(3) + ((-MDP(10) * t582 - MDP(11) * t587) * t588 * t653 + (-t588 * MDP(4) + (-MDP(10) * t587 + MDP(11) * t582 - MDP(3)) * t583) * t590) * t577; ((t549 * t660 - t737) * t565 + ((pkin(8) * t533 + t707) * qJD(3) + (t705 + (pkin(8) * t565 + t495) * t586) * qJD(4) + t625) * t587 + (-t533 * t646 + t494 * t659 + pkin(8) * t498 + t710 + ((-pkin(8) * t695 + t532) * qJD(2) + t442) * qJD(3)) * t582) * MDP(17) + (t738 * t565 + (t494 * t654 + (qJD(3) * t535 - t639) * pkin(8) + t599) * t587 + (-t535 * t646 - t494 * t660 + pkin(8) * t497 + t709 + (-pkin(8) * t702 - qJD(2) * t672 - t443) * qJD(3)) * t582) * MDP(18) + (-t591 * t587 + t438 * t557 + (-qJD(2) * t511 - t470) * t662) * MDP(22) + (-t402 * t511 - t437 * t470 + t438 * t607 - t512 * t591) * MDP(20) + (-t629 * t587 + t505 * t470 - t545 * t591 + t436 * t511 + t455 * t438 + t718 * t557 + (t391 * t587 + t557 * t680) * qJD(5) + (-t470 * t646 + (qJD(2) * t624 + t390) * qJD(3)) * t582) * MDP(24) + ((-t533 * t586 - t535 * t581) * t661 + (-t706 - t498 * t586 + (t533 * t581 - t535 * t586) * qJD(4)) * t582) * MDP(13) - MDP(8) * t693 + (pkin(8) * t693 + t587 * t595) * MDP(11) + (-pkin(8) * t688 + t582 * t595) * MDP(10) + (-t557 - t664) * t632 + (-t565 - t664) * t633 + (-t402 * t587 - t437 * t557 + (-qJD(2) * t512 - t607) * t662) * MDP(21) + (-t402 * t512 - t437 * t607) * MDP(19) + (-t619 * t587 - t505 * t607 + t545 * t402 - t436 * t512 + t455 * t437 + t717 * t557 + (t607 * t646 + (-qJD(2) * t680 - t391) * qJD(3)) * t582) * MDP(25) + ((t620 * t584 + t630) * t587 - t408 * t609 + t478 * t375 + t392 * t453 + t413 * t385 + ((qJD(6) * t406 + t686) * t584 + (-qJD(6) * t407 - t687) * t579) * t552 + (t609 * t646 + (-qJD(2) * t612 - t372) * qJD(3)) * t582) * MDP(32) + (-t375 * t587 - t385 * t552 + (qJD(2) * t453 - t609) * t662) * MDP(28) + (t375 * t453 - t385 * t609) * MDP(26) + (-t552 - t664) * MDP(30) * t662 + (-t631 * t587 - t408 * t424 - t478 * t592 + t392 * t452 + t413 * t386 + (t579 * t686 + t584 * t687) * t552 + (t372 * t587 + t552 * t612) * qJD(6) + (t424 * t646 + ((t406 * t584 - t407 * t579) * qJD(2) + t371) * qJD(3)) * t582) * MDP(31) + (t565 * t637 + t498 * t587 + (-t533 * t582 - t581 * t606) * qJD(3)) * MDP(15) + (t565 * t638 - t497 * t587 + (t535 * t582 + t586 * t606) * qJD(3)) * MDP(14) - 0.2e1 * t653 * t720 + 0.2e1 * t635 * t721 + (-t592 * t587 + t386 * t552 + (-qJD(2) * t452 + t424) * t662) * MDP(29) + (-t375 * t452 + t385 * t424 + t386 * t609 + t453 * t592) * MDP(27) + (t497 * t694 + t535 * t715) * MDP(12) + MDP(7) * t688; (-pkin(3) * t498 - t709 + t622 * t565 - t504 * t533 + (pkin(9) * t702 + t707) * qJD(4) + (-t442 * t582 + (-pkin(9) * t662 - t494 * t587) * t581) * qJD(2)) * MDP(17) + (-pkin(3) * t497 + t710 - t677 * t565 - t504 * t535 + (-pkin(9) * t565 * t581 + t494 * t586) * qJD(4) + (-t494 * t690 + (-pkin(9) * t654 + t443) * t582) * qJD(2)) * MDP(18) + ((t497 + t704) * t586 + (-t498 + t703) * t581) * MDP(13) + (-t535 * t702 + t706) * MDP(12) + (t639 + (-t565 * t695 + (t533 + t654) * t582) * qJD(2)) * MDP(15) + t615 * t664 * MDP(11) + (-t565 * t659 + (t565 * t690 + (-t535 + t663) * t582) * qJD(2)) * MDP(14) + (t571 * t402 + t436 * t537 + t679 * t455 - t607 * t617) * MDP(25) + (t392 * t474 + t684 * t413 - t424 * t681 - t510 * t592) * MDP(31) + (t510 * t375 + t392 * t475 + t685 * t413 - t609 * t681) * MDP(32) + (t375 * t475 - t609 * t685) * MDP(26) + (-t375 * t474 + t424 * t685 + t475 * t592 + t609 * t684) * MDP(27) + (t436 * t536 + t678 * t455 + t617 * t470 - t571 * t591) * MDP(24) + (t402 * t537 - t607 * t679) * MDP(19) + (-t402 * t536 - t470 * t679 + t537 * t591 + t607 * t678) * MDP(20) + (qJD(3) * t504 - t460) * MDP(10) + (-t587 * t721 + t720) * t590 + (t716 * MDP(25) + t678 * MDP(22) + (t554 * t657 + (-qJD(5) * t553 + t727) * t580 + t719) * MDP(24) - t679 * MDP(21)) * t557 + ((t579 * t613 + t584 * t614) * MDP(31) + t684 * MDP(29) + (-t579 * t614 + t584 * t613) * MDP(32) - t685 * MDP(28)) * t552 + ((-qJD(3) * t674 + t391) * MDP(25) + ((t456 * t584 - t457 * t579) * qJD(3) - t371) * MDP(31) + (-qJD(3) * t474 - t424) * MDP(29) + (-(t456 * t579 + t457 * t584) * qJD(3) + t372) * MDP(32) + (qJD(3) * t475 + t609) * MDP(28) + (-qJD(3) * t536 + t470) * MDP(22) + (qJD(3) * t621 - t390) * MDP(24) + (qJD(3) * t537 + t607) * MDP(21) - t548 * MDP(10) + t557 * MDP(23) + t552 * MDP(30) + t565 * MDP(16)) * t665; (t497 - t704) * MDP(14) + (t603 * t636 + (t383 * t584 - t384 * t579) * t552 + t449 * t424 - (-t579 * t585 - t697) * t650 + (t552 * t604 - t372) * qJD(6) + t726) * MDP(31) + (-t683 * t557 + (t535 * t607 + t557 * t657 - t580 * t636) * pkin(4) + t725) * MDP(25) + (-t604 * t636 - (t383 * t579 + t384 * t584) * t552 + t449 * t609 + (t584 * t585 - t698) * t650 + (t552 * t603 - t692) * qJD(6) + t742) * MDP(32) + (t627 * t557 + (-t470 * t535 + t557 * t658 + t585 * t636) * pkin(4) + t724) * MDP(24) + (-t442 * t565 + t494 * t533 - t599) * MDP(18) + (-t443 * t565 - t494 * t535 + t594) * MDP(17) + (-t498 - t703) * MDP(15) + (-t533 ^ 2 + t535 ^ 2) * MDP(13) + t535 * t533 * MDP(12) + qJD(2) * t633 + t741; (-t391 * t557 + t724) * MDP(24) + (-t390 * t557 + t725) * MDP(25) + ((-t381 * t579 - t691) * t552 - t372 * qJD(6) + (-t424 * t607 + t552 * t656 + t584 * t636) * pkin(5) + t726) * MDP(31) + (-t731 + t380 + (t382 * t552 - t369) * t579 + (-t381 * t552 - t620) * t584 + (t552 * t655 - t579 * t636 - t607 * t609) * pkin(5)) * MDP(32) + t741; (t649 + t730) * MDP(28) + (-t628 + t729) * MDP(29) + (-t372 * t552 + t726) * MDP(31) + (-t371 * t552 + t742) * MDP(32) + (MDP(28) * t708 + MDP(29) * t609 - MDP(31) * t372 - MDP(32) * t692) * qJD(6) + t743;];
tauc  = t1;
