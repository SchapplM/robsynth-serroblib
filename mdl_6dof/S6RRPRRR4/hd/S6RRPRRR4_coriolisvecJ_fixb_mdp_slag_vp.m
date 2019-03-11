% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRPRRR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:35:08
% EndTime: 2019-03-09 13:35:27
% DurationCPUTime: 11.74s
% Computational Cost: add. (12740->558), mult. (37349->758), div. (0->0), fcn. (31146->12), ass. (0->252)
t608 = cos(pkin(12));
t617 = cos(qJ(2));
t607 = sin(pkin(6));
t701 = qJD(1) * t607;
t679 = t617 * t701;
t589 = t608 * t679;
t606 = sin(pkin(12));
t613 = sin(qJ(2));
t680 = t613 * t701;
t573 = -t606 * t680 + t589;
t570 = qJD(4) - t573;
t567 = qJD(5) + t570;
t614 = cos(qJ(6));
t693 = qJD(6) * t614;
t634 = t606 * t617 + t608 * t613;
t576 = t634 * t701;
t609 = cos(pkin(6));
t700 = qJD(1) * t609;
t596 = qJD(2) + t700;
t612 = sin(qJ(4));
t616 = cos(qJ(4));
t535 = -t576 * t612 + t596 * t616;
t611 = sin(qJ(5));
t615 = cos(qJ(5));
t636 = -t576 * t616 - t596 * t612;
t472 = -t535 * t615 - t611 * t636;
t762 = t472 * t614;
t769 = t693 + t762;
t739 = pkin(1) * t613;
t689 = t609 * t739;
t716 = t607 * t617;
t737 = pkin(8) + qJ(3);
t571 = t716 * t737 + t689;
t560 = t571 * qJD(1);
t544 = t606 * t560;
t738 = pkin(1) * t617;
t688 = t609 * t738;
t594 = qJD(1) * t688;
t677 = t737 * t613;
t658 = t607 * t677;
t559 = -qJD(1) * t658 + t594;
t504 = t559 * t608 - t544;
t519 = pkin(2) * t680 + pkin(3) * t576 - pkin(9) * t573;
t508 = t616 * t519;
t598 = pkin(2) * t606 + pkin(9);
t736 = pkin(10) + t598;
t675 = qJD(4) * t736;
t767 = pkin(4) * t576 - t504 * t612 + t508 + (-pkin(10) * t573 + t675) * t616;
t706 = t504 * t616 + t519 * t612;
t720 = t573 * t612;
t766 = -pkin(10) * t720 + t612 * t675 + t706;
t698 = qJD(4) * t612;
t765 = t698 - t720;
t637 = t535 * t611 - t615 * t636;
t699 = qJD(2) * t607;
t678 = t613 * t699;
t656 = qJD(1) * t678;
t569 = qJD(2) * t589 - t606 * t656;
t697 = qJD(4) * t616;
t490 = t569 * t616 - t576 * t698 + t596 * t697;
t491 = -qJD(4) * t636 + t569 * t612;
t695 = qJD(5) * t615;
t696 = qJD(5) * t611;
t428 = t490 * t615 - t491 * t611 + t535 * t695 + t636 * t696;
t579 = t634 * t607;
t574 = qJD(2) * t579;
t568 = qJD(1) * t574;
t610 = sin(qJ(6));
t681 = t428 * t614 + t567 * t693 + t568 * t610;
t694 = qJD(6) * t610;
t411 = -t637 * t694 + t681;
t409 = t411 * t610;
t459 = t567 * t610 + t614 * t637;
t671 = t428 * t610 - t568 * t614;
t412 = qJD(6) * t459 + t671;
t429 = qJD(5) * t637 + t490 * t611 + t491 * t615;
t427 = t614 * t429;
t728 = t637 * t610;
t457 = -t567 * t614 + t728;
t690 = -qJD(6) - t472;
t425 = t610 * t429;
t711 = -t690 * t693 + t425;
t761 = t690 * t610;
t764 = t568 * MDP(24) - t429 * MDP(23) - t472 ^ 2 * MDP(21) + (t472 * t567 + t428) * MDP(22) + (MDP(20) * t472 + MDP(21) * t637 + MDP(23) * t567 + MDP(31) * t690) * t637 + (t459 * t769 + t409) * MDP(27) + (-t459 * t637 - t690 * t762 + t711) * MDP(29) + (t457 * t637 - t690 * t761 + t427) * MDP(30) + (t411 * t614 - t610 * t412 - t457 * t769 + t459 * t761) * MDP(28);
t539 = pkin(2) * t596 + t559;
t715 = t608 * t560;
t487 = t539 * t606 + t715;
t480 = pkin(9) * t596 + t487;
t652 = (-pkin(2) * t617 - pkin(1)) * t607;
t633 = qJD(1) * t652;
t582 = qJD(3) + t633;
t499 = -pkin(3) * t573 - pkin(9) * t576 + t582;
t450 = -t480 * t612 + t499 * t616;
t437 = pkin(10) * t636 + t450;
t431 = pkin(4) * t570 + t437;
t451 = t480 * t616 + t499 * t612;
t438 = pkin(10) * t535 + t451;
t734 = t438 * t611;
t400 = t431 * t615 - t734;
t396 = -pkin(5) * t567 - t400;
t763 = t396 * t472;
t585 = t611 * t612 - t615 * t616;
t760 = t567 * t585;
t586 = t611 * t616 + t612 * t615;
t704 = t567 * t586;
t759 = t412 * t585 + t457 * t704;
t444 = pkin(5) * t637 + pkin(11) * t472;
t486 = t539 * t608 - t544;
t479 = -pkin(3) * t596 - t486;
t456 = -pkin(4) * t535 + t479;
t592 = qJD(2) * t594;
t622 = (-qJD(2) * t677 + qJD(3) * t617) * t607;
t527 = qJD(1) * t622 + t592;
t717 = t607 * t613;
t541 = -qJD(2) * t571 - qJD(3) * t717;
t528 = t541 * qJD(1);
t465 = t527 * t608 + t528 * t606;
t591 = pkin(2) * t656;
t500 = pkin(3) * t568 - pkin(9) * t569 + t591;
t670 = -t465 * t612 + t500 * t616;
t620 = -qJD(4) * t451 + t670;
t402 = pkin(4) * t568 - pkin(10) * t490 + t620;
t629 = t465 * t616 - t480 * t698 + t499 * t697 + t500 * t612;
t407 = -pkin(10) * t491 + t629;
t659 = -t402 * t611 - t407 * t615 - t431 * t695 + t438 * t696;
t758 = t456 * t472 + t659;
t755 = pkin(2) * t678;
t753 = MDP(4) * t613;
t752 = MDP(5) * (t613 ^ 2 - t617 ^ 2);
t550 = t579 * t616 + t609 * t612;
t714 = t608 * t617;
t578 = t606 * t717 - t607 * t714;
t557 = (pkin(2) + t738) * t609 - t658;
t514 = t557 * t606 + t571 * t608;
t502 = pkin(9) * t609 + t514;
t522 = pkin(3) * t578 - pkin(9) * t579 + t652;
t665 = -t502 * t612 + t522 * t616;
t442 = pkin(4) * t578 - pkin(10) * t550 + t665;
t549 = t579 * t612 - t609 * t616;
t707 = t502 * t616 + t522 * t612;
t448 = -pkin(10) * t549 + t707;
t750 = t442 * t611 + t448 * t615;
t503 = t559 * t606 + t715;
t655 = pkin(4) * t765 - t503;
t583 = t736 * t612;
t584 = t736 * t616;
t635 = -t583 * t615 - t584 * t611;
t749 = -qJD(5) * t635 + t611 * t767 + t615 * t766;
t538 = -t583 * t611 + t584 * t615;
t748 = -qJD(5) * t538 + t611 * t766 - t615 * t767;
t733 = t438 * t615;
t401 = t431 * t611 + t733;
t672 = -t402 * t615 + t407 * t611;
t621 = -qJD(5) * t401 - t672;
t386 = -pkin(5) * t568 - t621;
t397 = pkin(11) * t567 + t401;
t423 = pkin(5) * t472 - pkin(11) * t637 + t456;
t645 = t397 * t610 - t423 * t614;
t746 = -t386 * t614 + t396 * t694 + t637 * t645;
t384 = t386 * t610;
t390 = t397 * t614 + t423 * t610;
t745 = t390 * t637 + t396 * t693 + t384;
t743 = -t456 * t637 + t621;
t741 = t567 * t760 - t568 * t586;
t732 = t472 * t576;
t729 = t637 * t576;
t727 = t535 * t570;
t726 = t535 * t576;
t725 = t636 * t570;
t724 = t636 * t576;
t722 = t568 * t612;
t719 = t586 * t614;
t603 = t607 ^ 2;
t618 = qJD(1) ^ 2;
t718 = t603 * t618;
t708 = pkin(5) * t576 - t748;
t687 = t603 * t739;
t685 = t586 * t425;
t684 = t586 * t427;
t683 = t617 * t718;
t599 = -pkin(2) * t608 - pkin(3);
t676 = qJD(1) * qJD(2) * t603;
t385 = pkin(11) * t568 - t659;
t464 = t527 * t606 - t528 * t608;
t446 = pkin(4) * t491 + t464;
t392 = pkin(5) * t429 - pkin(11) * t428 + t446;
t674 = -t385 * t610 + t392 * t614;
t669 = -t576 * t614 + t610 * t760;
t668 = t576 * t610 + t614 * t760;
t595 = qJD(2) * t688;
t540 = t595 + t622;
t482 = t540 * t608 + t541 * t606;
t575 = (-t606 * t613 + t714) * t699;
t520 = pkin(3) * t574 - pkin(9) * t575 + t755;
t667 = -t482 * t612 + t520 * t616;
t481 = t540 * t606 - t541 * t608;
t513 = t557 * t608 - t571 * t606;
t664 = t570 * t616;
t601 = pkin(4) * t611 + pkin(11);
t660 = -pkin(4) * t636 + qJD(6) * t601 + t444;
t657 = t617 * t676;
t404 = t437 * t611 + t733;
t654 = pkin(4) * t696 - t404;
t405 = t437 * t615 - t734;
t653 = -pkin(4) * t695 + t405;
t651 = t411 * t585 + t459 * t704;
t650 = -t567 * t704 - t585 * t568;
t590 = -pkin(4) * t616 + t599;
t532 = pkin(5) * t585 - pkin(11) * t586 + t590;
t649 = pkin(11) * t576 - qJD(6) * t532 + t749;
t648 = -pkin(5) * t704 - pkin(11) * t760 + qJD(6) * t538 - t655;
t647 = t385 * t614 + t392 * t610;
t646 = -t429 * t601 + t763;
t414 = pkin(11) * t578 + t750;
t501 = -pkin(3) * t609 - t513;
t460 = pkin(4) * t549 + t501;
t495 = t549 * t615 + t550 * t611;
t496 = -t549 * t611 + t550 * t615;
t424 = pkin(5) * t495 - pkin(11) * t496 + t460;
t644 = t414 * t614 + t424 * t610;
t643 = -t414 * t610 + t424 * t614;
t512 = -qJD(4) * t549 + t575 * t616;
t418 = pkin(4) * t574 - pkin(10) * t512 - qJD(4) * t707 + t667;
t511 = qJD(4) * t550 + t575 * t612;
t628 = t482 * t616 - t502 * t698 + t520 * t612 + t522 * t697;
t420 = -pkin(10) * t511 + t628;
t642 = t418 * t615 - t420 * t611;
t641 = t442 * t615 - t448 * t611;
t462 = t496 * t614 + t578 * t610;
t461 = t496 * t610 - t578 * t614;
t632 = t568 * t616 - t570 * t765;
t453 = pkin(4) * t511 + t481;
t631 = -pkin(8) * t716 - t689;
t630 = -pkin(8) * t656 + t592;
t627 = t418 * t611 + t420 * t615 + t442 * t695 - t448 * t696;
t626 = t479 * t570 - t568 * t598;
t625 = t586 * t693 - t669;
t624 = t586 * t694 + t668;
t623 = t631 * t596;
t602 = -pkin(4) * t615 - pkin(5);
t524 = t568 * t578;
t436 = qJD(5) * t496 + t511 * t615 + t512 * t611;
t435 = -qJD(5) * t495 - t511 * t611 + t512 * t615;
t422 = qJD(6) * t462 + t435 * t610 - t574 * t614;
t421 = -qJD(6) * t461 + t435 * t614 + t574 * t610;
t413 = -pkin(5) * t578 - t641;
t393 = pkin(5) * t436 - pkin(11) * t435 + t453;
t388 = -pkin(5) * t574 + qJD(5) * t750 - t642;
t387 = pkin(11) * t574 + t627;
t383 = -qJD(6) * t390 + t674;
t382 = -qJD(6) * t645 + t647;
t1 = [(t623 + (t609 * t631 - 0.2e1 * t687) * qJD(1)) * qJD(2) * MDP(9) + (MDP(6) * t617 * t699 - MDP(7) * t678) * (t596 + t700) + ((qJD(6) * t643 + t387 * t614 + t393 * t610) * t690 - t644 * t429 - t382 * t495 - t390 * t436 + t388 * t459 + t413 * t411 + t386 * t462 + t396 * t421) * MDP(33) + (-t412 * t495 + t422 * t690 - t429 * t461 - t436 * t457) * MDP(30) + (t411 * t495 - t421 * t690 + t429 * t462 + t436 * t459) * MDP(29) + (-t491 * t578 - t511 * t570 + t535 * t574 - t549 * t568) * MDP(16) + (-t429 * t578 - t436 * t567 - t472 * t574 - t495 * t568) * MDP(23) + (t464 * t579 - t465 * t578 + t481 * t576 + t482 * t573 - t486 * t575 - t487 * t574 - t513 * t569 - t514 * t568) * MDP(11) + (t570 * t574 + t524) * MDP(17) + (t567 * t574 + t524) * MDP(24) + (t428 * t578 + t435 * t567 + t496 * t568 + t574 * t637) * MDP(22) + (-t428 * t495 - t429 * t496 - t435 * t472 - t436 * t637) * MDP(21) + (t428 * t496 + t435 * t637) * MDP(20) + (-t411 * t461 - t412 * t462 - t421 * t457 - t422 * t459) * MDP(28) + (t411 * t462 + t421 * t459) * MDP(27) + (t429 * t495 - t436 * t690) * MDP(31) + (-(-qJD(6) * t644 - t387 * t610 + t393 * t614) * t690 + t643 * t429 + t383 * t495 - t645 * t436 + t388 * t457 + t413 * t412 + t386 * t461 + t396 * t422) * MDP(32) + (-t451 * t574 + t464 * t550 + t479 * t512 - t481 * t636 + t501 * t490 - t568 * t707 - t570 * t628 - t578 * t629) * MDP(19) + (t490 * t578 + t512 * t570 + t550 * t568 - t574 * t636) * MDP(15) + (-t490 * t549 - t491 * t550 + t511 * t636 + t512 * t535) * MDP(14) + (t490 * t550 - t512 * t636) * MDP(13) + (t642 * t567 + t641 * t568 - t672 * t578 + t400 * t574 + t453 * t472 + t460 * t429 + t446 * t495 + t456 * t436 + (-t401 * t578 - t567 * t750) * qJD(5)) * MDP(25) + (-t401 * t574 + t460 * t428 + t456 * t435 + t446 * t496 + t453 * t637 - t567 * t627 - t568 * t750 + t578 * t659) * MDP(26) + (-0.2e1 * pkin(1) * t657 - (-pkin(8) * t678 + t595) * t596 - t630 * t609) * MDP(10) - 0.2e1 * t676 * t752 + 0.2e1 * t657 * t753 + (-t464 * t513 + t465 * t514 - t486 * t481 + t487 * t482 + (t582 + t633) * t755) * MDP(12) + (t667 * t570 + t665 * t568 + t670 * t578 + t450 * t574 - t481 * t535 + t501 * t491 + t464 * t549 + t479 * t511 + (-t451 * t578 - t570 * t707) * qJD(4)) * MDP(18); (t590 * t428 + t446 * t586 - t456 * t760 - t538 * t568 + t567 * t749 + t655 * t637) * MDP(26) + (t428 * t586 - t637 * t760) * MDP(20) + (t669 * t459 + t668 * t457 + (-t409 - t412 * t614 + (t457 * t610 - t459 * t614) * qJD(6)) * t586) * MDP(28) + ((t490 + t727) * t616 + (-t491 + t725) * t612) * MDP(14) + (t632 - t726) * MDP(16) + (t650 + t732) * MDP(23) + (t625 * t690 - t685 - t759) * MDP(30) + (-t729 - t741) * MDP(22) + ((t486 - t504) * t573 + (-t568 * t606 - t569 * t608) * pkin(2)) * MDP(11) + (-t464 * t616 + t599 * t491 + t503 * t535 + (-t598 * t697 - t508) * t570 + (t504 * t570 + t626) * t612) * MDP(18) - ((-t487 + t503) * MDP(11) + t400 * MDP(25) - t451 * MDP(19) - t401 * MDP(26) + t450 * MDP(18) + t567 * MDP(24) + t570 * MDP(17)) * t576 + (t618 * t687 + (qJD(2) * t631 - t623) * qJD(1)) * MDP(9) + (t411 * t719 - t459 * t624) * MDP(27) + (t624 * t690 + t651 + t684) * MDP(29) + (pkin(1) * t683 + (-pkin(8) * t680 + t594) * t596 - t630) * MDP(10) + (t486 * t503 - t487 * t504 + (-t464 * t608 + t465 * t606 - t582 * t680) * pkin(2)) * MDP(12) + (MDP(6) * t679 - MDP(7) * t680) * (qJD(2) - t596) + (t464 * t612 + t599 * t490 + t503 * t636 + (t598 * t698 + t706) * t570 + t626 * t616) * MDP(19) + (t490 * t612 - t636 * t664) * MDP(13) + (t429 * t585 - t690 * t704) * MDP(31) + (t590 * t429 + t446 * t585 + t704 * t456 + t655 * t472 + t567 * t748 + t568 * t635) * MDP(25) + (-(t532 * t610 + t538 * t614) * t429 - t382 * t585 - t635 * t411 + t386 * t719 - (t610 * t648 + t614 * t649) * t690 + t708 * t459 - t704 * t390 - t624 * t396) * MDP(33) + ((t532 * t614 - t538 * t610) * t429 + t383 * t585 - t635 * t412 + t586 * t384 - (t610 * t649 - t614 * t648) * t690 + t708 * t457 - t704 * t645 + t625 * t396) * MDP(32) + (-t428 * t585 - t429 * t586 + t472 * t760 - t637 * t704) * MDP(21) + (t570 * t664 + t722 + t724) * MDP(15) - t683 * t753 + t718 * t752; (-t573 ^ 2 - t576 ^ 2) * MDP(11) + (t486 * t576 - t487 * t573 + t591) * MDP(12) + (t632 + t726) * MDP(18) + (-t570 ^ 2 * t616 - t722 + t724) * MDP(19) + (t650 - t732) * MDP(25) + (-t729 + t741) * MDP(26) + (-t685 + t759) * MDP(32) + (t651 - t684) * MDP(33) - (-MDP(32) * t625 + MDP(33) * t624) * t690; (t450 * t570 - t479 * t535 - t629) * MDP(19) + (-t491 - t725) * MDP(16) + (t451 * t570 + t479 * t636 + t620) * MDP(18) + (t404 * t567 + (t472 * t636 - t567 * t696 + t568 * t615) * pkin(4) + t743) * MDP(25) + (t405 * t567 + (-t567 * t695 - t568 * t611 + t636 * t637) * pkin(4) + t758) * MDP(26) + (t602 * t412 + t646 * t610 + t654 * t457 - (t610 * t653 - t614 * t660) * t690 + t746) * MDP(32) + (t602 * t411 + t646 * t614 + t654 * t459 - (t610 * t660 + t614 * t653) * t690 + t745) * MDP(33) + t568 * MDP(17) + (-t535 ^ 2 + t636 ^ 2) * MDP(14) + (t490 - t727) * MDP(15) + t636 * t535 * MDP(13) + t764; (t401 * t567 + t743) * MDP(25) + (t400 * t567 + t758) * MDP(26) + (-pkin(5) * t412 + (-t400 * t610 + t444 * t614) * t690 - t401 * t457 + t610 * t763 - t711 * pkin(11) + t746) * MDP(32) + (-pkin(5) * t411 - (t400 * t614 + t444 * t610) * t690 - t401 * t459 + t396 * t762 + (-t690 * t694 - t427) * pkin(11) + t745) * MDP(33) + t764; t459 * t457 * MDP(27) + (-t457 ^ 2 + t459 ^ 2) * MDP(28) + (-t457 * t690 + t681) * MDP(29) + (-t459 * t690 - t671) * MDP(30) + t429 * MDP(31) + (-t390 * t690 - t396 * t459 + t674) * MDP(32) + (t396 * t457 + t645 * t690 - t647) * MDP(33) + (-MDP(29) * t728 - MDP(30) * t459 - MDP(32) * t390 + MDP(33) * t645) * qJD(6);];
tauc  = t1;
