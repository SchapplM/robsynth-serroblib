% Calculate vector of inverse dynamics joint torques for
% S6PRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:40:09
% EndTime: 2019-03-09 00:40:21
% DurationCPUTime: 8.89s
% Computational Cost: add. (7613->526), mult. (17610->709), div. (0->0), fcn. (14381->18), ass. (0->237)
t589 = cos(qJ(6));
t678 = qJD(6) * t589;
t586 = sin(qJ(4));
t587 = sin(qJ(3));
t685 = qJD(2) * t587;
t663 = t586 * t685;
t591 = cos(qJ(4));
t592 = cos(qJ(3));
t684 = qJD(2) * t592;
t664 = t591 * t684;
t523 = -t663 + t664;
t524 = -t586 * t684 - t591 * t685;
t585 = sin(qJ(5));
t590 = cos(qJ(5));
t468 = t523 * t590 + t524 * t585;
t743 = t468 * t589;
t747 = t678 - t743;
t530 = t586 * t592 + t587 * t591;
t726 = pkin(8) + pkin(9);
t669 = qJD(3) * t726;
t535 = t587 * t669;
t536 = t592 * t669;
t593 = cos(qJ(2));
t582 = sin(pkin(6));
t687 = qJD(1) * t582;
t667 = t593 * t687;
t541 = t726 * t587;
t542 = t726 * t592;
t692 = -t541 * t586 + t542 * t591;
t746 = -qJD(4) * t692 + t530 * t667 + t535 * t586 - t536 * t591;
t529 = t586 * t587 - t591 * t592;
t682 = qJD(4) * t591;
t683 = qJD(4) * t586;
t745 = -t529 * t667 + t535 * t591 + t536 * t586 + t541 * t682 + t542 * t683;
t622 = t523 * t585 - t524 * t590;
t577 = qJD(3) + qJD(4);
t673 = qJDD(2) * t592;
t675 = qJD(2) * qJD(3);
t660 = t592 * t675;
t674 = qJDD(2) * t587;
t739 = t660 + t674;
t451 = qJD(4) * t664 - t577 * t663 + t586 * t673 + t591 * t739;
t493 = t577 * t530;
t631 = t586 * t674 - t591 * t673;
t452 = qJD(2) * t493 + t631;
t680 = qJD(5) * t590;
t681 = qJD(5) * t585;
t415 = t451 * t590 - t452 * t585 + t523 * t680 + t524 * t681;
t576 = qJDD(3) + qJDD(4);
t571 = qJDD(5) + t576;
t572 = qJD(5) + t577;
t584 = sin(qJ(6));
t670 = t415 * t589 + t571 * t584 + t572 * t678;
t679 = qJD(6) * t584;
t396 = -t622 * t679 + t670;
t394 = t396 * t584;
t395 = t396 * t589;
t455 = t572 * t584 + t589 * t622;
t651 = t415 * t584 - t571 * t589;
t397 = qJD(6) * t455 + t651;
t416 = qJD(5) * t622 + t451 * t585 + t452 * t590;
t414 = qJDD(6) + t416;
t412 = t589 * t414;
t709 = t622 * t584;
t453 = -t572 * t589 + t709;
t677 = -qJD(6) + t468;
t411 = t584 * t414;
t696 = -t677 * t678 + t411;
t742 = t677 * t584;
t744 = t571 * MDP(23) - t416 * MDP(22) - t468 ^ 2 * MDP(20) + (-t468 * t572 + t415) * MDP(21) + (-MDP(19) * t468 + MDP(20) * t622 + MDP(22) * t572 + MDP(30) * t677) * t622 + (t455 * t747 + t394) * MDP(26) + (-t455 * t622 + t677 * t743 + t696) * MDP(28) + (t453 * t622 - t677 * t742 + t412) * MDP(29) + (-t584 * t397 - t453 * t747 + t455 * t742 + t395) * MDP(27);
t515 = t524 * pkin(10);
t588 = sin(qJ(2));
t668 = t588 * t687;
t656 = qJD(2) * t726 + t668;
t583 = cos(pkin(6));
t686 = qJD(1) * t583;
t497 = t587 * t686 + t592 * t656;
t488 = t586 * t497;
t496 = -t587 * t656 + t592 * t686;
t719 = qJD(3) * pkin(3);
t491 = t496 + t719;
t648 = t491 * t591 - t488;
t434 = t515 + t648;
t426 = pkin(4) * t577 + t434;
t490 = t591 * t497;
t623 = -t491 * t586 - t490;
t724 = pkin(10) * t523;
t435 = -t623 + t724;
t713 = t435 * t585;
t404 = t426 * t590 - t713;
t402 = -pkin(5) * t572 - t404;
t716 = t402 * t468;
t672 = t583 * qJDD(1);
t550 = t592 * t672;
t676 = qJD(1) * qJD(2);
t505 = qJDD(2) * pkin(8) + (qJDD(1) * t588 + t593 * t676) * t582;
t655 = pkin(9) * qJDD(2) + t505;
t440 = qJDD(3) * pkin(3) - qJD(3) * t497 - t587 * t655 + t550;
t443 = qJD(3) * t496 + t587 * t672 + t592 * t655;
t607 = qJD(4) * t623 + t440 * t591 - t586 * t443;
t388 = pkin(4) * t576 - pkin(10) * t451 + t607;
t729 = -(qJD(4) * t491 + t443) * t591 - t586 * t440 + t497 * t683;
t391 = -pkin(10) * t452 - t729;
t712 = t435 * t590;
t405 = t426 * t585 + t712;
t728 = -qJD(5) * t405 + t388 * t590 - t585 * t391;
t376 = -pkin(5) * t571 - t728;
t718 = cos(pkin(12));
t658 = t718 * t588;
t581 = sin(pkin(12));
t706 = t581 * t593;
t517 = t583 * t658 + t706;
t657 = t718 * t593;
t707 = t581 * t588;
t519 = -t583 * t707 + t657;
t580 = qJ(3) + qJ(4);
t575 = qJ(5) + t580;
t563 = sin(t575);
t564 = cos(t575);
t659 = t582 * t718;
t705 = t582 * t588;
t708 = t581 * t582;
t615 = -g(3) * (-t563 * t705 + t564 * t583) - g(2) * (-t517 * t563 - t564 * t659) - g(1) * (-t519 * t563 + t564 * t708);
t611 = -t376 + t615;
t741 = -pkin(10) * t493 - t745;
t492 = t577 * t529;
t740 = -pkin(10) * t492 - t746;
t614 = -t587 * t719 + t668;
t486 = t529 * t590 + t530 * t585;
t487 = -t529 * t585 + t530 * t590;
t568 = -pkin(3) * t592 - pkin(2);
t507 = pkin(4) * t529 + t568;
t428 = pkin(5) * t486 - pkin(11) * t487 + t507;
t516 = -t583 * t657 + t707;
t518 = t583 * t706 + t658;
t636 = g(1) * t518 + g(2) * t516;
t703 = t582 * t593;
t612 = g(3) * t703 - t636;
t609 = t612 * t564;
t738 = t428 * t414 - t609;
t438 = pkin(5) * t622 - pkin(11) * t468;
t474 = t517 * t564 - t563 * t659;
t476 = t519 * t564 + t563 * t708;
t514 = qJD(2) * t568 - t667;
t478 = -pkin(4) * t523 + t514;
t503 = t563 * t583 + t564 * t705;
t730 = (qJD(5) * t426 + t391) * t590 + t585 * t388 - t435 * t681;
t601 = g(1) * t476 + g(2) * t474 + g(3) * t503 - t478 * t468 - t730;
t633 = pkin(4) * t493 - t614;
t403 = pkin(11) * t572 + t405;
t417 = -pkin(5) * t468 - pkin(11) * t622 + t478;
t626 = t403 * t584 - t417 * t589;
t654 = t402 * t679 + t622 * t626;
t384 = t403 * t589 + t417 * t584;
t613 = t384 * t622 + t402 * t678 - t584 * t611;
t599 = -t478 * t622 + t615 + t728;
t375 = pkin(11) * t571 + t730;
t421 = -qJD(5) * t486 - t492 * t590 - t493 * t585;
t646 = -t541 * t591 - t542 * t586;
t460 = -pkin(10) * t530 + t646;
t461 = -pkin(10) * t529 + t692;
t430 = t460 * t585 + t461 * t590;
t635 = g(1) * t519 + g(2) * t517;
t625 = t460 * t590 - t461 * t585;
t699 = -qJD(5) * t625 + t585 * t740 - t590 * t741;
t732 = -(qJD(6) * t417 + t375) * t486 + t376 * t487 + t402 * t421 - (-qJD(6) * t428 + t699) * t677 - t430 * t414 - g(3) * t705 - t635;
t594 = qJD(3) ^ 2;
t662 = t588 * t676;
t632 = -qJDD(1) * t703 + t582 * t662;
t731 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t594 + t582 * (-g(3) * t593 + t662) - t632 + t636;
t725 = pkin(4) * t524;
t720 = qJD(2) * pkin(2);
t715 = t402 * t487;
t704 = t582 * t592;
t702 = t585 * t586;
t701 = t586 * t590;
t700 = qJDD(1) - g(3);
t698 = qJD(5) * t430 + t585 * t741 + t590 * t740;
t647 = -t496 * t586 - t490;
t436 = t647 - t724;
t693 = t496 * t591 - t488;
t437 = t515 + t693;
t567 = pkin(3) * t591 + pkin(4);
t695 = t436 * t585 + t437 * t590 - t567 * t680 - (-t586 * t681 + (t590 * t591 - t702) * qJD(4)) * pkin(3);
t694 = t436 * t590 - t437 * t585 + t567 * t681 + (t586 * t680 + (t585 * t591 + t701) * qJD(4)) * pkin(3);
t691 = pkin(3) * t701 + t567 * t585;
t578 = t587 ^ 2;
t690 = -t592 ^ 2 + t578;
t665 = qJD(2) * t703;
t661 = t587 * t675;
t423 = t438 - t725;
t513 = pkin(11) + t691;
t569 = pkin(3) * t685;
t640 = qJD(6) * t513 + t423 + t569;
t565 = pkin(4) * t585 + pkin(11);
t639 = qJD(6) * t565 + t423;
t406 = t434 * t585 + t712;
t638 = pkin(4) * t681 - t406;
t407 = t434 * t590 - t713;
t637 = -pkin(4) * t680 + t407;
t422 = qJD(5) * t487 - t492 * t585 + t493 * t590;
t634 = pkin(5) * t422 - pkin(11) * t421 + t633;
t628 = -t414 * t513 - t716;
t627 = -t414 * t565 - t716;
t520 = t583 * t592 - t587 * t705;
t521 = t583 * t587 + t588 * t704;
t463 = t520 * t591 - t521 * t586;
t464 = t520 * t586 + t521 * t591;
t624 = t463 * t590 - t464 * t585;
t433 = t463 * t585 + t464 * t590;
t621 = -pkin(3) * t702 + t567 * t590;
t620 = -g(1) * t581 + g(2) * t718;
t618 = -t677 * t679 - t412;
t617 = t421 * t589 - t487 * t679;
t538 = -t667 - t720;
t610 = -qJD(2) * t538 - t505 + t635;
t479 = pkin(3) * t661 + qJDD(2) * t568 + t632;
t604 = -pkin(8) * qJDD(3) + (t538 + t667 - t720) * qJD(3);
t431 = pkin(4) * t452 + t479;
t573 = sin(t580);
t574 = cos(t580);
t600 = -g(1) * (-t519 * t574 - t573 * t708) - g(2) * (-t517 * t574 + t573 * t659) - g(3) * (-t573 * t583 - t574 * t705) - t514 * t523 + t729;
t598 = -g(1) * (-t519 * t573 + t574 * t708) - g(2) * (-t517 * t573 - t574 * t659) - g(3) * (-t573 * t705 + t574 * t583) + t514 * t524 + t607;
t596 = t524 * t523 * MDP(12) + (-t523 * t577 + t451) * MDP(14) + (-t631 + (-qJD(2) * t530 - t524) * t577) * MDP(15) + (-t523 ^ 2 + t524 ^ 2) * MDP(13) + t576 * MDP(16) + t744;
t595 = qJD(2) ^ 2;
t566 = -pkin(4) * t590 - pkin(5);
t512 = -pkin(5) - t621;
t501 = t569 - t725;
t495 = -qJD(3) * t521 - t587 * t665;
t494 = qJD(3) * t520 + t592 * t665;
t419 = -qJD(4) * t464 - t494 * t586 + t495 * t591;
t418 = qJD(4) * t463 + t494 * t591 + t495 * t586;
t386 = qJD(5) * t433 + t418 * t585 - t419 * t590;
t385 = qJD(5) * t624 + t418 * t590 + t419 * t585;
t379 = pkin(5) * t416 - pkin(11) * t415 + t431;
t378 = t589 * t379;
t1 = [t700 * MDP(1) + (qJD(3) * t495 + qJDD(3) * t520) * MDP(10) + (-qJD(3) * t494 - qJDD(3) * t521) * MDP(11) + (t419 * t577 + t463 * t576) * MDP(17) + (-t418 * t577 - t464 * t576) * MDP(18) + (-t386 * t572 + t571 * t624) * MDP(24) + (-t385 * t572 - t433 * t571) * MDP(25) + (-(-t385 * t584 - t433 * t678) * t677 - t433 * t411 + t386 * t453 - t624 * t397) * MDP(31) + ((t385 * t589 - t433 * t679) * t677 - t433 * t412 + t386 * t455 - t624 * t396) * MDP(32) + ((-qJDD(2) * MDP(4) + (-MDP(10) * t592 + MDP(11) * t587 - MDP(3)) * t595 + (-MDP(17) * t523 - MDP(18) * t524 - MDP(24) * t468 + MDP(25) * t622 - (MDP(31) * t589 - MDP(32) * t584) * t677) * qJD(2)) * t588 + (qJDD(2) * MDP(3) - t595 * MDP(4) + (-t661 + t673) * MDP(10) - t739 * MDP(11) - t452 * MDP(17) - t451 * MDP(18) - t416 * MDP(24) - t415 * MDP(25) + t618 * MDP(31) + t696 * MDP(32)) * t593) * t582; (-t384 * t422 - t625 * t396 + t698 * t455 + (-(-qJD(6) * t403 + t379) * t486 - qJD(6) * t715 - (qJD(6) * t430 - t634) * t677 - t738) * t584 + t732 * t589) * MDP(32) + (t378 * t486 - t626 * t422 - t625 * t397 + t698 * t453 + (-t634 * t677 + (-t403 * t486 + t430 * t677 + t715) * qJD(6) + t738) * t589 + t732 * t584) * MDP(31) + (t414 * t486 - t422 * t677) * MDP(30) + (t396 * t486 + t412 * t487 + t422 * t455 - t617 * t677) * MDP(28) + (-t487 * t411 - t397 * t486 - t422 * t453 - (-t421 * t584 - t487 * t678) * t677) * MDP(29) + (t416 * t507 + t422 * t478 + t431 * t486 - t468 * t633 + t571 * t625 - t572 * t698 - t609) * MDP(24) + (-t415 * t486 - t416 * t487 + t421 * t468 - t422 * t622) * MDP(20) + (-t492 * t577 + t530 * t576) * MDP(14) + (t451 * t530 + t492 * t524) * MDP(12) + (-t451 * t529 - t452 * t530 - t492 * t523 + t493 * t524) * MDP(13) + (t568 * t451 + t479 * t530 - t514 * t492 + t614 * t524 + t612 * t573 - t692 * t576 + t577 * t745) * MDP(18) + (t568 * t452 + t479 * t529 + t514 * t493 + t614 * t523 - t612 * t574 + t646 * t576 + t577 * t746) * MDP(17) + (t415 * t487 + t421 * t622) * MDP(19) + (t415 * t507 + t421 * t478 - t430 * t571 + t431 * t487 + t563 * t612 + t572 * t699 + t622 * t633) * MDP(25) + (-t422 * t572 - t486 * t571) * MDP(22) + (t421 * t572 + t487 * t571) * MDP(21) + qJDD(2) * MDP(2) + (qJDD(3) * t587 + t592 * t594) * MDP(7) + (qJDD(3) * t592 - t587 * t594) * MDP(8) + (-t493 * t577 - t529 * t576) * MDP(15) + (t604 * t587 + t592 * t731) * MDP(10) + (-t587 * t731 + t604 * t592) * MDP(11) + (qJDD(2) * t578 + 0.2e1 * t587 * t660) * MDP(5) + 0.2e1 * (t587 * t673 - t675 * t690) * MDP(6) + (t700 * t703 + t636) * MDP(3) + (-t700 * t705 + t635) * MDP(4) + (t395 * t487 + t455 * t617) * MDP(26) + ((-t453 * t589 - t455 * t584) * t421 + (-t394 - t397 * t589 + (t453 * t584 - t455 * t589) * qJD(6)) * t487) * MDP(27); t596 + (-g(3) * t520 + t587 * t610 + t620 * t704 + t550) * MDP(10) + qJDD(3) * MDP(9) + (t512 * t396 + t628 * t589 + t694 * t455 - (t584 * t640 + t589 * t695) * t677 + t613) * MDP(32) + (t512 * t397 + t694 * t453 + (-t677 * t695 + t628) * t584 + (t640 * t677 + t611) * t589 + t654) * MDP(31) + (g(3) * t521 + (-t582 * t620 - t672) * t587 + t610 * t592) * MDP(11) + (t693 * t577 + (t524 * t685 - t576 * t586 - t577 * t682) * pkin(3) + t600) * MDP(18) + (-t647 * t577 + (t523 * t685 + t576 * t591 - t577 * t683) * pkin(3) + t598) * MDP(17) + (t468 * t501 + t571 * t621 - t572 * t694 + t599) * MDP(24) + (-t501 * t622 - t571 * t691 + t572 * t695 + t601) * MDP(25) + MDP(7) * t674 + MDP(8) * t673 + (-MDP(5) * t587 * t592 + MDP(6) * t690) * t595; t596 + (t566 * t397 + t638 * t453 + (-t637 * t677 + t627) * t584 + (t639 * t677 + t611) * t589 + t654) * MDP(31) + (t406 * t572 + (-t468 * t524 + t571 * t590 - t572 * t681) * pkin(4) + t599) * MDP(24) + (t577 * t648 + t600) * MDP(18) + (t566 * t396 + t627 * t589 + t638 * t455 - (t584 * t639 + t589 * t637) * t677 + t613) * MDP(32) + (t407 * t572 + (t524 * t622 - t571 * t585 - t572 * t680) * pkin(4) + t601) * MDP(25) + (-t577 * t623 + t598) * MDP(17); (t405 * t572 + t599) * MDP(24) + (t404 * t572 + t601) * MDP(25) + (-pkin(5) * t397 - t405 * t453 + (-pkin(11) * t414 - t404 * t677 - t716) * t584 + (-(-pkin(11) * qJD(6) - t438) * t677 + t611) * t589 + t654) * MDP(31) + (-pkin(5) * t396 - (t404 * t589 + t438 * t584) * t677 - t405 * t455 - t402 * t743 + t618 * pkin(11) + t613) * MDP(32) + t744; t455 * t453 * MDP(26) + (-t453 ^ 2 + t455 ^ 2) * MDP(27) + (-t453 * t677 + t670) * MDP(28) + (-t455 * t677 - t651) * MDP(29) + t414 * MDP(30) + (-t584 * t375 + t378 - t384 * t677 - t402 * t455 - g(1) * (-t476 * t584 + t518 * t589) - g(2) * (-t474 * t584 + t516 * t589) - g(3) * (-t503 * t584 - t589 * t703)) * MDP(31) + (-t589 * t375 - t584 * t379 + t626 * t677 + t402 * t453 - g(1) * (-t476 * t589 - t518 * t584) - g(2) * (-t474 * t589 - t516 * t584) - g(3) * (-t503 * t589 + t584 * t703)) * MDP(32) + (-MDP(28) * t709 - MDP(29) * t455 - MDP(31) * t384 + MDP(32) * t626) * qJD(6);];
tau  = t1;
