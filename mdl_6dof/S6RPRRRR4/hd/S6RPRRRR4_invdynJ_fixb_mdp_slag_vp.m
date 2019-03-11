% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:06:15
% EndTime: 2019-03-09 07:06:27
% DurationCPUTime: 9.45s
% Computational Cost: add. (12230->507), mult. (30789->639), div. (0->0), fcn. (25752->18), ass. (0->238)
t586 = sin(qJ(5));
t591 = cos(qJ(5));
t583 = sin(pkin(11));
t588 = sin(qJ(3));
t584 = cos(pkin(11));
t593 = cos(qJ(3));
t687 = t593 * t584;
t530 = t583 * t588 - t687;
t520 = t530 * qJD(1);
t695 = t583 * t593;
t531 = t584 * t588 + t695;
t521 = t531 * qJD(1);
t587 = sin(qJ(4));
t592 = cos(qJ(4));
t620 = -t520 * t587 + t592 * t521;
t621 = -t520 * t592 - t587 * t521;
t452 = t586 * t620 - t591 * t621;
t590 = cos(qJ(6));
t667 = qJD(6) * t590;
t762 = t452 * t590 + t667;
t674 = qJD(1) * t588;
t655 = t583 * t674;
t661 = qJDD(1) * t593;
t662 = qJDD(1) * t588;
t673 = qJD(3) * t593;
t656 = t583 * t661 + (qJD(1) * t673 + t662) * t584;
t499 = -qJD(3) * t655 + t656;
t523 = t531 * qJD(3);
t556 = t584 * t661;
t628 = -t583 * t662 + t556;
t500 = qJD(1) * t523 - t628;
t671 = qJD(4) * t592;
t672 = qJD(4) * t587;
t440 = t592 * t499 - t587 * t500 - t520 * t671 - t521 * t672;
t441 = qJD(4) * t620 + t499 * t587 + t592 * t500;
t720 = t586 * t621 + t591 * t620;
t402 = qJD(5) * t720 + t440 * t586 + t591 * t441;
t400 = qJDD(6) + t402;
t398 = t590 * t400;
t582 = qJD(3) + qJD(4);
t576 = qJD(5) + t582;
t585 = sin(qJ(6));
t446 = t576 * t585 + t590 * t720;
t444 = -t590 * t576 + t585 * t720;
t706 = t444 * t720;
t744 = qJD(6) + t452;
t669 = qJD(5) * t591;
t670 = qJD(5) * t586;
t401 = t591 * t440 - t586 * t441 - t620 * t670 + t621 * t669;
t578 = qJDD(3) + qJDD(4);
t573 = qJDD(5) + t578;
t658 = t590 * t401 + t585 * t573 + t576 * t667;
t668 = qJD(6) * t585;
t395 = -t668 * t720 + t658;
t393 = t395 * t585;
t397 = t585 * t400;
t700 = t452 * t576;
t703 = t720 * t576;
t705 = t446 * t720;
t751 = (t401 + t700) * MDP(24) + (-t402 + t703) * MDP(25) - t452 ^ 2 * MDP(23) + (MDP(22) * t452 + MDP(23) * t720 - MDP(33) * t744) * t720 + (t762 * t446 + t393) * MDP(29) + (t762 * t744 + t397 - t705) * MDP(31);
t394 = t395 * t590;
t552 = t590 * t573;
t396 = qJD(6) * t446 + t401 * t585 - t552;
t752 = -t585 * t396 - t762 * t444 + t394;
t759 = t744 * t585;
t761 = t751 + (-t446 * t759 + t752) * MDP(30) + (-t744 * t759 + t398 + t706) * MDP(32);
t581 = pkin(11) + qJ(3);
t577 = qJ(4) + t581;
t567 = qJ(5) + t577;
t560 = sin(t567);
t589 = sin(qJ(1));
t594 = cos(qJ(1));
t630 = g(1) * t594 + g(2) * t589;
t760 = t630 * t560;
t712 = pkin(7) + qJ(2);
t545 = t712 * t583;
t532 = qJD(1) * t545;
t546 = t712 * t584;
t533 = qJD(1) * t546;
t619 = t532 * t588 - t533 * t593;
t477 = -pkin(8) * t520 - t619;
t472 = t587 * t477;
t732 = -t593 * t532 - t533 * t588;
t476 = -pkin(8) * t521 + t732;
t475 = qJD(3) * pkin(3) + t476;
t648 = t592 * t475 - t472;
t737 = pkin(9) * t620;
t428 = t648 - t737;
t426 = pkin(4) * t582 + t428;
t474 = t592 * t477;
t624 = -t475 * t587 - t474;
t736 = pkin(9) * t621;
t429 = -t624 + t736;
t707 = t429 * t591;
t408 = t426 * t586 + t707;
t406 = pkin(10) * t576 + t408;
t566 = -pkin(2) * t584 - pkin(1);
t540 = qJD(1) * t566 + qJD(2);
t503 = pkin(3) * t520 + t540;
t462 = -pkin(4) * t621 + t503;
t417 = pkin(5) * t452 - pkin(10) * t720 + t462;
t386 = -t406 * t585 + t417 * t590;
t708 = t429 * t586;
t407 = t426 * t591 - t708;
t405 = -pkin(5) * t576 - t407;
t758 = -t386 * t720 + t405 * t668 + t590 * t760;
t387 = t406 * t590 + t417 * t585;
t663 = qJD(1) * qJD(2);
t721 = qJDD(1) * t712 + t663;
t511 = t721 * t583;
t512 = t721 * t584;
t644 = -t593 * t511 - t588 * t512;
t435 = qJDD(3) * pkin(3) - pkin(8) * t499 + qJD(3) * t619 + t644;
t622 = -t588 * t511 + t593 * t512;
t439 = -pkin(8) * t500 + qJD(3) * t732 + t622;
t606 = qJD(4) * t624 + t592 * t435 - t587 * t439;
t390 = pkin(4) * t578 - pkin(9) * t440 + t606;
t723 = t592 * (qJD(4) * t475 + t439) + t587 * t435 - t477 * t672;
t391 = -pkin(9) * t441 + t723;
t722 = qJD(5) * t408 - t591 * t390 + t586 * t391;
t381 = -pkin(5) * t573 + t722;
t561 = cos(t567);
t713 = g(3) * t561;
t746 = t381 + t713;
t757 = t387 * t720 + t405 * t667 + t585 * t746;
t749 = t405 * t452;
t755 = -t585 * t760 + t757;
t754 = -t590 * t746 + t758;
t738 = pkin(5) * t720;
t753 = pkin(10) * t452 + t738;
t550 = g(3) * t560;
t724 = (qJD(5) * t426 + t391) * t591 + t586 * t390 - t429 * t670;
t743 = t452 * t462 + t561 * t630 + t550 - t724;
t727 = -t462 * t720 - t713 - t722 + t760;
t559 = t573 * MDP(26);
t697 = t621 * t582;
t698 = t620 * t582;
t742 = t578 * MDP(19) + t559 + (t620 ^ 2 - t621 ^ 2) * MDP(16) - t621 * MDP(15) * t620 + (t440 - t697) * MDP(17) + (-t441 + t698) * MDP(18);
t739 = pkin(4) * t620;
t733 = t744 * (pkin(10) * t744 + t738);
t643 = -t593 * t545 - t546 * t588;
t487 = -pkin(8) * t531 + t643;
t680 = -t588 * t545 + t593 * t546;
t488 = -pkin(8) * t530 + t680;
t683 = t587 * t487 + t592 * t488;
t731 = qJ(2) * qJDD(1);
t629 = g(1) * t589 - g(2) * t594;
t730 = qJDD(2) - t629;
t564 = sin(t577);
t565 = cos(t577);
t726 = -g(3) * t565 - t503 * t620 + t564 * t630 + t606;
t725 = g(3) * t564 - t503 * t621 + t565 * t630 - t723;
t502 = -t530 * t587 + t531 * t592;
t522 = t530 * qJD(3);
t459 = qJD(4) * t502 - t522 * t587 + t592 * t523;
t608 = -t545 * t673 + qJD(2) * t687 + (-qJD(2) * t583 - qJD(3) * t546) * t588;
t466 = -pkin(8) * t523 + t608;
t599 = -t531 * qJD(2) - qJD(3) * t680;
t467 = pkin(8) * t522 + t599;
t613 = t592 * t466 + t587 * t467 + t487 * t671 - t488 * t672;
t413 = -pkin(9) * t459 + t613;
t501 = t592 * t530 + t531 * t587;
t458 = -qJD(4) * t501 - t522 * t592 - t523 * t587;
t605 = -qJD(4) * t683 - t466 * t587 + t592 * t467;
t414 = -pkin(9) * t458 + t605;
t646 = t592 * t487 - t488 * t587;
t432 = -pkin(9) * t502 + t646;
t433 = -pkin(9) * t501 + t683;
t625 = t432 * t591 - t433 * t586;
t382 = qJD(5) * t625 + t413 * t591 + t414 * t586;
t416 = t432 * t586 + t433 * t591;
t460 = t591 * t501 + t502 * t586;
t418 = -qJD(5) * t460 + t458 * t591 - t459 * t586;
t461 = -t501 * t586 + t502 * t591;
t505 = pkin(3) * t530 + t566;
t470 = pkin(4) * t501 + t505;
t422 = pkin(5) * t460 - pkin(10) * t461 + t470;
t380 = pkin(10) * t573 + t724;
t638 = qJD(6) * t417 + t380;
t719 = t381 * t461 - t416 * t400 + t405 * t418 - (qJD(6) * t422 + t382) * t744 - t460 * t638;
t716 = pkin(3) * t523;
t711 = qJDD(1) * pkin(1);
t710 = t405 * t461;
t709 = t422 * t400;
t704 = t446 * t585;
t694 = t584 * MDP(4);
t693 = t585 * t589;
t692 = t585 * t594;
t691 = t586 * t587;
t690 = t587 * t591;
t689 = t589 * t590;
t688 = t590 * t594;
t684 = t592 * t476 - t472;
t647 = -t476 * t587 - t474;
t430 = t647 - t736;
t431 = t684 - t737;
t571 = pkin(3) * t592 + pkin(4);
t682 = t430 * t586 + t431 * t591 - t571 * t669 - (-t587 * t670 + (t591 * t592 - t691) * qJD(4)) * pkin(3);
t681 = t430 * t591 - t431 * t586 + t571 * t670 + (t587 * t669 + (t586 * t592 + t690) * qJD(4)) * pkin(3);
t678 = pkin(3) * t690 + t586 * t571;
t677 = t583 ^ 2 + t584 ^ 2;
t534 = qJDD(1) * t566 + qJDD(2);
t478 = pkin(3) * t500 + t534;
t425 = pkin(4) * t441 + t478;
t385 = pkin(5) * t402 - pkin(10) * t401 + t425;
t636 = qJD(6) * t406 - t385;
t468 = pkin(3) * t521 + t739;
t517 = pkin(10) + t678;
t635 = qJD(6) * t517 + t468 + t753;
t569 = pkin(4) * t586 + pkin(10);
t634 = qJD(6) * t569 + t739 + t753;
t633 = 0.2e1 * t677;
t409 = t428 * t586 + t707;
t632 = pkin(4) * t670 - t409;
t410 = t428 * t591 - t708;
t631 = -pkin(4) * t669 + t410;
t447 = pkin(4) * t459 + t716;
t627 = -t400 * t569 + t749;
t626 = -t400 * t517 + t749;
t618 = t711 - t730;
t617 = t398 - (t452 * t585 + t668) * t744;
t616 = -pkin(3) * t691 + t571 * t591;
t614 = t418 * t590 - t461 * t668;
t612 = -pkin(10) * t400 + t407 * t744 + t749;
t602 = t633 * t663 - t630;
t575 = cos(t581);
t574 = sin(t581);
t570 = -pkin(4) * t591 - pkin(5);
t516 = -pkin(5) - t616;
t510 = t561 * t688 + t693;
t509 = -t561 * t692 + t689;
t508 = -t561 * t689 + t692;
t507 = t561 * t693 + t688;
t419 = qJD(5) * t461 + t458 * t586 + t591 * t459;
t388 = pkin(5) * t419 - pkin(10) * t418 + t447;
t384 = t590 * t385;
t383 = qJD(5) * t416 + t413 * t586 - t414 * t591;
t1 = [(-t382 * t576 + t401 * t470 - t416 * t573 + t418 * t462 + t425 * t461 + t447 * t720 - t560 * t629) * MDP(28) + (-t401 * t460 - t402 * t461 - t418 * t452 - t419 * t720) * MDP(23) + (t401 * t461 + t418 * t720) * MDP(22) + (-t383 * t576 + t402 * t470 + t419 * t462 + t425 * t460 + t447 * t452 + t561 * t629 + t573 * t625) * MDP(27) + t629 * MDP(2) + t630 * MDP(3) + (t458 * t582 + t502 * t578) * MDP(17) + (-t459 * t582 - t501 * t578) * MDP(18) + (t418 * t576 + t461 * t573) * MDP(24) + (-t419 * t576 - t460 * t573) * MDP(25) + (t394 * t461 + t446 * t614) * MDP(29) + ((-t444 * t590 - t704) * t418 + (-t393 - t396 * t590 + (t444 * t585 - t446 * t590) * qJD(6)) * t461) * MDP(30) + (-qJD(3) * t523 - qJDD(3) * t530) * MDP(11) + (-t499 * t530 - t500 * t531 + t520 * t522 - t521 * t523) * MDP(9) + (t499 * t531 - t521 * t522) * MDP(8) + (-qJD(3) * t522 + qJDD(3) * t531) * MDP(10) + (-qJD(3) * t608 - qJDD(3) * t680 + t566 * t499 - t540 * t522 + t534 * t531 - t574 * t629) * MDP(14) + (t505 * t441 + t503 * t459 + t478 * t501 + t565 * t629 + t578 * t646 + t582 * t605 - t621 * t716) * MDP(20) + (t505 * t440 + t503 * t458 + t478 * t502 - t564 * t629 - t578 * t683 - t582 * t613 + t620 * t716) * MDP(21) + (t440 * t502 + t458 * t620) * MDP(15) + (-t440 * t501 - t441 * t502 + t458 * t621 - t459 * t620) * MDP(16) + (-MDP(5) * t583 + t694) * (t618 + t711) + (-g(1) * t508 - g(2) * t510 + t383 * t444 + t384 * t460 + t386 * t419 - t625 * t396 + (t388 * t744 + t709 + (-t406 * t460 - t416 * t744 + t710) * qJD(6)) * t590 + t719 * t585) * MDP(34) + (-g(1) * t507 - g(2) * t509 + t383 * t446 - t387 * t419 - t625 * t395 + (-(-qJD(6) * t416 + t388) * t744 - t709 + t636 * t460 - qJD(6) * t710) * t585 + t719 * t590) * MDP(35) + (t400 * t460 + t419 * t744) * MDP(33) + (-t461 * t397 - t396 * t460 - t419 * t444 + (-t418 * t585 - t461 * t667) * t744) * MDP(32) + (t395 * t460 + t398 * t461 + t419 * t446 + t614 * t744) * MDP(31) + (qJD(3) * t599 + qJDD(3) * t643 + t566 * t500 + t540 * t523 + t534 * t530 + t575 * t629) * MDP(13) + qJDD(1) * MDP(1) + (t618 * pkin(1) + (t677 * t731 + t602) * qJ(2)) * MDP(7) + (t633 * t731 + t602) * MDP(6); t730 * MDP(7) - t556 * MDP(13) + t656 * MDP(14) + (t441 + t698) * MDP(20) + (t440 + t697) * MDP(21) + (t402 + t703) * MDP(27) + (t401 - t700) * MDP(28) + (t617 - t706) * MDP(34) + (-t590 * t744 ^ 2 - t397 - t705) * MDP(35) + (-t694 - pkin(1) * MDP(7) + (MDP(13) * t588 + MDP(5)) * t583) * qJDD(1) + ((qJD(1) * t695 + t584 * t674 + t521) * MDP(13) + (-t520 - t655) * MDP(14)) * qJD(3) + (-MDP(7) * qJ(2) - MDP(6)) * qJD(1) ^ 2 * t677; t742 + t751 + (-t468 * t452 + t573 * t616 - t576 * t681 + t727) * MDP(27) + t628 * MDP(11) + (-t468 * t720 - t573 * t678 + t576 * t682 + t743) * MDP(28) + (-t520 ^ 2 + t521 ^ 2) * MDP(9) + (t617 + t706) * MDP(32) + (g(3) * t574 + t540 * t520 + t575 * t630 - t622) * MDP(14) + (t656 + (t520 - t655) * qJD(3)) * MDP(10) + (t516 * t395 + t626 * t590 + t681 * t446 + (t585 * t635 + t590 * t682) * t744 + t755) * MDP(35) + (t516 * t396 + t626 * t585 + t681 * t444 + (t585 * t682 - t590 * t635) * t744 + t754) * MDP(34) + (-g(3) * t575 - t540 * t521 + t574 * t630 + t644) * MDP(13) + qJDD(3) * MDP(12) + t521 * t520 * MDP(8) + (-t704 * t744 + t752) * MDP(30) + (t684 * t582 + (-t521 * t620 - t578 * t587 - t582 * t671) * pkin(3) + t725) * MDP(21) + (-t647 * t582 + (t521 * t621 + t578 * t592 - t582 * t672) * pkin(3) + t726) * MDP(20); (-t582 * t624 + t726) * MDP(20) + (t582 * t648 + t725) * MDP(21) + (t409 * t576 + (-t452 * t620 + t573 * t591 - t576 * t670) * pkin(4) + t727) * MDP(27) + (t410 * t576 + (-t573 * t586 - t576 * t669 - t620 * t720) * pkin(4) + t743) * MDP(28) + (t570 * t396 + t627 * t585 + t632 * t444 + (t585 * t631 - t590 * t634) * t744 + t754) * MDP(34) + (t570 * t395 + t627 * t590 + t632 * t446 + (t585 * t634 + t590 * t631) * t744 + t755) * MDP(35) + t742 + t761; t559 + (t408 * t576 + t727) * MDP(27) + (t407 * t576 + t743) * MDP(28) + (-pkin(5) * t396 - t408 * t444 + t612 * t585 + (-t746 - t733) * t590 + t758) * MDP(34) + (-pkin(5) * t395 - t408 * t446 + t612 * t590 + (-t760 + t733) * t585 + t757) * MDP(35) + t761; t446 * t444 * MDP(29) + (-t444 ^ 2 + t446 ^ 2) * MDP(30) + (t444 * t744 + t658) * MDP(31) + (t446 * t744 + t552) * MDP(32) + t400 * MDP(33) + (-g(1) * t509 + g(2) * t507 + t387 * t744 - t405 * t446 + t384) * MDP(34) + (g(1) * t510 - g(2) * t508 + t386 * t744 + t405 * t444) * MDP(35) + ((-t380 + t550) * MDP(35) + (-MDP(32) * t720 - MDP(34) * t406 - MDP(35) * t417) * qJD(6)) * t590 + (-qJD(6) * t720 * MDP(31) + (-qJD(6) * t576 - t401) * MDP(32) + (-t638 + t550) * MDP(34) + t636 * MDP(35)) * t585;];
tau  = t1;
