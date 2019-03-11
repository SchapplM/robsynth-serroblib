% Calculate vector of inverse dynamics joint torques for
% S6RRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:14:09
% EndTime: 2019-03-09 10:14:19
% DurationCPUTime: 7.83s
% Computational Cost: add. (6869->516), mult. (16240->620), div. (0->0), fcn. (12510->14), ass. (0->236)
t550 = qJD(2) + qJD(4);
t557 = sin(qJ(4));
t558 = sin(qJ(2));
t561 = cos(qJ(2));
t675 = sin(pkin(10));
t676 = cos(pkin(10));
t589 = -t558 * t676 - t561 * t675;
t502 = t589 * qJD(1);
t590 = t675 * t558 - t676 * t561;
t568 = qJD(2) * t502 - qJDD(1) * t590;
t634 = qJD(1) * qJD(2);
t571 = qJDD(1) * t589 + t590 * t634;
t686 = cos(qJ(4));
t569 = -t557 * t571 - t686 * t568;
t501 = t590 * qJD(1);
t600 = -t557 * t501 - t502 * t686;
t404 = qJD(4) * t600 + t569;
t456 = -t686 * t501 + t502 * t557;
t548 = qJDD(2) + qJDD(4);
t556 = sin(qJ(6));
t560 = cos(qJ(6));
t639 = qJD(6) * t560;
t631 = t556 * t404 - t456 * t639 + t560 * t548;
t640 = qJD(6) * t556;
t380 = -t550 * t640 + t631;
t379 = t380 * t560;
t399 = t560 * t404;
t440 = -t456 * t556 + t550 * t560;
t381 = qJD(6) * t440 + t548 * t556 - t399;
t627 = qJD(4) * t686;
t643 = qJD(4) * t557;
t567 = -t501 * t643 - t502 * t627 + t569;
t598 = -t501 * t627 + t502 * t643 + t557 * t568 - t571 * t686;
t664 = t456 * t550;
t584 = t598 - t664;
t666 = t600 * t550;
t694 = qJD(6) + t600;
t709 = t694 * t560;
t710 = t694 * t556;
t438 = t456 * t560 + t550 * t556;
t711 = t438 * t694;
t714 = t584 * MDP(15) + t548 * MDP(17) + (-t567 + t666) * MDP(16) + (-t440 * t710 + t379) * MDP(24) + (-t440 * t709 + (-t380 + t711) * t556 - t560 * t381) * MDP(25);
t402 = -qJDD(6) - t598;
t396 = t560 * t402;
t708 = -t694 * t710 - t396;
t671 = t402 * t556;
t591 = -t694 * t709 + t671;
t684 = t456 * pkin(5);
t674 = qJ(5) * t456;
t678 = qJ(3) + pkin(7);
t526 = t678 * t561;
t515 = qJD(1) * t526;
t504 = t675 * t515;
t525 = t678 * t558;
t514 = qJD(1) * t525;
t677 = qJD(2) * pkin(2);
t508 = -t514 + t677;
t464 = t676 * t508 - t504;
t680 = t502 * pkin(8);
t435 = qJD(2) * pkin(3) + t464 + t680;
t622 = t676 * t515;
t465 = t675 * t508 + t622;
t681 = t501 * pkin(8);
t437 = t465 - t681;
t409 = t557 * t435 + t686 * t437;
t397 = -qJ(5) * t550 - t409;
t384 = -t397 + t684;
t712 = t384 * t694;
t469 = t514 * t675 - t622;
t441 = t469 + t681;
t470 = -t676 * t514 - t504;
t442 = t470 + t680;
t629 = pkin(2) * t675;
t532 = t557 * t629;
t628 = t676 * pkin(2);
t538 = t628 + pkin(3);
t707 = qJD(4) * t532 + t557 * t441 + t442 * t686 - t538 * t627;
t546 = t561 * pkin(2);
t552 = qJ(2) + pkin(10);
t545 = qJ(4) + t552;
t536 = sin(t545);
t537 = cos(t545);
t648 = t537 * pkin(4) + t536 * qJ(5);
t706 = pkin(3) * cos(t552) + t546 + t648;
t408 = -t686 * t435 + t557 * t437;
t637 = qJD(5) + t408;
t703 = t600 ^ 2;
t702 = pkin(4) * t600;
t701 = pkin(5) * t600;
t540 = t546 + pkin(1);
t700 = t384 * t600;
t687 = pkin(4) + pkin(9);
t698 = t600 * t687;
t587 = t557 * t538 + t629 * t686;
t649 = t587 * qJD(4) + t441 * t686 - t557 * t442;
t697 = t649 * t550;
t650 = -qJD(5) + t707;
t541 = t548 * qJ(5);
t696 = -t550 * qJD(5) - t541;
t559 = sin(qJ(1));
t562 = cos(qJ(1));
t695 = g(1) * t559 - g(2) * t562;
t638 = t701 + t637;
t518 = -qJD(1) * t540 + qJD(3);
t471 = pkin(3) * t501 + t518;
t575 = -qJ(5) * t600 + t471;
t411 = -pkin(4) * t456 + t575;
t531 = g(3) * t537;
t623 = qJD(2) * t678;
t499 = -qJD(3) * t558 - t561 * t623;
t463 = qJDD(2) * pkin(2) + qJD(1) * t499 - qJDD(1) * t525;
t498 = qJD(3) * t561 - t558 * t623;
t468 = qJD(1) * t498 + qJDD(1) * t526;
t419 = t676 * t463 - t468 * t675;
t407 = qJDD(2) * pkin(3) + pkin(8) * t571 + t419;
t420 = t675 * t463 + t676 * t468;
t410 = pkin(8) * t568 + t420;
t614 = -t686 * t407 + t557 * t410 + t435 * t643 + t437 * t627;
t660 = t536 * t562;
t661 = t536 * t559;
t595 = -g(1) * t660 - g(2) * t661 + t531 + t614;
t576 = t411 * t600 + qJDD(5) + t595;
t609 = t538 * t686 - t532;
t496 = -pkin(4) - t609;
t488 = -pkin(9) + t496;
t693 = -t694 * (-t649 + t684) - t488 * t402;
t692 = -t471 * t600 - t595;
t446 = -t498 * t675 + t676 * t499;
t579 = qJD(2) * t590;
t427 = pkin(8) * t579 + t446;
t447 = t676 * t498 + t675 * t499;
t578 = qJD(2) * t589;
t428 = pkin(8) * t578 + t447;
t472 = -t676 * t525 - t526 * t675;
t448 = pkin(8) * t589 + t472;
t473 = -t675 * t525 + t676 * t526;
t449 = -pkin(8) * t590 + t473;
t375 = -t557 * t427 - t686 * t428 - t448 * t627 + t449 * t643;
t417 = t557 * t448 + t449 * t686;
t691 = -t375 * t550 + t417 * t548 + t536 * t695;
t376 = qJD(4) * t417 - t427 * t686 + t557 * t428;
t416 = -t448 * t686 + t557 * t449;
t690 = t376 * t550 + t416 * t548 - t537 * t695;
t382 = -t550 * t687 + t638;
t387 = -t456 * t687 + t575;
t370 = t382 * t560 - t387 * t556;
t371 = t382 * t556 + t387 * t560;
t689 = MDP(13) * t600 + t471 * MDP(19) - t411 * MDP(22) + MDP(28) * t694 + t370 * MDP(29) - t371 * MDP(30);
t685 = pkin(4) * t548;
t530 = g(3) * t536;
t682 = g(3) * t561;
t679 = t558 * pkin(2);
t673 = qJDD(1) * pkin(1);
t577 = t686 * t590;
t466 = -t557 * t589 + t577;
t580 = t557 * t590;
t467 = -t589 * t686 - t580;
t478 = pkin(3) * t590 - t540;
t574 = -t467 * qJ(5) + t478;
t391 = t466 * t687 + t574;
t672 = t391 * t402;
t670 = t409 * t550;
t669 = t438 * t456;
t668 = t440 * t456;
t665 = t456 ^ 2;
t663 = t466 * t556;
t659 = t537 * t559;
t658 = t537 * t562;
t657 = t556 * t559;
t656 = t556 * t562;
t655 = t559 * t560;
t654 = t560 * t562;
t613 = -t557 * t407 - t686 * t410 - t435 * t627 + t437 * t643;
t367 = t613 + t696;
t365 = -pkin(5) * t404 - t367;
t653 = t365 * t556 + t384 * t639;
t651 = -t650 + t701;
t553 = t558 ^ 2;
t646 = -t561 ^ 2 + t553;
t642 = qJD(6) * t387;
t641 = qJD(6) * t550;
t633 = qJDD(1) * t561;
t626 = t558 * t634;
t632 = pkin(2) * t626 + qJDD(3);
t543 = t558 * t677;
t474 = -pkin(3) * t502 + qJD(1) * t679;
t625 = -pkin(4) * t536 - pkin(3) * sin(t552) - t679;
t436 = -pkin(2) * t633 - pkin(3) * t568 + t632 - t673;
t369 = t404 * pkin(4) - qJ(5) * t598 - qJD(5) * t600 + t436;
t366 = t404 * pkin(9) + t369;
t615 = qJD(6) * t382 + t366;
t611 = g(1) * t562 + g(2) * t559;
t608 = -t642 + t531;
t605 = (t409 + t684) * t694 - t687 * t402;
t604 = qJDD(5) + t614;
t603 = t474 - t674;
t602 = pkin(1) + t706;
t601 = -0.2e1 * pkin(1) * t634 - pkin(7) * qJDD(2);
t422 = -qJD(4) * t580 - t557 * t579 - t578 * t686 - t589 * t627;
t599 = t422 * t556 + t466 * t639;
t596 = -g(1) * t658 - g(2) * t659 - t530 - t613;
t593 = -qJDD(1) * t540 + t632;
t392 = t467 * pkin(5) + t416;
t592 = t365 * t466 + t384 * t422 + t392 * t402;
t586 = -t537 * t611 - t530;
t585 = t596 - t696;
t564 = qJD(2) ^ 2;
t582 = -pkin(7) * t564 + 0.2e1 * t673 + t695;
t565 = qJD(1) ^ 2;
t581 = pkin(1) * t565 - pkin(7) * qJDD(1) + t611;
t475 = -pkin(3) * t578 + t543;
t573 = (-qJD(6) * t488 + t603 + t698) * t694 + t586;
t572 = (qJD(6) * t687 - t674 + t698) * t694 + t586;
t421 = t550 * t577 - t557 * t578 - t589 * t643;
t377 = t422 * pkin(4) + t421 * qJ(5) - t467 * qJD(5) + t475;
t549 = -pkin(8) - t678;
t520 = qJ(5) * t658;
t519 = qJ(5) * t659;
t495 = qJ(5) + t587;
t494 = -t536 * t657 + t654;
t493 = t536 * t655 + t656;
t492 = t536 * t656 + t655;
t491 = t536 * t654 - t657;
t418 = -t674 + t702;
t415 = t466 * pkin(4) + t574;
t414 = t603 + t702;
t395 = -pkin(4) * t550 + t637;
t393 = -t466 * pkin(5) + t417;
t374 = t422 * pkin(9) + t377;
t373 = -t421 * pkin(5) + t376;
t372 = -pkin(5) * t422 - t375;
t368 = t604 - t685;
t364 = pkin(5) * t598 - t548 * t687 + t604;
t363 = t365 * t560;
t361 = t560 * t364;
t1 = [qJDD(1) * MDP(1) + ((-t438 * t556 + t440 * t560) * t422 + (t379 - t381 * t556 + (-t438 * t560 - t440 * t556) * qJD(6)) * t466) * MDP(25) + (t380 * t663 + t440 * t599) * MDP(24) + 0.2e1 * (t558 * t633 - t634 * t646) * MDP(5) + (-t421 * t471 + t436 * t467 + t475 * t600 + t478 * t598 - t691) * MDP(19) + (-t369 * t467 - t377 * t600 + t411 * t421 - t415 * t598 + t691) * MDP(22) + (-t421 * t600 + t467 * t598) * MDP(13) + (-g(1) * t494 - g(2) * t492 + t361 * t467 - t370 * t421 + t372 * t438 + t393 * t381 + (-t366 * t467 - t374 * t694 + t672) * t556 + (t373 * t694 - t592) * t560 + ((-t391 * t560 - t392 * t556) * t694 - t371 * t467 + t384 * t663) * qJD(6)) * MDP(29) + (g(1) * t493 - g(2) * t491 + t371 * t421 + t372 * t440 + t393 * t380 + (-(qJD(6) * t392 + t374) * t694 + t672 - t615 * t467 + t384 * qJD(6) * t466) * t560 + (-(-qJD(6) * t391 + t373) * t694 - (t364 - t642) * t467 + t592) * t556) * MDP(30) + (t380 * t467 - t402 * t663 - t421 * t440 + t599 * t694) * MDP(26) + (-t466 * t396 - t381 * t467 + t421 * t438 + (t422 * t560 - t466 * t640) * t694) * MDP(27) + (-t402 * t467 - t421 * t694) * MDP(28) + (-t369 * t466 + t377 * t456 - t404 * t415 - t411 * t422 + t690) * MDP(21) + (t404 * t478 + t422 * t471 + t436 * t466 - t456 * t475 - t690) * MDP(18) + (t367 * t466 + t368 * t467 - t375 * t456 + t376 * t600 - t395 * t421 + t397 * t422 - t404 * t417 + t416 * t598 - t611) * MDP(20) + (-t404 * t467 - t421 * t456 - t422 * t600 - t466 * t598) * MDP(14) + (t419 * t589 + t446 * t502 - t447 * t501 + t465 * t578 + t472 * t571 + t473 * t568 - t611 + (qJD(2) * t464 - t420) * t590) * MDP(11) + (t420 * t473 + t465 * t447 + t419 * t472 + t464 * t446 - t593 * t540 + t518 * t543 - g(1) * (-t540 * t559 + t562 * t678) - g(2) * (t540 * t562 + t559 * t678)) * MDP(12) + t695 * MDP(2) + (qJDD(1) * t553 + 0.2e1 * t561 * t626) * MDP(4) + t611 * MDP(3) + (-t367 * t417 + t368 * t416 + t369 * t415 + t397 * t375 + t395 * t376 + t411 * t377 + (g(1) * t549 - g(2) * t602) * t562 + (g(1) * t602 + g(2) * t549) * t559) * MDP(23) + (t558 * t601 + t561 * t582) * MDP(9) + (-t558 * t582 + t561 * t601) * MDP(10) + (-t422 * t550 - t466 * t548) * MDP(16) + (-t421 * t550 + t467 * t548) * MDP(15) + (qJDD(2) * t558 + t561 * t564) * MDP(6) + (qJDD(2) * t561 - t558 * t564) * MDP(7); (-MDP(20) * t395 - t689) * t456 + (-t558 * t561 * MDP(4) + MDP(5) * t646) * t565 + qJDD(2) * MDP(8) + (t495 * t380 + t363 + t651 * t440 + t573 * t560 + (-t693 - t712) * t556) * MDP(30) + (t591 + t669) * MDP(27) + (-t367 * t495 + t368 * t496 - t411 * t414 - g(1) * (t562 * t625 + t520) - g(2) * (t559 * t625 + t519) - g(3) * t706 + t650 * t397 + t649 * t395) * MDP(23) + t558 * qJDD(1) * MDP(6) + (-t474 * t600 - t587 * t548 + t550 * t707 - t596) * MDP(19) + (-t668 + t708) * MDP(26) + (t456 * t474 + t548 * t609 + t692 - t697) * MDP(18) + (-t414 * t456 + t697 + (-pkin(4) + t496) * t548 + t576) * MDP(21) + (t598 * t496 - t404 * t495 - t456 * t650 + (-t397 + t649) * t600) * MDP(20) + (-t665 + t703) * MDP(14) + (t414 * t600 + t495 * t548 - t550 * t650 + t585) * MDP(22) + (t568 * t629 + t571 * t628 + (-t465 - t469) * t502 - (t464 - t470) * t501) * MDP(11) + (t495 * t381 + t651 * t438 + (t693 + t700) * t560 + t573 * t556 + t653) * MDP(29) + (g(3) * t558 + t561 * t581) * MDP(10) + (t558 * t581 - t682) * MDP(9) + (-t464 * t469 - t465 * t470 + (t675 * t420 + t676 * t419 - t682 + (-qJD(1) * t518 + t611) * t558) * pkin(2)) * MDP(12) + t714 + MDP(7) * t633; (-t501 ^ 2 - t502 ^ 2) * MDP(11) + (-t464 * t502 + t465 * t501 + t593 - t695) * MDP(12) + (-t665 - t703) * MDP(20) + (-t395 * t600 + t397 * t456 + t369 - t695) * MDP(23) + (t591 - t669) * MDP(29) + (-t668 - t708) * MDP(30) + (MDP(19) - MDP(22)) * (t598 + t664) + (-MDP(18) + MDP(21)) * (-t567 - t666); t703 * MDP(14) + (t670 + t692) * MDP(18) + (-t408 * t550 - t596) * MDP(19) + (-pkin(4) * t598 - qJ(5) * t404 + (-t397 - t409) * t600) * MDP(20) + (t576 - t670 - 0.2e1 * t685) * MDP(21) + (t418 * t600 + t550 * t637 + t541 + t585) * MDP(22) + (-t367 * qJ(5) - t368 * pkin(4) - t411 * t418 - t395 * t409 - g(1) * (-pkin(4) * t660 + t520) - g(2) * (-pkin(4) * t661 + t519) - g(3) * t648 - t637 * t397) * MDP(23) + t708 * MDP(26) + t591 * MDP(27) + (qJ(5) * t381 + t638 * t438 + (-t605 + t700) * t560 + t572 * t556 + t653) * MDP(29) + (qJ(5) * t380 + t363 + t638 * t440 + (t605 - t712) * t556 + t572 * t560) * MDP(30) - ((t395 - t637) * MDP(20) + t418 * MDP(21) + t440 * MDP(26) - t438 * MDP(27) + MDP(14) * t456 + t689) * t456 + t714; t584 * MDP(20) + (t456 * t600 + t548) * MDP(21) + (-t550 ^ 2 - t703) * MDP(22) + (t397 * t550 + t576 - t685) * MDP(23) + (-t438 * t550 - t396) * MDP(29) + (-t440 * t550 + t671) * MDP(30) + (-MDP(29) * t710 - MDP(30) * t709) * t694; t440 * t438 * MDP(24) + (-t438 ^ 2 + t440 ^ 2) * MDP(25) + (t631 + t711) * MDP(26) + (t440 * t694 + t399) * MDP(27) - t402 * MDP(28) + (-g(1) * t491 - g(2) * t493 + t371 * t694 - t384 * t440 + t361) * MDP(29) + (g(1) * t492 - g(2) * t494 + t370 * t694 + t384 * t438) * MDP(30) + (-MDP(27) * t641 + MDP(29) * t608 - MDP(30) * t615) * t560 + (-MDP(26) * t641 + (qJD(6) * t456 - t548) * MDP(27) - t615 * MDP(29) + (-t364 - t608) * MDP(30)) * t556;];
tau  = t1;
