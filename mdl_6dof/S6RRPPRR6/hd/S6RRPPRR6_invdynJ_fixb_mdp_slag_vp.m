% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:15:33
% EndTime: 2019-03-09 09:15:45
% DurationCPUTime: 9.91s
% Computational Cost: add. (5263->552), mult. (11626->698), div. (0->0), fcn. (8647->12), ass. (0->245)
t584 = sin(pkin(10));
t585 = cos(pkin(10));
t593 = cos(qJ(2));
t674 = qJD(1) * t593;
t589 = sin(qJ(2));
t675 = qJD(1) * t589;
t491 = t584 * t675 + t585 * t674;
t494 = -t584 * t674 + t585 * t675;
t588 = sin(qJ(5));
t592 = cos(qJ(5));
t437 = t491 * t588 - t494 * t592;
t576 = qJDD(2) - qJDD(5);
t577 = qJD(2) - qJD(5);
t587 = sin(qJ(6));
t591 = cos(qJ(6));
t510 = t584 * t589 + t585 * t593;
t660 = qJD(1) * qJD(2);
t650 = t589 * t660;
t529 = t585 * t650;
t649 = t593 * t660;
t449 = qJDD(1) * t510 + t584 * t649 - t529;
t500 = t510 * qJD(2);
t658 = qJDD(1) * t593;
t659 = qJDD(1) * t589;
t626 = -t584 * t658 + t585 * t659;
t450 = qJD(1) * t500 + t626;
t669 = qJD(5) * t592;
t670 = qJD(5) * t588;
t611 = t588 * t449 - t592 * t450 + t491 * t669 + t494 * t670;
t665 = qJD(6) * t591;
t655 = -t587 * t576 - t577 * t665 - t591 * t611;
t666 = qJD(6) * t587;
t379 = t437 * t666 + t655;
t418 = -t437 * t591 - t577 * t587;
t543 = t591 * t576;
t380 = qJD(6) * t418 - t587 * t611 + t543;
t395 = -qJD(5) * t437 + t592 * t449 + t450 * t588;
t393 = qJDD(6) + t395;
t390 = t591 * t393;
t712 = t592 * t491 + t494 * t588;
t722 = qJD(6) + t712;
t693 = t418 * t722;
t416 = -t437 * t587 + t591 * t577;
t694 = t416 * t722;
t695 = t393 * t587;
t697 = t379 * t587;
t724 = t577 * t712;
t725 = t437 * t577;
t726 = t418 * t437;
t727 = t416 * t437;
t733 = t722 * t591;
t734 = t722 ^ 2;
t742 = -((t380 + t693) * t587 - (t379 - t694) * t591) * MDP(27) + (t418 * t733 + t697) * MDP(26) + (t722 * t733 + t695 + t726) * MDP(28) - (t587 * t734 - t390 + t727) * MDP(29) + (t437 ^ 2 - t712 ^ 2) * MDP(20) - (t611 + t724) * MDP(21) + (-t395 + t725) * MDP(22) - t576 * MDP(23) + (-MDP(19) * t712 + t722 * MDP(30)) * t437;
t730 = pkin(5) * t437;
t560 = pkin(7) * t675;
t517 = qJ(4) * t675 - t560;
t705 = pkin(2) + pkin(3);
t651 = t705 * qJD(2);
t475 = qJD(3) - t651 - t517;
t561 = pkin(7) * t674;
t519 = -qJ(4) * t674 + t561;
t580 = qJD(2) * qJ(3);
t505 = t519 + t580;
t426 = t585 * t475 - t505 * t584;
t703 = pkin(8) * t494;
t411 = -qJD(2) * pkin(4) + t426 - t703;
t427 = t584 * t475 + t585 * t505;
t704 = pkin(8) * t491;
t412 = t427 - t704;
t384 = t411 * t588 + t412 * t592;
t382 = -pkin(9) * t577 + t384;
t506 = -qJD(1) * pkin(1) - pkin(2) * t674 - qJ(3) * t675;
t469 = pkin(3) * t674 + qJD(4) - t506;
t436 = pkin(4) * t491 + t469;
t391 = pkin(5) * t712 + pkin(9) * t437 + t436;
t374 = -t382 * t587 + t391 * t591;
t729 = t374 * t437;
t375 = t382 * t591 + t391 * t587;
t728 = t375 * t437;
t723 = t649 + t659;
t590 = sin(qJ(1));
t657 = pkin(10) + qJ(5);
t564 = sin(t657);
t645 = cos(t657);
t625 = t589 * t645;
t710 = -t593 * t564 + t625;
t465 = t710 * t590;
t594 = cos(qJ(1));
t686 = t593 * t594;
t467 = t564 * t686 - t594 * t625;
t606 = t589 * t564 + t593 * t645;
t610 = g(1) * t467 - g(2) * t465 + g(3) * t606;
t556 = pkin(7) * t659;
t648 = pkin(7) * t649 + qJDD(3) + t556;
t671 = qJD(4) * t589;
t432 = -qJ(4) * t723 - qJD(1) * t671 - t705 * qJDD(2) + t648;
t673 = qJD(2) * t589;
t700 = pkin(7) - qJ(4);
t487 = -qJD(4) * t593 - t673 * t700;
t557 = pkin(7) * t658;
t578 = qJDD(2) * qJ(3);
t579 = qJD(2) * qJD(3);
t654 = t557 + t578 + t579;
t435 = -qJ(4) * t658 + qJD(1) * t487 + t654;
t643 = -t585 * t432 + t435 * t584;
t388 = -qJDD(2) * pkin(4) - pkin(8) * t450 - t643;
t403 = t584 * t432 + t585 * t435;
t389 = -pkin(8) * t449 + t403;
t615 = -t592 * t388 + t588 * t389 + t411 * t670 + t412 * t669;
t717 = t436 * t437 + t610 - t615;
t468 = t606 * t594;
t612 = t588 * t388 + t592 * t389 + t411 * t669 - t412 * t670;
t466 = t606 * t590;
t628 = g(2) * t466 + g(3) * t710;
t716 = g(1) * t468 + t436 * t712 - t612 + t628;
t679 = t593 * pkin(2) + t589 * qJ(3);
t713 = -pkin(1) - t679;
t520 = -qJ(3) * t584 - t585 * t705;
t516 = -pkin(4) + t520;
t521 = qJ(3) * t585 - t584 * t705;
t681 = t588 * t516 + t592 * t521;
t575 = g(1) * t594;
t711 = g(2) * t590 + t575;
t554 = qJ(3) * t674;
t484 = -t675 * t705 + t554;
t448 = -pkin(4) * t494 + t484;
t454 = -pkin(9) + t681;
t371 = pkin(5) * t576 + t615;
t607 = -t371 + t610;
t709 = t722 * (-pkin(9) * t712 + qJD(6) * t454 + t448 + t730) + t607;
t708 = t722 * (pkin(9) * t722 - t730) - t607;
t699 = pkin(7) * qJDD(2);
t707 = qJD(2) * (qJD(1) * t713 + t506) - t699;
t528 = t700 * t593;
t488 = qJD(2) * t528 - t671;
t428 = -t487 * t584 + t585 * t488;
t413 = -pkin(8) * t500 + t428;
t429 = t585 * t487 + t584 * t488;
t672 = qJD(2) * t593;
t499 = t584 * t672 - t585 * t673;
t414 = -pkin(8) * t499 + t429;
t527 = t700 * t589;
t460 = t585 * t527 - t528 * t584;
t618 = t584 * t593 - t585 * t589;
t430 = pkin(8) * t618 + t460;
t461 = t584 * t527 + t585 * t528;
t431 = -pkin(8) * t510 + t461;
t621 = t430 * t592 - t431 * t588;
t376 = qJD(5) * t621 + t413 * t588 + t414 * t592;
t383 = t411 * t592 - t412 * t588;
t381 = pkin(5) * t577 - t383;
t451 = t592 * t510 - t588 * t618;
t452 = -t510 * t588 - t592 * t618;
t568 = t593 * pkin(3);
t653 = t568 + t679;
t508 = pkin(1) + t653;
t459 = pkin(4) * t510 + t508;
t399 = pkin(5) * t451 - pkin(9) * t452 + t459;
t401 = t430 * t588 + t431 * t592;
t406 = -qJD(5) * t451 - t499 * t588 + t500 * t592;
t636 = -pkin(9) * t576 + qJD(6) * t391 + t612;
t706 = t371 * t452 + t381 * t406 - t401 * t393 - (qJD(6) * t399 + t376) * t722 - t451 * t636 + t575;
t702 = g(1) * t466;
t574 = g(1) * t590;
t701 = g(2) * t594;
t583 = qJDD(1) * pkin(1);
t698 = qJDD(2) * pkin(2);
t696 = t381 * t452;
t691 = t589 * t590;
t690 = t589 * t594;
t597 = qJD(1) ^ 2;
t689 = t589 * t597;
t688 = t590 * t593;
t455 = -t517 * t584 + t585 * t519;
t419 = t455 - t704;
t456 = t585 * t517 + t584 * t519;
t420 = t456 + t703;
t509 = t584 * t588 - t585 * t592;
t620 = t516 * t592 - t521 * t588;
t685 = -qJD(3) * t509 + qJD(5) * t620 - t419 * t588 - t420 * t592;
t512 = t584 * t592 + t585 * t588;
t684 = t512 * qJD(3) + qJD(5) * t681 + t419 * t592 - t420 * t588;
t683 = t577 * t509;
t682 = t577 * t512;
t565 = t589 * qJD(3);
t680 = qJ(3) * t672 + t565;
t581 = t589 ^ 2;
t582 = t593 ^ 2;
t677 = t581 - t582;
t668 = qJD(6) * t382;
t667 = qJD(6) * t437;
t656 = t593 * t689;
t652 = -g(1) * t690 - g(2) * t691 + g(3) * t593;
t647 = t574 - t701;
t644 = -qJD(2) * pkin(2) + qJD(3);
t638 = qJD(3) * t584 + t455;
t637 = qJD(3) * t585 - t456;
t633 = t594 * pkin(1) + pkin(2) * t686 + t590 * pkin(7) + qJ(3) * t690;
t632 = -t556 - t652;
t631 = t589 * t651;
t596 = qJD(2) ^ 2;
t630 = pkin(7) * t596 + t701;
t627 = t399 * t393 + t702;
t624 = pkin(2) * t658 + qJ(3) * t723 + qJD(1) * t565 + t583;
t623 = -t668 - t701;
t622 = t426 * t584 - t427 * t585;
t524 = t560 + t644;
t526 = t561 + t580;
t619 = t524 * t593 - t526 * t589;
t617 = qJD(6) * t512 + t675;
t483 = t648 - t698;
t614 = -0.2e1 * pkin(1) * t660 - t699;
t613 = t406 * t591 - t452 * t666;
t464 = -t631 + t680;
t609 = -t630 + 0.2e1 * t583;
t608 = pkin(3) * t658 + qJDD(4) + t624;
t605 = t628 - t636;
t434 = pkin(4) * t499 + t464;
t604 = -pkin(9) * t393 + (t381 + t383) * t722;
t602 = -t454 * t393 + (-t381 - t685) * t722;
t447 = pkin(2) * t650 - t624;
t489 = pkin(2) * t673 - t680;
t601 = -qJD(1) * t489 - qJDD(1) * t713 - t447 - t630;
t474 = -pkin(7) * t650 + t654;
t599 = qJD(2) * t619 + t474 * t593 + t483 * t589;
t415 = -qJD(1) * t631 + t608;
t404 = pkin(4) * t449 + t415;
t570 = t594 * pkin(7);
t547 = g(1) * t688;
t541 = qJ(3) * t686;
t539 = qJ(3) * t688;
t518 = pkin(2) * t675 - t554;
t482 = t510 * t594;
t481 = t618 * t594;
t480 = t510 * t590;
t479 = t618 * t590;
t458 = t468 * t591 - t587 * t590;
t457 = -t468 * t587 - t590 * t591;
t453 = pkin(5) - t620;
t407 = qJD(5) * t452 + t592 * t499 + t500 * t588;
t378 = pkin(5) * t407 - pkin(9) * t406 + t434;
t377 = qJD(5) * t401 - t413 * t592 + t414 * t588;
t373 = pkin(5) * t395 + pkin(9) * t611 + t404;
t372 = t591 * t373;
t1 = [t647 * MDP(2) + (t407 * t577 + t451 * t576) * MDP(22) + (-t406 * t577 - t452 * t576) * MDP(21) + (t379 * t452 * t591 + t418 * t613) * MDP(26) + (t589 * t614 + t593 * t609 + t547) * MDP(9) + (g(1) * t480 - g(2) * t482 - qJD(2) * t428 - qJDD(2) * t460 + t415 * t510 + t449 * t508 + t464 * t491 + t469 * t499) * MDP(15) + (-t707 * t593 + (t601 + t574) * t589) * MDP(13) + (qJDD(1) * t581 + 0.2e1 * t589 * t649) * MDP(4) + (t589 * t707 + t601 * t593 + t547) * MDP(11) + (-t406 * t437 - t452 * t611) * MDP(19) + (-t395 * t452 - t406 * t712 + t407 * t437 + t451 * t611) * MDP(20) + (-g(1) * t479 + g(2) * t481 + qJD(2) * t429 + qJDD(2) * t461 - t415 * t618 + t450 * t508 + t464 * t494 + t469 * t500) * MDP(16) + (t403 * t461 + t427 * t429 - t643 * t460 + t426 * t428 + t415 * t508 + t469 * t464 - g(1) * (-qJ(4) * t594 + t570) - g(2) * (pkin(3) * t686 + t633) + (-g(1) * (t713 - t568) + g(2) * qJ(4)) * t590) * MDP(18) + (-t403 * t510 - t426 * t500 - t427 * t499 - t428 * t494 - t429 * t491 - t449 * t461 - t450 * t460 - t618 * t643 + t711) * MDP(17) + ((t581 + t582) * qJDD(1) * pkin(7) + t599 - t711) * MDP(12) + t711 * MDP(3) + (g(1) * t465 + g(2) * t467 + t376 * t577 + t401 * t576 + t404 * t452 + t406 * t436 - t434 * t437 - t459 * t611) * MDP(25) + (t393 * t451 + t407 * t722) * MDP(30) + (-g(2) * t457 - t375 * t407 + t377 * t418 - t621 * t379 + (-(-qJD(6) * t401 + t378) * t722 - (t373 - t668) * t451 - qJD(6) * t696 - t627) * t587 + t706 * t591) * MDP(32) + (-g(2) * t458 + t372 * t451 + t374 * t407 + t377 * t416 - t621 * t380 + (t378 * t722 + (-t382 * t451 - t401 * t722 + t696) * qJD(6) + t627) * t591 + t706 * t587) * MDP(31) + (t379 * t451 + t390 * t452 + t407 * t418 + t613 * t722) * MDP(28) + (-t452 * t695 - t380 * t451 - t407 * t416 + (-t406 * t587 - t452 * t665) * t722) * MDP(29) + (-g(2) * t468 + t377 * t577 + t395 * t459 + t404 * t451 + t407 * t436 + t434 * t712 - t576 * t621 + t702) * MDP(24) + (qJDD(2) * t589 + t593 * t596) * MDP(6) + (qJDD(2) * t593 - t589 * t596) * MDP(7) + qJDD(1) * MDP(1) + (pkin(7) * t599 - g(1) * t570 - g(2) * t633 + t506 * t489 + (t447 - t574) * t713) * MDP(14) + 0.2e1 * (t589 * t658 - t660 * t677) * MDP(5) + ((-t416 * t591 - t418 * t587) * t406 + (-t697 - t380 * t591 + (t416 * t587 - t418 * t591) * qJD(6)) * t452) * MDP(27) + (t614 * t593 + (-t609 - t574) * t589) * MDP(10); (-g(1) * t481 - g(2) * t479 - g(3) * t510 + qJD(2) * t638 - qJDD(2) * t520 + t469 * t494 - t484 * t491 + t643) * MDP(15) + MDP(7) * t658 + MDP(6) * t659 + (-t449 * t521 - t450 * t520 + (-t427 + t638) * t494 + (t426 - t637) * t491) * MDP(17) + (-g(1) * t482 - g(2) * t480 + g(3) * t618 + qJD(2) * t637 + qJDD(2) * t521 - t469 * t491 - t484 * t494 + t403) * MDP(16) + (t557 + 0.2e1 * t578 + 0.2e1 * t579 + (qJD(1) * t518 - g(3)) * t589 + (qJD(1) * t506 - t711) * t593) * MDP(13) + (g(3) * t589 - t557 + (pkin(1) * t597 + t711) * t593) * MDP(10) + (t453 * t379 + t684 * t418 + t587 * t709 + t602 * t591 + t728) * MDP(32) + (t453 * t380 + t684 * t416 + t602 * t587 - t591 * t709 - t729) * MDP(31) + (t437 * t448 + t576 * t681 + t577 * t685 - t716) * MDP(25) + (-t448 * t712 - t576 * t620 + t577 * t684 - t717) * MDP(24) + (t589 * t705 * t711 - g(1) * t541 - g(2) * t539 - g(3) * t653 - t622 * qJD(3) + t403 * t521 - t426 * t455 - t427 * t456 - t469 * t484 - t520 * t643) * MDP(18) + qJDD(2) * MDP(8) + t677 * MDP(5) * t597 - MDP(4) * t656 + ((-pkin(2) * t589 + qJ(3) * t593) * qJDD(1) + ((t526 - t580) * t589 + (-t524 + t644) * t593) * qJD(1)) * MDP(12) + (pkin(1) * t689 + t632) * MDP(9) + (t474 * qJ(3) + t526 * qJD(3) - t483 * pkin(2) - t506 * t518 - g(1) * (-pkin(2) * t690 + t541) - g(2) * (-pkin(2) * t691 + t539) - g(3) * t679 - t619 * qJD(1) * pkin(7)) * MDP(14) + (0.2e1 * t698 - qJDD(3) + (-t506 * t589 + t518 * t593) * qJD(1) + t632) * MDP(11) - t742; (-qJDD(2) - t656) * MDP(11) + MDP(12) * t659 + (-t581 * t597 - t596) * MDP(13) + (-qJD(2) * t526 + t506 * t675 + t483 + t652) * MDP(14) + (-qJDD(2) * t585 - t491 * t675 - t584 * t596) * MDP(15) + (qJDD(2) * t584 - t494 * t675 - t585 * t596) * MDP(16) + (-t449 * t584 - t450 * t585 + (t491 * t585 - t494 * t584) * qJD(2)) * MDP(17) + (qJD(2) * t622 + t403 * t584 - t469 * t675 - t585 * t643 + t652) * MDP(18) + (t509 * t576 - t577 * t682 - t675 * t712) * MDP(24) + (t437 * t675 + t512 * t576 + t577 * t683) * MDP(25) + (-t512 * t695 + t509 * t380 - t682 * t416 + (-t587 * t683 - t591 * t617) * t722) * MDP(31) + (-t512 * t390 + t509 * t379 - t682 * t418 + (t587 * t617 - t591 * t683) * t722) * MDP(32); (t584 * t659 + t585 * t658 - t529) * MDP(15) + t626 * MDP(16) + (-t491 ^ 2 - t494 ^ 2) * MDP(17) + (t426 * t494 + t427 * t491 + t608 + t647) * MDP(18) + (t395 + t725) * MDP(24) + (-t611 + t724) * MDP(25) + (t390 + t727) * MDP(31) + (-t695 + t726) * MDP(32) + (-t494 * MDP(15) + t491 * MDP(16) + ((MDP(15) * t584 + MDP(16) * t585) * t593 + (t584 * MDP(16) - MDP(18) * t705) * t589) * qJD(1)) * qJD(2) - (MDP(31) * t587 + MDP(32) * t591) * t734; (-t384 * t577 + t717) * MDP(24) + (-t383 * t577 + t716) * MDP(25) + (-pkin(5) * t380 - t384 * t416 + t604 * t587 - t591 * t708 + t729) * MDP(31) + (-pkin(5) * t379 - t384 * t418 + t587 * t708 + t604 * t591 - t728) * MDP(32) + t742; t418 * t416 * MDP(26) + (-t416 ^ 2 + t418 ^ 2) * MDP(27) + (t655 + t694) * MDP(28) + (-t543 + t693) * MDP(29) + t393 * MDP(30) + (-g(1) * t457 + t375 * t722 - t381 * t418 + t372) * MDP(31) + (g(1) * t458 + t374 * t722 + t381 * t416) * MDP(32) + (MDP(29) * t667 + MDP(31) * t623 + MDP(32) * t605) * t591 + (MDP(28) * t667 + (qJD(6) * t577 + t611) * MDP(29) + t605 * MDP(31) + (-t373 - t623) * MDP(32)) * t587;];
tau  = t1;
