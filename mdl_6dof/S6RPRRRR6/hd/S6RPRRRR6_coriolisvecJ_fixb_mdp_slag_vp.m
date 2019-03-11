% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:14:51
% EndTime: 2019-03-09 07:15:09
% DurationCPUTime: 11.17s
% Computational Cost: add. (9199->523), mult. (23678->685), div. (0->0), fcn. (19062->10), ass. (0->218)
t578 = sin(pkin(11));
t579 = cos(pkin(11));
t583 = sin(qJ(3));
t587 = cos(qJ(3));
t548 = t578 * t587 + t579 * t583;
t536 = t548 * qJD(1);
t582 = sin(qJ(4));
t586 = cos(qJ(4));
t631 = t586 * qJD(3);
t511 = t536 * t582 - t631;
t513 = qJD(3) * t582 + t536 * t586;
t581 = sin(qJ(5));
t585 = cos(qJ(5));
t449 = t585 * t511 + t513 * t581;
t584 = cos(qJ(6));
t580 = sin(qJ(6));
t604 = t511 * t581 - t585 * t513;
t680 = t604 * t580;
t407 = -t584 * t449 + t680;
t538 = t548 * qJD(3);
t526 = qJD(1) * t538;
t605 = t449 * t580 + t584 * t604;
t718 = t526 * MDP(33) + (-t407 ^ 2 + t605 ^ 2) * MDP(30) + t407 * MDP(29) * t605;
t640 = qJD(1) * t587;
t568 = t579 * t640;
t641 = qJD(1) * t583;
t625 = t578 * t641;
t535 = t568 - t625;
t691 = qJD(4) + qJD(5);
t717 = t535 - t691;
t565 = qJD(3) * t568;
t525 = -qJD(3) * t625 + t565;
t637 = qJD(4) * t582;
t461 = qJD(4) * t631 + t586 * t525 - t536 * t637;
t462 = qJD(4) * t513 + t525 * t582;
t634 = qJD(5) * t585;
t635 = qJD(5) * t581;
t393 = t585 * t461 - t581 * t462 - t511 * t634 - t513 * t635;
t589 = qJD(5) * t604 - t461 * t581 - t585 * t462;
t632 = qJD(6) * t584;
t627 = t584 * t393 - t449 * t632 + t580 * t589;
t633 = qJD(6) * t580;
t362 = t604 * t633 + t627;
t530 = qJD(4) - t535;
t524 = qJD(5) + t530;
t619 = t393 * t580 - t584 * t589;
t590 = qJD(6) * t605 - t619;
t520 = qJD(6) + t524;
t707 = t520 * t605;
t708 = t407 * t520;
t716 = t526 * MDP(26) + (-t449 ^ 2 + t604 ^ 2) * MDP(23) + (t449 * t524 + t393) * MDP(24) + (-t524 * t604 + t589) * MDP(25) - t449 * MDP(22) * t604 + (t590 - t707) * MDP(32) + (t362 - t708) * MDP(31) + t718;
t570 = -pkin(2) * t579 - pkin(1);
t558 = qJD(1) * t570 + qJD(2);
t467 = -pkin(3) * t535 - pkin(8) * t536 + t558;
t687 = pkin(7) + qJ(2);
t559 = t687 * t578;
t549 = qJD(1) * t559;
t560 = t687 * t579;
t550 = qJD(1) * t560;
t498 = -t583 * t549 + t587 * t550;
t492 = qJD(3) * pkin(8) + t498;
t434 = t586 * t467 - t492 * t582;
t416 = -pkin(9) * t513 + t434;
t397 = pkin(4) * t530 + t416;
t435 = t467 * t582 + t492 * t586;
t417 = -pkin(9) * t511 + t435;
t413 = t585 * t417;
t379 = t397 * t581 + t413;
t710 = pkin(10) * t449;
t369 = t379 - t710;
t367 = t369 * t633;
t695 = -t549 * t587 - t583 * t550;
t491 = -qJD(3) * pkin(3) - t695;
t442 = pkin(4) * t511 + t491;
t409 = pkin(5) * t449 + t442;
t701 = -t409 * t407 + t367;
t551 = t581 * t582 - t585 * t586;
t648 = t717 * t551;
t660 = t581 * t586;
t552 = t582 * t585 + t660;
t647 = t717 * t552;
t547 = t578 * t583 - t587 * t579;
t593 = t547 * qJD(2);
t453 = -qJD(1) * t593 + qJD(3) * t695;
t476 = pkin(3) * t526 - pkin(8) * t525;
t472 = t586 * t476;
t592 = -qJD(4) * t435 - t453 * t582 + t472;
t376 = pkin(4) * t526 - pkin(9) * t461 + t592;
t636 = qJD(4) * t586;
t597 = t586 * t453 + t467 * t636 + t582 * t476 - t492 * t637;
t382 = -pkin(9) * t462 + t597;
t621 = t585 * t376 - t581 * t382;
t591 = -qJD(5) * t379 + t621;
t356 = pkin(5) * t526 - pkin(10) * t393 + t591;
t612 = -t581 * t376 - t585 * t382 - t397 * t634 + t417 * t635;
t357 = pkin(10) * t589 - t612;
t622 = t584 * t356 - t580 * t357;
t700 = t409 * t605 + t622;
t493 = pkin(3) * t536 - pkin(8) * t535;
t482 = t586 * t493;
t688 = pkin(8) + pkin(9);
t626 = qJD(4) * t688;
t713 = pkin(4) * t536 - t582 * t695 + t482 + (-pkin(9) * t535 + t626) * t586;
t650 = t582 * t493 + t586 * t695;
t671 = t535 * t582;
t712 = -pkin(9) * t671 + t582 * t626 + t650;
t711 = t637 - t671;
t709 = pkin(10) * t604;
t705 = (t578 ^ 2 + t579 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t501 = -t551 * t580 + t552 * t584;
t654 = qJD(6) * t501 + t580 * t648 - t584 * t647;
t703 = t548 * qJD(2);
t623 = t548 * t636;
t537 = t547 * qJD(3);
t670 = t537 * t582;
t702 = t623 - t670;
t699 = t442 * t449 + t612;
t698 = t442 * t604 + t591;
t484 = t552 * t548;
t495 = pkin(3) * t547 - pkin(8) * t548 + t570;
t487 = t586 * t495;
t510 = -t559 * t583 + t560 * t587;
t667 = t548 * t586;
t423 = pkin(4) * t547 - pkin(9) * t667 - t510 * t582 + t487;
t503 = t586 * t510;
t649 = t582 * t495 + t503;
t668 = t548 * t582;
t436 = -pkin(9) * t668 + t649;
t652 = t581 * t423 + t585 * t436;
t611 = pkin(4) * t711 - t498;
t694 = t713 * t585;
t509 = t559 * t587 + t583 * t560;
t561 = t688 * t582;
t562 = t688 * t586;
t645 = -t581 * t561 + t585 * t562;
t693 = t561 * t634 + t562 * t635 + t581 * t713 + t712 * t585;
t500 = t584 * t551 + t552 * t580;
t655 = -qJD(6) * t500 + t580 * t647 + t584 * t648;
t690 = -t501 * t526 - t520 * t655;
t689 = -t524 * t648 - t526 * t552;
t411 = t581 * t417;
t378 = t585 * t397 - t411;
t368 = t378 + t709;
t366 = pkin(5) * t524 + t368;
t685 = t366 * t584;
t684 = t407 * t536;
t683 = t605 * t536;
t682 = t449 * t536;
t681 = t604 * t536;
t679 = t461 * t582;
t677 = t511 * t530;
t676 = t511 * t536;
t675 = t513 * t530;
t674 = t513 * t536;
t669 = t537 * t586;
t662 = t581 * t526;
t661 = t581 * t584;
t659 = t582 * t526;
t658 = t584 * t369;
t657 = t584 * t526;
t515 = t586 * t526;
t656 = t585 * t416 - t411;
t651 = -pkin(5) * t647 + t611;
t639 = qJD(3) * t583;
t638 = qJD(3) * t587;
t630 = qJD(1) * qJD(2);
t574 = -pkin(4) * t586 - pkin(3);
t624 = t548 * t637;
t468 = -qJD(3) * t509 - t593;
t494 = pkin(3) * t538 + pkin(8) * t537;
t483 = t586 * t494;
t386 = pkin(9) * t669 + pkin(4) * t538 - t468 * t582 + t483 + (-t503 + (pkin(9) * t548 - t495) * t582) * qJD(4);
t596 = t586 * t468 + t582 * t494 + t495 * t636 - t510 * t637;
t389 = -pkin(9) * t702 + t596;
t620 = t585 * t386 - t389 * t581;
t618 = -t416 * t581 - t413;
t617 = t585 * t423 - t436 * t581;
t615 = -t585 * t561 - t562 * t581;
t614 = t530 * t586;
t613 = qJD(6) * t366 + t357;
t454 = qJD(1) * t703 - t549 * t639 + t550 * t638;
t469 = -t559 * t639 + t560 * t638 + t703;
t610 = -t500 * t526 - t520 * t654;
t609 = t524 * t647 - t551 * t526;
t470 = pkin(4) * t668 + t509;
t479 = -pkin(10) * t551 + t645;
t608 = pkin(5) * t536 + pkin(10) * t648 + t645 * qJD(5) + qJD(6) * t479 - t581 * t712 + t694;
t478 = -pkin(10) * t552 + t615;
t607 = -pkin(10) * t647 - qJD(6) * t478 + t693;
t361 = t580 * t366 + t658;
t485 = t551 * t548;
t438 = t584 * t484 - t485 * t580;
t439 = -t484 * t580 - t485 * t584;
t601 = -t530 * t711 + t515;
t441 = pkin(4) * t702 + t469;
t599 = -t624 - t669;
t422 = pkin(4) * t462 + t454;
t598 = -pkin(8) * t526 + t491 * t530;
t595 = t581 * t386 + t585 * t389 + t423 * t634 - t436 * t635;
t573 = pkin(4) * t585 + pkin(5);
t523 = pkin(5) * t551 + t574;
t499 = t526 * t547;
t440 = pkin(4) * t513 - pkin(5) * t604;
t437 = pkin(5) * t484 + t470;
t400 = -t537 * t660 - t581 * t624 - t635 * t668 + (t667 * t691 - t670) * t585;
t399 = -t484 * t691 + t551 * t537;
t387 = pkin(5) * t400 + t441;
t383 = -pkin(10) * t484 + t652;
t380 = pkin(5) * t547 + pkin(10) * t485 + t617;
t373 = -pkin(5) * t589 + t422;
t371 = t656 + t709;
t370 = t618 + t710;
t365 = qJD(6) * t439 + t399 * t580 + t584 * t400;
t364 = -qJD(6) * t438 + t399 * t584 - t400 * t580;
t360 = -t369 * t580 + t685;
t359 = -pkin(10) * t400 + t595;
t358 = pkin(5) * t538 - pkin(10) * t399 - qJD(5) * t652 + t620;
t1 = [(-MDP(10) * t537 - MDP(11) * t538 - MDP(13) * t469 - MDP(14) * t468) * qJD(3) + 0.2e1 * t630 * t705 + (t525 * t570 - t537 * t558) * MDP(14) + (-t379 * t538 + t470 * t393 + t442 * t399 - t422 * t485 - t441 * t604 - t524 * t595 - t526 * t652 + t547 * t612) * MDP(28) + (t393 * t547 + t399 * t524 - t485 * t526 - t538 * t604) * MDP(24) + (-t393 * t485 - t399 * t604) * MDP(22) + (t526 * t570 + t538 * t558) * MDP(13) + (-t462 * t547 - t511 * t538 - t530 * t702 - t548 * t659) * MDP(18) + (-t361 * t538 + t437 * t362 + t409 * t364 + t367 * t547 + t373 * t439 - t387 * t605 + (-(-qJD(6) * t383 + t358) * t520 - t380 * t526 - t356 * t547) * t580 + (-(qJD(6) * t380 + t359) * t520 - t383 * t526 - t613 * t547) * t584) * MDP(35) + (t362 * t547 + t364 * t520 + t439 * t526 - t538 * t605) * MDP(31) + (t362 * t439 - t364 * t605) * MDP(29) + (-t393 * t484 - t399 * t449 + t400 * t604 - t485 * t589) * MDP(23) + (-t400 * t524 - t449 * t538 - t484 * t526 + t547 * t589) * MDP(25) + (t620 * t524 + t617 * t526 + t621 * t547 + t378 * t538 + t441 * t449 - t470 * t589 + t422 * t484 + t442 * t400 + (-t379 * t547 - t524 * t652) * qJD(5)) * MDP(27) + (-(-t511 * t586 - t513 * t582) * t537 + (-t679 - t462 * t586 + (t511 * t582 - t513 * t586) * qJD(4)) * t548) * MDP(16) + ((-t510 * t636 + t483) * t530 + t487 * t526 + (-t492 * t636 + t472) * t547 + t434 * t538 + t469 * t511 + t509 * t462 + t491 * t623 + ((-qJD(4) * t495 - t468) * t530 - t510 * t526 + (-qJD(4) * t467 - t453) * t547 + t454 * t548 - t491 * t537) * t582) * MDP(20) + (-t525 * t547 - t526 * t548 - t535 * t537 - t536 * t538) * MDP(9) + (t525 * t548 - t536 * t537) * MDP(8) + (t461 * t667 + t513 * t599) * MDP(15) + (-t435 * t538 + t454 * t667 + t509 * t461 + t469 * t513 + t491 * t599 - t526 * t649 - t530 * t596 - t547 * t597) * MDP(21) + (t461 * t547 + t513 * t538 + t515 * t548 + t530 * t599) * MDP(17) + (t524 * t538 + t499) * MDP(26) + (t520 * t538 + t499) * MDP(33) + (t530 * t538 + t499) * MDP(19) + (-t362 * t438 + t364 * t407 + t365 * t605 + t439 * t590) * MDP(30) + ((t358 * t584 - t359 * t580) * t520 + (t380 * t584 - t383 * t580) * t526 + t622 * t547 + t360 * t538 - t387 * t407 - t437 * t590 + t373 * t438 + t409 * t365 + ((-t380 * t580 - t383 * t584) * t520 - t361 * t547) * qJD(6)) * MDP(34) + (-t365 * t520 + t407 * t538 - t438 * t526 + t547 * t590) * MDP(32); t565 * MDP(14) + (t601 - t676) * MDP(20) + (-t530 ^ 2 * t586 - t659 - t674) * MDP(21) + (t609 - t682) * MDP(27) + (t681 + t689) * MDP(28) + (t610 + t684) * MDP(34) + (t683 + t690) * MDP(35) + ((t578 * t640 + t579 * t641 + t536) * MDP(13) + (t535 - t625) * MDP(14)) * qJD(3) - qJD(1) ^ 2 * t705; (t609 + t682) * MDP(25) + (t683 - t690) * MDP(31) + (t610 - t684) * MDP(32) + (t601 + t676) * MDP(18) + ((t461 - t677) * t586 + (-t462 - t675) * t582) * MDP(16) + (t513 * t614 + t679) * MDP(15) + (t681 - t689) * MDP(24) + (t530 * t614 + t659 - t674) * MDP(17) + (t362 * t501 - t605 * t655) * MDP(29) + (-(t478 * t580 + t479 * t584) * t526 + t523 * t362 + t373 * t501 + (t580 * t608 + t584 * t607) * t520 + t655 * t409 - t651 * t605) * MDP(35) + (-t362 * t500 + t407 * t655 + t501 * t590 + t605 * t654) * MDP(30) + (t393 * t552 - t604 * t648) * MDP(22) + (-t393 * t551 - t449 * t648 + t552 * t589 - t604 * t647) * MDP(23) + (-pkin(3) * t461 + t454 * t582 - t498 * t513 + (pkin(8) * t637 + t650) * t530 + t598 * t586) * MDP(21) + (t574 * t393 + t422 * t552 + t648 * t442 + t524 * t693 - t645 * t526 - t604 * t611) * MDP(28) + ((t478 * t584 - t479 * t580) * t526 - t523 * t590 + t373 * t500 + (t580 * t607 - t584 * t608) * t520 + t654 * t409 - t651 * t407) * MDP(34) + (t615 * t526 - t574 * t589 + t422 * t551 + (-t562 * t634 + (qJD(5) * t561 + t712) * t581 - t694) * t524 + t611 * t449 - t647 * t442) * MDP(27) + (-pkin(3) * t462 - t454 * t586 - t498 * t511 + (-pkin(8) * t636 - t482) * t530 + (t530 * t695 + t598) * t582) * MDP(20) + (-t535 * t558 + t547 * t630) * MDP(14) + (t565 + (-t535 - t625) * qJD(3)) * MDP(10) - t535 ^ 2 * MDP(9) + (qJD(3) * t498 - t454) * MDP(13) + (-t558 * MDP(13) - t530 * MDP(19) - t434 * MDP(20) + t435 * MDP(21) - t524 * MDP(26) - t378 * MDP(27) + t379 * MDP(28) - t520 * MDP(33) - t360 * MDP(34) + t361 * MDP(35) - MDP(8) * t535 + MDP(9) * t536) * t536; (-t618 * t524 + (-t449 * t513 - t524 * t635 + t526 * t585) * pkin(4) + t698) * MDP(27) + (t656 * t524 + (t513 * t604 - t524 * t634 - t662) * pkin(4) + t699) * MDP(28) + (t440 * t605 + (-t573 * t526 - t356 + (t370 - (-qJD(5) - qJD(6)) * t581 * pkin(4)) * t520) * t580 + (-pkin(4) * t662 + (-pkin(4) * t634 - qJD(6) * t573 + t371) * t520 - t613) * t584 + t701) * MDP(35) + (-t462 + t675) * MDP(18) + (t461 + t677) * MDP(17) + (t573 * t657 - (t370 * t584 - t371 * t580) * t520 + t440 * t407 + (-t580 * t662 + (-t580 * t585 - t661) * t520 * qJD(5)) * pkin(4) + ((-pkin(4) * t661 - t573 * t580) * t520 - t361) * qJD(6) + t700) * MDP(34) + t526 * MDP(19) + (-t511 ^ 2 + t513 ^ 2) * MDP(16) + (t434 * t530 + t491 * t511 - t597) * MDP(21) + (t435 * t530 - t491 * t513 + t592) * MDP(20) + t513 * t511 * MDP(15) + t716; (t379 * t524 + t698) * MDP(27) + (t378 * t524 + t699) * MDP(28) + (-(-t368 * t580 - t658) * t520 - t361 * qJD(6) + (-t407 * t604 - t520 * t633 + t657) * pkin(5) + t700) * MDP(34) + ((-t369 * t520 - t356) * t580 + (t368 * t520 - t613) * t584 + (-t520 * t632 - t580 * t526 - t604 * t605) * pkin(5) + t701) * MDP(35) + t716; (t627 - t708) * MDP(31) + (-t619 - t707) * MDP(32) + (t361 * t520 + t700) * MDP(34) + (-t580 * t356 - t584 * t357 + t360 * t520 + t701) * MDP(35) + (MDP(31) * t680 + MDP(32) * t605 - MDP(34) * t361 - MDP(35) * t685) * qJD(6) + t718;];
tauc  = t1;
