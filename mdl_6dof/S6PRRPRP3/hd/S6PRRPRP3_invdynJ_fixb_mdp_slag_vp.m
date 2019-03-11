% Calculate vector of inverse dynamics joint torques for
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:38:37
% EndTime: 2019-03-08 21:38:50
% DurationCPUTime: 10.17s
% Computational Cost: add. (5791->591), mult. (13274->777), div. (0->0), fcn. (10410->14), ass. (0->246)
t576 = sin(qJ(3));
t578 = cos(qJ(3));
t626 = pkin(3) * t576 - qJ(4) * t578;
t505 = qJD(3) * t626 - qJD(4) * t576;
t571 = sin(pkin(11));
t573 = cos(pkin(11));
t577 = sin(qJ(2));
t668 = qJD(3) * t576;
t660 = pkin(8) * t668;
t572 = sin(pkin(6));
t674 = qJD(1) * t572;
t579 = cos(qJ(2));
t690 = t578 * t579;
t681 = t573 * t505 + t571 * t660 - (-t571 * t690 + t573 * t577) * t674;
t696 = t571 * t577;
t738 = t571 * t505 - (t573 * t690 + t696) * t674;
t691 = t573 * t578;
t621 = pkin(4) * t576 - pkin(9) * t691;
t737 = qJD(3) * t621 + t681;
t692 = t573 * t576;
t695 = t571 * t578;
t736 = (-pkin(8) * t692 - pkin(9) * t695) * qJD(3) + t738;
t669 = qJD(3) * t573;
t671 = qJD(2) * t576;
t523 = -t571 * t671 + t669;
t524 = qJD(3) * t571 + t573 * t671;
t575 = sin(qJ(5));
t717 = cos(qJ(5));
t448 = -t717 * t523 + t524 * t575;
t735 = t448 ^ 2;
t612 = -t575 * t523 - t524 * t717;
t718 = t612 ^ 2;
t670 = qJD(2) * t578;
t558 = -qJD(5) + t670;
t734 = t448 * t558;
t733 = t558 * t612;
t657 = t717 * t573;
t610 = -t575 * t571 + t657;
t649 = qJD(5) * t717;
t666 = qJD(5) * t575;
t726 = -t571 * t666 + t573 * t649;
t679 = -t610 * t670 + t726;
t528 = t571 * t717 + t575 * t573;
t512 = t528 * qJD(5);
t599 = t578 * t528;
t678 = -qJD(2) * t599 + t512;
t663 = qJD(2) * qJD(3);
t646 = t578 * t663;
t662 = qJDD(2) * t576;
t732 = t646 + t662;
t627 = pkin(3) * t578 + qJ(4) * t576;
t535 = -pkin(2) - t627;
t522 = t573 * t535;
t455 = -pkin(9) * t692 + t522 + (-pkin(8) * t571 - pkin(4)) * t578;
t489 = pkin(8) * t691 + t571 * t535;
t697 = t571 * t576;
t469 = -pkin(9) * t697 + t489;
t731 = t455 * t649 - t469 * t666 + t737 * t575 + t736 * t717;
t529 = t626 * qJD(2);
t531 = qJD(2) * pkin(8) + t577 * t674;
t709 = cos(pkin(6));
t642 = qJD(1) * t709;
t600 = t576 * t531 - t578 * t642;
t437 = t573 * t529 + t571 * t600;
t415 = qJD(2) * t621 + t437;
t438 = t571 * t529 - t573 * t600;
t654 = t571 * t670;
t425 = -pkin(9) * t654 + t438;
t711 = pkin(9) + qJ(4);
t538 = t711 * t571;
t539 = t711 * t573;
t611 = -t538 * t717 - t575 * t539;
t730 = qJD(4) * t610 + qJD(5) * t611 - t575 * t415 - t717 * t425;
t471 = -t575 * t538 + t539 * t717;
t729 = qJD(4) * t528 + qJD(5) * t471 + t415 * t717 - t575 * t425;
t728 = t575 * t455 + t717 * t469;
t667 = qJD(3) * t578;
t651 = t571 * t667;
t518 = pkin(4) * t651 + pkin(8) * t667;
t655 = t579 * t674;
t635 = t576 * t655;
t727 = t518 - t635;
t566 = t578 * qJDD(2);
t725 = t576 * t663 - t566;
t724 = MDP(21) + MDP(23);
t723 = MDP(22) - MDP(25);
t703 = qJDD(3) * pkin(3);
t722 = qJDD(4) - t703;
t708 = cos(pkin(10));
t629 = t709 * t708;
t707 = sin(pkin(10));
t508 = t577 * t629 + t579 * t707;
t628 = t709 * t707;
t510 = -t577 * t628 + t579 * t708;
t694 = t572 * t577;
t514 = t576 * t694 - t578 * t709;
t643 = t572 * t707;
t644 = t572 * t708;
t603 = g(3) * t514 + g(2) * (t508 * t576 + t578 * t644) + g(1) * (t510 * t576 - t578 * t643);
t526 = qJDD(5) + t725;
t568 = pkin(11) + qJ(5);
t564 = sin(t568);
t721 = -t471 * t526 - t603 * t564;
t580 = qJD(3) ^ 2;
t664 = qJD(1) * qJD(2);
t648 = t577 * t664;
t693 = t572 * t579;
t625 = -qJDD(1) * t693 + t572 * t648;
t507 = t577 * t707 - t579 * t629;
t509 = t577 * t708 + t579 * t628;
t632 = g(1) * t509 + g(2) * t507;
t720 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t580 + t572 * (-g(3) * t579 + t648) - t625 + t632;
t661 = t571 * qJDD(3);
t588 = t573 * t732 + t661;
t677 = t732 * t571;
t624 = qJDD(3) * t573 - t677;
t393 = -qJD(5) * t612 + t575 * t588 - t717 * t624;
t719 = qJD(5) * t728 + t736 * t575 - t737 * t717;
t715 = pkin(5) * t526;
t710 = qJD(2) * pkin(2);
t705 = qJ(6) * t526;
t552 = t576 * t642;
t486 = t578 * t531 + t552;
t476 = qJD(3) * qJ(4) + t486;
t487 = qJD(2) * t535 - t655;
t416 = -t476 * t571 + t573 * t487;
t397 = -pkin(4) * t670 - pkin(9) * t524 + t416;
t417 = t573 * t476 + t571 * t487;
t403 = pkin(9) * t523 + t417;
t380 = t575 * t397 + t403 * t717;
t702 = t380 * t558;
t701 = t448 * t612;
t699 = t564 * t578;
t565 = cos(t568);
t698 = t565 * t578;
t689 = qJDD(1) - g(3);
t688 = qJ(6) * t668 - qJD(6) * t578 + t731;
t687 = -pkin(5) * t668 + t719;
t497 = qJDD(2) * pkin(8) + (qJDD(1) * t577 + t579 * t664) * t572;
t640 = qJDD(1) * t709;
t630 = t576 * t640;
t407 = t630 + qJDD(3) * qJ(4) + t578 * t497 + (qJD(4) - t600) * qJD(3);
t423 = qJD(2) * t505 + qJDD(2) * t535 + t625;
t386 = t573 * t407 + t571 * t423;
t457 = pkin(4) * t654 + t486;
t685 = pkin(5) * t678 - qJ(6) * t679 - qJD(6) * t528 - t457;
t684 = -qJ(6) * t671 + t730;
t683 = pkin(5) * t671 + t729;
t638 = t573 * t660;
t680 = -t638 + t738;
t676 = pkin(2) * t693 + pkin(8) * t694;
t530 = pkin(4) * t697 + t576 * pkin(8);
t569 = t576 ^ 2;
t675 = -t578 ^ 2 + t569;
t673 = qJD(1) * t576;
t672 = qJD(2) * t572;
t379 = t397 * t717 - t575 * t403;
t665 = qJD(6) - t379;
t659 = t564 * t693;
t562 = pkin(4) * t573 + pkin(3);
t658 = pkin(4) * t571 + pkin(8);
t656 = qJ(4) * t668;
t653 = t577 * t672;
t652 = t579 * t672;
t645 = t579 * t663;
t641 = -qJD(3) * pkin(3) + qJD(4);
t385 = -t571 * t407 + t573 * t423;
t488 = -pkin(8) * t695 + t522;
t639 = qJD(2) * t488 + t416;
t637 = -qJDD(2) * t488 - t385;
t378 = pkin(4) * t725 - t588 * pkin(9) + t385;
t383 = pkin(9) * t624 + t386;
t636 = -t717 * t378 + t575 * t383 + t397 * t666 + t403 * t649;
t631 = g(1) * t510 + g(2) * t508;
t623 = t562 * t578 + t576 * t711;
t581 = qJD(2) ^ 2;
t622 = qJDD(2) * t579 - t577 * t581;
t440 = t512 * t576 + t575 * t651 - t657 * t667;
t441 = qJD(3) * t599 + t576 * t726;
t501 = t610 * t576;
t387 = pkin(5) * t441 + qJ(6) * t440 - qJD(6) * t501 + t518;
t620 = -t387 + t635;
t616 = t455 * t717 - t575 * t469;
t515 = t576 * t709 + t578 * t694;
t461 = -t515 * t571 - t573 * t693;
t462 = t515 * t573 - t571 * t693;
t614 = t461 * t717 - t575 * t462;
t406 = t575 * t461 + t462 * t717;
t608 = t573 * t662 + t661;
t607 = t575 * t378 + t717 * t383 + t397 * t649 - t403 * t666;
t392 = -t523 * t649 + t524 * t666 - t575 * t624 - t717 * t588;
t428 = -t507 * t699 - t508 * t565;
t430 = -t509 * t699 - t510 * t565;
t480 = -t565 * t694 + t578 * t659;
t605 = g(1) * t430 + g(2) * t428 + g(3) * t480;
t429 = -t507 * t698 + t508 * t564;
t431 = -t509 * t698 + t510 * t564;
t481 = (t564 * t577 + t565 * t690) * t572;
t604 = -g(1) * t431 - g(2) * t429 - g(3) * t481;
t464 = t508 * t578 - t576 * t644;
t466 = t510 * t578 + t576 * t643;
t602 = g(1) * t466 + g(2) * t464 + g(3) * t515;
t598 = g(3) * t693 - t632;
t597 = -g(3) * t694 - t631;
t595 = -qJD(3) * t552 - t576 * t497 - t531 * t667 + t578 * t640;
t409 = -t595 + t722;
t596 = -t409 + t603;
t472 = t600 + t641;
t591 = t596 + t703;
t411 = t464 * t564 - t507 * t565;
t413 = t466 * t564 - t509 * t565;
t453 = t515 * t564 + t565 * t693;
t590 = g(1) * t413 + g(2) * t411 + g(3) * t453 - t636;
t589 = t526 * t611 + t565 * t603;
t532 = -t655 - t710;
t587 = -pkin(8) * qJDD(3) + (t532 + t655 - t710) * qJD(3);
t439 = -t523 * pkin(4) + t472;
t412 = t464 * t565 + t507 * t564;
t414 = t466 * t565 + t509 * t564;
t454 = t515 * t565 - t659;
t586 = -g(1) * t414 - g(2) * t412 - g(3) * t454 + t607;
t390 = t448 * pkin(5) + qJ(6) * t612 + t439;
t585 = -t390 * t612 + qJDD(6) - t590;
t583 = t595 + t603;
t394 = -pkin(4) * t624 + t409;
t369 = t393 * pkin(5) + t392 * qJ(6) + qJD(6) * t612 + t394;
t504 = t509 * pkin(2);
t503 = t507 * pkin(2);
t500 = t528 * t576;
t468 = qJD(3) * t515 + t576 * t652;
t467 = -qJD(3) * t514 + t578 * t652;
t446 = -pkin(5) * t610 - qJ(6) * t528 - t562;
t436 = t467 * t573 + t571 * t653;
t435 = -t467 * t571 + t573 * t653;
t424 = pkin(5) * t500 - qJ(6) * t501 + t530;
t402 = -pkin(5) * t612 + qJ(6) * t448;
t401 = t578 * pkin(5) - t616;
t400 = -qJ(6) * t578 + t728;
t384 = -t392 - t734;
t377 = qJD(5) * t406 - t435 * t717 + t575 * t436;
t376 = qJD(5) * t614 + t575 * t435 + t436 * t717;
t375 = -t558 * qJ(6) + t380;
t371 = t558 * pkin(5) + t665;
t368 = qJDD(6) + t636 - t715;
t367 = -qJD(6) * t558 + t607 + t705;
t1 = [t689 * MDP(1) + (-qJD(3) * t468 - qJDD(3) * t514) * MDP(10) + (-qJD(3) * t467 - qJDD(3) * t515) * MDP(11) + (-t461 * t566 - t468 * t523 - t514 * t624 + (-t435 * t578 + t461 * t668) * qJD(2)) * MDP(12) + (t462 * t566 + t468 * t524 + t514 * t608 + (t436 * t578 + (-t462 * t576 + t514 * t691) * qJD(3)) * qJD(2)) * MDP(13) + (-t435 * t524 + t436 * t523 - t461 * t588 + t462 * t624) * MDP(14) + (t385 * t461 + t386 * t462 + t409 * t514 + t416 * t435 + t417 * t436 + t468 * t472 - g(3)) * MDP(15) + (-t376 * t448 - t377 * t612 + t392 * t614 - t393 * t406) * MDP(24) + (t367 * t406 - t368 * t614 + t369 * t514 + t371 * t377 + t375 * t376 + t390 * t468 - g(3)) * MDP(26) + t724 * (t377 * t558 + t514 * t393 + t468 * t448 + t526 * t614) + t723 * (t376 * t558 - t392 * t514 - t406 * t526 - t468 * t612) + (t622 * MDP(3) + (-qJDD(2) * t577 - t579 * t581) * MDP(4) + (-t576 * t645 + t578 * t622) * MDP(10) + (-t576 * t622 - t578 * t645) * MDP(11)) * t572; (-t576 * t720 + t587 * t578) * MDP(11) + (t597 * t571 + (-pkin(8) * t624 + qJD(3) * t639 + t409 * t571 + t523 * t655) * t576 + ((-pkin(8) * t523 + t472 * t571) * qJD(3) - t681 * qJD(2) - t598 * t573 + t637) * t578) * MDP(12) + (t597 * t573 + (-t524 * t655 + t409 * t573 + (-qJD(2) * t489 - t417) * qJD(3) + t608 * pkin(8)) * t576 + (t489 * qJDD(2) + t386 + (pkin(8) * t524 + t472 * t573) * qJD(3) + t598 * t571 + (t638 + t680) * qJD(2)) * t578) * MDP(13) + (t489 * t624 - t488 * t661 - t681 * t524 + t680 * t523 + (-t417 * t571 - t573 * t639) * t667 + (-t386 * t571 + t573 * t637 - t598) * t576) * MDP(14) + (t386 * t489 + t385 * t488 - g(1) * (-t509 * t627 - t504) - g(2) * (-t507 * t627 - t503) - g(3) * t676 + (-g(3) * t627 - t472 * t673) * t693 + t680 * t417 + t681 * t416 + (t409 * t576 + t472 * t667 - t631) * pkin(8)) * MDP(15) + (-t392 * t501 + t440 * t612) * MDP(16) + (t392 * t500 - t393 * t501 + t440 * t448 + t441 * t612) * MDP(17) + (t392 * t578 + t440 * t558 + t501 * t526 - t612 * t668) * MDP(18) + (t393 * t578 + t441 * t558 - t448 * t668 - t500 * t526) * MDP(19) + (-t526 * t578 - t558 * t668) * MDP(20) + (t379 * t668 + t530 * t393 + t394 * t500 + t439 * t441 + t448 * t727 + t616 * t526 + t558 * t719 + t636 * t578 + t604) * MDP(21) + (-t380 * t668 - t530 * t392 + t394 * t501 - t439 * t440 - t526 * t728 + t558 * t731 + t607 * t578 - t612 * t727 + t605) * MDP(22) + (t368 * t578 + t369 * t500 - t371 * t668 + t390 * t441 + t393 * t424 - t401 * t526 - t448 * t620 + t558 * t687 + t604) * MDP(23) + (-t367 * t500 + t368 * t501 - t371 * t440 - t375 * t441 - t392 * t401 - t393 * t400 - t448 * t688 - t576 * t598 - t612 * t687) * MDP(24) + (-t367 * t578 - t369 * t501 + t375 * t668 + t390 * t440 + t392 * t424 + t400 * t526 - t558 * t688 - t612 * t620 - t605) * MDP(25) + (t367 * t400 + t369 * t424 + t390 * t387 + t368 * t401 - g(1) * (pkin(5) * t431 + qJ(6) * t430 - t509 * t623 + t510 * t658 - t504) - g(2) * (pkin(5) * t429 + qJ(6) * t428 - t507 * t623 + t508 * t658 - t503) - g(3) * (pkin(5) * t481 + qJ(6) * t480 + t676) + t688 * t375 + t687 * t371 + (-g(3) * pkin(4) * t696 + (-g(3) * t623 - t390 * t673) * t579) * t572) * MDP(26) + qJDD(2) * MDP(2) + (t689 * t693 + t632) * MDP(3) + (-t689 * t694 + t631) * MDP(4) + (qJDD(2) * t569 + 0.2e1 * t576 * t646) * MDP(5) + 0.2e1 * (t566 * t576 - t663 * t675) * MDP(6) + (qJDD(3) * t576 + t578 * t580) * MDP(7) + (qJDD(3) * t578 - t576 * t580) * MDP(8) + (t587 * t576 + t578 * t720) * MDP(10); MDP(7) * t662 + MDP(8) * t566 + qJDD(3) * MDP(9) + (t486 * qJD(3) - t532 * t671 + t583) * MDP(10) + (-t630 + (-qJD(2) * t532 - t497) * t578 + t602) * MDP(11) + (t571 * qJ(4) * t566 - pkin(3) * t677 + t486 * t523 + t591 * t573 + (-t416 * t576 + t437 * t578 + (-t656 + (qJD(4) - t472) * t578) * t571) * qJD(2)) * MDP(12) + (-t486 * t524 - t626 * t573 * qJDD(2) - t591 * t571 + (t417 * t576 - t438 * t578 + (-t656 + (-t472 + t641) * t578) * t573) * qJD(2)) * MDP(13) + (t437 * t524 - t438 * t523 + (qJ(4) * t624 + qJD(4) * t523 + t416 * t670 + t386) * t573 + (qJ(4) * t588 + qJD(4) * t524 + t417 * t670 - t385) * t571 - t602) * MDP(14) + (-t416 * t437 - t417 * t438 - t472 * t486 + (-t416 * t571 + t417 * t573) * qJD(4) + t596 * pkin(3) + (-t385 * t571 + t386 * t573 - t602) * qJ(4)) * MDP(15) + (-t392 * t528 - t612 * t679) * MDP(16) + (-t392 * t610 - t393 * t528 - t448 * t679 + t612 * t678) * MDP(17) + (t526 * t528 - t558 * t679 + t612 * t671) * MDP(18) + (t448 * t671 + t526 * t610 + t558 * t678) * MDP(19) + t558 * MDP(20) * t671 + (-t379 * t671 - t562 * t393 - t394 * t610 + t678 * t439 - t457 * t448 + t558 * t729 + t589) * MDP(21) + (t380 * t671 + t562 * t392 + t394 * t528 + t679 * t439 + t457 * t612 + t558 * t730 + t721) * MDP(22) + (-t369 * t610 + t371 * t671 + t390 * t678 + t393 * t446 + t448 * t685 + t558 * t683 + t589) * MDP(23) + (t367 * t610 + t368 * t528 + t371 * t679 - t375 * t678 + t392 * t611 - t393 * t471 - t448 * t684 - t612 * t683 - t602) * MDP(24) + (-t369 * t528 - t375 * t671 - t390 * t679 + t392 * t446 - t558 * t684 + t612 * t685 - t721) * MDP(25) + (-MDP(5) * t576 * t578 + MDP(6) * t675) * t581 + (t367 * t471 - t368 * t611 + t369 * t446 + t371 * t683 + t375 * t684 + t390 * t685 - t711 * t602 + t603 * (pkin(5) * t565 + qJ(6) * t564 + t562)) * MDP(26); (-t524 * t670 - t624) * MDP(12) + ((-t523 + t669) * t670 + t608) * MDP(13) + (-t523 ^ 2 - t524 ^ 2) * MDP(14) + (t416 * t524 - t417 * t523 - t583 + t722) * MDP(15) + (-t718 - t735) * MDP(24) + (t371 * t612 + t375 * t448 + t369 - t603) * MDP(26) - t723 * (t392 - t734) + t724 * (t393 + t733); -MDP(16) * t701 + (t718 - t735) * MDP(17) + t384 * MDP(18) + (-t393 + t733) * MDP(19) + t526 * MDP(20) + (t439 * t612 + t590 - t702) * MDP(21) + (-t379 * t558 + t439 * t448 - t586) * MDP(22) + (-t402 * t448 - t585 - t702 + 0.2e1 * t715) * MDP(23) + (pkin(5) * t392 - qJ(6) * t393 - (t375 - t380) * t612 + (t371 - t665) * t448) * MDP(24) + (0.2e1 * t705 - t390 * t448 - t402 * t612 + (-0.2e1 * qJD(6) + t379) * t558 + t586) * MDP(25) + (t367 * qJ(6) - t368 * pkin(5) - t390 * t402 - t371 * t380 - g(1) * (-pkin(5) * t413 + qJ(6) * t414) - g(2) * (-pkin(5) * t411 + qJ(6) * t412) - g(3) * (-pkin(5) * t453 + qJ(6) * t454) + t665 * t375) * MDP(26); (-t526 - t701) * MDP(23) + t384 * MDP(24) + (-t558 ^ 2 - t718) * MDP(25) + (t375 * t558 + t585 - t715) * MDP(26);];
tau  = t1;
