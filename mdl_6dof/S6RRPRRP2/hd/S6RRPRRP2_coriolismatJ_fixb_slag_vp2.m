% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:43:13
% EndTime: 2019-03-09 11:43:36
% DurationCPUTime: 15.07s
% Computational Cost: add. (32520->688), mult. (63128->902), div. (0->0), fcn. (74143->8), ass. (0->356)
t415 = sin(qJ(4));
t412 = sin(pkin(10));
t413 = cos(pkin(10));
t416 = sin(qJ(2));
t599 = -qJ(3) - pkin(7);
t486 = t599 * t416;
t418 = cos(qJ(2));
t672 = t599 * t418;
t344 = -t412 * t486 + t413 * t672;
t378 = -t412 * t416 + t413 * t418;
t445 = -t378 * pkin(8) + t344;
t617 = cos(qJ(4));
t379 = -t412 * t418 - t413 * t416;
t695 = t412 * t672 + t413 * t486;
t725 = t379 * pkin(8) + t695;
t174 = t415 * t445 + t617 * t725;
t748 = t174 * mrSges(5,2);
t613 = pkin(2) * t413;
t401 = pkin(3) + t613;
t614 = pkin(2) * t412;
t365 = t415 * t401 + t617 * t614;
t747 = t174 * t365;
t414 = sin(qJ(5));
t746 = t174 * t414;
t740 = t415 * t725 - t617 * t445;
t562 = t174 * t740;
t417 = cos(qJ(5));
t745 = t417 * t174;
t612 = pkin(4) * t740;
t364 = t617 * t401 - t415 * t614;
t362 = -pkin(4) - t364;
t744 = t362 * t740;
t337 = t617 * t378 + t379 * t415;
t708 = t417 * t337;
t525 = mrSges(7,2) * t708;
t455 = t415 * t378 - t617 * t379;
t679 = t455 * mrSges(7,1);
t730 = t525 - t679;
t643 = t730 / 0.2e1;
t720 = mrSges(6,1) * t455 - mrSges(6,3) * t708;
t741 = t643 - t720 / 0.2e1;
t743 = t741 * t417;
t465 = pkin(5) * t414 - qJ(6) * t417;
t101 = t455 * t465 - t174;
t580 = t417 * mrSges(6,2);
t582 = t414 * mrSges(6,1);
t390 = t580 + t582;
t223 = t390 * t455;
t557 = t337 * qJ(6);
t538 = t417 * t740;
t402 = -pkin(2) * t418 - pkin(1);
t358 = -t378 * pkin(3) + t402;
t606 = t455 * pkin(9);
t171 = -pkin(4) * t337 + t358 - t606;
t546 = t414 * t171;
t84 = t538 + t546;
t63 = t84 - t557;
t83 = t171 * t417 - t414 * t740;
t64 = pkin(5) * t337 - t83;
t710 = t390 * t337;
t579 = t417 * mrSges(7,3);
t581 = t414 * mrSges(7,1);
t389 = -t579 + t581;
t711 = t389 * t337;
t709 = t414 * t337;
t721 = -mrSges(6,2) * t455 - mrSges(6,3) * t709;
t722 = -mrSges(7,2) * t709 + mrSges(7,3) * t455;
t742 = t101 * t711 - t174 * t710 + t740 * t223 + t63 * t722 + t64 * t730 + t83 * t720 + t84 * t721;
t739 = -t720 + t730;
t738 = t379 * mrSges(4,1) - t378 * mrSges(4,2);
t619 = t417 / 0.2e1;
t621 = -t417 / 0.2e1;
t623 = t414 / 0.2e1;
t729 = t722 + t721;
t736 = t619 * t720 + t621 * t730 + t623 * t729;
t410 = t414 ^ 2;
t411 = t417 ^ 2;
t731 = t410 + t411;
t637 = -t337 / 0.2e1;
t693 = t455 / 0.2e1;
t694 = -t455 / 0.2e1;
t666 = Ifges(7,4) * t417 + Ifges(7,6) * t414;
t667 = Ifges(6,5) * t417 - Ifges(6,6) * t414;
t701 = t666 + t667;
t712 = Ifges(7,2) + Ifges(6,3);
t594 = Ifges(6,4) * t417;
t470 = -Ifges(6,2) * t414 + t594;
t698 = Ifges(6,6) * t455 + t470 * t337;
t592 = Ifges(7,5) * t417;
t467 = Ifges(7,3) * t414 + t592;
t699 = Ifges(7,6) * t455 + t467 * t337;
t723 = -t698 / 0.2e1 + t699 / 0.2e1;
t595 = Ifges(6,4) * t414;
t472 = Ifges(6,1) * t417 - t595;
t696 = Ifges(6,5) * t455 + t472 * t337;
t593 = Ifges(7,5) * t414;
t471 = Ifges(7,1) * t417 + t593;
t697 = Ifges(7,4) * t455 + t471 * t337;
t724 = t697 / 0.2e1 + t696 / 0.2e1;
t728 = t723 * t414 + t724 * t417 + t358 * mrSges(5,1) + 0.2e1 * Ifges(5,4) * t694 + t701 * t693 + (Ifges(5,2) + t712) * t637;
t726 = pkin(4) * t710;
t554 = t455 * t414;
t526 = mrSges(6,3) * t554;
t231 = mrSges(6,2) * t337 - t526;
t585 = t337 * mrSges(7,3);
t239 = -mrSges(7,2) * t554 - t585;
t671 = t231 + t239;
t609 = pkin(9) * t337;
t611 = pkin(4) * t455;
t243 = -t609 + t611;
t688 = mrSges(7,2) + mrSges(6,3);
t719 = t688 * t731;
t391 = -Ifges(7,3) * t417 + t593;
t396 = Ifges(6,2) * t417 + t595;
t491 = t396 / 0.4e1 - t391 / 0.4e1;
t718 = t696 / 0.4e1 + t697 / 0.4e1 - t491 * t337;
t397 = Ifges(7,1) * t414 - t592;
t398 = Ifges(6,1) * t414 + t594;
t490 = -t398 / 0.4e1 - t397 / 0.4e1;
t717 = t698 / 0.4e1 - t699 / 0.4e1 - t490 * t337;
t653 = m(7) / 0.2e1;
t715 = 0.2e1 * t653;
t714 = t337 / 0.2e1;
t713 = -t679 / 0.2e1;
t597 = t64 + t83;
t705 = t597 * t417;
t443 = t688 * (t410 / 0.2e1 + t411 / 0.2e1);
t668 = t398 + t397;
t451 = t471 * t455;
t138 = -Ifges(7,4) * t337 + t451;
t452 = t472 * t455;
t141 = -Ifges(6,5) * t337 + t452;
t522 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t479 = t522 * t337;
t704 = t479 - t138 / 0.2e1 - t141 / 0.2e1;
t703 = t731 * t364;
t702 = t731 * t337;
t520 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t521 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t700 = -(t520 * t414 - t522 * t417) * t337 + t521 * t455;
t691 = pkin(5) * t455;
t690 = t741 * pkin(9);
t689 = -mrSges(6,1) - mrSges(7,1);
t687 = Ifges(7,4) + Ifges(6,5);
t474 = t417 * mrSges(6,1) - t414 * mrSges(6,2);
t678 = t474 + mrSges(5,1);
t677 = qJ(6) * t455;
t549 = t414 * qJ(6);
t466 = t417 * pkin(5) + t549;
t384 = -pkin(4) - t466;
t343 = t384 - t364;
t676 = t343 * t455;
t675 = t384 * t455;
t392 = Ifges(6,5) * t414 + Ifges(6,6) * t417;
t394 = Ifges(7,4) * t414 - Ifges(7,6) * t417;
t663 = t394 / 0.4e1 + t392 / 0.4e1;
t674 = (-Ifges(5,6) / 0.2e1 + t663) * t455;
t473 = t417 * mrSges(7,1) + t414 * mrSges(7,3);
t665 = (-t473 - t678) * t365 + (-mrSges(5,2) + t719) * t364;
t494 = 0.2e1 * t714;
t662 = Ifges(5,1) * t694 - Ifges(5,4) * t337;
t403 = m(7) * qJ(6) + mrSges(7,3);
t553 = t455 * t417;
t236 = -mrSges(6,1) * t337 - mrSges(6,3) * t553;
t237 = mrSges(7,1) * t337 + mrSges(7,2) * t553;
t624 = -t414 / 0.2e1;
t661 = t236 * t621 + t237 * t619 + t624 * t671;
t660 = (t396 / 0.2e1 - t391 / 0.2e1) * t414 - (t398 / 0.2e1 + t397 / 0.2e1) * t417 - Ifges(5,5);
t217 = t466 * t455;
t222 = t389 * t455;
t659 = -t465 * t222 / 0.2e1 + t217 * t473 / 0.2e1 + t174 * t390 / 0.2e1 - t101 * t389 / 0.2e1;
t658 = 0.2e1 * pkin(9);
t657 = 2 * qJD(4);
t656 = m(5) / 0.2e1;
t655 = m(6) / 0.2e1;
t654 = -m(7) / 0.2e1;
t652 = -mrSges(5,1) / 0.2e1;
t651 = mrSges(6,1) / 0.2e1;
t650 = -mrSges(7,1) / 0.2e1;
t649 = -mrSges(6,2) / 0.2e1;
t409 = t416 * pkin(2);
t360 = -pkin(3) * t379 + t409;
t176 = t360 + t243;
t86 = t414 * t176 + t745;
t65 = t86 + t677;
t648 = -t65 / 0.2e1;
t97 = t414 * t243 + t745;
t72 = t97 + t677;
t647 = t72 / 0.2e1;
t96 = t243 * t417 - t746;
t73 = -t96 - t691;
t646 = t73 / 0.2e1;
t85 = t176 * t417 - t746;
t66 = -t85 - t691;
t645 = m(7) * t66;
t642 = -t722 / 0.2e1;
t641 = t722 / 0.2e1;
t633 = t343 / 0.2e1;
t632 = t362 / 0.2e1;
t631 = t364 / 0.2e1;
t630 = -t384 / 0.2e1;
t629 = t384 / 0.2e1;
t628 = -t473 / 0.2e1;
t627 = -t474 / 0.2e1;
t622 = t414 / 0.4e1;
t620 = -t417 / 0.4e1;
t618 = t417 / 0.4e1;
t616 = m(7) * t465;
t615 = m(7) * t414;
t610 = pkin(4) * t390;
t607 = pkin(9) * t417;
t603 = t72 * mrSges(7,2);
t602 = t73 * mrSges(7,2);
t601 = t96 * mrSges(6,3);
t600 = t97 * mrSges(6,3);
t598 = t63 - t84;
t590 = t740 * mrSges(5,1);
t100 = t337 * t465 + t740;
t333 = t337 * mrSges(5,2);
t480 = Ifges(5,2) / 0.2e1 + t521;
t539 = t417 * t141;
t540 = t417 * t138;
t450 = t470 * t455;
t135 = -Ifges(6,6) * t337 + t450;
t547 = t414 * t135;
t314 = Ifges(7,5) * t553;
t130 = -Ifges(7,6) * t337 + Ifges(7,3) * t554 + t314;
t548 = t414 * t130;
t3 = m(6) * (t83 * t85 + t84 * t86 - t562) + m(7) * (t100 * t101 + t63 * t65 + t64 * t66) + t65 * t239 + t86 * t231 + t85 * t236 + t66 * t237 + t100 * t222 + (-mrSges(4,2) * t409 - Ifges(4,4) * t379) * t379 + (-mrSges(4,1) * t409 + Ifges(4,4) * t378 + (-Ifges(4,1) + Ifges(4,2)) * t379) * t378 + (t360 * mrSges(5,2) + Ifges(5,1) * t714 + t728) * t455 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t418 + (Ifges(3,1) - Ifges(3,2)) * t416) * t418 + (-t480 * t455 - t360 * mrSges(5,1) + t540 / 0.2e1 + t539 / 0.2e1 + t548 / 0.2e1 - t547 / 0.2e1 + t701 * t637 - t662) * t337 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t416) * t416 + (m(4) * t409 - t738) * t402 + t742 + (m(5) * t360 + t333) * t358;
t587 = t3 * qJD(1);
t102 = pkin(5) * t709 - qJ(6) * t708 + t740;
t478 = t520 * t337;
t446 = -t130 / 0.2e1 + t135 / 0.2e1 - t478;
t4 = t102 * t222 + t97 * t231 + t96 * t236 + t73 * t237 + t72 * t239 + m(6) * (t83 * t96 + t84 * t97 - t562) + m(7) * (t101 * t102 + t63 * t72 + t64 * t73) + t728 * t455 - (-t358 * mrSges(5,2) + t704 * t417 + t446 * t414 + (-Ifges(5,1) / 0.2e1 + t480) * t455 + t662) * t337 + t742;
t583 = t4 * qJD(1);
t404 = t417 * mrSges(7,2);
t578 = t63 * t414;
t577 = t66 * t414;
t576 = t84 * t414;
t575 = t86 * t417;
t224 = t391 * t455;
t225 = t396 * t455;
t226 = -Ifges(7,1) * t554 + t314;
t227 = t398 * t455;
t453 = t473 * t455;
t454 = t474 * t455;
t9 = t174 * t454 - t84 * t237 - m(7) * (t101 * t217 + t64 * t84) - t217 * t222 - t101 * t453 + t84 * t236 + ((-t225 / 0.2e1 + t224 / 0.2e1 + t64 * mrSges(7,2) - t704) * t414 + (-t226 / 0.2e1 + t227 / 0.2e1 + t63 * mrSges(7,2) + t84 * mrSges(6,3) + t446) * t417) * t455 + (-m(7) * t63 - t526 - t671) * t83;
t574 = t9 * qJD(1);
t560 = t174 * t455;
t10 = (t378 ^ 2 + t379 ^ 2) * mrSges(4,3) + (mrSges(5,3) * t455 + t222 + t223) * t455 + (mrSges(5,3) * t337 + t671 * t417 + (-t236 + t237) * t414) * t337 + m(7) * (t101 * t455 + (t414 * t64 + t417 * t63) * t337) + m(6) * (-t560 + (-t414 * t83 + t417 * t84) * t337) + m(5) * (t337 * t740 - t560) + m(4) * (-t344 * t378 + t379 * t695);
t573 = qJD(1) * t10;
t28 = m(7) * (-t101 * t553 - t63 * t337) - t337 * t239 - t222 * t553;
t572 = qJD(1) * t28;
t571 = t100 * t384;
t570 = t100 * t473;
t569 = t101 * t465;
t567 = t101 * t414;
t566 = t102 * t473;
t439 = t443 * t337;
t528 = m(4) * pkin(2) / 0.2e1;
t363 = pkin(9) + t365;
t534 = t702 * t363;
t422 = t439 + (t337 * t365 - t364 * t455) * t656 + (t362 * t455 + t534) * t655 + (t534 + t676) * t653 + (t378 * t412 + t379 * t413) * t528;
t431 = t360 * t656 + (t414 * t86 + t417 * t85) * t655 + (t414 * t65 - t417 * t66) * t653 + t416 * t528;
t493 = t627 + t628;
t497 = t642 - t721 / 0.2e1;
t11 = -t333 + t743 + t497 * t414 + (-mrSges(5,1) + t493) * t455 + t422 - t431 + t738;
t565 = t11 * qJD(1);
t430 = (t597 * t414 + t598 * t417) * t654 + t236 * t623 + t237 * t624 + t671 * t621;
t435 = -t582 / 0.2e1 - t581 / 0.2e1 + t579 / 0.2e1 - t580 / 0.2e1 + t465 * t654;
t433 = t435 * t337;
t13 = t433 + t430;
t564 = t13 * qJD(1);
t532 = t702 * pkin(9);
t423 = t439 + (t652 + t493) * t455 - t333 / 0.2e1 + (t532 - t611) * t655 + (t532 + t675) * t653;
t434 = -m(6) * (t414 * t97 + t417 * t96) / 0.2e1 + (t414 * t72 - t417 * t73) * t654 + mrSges(5,2) * t637 + mrSges(5,1) * t694;
t495 = t641 + t721 / 0.2e1;
t15 = -t495 * t414 + t423 + t434 + t743;
t563 = t15 * qJD(1);
t561 = t740 * t474;
t545 = t414 * t222;
t542 = t414 * t391;
t541 = t414 * t396;
t535 = t417 * t473;
t533 = t703 * t363;
t531 = t703 * pkin(9);
t212 = t494 * t615;
t529 = t212 * qJD(1);
t527 = m(7) * t646;
t519 = mrSges(6,3) * t624;
t518 = t404 / 0.2e1;
t517 = t576 / 0.2e1;
t516 = -t545 / 0.2e1 - t535 * t694 - t337 * t518;
t514 = t545 / 0.2e1;
t511 = t709 / 0.2e1;
t508 = t708 / 0.2e1;
t498 = t237 / 0.2e1 - t236 / 0.2e1;
t496 = t239 / 0.2e1 + t231 / 0.2e1;
t482 = t394 / 0.2e1 + t392 / 0.2e1 - Ifges(5,6);
t476 = (-t343 - t384) * t653;
t464 = t65 * t417 + t577;
t463 = -t85 * t414 + t575;
t424 = t496 * t364 + t495 * t363 + (t363 * t97 + t364 * t84) * t655 + (t363 * t72 + t364 * t63) * t653 + t717;
t425 = t741 * t363 + t498 * t364 + (-t363 * t96 - t364 * t83) * t655 + (t363 * t73 + t364 * t64) * t653 + t718;
t426 = (t223 / 0.2e1 + t222 / 0.2e1) * t365 + t674 + (t744 - t747) * t655 + (t101 * t365 + t102 * t343) * t653 + t711 * t633 + t710 * t632;
t456 = t726 / 0.2e1 + t711 * t630;
t1 = (t497 * pkin(9) + (t97 / 0.2e1 - t86 / 0.2e1) * mrSges(6,3) + (t647 + t648) * mrSges(7,2) + (-m(6) * t86 / 0.4e1 - m(7) * t65 / 0.4e1) * t658 + t424 - t717) * t417 + (-t690 + (-t96 / 0.2e1 + t85 / 0.2e1) * mrSges(6,3) + (t646 - t66 / 0.2e1) * mrSges(7,2) + (m(6) * t85 / 0.4e1 - t645 / 0.4e1) * t658 + t425 - t718) * t414 + t426 - (t102 / 0.2e1 - t100 / 0.2e1) * t473 - t674 + (t714 + t637) * Ifges(5,5) + t612 * t655 + t571 * t654 + t456;
t30 = m(7) * (t343 * t365 + t533) + m(6) * (t362 * t365 + t533) + t665;
t462 = t1 * qJD(1) + t30 * qJD(2);
t427 = -t659 + t539 / 0.4e1 + t540 / 0.4e1 - t547 / 0.4e1 + t548 / 0.4e1 - t224 * t620 + mrSges(6,3) * t517 + t84 * t519 - t414 * t450 / 0.4e1 - t701 * t337 / 0.4e1 + (-t227 + t226) * t622 - t491 * t553 + t597 * t518 + (-t578 / 0.2e1 + t517) * mrSges(7,2) + (t452 + t451 - t225) * t618 + (t467 / 0.4e1 - t668 / 0.4e1) * t554;
t419 = t427 + ((t576 - t578 + t705) * t653 - t443 * t455 + t661) * t363 + (t473 * t633 + t474 * t632) * t455 + (t217 * t343 + t569) * t653;
t429 = (-pkin(5) * t66 + qJ(6) * t65) * t654 + pkin(5) * t643 + qJ(6) * t642 + mrSges(7,3) * t648 + t66 * mrSges(7,1) / 0.2e1 - t85 * mrSges(6,1) / 0.2e1 + t86 * mrSges(6,2) / 0.2e1;
t5 = t419 + t429 - t700;
t332 = t343 * t389;
t340 = t362 * t390;
t421 = t465 * t473 - t542 / 0.2e1 + t541 / 0.2e1 + t467 * t619 + (t472 + t471) * t624 + (t470 + t668) * t621;
t54 = -t343 * t616 - t332 - t340 + t421;
t461 = t5 * qJD(1) - t54 * qJD(2);
t460 = t217 * t384 + t569;
t440 = m(7) * (-t567 + (-t337 * t363 - t676) * t417);
t457 = t645 / 0.2e1 + t713;
t23 = t514 + (t494 * mrSges(7,2) + t455 * t628) * t417 - t440 / 0.2e1 + t457;
t282 = (m(7) * t343 - t473) * t414;
t459 = qJD(1) * t23 + qJD(2) * t282;
t34 = -t585 + 0.2e1 * (-t557 / 0.2e1 + t546 / 0.4e1 + t538 / 0.4e1 - t84 / 0.4e1) * m(7);
t458 = qJD(1) * t34 + qJD(5) * t403;
t359 = t384 * t389;
t420 = t332 / 0.2e1 + t340 / 0.2e1 + t359 / 0.2e1 - t610 / 0.2e1 - t465 * t476 - t421;
t432 = t435 * t364;
t38 = t432 - t420;
t428 = (-pkin(5) * t73 + qJ(6) * t72) * t653 - pkin(5) * t730 / 0.2e1 + qJ(6) * t641 + mrSges(7,3) * t647 + t73 * t650 + t96 * t651 + t97 * t649;
t7 = -(-t667 / 0.4e1 - t666 / 0.4e1) * t337 + t460 * t654 + (t225 / 0.4e1 - t224 / 0.4e1 - t141 / 0.4e1 - t138 / 0.4e1 + t479 + (-t83 / 0.2e1 - t64 / 0.2e1) * mrSges(7,2) + (t597 * t654 - t498) * pkin(9)) * t417 + (t227 / 0.4e1 - t226 / 0.4e1 + t135 / 0.4e1 - t130 / 0.4e1 - t478 + (-t84 / 0.2e1 + t63 / 0.2e1) * mrSges(7,2) + (t598 * t653 + t496) * pkin(9)) * t414 + ((pkin(4) * t651 + mrSges(7,1) * t630 + (-Ifges(6,1) / 0.4e1 - Ifges(7,1) / 0.4e1) * t417 + t491) * t417 + (pkin(4) * t649 + mrSges(7,3) * t630 + (-Ifges(6,2) / 0.4e1 - Ifges(7,3) / 0.4e1) * t414 + (-Ifges(7,5) / 0.2e1 + Ifges(6,4) / 0.2e1) * t417 - t490) * t414 + t443 * pkin(9) + t521) * t455 + t428 + t659;
t87 = -t384 * t616 - t359 + t421 + t610;
t448 = t7 * qJD(1) + t38 * qJD(2) + t87 * qJD(4);
t181 = (-t473 + (t631 + t633 + t629) * m(7)) * t414;
t441 = m(7) * (-t567 + (-t609 - t675) * t417);
t26 = t525 + t514 + (t650 - t535 / 0.2e1) * t455 + t527 - t441 / 0.2e1;
t341 = (m(7) * t384 - t473) * t414;
t447 = qJD(1) * t26 + qJD(2) * t181 + qJD(4) * t341;
t437 = (-mrSges(7,2) * t549 - pkin(5) * t404 + t701) * qJD(5);
t436 = qJD(5) * (-m(7) * t466 - t473 - t474);
t385 = m(7) * t607 + t404;
t342 = m(7) * t363 * t417 + t404;
t211 = m(7) * t511 - t653 * t709;
t182 = (m(7) * t631 + t473 + t476) * t414;
t39 = t432 + t420;
t31 = t63 * t715 + t239;
t27 = t441 / 0.2e1 + t525 / 0.2e1 + t713 + t527 + t516;
t24 = t440 / 0.2e1 + mrSges(7,2) * t508 + t457 + t516;
t16 = t423 - t434 + t736;
t14 = t433 - t430;
t12 = t493 * t455 + t422 + t431 + t736;
t8 = t427 + t428 + t460 * t653 + t453 * t629 - pkin(4) * t454 / 0.2e1 + ((-t598 * t414 + t705) * t653 + t661) * pkin(9) + t700 - t606 * t719 / 0.2e1;
t6 = t419 - Ifges(6,6) * t709 / 0.2e1 + Ifges(7,6) * t511 - t429 + t712 * t693 + t687 * t508;
t2 = (-t541 / 0.4e1 + t542 / 0.4e1) * t337 + t729 * t607 / 0.2e1 - t456 - t748 + t426 - t590 / 0.2e1 + (t627 + t652) * t740 + (t600 / 0.2e1 + t603 / 0.2e1 + t424) * t417 + mrSges(7,2) * t577 / 0.2e1 + mrSges(6,3) * t575 / 0.2e1 + t494 * Ifges(5,5) + t663 * t455 + Ifges(5,6) * t694 + (-t601 / 0.2e1 + t602 / 0.2e1 + t425 + t690) * t414 + (pkin(9) * t464 + t571) * t653 + (pkin(9) * t463 - t612) * t655 + t65 * t518 + t85 * t519 + t699 * t620 + t698 * t618 + (t696 + t697) * t622 + t668 * t708 / 0.4e1 - t561 / 0.2e1 - t566 / 0.2e1 - t570 / 0.2e1;
t17 = [qJD(2) * t3 + qJD(3) * t10 + qJD(4) * t4 - qJD(5) * t9 + qJD(6) * t28, t12 * qJD(3) + t2 * qJD(4) + t6 * qJD(5) + t24 * qJD(6) + t587 + (t344 * mrSges(4,1) - t695 * mrSges(4,2) + Ifges(3,5) * t418 + Ifges(4,5) * t378 - Ifges(3,6) * t416 + Ifges(4,6) * t379 + t343 * t711 + t362 * t710 - t561 - t570 - t748 - t590 + (-t365 * mrSges(5,3) + t482) * t455 + (t65 * mrSges(7,2) + t86 * mrSges(6,3) + t363 * t729 - t723) * t417 + (t66 * mrSges(7,2) - t85 * mrSges(6,3) + t363 * t739 + t724) * t414 + (t100 * t343 + t363 * t464) * t715 + 0.2e1 * (t363 * t463 + t744) * t655 + 0.2e1 * (-t364 * t740 + t747) * t656 + 0.2e1 * (t344 * t413 + t412 * t695) * t528 + (-t364 * mrSges(5,3) - t660) * t337 + (-mrSges(3,1) * t418 + mrSges(3,2) * t416) * pkin(7) + (-t378 * t613 + t379 * t614) * mrSges(4,3)) * qJD(2), qJD(2) * t12 + qJD(4) * t16 + qJD(5) * t14 + qJD(6) * t211 + t573, t583 + t2 * qJD(2) + t16 * qJD(3) + t8 * qJD(5) + t27 * qJD(6) + ((t102 * t384 + (t73 * t414 + t72 * t417) * pkin(9)) * t653 + (-t612 + (-t96 * t414 + t97 * t417) * pkin(9)) * t655) * t657 + (-t726 + t384 * t711 - t566 - t748 + (pkin(9) * t729 + t600 + t603 - t723) * t417 + (pkin(9) * t739 - t601 + t602 + t724) * t414 + t482 * t455 - t660 * t337 - t678 * t740) * qJD(4), t6 * qJD(2) + t14 * qJD(3) + t8 * qJD(4) + t31 * qJD(6) - t574 + ((-m(7) * pkin(5) + t689) * t84 + (-mrSges(6,2) + t403) * t83 + ((-mrSges(7,2) * qJ(6) - Ifges(6,6) + Ifges(7,6)) * t417 + (mrSges(7,2) * pkin(5) - t687) * t414) * t455) * qJD(5), qJD(2) * t24 + qJD(3) * t211 + qJD(4) * t27 + qJD(5) * t31 + t572; qJD(3) * t11 + qJD(4) * t1 + qJD(5) * t5 - qJD(6) * t23 - t587, qJD(4) * t30 - qJD(5) * t54 - qJD(6) * t282, t565, t39 * qJD(5) + t182 * qJD(6) + ((t365 * t384 + t531) * t653 + (-pkin(4) * t365 + t531) * t655) * t657 + t462 + t665 * qJD(4), t39 * qJD(4) + t342 * qJD(6) + t363 * t436 + t437 + t461, qJD(4) * t182 + qJD(5) * t342 - t459; -qJD(2) * t11 - qJD(4) * t15 - qJD(5) * t13 - qJD(6) * t212 - t573, -t565, 0, -t563, -t564 + (t579 - t580 - t616) * qJD(5) + (m(7) * qJD(6) + t689 * qJD(5)) * t414, qJD(5) * t615 - t529; -qJD(2) * t1 + qJD(3) * t15 - qJD(5) * t7 - qJD(6) * t26 - t583, -qJD(5) * t38 - qJD(6) * t181 - t462, t563, -qJD(5) * t87 - qJD(6) * t341, pkin(9) * t436 + t385 * qJD(6) + t437 - t448, qJD(5) * t385 - t447; -qJD(2) * t5 + qJD(3) * t13 + qJD(4) * t7 + qJD(6) * t34 + t574, qJD(4) * t38 - t461, t564, t448, t403 * qJD(6), t458; qJD(2) * t23 + qJD(3) * t212 + qJD(4) * t26 - qJD(5) * t34 - t572, qJD(4) * t181 + t459, t529, t447, -t458, 0;];
Cq  = t17;
