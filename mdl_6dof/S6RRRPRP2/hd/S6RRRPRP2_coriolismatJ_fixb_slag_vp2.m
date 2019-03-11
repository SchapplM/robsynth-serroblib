% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRRPRP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:35:14
% EndTime: 2019-03-09 16:35:38
% DurationCPUTime: 13.58s
% Computational Cost: add. (32716->676), mult. (63471->877), div. (0->0), fcn. (74113->8), ass. (0->376)
t441 = sin(qJ(3));
t442 = sin(qJ(2));
t444 = cos(qJ(3));
t445 = cos(qJ(2));
t401 = -t441 * t442 + t444 * t445;
t402 = -t441 * t445 - t444 * t442;
t439 = sin(pkin(10));
t619 = cos(pkin(10));
t346 = t619 * t401 + t402 * t439;
t440 = sin(qJ(5));
t443 = cos(qJ(5));
t653 = Ifges(7,5) * t443;
t411 = Ifges(7,3) * t440 + t653;
t471 = t439 * t401 - t402 * t619;
t155 = Ifges(7,6) * t471 + t346 * t411;
t435 = Ifges(6,4) * t443;
t718 = -Ifges(6,2) * t440 + t435;
t157 = Ifges(6,6) * t471 + t346 * t718;
t672 = t443 / 0.2e1;
t673 = -t443 / 0.2e1;
t675 = t440 / 0.2e1;
t691 = t471 / 0.2e1;
t654 = Ifges(6,4) * t440;
t419 = Ifges(6,1) * t443 - t654;
t432 = Ifges(7,5) * t440;
t721 = Ifges(7,1) * t443 + t432;
t742 = t721 + t419;
t756 = Ifges(6,5) + Ifges(7,4);
t727 = t742 * t346 + t471 * t756;
t414 = Ifges(6,2) * t443 + t654;
t676 = -t440 / 0.2e1;
t714 = Ifges(7,3) * t443 - t432;
t416 = Ifges(7,1) * t440 - t653;
t418 = Ifges(6,1) * t440 + t435;
t771 = t416 + t418;
t736 = t414 * t676 + t672 * t771 - t675 * t714;
t754 = Ifges(6,6) - Ifges(7,6);
t459 = Ifges(4,5) * t401 + Ifges(4,6) * t402 - Ifges(5,6) * t471 + t155 * t673 + t157 * t672 + (t756 * t440 + t754 * t443) * t691 + t727 * t675 + (Ifges(5,5) + t736) * t346;
t695 = -pkin(8) - pkin(7);
t421 = t695 * t442;
t422 = t695 * t445;
t722 = t444 * t421 + t441 * t422;
t750 = t722 * mrSges(4,2);
t360 = t421 * t441 - t422 * t444;
t751 = t360 * mrSges(4,1);
t494 = t443 * mrSges(6,1) - t440 * mrSges(6,2);
t311 = qJ(4) * t401 + t360;
t741 = t402 * qJ(4) + t722;
t763 = t619 * t311 + t439 * t741;
t773 = t763 * t494;
t779 = t763 * mrSges(5,1);
t762 = -t439 * t311 + t619 * t741;
t780 = t762 * mrSges(5,2);
t493 = t443 * mrSges(7,1) + t440 * mrSges(7,3);
t594 = t346 * t440;
t332 = pkin(5) * t594;
t593 = t346 * t443;
t496 = qJ(6) * t593 - t332;
t769 = t763 - t496;
t784 = t769 * t493;
t790 = t459 - t750 - t751 - t773 - t779 - t780 - t784;
t789 = -t750 / 0.2e1 - t751 / 0.2e1 - t773 / 0.2e1 - t779 / 0.2e1 - t780 / 0.2e1 - t784 / 0.2e1;
t533 = -pkin(2) * t445 - pkin(1);
t369 = -t401 * pkin(3) + t533;
t788 = m(5) * t369;
t487 = t440 * pkin(5) - qJ(6) * t443;
t118 = t471 * t487 - t762;
t787 = t118 * t769;
t668 = pkin(2) * t441;
t423 = t439 * t668;
t666 = pkin(2) * t444;
t428 = pkin(3) + t666;
t376 = t428 * t619 - t423;
t373 = -pkin(4) - t376;
t574 = t440 * qJ(6);
t488 = t443 * pkin(5) + t574;
t355 = t373 - t488;
t786 = t355 * t769;
t532 = t619 * pkin(3);
t425 = -t532 - pkin(4);
t384 = -t488 + t425;
t785 = t384 * t769;
t200 = -pkin(4) * t346 - pkin(9) * t471 + t369;
t100 = t200 * t443 - t440 * t763;
t559 = t443 * t763;
t568 = t440 * t200;
t101 = t559 + t568;
t627 = t443 * mrSges(7,3);
t629 = t440 * mrSges(7,1);
t408 = -t627 + t629;
t244 = t408 * t346;
t628 = t443 * mrSges(6,2);
t630 = t440 * mrSges(6,1);
t409 = t628 + t630;
t245 = t409 * t346;
t246 = t408 * t471;
t252 = -mrSges(6,2) * t471 - mrSges(6,3) * t594;
t254 = mrSges(6,1) * t471 - mrSges(6,3) * t593;
t642 = t471 * mrSges(7,1);
t255 = mrSges(7,2) * t593 - t642;
t258 = -mrSges(7,2) * t594 + mrSges(7,3) * t471;
t340 = t346 * mrSges(5,2);
t504 = mrSges(5,1) * t471 + t340;
t631 = t402 * mrSges(4,3);
t729 = t471 * mrSges(5,3);
t595 = t346 * qJ(6);
t78 = t101 - t595;
t79 = pkin(5) * t346 - t100;
t783 = t769 * t246 + t100 * t254 + t101 * t252 + t118 * t244 - t762 * t245 - t763 * t729 + t79 * t255 + t78 * t258 + t360 * t631 + t369 * t504 + t533 * (-mrSges(4,1) * t402 + mrSges(4,2) * t401);
t778 = t373 * t763;
t777 = t425 * t763;
t776 = t440 * t762;
t775 = t443 * t762;
t774 = t762 * t763;
t247 = t409 * t471;
t772 = t247 + t729;
t437 = t440 ^ 2;
t438 = t443 ^ 2;
t770 = t437 + t438;
t505 = t619 * t441;
t377 = pkin(2) * t505 + t439 * t428;
t768 = -t376 * t763 + t377 * t762;
t767 = t439 * t762 - t619 * t763;
t732 = mrSges(7,2) + mrSges(6,3);
t726 = t252 + t258;
t592 = t471 * t440;
t547 = mrSges(6,3) * t592;
t253 = mrSges(6,2) * t346 - t547;
t639 = t346 * mrSges(7,3);
t259 = -mrSges(7,2) * t592 - t639;
t725 = t253 + t259;
t764 = t732 * t770;
t701 = m(7) / 0.2e1;
t760 = 0.2e1 * t701;
t702 = -m(7) / 0.2e1;
t759 = t255 / 0.2e1;
t758 = -t258 / 0.2e1;
t689 = -t346 / 0.2e1;
t690 = t346 / 0.2e1;
t757 = t493 / 0.2e1;
t755 = Ifges(4,2) - Ifges(4,1);
t753 = Ifges(5,4) * t346;
t730 = Ifges(5,4) * t471;
t752 = t346 * mrSges(5,3);
t671 = m(7) * qJ(6);
t429 = mrSges(7,3) + t671;
t621 = t100 + t79;
t719 = Ifges(7,4) * t443 + Ifges(7,6) * t440;
t720 = Ifges(6,5) * t443 - Ifges(6,6) * t440;
t717 = t719 + t720;
t749 = t346 * t717;
t748 = t443 * t621;
t723 = t384 + t355;
t498 = t723 * t702;
t747 = t487 * t498;
t464 = t732 * (t438 / 0.2e1 + t437 / 0.2e1);
t665 = pkin(3) * t402;
t213 = pkin(4) * t471 - pkin(9) * t346 - t665;
t106 = t213 * t443 - t776;
t107 = t440 * t213 + t775;
t502 = t254 * t672 + t255 * t673 + t675 * t726;
t703 = m(6) / 0.2e1;
t706 = -m(5) / 0.2e1;
t339 = t471 * qJ(6);
t85 = t339 + t107;
t663 = pkin(5) * t471;
t86 = -t106 - t663;
t746 = (t106 * t443 + t107 * t440) * t703 + (t440 * t85 - t443 * t86) * t701 + t665 * t706 + t502;
t745 = -t493 - t494;
t662 = t439 * pkin(3);
t424 = pkin(9) + t662;
t744 = t770 * t424;
t383 = t619 * t666 - t423;
t743 = t770 * t383;
t740 = (-mrSges(4,1) * t441 - mrSges(4,2) * t444) * pkin(2);
t739 = t360 * mrSges(4,3) + Ifges(4,4) * t402;
t738 = Ifges(5,1) * t471 + t753;
t734 = -t471 / 0.2e1;
t733 = -mrSges(6,1) - mrSges(7,1);
t731 = Ifges(7,2) + Ifges(6,3);
t728 = t355 * t471;
t509 = -t493 / 0.2e1 - t494 / 0.2e1;
t712 = t246 + t772;
t711 = pkin(5) * t759 + qJ(6) * t758;
t710 = m(7) * t384 - t493;
t582 = t424 * t443;
t709 = t726 * t582;
t591 = t471 * t443;
t256 = -mrSges(6,1) * t346 - mrSges(6,3) * t591;
t257 = mrSges(7,1) * t346 + mrSges(7,2) * t591;
t708 = t256 * t673 + t257 * t672 + t676 * t725;
t243 = t488 * t471;
t707 = -t487 * t246 / 0.2e1 + t243 * t757 + t762 * t409 / 0.2e1 - t118 * t408 / 0.2e1;
t705 = m(5) / 0.2e1;
t704 = -m(6) / 0.2e1;
t700 = m(5) * pkin(3);
t699 = -mrSges(6,1) / 0.2e1;
t698 = -mrSges(7,1) / 0.2e1;
t697 = -mrSges(6,2) / 0.2e1;
t696 = mrSges(7,3) / 0.2e1;
t667 = pkin(2) * t442;
t201 = t213 + t667;
t102 = t201 * t443 - t776;
t82 = -t102 - t663;
t694 = m(7) * t82;
t693 = m(7) * t86;
t684 = t383 / 0.2e1;
t683 = -t384 / 0.2e1;
t678 = -t425 / 0.2e1;
t677 = t425 / 0.2e1;
t670 = m(7) * t487;
t669 = m(7) * t440;
t659 = mrSges(4,3) * t401;
t651 = Ifges(7,2) * t471;
t649 = Ifges(6,3) * t471;
t103 = t440 * t201 + t775;
t388 = Ifges(4,4) * t401;
t495 = -mrSges(5,1) * t346 + mrSges(5,2) * t471;
t497 = -t665 + t667;
t162 = -Ifges(6,5) * t346 + t419 * t471;
t560 = t443 * t162;
t160 = -Ifges(7,4) * t346 + t471 * t721;
t561 = t443 * t160;
t158 = -Ifges(6,6) * t346 + t471 * t718;
t569 = t440 * t158;
t570 = t440 * t157;
t331 = Ifges(7,5) * t591;
t156 = -Ifges(7,6) * t346 + Ifges(7,3) * t592 + t331;
t571 = t440 * t156;
t572 = t440 * t155;
t81 = t339 + t103;
t2 = -m(6) * (t100 * t102 + t101 * t103 - t774) - m(7) * (t78 * t81 + t79 * t82 + t787) - t81 * t259 - t763 * t247 - t103 * t253 - t102 * t256 - t82 * t257 + (t442 ^ 2 - t445 ^ 2) * Ifges(3,4) - t722 * t659 - t762 * t752 + (mrSges(4,3) * t722 - t402 * t755 - t388) * t401 + (pkin(1) * mrSges(3,1) + (-Ifges(3,1) + Ifges(3,2)) * t445) * t442 + t739 * t402 + pkin(1) * mrSges(3,2) * t445 + (-t560 / 0.2e1 - t561 / 0.2e1 + t569 / 0.2e1 - t571 / 0.2e1 + t762 * mrSges(5,3) + t717 * t690 - t738) * t346 + (-(-Ifges(5,2) - t731) * t346 + t570 / 0.2e1 - t572 / 0.2e1 - t763 * mrSges(5,3) + t717 * t734 + t727 * t673 + t730) * t471 + (-m(4) * t533 + mrSges(4,1) * t401 + mrSges(4,2) * t402) * t667 - t783 + (-t495 - t788) * t497;
t647 = t2 * qJD(1);
t528 = t593 / 0.2e1;
t529 = t594 / 0.2e1;
t530 = -t594 / 0.2e1;
t4 = t156 * t529 + t401 * t388 + t86 * t257 + t85 * t259 + t107 * t253 + t106 * t256 + m(7) * (t78 * t85 + t79 * t86 + t787) + t158 * t530 + m(6) * (t100 * t106 + t101 * t107 - t774) - t665 * t788 + (t649 + t651 + t749) * t689 + (t162 + t160) * t528 + t727 * t591 / 0.2e1 + (-pkin(3) * t495 + t401 * t755 - t739) * t402 + (t471 * t717 + t572 - t730) * t691 + (t570 + t730) * t734 + (-Ifges(5,2) * t471 + t738 + t753) * t690 + t772 * t763 + (Ifges(5,2) * t734 + (Ifges(5,1) - t731) * t691) * t346 + t783;
t632 = t4 * qJD(1);
t430 = t443 * mrSges(7,2);
t626 = t81 * t443;
t625 = t82 * t440;
t624 = t85 * t443;
t623 = t86 * t440;
t248 = t714 * t471;
t249 = t414 * t471;
t250 = -Ifges(7,1) * t592 + t331;
t251 = t418 * t471;
t472 = t493 * t471;
t473 = t494 * t471;
t540 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t541 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t9 = -m(7) * (t101 * t79 + t118 * t243) - t243 * t246 - t118 * t472 + t762 * t473 + t101 * t256 - t101 * t257 + ((t162 / 0.2e1 - t249 / 0.2e1 + t160 / 0.2e1 - t248 / 0.2e1 + t79 * mrSges(7,2) - t541 * t346) * t440 + (t251 / 0.2e1 + t158 / 0.2e1 - t250 / 0.2e1 - t156 / 0.2e1 + t78 * mrSges(7,2) + t101 * mrSges(6,3) - t540 * t346) * t443) * t471 + (-m(7) * t78 - t547 - t725) * t100;
t622 = t9 * qJD(1);
t620 = t101 - t78;
t480 = -t100 * t440 + t101 * t443;
t486 = t440 * t79 + t443 * t78;
t599 = t762 * t471;
t10 = t712 * t471 + (t752 + t725 * t443 + (-t256 + t257) * t440) * t346 + m(7) * (t118 * t471 + t346 * t486) + m(6) * (t346 * t480 - t599) + m(5) * (t346 * t763 - t599);
t617 = qJD(1) * t10;
t30 = m(7) * (-t118 * t591 - t346 * t78) - t246 * t591 - t346 * t259;
t616 = qJD(1) * t30;
t615 = t101 * t440;
t614 = t102 * t440;
t613 = t103 * t443;
t612 = t106 * t440;
t611 = t107 * t443;
t374 = pkin(9) + t377;
t577 = t438 * t346;
t579 = t437 * t346;
t556 = (t577 + t579) * t374;
t450 = t464 * t346 + (t346 * t377 - t376 * t471) * t705 + (t373 * t471 + t556) * t703 + (t556 + t728) * t701;
t453 = t497 * t706 + (t102 * t443 + t103 * t440) * t704 + (t440 * t81 - t443 * t82) * t702;
t11 = -t340 + (t759 - t254 / 0.2e1) * t443 + (t758 - t252 / 0.2e1) * t440 + (-mrSges(5,1) + t509) * t471 + t450 + t453;
t610 = t11 * qJD(1);
t607 = t118 * t487;
t550 = t700 / 0.2e1;
t449 = (t384 * t701 + t425 * t703 - t550 * t619 + t509) * t471 + t732 * (t579 / 0.2e1 + t577 / 0.2e1) + (t439 * t550 + (t703 + t701) * t744) * t346;
t13 = t449 - t504 - t746;
t605 = t13 * qJD(1);
t470 = t621 * t702 - t257 / 0.2e1 + t256 / 0.2e1;
t499 = (t699 + t698) * t346;
t512 = t253 / 0.2e1 + t259 / 0.2e1;
t14 = t332 * t702 + (t620 * t701 + (t697 + t696 + t671 / 0.2e1) * t346 - t512) * t443 + (t499 + t470) * t440;
t604 = t14 * qJD(1);
t382 = (t439 * t444 + t505) * pkin(2);
t598 = t762 * t382;
t590 = t355 * t244;
t589 = t373 * t245;
t588 = t373 * t409;
t587 = t374 * t443;
t584 = t384 * t244;
t581 = t425 * t245;
t580 = t425 * t409;
t573 = t440 * t118;
t566 = t440 * t254;
t565 = t440 * t255;
t564 = t440 * t256;
t563 = t440 * t257;
t554 = t743 * t374;
t511 = 0.2e1 * t689;
t240 = t511 * t669;
t551 = t240 * qJD(1);
t549 = mrSges(7,2) * t625;
t548 = mrSges(7,2) * t624;
t546 = mrSges(6,3) * t613;
t545 = mrSges(6,3) * t612;
t544 = mrSges(6,3) * t611;
t543 = t376 * t752;
t542 = t377 * t729;
t539 = t424 * t566;
t538 = t424 * t565;
t537 = mrSges(7,2) * t676;
t536 = mrSges(6,3) * t676;
t535 = -t628 / 0.2e1;
t534 = t430 / 0.2e1;
t527 = t493 * t734;
t517 = t246 * t676;
t510 = -t355 / 0.2e1 + t683;
t506 = t408 + t670;
t503 = t729 * t662;
t485 = t625 + t626;
t484 = t623 + t624;
t483 = t532 * t752;
t478 = t611 - t612;
t446 = (t480 * t704 + t763 * t706 + t486 * t702 - t563 / 0.2e1 + t725 * t673) * t383 + t545 / 0.2e1 + t542 / 0.2e1 + t543 / 0.2e1 + (t118 * t382 + t786) * t702 - t544 / 0.2e1 - t726 * t587 / 0.2e1 - t712 * t382 / 0.2e1 + (t478 * t704 + t484 * t702 - t565 / 0.2e1 + t566 / 0.2e1) * t374 + (-t598 + t778) * t704 + (-t598 + t768) * t706 - t548 / 0.2e1 + (-t752 + t564) * t684 + t86 * t537 - t590 / 0.2e1 - t589 / 0.2e1 - t789;
t479 = t613 - t614;
t447 = (t424 * t479 + t777) * t703 + (t424 * t485 + t785) * t701 + t584 / 0.2e1 + t581 / 0.2e1 + t767 * t550 + t102 * t536 + t546 / 0.2e1 - t539 / 0.2e1 + t538 / 0.2e1 + t81 * t534 + t549 / 0.2e1 - t503 / 0.2e1 - t483 / 0.2e1 + t709 / 0.2e1 + t789;
t3 = t447 + t446;
t34 = t740 + (-mrSges(5,1) + t745) * t382 + (-mrSges(5,2) + t764) * t383 + m(7) * (t355 * t382 + t554) + m(6) * (t373 * t382 + t554) + m(5) * (-t376 * t382 + t377 * t383);
t482 = -t3 * qJD(1) + t34 * qJD(2);
t451 = -t707 + t571 / 0.4e1 - t569 / 0.4e1 + t561 / 0.4e1 + t101 * t536 + t78 * t537 + t560 / 0.4e1 - t414 * t591 / 0.4e1 + t411 * t592 / 0.4e1 - (t249 + t248) * t443 / 0.4e1 - t749 / 0.4e1 + (-t251 + t250) * t440 / 0.4e1 + t621 * t534 + t732 * t615 / 0.2e1 - (t718 + t771) * t592 / 0.4e1 + (-t714 + t742) * t591 / 0.4e1;
t448 = t451 + (t243 * t355 + t607) * t701 + (t373 * t494 / 0.2e1 + t355 * t757) * t471 + ((-t440 * t78 + t615 + t748) * t701 - t464 * t471 + t708) * t374;
t455 = (-pkin(5) * t82 + qJ(6) * t81) * t702 + t102 * t699 + t103 * mrSges(6,2) / 0.2e1 - t81 * mrSges(7,3) / 0.2e1 + t82 * mrSges(7,1) / 0.2e1;
t6 = t448 + (t440 * t540 - t443 * t541) * t346 + (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t471 + t455 + t711;
t370 = t487 * t493;
t466 = t411 * t673 + t718 * t672 + t675 * t742 - t370 + t736;
t66 = t355 * t506 + t466 + t588;
t481 = t6 * qJD(1) + t66 * qJD(2);
t477 = t243 * t384 + t607;
t458 = (mrSges(7,2) * t511 - t527) * t443 + t642 / 0.2e1 + t517;
t463 = (-t573 + (-t346 * t374 - t728) * t443) * t701;
t26 = t463 - t694 / 0.2e1 + t458;
t296 = (m(7) * t355 - t493) * t440;
t476 = -qJD(1) * t26 + qJD(2) * t296;
t36 = -t639 + 0.2e1 * (-t595 / 0.2e1 + t568 / 0.4e1 + t559 / 0.4e1 - t101 / 0.4e1) * m(7);
t475 = qJD(1) * t36 + qJD(5) * t429;
t454 = (t487 * t702 + t627 / 0.2e1 + t535 - t630 / 0.2e1 - t629 / 0.2e1) * t383;
t40 = t370 + (-t373 / 0.2e1 + t678) * t409 + t510 * t408 + t747 + (-t416 / 0.2e1 + t411 / 0.2e1 - t418 / 0.2e1 - t718 / 0.2e1) * t443 + (-t721 / 0.2e1 + t714 / 0.2e1 - t419 / 0.2e1 + t414 / 0.2e1) * t440 + t454;
t468 = Ifges(7,6) * t529 + Ifges(6,6) * t530 + t649 / 0.2e1 + t651 / 0.2e1 - t711 + t756 * t528;
t452 = (-pkin(5) * t86 + qJ(6) * t85) * t701 + t106 * mrSges(6,1) / 0.2e1 + t107 * t697 + t85 * t696 + t86 * t698 + t468;
t7 = t452 + ((mrSges(6,1) * t678 + mrSges(7,1) * t683 - t419 / 0.4e1 - t721 / 0.4e1 + t414 / 0.4e1 + t714 / 0.4e1) * t443 + (mrSges(6,2) * t677 + mrSges(7,3) * t683 + t418 / 0.4e1 + t416 / 0.4e1 + t718 / 0.4e1 - t411 / 0.4e1) * t440 + t464 * t424) * t471 + t477 * t702 - (-t720 / 0.4e1 - t719 / 0.4e1) * t346 + (t249 / 0.4e1 + t248 / 0.4e1 - t162 / 0.4e1 - t160 / 0.4e1 + (-t79 / 0.2e1 - t100 / 0.2e1) * mrSges(7,2) + t470 * t424) * t443 + (t158 / 0.4e1 - t156 / 0.4e1 + t251 / 0.4e1 - t250 / 0.4e1 + (-t101 / 0.2e1 + t78 / 0.2e1) * mrSges(7,2) + (t620 * t702 + t512) * t424) * t440 + t707;
t87 = t384 * t506 + t466 + t580;
t469 = t7 * qJD(1) + t40 * qJD(2) - t87 * qJD(3);
t211 = (-t493 + (t684 - t510) * m(7)) * t440;
t462 = (-t573 + (-t346 * t424 - t384 * t471) * t443) * t701;
t29 = t462 - t693 / 0.2e1 + t458;
t348 = t710 * t440;
t467 = -qJD(1) * t29 + qJD(2) * t211 + qJD(3) * t348;
t461 = mrSges(7,2) * t528 - t346 * t534 - t443 * t527 - t642 / 0.2e1 + t517;
t460 = (-mrSges(7,2) * t574 - pkin(5) * t430 + t717) * qJD(5);
t457 = qJD(5) * (-m(7) * t488 + t745);
t381 = m(7) * t582 + t430;
t354 = m(7) * t587 + t430;
t239 = m(7) * t529 + t669 * t689;
t212 = (m(7) * t684 + t493 + t498) * t440;
t41 = -t747 + t580 / 0.2e1 + t588 / 0.2e1 + t466 + t454 + t723 * t408 / 0.2e1;
t33 = t760 * t78 + t259;
t28 = t462 + t693 / 0.2e1 + t461;
t25 = t463 + t694 / 0.2e1 + t461;
t16 = t449 + t746;
t15 = t563 / 0.2e1 - t564 / 0.2e1 + t346 * t535 + mrSges(7,3) * t528 + t440 * t499 + (t440 * t621 - t443 * t620 + t496) * t701 + t725 * t672;
t12 = t471 * t509 + t450 - t453 + t502;
t8 = t451 + t452 + t384 * t472 / 0.2e1 + t473 * t677 + t477 * t701 + ((t440 * t620 + t748) * t701 + t708 + t734 * t764) * t424;
t5 = t448 - t455 + t468;
t1 = t447 - t446 + t459;
t17 = [-qJD(2) * t2 + qJD(3) * t4 + qJD(4) * t10 - qJD(5) * t9 + qJD(6) * t30, t1 * qJD(3) + t12 * qJD(4) + t5 * qJD(5) + t25 * qJD(6) - t647 + (m(4) * (-t360 * t444 + t441 * t722) * pkin(2) + t546 - t542 - t543 + (t374 * t485 + t786) * t760 + Ifges(3,5) * t445 - Ifges(3,6) * t442 + (t726 * t443 + (-t254 + t255) * t440) * t374 + 0.2e1 * (t374 * t479 + t778) * t703 + 0.2e1 * t768 * t705 + t549 + (-mrSges(3,1) * t445 + mrSges(3,2) * t442) * pkin(7) + t590 + t589 + mrSges(7,2) * t626 - mrSges(6,3) * t614 + t631 * t668 - t659 * t666 + t790) * qJD(2), t632 + t1 * qJD(2) + (t538 - t545 - t503 + t544 + t548 + t767 * t700 + m(6) * (t424 * t478 + t777) - t539 + t581 - t483 + m(7) * (t424 * t484 + t785) + t584 + t709 + mrSges(7,2) * t623 + t790) * qJD(3) + t16 * qJD(4) + t8 * qJD(5) + t28 * qJD(6), qJD(2) * t12 + qJD(3) * t16 + qJD(5) * t15 + qJD(6) * t239 + t617, t5 * qJD(2) + t8 * qJD(3) + t15 * qJD(4) + t33 * qJD(6) - t622 + ((-m(7) * pkin(5) + t733) * t101 + (-mrSges(6,2) + t429) * t100 + ((-mrSges(7,2) * qJ(6) - t754) * t443 + (mrSges(7,2) * pkin(5) - t756) * t440) * t471) * qJD(5), qJD(2) * t25 + qJD(3) * t28 + qJD(4) * t239 + qJD(5) * t33 + t616; -qJD(3) * t3 + qJD(4) * t11 + qJD(5) * t6 + qJD(6) * t26 + t647, qJD(3) * t34 + qJD(5) * t66 - qJD(6) * t296, t41 * qJD(5) + t212 * qJD(6) + t482 + ((m(6) * t425 - t619 * t700 - mrSges(5,1) - t494 + t710) * t382 + t740 + t732 * t743 + ((m(7) + m(6)) * t744 + t439 * t700 - mrSges(5,2)) * t383) * qJD(3), t610, t41 * qJD(3) + t354 * qJD(6) + t374 * t457 + t460 + t481, qJD(3) * t212 + qJD(5) * t354 - t476; qJD(2) * t3 + qJD(4) * t13 - qJD(5) * t7 + qJD(6) * t29 - t632, -qJD(5) * t40 - qJD(6) * t211 - t482, qJD(5) * t87 - qJD(6) * t348, t605, t381 * qJD(6) + t424 * t457 + t460 - t469, qJD(5) * t381 - t467; -qJD(2) * t11 - qJD(3) * t13 - qJD(5) * t14 + qJD(6) * t240 - t617, -t610, -t605, 0, -t604 + (t627 - t628 - t670) * qJD(5) + (m(7) * qJD(6) + qJD(5) * t733) * t440, qJD(5) * t669 + t551; -qJD(2) * t6 + qJD(3) * t7 + qJD(4) * t14 + qJD(6) * t36 + t622, qJD(3) * t40 - t481, t469, t604, t429 * qJD(6), t475; -qJD(2) * t26 - qJD(3) * t29 - qJD(4) * t240 - qJD(5) * t36 - t616, qJD(3) * t211 + t476, t467, -t551, -t475, 0;];
Cq  = t17;
