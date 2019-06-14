% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPRP2
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 03:50
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:44:01
% EndTime: 2019-05-05 03:44:14
% DurationCPUTime: 9.41s
% Computational Cost: add. (142710->317), mult. (304540->400), div. (0->0), fcn. (213994->12), ass. (0->132)
t712 = -2 * qJD(4);
t711 = Ifges(6,1) + Ifges(7,1);
t705 = Ifges(6,4) - Ifges(7,5);
t710 = -Ifges(6,5) - Ifges(7,4);
t709 = Ifges(6,2) + Ifges(7,3);
t703 = Ifges(6,6) - Ifges(7,6);
t708 = -Ifges(6,3) - Ifges(7,2);
t664 = sin(pkin(10));
t667 = cos(pkin(10));
t653 = g(1) * t664 - g(2) * t667;
t654 = -g(1) * t667 - g(2) * t664;
t662 = -g(3) + qJDD(1);
t673 = cos(qJ(2));
t668 = cos(pkin(6));
t671 = sin(qJ(2));
t700 = t668 * t671;
t665 = sin(pkin(6));
t701 = t665 * t671;
t618 = t653 * t700 + t673 * t654 + t662 * t701;
t675 = qJD(2) ^ 2;
t613 = -pkin(2) * t675 + qJDD(2) * pkin(8) + t618;
t634 = -t653 * t665 + t662 * t668;
t670 = sin(qJ(3));
t672 = cos(qJ(3));
t587 = -t613 * t670 + t672 * t634;
t691 = qJD(2) * qJD(3);
t689 = t672 * t691;
t651 = qJDD(2) * t670 + t689;
t584 = (-t651 + t689) * qJ(4) + (t670 * t672 * t675 + qJDD(3)) * pkin(3) + t587;
t588 = t672 * t613 + t670 * t634;
t652 = qJDD(2) * t672 - t670 * t691;
t694 = qJD(2) * t670;
t655 = qJD(3) * pkin(3) - qJ(4) * t694;
t661 = t672 ^ 2;
t585 = -pkin(3) * t661 * t675 + qJ(4) * t652 - qJD(3) * t655 + t588;
t663 = sin(pkin(11));
t666 = cos(pkin(11));
t640 = (t663 * t672 + t666 * t670) * qJD(2);
t577 = t584 * t666 - t663 * t585 + t640 * t712;
t639 = (t663 * t670 - t666 * t672) * qJD(2);
t617 = -t671 * t654 + (t653 * t668 + t662 * t665) * t673;
t707 = cos(qJ(5));
t706 = -mrSges(6,3) - mrSges(7,2);
t679 = -qJDD(2) * pkin(2) - t617;
t612 = -pkin(8) * t675 + t679;
t656 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t694;
t693 = qJD(2) * t672;
t657 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t693;
t578 = t663 * t584 + t666 * t585 + t639 * t712;
t621 = pkin(4) * t639 - pkin(9) * t640;
t674 = qJD(3) ^ 2;
t576 = -pkin(4) * t674 + qJDD(3) * pkin(9) - t621 * t639 + t578;
t586 = -pkin(3) * t652 + qJDD(4) + t655 * t694 + (-qJ(4) * t661 - pkin(8)) * t675 + t679;
t626 = -t651 * t663 + t652 * t666;
t627 = t651 * t666 + t652 * t663;
t580 = (qJD(3) * t639 - t627) * pkin(9) + (qJD(3) * t640 - t626) * pkin(4) + t586;
t669 = sin(qJ(5));
t573 = t707 * t576 + t669 * t580;
t629 = t669 * qJD(3) + t640 * t707;
t599 = qJD(5) * t629 - qJDD(3) * t707 + t627 * t669;
t638 = qJD(5) + t639;
t610 = mrSges(6,1) * t638 - mrSges(6,3) * t629;
t625 = qJDD(5) - t626;
t628 = -qJD(3) * t707 + t640 * t669;
t603 = pkin(5) * t628 - qJ(6) * t629;
t637 = t638 ^ 2;
t569 = -pkin(5) * t637 + qJ(6) * t625 + 0.2e1 * qJD(6) * t638 - t603 * t628 + t573;
t611 = -mrSges(7,1) * t638 + mrSges(7,2) * t629;
t690 = m(7) * t569 + t625 * mrSges(7,3) + t638 * t611;
t604 = mrSges(7,1) * t628 - mrSges(7,3) * t629;
t695 = -mrSges(6,1) * t628 - mrSges(6,2) * t629 - t604;
t564 = m(6) * t573 - mrSges(6,2) * t625 + t599 * t706 - t610 * t638 + t628 * t695 + t690;
t572 = -t669 * t576 + t580 * t707;
t600 = -t628 * qJD(5) + t669 * qJDD(3) + t627 * t707;
t609 = -mrSges(6,2) * t638 - mrSges(6,3) * t628;
t570 = -t625 * pkin(5) - t637 * qJ(6) + t629 * t603 + qJDD(6) - t572;
t608 = -mrSges(7,2) * t628 + mrSges(7,3) * t638;
t683 = -m(7) * t570 + t625 * mrSges(7,1) + t638 * t608;
t566 = m(6) * t572 + mrSges(6,1) * t625 + t600 * t706 + t609 * t638 + t629 * t695 + t683;
t559 = t669 * t564 + t707 * t566;
t632 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t639;
t633 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t640;
t678 = m(5) * t586 - t626 * mrSges(5,1) + mrSges(5,2) * t627 + t639 * t632 + t633 * t640 + t559;
t676 = -m(4) * t612 + t652 * mrSges(4,1) - mrSges(4,2) * t651 - t656 * t694 + t657 * t693 - t678;
t553 = m(3) * t617 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t675 + t676;
t702 = t553 * t673;
t620 = mrSges(5,1) * t639 + mrSges(5,2) * t640;
t685 = t707 * t564 - t566 * t669;
t556 = m(5) * t578 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t626 - qJD(3) * t633 - t620 * t639 + t685;
t575 = -qJDD(3) * pkin(4) - pkin(9) * t674 + t640 * t621 - t577;
t571 = -0.2e1 * qJD(6) * t629 + (t628 * t638 - t600) * qJ(6) + (t629 * t638 + t599) * pkin(5) + t575;
t567 = m(7) * t571 + mrSges(7,1) * t599 - t600 * mrSges(7,3) + t608 * t628 - t629 * t611;
t677 = -m(6) * t575 - t599 * mrSges(6,1) - mrSges(6,2) * t600 - t628 * t609 - t610 * t629 - t567;
t561 = m(5) * t577 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t627 + qJD(3) * t632 - t620 * t640 + t677;
t550 = t663 * t556 + t666 * t561;
t650 = (-mrSges(4,1) * t672 + mrSges(4,2) * t670) * qJD(2);
t548 = m(4) * t587 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t651 + qJD(3) * t657 - t650 * t694 + t550;
t686 = t666 * t556 - t561 * t663;
t549 = m(4) * t588 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t652 - qJD(3) * t656 + t650 * t693 + t686;
t687 = -t548 * t670 + t672 * t549;
t539 = m(3) * t618 - mrSges(3,1) * t675 - qJDD(2) * mrSges(3,2) + t687;
t542 = t672 * t548 + t670 * t549;
t541 = m(3) * t634 + t542;
t530 = t539 * t700 - t541 * t665 + t668 * t702;
t528 = m(2) * t653 + t530;
t535 = t673 * t539 - t553 * t671;
t534 = m(2) * t654 + t535;
t699 = t667 * t528 + t664 * t534;
t698 = t628 * t709 - t629 * t705 - t638 * t703;
t697 = t628 * t703 + t629 * t710 + t638 * t708;
t696 = -t628 * t705 + t629 * t711 - t638 * t710;
t529 = t539 * t701 + t668 * t541 + t665 * t702;
t688 = -t528 * t664 + t667 * t534;
t557 = -mrSges(6,1) * t575 - mrSges(7,1) * t571 + mrSges(7,2) * t569 + mrSges(6,3) * t573 - pkin(5) * t567 - t599 * t709 + t705 * t600 + t703 * t625 + t697 * t629 + t696 * t638;
t558 = mrSges(6,2) * t575 + mrSges(7,2) * t570 - mrSges(6,3) * t572 - mrSges(7,3) * t571 - qJ(6) * t567 - t705 * t599 + t600 * t711 - t625 * t710 + t697 * t628 + t698 * t638;
t614 = Ifges(5,5) * t640 - Ifges(5,6) * t639 + Ifges(5,3) * qJD(3);
t615 = Ifges(5,4) * t640 - Ifges(5,2) * t639 + Ifges(5,6) * qJD(3);
t543 = mrSges(5,2) * t586 - mrSges(5,3) * t577 + Ifges(5,1) * t627 + Ifges(5,4) * t626 + Ifges(5,5) * qJDD(3) - pkin(9) * t559 - qJD(3) * t615 - t669 * t557 + t558 * t707 - t639 * t614;
t616 = Ifges(5,1) * t640 - Ifges(5,4) * t639 + Ifges(5,5) * qJD(3);
t544 = Ifges(5,4) * t627 + Ifges(5,2) * t626 + Ifges(5,6) * qJDD(3) - t640 * t614 + qJD(3) * t616 - mrSges(5,1) * t586 + mrSges(5,3) * t578 - mrSges(6,1) * t572 + mrSges(6,2) * t573 + mrSges(7,1) * t570 - mrSges(7,3) * t569 - pkin(5) * t683 - qJ(6) * t690 - pkin(4) * t559 + (pkin(5) * t604 + t698) * t629 + (qJ(6) * t604 - t696) * t628 + t708 * t625 + (mrSges(7,2) * pkin(5) + t710) * t600 + (mrSges(7,2) * qJ(6) + t703) * t599;
t642 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t670 + Ifges(4,6) * t672) * qJD(2);
t644 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t670 + Ifges(4,4) * t672) * qJD(2);
t526 = -mrSges(4,1) * t612 + mrSges(4,3) * t588 + Ifges(4,4) * t651 + Ifges(4,2) * t652 + Ifges(4,6) * qJDD(3) - pkin(3) * t678 + qJ(4) * t686 + qJD(3) * t644 + t663 * t543 + t666 * t544 - t642 * t694;
t643 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t670 + Ifges(4,2) * t672) * qJD(2);
t531 = mrSges(4,2) * t612 - mrSges(4,3) * t587 + Ifges(4,1) * t651 + Ifges(4,4) * t652 + Ifges(4,5) * qJDD(3) - qJ(4) * t550 - qJD(3) * t643 + t543 * t666 - t544 * t663 + t642 * t693;
t524 = mrSges(3,2) * t634 - mrSges(3,3) * t617 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t675 - pkin(8) * t542 - t526 * t670 + t531 * t672;
t525 = t675 * Ifges(3,5) - t640 * t615 - t639 * t616 - pkin(2) * t542 + mrSges(3,3) * t618 - mrSges(3,1) * t634 + Ifges(3,6) * qJDD(2) - pkin(3) * t550 - Ifges(4,5) * t651 - Ifges(4,6) * t652 - mrSges(4,1) * t587 + mrSges(4,2) * t588 - t707 * t557 - pkin(4) * t677 - pkin(9) * t685 - Ifges(5,5) * t627 - Ifges(5,6) * t626 - mrSges(5,1) * t577 + mrSges(5,2) * t578 - t669 * t558 + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) + (-t643 * t670 + t644 * t672) * qJD(2);
t680 = pkin(7) * t535 + t524 * t671 + t525 * t673;
t523 = mrSges(3,1) * t617 - mrSges(3,2) * t618 + Ifges(3,3) * qJDD(2) + pkin(2) * t676 + pkin(8) * t687 + t672 * t526 + t670 * t531;
t522 = mrSges(2,2) * t662 - mrSges(2,3) * t653 + t524 * t673 - t525 * t671 + (-t529 * t665 - t530 * t668) * pkin(7);
t521 = -mrSges(2,1) * t662 + mrSges(2,3) * t654 - pkin(1) * t529 - t523 * t665 + t668 * t680;
t1 = [-m(1) * g(1) + t688; -m(1) * g(2) + t699; -m(1) * g(3) + m(2) * t662 + t529; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t699 - t664 * t521 + t667 * t522; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t688 + t667 * t521 + t664 * t522; -mrSges(1,1) * g(2) + mrSges(2,1) * t653 + mrSges(1,2) * g(1) - mrSges(2,2) * t654 + pkin(1) * t530 + t523 * t668 + t665 * t680;];
tauB  = t1;
