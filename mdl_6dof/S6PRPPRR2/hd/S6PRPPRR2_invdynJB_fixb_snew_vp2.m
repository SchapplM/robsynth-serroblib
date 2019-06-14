% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 21:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPPRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:50:34
% EndTime: 2019-05-04 21:50:39
% DurationCPUTime: 4.77s
% Computational Cost: add. (62187->250), mult. (109338->311), div. (0->0), fcn. (72084->12), ass. (0->118)
t689 = sin(pkin(10));
t692 = cos(pkin(10));
t672 = t689 * g(1) - t692 * g(2);
t673 = -t692 * g(1) - t689 * g(2);
t685 = -g(3) + qJDD(1);
t696 = sin(qJ(2));
t693 = cos(pkin(6));
t699 = cos(qJ(2));
t724 = t693 * t699;
t690 = sin(pkin(6));
t726 = t690 * t699;
t635 = t672 * t724 - t673 * t696 + t685 * t726;
t633 = qJDD(2) * pkin(2) + t635;
t725 = t693 * t696;
t727 = t690 * t696;
t636 = t672 * t725 + t699 * t673 + t685 * t727;
t701 = qJD(2) ^ 2;
t634 = -pkin(2) * t701 + t636;
t688 = sin(pkin(11));
t691 = cos(pkin(11));
t629 = t688 * t633 + t691 * t634;
t732 = -qJDD(2) * qJ(4) - (2 * qJD(4) * qJD(2)) - t629;
t731 = -pkin(3) - pkin(8);
t730 = mrSges(4,1) - mrSges(5,2);
t729 = -Ifges(5,4) + Ifges(4,5);
t728 = Ifges(5,5) - Ifges(4,6);
t651 = -t690 * t672 + t693 * t685;
t650 = qJDD(3) + t651;
t695 = sin(qJ(5));
t723 = t695 * t650;
t628 = t691 * t633 - t688 * t634;
t708 = -t701 * qJ(4) + qJDD(4) - t628;
t625 = qJDD(2) * t731 + t708;
t698 = cos(qJ(5));
t621 = t695 * t625 + t698 * t650;
t668 = (mrSges(6,1) * t695 + mrSges(6,2) * t698) * qJD(2);
t719 = qJD(2) * qJD(5);
t715 = t698 * t719;
t670 = -t695 * qJDD(2) - t715;
t721 = qJD(2) * t698;
t675 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t721;
t669 = (pkin(5) * t695 - pkin(9) * t698) * qJD(2);
t700 = qJD(5) ^ 2;
t720 = t695 * qJD(2);
t617 = -pkin(5) * t700 + qJDD(5) * pkin(9) - t669 * t720 + t621;
t624 = t701 * t731 - t732;
t716 = t695 * t719;
t671 = t698 * qJDD(2) - t716;
t618 = (-t671 + t716) * pkin(9) + (-t670 + t715) * pkin(5) + t624;
t694 = sin(qJ(6));
t697 = cos(qJ(6));
t613 = -t617 * t694 + t618 * t697;
t666 = t697 * qJD(5) - t694 * t721;
t643 = t666 * qJD(6) + t694 * qJDD(5) + t697 * t671;
t667 = t694 * qJD(5) + t697 * t721;
t644 = -t666 * mrSges(7,1) + t667 * mrSges(7,2);
t678 = qJD(6) + t720;
t648 = -t678 * mrSges(7,2) + t666 * mrSges(7,3);
t664 = qJDD(6) - t670;
t611 = m(7) * t613 + mrSges(7,1) * t664 - mrSges(7,3) * t643 - t644 * t667 + t648 * t678;
t614 = t617 * t697 + t618 * t694;
t642 = -t667 * qJD(6) + t697 * qJDD(5) - t694 * t671;
t649 = t678 * mrSges(7,1) - t667 * mrSges(7,3);
t612 = m(7) * t614 - mrSges(7,2) * t664 + mrSges(7,3) * t642 + t644 * t666 - t649 * t678;
t711 = -t694 * t611 + t697 * t612;
t600 = m(6) * t621 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t670 - qJD(5) * t675 - t668 * t720 + t711;
t620 = t625 * t698 - t723;
t674 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t720;
t616 = -qJDD(5) * pkin(5) - t700 * pkin(9) + t723 + (qJD(2) * t669 - t625) * t698;
t705 = -m(7) * t616 + t642 * mrSges(7,1) - mrSges(7,2) * t643 + t666 * t648 - t649 * t667;
t607 = m(6) * t620 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t671 + qJD(5) * t674 - t668 * t721 + t705;
t595 = t695 * t600 + t698 * t607;
t627 = -qJDD(2) * pkin(3) + t708;
t707 = -m(5) * t627 + t701 * mrSges(5,3) - t595;
t589 = m(4) * t628 - t701 * mrSges(4,2) + qJDD(2) * t730 + t707;
t626 = t701 * pkin(3) + t732;
t602 = t697 * t611 + t694 * t612;
t706 = -m(6) * t624 + mrSges(6,1) * t670 - t671 * mrSges(6,2) - t674 * t720 - t675 * t721 - t602;
t703 = -m(5) * t626 + t701 * mrSges(5,2) + qJDD(2) * mrSges(5,3) - t706;
t598 = m(4) * t629 - mrSges(4,1) * t701 - qJDD(2) * mrSges(4,2) + t703;
t585 = t691 * t589 + t688 * t598;
t583 = m(3) * t635 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t701 + t585;
t713 = -t589 * t688 + t691 * t598;
t584 = m(3) * t636 - mrSges(3,1) * t701 - qJDD(2) * mrSges(3,2) + t713;
t712 = t698 * t600 - t607 * t695;
t594 = m(5) * t650 + t712;
t593 = m(4) * t650 + t594;
t592 = m(3) * t651 + t593;
t572 = t583 * t724 + t584 * t725 - t592 * t690;
t570 = m(2) * t672 + t572;
t577 = -t583 * t696 + t699 * t584;
t576 = m(2) * t673 + t577;
t722 = t692 * t570 + t689 * t576;
t571 = t583 * t726 + t584 * t727 + t693 * t592;
t714 = -t689 * t570 + t692 * t576;
t710 = m(2) * t685 + t571;
t637 = Ifges(7,5) * t667 + Ifges(7,6) * t666 + Ifges(7,3) * t678;
t639 = Ifges(7,1) * t667 + Ifges(7,4) * t666 + Ifges(7,5) * t678;
t605 = -mrSges(7,1) * t616 + mrSges(7,3) * t614 + Ifges(7,4) * t643 + Ifges(7,2) * t642 + Ifges(7,6) * t664 - t637 * t667 + t639 * t678;
t638 = Ifges(7,4) * t667 + Ifges(7,2) * t666 + Ifges(7,6) * t678;
t606 = mrSges(7,2) * t616 - mrSges(7,3) * t613 + Ifges(7,1) * t643 + Ifges(7,4) * t642 + Ifges(7,5) * t664 + t637 * t666 - t638 * t678;
t656 = (Ifges(6,3) * qJD(5)) + (Ifges(6,5) * t698 - Ifges(6,6) * t695) * qJD(2);
t657 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t698 - Ifges(6,2) * t695) * qJD(2);
t586 = mrSges(6,2) * t624 - mrSges(6,3) * t620 + Ifges(6,1) * t671 + Ifges(6,4) * t670 + Ifges(6,5) * qJDD(5) - pkin(9) * t602 - qJD(5) * t657 - t605 * t694 + t606 * t697 - t656 * t720;
t658 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t698 - Ifges(6,4) * t695) * qJD(2);
t702 = mrSges(7,1) * t613 - mrSges(7,2) * t614 + Ifges(7,5) * t643 + Ifges(7,6) * t642 + Ifges(7,3) * t664 + t638 * t667 - t639 * t666;
t587 = -mrSges(6,1) * t624 + mrSges(6,3) * t621 + Ifges(6,4) * t671 + Ifges(6,2) * t670 + Ifges(6,6) * qJDD(5) - pkin(5) * t602 + qJD(5) * t658 - t656 * t721 - t702;
t568 = -mrSges(5,1) * t626 + mrSges(4,3) * t629 - pkin(3) * t594 - pkin(4) * t706 - pkin(8) * t712 - qJDD(2) * t728 - t695 * t586 - t698 * t587 - t650 * t730 + t701 * t729;
t704 = mrSges(6,1) * t620 - mrSges(6,2) * t621 + Ifges(6,5) * t671 + Ifges(6,6) * t670 + Ifges(6,3) * qJDD(5) + pkin(5) * t705 + pkin(9) * t711 + t697 * t605 + t694 * t606 + t657 * t721 + t658 * t720;
t573 = t704 + t728 * t701 + (mrSges(4,2) - mrSges(5,3)) * t650 + t729 * qJDD(2) + mrSges(5,1) * t627 - mrSges(4,3) * t628 - qJ(4) * t594 + pkin(4) * t595;
t565 = -mrSges(3,1) * t651 + mrSges(3,3) * t636 + t701 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t593 + qJ(3) * t713 + t691 * t568 + t688 * t573;
t566 = mrSges(3,2) * t651 - mrSges(3,3) * t635 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t701 - qJ(3) * t585 - t568 * t688 + t691 * t573;
t709 = pkin(7) * t577 + t565 * t699 + t566 * t696;
t590 = qJDD(2) * mrSges(5,2) - t707;
t567 = pkin(2) * t585 + mrSges(3,1) * t635 - mrSges(3,2) * t636 - pkin(3) * t590 + qJ(4) * t703 + t698 * t586 - t695 * t587 - pkin(8) * t595 + mrSges(4,1) * t628 - mrSges(4,2) * t629 + mrSges(5,2) * t627 - mrSges(5,3) * t626 + (Ifges(3,3) + Ifges(4,3) + Ifges(5,1)) * qJDD(2);
t564 = mrSges(2,2) * t685 - mrSges(2,3) * t672 - t696 * t565 + t699 * t566 + (-t571 * t690 - t572 * t693) * pkin(7);
t563 = -mrSges(2,1) * t685 + mrSges(2,3) * t673 - pkin(1) * t571 - t690 * t567 + t693 * t709;
t1 = [-m(1) * g(1) + t714; -m(1) * g(2) + t722; -m(1) * g(3) + t710; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t722 - t689 * t563 + t692 * t564; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t714 + t692 * t563 + t689 * t564; -mrSges(1,1) * g(2) + mrSges(2,1) * t672 + mrSges(1,2) * g(1) - mrSges(2,2) * t673 + pkin(1) * t572 + t693 * t567 + t690 * t709; t710; t567; t593; t590; t704; t702;];
tauJB  = t1;
