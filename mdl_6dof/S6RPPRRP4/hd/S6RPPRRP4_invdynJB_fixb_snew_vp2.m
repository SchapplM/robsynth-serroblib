% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-05-05 14:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:54:15
% EndTime: 2019-05-05 14:54:18
% DurationCPUTime: 3.01s
% Computational Cost: add. (35502->273), mult. (62831->316), div. (0->0), fcn. (30626->8), ass. (0->114)
t734 = Ifges(6,1) + Ifges(7,1);
t723 = Ifges(6,4) - Ifges(7,5);
t722 = -Ifges(6,5) - Ifges(7,4);
t733 = Ifges(6,2) + Ifges(7,3);
t720 = Ifges(6,6) - Ifges(7,6);
t732 = -Ifges(6,3) - Ifges(7,2);
t697 = qJD(1) ^ 2;
t692 = sin(qJ(1));
t695 = cos(qJ(1));
t672 = -t695 * g(1) - t692 * g(2);
t702 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t672;
t728 = -pkin(1) - pkin(2);
t648 = t728 * t697 + t702;
t671 = t692 * g(1) - t695 * g(2);
t701 = -t697 * qJ(2) + qJDD(2) - t671;
t651 = t728 * qJDD(1) + t701;
t689 = sin(pkin(9));
t690 = cos(pkin(9));
t623 = t690 * t648 + t689 * t651;
t620 = -t697 * pkin(3) - qJDD(1) * pkin(7) + t623;
t686 = g(3) + qJDD(3);
t691 = sin(qJ(4));
t694 = cos(qJ(4));
t615 = -t691 * t620 + t694 * t686;
t665 = (pkin(4) * t694 + pkin(8) * t691) * qJD(1);
t696 = qJD(4) ^ 2;
t713 = qJD(1) * t691;
t613 = -qJDD(4) * pkin(4) - t696 * pkin(8) - t665 * t713 - t615;
t693 = cos(qJ(5));
t727 = sin(qJ(5));
t663 = -t727 * qJD(4) + t693 * t713;
t711 = qJD(1) * qJD(4);
t708 = t694 * t711;
t666 = -t691 * qJDD(1) - t708;
t634 = -t663 * qJD(5) - t693 * qJDD(4) + t727 * t666;
t662 = t693 * qJD(4) + t727 * t713;
t635 = t662 * qJD(5) + t727 * qJDD(4) + t693 * t666;
t712 = t694 * qJD(1);
t674 = qJD(5) + t712;
t729 = 2 * qJD(6);
t607 = t663 * t729 + (-t662 * t674 - t635) * qJ(6) + (-t663 * t674 + t634) * pkin(5) + t613;
t646 = -t674 * mrSges(7,1) - t663 * mrSges(7,2);
t647 = t662 * mrSges(7,2) + t674 * mrSges(7,3);
t603 = m(7) * t607 + t634 * mrSges(7,1) - t635 * mrSges(7,3) + t663 * t646 - t662 * t647;
t622 = -t689 * t648 + t690 * t651;
t619 = qJDD(1) * pkin(3) - t697 * pkin(7) - t622;
t709 = t691 * t711;
t667 = -t694 * qJDD(1) + t709;
t611 = (-t666 + t708) * pkin(8) + (-t667 - t709) * pkin(4) + t619;
t616 = t694 * t620 + t691 * t686;
t614 = -t696 * pkin(4) + qJDD(4) * pkin(8) - t665 * t712 + t616;
t609 = t727 * t611 + t693 * t614;
t638 = -t662 * pkin(5) + t663 * qJ(6);
t661 = qJDD(5) - t667;
t673 = t674 ^ 2;
t605 = -t673 * pkin(5) + t661 * qJ(6) + t662 * t638 + t674 * t729 + t609;
t715 = t723 * t662 - t734 * t663 - t722 * t674;
t717 = t720 * t662 + t722 * t663 - t732 * t674;
t592 = -mrSges(6,1) * t613 - mrSges(7,1) * t607 + mrSges(7,2) * t605 + mrSges(6,3) * t609 - pkin(5) * t603 - t733 * t634 + t723 * t635 + t720 * t661 + t717 * t663 + t715 * t674;
t608 = t693 * t611 - t727 * t614;
t606 = -t661 * pkin(5) - t673 * qJ(6) - t663 * t638 + qJDD(6) - t608;
t716 = t733 * t662 - t723 * t663 + t720 * t674;
t593 = mrSges(6,2) * t613 + mrSges(7,2) * t606 - mrSges(6,3) * t608 - mrSges(7,3) * t607 - qJ(6) * t603 - t723 * t634 + t734 * t635 - t722 * t661 + t717 * t662 - t716 * t674;
t645 = t674 * mrSges(6,1) + t663 * mrSges(6,3);
t710 = m(7) * t605 + t661 * mrSges(7,3) + t674 * t646;
t639 = -t662 * mrSges(7,1) + t663 * mrSges(7,3);
t714 = -t662 * mrSges(6,1) - t663 * mrSges(6,2) + t639;
t725 = -mrSges(6,3) - mrSges(7,2);
t598 = m(6) * t609 - t661 * mrSges(6,2) + t725 * t634 - t674 * t645 + t714 * t662 + t710;
t644 = -t674 * mrSges(6,2) + t662 * mrSges(6,3);
t705 = -m(7) * t606 + t661 * mrSges(7,1) + t674 * t647;
t599 = m(6) * t608 + t661 * mrSges(6,1) + t725 * t635 + t674 * t644 + t714 * t663 + t705;
t595 = t693 * t598 - t727 * t599;
t600 = -m(6) * t613 - t634 * mrSges(6,1) - t635 * mrSges(6,2) + t662 * t644 + t663 * t645 - t603;
t655 = Ifges(5,6) * qJD(4) + (-Ifges(5,4) * t691 - Ifges(5,2) * t694) * qJD(1);
t656 = Ifges(5,5) * qJD(4) + (-Ifges(5,1) * t691 - Ifges(5,4) * t694) * qJD(1);
t731 = mrSges(5,1) * t615 - mrSges(5,2) * t616 + Ifges(5,5) * t666 + Ifges(5,6) * t667 + Ifges(5,3) * qJDD(4) + pkin(4) * t600 + pkin(8) * t595 - (t655 * t691 - t656 * t694) * qJD(1) + t693 * t592 + t727 * t593;
t602 = t635 * mrSges(7,2) - t663 * t639 - t705;
t730 = -t720 * t634 - t722 * t635 - t732 * t661 - t715 * t662 - t716 * t663 + mrSges(6,1) * t608 - mrSges(7,1) * t606 - mrSges(6,2) * t609 + mrSges(7,3) * t605 - pkin(5) * t602 + qJ(6) * (-t634 * mrSges(7,2) + t662 * t639 + t710);
t726 = mrSges(2,1) + mrSges(3,1);
t724 = Ifges(3,4) + Ifges(2,5);
t721 = Ifges(2,6) - Ifges(3,6);
t652 = -t697 * pkin(1) + t702;
t664 = (mrSges(5,1) * t694 - mrSges(5,2) * t691) * qJD(1);
t669 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t713;
t590 = m(5) * t616 - qJDD(4) * mrSges(5,2) + t667 * mrSges(5,3) - qJD(4) * t669 - t664 * t712 + t595;
t670 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t712;
t596 = m(5) * t615 + qJDD(4) * mrSges(5,1) - t666 * mrSges(5,3) + qJD(4) * t670 + t664 * t713 + t600;
t587 = t694 * t590 - t691 * t596;
t583 = m(4) * t623 - t697 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t587;
t594 = t727 * t598 + t693 * t599;
t591 = -m(5) * t619 + t667 * mrSges(5,1) - t666 * mrSges(5,2) + t669 * t713 - t670 * t712 - t594;
t588 = m(4) * t622 - qJDD(1) * mrSges(4,1) - t697 * mrSges(4,2) + t591;
t706 = t690 * t583 - t689 * t588;
t703 = m(3) * t652 + qJDD(1) * mrSges(3,3) + t706;
t575 = m(2) * t672 - qJDD(1) * mrSges(2,2) - t726 * t697 + t703;
t580 = t689 * t583 + t690 * t588;
t653 = -qJDD(1) * pkin(1) + t701;
t579 = m(3) * t653 - qJDD(1) * mrSges(3,1) - t697 * mrSges(3,3) + t580;
t576 = m(2) * t671 + qJDD(1) * mrSges(2,1) - t697 * mrSges(2,2) - t579;
t718 = t692 * t575 + t695 * t576;
t707 = t695 * t575 - t692 * t576;
t586 = t691 * t590 + t694 * t596;
t585 = m(4) * t686 + t586;
t654 = Ifges(5,3) * qJD(4) + (-Ifges(5,5) * t691 - Ifges(5,6) * t694) * qJD(1);
t571 = mrSges(5,2) * t619 - mrSges(5,3) * t615 + Ifges(5,1) * t666 + Ifges(5,4) * t667 + Ifges(5,5) * qJDD(4) - pkin(8) * t594 - qJD(4) * t655 - t727 * t592 + t693 * t593 - t654 * t712;
t581 = -mrSges(5,1) * t619 + mrSges(5,3) * t616 + Ifges(5,4) * t666 + Ifges(5,2) * t667 + Ifges(5,6) * qJDD(4) - pkin(4) * t594 + qJD(4) * t656 + t654 * t713 - t730;
t698 = -mrSges(3,1) * t653 - mrSges(4,1) * t622 - mrSges(2,2) * t672 - pkin(2) * t580 - pkin(3) * t591 - pkin(7) * t587 - t691 * t571 - t694 * t581 + qJ(2) * (-t697 * mrSges(3,1) + t703) - pkin(1) * t579 + mrSges(4,2) * t623 + mrSges(3,3) * t652 + mrSges(2,1) * t671 + (Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);
t584 = -m(3) * g(3) - t585;
t570 = -mrSges(4,1) * t686 + mrSges(4,3) * t623 + t697 * Ifges(4,5) - Ifges(4,6) * qJDD(1) - pkin(3) * t586 - t731;
t569 = mrSges(4,2) * t686 - mrSges(4,3) * t622 - Ifges(4,5) * qJDD(1) - t697 * Ifges(4,6) - pkin(7) * t586 + t694 * t571 - t691 * t581;
t568 = mrSges(3,2) * t653 - mrSges(2,3) * t671 - qJ(2) * t584 - qJ(3) * t580 + t690 * t569 - t689 * t570 - t721 * t697 + t724 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t567 = mrSges(3,2) * t652 + mrSges(2,3) * t672 - pkin(1) * t584 + pkin(2) * t585 + t726 * g(3) - qJ(3) * t706 + t721 * qJDD(1) - t689 * t569 - t690 * t570 + t724 * t697;
t1 = [-m(1) * g(1) + t707; -m(1) * g(2) + t718; (-m(1) - m(2) - m(3)) * g(3) - t585; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t718 - t692 * t567 + t695 * t568; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t707 + t695 * t567 + t692 * t568; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t698; t698; t579; t585; t731; t730; t602;];
tauJB  = t1;
