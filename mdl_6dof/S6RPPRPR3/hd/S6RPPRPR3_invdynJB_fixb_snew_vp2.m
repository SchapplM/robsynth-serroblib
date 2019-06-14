% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 14:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:09:18
% EndTime: 2019-05-05 14:09:23
% DurationCPUTime: 4.96s
% Computational Cost: add. (61150->296), mult. (124194->362), div. (0->0), fcn. (74973->10), ass. (0->123)
t755 = sin(qJ(1));
t758 = cos(qJ(1));
t730 = t755 * g(1) - t758 * g(2);
t721 = qJDD(1) * pkin(1) + t730;
t731 = -t758 * g(1) - t755 * g(2);
t760 = qJD(1) ^ 2;
t723 = -t760 * pkin(1) + t731;
t750 = sin(pkin(9));
t752 = cos(pkin(9));
t696 = t750 * t721 + t752 * t723;
t771 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t696;
t749 = sin(pkin(10));
t751 = cos(pkin(10));
t754 = sin(qJ(4));
t757 = cos(qJ(4));
t707 = (t749 * t757 + t751 * t754) * qJD(1);
t788 = 2 * qJD(5);
t787 = -pkin(2) - pkin(7);
t786 = pkin(4) * t760;
t785 = mrSges(3,1) - mrSges(4,2);
t784 = Ifges(3,5) - Ifges(4,4);
t783 = -Ifges(3,6) + Ifges(4,5);
t695 = t752 * t721 - t750 * t723;
t768 = -t760 * qJ(3) + qJDD(3) - t695;
t680 = t787 * qJDD(1) + t768;
t675 = t757 * t680;
t779 = qJD(1) * qJD(4);
t725 = t757 * qJDD(1) - t754 * t779;
t746 = -g(3) + qJDD(2);
t658 = qJDD(4) * pkin(4) - t725 * qJ(5) + t675 + (-qJ(5) * t779 - t757 * t786 - t746) * t754;
t672 = t754 * t680 + t757 * t746;
t724 = -t754 * qJDD(1) - t757 * t779;
t780 = qJD(1) * t757;
t728 = qJD(4) * pkin(4) - qJ(5) * t780;
t745 = t754 ^ 2;
t659 = t724 * qJ(5) - qJD(4) * t728 - t745 * t786 + t672;
t654 = t749 * t658 + t751 * t659 - t707 * t788;
t781 = qJD(1) * t754;
t708 = -t749 * t781 + t751 * t780;
t688 = t707 * mrSges(6,1) + t708 * mrSges(6,2);
t697 = t751 * t724 - t749 * t725;
t702 = qJD(4) * mrSges(6,1) - t708 * mrSges(6,3);
t689 = t707 * pkin(5) - t708 * pkin(8);
t759 = qJD(4) ^ 2;
t651 = -t759 * pkin(5) + qJDD(4) * pkin(8) - t707 * t689 + t654;
t661 = -t724 * pkin(4) + qJDD(5) + t728 * t780 + (-qJ(5) * t745 + t787) * t760 + t771;
t698 = t749 * t724 + t751 * t725;
t655 = (qJD(4) * t707 - t698) * pkin(8) + (qJD(4) * t708 - t697) * pkin(5) + t661;
t753 = sin(qJ(6));
t756 = cos(qJ(6));
t648 = -t753 * t651 + t756 * t655;
t699 = t756 * qJD(4) - t753 * t708;
t668 = t699 * qJD(6) + t753 * qJDD(4) + t756 * t698;
t700 = t753 * qJD(4) + t756 * t708;
t673 = -t699 * mrSges(7,1) + t700 * mrSges(7,2);
t706 = qJD(6) + t707;
t677 = -t706 * mrSges(7,2) + t699 * mrSges(7,3);
t694 = qJDD(6) - t697;
t645 = m(7) * t648 + t694 * mrSges(7,1) - t668 * mrSges(7,3) - t700 * t673 + t706 * t677;
t649 = t756 * t651 + t753 * t655;
t667 = -t700 * qJD(6) + t756 * qJDD(4) - t753 * t698;
t678 = t706 * mrSges(7,1) - t700 * mrSges(7,3);
t646 = m(7) * t649 - t694 * mrSges(7,2) + t667 * mrSges(7,3) + t699 * t673 - t706 * t678;
t772 = -t753 * t645 + t756 * t646;
t632 = m(6) * t654 - qJDD(4) * mrSges(6,2) + t697 * mrSges(6,3) - qJD(4) * t702 - t707 * t688 + t772;
t770 = -t751 * t658 + t749 * t659;
t653 = -0.2e1 * qJD(5) * t708 - t770;
t701 = -qJD(4) * mrSges(6,2) - t707 * mrSges(6,3);
t650 = -qJDD(4) * pkin(5) - t759 * pkin(8) + (t788 + t689) * t708 + t770;
t766 = -m(7) * t650 + t667 * mrSges(7,1) - t668 * mrSges(7,2) + t699 * t677 - t700 * t678;
t641 = m(6) * t653 + qJDD(4) * mrSges(6,1) - t698 * mrSges(6,3) + qJD(4) * t701 - t708 * t688 + t766;
t625 = t749 * t632 + t751 * t641;
t671 = -t754 * t746 + t675;
t722 = (mrSges(5,1) * t754 + mrSges(5,2) * t757) * qJD(1);
t727 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t781;
t622 = m(5) * t671 + qJDD(4) * mrSges(5,1) - t725 * mrSges(5,3) + qJD(4) * t727 - t722 * t780 + t625;
t729 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t780;
t773 = t751 * t632 - t749 * t641;
t623 = m(5) * t672 - qJDD(4) * mrSges(5,2) + t724 * mrSges(5,3) - qJD(4) * t729 - t722 * t781 + t773;
t618 = t757 * t622 + t754 * t623;
t683 = -qJDD(1) * pkin(2) + t768;
t767 = -m(4) * t683 + t760 * mrSges(4,3) - t618;
t613 = m(3) * t695 - t760 * mrSges(3,2) + t785 * qJDD(1) + t767;
t681 = t760 * pkin(2) - t771;
t635 = t756 * t645 + t753 * t646;
t633 = m(6) * t661 - t697 * mrSges(6,1) + t698 * mrSges(6,2) + t707 * t701 + t708 * t702 + t635;
t679 = t787 * t760 + t771;
t765 = -m(5) * t679 + t724 * mrSges(5,1) - t725 * mrSges(5,2) - t727 * t781 - t729 * t780 - t633;
t763 = -m(4) * t681 + t760 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t765;
t628 = m(3) * t696 - t760 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t763;
t611 = t752 * t613 + t750 * t628;
t608 = m(2) * t730 + qJDD(1) * mrSges(2,1) - t760 * mrSges(2,2) + t611;
t775 = -t750 * t613 + t752 * t628;
t609 = m(2) * t731 - t760 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t775;
t782 = t758 * t608 + t755 * t609;
t776 = -t755 * t608 + t758 * t609;
t774 = -t754 * t622 + t757 * t623;
t617 = m(4) * t746 + t774;
t616 = m(3) * t746 + t617;
t663 = Ifges(7,4) * t700 + Ifges(7,2) * t699 + Ifges(7,6) * t706;
t664 = Ifges(7,1) * t700 + Ifges(7,4) * t699 + Ifges(7,5) * t706;
t764 = mrSges(7,1) * t648 - mrSges(7,2) * t649 + Ifges(7,5) * t668 + Ifges(7,6) * t667 + Ifges(7,3) * t694 + t700 * t663 - t699 * t664;
t662 = Ifges(7,5) * t700 + Ifges(7,6) * t699 + Ifges(7,3) * t706;
t638 = -mrSges(7,1) * t650 + mrSges(7,3) * t649 + Ifges(7,4) * t668 + Ifges(7,2) * t667 + Ifges(7,6) * t694 - t700 * t662 + t706 * t664;
t639 = mrSges(7,2) * t650 - mrSges(7,3) * t648 + Ifges(7,1) * t668 + Ifges(7,4) * t667 + Ifges(7,5) * t694 + t699 * t662 - t706 * t663;
t684 = Ifges(6,5) * t708 - Ifges(6,6) * t707 + (Ifges(6,3) * qJD(4));
t685 = Ifges(6,4) * t708 - Ifges(6,2) * t707 + Ifges(6,6) * qJD(4);
t619 = mrSges(6,2) * t661 - mrSges(6,3) * t653 + Ifges(6,1) * t698 + Ifges(6,4) * t697 + Ifges(6,5) * qJDD(4) - pkin(8) * t635 - qJD(4) * t685 - t753 * t638 + t756 * t639 - t707 * t684;
t686 = Ifges(6,1) * t708 - Ifges(6,4) * t707 + Ifges(6,5) * qJD(4);
t620 = -mrSges(6,1) * t661 + mrSges(6,3) * t654 + Ifges(6,4) * t698 + Ifges(6,2) * t697 + Ifges(6,6) * qJDD(4) - pkin(5) * t635 + qJD(4) * t686 - t708 * t684 - t764;
t712 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t757 - Ifges(5,6) * t754) * qJD(1);
t714 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t757 - Ifges(5,4) * t754) * qJD(1);
t602 = -mrSges(5,1) * t679 + mrSges(5,3) * t672 + Ifges(5,4) * t725 + Ifges(5,2) * t724 + Ifges(5,6) * qJDD(4) - pkin(4) * t633 + qJ(5) * t773 + qJD(4) * t714 + t749 * t619 + t751 * t620 - t712 * t780;
t713 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t757 - Ifges(5,2) * t754) * qJD(1);
t604 = mrSges(5,2) * t679 - mrSges(5,3) * t671 + Ifges(5,1) * t725 + Ifges(5,4) * t724 + Ifges(5,5) * qJDD(4) - qJ(5) * t625 - qJD(4) * t713 + t751 * t619 - t749 * t620 - t712 * t781;
t615 = qJDD(1) * mrSges(4,2) - t767;
t762 = mrSges(2,1) * t730 + mrSges(3,1) * t695 - mrSges(2,2) * t731 - mrSges(3,2) * t696 + mrSges(4,2) * t683 - mrSges(4,3) * t681 + pkin(1) * t611 - pkin(2) * t615 - pkin(7) * t618 + qJ(3) * t763 - t754 * t602 + t757 * t604 + (Ifges(2,3) + Ifges(3,3) + Ifges(4,1)) * qJDD(1);
t761 = mrSges(5,1) * t671 + mrSges(6,1) * t653 - mrSges(5,2) * t672 - mrSges(6,2) * t654 + Ifges(6,6) * t697 + pkin(4) * t625 + pkin(5) * t766 + pkin(8) * t772 + t756 * t638 + t753 * t639 + t708 * t685 + t707 * t686 + Ifges(6,5) * t698 + t714 * t781 + t713 * t780 + Ifges(5,6) * t724 + Ifges(5,5) * t725 + (Ifges(5,3) + Ifges(6,3)) * qJDD(4);
t601 = t761 + pkin(3) * t618 - qJ(3) * t617 + t784 * qJDD(1) + t783 * t760 + (mrSges(3,2) - mrSges(4,3)) * t746 + mrSges(4,1) * t683 - mrSges(3,3) * t695;
t600 = -mrSges(4,1) * t681 + mrSges(3,3) * t696 - pkin(2) * t617 - pkin(3) * t765 - pkin(7) * t774 - t783 * qJDD(1) - t757 * t602 - t754 * t604 - t785 * t746 + t784 * t760;
t599 = -mrSges(2,2) * g(3) - mrSges(2,3) * t730 + Ifges(2,5) * qJDD(1) - t760 * Ifges(2,6) - qJ(2) * t611 - t750 * t600 + t752 * t601;
t598 = mrSges(2,1) * g(3) + mrSges(2,3) * t731 + t760 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t616 + qJ(2) * t775 + t752 * t600 + t750 * t601;
t1 = [-m(1) * g(1) + t776; -m(1) * g(2) + t782; (-m(1) - m(2)) * g(3) + t616; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t782 - t755 * t598 + t758 * t599; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t776 + t758 * t598 + t755 * t599; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t762; t762; t616; t615; t761; t633; t764;];
tauJB  = t1;
