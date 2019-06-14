% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-05-05 00:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:59:40
% EndTime: 2019-05-04 23:59:47
% DurationCPUTime: 4.13s
% Computational Cost: add. (49153->268), mult. (89566->323), div. (0->0), fcn. (55663->10), ass. (0->121)
t823 = Ifges(6,1) + Ifges(7,1);
t813 = Ifges(6,4) - Ifges(7,5);
t811 = -Ifges(6,5) - Ifges(7,4);
t822 = Ifges(6,2) + Ifges(7,3);
t810 = Ifges(6,6) - Ifges(7,6);
t821 = -Ifges(6,3) - Ifges(7,2);
t768 = sin(pkin(10));
t770 = cos(pkin(10));
t751 = g(1) * t768 - g(2) * t770;
t752 = -g(1) * t770 - g(2) * t768;
t765 = -g(3) + qJDD(1);
t776 = cos(qJ(2));
t771 = cos(pkin(6));
t774 = sin(qJ(2));
t806 = t771 * t774;
t769 = sin(pkin(6));
t807 = t769 * t774;
t702 = t751 * t806 + t776 * t752 + t765 * t807;
t820 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t702;
t701 = -t774 * t752 + (t751 * t771 + t765 * t769) * t776;
t772 = sin(qJ(5));
t775 = cos(qJ(4));
t799 = qJD(2) * t775;
t817 = cos(qJ(5));
t745 = -qJD(4) * t817 + t772 * t799;
t773 = sin(qJ(4));
t798 = qJD(2) * qJD(4);
t795 = t773 * t798;
t750 = qJDD(2) * t775 - t795;
t714 = -t745 * qJD(5) + t772 * qJDD(4) + t750 * t817;
t746 = t772 * qJD(4) + t799 * t817;
t718 = mrSges(7,1) * t745 - mrSges(7,3) * t746;
t778 = qJD(2) ^ 2;
t783 = -t778 * qJ(3) + qJDD(3) - t701;
t818 = -pkin(2) - pkin(8);
t698 = qJDD(2) * t818 + t783;
t727 = -t751 * t769 + t765 * t771;
t694 = t773 * t698 + t775 * t727;
t748 = (pkin(4) * t773 - pkin(9) * t775) * qJD(2);
t777 = qJD(4) ^ 2;
t800 = qJD(2) * t773;
t689 = -pkin(4) * t777 + qJDD(4) * pkin(9) - t748 * t800 + t694;
t697 = t778 * t818 - t820;
t794 = t775 * t798;
t749 = -qJDD(2) * t773 - t794;
t691 = (-t750 + t795) * pkin(9) + (-t749 + t794) * pkin(4) + t697;
t685 = -t772 * t689 + t691 * t817;
t717 = pkin(5) * t745 - qJ(6) * t746;
t742 = qJDD(5) - t749;
t757 = qJD(5) + t800;
t756 = t757 ^ 2;
t683 = -t742 * pkin(5) - t756 * qJ(6) + t746 * t717 + qJDD(6) - t685;
t725 = -mrSges(7,2) * t745 + mrSges(7,3) * t757;
t788 = -m(7) * t683 + t742 * mrSges(7,1) + t757 * t725;
t679 = t714 * mrSges(7,2) + t746 * t718 - t788;
t686 = t817 * t689 + t772 * t691;
t682 = -pkin(5) * t756 + qJ(6) * t742 + 0.2e1 * qJD(6) * t757 - t717 * t745 + t686;
t713 = t746 * qJD(5) - qJDD(4) * t817 + t772 * t750;
t724 = -mrSges(7,1) * t757 + mrSges(7,2) * t746;
t796 = m(7) * t682 + t742 * mrSges(7,3) + t757 * t724;
t802 = t813 * t745 - t823 * t746 + t811 * t757;
t803 = t822 * t745 - t813 * t746 - t810 * t757;
t819 = -t713 * t810 - t714 * t811 - t821 * t742 - t745 * t802 - t803 * t746 + mrSges(6,1) * t685 - mrSges(7,1) * t683 - mrSges(6,2) * t686 + mrSges(7,3) * t682 - pkin(5) * t679 + qJ(6) * (-t713 * mrSges(7,2) - t745 * t718 + t796);
t816 = mrSges(3,1) - mrSges(4,2);
t815 = -mrSges(6,3) - mrSges(7,2);
t814 = (-Ifges(4,4) + Ifges(3,5));
t812 = Ifges(4,5) - Ifges(3,6);
t747 = (mrSges(5,1) * t773 + mrSges(5,2) * t775) * qJD(2);
t754 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t799;
t723 = mrSges(6,1) * t757 - mrSges(6,3) * t746;
t801 = -mrSges(6,1) * t745 - mrSges(6,2) * t746 - t718;
t674 = m(6) * t686 - t742 * mrSges(6,2) + t713 * t815 - t757 * t723 + t745 * t801 + t796;
t722 = -mrSges(6,2) * t757 - mrSges(6,3) * t745;
t676 = m(6) * t685 + t742 * mrSges(6,1) + t714 * t815 + t757 * t722 + t746 * t801 + t788;
t791 = t817 * t674 - t676 * t772;
t664 = m(5) * t694 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t749 - qJD(4) * t754 - t747 * t800 + t791;
t693 = t775 * t698 - t773 * t727;
t753 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t800;
t688 = -qJDD(4) * pkin(4) - t777 * pkin(9) + t748 * t799 - t693;
t684 = -0.2e1 * qJD(6) * t746 + (t745 * t757 - t714) * qJ(6) + (t746 * t757 + t713) * pkin(5) + t688;
t680 = m(7) * t684 + mrSges(7,1) * t713 - t714 * mrSges(7,3) - t746 * t724 + t725 * t745;
t779 = -m(6) * t688 - t713 * mrSges(6,1) - mrSges(6,2) * t714 - t745 * t722 - t723 * t746 - t680;
t671 = m(5) * t693 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t750 + qJD(4) * t753 - t747 * t799 + t779;
t658 = t773 * t664 + t775 * t671;
t700 = -qJDD(2) * pkin(2) + t783;
t785 = -m(4) * t700 + (t778 * mrSges(4,3)) - t658;
t653 = m(3) * t701 - (t778 * mrSges(3,2)) + qJDD(2) * t816 + t785;
t808 = t653 * t776;
t792 = t775 * t664 - t671 * t773;
t657 = m(4) * t727 + t792;
t656 = m(3) * t727 + t657;
t699 = t778 * pkin(2) + t820;
t670 = t772 * t674 + t817 * t676;
t784 = -m(5) * t697 + mrSges(5,1) * t749 - t750 * mrSges(5,2) - t753 * t800 - t754 * t799 - t670;
t781 = -m(4) * t699 + (t778 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t784;
t662 = m(3) * t702 - (mrSges(3,1) * t778) - qJDD(2) * mrSges(3,2) + t781;
t644 = -t656 * t769 + t662 * t806 + t771 * t808;
t642 = m(2) * t751 + t644;
t649 = -t653 * t774 + t776 * t662;
t648 = m(2) * t752 + t649;
t805 = t770 * t642 + t768 * t648;
t804 = t810 * t745 + t811 * t746 + t821 * t757;
t643 = t771 * t656 + t662 * t807 + t769 * t808;
t793 = -t642 * t768 + t770 * t648;
t789 = m(2) * t765 + t643;
t666 = -mrSges(6,1) * t688 - mrSges(7,1) * t684 + mrSges(7,2) * t682 + mrSges(6,3) * t686 - pkin(5) * t680 - t822 * t713 + t813 * t714 + t810 * t742 + t804 * t746 - t802 * t757;
t668 = mrSges(6,2) * t688 + mrSges(7,2) * t683 - mrSges(6,3) * t685 - mrSges(7,3) * t684 - qJ(6) * t680 - t813 * t713 + t823 * t714 - t811 * t742 + t804 * t745 + t803 * t757;
t732 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t775 - Ifges(5,6) * t773) * qJD(2);
t733 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t775 - Ifges(5,2) * t773) * qJD(2);
t645 = mrSges(5,2) * t697 - mrSges(5,3) * t693 + Ifges(5,1) * t750 + Ifges(5,4) * t749 + Ifges(5,5) * qJDD(4) - pkin(9) * t670 - qJD(4) * t733 - t772 * t666 + t668 * t817 - t732 * t800;
t734 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t775 - Ifges(5,4) * t773) * qJD(2);
t650 = -mrSges(5,1) * t697 + mrSges(5,3) * t694 + Ifges(5,4) * t750 + Ifges(5,2) * t749 + Ifges(5,6) * qJDD(4) - pkin(4) * t670 + qJD(4) * t734 - t732 * t799 - t819;
t639 = -mrSges(4,1) * t699 + mrSges(3,3) * t702 - pkin(2) * t657 - pkin(3) * t784 - pkin(8) * t792 - qJDD(2) * t812 - t773 * t645 - t775 * t650 - t727 * t816 + (t778 * t814);
t782 = mrSges(5,1) * t693 - mrSges(5,2) * t694 + Ifges(5,5) * t750 + Ifges(5,6) * t749 + Ifges(5,3) * qJDD(4) + pkin(4) * t779 + pkin(9) * t791 + t817 * t666 + t772 * t668 + t733 * t799 + t734 * t800;
t640 = mrSges(4,1) * t700 - mrSges(3,3) * t701 + (mrSges(3,2) - mrSges(4,3)) * t727 + t812 * t778 + t782 - qJ(3) * t657 + pkin(3) * t658 + t814 * qJDD(2);
t786 = pkin(7) * t649 + t639 * t776 + t640 * t774;
t654 = qJDD(2) * mrSges(4,2) - t785;
t638 = mrSges(3,1) * t701 - mrSges(3,2) * t702 + mrSges(4,2) * t700 - mrSges(4,3) * t699 + t775 * t645 - t773 * t650 - pkin(8) * t658 - pkin(2) * t654 + qJ(3) * t781 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2);
t637 = mrSges(2,2) * t765 - mrSges(2,3) * t751 - t774 * t639 + t776 * t640 + (-t643 * t769 - t644 * t771) * pkin(7);
t636 = -mrSges(2,1) * t765 + mrSges(2,3) * t752 - pkin(1) * t643 - t769 * t638 + t786 * t771;
t1 = [-m(1) * g(1) + t793; -m(1) * g(2) + t805; -m(1) * g(3) + t789; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t805 - t768 * t636 + t770 * t637; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t793 + t770 * t636 + t768 * t637; -mrSges(1,1) * g(2) + mrSges(2,1) * t751 + mrSges(1,2) * g(1) - mrSges(2,2) * t752 + pkin(1) * t644 + t771 * t638 + t769 * t786; t789; t638; t654; t782; t819; t679;];
tauJB  = t1;
