% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPPRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:31:13
% EndTime: 2020-01-03 11:31:20
% DurationCPUTime: 8.02s
% Computational Cost: add. (73722->278), mult. (201875->377), div. (0->0), fcn. (138541->10), ass. (0->134)
t797 = sin(qJ(1));
t800 = cos(qJ(1));
t773 = -g(2) * t797 + t800 * g(3);
t801 = qJD(1) ^ 2;
t845 = -pkin(1) * t801 + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t773;
t774 = -t800 * g(2) - t797 * g(3);
t812 = -t801 * qJ(2) + qJDD(2) - t774;
t792 = sin(pkin(8));
t794 = cos(pkin(8));
t818 = -pkin(2) * t794 - qJ(3) * t792;
t836 = t792 * qJD(1);
t844 = (-pkin(1) + t818) * qJDD(1) + t812 - 0.2e1 * qJD(3) * t836;
t746 = -t794 * g(1) - t792 * t845;
t843 = mrSges(3,2) * t792;
t789 = t792 ^ 2;
t842 = t789 * t801;
t791 = sin(pkin(9));
t841 = t791 * t792;
t793 = cos(pkin(9));
t840 = t792 * t793;
t747 = -t792 * g(1) + t794 * t845;
t766 = (-mrSges(3,1) * t794 + t843) * qJD(1);
t765 = t818 * qJD(1);
t835 = t794 * qJD(1);
t736 = t765 * t835 + t747;
t815 = -pkin(3) * t794 - pkin(6) * t840;
t838 = t844 * t793;
t715 = t815 * qJDD(1) + (-t736 + (-pkin(3) * t789 * t793 + pkin(6) * t792 * t794) * t801) * t791 + t838;
t724 = t793 * t736 + t791 * t844;
t764 = t815 * qJD(1);
t833 = qJDD(1) * t792;
t829 = t791 * t833;
t831 = t791 ^ 2 * t842;
t716 = -pkin(3) * t831 - pkin(6) * t829 + t764 * t835 + t724;
t796 = sin(qJ(4));
t799 = cos(qJ(4));
t701 = t715 * t799 - t796 * t716;
t811 = (-t791 * t799 - t793 * t796) * t792;
t754 = qJD(1) * t811;
t810 = (-t791 * t796 + t793 * t799) * t792;
t740 = t754 * qJD(4) + qJDD(1) * t810;
t755 = qJD(1) * t810;
t832 = t794 * qJDD(1);
t777 = qJDD(4) - t832;
t778 = qJD(4) - t835;
t699 = (t754 * t778 - t740) * pkin(7) + (t754 * t755 + t777) * pkin(4) + t701;
t702 = t715 * t796 + t716 * t799;
t739 = -t755 * qJD(4) + qJDD(1) * t811;
t745 = pkin(4) * t778 - pkin(7) * t755;
t753 = t754 ^ 2;
t700 = -pkin(4) * t753 + pkin(7) * t739 - t745 * t778 + t702;
t795 = sin(qJ(5));
t798 = cos(qJ(5));
t697 = t699 * t798 - t700 * t795;
t733 = t754 * t798 - t755 * t795;
t711 = qJD(5) * t733 + t739 * t795 + t740 * t798;
t734 = t754 * t795 + t755 * t798;
t722 = -mrSges(6,1) * t733 + mrSges(6,2) * t734;
t776 = qJD(5) + t778;
t726 = -mrSges(6,2) * t776 + mrSges(6,3) * t733;
t772 = qJDD(5) + t777;
t693 = m(6) * t697 + mrSges(6,1) * t772 - mrSges(6,3) * t711 - t722 * t734 + t726 * t776;
t698 = t699 * t795 + t700 * t798;
t710 = -qJD(5) * t734 + t739 * t798 - t740 * t795;
t727 = mrSges(6,1) * t776 - mrSges(6,3) * t734;
t694 = m(6) * t698 - mrSges(6,2) * t772 + mrSges(6,3) * t710 + t722 * t733 - t727 * t776;
t685 = t693 * t798 + t694 * t795;
t737 = -mrSges(5,1) * t754 + mrSges(5,2) * t755;
t741 = -mrSges(5,2) * t778 + mrSges(5,3) * t754;
t683 = m(5) * t701 + mrSges(5,1) * t777 - mrSges(5,3) * t740 - t737 * t755 + t741 * t778 + t685;
t742 = mrSges(5,1) * t778 - mrSges(5,3) * t755;
t824 = -t693 * t795 + t694 * t798;
t684 = m(5) * t702 - mrSges(5,2) * t777 + mrSges(5,3) * t739 + t737 * t754 - t742 * t778 + t824;
t679 = t683 * t799 + t684 * t796;
t723 = -t791 * t736 + t838;
t822 = mrSges(4,1) * t791 + mrSges(4,2) * t793;
t758 = t822 * t836;
t813 = mrSges(4,2) * t794 - mrSges(4,3) * t841;
t761 = t813 * qJD(1);
t814 = -mrSges(4,1) * t794 - mrSges(4,3) * t840;
t677 = m(4) * t723 + t814 * qJDD(1) + (-t758 * t840 - t761 * t794) * qJD(1) + t679;
t762 = t814 * qJD(1);
t825 = -t796 * t683 + t684 * t799;
t678 = m(4) * t724 + t813 * qJDD(1) + (-t758 * t841 + t762 * t794) * qJD(1) + t825;
t826 = -t791 * t677 + t678 * t793;
t670 = m(3) * t747 + (qJDD(1) * mrSges(3,3) + qJD(1) * t766) * t794 + t826;
t735 = t765 * t836 + qJDD(3) - t746;
t725 = t764 * t793 * t836 + pkin(3) * t829 - pkin(6) * t831 + t735;
t704 = -pkin(4) * t739 - pkin(7) * t753 + t745 * t755 + t725;
t817 = m(6) * t704 - mrSges(6,1) * t710 + mrSges(6,2) * t711 - t726 * t733 + t727 * t734;
t804 = m(5) * t725 - mrSges(5,1) * t739 + t740 * mrSges(5,2) - t741 * t754 + t755 * t742 + t817;
t803 = m(4) * t735 + t804;
t816 = t761 * t791 + t762 * t793;
t689 = -t803 + m(3) * t746 + ((-mrSges(3,3) - t822) * qJDD(1) + (-t766 - t816) * qJD(1)) * t792;
t827 = t670 * t794 - t689 * t792;
t662 = m(2) * t773 - mrSges(2,1) * t801 - qJDD(1) * mrSges(2,2) + t827;
t673 = t793 * t677 + t791 * t678;
t760 = -qJDD(1) * pkin(1) + t812;
t805 = -m(3) * t760 + mrSges(3,1) * t832 - t673 + (t794 ^ 2 * t801 + t842) * mrSges(3,3);
t667 = m(2) * t774 - t801 * mrSges(2,2) + (mrSges(2,1) - t843) * qJDD(1) + t805;
t839 = t662 * t797 + t667 * t800;
t664 = t670 * t792 + t689 * t794;
t828 = -t662 * t800 + t667 * t797;
t821 = Ifges(3,1) * t792 + Ifges(3,4) * t794;
t820 = Ifges(3,5) * t792 + Ifges(3,6) * t794;
t819 = Ifges(4,5) * t793 - Ifges(4,6) * t791;
t717 = Ifges(6,5) * t734 + Ifges(6,6) * t733 + Ifges(6,3) * t776;
t719 = Ifges(6,1) * t734 + Ifges(6,4) * t733 + Ifges(6,5) * t776;
t686 = -mrSges(6,1) * t704 + mrSges(6,3) * t698 + Ifges(6,4) * t711 + Ifges(6,2) * t710 + Ifges(6,6) * t772 - t717 * t734 + t719 * t776;
t718 = Ifges(6,4) * t734 + Ifges(6,2) * t733 + Ifges(6,6) * t776;
t687 = mrSges(6,2) * t704 - mrSges(6,3) * t697 + Ifges(6,1) * t711 + Ifges(6,4) * t710 + Ifges(6,5) * t772 + t717 * t733 - t718 * t776;
t728 = Ifges(5,5) * t755 + Ifges(5,6) * t754 + Ifges(5,3) * t778;
t730 = Ifges(5,1) * t755 + Ifges(5,4) * t754 + Ifges(5,5) * t778;
t674 = -mrSges(5,1) * t725 + mrSges(5,3) * t702 + Ifges(5,4) * t740 + Ifges(5,2) * t739 + Ifges(5,6) * t777 - pkin(4) * t817 + pkin(7) * t824 + t798 * t686 + t795 * t687 - t755 * t728 + t778 * t730;
t729 = Ifges(5,4) * t755 + Ifges(5,2) * t754 + Ifges(5,6) * t778;
t675 = mrSges(5,2) * t725 - mrSges(5,3) * t701 + Ifges(5,1) * t740 + Ifges(5,4) * t739 + Ifges(5,5) * t777 - pkin(7) * t685 - t686 * t795 + t687 * t798 + t728 * t754 - t729 * t778;
t749 = (-Ifges(4,3) * t794 + t792 * t819) * qJD(1);
t808 = -Ifges(4,5) * t794 + (Ifges(4,1) * t793 - Ifges(4,4) * t791) * t792;
t751 = t808 * qJD(1);
t807 = -Ifges(4,6) * t794 + (Ifges(4,4) * t793 - Ifges(4,2) * t791) * t792;
t659 = -mrSges(4,1) * t735 + mrSges(4,3) * t724 + t796 * t675 + t799 * t674 - pkin(3) * t804 + pkin(6) * t825 + (-t749 * t840 - t751 * t794) * qJD(1) + t807 * qJDD(1);
t750 = t807 * qJD(1);
t660 = mrSges(4,2) * t735 - mrSges(4,3) * t723 - pkin(6) * t679 - t796 * t674 + t799 * t675 + (-t749 * t841 + t750 * t794) * qJD(1) + t808 * qJDD(1);
t767 = t820 * qJD(1);
t656 = mrSges(3,2) * t760 - mrSges(3,3) * t746 - qJ(3) * t673 + qJDD(1) * t821 - t791 * t659 + t793 * t660 + t767 * t835;
t806 = -mrSges(6,1) * t697 + mrSges(6,2) * t698 - Ifges(6,5) * t711 - Ifges(6,6) * t710 - Ifges(6,3) * t772 - t718 * t734 + t733 * t719;
t802 = mrSges(5,1) * t701 - mrSges(5,2) * t702 + Ifges(5,5) * t740 + Ifges(5,6) * t739 + Ifges(5,3) * t777 + pkin(4) * t685 + t755 * t729 - t754 * t730 - t806;
t658 = (Ifges(3,2) + Ifges(4,3)) * t832 - t802 - mrSges(4,1) * t723 + mrSges(4,2) * t724 - pkin(3) * t679 - pkin(2) * t673 + ((Ifges(3,4) - t819) * qJDD(1) + (-t750 * t793 - t751 * t791 - t767) * qJD(1)) * t792 + mrSges(3,3) * t747 - mrSges(3,1) * t760;
t672 = mrSges(3,2) * t833 - t805;
t809 = mrSges(2,1) * t774 - mrSges(2,2) * t773 + Ifges(2,3) * qJDD(1) - pkin(1) * t672 + qJ(2) * t827 + t656 * t792 + t658 * t794;
t695 = (qJD(1) * t816 + qJDD(1) * t822) * t792 + t803;
t654 = mrSges(2,1) * g(1) + mrSges(2,3) * t773 - mrSges(3,1) * t746 + mrSges(3,2) * t747 - t791 * t660 - t793 * t659 + pkin(2) * t695 - qJ(3) * t826 - pkin(1) * t664 + (Ifges(2,6) - t820) * qJDD(1) + (Ifges(2,5) - t792 * (Ifges(3,4) * t792 + Ifges(3,2) * t794) + t794 * t821) * t801;
t653 = -mrSges(2,2) * g(1) - mrSges(2,3) * t774 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t801 - qJ(2) * t664 + t656 * t794 - t658 * t792;
t1 = [(-m(1) - m(2)) * g(1) + t664; -m(1) * g(2) + t839; -m(1) * g(3) + t828; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t809; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t828 + t797 * t653 + t800 * t654; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t839 - t653 * t800 + t654 * t797; t809; t672; t695; t802; -t806;];
tauJB = t1;
