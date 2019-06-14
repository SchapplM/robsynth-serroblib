% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-05-05 15:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:10:32
% EndTime: 2019-05-05 15:10:45
% DurationCPUTime: 13.44s
% Computational Cost: add. (221441->323), mult. (503067->401), div. (0->0), fcn. (371304->12), ass. (0->141)
t865 = qJD(1) ^ 2;
t855 = cos(pkin(11));
t895 = pkin(3) * t855;
t853 = sin(pkin(11));
t894 = mrSges(4,2) * t853;
t848 = t855 ^ 2;
t893 = t848 * t865;
t860 = sin(qJ(1));
t864 = cos(qJ(1));
t833 = t860 * g(1) - g(2) * t864;
t830 = qJDD(1) * pkin(1) + t833;
t834 = -g(1) * t864 - g(2) * t860;
t831 = -pkin(1) * t865 + t834;
t854 = sin(pkin(10));
t856 = cos(pkin(10));
t818 = t854 * t830 + t856 * t831;
t809 = -pkin(2) * t865 + qJDD(1) * qJ(3) + t818;
t852 = -g(3) + qJDD(2);
t888 = qJD(1) * qJD(3);
t891 = t855 * t852 - 0.2e1 * t853 * t888;
t792 = (-pkin(7) * qJDD(1) + t865 * t895 - t809) * t853 + t891;
t796 = t853 * t852 + (t809 + 0.2e1 * t888) * t855;
t887 = qJDD(1) * t855;
t793 = -pkin(3) * t893 + pkin(7) * t887 + t796;
t859 = sin(qJ(4));
t863 = cos(qJ(4));
t768 = t863 * t792 - t859 * t793;
t876 = t853 * t863 + t855 * t859;
t875 = -t853 * t859 + t855 * t863;
t823 = t875 * qJD(1);
t889 = t823 * qJD(4);
t815 = qJDD(1) * t876 + t889;
t824 = t876 * qJD(1);
t757 = (-t815 + t889) * pkin(8) + (t823 * t824 + qJDD(4)) * pkin(4) + t768;
t769 = t859 * t792 + t863 * t793;
t814 = -t824 * qJD(4) + qJDD(1) * t875;
t821 = qJD(4) * pkin(4) - pkin(8) * t824;
t822 = t823 ^ 2;
t759 = -pkin(4) * t822 + pkin(8) * t814 - qJD(4) * t821 + t769;
t858 = sin(qJ(5));
t862 = cos(qJ(5));
t755 = t858 * t757 + t862 * t759;
t808 = t823 * t858 + t824 * t862;
t777 = -qJD(5) * t808 + t814 * t862 - t815 * t858;
t807 = t823 * t862 - t824 * t858;
t787 = -mrSges(6,1) * t807 + mrSges(6,2) * t808;
t849 = qJD(4) + qJD(5);
t800 = mrSges(6,1) * t849 - mrSges(6,3) * t808;
t846 = qJDD(4) + qJDD(5);
t788 = -pkin(5) * t807 - pkin(9) * t808;
t845 = t849 ^ 2;
t751 = -pkin(5) * t845 + pkin(9) * t846 + t788 * t807 + t755;
t847 = t853 ^ 2;
t817 = t856 * t830 - t854 * t831;
t877 = qJDD(3) - t817;
t794 = (-pkin(2) - t895) * qJDD(1) + (-qJ(3) + (-t847 - t848) * pkin(7)) * t865 + t877;
t764 = -t814 * pkin(4) - t822 * pkin(8) + t824 * t821 + t794;
t778 = qJD(5) * t807 + t814 * t858 + t815 * t862;
t752 = (-t807 * t849 - t778) * pkin(9) + (t808 * t849 - t777) * pkin(5) + t764;
t857 = sin(qJ(6));
t861 = cos(qJ(6));
t748 = -t751 * t857 + t752 * t861;
t797 = -t808 * t857 + t849 * t861;
t762 = qJD(6) * t797 + t778 * t861 + t846 * t857;
t776 = qJDD(6) - t777;
t798 = t808 * t861 + t849 * t857;
t779 = -mrSges(7,1) * t797 + mrSges(7,2) * t798;
t802 = qJD(6) - t807;
t780 = -mrSges(7,2) * t802 + mrSges(7,3) * t797;
t744 = m(7) * t748 + mrSges(7,1) * t776 - mrSges(7,3) * t762 - t779 * t798 + t780 * t802;
t749 = t751 * t861 + t752 * t857;
t761 = -qJD(6) * t798 - t778 * t857 + t846 * t861;
t781 = mrSges(7,1) * t802 - mrSges(7,3) * t798;
t745 = m(7) * t749 - mrSges(7,2) * t776 + mrSges(7,3) * t761 + t779 * t797 - t781 * t802;
t881 = -t744 * t857 + t861 * t745;
t731 = m(6) * t755 - mrSges(6,2) * t846 + mrSges(6,3) * t777 + t787 * t807 - t800 * t849 + t881;
t754 = t757 * t862 - t759 * t858;
t799 = -mrSges(6,2) * t849 + mrSges(6,3) * t807;
t750 = -pkin(5) * t846 - pkin(9) * t845 + t788 * t808 - t754;
t872 = -m(7) * t750 + t761 * mrSges(7,1) - mrSges(7,2) * t762 + t797 * t780 - t781 * t798;
t740 = m(6) * t754 + mrSges(6,1) * t846 - mrSges(6,3) * t778 - t787 * t808 + t799 * t849 + t872;
t724 = t858 * t731 + t862 * t740;
t812 = -mrSges(5,1) * t823 + mrSges(5,2) * t824;
t819 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t823;
t722 = m(5) * t768 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t815 + qJD(4) * t819 - t812 * t824 + t724;
t820 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t824;
t882 = t862 * t731 - t740 * t858;
t723 = m(5) * t769 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t814 - qJD(4) * t820 + t812 * t823 + t882;
t716 = t863 * t722 + t859 * t723;
t795 = -t809 * t853 + t891;
t874 = mrSges(4,3) * qJDD(1) + t865 * (-mrSges(4,1) * t855 + t894);
t714 = m(4) * t795 - t853 * t874 + t716;
t883 = -t859 * t722 + t863 * t723;
t715 = m(4) * t796 + t855 * t874 + t883;
t884 = -t714 * t853 + t855 * t715;
t706 = m(3) * t818 - mrSges(3,1) * t865 - qJDD(1) * mrSges(3,2) + t884;
t803 = -qJDD(1) * pkin(2) - t865 * qJ(3) + t877;
t733 = t861 * t744 + t857 * t745;
t873 = m(6) * t764 - t777 * mrSges(6,1) + t778 * mrSges(6,2) - t807 * t799 + t808 * t800 + t733;
t869 = m(5) * t794 - t814 * mrSges(5,1) + t815 * mrSges(5,2) - t823 * t819 + t824 * t820 + t873;
t867 = -m(4) * t803 + mrSges(4,1) * t887 - t869 + (t847 * t865 + t893) * mrSges(4,3);
t726 = t867 + (mrSges(3,1) - t894) * qJDD(1) - t865 * mrSges(3,2) + m(3) * t817;
t702 = t854 * t706 + t856 * t726;
t699 = m(2) * t833 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t865 + t702;
t885 = t856 * t706 - t726 * t854;
t700 = m(2) * t834 - mrSges(2,1) * t865 - qJDD(1) * mrSges(2,2) + t885;
t892 = t864 * t699 + t860 * t700;
t709 = t855 * t714 + t853 * t715;
t878 = Ifges(4,5) * t853 + Ifges(4,6) * t855;
t890 = t865 * t878;
t707 = m(3) * t852 + t709;
t886 = -t699 * t860 + t864 * t700;
t880 = Ifges(4,1) * t853 + Ifges(4,4) * t855;
t879 = Ifges(4,4) * t853 + Ifges(4,2) * t855;
t765 = Ifges(7,5) * t798 + Ifges(7,6) * t797 + Ifges(7,3) * t802;
t767 = Ifges(7,1) * t798 + Ifges(7,4) * t797 + Ifges(7,5) * t802;
t737 = -mrSges(7,1) * t750 + mrSges(7,3) * t749 + Ifges(7,4) * t762 + Ifges(7,2) * t761 + Ifges(7,6) * t776 - t765 * t798 + t767 * t802;
t766 = Ifges(7,4) * t798 + Ifges(7,2) * t797 + Ifges(7,6) * t802;
t738 = mrSges(7,2) * t750 - mrSges(7,3) * t748 + Ifges(7,1) * t762 + Ifges(7,4) * t761 + Ifges(7,5) * t776 + t765 * t797 - t766 * t802;
t782 = Ifges(6,5) * t808 + Ifges(6,6) * t807 + Ifges(6,3) * t849;
t783 = Ifges(6,4) * t808 + Ifges(6,2) * t807 + Ifges(6,6) * t849;
t717 = mrSges(6,2) * t764 - mrSges(6,3) * t754 + Ifges(6,1) * t778 + Ifges(6,4) * t777 + Ifges(6,5) * t846 - pkin(9) * t733 - t737 * t857 + t738 * t861 + t782 * t807 - t783 * t849;
t784 = Ifges(6,1) * t808 + Ifges(6,4) * t807 + Ifges(6,5) * t849;
t868 = mrSges(7,1) * t748 - mrSges(7,2) * t749 + Ifges(7,5) * t762 + Ifges(7,6) * t761 + Ifges(7,3) * t776 + t766 * t798 - t767 * t797;
t718 = -mrSges(6,1) * t764 + mrSges(6,3) * t755 + Ifges(6,4) * t778 + Ifges(6,2) * t777 + Ifges(6,6) * t846 - pkin(5) * t733 - t782 * t808 + t784 * t849 - t868;
t804 = Ifges(5,5) * t824 + Ifges(5,6) * t823 + Ifges(5,3) * qJD(4);
t806 = Ifges(5,1) * t824 + Ifges(5,4) * t823 + Ifges(5,5) * qJD(4);
t703 = -mrSges(5,1) * t794 + mrSges(5,3) * t769 + Ifges(5,4) * t815 + Ifges(5,2) * t814 + Ifges(5,6) * qJDD(4) - pkin(4) * t873 + pkin(8) * t882 + qJD(4) * t806 + t858 * t717 + t862 * t718 - t824 * t804;
t805 = Ifges(5,4) * t824 + Ifges(5,2) * t823 + Ifges(5,6) * qJD(4);
t710 = mrSges(5,2) * t794 - mrSges(5,3) * t768 + Ifges(5,1) * t815 + Ifges(5,4) * t814 + Ifges(5,5) * qJDD(4) - pkin(8) * t724 - qJD(4) * t805 + t717 * t862 - t718 * t858 + t804 * t823;
t692 = -mrSges(4,1) * t803 + mrSges(4,3) * t796 - pkin(3) * t869 + pkin(7) * t883 + qJDD(1) * t879 + t863 * t703 + t859 * t710 - t853 * t890;
t694 = mrSges(4,2) * t803 - mrSges(4,3) * t795 - pkin(7) * t716 + qJDD(1) * t880 - t859 * t703 + t863 * t710 + t855 * t890;
t728 = qJDD(1) * t894 - t867;
t871 = mrSges(2,1) * t833 + mrSges(3,1) * t817 - mrSges(2,2) * t834 - mrSges(3,2) * t818 + pkin(1) * t702 - pkin(2) * t728 + qJ(3) * t884 + t855 * t692 + t853 * t694 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t870 = -mrSges(6,1) * t754 + mrSges(6,2) * t755 - Ifges(6,5) * t778 - Ifges(6,6) * t777 - Ifges(6,3) * t846 - pkin(5) * t872 - pkin(9) * t881 - t861 * t737 - t857 * t738 - t808 * t783 + t807 * t784;
t866 = mrSges(5,1) * t768 - mrSges(5,2) * t769 + Ifges(5,5) * t815 + Ifges(5,6) * t814 + Ifges(5,3) * qJDD(4) + pkin(4) * t724 + t824 * t805 - t823 * t806 - t870;
t695 = -t866 + (Ifges(3,6) - t878) * qJDD(1) - mrSges(3,1) * t852 + mrSges(3,3) * t818 - mrSges(4,1) * t795 + mrSges(4,2) * t796 - pkin(3) * t716 - pkin(2) * t709 + (-t853 * t879 + t855 * t880 + Ifges(3,5)) * t865;
t690 = mrSges(3,2) * t852 - mrSges(3,3) * t817 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t865 - qJ(3) * t709 - t692 * t853 + t694 * t855;
t689 = -mrSges(2,2) * g(3) - mrSges(2,3) * t833 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t865 - qJ(2) * t702 + t690 * t856 - t695 * t854;
t688 = mrSges(2,1) * g(3) + mrSges(2,3) * t834 + t865 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t707 + qJ(2) * t885 + t854 * t690 + t856 * t695;
t1 = [-m(1) * g(1) + t886; -m(1) * g(2) + t892; (-m(1) - m(2)) * g(3) + t707; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t892 - t860 * t688 + t864 * t689; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t886 + t864 * t688 + t860 * t689; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t871; t871; t707; t728; t866; -t870; t868;];
tauJB  = t1;
