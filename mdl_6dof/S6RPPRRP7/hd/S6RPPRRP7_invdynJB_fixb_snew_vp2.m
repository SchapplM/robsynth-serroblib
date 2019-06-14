% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRP7
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
% Datum: 2019-05-05 15:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRP7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:03:40
% EndTime: 2019-05-05 15:03:45
% DurationCPUTime: 4.88s
% Computational Cost: add. (53402->292), mult. (119850->339), div. (0->0), fcn. (81411->8), ass. (0->129)
t901 = Ifges(6,1) + Ifges(7,1);
t889 = Ifges(6,4) + Ifges(7,4);
t888 = Ifges(6,5) + Ifges(7,5);
t900 = Ifges(6,2) + Ifges(7,2);
t886 = Ifges(6,6) + Ifges(7,6);
t899 = Ifges(6,3) + Ifges(7,3);
t841 = sin(qJ(1));
t844 = cos(qJ(1));
t818 = t841 * g(1) - t844 * g(2);
t846 = qJD(1) ^ 2;
t855 = -t846 * qJ(2) + qJDD(2) - t818;
t884 = -pkin(1) - qJ(3);
t898 = -(2 * qJD(1) * qJD(3)) + qJDD(1) * t884 + t855;
t819 = -t844 * g(1) - t841 * g(2);
t897 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t819;
t837 = sin(pkin(9));
t838 = cos(pkin(9));
t840 = sin(qJ(4));
t843 = cos(qJ(4));
t859 = t837 * t843 + t838 * t840;
t813 = t859 * qJD(1);
t809 = t846 * pkin(1) - t897;
t896 = -m(3) * t809 + t846 * mrSges(3,2) + qJDD(1) * mrSges(3,3);
t858 = -t837 * t840 + t838 * t843;
t814 = t858 * qJD(1);
t874 = t814 * qJD(4);
t797 = -qJDD(1) * t859 - t874;
t875 = t813 * qJD(4);
t798 = qJDD(1) * t858 - t875;
t839 = sin(qJ(5));
t842 = cos(qJ(5));
t801 = qJD(4) * t842 - t814 * t839;
t766 = qJD(5) * t801 + qJDD(4) * t839 + t798 * t842;
t802 = qJD(4) * t839 + t814 * t842;
t770 = -mrSges(7,1) * t801 + mrSges(7,2) * t802;
t794 = t837 * g(3) + t898 * t838;
t893 = pkin(3) * t846;
t775 = (-pkin(7) * qJDD(1) - t837 * t893) * t838 + t794;
t795 = -g(3) * t838 + t898 * t837;
t828 = t837 ^ 2;
t872 = qJDD(1) * t837;
t776 = -pkin(7) * t872 - t828 * t893 + t795;
t752 = t840 * t775 + t843 * t776;
t796 = pkin(4) * t813 - pkin(8) * t814;
t845 = qJD(4) ^ 2;
t746 = -pkin(4) * t845 + qJDD(4) * pkin(8) - t796 * t813 + t752;
t854 = qJDD(3) + t897;
t877 = -t838 ^ 2 - t828;
t783 = pkin(3) * t872 + (pkin(7) * t877 + t884) * t846 + t854;
t749 = (-t798 + t875) * pkin(8) + (-t797 + t874) * pkin(4) + t783;
t742 = -t839 * t746 + t842 * t749;
t793 = qJDD(5) - t797;
t811 = qJD(5) + t813;
t738 = -0.2e1 * qJD(6) * t802 + (t801 * t811 - t766) * qJ(6) + (t801 * t802 + t793) * pkin(5) + t742;
t777 = -mrSges(7,2) * t811 + mrSges(7,3) * t801;
t870 = m(7) * t738 + t793 * mrSges(7,1) + t811 * t777;
t735 = -t766 * mrSges(7,3) - t802 * t770 + t870;
t743 = t842 * t746 + t839 * t749;
t765 = -qJD(5) * t802 + qJDD(4) * t842 - t798 * t839;
t779 = pkin(5) * t811 - qJ(6) * t802;
t800 = t801 ^ 2;
t740 = -pkin(5) * t800 + t765 * qJ(6) + 0.2e1 * qJD(6) * t801 - t779 * t811 + t743;
t879 = t889 * t801 + t901 * t802 + t888 * t811;
t880 = -t900 * t801 - t889 * t802 - t886 * t811;
t895 = mrSges(6,1) * t742 + mrSges(7,1) * t738 - mrSges(6,2) * t743 - mrSges(7,2) * t740 + pkin(5) * t735 + t765 * t886 + t766 * t888 + t899 * t793 - t801 * t879 - t802 * t880;
t892 = mrSges(2,1) - mrSges(3,2);
t891 = -mrSges(6,2) - mrSges(7,2);
t890 = -Ifges(3,4) + Ifges(2,5);
t887 = -Ifges(2,6) + Ifges(3,5);
t883 = mrSges(4,2) * t838;
t789 = mrSges(5,1) * t813 + mrSges(5,2) * t814;
t807 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t814;
t771 = -mrSges(6,1) * t801 + mrSges(6,2) * t802;
t778 = -mrSges(6,2) * t811 + mrSges(6,3) * t801;
t728 = m(6) * t742 + t793 * mrSges(6,1) + t811 * t778 + (-t770 - t771) * t802 + (-mrSges(6,3) - mrSges(7,3)) * t766 + t870;
t869 = m(7) * t740 + t765 * mrSges(7,3) + t801 * t770;
t780 = mrSges(7,1) * t811 - mrSges(7,3) * t802;
t878 = -mrSges(6,1) * t811 + mrSges(6,3) * t802 - t780;
t731 = m(6) * t743 + t765 * mrSges(6,3) + t801 * t771 + t793 * t891 + t811 * t878 + t869;
t864 = -t728 * t839 + t842 * t731;
t722 = m(5) * t752 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t797 - qJD(4) * t807 - t789 * t813 + t864;
t751 = t775 * t843 - t840 * t776;
t806 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t813;
t745 = -qJDD(4) * pkin(4) - pkin(8) * t845 + t814 * t796 - t751;
t741 = -t765 * pkin(5) - qJ(6) * t800 + t779 * t802 + qJDD(6) + t745;
t863 = -m(7) * t741 + t765 * mrSges(7,1) + t801 * t777;
t848 = -m(6) * t745 + t765 * mrSges(6,1) + t766 * t891 + t801 * t778 + t802 * t878 + t863;
t733 = m(5) * t751 + qJDD(4) * mrSges(5,1) - t798 * mrSges(5,3) + qJD(4) * t806 - t814 * t789 + t848;
t712 = t840 * t722 + t843 * t733;
t857 = -qJDD(1) * mrSges(4,3) - t846 * (mrSges(4,1) * t837 + t883);
t710 = m(4) * t794 + t838 * t857 + t712;
t865 = t843 * t722 - t840 * t733;
t711 = m(4) * t795 + t837 * t857 + t865;
t707 = t838 * t710 + t837 * t711;
t812 = -qJDD(1) * pkin(1) + t855;
t853 = -m(3) * t812 + t846 * mrSges(3,3) - t707;
t703 = m(2) * t818 - t846 * mrSges(2,2) + qJDD(1) * t892 + t853;
t805 = t846 * t884 + t854;
t726 = t842 * t728 + t839 * t731;
t852 = m(5) * t783 - t797 * mrSges(5,1) + t798 * mrSges(5,2) + t813 * t806 + t814 * t807 + t726;
t850 = m(4) * t805 + mrSges(4,1) * t872 + qJDD(1) * t883 + t852;
t868 = t877 * mrSges(4,3);
t717 = (-mrSges(2,1) + t868) * t846 + m(2) * t819 - qJDD(1) * mrSges(2,2) + t850 + t896;
t882 = t844 * t703 + t841 * t717;
t881 = -t886 * t801 - t888 * t802 - t899 * t811;
t860 = Ifges(4,5) * t838 - Ifges(4,6) * t837;
t876 = t846 * t860;
t867 = -t703 * t841 + t844 * t717;
t866 = -t837 * t710 + t838 * t711;
t862 = Ifges(4,1) * t838 - Ifges(4,4) * t837;
t861 = Ifges(4,4) * t838 - Ifges(4,2) * t837;
t736 = t766 * mrSges(7,2) + t802 * t780 - t863;
t714 = -mrSges(6,1) * t745 + mrSges(6,3) * t743 - mrSges(7,1) * t741 + mrSges(7,3) * t740 - pkin(5) * t736 + qJ(6) * t869 + (-qJ(6) * t780 + t879) * t811 + t881 * t802 + (-mrSges(7,2) * qJ(6) + t886) * t793 + t889 * t766 + t900 * t765;
t724 = mrSges(6,2) * t745 + mrSges(7,2) * t741 - mrSges(6,3) * t742 - mrSges(7,3) * t738 - qJ(6) * t735 + t889 * t765 + t901 * t766 + t888 * t793 - t881 * t801 + t880 * t811;
t785 = Ifges(5,4) * t814 - Ifges(5,2) * t813 + Ifges(5,6) * qJD(4);
t786 = Ifges(5,1) * t814 - Ifges(5,4) * t813 + Ifges(5,5) * qJD(4);
t849 = mrSges(5,1) * t751 - mrSges(5,2) * t752 + Ifges(5,5) * t798 + Ifges(5,6) * t797 + Ifges(5,3) * qJDD(4) + pkin(4) * t848 + pkin(8) * t864 + t842 * t714 + t839 * t724 + t814 * t785 + t813 * t786;
t784 = Ifges(5,5) * t814 - Ifges(5,6) * t813 + Ifges(5,3) * qJD(4);
t701 = mrSges(5,2) * t783 - mrSges(5,3) * t751 + Ifges(5,1) * t798 + Ifges(5,4) * t797 + Ifges(5,5) * qJDD(4) - pkin(8) * t726 - qJD(4) * t785 - t714 * t839 + t724 * t842 - t784 * t813;
t708 = -mrSges(5,1) * t783 + mrSges(5,3) * t752 + Ifges(5,4) * t798 + Ifges(5,2) * t797 + Ifges(5,6) * qJDD(4) - pkin(4) * t726 + qJD(4) * t786 - t814 * t784 - t895;
t698 = -mrSges(4,1) * t805 + mrSges(4,3) * t795 - pkin(3) * t852 + pkin(7) * t865 + qJDD(1) * t861 + t840 * t701 + t843 * t708 - t838 * t876;
t700 = mrSges(4,2) * t805 - mrSges(4,3) * t794 - pkin(7) * t712 + qJDD(1) * t862 + t843 * t701 - t840 * t708 - t837 * t876;
t705 = qJDD(1) * mrSges(3,2) - t853;
t719 = t846 * t868 + t850;
t847 = -mrSges(2,2) * t819 - mrSges(3,3) * t809 - pkin(1) * t705 - qJ(3) * t707 - t698 * t837 + t838 * t700 + qJ(2) * (t719 + t896) + mrSges(3,2) * t812 + mrSges(2,1) * t818 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t706 = -m(3) * g(3) + t866;
t697 = (t860 + t890) * qJDD(1) - mrSges(2,3) * t818 + mrSges(3,1) * t812 + mrSges(4,1) * t794 - mrSges(4,2) * t795 + pkin(3) * t712 + t849 + pkin(2) * t707 - qJ(2) * t706 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t837 * t862 + t838 * t861 + t887) * t846;
t696 = -mrSges(3,1) * t809 + mrSges(2,3) * t819 - pkin(1) * t706 + pkin(2) * t719 + g(3) * t892 - qJ(3) * t866 - qJDD(1) * t887 - t838 * t698 - t837 * t700 + t846 * t890;
t1 = [-m(1) * g(1) + t867; -m(1) * g(2) + t882; (-m(1) - m(2) - m(3)) * g(3) + t866; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t882 - t841 * t696 + t844 * t697; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t867 + t844 * t696 + t841 * t697; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t847; t847; t705; t719; t849; t895; t736;];
tauJB  = t1;
