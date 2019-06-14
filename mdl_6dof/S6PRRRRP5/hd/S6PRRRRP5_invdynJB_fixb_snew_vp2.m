% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRRP5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 10:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:02:18
% EndTime: 2019-05-05 10:02:42
% DurationCPUTime: 23.50s
% Computational Cost: add. (408903->331), mult. (839332->425), div. (0->0), fcn. (666136->14), ass. (0->151)
t940 = Ifges(6,1) + Ifges(7,1);
t934 = Ifges(6,4) + Ifges(7,4);
t933 = Ifges(6,5) + Ifges(7,5);
t939 = Ifges(6,2) + Ifges(7,2);
t932 = Ifges(6,6) + Ifges(7,6);
t938 = Ifges(6,3) + Ifges(7,3);
t888 = sin(pkin(7));
t895 = sin(qJ(3));
t899 = cos(qJ(3));
t917 = qJD(2) * qJD(3);
t872 = (-qJDD(2) * t899 + t895 * t917) * t888;
t887 = sin(pkin(12));
t890 = cos(pkin(12));
t878 = g(1) * t887 - g(2) * t890;
t879 = -g(1) * t890 - g(2) * t887;
t886 = -g(3) + qJDD(1);
t896 = sin(qJ(2));
t892 = cos(pkin(6));
t900 = cos(qJ(2));
t924 = t892 * t900;
t889 = sin(pkin(6));
t927 = t889 * t900;
t849 = t878 * t924 - t879 * t896 + t886 * t927;
t901 = qJD(2) ^ 2;
t936 = pkin(9) * t888;
t845 = qJDD(2) * pkin(2) + t901 * t936 + t849;
t925 = t892 * t896;
t928 = t889 * t896;
t850 = t878 * t925 + t900 * t879 + t886 * t928;
t846 = -pkin(2) * t901 + qJDD(2) * t936 + t850;
t865 = -t878 * t889 + t886 * t892;
t891 = cos(pkin(7));
t807 = -t895 * t846 + (t845 * t891 + t865 * t888) * t899;
t884 = qJD(2) * t891 + qJD(3);
t894 = sin(qJ(4));
t898 = cos(qJ(4));
t918 = qJD(2) * t888;
t914 = t895 * t918;
t863 = t884 * t898 - t894 * t914;
t871 = (qJDD(2) * t895 + t899 * t917) * t888;
t883 = qJDD(2) * t891 + qJDD(3);
t841 = qJD(4) * t863 + t871 * t898 + t883 * t894;
t864 = t884 * t894 + t898 * t914;
t913 = t899 * t918;
t877 = qJD(4) - t913;
t893 = sin(qJ(5));
t897 = cos(qJ(5));
t852 = -t864 * t893 + t877 * t897;
t866 = qJDD(4) + t872;
t813 = qJD(5) * t852 + t841 * t897 + t866 * t893;
t853 = t864 * t897 + t877 * t893;
t824 = -mrSges(7,1) * t852 + mrSges(7,2) * t853;
t926 = t891 * t895;
t929 = t888 * t895;
t808 = t845 * t926 + t899 * t846 + t865 * t929;
t870 = (-pkin(3) * t899 - pkin(10) * t895) * t918;
t882 = t884 ^ 2;
t804 = -pkin(3) * t882 + pkin(10) * t883 + t870 * t913 + t808;
t861 = t891 * t865;
t806 = t872 * pkin(3) - t871 * pkin(10) + t861 + (-t845 + (pkin(3) * t895 - pkin(10) * t899) * t884 * qJD(2)) * t888;
t800 = t898 * t804 + t894 * t806;
t848 = -pkin(4) * t863 - pkin(11) * t864;
t876 = t877 ^ 2;
t795 = -pkin(4) * t876 + pkin(11) * t866 + t848 * t863 + t800;
t803 = -t883 * pkin(3) - t882 * pkin(10) + t870 * t914 - t807;
t840 = -qJD(4) * t864 - t871 * t894 + t883 * t898;
t798 = (-t863 * t877 - t841) * pkin(11) + (t864 * t877 - t840) * pkin(4) + t803;
t790 = -t893 * t795 + t897 * t798;
t838 = qJDD(5) - t840;
t862 = qJD(5) - t863;
t787 = -0.2e1 * qJD(6) * t853 + (t852 * t862 - t813) * qJ(6) + (t852 * t853 + t838) * pkin(5) + t790;
t828 = -mrSges(7,2) * t862 + mrSges(7,3) * t852;
t916 = m(7) * t787 + t838 * mrSges(7,1) + t862 * t828;
t784 = -t813 * mrSges(7,3) - t853 * t824 + t916;
t791 = t897 * t795 + t893 * t798;
t812 = -qJD(5) * t853 - t841 * t893 + t866 * t897;
t830 = pkin(5) * t862 - qJ(6) * t853;
t851 = t852 ^ 2;
t789 = -pkin(5) * t851 + qJ(6) * t812 + 0.2e1 * qJD(6) * t852 - t830 * t862 + t791;
t920 = t934 * t852 + t940 * t853 + t933 * t862;
t921 = -t939 * t852 - t934 * t853 - t932 * t862;
t937 = mrSges(6,1) * t790 + mrSges(7,1) * t787 - mrSges(6,2) * t791 - mrSges(7,2) * t789 + pkin(5) * t784 + t812 * t932 + t813 * t933 + t938 * t838 - t852 * t920 - t853 * t921;
t935 = -mrSges(6,2) - mrSges(7,2);
t868 = -mrSges(4,2) * t884 + mrSges(4,3) * t913;
t869 = (-mrSges(4,1) * t899 + mrSges(4,2) * t895) * t918;
t825 = -mrSges(6,1) * t852 + mrSges(6,2) * t853;
t829 = -mrSges(6,2) * t862 + mrSges(6,3) * t852;
t778 = m(6) * t790 + t838 * mrSges(6,1) + t862 * t829 + (-t824 - t825) * t853 + (-mrSges(6,3) - mrSges(7,3)) * t813 + t916;
t915 = m(7) * t789 + t812 * mrSges(7,3) + t852 * t824;
t831 = mrSges(7,1) * t862 - mrSges(7,3) * t853;
t919 = -mrSges(6,1) * t862 + mrSges(6,3) * t853 - t831;
t780 = m(6) * t791 + t812 * mrSges(6,3) + t852 * t825 + t838 * t935 + t862 * t919 + t915;
t776 = t778 * t897 + t780 * t893;
t854 = -mrSges(5,2) * t877 + mrSges(5,3) * t863;
t855 = mrSges(5,1) * t877 - mrSges(5,3) * t864;
t903 = -m(5) * t803 + t840 * mrSges(5,1) - mrSges(5,2) * t841 + t863 * t854 - t855 * t864 - t776;
t771 = m(4) * t807 + mrSges(4,1) * t883 - mrSges(4,3) * t871 + t868 * t884 - t869 * t914 + t903;
t930 = t771 * t899;
t867 = mrSges(4,1) * t884 - mrSges(4,3) * t914;
t777 = -t778 * t893 + t897 * t780;
t847 = -mrSges(5,1) * t863 + mrSges(5,2) * t864;
t774 = m(5) * t800 - mrSges(5,2) * t866 + mrSges(5,3) * t840 + t847 * t863 - t855 * t877 + t777;
t799 = -t894 * t804 + t806 * t898;
t794 = -pkin(4) * t866 - pkin(11) * t876 + t864 * t848 - t799;
t792 = -pkin(5) * t812 - qJ(6) * t851 + t830 * t853 + qJDD(6) + t794;
t909 = -m(7) * t792 + t812 * mrSges(7,1) + t852 * t828;
t783 = -m(6) * t794 + t812 * mrSges(6,1) + t813 * t935 + t852 * t829 + t853 * t919 + t909;
t782 = m(5) * t799 + t866 * mrSges(5,1) - t841 * mrSges(5,3) - t864 * t847 + t877 * t854 + t783;
t911 = t898 * t774 - t782 * t894;
t764 = m(4) * t808 - mrSges(4,2) * t883 - mrSges(4,3) * t872 - t867 * t884 + t869 * t913 + t911;
t767 = t894 * t774 + t898 * t782;
t826 = -t888 * t845 + t861;
t766 = m(4) * t826 + t872 * mrSges(4,1) + t871 * mrSges(4,2) + (t867 * t895 - t868 * t899) * t918 + t767;
t753 = t764 * t926 - t766 * t888 + t891 * t930;
t749 = m(3) * t849 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t901 + t753;
t752 = t764 * t929 + t891 * t766 + t888 * t930;
t751 = m(3) * t865 + t752;
t759 = t899 * t764 - t771 * t895;
t758 = m(3) * t850 - mrSges(3,1) * t901 - qJDD(2) * mrSges(3,2) + t759;
t739 = t749 * t924 - t751 * t889 + t758 * t925;
t737 = m(2) * t878 + t739;
t745 = -t749 * t896 + t900 * t758;
t744 = m(2) * t879 + t745;
t923 = t890 * t737 + t887 * t744;
t922 = -t932 * t852 - t933 * t853 - t938 * t862;
t738 = t749 * t927 + t892 * t751 + t758 * t928;
t912 = -t737 * t887 + t890 * t744;
t910 = m(2) * t886 + t738;
t786 = t813 * mrSges(7,2) + t853 * t831 - t909;
t768 = -mrSges(6,1) * t794 + mrSges(6,3) * t791 - mrSges(7,1) * t792 + mrSges(7,3) * t789 - pkin(5) * t786 + qJ(6) * t915 + (-qJ(6) * t831 + t920) * t862 + t922 * t853 + (-mrSges(7,2) * qJ(6) + t932) * t838 + t934 * t813 + t939 * t812;
t775 = mrSges(6,2) * t794 + mrSges(7,2) * t792 - mrSges(6,3) * t790 - mrSges(7,3) * t787 - qJ(6) * t784 + t934 * t812 + t940 * t813 + t933 * t838 - t922 * t852 + t921 * t862;
t834 = Ifges(5,5) * t864 + Ifges(5,6) * t863 + Ifges(5,3) * t877;
t835 = Ifges(5,4) * t864 + Ifges(5,2) * t863 + Ifges(5,6) * t877;
t754 = mrSges(5,2) * t803 - mrSges(5,3) * t799 + Ifges(5,1) * t841 + Ifges(5,4) * t840 + Ifges(5,5) * t866 - pkin(11) * t776 - t768 * t893 + t775 * t897 + t834 * t863 - t835 * t877;
t836 = Ifges(5,1) * t864 + Ifges(5,4) * t863 + Ifges(5,5) * t877;
t760 = -mrSges(5,1) * t803 + mrSges(5,3) * t800 + Ifges(5,4) * t841 + Ifges(5,2) * t840 + Ifges(5,6) * t866 - pkin(4) * t776 - t864 * t834 + t877 * t836 - t937;
t858 = Ifges(4,6) * t884 + (Ifges(4,4) * t895 + Ifges(4,2) * t899) * t918;
t859 = Ifges(4,5) * t884 + (Ifges(4,1) * t895 + Ifges(4,4) * t899) * t918;
t740 = Ifges(4,5) * t871 - Ifges(4,6) * t872 + Ifges(4,3) * t883 + mrSges(4,1) * t807 - mrSges(4,2) * t808 + t894 * t754 + t898 * t760 + pkin(3) * t903 + pkin(10) * t911 + (t858 * t895 - t859 * t899) * t918;
t857 = Ifges(4,3) * t884 + (Ifges(4,5) * t895 + Ifges(4,6) * t899) * t918;
t741 = mrSges(4,2) * t826 - mrSges(4,3) * t807 + Ifges(4,1) * t871 - Ifges(4,4) * t872 + Ifges(4,5) * t883 - pkin(10) * t767 + t754 * t898 - t760 * t894 + t857 * t913 - t858 * t884;
t902 = mrSges(5,1) * t799 - mrSges(5,2) * t800 + Ifges(5,5) * t841 + Ifges(5,6) * t840 + Ifges(5,3) * t866 + pkin(4) * t783 + pkin(11) * t777 + t897 * t768 + t893 * t775 + t864 * t835 - t863 * t836;
t746 = -mrSges(4,1) * t826 + mrSges(4,3) * t808 + Ifges(4,4) * t871 - Ifges(4,2) * t872 + Ifges(4,6) * t883 - pkin(3) * t767 - t857 * t914 + t884 * t859 - t902;
t905 = pkin(9) * t759 + t741 * t895 + t746 * t899;
t734 = -mrSges(3,1) * t865 + mrSges(3,3) * t850 + t901 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t752 - t888 * t740 + t891 * t905;
t735 = mrSges(3,2) * t865 - mrSges(3,3) * t849 + Ifges(3,5) * qJDD(2) - t901 * Ifges(3,6) + t899 * t741 - t895 * t746 + (-t752 * t888 - t753 * t891) * pkin(9);
t906 = pkin(8) * t745 + t734 * t900 + t735 * t896;
t733 = mrSges(3,1) * t849 - mrSges(3,2) * t850 + Ifges(3,3) * qJDD(2) + pkin(2) * t753 + t891 * t740 + t888 * t905;
t732 = mrSges(2,2) * t886 - mrSges(2,3) * t878 - t896 * t734 + t900 * t735 + (-t738 * t889 - t739 * t892) * pkin(8);
t731 = -mrSges(2,1) * t886 + mrSges(2,3) * t879 - pkin(1) * t738 - t889 * t733 + t892 * t906;
t1 = [-m(1) * g(1) + t912; -m(1) * g(2) + t923; -m(1) * g(3) + t910; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t923 - t887 * t731 + t890 * t732; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t912 + t890 * t731 + t887 * t732; -mrSges(1,1) * g(2) + mrSges(2,1) * t878 + mrSges(1,2) * g(1) - mrSges(2,2) * t879 + pkin(1) * t739 + t892 * t733 + t889 * t906; t910; t733; t740; t902; t937; t786;];
tauJB  = t1;
