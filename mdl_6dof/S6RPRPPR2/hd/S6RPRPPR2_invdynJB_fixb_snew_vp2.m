% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-05-05 16:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:33:10
% EndTime: 2019-05-05 16:33:18
% DurationCPUTime: 6.80s
% Computational Cost: add. (76911->326), mult. (167411->391), div. (0->0), fcn. (104739->10), ass. (0->134)
t907 = -2 * qJD(4);
t906 = Ifges(5,1) + Ifges(6,2);
t905 = -Ifges(6,1) - Ifges(5,3);
t901 = Ifges(5,4) + Ifges(6,6);
t900 = Ifges(5,5) - Ifges(6,4);
t904 = -Ifges(5,2) - Ifges(6,3);
t899 = Ifges(5,6) - Ifges(6,5);
t863 = sin(qJ(1));
t866 = cos(qJ(1));
t846 = t863 * g(1) - t866 * g(2);
t837 = qJDD(1) * pkin(1) + t846;
t847 = -t866 * g(1) - t863 * g(2);
t868 = qJD(1) ^ 2;
t839 = -t868 * pkin(1) + t847;
t859 = sin(pkin(9));
t860 = cos(pkin(9));
t808 = t859 * t837 + t860 * t839;
t793 = -t868 * pkin(2) + qJDD(1) * pkin(7) + t808;
t857 = -g(3) + qJDD(2);
t862 = sin(qJ(3));
t865 = cos(qJ(3));
t781 = -t862 * t793 + t865 * t857;
t887 = qJD(1) * qJD(3);
t885 = t865 * t887;
t840 = t862 * qJDD(1) + t885;
t766 = (-t840 + t885) * qJ(4) + (t862 * t865 * t868 + qJDD(3)) * pkin(3) + t781;
t782 = t865 * t793 + t862 * t857;
t841 = t865 * qJDD(1) - t862 * t887;
t891 = qJD(1) * t862;
t843 = qJD(3) * pkin(3) - qJ(4) * t891;
t856 = t865 ^ 2;
t767 = -t856 * t868 * pkin(3) + t841 * qJ(4) - qJD(3) * t843 + t782;
t858 = sin(pkin(10));
t898 = cos(pkin(10));
t826 = (t858 * t865 + t898 * t862) * qJD(1);
t760 = t898 * t766 - t858 * t767 + t826 * t907;
t890 = qJD(1) * t865;
t825 = t858 * t891 - t898 * t890;
t798 = t825 * mrSges(5,1) + t826 * mrSges(5,2);
t810 = t898 * t840 + t858 * t841;
t814 = -qJD(3) * mrSges(5,2) - t825 * mrSges(5,3);
t816 = t825 * mrSges(6,1) - qJD(3) * mrSges(6,3);
t797 = t825 * pkin(4) - t826 * qJ(5);
t867 = qJD(3) ^ 2;
t757 = -qJDD(3) * pkin(4) - t867 * qJ(5) + t826 * t797 + qJDD(5) - t760;
t889 = qJD(3) * t825;
t752 = (t825 * t826 - qJDD(3)) * pkin(8) + (t810 + t889) * pkin(5) + t757;
t809 = t858 * t840 - t898 * t841;
t818 = t826 * pkin(5) - qJD(3) * pkin(8);
t824 = t825 ^ 2;
t807 = t860 * t837 - t859 * t839;
t878 = -qJDD(1) * pkin(2) - t807;
t768 = -t841 * pkin(3) + qJDD(4) + t843 * t891 + (-qJ(4) * t856 - pkin(7)) * t868 + t878;
t902 = -2 * qJD(5);
t871 = (-t810 + t889) * qJ(5) + t768 + (pkin(4) * qJD(3) + t902) * t826;
t755 = -t824 * pkin(5) - t826 * t818 + (pkin(4) + pkin(8)) * t809 + t871;
t861 = sin(qJ(6));
t864 = cos(qJ(6));
t750 = t864 * t752 - t861 * t755;
t811 = -t861 * qJD(3) + t864 * t825;
t777 = t811 * qJD(6) + t864 * qJDD(3) + t861 * t809;
t812 = t864 * qJD(3) + t861 * t825;
t778 = -t811 * mrSges(7,1) + t812 * mrSges(7,2);
t823 = qJD(6) + t826;
t783 = -t823 * mrSges(7,2) + t811 * mrSges(7,3);
t806 = qJDD(6) + t810;
t747 = m(7) * t750 + t806 * mrSges(7,1) - t777 * mrSges(7,3) - t812 * t778 + t823 * t783;
t751 = t861 * t752 + t864 * t755;
t776 = -t812 * qJD(6) - t861 * qJDD(3) + t864 * t809;
t784 = t823 * mrSges(7,1) - t812 * mrSges(7,3);
t748 = m(7) * t751 - t806 * mrSges(7,2) + t776 * mrSges(7,3) + t811 * t778 - t823 * t784;
t738 = t864 * t747 + t861 * t748;
t799 = -t825 * mrSges(6,2) - t826 * mrSges(6,3);
t875 = -m(6) * t757 - t810 * mrSges(6,1) - t826 * t799 - t738;
t733 = m(5) * t760 - t810 * mrSges(5,3) - t826 * t798 + (mrSges(5,1) - mrSges(6,2)) * qJDD(3) + (t814 - t816) * qJD(3) + t875;
t821 = t825 * t907;
t895 = t858 * t766 + t898 * t767;
t761 = t821 + t895;
t815 = qJD(3) * mrSges(5,1) - t826 * mrSges(5,3);
t877 = t867 * pkin(4) - qJDD(3) * qJ(5) - t895;
t756 = qJD(3) * t902 + ((2 * qJD(4)) + t797) * t825 + t877;
t817 = t826 * mrSges(6,1) + qJD(3) * mrSges(6,2);
t754 = -t809 * pkin(5) - t824 * pkin(8) - t825 * t797 + t821 + ((2 * qJD(5)) + t818) * qJD(3) - t877;
t876 = -m(7) * t754 + t776 * mrSges(7,1) - t777 * mrSges(7,2) + t811 * t783 - t812 * t784;
t873 = -m(6) * t756 + qJDD(3) * mrSges(6,3) + qJD(3) * t817 - t876;
t743 = m(5) * t761 - qJDD(3) * mrSges(5,2) - qJD(3) * t815 + (-t798 - t799) * t825 + (-mrSges(5,3) - mrSges(6,1)) * t809 + t873;
t729 = t898 * t733 + t858 * t743;
t736 = qJDD(3) * mrSges(6,2) + qJD(3) * t816 - t875;
t769 = Ifges(7,5) * t812 + Ifges(7,6) * t811 + Ifges(7,3) * t823;
t771 = Ifges(7,1) * t812 + Ifges(7,4) * t811 + Ifges(7,5) * t823;
t739 = -mrSges(7,1) * t754 + mrSges(7,3) * t751 + Ifges(7,4) * t777 + Ifges(7,2) * t776 + Ifges(7,6) * t806 - t812 * t769 + t823 * t771;
t770 = Ifges(7,4) * t812 + Ifges(7,2) * t811 + Ifges(7,6) * t823;
t740 = mrSges(7,2) * t754 - mrSges(7,3) * t750 + Ifges(7,1) * t777 + Ifges(7,4) * t776 + Ifges(7,5) * t806 + t811 * t769 - t823 * t770;
t831 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t862 + Ifges(4,2) * t865) * qJD(1);
t832 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t862 + Ifges(4,4) * t865) * qJD(1);
t892 = qJD(3) * t900 - t825 * t901 + t826 * t906;
t893 = qJD(3) * t899 + t825 * t904 + t826 * t901;
t903 = (t862 * t831 - t865 * t832) * qJD(1) + (Ifges(4,3) - t905) * qJDD(3) - t899 * t809 + t900 * t810 + t892 * t825 + t893 * t826 + mrSges(4,1) * t781 + mrSges(5,1) * t760 - mrSges(4,2) * t782 - mrSges(5,2) * t761 + mrSges(6,2) * t757 - mrSges(6,3) * t756 + Ifges(4,5) * t840 + Ifges(4,6) * t841 + pkin(3) * t729 - pkin(4) * t736 - pkin(8) * t738 + qJ(5) * (-t809 * mrSges(6,1) - t825 * t799 + t873) - t861 * t739 + t864 * t740;
t838 = (-mrSges(4,1) * t865 + mrSges(4,2) * t862) * qJD(1);
t845 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t890;
t727 = m(4) * t781 + qJDD(3) * mrSges(4,1) - t840 * mrSges(4,3) + qJD(3) * t845 - t838 * t891 + t729;
t844 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t891;
t881 = -t858 * t733 + t898 * t743;
t728 = m(4) * t782 - qJDD(3) * mrSges(4,2) + t841 * mrSges(4,3) - qJD(3) * t844 + t838 * t890 + t881;
t882 = -t862 * t727 + t865 * t728;
t719 = m(3) * t808 - t868 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t882;
t759 = t809 * pkin(4) + t871;
t896 = -t861 * t747 + t864 * t748;
t737 = m(6) * t759 - t809 * mrSges(6,2) - t810 * mrSges(6,3) - t825 * t816 - t826 * t817 + t896;
t735 = m(5) * t768 + t809 * mrSges(5,1) + t810 * mrSges(5,2) + t825 * t814 + t826 * t815 + t737;
t792 = -t868 * pkin(7) + t878;
t870 = -m(4) * t792 + t841 * mrSges(4,1) - t840 * mrSges(4,2) - t844 * t891 + t845 * t890 - t735;
t731 = m(3) * t807 + qJDD(1) * mrSges(3,1) - t868 * mrSges(3,2) + t870;
t715 = t859 * t719 + t860 * t731;
t712 = m(2) * t846 + qJDD(1) * mrSges(2,1) - t868 * mrSges(2,2) + t715;
t883 = t860 * t719 - t859 * t731;
t713 = m(2) * t847 - t868 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t883;
t897 = t866 * t712 + t863 * t713;
t722 = t865 * t727 + t862 * t728;
t894 = qJD(3) * t905 + t825 * t899 - t826 * t900;
t720 = m(3) * t857 + t722;
t884 = -t863 * t712 + t866 * t713;
t874 = mrSges(7,1) * t750 - mrSges(7,2) * t751 + Ifges(7,5) * t777 + Ifges(7,6) * t776 + Ifges(7,3) * t806 + t812 * t770 - t811 * t771;
t716 = -mrSges(5,1) * t768 - mrSges(6,1) * t756 + mrSges(6,2) * t759 + mrSges(5,3) * t761 - pkin(4) * t737 - pkin(5) * t876 - pkin(8) * t896 + t892 * qJD(3) + t899 * qJDD(3) - t864 * t739 - t861 * t740 + t809 * t904 + t901 * t810 + t894 * t826;
t723 = mrSges(6,1) * t757 + mrSges(5,2) * t768 - mrSges(5,3) * t760 - mrSges(6,3) * t759 + pkin(5) * t738 - qJ(5) * t737 - t893 * qJD(3) + t900 * qJDD(3) - t901 * t809 + t810 * t906 + t894 * t825 + t874;
t830 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t862 + Ifges(4,6) * t865) * qJD(1);
t705 = -mrSges(4,1) * t792 + mrSges(4,3) * t782 + Ifges(4,4) * t840 + Ifges(4,2) * t841 + Ifges(4,6) * qJDD(3) - pkin(3) * t735 + qJ(4) * t881 + qJD(3) * t832 + t898 * t716 + t858 * t723 - t830 * t891;
t708 = mrSges(4,2) * t792 - mrSges(4,3) * t781 + Ifges(4,1) * t840 + Ifges(4,4) * t841 + Ifges(4,5) * qJDD(3) - qJ(4) * t729 - qJD(3) * t831 - t858 * t716 + t898 * t723 + t830 * t890;
t872 = mrSges(2,1) * t846 + mrSges(3,1) * t807 - mrSges(2,2) * t847 - mrSges(3,2) * t808 + pkin(1) * t715 + pkin(2) * t870 + pkin(7) * t882 + t865 * t705 + t862 * t708 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t706 = -mrSges(3,1) * t857 + mrSges(3,3) * t808 + t868 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t722 - t903;
t703 = mrSges(3,2) * t857 - mrSges(3,3) * t807 + Ifges(3,5) * qJDD(1) - t868 * Ifges(3,6) - pkin(7) * t722 - t862 * t705 + t865 * t708;
t702 = -mrSges(2,2) * g(3) - mrSges(2,3) * t846 + Ifges(2,5) * qJDD(1) - t868 * Ifges(2,6) - qJ(2) * t715 + t860 * t703 - t859 * t706;
t701 = mrSges(2,1) * g(3) + mrSges(2,3) * t847 + t868 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t720 + qJ(2) * t883 + t859 * t703 + t860 * t706;
t1 = [-m(1) * g(1) + t884; -m(1) * g(2) + t897; (-m(1) - m(2)) * g(3) + t720; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t897 - t863 * t701 + t866 * t702; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t884 + t866 * t701 + t863 * t702; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t872; t872; t720; t903; t735; t736; t874;];
tauJB  = t1;
