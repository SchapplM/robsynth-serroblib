% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 09:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRRRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:51:48
% EndTime: 2019-05-05 09:52:00
% DurationCPUTime: 11.58s
% Computational Cost: add. (194870->318), mult. (375228->395), div. (0->0), fcn. (265818->12), ass. (0->138)
t907 = Ifges(6,1) + Ifges(7,1);
t897 = Ifges(6,4) - Ifges(7,5);
t905 = Ifges(7,4) + Ifges(6,5);
t906 = Ifges(6,2) + Ifges(7,3);
t904 = Ifges(6,6) - Ifges(7,6);
t903 = -Ifges(6,3) - Ifges(7,2);
t865 = sin(qJ(4));
t868 = cos(qJ(4));
t866 = sin(qJ(3));
t890 = qJD(2) * t866;
t843 = qJD(3) * t868 - t865 * t890;
t844 = qJD(3) * t865 + t868 * t890;
t864 = sin(qJ(5));
t899 = cos(qJ(5));
t818 = -t843 * t899 + t844 * t864;
t819 = t864 * t843 + t844 * t899;
t869 = cos(qJ(3));
t889 = qJD(2) * t869;
t856 = qJD(4) - t889;
t855 = qJD(5) + t856;
t902 = t818 * t906 - t819 * t897 - t855 * t904;
t901 = -t897 * t818 + t819 * t907 + t905 * t855;
t860 = sin(pkin(11));
t862 = cos(pkin(11));
t849 = g(1) * t860 - g(2) * t862;
t850 = -g(1) * t862 - g(2) * t860;
t859 = -g(3) + qJDD(1);
t861 = sin(pkin(6));
t863 = cos(pkin(6));
t867 = sin(qJ(2));
t870 = cos(qJ(2));
t807 = -t867 * t850 + (t849 * t863 + t859 * t861) * t870;
t894 = t863 * t867;
t895 = t861 * t867;
t808 = t849 * t894 + t870 * t850 + t859 * t895;
t872 = qJD(2) ^ 2;
t802 = -pkin(2) * t872 + qJDD(2) * pkin(8) + t808;
t825 = -t849 * t861 + t859 * t863;
t792 = -t866 * t802 + t869 * t825;
t846 = (-pkin(3) * t869 - pkin(9) * t866) * qJD(2);
t871 = qJD(3) ^ 2;
t778 = -qJDD(3) * pkin(3) - t871 * pkin(9) + t846 * t890 - t792;
t888 = qJD(2) * qJD(3);
t886 = t869 * t888;
t847 = qJDD(2) * t866 + t886;
t816 = -qJD(4) * t844 + qJDD(3) * t868 - t847 * t865;
t824 = pkin(4) * t856 - pkin(10) * t844;
t839 = t843 ^ 2;
t765 = -t816 * pkin(4) - t839 * pkin(10) + t844 * t824 + t778;
t817 = qJD(4) * t843 + qJDD(3) * t865 + t847 * t868;
t775 = qJD(5) * t819 - t816 * t899 + t817 * t864;
t776 = -t818 * qJD(5) + t864 * t816 + t817 * t899;
t758 = -0.2e1 * qJD(6) * t819 + (t818 * t855 - t776) * qJ(6) + (t819 * t855 + t775) * pkin(5) + t765;
t803 = -mrSges(7,2) * t818 + mrSges(7,3) * t855;
t806 = -mrSges(7,1) * t855 + mrSges(7,2) * t819;
t749 = m(7) * t758 + t775 * mrSges(7,1) - t776 * mrSges(7,3) + t818 * t803 - t819 * t806;
t793 = t869 * t802 + t866 * t825;
t779 = -pkin(3) * t871 + qJDD(3) * pkin(9) + t846 * t889 + t793;
t801 = -qJDD(2) * pkin(2) - t872 * pkin(8) - t807;
t857 = t866 * t888;
t848 = qJDD(2) * t869 - t857;
t782 = (-t847 - t886) * pkin(9) + (-t848 + t857) * pkin(3) + t801;
t763 = -t865 * t779 + t868 * t782;
t840 = qJDD(4) - t848;
t760 = (t843 * t856 - t817) * pkin(10) + (t843 * t844 + t840) * pkin(4) + t763;
t764 = t868 * t779 + t865 * t782;
t762 = -pkin(4) * t839 + pkin(10) * t816 - t824 * t856 + t764;
t756 = t864 * t760 + t899 * t762;
t794 = pkin(5) * t818 - qJ(6) * t819;
t836 = qJDD(5) + t840;
t854 = t855 ^ 2;
t752 = -pkin(5) * t854 + qJ(6) * t836 + 0.2e1 * qJD(6) * t855 - t794 * t818 + t756;
t892 = t904 * t818 - t905 * t819 + t903 * t855;
t734 = -mrSges(6,1) * t765 - mrSges(7,1) * t758 + mrSges(7,2) * t752 + mrSges(6,3) * t756 - pkin(5) * t749 - t775 * t906 + t897 * t776 + t892 * t819 + t904 * t836 + t901 * t855;
t755 = t760 * t899 - t864 * t762;
t753 = -t836 * pkin(5) - t854 * qJ(6) + t819 * t794 + qJDD(6) - t755;
t735 = mrSges(6,2) * t765 + mrSges(7,2) * t753 - mrSges(6,3) * t755 - mrSges(7,3) * t758 - qJ(6) * t749 - t897 * t775 + t776 * t907 + t892 * t818 + t905 * t836 + t902 * t855;
t810 = Ifges(5,5) * t844 + Ifges(5,6) * t843 + Ifges(5,3) * t856;
t812 = Ifges(5,1) * t844 + Ifges(5,4) * t843 + Ifges(5,5) * t856;
t804 = -mrSges(6,2) * t855 - mrSges(6,3) * t818;
t805 = mrSges(6,1) * t855 - mrSges(6,3) * t819;
t877 = m(6) * t765 + t775 * mrSges(6,1) + t776 * mrSges(6,2) + t818 * t804 + t819 * t805 + t749;
t887 = m(7) * t752 + t836 * mrSges(7,3) + t855 * t806;
t795 = mrSges(7,1) * t818 - mrSges(7,3) * t819;
t891 = -mrSges(6,1) * t818 - mrSges(6,2) * t819 - t795;
t898 = -mrSges(6,3) - mrSges(7,2);
t741 = m(6) * t756 - t836 * mrSges(6,2) + t775 * t898 - t855 * t805 + t818 * t891 + t887;
t881 = -m(7) * t753 + t836 * mrSges(7,1) + t855 * t803;
t743 = m(6) * t755 + t836 * mrSges(6,1) + t776 * t898 + t855 * t804 + t819 * t891 + t881;
t883 = t899 * t741 - t743 * t864;
t714 = -mrSges(5,1) * t778 + mrSges(5,3) * t764 + Ifges(5,4) * t817 + Ifges(5,2) * t816 + Ifges(5,6) * t840 - pkin(4) * t877 + pkin(10) * t883 + t734 * t899 + t864 * t735 - t844 * t810 + t856 * t812;
t736 = t864 * t741 + t899 * t743;
t811 = Ifges(5,4) * t844 + Ifges(5,2) * t843 + Ifges(5,6) * t856;
t715 = mrSges(5,2) * t778 - mrSges(5,3) * t763 + Ifges(5,1) * t817 + Ifges(5,4) * t816 + Ifges(5,5) * t840 - pkin(10) * t736 - t864 * t734 + t735 * t899 + t843 * t810 - t856 * t811;
t820 = -mrSges(5,1) * t843 + mrSges(5,2) * t844;
t822 = -mrSges(5,2) * t856 + mrSges(5,3) * t843;
t732 = m(5) * t763 + mrSges(5,1) * t840 - mrSges(5,3) * t817 - t820 * t844 + t822 * t856 + t736;
t823 = mrSges(5,1) * t856 - mrSges(5,3) * t844;
t733 = m(5) * t764 - mrSges(5,2) * t840 + mrSges(5,3) * t816 + t820 * t843 - t823 * t856 + t883;
t730 = -t732 * t865 + t868 * t733;
t744 = -m(5) * t778 + t816 * mrSges(5,1) - t817 * mrSges(5,2) + t843 * t822 - t844 * t823 - t877;
t834 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t866 + Ifges(4,2) * t869) * qJD(2);
t835 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t866 + Ifges(4,4) * t869) * qJD(2);
t900 = mrSges(4,1) * t792 - mrSges(4,2) * t793 + Ifges(4,5) * t847 + Ifges(4,6) * t848 + Ifges(4,3) * qJDD(3) + pkin(3) * t744 + pkin(9) * t730 + t868 * t714 + t865 * t715 + (t834 * t866 - t835 * t869) * qJD(2);
t729 = t732 * t868 + t733 * t865;
t851 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t890;
t852 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t889;
t876 = -m(4) * t801 + t848 * mrSges(4,1) - mrSges(4,2) * t847 - t851 * t890 + t852 * t889 - t729;
t725 = m(3) * t807 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t872 + t876;
t896 = t725 * t870;
t845 = (-mrSges(4,1) * t869 + mrSges(4,2) * t866) * qJD(2);
t728 = m(4) * t793 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t848 - qJD(3) * t851 + t845 * t889 + t730;
t738 = m(4) * t792 + qJDD(3) * mrSges(4,1) - t847 * mrSges(4,3) + qJD(3) * t852 - t845 * t890 + t744;
t884 = t869 * t728 - t738 * t866;
t719 = m(3) * t808 - mrSges(3,1) * t872 - qJDD(2) * mrSges(3,2) + t884;
t722 = t866 * t728 + t869 * t738;
t721 = m(3) * t825 + t722;
t708 = t719 * t894 - t721 * t861 + t863 * t896;
t706 = m(2) * t849 + t708;
t712 = t870 * t719 - t725 * t867;
t711 = m(2) * t850 + t712;
t893 = t862 * t706 + t860 * t711;
t707 = t719 * t895 + t863 * t721 + t861 * t896;
t885 = -t706 * t860 + t862 * t711;
t882 = m(2) * t859 + t707;
t833 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t866 + Ifges(4,6) * t869) * qJD(2);
t704 = mrSges(4,2) * t801 - mrSges(4,3) * t792 + Ifges(4,1) * t847 + Ifges(4,4) * t848 + Ifges(4,5) * qJDD(3) - pkin(9) * t729 - qJD(3) * t834 - t714 * t865 + t715 * t868 + t833 * t889;
t748 = t776 * mrSges(7,2) + t819 * t795 - t881;
t875 = -mrSges(6,1) * t755 + mrSges(7,1) * t753 + mrSges(6,2) * t756 - mrSges(7,3) * t752 + pkin(5) * t748 - qJ(6) * t887 + t903 * t836 + t902 * t819 + (qJ(6) * t795 - t901) * t818 - t905 * t776 + (mrSges(7,2) * qJ(6) + t904) * t775;
t873 = mrSges(5,1) * t763 - mrSges(5,2) * t764 + Ifges(5,5) * t817 + Ifges(5,6) * t816 + Ifges(5,3) * t840 + pkin(4) * t736 + t844 * t811 - t843 * t812 - t875;
t713 = -mrSges(4,1) * t801 + mrSges(4,3) * t793 + Ifges(4,4) * t847 + Ifges(4,2) * t848 + Ifges(4,6) * qJDD(3) - pkin(3) * t729 + qJD(3) * t835 - t833 * t890 - t873;
t702 = mrSges(3,2) * t825 - mrSges(3,3) * t807 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t872 - pkin(8) * t722 + t704 * t869 - t713 * t866;
t703 = -mrSges(3,1) * t825 + mrSges(3,3) * t808 + t872 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t722 - t900;
t878 = pkin(7) * t712 + t702 * t867 + t703 * t870;
t701 = mrSges(3,1) * t807 - mrSges(3,2) * t808 + Ifges(3,3) * qJDD(2) + pkin(2) * t876 + pkin(8) * t884 + t866 * t704 + t869 * t713;
t700 = mrSges(2,2) * t859 - mrSges(2,3) * t849 + t870 * t702 - t867 * t703 + (-t707 * t861 - t708 * t863) * pkin(7);
t699 = -mrSges(2,1) * t859 + mrSges(2,3) * t850 - pkin(1) * t707 - t861 * t701 + t863 * t878;
t1 = [-m(1) * g(1) + t885; -m(1) * g(2) + t893; -m(1) * g(3) + t882; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t893 - t860 * t699 + t862 * t700; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t885 + t862 * t699 + t860 * t700; -mrSges(1,1) * g(2) + mrSges(2,1) * t849 + mrSges(1,2) * g(1) - mrSges(2,2) * t850 + pkin(1) * t708 + t863 * t701 + t861 * t878; t882; t701; t900; t873; -t875; t748;];
tauJB  = t1;
