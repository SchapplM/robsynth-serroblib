% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 02:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:44:45
% EndTime: 2019-05-06 02:45:00
% DurationCPUTime: 15.57s
% Computational Cost: add. (270102->346), mult. (527649->432), div. (0->0), fcn. (364806->12), ass. (0->141)
t870 = sin(qJ(1));
t875 = cos(qJ(1));
t849 = t870 * g(1) - g(2) * t875;
t840 = qJDD(1) * pkin(1) + t849;
t850 = -g(1) * t875 - g(2) * t870;
t876 = qJD(1) ^ 2;
t842 = -pkin(1) * t876 + t850;
t864 = sin(pkin(11));
t865 = cos(pkin(11));
t819 = t864 * t840 + t865 * t842;
t814 = -pkin(2) * t876 + qJDD(1) * pkin(7) + t819;
t863 = -g(3) + qJDD(2);
t869 = sin(qJ(3));
t874 = cos(qJ(3));
t801 = -t869 * t814 + t874 * t863;
t895 = qJD(1) * qJD(3);
t894 = t874 * t895;
t843 = qJDD(1) * t869 + t894;
t783 = (-t843 + t894) * pkin(8) + (t869 * t874 * t876 + qJDD(3)) * pkin(3) + t801;
t802 = t874 * t814 + t869 * t863;
t844 = qJDD(1) * t874 - t869 * t895;
t897 = qJD(1) * t869;
t848 = qJD(3) * pkin(3) - pkin(8) * t897;
t862 = t874 ^ 2;
t784 = -pkin(3) * t862 * t876 + pkin(8) * t844 - qJD(3) * t848 + t802;
t868 = sin(qJ(4));
t873 = cos(qJ(4));
t769 = t868 * t783 + t873 * t784;
t835 = (t868 * t874 + t869 * t873) * qJD(1);
t803 = -t835 * qJD(4) - t868 * t843 + t844 * t873;
t896 = qJD(1) * t874;
t834 = -t868 * t897 + t873 * t896;
t815 = -mrSges(5,1) * t834 + mrSges(5,2) * t835;
t859 = qJD(3) + qJD(4);
t824 = mrSges(5,1) * t859 - mrSges(5,3) * t835;
t858 = qJDD(3) + qJDD(4);
t816 = -pkin(4) * t834 - pkin(9) * t835;
t857 = t859 ^ 2;
t762 = -pkin(4) * t857 + pkin(9) * t858 + t816 * t834 + t769;
t818 = t865 * t840 - t864 * t842;
t886 = -qJDD(1) * pkin(2) - t818;
t789 = -t844 * pkin(3) + t848 * t897 + (-pkin(8) * t862 - pkin(7)) * t876 + t886;
t804 = qJD(4) * t834 + t843 * t873 + t844 * t868;
t765 = (-t834 * t859 - t804) * pkin(9) + (t835 * t859 - t803) * pkin(4) + t789;
t867 = sin(qJ(5));
t872 = cos(qJ(5));
t752 = -t867 * t762 + t872 * t765;
t821 = -t835 * t867 + t859 * t872;
t777 = qJD(5) * t821 + t804 * t872 + t858 * t867;
t800 = qJDD(5) - t803;
t822 = t835 * t872 + t859 * t867;
t827 = qJD(5) - t834;
t750 = (t821 * t827 - t777) * pkin(10) + (t821 * t822 + t800) * pkin(5) + t752;
t753 = t872 * t762 + t867 * t765;
t776 = -qJD(5) * t822 - t804 * t867 + t858 * t872;
t807 = pkin(5) * t827 - pkin(10) * t822;
t820 = t821 ^ 2;
t751 = -pkin(5) * t820 + pkin(10) * t776 - t807 * t827 + t753;
t866 = sin(qJ(6));
t871 = cos(qJ(6));
t748 = t750 * t871 - t751 * t866;
t790 = t821 * t871 - t822 * t866;
t759 = qJD(6) * t790 + t776 * t866 + t777 * t871;
t791 = t821 * t866 + t822 * t871;
t774 = -mrSges(7,1) * t790 + mrSges(7,2) * t791;
t825 = qJD(6) + t827;
t781 = -mrSges(7,2) * t825 + mrSges(7,3) * t790;
t795 = qJDD(6) + t800;
t743 = m(7) * t748 + mrSges(7,1) * t795 - mrSges(7,3) * t759 - t774 * t791 + t781 * t825;
t749 = t750 * t866 + t751 * t871;
t758 = -qJD(6) * t791 + t776 * t871 - t777 * t866;
t782 = mrSges(7,1) * t825 - mrSges(7,3) * t791;
t744 = m(7) * t749 - mrSges(7,2) * t795 + mrSges(7,3) * t758 + t774 * t790 - t782 * t825;
t735 = t871 * t743 + t866 * t744;
t792 = -mrSges(6,1) * t821 + mrSges(6,2) * t822;
t805 = -mrSges(6,2) * t827 + mrSges(6,3) * t821;
t733 = m(6) * t752 + mrSges(6,1) * t800 - mrSges(6,3) * t777 - t792 * t822 + t805 * t827 + t735;
t806 = mrSges(6,1) * t827 - mrSges(6,3) * t822;
t888 = -t743 * t866 + t871 * t744;
t734 = m(6) * t753 - mrSges(6,2) * t800 + mrSges(6,3) * t776 + t792 * t821 - t806 * t827 + t888;
t889 = -t733 * t867 + t872 * t734;
t726 = m(5) * t769 - mrSges(5,2) * t858 + mrSges(5,3) * t803 + t815 * t834 - t824 * t859 + t889;
t768 = t783 * t873 - t868 * t784;
t823 = -mrSges(5,2) * t859 + mrSges(5,3) * t834;
t761 = -pkin(4) * t858 - pkin(9) * t857 + t835 * t816 - t768;
t754 = -pkin(5) * t776 - pkin(10) * t820 + t807 * t822 + t761;
t885 = m(7) * t754 - t758 * mrSges(7,1) + mrSges(7,2) * t759 - t790 * t781 + t782 * t791;
t880 = -m(6) * t761 + t776 * mrSges(6,1) - mrSges(6,2) * t777 + t821 * t805 - t806 * t822 - t885;
t739 = m(5) * t768 + mrSges(5,1) * t858 - mrSges(5,3) * t804 - t815 * t835 + t823 * t859 + t880;
t716 = t868 * t726 + t873 * t739;
t832 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t869 + Ifges(4,2) * t874) * qJD(1);
t833 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t869 + Ifges(4,4) * t874) * qJD(1);
t770 = Ifges(7,5) * t791 + Ifges(7,6) * t790 + Ifges(7,3) * t825;
t772 = Ifges(7,1) * t791 + Ifges(7,4) * t790 + Ifges(7,5) * t825;
t736 = -mrSges(7,1) * t754 + mrSges(7,3) * t749 + Ifges(7,4) * t759 + Ifges(7,2) * t758 + Ifges(7,6) * t795 - t770 * t791 + t772 * t825;
t771 = Ifges(7,4) * t791 + Ifges(7,2) * t790 + Ifges(7,6) * t825;
t737 = mrSges(7,2) * t754 - mrSges(7,3) * t748 + Ifges(7,1) * t759 + Ifges(7,4) * t758 + Ifges(7,5) * t795 + t770 * t790 - t771 * t825;
t785 = Ifges(6,5) * t822 + Ifges(6,6) * t821 + Ifges(6,3) * t827;
t787 = Ifges(6,1) * t822 + Ifges(6,4) * t821 + Ifges(6,5) * t827;
t718 = -mrSges(6,1) * t761 + mrSges(6,3) * t753 + Ifges(6,4) * t777 + Ifges(6,2) * t776 + Ifges(6,6) * t800 - pkin(5) * t885 + pkin(10) * t888 + t871 * t736 + t866 * t737 - t822 * t785 + t827 * t787;
t786 = Ifges(6,4) * t822 + Ifges(6,2) * t821 + Ifges(6,6) * t827;
t720 = mrSges(6,2) * t761 - mrSges(6,3) * t752 + Ifges(6,1) * t777 + Ifges(6,4) * t776 + Ifges(6,5) * t800 - pkin(10) * t735 - t736 * t866 + t737 * t871 + t785 * t821 - t786 * t827;
t810 = Ifges(5,4) * t835 + Ifges(5,2) * t834 + Ifges(5,6) * t859;
t811 = Ifges(5,1) * t835 + Ifges(5,4) * t834 + Ifges(5,5) * t859;
t881 = -mrSges(5,1) * t768 + mrSges(5,2) * t769 - Ifges(5,5) * t804 - Ifges(5,6) * t803 - Ifges(5,3) * t858 - pkin(4) * t880 - pkin(9) * t889 - t872 * t718 - t867 * t720 - t835 * t810 + t834 * t811;
t899 = mrSges(4,1) * t801 - mrSges(4,2) * t802 + Ifges(4,5) * t843 + Ifges(4,6) * t844 + Ifges(4,3) * qJDD(3) + pkin(3) * t716 + (t832 * t869 - t833 * t874) * qJD(1) - t881;
t841 = (-mrSges(4,1) * t874 + mrSges(4,2) * t869) * qJD(1);
t847 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t896;
t714 = m(4) * t801 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t843 + qJD(3) * t847 - t841 * t897 + t716;
t846 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t897;
t890 = t873 * t726 - t739 * t868;
t715 = m(4) * t802 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t844 - qJD(3) * t846 + t841 * t896 + t890;
t891 = -t714 * t869 + t874 * t715;
t706 = m(3) * t819 - mrSges(3,1) * t876 - qJDD(1) * mrSges(3,2) + t891;
t813 = -t876 * pkin(7) + t886;
t728 = t872 * t733 + t867 * t734;
t883 = m(5) * t789 - t803 * mrSges(5,1) + mrSges(5,2) * t804 - t834 * t823 + t824 * t835 + t728;
t879 = -m(4) * t813 + t844 * mrSges(4,1) - mrSges(4,2) * t843 - t846 * t897 + t847 * t896 - t883;
t722 = m(3) * t818 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t876 + t879;
t702 = t864 * t706 + t865 * t722;
t699 = m(2) * t849 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t876 + t702;
t892 = t865 * t706 - t722 * t864;
t700 = m(2) * t850 - mrSges(2,1) * t876 - qJDD(1) * mrSges(2,2) + t892;
t898 = t875 * t699 + t870 * t700;
t709 = t874 * t714 + t869 * t715;
t707 = m(3) * t863 + t709;
t893 = -t699 * t870 + t875 * t700;
t884 = -mrSges(7,1) * t748 + mrSges(7,2) * t749 - Ifges(7,5) * t759 - Ifges(7,6) * t758 - Ifges(7,3) * t795 - t791 * t771 + t790 * t772;
t809 = Ifges(5,5) * t835 + Ifges(5,6) * t834 + Ifges(5,3) * t859;
t703 = mrSges(5,2) * t789 - mrSges(5,3) * t768 + Ifges(5,1) * t804 + Ifges(5,4) * t803 + Ifges(5,5) * t858 - pkin(9) * t728 - t718 * t867 + t720 * t872 + t809 * t834 - t810 * t859;
t877 = mrSges(6,1) * t752 - mrSges(6,2) * t753 + Ifges(6,5) * t777 + Ifges(6,6) * t776 + Ifges(6,3) * t800 + pkin(5) * t735 + t822 * t786 - t821 * t787 - t884;
t710 = -mrSges(5,1) * t789 + mrSges(5,3) * t769 + Ifges(5,4) * t804 + Ifges(5,2) * t803 + Ifges(5,6) * t858 - pkin(4) * t728 - t835 * t809 + t859 * t811 - t877;
t831 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t869 + Ifges(4,6) * t874) * qJD(1);
t692 = -mrSges(4,1) * t813 + mrSges(4,3) * t802 + Ifges(4,4) * t843 + Ifges(4,2) * t844 + Ifges(4,6) * qJDD(3) - pkin(3) * t883 + pkin(8) * t890 + qJD(3) * t833 + t868 * t703 + t873 * t710 - t831 * t897;
t695 = mrSges(4,2) * t813 - mrSges(4,3) * t801 + Ifges(4,1) * t843 + Ifges(4,4) * t844 + Ifges(4,5) * qJDD(3) - pkin(8) * t716 - qJD(3) * t832 + t703 * t873 - t710 * t868 + t831 * t896;
t882 = mrSges(2,1) * t849 + mrSges(3,1) * t818 - mrSges(2,2) * t850 - mrSges(3,2) * t819 + pkin(1) * t702 + pkin(2) * t879 + pkin(7) * t891 + t874 * t692 + t869 * t695 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t693 = -mrSges(3,1) * t863 + mrSges(3,3) * t819 + t876 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t709 - t899;
t690 = mrSges(3,2) * t863 - mrSges(3,3) * t818 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t876 - pkin(7) * t709 - t692 * t869 + t695 * t874;
t689 = -mrSges(2,2) * g(3) - mrSges(2,3) * t849 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t876 - qJ(2) * t702 + t690 * t865 - t693 * t864;
t688 = mrSges(2,1) * g(3) + mrSges(2,3) * t850 + t876 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t707 + qJ(2) * t892 + t864 * t690 + t865 * t693;
t1 = [-m(1) * g(1) + t893; -m(1) * g(2) + t898; (-m(1) - m(2)) * g(3) + t707; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t898 - t870 * t688 + t875 * t689; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t893 + t875 * t688 + t870 * t689; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t882; t882; t707; t899; -t881; t877; -t884;];
tauJB  = t1;
