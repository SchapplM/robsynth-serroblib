% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 18:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:14:14
% EndTime: 2019-05-05 18:14:29
% DurationCPUTime: 15.34s
% Computational Cost: add. (241465->344), mult. (532125->434), div. (0->0), fcn. (371278->12), ass. (0->138)
t855 = sin(qJ(1));
t859 = cos(qJ(1));
t834 = t855 * g(1) - g(2) * t859;
t825 = qJDD(1) * pkin(1) + t834;
t835 = -g(1) * t859 - g(2) * t855;
t860 = qJD(1) ^ 2;
t827 = -pkin(1) * t860 + t835;
t849 = sin(pkin(10));
t851 = cos(pkin(10));
t805 = t849 * t825 + t851 * t827;
t798 = -pkin(2) * t860 + qJDD(1) * pkin(7) + t805;
t847 = -g(3) + qJDD(2);
t854 = sin(qJ(3));
t858 = cos(qJ(3));
t786 = -t798 * t854 + t858 * t847;
t877 = qJD(1) * qJD(3);
t876 = t858 * t877;
t828 = qJDD(1) * t854 + t876;
t781 = (-t828 + t876) * qJ(4) + (t854 * t858 * t860 + qJDD(3)) * pkin(3) + t786;
t787 = t858 * t798 + t854 * t847;
t829 = qJDD(1) * t858 - t854 * t877;
t879 = qJD(1) * t854;
t831 = qJD(3) * pkin(3) - qJ(4) * t879;
t846 = t858 ^ 2;
t782 = -pkin(3) * t846 * t860 + qJ(4) * t829 - qJD(3) * t831 + t787;
t848 = sin(pkin(11));
t850 = cos(pkin(11));
t815 = (t848 * t858 + t850 * t854) * qJD(1);
t749 = -0.2e1 * qJD(4) * t815 + t850 * t781 - t782 * t848;
t807 = t828 * t850 + t829 * t848;
t814 = (-t848 * t854 + t850 * t858) * qJD(1);
t746 = (qJD(3) * t814 - t807) * pkin(8) + (t814 * t815 + qJDD(3)) * pkin(4) + t749;
t750 = 0.2e1 * qJD(4) * t814 + t848 * t781 + t850 * t782;
t806 = -t828 * t848 + t829 * t850;
t810 = qJD(3) * pkin(4) - pkin(8) * t815;
t813 = t814 ^ 2;
t748 = -pkin(4) * t813 + pkin(8) * t806 - qJD(3) * t810 + t750;
t853 = sin(qJ(5));
t857 = cos(qJ(5));
t743 = t853 * t746 + t857 * t748;
t796 = t814 * t853 + t815 * t857;
t766 = -qJD(5) * t796 + t806 * t857 - t807 * t853;
t795 = t814 * t857 - t815 * t853;
t779 = -mrSges(6,1) * t795 + mrSges(6,2) * t796;
t843 = qJD(3) + qJD(5);
t789 = mrSges(6,1) * t843 - mrSges(6,3) * t796;
t842 = qJDD(3) + qJDD(5);
t780 = -pkin(5) * t795 - pkin(9) * t796;
t841 = t843 ^ 2;
t740 = -pkin(5) * t841 + pkin(9) * t842 + t780 * t795 + t743;
t804 = t825 * t851 - t849 * t827;
t868 = -qJDD(1) * pkin(2) - t804;
t783 = -pkin(3) * t829 + qJDD(4) + t831 * t879 + (-qJ(4) * t846 - pkin(7)) * t860 + t868;
t755 = -pkin(4) * t806 - pkin(8) * t813 + t815 * t810 + t783;
t767 = qJD(5) * t795 + t806 * t853 + t807 * t857;
t744 = (-t795 * t843 - t767) * pkin(9) + (t796 * t843 - t766) * pkin(5) + t755;
t852 = sin(qJ(6));
t856 = cos(qJ(6));
t737 = -t740 * t852 + t744 * t856;
t784 = -t796 * t852 + t843 * t856;
t753 = qJD(6) * t784 + t767 * t856 + t842 * t852;
t765 = qJDD(6) - t766;
t785 = t796 * t856 + t843 * t852;
t768 = -mrSges(7,1) * t784 + mrSges(7,2) * t785;
t791 = qJD(6) - t795;
t769 = -mrSges(7,2) * t791 + mrSges(7,3) * t784;
t733 = m(7) * t737 + mrSges(7,1) * t765 - mrSges(7,3) * t753 - t768 * t785 + t769 * t791;
t738 = t740 * t856 + t744 * t852;
t752 = -qJD(6) * t785 - t767 * t852 + t842 * t856;
t770 = mrSges(7,1) * t791 - mrSges(7,3) * t785;
t734 = m(7) * t738 - mrSges(7,2) * t765 + mrSges(7,3) * t752 + t768 * t784 - t770 * t791;
t870 = -t733 * t852 + t856 * t734;
t719 = m(6) * t743 - mrSges(6,2) * t842 + mrSges(6,3) * t766 + t779 * t795 - t789 * t843 + t870;
t742 = t746 * t857 - t748 * t853;
t788 = -mrSges(6,2) * t843 + mrSges(6,3) * t795;
t739 = -pkin(5) * t842 - pkin(9) * t841 + t780 * t796 - t742;
t866 = -m(7) * t739 + t752 * mrSges(7,1) - mrSges(7,2) * t753 + t784 * t769 - t770 * t785;
t729 = m(6) * t742 + mrSges(6,1) * t842 - mrSges(6,3) * t767 - t779 * t796 + t788 * t843 + t866;
t713 = t853 * t719 + t857 * t729;
t801 = -mrSges(5,1) * t814 + mrSges(5,2) * t815;
t808 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t814;
t711 = m(5) * t749 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t807 + qJD(3) * t808 - t801 * t815 + t713;
t809 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t815;
t871 = t857 * t719 - t729 * t853;
t712 = m(5) * t750 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t806 - qJD(3) * t809 + t801 * t814 + t871;
t705 = t850 * t711 + t848 * t712;
t793 = Ifges(5,4) * t815 + Ifges(5,2) * t814 + Ifges(5,6) * qJD(3);
t794 = Ifges(5,1) * t815 + Ifges(5,4) * t814 + Ifges(5,5) * qJD(3);
t820 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t854 + Ifges(4,2) * t858) * qJD(1);
t821 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t854 + Ifges(4,4) * t858) * qJD(1);
t756 = Ifges(7,5) * t785 + Ifges(7,6) * t784 + Ifges(7,3) * t791;
t758 = Ifges(7,1) * t785 + Ifges(7,4) * t784 + Ifges(7,5) * t791;
t726 = -mrSges(7,1) * t739 + mrSges(7,3) * t738 + Ifges(7,4) * t753 + Ifges(7,2) * t752 + Ifges(7,6) * t765 - t756 * t785 + t758 * t791;
t757 = Ifges(7,4) * t785 + Ifges(7,2) * t784 + Ifges(7,6) * t791;
t727 = mrSges(7,2) * t739 - mrSges(7,3) * t737 + Ifges(7,1) * t753 + Ifges(7,4) * t752 + Ifges(7,5) * t765 + t756 * t784 - t757 * t791;
t772 = Ifges(6,4) * t796 + Ifges(6,2) * t795 + Ifges(6,6) * t843;
t773 = Ifges(6,1) * t796 + Ifges(6,4) * t795 + Ifges(6,5) * t843;
t864 = -mrSges(6,1) * t742 + mrSges(6,2) * t743 - Ifges(6,5) * t767 - Ifges(6,6) * t766 - Ifges(6,3) * t842 - pkin(5) * t866 - pkin(9) * t870 - t856 * t726 - t852 * t727 - t796 * t772 + t795 * t773;
t882 = mrSges(4,1) * t786 + mrSges(5,1) * t749 - mrSges(4,2) * t787 - mrSges(5,2) * t750 + Ifges(4,5) * t828 + Ifges(5,5) * t807 + Ifges(4,6) * t829 + Ifges(5,6) * t806 + pkin(3) * t705 + pkin(4) * t713 + (t820 * t854 - t821 * t858) * qJD(1) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t815 * t793 - t814 * t794 - t864;
t826 = (-mrSges(4,1) * t858 + mrSges(4,2) * t854) * qJD(1);
t878 = qJD(1) * t858;
t833 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t878;
t703 = m(4) * t786 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t828 + qJD(3) * t833 - t826 * t879 + t705;
t832 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t879;
t872 = -t711 * t848 + t850 * t712;
t704 = m(4) * t787 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t829 - qJD(3) * t832 + t826 * t878 + t872;
t873 = -t703 * t854 + t858 * t704;
t695 = m(3) * t805 - mrSges(3,1) * t860 - qJDD(1) * mrSges(3,2) + t873;
t722 = t856 * t733 + t852 * t734;
t867 = m(6) * t755 - t766 * mrSges(6,1) + t767 * mrSges(6,2) - t795 * t788 + t796 * t789 + t722;
t720 = m(5) * t783 - t806 * mrSges(5,1) + mrSges(5,2) * t807 - t814 * t808 + t809 * t815 + t867;
t797 = -pkin(7) * t860 + t868;
t862 = -m(4) * t797 + t829 * mrSges(4,1) - mrSges(4,2) * t828 - t832 * t879 + t833 * t878 - t720;
t715 = m(3) * t804 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t860 + t862;
t691 = t849 * t695 + t851 * t715;
t688 = m(2) * t834 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t860 + t691;
t874 = t851 * t695 - t715 * t849;
t689 = m(2) * t835 - mrSges(2,1) * t860 - qJDD(1) * mrSges(2,2) + t874;
t880 = t859 * t688 + t855 * t689;
t698 = t858 * t703 + t854 * t704;
t696 = m(3) * t847 + t698;
t875 = -t688 * t855 + t859 * t689;
t771 = Ifges(6,5) * t796 + Ifges(6,6) * t795 + Ifges(6,3) * t843;
t706 = mrSges(6,2) * t755 - mrSges(6,3) * t742 + Ifges(6,1) * t767 + Ifges(6,4) * t766 + Ifges(6,5) * t842 - pkin(9) * t722 - t726 * t852 + t727 * t856 + t771 * t795 - t772 * t843;
t863 = mrSges(7,1) * t737 - mrSges(7,2) * t738 + Ifges(7,5) * t753 + Ifges(7,6) * t752 + Ifges(7,3) * t765 + t757 * t785 - t758 * t784;
t707 = -mrSges(6,1) * t755 + mrSges(6,3) * t743 + Ifges(6,4) * t767 + Ifges(6,2) * t766 + Ifges(6,6) * t842 - pkin(5) * t722 - t771 * t796 + t773 * t843 - t863;
t792 = Ifges(5,5) * t815 + Ifges(5,6) * t814 + Ifges(5,3) * qJD(3);
t692 = -mrSges(5,1) * t783 + mrSges(5,3) * t750 + Ifges(5,4) * t807 + Ifges(5,2) * t806 + Ifges(5,6) * qJDD(3) - pkin(4) * t867 + pkin(8) * t871 + qJD(3) * t794 + t853 * t706 + t857 * t707 - t815 * t792;
t699 = mrSges(5,2) * t783 - mrSges(5,3) * t749 + Ifges(5,1) * t807 + Ifges(5,4) * t806 + Ifges(5,5) * qJDD(3) - pkin(8) * t713 - qJD(3) * t793 + t706 * t857 - t707 * t853 + t792 * t814;
t819 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t854 + Ifges(4,6) * t858) * qJD(1);
t681 = -mrSges(4,1) * t797 + mrSges(4,3) * t787 + Ifges(4,4) * t828 + Ifges(4,2) * t829 + Ifges(4,6) * qJDD(3) - pkin(3) * t720 + qJ(4) * t872 + qJD(3) * t821 + t850 * t692 + t848 * t699 - t819 * t879;
t683 = mrSges(4,2) * t797 - mrSges(4,3) * t786 + Ifges(4,1) * t828 + Ifges(4,4) * t829 + Ifges(4,5) * qJDD(3) - qJ(4) * t705 - qJD(3) * t820 - t692 * t848 + t699 * t850 + t819 * t878;
t865 = mrSges(2,1) * t834 + mrSges(3,1) * t804 - mrSges(2,2) * t835 - mrSges(3,2) * t805 + pkin(1) * t691 + pkin(2) * t862 + pkin(7) * t873 + t858 * t681 + t854 * t683 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t684 = -mrSges(3,1) * t847 + mrSges(3,3) * t805 + t860 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t698 - t882;
t679 = mrSges(3,2) * t847 - mrSges(3,3) * t804 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t860 - pkin(7) * t698 - t681 * t854 + t683 * t858;
t678 = -mrSges(2,2) * g(3) - mrSges(2,3) * t834 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t860 - qJ(2) * t691 + t679 * t851 - t684 * t849;
t677 = mrSges(2,1) * g(3) + mrSges(2,3) * t835 + t860 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t696 + qJ(2) * t874 + t849 * t679 + t851 * t684;
t1 = [-m(1) * g(1) + t875; -m(1) * g(2) + t880; (-m(1) - m(2)) * g(3) + t696; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t880 - t855 * t677 + t859 * t678; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t875 + t859 * t677 + t855 * t678; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t865; t865; t696; t882; t720; -t864; t863;];
tauJB  = t1;
