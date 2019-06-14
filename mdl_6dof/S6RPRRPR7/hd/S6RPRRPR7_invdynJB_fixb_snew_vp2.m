% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-05-05 23:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:01:18
% EndTime: 2019-05-05 23:01:28
% DurationCPUTime: 9.72s
% Computational Cost: add. (153947->339), mult. (331863->417), div. (0->0), fcn. (231141->10), ass. (0->137)
t852 = sin(qJ(1));
t856 = cos(qJ(1));
t830 = -t856 * g(1) - t852 * g(2);
t868 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t830;
t886 = 2 * qJD(5);
t885 = -pkin(1) - pkin(7);
t884 = mrSges(2,1) - mrSges(3,2);
t883 = Ifges(2,5) - Ifges(3,4);
t882 = (-Ifges(2,6) + Ifges(3,5));
t829 = t852 * g(1) - t856 * g(2);
t857 = qJD(1) ^ 2;
t867 = -t857 * qJ(2) + qJDD(2) - t829;
t804 = t885 * qJDD(1) + t867;
t851 = sin(qJ(3));
t855 = cos(qJ(3));
t795 = t851 * g(3) + t855 * t804;
t878 = qJD(1) * qJD(3);
t876 = t851 * t878;
t824 = t855 * qJDD(1) - t876;
t770 = (-t824 - t876) * pkin(8) + (-t851 * t855 * t857 + qJDD(3)) * pkin(3) + t795;
t796 = -t855 * g(3) + t851 * t804;
t823 = -t851 * qJDD(1) - t855 * t878;
t879 = qJD(1) * t855;
t828 = qJD(3) * pkin(3) - pkin(8) * t879;
t844 = t851 ^ 2;
t771 = -t844 * t857 * pkin(3) + t823 * pkin(8) - qJD(3) * t828 + t796;
t850 = sin(qJ(4));
t854 = cos(qJ(4));
t750 = t854 * t770 - t850 * t771;
t814 = (-t850 * t855 - t851 * t854) * qJD(1);
t781 = t814 * qJD(4) + t850 * t823 + t854 * t824;
t815 = (-t850 * t851 + t854 * t855) * qJD(1);
t837 = qJDD(3) + qJDD(4);
t838 = qJD(3) + qJD(4);
t736 = (t814 * t838 - t781) * qJ(5) + (t814 * t815 + t837) * pkin(4) + t750;
t751 = t850 * t770 + t854 * t771;
t780 = -t815 * qJD(4) + t854 * t823 - t850 * t824;
t801 = t838 * pkin(4) - t815 * qJ(5);
t810 = t814 ^ 2;
t738 = -t810 * pkin(4) + t780 * qJ(5) - t838 * t801 + t751;
t847 = sin(pkin(10));
t848 = cos(pkin(10));
t791 = t848 * t814 - t847 * t815;
t733 = t847 * t736 + t848 * t738 + t791 * t886;
t756 = t848 * t780 - t847 * t781;
t792 = t847 * t814 + t848 * t815;
t765 = -t791 * mrSges(6,1) + t792 * mrSges(6,2);
t783 = t838 * mrSges(6,1) - t792 * mrSges(6,3);
t766 = -t791 * pkin(5) - t792 * pkin(9);
t836 = t838 ^ 2;
t730 = -t836 * pkin(5) + t837 * pkin(9) + t791 * t766 + t733;
t774 = -t823 * pkin(3) + t828 * t879 + (-pkin(8) * t844 + t885) * t857 + t868;
t743 = -t780 * pkin(4) - t810 * qJ(5) + t815 * t801 + qJDD(5) + t774;
t757 = t847 * t780 + t848 * t781;
t734 = (-t791 * t838 - t757) * pkin(9) + (t792 * t838 - t756) * pkin(5) + t743;
t849 = sin(qJ(6));
t853 = cos(qJ(6));
t727 = -t849 * t730 + t853 * t734;
t778 = -t849 * t792 + t853 * t838;
t741 = t778 * qJD(6) + t853 * t757 + t849 * t837;
t755 = qJDD(6) - t756;
t779 = t853 * t792 + t849 * t838;
t758 = -t778 * mrSges(7,1) + t779 * mrSges(7,2);
t785 = qJD(6) - t791;
t759 = -t785 * mrSges(7,2) + t778 * mrSges(7,3);
t724 = m(7) * t727 + t755 * mrSges(7,1) - t741 * mrSges(7,3) - t779 * t758 + t785 * t759;
t728 = t853 * t730 + t849 * t734;
t740 = -t779 * qJD(6) - t849 * t757 + t853 * t837;
t760 = t785 * mrSges(7,1) - t779 * mrSges(7,3);
t725 = m(7) * t728 - t755 * mrSges(7,2) + t740 * mrSges(7,3) + t778 * t758 - t785 * t760;
t871 = -t849 * t724 + t853 * t725;
t711 = m(6) * t733 - t837 * mrSges(6,2) + t756 * mrSges(6,3) + t791 * t765 - t838 * t783 + t871;
t870 = -t848 * t736 + t847 * t738;
t732 = -0.2e1 * qJD(5) * t792 - t870;
t782 = -t838 * mrSges(6,2) + t791 * mrSges(6,3);
t729 = -t837 * pkin(5) - t836 * pkin(9) + (t886 + t766) * t792 + t870;
t865 = -m(7) * t729 + t740 * mrSges(7,1) - t741 * mrSges(7,2) + t778 * t759 - t779 * t760;
t720 = m(6) * t732 + t837 * mrSges(6,1) - t757 * mrSges(6,3) - t792 * t765 + t838 * t782 + t865;
t704 = t847 * t711 + t848 * t720;
t793 = -t814 * mrSges(5,1) + t815 * mrSges(5,2);
t800 = -t838 * mrSges(5,2) + t814 * mrSges(5,3);
t701 = m(5) * t750 + t837 * mrSges(5,1) - t781 * mrSges(5,3) - t815 * t793 + t838 * t800 + t704;
t802 = t838 * mrSges(5,1) - t815 * mrSges(5,3);
t872 = t848 * t711 - t847 * t720;
t702 = m(5) * t751 - t837 * mrSges(5,2) + t780 * mrSges(5,3) + t814 * t793 - t838 * t802 + t872;
t695 = t854 * t701 + t850 * t702;
t822 = (mrSges(4,1) * t851 + mrSges(4,2) * t855) * qJD(1);
t880 = qJD(1) * t851;
t826 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t880;
t692 = m(4) * t795 + qJDD(3) * mrSges(4,1) - t824 * mrSges(4,3) + qJD(3) * t826 - t822 * t879 + t695;
t827 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t879;
t873 = -t850 * t701 + t854 * t702;
t693 = m(4) * t796 - qJDD(3) * mrSges(4,2) + t823 * mrSges(4,3) - qJD(3) * t827 - t822 * t880 + t873;
t689 = t855 * t692 + t851 * t693;
t809 = -qJDD(1) * pkin(1) + t867;
t866 = -m(3) * t809 + (t857 * mrSges(3,3)) - t689;
t685 = m(2) * t829 - (t857 * mrSges(2,2)) + t884 * qJDD(1) + t866;
t807 = t857 * pkin(1) - t868;
t803 = t885 * t857 + t868;
t714 = t853 * t724 + t849 * t725;
t712 = m(6) * t743 - t756 * mrSges(6,1) + t757 * mrSges(6,2) - t791 * t782 + t792 * t783 + t714;
t862 = m(5) * t774 - t780 * mrSges(5,1) + t781 * mrSges(5,2) - t814 * t800 + t815 * t802 + t712;
t861 = -m(4) * t803 + t823 * mrSges(4,1) - t824 * mrSges(4,2) - t826 * t880 - t827 * t879 - t862;
t859 = -m(3) * t807 + (t857 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t861;
t707 = m(2) * t830 - (t857 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t859;
t881 = t856 * t685 + t852 * t707;
t875 = -t852 * t685 + t856 * t707;
t874 = -t851 * t692 + t855 * t693;
t744 = Ifges(7,5) * t779 + Ifges(7,6) * t778 + Ifges(7,3) * t785;
t746 = Ifges(7,1) * t779 + Ifges(7,4) * t778 + Ifges(7,5) * t785;
t717 = -mrSges(7,1) * t729 + mrSges(7,3) * t728 + Ifges(7,4) * t741 + Ifges(7,2) * t740 + Ifges(7,6) * t755 - t779 * t744 + t785 * t746;
t745 = Ifges(7,4) * t779 + Ifges(7,2) * t778 + Ifges(7,6) * t785;
t718 = mrSges(7,2) * t729 - mrSges(7,3) * t727 + Ifges(7,1) * t741 + Ifges(7,4) * t740 + Ifges(7,5) * t755 + t778 * t744 - t785 * t745;
t761 = Ifges(6,5) * t792 + Ifges(6,6) * t791 + Ifges(6,3) * t838;
t762 = Ifges(6,4) * t792 + Ifges(6,2) * t791 + Ifges(6,6) * t838;
t696 = mrSges(6,2) * t743 - mrSges(6,3) * t732 + Ifges(6,1) * t757 + Ifges(6,4) * t756 + Ifges(6,5) * t837 - pkin(9) * t714 - t849 * t717 + t853 * t718 + t791 * t761 - t838 * t762;
t763 = Ifges(6,1) * t792 + Ifges(6,4) * t791 + Ifges(6,5) * t838;
t863 = mrSges(7,1) * t727 - mrSges(7,2) * t728 + Ifges(7,5) * t741 + Ifges(7,6) * t740 + Ifges(7,3) * t755 + t779 * t745 - t778 * t746;
t697 = -mrSges(6,1) * t743 + mrSges(6,3) * t733 + Ifges(6,4) * t757 + Ifges(6,2) * t756 + Ifges(6,6) * t837 - pkin(5) * t714 - t792 * t761 + t838 * t763 - t863;
t786 = Ifges(5,5) * t815 + Ifges(5,6) * t814 + Ifges(5,3) * t838;
t788 = Ifges(5,1) * t815 + Ifges(5,4) * t814 + Ifges(5,5) * t838;
t683 = -mrSges(5,1) * t774 + mrSges(5,3) * t751 + Ifges(5,4) * t781 + Ifges(5,2) * t780 + Ifges(5,6) * t837 - pkin(4) * t712 + qJ(5) * t872 + t847 * t696 + t848 * t697 - t815 * t786 + t838 * t788;
t787 = Ifges(5,4) * t815 + Ifges(5,2) * t814 + Ifges(5,6) * t838;
t690 = mrSges(5,2) * t774 - mrSges(5,3) * t750 + Ifges(5,1) * t781 + Ifges(5,4) * t780 + Ifges(5,5) * t837 - qJ(5) * t704 + t848 * t696 - t847 * t697 + t814 * t786 - t838 * t787;
t811 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t855 - Ifges(4,6) * t851) * qJD(1);
t813 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t855 - Ifges(4,4) * t851) * qJD(1);
t680 = -mrSges(4,1) * t803 + mrSges(4,3) * t796 + Ifges(4,4) * t824 + Ifges(4,2) * t823 + Ifges(4,6) * qJDD(3) - pkin(3) * t862 + pkin(8) * t873 + qJD(3) * t813 + t854 * t683 + t850 * t690 - t811 * t879;
t812 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t855 - Ifges(4,2) * t851) * qJD(1);
t682 = mrSges(4,2) * t803 - mrSges(4,3) * t795 + Ifges(4,1) * t824 + Ifges(4,4) * t823 + Ifges(4,5) * qJDD(3) - pkin(8) * t695 - qJD(3) * t812 - t850 * t683 + t854 * t690 - t811 * t880;
t687 = qJDD(1) * mrSges(3,2) - t866;
t864 = mrSges(2,1) * t829 - mrSges(2,2) * t830 + mrSges(3,2) * t809 - mrSges(3,3) * t807 - pkin(1) * t687 - pkin(7) * t689 + qJ(2) * t859 - t851 * t680 + t855 * t682 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t860 = mrSges(5,1) * t750 + mrSges(6,1) * t732 - mrSges(5,2) * t751 - mrSges(6,2) * t733 + pkin(4) * t704 + pkin(5) * t865 + pkin(9) * t871 + t853 * t717 + t849 * t718 + t792 * t762 - t791 * t763 - t814 * t788 + Ifges(6,6) * t756 + Ifges(6,5) * t757 + t815 * t787 + Ifges(5,6) * t780 + Ifges(5,5) * t781 + (Ifges(6,3) + Ifges(5,3)) * t837;
t858 = mrSges(4,1) * t795 - mrSges(4,2) * t796 + Ifges(4,5) * t824 + Ifges(4,6) * t823 + Ifges(4,3) * qJDD(3) + pkin(3) * t695 + t812 * t879 + t813 * t880 + t860;
t688 = -m(3) * g(3) + t874;
t679 = pkin(2) * t689 - qJ(2) * t688 + t858 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t882 * t857) + t883 * qJDD(1) - mrSges(2,3) * t829 + mrSges(3,1) * t809;
t678 = -mrSges(3,1) * t807 + mrSges(2,3) * t830 - pkin(1) * t688 - pkin(2) * t861 - pkin(7) * t874 + t884 * g(3) - t882 * qJDD(1) - t855 * t680 - t851 * t682 + t883 * t857;
t1 = [-m(1) * g(1) + t875; -m(1) * g(2) + t881; (-m(1) - m(2) - m(3)) * g(3) + t874; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t881 - t852 * t678 + t856 * t679; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t875 + t856 * t678 + t852 * t679; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t864; t864; t687; t858; t860; t712; t863;];
tauJB  = t1;
