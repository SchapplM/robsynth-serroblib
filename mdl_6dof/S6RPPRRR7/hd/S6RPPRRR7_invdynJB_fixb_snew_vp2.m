% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-05-05 16:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:05:52
% EndTime: 2019-05-05 16:06:00
% DurationCPUTime: 8.15s
% Computational Cost: add. (131129->314), mult. (301590->382), div. (0->0), fcn. (223194->10), ass. (0->137)
t853 = sin(qJ(1));
t857 = cos(qJ(1));
t825 = t853 * g(1) - t857 * g(2);
t858 = qJD(1) ^ 2;
t869 = -t858 * qJ(2) + qJDD(2) - t825;
t891 = -pkin(1) - qJ(3);
t899 = -(2 * qJD(1) * qJD(3)) + t891 * qJDD(1) + t869;
t826 = -t857 * g(1) - t853 * g(2);
t898 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t826;
t816 = t858 * pkin(1) - t898;
t897 = -m(3) * t816 + t858 * mrSges(3,2) + qJDD(1) * mrSges(3,3);
t895 = pkin(3) * t858;
t894 = mrSges(2,1) - mrSges(3,2);
t893 = -Ifges(3,4) + Ifges(2,5);
t892 = -Ifges(2,6) + Ifges(3,5);
t849 = cos(pkin(10));
t890 = mrSges(4,2) * t849;
t848 = sin(pkin(10));
t804 = t848 * g(3) + t899 * t849;
t784 = (-pkin(7) * qJDD(1) - t848 * t895) * t849 + t804;
t805 = -t849 * g(3) + t899 * t848;
t838 = t848 ^ 2;
t884 = qJDD(1) * t848;
t787 = -pkin(7) * t884 - t838 * t895 + t805;
t852 = sin(qJ(4));
t856 = cos(qJ(4));
t769 = t856 * t784 - t852 * t787;
t872 = -t848 * t852 + t849 * t856;
t873 = -t848 * t856 - t849 * t852;
t820 = t873 * qJD(1);
t886 = t820 * qJD(4);
t807 = t872 * qJDD(1) + t886;
t821 = t872 * qJD(1);
t749 = (-t807 + t886) * pkin(8) + (t820 * t821 + qJDD(4)) * pkin(4) + t769;
t770 = t852 * t784 + t856 * t787;
t806 = -t821 * qJD(4) + t873 * qJDD(1);
t814 = qJD(4) * pkin(4) - t821 * pkin(8);
t819 = t820 ^ 2;
t751 = -t819 * pkin(4) + t806 * pkin(8) - qJD(4) * t814 + t770;
t851 = sin(qJ(5));
t855 = cos(qJ(5));
t747 = t851 * t749 + t855 * t751;
t797 = t851 * t820 + t855 * t821;
t766 = -t797 * qJD(5) + t855 * t806 - t851 * t807;
t796 = t855 * t820 - t851 * t821;
t778 = -t796 * mrSges(6,1) + t797 * mrSges(6,2);
t840 = qJD(4) + qJD(5);
t789 = t840 * mrSges(6,1) - t797 * mrSges(6,3);
t837 = qJDD(4) + qJDD(5);
t779 = -t796 * pkin(5) - t797 * pkin(9);
t836 = t840 ^ 2;
t743 = -t836 * pkin(5) + t837 * pkin(9) + t796 * t779 + t747;
t868 = qJDD(3) + t898;
t888 = -t849 ^ 2 - t838;
t791 = pkin(3) * t884 + (t888 * pkin(7) + t891) * t858 + t868;
t760 = -t806 * pkin(4) - t819 * pkin(8) + t821 * t814 + t791;
t767 = t796 * qJD(5) + t851 * t806 + t855 * t807;
t744 = (-t796 * t840 - t767) * pkin(9) + (t797 * t840 - t766) * pkin(5) + t760;
t850 = sin(qJ(6));
t854 = cos(qJ(6));
t740 = -t850 * t743 + t854 * t744;
t785 = -t850 * t797 + t854 * t840;
t754 = t785 * qJD(6) + t854 * t767 + t850 * t837;
t765 = qJDD(6) - t766;
t786 = t854 * t797 + t850 * t840;
t771 = -t785 * mrSges(7,1) + t786 * mrSges(7,2);
t792 = qJD(6) - t796;
t772 = -t792 * mrSges(7,2) + t785 * mrSges(7,3);
t737 = m(7) * t740 + t765 * mrSges(7,1) - t754 * mrSges(7,3) - t786 * t771 + t792 * t772;
t741 = t854 * t743 + t850 * t744;
t753 = -t786 * qJD(6) - t850 * t767 + t854 * t837;
t773 = t792 * mrSges(7,1) - t786 * mrSges(7,3);
t738 = m(7) * t741 - t765 * mrSges(7,2) + t753 * mrSges(7,3) + t785 * t771 - t792 * t773;
t877 = -t850 * t737 + t854 * t738;
t725 = m(6) * t747 - t837 * mrSges(6,2) + t766 * mrSges(6,3) + t796 * t778 - t840 * t789 + t877;
t746 = t855 * t749 - t851 * t751;
t788 = -t840 * mrSges(6,2) + t796 * mrSges(6,3);
t742 = -t837 * pkin(5) - t836 * pkin(9) + t797 * t779 - t746;
t866 = -m(7) * t742 + t753 * mrSges(7,1) - t754 * mrSges(7,2) + t785 * t772 - t786 * t773;
t733 = m(6) * t746 + t837 * mrSges(6,1) - t767 * mrSges(6,3) - t797 * t778 + t840 * t788 + t866;
t717 = t851 * t725 + t855 * t733;
t800 = -t820 * mrSges(5,1) + t821 * mrSges(5,2);
t812 = -qJD(4) * mrSges(5,2) + t820 * mrSges(5,3);
t714 = m(5) * t769 + qJDD(4) * mrSges(5,1) - t807 * mrSges(5,3) + qJD(4) * t812 - t821 * t800 + t717;
t813 = qJD(4) * mrSges(5,1) - t821 * mrSges(5,3);
t878 = t855 * t725 - t851 * t733;
t715 = m(5) * t770 - qJDD(4) * mrSges(5,2) + t806 * mrSges(5,3) - qJD(4) * t813 + t820 * t800 + t878;
t708 = t856 * t714 + t852 * t715;
t871 = -mrSges(4,3) * qJDD(1) - t858 * (mrSges(4,1) * t848 + t890);
t706 = m(4) * t804 + t871 * t849 + t708;
t879 = -t852 * t714 + t856 * t715;
t707 = m(4) * t805 + t871 * t848 + t879;
t703 = t849 * t706 + t848 * t707;
t818 = -qJDD(1) * pkin(1) + t869;
t867 = -m(3) * t818 + t858 * mrSges(3,3) - t703;
t699 = m(2) * t825 - t858 * mrSges(2,2) + t894 * qJDD(1) + t867;
t811 = t891 * t858 + t868;
t727 = t854 * t737 + t850 * t738;
t865 = m(6) * t760 - t766 * mrSges(6,1) + t767 * mrSges(6,2) - t796 * t788 + t797 * t789 + t727;
t861 = m(5) * t791 - t806 * mrSges(5,1) + t807 * mrSges(5,2) - t820 * t812 + t821 * t813 + t865;
t860 = m(4) * t811 + mrSges(4,1) * t884 + qJDD(1) * t890 + t861;
t882 = t888 * mrSges(4,3);
t720 = t860 + m(2) * t826 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) + t882) * t858 + t897;
t889 = t857 * t699 + t853 * t720;
t874 = Ifges(4,5) * t849 - Ifges(4,6) * t848;
t887 = t858 * t874;
t881 = -t853 * t699 + t857 * t720;
t880 = -t848 * t706 + t849 * t707;
t876 = Ifges(4,1) * t849 - Ifges(4,4) * t848;
t875 = Ifges(4,4) * t849 - Ifges(4,2) * t848;
t755 = Ifges(7,5) * t786 + Ifges(7,6) * t785 + Ifges(7,3) * t792;
t757 = Ifges(7,1) * t786 + Ifges(7,4) * t785 + Ifges(7,5) * t792;
t730 = -mrSges(7,1) * t742 + mrSges(7,3) * t741 + Ifges(7,4) * t754 + Ifges(7,2) * t753 + Ifges(7,6) * t765 - t786 * t755 + t792 * t757;
t756 = Ifges(7,4) * t786 + Ifges(7,2) * t785 + Ifges(7,6) * t792;
t731 = mrSges(7,2) * t742 - mrSges(7,3) * t740 + Ifges(7,1) * t754 + Ifges(7,4) * t753 + Ifges(7,5) * t765 + t785 * t755 - t792 * t756;
t775 = Ifges(6,4) * t797 + Ifges(6,2) * t796 + Ifges(6,6) * t840;
t776 = Ifges(6,1) * t797 + Ifges(6,4) * t796 + Ifges(6,5) * t840;
t864 = mrSges(6,1) * t746 - mrSges(6,2) * t747 + Ifges(6,5) * t767 + Ifges(6,6) * t766 + Ifges(6,3) * t837 + pkin(5) * t866 + pkin(9) * t877 + t854 * t730 + t850 * t731 + t797 * t775 - t796 * t776;
t774 = Ifges(6,5) * t797 + Ifges(6,6) * t796 + Ifges(6,3) * t840;
t709 = mrSges(6,2) * t760 - mrSges(6,3) * t746 + Ifges(6,1) * t767 + Ifges(6,4) * t766 + Ifges(6,5) * t837 - pkin(9) * t727 - t850 * t730 + t854 * t731 + t796 * t774 - t840 * t775;
t862 = mrSges(7,1) * t740 - mrSges(7,2) * t741 + Ifges(7,5) * t754 + Ifges(7,6) * t753 + Ifges(7,3) * t765 + t786 * t756 - t785 * t757;
t710 = -mrSges(6,1) * t760 + mrSges(6,3) * t747 + Ifges(6,4) * t767 + Ifges(6,2) * t766 + Ifges(6,6) * t837 - pkin(5) * t727 - t797 * t774 + t840 * t776 - t862;
t793 = Ifges(5,5) * t821 + Ifges(5,6) * t820 + Ifges(5,3) * qJD(4);
t795 = Ifges(5,1) * t821 + Ifges(5,4) * t820 + Ifges(5,5) * qJD(4);
t697 = -mrSges(5,1) * t791 + mrSges(5,3) * t770 + Ifges(5,4) * t807 + Ifges(5,2) * t806 + Ifges(5,6) * qJDD(4) - pkin(4) * t865 + pkin(8) * t878 + qJD(4) * t795 + t851 * t709 + t855 * t710 - t821 * t793;
t794 = Ifges(5,4) * t821 + Ifges(5,2) * t820 + Ifges(5,6) * qJD(4);
t704 = mrSges(5,2) * t791 - mrSges(5,3) * t769 + Ifges(5,1) * t807 + Ifges(5,4) * t806 + Ifges(5,5) * qJDD(4) - pkin(8) * t717 - qJD(4) * t794 + t855 * t709 - t851 * t710 + t820 * t793;
t694 = -mrSges(4,1) * t811 + mrSges(4,3) * t805 - pkin(3) * t861 + pkin(7) * t879 + t875 * qJDD(1) + t856 * t697 + t852 * t704 - t849 * t887;
t696 = mrSges(4,2) * t811 - mrSges(4,3) * t804 - pkin(7) * t708 + t876 * qJDD(1) - t852 * t697 + t856 * t704 - t848 * t887;
t701 = qJDD(1) * mrSges(3,2) - t867;
t722 = t858 * t882 + t860;
t863 = -mrSges(2,2) * t826 - mrSges(3,3) * t816 - pkin(1) * t701 - qJ(3) * t703 - t848 * t694 + t849 * t696 + qJ(2) * (t722 + t897) + mrSges(3,2) * t818 + mrSges(2,1) * t825 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t859 = mrSges(5,1) * t769 - mrSges(5,2) * t770 + Ifges(5,5) * t807 + Ifges(5,6) * t806 + Ifges(5,3) * qJDD(4) + pkin(4) * t717 + t821 * t794 - t820 * t795 + t864;
t702 = -m(3) * g(3) + t880;
t693 = t859 - mrSges(2,3) * t825 + mrSges(3,1) * t818 + mrSges(4,1) * t804 - mrSges(4,2) * t805 + pkin(3) * t708 + pkin(2) * t703 - qJ(2) * t702 + (t874 + t893) * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t848 * t876 + t849 * t875 + t892) * t858;
t692 = -mrSges(3,1) * t816 + mrSges(2,3) * t826 - pkin(1) * t702 + pkin(2) * t722 + t894 * g(3) - qJ(3) * t880 - t892 * qJDD(1) - t849 * t694 - t848 * t696 + t893 * t858;
t1 = [-m(1) * g(1) + t881; -m(1) * g(2) + t889; (-m(1) - m(2) - m(3)) * g(3) + t880; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t889 - t853 * t692 + t857 * t693; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t881 + t857 * t692 + t853 * t693; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t863; t863; t701; t722; t859; t864; t862;];
tauJB  = t1;
