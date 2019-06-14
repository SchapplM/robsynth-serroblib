% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 04:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPPR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:31:15
% EndTime: 2019-05-07 04:31:23
% DurationCPUTime: 5.03s
% Computational Cost: add. (50299->349), mult. (105192->406), div. (0->0), fcn. (68496->8), ass. (0->135)
t852 = Ifges(5,1) + Ifges(4,1) + Ifges(6,2);
t851 = Ifges(6,1) + Ifges(5,3) + Ifges(4,2);
t832 = Ifges(4,5) + Ifges(5,4) + Ifges(6,6);
t831 = -Ifges(5,5) + Ifges(4,4) + Ifges(6,4);
t830 = Ifges(6,5) + Ifges(4,6) - Ifges(5,6);
t850 = -Ifges(5,2) - Ifges(4,3) - Ifges(6,3);
t797 = qJD(2) + qJD(3);
t849 = t797 ^ 2;
t848 = -2 * qJD(4);
t847 = 2 * qJD(4);
t846 = -2 * qJD(5);
t845 = -pkin(4) - pkin(9);
t844 = cos(qJ(3));
t843 = (pkin(3) * t797);
t842 = -mrSges(4,3) - mrSges(5,2);
t803 = sin(qJ(3));
t807 = cos(qJ(2));
t834 = qJD(1) * t807;
t804 = sin(qJ(2));
t835 = qJD(1) * t804;
t768 = t803 * t835 - t844 * t834;
t841 = t768 * t797;
t809 = qJD(1) ^ 2;
t840 = t807 ^ 2 * t809;
t805 = sin(qJ(1));
t808 = cos(qJ(1));
t784 = -g(1) * t808 - g(2) * t805;
t771 = -pkin(1) * t809 + qJDD(1) * pkin(7) + t784;
t839 = t804 * t771;
t833 = qJD(1) * qJD(2);
t778 = qJDD(1) * t804 + t807 * t833;
t705 = qJDD(2) * pkin(2) - t778 * pkin(8) - t839 + (pkin(2) * t804 * t809 + pkin(8) * t833 - g(3)) * t807;
t751 = -g(3) * t804 + t807 * t771;
t779 = qJDD(1) * t807 - t804 * t833;
t782 = qJD(2) * pkin(2) - pkin(8) * t835;
t706 = -pkin(2) * t840 + pkin(8) * t779 - qJD(2) * t782 + t751;
t692 = t844 * t705 - t803 * t706;
t720 = -t768 * qJD(3) + t844 * t778 + t803 * t779;
t769 = (t803 * t807 + t844 * t804) * qJD(1);
t739 = mrSges(6,1) * t769 + mrSges(6,2) * t768;
t753 = -mrSges(4,2) * t797 - mrSges(4,3) * t768;
t796 = qJDD(2) + qJDD(3);
t740 = pkin(3) * t768 - qJ(4) * t769;
t819 = -qJ(4) * t849 + t769 * t740 + qJDD(4) - t692;
t691 = -t796 * pkin(3) + t819;
t758 = -mrSges(5,2) * t768 + mrSges(5,3) * t797;
t812 = (-t720 - t841) * qJ(5) + t819 + (t768 * pkin(4) + t846) * t769;
t685 = (-pkin(3) - pkin(4)) * t796 + t812;
t752 = -mrSges(6,1) * t797 - mrSges(6,3) * t768;
t719 = t769 * qJD(3) + t803 * t778 - t844 * t779;
t754 = -pkin(4) * t797 - qJ(5) * t769;
t764 = t768 ^ 2;
t783 = t805 * g(1) - t808 * g(2);
t770 = -qJDD(1) * pkin(1) - t809 * pkin(7) - t783;
t721 = -t779 * pkin(2) - pkin(8) * t840 + t782 * t835 + t770;
t818 = t719 * pkin(3) + t721 + (-t720 + t841) * qJ(4);
t813 = -t764 * qJ(5) + qJDD(5) - t818 + (t754 + t847) * t769;
t679 = t813 + t845 * t719 + (-pkin(5) * t768 + (-pkin(3) - pkin(9)) * t769) * t797 + t720 * pkin(5);
t743 = pkin(5) * t769 - pkin(9) * t768;
t680 = -t849 * pkin(5) - t769 * t743 + (-pkin(3) + t845) * t796 + t812;
t802 = sin(qJ(6));
t806 = cos(qJ(6));
t677 = t679 * t806 - t680 * t802;
t748 = -t768 * t802 - t797 * t806;
t696 = qJD(6) * t748 + t719 * t806 - t796 * t802;
t749 = t768 * t806 - t797 * t802;
t704 = -mrSges(7,1) * t748 + mrSges(7,2) * t749;
t717 = qJDD(6) + t720;
t763 = qJD(6) + t769;
t722 = -mrSges(7,2) * t763 + mrSges(7,3) * t748;
t675 = m(7) * t677 + mrSges(7,1) * t717 - mrSges(7,3) * t696 - t704 * t749 + t722 * t763;
t678 = t679 * t802 + t680 * t806;
t695 = -qJD(6) * t749 - t719 * t802 - t796 * t806;
t723 = mrSges(7,1) * t763 - mrSges(7,3) * t749;
t676 = m(7) * t678 - mrSges(7,2) * t717 + mrSges(7,3) * t695 + t704 * t748 - t723 * t763;
t837 = -t802 * t675 + t806 * t676;
t823 = m(6) * t685 + t796 * mrSges(6,2) + t797 * t752 + t837;
t817 = -m(5) * t691 + t796 * mrSges(5,1) + t797 * t758 - t823;
t741 = mrSges(5,1) * t768 - mrSges(5,3) * t769;
t836 = -mrSges(4,1) * t768 - mrSges(4,2) * t769 - t741;
t662 = m(4) * t692 + t796 * mrSges(4,1) + t797 * t753 + (t739 + t836) * t769 + (mrSges(6,3) + t842) * t720 + t817;
t693 = t803 * t705 + t844 * t706;
t755 = mrSges(4,1) * t797 - mrSges(4,3) * t769;
t757 = mrSges(6,2) * t797 - mrSges(6,3) * t769;
t786 = t797 * t847;
t822 = pkin(3) * t849 - t796 * qJ(4) + t768 * t740 - t693;
t690 = t786 - t822;
t756 = -mrSges(5,1) * t797 + mrSges(5,2) * t769;
t816 = t764 * pkin(4) - t719 * qJ(5) + t822;
t686 = t768 * t846 + (t848 - t754) * t797 + t816;
t682 = t796 * pkin(5) - t849 * pkin(9) + t797 * t754 + t786 + ((2 * qJD(5)) + t743) * t768 - t816;
t820 = -m(7) * t682 + t695 * mrSges(7,1) - t696 * mrSges(7,2) + t748 * t722 - t749 * t723;
t815 = -m(6) * t686 + t719 * mrSges(6,3) + t768 * t739 - t820;
t814 = m(5) * t690 + t796 * mrSges(5,3) + t797 * t756 + t815;
t668 = t814 + t836 * t768 + m(4) * t693 + (-mrSges(4,2) + mrSges(6,1)) * t796 + (-t755 + t757) * t797 + t842 * t719;
t658 = t844 * t662 + t803 * t668;
t750 = -t807 * g(3) - t839;
t777 = (-mrSges(3,1) * t807 + mrSges(3,2) * t804) * qJD(1);
t781 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t834;
t656 = m(3) * t750 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t778 + qJD(2) * t781 - t777 * t835 + t658;
t780 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t835;
t824 = -t662 * t803 + t844 * t668;
t657 = m(3) * t751 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t779 - qJD(2) * t780 + t777 * t834 + t824;
t825 = -t656 * t804 + t807 * t657;
t651 = m(2) * t784 - mrSges(2,1) * t809 - qJDD(1) * mrSges(2,2) + t825;
t688 = (t848 + t843) * t769 + t818;
t665 = t806 * t675 + t802 * t676;
t684 = -t719 * pkin(4) - t769 * t843 + t813;
t821 = -m(6) * t684 - t720 * mrSges(6,1) - t719 * mrSges(6,2) - t768 * t752 - t769 * t757 - t665;
t663 = m(5) * t688 + t719 * mrSges(5,1) - t720 * mrSges(5,3) - t769 * t756 + t768 * t758 + t821;
t811 = m(4) * t721 + t719 * mrSges(4,1) + t720 * mrSges(4,2) + t768 * t753 + t769 * t755 + t663;
t810 = -m(3) * t770 + t779 * mrSges(3,1) - t778 * mrSges(3,2) - t780 * t835 + t781 * t834 - t811;
t660 = m(2) * t783 + qJDD(1) * mrSges(2,1) - t809 * mrSges(2,2) + t810;
t838 = t805 * t651 + t808 * t660;
t652 = t807 * t656 + t804 * t657;
t829 = t851 * t768 - t831 * t769 - t830 * t797;
t828 = t830 * t768 - t832 * t769 + t850 * t797;
t827 = t831 * t768 - t852 * t769 - t832 * t797;
t826 = t808 * t651 - t660 * t805;
t767 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t804 + Ifges(3,4) * t807) * qJD(1);
t766 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t804 + Ifges(3,2) * t807) * qJD(1);
t765 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t804 + Ifges(3,6) * t807) * qJD(1);
t699 = Ifges(7,1) * t749 + Ifges(7,4) * t748 + Ifges(7,5) * t763;
t698 = Ifges(7,4) * t749 + Ifges(7,2) * t748 + Ifges(7,6) * t763;
t697 = Ifges(7,5) * t749 + Ifges(7,6) * t748 + Ifges(7,3) * t763;
t670 = mrSges(7,2) * t682 - mrSges(7,3) * t677 + Ifges(7,1) * t696 + Ifges(7,4) * t695 + Ifges(7,5) * t717 + t697 * t748 - t698 * t763;
t669 = -mrSges(7,1) * t682 + mrSges(7,3) * t678 + Ifges(7,4) * t696 + Ifges(7,2) * t695 + Ifges(7,6) * t717 - t697 * t749 + t699 * t763;
t664 = -t720 * mrSges(6,3) - t769 * t739 + t823;
t648 = mrSges(6,1) * t684 - mrSges(6,3) * t685 - mrSges(5,3) * t688 + Ifges(7,3) * t717 + pkin(5) * t665 + mrSges(5,2) * t691 - mrSges(4,3) * t692 + mrSges(7,1) * t677 - mrSges(7,2) * t678 - qJ(4) * t663 - qJ(5) * t664 + Ifges(7,6) * t695 + Ifges(7,5) * t696 + mrSges(4,2) * t721 - t748 * t699 + t749 * t698 + t829 * t797 + t832 * t796 + t828 * t768 + t852 * t720 - t831 * t719;
t647 = -mrSges(6,2) * t684 + mrSges(6,3) * t686 - mrSges(5,1) * t688 + pkin(9) * t665 + mrSges(5,2) * t690 + mrSges(4,3) * t693 - pkin(3) * t663 - t806 * t670 + t802 * t669 - mrSges(4,1) * t721 - pkin(4) * t821 - qJ(5) * t815 + (-qJ(5) * t757 - t827) * t797 + (-mrSges(6,1) * qJ(5) + t830) * t796 + t828 * t769 + t831 * t720 - t851 * t719;
t646 = -mrSges(6,2) * t685 + mrSges(6,1) * t686 + pkin(9) * t837 + (-mrSges(6,1) * qJ(4) + t850) * t796 - pkin(1) * t652 + (mrSges(5,2) * qJ(4) + t830) * t719 + (-pkin(3) * (-mrSges(5,2) + mrSges(6,3)) - t832) * t720 - pkin(2) * t658 + (qJ(4) * t741 + t827) * t768 + (-pkin(3) * (t739 - t741) + t829) * t769 - mrSges(5,3) * t690 + mrSges(5,1) * t691 - mrSges(4,1) * t692 + mrSges(4,2) * t693 - qJ(4) * (t797 * t757 + t814) - pkin(3) * t817 + pkin(5) * t820 + pkin(4) * t664 - Ifges(3,3) * qJDD(2) + mrSges(2,1) * g(3) + (-t766 * t804 + t767 * t807) * qJD(1) + t806 * t669 + t809 * Ifges(2,5) + t802 * t670 - mrSges(3,1) * t750 + mrSges(3,2) * t751 - Ifges(3,5) * t778 - Ifges(3,6) * t779 + mrSges(2,3) * t784 + Ifges(2,6) * qJDD(1);
t645 = mrSges(3,2) * t770 - mrSges(3,3) * t750 + Ifges(3,1) * t778 + Ifges(3,4) * t779 + Ifges(3,5) * qJDD(2) - pkin(8) * t658 - qJD(2) * t766 - t803 * t647 + t844 * t648 + t765 * t834;
t644 = -mrSges(3,1) * t770 + mrSges(3,3) * t751 + Ifges(3,4) * t778 + Ifges(3,2) * t779 + Ifges(3,6) * qJDD(2) - pkin(2) * t811 + pkin(8) * t824 + qJD(2) * t767 + t844 * t647 + t803 * t648 - t765 * t835;
t643 = -mrSges(2,2) * g(3) - mrSges(2,3) * t783 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t809 - pkin(7) * t652 - t644 * t804 + t645 * t807;
t1 = [-m(1) * g(1) + t826; -m(1) * g(2) + t838; (-m(1) - m(2)) * g(3) + t652; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t838 + t808 * t643 - t805 * t646; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t826 + t805 * t643 + t808 * t646; -mrSges(1,1) * g(2) + mrSges(2,1) * t783 + mrSges(1,2) * g(1) - mrSges(2,2) * t784 + Ifges(2,3) * qJDD(1) + pkin(1) * t810 + pkin(7) * t825 + t807 * t644 + t804 * t645;];
tauB  = t1;
