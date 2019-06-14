% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-05-04 23:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPRRP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:54:40
% EndTime: 2019-05-04 23:54:46
% DurationCPUTime: 4.16s
% Computational Cost: add. (50136->270), mult. (92140->323), div. (0->0), fcn. (57517->10), ass. (0->122)
t825 = Ifges(6,1) + Ifges(7,1);
t816 = Ifges(6,4) + Ifges(7,4);
t814 = Ifges(6,5) + Ifges(7,5);
t824 = Ifges(6,2) + Ifges(7,2);
t813 = Ifges(6,6) + Ifges(7,6);
t823 = Ifges(6,3) + Ifges(7,3);
t769 = sin(pkin(10));
t771 = cos(pkin(10));
t755 = t769 * g(1) - t771 * g(2);
t756 = -t771 * g(1) - t769 * g(2);
t766 = -g(3) + qJDD(1);
t778 = cos(qJ(2));
t772 = cos(pkin(6));
t775 = sin(qJ(2));
t809 = t772 * t775;
t770 = sin(pkin(6));
t810 = t770 * t775;
t706 = t755 * t809 + t778 * t756 + t766 * t810;
t822 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t706;
t705 = -t775 * t756 + (t755 * t772 + t766 * t770) * t778;
t773 = sin(qJ(5));
t776 = cos(qJ(5));
t777 = cos(qJ(4));
t803 = qJD(2) * t777;
t749 = t776 * qJD(4) - t773 * t803;
t774 = sin(qJ(4));
t801 = qJD(2) * qJD(4);
t797 = t774 * t801;
t754 = t777 * qJDD(2) - t797;
t719 = t749 * qJD(5) + t773 * qJDD(4) + t776 * t754;
t750 = t773 * qJD(4) + t776 * t803;
t721 = -t749 * mrSges(7,1) + t750 * mrSges(7,2);
t780 = qJD(2) ^ 2;
t784 = -t780 * qJ(3) + qJDD(3) - t705;
t820 = -pkin(2) - pkin(8);
t701 = t820 * qJDD(2) + t784;
t731 = -t770 * t755 + t772 * t766;
t697 = t774 * t701 + t777 * t731;
t752 = (pkin(4) * t774 - pkin(9) * t777) * qJD(2);
t779 = qJD(4) ^ 2;
t802 = t774 * qJD(2);
t691 = -t779 * pkin(4) + qJDD(4) * pkin(9) - t752 * t802 + t697;
t700 = t820 * t780 - t822;
t796 = t777 * t801;
t753 = -t774 * qJDD(2) - t796;
t694 = (-t754 + t797) * pkin(9) + (-t753 + t796) * pkin(4) + t700;
t686 = -t773 * t691 + t776 * t694;
t746 = qJDD(5) - t753;
t760 = qJD(5) + t802;
t683 = -0.2e1 * qJD(6) * t750 + (t749 * t760 - t719) * qJ(6) + (t749 * t750 + t746) * pkin(5) + t686;
t725 = -t760 * mrSges(7,2) + t749 * mrSges(7,3);
t799 = m(7) * t683 + t746 * mrSges(7,1) + t760 * t725;
t680 = -t719 * mrSges(7,3) - t750 * t721 + t799;
t687 = t776 * t691 + t773 * t694;
t718 = -t750 * qJD(5) + t776 * qJDD(4) - t773 * t754;
t727 = t760 * pkin(5) - t750 * qJ(6);
t745 = t749 ^ 2;
t685 = -t745 * pkin(5) + t718 * qJ(6) + 0.2e1 * qJD(6) * t749 - t760 * t727 + t687;
t805 = t816 * t749 + t825 * t750 + t814 * t760;
t806 = -t824 * t749 - t816 * t750 - t813 * t760;
t821 = mrSges(6,1) * t686 + mrSges(7,1) * t683 - mrSges(6,2) * t687 - mrSges(7,2) * t685 + pkin(5) * t680 + t813 * t718 + t814 * t719 + t823 * t746 - t805 * t749 - t806 * t750;
t819 = mrSges(3,1) - mrSges(4,2);
t818 = -mrSges(6,2) - mrSges(7,2);
t817 = (-Ifges(4,4) + Ifges(3,5));
t815 = Ifges(4,5) - Ifges(3,6);
t751 = (mrSges(5,1) * t774 + mrSges(5,2) * t777) * qJD(2);
t758 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t803;
t722 = -t749 * mrSges(6,1) + t750 * mrSges(6,2);
t726 = -t760 * mrSges(6,2) + t749 * mrSges(6,3);
t674 = m(6) * t686 + t746 * mrSges(6,1) + t760 * t726 + (-t721 - t722) * t750 + (-mrSges(6,3) - mrSges(7,3)) * t719 + t799;
t798 = m(7) * t685 + t718 * mrSges(7,3) + t749 * t721;
t728 = t760 * mrSges(7,1) - t750 * mrSges(7,3);
t804 = -t760 * mrSges(6,1) + t750 * mrSges(6,3) - t728;
t677 = m(6) * t687 + t718 * mrSges(6,3) + t749 * t722 + t818 * t746 + t804 * t760 + t798;
t793 = -t773 * t674 + t776 * t677;
t668 = m(5) * t697 - qJDD(4) * mrSges(5,2) + t753 * mrSges(5,3) - qJD(4) * t758 - t751 * t802 + t793;
t696 = t777 * t701 - t774 * t731;
t757 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t802;
t690 = -qJDD(4) * pkin(4) - t779 * pkin(9) + t752 * t803 - t696;
t688 = -t718 * pkin(5) - t745 * qJ(6) + t750 * t727 + qJDD(6) + t690;
t790 = -m(7) * t688 + t718 * mrSges(7,1) + t749 * t725;
t781 = -m(6) * t690 + t718 * mrSges(6,1) + t818 * t719 + t749 * t726 + t804 * t750 + t790;
t678 = m(5) * t696 + qJDD(4) * mrSges(5,1) - t754 * mrSges(5,3) + qJD(4) * t757 - t751 * t803 + t781;
t660 = t774 * t668 + t777 * t678;
t703 = -qJDD(2) * pkin(2) + t784;
t787 = -m(4) * t703 + (t780 * mrSges(4,3)) - t660;
t655 = m(3) * t705 - (t780 * mrSges(3,2)) + t819 * qJDD(2) + t787;
t811 = t655 * t778;
t794 = t777 * t668 - t774 * t678;
t659 = m(4) * t731 + t794;
t658 = m(3) * t731 + t659;
t702 = t780 * pkin(2) + t822;
t672 = t776 * t674 + t773 * t677;
t786 = -m(5) * t700 + t753 * mrSges(5,1) - t754 * mrSges(5,2) - t757 * t802 - t758 * t803 - t672;
t782 = -m(4) * t702 + (t780 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t786;
t666 = m(3) * t706 - (t780 * mrSges(3,1)) - qJDD(2) * mrSges(3,2) + t782;
t646 = -t770 * t658 + t666 * t809 + t772 * t811;
t644 = m(2) * t755 + t646;
t651 = -t775 * t655 + t778 * t666;
t650 = m(2) * t756 + t651;
t808 = t771 * t644 + t769 * t650;
t807 = -t813 * t749 - t814 * t750 - t823 * t760;
t645 = t772 * t658 + t666 * t810 + t770 * t811;
t795 = -t769 * t644 + t771 * t650;
t791 = m(2) * t766 + t645;
t681 = t719 * mrSges(7,2) + t750 * t728 - t790;
t662 = -mrSges(6,1) * t690 + mrSges(6,3) * t687 - mrSges(7,1) * t688 + mrSges(7,3) * t685 - pkin(5) * t681 + qJ(6) * t798 + (-qJ(6) * t728 + t805) * t760 + t807 * t750 + (-qJ(6) * mrSges(7,2) + t813) * t746 + t816 * t719 + t824 * t718;
t670 = mrSges(6,2) * t690 + mrSges(7,2) * t688 - mrSges(6,3) * t686 - mrSges(7,3) * t683 - qJ(6) * t680 + t816 * t718 + t825 * t719 + t814 * t746 - t807 * t749 + t806 * t760;
t736 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t777 - Ifges(5,6) * t774) * qJD(2);
t737 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t777 - Ifges(5,2) * t774) * qJD(2);
t647 = mrSges(5,2) * t700 - mrSges(5,3) * t696 + Ifges(5,1) * t754 + Ifges(5,4) * t753 + Ifges(5,5) * qJDD(4) - pkin(9) * t672 - qJD(4) * t737 - t773 * t662 + t776 * t670 - t736 * t802;
t738 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t777 - Ifges(5,4) * t774) * qJD(2);
t652 = -mrSges(5,1) * t700 + mrSges(5,3) * t697 + Ifges(5,4) * t754 + Ifges(5,2) * t753 + Ifges(5,6) * qJDD(4) - pkin(4) * t672 + qJD(4) * t738 - t736 * t803 - t821;
t641 = -mrSges(4,1) * t702 + mrSges(3,3) * t706 - pkin(2) * t659 - pkin(3) * t786 - pkin(8) * t794 - t815 * qJDD(2) - t774 * t647 - t777 * t652 - t819 * t731 + (t817 * t780);
t783 = mrSges(5,1) * t696 - mrSges(5,2) * t697 + Ifges(5,5) * t754 + Ifges(5,6) * t753 + Ifges(5,3) * qJDD(4) + pkin(4) * t781 + pkin(9) * t793 + t776 * t662 + t773 * t670 + t737 * t803 + t738 * t802;
t642 = mrSges(4,1) * t703 - mrSges(3,3) * t705 + t783 + pkin(3) * t660 + (mrSges(3,2) - mrSges(4,3)) * t731 + t815 * t780 + t817 * qJDD(2) - qJ(3) * t659;
t788 = pkin(7) * t651 + t641 * t778 + t642 * t775;
t656 = qJDD(2) * mrSges(4,2) - t787;
t640 = mrSges(3,1) * t705 - mrSges(3,2) * t706 + mrSges(4,2) * t703 - mrSges(4,3) * t702 + t777 * t647 - t774 * t652 - pkin(8) * t660 - pkin(2) * t656 + qJ(3) * t782 + (Ifges(3,3) + Ifges(4,1)) * qJDD(2);
t639 = mrSges(2,2) * t766 - mrSges(2,3) * t755 - t775 * t641 + t778 * t642 + (-t645 * t770 - t646 * t772) * pkin(7);
t638 = -mrSges(2,1) * t766 + mrSges(2,3) * t756 - pkin(1) * t645 - t770 * t640 + t788 * t772;
t1 = [-m(1) * g(1) + t795; -m(1) * g(2) + t808; -m(1) * g(3) + t791; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t808 - t769 * t638 + t771 * t639; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t795 + t771 * t638 + t769 * t639; -mrSges(1,1) * g(2) + mrSges(2,1) * t755 + mrSges(1,2) * g(1) - mrSges(2,2) * t756 + pkin(1) * t646 + t772 * t640 + t788 * t770; t791; t640; t656; t783; t821; t681;];
tauJB  = t1;
