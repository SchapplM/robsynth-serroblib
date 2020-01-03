% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRR12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR12_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR12_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR12_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:24
% EndTime: 2019-12-31 20:29:30
% DurationCPUTime: 3.65s
% Computational Cost: add. (34837->291), mult. (71635->356), div. (0->0), fcn. (42139->8), ass. (0->119)
t828 = Ifges(3,1) + Ifges(4,1);
t820 = Ifges(3,4) - Ifges(4,5);
t819 = Ifges(3,5) + Ifges(4,4);
t827 = Ifges(3,2) + Ifges(4,3);
t818 = Ifges(3,6) - Ifges(4,6);
t826 = Ifges(3,3) + Ifges(4,2);
t784 = sin(qJ(2));
t788 = cos(qJ(2));
t750 = (-mrSges(4,1) * t788 - mrSges(4,3) * t784) * qJD(1);
t809 = qJD(1) * qJD(2);
t807 = t788 * t809;
t752 = t784 * qJDD(1) + t807;
t785 = sin(qJ(1));
t789 = cos(qJ(1));
t763 = -t789 * g(1) - t785 * g(2);
t791 = qJD(1) ^ 2;
t740 = -t791 * pkin(1) + qJDD(1) * pkin(6) + t763;
t722 = -t784 * g(3) + t788 * t740;
t749 = (-pkin(2) * t788 - qJ(3) * t784) * qJD(1);
t790 = qJD(2) ^ 2;
t810 = qJD(1) * t788;
t822 = 2 * qJD(3);
t701 = -t790 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t822 + t749 * t810 + t722;
t808 = t784 * t809;
t753 = t788 * qJDD(1) - t808;
t811 = qJD(1) * t784;
t761 = -qJD(2) * pkin(3) - pkin(7) * t811;
t816 = t788 ^ 2 * t791;
t691 = -pkin(3) * t816 - t753 * pkin(7) + qJD(2) * t761 + t701;
t721 = -t788 * g(3) - t784 * t740;
t707 = -qJDD(2) * pkin(2) - t790 * qJ(3) + t749 * t811 + qJDD(3) - t721;
t692 = (-t752 + t807) * pkin(7) + (-t784 * t788 * t791 - qJDD(2)) * pkin(3) + t707;
t783 = sin(qJ(4));
t787 = cos(qJ(4));
t683 = t787 * t691 + t783 * t692;
t738 = (-t783 * t788 + t784 * t787) * qJD(1);
t708 = -t738 * qJD(4) - t783 * t752 - t787 * t753;
t737 = (t783 * t784 + t787 * t788) * qJD(1);
t717 = t737 * mrSges(5,1) + t738 * mrSges(5,2);
t774 = -qJD(2) + qJD(4);
t724 = t774 * mrSges(5,1) - t738 * mrSges(5,3);
t773 = -qJDD(2) + qJDD(4);
t762 = t785 * g(1) - t789 * g(2);
t739 = -qJDD(1) * pkin(1) - t791 * pkin(6) - t762;
t800 = -t753 * pkin(2) + t739 + (-t752 - t807) * qJ(3);
t685 = -pkin(2) * t808 + t753 * pkin(3) - pkin(7) * t816 - t800 + (t761 + t822) * t811;
t709 = -t737 * qJD(4) + t787 * t752 - t783 * t753;
t678 = (t738 * t774 - t708) * pkin(4) + (t737 * t774 - t709) * pkin(8) + t685;
t718 = t737 * pkin(4) - t738 * pkin(8);
t772 = t774 ^ 2;
t680 = -t772 * pkin(4) + t773 * pkin(8) - t737 * t718 + t683;
t782 = sin(qJ(5));
t786 = cos(qJ(5));
t676 = t786 * t678 - t782 * t680;
t719 = -t782 * t738 + t786 * t774;
t688 = t719 * qJD(5) + t786 * t709 + t782 * t773;
t720 = t786 * t738 + t782 * t774;
t698 = -t719 * mrSges(6,1) + t720 * mrSges(6,2);
t706 = qJDD(5) - t708;
t730 = qJD(5) + t737;
t710 = -t730 * mrSges(6,2) + t719 * mrSges(6,3);
t672 = m(6) * t676 + t706 * mrSges(6,1) - t688 * mrSges(6,3) - t720 * t698 + t730 * t710;
t677 = t782 * t678 + t786 * t680;
t687 = -t720 * qJD(5) - t782 * t709 + t786 * t773;
t711 = t730 * mrSges(6,1) - t720 * mrSges(6,3);
t673 = m(6) * t677 - t706 * mrSges(6,2) + t687 * mrSges(6,3) + t719 * t698 - t730 * t711;
t803 = -t782 * t672 + t786 * t673;
t660 = m(5) * t683 - t773 * mrSges(5,2) + t708 * mrSges(5,3) - t737 * t717 - t774 * t724 + t803;
t682 = -t783 * t691 + t787 * t692;
t723 = -t774 * mrSges(5,2) - t737 * mrSges(5,3);
t679 = -t773 * pkin(4) - t772 * pkin(8) + t738 * t718 - t682;
t797 = -m(6) * t679 + t687 * mrSges(6,1) - t688 * mrSges(6,2) + t719 * t710 - t720 * t711;
t668 = m(5) * t682 + t773 * mrSges(5,1) - t709 * mrSges(5,3) - t738 * t717 + t774 * t723 + t797;
t654 = t783 * t660 + t787 * t668;
t760 = mrSges(4,2) * t810 + qJD(2) * mrSges(4,3);
t796 = -m(4) * t707 + qJDD(2) * mrSges(4,1) + qJD(2) * t760 - t654;
t653 = t752 * mrSges(4,2) + t750 * t811 - t796;
t693 = Ifges(6,5) * t720 + Ifges(6,6) * t719 + Ifges(6,3) * t730;
t695 = Ifges(6,1) * t720 + Ifges(6,4) * t719 + Ifges(6,5) * t730;
t666 = -mrSges(6,1) * t679 + mrSges(6,3) * t677 + Ifges(6,4) * t688 + Ifges(6,2) * t687 + Ifges(6,6) * t706 - t720 * t693 + t730 * t695;
t694 = Ifges(6,4) * t720 + Ifges(6,2) * t719 + Ifges(6,6) * t730;
t667 = mrSges(6,2) * t679 - mrSges(6,3) * t676 + Ifges(6,1) * t688 + Ifges(6,4) * t687 + Ifges(6,5) * t706 + t719 * t693 - t730 * t694;
t713 = Ifges(5,4) * t738 - Ifges(5,2) * t737 + Ifges(5,6) * t774;
t714 = Ifges(5,1) * t738 - Ifges(5,4) * t737 + Ifges(5,5) * t774;
t795 = -mrSges(5,1) * t682 + mrSges(5,2) * t683 - Ifges(5,5) * t709 - Ifges(5,6) * t708 - Ifges(5,3) * t773 - pkin(4) * t797 - pkin(8) * t803 - t786 * t666 - t782 * t667 - t738 * t713 - t737 * t714;
t758 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t811;
t804 = t787 * t660 - t783 * t668;
t799 = m(4) * t701 + qJDD(2) * mrSges(4,3) + qJD(2) * t758 + t750 * t810 + t804;
t812 = t819 * qJD(2) + (t828 * t784 + t820 * t788) * qJD(1);
t813 = -t818 * qJD(2) + (-t820 * t784 - t827 * t788) * qJD(1);
t825 = -(t813 * t784 + t812 * t788) * qJD(1) + t826 * qJDD(2) + t819 * t752 + t818 * t753 + mrSges(3,1) * t721 - mrSges(4,1) * t707 - mrSges(3,2) * t722 + mrSges(4,3) * t701 - pkin(2) * t653 - pkin(3) * t654 + qJ(3) * (t753 * mrSges(4,2) + t799) + t795;
t821 = mrSges(3,3) + mrSges(4,2);
t751 = (-mrSges(3,1) * t788 + mrSges(3,2) * t784) * qJD(1);
t757 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t811;
t650 = m(3) * t722 - qJDD(2) * mrSges(3,2) - qJD(2) * t757 + t751 * t810 + t821 * t753 + t799;
t759 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t810;
t651 = m(3) * t721 + qJDD(2) * mrSges(3,1) + qJD(2) * t759 - t821 * t752 + (-t750 - t751) * t811 + t796;
t805 = t788 * t650 - t784 * t651;
t642 = m(2) * t763 - t791 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t805;
t697 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t811 + t800;
t662 = t786 * t672 + t782 * t673;
t801 = -m(5) * t685 + t708 * mrSges(5,1) - t709 * mrSges(5,2) - t737 * t723 - t738 * t724 - t662;
t658 = m(4) * t697 - t753 * mrSges(4,1) - t752 * mrSges(4,3) - t758 * t811 - t760 * t810 + t801;
t793 = -m(3) * t739 + t753 * mrSges(3,1) - t752 * mrSges(3,2) - t757 * t811 + t759 * t810 - t658;
t656 = m(2) * t762 + qJDD(1) * mrSges(2,1) - t791 * mrSges(2,2) + t793;
t815 = t785 * t642 + t789 * t656;
t644 = t784 * t650 + t788 * t651;
t814 = t826 * qJD(2) + (t819 * t784 + t818 * t788) * qJD(1);
t806 = t789 * t642 - t785 * t656;
t712 = Ifges(5,5) * t738 - Ifges(5,6) * t737 + Ifges(5,3) * t774;
t645 = mrSges(5,2) * t685 - mrSges(5,3) * t682 + Ifges(5,1) * t709 + Ifges(5,4) * t708 + Ifges(5,5) * t773 - pkin(8) * t662 - t782 * t666 + t786 * t667 - t737 * t712 - t774 * t713;
t794 = mrSges(6,1) * t676 - mrSges(6,2) * t677 + Ifges(6,5) * t688 + Ifges(6,6) * t687 + Ifges(6,3) * t706 + t720 * t694 - t719 * t695;
t646 = -mrSges(5,1) * t685 + mrSges(5,3) * t683 + Ifges(5,4) * t709 + Ifges(5,2) * t708 + Ifges(5,6) * t773 - pkin(4) * t662 - t738 * t712 + t774 * t714 - t794;
t637 = -mrSges(3,1) * t739 - mrSges(4,1) * t697 + mrSges(4,2) * t701 + mrSges(3,3) * t722 - pkin(2) * t658 - pkin(3) * t801 - pkin(7) * t804 + t812 * qJD(2) + t818 * qJDD(2) - t783 * t645 - t787 * t646 + t820 * t752 + t827 * t753 - t814 * t811;
t639 = mrSges(3,2) * t739 + mrSges(4,2) * t707 - mrSges(3,3) * t721 - mrSges(4,3) * t697 - pkin(7) * t654 - qJ(3) * t658 + t813 * qJD(2) + t819 * qJDD(2) + t787 * t645 - t783 * t646 + t828 * t752 + t820 * t753 + t814 * t810;
t798 = mrSges(2,1) * t762 - mrSges(2,2) * t763 + Ifges(2,3) * qJDD(1) + pkin(1) * t793 + pkin(6) * t805 + t788 * t637 + t784 * t639;
t635 = mrSges(2,1) * g(3) + mrSges(2,3) * t763 + t791 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t644 - t825;
t634 = -mrSges(2,2) * g(3) - mrSges(2,3) * t762 + Ifges(2,5) * qJDD(1) - t791 * Ifges(2,6) - pkin(6) * t644 - t784 * t637 + t788 * t639;
t1 = [-m(1) * g(1) + t806; -m(1) * g(2) + t815; (-m(1) - m(2)) * g(3) + t644; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t815 + t789 * t634 - t785 * t635; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t806 + t785 * t634 + t789 * t635; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t798; t798; t825; t653; -t795; t794;];
tauJB = t1;
