% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPPR10
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPPR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR10_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR10_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR10_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:24
% EndTime: 2019-12-31 19:43:30
% DurationCPUTime: 3.89s
% Computational Cost: add. (34956->291), mult. (76395->351), div. (0->0), fcn. (47044->8), ass. (0->118)
t799 = -2 * qJD(3);
t798 = Ifges(4,1) + Ifges(5,1);
t789 = Ifges(4,4) - Ifges(5,5);
t797 = Ifges(4,5) + Ifges(5,4);
t796 = Ifges(4,2) + Ifges(5,3);
t795 = Ifges(5,2) + Ifges(4,3);
t794 = Ifges(4,6) - Ifges(5,6);
t753 = sin(qJ(1));
t756 = cos(qJ(1));
t739 = t753 * g(1) - t756 * g(2);
t758 = qJD(1) ^ 2;
t720 = -qJDD(1) * pkin(1) - t758 * pkin(6) - t739;
t752 = sin(qJ(2));
t755 = cos(qJ(2));
t776 = qJD(1) * qJD(2);
t774 = t755 * t776;
t734 = t752 * qJDD(1) + t774;
t743 = t752 * t776;
t735 = t755 * qJDD(1) - t743;
t680 = (-t734 - t774) * qJ(3) + (-t735 + t743) * pkin(2) + t720;
t740 = -t756 * g(1) - t753 * g(2);
t721 = -t758 * pkin(1) + qJDD(1) * pkin(6) + t740;
t703 = -t752 * g(3) + t755 * t721;
t732 = (-pkin(2) * t755 - qJ(3) * t752) * qJD(1);
t757 = qJD(2) ^ 2;
t777 = t755 * qJD(1);
t684 = -t757 * pkin(2) + qJDD(2) * qJ(3) + t732 * t777 + t703;
t750 = sin(pkin(8));
t779 = qJD(1) * t752;
t786 = cos(pkin(8));
t727 = t750 * qJD(2) + t786 * t779;
t666 = t786 * t680 - t750 * t684 + t727 * t799;
t710 = t750 * qJDD(2) + t786 * t734;
t702 = -t755 * g(3) - t752 * t721;
t763 = qJDD(2) * pkin(2) + t757 * qJ(3) - t732 * t779 - qJDD(3) + t702;
t726 = -t786 * qJD(2) + t750 * t779;
t775 = t726 * t777;
t793 = -(t710 + t775) * qJ(4) - t763;
t699 = t726 * pkin(3) - t727 * qJ(4);
t785 = t755 ^ 2 * t758;
t664 = t735 * pkin(3) - qJ(4) * t785 + t727 * t699 + qJDD(4) - t666;
t658 = (-t710 + t775) * pkin(7) + (t726 * t727 + t735) * pkin(4) + t664;
t667 = t750 * t680 + t786 * t684 + t726 * t799;
t791 = -2 * qJD(4);
t663 = -pkin(3) * t785 - t735 * qJ(4) - t726 * t699 + t777 * t791 + t667;
t709 = -t786 * qJDD(2) + t750 * t734;
t711 = pkin(4) * t777 - t727 * pkin(7);
t724 = t726 ^ 2;
t659 = -t724 * pkin(4) + t709 * pkin(7) - t711 * t777 + t663;
t751 = sin(qJ(5));
t754 = cos(qJ(5));
t657 = t751 * t658 + t754 * t659;
t661 = -t724 * pkin(7) + (-pkin(3) - pkin(4)) * t709 + (pkin(3) * t777 + (2 * qJD(4)) + t711) * t727 - t793;
t698 = t751 * t726 + t754 * t727;
t672 = -t698 * qJD(5) + t754 * t709 - t751 * t710;
t697 = t754 * t726 - t751 * t727;
t673 = t697 * qJD(5) + t751 * t709 + t754 * t710;
t741 = qJD(5) + t777;
t674 = Ifges(6,5) * t698 + Ifges(6,6) * t697 + Ifges(6,3) * t741;
t676 = Ifges(6,1) * t698 + Ifges(6,4) * t697 + Ifges(6,5) * t741;
t731 = qJDD(5) + t735;
t646 = -mrSges(6,1) * t661 + mrSges(6,3) * t657 + Ifges(6,4) * t673 + Ifges(6,2) * t672 + Ifges(6,6) * t731 - t698 * t674 + t741 * t676;
t656 = t754 * t658 - t751 * t659;
t675 = Ifges(6,4) * t698 + Ifges(6,2) * t697 + Ifges(6,6) * t741;
t647 = mrSges(6,2) * t661 - mrSges(6,3) * t656 + Ifges(6,1) * t673 + Ifges(6,4) * t672 + Ifges(6,5) * t731 + t697 * t674 - t741 * t675;
t665 = t727 * t791 + (-t727 * t777 + t709) * pkin(3) + t793;
t706 = -t726 * mrSges(5,2) - mrSges(5,3) * t777;
t708 = mrSges(5,1) * t777 + t727 * mrSges(5,2);
t685 = -t741 * mrSges(6,2) + t697 * mrSges(6,3);
t686 = t741 * mrSges(6,1) - t698 * mrSges(6,3);
t765 = -m(6) * t661 + t672 * mrSges(6,1) - t673 * mrSges(6,2) + t697 * t685 - t698 * t686;
t651 = m(5) * t665 + t709 * mrSges(5,1) - t710 * mrSges(5,3) + t726 * t706 - t727 * t708 + t765;
t678 = -t697 * mrSges(6,1) + t698 * mrSges(6,2);
t653 = m(6) * t656 + t731 * mrSges(6,1) - t673 * mrSges(6,3) - t698 * t678 + t741 * t685;
t654 = m(6) * t657 - t731 * mrSges(6,2) + t672 * mrSges(6,3) + t697 * t678 - t741 * t686;
t771 = -t751 * t653 + t754 * t654;
t781 = t789 * t726 - t798 * t727 + t797 * t777;
t782 = t794 * t726 - t797 * t727 + t795 * t777;
t626 = mrSges(4,1) * t763 - mrSges(5,1) * t665 + mrSges(5,2) * t663 + mrSges(4,3) * t667 - pkin(3) * t651 - pkin(4) * t765 - pkin(7) * t771 - t754 * t646 - t751 * t647 - t796 * t709 + t789 * t710 + t782 * t727 - t735 * t794 + t781 * t777;
t645 = t754 * t653 + t751 * t654;
t783 = t796 * t726 - t789 * t727 + t794 * t777;
t627 = -mrSges(4,2) * t763 + mrSges(5,2) * t664 - mrSges(4,3) * t666 - mrSges(5,3) * t665 - pkin(7) * t645 - qJ(4) * t651 - t751 * t646 + t754 * t647 - t789 * t709 + t798 * t710 + t782 * t726 - t735 * t797 - t783 * t777;
t707 = -mrSges(4,1) * t777 - t727 * mrSges(4,3);
t767 = m(5) * t663 - t735 * mrSges(5,3) + t771;
t700 = t726 * mrSges(5,1) - t727 * mrSges(5,3);
t780 = -t726 * mrSges(4,1) - t727 * mrSges(4,2) - t700;
t790 = -mrSges(4,3) - mrSges(5,2);
t642 = m(4) * t667 + t735 * mrSges(4,2) + t780 * t726 + t790 * t709 + (t707 - t708) * t777 + t767;
t761 = -m(5) * t664 - t735 * mrSges(5,1) - t706 * t777 - t645;
t766 = mrSges(4,2) * t777 - t726 * mrSges(4,3);
t643 = m(4) * t666 - t735 * mrSges(4,1) + t790 * t710 + t780 * t727 - t766 * t777 + t761;
t640 = t786 * t642 - t750 * t643;
t650 = -m(4) * t763 + t709 * mrSges(4,1) + t710 * mrSges(4,2) + t727 * t707 + t726 * t766 + t651;
t718 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t752 + Ifges(3,2) * t755) * qJD(1);
t719 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t752 + Ifges(3,4) * t755) * qJD(1);
t792 = mrSges(3,1) * t702 - mrSges(3,2) * t703 + Ifges(3,5) * t734 + Ifges(3,6) * t735 + Ifges(3,3) * qJDD(2) - pkin(2) * t650 + qJ(3) * t640 + (t718 * t752 - t719 * t755) * qJD(1) + t786 * t626 + t750 * t627;
t733 = (-mrSges(3,1) * t755 + mrSges(3,2) * t752) * qJD(1);
t737 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t779;
t638 = m(3) * t703 - qJDD(2) * mrSges(3,2) + t735 * mrSges(3,3) - qJD(2) * t737 + t733 * t777 + t640;
t738 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t777;
t649 = m(3) * t702 + qJDD(2) * mrSges(3,1) - t734 * mrSges(3,3) + qJD(2) * t738 - t733 * t779 - t650;
t772 = t755 * t638 - t752 * t649;
t630 = m(2) * t740 - t758 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t772;
t639 = t750 * t642 + t786 * t643;
t760 = -m(3) * t720 + t735 * mrSges(3,1) - t734 * mrSges(3,2) - t737 * t779 + t738 * t777 - t639;
t634 = m(2) * t739 + qJDD(1) * mrSges(2,1) - t758 * mrSges(2,2) + t760;
t784 = t753 * t630 + t756 * t634;
t632 = t752 * t638 + t755 * t649;
t773 = t756 * t630 - t753 * t634;
t717 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t752 + Ifges(3,6) * t755) * qJD(1);
t623 = mrSges(3,2) * t720 - mrSges(3,3) * t702 + Ifges(3,1) * t734 + Ifges(3,4) * t735 + Ifges(3,5) * qJDD(2) - qJ(3) * t639 - qJD(2) * t718 - t750 * t626 + t786 * t627 + t717 * t777;
t644 = t710 * mrSges(5,2) + t727 * t700 - t761;
t762 = mrSges(6,1) * t656 - mrSges(6,2) * t657 + Ifges(6,5) * t673 + Ifges(6,6) * t672 + Ifges(6,3) * t731 + t698 * t675 - t697 * t676;
t625 = t762 + (Ifges(3,2) + t795) * t735 + t783 * t727 + Ifges(3,4) * t734 + qJD(2) * t719 - mrSges(3,1) * t720 + mrSges(3,3) * t703 - mrSges(5,3) * t663 + mrSges(5,1) * t664 - mrSges(4,1) * t666 + mrSges(4,2) * t667 + pkin(4) * t645 + pkin(3) * t644 + (qJ(4) * t700 + t781) * t726 - pkin(2) * t639 + Ifges(3,6) * qJDD(2) + (qJ(4) * mrSges(5,2) + t794) * t709 - qJ(4) * (-t708 * t777 + t767) - t797 * t710 - t717 * t779;
t764 = mrSges(2,1) * t739 - mrSges(2,2) * t740 + Ifges(2,3) * qJDD(1) + pkin(1) * t760 + pkin(6) * t772 + t752 * t623 + t755 * t625;
t621 = mrSges(2,1) * g(3) + mrSges(2,3) * t740 + t758 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t632 - t792;
t620 = -mrSges(2,2) * g(3) - mrSges(2,3) * t739 + Ifges(2,5) * qJDD(1) - t758 * Ifges(2,6) - pkin(6) * t632 + t755 * t623 - t752 * t625;
t1 = [-m(1) * g(1) + t773; -m(1) * g(2) + t784; (-m(1) - m(2)) * g(3) + t632; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t784 + t756 * t620 - t753 * t621; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t773 + t753 * t620 + t756 * t621; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t764; t764; t792; t650; t644; t762;];
tauJB = t1;
