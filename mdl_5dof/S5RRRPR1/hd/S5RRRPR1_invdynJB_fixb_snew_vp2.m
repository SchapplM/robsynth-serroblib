% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRPR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:32
% EndTime: 2019-12-05 18:38:44
% DurationCPUTime: 11.68s
% Computational Cost: add. (168226->315), mult. (386471->400), div. (0->0), fcn. (275400->10), ass. (0->125)
t775 = sin(qJ(2));
t779 = cos(qJ(2));
t796 = qJD(1) * qJD(2);
t753 = t775 * qJDD(1) + t779 * t796;
t776 = sin(qJ(1));
t780 = cos(qJ(1));
t760 = -t780 * g(1) - t776 * g(2);
t781 = qJD(1) ^ 2;
t748 = -t781 * pkin(1) + qJDD(1) * pkin(6) + t760;
t800 = t775 * t748;
t801 = pkin(2) * t781;
t713 = qJDD(2) * pkin(2) - t753 * pkin(7) - t800 + (pkin(7) * t796 + t775 * t801 - g(3)) * t779;
t736 = -t775 * g(3) + t779 * t748;
t754 = t779 * qJDD(1) - t775 * t796;
t798 = qJD(1) * t775;
t758 = qJD(2) * pkin(2) - pkin(7) * t798;
t770 = t779 ^ 2;
t714 = t754 * pkin(7) - qJD(2) * t758 - t770 * t801 + t736;
t774 = sin(qJ(3));
t778 = cos(qJ(3));
t691 = t778 * t713 - t774 * t714;
t745 = (-t774 * t775 + t778 * t779) * qJD(1);
t719 = t745 * qJD(3) + t778 * t753 + t774 * t754;
t746 = (t774 * t779 + t775 * t778) * qJD(1);
t767 = qJDD(2) + qJDD(3);
t768 = qJD(2) + qJD(3);
t678 = (t745 * t768 - t719) * qJ(4) + (t745 * t746 + t767) * pkin(3) + t691;
t692 = t774 * t713 + t778 * t714;
t718 = -t746 * qJD(3) - t774 * t753 + t778 * t754;
t738 = t768 * pkin(3) - t746 * qJ(4);
t741 = t745 ^ 2;
t680 = -t741 * pkin(3) + t718 * qJ(4) - t768 * t738 + t692;
t771 = sin(pkin(9));
t772 = cos(pkin(9));
t733 = t771 * t745 + t772 * t746;
t665 = -0.2e1 * qJD(4) * t733 + t772 * t678 - t771 * t680;
t698 = t771 * t718 + t772 * t719;
t732 = t772 * t745 - t771 * t746;
t662 = (t732 * t768 - t698) * pkin(8) + (t732 * t733 + t767) * pkin(4) + t665;
t666 = 0.2e1 * qJD(4) * t732 + t771 * t678 + t772 * t680;
t697 = t772 * t718 - t771 * t719;
t723 = t768 * pkin(4) - t733 * pkin(8);
t729 = t732 ^ 2;
t663 = -t729 * pkin(4) + t697 * pkin(8) - t768 * t723 + t666;
t773 = sin(qJ(5));
t777 = cos(qJ(5));
t660 = t777 * t662 - t773 * t663;
t706 = t777 * t732 - t773 * t733;
t674 = t706 * qJD(5) + t773 * t697 + t777 * t698;
t707 = t773 * t732 + t777 * t733;
t686 = -t706 * mrSges(6,1) + t707 * mrSges(6,2);
t765 = qJD(5) + t768;
t699 = -t765 * mrSges(6,2) + t706 * mrSges(6,3);
t764 = qJDD(5) + t767;
t656 = m(6) * t660 + t764 * mrSges(6,1) - t674 * mrSges(6,3) - t707 * t686 + t765 * t699;
t661 = t773 * t662 + t777 * t663;
t673 = -t707 * qJD(5) + t777 * t697 - t773 * t698;
t700 = t765 * mrSges(6,1) - t707 * mrSges(6,3);
t657 = m(6) * t661 - t764 * mrSges(6,2) + t673 * mrSges(6,3) + t706 * t686 - t765 * t700;
t647 = t777 * t656 + t773 * t657;
t708 = -t732 * mrSges(5,1) + t733 * mrSges(5,2);
t721 = -t768 * mrSges(5,2) + t732 * mrSges(5,3);
t644 = m(5) * t665 + t767 * mrSges(5,1) - t698 * mrSges(5,3) - t733 * t708 + t768 * t721 + t647;
t722 = t768 * mrSges(5,1) - t733 * mrSges(5,3);
t791 = -t773 * t656 + t777 * t657;
t645 = m(5) * t666 - t767 * mrSges(5,2) + t697 * mrSges(5,3) + t732 * t708 - t768 * t722 + t791;
t640 = t772 * t644 + t771 * t645;
t734 = -t745 * mrSges(4,1) + t746 * mrSges(4,2);
t737 = -t768 * mrSges(4,2) + t745 * mrSges(4,3);
t637 = m(4) * t691 + t767 * mrSges(4,1) - t719 * mrSges(4,3) - t746 * t734 + t768 * t737 + t640;
t739 = t768 * mrSges(4,1) - t746 * mrSges(4,3);
t792 = -t771 * t644 + t772 * t645;
t638 = m(4) * t692 - t767 * mrSges(4,2) + t718 * mrSges(4,3) + t745 * t734 - t768 * t739 + t792;
t631 = t778 * t637 + t774 * t638;
t735 = -t779 * g(3) - t800;
t743 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t775 + Ifges(3,2) * t779) * qJD(1);
t744 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t775 + Ifges(3,4) * t779) * qJD(1);
t702 = Ifges(5,4) * t733 + Ifges(5,2) * t732 + Ifges(5,6) * t768;
t703 = Ifges(5,1) * t733 + Ifges(5,4) * t732 + Ifges(5,5) * t768;
t727 = Ifges(4,4) * t746 + Ifges(4,2) * t745 + Ifges(4,6) * t768;
t728 = Ifges(4,1) * t746 + Ifges(4,4) * t745 + Ifges(4,5) * t768;
t682 = Ifges(6,4) * t707 + Ifges(6,2) * t706 + Ifges(6,6) * t765;
t683 = Ifges(6,1) * t707 + Ifges(6,4) * t706 + Ifges(6,5) * t765;
t786 = -mrSges(6,1) * t660 + mrSges(6,2) * t661 - Ifges(6,5) * t674 - Ifges(6,6) * t673 - Ifges(6,3) * t764 - t707 * t682 + t706 * t683;
t783 = -mrSges(4,1) * t691 - mrSges(5,1) * t665 + mrSges(4,2) * t692 + mrSges(5,2) * t666 - Ifges(4,5) * t719 - Ifges(5,5) * t698 - Ifges(4,6) * t718 - Ifges(5,6) * t697 - pkin(3) * t640 - pkin(4) * t647 - t733 * t702 + t732 * t703 - t746 * t727 + t745 * t728 + t786 + (-Ifges(4,3) - Ifges(5,3)) * t767;
t802 = mrSges(3,1) * t735 - mrSges(3,2) * t736 + Ifges(3,5) * t753 + Ifges(3,6) * t754 + Ifges(3,3) * qJDD(2) + pkin(2) * t631 + (t775 * t743 - t779 * t744) * qJD(1) - t783;
t752 = (-mrSges(3,1) * t779 + mrSges(3,2) * t775) * qJD(1);
t797 = qJD(1) * t779;
t757 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t797;
t629 = m(3) * t735 + qJDD(2) * mrSges(3,1) - t753 * mrSges(3,3) + qJD(2) * t757 - t752 * t798 + t631;
t756 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t798;
t793 = -t774 * t637 + t778 * t638;
t630 = m(3) * t736 - qJDD(2) * mrSges(3,2) + t754 * mrSges(3,3) - qJD(2) * t756 + t752 * t797 + t793;
t794 = -t775 * t629 + t779 * t630;
t621 = m(2) * t760 - t781 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t794;
t759 = t776 * g(1) - t780 * g(2);
t788 = -qJDD(1) * pkin(1) - t759;
t747 = -t781 * pkin(6) + t788;
t720 = -t754 * pkin(2) + t758 * t798 + (-pkin(7) * t770 - pkin(6)) * t781 + t788;
t688 = -t718 * pkin(3) - t741 * qJ(4) + t746 * t738 + qJDD(4) + t720;
t668 = -t697 * pkin(4) - t729 * pkin(8) + t733 * t723 + t688;
t790 = m(6) * t668 - t673 * mrSges(6,1) + t674 * mrSges(6,2) - t706 * t699 + t707 * t700;
t658 = m(5) * t688 - t697 * mrSges(5,1) + t698 * mrSges(5,2) - t732 * t721 + t733 * t722 + t790;
t785 = m(4) * t720 - t718 * mrSges(4,1) + t719 * mrSges(4,2) - t745 * t737 + t746 * t739 + t658;
t784 = -m(3) * t747 + t754 * mrSges(3,1) - t753 * mrSges(3,2) - t756 * t798 + t757 * t797 - t785;
t651 = m(2) * t759 + qJDD(1) * mrSges(2,1) - t781 * mrSges(2,2) + t784;
t799 = t776 * t621 + t780 * t651;
t623 = t779 * t629 + t775 * t630;
t795 = t780 * t621 - t776 * t651;
t681 = Ifges(6,5) * t707 + Ifges(6,6) * t706 + Ifges(6,3) * t765;
t648 = -mrSges(6,1) * t668 + mrSges(6,3) * t661 + Ifges(6,4) * t674 + Ifges(6,2) * t673 + Ifges(6,6) * t764 - t707 * t681 + t765 * t683;
t649 = mrSges(6,2) * t668 - mrSges(6,3) * t660 + Ifges(6,1) * t674 + Ifges(6,4) * t673 + Ifges(6,5) * t764 + t706 * t681 - t765 * t682;
t701 = Ifges(5,5) * t733 + Ifges(5,6) * t732 + Ifges(5,3) * t768;
t632 = -mrSges(5,1) * t688 + mrSges(5,3) * t666 + Ifges(5,4) * t698 + Ifges(5,2) * t697 + Ifges(5,6) * t767 - pkin(4) * t790 + pkin(8) * t791 + t777 * t648 + t773 * t649 - t733 * t701 + t768 * t703;
t633 = mrSges(5,2) * t688 - mrSges(5,3) * t665 + Ifges(5,1) * t698 + Ifges(5,4) * t697 + Ifges(5,5) * t767 - pkin(8) * t647 - t773 * t648 + t777 * t649 + t732 * t701 - t768 * t702;
t726 = Ifges(4,5) * t746 + Ifges(4,6) * t745 + Ifges(4,3) * t768;
t624 = -mrSges(4,1) * t720 + mrSges(4,3) * t692 + Ifges(4,4) * t719 + Ifges(4,2) * t718 + Ifges(4,6) * t767 - pkin(3) * t658 + qJ(4) * t792 + t772 * t632 + t771 * t633 - t746 * t726 + t768 * t728;
t625 = mrSges(4,2) * t720 - mrSges(4,3) * t691 + Ifges(4,1) * t719 + Ifges(4,4) * t718 + Ifges(4,5) * t767 - qJ(4) * t640 - t771 * t632 + t772 * t633 + t745 * t726 - t768 * t727;
t742 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t775 + Ifges(3,6) * t779) * qJD(1);
t615 = -mrSges(3,1) * t747 + mrSges(3,3) * t736 + Ifges(3,4) * t753 + Ifges(3,2) * t754 + Ifges(3,6) * qJDD(2) - pkin(2) * t785 + pkin(7) * t793 + qJD(2) * t744 + t778 * t624 + t774 * t625 - t742 * t798;
t617 = mrSges(3,2) * t747 - mrSges(3,3) * t735 + Ifges(3,1) * t753 + Ifges(3,4) * t754 + Ifges(3,5) * qJDD(2) - pkin(7) * t631 - qJD(2) * t743 - t774 * t624 + t778 * t625 + t742 * t797;
t787 = mrSges(2,1) * t759 - mrSges(2,2) * t760 + Ifges(2,3) * qJDD(1) + pkin(1) * t784 + pkin(6) * t794 + t779 * t615 + t775 * t617;
t618 = mrSges(2,1) * g(3) + mrSges(2,3) * t760 + t781 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t623 - t802;
t613 = -mrSges(2,2) * g(3) - mrSges(2,3) * t759 + Ifges(2,5) * qJDD(1) - t781 * Ifges(2,6) - pkin(6) * t623 - t775 * t615 + t779 * t617;
t1 = [-m(1) * g(1) + t795; -m(1) * g(2) + t799; (-m(1) - m(2)) * g(3) + t623; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t799 + t780 * t613 - t776 * t618; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t795 + t776 * t613 + t780 * t618; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t787; t787; t802; -t783; t658; -t786;];
tauJB = t1;
