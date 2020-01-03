% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRP9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRP9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP9_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP9_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP9_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:58
% EndTime: 2019-12-31 18:49:03
% DurationCPUTime: 5.28s
% Computational Cost: add. (51513->268), mult. (124838->324), div. (0->0), fcn. (89775->8), ass. (0->117)
t800 = Ifges(5,1) + Ifges(6,1);
t790 = Ifges(5,4) - Ifges(6,5);
t798 = Ifges(6,4) + Ifges(5,5);
t799 = Ifges(5,2) + Ifges(6,3);
t796 = Ifges(5,6) - Ifges(6,6);
t797 = -Ifges(6,2) - Ifges(5,3);
t753 = sin(pkin(8));
t754 = cos(pkin(8));
t756 = sin(qJ(3));
t758 = cos(qJ(3));
t768 = -t753 * t756 + t754 * t758;
t729 = t768 * qJD(1);
t769 = t753 * t758 + t754 * t756;
t730 = t769 * qJD(1);
t755 = sin(qJ(4));
t793 = cos(qJ(4));
t710 = -t793 * t729 + t755 * t730;
t711 = t755 * t729 + t793 * t730;
t751 = qJD(3) + qJD(4);
t795 = t799 * t710 - t790 * t711 - t796 * t751;
t794 = -t790 * t710 + t800 * t711 + t798 * t751;
t760 = qJD(1) ^ 2;
t792 = pkin(2) * t754;
t791 = -mrSges(5,3) - mrSges(6,2);
t789 = mrSges(3,2) * t753;
t750 = t754 ^ 2;
t788 = t750 * t760;
t757 = sin(qJ(1));
t759 = cos(qJ(1));
t736 = -t759 * g(1) - t757 * g(2);
t731 = -t760 * pkin(1) + qJDD(1) * qJ(2) + t736;
t782 = qJD(1) * qJD(2);
t779 = -t754 * g(3) - 0.2e1 * t753 * t782;
t705 = (-pkin(6) * qJDD(1) + t760 * t792 - t731) * t753 + t779;
t721 = -t753 * g(3) + (t731 + 0.2e1 * t782) * t754;
t781 = qJDD(1) * t754;
t706 = -pkin(2) * t788 + pkin(6) * t781 + t721;
t681 = t758 * t705 - t756 * t706;
t783 = t729 * qJD(3);
t719 = t769 * qJDD(1) + t783;
t664 = (-t719 + t783) * pkin(7) + (t729 * t730 + qJDD(3)) * pkin(3) + t681;
t682 = t756 * t705 + t758 * t706;
t718 = -t730 * qJD(3) + t768 * qJDD(1);
t724 = qJD(3) * pkin(3) - t730 * pkin(7);
t728 = t729 ^ 2;
t666 = -t728 * pkin(3) + t718 * pkin(7) - qJD(3) * t724 + t682;
t662 = t755 * t664 + t793 * t666;
t679 = t711 * qJD(4) - t793 * t718 + t755 * t719;
t701 = t751 * mrSges(5,1) - t711 * mrSges(5,3);
t748 = qJDD(3) + qJDD(4);
t693 = t710 * pkin(4) - t711 * qJ(5);
t747 = t751 ^ 2;
t656 = -t747 * pkin(4) + t748 * qJ(5) + 0.2e1 * qJD(5) * t751 - t710 * t693 + t662;
t702 = -t751 * mrSges(6,1) + t711 * mrSges(6,2);
t780 = m(6) * t656 + t748 * mrSges(6,3) + t751 * t702;
t694 = t710 * mrSges(6,1) - t711 * mrSges(6,3);
t785 = -t710 * mrSges(5,1) - t711 * mrSges(5,2) - t694;
t646 = m(5) * t662 - t748 * mrSges(5,2) + t791 * t679 - t751 * t701 + t785 * t710 + t780;
t661 = t793 * t664 - t755 * t666;
t680 = -t710 * qJD(4) + t755 * t718 + t793 * t719;
t700 = -t751 * mrSges(5,2) - t710 * mrSges(5,3);
t657 = -t748 * pkin(4) - t747 * qJ(5) + t711 * t693 + qJDD(5) - t661;
t703 = -t710 * mrSges(6,2) + t751 * mrSges(6,3);
t774 = -m(6) * t657 + t748 * mrSges(6,1) + t751 * t703;
t648 = m(5) * t661 + t748 * mrSges(5,1) + t791 * t680 + t751 * t700 + t785 * t711 + t774;
t639 = t755 * t646 + t793 * t648;
t715 = -t729 * mrSges(4,1) + t730 * mrSges(4,2);
t722 = -qJD(3) * mrSges(4,2) + t729 * mrSges(4,3);
t637 = m(4) * t681 + qJDD(3) * mrSges(4,1) - t719 * mrSges(4,3) + qJD(3) * t722 - t730 * t715 + t639;
t723 = qJD(3) * mrSges(4,1) - t730 * mrSges(4,3);
t775 = t793 * t646 - t755 * t648;
t638 = m(4) * t682 - qJDD(3) * mrSges(4,2) + t718 * mrSges(4,3) - qJD(3) * t723 + t729 * t715 + t775;
t631 = t758 * t637 + t756 * t638;
t720 = -t753 * t731 + t779;
t767 = mrSges(3,3) * qJDD(1) + t760 * (-mrSges(3,1) * t754 + t789);
t629 = m(3) * t720 - t767 * t753 + t631;
t776 = -t756 * t637 + t758 * t638;
t630 = m(3) * t721 + t767 * t754 + t776;
t777 = -t753 * t629 + t754 * t630;
t621 = m(2) * t736 - t760 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t777;
t735 = t757 * g(1) - t759 * g(2);
t773 = qJDD(2) - t735;
t727 = -qJDD(1) * pkin(1) - t760 * qJ(2) + t773;
t749 = t753 ^ 2;
t717 = (-pkin(1) - t792) * qJDD(1) + (-qJ(2) + (-t749 - t750) * pkin(6)) * t760 + t773;
t670 = -t718 * pkin(3) - t728 * pkin(7) + t730 * t724 + t717;
t659 = -0.2e1 * qJD(5) * t711 + (t710 * t751 - t680) * qJ(5) + (t711 * t751 + t679) * pkin(4) + t670;
t649 = m(6) * t659 + t679 * mrSges(6,1) - t680 * mrSges(6,3) - t711 * t702 + t710 * t703;
t765 = m(5) * t670 + t679 * mrSges(5,1) + t680 * mrSges(5,2) + t710 * t700 + t711 * t701 + t649;
t764 = m(4) * t717 - t718 * mrSges(4,1) + t719 * mrSges(4,2) - t729 * t722 + t730 * t723 + t765;
t762 = -m(3) * t727 + mrSges(3,1) * t781 - t764 + (t749 * t760 + t788) * mrSges(3,3);
t641 = t762 - t760 * mrSges(2,2) + m(2) * t735 + (mrSges(2,1) - t789) * qJDD(1);
t787 = t757 * t621 + t759 * t641;
t623 = t754 * t629 + t753 * t630;
t786 = t796 * t710 - t798 * t711 + t797 * t751;
t770 = Ifges(3,5) * t753 + Ifges(3,6) * t754;
t784 = t760 * t770;
t778 = t759 * t621 - t757 * t641;
t772 = Ifges(3,1) * t753 + Ifges(3,4) * t754;
t771 = Ifges(3,4) * t753 + Ifges(3,2) * t754;
t632 = -mrSges(5,1) * t670 - mrSges(6,1) * t659 + mrSges(6,2) * t656 + mrSges(5,3) * t662 - pkin(4) * t649 - t799 * t679 + t790 * t680 + t786 * t711 + t796 * t748 + t794 * t751;
t633 = mrSges(5,2) * t670 + mrSges(6,2) * t657 - mrSges(5,3) * t661 - mrSges(6,3) * t659 - qJ(5) * t649 - t790 * t679 + t800 * t680 + t786 * t710 + t798 * t748 + t795 * t751;
t707 = Ifges(4,5) * t730 + Ifges(4,6) * t729 + Ifges(4,3) * qJD(3);
t709 = Ifges(4,1) * t730 + Ifges(4,4) * t729 + Ifges(4,5) * qJD(3);
t624 = -mrSges(4,1) * t717 + mrSges(4,3) * t682 + Ifges(4,4) * t719 + Ifges(4,2) * t718 + Ifges(4,6) * qJDD(3) - pkin(3) * t765 + pkin(7) * t775 + qJD(3) * t709 + t793 * t632 + t755 * t633 - t730 * t707;
t708 = Ifges(4,4) * t730 + Ifges(4,2) * t729 + Ifges(4,6) * qJD(3);
t625 = mrSges(4,2) * t717 - mrSges(4,3) * t681 + Ifges(4,1) * t719 + Ifges(4,4) * t718 + Ifges(4,5) * qJDD(3) - pkin(7) * t639 - qJD(3) * t708 - t755 * t632 + t793 * t633 + t729 * t707;
t615 = -mrSges(3,1) * t727 + mrSges(3,3) * t721 - pkin(2) * t764 + pkin(6) * t776 + t771 * qJDD(1) + t758 * t624 + t756 * t625 - t753 * t784;
t617 = mrSges(3,2) * t727 - mrSges(3,3) * t720 - pkin(6) * t631 + t772 * qJDD(1) - t756 * t624 + t758 * t625 + t754 * t784;
t643 = qJDD(1) * t789 - t762;
t766 = mrSges(2,1) * t735 - mrSges(2,2) * t736 + Ifges(2,3) * qJDD(1) - pkin(1) * t643 + qJ(2) * t777 + t754 * t615 + t753 * t617;
t653 = t680 * mrSges(6,2) + t711 * t694 - t774;
t763 = -mrSges(5,1) * t661 + mrSges(6,1) * t657 + mrSges(5,2) * t662 - mrSges(6,3) * t656 + pkin(4) * t653 - qJ(5) * t780 + t797 * t748 + t795 * t711 + (qJ(5) * t694 - t794) * t710 - t798 * t680 + (qJ(5) * mrSges(6,2) + t796) * t679;
t761 = mrSges(4,1) * t681 - mrSges(4,2) * t682 + Ifges(4,5) * t719 + Ifges(4,6) * t718 + Ifges(4,3) * qJDD(3) + pkin(3) * t639 + t730 * t708 - t729 * t709 - t763;
t618 = -t761 + mrSges(2,1) * g(3) + (Ifges(2,6) - t770) * qJDD(1) + mrSges(2,3) * t736 - mrSges(3,1) * t720 + mrSges(3,2) * t721 - pkin(2) * t631 - pkin(1) * t623 + (-t753 * t771 + t754 * t772 + Ifges(2,5)) * t760;
t613 = -mrSges(2,2) * g(3) - mrSges(2,3) * t735 + Ifges(2,5) * qJDD(1) - t760 * Ifges(2,6) - qJ(2) * t623 - t753 * t615 + t754 * t617;
t1 = [-m(1) * g(1) + t778; -m(1) * g(2) + t787; (-m(1) - m(2)) * g(3) + t623; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t787 + t759 * t613 - t757 * t618; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t778 + t757 * t613 + t759 * t618; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t766; t766; t643; t761; -t763; t653;];
tauJB = t1;
