% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRR2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:03
% EndTime: 2019-12-05 18:12:12
% DurationCPUTime: 9.67s
% Computational Cost: add. (139313->293), mult. (345740->366), div. (0->0), fcn. (262413->10), ass. (0->127)
t769 = qJD(1) ^ 2;
t760 = cos(pkin(9));
t798 = pkin(2) * t760;
t759 = sin(pkin(9));
t797 = mrSges(3,2) * t759;
t756 = t760 ^ 2;
t796 = t756 * t769;
t764 = sin(qJ(1));
t768 = cos(qJ(1));
t743 = -t768 * g(1) - t764 * g(2);
t738 = -t769 * pkin(1) + qJDD(1) * qJ(2) + t743;
t792 = qJD(1) * qJD(2);
t790 = -t760 * g(3) - 0.2e1 * t759 * t792;
t713 = (-pkin(6) * qJDD(1) + t769 * t798 - t738) * t759 + t790;
t729 = -t759 * g(3) + (t738 + 0.2e1 * t792) * t760;
t791 = qJDD(1) * t760;
t714 = -pkin(2) * t796 + pkin(6) * t791 + t729;
t763 = sin(qJ(3));
t767 = cos(qJ(3));
t695 = t767 * t713 - t763 * t714;
t779 = t759 * t767 + t760 * t763;
t778 = -t759 * t763 + t760 * t767;
t736 = t778 * qJD(1);
t793 = t736 * qJD(3);
t727 = t779 * qJDD(1) + t793;
t737 = t779 * qJD(1);
t677 = (-t727 + t793) * pkin(7) + (t736 * t737 + qJDD(3)) * pkin(3) + t695;
t696 = t763 * t713 + t767 * t714;
t726 = -t737 * qJD(3) + t778 * qJDD(1);
t732 = qJD(3) * pkin(3) - t737 * pkin(7);
t735 = t736 ^ 2;
t685 = -t735 * pkin(3) + t726 * pkin(7) - qJD(3) * t732 + t696;
t762 = sin(qJ(4));
t766 = cos(qJ(4));
t666 = t766 * t677 - t762 * t685;
t719 = t766 * t736 - t762 * t737;
t694 = t719 * qJD(4) + t762 * t726 + t766 * t727;
t720 = t762 * t736 + t766 * t737;
t754 = qJDD(3) + qJDD(4);
t757 = qJD(3) + qJD(4);
t661 = (t719 * t757 - t694) * pkin(8) + (t719 * t720 + t754) * pkin(4) + t666;
t667 = t762 * t677 + t766 * t685;
t693 = -t720 * qJD(4) + t766 * t726 - t762 * t727;
t711 = t757 * pkin(4) - t720 * pkin(8);
t715 = t719 ^ 2;
t662 = -t715 * pkin(4) + t693 * pkin(8) - t757 * t711 + t667;
t761 = sin(qJ(5));
t765 = cos(qJ(5));
t659 = t765 * t661 - t761 * t662;
t704 = t765 * t719 - t761 * t720;
t673 = t704 * qJD(5) + t761 * t693 + t765 * t694;
t705 = t761 * t719 + t765 * t720;
t683 = -t704 * mrSges(6,1) + t705 * mrSges(6,2);
t752 = qJD(5) + t757;
t697 = -t752 * mrSges(6,2) + t704 * mrSges(6,3);
t751 = qJDD(5) + t754;
t656 = m(6) * t659 + t751 * mrSges(6,1) - t673 * mrSges(6,3) - t705 * t683 + t752 * t697;
t660 = t761 * t661 + t765 * t662;
t672 = -t705 * qJD(5) + t765 * t693 - t761 * t694;
t698 = t752 * mrSges(6,1) - t705 * mrSges(6,3);
t657 = m(6) * t660 - t751 * mrSges(6,2) + t672 * mrSges(6,3) + t704 * t683 - t752 * t698;
t646 = t765 * t656 + t761 * t657;
t706 = -t719 * mrSges(5,1) + t720 * mrSges(5,2);
t709 = -t757 * mrSges(5,2) + t719 * mrSges(5,3);
t643 = m(5) * t666 + t754 * mrSges(5,1) - t694 * mrSges(5,3) - t720 * t706 + t757 * t709 + t646;
t710 = t757 * mrSges(5,1) - t720 * mrSges(5,3);
t785 = -t761 * t656 + t765 * t657;
t644 = m(5) * t667 - t754 * mrSges(5,2) + t693 * mrSges(5,3) + t719 * t706 - t757 * t710 + t785;
t639 = t766 * t643 + t762 * t644;
t723 = -t736 * mrSges(4,1) + t737 * mrSges(4,2);
t730 = -qJD(3) * mrSges(4,2) + t736 * mrSges(4,3);
t637 = m(4) * t695 + qJDD(3) * mrSges(4,1) - t727 * mrSges(4,3) + qJD(3) * t730 - t737 * t723 + t639;
t731 = qJD(3) * mrSges(4,1) - t737 * mrSges(4,3);
t786 = -t762 * t643 + t766 * t644;
t638 = m(4) * t696 - qJDD(3) * mrSges(4,2) + t726 * mrSges(4,3) - qJD(3) * t731 + t736 * t723 + t786;
t631 = t767 * t637 + t763 * t638;
t728 = -t759 * t738 + t790;
t777 = mrSges(3,3) * qJDD(1) + t769 * (-mrSges(3,1) * t760 + t797);
t629 = m(3) * t728 - t777 * t759 + t631;
t787 = -t763 * t637 + t767 * t638;
t630 = m(3) * t729 + t777 * t760 + t787;
t788 = -t759 * t629 + t760 * t630;
t621 = m(2) * t743 - t769 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t788;
t742 = t764 * g(1) - t768 * g(2);
t784 = qJDD(2) - t742;
t734 = -qJDD(1) * pkin(1) - t769 * qJ(2) + t784;
t755 = t759 ^ 2;
t725 = (-pkin(1) - t798) * qJDD(1) + (-qJ(2) + (-t755 - t756) * pkin(6)) * t769 + t784;
t688 = -t726 * pkin(3) - t735 * pkin(7) + t737 * t732 + t725;
t664 = -t693 * pkin(4) - t715 * pkin(8) + t720 * t711 + t688;
t780 = m(6) * t664 - t672 * mrSges(6,1) + t673 * mrSges(6,2) - t704 * t697 + t705 * t698;
t775 = m(5) * t688 - t693 * mrSges(5,1) + t694 * mrSges(5,2) - t719 * t709 + t720 * t710 + t780;
t773 = m(4) * t725 - t726 * mrSges(4,1) + t727 * mrSges(4,2) - t736 * t730 + t737 * t731 + t775;
t771 = -m(3) * t734 + mrSges(3,1) * t791 - t773 + (t755 * t769 + t796) * mrSges(3,3);
t650 = t771 + (mrSges(2,1) - t797) * qJDD(1) - t769 * mrSges(2,2) + m(2) * t742;
t795 = t764 * t621 + t768 * t650;
t623 = t760 * t629 + t759 * t630;
t781 = Ifges(3,5) * t759 + Ifges(3,6) * t760;
t794 = t769 * t781;
t789 = t768 * t621 - t764 * t650;
t783 = Ifges(3,1) * t759 + Ifges(3,4) * t760;
t782 = Ifges(3,4) * t759 + Ifges(3,2) * t760;
t678 = Ifges(6,5) * t705 + Ifges(6,6) * t704 + Ifges(6,3) * t752;
t680 = Ifges(6,1) * t705 + Ifges(6,4) * t704 + Ifges(6,5) * t752;
t647 = -mrSges(6,1) * t664 + mrSges(6,3) * t660 + Ifges(6,4) * t673 + Ifges(6,2) * t672 + Ifges(6,6) * t751 - t705 * t678 + t752 * t680;
t679 = Ifges(6,4) * t705 + Ifges(6,2) * t704 + Ifges(6,6) * t752;
t648 = mrSges(6,2) * t664 - mrSges(6,3) * t659 + Ifges(6,1) * t673 + Ifges(6,4) * t672 + Ifges(6,5) * t751 + t704 * t678 - t752 * t679;
t699 = Ifges(5,5) * t720 + Ifges(5,6) * t719 + Ifges(5,3) * t757;
t701 = Ifges(5,1) * t720 + Ifges(5,4) * t719 + Ifges(5,5) * t757;
t632 = -mrSges(5,1) * t688 + mrSges(5,3) * t667 + Ifges(5,4) * t694 + Ifges(5,2) * t693 + Ifges(5,6) * t754 - pkin(4) * t780 + pkin(8) * t785 + t765 * t647 + t761 * t648 - t720 * t699 + t757 * t701;
t700 = Ifges(5,4) * t720 + Ifges(5,2) * t719 + Ifges(5,6) * t757;
t633 = mrSges(5,2) * t688 - mrSges(5,3) * t666 + Ifges(5,1) * t694 + Ifges(5,4) * t693 + Ifges(5,5) * t754 - pkin(8) * t646 - t761 * t647 + t765 * t648 + t719 * t699 - t757 * t700;
t716 = Ifges(4,5) * t737 + Ifges(4,6) * t736 + Ifges(4,3) * qJD(3);
t718 = Ifges(4,1) * t737 + Ifges(4,4) * t736 + Ifges(4,5) * qJD(3);
t624 = -mrSges(4,1) * t725 + mrSges(4,3) * t696 + Ifges(4,4) * t727 + Ifges(4,2) * t726 + Ifges(4,6) * qJDD(3) - pkin(3) * t775 + pkin(7) * t786 + qJD(3) * t718 + t766 * t632 + t762 * t633 - t737 * t716;
t717 = Ifges(4,4) * t737 + Ifges(4,2) * t736 + Ifges(4,6) * qJD(3);
t625 = mrSges(4,2) * t725 - mrSges(4,3) * t695 + Ifges(4,1) * t727 + Ifges(4,4) * t726 + Ifges(4,5) * qJDD(3) - pkin(7) * t639 - qJD(3) * t717 - t762 * t632 + t766 * t633 + t736 * t716;
t615 = -mrSges(3,1) * t734 + mrSges(3,3) * t729 - pkin(2) * t773 + pkin(6) * t787 + t782 * qJDD(1) + t767 * t624 + t763 * t625 - t759 * t794;
t617 = mrSges(3,2) * t734 - mrSges(3,3) * t728 - pkin(6) * t631 + t783 * qJDD(1) - t763 * t624 + t767 * t625 + t760 * t794;
t652 = qJDD(1) * t797 - t771;
t776 = mrSges(2,1) * t742 - mrSges(2,2) * t743 + Ifges(2,3) * qJDD(1) - pkin(1) * t652 + qJ(2) * t788 + t760 * t615 + t759 * t617;
t774 = -mrSges(6,1) * t659 + mrSges(6,2) * t660 - Ifges(6,5) * t673 - Ifges(6,6) * t672 - Ifges(6,3) * t751 - t705 * t679 + t704 * t680;
t772 = -mrSges(5,1) * t666 + mrSges(5,2) * t667 - Ifges(5,5) * t694 - Ifges(5,6) * t693 - Ifges(5,3) * t754 - pkin(4) * t646 - t720 * t700 + t719 * t701 + t774;
t770 = mrSges(4,1) * t695 - mrSges(4,2) * t696 + Ifges(4,5) * t727 + Ifges(4,6) * t726 + Ifges(4,3) * qJDD(3) + pkin(3) * t639 + t737 * t717 - t736 * t718 - t772;
t618 = -t770 - pkin(1) * t623 + mrSges(2,1) * g(3) + (Ifges(2,6) - t781) * qJDD(1) + mrSges(2,3) * t743 - mrSges(3,1) * t728 + mrSges(3,2) * t729 - pkin(2) * t631 + (-t759 * t782 + t760 * t783 + Ifges(2,5)) * t769;
t613 = -mrSges(2,2) * g(3) - mrSges(2,3) * t742 + Ifges(2,5) * qJDD(1) - t769 * Ifges(2,6) - qJ(2) * t623 - t759 * t615 + t760 * t617;
t1 = [-m(1) * g(1) + t789; -m(1) * g(2) + t795; (-m(1) - m(2)) * g(3) + t623; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t795 + t768 * t613 - t764 * t618; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t789 + t764 * t613 + t768 * t618; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t776; t776; t652; t770; -t772; -t774;];
tauJB = t1;
