% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PPRRRP2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-05-04 20:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PPRRRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:32:37
% EndTime: 2019-05-04 20:32:48
% DurationCPUTime: 10.48s
% Computational Cost: add. (186281->270), mult. (332625->344), div. (0->0), fcn. (253666->14), ass. (0->128)
t807 = Ifges(6,1) + Ifges(7,1);
t799 = Ifges(6,4) - Ifges(7,5);
t798 = -Ifges(6,5) - Ifges(7,4);
t806 = Ifges(6,2) + Ifges(7,3);
t797 = Ifges(6,6) - Ifges(7,6);
t805 = -Ifges(6,3) - Ifges(7,2);
t754 = sin(pkin(11));
t758 = cos(pkin(11));
t744 = -g(1) * t758 - g(2) * t754;
t753 = sin(pkin(12));
t757 = cos(pkin(12));
t743 = g(1) * t754 - g(2) * t758;
t752 = -g(3) + qJDD(1);
t756 = sin(pkin(6));
t760 = cos(pkin(6));
t773 = t743 * t760 + t752 * t756;
t703 = -t753 * t744 + t757 * t773;
t704 = t757 * t744 + t753 * t773;
t727 = -t743 * t756 + t752 * t760 + qJDD(2);
t765 = cos(qJ(3));
t759 = cos(pkin(7));
t763 = sin(qJ(3));
t791 = t759 * t763;
t755 = sin(pkin(7));
t792 = t755 * t763;
t697 = t703 * t791 + t765 * t704 + t727 * t792;
t767 = qJD(3) ^ 2;
t695 = -pkin(3) * t767 + qJDD(3) * pkin(9) + t697;
t699 = -t703 * t755 + t727 * t759;
t762 = sin(qJ(4));
t764 = cos(qJ(4));
t689 = t764 * t695 + t762 * t699;
t740 = (-pkin(4) * t764 - pkin(10) * t762) * qJD(3);
t766 = qJD(4) ^ 2;
t784 = qJD(3) * t764;
t687 = -pkin(4) * t766 + qJDD(4) * pkin(10) + t740 * t784 + t689;
t696 = -t763 * t704 + (t703 * t759 + t727 * t755) * t765;
t694 = -qJDD(3) * pkin(3) - t767 * pkin(9) - t696;
t783 = qJD(3) * qJD(4);
t780 = t764 * t783;
t741 = qJDD(3) * t762 + t780;
t781 = t762 * t783;
t742 = qJDD(3) * t764 - t781;
t691 = (-t741 - t780) * pkin(10) + (-t742 + t781) * pkin(4) + t694;
t761 = sin(qJ(5));
t801 = cos(qJ(5));
t683 = t687 * t801 + t761 * t691;
t785 = qJD(3) * t762;
t738 = t761 * qJD(4) + t785 * t801;
t715 = qJD(5) * t738 - qJDD(4) * t801 + t741 * t761;
t748 = qJD(5) - t784;
t724 = mrSges(6,1) * t748 - mrSges(6,3) * t738;
t736 = qJDD(5) - t742;
t737 = -qJD(4) * t801 + t761 * t785;
t719 = pkin(5) * t737 - qJ(6) * t738;
t747 = t748 ^ 2;
t680 = -pkin(5) * t747 + qJ(6) * t736 + 0.2e1 * qJD(6) * t748 - t719 * t737 + t683;
t725 = -mrSges(7,1) * t748 + mrSges(7,2) * t738;
t782 = m(7) * t680 + t736 * mrSges(7,3) + t748 * t725;
t720 = mrSges(7,1) * t737 - mrSges(7,3) * t738;
t786 = -mrSges(6,1) * t737 - mrSges(6,2) * t738 - t720;
t800 = -mrSges(6,3) - mrSges(7,2);
t673 = m(6) * t683 - t736 * mrSges(6,2) + t715 * t800 - t748 * t724 + t737 * t786 + t782;
t682 = -t761 * t687 + t691 * t801;
t716 = -t737 * qJD(5) + t761 * qJDD(4) + t741 * t801;
t723 = -mrSges(6,2) * t748 - mrSges(6,3) * t737;
t681 = -t736 * pkin(5) - t747 * qJ(6) + t738 * t719 + qJDD(6) - t682;
t726 = -mrSges(7,2) * t737 + mrSges(7,3) * t748;
t776 = -m(7) * t681 + t736 * mrSges(7,1) + t748 * t726;
t674 = m(6) * t682 + t736 * mrSges(6,1) + t716 * t800 + t748 * t723 + t738 * t786 + t776;
t669 = t673 * t801 - t674 * t761;
t739 = (-mrSges(5,1) * t764 + mrSges(5,2) * t762) * qJD(3);
t745 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t785;
t665 = m(5) * t689 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t742 - qJD(4) * t745 + t739 * t784 + t669;
t688 = -t762 * t695 + t764 * t699;
t686 = -qJDD(4) * pkin(4) - t766 * pkin(10) + t740 * t785 - t688;
t684 = -0.2e1 * qJD(6) * t738 + (t737 * t748 - t716) * qJ(6) + (t738 * t748 + t715) * pkin(5) + t686;
t678 = m(7) * t684 + mrSges(7,1) * t715 - t716 * mrSges(7,3) - t738 * t725 + t726 * t737;
t675 = -m(6) * t686 - t715 * mrSges(6,1) - mrSges(6,2) * t716 - t737 * t723 - t724 * t738 - t678;
t746 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t784;
t671 = m(5) * t688 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t741 + qJD(4) * t746 - t739 * t785 + t675;
t778 = t764 * t665 - t671 * t762;
t656 = m(4) * t697 - mrSges(4,1) * t767 - qJDD(3) * mrSges(4,2) + t778;
t659 = t762 * t665 + t764 * t671;
t658 = m(4) * t699 + t659;
t668 = t761 * t673 + t674 * t801;
t769 = -m(5) * t694 + t742 * mrSges(5,1) - t741 * mrSges(5,2) - t745 * t785 + t746 * t784 - t668;
t662 = m(4) * t696 + qJDD(3) * mrSges(4,1) - t767 * mrSges(4,2) + t769;
t793 = t662 * t765;
t645 = t656 * t791 - t658 * t755 + t759 * t793;
t641 = m(3) * t703 + t645;
t651 = t765 * t656 - t662 * t763;
t650 = m(3) * t704 + t651;
t804 = t641 * t757 + t650 * t753;
t787 = t799 * t737 - t738 * t807 + t798 * t748;
t789 = t737 * t797 + t738 * t798 + t748 * t805;
t666 = -mrSges(6,1) * t686 - mrSges(7,1) * t684 + mrSges(7,2) * t680 + mrSges(6,3) * t683 - pkin(5) * t678 - t715 * t806 + t799 * t716 + t797 * t736 + t789 * t738 - t787 * t748;
t788 = t737 * t806 - t799 * t738 - t748 * t797;
t667 = mrSges(6,2) * t686 + mrSges(7,2) * t681 - mrSges(6,3) * t682 - mrSges(7,3) * t684 - qJ(6) * t678 - t799 * t715 + t716 * t807 - t798 * t736 + t789 * t737 + t788 * t748;
t730 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t762 + Ifges(5,2) * t764) * qJD(3);
t731 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t762 + Ifges(5,4) * t764) * qJD(3);
t803 = mrSges(5,1) * t688 - mrSges(5,2) * t689 + Ifges(5,5) * t741 + Ifges(5,6) * t742 + Ifges(5,3) * qJDD(4) + pkin(4) * t675 + pkin(10) * t669 + (t730 * t762 - t731 * t764) * qJD(3) + t666 * t801 + t761 * t667;
t677 = t716 * mrSges(7,2) + t738 * t720 - t776;
t802 = -t715 * t797 - t716 * t798 - t805 * t736 - t737 * t787 - t738 * t788 + mrSges(6,1) * t682 - mrSges(7,1) * t681 - mrSges(6,2) * t683 + mrSges(7,3) * t680 - pkin(5) * t677 + qJ(6) * (-t715 * mrSges(7,2) - t737 * t720 + t782);
t644 = t656 * t792 + t759 * t658 + t755 * t793;
t643 = m(3) * t727 + t644;
t631 = -t643 * t756 + t760 * t804;
t629 = m(2) * t743 + t631;
t637 = -t641 * t753 + t757 * t650;
t636 = m(2) * t744 + t637;
t790 = t758 * t629 + t754 * t636;
t630 = t760 * t643 + t756 * t804;
t779 = -t629 * t754 + t758 * t636;
t777 = m(2) * t752 + t630;
t729 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t762 + Ifges(5,6) * t764) * qJD(3);
t646 = mrSges(5,2) * t694 - mrSges(5,3) * t688 + Ifges(5,1) * t741 + Ifges(5,4) * t742 + Ifges(5,5) * qJDD(4) - pkin(10) * t668 - qJD(4) * t730 - t761 * t666 + t667 * t801 + t729 * t784;
t652 = -mrSges(5,1) * t694 + mrSges(5,3) * t689 + Ifges(5,4) * t741 + Ifges(5,2) * t742 + Ifges(5,6) * qJDD(4) - pkin(4) * t668 + qJD(4) * t731 - t729 * t785 - t802;
t633 = mrSges(4,2) * t699 - mrSges(4,3) * t696 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t767 - pkin(9) * t659 + t646 * t764 - t652 * t762;
t638 = -mrSges(4,1) * t699 + mrSges(4,3) * t697 + t767 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t659 - t803;
t772 = pkin(8) * t651 + t633 * t763 + t638 * t765;
t632 = mrSges(4,1) * t696 - mrSges(4,2) * t697 + Ifges(4,3) * qJDD(3) + pkin(3) * t769 + pkin(9) * t778 + t762 * t646 + t764 * t652;
t626 = -mrSges(3,1) * t727 + mrSges(3,3) * t704 - pkin(2) * t644 - t755 * t632 + t759 * t772;
t627 = mrSges(3,2) * t727 - mrSges(3,3) * t703 + t765 * t633 - t763 * t638 + (-t644 * t755 - t645 * t759) * pkin(8);
t771 = qJ(2) * t637 + t626 * t757 + t627 * t753;
t625 = mrSges(3,1) * t703 - mrSges(3,2) * t704 + pkin(2) * t645 + t759 * t632 + t755 * t772;
t624 = mrSges(2,2) * t752 - mrSges(2,3) * t743 - t753 * t626 + t757 * t627 + (-t630 * t756 - t631 * t760) * qJ(2);
t623 = -mrSges(2,1) * t752 + mrSges(2,3) * t744 - pkin(1) * t630 - t756 * t625 + t760 * t771;
t1 = [-m(1) * g(1) + t779; -m(1) * g(2) + t790; -m(1) * g(3) + t777; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t790 - t754 * t623 + t758 * t624; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t779 + t758 * t623 + t754 * t624; -mrSges(1,1) * g(2) + mrSges(2,1) * t743 + mrSges(1,2) * g(1) - mrSges(2,2) * t744 + pkin(1) * t631 + t760 * t625 + t756 * t771; t777; t643; t632; t803; t802; t677;];
tauJB  = t1;
