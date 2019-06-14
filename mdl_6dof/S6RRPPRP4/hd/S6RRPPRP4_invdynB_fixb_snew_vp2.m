% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-05-06 09:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:19:09
% EndTime: 2019-05-06 09:19:19
% DurationCPUTime: 5.70s
% Computational Cost: add. (53722->339), mult. (117182->397), div. (0->0), fcn. (75223->8), ass. (0->131)
t807 = -2 * qJD(3);
t806 = -2 * qJD(4);
t805 = Ifges(4,1) + Ifges(5,1);
t804 = Ifges(6,1) + Ifges(7,1);
t794 = Ifges(4,4) - Ifges(5,5);
t793 = Ifges(6,4) - Ifges(7,5);
t792 = Ifges(7,4) + Ifges(6,5);
t803 = Ifges(4,5) + Ifges(5,4);
t802 = -Ifges(4,2) - Ifges(5,3);
t801 = Ifges(5,2) + Ifges(4,3);
t800 = Ifges(6,2) + Ifges(7,3);
t799 = Ifges(4,6) - Ifges(5,6);
t789 = Ifges(6,6) - Ifges(7,6);
t798 = Ifges(6,3) + Ifges(7,2);
t753 = sin(qJ(1));
t755 = cos(qJ(1));
t738 = -g(1) * t755 - g(2) * t753;
t757 = qJD(1) ^ 2;
t719 = -pkin(1) * t757 + qJDD(1) * pkin(7) + t738;
t752 = sin(qJ(2));
t754 = cos(qJ(2));
t694 = -t754 * g(3) - t752 * t719;
t731 = (-pkin(2) * t754 - qJ(3) * t752) * qJD(1);
t756 = qJD(2) ^ 2;
t777 = qJD(1) * t752;
t671 = -qJDD(2) * pkin(2) - t756 * qJ(3) + t731 * t777 + qJDD(3) - t694;
t773 = qJD(1) * qJD(2);
t770 = t754 * t773;
t733 = qJDD(1) * t752 + t770;
t750 = sin(pkin(9));
t788 = cos(pkin(9));
t704 = -qJDD(2) * t788 + t733 * t750;
t705 = t750 * qJDD(2) + t733 * t788;
t726 = t750 * qJD(2) + t777 * t788;
t725 = -qJD(2) * t788 + t750 * t777;
t776 = qJD(1) * t754;
t771 = t725 * t776;
t643 = t671 - (t705 + t771) * qJ(4) - (t726 * t776 - t704) * pkin(3) + t726 * t806;
t737 = t753 * g(1) - t755 * g(2);
t718 = -qJDD(1) * pkin(1) - t757 * pkin(7) - t737;
t742 = t752 * t773;
t734 = qJDD(1) * t754 - t742;
t666 = (-t733 - t770) * qJ(3) + (-t734 + t742) * pkin(2) + t718;
t695 = -g(3) * t752 + t754 * t719;
t672 = -pkin(2) * t756 + qJDD(2) * qJ(3) + t731 * t776 + t695;
t644 = t666 * t788 - t750 * t672 + t726 * t807;
t797 = cos(qJ(5));
t796 = -mrSges(4,3) - mrSges(5,2);
t795 = -mrSges(6,3) - mrSges(7,2);
t787 = t754 ^ 2 * t757;
t732 = (-mrSges(3,1) * t754 + mrSges(3,2) * t752) * qJD(1);
t735 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t777;
t645 = t750 * t666 + t788 * t672 + t725 * t807;
t702 = -mrSges(4,1) * t776 - mrSges(4,3) * t726;
t703 = mrSges(5,1) * t776 + mrSges(5,2) * t726;
t690 = pkin(3) * t725 - qJ(4) * t726;
t641 = -pkin(3) * t787 - t734 * qJ(4) - t725 * t690 + t776 * t806 + t645;
t642 = t734 * pkin(3) - qJ(4) * t787 + t726 * t690 + qJDD(4) - t644;
t635 = (-t705 + t771) * pkin(8) + (t725 * t726 + t734) * pkin(4) + t642;
t706 = pkin(4) * t776 - pkin(8) * t726;
t723 = t725 ^ 2;
t637 = -pkin(4) * t723 + pkin(8) * t704 - t706 * t776 + t641;
t751 = sin(qJ(5));
t631 = t751 * t635 + t637 * t797;
t689 = t751 * t725 + t726 * t797;
t650 = qJD(5) * t689 - t704 * t797 + t705 * t751;
t740 = qJD(5) + t776;
t674 = mrSges(6,1) * t740 - mrSges(6,3) * t689;
t688 = -t725 * t797 + t726 * t751;
t730 = qJDD(5) + t734;
t662 = pkin(5) * t688 - qJ(6) * t689;
t739 = t740 ^ 2;
t628 = -pkin(5) * t739 + qJ(6) * t730 + 0.2e1 * qJD(6) * t740 - t662 * t688 + t631;
t675 = -mrSges(7,1) * t740 + mrSges(7,2) * t689;
t772 = m(7) * t628 + t730 * mrSges(7,3) + t740 * t675;
t663 = mrSges(7,1) * t688 - mrSges(7,3) * t689;
t782 = -mrSges(6,1) * t688 - mrSges(6,2) * t689 - t663;
t623 = m(6) * t631 - t730 * mrSges(6,2) + t650 * t795 - t740 * t674 + t688 * t782 + t772;
t630 = t635 * t797 - t751 * t637;
t651 = -t688 * qJD(5) + t751 * t704 + t705 * t797;
t673 = -mrSges(6,2) * t740 - mrSges(6,3) * t688;
t629 = -t730 * pkin(5) - t739 * qJ(6) + t689 * t662 + qJDD(6) - t630;
t676 = -mrSges(7,2) * t688 + mrSges(7,3) * t740;
t765 = -m(7) * t629 + t730 * mrSges(7,1) + t740 * t676;
t624 = m(6) * t630 + t730 * mrSges(6,1) + t651 * t795 + t740 * t673 + t689 * t782 + t765;
t766 = t623 * t797 - t751 * t624;
t763 = m(5) * t641 - t734 * mrSges(5,3) + t766;
t691 = mrSges(5,1) * t725 - mrSges(5,3) * t726;
t778 = -mrSges(4,1) * t725 - mrSges(4,2) * t726 - t691;
t615 = m(4) * t645 + t734 * mrSges(4,2) + t778 * t725 + t796 * t704 + (t702 - t703) * t776 + t763;
t700 = -mrSges(5,2) * t725 - mrSges(5,3) * t776;
t701 = mrSges(4,2) * t776 - mrSges(4,3) * t725;
t619 = t751 * t623 + t624 * t797;
t762 = -m(5) * t642 - t734 * mrSges(5,1) - t619;
t616 = m(4) * t644 - t734 * mrSges(4,1) + t778 * t726 + t796 * t705 + (-t700 - t701) * t776 + t762;
t767 = t788 * t615 - t616 * t750;
t612 = m(3) * t695 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t734 - qJD(2) * t735 + t732 * t776 + t767;
t736 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t776;
t639 = -t704 * pkin(4) - t723 * pkin(8) + t726 * t706 - t643;
t633 = t639 - 0.2e1 * qJD(6) * t689 + (t689 * t740 + t650) * pkin(5) + (t688 * t740 - t651) * qJ(6);
t626 = m(7) * t633 + t650 * mrSges(7,1) - t651 * mrSges(7,3) - t689 * t675 + t688 * t676;
t761 = -m(6) * t639 - t650 * mrSges(6,1) - t651 * mrSges(6,2) - t688 * t673 - t689 * t674 - t626;
t625 = m(5) * t643 + t704 * mrSges(5,1) - t705 * mrSges(5,3) + t725 * t700 - t726 * t703 + t761;
t758 = -m(4) * t671 - t704 * mrSges(4,1) - t705 * mrSges(4,2) - t725 * t701 - t726 * t702 - t625;
t621 = m(3) * t694 + qJDD(2) * mrSges(3,1) - t733 * mrSges(3,3) + qJD(2) * t736 - t732 * t777 + t758;
t768 = t754 * t612 - t621 * t752;
t606 = m(2) * t738 - mrSges(2,1) * t757 - qJDD(1) * mrSges(2,2) + t768;
t613 = t750 * t615 + t616 * t788;
t759 = -m(3) * t718 + t734 * mrSges(3,1) - t733 * mrSges(3,2) - t735 * t777 + t736 * t776 - t613;
t609 = m(2) * t737 + qJDD(1) * mrSges(2,1) - t757 * mrSges(2,2) + t759;
t786 = t753 * t606 + t755 * t609;
t607 = t752 * t612 + t754 * t621;
t785 = t800 * t688 - t793 * t689 - t789 * t740;
t784 = t789 * t688 - t792 * t689 - t798 * t740;
t783 = -t793 * t688 + t804 * t689 + t792 * t740;
t781 = t802 * t725 + t794 * t726 - t799 * t776;
t780 = t799 * t725 - t803 * t726 + t801 * t776;
t779 = t794 * t725 - t805 * t726 + t803 * t776;
t769 = t755 * t606 - t609 * t753;
t717 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t752 + Ifges(3,4) * t754) * qJD(1);
t716 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t752 + Ifges(3,2) * t754) * qJD(1);
t715 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t752 + Ifges(3,6) * t754) * qJD(1);
t618 = mrSges(6,2) * t639 + mrSges(7,2) * t629 - mrSges(6,3) * t630 - mrSges(7,3) * t633 - qJ(6) * t626 - t793 * t650 + t804 * t651 + t784 * t688 + t792 * t730 + t785 * t740;
t617 = -mrSges(6,1) * t639 - mrSges(7,1) * t633 + mrSges(7,2) * t628 + mrSges(6,3) * t631 - pkin(5) * t626 - t800 * t650 + t793 * t651 + t784 * t689 + t789 * t730 + t783 * t740;
t603 = mrSges(4,2) * t671 + mrSges(5,2) * t642 - mrSges(4,3) * t644 - mrSges(5,3) * t643 - pkin(8) * t619 - qJ(4) * t625 - t751 * t617 + t797 * t618 - t794 * t704 + t805 * t705 + t780 * t725 - t734 * t803 + t781 * t776;
t602 = -mrSges(4,1) * t671 - mrSges(5,1) * t643 + mrSges(5,2) * t641 + mrSges(4,3) * t645 - pkin(3) * t625 - pkin(4) * t761 - pkin(8) * t766 - t797 * t617 - t751 * t618 + t802 * t704 + t794 * t705 + t780 * t726 - t734 * t799 + t779 * t776;
t601 = pkin(5) * t765 - qJ(4) * (-t703 * t776 + t763) - pkin(3) * (-t700 * t776 + t762) + qJ(6) * t772 + (-mrSges(7,2) * pkin(5) + t792) * t651 + pkin(4) * t619 + qJD(2) * t717 - mrSges(3,1) * t718 - pkin(2) * t613 + (Ifges(3,2) + t801) * t734 + t798 * t730 + (qJ(4) * t691 + t779) * t725 + (pkin(3) * t691 - t781) * t726 + (-qJ(6) * t663 + t783) * t688 + (-pkin(5) * t663 - t785) * t689 + Ifges(3,6) * qJDD(2) + (-mrSges(7,2) * qJ(6) - t789) * t650 - mrSges(5,3) * t641 - t715 * t777 + mrSges(7,3) * t628 - mrSges(7,1) * t629 + mrSges(5,1) * t642 - mrSges(4,1) * t644 + mrSges(4,2) * t645 + Ifges(3,4) * t733 + mrSges(6,1) * t630 - mrSges(6,2) * t631 + mrSges(3,3) * t695 + (mrSges(5,2) * qJ(4) + t799) * t704 + (mrSges(5,2) * pkin(3) - t803) * t705;
t600 = mrSges(3,2) * t718 - mrSges(3,3) * t694 + Ifges(3,1) * t733 + Ifges(3,4) * t734 + Ifges(3,5) * qJDD(2) - qJ(3) * t613 - qJD(2) * t716 - t750 * t602 + t603 * t788 + t715 * t776;
t599 = Ifges(2,6) * qJDD(1) + t757 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t738 - Ifges(3,5) * t733 - Ifges(3,6) * t734 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t694 + mrSges(3,2) * t695 - t750 * t603 - t788 * t602 - pkin(2) * t758 - qJ(3) * t767 - pkin(1) * t607 + (-t716 * t752 + t717 * t754) * qJD(1);
t598 = -mrSges(2,2) * g(3) - mrSges(2,3) * t737 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t757 - pkin(7) * t607 + t600 * t754 - t601 * t752;
t1 = [-m(1) * g(1) + t769; -m(1) * g(2) + t786; (-m(1) - m(2)) * g(3) + t607; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t786 + t755 * t598 - t753 * t599; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t769 + t753 * t598 + t755 * t599; -mrSges(1,1) * g(2) + mrSges(2,1) * t737 + mrSges(1,2) * g(1) - mrSges(2,2) * t738 + Ifges(2,3) * qJDD(1) + pkin(1) * t759 + pkin(7) * t768 + t752 * t600 + t754 * t601;];
tauB  = t1;
