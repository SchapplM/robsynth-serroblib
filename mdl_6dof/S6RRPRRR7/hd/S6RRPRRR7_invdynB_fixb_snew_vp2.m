% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-05-06 22:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 22:17:07
% EndTime: 2019-05-06 22:17:23
% DurationCPUTime: 8.60s
% Computational Cost: add. (128764->364), mult. (261317->445), div. (0->0), fcn. (170767->10), ass. (0->142)
t800 = Ifges(3,1) + Ifges(4,1);
t795 = Ifges(3,4) - Ifges(4,5);
t794 = Ifges(3,5) + Ifges(4,4);
t799 = Ifges(3,2) + Ifges(4,3);
t793 = Ifges(3,6) - Ifges(4,6);
t798 = Ifges(3,3) + Ifges(4,2);
t797 = 2 * qJD(3);
t796 = mrSges(3,3) + mrSges(4,2);
t767 = cos(qJ(2));
t770 = qJD(1) ^ 2;
t792 = t767 ^ 2 * t770;
t763 = sin(qJ(1));
t768 = cos(qJ(1));
t740 = -t768 * g(1) - t763 * g(2);
t717 = -t770 * pkin(1) + qJDD(1) * pkin(7) + t740;
t762 = sin(qJ(2));
t697 = -t762 * g(3) + t767 * t717;
t729 = (-mrSges(3,1) * t767 + mrSges(3,2) * t762) * qJD(1);
t785 = qJD(1) * qJD(2);
t784 = t762 * t785;
t731 = t767 * qJDD(1) - t784;
t787 = qJD(1) * t762;
t734 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t787;
t727 = (-pkin(2) * t767 - qJ(3) * t762) * qJD(1);
t769 = qJD(2) ^ 2;
t786 = qJD(1) * t767;
t676 = -t769 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t797 + t727 * t786 + t697;
t728 = (-mrSges(4,1) * t767 - mrSges(4,3) * t762) * qJD(1);
t735 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t787;
t738 = -qJD(2) * pkin(3) - pkin(8) * t787;
t663 = -pkin(3) * t792 - t731 * pkin(8) + qJD(2) * t738 + t676;
t696 = -t767 * g(3) - t762 * t717;
t680 = -qJDD(2) * pkin(2) - t769 * qJ(3) + t727 * t787 + qJDD(3) - t696;
t783 = t767 * t785;
t730 = t762 * qJDD(1) + t783;
t664 = (-t730 + t783) * pkin(8) + (-t762 * t767 * t770 - qJDD(2)) * pkin(3) + t680;
t761 = sin(qJ(4));
t766 = cos(qJ(4));
t647 = t766 * t663 + t761 * t664;
t715 = (-t761 * t767 + t762 * t766) * qJD(1);
t681 = -t715 * qJD(4) - t761 * t730 - t766 * t731;
t714 = -t761 * t787 - t766 * t786;
t691 = -t714 * mrSges(5,1) + t715 * mrSges(5,2);
t752 = -qJD(2) + qJD(4);
t699 = t752 * mrSges(5,1) - t715 * mrSges(5,3);
t751 = -qJDD(2) + qJDD(4);
t739 = t763 * g(1) - t768 * g(2);
t716 = -qJDD(1) * pkin(1) - t770 * pkin(7) - t739;
t776 = -t731 * pkin(2) + t716 + (-t730 - t783) * qJ(3);
t654 = -pkin(2) * t784 + t731 * pkin(3) - pkin(8) * t792 - t776 + (t738 + t797) * t787;
t682 = t714 * qJD(4) + t766 * t730 - t761 * t731;
t642 = t654 + (-t714 * t752 - t682) * pkin(9) + (t715 * t752 - t681) * pkin(4);
t692 = -t714 * pkin(4) - t715 * pkin(9);
t750 = t752 ^ 2;
t645 = -t750 * pkin(4) + t751 * pkin(9) + t714 * t692 + t647;
t760 = sin(qJ(5));
t765 = cos(qJ(5));
t634 = t765 * t642 - t760 * t645;
t694 = -t760 * t715 + t765 * t752;
t657 = t694 * qJD(5) + t765 * t682 + t760 * t751;
t679 = qJDD(5) - t681;
t695 = t765 * t715 + t760 * t752;
t707 = qJD(5) - t714;
t632 = (t694 * t707 - t657) * pkin(10) + (t694 * t695 + t679) * pkin(5) + t634;
t635 = t760 * t642 + t765 * t645;
t656 = -t695 * qJD(5) - t760 * t682 + t765 * t751;
t685 = t707 * pkin(5) - t695 * pkin(10);
t693 = t694 ^ 2;
t633 = -t693 * pkin(5) + t656 * pkin(10) - t707 * t685 + t635;
t759 = sin(qJ(6));
t764 = cos(qJ(6));
t630 = t764 * t632 - t759 * t633;
t670 = t764 * t694 - t759 * t695;
t639 = t670 * qJD(6) + t759 * t656 + t764 * t657;
t671 = t759 * t694 + t764 * t695;
t652 = -t670 * mrSges(7,1) + t671 * mrSges(7,2);
t702 = qJD(6) + t707;
t661 = -t702 * mrSges(7,2) + t670 * mrSges(7,3);
t674 = qJDD(6) + t679;
t628 = m(7) * t630 + t674 * mrSges(7,1) - t639 * mrSges(7,3) - t671 * t652 + t702 * t661;
t631 = t759 * t632 + t764 * t633;
t638 = -t671 * qJD(6) + t764 * t656 - t759 * t657;
t662 = t702 * mrSges(7,1) - t671 * mrSges(7,3);
t629 = m(7) * t631 - t674 * mrSges(7,2) + t638 * mrSges(7,3) + t670 * t652 - t702 * t662;
t621 = t764 * t628 + t759 * t629;
t672 = -t694 * mrSges(6,1) + t695 * mrSges(6,2);
t683 = -t707 * mrSges(6,2) + t694 * mrSges(6,3);
t619 = m(6) * t634 + t679 * mrSges(6,1) - t657 * mrSges(6,3) - t695 * t672 + t707 * t683 + t621;
t684 = t707 * mrSges(6,1) - t695 * mrSges(6,3);
t778 = -t759 * t628 + t764 * t629;
t620 = m(6) * t635 - t679 * mrSges(6,2) + t656 * mrSges(6,3) + t694 * t672 - t707 * t684 + t778;
t779 = -t760 * t619 + t765 * t620;
t614 = m(5) * t647 - t751 * mrSges(5,2) + t681 * mrSges(5,3) + t714 * t691 - t752 * t699 + t779;
t646 = -t761 * t663 + t766 * t664;
t698 = -t752 * mrSges(5,2) + t714 * mrSges(5,3);
t644 = -t751 * pkin(4) - t750 * pkin(9) + t715 * t692 - t646;
t636 = -t656 * pkin(5) - t693 * pkin(10) + t695 * t685 + t644;
t774 = m(7) * t636 - t638 * mrSges(7,1) + t639 * mrSges(7,2) - t670 * t661 + t671 * t662;
t772 = -m(6) * t644 + t656 * mrSges(6,1) - t657 * mrSges(6,2) + t694 * t683 - t695 * t684 - t774;
t624 = m(5) * t646 + t751 * mrSges(5,1) - t682 * mrSges(5,3) - t715 * t691 + t752 * t698 + t772;
t780 = t766 * t614 - t761 * t624;
t775 = m(4) * t676 + qJDD(2) * mrSges(4,3) + qJD(2) * t735 + t728 * t786 + t780;
t605 = m(3) * t697 - qJDD(2) * mrSges(3,2) - qJD(2) * t734 + t729 * t786 + t796 * t731 + t775;
t736 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t786;
t607 = t761 * t614 + t766 * t624;
t737 = mrSges(4,2) * t786 + qJD(2) * mrSges(4,3);
t773 = -m(4) * t680 + qJDD(2) * mrSges(4,1) + qJD(2) * t737 - t607;
t606 = m(3) * t696 + qJDD(2) * mrSges(3,1) + qJD(2) * t736 - t796 * t730 + (-t728 - t729) * t787 + t773;
t781 = t767 * t605 - t762 * t606;
t599 = m(2) * t740 - t770 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t781;
t669 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t787 + t776;
t615 = t765 * t619 + t760 * t620;
t777 = -m(5) * t654 + t681 * mrSges(5,1) - t682 * mrSges(5,2) + t714 * t698 - t715 * t699 - t615;
t612 = m(4) * t669 - t731 * mrSges(4,1) - t730 * mrSges(4,3) - t735 * t787 - t737 * t786 + t777;
t771 = -m(3) * t716 + t731 * mrSges(3,1) - t730 * mrSges(3,2) - t734 * t787 + t736 * t786 - t612;
t611 = m(2) * t739 + qJDD(1) * mrSges(2,1) - t770 * mrSges(2,2) + t771;
t791 = t763 * t599 + t768 * t611;
t600 = t762 * t605 + t767 * t606;
t790 = t798 * qJD(2) + (t794 * t762 + t793 * t767) * qJD(1);
t789 = -t793 * qJD(2) + (-t795 * t762 - t799 * t767) * qJD(1);
t788 = t794 * qJD(2) + (t800 * t762 + t795 * t767) * qJD(1);
t782 = t768 * t599 - t611 * t763;
t688 = Ifges(5,1) * t715 + Ifges(5,4) * t714 + Ifges(5,5) * t752;
t687 = Ifges(5,4) * t715 + Ifges(5,2) * t714 + Ifges(5,6) * t752;
t686 = Ifges(5,5) * t715 + Ifges(5,6) * t714 + Ifges(5,3) * t752;
t667 = Ifges(6,1) * t695 + Ifges(6,4) * t694 + Ifges(6,5) * t707;
t666 = Ifges(6,4) * t695 + Ifges(6,2) * t694 + Ifges(6,6) * t707;
t665 = Ifges(6,5) * t695 + Ifges(6,6) * t694 + Ifges(6,3) * t707;
t650 = Ifges(7,1) * t671 + Ifges(7,4) * t670 + Ifges(7,5) * t702;
t649 = Ifges(7,4) * t671 + Ifges(7,2) * t670 + Ifges(7,6) * t702;
t648 = Ifges(7,5) * t671 + Ifges(7,6) * t670 + Ifges(7,3) * t702;
t623 = mrSges(7,2) * t636 - mrSges(7,3) * t630 + Ifges(7,1) * t639 + Ifges(7,4) * t638 + Ifges(7,5) * t674 + t670 * t648 - t702 * t649;
t622 = -mrSges(7,1) * t636 + mrSges(7,3) * t631 + Ifges(7,4) * t639 + Ifges(7,2) * t638 + Ifges(7,6) * t674 - t671 * t648 + t702 * t650;
t609 = mrSges(6,2) * t644 - mrSges(6,3) * t634 + Ifges(6,1) * t657 + Ifges(6,4) * t656 + Ifges(6,5) * t679 - pkin(10) * t621 - t759 * t622 + t764 * t623 + t694 * t665 - t707 * t666;
t608 = -mrSges(6,1) * t644 + mrSges(6,3) * t635 + Ifges(6,4) * t657 + Ifges(6,2) * t656 + Ifges(6,6) * t679 - pkin(5) * t774 + pkin(10) * t778 + t764 * t622 + t759 * t623 - t695 * t665 + t707 * t667;
t601 = Ifges(5,4) * t682 + Ifges(5,2) * t681 + Ifges(5,6) * t751 - t715 * t686 + t752 * t688 - mrSges(5,1) * t654 + mrSges(5,3) * t647 - Ifges(6,5) * t657 - Ifges(6,6) * t656 - Ifges(6,3) * t679 - t695 * t666 + t694 * t667 - mrSges(6,1) * t634 + mrSges(6,2) * t635 - Ifges(7,5) * t639 - Ifges(7,6) * t638 - Ifges(7,3) * t674 - t671 * t649 + t670 * t650 - mrSges(7,1) * t630 + mrSges(7,2) * t631 - pkin(5) * t621 - pkin(4) * t615;
t596 = mrSges(5,2) * t654 - mrSges(5,3) * t646 + Ifges(5,1) * t682 + Ifges(5,4) * t681 + Ifges(5,5) * t751 - pkin(9) * t615 - t760 * t608 + t765 * t609 + t714 * t686 - t752 * t687;
t595 = mrSges(3,2) * t716 + mrSges(4,2) * t680 - mrSges(3,3) * t696 - mrSges(4,3) * t669 - pkin(8) * t607 - qJ(3) * t612 + t789 * qJD(2) + t794 * qJDD(2) + t766 * t596 - t761 * t601 + t800 * t730 + t795 * t731 + t790 * t786;
t594 = -mrSges(3,1) * t716 - mrSges(4,1) * t669 + mrSges(4,2) * t676 + mrSges(3,3) * t697 - pkin(2) * t612 - pkin(3) * t777 - pkin(8) * t780 + t788 * qJD(2) + t793 * qJDD(2) - t761 * t596 - t766 * t601 + t795 * t730 + t799 * t731 - t790 * t787;
t593 = (-mrSges(4,2) * qJ(3) - t793) * t731 + (mrSges(4,2) * pkin(2) - t794) * t730 + (t788 * t767 + (pkin(2) * t728 + t789) * t762) * qJD(1) - t798 * qJDD(2) + pkin(4) * t772 - pkin(2) * t773 - qJ(3) * t775 + mrSges(2,1) * g(3) - pkin(1) * t600 + Ifges(5,3) * t751 + mrSges(4,1) * t680 + Ifges(5,6) * t681 + Ifges(5,5) * t682 + Ifges(2,6) * qJDD(1) - t714 * t688 + t715 * t687 + mrSges(5,1) * t646 - mrSges(5,2) * t647 + pkin(3) * t607 + mrSges(2,3) * t740 - mrSges(4,3) * t676 + pkin(9) * t779 - mrSges(3,1) * t696 + mrSges(3,2) * t697 + t760 * t609 + t765 * t608 + t770 * Ifges(2,5);
t592 = -mrSges(2,2) * g(3) - mrSges(2,3) * t739 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t770 - pkin(7) * t600 - t594 * t762 + t595 * t767;
t1 = [-m(1) * g(1) + t782; -m(1) * g(2) + t791; (-m(1) - m(2)) * g(3) + t600; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t791 + t768 * t592 - t763 * t593; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t782 + t763 * t592 + t768 * t593; -mrSges(1,1) * g(2) + mrSges(2,1) * t739 + mrSges(1,2) * g(1) - mrSges(2,2) * t740 + Ifges(2,3) * qJDD(1) + pkin(1) * t771 + pkin(7) * t781 + t767 * t594 + t762 * t595;];
tauB  = t1;
