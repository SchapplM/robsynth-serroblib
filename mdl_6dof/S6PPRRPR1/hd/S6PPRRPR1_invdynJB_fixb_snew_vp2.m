% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PPRRPR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-05-04 20:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PPRRPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:04:31
% EndTime: 2019-05-04 20:04:55
% DurationCPUTime: 22.95s
% Computational Cost: add. (396020->293), mult. (750100->387), div. (0->0), fcn. (585404->16), ass. (0->135)
t751 = sin(pkin(11));
t756 = cos(pkin(11));
t742 = -g(1) * t756 - g(2) * t751;
t750 = sin(pkin(12));
t755 = cos(pkin(12));
t741 = g(1) * t751 - g(2) * t756;
t748 = -g(3) + qJDD(1);
t753 = sin(pkin(6));
t758 = cos(pkin(6));
t774 = t741 * t758 + t748 * t753;
t709 = -t750 * t742 + t755 * t774;
t710 = t755 * t742 + t750 * t774;
t720 = -t741 * t753 + t748 * t758 + qJDD(2);
t764 = cos(qJ(3));
t757 = cos(pkin(7));
t761 = sin(qJ(3));
t786 = t757 * t761;
t752 = sin(pkin(7));
t787 = t752 * t761;
t693 = t709 * t786 + t764 * t710 + t720 * t787;
t766 = qJD(3) ^ 2;
t691 = -pkin(3) * t766 + qJDD(3) * pkin(9) + t693;
t703 = -t709 * t752 + t720 * t757;
t760 = sin(qJ(4));
t763 = cos(qJ(4));
t684 = t763 * t691 + t760 * t703;
t737 = (-pkin(4) * t763 - qJ(5) * t760) * qJD(3);
t765 = qJD(4) ^ 2;
t783 = qJD(3) * t763;
t682 = -pkin(4) * t765 + qJDD(4) * qJ(5) + t737 * t783 + t684;
t692 = -t761 * t710 + (t709 * t757 + t720 * t752) * t764;
t690 = -qJDD(3) * pkin(3) - t766 * pkin(9) - t692;
t782 = qJD(3) * qJD(4);
t781 = t763 * t782;
t739 = qJDD(3) * t760 + t781;
t746 = t760 * t782;
t740 = qJDD(3) * t763 - t746;
t687 = (-t739 - t781) * qJ(5) + (-t740 + t746) * pkin(4) + t690;
t749 = sin(pkin(13));
t754 = cos(pkin(13));
t784 = qJD(3) * t760;
t734 = qJD(4) * t749 + t754 * t784;
t677 = -0.2e1 * qJD(5) * t734 - t749 * t682 + t754 * t687;
t723 = qJDD(4) * t749 + t739 * t754;
t733 = qJD(4) * t754 - t749 * t784;
t675 = (-t733 * t783 - t723) * pkin(10) + (t733 * t734 - t740) * pkin(5) + t677;
t678 = 0.2e1 * qJD(5) * t733 + t754 * t682 + t749 * t687;
t722 = qJDD(4) * t754 - t739 * t749;
t724 = -pkin(5) * t783 - pkin(10) * t734;
t732 = t733 ^ 2;
t676 = -pkin(5) * t732 + pkin(10) * t722 + t724 * t783 + t678;
t759 = sin(qJ(6));
t762 = cos(qJ(6));
t672 = t675 * t762 - t676 * t759;
t715 = t733 * t762 - t734 * t759;
t696 = qJD(6) * t715 + t722 * t759 + t723 * t762;
t716 = t733 * t759 + t734 * t762;
t702 = -mrSges(7,1) * t715 + mrSges(7,2) * t716;
t745 = qJD(6) - t783;
t705 = -mrSges(7,2) * t745 + mrSges(7,3) * t715;
t736 = qJDD(6) - t740;
t668 = m(7) * t672 + mrSges(7,1) * t736 - mrSges(7,3) * t696 - t702 * t716 + t705 * t745;
t673 = t675 * t759 + t676 * t762;
t695 = -qJD(6) * t716 + t722 * t762 - t723 * t759;
t706 = mrSges(7,1) * t745 - mrSges(7,3) * t716;
t669 = m(7) * t673 - mrSges(7,2) * t736 + mrSges(7,3) * t695 + t702 * t715 - t706 * t745;
t662 = t762 * t668 + t759 * t669;
t717 = -mrSges(6,1) * t733 + mrSges(6,2) * t734;
t773 = mrSges(6,2) * t783 + mrSges(6,3) * t733;
t660 = m(6) * t677 - t740 * mrSges(6,1) - t723 * mrSges(6,3) - t734 * t717 - t773 * t783 + t662;
t721 = -mrSges(6,1) * t783 - mrSges(6,3) * t734;
t778 = -t668 * t759 + t762 * t669;
t661 = m(6) * t678 + mrSges(6,2) * t740 + mrSges(6,3) * t722 + t717 * t733 + t721 * t783 + t778;
t658 = -t660 * t749 + t754 * t661;
t738 = (-mrSges(5,1) * t763 + mrSges(5,2) * t760) * qJD(3);
t743 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t784;
t656 = m(5) * t684 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t740 - qJD(4) * t743 + t738 * t783 + t658;
t683 = -t760 * t691 + t703 * t763;
t681 = -qJDD(4) * pkin(4) - qJ(5) * t765 + t737 * t784 + qJDD(5) - t683;
t679 = -pkin(5) * t722 - pkin(10) * t732 + t724 * t734 + t681;
t770 = m(7) * t679 - t695 * mrSges(7,1) + mrSges(7,2) * t696 - t715 * t705 + t706 * t716;
t674 = m(6) * t681 - t722 * mrSges(6,1) + mrSges(6,2) * t723 + t721 * t734 - t733 * t773 + t770;
t744 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t783;
t671 = m(5) * t683 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t739 + qJD(4) * t744 - t738 * t784 - t674;
t779 = t763 * t656 - t671 * t760;
t645 = m(4) * t693 - mrSges(4,1) * t766 - qJDD(3) * mrSges(4,2) + t779;
t648 = t760 * t656 + t763 * t671;
t647 = m(4) * t703 + t648;
t657 = t660 * t754 + t661 * t749;
t769 = -m(5) * t690 + t740 * mrSges(5,1) - mrSges(5,2) * t739 - t743 * t784 + t744 * t783 - t657;
t653 = m(4) * t692 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t766 + t769;
t788 = t653 * t764;
t634 = t645 * t786 - t647 * t752 + t757 * t788;
t630 = m(3) * t709 + t634;
t640 = t764 * t645 - t653 * t761;
t639 = m(3) * t710 + t640;
t792 = t630 * t755 + t639 * t750;
t697 = Ifges(7,5) * t716 + Ifges(7,6) * t715 + Ifges(7,3) * t745;
t699 = Ifges(7,1) * t716 + Ifges(7,4) * t715 + Ifges(7,5) * t745;
t663 = -mrSges(7,1) * t679 + mrSges(7,3) * t673 + Ifges(7,4) * t696 + Ifges(7,2) * t695 + Ifges(7,6) * t736 - t697 * t716 + t699 * t745;
t698 = Ifges(7,4) * t716 + Ifges(7,2) * t715 + Ifges(7,6) * t745;
t664 = mrSges(7,2) * t679 - mrSges(7,3) * t672 + Ifges(7,1) * t696 + Ifges(7,4) * t695 + Ifges(7,5) * t736 + t697 * t715 - t698 * t745;
t711 = Ifges(6,5) * t734 + Ifges(6,6) * t733 - Ifges(6,3) * t783;
t713 = Ifges(6,1) * t734 + Ifges(6,4) * t733 - Ifges(6,5) * t783;
t649 = -mrSges(6,1) * t681 + mrSges(6,3) * t678 + Ifges(6,4) * t723 + Ifges(6,2) * t722 - Ifges(6,6) * t740 - pkin(5) * t770 + pkin(10) * t778 + t762 * t663 + t759 * t664 - t734 * t711 - t713 * t783;
t712 = Ifges(6,4) * t734 + Ifges(6,2) * t733 - Ifges(6,6) * t783;
t650 = mrSges(6,2) * t681 - mrSges(6,3) * t677 + Ifges(6,1) * t723 + Ifges(6,4) * t722 - Ifges(6,5) * t740 - pkin(10) * t662 - t663 * t759 + t664 * t762 + t711 * t733 + t712 * t783;
t729 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t760 + Ifges(5,2) * t763) * qJD(3);
t730 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t760 + Ifges(5,4) * t763) * qJD(3);
t791 = mrSges(5,1) * t683 - mrSges(5,2) * t684 + Ifges(5,5) * t739 + Ifges(5,6) * t740 + Ifges(5,3) * qJDD(4) - pkin(4) * t674 + qJ(5) * t658 + t754 * t649 + t749 * t650 + (t729 * t760 - t730 * t763) * qJD(3);
t633 = t645 * t787 + t757 * t647 + t752 * t788;
t632 = m(3) * t720 + t633;
t620 = -t632 * t753 + t758 * t792;
t618 = m(2) * t741 + t620;
t626 = -t630 * t750 + t755 * t639;
t625 = m(2) * t742 + t626;
t785 = t756 * t618 + t751 * t625;
t619 = t758 * t632 + t753 * t792;
t780 = -t618 * t751 + t756 * t625;
t777 = m(2) * t748 + t619;
t728 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t760 + Ifges(5,6) * t763) * qJD(3);
t635 = mrSges(5,2) * t690 - mrSges(5,3) * t683 + Ifges(5,1) * t739 + Ifges(5,4) * t740 + Ifges(5,5) * qJDD(4) - qJ(5) * t657 - qJD(4) * t729 - t649 * t749 + t650 * t754 + t728 * t783;
t768 = mrSges(7,1) * t672 - mrSges(7,2) * t673 + Ifges(7,5) * t696 + Ifges(7,6) * t695 + Ifges(7,3) * t736 + t716 * t698 - t715 * t699;
t641 = Ifges(5,6) * qJDD(4) - t728 * t784 - pkin(4) * t657 - t768 + (Ifges(5,2) + Ifges(6,3)) * t740 + t733 * t713 - t734 * t712 + Ifges(5,4) * t739 + qJD(4) * t730 - Ifges(6,6) * t722 - Ifges(6,5) * t723 - pkin(5) * t662 - mrSges(6,1) * t677 + mrSges(6,2) * t678 + mrSges(5,3) * t684 - mrSges(5,1) * t690;
t622 = mrSges(4,2) * t703 - mrSges(4,3) * t692 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t766 - pkin(9) * t648 + t635 * t763 - t641 * t760;
t627 = -mrSges(4,1) * t703 + mrSges(4,3) * t693 + t766 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t648 - t791;
t772 = pkin(8) * t640 + t622 * t761 + t627 * t764;
t621 = mrSges(4,1) * t692 - mrSges(4,2) * t693 + Ifges(4,3) * qJDD(3) + pkin(3) * t769 + pkin(9) * t779 + t760 * t635 + t763 * t641;
t615 = -mrSges(3,1) * t720 + mrSges(3,3) * t710 - pkin(2) * t633 - t752 * t621 + t757 * t772;
t616 = mrSges(3,2) * t720 - mrSges(3,3) * t709 + t764 * t622 - t761 * t627 + (-t633 * t752 - t634 * t757) * pkin(8);
t771 = qJ(2) * t626 + t615 * t755 + t616 * t750;
t614 = mrSges(3,1) * t709 - mrSges(3,2) * t710 + pkin(2) * t634 + t757 * t621 + t752 * t772;
t613 = mrSges(2,2) * t748 - mrSges(2,3) * t741 - t750 * t615 + t755 * t616 + (-t619 * t753 - t620 * t758) * qJ(2);
t612 = -mrSges(2,1) * t748 + mrSges(2,3) * t742 - pkin(1) * t619 - t753 * t614 + t758 * t771;
t1 = [-m(1) * g(1) + t780; -m(1) * g(2) + t785; -m(1) * g(3) + t777; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t785 - t751 * t612 + t756 * t613; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t780 + t756 * t612 + t751 * t613; -mrSges(1,1) * g(2) + mrSges(2,1) * t741 + mrSges(1,2) * g(1) - mrSges(2,2) * t742 + pkin(1) * t620 + t758 * t614 + t753 * t771; t777; t632; t621; t791; t674; t768;];
tauJB  = t1;
