% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-05-05 18:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:52:34
% EndTime: 2019-05-05 18:52:42
% DurationCPUTime: 7.92s
% Computational Cost: add. (106520->341), mult. (257858->416), div. (0->0), fcn. (190738->10), ass. (0->142)
t798 = Ifges(4,1) + Ifges(5,1);
t791 = Ifges(4,4) - Ifges(5,5);
t790 = Ifges(4,5) + Ifges(5,4);
t797 = Ifges(4,2) + Ifges(5,3);
t789 = Ifges(4,6) - Ifges(5,6);
t796 = -Ifges(4,3) - Ifges(5,2);
t745 = sin(pkin(10));
t755 = qJD(1) ^ 2;
t746 = cos(pkin(10));
t786 = t746 ^ 2 * t755;
t795 = t745 ^ 2 * t755 + t786;
t794 = 2 * qJD(4);
t793 = cos(qJ(3));
t792 = -mrSges(4,3) - mrSges(5,2);
t788 = mrSges(3,2) * t745;
t750 = sin(qJ(1));
t753 = cos(qJ(1));
t721 = -g(1) * t753 - g(2) * t750;
t717 = -pkin(1) * t755 + qJDD(1) * qJ(2) + t721;
t777 = qJD(1) * qJD(2);
t773 = -g(3) * t746 - 0.2e1 * t745 * t777;
t675 = (pkin(2) * t746 * t755 - pkin(7) * qJDD(1) - t717) * t745 + t773;
t702 = -g(3) * t745 + (t717 + 0.2e1 * t777) * t746;
t775 = qJDD(1) * t746;
t676 = -pkin(2) * t786 + pkin(7) * t775 + t702;
t749 = sin(qJ(3));
t656 = t749 * t675 + t793 * t676;
t774 = t746 * t793;
t776 = qJDD(1) * t745;
t763 = t793 * t745 + t746 * t749;
t716 = t763 * qJD(1);
t779 = qJD(3) * t716;
t699 = -qJDD(1) * t774 + t749 * t776 + t779;
t706 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t716;
t780 = qJD(1) * t745;
t715 = -qJD(1) * t774 + t749 * t780;
t691 = pkin(3) * t715 - qJ(4) * t716;
t754 = qJD(3) ^ 2;
t643 = -pkin(3) * t754 + qJDD(3) * qJ(4) + qJD(3) * t794 - t715 * t691 + t656;
t707 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t716;
t655 = t793 * t675 - t749 * t676;
t644 = -qJDD(3) * pkin(3) - t754 * qJ(4) + t716 * t691 + qJDD(4) - t655;
t778 = t715 * qJD(3);
t700 = t763 * qJDD(1) - t778;
t634 = (-t700 - t778) * pkin(8) + (t715 * t716 - qJDD(3)) * pkin(4) + t644;
t709 = -qJD(3) * pkin(4) - pkin(8) * t716;
t714 = t715 ^ 2;
t636 = -pkin(4) * t714 + pkin(8) * t699 + qJD(3) * t709 + t643;
t748 = sin(qJ(5));
t752 = cos(qJ(5));
t630 = t748 * t634 + t752 * t636;
t685 = t715 * t748 + t716 * t752;
t652 = -qJD(5) * t685 + t699 * t752 - t700 * t748;
t684 = t715 * t752 - t716 * t748;
t664 = -mrSges(6,1) * t684 + mrSges(6,2) * t685;
t739 = -qJD(3) + qJD(5);
t672 = mrSges(6,1) * t739 - mrSges(6,3) * t685;
t736 = -qJDD(3) + qJDD(5);
t665 = -pkin(5) * t684 - pkin(9) * t685;
t735 = t739 ^ 2;
t627 = -pkin(5) * t735 + pkin(9) * t736 + t665 * t684 + t630;
t720 = t750 * g(1) - t753 * g(2);
t713 = -qJDD(1) * pkin(1) - t755 * qJ(2) + qJDD(2) - t720;
t698 = -pkin(2) * t775 - t795 * pkin(7) + t713;
t758 = t699 * pkin(3) + t698 + (-t700 + t778) * qJ(4);
t632 = -pkin(3) * t779 - pkin(4) * t699 - pkin(8) * t714 - t758 + (t709 + t794) * t716;
t653 = qJD(5) * t684 + t699 * t748 + t700 * t752;
t628 = t632 + (t685 * t739 - t652) * pkin(5) + (-t684 * t739 - t653) * pkin(9);
t747 = sin(qJ(6));
t751 = cos(qJ(6));
t624 = -t627 * t747 + t628 * t751;
t668 = -t685 * t747 + t739 * t751;
t639 = qJD(6) * t668 + t653 * t751 + t736 * t747;
t651 = qJDD(6) - t652;
t669 = t685 * t751 + t739 * t747;
t654 = -mrSges(7,1) * t668 + mrSges(7,2) * t669;
t677 = qJD(6) - t684;
t657 = -mrSges(7,2) * t677 + mrSges(7,3) * t668;
t622 = m(7) * t624 + mrSges(7,1) * t651 - mrSges(7,3) * t639 - t654 * t669 + t657 * t677;
t625 = t627 * t751 + t628 * t747;
t638 = -qJD(6) * t669 - t653 * t747 + t736 * t751;
t658 = mrSges(7,1) * t677 - mrSges(7,3) * t669;
t623 = m(7) * t625 - mrSges(7,2) * t651 + mrSges(7,3) * t638 + t654 * t668 - t658 * t677;
t768 = -t622 * t747 + t751 * t623;
t614 = m(6) * t630 - mrSges(6,2) * t736 + mrSges(6,3) * t652 + t664 * t684 - t672 * t739 + t768;
t629 = t634 * t752 - t636 * t748;
t671 = -mrSges(6,2) * t739 + mrSges(6,3) * t684;
t626 = -pkin(5) * t736 - pkin(9) * t735 + t665 * t685 - t629;
t760 = -m(7) * t626 + t638 * mrSges(7,1) - mrSges(7,2) * t639 + t668 * t657 - t658 * t669;
t618 = m(6) * t629 + mrSges(6,1) * t736 - mrSges(6,3) * t653 - t664 * t685 + t671 * t739 + t760;
t769 = t752 * t614 - t618 * t748;
t762 = m(5) * t643 + qJDD(3) * mrSges(5,3) + qJD(3) * t707 + t769;
t692 = mrSges(5,1) * t715 - mrSges(5,3) * t716;
t781 = -mrSges(4,1) * t715 - mrSges(4,2) * t716 - t692;
t607 = m(4) * t656 - qJDD(3) * mrSges(4,2) - qJD(3) * t706 + t792 * t699 + t781 * t715 + t762;
t705 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t715;
t609 = t614 * t748 + t618 * t752;
t708 = -mrSges(5,2) * t715 + qJD(3) * mrSges(5,3);
t759 = -m(5) * t644 + qJDD(3) * mrSges(5,1) + qJD(3) * t708 - t609;
t608 = m(4) * t655 + qJDD(3) * mrSges(4,1) + qJD(3) * t705 + t792 * t700 + t781 * t716 + t759;
t601 = t749 * t607 + t793 * t608;
t701 = -t717 * t745 + t773;
t764 = mrSges(3,3) * qJDD(1) + t755 * (-mrSges(3,1) * t746 + t788);
t599 = m(3) * t701 - t764 * t745 + t601;
t770 = t793 * t607 - t608 * t749;
t600 = m(3) * t702 + t764 * t746 + t770;
t771 = -t599 * t745 + t746 * t600;
t594 = m(2) * t721 - mrSges(2,1) * t755 - qJDD(1) * mrSges(2,2) + t771;
t641 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t716 + t758;
t615 = t751 * t622 + t747 * t623;
t761 = -m(6) * t632 + t652 * mrSges(6,1) - t653 * mrSges(6,2) + t684 * t671 - t685 * t672 - t615;
t612 = m(5) * t641 + t699 * mrSges(5,1) - t700 * mrSges(5,3) - t716 * t707 + t715 * t708 + t761;
t757 = m(4) * t698 + t699 * mrSges(4,1) + t700 * mrSges(4,2) + t715 * t705 + t716 * t706 + t612;
t756 = -m(3) * t713 + mrSges(3,1) * t775 + t795 * mrSges(3,3) - t757;
t611 = t756 + (mrSges(2,1) - t788) * qJDD(1) - t755 * mrSges(2,2) + m(2) * t720;
t785 = t750 * t594 + t753 * t611;
t595 = t746 * t599 + t745 * t600;
t784 = -t789 * qJD(3) + t797 * t715 - t791 * t716;
t783 = t796 * qJD(3) + t789 * t715 - t790 * t716;
t782 = t790 * qJD(3) - t791 * t715 + t798 * t716;
t772 = t753 * t594 - t611 * t750;
t767 = Ifges(3,1) * t745 + Ifges(3,4) * t746;
t766 = Ifges(3,4) * t745 + Ifges(3,2) * t746;
t765 = Ifges(3,5) * t745 + Ifges(3,6) * t746;
t719 = t765 * qJD(1);
t661 = Ifges(6,1) * t685 + Ifges(6,4) * t684 + Ifges(6,5) * t739;
t660 = Ifges(6,4) * t685 + Ifges(6,2) * t684 + Ifges(6,6) * t739;
t659 = Ifges(6,5) * t685 + Ifges(6,6) * t684 + Ifges(6,3) * t739;
t647 = Ifges(7,1) * t669 + Ifges(7,4) * t668 + Ifges(7,5) * t677;
t646 = Ifges(7,4) * t669 + Ifges(7,2) * t668 + Ifges(7,6) * t677;
t645 = Ifges(7,5) * t669 + Ifges(7,6) * t668 + Ifges(7,3) * t677;
t617 = mrSges(7,2) * t626 - mrSges(7,3) * t624 + Ifges(7,1) * t639 + Ifges(7,4) * t638 + Ifges(7,5) * t651 + t645 * t668 - t646 * t677;
t616 = -mrSges(7,1) * t626 + mrSges(7,3) * t625 + Ifges(7,4) * t639 + Ifges(7,2) * t638 + Ifges(7,6) * t651 - t645 * t669 + t647 * t677;
t603 = -mrSges(6,1) * t632 - mrSges(7,1) * t624 + mrSges(7,2) * t625 + mrSges(6,3) * t630 + Ifges(6,4) * t653 - Ifges(7,5) * t639 + Ifges(6,2) * t652 + Ifges(6,6) * t736 - Ifges(7,6) * t638 - Ifges(7,3) * t651 - pkin(5) * t615 - t646 * t669 + t647 * t668 - t659 * t685 + t661 * t739;
t602 = mrSges(6,2) * t632 - mrSges(6,3) * t629 + Ifges(6,1) * t653 + Ifges(6,4) * t652 + Ifges(6,5) * t736 - pkin(9) * t615 - t616 * t747 + t617 * t751 + t659 * t684 - t660 * t739;
t591 = mrSges(4,2) * t698 + mrSges(5,2) * t644 - mrSges(4,3) * t655 - mrSges(5,3) * t641 - pkin(8) * t609 - qJ(4) * t612 + t784 * qJD(3) + t790 * qJDD(3) + t602 * t752 - t603 * t748 - t791 * t699 + t798 * t700 + t783 * t715;
t590 = -mrSges(4,1) * t698 - mrSges(5,1) * t641 + mrSges(5,2) * t643 + mrSges(4,3) * t656 - pkin(3) * t612 - pkin(4) * t761 - pkin(8) * t769 + t782 * qJD(3) + t789 * qJDD(3) - t748 * t602 - t752 * t603 - t797 * t699 + t791 * t700 + t783 * t716;
t589 = (Ifges(2,6) - t765) * qJDD(1) - qJ(4) * t762 - pkin(3) * t759 + pkin(5) * t760 + (qJ(4) * t692 - t782) * t715 + (pkin(3) * t692 + t784) * t716 + t796 * qJDD(3) + t751 * t616 + t747 * t617 + Ifges(6,3) * t736 + mrSges(2,3) * t721 - mrSges(3,1) * t701 + mrSges(3,2) * t702 - t684 * t661 + t685 * t660 + mrSges(4,2) * t656 + Ifges(6,6) * t652 + Ifges(6,5) * t653 + mrSges(2,1) * g(3) - mrSges(4,1) * t655 - mrSges(5,3) * t643 + mrSges(5,1) * t644 + mrSges(6,1) * t629 - mrSges(6,2) * t630 + pkin(4) * t609 - pkin(2) * t601 + pkin(9) * t768 + (qJ(4) * mrSges(5,2) + t789) * t699 + (pkin(3) * mrSges(5,2) - t790) * t700 - pkin(1) * t595 + (-t745 * t766 + t746 * t767 + Ifges(2,5)) * t755;
t588 = t746 * qJD(1) * t719 + mrSges(3,2) * t713 - mrSges(3,3) * t701 - pkin(7) * t601 + t767 * qJDD(1) - t749 * t590 + t793 * t591;
t587 = -mrSges(3,1) * t713 + mrSges(3,3) * t702 - pkin(2) * t757 + pkin(7) * t770 + t766 * qJDD(1) + t793 * t590 + t749 * t591 - t719 * t780;
t586 = -mrSges(2,2) * g(3) - mrSges(2,3) * t720 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t755 - qJ(2) * t595 - t587 * t745 + t588 * t746;
t1 = [-m(1) * g(1) + t772; -m(1) * g(2) + t785; (-m(1) - m(2)) * g(3) + t595; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t785 + t753 * t586 - t750 * t589; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t772 + t750 * t586 + t753 * t589; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t720 - mrSges(2,2) * t721 + t745 * t588 + t746 * t587 + pkin(1) * (-mrSges(3,2) * t776 + t756) + qJ(2) * t771;];
tauB  = t1;
