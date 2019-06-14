% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 10:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:40:24
% EndTime: 2019-05-07 10:40:45
% DurationCPUTime: 10.14s
% Computational Cost: add. (149083->367), mult. (306446->444), div. (0->0), fcn. (214210->10), ass. (0->143)
t791 = -2 * qJD(4);
t790 = Ifges(4,1) + Ifges(5,2);
t789 = -Ifges(5,1) - Ifges(4,3);
t785 = Ifges(4,4) + Ifges(5,6);
t784 = Ifges(4,5) - Ifges(5,4);
t788 = Ifges(4,2) + Ifges(5,3);
t783 = Ifges(4,6) - Ifges(5,5);
t752 = sin(qJ(2));
t756 = cos(qJ(2));
t772 = qJD(1) * qJD(2);
t733 = qJDD(1) * t752 + t756 * t772;
t753 = sin(qJ(1));
t757 = cos(qJ(1));
t739 = -g(1) * t757 - g(2) * t753;
t758 = qJD(1) ^ 2;
t727 = -pkin(1) * t758 + qJDD(1) * pkin(7) + t739;
t781 = t727 * t752;
t786 = pkin(2) * t758;
t674 = qJDD(2) * pkin(2) - pkin(8) * t733 - t781 + (pkin(8) * t772 + t752 * t786 - g(3)) * t756;
t709 = -g(3) * t752 + t756 * t727;
t734 = qJDD(1) * t756 - t752 * t772;
t775 = qJD(1) * t752;
t737 = qJD(2) * pkin(2) - pkin(8) * t775;
t748 = t756 ^ 2;
t675 = pkin(8) * t734 - qJD(2) * t737 - t748 * t786 + t709;
t751 = sin(qJ(3));
t787 = cos(qJ(3));
t658 = t751 * t674 + t787 * t675;
t774 = qJD(1) * t756;
t724 = t751 * t775 - t774 * t787;
t725 = (t751 * t756 + t752 * t787) * qJD(1);
t700 = pkin(3) * t724 - qJ(4) * t725;
t747 = qJD(2) + qJD(3);
t745 = t747 ^ 2;
t746 = qJDD(2) + qJDD(3);
t650 = pkin(3) * t745 - t746 * qJ(4) + t724 * t700 + t747 * t791 - t658;
t782 = t724 * t747;
t657 = t674 * t787 - t751 * t675;
t686 = -t724 * qJD(3) + t733 * t787 + t751 * t734;
t701 = mrSges(4,1) * t724 + mrSges(4,2) * t725;
t710 = -mrSges(4,2) * t747 - mrSges(4,3) * t724;
t712 = mrSges(5,1) * t724 - mrSges(5,3) * t747;
t685 = qJD(3) * t725 + t733 * t751 - t734 * t787;
t714 = pkin(4) * t725 - pkin(9) * t747;
t720 = t724 ^ 2;
t738 = g(1) * t753 - t757 * g(2);
t767 = -qJDD(1) * pkin(1) - t738;
t687 = -pkin(2) * t734 + t737 * t775 + (-pkin(8) * t748 - pkin(7)) * t758 + t767;
t760 = (-t686 + t782) * qJ(4) + t687 + (pkin(3) * t747 + t791) * t725;
t639 = -pkin(4) * t720 - t714 * t725 + (pkin(3) + pkin(9)) * t685 + t760;
t651 = -t746 * pkin(3) - t745 * qJ(4) + t725 * t700 + qJDD(4) - t657;
t642 = (t724 * t725 - t746) * pkin(9) + (t686 + t782) * pkin(4) + t651;
t750 = sin(qJ(5));
t755 = cos(qJ(5));
t634 = -t639 * t750 + t755 * t642;
t706 = t724 * t755 - t747 * t750;
t661 = qJD(5) * t706 + t685 * t750 + t746 * t755;
t684 = qJDD(5) + t686;
t707 = t724 * t750 + t747 * t755;
t719 = qJD(5) + t725;
t632 = (t706 * t719 - t661) * pkin(10) + (t706 * t707 + t684) * pkin(5) + t634;
t635 = t755 * t639 + t750 * t642;
t660 = -qJD(5) * t707 + t685 * t755 - t746 * t750;
t690 = pkin(5) * t719 - pkin(10) * t707;
t705 = t706 ^ 2;
t633 = -pkin(5) * t705 + pkin(10) * t660 - t690 * t719 + t635;
t749 = sin(qJ(6));
t754 = cos(qJ(6));
t630 = t632 * t754 - t633 * t749;
t668 = t706 * t754 - t707 * t749;
t647 = qJD(6) * t668 + t660 * t749 + t661 * t754;
t669 = t706 * t749 + t707 * t754;
t656 = -mrSges(7,1) * t668 + mrSges(7,2) * t669;
t717 = qJD(6) + t719;
t662 = -mrSges(7,2) * t717 + mrSges(7,3) * t668;
t676 = qJDD(6) + t684;
t628 = m(7) * t630 + mrSges(7,1) * t676 - mrSges(7,3) * t647 - t656 * t669 + t662 * t717;
t631 = t632 * t749 + t633 * t754;
t646 = -qJD(6) * t669 + t660 * t754 - t661 * t749;
t663 = mrSges(7,1) * t717 - mrSges(7,3) * t669;
t629 = m(7) * t631 - mrSges(7,2) * t676 + mrSges(7,3) * t646 + t656 * t668 - t663 * t717;
t619 = t754 * t628 + t749 * t629;
t673 = -mrSges(6,1) * t706 + mrSges(6,2) * t707;
t688 = -mrSges(6,2) * t719 + mrSges(6,3) * t706;
t617 = m(6) * t634 + mrSges(6,1) * t684 - mrSges(6,3) * t661 - t673 * t707 + t688 * t719 + t619;
t689 = mrSges(6,1) * t719 - mrSges(6,3) * t707;
t768 = -t628 * t749 + t754 * t629;
t618 = m(6) * t635 - mrSges(6,2) * t684 + mrSges(6,3) * t660 + t673 * t706 - t689 * t719 + t768;
t614 = t617 * t755 + t618 * t750;
t702 = -mrSges(5,2) * t724 - mrSges(5,3) * t725;
t764 = -m(5) * t651 - t686 * mrSges(5,1) - t725 * t702 - t614;
t612 = m(4) * t657 - mrSges(4,3) * t686 - t701 * t725 + (t710 - t712) * t747 + (mrSges(4,1) - mrSges(5,2)) * t746 + t764;
t711 = mrSges(4,1) * t747 - mrSges(4,3) * t725;
t713 = mrSges(5,1) * t725 + mrSges(5,2) * t747;
t644 = -pkin(4) * t685 - pkin(9) * t720 + t747 * t714 - t650;
t637 = -pkin(5) * t660 - pkin(10) * t705 + t690 * t707 + t644;
t765 = m(7) * t637 - mrSges(7,1) * t646 + t647 * mrSges(7,2) - t662 * t668 + t669 * t663;
t763 = -m(6) * t644 + mrSges(6,1) * t660 - t661 * mrSges(6,2) + t688 * t706 - t707 * t689 - t765;
t761 = -m(5) * t650 + t746 * mrSges(5,3) + t747 * t713 - t763;
t624 = (-t701 - t702) * t724 - t711 * t747 + (-mrSges(4,3) - mrSges(5,1)) * t685 + m(4) * t658 - mrSges(4,2) * t746 + t761;
t606 = t787 * t612 + t751 * t624;
t708 = -g(3) * t756 - t781;
t732 = (-mrSges(3,1) * t756 + mrSges(3,2) * t752) * qJD(1);
t736 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t774;
t604 = m(3) * t708 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t733 + qJD(2) * t736 - t732 * t775 + t606;
t735 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t775;
t769 = -t612 * t751 + t787 * t624;
t605 = m(3) * t709 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t734 - qJD(2) * t735 + t732 * t774 + t769;
t770 = -t604 * t752 + t756 * t605;
t599 = m(2) * t739 - mrSges(2,1) * t758 - qJDD(1) * mrSges(2,2) + t770;
t726 = -pkin(7) * t758 + t767;
t649 = pkin(3) * t685 + t760;
t779 = -t750 * t617 + t755 * t618;
t613 = m(5) * t649 - t685 * mrSges(5,2) - t686 * mrSges(5,3) - t724 * t712 - t725 * t713 + t779;
t762 = m(4) * t687 + t685 * mrSges(4,1) + mrSges(4,2) * t686 + t724 * t710 + t711 * t725 + t613;
t759 = -m(3) * t726 + t734 * mrSges(3,1) - mrSges(3,2) * t733 - t735 * t775 + t736 * t774 - t762;
t610 = m(2) * t738 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t758 + t759;
t780 = t753 * t599 + t757 * t610;
t600 = t756 * t604 + t752 * t605;
t778 = t788 * t724 - t785 * t725 - t783 * t747;
t777 = t783 * t724 - t784 * t725 + t789 * t747;
t776 = -t785 * t724 + t790 * t725 + t784 * t747;
t771 = t757 * t599 - t610 * t753;
t723 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t752 + Ifges(3,4) * t756) * qJD(1);
t722 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t752 + Ifges(3,2) * t756) * qJD(1);
t721 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t752 + Ifges(3,6) * t756) * qJD(1);
t666 = Ifges(6,1) * t707 + Ifges(6,4) * t706 + Ifges(6,5) * t719;
t665 = Ifges(6,4) * t707 + Ifges(6,2) * t706 + Ifges(6,6) * t719;
t664 = Ifges(6,5) * t707 + Ifges(6,6) * t706 + Ifges(6,3) * t719;
t654 = Ifges(7,1) * t669 + Ifges(7,4) * t668 + Ifges(7,5) * t717;
t653 = Ifges(7,4) * t669 + Ifges(7,2) * t668 + Ifges(7,6) * t717;
t652 = Ifges(7,5) * t669 + Ifges(7,6) * t668 + Ifges(7,3) * t717;
t621 = mrSges(7,2) * t637 - mrSges(7,3) * t630 + Ifges(7,1) * t647 + Ifges(7,4) * t646 + Ifges(7,5) * t676 + t652 * t668 - t653 * t717;
t620 = -mrSges(7,1) * t637 + mrSges(7,3) * t631 + Ifges(7,4) * t647 + Ifges(7,2) * t646 + Ifges(7,6) * t676 - t652 * t669 + t654 * t717;
t608 = mrSges(6,2) * t644 - mrSges(6,3) * t634 + Ifges(6,1) * t661 + Ifges(6,4) * t660 + Ifges(6,5) * t684 - pkin(10) * t619 - t620 * t749 + t621 * t754 + t664 * t706 - t665 * t719;
t607 = -mrSges(6,1) * t644 + mrSges(6,3) * t635 + Ifges(6,4) * t661 + Ifges(6,2) * t660 + Ifges(6,6) * t684 - pkin(5) * t765 + pkin(10) * t768 + t754 * t620 + t749 * t621 - t707 * t664 + t719 * t666;
t596 = pkin(4) * t614 + t790 * t686 + t777 * t724 + t778 * t747 - qJ(4) * t613 + t784 * t746 - t785 * t685 + t707 * t665 - t706 * t666 + Ifges(6,3) * t684 + mrSges(4,2) * t687 + Ifges(7,3) * t676 - t668 * t654 + t669 * t653 + Ifges(6,6) * t660 + Ifges(6,5) * t661 - mrSges(4,3) * t657 + Ifges(7,5) * t647 - mrSges(5,3) * t649 + mrSges(5,1) * t651 + Ifges(7,6) * t646 - mrSges(6,2) * t635 + mrSges(6,1) * t634 - mrSges(7,2) * t631 + pkin(5) * t619 + mrSges(7,1) * t630;
t595 = -mrSges(4,1) * t687 - mrSges(5,1) * t650 + mrSges(5,2) * t649 + mrSges(4,3) * t658 - pkin(3) * t613 - pkin(4) * t763 - pkin(9) * t779 - t755 * t607 - t750 * t608 - t788 * t685 + t785 * t686 + t777 * t725 + t783 * t746 + t776 * t747;
t594 = pkin(9) * t614 + mrSges(2,1) * g(3) + Ifges(2,6) * qJDD(1) + (pkin(3) * mrSges(5,2) + t789) * t746 + (qJ(4) * t702 - t776) * t724 + t778 * t725 + (mrSges(5,1) * qJ(4) + t783) * t685 - t784 * t686 - qJ(4) * t761 - pkin(3) * (-t712 * t747 + t764) + (-t722 * t752 + t723 * t756) * qJD(1) - Ifges(3,3) * qJDD(2) - pkin(2) * t606 - pkin(1) * t600 - t755 * t608 + t758 * Ifges(2,5) + t750 * t607 + mrSges(2,3) * t739 - Ifges(3,5) * t733 - Ifges(3,6) * t734 - mrSges(3,1) * t708 + mrSges(3,2) * t709 - mrSges(4,1) * t657 + mrSges(4,2) * t658 + mrSges(5,3) * t650 - mrSges(5,2) * t651;
t593 = mrSges(3,2) * t726 - mrSges(3,3) * t708 + Ifges(3,1) * t733 + Ifges(3,4) * t734 + Ifges(3,5) * qJDD(2) - pkin(8) * t606 - qJD(2) * t722 - t751 * t595 + t596 * t787 + t721 * t774;
t592 = -mrSges(3,1) * t726 + mrSges(3,3) * t709 + Ifges(3,4) * t733 + Ifges(3,2) * t734 + Ifges(3,6) * qJDD(2) - pkin(2) * t762 + pkin(8) * t769 + qJD(2) * t723 + t595 * t787 + t751 * t596 - t721 * t775;
t591 = -mrSges(2,2) * g(3) - mrSges(2,3) * t738 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t758 - pkin(7) * t600 - t592 * t752 + t593 * t756;
t1 = [-m(1) * g(1) + t771; -m(1) * g(2) + t780; (-m(1) - m(2)) * g(3) + t600; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t780 + t757 * t591 - t753 * t594; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t771 + t753 * t591 + t757 * t594; -mrSges(1,1) * g(2) + mrSges(2,1) * t738 + mrSges(1,2) * g(1) - mrSges(2,2) * t739 + Ifges(2,3) * qJDD(1) + pkin(1) * t759 + pkin(7) * t770 + t756 * t592 + t752 * t593;];
tauB  = t1;
