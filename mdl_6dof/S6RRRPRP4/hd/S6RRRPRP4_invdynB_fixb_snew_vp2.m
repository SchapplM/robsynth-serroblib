% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 07:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP4_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP4_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP4_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:44:58
% EndTime: 2019-05-07 07:45:11
% DurationCPUTime: 6.02s
% Computational Cost: add. (65493->343), mult. (134776->402), div. (0->0), fcn. (90278->8), ass. (0->133)
t792 = -2 * qJD(4);
t791 = Ifges(4,1) + Ifges(5,2);
t790 = -Ifges(5,1) - Ifges(4,3);
t789 = Ifges(6,1) + Ifges(7,1);
t781 = Ifges(6,4) - Ifges(7,5);
t780 = Ifges(7,4) + Ifges(6,5);
t779 = Ifges(4,5) - Ifges(5,4);
t788 = Ifges(4,2) + Ifges(5,3);
t787 = Ifges(6,2) + Ifges(7,3);
t778 = Ifges(4,6) - Ifges(5,5);
t777 = -Ifges(5,6) - Ifges(4,4);
t776 = Ifges(6,6) - Ifges(7,6);
t786 = Ifges(6,3) + Ifges(7,2);
t743 = sin(qJ(2));
t745 = cos(qJ(2));
t761 = qJD(1) * qJD(2);
t725 = qJDD(1) * t743 + t745 * t761;
t744 = sin(qJ(1));
t746 = cos(qJ(1));
t731 = -g(1) * t746 - g(2) * t744;
t747 = qJD(1) ^ 2;
t719 = -pkin(1) * t747 + qJDD(1) * pkin(7) + t731;
t774 = t719 * t743;
t783 = pkin(2) * t747;
t663 = qJDD(2) * pkin(2) - pkin(8) * t725 - t774 + (pkin(8) * t761 + t743 * t783 - g(3)) * t745;
t700 = -g(3) * t743 + t745 * t719;
t726 = qJDD(1) * t745 - t743 * t761;
t764 = qJD(1) * t743;
t729 = qJD(2) * pkin(2) - pkin(8) * t764;
t740 = t745 ^ 2;
t664 = pkin(8) * t726 - qJD(2) * t729 - t740 * t783 + t700;
t742 = sin(qJ(3));
t785 = cos(qJ(3));
t642 = t742 * t663 + t785 * t664;
t763 = qJD(1) * t745;
t716 = t742 * t764 - t785 * t763;
t717 = (t742 * t745 + t785 * t743) * qJD(1);
t692 = pkin(3) * t716 - qJ(4) * t717;
t739 = qJD(2) + qJD(3);
t737 = t739 ^ 2;
t738 = qJDD(2) + qJDD(3);
t639 = pkin(3) * t737 - t738 * qJ(4) + t716 * t692 + t739 * t792 - t642;
t784 = cos(qJ(5));
t782 = -mrSges(6,3) - mrSges(7,2);
t775 = t716 * t739;
t641 = t785 * t663 - t742 * t664;
t677 = -t716 * qJD(3) + t785 * t725 + t742 * t726;
t693 = mrSges(4,1) * t716 + mrSges(4,2) * t717;
t701 = -mrSges(4,2) * t739 - mrSges(4,3) * t716;
t703 = mrSges(5,1) * t716 - mrSges(5,3) * t739;
t676 = qJD(3) * t717 + t725 * t742 - t785 * t726;
t705 = pkin(4) * t717 - pkin(9) * t739;
t712 = t716 ^ 2;
t730 = g(1) * t744 - t746 * g(2);
t755 = -qJDD(1) * pkin(1) - t730;
t678 = -pkin(2) * t726 + t729 * t764 + (-pkin(8) * t740 - pkin(7)) * t747 + t755;
t749 = (-t677 + t775) * qJ(4) + t678 + (pkin(3) * t739 + t792) * t717;
t632 = -pkin(4) * t712 - t705 * t717 + (pkin(3) + pkin(9)) * t676 + t749;
t640 = -t738 * pkin(3) - t737 * qJ(4) + t717 * t692 + qJDD(4) - t641;
t634 = (t716 * t717 - t738) * pkin(9) + (t677 + t775) * pkin(4) + t640;
t741 = sin(qJ(5));
t628 = t784 * t632 + t741 * t634;
t698 = t741 * t716 + t784 * t739;
t645 = qJD(5) * t698 - t784 * t676 + t738 * t741;
t675 = qJDD(5) + t677;
t711 = qJD(5) + t717;
t681 = mrSges(6,1) * t711 - mrSges(6,3) * t698;
t697 = -t784 * t716 + t739 * t741;
t660 = pkin(5) * t697 - qJ(6) * t698;
t710 = t711 ^ 2;
t625 = -pkin(5) * t710 + qJ(6) * t675 + 0.2e1 * qJD(6) * t711 - t660 * t697 + t628;
t682 = -mrSges(7,1) * t711 + mrSges(7,2) * t698;
t760 = m(7) * t625 + t675 * mrSges(7,3) + t711 * t682;
t661 = mrSges(7,1) * t697 - mrSges(7,3) * t698;
t768 = -mrSges(6,1) * t697 - mrSges(6,2) * t698 - t661;
t620 = m(6) * t628 - mrSges(6,2) * t675 + t782 * t645 - t681 * t711 + t768 * t697 + t760;
t627 = -t741 * t632 + t784 * t634;
t646 = -t697 * qJD(5) + t741 * t676 + t784 * t738;
t680 = -mrSges(6,2) * t711 - mrSges(6,3) * t697;
t626 = -t675 * pkin(5) - t710 * qJ(6) + t698 * t660 + qJDD(6) - t627;
t679 = -mrSges(7,2) * t697 + mrSges(7,3) * t711;
t756 = -m(7) * t626 + t675 * mrSges(7,1) + t711 * t679;
t622 = m(6) * t627 + mrSges(6,1) * t675 + t782 * t646 + t680 * t711 + t768 * t698 + t756;
t615 = t741 * t620 + t784 * t622;
t694 = -mrSges(5,2) * t716 - mrSges(5,3) * t717;
t753 = -m(5) * t640 - t677 * mrSges(5,1) - t717 * t694 - t615;
t611 = m(4) * t641 - t677 * mrSges(4,3) - t717 * t693 + (t701 - t703) * t739 + (mrSges(4,1) - mrSges(5,2)) * t738 + t753;
t702 = mrSges(4,1) * t739 - mrSges(4,3) * t717;
t704 = mrSges(5,1) * t717 + mrSges(5,2) * t739;
t636 = -pkin(4) * t676 - pkin(9) * t712 + t739 * t705 - t639;
t630 = -0.2e1 * qJD(6) * t698 + (t697 * t711 - t646) * qJ(6) + (t698 * t711 + t645) * pkin(5) + t636;
t623 = m(7) * t630 + t645 * mrSges(7,1) - mrSges(7,3) * t646 + t697 * t679 - t682 * t698;
t752 = m(6) * t636 + mrSges(6,1) * t645 + t646 * mrSges(6,2) + t680 * t697 + t698 * t681 + t623;
t750 = -m(5) * t639 + t738 * mrSges(5,3) + t739 * t704 + t752;
t618 = (-t693 - t694) * t716 + (-mrSges(4,3) - mrSges(5,1)) * t676 + m(4) * t642 - mrSges(4,2) * t738 - t702 * t739 + t750;
t607 = t785 * t611 + t742 * t618;
t699 = -g(3) * t745 - t774;
t724 = (-mrSges(3,1) * t745 + mrSges(3,2) * t743) * qJD(1);
t728 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t763;
t605 = m(3) * t699 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t725 + qJD(2) * t728 - t724 * t764 + t607;
t727 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t764;
t757 = -t611 * t742 + t785 * t618;
t606 = m(3) * t700 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t726 - qJD(2) * t727 + t724 * t763 + t757;
t758 = -t605 * t743 + t745 * t606;
t600 = m(2) * t731 - mrSges(2,1) * t747 - qJDD(1) * mrSges(2,2) + t758;
t718 = -pkin(7) * t747 + t755;
t638 = pkin(3) * t676 + t749;
t772 = t784 * t620 - t741 * t622;
t612 = m(5) * t638 - t676 * mrSges(5,2) - t677 * mrSges(5,3) - t716 * t703 - t717 * t704 + t772;
t751 = m(4) * t678 + t676 * mrSges(4,1) + mrSges(4,2) * t677 + t716 * t701 + t702 * t717 + t612;
t748 = -m(3) * t718 + t726 * mrSges(3,1) - mrSges(3,2) * t725 - t727 * t764 + t728 * t763 - t751;
t609 = m(2) * t730 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t747 + t748;
t773 = t744 * t600 + t746 * t609;
t601 = t745 * t605 + t743 * t606;
t771 = t787 * t697 - t781 * t698 - t776 * t711;
t770 = t776 * t697 - t780 * t698 - t786 * t711;
t769 = -t781 * t697 + t789 * t698 + t780 * t711;
t767 = t788 * t716 + t777 * t717 - t778 * t739;
t766 = t778 * t716 - t779 * t717 + t790 * t739;
t765 = t777 * t716 + t791 * t717 + t779 * t739;
t759 = t746 * t600 - t609 * t744;
t715 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t743 + Ifges(3,4) * t745) * qJD(1);
t714 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t743 + Ifges(3,2) * t745) * qJD(1);
t713 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t743 + Ifges(3,6) * t745) * qJD(1);
t614 = mrSges(6,2) * t636 + mrSges(7,2) * t626 - mrSges(6,3) * t627 - mrSges(7,3) * t630 - qJ(6) * t623 - t781 * t645 + t789 * t646 + t780 * t675 + t770 * t697 + t771 * t711;
t613 = -mrSges(6,1) * t636 - mrSges(7,1) * t630 + mrSges(7,2) * t625 + mrSges(6,3) * t628 - pkin(5) * t623 - t787 * t645 + t781 * t646 + t776 * t675 + t770 * t698 + t769 * t711;
t597 = pkin(5) * t756 + qJ(6) * t760 + mrSges(4,2) * t678 - mrSges(5,3) * t638 + mrSges(5,1) * t640 - mrSges(4,3) * t641 - mrSges(6,2) * t628 - mrSges(7,1) * t626 + mrSges(6,1) * t627 + mrSges(7,3) * t625 + pkin(4) * t615 - qJ(4) * t612 + t767 * t739 + t779 * t738 + t766 * t716 + (-pkin(5) * t661 - t771) * t698 + (-qJ(6) * t661 + t769) * t697 + t791 * t677 + t777 * t676 + t786 * t675 + (-pkin(5) * mrSges(7,2) + t780) * t646 + (-mrSges(7,2) * qJ(6) - t776) * t645;
t596 = -mrSges(4,1) * t678 - mrSges(5,1) * t639 + mrSges(5,2) * t638 + mrSges(4,3) * t642 - pkin(3) * t612 + pkin(4) * t752 - pkin(9) * t772 - t784 * t613 - t741 * t614 - t788 * t676 - t777 * t677 + t766 * t717 + t778 * t738 + t765 * t739;
t595 = -t784 * t614 + (qJ(4) * mrSges(5,1) + t778) * t676 - t779 * t677 + (qJ(4) * t694 - t765) * t716 + t767 * t717 - Ifges(3,3) * qJDD(2) + mrSges(2,1) * g(3) - qJ(4) * t750 - pkin(3) * (-t739 * t703 + t753) + t747 * Ifges(2,5) + t741 * t613 - Ifges(3,5) * t725 - Ifges(3,6) * t726 + mrSges(2,3) * t731 - mrSges(3,1) * t699 + mrSges(3,2) * t700 + (pkin(3) * mrSges(5,2) + t790) * t738 + mrSges(5,3) * t639 - mrSges(5,2) * t640 - mrSges(4,1) * t641 + mrSges(4,2) * t642 + pkin(9) * t615 - pkin(2) * t607 + (-t743 * t714 + t745 * t715) * qJD(1) - pkin(1) * t601 + Ifges(2,6) * qJDD(1);
t594 = mrSges(3,2) * t718 - mrSges(3,3) * t699 + Ifges(3,1) * t725 + Ifges(3,4) * t726 + Ifges(3,5) * qJDD(2) - pkin(8) * t607 - qJD(2) * t714 - t742 * t596 + t785 * t597 + t713 * t763;
t593 = -mrSges(3,1) * t718 + mrSges(3,3) * t700 + Ifges(3,4) * t725 + Ifges(3,2) * t726 + Ifges(3,6) * qJDD(2) - pkin(2) * t751 + pkin(8) * t757 + qJD(2) * t715 + t785 * t596 + t742 * t597 - t713 * t764;
t592 = -mrSges(2,2) * g(3) - mrSges(2,3) * t730 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t747 - pkin(7) * t601 - t593 * t743 + t594 * t745;
t1 = [-m(1) * g(1) + t759; -m(1) * g(2) + t773; (-m(1) - m(2)) * g(3) + t601; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t773 + t746 * t592 - t744 * t595; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t759 + t744 * t592 + t746 * t595; -mrSges(1,1) * g(2) + mrSges(2,1) * t730 + mrSges(1,2) * g(1) - mrSges(2,2) * t731 + Ifges(2,3) * qJDD(1) + pkin(1) * t748 + pkin(7) * t758 + t745 * t593 + t743 * t594;];
tauB  = t1;
