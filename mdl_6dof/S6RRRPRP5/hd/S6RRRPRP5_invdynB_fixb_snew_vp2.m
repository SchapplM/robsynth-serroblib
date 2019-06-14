% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:52:12
% EndTime: 2019-05-07 07:52:32
% DurationCPUTime: 14.83s
% Computational Cost: add. (223945->361), mult. (468211->442), div. (0->0), fcn. (335291->10), ass. (0->138)
t751 = Ifges(6,1) + Ifges(7,1);
t746 = Ifges(6,4) - Ifges(7,5);
t745 = Ifges(7,4) + Ifges(6,5);
t750 = Ifges(6,2) + Ifges(7,3);
t744 = Ifges(6,6) - Ifges(7,6);
t749 = -Ifges(6,3) - Ifges(7,2);
t748 = cos(qJ(5));
t747 = -mrSges(6,3) - mrSges(7,2);
t718 = sin(qJ(1));
t721 = cos(qJ(1));
t706 = -g(1) * t721 - g(2) * t718;
t723 = qJD(1) ^ 2;
t691 = -pkin(1) * t723 + qJDD(1) * pkin(7) + t706;
t717 = sin(qJ(2));
t720 = cos(qJ(2));
t681 = -g(3) * t717 + t720 * t691;
t699 = (-mrSges(3,1) * t720 + mrSges(3,2) * t717) * qJD(1);
t736 = qJD(1) * qJD(2);
t710 = t717 * t736;
t702 = qJDD(1) * t720 - t710;
t738 = qJD(1) * t717;
t703 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t738;
t705 = t718 * g(1) - t721 * g(2);
t690 = -qJDD(1) * pkin(1) - t723 * pkin(7) - t705;
t734 = t720 * t736;
t701 = qJDD(1) * t717 + t734;
t655 = (-t701 - t734) * pkin(8) + (-t702 + t710) * pkin(2) + t690;
t700 = (-pkin(2) * t720 - pkin(8) * t717) * qJD(1);
t722 = qJD(2) ^ 2;
t737 = qJD(1) * t720;
t658 = -pkin(2) * t722 + qJDD(2) * pkin(8) + t700 * t737 + t681;
t716 = sin(qJ(3));
t719 = cos(qJ(3));
t633 = t719 * t655 - t716 * t658;
t697 = qJD(2) * t719 - t716 * t738;
t672 = qJD(3) * t697 + qJDD(2) * t716 + t701 * t719;
t696 = qJDD(3) - t702;
t698 = qJD(2) * t716 + t719 * t738;
t709 = qJD(3) - t737;
t616 = (t697 * t709 - t672) * qJ(4) + (t697 * t698 + t696) * pkin(3) + t633;
t634 = t716 * t655 + t719 * t658;
t671 = -qJD(3) * t698 + qJDD(2) * t719 - t701 * t716;
t678 = pkin(3) * t709 - qJ(4) * t698;
t695 = t697 ^ 2;
t618 = -pkin(3) * t695 + qJ(4) * t671 - t678 * t709 + t634;
t713 = sin(pkin(10));
t714 = cos(pkin(10));
t675 = t697 * t713 + t698 * t714;
t604 = -0.2e1 * qJD(4) * t675 + t714 * t616 - t713 * t618;
t647 = t671 * t713 + t672 * t714;
t674 = t697 * t714 - t698 * t713;
t601 = (t674 * t709 - t647) * pkin(9) + (t674 * t675 + t696) * pkin(4) + t604;
t605 = 0.2e1 * qJD(4) * t674 + t713 * t616 + t714 * t618;
t646 = t671 * t714 - t672 * t713;
t661 = pkin(4) * t709 - pkin(9) * t675;
t673 = t674 ^ 2;
t603 = -pkin(4) * t673 + pkin(9) * t646 - t661 * t709 + t605;
t715 = sin(qJ(5));
t597 = t715 * t601 + t603 * t748;
t651 = t715 * t674 + t675 * t748;
t612 = qJD(5) * t651 - t646 * t748 + t647 * t715;
t708 = qJD(5) + t709;
t639 = mrSges(6,1) * t708 - mrSges(6,3) * t651;
t650 = -t674 * t748 + t675 * t715;
t692 = qJDD(5) + t696;
t629 = pkin(5) * t650 - qJ(6) * t651;
t707 = t708 ^ 2;
t594 = -pkin(5) * t707 + qJ(6) * t692 + 0.2e1 * qJD(6) * t708 - t629 * t650 + t597;
t640 = -mrSges(7,1) * t708 + mrSges(7,2) * t651;
t735 = m(7) * t594 + t692 * mrSges(7,3) + t708 * t640;
t630 = mrSges(7,1) * t650 - mrSges(7,3) * t651;
t739 = -mrSges(6,1) * t650 - mrSges(6,2) * t651 - t630;
t587 = m(6) * t597 - t692 * mrSges(6,2) + t612 * t747 - t708 * t639 + t650 * t739 + t735;
t596 = t601 * t748 - t715 * t603;
t613 = -t650 * qJD(5) + t715 * t646 + t647 * t748;
t638 = -mrSges(6,2) * t708 - mrSges(6,3) * t650;
t595 = -t692 * pkin(5) - t707 * qJ(6) + t651 * t629 + qJDD(6) - t596;
t637 = -mrSges(7,2) * t650 + mrSges(7,3) * t708;
t728 = -m(7) * t595 + t692 * mrSges(7,1) + t708 * t637;
t589 = m(6) * t596 + t692 * mrSges(6,1) + t613 * t747 + t708 * t638 + t651 * t739 + t728;
t582 = t715 * t587 + t589 * t748;
t652 = -mrSges(5,1) * t674 + mrSges(5,2) * t675;
t659 = -mrSges(5,2) * t709 + mrSges(5,3) * t674;
t580 = m(5) * t604 + mrSges(5,1) * t696 - mrSges(5,3) * t647 - t652 * t675 + t659 * t709 + t582;
t660 = mrSges(5,1) * t709 - mrSges(5,3) * t675;
t729 = t587 * t748 - t589 * t715;
t581 = m(5) * t605 - mrSges(5,2) * t696 + mrSges(5,3) * t646 + t652 * t674 - t660 * t709 + t729;
t576 = t714 * t580 + t713 * t581;
t676 = -mrSges(4,1) * t697 + mrSges(4,2) * t698;
t677 = -mrSges(4,2) * t709 + mrSges(4,3) * t697;
t574 = m(4) * t633 + mrSges(4,1) * t696 - mrSges(4,3) * t672 - t676 * t698 + t677 * t709 + t576;
t679 = mrSges(4,1) * t709 - mrSges(4,3) * t698;
t730 = -t580 * t713 + t714 * t581;
t575 = m(4) * t634 - mrSges(4,2) * t696 + mrSges(4,3) * t671 + t676 * t697 - t679 * t709 + t730;
t731 = -t574 * t716 + t719 * t575;
t569 = m(3) * t681 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t702 - qJD(2) * t703 + t699 * t737 + t731;
t680 = -t720 * g(3) - t717 * t691;
t704 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t737;
t657 = -qJDD(2) * pkin(2) - t722 * pkin(8) + t700 * t738 - t680;
t632 = -t671 * pkin(3) - t695 * qJ(4) + t698 * t678 + qJDD(4) + t657;
t607 = -t646 * pkin(4) - t673 * pkin(9) + t675 * t661 + t632;
t599 = t607 - 0.2e1 * qJD(6) * t651 + (t651 * t708 + t612) * pkin(5) + (t650 * t708 - t613) * qJ(6);
t592 = m(7) * t599 + t612 * mrSges(7,1) - t613 * mrSges(7,3) + t650 * t637 - t651 * t640;
t727 = m(6) * t607 + t612 * mrSges(6,1) + t613 * mrSges(6,2) + t650 * t638 + t651 * t639 + t592;
t725 = m(5) * t632 - t646 * mrSges(5,1) + t647 * mrSges(5,2) - t674 * t659 + t675 * t660 + t727;
t724 = -m(4) * t657 + t671 * mrSges(4,1) - t672 * mrSges(4,2) + t697 * t677 - t698 * t679 - t725;
t591 = m(3) * t680 + qJDD(2) * mrSges(3,1) - t701 * mrSges(3,3) + qJD(2) * t704 - t699 * t738 + t724;
t732 = t720 * t569 - t591 * t717;
t563 = m(2) * t706 - mrSges(2,1) * t723 - qJDD(1) * mrSges(2,2) + t732;
t570 = t574 * t719 + t575 * t716;
t726 = -m(3) * t690 + t702 * mrSges(3,1) - mrSges(3,2) * t701 - t703 * t738 + t704 * t737 - t570;
t566 = m(2) * t705 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t723 + t726;
t743 = t718 * t563 + t721 * t566;
t564 = t717 * t569 + t720 * t591;
t742 = t750 * t650 - t746 * t651 - t744 * t708;
t741 = t744 * t650 - t745 * t651 + t749 * t708;
t740 = -t746 * t650 + t751 * t651 + t745 * t708;
t733 = t721 * t563 - t566 * t718;
t689 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t717 + Ifges(3,4) * t720) * qJD(1);
t688 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t717 + Ifges(3,2) * t720) * qJD(1);
t687 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t717 + Ifges(3,6) * t720) * qJD(1);
t664 = Ifges(4,1) * t698 + Ifges(4,4) * t697 + Ifges(4,5) * t709;
t663 = Ifges(4,4) * t698 + Ifges(4,2) * t697 + Ifges(4,6) * t709;
t662 = Ifges(4,5) * t698 + Ifges(4,6) * t697 + Ifges(4,3) * t709;
t645 = Ifges(5,1) * t675 + Ifges(5,4) * t674 + Ifges(5,5) * t709;
t644 = Ifges(5,4) * t675 + Ifges(5,2) * t674 + Ifges(5,6) * t709;
t643 = Ifges(5,5) * t675 + Ifges(5,6) * t674 + Ifges(5,3) * t709;
t584 = mrSges(6,2) * t607 + mrSges(7,2) * t595 - mrSges(6,3) * t596 - mrSges(7,3) * t599 - qJ(6) * t592 - t746 * t612 + t751 * t613 + t741 * t650 + t745 * t692 + t742 * t708;
t583 = -mrSges(6,1) * t607 - mrSges(7,1) * t599 + mrSges(7,2) * t594 + mrSges(6,3) * t597 - pkin(5) * t592 - t750 * t612 + t746 * t613 + t741 * t651 + t744 * t692 + t740 * t708;
t572 = mrSges(5,2) * t632 - mrSges(5,3) * t604 + Ifges(5,1) * t647 + Ifges(5,4) * t646 + Ifges(5,5) * t696 - pkin(9) * t582 - t715 * t583 + t584 * t748 + t674 * t643 - t709 * t644;
t571 = -mrSges(5,1) * t632 + mrSges(5,3) * t605 + Ifges(5,4) * t647 + Ifges(5,2) * t646 + Ifges(5,6) * t696 - pkin(4) * t727 + pkin(9) * t729 + t583 * t748 + t715 * t584 - t675 * t643 + t709 * t645;
t560 = mrSges(4,2) * t657 - mrSges(4,3) * t633 + Ifges(4,1) * t672 + Ifges(4,4) * t671 + Ifges(4,5) * t696 - qJ(4) * t576 - t571 * t713 + t572 * t714 + t662 * t697 - t663 * t709;
t559 = -mrSges(4,1) * t657 + mrSges(4,3) * t634 + Ifges(4,4) * t672 + Ifges(4,2) * t671 + Ifges(4,6) * t696 - pkin(3) * t725 + qJ(4) * t730 + t714 * t571 + t713 * t572 - t698 * t662 + t709 * t664;
t558 = (mrSges(7,2) * qJ(6) + t744) * t612 + (mrSges(7,2) * pkin(5) - t745) * t613 - t687 * t738 + (qJ(6) * t630 - t740) * t650 + (pkin(5) * t630 + t742) * t651 - qJ(6) * t735 - pkin(5) * t728 + Ifges(3,6) * qJDD(2) + (-Ifges(4,3) - Ifges(5,3)) * t696 - t698 * t663 + Ifges(3,4) * t701 + Ifges(3,2) * t702 + t697 * t664 + mrSges(3,3) * t681 + qJD(2) * t689 - mrSges(3,1) * t690 - Ifges(4,6) * t671 - Ifges(4,5) * t672 + t674 * t645 - t675 * t644 - Ifges(5,6) * t646 - Ifges(5,5) * t647 - mrSges(4,1) * t633 + mrSges(4,2) * t634 - mrSges(5,1) * t604 + mrSges(5,2) * t605 - mrSges(6,1) * t596 + mrSges(6,2) * t597 - mrSges(7,3) * t594 + mrSges(7,1) * t595 - pkin(4) * t582 + t749 * t692 - pkin(3) * t576 - pkin(2) * t570;
t557 = mrSges(3,2) * t690 - mrSges(3,3) * t680 + Ifges(3,1) * t701 + Ifges(3,4) * t702 + Ifges(3,5) * qJDD(2) - pkin(8) * t570 - qJD(2) * t688 - t559 * t716 + t560 * t719 + t687 * t737;
t556 = Ifges(2,6) * qJDD(1) + t723 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t706 - Ifges(3,5) * t701 - Ifges(3,6) * t702 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t680 + mrSges(3,2) * t681 - t716 * t560 - t719 * t559 - pkin(2) * t724 - pkin(8) * t731 - pkin(1) * t564 + (-t688 * t717 + t689 * t720) * qJD(1);
t555 = -mrSges(2,2) * g(3) - mrSges(2,3) * t705 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t723 - pkin(7) * t564 + t557 * t720 - t558 * t717;
t1 = [-m(1) * g(1) + t733; -m(1) * g(2) + t743; (-m(1) - m(2)) * g(3) + t564; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t743 + t721 * t555 - t718 * t556; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t733 + t718 * t555 + t721 * t556; -mrSges(1,1) * g(2) + mrSges(2,1) * t705 + mrSges(1,2) * g(1) - mrSges(2,2) * t706 + Ifges(2,3) * qJDD(1) + pkin(1) * t726 + pkin(7) * t732 + t717 * t557 + t720 * t558;];
tauB  = t1;
