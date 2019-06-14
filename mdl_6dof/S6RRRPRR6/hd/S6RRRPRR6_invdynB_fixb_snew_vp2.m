% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 11:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 11:00:35
% EndTime: 2019-05-07 11:01:57
% DurationCPUTime: 37.48s
% Computational Cost: add. (592276->385), mult. (1256171->484), div. (0->0), fcn. (928527->12), ass. (0->148)
t727 = sin(qJ(1));
t732 = cos(qJ(1));
t715 = -g(1) * t732 - g(2) * t727;
t734 = qJD(1) ^ 2;
t699 = -pkin(1) * t734 + qJDD(1) * pkin(7) + t715;
t726 = sin(qJ(2));
t731 = cos(qJ(2));
t690 = -g(3) * t726 + t731 * t699;
t707 = (-mrSges(3,1) * t731 + mrSges(3,2) * t726) * qJD(1);
t747 = qJD(1) * qJD(2);
t718 = t726 * t747;
t710 = qJDD(1) * t731 - t718;
t749 = qJD(1) * t726;
t712 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t749;
t714 = g(1) * t727 - t732 * g(2);
t698 = -qJDD(1) * pkin(1) - pkin(7) * t734 - t714;
t746 = t731 * t747;
t709 = qJDD(1) * t726 + t746;
t665 = (-t709 - t746) * pkin(8) + (-t710 + t718) * pkin(2) + t698;
t708 = (-pkin(2) * t731 - pkin(8) * t726) * qJD(1);
t733 = qJD(2) ^ 2;
t748 = qJD(1) * t731;
t668 = -pkin(2) * t733 + qJDD(2) * pkin(8) + t708 * t748 + t690;
t725 = sin(qJ(3));
t730 = cos(qJ(3));
t646 = t730 * t665 - t668 * t725;
t705 = qJD(2) * t730 - t725 * t749;
t681 = qJD(3) * t705 + qJDD(2) * t725 + t709 * t730;
t704 = qJDD(3) - t710;
t706 = qJD(2) * t725 + t730 * t749;
t717 = qJD(3) - t748;
t632 = (t705 * t717 - t681) * qJ(4) + (t705 * t706 + t704) * pkin(3) + t646;
t647 = t725 * t665 + t730 * t668;
t680 = -qJD(3) * t706 + qJDD(2) * t730 - t709 * t725;
t687 = pkin(3) * t717 - qJ(4) * t706;
t703 = t705 ^ 2;
t634 = -pkin(3) * t703 + qJ(4) * t680 - t687 * t717 + t647;
t721 = sin(pkin(11));
t722 = cos(pkin(11));
t684 = t705 * t721 + t706 * t722;
t616 = -0.2e1 * qJD(4) * t684 + t722 * t632 - t634 * t721;
t656 = t680 * t721 + t681 * t722;
t683 = t705 * t722 - t706 * t721;
t609 = (t683 * t717 - t656) * pkin(9) + (t683 * t684 + t704) * pkin(4) + t616;
t617 = 0.2e1 * qJD(4) * t683 + t721 * t632 + t722 * t634;
t655 = t680 * t722 - t681 * t721;
t671 = pkin(4) * t717 - pkin(9) * t684;
t682 = t683 ^ 2;
t615 = -pkin(4) * t682 + pkin(9) * t655 - t671 * t717 + t617;
t724 = sin(qJ(5));
t729 = cos(qJ(5));
t603 = t729 * t609 - t615 * t724;
t660 = t683 * t729 - t684 * t724;
t629 = qJD(5) * t660 + t655 * t724 + t656 * t729;
t661 = t683 * t724 + t684 * t729;
t700 = qJDD(5) + t704;
t716 = qJD(5) + t717;
t601 = (t660 * t716 - t629) * pkin(10) + (t660 * t661 + t700) * pkin(5) + t603;
t604 = t724 * t609 + t729 * t615;
t628 = -qJD(5) * t661 + t655 * t729 - t656 * t724;
t650 = pkin(5) * t716 - pkin(10) * t661;
t659 = t660 ^ 2;
t602 = -pkin(5) * t659 + pkin(10) * t628 - t650 * t716 + t604;
t723 = sin(qJ(6));
t728 = cos(qJ(6));
t599 = t601 * t728 - t602 * t723;
t642 = t660 * t728 - t661 * t723;
t613 = qJD(6) * t642 + t628 * t723 + t629 * t728;
t643 = t660 * t723 + t661 * t728;
t625 = -mrSges(7,1) * t642 + mrSges(7,2) * t643;
t711 = qJD(6) + t716;
t635 = -mrSges(7,2) * t711 + mrSges(7,3) * t642;
t694 = qJDD(6) + t700;
t595 = m(7) * t599 + mrSges(7,1) * t694 - mrSges(7,3) * t613 - t625 * t643 + t635 * t711;
t600 = t601 * t723 + t602 * t728;
t612 = -qJD(6) * t643 + t628 * t728 - t629 * t723;
t636 = mrSges(7,1) * t711 - mrSges(7,3) * t643;
t596 = m(7) * t600 - mrSges(7,2) * t694 + mrSges(7,3) * t612 + t625 * t642 - t636 * t711;
t589 = t728 * t595 + t723 * t596;
t644 = -mrSges(6,1) * t660 + mrSges(6,2) * t661;
t648 = -mrSges(6,2) * t716 + mrSges(6,3) * t660;
t587 = m(6) * t603 + mrSges(6,1) * t700 - mrSges(6,3) * t629 - t644 * t661 + t648 * t716 + t589;
t649 = mrSges(6,1) * t716 - mrSges(6,3) * t661;
t740 = -t595 * t723 + t728 * t596;
t588 = m(6) * t604 - mrSges(6,2) * t700 + mrSges(6,3) * t628 + t644 * t660 - t649 * t716 + t740;
t583 = t729 * t587 + t724 * t588;
t662 = -mrSges(5,1) * t683 + mrSges(5,2) * t684;
t669 = -mrSges(5,2) * t717 + mrSges(5,3) * t683;
t581 = m(5) * t616 + mrSges(5,1) * t704 - mrSges(5,3) * t656 - t662 * t684 + t669 * t717 + t583;
t670 = mrSges(5,1) * t717 - mrSges(5,3) * t684;
t741 = -t587 * t724 + t729 * t588;
t582 = m(5) * t617 - mrSges(5,2) * t704 + mrSges(5,3) * t655 + t662 * t683 - t670 * t717 + t741;
t575 = t722 * t581 + t721 * t582;
t685 = -mrSges(4,1) * t705 + mrSges(4,2) * t706;
t686 = -mrSges(4,2) * t717 + mrSges(4,3) * t705;
t573 = m(4) * t646 + mrSges(4,1) * t704 - mrSges(4,3) * t681 - t685 * t706 + t686 * t717 + t575;
t688 = mrSges(4,1) * t717 - mrSges(4,3) * t706;
t742 = -t581 * t721 + t722 * t582;
t574 = m(4) * t647 - mrSges(4,2) * t704 + mrSges(4,3) * t680 + t685 * t705 - t688 * t717 + t742;
t743 = -t573 * t725 + t730 * t574;
t568 = m(3) * t690 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t710 - qJD(2) * t712 + t707 * t748 + t743;
t689 = -t731 * g(3) - t726 * t699;
t713 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t748;
t667 = -qJDD(2) * pkin(2) - pkin(8) * t733 + t708 * t749 - t689;
t645 = -pkin(3) * t680 - qJ(4) * t703 + t706 * t687 + qJDD(4) + t667;
t619 = -pkin(4) * t655 - pkin(9) * t682 + t684 * t671 + t645;
t606 = -pkin(5) * t628 - pkin(10) * t659 + t650 * t661 + t619;
t739 = m(7) * t606 - t612 * mrSges(7,1) + t613 * mrSges(7,2) - t642 * t635 + t643 * t636;
t738 = m(6) * t619 - t628 * mrSges(6,1) + t629 * mrSges(6,2) - t660 * t648 + t661 * t649 + t739;
t736 = m(5) * t645 - t655 * mrSges(5,1) + t656 * mrSges(5,2) - t683 * t669 + t684 * t670 + t738;
t735 = -m(4) * t667 + t680 * mrSges(4,1) - t681 * mrSges(4,2) + t705 * t686 - t706 * t688 - t736;
t598 = m(3) * t689 + qJDD(2) * mrSges(3,1) - t709 * mrSges(3,3) + qJD(2) * t713 - t707 * t749 + t735;
t744 = t731 * t568 - t598 * t726;
t562 = m(2) * t715 - mrSges(2,1) * t734 - qJDD(1) * mrSges(2,2) + t744;
t569 = t573 * t730 + t574 * t725;
t737 = -m(3) * t698 + t710 * mrSges(3,1) - mrSges(3,2) * t709 - t712 * t749 + t713 * t748 - t569;
t565 = m(2) * t714 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t734 + t737;
t750 = t727 * t562 + t732 * t565;
t563 = t726 * t568 + t731 * t598;
t745 = t732 * t562 - t565 * t727;
t697 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t726 + Ifges(3,4) * t731) * qJD(1);
t696 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t726 + Ifges(3,2) * t731) * qJD(1);
t695 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t726 + Ifges(3,6) * t731) * qJD(1);
t674 = Ifges(4,1) * t706 + Ifges(4,4) * t705 + Ifges(4,5) * t717;
t673 = Ifges(4,4) * t706 + Ifges(4,2) * t705 + Ifges(4,6) * t717;
t672 = Ifges(4,5) * t706 + Ifges(4,6) * t705 + Ifges(4,3) * t717;
t654 = Ifges(5,1) * t684 + Ifges(5,4) * t683 + Ifges(5,5) * t717;
t653 = Ifges(5,4) * t684 + Ifges(5,2) * t683 + Ifges(5,6) * t717;
t652 = Ifges(5,5) * t684 + Ifges(5,6) * t683 + Ifges(5,3) * t717;
t639 = Ifges(6,1) * t661 + Ifges(6,4) * t660 + Ifges(6,5) * t716;
t638 = Ifges(6,4) * t661 + Ifges(6,2) * t660 + Ifges(6,6) * t716;
t637 = Ifges(6,5) * t661 + Ifges(6,6) * t660 + Ifges(6,3) * t716;
t622 = Ifges(7,1) * t643 + Ifges(7,4) * t642 + Ifges(7,5) * t711;
t621 = Ifges(7,4) * t643 + Ifges(7,2) * t642 + Ifges(7,6) * t711;
t620 = Ifges(7,5) * t643 + Ifges(7,6) * t642 + Ifges(7,3) * t711;
t591 = mrSges(7,2) * t606 - mrSges(7,3) * t599 + Ifges(7,1) * t613 + Ifges(7,4) * t612 + Ifges(7,5) * t694 + t620 * t642 - t621 * t711;
t590 = -mrSges(7,1) * t606 + mrSges(7,3) * t600 + Ifges(7,4) * t613 + Ifges(7,2) * t612 + Ifges(7,6) * t694 - t620 * t643 + t622 * t711;
t577 = mrSges(6,2) * t619 - mrSges(6,3) * t603 + Ifges(6,1) * t629 + Ifges(6,4) * t628 + Ifges(6,5) * t700 - pkin(10) * t589 - t590 * t723 + t591 * t728 + t637 * t660 - t638 * t716;
t576 = -mrSges(6,1) * t619 + mrSges(6,3) * t604 + Ifges(6,4) * t629 + Ifges(6,2) * t628 + Ifges(6,6) * t700 - pkin(5) * t739 + pkin(10) * t740 + t728 * t590 + t723 * t591 - t661 * t637 + t716 * t639;
t571 = mrSges(5,2) * t645 - mrSges(5,3) * t616 + Ifges(5,1) * t656 + Ifges(5,4) * t655 + Ifges(5,5) * t704 - pkin(9) * t583 - t576 * t724 + t577 * t729 + t652 * t683 - t653 * t717;
t570 = -mrSges(5,1) * t645 + mrSges(5,3) * t617 + Ifges(5,4) * t656 + Ifges(5,2) * t655 + Ifges(5,6) * t704 - pkin(4) * t738 + pkin(9) * t741 + t729 * t576 + t724 * t577 - t684 * t652 + t717 * t654;
t559 = -pkin(2) * t569 + (-Ifges(4,3) - Ifges(5,3)) * t704 - t695 * t749 + Ifges(3,4) * t709 + Ifges(3,2) * t710 + t705 * t674 - t706 * t673 + mrSges(3,3) * t690 - Ifges(7,3) * t694 + qJD(2) * t697 - mrSges(3,1) * t698 - Ifges(6,3) * t700 - Ifges(4,6) * t680 - Ifges(4,5) * t681 + t683 * t654 - t684 * t653 - Ifges(5,5) * t656 + t660 * t639 - t661 * t638 - Ifges(5,6) * t655 + mrSges(4,2) * t647 - mrSges(4,1) * t646 + t642 * t622 - t643 * t621 - Ifges(6,6) * t628 - Ifges(6,5) * t629 - mrSges(5,1) * t616 + mrSges(5,2) * t617 - Ifges(7,6) * t612 - Ifges(7,5) * t613 + mrSges(6,2) * t604 - mrSges(6,1) * t603 + mrSges(7,2) * t600 - mrSges(7,1) * t599 - pkin(5) * t589 - pkin(4) * t583 + Ifges(3,6) * qJDD(2) - pkin(3) * t575;
t558 = mrSges(4,2) * t667 - mrSges(4,3) * t646 + Ifges(4,1) * t681 + Ifges(4,4) * t680 + Ifges(4,5) * t704 - qJ(4) * t575 - t570 * t721 + t571 * t722 + t672 * t705 - t673 * t717;
t557 = -mrSges(4,1) * t667 + mrSges(4,3) * t647 + Ifges(4,4) * t681 + Ifges(4,2) * t680 + Ifges(4,6) * t704 - pkin(3) * t736 + qJ(4) * t742 + t722 * t570 + t721 * t571 - t706 * t672 + t717 * t674;
t556 = mrSges(3,2) * t698 - mrSges(3,3) * t689 + Ifges(3,1) * t709 + Ifges(3,4) * t710 + Ifges(3,5) * qJDD(2) - pkin(8) * t569 - qJD(2) * t696 - t557 * t725 + t558 * t730 + t695 * t748;
t555 = Ifges(2,6) * qJDD(1) + t734 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t715 - Ifges(3,5) * t709 - Ifges(3,6) * t710 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t689 + mrSges(3,2) * t690 - t725 * t558 - t730 * t557 - pkin(2) * t735 - pkin(8) * t743 - pkin(1) * t563 + (-t696 * t726 + t697 * t731) * qJD(1);
t554 = -mrSges(2,2) * g(3) - mrSges(2,3) * t714 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t734 - pkin(7) * t563 + t556 * t731 - t559 * t726;
t1 = [-m(1) * g(1) + t745; -m(1) * g(2) + t750; (-m(1) - m(2)) * g(3) + t563; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t750 + t732 * t554 - t727 * t555; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t745 + t727 * t554 + t732 * t555; -mrSges(1,1) * g(2) + mrSges(2,1) * t714 + mrSges(1,2) * g(1) - mrSges(2,2) * t715 + Ifges(2,3) * qJDD(1) + pkin(1) * t737 + pkin(7) * t744 + t726 * t556 + t731 * t559;];
tauB  = t1;
