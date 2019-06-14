% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-05-05 04:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:07:15
% EndTime: 2019-05-05 04:07:23
% DurationCPUTime: 4.65s
% Computational Cost: add. (51767->302), mult. (100782->359), div. (0->0), fcn. (59606->10), ass. (0->132)
t746 = -2 * qJD(4);
t745 = Ifges(4,1) + Ifges(5,2);
t744 = Ifges(6,1) + Ifges(7,1);
t735 = Ifges(4,4) + Ifges(5,6);
t734 = Ifges(6,4) - Ifges(7,5);
t733 = Ifges(7,4) + Ifges(6,5);
t732 = Ifges(4,5) - Ifges(5,4);
t743 = Ifges(4,2) + Ifges(5,3);
t742 = Ifges(6,2) + Ifges(7,3);
t731 = Ifges(4,6) - Ifges(5,5);
t730 = Ifges(6,6) - Ifges(7,6);
t741 = (Ifges(4,3) + Ifges(5,1));
t740 = Ifges(6,3) + Ifges(7,2);
t684 = sin(pkin(10));
t686 = cos(pkin(10));
t666 = g(1) * t684 - g(2) * t686;
t667 = -g(1) * t686 - g(2) * t684;
t683 = -g(3) + qJDD(1);
t693 = cos(qJ(2));
t687 = cos(pkin(6));
t690 = sin(qJ(2));
t726 = t687 * t690;
t685 = sin(pkin(6));
t727 = t685 * t690;
t612 = t666 * t726 + t693 * t667 + t683 * t727;
t695 = qJD(2) ^ 2;
t610 = -pkin(2) * t695 + qJDD(2) * pkin(8) + t612;
t635 = -t666 * t685 + t683 * t687;
t689 = sin(qJ(3));
t692 = cos(qJ(3));
t606 = t692 * t610 + t689 * t635;
t659 = (-pkin(3) * t692 - qJ(4) * t689) * qJD(2);
t694 = qJD(3) ^ 2;
t715 = qJD(2) * t692;
t602 = pkin(3) * t694 - qJDD(3) * qJ(4) + (qJD(3) * t746) - t659 * t715 - t606;
t611 = -t690 * t667 + (t666 * t687 + t683 * t685) * t693;
t739 = -pkin(3) - pkin(9);
t738 = pkin(8) * t695;
t737 = pkin(9) * t695;
t736 = -mrSges(6,3) - mrSges(7,2);
t700 = -qJDD(2) * pkin(2) - t611;
t609 = t700 - t738;
t714 = qJD(2) * qJD(3);
t711 = t692 * t714;
t662 = qJDD(2) * t689 + t711;
t710 = t689 * t714;
t663 = qJDD(2) * t692 - t710;
t716 = qJD(2) * t689;
t668 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t716;
t669 = -(qJD(3) * mrSges(4,2)) + mrSges(4,3) * t715;
t670 = -mrSges(5,1) * t715 - (qJD(3) * mrSges(5,3));
t698 = pkin(3) * t710 + t716 * t746 + (-t662 - t711) * qJ(4) + t700;
t604 = -pkin(3) * t663 + t698 - t738;
t671 = mrSges(5,1) * t716 + (qJD(3) * mrSges(5,2));
t607 = t689 * t610;
t705 = -qJ(4) * t694 + t659 * t716 + qJDD(4) + t607;
t599 = pkin(4) * t662 + t739 * qJDD(3) + (-pkin(4) * t714 - t689 * t737 - t635) * t692 + t705;
t673 = pkin(4) * t716 - qJD(3) * pkin(9);
t682 = t692 ^ 2;
t601 = -t673 * t716 + (-pkin(4) * t682 - pkin(8)) * t695 + t739 * t663 + t698;
t688 = sin(qJ(5));
t691 = cos(qJ(5));
t594 = t688 * t599 + t691 * t601;
t658 = qJD(3) * t691 - t688 * t715;
t623 = qJD(5) * t658 + qJDD(3) * t688 + t691 * t663;
t676 = qJD(5) + t716;
t633 = mrSges(6,1) * t676 - mrSges(6,3) * t658;
t654 = qJDD(5) + t662;
t657 = qJD(3) * t688 + t691 * t715;
t627 = pkin(5) * t657 - qJ(6) * t658;
t674 = t676 ^ 2;
t590 = -pkin(5) * t674 + qJ(6) * t654 + 0.2e1 * qJD(6) * t676 - t627 * t657 + t594;
t634 = -mrSges(7,1) * t676 + mrSges(7,2) * t658;
t712 = m(7) * t590 + t654 * mrSges(7,3) + t676 * t634;
t628 = mrSges(7,1) * t657 - mrSges(7,3) * t658;
t720 = -mrSges(6,1) * t657 - mrSges(6,2) * t658 - t628;
t585 = m(6) * t594 - mrSges(6,2) * t654 + t623 * t736 - t633 * t676 + t657 * t720 + t712;
t593 = t599 * t691 - t601 * t688;
t624 = -qJD(5) * t657 + qJDD(3) * t691 - t663 * t688;
t631 = -mrSges(6,2) * t676 - mrSges(6,3) * t657;
t591 = -pkin(5) * t654 - qJ(6) * t674 + t627 * t658 + qJDD(6) - t593;
t632 = -mrSges(7,2) * t657 + mrSges(7,3) * t676;
t707 = -m(7) * t591 + t654 * mrSges(7,1) + t676 * t632;
t587 = m(6) * t593 + mrSges(6,1) * t654 + t624 * t736 + t631 * t676 + t658 * t720 + t707;
t724 = t691 * t585 - t688 * t587;
t704 = -m(5) * t604 - t663 * mrSges(5,2) + t671 * t716 - t724;
t696 = -m(4) * t609 + t669 * t715 + t663 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t662 + (-t668 * t689 - t670 * t692) * qJD(2) + t704;
t574 = m(3) * t611 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t695 + t696;
t729 = t574 * t693;
t728 = t635 * t692;
t605 = -t607 + t728;
t660 = (mrSges(5,2) * t692 - mrSges(5,3) * t689) * qJD(2);
t661 = (-mrSges(4,1) * t692 + mrSges(4,2) * t689) * qJD(2);
t580 = t585 * t688 + t587 * t691;
t603 = -qJDD(3) * pkin(3) + t705 - t728;
t701 = -m(5) * t603 - t662 * mrSges(5,1) - t580;
t576 = m(4) * t605 - mrSges(4,3) * t662 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t669 - t670) * qJD(3) + (-t660 - t661) * t716 + t701;
t598 = pkin(4) * t663 + qJD(3) * t673 - t682 * t737 - t602;
t595 = -0.2e1 * qJD(6) * t658 + (t657 * t676 - t624) * qJ(6) + (t658 * t676 + t623) * pkin(5) + t598;
t588 = m(7) * t595 + t623 * mrSges(7,1) - mrSges(7,3) * t624 + t657 * t632 - t634 * t658;
t699 = m(6) * t598 + mrSges(6,1) * t623 + t624 * mrSges(6,2) + t631 * t657 + t658 * t633 + t588;
t697 = -m(5) * t602 + qJDD(3) * mrSges(5,3) + qJD(3) * t671 + t660 * t715 + t699;
t583 = t697 + (mrSges(4,3) + mrSges(5,1)) * t663 + m(4) * t606 - qJD(3) * t668 - qJDD(3) * mrSges(4,2) + t661 * t715;
t708 = -t576 * t689 + t692 * t583;
t568 = m(3) * t612 - mrSges(3,1) * t695 - qJDD(2) * mrSges(3,2) + t708;
t571 = t692 * t576 + t689 * t583;
t570 = m(3) * t635 + t571;
t559 = t568 * t726 - t570 * t685 + t687 * t729;
t557 = m(2) * t666 + t559;
t564 = t693 * t568 - t574 * t690;
t563 = m(2) * t667 + t564;
t725 = t686 * t557 + t684 * t563;
t723 = t742 * t657 - t734 * t658 - t730 * t676;
t722 = t730 * t657 - t733 * t658 - t740 * t676;
t721 = -t734 * t657 + t744 * t658 + t733 * t676;
t719 = (t741 * qJD(3)) + (t732 * t689 + t731 * t692) * qJD(2);
t718 = -t731 * qJD(3) + (-t735 * t689 - t743 * t692) * qJD(2);
t717 = t732 * qJD(3) + (t745 * t689 + t735 * t692) * qJD(2);
t558 = t568 * t727 + t687 * t570 + t685 * t729;
t709 = -t557 * t684 + t686 * t563;
t577 = -mrSges(5,3) * t662 + t670 * t715 - t704;
t578 = -mrSges(6,1) * t598 - mrSges(7,1) * t595 + mrSges(7,2) * t590 + mrSges(6,3) * t594 - pkin(5) * t588 - t742 * t623 + t734 * t624 + t730 * t654 + t722 * t658 + t721 * t676;
t579 = mrSges(6,2) * t598 + mrSges(7,2) * t591 - mrSges(6,3) * t593 - mrSges(7,3) * t595 - qJ(6) * t588 - t734 * t623 + t744 * t624 + t733 * t654 + t722 * t657 + t723 * t676;
t555 = -mrSges(4,1) * t609 - mrSges(5,1) * t602 + mrSges(5,2) * t604 + mrSges(4,3) * t606 - pkin(3) * t577 + pkin(4) * t699 - pkin(9) * t724 + t717 * qJD(3) + t731 * qJDD(3) - t691 * t578 - t688 * t579 + t735 * t662 + t743 * t663 - t719 * t716;
t560 = pkin(4) * t580 - qJ(4) * t577 + mrSges(7,3) * t590 - mrSges(7,1) * t591 + mrSges(6,1) * t593 - mrSges(6,2) * t594 + mrSges(5,1) * t603 - mrSges(5,3) * t604 - mrSges(4,3) * t605 + mrSges(4,2) * t609 + qJ(6) * t712 + pkin(5) * t707 + t735 * t663 + t745 * t662 + (-pkin(5) * t628 - t723) * t658 + (-qJ(6) * t628 + t721) * t657 + t740 * t654 + (-mrSges(7,2) * pkin(5) + t733) * t624 + (-mrSges(7,2) * qJ(6) - t730) * t623 + t732 * qJDD(3) + t718 * qJD(3) + t719 * t715;
t553 = mrSges(3,2) * t635 - mrSges(3,3) * t611 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t695 - pkin(8) * t571 - t555 * t689 + t560 * t692;
t554 = -pkin(2) * t571 + mrSges(3,3) * t612 - mrSges(3,1) * t635 + Ifges(3,6) * qJDD(2) - pkin(3) * (-qJD(3) * t670 + t701) - qJ(4) * t697 - t691 * t579 + t688 * t578 + pkin(9) * t580 - mrSges(4,1) * t605 + mrSges(4,2) * t606 - mrSges(5,2) * t603 + mrSges(5,3) * t602 + t695 * Ifges(3,5) + (-mrSges(5,1) * qJ(4) - t731) * t663 - t732 * t662 + (mrSges(5,2) * pkin(3) - t741) * qJDD(3) + (t717 * t692 + (pkin(3) * t660 + t718) * t689) * qJD(2);
t702 = pkin(7) * t564 + t553 * t690 + t554 * t693;
t552 = mrSges(3,1) * t611 - mrSges(3,2) * t612 + Ifges(3,3) * qJDD(2) + pkin(2) * t696 + pkin(8) * t708 + t692 * t555 + t689 * t560;
t551 = mrSges(2,2) * t683 - mrSges(2,3) * t666 + t553 * t693 - t554 * t690 + (-t558 * t685 - t559 * t687) * pkin(7);
t550 = -mrSges(2,1) * t683 + mrSges(2,3) * t667 - pkin(1) * t558 - t552 * t685 + t687 * t702;
t1 = [-m(1) * g(1) + t709; -m(1) * g(2) + t725; -m(1) * g(3) + m(2) * t683 + t558; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t725 - t684 * t550 + t686 * t551; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t709 + t686 * t550 + t684 * t551; -mrSges(1,1) * g(2) + mrSges(2,1) * t666 + mrSges(1,2) * g(1) - mrSges(2,2) * t667 + pkin(1) * t559 + t552 * t687 + t685 * t702;];
tauB  = t1;
