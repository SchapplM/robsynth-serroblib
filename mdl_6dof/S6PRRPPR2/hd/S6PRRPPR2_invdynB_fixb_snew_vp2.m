% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-05-05 02:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:44:32
% EndTime: 2019-05-05 02:44:43
% DurationCPUTime: 7.51s
% Computational Cost: add. (110871->324), mult. (240336->400), div. (0->0), fcn. (163916->12), ass. (0->137)
t729 = -2 * qJD(4);
t728 = Ifges(5,1) + Ifges(6,2);
t727 = -Ifges(6,1) - Ifges(5,3);
t724 = Ifges(5,4) + Ifges(6,6);
t723 = Ifges(5,5) - Ifges(6,4);
t726 = Ifges(5,2) + Ifges(6,3);
t722 = Ifges(5,6) - Ifges(6,5);
t680 = sin(pkin(10));
t682 = cos(pkin(10));
t669 = g(1) * t680 - g(2) * t682;
t670 = -g(1) * t682 - g(2) * t680;
t678 = -g(3) + qJDD(1);
t689 = cos(qJ(2));
t683 = cos(pkin(6));
t686 = sin(qJ(2));
t718 = t683 * t686;
t681 = sin(pkin(6));
t719 = t681 * t686;
t623 = t669 * t718 + t689 * t670 + t678 * t719;
t691 = qJD(2) ^ 2;
t615 = -pkin(2) * t691 + qJDD(2) * pkin(8) + t623;
t647 = -t669 * t681 + t678 * t683;
t685 = sin(qJ(3));
t688 = cos(qJ(3));
t599 = -t615 * t685 + t688 * t647;
t707 = qJD(2) * qJD(3);
t706 = t688 * t707;
t667 = qJDD(2) * t685 + t706;
t596 = (-t667 + t706) * qJ(4) + (t685 * t688 * t691 + qJDD(3)) * pkin(3) + t599;
t600 = t688 * t615 + t685 * t647;
t668 = qJDD(2) * t688 - t685 * t707;
t711 = qJD(2) * t685;
t671 = qJD(3) * pkin(3) - qJ(4) * t711;
t677 = t688 ^ 2;
t597 = -pkin(3) * t677 * t691 + qJ(4) * t668 - qJD(3) * t671 + t600;
t679 = sin(pkin(11));
t721 = cos(pkin(11));
t655 = (t679 * t688 + t685 * t721) * qJD(2);
t589 = t596 * t721 - t679 * t597 + t655 * t729;
t622 = -t686 * t670 + (t669 * t683 + t678 * t681) * t689;
t725 = -2 * qJD(5);
t696 = -qJDD(2) * pkin(2) - t622;
t614 = -pkin(8) * t691 + t696;
t672 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t711;
t710 = qJD(2) * t688;
t673 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t710;
t635 = t667 * t679 - t668 * t721;
t598 = -pkin(3) * t668 + qJDD(4) + t671 * t711 + (-qJ(4) * t677 - pkin(8)) * t691 + t696;
t636 = t667 * t721 + t679 * t668;
t654 = t679 * t711 - t710 * t721;
t709 = qJD(3) * t654;
t692 = (-t636 + t709) * qJ(5) + t598 + (qJD(3) * pkin(4) + t725) * t655;
t592 = pkin(4) * t635 + t692;
t644 = mrSges(6,1) * t654 - qJD(3) * mrSges(6,3);
t645 = mrSges(6,1) * t655 + qJD(3) * mrSges(6,2);
t627 = pkin(4) * t654 - qJ(5) * t655;
t690 = qJD(3) ^ 2;
t587 = -qJDD(3) * pkin(4) - t690 * qJ(5) + t655 * t627 + qJDD(5) - t589;
t583 = (t654 * t655 - qJDD(3)) * pkin(9) + (t636 + t709) * pkin(5) + t587;
t646 = pkin(5) * t655 - qJD(3) * pkin(9);
t653 = t654 ^ 2;
t588 = -pkin(5) * t653 - t646 * t655 + (pkin(4) + pkin(9)) * t635 + t692;
t684 = sin(qJ(6));
t687 = cos(qJ(6));
t581 = t583 * t687 - t588 * t684;
t637 = -qJD(3) * t684 + t654 * t687;
t607 = qJD(6) * t637 + qJDD(3) * t687 + t635 * t684;
t638 = qJD(3) * t687 + t654 * t684;
t608 = -mrSges(7,1) * t637 + mrSges(7,2) * t638;
t652 = qJD(6) + t655;
t612 = -mrSges(7,2) * t652 + mrSges(7,3) * t637;
t634 = qJDD(6) + t636;
t579 = m(7) * t581 + mrSges(7,1) * t634 - mrSges(7,3) * t607 - t608 * t638 + t612 * t652;
t582 = t583 * t684 + t588 * t687;
t606 = -qJD(6) * t638 - qJDD(3) * t684 + t635 * t687;
t613 = mrSges(7,1) * t652 - mrSges(7,3) * t638;
t580 = m(7) * t582 - mrSges(7,2) * t634 + mrSges(7,3) * t606 + t608 * t637 - t613 * t652;
t716 = -t684 * t579 + t687 * t580;
t570 = m(6) * t592 - t635 * mrSges(6,2) - t636 * mrSges(6,3) - t654 * t644 - t655 * t645 + t716;
t642 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t654;
t643 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t655;
t694 = m(5) * t598 + t635 * mrSges(5,1) + mrSges(5,2) * t636 + t654 * t642 + t643 * t655 + t570;
t693 = -m(4) * t614 + t668 * mrSges(4,1) - mrSges(4,2) * t667 - t672 * t711 + t673 * t710 - t694;
t569 = m(3) * t622 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t691 + t693;
t720 = t569 * t689;
t628 = mrSges(5,1) * t654 + mrSges(5,2) * t655;
t571 = t579 * t687 + t580 * t684;
t629 = -mrSges(6,2) * t654 - mrSges(6,3) * t655;
t697 = -m(6) * t587 - t636 * mrSges(6,1) - t655 * t629 - t571;
t566 = m(5) * t589 - mrSges(5,3) * t636 - t628 * t655 + (mrSges(5,1) - mrSges(6,2)) * qJDD(3) + (t642 - t644) * qJD(3) + t697;
t650 = t654 * t729;
t715 = t679 * t596 + t721 * t597;
t590 = t650 + t715;
t700 = pkin(4) * t690 - qJDD(3) * qJ(5) - t715;
t586 = qJD(3) * t725 + ((2 * qJD(4)) + t627) * t654 + t700;
t585 = -pkin(5) * t635 - pkin(9) * t653 - t627 * t654 + t650 + ((2 * qJD(5)) + t646) * qJD(3) - t700;
t698 = -m(7) * t585 + mrSges(7,1) * t606 - t607 * mrSges(7,2) + t612 * t637 - t638 * t613;
t695 = -m(6) * t586 + qJDD(3) * mrSges(6,3) + qJD(3) * t645 - t698;
t576 = m(5) * t590 - qJDD(3) * mrSges(5,2) - qJD(3) * t643 + (-t628 - t629) * t654 + (-mrSges(5,3) - mrSges(6,1)) * t635 + t695;
t564 = t721 * t566 + t679 * t576;
t666 = (-mrSges(4,1) * t688 + mrSges(4,2) * t685) * qJD(2);
t562 = m(4) * t599 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t667 + qJD(3) * t673 - t666 * t711 + t564;
t703 = -t566 * t679 + t721 * t576;
t563 = m(4) * t600 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t668 - qJD(3) * t672 + t666 * t710 + t703;
t704 = -t562 * t685 + t688 * t563;
t554 = m(3) * t623 - mrSges(3,1) * t691 - qJDD(2) * mrSges(3,2) + t704;
t557 = t688 * t562 + t685 * t563;
t556 = m(3) * t647 + t557;
t545 = t554 * t718 - t556 * t681 + t683 * t720;
t543 = m(2) * t669 + t545;
t549 = t689 * t554 - t569 * t686;
t548 = m(2) * t670 + t549;
t717 = t682 * t543 + t680 * t548;
t714 = t727 * qJD(3) + t722 * t654 - t723 * t655;
t713 = -t722 * qJD(3) + t726 * t654 - t724 * t655;
t712 = t723 * qJD(3) - t724 * t654 + t728 * t655;
t544 = t554 * t719 + t683 * t556 + t681 * t720;
t705 = -t543 * t680 + t682 * t548;
t601 = Ifges(7,5) * t638 + Ifges(7,6) * t637 + Ifges(7,3) * t652;
t603 = Ifges(7,1) * t638 + Ifges(7,4) * t637 + Ifges(7,5) * t652;
t574 = -mrSges(7,1) * t585 + mrSges(7,3) * t582 + Ifges(7,4) * t607 + Ifges(7,2) * t606 + Ifges(7,6) * t634 - t601 * t638 + t603 * t652;
t602 = Ifges(7,4) * t638 + Ifges(7,2) * t637 + Ifges(7,6) * t652;
t575 = mrSges(7,2) * t585 - mrSges(7,3) * t581 + Ifges(7,1) * t607 + Ifges(7,4) * t606 + Ifges(7,5) * t634 + t601 * t637 - t602 * t652;
t550 = -mrSges(5,1) * t598 - mrSges(6,1) * t586 + mrSges(6,2) * t592 + mrSges(5,3) * t590 - pkin(4) * t570 - pkin(5) * t698 - pkin(9) * t716 + t712 * qJD(3) + t722 * qJDD(3) - t687 * t574 - t684 * t575 - t726 * t635 + t724 * t636 + t714 * t655;
t558 = mrSges(6,1) * t587 + mrSges(7,1) * t581 + mrSges(5,2) * t598 - mrSges(7,2) * t582 - mrSges(5,3) * t589 - mrSges(6,3) * t592 + Ifges(7,5) * t607 + Ifges(7,6) * t606 + Ifges(7,3) * t634 + pkin(5) * t571 - qJ(5) * t570 + t638 * t602 - t637 * t603 + t714 * t654 + t728 * t636 - t724 * t635 + t723 * qJDD(3) + t713 * qJD(3);
t657 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t685 + Ifges(4,6) * t688) * qJD(2);
t659 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t685 + Ifges(4,4) * t688) * qJD(2);
t539 = -mrSges(4,1) * t614 + mrSges(4,3) * t600 + Ifges(4,4) * t667 + Ifges(4,2) * t668 + Ifges(4,6) * qJDD(3) - pkin(3) * t694 + qJ(4) * t703 + qJD(3) * t659 + t550 * t721 + t679 * t558 - t657 * t711;
t658 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t685 + Ifges(4,2) * t688) * qJD(2);
t541 = mrSges(4,2) * t614 - mrSges(4,3) * t599 + Ifges(4,1) * t667 + Ifges(4,4) * t668 + Ifges(4,5) * qJDD(3) - qJ(4) * t564 - qJD(3) * t658 - t679 * t550 + t558 * t721 + t657 * t710;
t538 = mrSges(3,2) * t647 - mrSges(3,3) * t622 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t691 - pkin(8) * t557 - t539 * t685 + t541 * t688;
t540 = (qJ(5) * mrSges(6,1) + t722) * t635 - t723 * t636 + (qJ(5) * t629 - t712) * t654 + t713 * t655 + Ifges(3,6) * qJDD(2) + (-t685 * t658 + t688 * t659) * qJD(2) - qJ(5) * t695 - pkin(4) * (-qJD(3) * t644 + t697) + (pkin(4) * mrSges(6,2) - Ifges(4,3) + t727) * qJDD(3) + t691 * Ifges(3,5) + t684 * t574 - t687 * t575 - Ifges(4,5) * t667 - Ifges(4,6) * t668 - mrSges(3,1) * t647 + mrSges(3,3) * t623 - mrSges(4,1) * t599 + mrSges(4,2) * t600 + mrSges(5,2) * t590 + mrSges(6,3) * t586 - mrSges(6,2) * t587 - mrSges(5,1) * t589 + pkin(9) * t571 - pkin(3) * t564 - pkin(2) * t557;
t699 = pkin(7) * t549 + t538 * t686 + t540 * t689;
t537 = mrSges(3,1) * t622 - mrSges(3,2) * t623 + Ifges(3,3) * qJDD(2) + pkin(2) * t693 + pkin(8) * t704 + t688 * t539 + t685 * t541;
t536 = mrSges(2,2) * t678 - mrSges(2,3) * t669 + t538 * t689 - t540 * t686 + (-t544 * t681 - t545 * t683) * pkin(7);
t535 = -mrSges(2,1) * t678 + mrSges(2,3) * t670 - pkin(1) * t544 - t537 * t681 + t683 * t699;
t1 = [-m(1) * g(1) + t705; -m(1) * g(2) + t717; -m(1) * g(3) + m(2) * t678 + t544; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t717 - t680 * t535 + t682 * t536; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t705 + t682 * t535 + t680 * t536; -mrSges(1,1) * g(2) + mrSges(2,1) * t669 + mrSges(1,2) * g(1) - mrSges(2,2) * t670 + pkin(1) * t545 + t537 * t683 + t681 * t699;];
tauB  = t1;
