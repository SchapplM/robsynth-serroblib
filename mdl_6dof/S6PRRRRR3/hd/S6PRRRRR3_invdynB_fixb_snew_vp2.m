% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRRR3
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 11:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:02:27
% EndTime: 2019-05-05 11:02:53
% DurationCPUTime: 25.40s
% Computational Cost: add. (439561->341), mult. (866769->437), div. (0->0), fcn. (638730->14), ass. (0->143)
t663 = sin(pkin(12));
t665 = cos(pkin(12));
t653 = g(1) * t663 - g(2) * t665;
t654 = -g(1) * t665 - g(2) * t663;
t662 = -g(3) + qJDD(1);
t664 = sin(pkin(6));
t666 = cos(pkin(6));
t671 = sin(qJ(2));
t676 = cos(qJ(2));
t615 = -t671 * t654 + (t653 * t666 + t662 * t664) * t676;
t678 = qJD(2) ^ 2;
t695 = t666 * t671;
t696 = t664 * t671;
t616 = t653 * t695 + t676 * t654 + t662 * t696;
t611 = -pkin(2) * t678 + qJDD(2) * pkin(8) + t616;
t633 = -t653 * t664 + t662 * t666;
t670 = sin(qJ(3));
t675 = cos(qJ(3));
t604 = t675 * t611 + t670 * t633;
t650 = (-pkin(3) * t675 - pkin(9) * t670) * qJD(2);
t677 = qJD(3) ^ 2;
t692 = qJD(2) * t675;
t593 = -pkin(3) * t677 + qJDD(3) * pkin(9) + t650 * t692 + t604;
t610 = -qJDD(2) * pkin(2) - t678 * pkin(8) - t615;
t691 = qJD(2) * qJD(3);
t690 = t675 * t691;
t651 = qJDD(2) * t670 + t690;
t661 = t670 * t691;
t652 = qJDD(2) * t675 - t661;
t598 = (-t651 - t690) * pkin(9) + (-t652 + t661) * pkin(3) + t610;
t669 = sin(qJ(4));
t674 = cos(qJ(4));
t579 = -t669 * t593 + t674 * t598;
t693 = qJD(2) * t670;
t647 = qJD(3) * t674 - t669 * t693;
t624 = qJD(4) * t647 + qJDD(3) * t669 + t651 * t674;
t644 = qJDD(4) - t652;
t648 = qJD(3) * t669 + t674 * t693;
t660 = qJD(4) - t692;
t576 = (t647 * t660 - t624) * pkin(10) + (t647 * t648 + t644) * pkin(4) + t579;
t580 = t674 * t593 + t669 * t598;
t623 = -qJD(4) * t648 + qJDD(3) * t674 - t651 * t669;
t632 = pkin(4) * t660 - pkin(10) * t648;
t643 = t647 ^ 2;
t578 = -pkin(4) * t643 + pkin(10) * t623 - t632 * t660 + t580;
t668 = sin(qJ(5));
t673 = cos(qJ(5));
t566 = t673 * t576 - t668 * t578;
t626 = t647 * t673 - t648 * t668;
t590 = qJD(5) * t626 + t623 * t668 + t624 * t673;
t627 = t647 * t668 + t648 * t673;
t640 = qJDD(5) + t644;
t659 = qJD(5) + t660;
t564 = (t626 * t659 - t590) * pkin(11) + (t626 * t627 + t640) * pkin(5) + t566;
t567 = t668 * t576 + t673 * t578;
t589 = -qJD(5) * t627 + t623 * t673 - t624 * t668;
t614 = pkin(5) * t659 - pkin(11) * t627;
t625 = t626 ^ 2;
t565 = -pkin(5) * t625 + pkin(11) * t589 - t614 * t659 + t567;
t667 = sin(qJ(6));
t672 = cos(qJ(6));
t562 = t564 * t672 - t565 * t667;
t605 = t626 * t672 - t627 * t667;
t573 = qJD(6) * t605 + t589 * t667 + t590 * t672;
t606 = t626 * t667 + t627 * t672;
t587 = -mrSges(7,1) * t605 + mrSges(7,2) * t606;
t655 = qJD(6) + t659;
t596 = -mrSges(7,2) * t655 + mrSges(7,3) * t605;
t635 = qJDD(6) + t640;
t560 = m(7) * t562 + mrSges(7,1) * t635 - mrSges(7,3) * t573 - t587 * t606 + t596 * t655;
t563 = t564 * t667 + t565 * t672;
t572 = -qJD(6) * t606 + t589 * t672 - t590 * t667;
t597 = mrSges(7,1) * t655 - mrSges(7,3) * t606;
t561 = m(7) * t563 - mrSges(7,2) * t635 + mrSges(7,3) * t572 + t587 * t605 - t597 * t655;
t552 = t672 * t560 + t667 * t561;
t607 = -mrSges(6,1) * t626 + mrSges(6,2) * t627;
t612 = -mrSges(6,2) * t659 + mrSges(6,3) * t626;
t550 = m(6) * t566 + mrSges(6,1) * t640 - mrSges(6,3) * t590 - t607 * t627 + t612 * t659 + t552;
t613 = mrSges(6,1) * t659 - mrSges(6,3) * t627;
t685 = -t560 * t667 + t672 * t561;
t551 = m(6) * t567 - mrSges(6,2) * t640 + mrSges(6,3) * t589 + t607 * t626 - t613 * t659 + t685;
t546 = t673 * t550 + t668 * t551;
t628 = -mrSges(5,1) * t647 + mrSges(5,2) * t648;
t630 = -mrSges(5,2) * t660 + mrSges(5,3) * t647;
t544 = m(5) * t579 + mrSges(5,1) * t644 - mrSges(5,3) * t624 - t628 * t648 + t630 * t660 + t546;
t631 = mrSges(5,1) * t660 - mrSges(5,3) * t648;
t686 = -t550 * t668 + t673 * t551;
t545 = m(5) * t580 - mrSges(5,2) * t644 + mrSges(5,3) * t623 + t628 * t647 - t631 * t660 + t686;
t540 = t544 * t674 + t545 * t669;
t656 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t693;
t657 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t692;
t680 = -m(4) * t610 + t652 * mrSges(4,1) - mrSges(4,2) * t651 - t656 * t693 + t657 * t692 - t540;
t536 = m(3) * t615 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t678 + t680;
t697 = t536 * t676;
t649 = (-mrSges(4,1) * t675 + mrSges(4,2) * t670) * qJD(2);
t687 = -t544 * t669 + t674 * t545;
t539 = m(4) * t604 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t652 - qJD(3) * t656 + t649 * t692 + t687;
t603 = -t670 * t611 + t633 * t675;
t592 = -qJDD(3) * pkin(3) - pkin(9) * t677 + t650 * t693 - t603;
t581 = -pkin(4) * t623 - pkin(10) * t643 + t648 * t632 + t592;
t569 = -pkin(5) * t589 - pkin(11) * t625 + t614 * t627 + t581;
t684 = m(7) * t569 - t572 * mrSges(7,1) + t573 * mrSges(7,2) - t605 * t596 + t606 * t597;
t681 = m(6) * t581 - t589 * mrSges(6,1) + t590 * mrSges(6,2) - t626 * t612 + t627 * t613 + t684;
t679 = -m(5) * t592 + t623 * mrSges(5,1) - t624 * mrSges(5,2) + t647 * t630 - t648 * t631 - t681;
t556 = m(4) * t603 + qJDD(3) * mrSges(4,1) - t651 * mrSges(4,3) + qJD(3) * t657 - t649 * t693 + t679;
t688 = t675 * t539 - t556 * t670;
t530 = m(3) * t616 - mrSges(3,1) * t678 - qJDD(2) * mrSges(3,2) + t688;
t533 = t670 * t539 + t675 * t556;
t532 = m(3) * t633 + t533;
t519 = t530 * t695 - t532 * t664 + t666 * t697;
t517 = m(2) * t653 + t519;
t523 = t676 * t530 - t536 * t671;
t522 = m(2) * t654 + t523;
t694 = t665 * t517 + t663 * t522;
t518 = t530 * t696 + t666 * t532 + t664 * t697;
t689 = -t517 * t663 + t665 * t522;
t582 = Ifges(7,5) * t606 + Ifges(7,6) * t605 + Ifges(7,3) * t655;
t584 = Ifges(7,1) * t606 + Ifges(7,4) * t605 + Ifges(7,5) * t655;
t553 = -mrSges(7,1) * t569 + mrSges(7,3) * t563 + Ifges(7,4) * t573 + Ifges(7,2) * t572 + Ifges(7,6) * t635 - t582 * t606 + t584 * t655;
t583 = Ifges(7,4) * t606 + Ifges(7,2) * t605 + Ifges(7,6) * t655;
t554 = mrSges(7,2) * t569 - mrSges(7,3) * t562 + Ifges(7,1) * t573 + Ifges(7,4) * t572 + Ifges(7,5) * t635 + t582 * t605 - t583 * t655;
t599 = Ifges(6,5) * t627 + Ifges(6,6) * t626 + Ifges(6,3) * t659;
t601 = Ifges(6,1) * t627 + Ifges(6,4) * t626 + Ifges(6,5) * t659;
t541 = -mrSges(6,1) * t581 + mrSges(6,3) * t567 + Ifges(6,4) * t590 + Ifges(6,2) * t589 + Ifges(6,6) * t640 - pkin(5) * t684 + pkin(11) * t685 + t672 * t553 + t667 * t554 - t627 * t599 + t659 * t601;
t600 = Ifges(6,4) * t627 + Ifges(6,2) * t626 + Ifges(6,6) * t659;
t542 = mrSges(6,2) * t581 - mrSges(6,3) * t566 + Ifges(6,1) * t590 + Ifges(6,4) * t589 + Ifges(6,5) * t640 - pkin(11) * t552 - t553 * t667 + t554 * t672 + t599 * t626 - t600 * t659;
t617 = Ifges(5,5) * t648 + Ifges(5,6) * t647 + Ifges(5,3) * t660;
t619 = Ifges(5,1) * t648 + Ifges(5,4) * t647 + Ifges(5,5) * t660;
t525 = -mrSges(5,1) * t592 + mrSges(5,3) * t580 + Ifges(5,4) * t624 + Ifges(5,2) * t623 + Ifges(5,6) * t644 - pkin(4) * t681 + pkin(10) * t686 + t673 * t541 + t668 * t542 - t648 * t617 + t660 * t619;
t618 = Ifges(5,4) * t648 + Ifges(5,2) * t647 + Ifges(5,6) * t660;
t526 = mrSges(5,2) * t592 - mrSges(5,3) * t579 + Ifges(5,1) * t624 + Ifges(5,4) * t623 + Ifges(5,5) * t644 - pkin(10) * t546 - t541 * t668 + t542 * t673 + t617 * t647 - t618 * t660;
t637 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t670 + Ifges(4,6) * t675) * qJD(2);
t638 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t670 + Ifges(4,2) * t675) * qJD(2);
t515 = mrSges(4,2) * t610 - mrSges(4,3) * t603 + Ifges(4,1) * t651 + Ifges(4,4) * t652 + Ifges(4,5) * qJDD(3) - pkin(9) * t540 - qJD(3) * t638 - t525 * t669 + t526 * t674 + t637 * t692;
t639 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t670 + Ifges(4,4) * t675) * qJD(2);
t524 = -pkin(5) * t552 - pkin(4) * t546 - pkin(3) * t540 - t637 * t693 + Ifges(4,6) * qJDD(3) - Ifges(5,3) * t644 + t647 * t619 - t648 * t618 + Ifges(4,4) * t651 + Ifges(4,2) * t652 - Ifges(7,3) * t635 + qJD(3) * t639 - Ifges(6,3) * t640 - Ifges(5,6) * t623 - Ifges(5,5) * t624 + t626 * t601 - t627 * t600 - t606 * t583 - mrSges(4,1) * t610 + mrSges(4,3) * t604 + t605 * t584 - Ifges(6,6) * t589 - Ifges(6,5) * t590 - mrSges(5,1) * t579 + mrSges(5,2) * t580 - Ifges(7,6) * t572 - Ifges(7,5) * t573 - mrSges(6,1) * t566 + mrSges(6,2) * t567 + mrSges(7,2) * t563 - mrSges(7,1) * t562;
t513 = mrSges(3,2) * t633 - mrSges(3,3) * t615 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t678 - pkin(8) * t533 + t515 * t675 - t524 * t670;
t514 = Ifges(3,6) * qJDD(2) + t678 * Ifges(3,5) - mrSges(3,1) * t633 + mrSges(3,3) * t616 - Ifges(4,5) * t651 - Ifges(4,6) * t652 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t603 + mrSges(4,2) * t604 - t669 * t526 - t674 * t525 - pkin(3) * t679 - pkin(9) * t687 - pkin(2) * t533 + (-t638 * t670 + t639 * t675) * qJD(2);
t682 = pkin(7) * t523 + t513 * t671 + t514 * t676;
t512 = mrSges(3,1) * t615 - mrSges(3,2) * t616 + Ifges(3,3) * qJDD(2) + pkin(2) * t680 + pkin(8) * t688 + t670 * t515 + t675 * t524;
t511 = mrSges(2,2) * t662 - mrSges(2,3) * t653 + t676 * t513 - t671 * t514 + (-t518 * t664 - t519 * t666) * pkin(7);
t510 = -mrSges(2,1) * t662 + mrSges(2,3) * t654 - pkin(1) * t518 - t664 * t512 + t666 * t682;
t1 = [-m(1) * g(1) + t689; -m(1) * g(2) + t694; -m(1) * g(3) + m(2) * t662 + t518; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t694 - t663 * t510 + t665 * t511; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t689 + t665 * t510 + t663 * t511; -mrSges(1,1) * g(2) + mrSges(2,1) * t653 + mrSges(1,2) * g(1) - mrSges(2,2) * t654 + pkin(1) * t519 + t666 * t512 + t664 * t682;];
tauB  = t1;
