% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-05-05 06:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRPP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:50:07
% EndTime: 2019-05-05 06:50:15
% DurationCPUTime: 5.52s
% Computational Cost: add. (61943->297), mult. (116765->353), div. (0->0), fcn. (75705->10), ass. (0->124)
t694 = Ifges(5,1) + Ifges(6,2) + Ifges(7,3);
t677 = Ifges(5,4) + Ifges(6,6) - Ifges(7,6);
t676 = Ifges(5,5) - Ifges(6,4) + Ifges(7,5);
t693 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t692 = -Ifges(5,6) + Ifges(6,5) + Ifges(7,4);
t691 = Ifges(5,3) + Ifges(6,1) + Ifges(7,1);
t644 = sin(pkin(10));
t646 = cos(pkin(10));
t633 = t644 * g(1) - t646 * g(2);
t634 = -t646 * g(1) - t644 * g(2);
t643 = -g(3) + qJDD(1);
t645 = sin(pkin(6));
t647 = cos(pkin(6));
t650 = sin(qJ(2));
t652 = cos(qJ(2));
t573 = -t650 * t634 + (t633 * t647 + t643 * t645) * t652;
t690 = -2 * qJD(5);
t689 = 2 * qJD(6);
t688 = cos(qJ(4));
t687 = -mrSges(7,1) - mrSges(5,3);
t654 = qJD(2) ^ 2;
t683 = t647 * t650;
t684 = t645 * t650;
t574 = t633 * t683 + t652 * t634 + t643 * t684;
t570 = -t654 * pkin(2) + qJDD(2) * pkin(8) + t574;
t611 = -t645 * t633 + t647 * t643;
t649 = sin(qJ(3));
t651 = cos(qJ(3));
t566 = t651 * t570 + t649 * t611;
t630 = (-pkin(3) * t651 - pkin(9) * t649) * qJD(2);
t653 = qJD(3) ^ 2;
t679 = t651 * qJD(2);
t562 = -t653 * pkin(3) + qJDD(3) * pkin(9) + t630 * t679 + t566;
t569 = -qJDD(2) * pkin(2) - t654 * pkin(8) - t573;
t678 = qJD(2) * qJD(3);
t669 = t651 * t678;
t631 = t649 * qJDD(2) + t669;
t670 = t649 * t678;
t632 = t651 * qJDD(2) - t670;
t564 = (-t631 - t669) * pkin(9) + (-t632 + t670) * pkin(3) + t569;
t648 = sin(qJ(4));
t557 = -t648 * t562 + t688 * t564;
t680 = qJD(2) * t649;
t627 = -t688 * qJD(3) + t648 * t680;
t594 = -t627 * qJD(4) + t648 * qJDD(3) + t688 * t631;
t628 = t648 * qJD(3) + t688 * t680;
t598 = -t628 * mrSges(7,2) + t627 * mrSges(7,3);
t600 = t627 * mrSges(5,1) + t628 * mrSges(5,2);
t639 = -qJD(4) + t679;
t603 = t639 * mrSges(5,2) - t627 * mrSges(5,3);
t607 = t627 * mrSges(6,1) + t639 * mrSges(6,3);
t624 = qJDD(4) - t632;
t599 = t627 * pkin(4) - t628 * qJ(5);
t638 = t639 ^ 2;
t555 = -t624 * pkin(4) - t638 * qJ(5) + t628 * t599 + qJDD(5) - t557;
t601 = -t627 * mrSges(6,2) - t628 * mrSges(6,3);
t685 = t627 * t639;
t549 = t639 * t689 + (t627 * t628 - t624) * qJ(6) + (t594 - t685) * pkin(5) + t555;
t608 = -t627 * mrSges(7,1) - t639 * mrSges(7,2);
t665 = m(7) * t549 - t624 * mrSges(7,3) + t639 * t608;
t660 = -m(6) * t555 - t594 * mrSges(6,1) - t628 * t601 - t665;
t544 = m(5) * t557 + (-t603 + t607) * t639 + (-t598 - t600) * t628 + (mrSges(5,1) - mrSges(6,2)) * t624 + t687 * t594 + t660;
t558 = t688 * t562 + t648 * t564;
t593 = t628 * qJD(4) - t688 * qJDD(3) + t648 * t631;
t604 = -t639 * mrSges(5,1) - t628 * mrSges(5,3);
t658 = -t638 * pkin(4) + t624 * qJ(5) - t627 * t599 + t558;
t554 = 0.2e1 * qJD(5) * t639 - t658;
t609 = t628 * mrSges(6,1) - t639 * mrSges(6,2);
t605 = t628 * pkin(5) + t639 * qJ(6);
t623 = t627 ^ 2;
t551 = -t593 * pkin(5) - t623 * qJ(6) + qJDD(6) + (t690 - t605) * t639 + t658;
t606 = t628 * mrSges(7,1) + t639 * mrSges(7,3);
t671 = -m(7) * t551 - t624 * mrSges(7,2) + t639 * t606;
t662 = -m(6) * t554 + t624 * mrSges(6,3) - t639 * t609 - t671;
t681 = -t598 - t601;
t546 = m(5) * t558 - t624 * mrSges(5,2) + t639 * t604 + (-t600 + t681) * t627 + (-mrSges(6,1) + t687) * t593 + t662;
t541 = t688 * t544 + t648 * t546;
t635 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t680;
t636 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t679;
t657 = -m(4) * t569 + t632 * mrSges(4,1) - t631 * mrSges(4,2) - t635 * t680 + t636 * t679 - t541;
t537 = m(3) * t573 + qJDD(2) * mrSges(3,1) - t654 * mrSges(3,2) + t657;
t686 = t537 * t652;
t629 = (-mrSges(4,1) * t651 + mrSges(4,2) * t649) * qJD(2);
t666 = -t648 * t544 + t688 * t546;
t540 = m(4) * t566 - qJDD(3) * mrSges(4,2) + t632 * mrSges(4,3) - qJD(3) * t635 + t629 * t679 + t666;
t565 = -t649 * t570 + t651 * t611;
t561 = -qJDD(3) * pkin(3) - t653 * pkin(9) + t630 * t680 - t565;
t656 = (-t594 - t685) * qJ(5) + t561 + (-t639 * pkin(4) + t690) * t628;
t556 = t593 * pkin(4) + t656;
t553 = -t623 * pkin(5) + t627 * t689 - t628 * t605 + (pkin(4) + qJ(6)) * t593 + t656;
t664 = m(7) * t553 - t594 * mrSges(7,2) + t593 * mrSges(7,3) - t628 * t606 + t627 * t608;
t659 = -m(6) * t556 + t593 * mrSges(6,2) + t627 * t607 - t664;
t655 = -m(5) * t561 - t593 * mrSges(5,1) - t627 * t603 + (-t604 + t609) * t628 + (-mrSges(5,2) + mrSges(6,3)) * t594 + t659;
t543 = m(4) * t565 + qJDD(3) * mrSges(4,1) - t631 * mrSges(4,3) + qJD(3) * t636 - t629 * t680 + t655;
t667 = t651 * t540 - t649 * t543;
t530 = m(3) * t574 - t654 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t667;
t533 = t649 * t540 + t651 * t543;
t532 = m(3) * t611 + t533;
t519 = t530 * t683 - t645 * t532 + t647 * t686;
t517 = m(2) * t633 + t519;
t525 = t652 * t530 - t650 * t537;
t524 = m(2) * t634 + t525;
t682 = t646 * t517 + t644 * t524;
t518 = t530 * t684 + t647 * t532 + t645 * t686;
t674 = -t692 * t627 - t676 * t628 + t691 * t639;
t673 = t693 * t627 + t677 * t628 + t692 * t639;
t672 = t677 * t627 - t694 * t628 + t676 * t639;
t668 = -t644 * t517 + t646 * t524;
t547 = -t594 * mrSges(6,3) - t628 * t609 - t659;
t526 = -mrSges(5,1) * t561 + mrSges(5,3) * t558 - mrSges(6,1) * t554 + mrSges(6,2) * t556 + mrSges(7,1) * t551 - mrSges(7,3) * t553 - pkin(5) * (t627 * t598 + t671) - qJ(6) * t664 - pkin(4) * t547 + t672 * t639 + t674 * t628 - t692 * t624 + t677 * t594 + (-pkin(5) * mrSges(7,1) + t693) * t593;
t548 = t594 * mrSges(7,1) + t628 * t598 + t665;
t534 = mrSges(6,1) * t555 + mrSges(7,1) * t549 + mrSges(5,2) * t561 - mrSges(7,2) * t553 - mrSges(5,3) * t557 - mrSges(6,3) * t556 + pkin(5) * t548 - qJ(5) * t547 - t677 * t593 + t694 * t594 + t676 * t624 + t674 * t627 + t673 * t639;
t615 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t649 + Ifges(4,6) * t651) * qJD(2);
t616 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t649 + Ifges(4,2) * t651) * qJD(2);
t520 = mrSges(4,2) * t569 - mrSges(4,3) * t565 + Ifges(4,1) * t631 + Ifges(4,4) * t632 + Ifges(4,5) * qJDD(3) - pkin(9) * t541 - qJD(3) * t616 - t648 * t526 + t688 * t534 + t615 * t679;
t617 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t649 + Ifges(4,4) * t651) * qJD(2);
t521 = -t615 * t680 + mrSges(6,3) * t554 - mrSges(6,2) * t555 - mrSges(5,1) * t557 + mrSges(5,2) * t558 + qJ(6) * t548 + mrSges(7,3) * t549 - mrSges(7,2) * t551 - pkin(3) * t541 + Ifges(4,4) * t631 - pkin(4) * (t639 * t607 + t660) + mrSges(4,3) * t566 - mrSges(4,1) * t569 - qJ(5) * t662 + qJD(3) * t617 + Ifges(4,6) * qJDD(3) + Ifges(4,2) * t632 + (pkin(4) * t598 - t673) * t628 + (pkin(4) * mrSges(6,2) - t691) * t624 + (pkin(4) * mrSges(7,1) - t676) * t594 + (-qJ(5) * t681 + t672) * t627 + (-qJ(5) * (-mrSges(6,1) - mrSges(7,1)) - t692) * t593;
t514 = mrSges(3,2) * t611 - mrSges(3,3) * t573 + Ifges(3,5) * qJDD(2) - t654 * Ifges(3,6) - pkin(8) * t533 + t651 * t520 - t649 * t521;
t515 = Ifges(3,6) * qJDD(2) + t654 * Ifges(3,5) - mrSges(3,1) * t611 + mrSges(3,3) * t574 - Ifges(4,5) * t631 - Ifges(4,6) * t632 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t565 + mrSges(4,2) * t566 - t648 * t534 - t688 * t526 - pkin(3) * t655 - pkin(9) * t666 - pkin(2) * t533 + (-t649 * t616 + t651 * t617) * qJD(2);
t661 = pkin(7) * t525 + t514 * t650 + t515 * t652;
t513 = mrSges(3,1) * t573 - mrSges(3,2) * t574 + Ifges(3,3) * qJDD(2) + pkin(2) * t657 + pkin(8) * t667 + t649 * t520 + t651 * t521;
t512 = mrSges(2,2) * t643 - mrSges(2,3) * t633 + t652 * t514 - t650 * t515 + (-t518 * t645 - t519 * t647) * pkin(7);
t511 = -mrSges(2,1) * t643 + mrSges(2,3) * t634 - pkin(1) * t518 - t645 * t513 + t661 * t647;
t1 = [-m(1) * g(1) + t668; -m(1) * g(2) + t682; -m(1) * g(3) + m(2) * t643 + t518; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t682 - t644 * t511 + t646 * t512; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t668 + t646 * t511 + t644 * t512; -mrSges(1,1) * g(2) + mrSges(2,1) * t633 + mrSges(1,2) * g(1) - mrSges(2,2) * t634 + pkin(1) * t519 + t647 * t513 + t661 * t645;];
tauB  = t1;
