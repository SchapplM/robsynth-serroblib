% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-05-05 17:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:32:58
% EndTime: 2019-05-05 17:33:08
% DurationCPUTime: 7.51s
% Computational Cost: add. (83966->319), mult. (179663->391), div. (0->0), fcn. (116279->10), ass. (0->125)
t683 = -2 * qJD(4);
t682 = Ifges(6,1) + Ifges(7,1);
t676 = Ifges(6,4) - Ifges(7,5);
t681 = -Ifges(6,5) - Ifges(7,4);
t680 = Ifges(6,2) + Ifges(7,3);
t674 = Ifges(6,6) - Ifges(7,6);
t679 = -Ifges(6,3) - Ifges(7,2);
t645 = sin(qJ(1));
t647 = cos(qJ(1));
t629 = t645 * g(1) - g(2) * t647;
t621 = qJDD(1) * pkin(1) + t629;
t630 = -g(1) * t647 - g(2) * t645;
t649 = qJD(1) ^ 2;
t623 = -pkin(1) * t649 + t630;
t640 = sin(pkin(9));
t642 = cos(pkin(9));
t599 = t640 * t621 + t642 * t623;
t590 = -pkin(2) * t649 + qJDD(1) * pkin(7) + t599;
t638 = -g(3) + qJDD(2);
t644 = sin(qJ(3));
t646 = cos(qJ(3));
t579 = -t644 * t590 + t646 * t638;
t665 = qJD(1) * qJD(3);
t662 = t646 * t665;
t624 = qJDD(1) * t644 + t662;
t558 = (-t624 + t662) * qJ(4) + (t644 * t646 * t649 + qJDD(3)) * pkin(3) + t579;
t580 = t646 * t590 + t644 * t638;
t625 = qJDD(1) * t646 - t644 * t665;
t668 = qJD(1) * t644;
t626 = qJD(3) * pkin(3) - qJ(4) * t668;
t637 = t646 ^ 2;
t559 = -pkin(3) * t637 * t649 + qJ(4) * t625 - qJD(3) * t626 + t580;
t639 = sin(pkin(10));
t641 = cos(pkin(10));
t611 = (t639 * t646 + t641 * t644) * qJD(1);
t551 = t641 * t558 - t639 * t559 + t611 * t683;
t610 = (t639 * t644 - t641 * t646) * qJD(1);
t678 = cos(qJ(5));
t677 = -mrSges(6,3) - mrSges(7,2);
t552 = t639 * t558 + t641 * t559 + t610 * t683;
t592 = mrSges(5,1) * t610 + mrSges(5,2) * t611;
t600 = -t624 * t639 + t625 * t641;
t605 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t611;
t593 = pkin(4) * t610 - pkin(8) * t611;
t648 = qJD(3) ^ 2;
t550 = -pkin(4) * t648 + qJDD(3) * pkin(8) - t593 * t610 + t552;
t598 = t642 * t621 - t640 * t623;
t653 = -qJDD(1) * pkin(2) - t598;
t560 = -t625 * pkin(3) + qJDD(4) + t626 * t668 + (-qJ(4) * t637 - pkin(7)) * t649 + t653;
t601 = t624 * t641 + t625 * t639;
t554 = (qJD(3) * t610 - t601) * pkin(8) + (qJD(3) * t611 - t600) * pkin(4) + t560;
t643 = sin(qJ(5));
t547 = t678 * t550 + t643 * t554;
t603 = t643 * qJD(3) + t611 * t678;
t571 = t603 * qJD(5) - qJDD(3) * t678 + t643 * t601;
t609 = qJD(5) + t610;
t583 = mrSges(6,1) * t609 - mrSges(6,3) * t603;
t597 = qJDD(5) - t600;
t602 = -qJD(3) * t678 + t643 * t611;
t575 = pkin(5) * t602 - qJ(6) * t603;
t608 = t609 ^ 2;
t543 = -pkin(5) * t608 + qJ(6) * t597 + 0.2e1 * qJD(6) * t609 - t575 * t602 + t547;
t584 = -mrSges(7,1) * t609 + mrSges(7,2) * t603;
t664 = m(7) * t543 + t597 * mrSges(7,3) + t609 * t584;
t576 = mrSges(7,1) * t602 - mrSges(7,3) * t603;
t669 = -mrSges(6,1) * t602 - mrSges(6,2) * t603 - t576;
t538 = m(6) * t547 - t597 * mrSges(6,2) + t571 * t677 - t609 * t583 + t602 * t669 + t664;
t546 = -t643 * t550 + t554 * t678;
t572 = -t602 * qJD(5) + t643 * qJDD(3) + t601 * t678;
t582 = -mrSges(6,2) * t609 - mrSges(6,3) * t602;
t544 = -t597 * pkin(5) - t608 * qJ(6) + t603 * t575 + qJDD(6) - t546;
t581 = -mrSges(7,2) * t602 + mrSges(7,3) * t609;
t655 = -m(7) * t544 + t597 * mrSges(7,1) + t609 * t581;
t540 = m(6) * t546 + t597 * mrSges(6,1) + t572 * t677 + t609 * t582 + t603 * t669 + t655;
t657 = t678 * t538 - t540 * t643;
t530 = m(5) * t552 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t600 - qJD(3) * t605 - t592 * t610 + t657;
t604 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t610;
t549 = -qJDD(3) * pkin(4) - t648 * pkin(8) + t611 * t593 - t551;
t545 = -0.2e1 * qJD(6) * t603 + (t602 * t609 - t572) * qJ(6) + (t603 * t609 + t571) * pkin(5) + t549;
t541 = m(7) * t545 + mrSges(7,1) * t571 - t572 * mrSges(7,3) + t581 * t602 - t603 * t584;
t651 = -m(6) * t549 - t571 * mrSges(6,1) - mrSges(6,2) * t572 - t602 * t582 - t583 * t603 - t541;
t535 = m(5) * t551 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t601 + qJD(3) * t604 - t592 * t611 + t651;
t525 = t639 * t530 + t641 * t535;
t622 = (-mrSges(4,1) * t646 + mrSges(4,2) * t644) * qJD(1);
t667 = qJD(1) * t646;
t628 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t667;
t523 = m(4) * t579 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t624 + qJD(3) * t628 - t622 * t668 + t525;
t627 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t668;
t658 = t641 * t530 - t535 * t639;
t524 = m(4) * t580 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t625 - qJD(3) * t627 + t622 * t667 + t658;
t659 = -t523 * t644 + t646 * t524;
t516 = m(3) * t599 - mrSges(3,1) * t649 - qJDD(1) * mrSges(3,2) + t659;
t589 = -t649 * pkin(7) + t653;
t533 = t643 * t538 + t678 * t540;
t652 = m(5) * t560 - t600 * mrSges(5,1) + mrSges(5,2) * t601 + t610 * t604 + t605 * t611 + t533;
t650 = -m(4) * t589 + t625 * mrSges(4,1) - mrSges(4,2) * t624 - t627 * t668 + t628 * t667 - t652;
t527 = m(3) * t598 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t649 + t650;
t513 = t640 * t516 + t642 * t527;
t511 = m(2) * t629 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t649 + t513;
t660 = t642 * t516 - t527 * t640;
t512 = m(2) * t630 - mrSges(2,1) * t649 - qJDD(1) * mrSges(2,2) + t660;
t673 = t647 * t511 + t645 * t512;
t517 = t646 * t523 + t644 * t524;
t672 = t680 * t602 - t676 * t603 - t674 * t609;
t671 = t674 * t602 + t681 * t603 + t679 * t609;
t670 = -t676 * t602 + t682 * t603 - t681 * t609;
t663 = m(3) * t638 + t517;
t661 = -t511 * t645 + t647 * t512;
t617 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t644 + Ifges(4,4) * t646) * qJD(1);
t616 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t644 + Ifges(4,2) * t646) * qJD(1);
t615 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t644 + Ifges(4,6) * t646) * qJD(1);
t588 = Ifges(5,1) * t611 - Ifges(5,4) * t610 + Ifges(5,5) * qJD(3);
t587 = Ifges(5,4) * t611 - Ifges(5,2) * t610 + Ifges(5,6) * qJD(3);
t586 = Ifges(5,5) * t611 - Ifges(5,6) * t610 + Ifges(5,3) * qJD(3);
t532 = mrSges(6,2) * t549 + mrSges(7,2) * t544 - mrSges(6,3) * t546 - mrSges(7,3) * t545 - qJ(6) * t541 - t676 * t571 + t682 * t572 - t597 * t681 + t671 * t602 + t672 * t609;
t531 = -mrSges(6,1) * t549 - mrSges(7,1) * t545 + mrSges(7,2) * t543 + mrSges(6,3) * t547 - pkin(5) * t541 - t680 * t571 + t676 * t572 + t674 * t597 + t671 * t603 + t670 * t609;
t519 = Ifges(5,4) * t601 + Ifges(5,2) * t600 + Ifges(5,6) * qJDD(3) - t611 * t586 + qJD(3) * t588 - mrSges(5,1) * t560 + mrSges(5,3) * t552 - mrSges(6,1) * t546 + mrSges(6,2) * t547 + mrSges(7,1) * t544 - mrSges(7,3) * t543 - pkin(5) * t655 - qJ(6) * t664 - pkin(4) * t533 + (pkin(5) * t576 + t672) * t603 + (qJ(6) * t576 - t670) * t602 + t679 * t597 + (mrSges(7,2) * pkin(5) + t681) * t572 + (mrSges(7,2) * qJ(6) + t674) * t571;
t518 = mrSges(5,2) * t560 - mrSges(5,3) * t551 + Ifges(5,1) * t601 + Ifges(5,4) * t600 + Ifges(5,5) * qJDD(3) - pkin(8) * t533 - qJD(3) * t587 - t643 * t531 + t532 * t678 - t610 * t586;
t507 = mrSges(4,2) * t589 - mrSges(4,3) * t579 + Ifges(4,1) * t624 + Ifges(4,4) * t625 + Ifges(4,5) * qJDD(3) - qJ(4) * t525 - qJD(3) * t616 + t518 * t641 - t519 * t639 + t615 * t667;
t506 = -mrSges(4,1) * t589 + mrSges(4,3) * t580 + Ifges(4,4) * t624 + Ifges(4,2) * t625 + Ifges(4,6) * qJDD(3) - pkin(3) * t652 + qJ(4) * t658 + qJD(3) * t617 + t639 * t518 + t641 * t519 - t615 * t668;
t505 = -pkin(2) * t517 - mrSges(3,1) * t638 + mrSges(3,3) * t599 - pkin(3) * t525 - Ifges(4,5) * t624 - Ifges(4,6) * t625 - mrSges(4,1) * t579 + mrSges(4,2) * t580 - mrSges(5,1) * t551 + mrSges(5,2) * t552 - t643 * t532 - t678 * t531 - pkin(4) * t651 - pkin(8) * t657 - Ifges(5,5) * t601 - Ifges(5,6) * t600 - t610 * t588 - t611 * t587 + t649 * Ifges(3,5) + Ifges(3,6) * qJDD(1) + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) + (-t616 * t644 + t617 * t646) * qJD(1);
t504 = mrSges(3,2) * t638 - mrSges(3,3) * t598 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t649 - pkin(7) * t517 - t506 * t644 + t507 * t646;
t503 = -mrSges(2,2) * g(3) - mrSges(2,3) * t629 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t649 - qJ(2) * t513 + t504 * t642 - t505 * t640;
t502 = mrSges(2,1) * g(3) + mrSges(2,3) * t630 + t649 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t663 + qJ(2) * t660 + t640 * t504 + t642 * t505;
t1 = [-m(1) * g(1) + t661; -m(1) * g(2) + t673; (-m(1) - m(2)) * g(3) + t663; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t673 - t645 * t502 + t647 * t503; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t661 + t647 * t502 + t645 * t503; pkin(1) * t513 - mrSges(2,2) * t630 + mrSges(2,1) * t629 + t644 * t507 + t646 * t506 + pkin(2) * t650 + pkin(7) * t659 + mrSges(3,1) * t598 - mrSges(3,2) * t599 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
