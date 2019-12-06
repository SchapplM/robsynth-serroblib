% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP2_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP2_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP2_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:48
% EndTime: 2019-12-05 15:30:50
% DurationCPUTime: 2.19s
% Computational Cost: add. (15329->218), mult. (32871->277), div. (0->0), fcn. (19734->8), ass. (0->106)
t624 = sin(pkin(7));
t626 = cos(pkin(7));
t606 = g(1) * t624 - g(2) * t626;
t607 = -g(1) * t626 - g(2) * t624;
t628 = sin(qJ(2));
t630 = cos(qJ(2));
t585 = t628 * t606 + t630 * t607;
t631 = qJD(2) ^ 2;
t681 = -pkin(2) * t631 + qJDD(2) * qJ(3) + (2 * qJD(2) * qJD(3)) + t585;
t680 = Ifges(5,1) + Ifges(6,1);
t673 = Ifges(5,4) + Ifges(6,4);
t672 = Ifges(5,5) + Ifges(6,5);
t679 = Ifges(5,2) + Ifges(6,2);
t671 = Ifges(5,6) + Ifges(6,6);
t678 = Ifges(5,3) + Ifges(6,3);
t622 = -g(3) + qJDD(1);
t623 = sin(pkin(8));
t625 = cos(pkin(8));
t568 = t622 * t625 - t681 * t623;
t627 = sin(qJ(4));
t629 = cos(qJ(4));
t659 = qJD(2) * t625;
t610 = qJD(4) - t659;
t660 = qJD(2) * t623;
t662 = (-t673 * t627 + t680 * t629) * t660 + t672 * t610;
t663 = (-t679 * t627 + t673 * t629) * t660 + t671 * t610;
t677 = t662 * t627 + t663 * t629;
t595 = (mrSges(6,1) * t627 + mrSges(6,2) * t629) * t660;
t656 = qJD(2) * qJD(4);
t598 = (qJDD(2) * t629 - t627 * t656) * t623;
t650 = t629 * t660;
t569 = t623 * t622 + t681 * t625;
t640 = -pkin(3) * t625 - pkin(6) * t623;
t605 = t640 * qJD(2);
t567 = t605 * t659 + t569;
t584 = t606 * t630 - t628 * t607;
t636 = -qJ(3) * t631 + qJDD(3) - t584;
t572 = (-pkin(2) + t640) * qJDD(2) + t636;
t571 = t629 * t572;
t655 = qJDD(2) * t625;
t609 = qJDD(4) - t655;
t642 = -0.2e1 * qJD(5) * t660;
t667 = t623 ^ 2 * t631;
t559 = t629 * t642 + pkin(4) * t609 - t598 * qJ(5) + t571 + (-pkin(4) * t629 * t667 - qJ(5) * t610 * t660 - t567) * t627;
t651 = t627 * t660;
t590 = -mrSges(6,2) * t610 - mrSges(6,3) * t651;
t653 = m(6) * t559 + t609 * mrSges(6,1) + t610 * t590;
t556 = -t598 * mrSges(6,3) - t595 * t650 + t653;
t564 = t629 * t567 + t627 * t572;
t592 = pkin(4) * t610 - qJ(5) * t650;
t597 = (-qJDD(2) * t627 - t629 * t656) * t623;
t654 = t627 ^ 2 * t667;
t561 = -pkin(4) * t654 + t597 * qJ(5) - t592 * t610 + t627 * t642 + t564;
t563 = -t567 * t627 + t571;
t676 = mrSges(5,1) * t563 + mrSges(6,1) * t559 - mrSges(5,2) * t564 - mrSges(6,2) * t561 + pkin(4) * t556 + t671 * t597 + t672 * t598 + t678 * t609;
t675 = t625 ^ 2;
t674 = -mrSges(5,2) - mrSges(6,2);
t669 = mrSges(4,2) * t623;
t668 = t598 * mrSges(6,2);
t603 = (-mrSges(4,1) * t625 + t669) * qJD(2);
t591 = -mrSges(5,2) * t610 - mrSges(5,3) * t651;
t639 = (-t595 - (mrSges(5,1) * t627 + mrSges(5,2) * t629) * t660) * t660;
t552 = m(5) * t563 + mrSges(5,1) * t609 + t591 * t610 + (-mrSges(5,3) - mrSges(6,3)) * t598 + t629 * t639 + t653;
t593 = mrSges(6,1) * t610 - mrSges(6,3) * t650;
t661 = -mrSges(5,1) * t610 + mrSges(5,3) * t650 - t593;
t665 = m(6) * t561 + t597 * mrSges(6,3);
t553 = m(5) * t564 + t597 * mrSges(5,3) + t674 * t609 + t661 * t610 + t627 * t639 + t665;
t644 = -t552 * t627 + t629 * t553;
t658 = qJDD(2) * mrSges(4,3);
t546 = m(4) * t569 + (qJD(2) * t603 + t658) * t625 + t644;
t635 = t661 * t629 + (-t590 - t591) * t627;
t566 = t605 * t660 - t568;
t562 = -t597 * pkin(4) - qJ(5) * t654 + t592 * t650 + qJDD(5) + t566;
t648 = m(6) * t562 - t597 * mrSges(6,1);
t637 = -m(5) * t566 + t597 * mrSges(5,1) - t648;
t555 = m(4) * t568 + t674 * t598 + (-t658 + (-t603 + t635) * qJD(2)) * t623 + t637;
t645 = t625 * t546 - t555 * t623;
t538 = m(3) * t585 - mrSges(3,1) * t631 - qJDD(2) * mrSges(3,2) + t645;
t550 = t552 * t629 + t553 * t627;
t575 = -qJDD(2) * pkin(2) + t636;
t633 = -m(4) * t575 + mrSges(4,1) * t655 - t550 + (t631 * t675 + t667) * mrSges(4,3);
t543 = m(3) * t584 - mrSges(3,2) * t631 + (mrSges(3,1) - t669) * qJDD(2) + t633;
t533 = t628 * t538 + t630 * t543;
t531 = m(2) * t606 + t533;
t646 = t630 * t538 - t628 * t543;
t532 = m(2) * t607 + t646;
t666 = t626 * t531 + t624 * t532;
t540 = t623 * t546 + t625 * t555;
t664 = (t671 * t627 - t672 * t629) * t660 - t678 * t610;
t652 = m(3) * t622 + t540;
t647 = -t531 * t624 + t626 * t532;
t641 = m(2) * t622 + t652;
t638 = Ifges(4,5) * t623 + Ifges(4,6) * t625;
t557 = t668 + (t590 * t627 + t593 * t629) * t660 + t648;
t541 = -mrSges(5,1) * t566 + mrSges(5,3) * t564 - mrSges(6,1) * t562 + mrSges(6,3) * t561 - pkin(4) * t557 + qJ(5) * t665 + (-qJ(5) * t593 + t662) * t610 + (-qJ(5) * mrSges(6,2) + t671) * t609 + t673 * t598 + t679 * t597 + (-qJ(5) * t595 * t627 + t664 * t629) * t660;
t548 = mrSges(5,2) * t566 + mrSges(6,2) * t562 - mrSges(5,3) * t563 - mrSges(6,3) * t559 - qJ(5) * t556 + t673 * t597 + t680 * t598 + t672 * t609 - t663 * t610 + t664 * t651;
t604 = t638 * qJD(2);
t527 = t604 * t659 + mrSges(4,2) * t575 - mrSges(4,3) * t568 - pkin(6) * t550 - t541 * t627 + t548 * t629 + (Ifges(4,1) * t623 + Ifges(4,4) * t625) * qJDD(2);
t535 = Ifges(4,2) * t655 - mrSges(4,1) * t575 + mrSges(4,3) * t569 - pkin(3) * t550 + (Ifges(4,4) * qJDD(2) + (-t604 - t677) * qJD(2)) * t623 - t676;
t549 = qJDD(2) * t669 - t633;
t634 = mrSges(3,1) * t584 - mrSges(3,2) * t585 + Ifges(3,3) * qJDD(2) - pkin(2) * t549 + qJ(3) * t645 + t623 * t527 + t625 * t535;
t525 = t631 * Ifges(3,5) - mrSges(3,1) * t622 + mrSges(3,3) * t585 - mrSges(4,1) * t568 + mrSges(4,2) * t569 - t627 * t548 - t629 * t541 - pkin(3) * (-t598 * mrSges(5,2) + t637 - t668) - pkin(6) * t644 - pkin(2) * t540 + (Ifges(3,6) - t638) * qJDD(2) + (Ifges(4,4) * t675 * qJD(2) + (-pkin(3) * t635 + (-Ifges(4,4) * t623 + (Ifges(4,1) - Ifges(4,2)) * t625) * qJD(2)) * t623) * qJD(2);
t524 = mrSges(3,2) * t622 - mrSges(3,3) * t584 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t631 - qJ(3) * t540 + t527 * t625 - t535 * t623;
t523 = mrSges(2,2) * t622 - mrSges(2,3) * t606 - pkin(5) * t533 + t524 * t630 - t525 * t628;
t522 = -mrSges(2,1) * t622 + mrSges(2,3) * t607 - pkin(1) * t652 + pkin(5) * t646 + t628 * t524 + t630 * t525;
t1 = [-m(1) * g(1) + t647; -m(1) * g(2) + t666; -m(1) * g(3) + t641; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t666 - t624 * t522 + t626 * t523; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t647 + t626 * t522 + t624 * t523; -mrSges(1,1) * g(2) + mrSges(2,1) * t606 + mrSges(1,2) * g(1) - mrSges(2,2) * t607 + pkin(1) * t533 + t634; t641; t634; t549; t677 * t660 + t676; t557;];
tauJB = t1;
