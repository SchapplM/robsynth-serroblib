% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-05-05 14:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:47:17
% EndTime: 2019-05-05 14:47:24
% DurationCPUTime: 6.91s
% Computational Cost: add. (77026->298), mult. (169780->358), div. (0->0), fcn. (116025->10), ass. (0->127)
t682 = Ifges(6,1) + Ifges(7,1);
t675 = Ifges(6,4) - Ifges(7,5);
t681 = -Ifges(6,5) - Ifges(7,4);
t680 = Ifges(6,2) + Ifges(7,3);
t673 = Ifges(6,6) - Ifges(7,6);
t679 = -Ifges(6,3) - Ifges(7,2);
t641 = qJD(1) ^ 2;
t631 = sin(pkin(10));
t633 = cos(pkin(10));
t636 = sin(qJ(4));
t638 = cos(qJ(4));
t646 = t631 * t636 - t633 * t638;
t606 = t646 * qJD(1);
t647 = t631 * t638 + t633 * t636;
t607 = t647 * qJD(1);
t662 = qJD(4) * t607;
t596 = -qJDD(1) * t646 - t662;
t678 = cos(qJ(5));
t677 = pkin(3) * t633;
t676 = -mrSges(6,3) - mrSges(7,2);
t672 = mrSges(4,2) * t631;
t629 = t633 ^ 2;
t671 = t629 * t641;
t637 = sin(qJ(1));
t639 = cos(qJ(1));
t615 = t637 * g(1) - g(2) * t639;
t613 = qJDD(1) * pkin(1) + t615;
t616 = -g(1) * t639 - g(2) * t637;
t614 = -pkin(1) * t641 + t616;
t632 = sin(pkin(9));
t634 = cos(pkin(9));
t599 = t632 * t613 + t634 * t614;
t588 = -pkin(2) * t641 + qJDD(1) * qJ(3) + t599;
t630 = -g(3) + qJDD(2);
t661 = qJD(1) * qJD(3);
t665 = t633 * t630 - 0.2e1 * t631 * t661;
t566 = (-pkin(7) * qJDD(1) + t641 * t677 - t588) * t631 + t665;
t574 = t631 * t630 + (t588 + 0.2e1 * t661) * t633;
t660 = qJDD(1) * t633;
t569 = -pkin(3) * t671 + pkin(7) * t660 + t574;
t552 = t636 * t566 + t638 * t569;
t592 = mrSges(5,1) * t606 + mrSges(5,2) * t607;
t603 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t607;
t595 = pkin(4) * t606 - pkin(8) * t607;
t640 = qJD(4) ^ 2;
t548 = -pkin(4) * t640 + qJDD(4) * pkin(8) - t595 * t606 + t552;
t628 = t631 ^ 2;
t598 = t613 * t634 - t632 * t614;
t648 = qJDD(3) - t598;
t572 = (-pkin(2) - t677) * qJDD(1) + (-qJ(3) + (-t628 - t629) * pkin(7)) * t641 + t648;
t663 = qJD(4) * t606;
t597 = qJDD(1) * t647 - t663;
t550 = (-t597 + t663) * pkin(8) + (-t596 + t662) * pkin(4) + t572;
t635 = sin(qJ(5));
t545 = t678 * t548 + t635 * t550;
t601 = t635 * qJD(4) + t607 * t678;
t567 = qJD(5) * t601 - qJDD(4) * t678 + t597 * t635;
t605 = qJD(5) + t606;
t581 = mrSges(6,1) * t605 - mrSges(6,3) * t601;
t594 = qJDD(5) - t596;
t600 = -qJD(4) * t678 + t607 * t635;
t575 = pkin(5) * t600 - qJ(6) * t601;
t604 = t605 ^ 2;
t541 = -pkin(5) * t604 + qJ(6) * t594 + 0.2e1 * qJD(6) * t605 - t575 * t600 + t545;
t582 = -mrSges(7,1) * t605 + mrSges(7,2) * t601;
t659 = m(7) * t541 + t594 * mrSges(7,3) + t605 * t582;
t576 = mrSges(7,1) * t600 - mrSges(7,3) * t601;
t666 = -mrSges(6,1) * t600 - mrSges(6,2) * t601 - t576;
t536 = m(6) * t545 - mrSges(6,2) * t594 + t567 * t676 - t581 * t605 + t600 * t666 + t659;
t544 = -t635 * t548 + t550 * t678;
t568 = -t600 * qJD(5) + t635 * qJDD(4) + t597 * t678;
t580 = -mrSges(6,2) * t605 - mrSges(6,3) * t600;
t542 = -t594 * pkin(5) - t604 * qJ(6) + t601 * t575 + qJDD(6) - t544;
t579 = -mrSges(7,2) * t600 + mrSges(7,3) * t605;
t652 = -m(7) * t542 + t594 * mrSges(7,1) + t605 * t579;
t538 = m(6) * t544 + mrSges(6,1) * t594 + t568 * t676 + t580 * t605 + t601 * t666 + t652;
t653 = t678 * t536 - t538 * t635;
t528 = m(5) * t552 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t596 - qJD(4) * t603 - t592 * t606 + t653;
t551 = t566 * t638 - t636 * t569;
t602 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t606;
t547 = -qJDD(4) * pkin(4) - pkin(8) * t640 + t607 * t595 - t551;
t543 = -0.2e1 * qJD(6) * t601 + (t600 * t605 - t568) * qJ(6) + (t601 * t605 + t567) * pkin(5) + t547;
t539 = m(7) * t543 + mrSges(7,1) * t567 - t568 * mrSges(7,3) + t579 * t600 - t601 * t582;
t642 = -m(6) * t547 - t567 * mrSges(6,1) - mrSges(6,2) * t568 - t600 * t580 - t581 * t601 - t539;
t533 = m(5) * t551 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t597 + qJD(4) * t602 - t592 * t607 + t642;
t523 = t636 * t528 + t638 * t533;
t573 = -t588 * t631 + t665;
t645 = mrSges(4,3) * qJDD(1) + t641 * (-mrSges(4,1) * t633 + t672);
t521 = m(4) * t573 - t631 * t645 + t523;
t654 = t638 * t528 - t533 * t636;
t522 = m(4) * t574 + t633 * t645 + t654;
t655 = -t521 * t631 + t633 * t522;
t514 = m(3) * t599 - mrSges(3,1) * t641 - qJDD(1) * mrSges(3,2) + t655;
t584 = -qJDD(1) * pkin(2) - qJ(3) * t641 + t648;
t531 = t635 * t536 + t678 * t538;
t644 = m(5) * t572 - t596 * mrSges(5,1) + mrSges(5,2) * t597 + t606 * t602 + t603 * t607 + t531;
t643 = -m(4) * t584 + mrSges(4,1) * t660 - t644 + (t628 * t641 + t671) * mrSges(4,3);
t525 = t643 + (mrSges(3,1) - t672) * qJDD(1) + m(3) * t598 - mrSges(3,2) * t641;
t511 = t632 * t514 + t634 * t525;
t509 = m(2) * t615 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t641 + t511;
t656 = t634 * t514 - t525 * t632;
t510 = m(2) * t616 - mrSges(2,1) * t641 - qJDD(1) * mrSges(2,2) + t656;
t670 = t639 * t509 + t637 * t510;
t515 = t633 * t521 + t631 * t522;
t669 = t600 * t680 - t601 * t675 - t605 * t673;
t668 = t600 * t673 + t601 * t681 + t605 * t679;
t667 = -t675 * t600 + t601 * t682 - t681 * t605;
t649 = Ifges(4,5) * t631 + Ifges(4,6) * t633;
t664 = t641 * t649;
t658 = m(3) * t630 + t515;
t657 = -t509 * t637 + t639 * t510;
t651 = Ifges(4,1) * t631 + Ifges(4,4) * t633;
t650 = Ifges(4,4) * t631 + Ifges(4,2) * t633;
t587 = Ifges(5,1) * t607 - Ifges(5,4) * t606 + Ifges(5,5) * qJD(4);
t586 = Ifges(5,4) * t607 - Ifges(5,2) * t606 + Ifges(5,6) * qJD(4);
t585 = Ifges(5,5) * t607 - Ifges(5,6) * t606 + Ifges(5,3) * qJD(4);
t530 = mrSges(6,2) * t547 + mrSges(7,2) * t542 - mrSges(6,3) * t544 - mrSges(7,3) * t543 - qJ(6) * t539 - t675 * t567 + t568 * t682 - t594 * t681 + t668 * t600 + t669 * t605;
t529 = -mrSges(6,1) * t547 - mrSges(7,1) * t543 + mrSges(7,2) * t541 + mrSges(6,3) * t545 - pkin(5) * t539 - t567 * t680 + t675 * t568 + t673 * t594 + t668 * t601 + t667 * t605;
t517 = Ifges(5,4) * t597 + Ifges(5,2) * t596 + Ifges(5,6) * qJDD(4) - t607 * t585 + qJD(4) * t587 - mrSges(5,1) * t572 + mrSges(5,3) * t552 - mrSges(6,1) * t544 + mrSges(6,2) * t545 + mrSges(7,1) * t542 - mrSges(7,3) * t541 - pkin(5) * t652 - qJ(6) * t659 - pkin(4) * t531 + (pkin(5) * t576 + t669) * t601 + (qJ(6) * t576 - t667) * t600 + t679 * t594 + (mrSges(7,2) * pkin(5) + t681) * t568 + (mrSges(7,2) * qJ(6) + t673) * t567;
t516 = mrSges(5,2) * t572 - mrSges(5,3) * t551 + Ifges(5,1) * t597 + Ifges(5,4) * t596 + Ifges(5,5) * qJDD(4) - pkin(8) * t531 - qJD(4) * t586 - t635 * t529 + t530 * t678 - t606 * t585;
t505 = mrSges(4,2) * t584 - mrSges(4,3) * t573 - pkin(7) * t523 + qJDD(1) * t651 + t516 * t638 - t517 * t636 + t633 * t664;
t504 = -mrSges(4,1) * t584 + mrSges(4,3) * t574 - pkin(3) * t644 + pkin(7) * t654 + qJDD(1) * t650 + t636 * t516 + t638 * t517 - t631 * t664;
t503 = -pkin(2) * t515 - mrSges(3,1) * t630 + mrSges(3,3) * t599 - pkin(3) * t523 + mrSges(4,2) * t574 - mrSges(4,1) * t573 - pkin(8) * t653 - t635 * t530 - t678 * t529 - pkin(4) * t642 - t606 * t587 - mrSges(5,1) * t551 + mrSges(5,2) * t552 - Ifges(5,5) * t597 - Ifges(5,6) * t596 - Ifges(5,3) * qJDD(4) - t607 * t586 + (Ifges(3,6) - t649) * qJDD(1) + (-t631 * t650 + t633 * t651 + Ifges(3,5)) * t641;
t502 = mrSges(3,2) * t630 - mrSges(3,3) * t598 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t641 - qJ(3) * t515 - t504 * t631 + t505 * t633;
t501 = -mrSges(2,2) * g(3) - mrSges(2,3) * t615 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t641 - qJ(2) * t511 + t502 * t634 - t503 * t632;
t500 = mrSges(2,1) * g(3) + mrSges(2,3) * t616 + t641 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t658 + qJ(2) * t656 + t632 * t502 + t634 * t503;
t1 = [-m(1) * g(1) + t657; -m(1) * g(2) + t670; (-m(1) - m(2)) * g(3) + t658; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t670 - t637 * t500 + t639 * t501; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t657 + t639 * t500 + t637 * t501; pkin(1) * t511 + mrSges(2,1) * t615 - mrSges(2,2) * t616 + t631 * t505 + t633 * t504 + pkin(2) * t643 + qJ(3) * t655 + mrSges(3,1) * t598 - mrSges(3,2) * t599 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * t672 + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
