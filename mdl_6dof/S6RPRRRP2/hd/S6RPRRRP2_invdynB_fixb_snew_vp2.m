% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:13:08
% EndTime: 2019-05-06 01:13:20
% DurationCPUTime: 7.19s
% Computational Cost: add. (111498->321), mult. (216386->386), div. (0->0), fcn. (142526->10), ass. (0->126)
t678 = Ifges(6,1) + Ifges(7,1);
t675 = Ifges(6,4) + Ifges(7,4);
t674 = Ifges(6,5) + Ifges(7,5);
t677 = Ifges(6,2) + Ifges(7,2);
t673 = -Ifges(6,6) - Ifges(7,6);
t676 = -Ifges(6,3) - Ifges(7,3);
t646 = sin(qJ(1));
t650 = cos(qJ(1));
t632 = t646 * g(1) - g(2) * t650;
t624 = qJDD(1) * pkin(1) + t632;
t633 = -g(1) * t650 - g(2) * t646;
t652 = qJD(1) ^ 2;
t626 = -pkin(1) * t652 + t633;
t641 = sin(pkin(10));
t642 = cos(pkin(10));
t604 = t641 * t624 + t642 * t626;
t591 = -pkin(2) * t652 + qJDD(1) * pkin(7) + t604;
t640 = -g(3) + qJDD(2);
t645 = sin(qJ(3));
t649 = cos(qJ(3));
t582 = t649 * t591 + t645 * t640;
t625 = (-mrSges(4,1) * t649 + mrSges(4,2) * t645) * qJD(1);
t666 = qJD(1) * qJD(3);
t637 = t645 * t666;
t629 = qJDD(1) * t649 - t637;
t668 = qJD(1) * t645;
t630 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t668;
t603 = t642 * t624 - t641 * t626;
t590 = -qJDD(1) * pkin(2) - t652 * pkin(7) - t603;
t662 = t649 * t666;
t628 = qJDD(1) * t645 + t662;
t571 = (-t628 - t662) * pkin(8) + (-t629 + t637) * pkin(3) + t590;
t627 = (-pkin(3) * t649 - pkin(8) * t645) * qJD(1);
t651 = qJD(3) ^ 2;
t667 = qJD(1) * t649;
t577 = -pkin(3) * t651 + qJDD(3) * pkin(8) + t627 * t667 + t582;
t644 = sin(qJ(4));
t648 = cos(qJ(4));
t554 = t648 * t571 - t644 * t577;
t622 = qJD(3) * t648 - t644 * t668;
t599 = qJD(4) * t622 + qJDD(3) * t644 + t628 * t648;
t621 = qJDD(4) - t629;
t623 = qJD(3) * t644 + t648 * t668;
t635 = qJD(4) - t667;
t550 = (t622 * t635 - t599) * pkin(9) + (t622 * t623 + t621) * pkin(4) + t554;
t555 = t644 * t571 + t648 * t577;
t598 = -qJD(4) * t623 + qJDD(3) * t648 - t628 * t644;
t608 = pkin(4) * t635 - pkin(9) * t623;
t620 = t622 ^ 2;
t552 = -pkin(4) * t620 + pkin(9) * t598 - t608 * t635 + t555;
t643 = sin(qJ(5));
t647 = cos(qJ(5));
t544 = t647 * t550 - t643 * t552;
t601 = t622 * t647 - t623 * t643;
t561 = qJD(5) * t601 + t598 * t643 + t599 * t647;
t602 = t622 * t643 + t623 * t647;
t578 = -mrSges(7,1) * t601 + mrSges(7,2) * t602;
t579 = -mrSges(6,1) * t601 + mrSges(6,2) * t602;
t634 = qJD(5) + t635;
t584 = -mrSges(6,2) * t634 + mrSges(6,3) * t601;
t617 = qJDD(5) + t621;
t541 = -0.2e1 * qJD(6) * t602 + (t601 * t634 - t561) * qJ(6) + (t601 * t602 + t617) * pkin(5) + t544;
t583 = -mrSges(7,2) * t634 + mrSges(7,3) * t601;
t665 = m(7) * t541 + t617 * mrSges(7,1) + t634 * t583;
t533 = m(6) * t544 + t617 * mrSges(6,1) + t634 * t584 + (-t578 - t579) * t602 + (-mrSges(6,3) - mrSges(7,3)) * t561 + t665;
t545 = t643 * t550 + t647 * t552;
t560 = -qJD(5) * t602 + t598 * t647 - t599 * t643;
t586 = mrSges(7,1) * t634 - mrSges(7,3) * t602;
t587 = mrSges(6,1) * t634 - mrSges(6,3) * t602;
t585 = pkin(5) * t634 - qJ(6) * t602;
t600 = t601 ^ 2;
t543 = -pkin(5) * t600 + qJ(6) * t560 + 0.2e1 * qJD(6) * t601 - t585 * t634 + t545;
t664 = m(7) * t543 + t560 * mrSges(7,3) + t601 * t578;
t536 = m(6) * t545 + t560 * mrSges(6,3) + t601 * t579 + (-t586 - t587) * t634 + (-mrSges(6,2) - mrSges(7,2)) * t617 + t664;
t531 = t647 * t533 + t643 * t536;
t605 = -mrSges(5,1) * t622 + mrSges(5,2) * t623;
t606 = -mrSges(5,2) * t635 + mrSges(5,3) * t622;
t528 = m(5) * t554 + mrSges(5,1) * t621 - mrSges(5,3) * t599 - t605 * t623 + t606 * t635 + t531;
t607 = mrSges(5,1) * t635 - mrSges(5,3) * t623;
t657 = -t533 * t643 + t647 * t536;
t529 = m(5) * t555 - mrSges(5,2) * t621 + mrSges(5,3) * t598 + t605 * t622 - t607 * t635 + t657;
t658 = -t528 * t644 + t648 * t529;
t524 = m(4) * t582 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t629 - qJD(3) * t630 + t625 * t667 + t658;
t581 = -t645 * t591 + t640 * t649;
t631 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t667;
t576 = -qJDD(3) * pkin(3) - pkin(8) * t651 + t627 * t668 - t581;
t553 = -pkin(4) * t598 - pkin(9) * t620 + t623 * t608 + t576;
t547 = -pkin(5) * t560 - qJ(6) * t600 + t585 * t602 + qJDD(6) + t553;
t656 = m(7) * t547 - t560 * mrSges(7,1) + t561 * mrSges(7,2) - t601 * t583 + t602 * t586;
t655 = m(6) * t553 - t560 * mrSges(6,1) + t561 * mrSges(6,2) - t601 * t584 + t602 * t587 + t656;
t653 = -m(5) * t576 + t598 * mrSges(5,1) - t599 * mrSges(5,2) + t622 * t606 - t623 * t607 - t655;
t538 = m(4) * t581 + qJDD(3) * mrSges(4,1) - t628 * mrSges(4,3) + qJD(3) * t631 - t625 * t668 + t653;
t659 = t649 * t524 - t538 * t645;
t518 = m(3) * t604 - mrSges(3,1) * t652 - qJDD(1) * mrSges(3,2) + t659;
t525 = t528 * t648 + t529 * t644;
t654 = -m(4) * t590 + t629 * mrSges(4,1) - mrSges(4,2) * t628 - t630 * t668 + t631 * t667 - t525;
t521 = m(3) * t603 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t652 + t654;
t512 = t641 * t518 + t642 * t521;
t510 = m(2) * t632 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t652 + t512;
t660 = t642 * t518 - t521 * t641;
t511 = m(2) * t633 - mrSges(2,1) * t652 - qJDD(1) * mrSges(2,2) + t660;
t672 = t650 * t510 + t646 * t511;
t519 = t645 * t524 + t649 * t538;
t671 = t673 * t601 - t674 * t602 + t676 * t634;
t670 = -t677 * t601 - t675 * t602 + t673 * t634;
t669 = t675 * t601 + t678 * t602 + t674 * t634;
t663 = m(3) * t640 + t519;
t661 = -t510 * t646 + t650 * t511;
t616 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t645 + Ifges(4,4) * t649) * qJD(1);
t615 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t645 + Ifges(4,2) * t649) * qJD(1);
t614 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t645 + Ifges(4,6) * t649) * qJD(1);
t594 = Ifges(5,1) * t623 + Ifges(5,4) * t622 + Ifges(5,5) * t635;
t593 = Ifges(5,4) * t623 + Ifges(5,2) * t622 + Ifges(5,6) * t635;
t592 = Ifges(5,5) * t623 + Ifges(5,6) * t622 + Ifges(5,3) * t635;
t539 = -t561 * mrSges(7,3) - t602 * t578 + t665;
t530 = mrSges(6,2) * t553 + mrSges(7,2) * t547 - mrSges(6,3) * t544 - mrSges(7,3) * t541 - qJ(6) * t539 + t675 * t560 + t678 * t561 - t671 * t601 + t674 * t617 + t670 * t634;
t526 = -mrSges(6,1) * t553 + mrSges(6,3) * t545 - mrSges(7,1) * t547 + mrSges(7,3) * t543 - pkin(5) * t656 + qJ(6) * t664 + (-qJ(6) * t586 + t669) * t634 + (-mrSges(7,2) * qJ(6) - t673) * t617 + t671 * t602 + t675 * t561 + t677 * t560;
t515 = mrSges(5,2) * t576 - mrSges(5,3) * t554 + Ifges(5,1) * t599 + Ifges(5,4) * t598 + Ifges(5,5) * t621 - pkin(9) * t531 - t526 * t643 + t530 * t647 + t592 * t622 - t593 * t635;
t514 = -mrSges(5,1) * t576 + mrSges(5,3) * t555 + Ifges(5,4) * t599 + Ifges(5,2) * t598 + Ifges(5,6) * t621 - pkin(4) * t655 + pkin(9) * t657 + t647 * t526 + t643 * t530 - t623 * t592 + t635 * t594;
t513 = t676 * t617 + Ifges(4,6) * qJDD(3) + t673 * t560 - t674 * t561 + t669 * t601 + t670 * t602 - t623 * t593 + Ifges(4,4) * t628 + Ifges(4,2) * t629 - Ifges(5,3) * t621 + t622 * t594 + qJD(3) * t616 - Ifges(5,6) * t598 - Ifges(5,5) * t599 - mrSges(4,1) * t590 + mrSges(4,3) * t582 - mrSges(5,1) * t554 + mrSges(5,2) * t555 + mrSges(7,2) * t543 - mrSges(6,1) * t544 + mrSges(6,2) * t545 - pkin(5) * t539 - mrSges(7,1) * t541 - pkin(4) * t531 - pkin(3) * t525 - t614 * t668;
t506 = mrSges(4,2) * t590 - mrSges(4,3) * t581 + Ifges(4,1) * t628 + Ifges(4,4) * t629 + Ifges(4,5) * qJDD(3) - pkin(8) * t525 - qJD(3) * t615 - t514 * t644 + t515 * t648 + t614 * t667;
t505 = Ifges(3,6) * qJDD(1) + t652 * Ifges(3,5) - mrSges(3,1) * t640 + mrSges(3,3) * t604 - Ifges(4,5) * t628 - Ifges(4,6) * t629 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t581 + mrSges(4,2) * t582 - t644 * t515 - t648 * t514 - pkin(3) * t653 - pkin(8) * t658 - pkin(2) * t519 + (-t615 * t645 + t616 * t649) * qJD(1);
t504 = mrSges(3,2) * t640 - mrSges(3,3) * t603 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t652 - pkin(7) * t519 + t506 * t649 - t513 * t645;
t503 = -mrSges(2,2) * g(3) - mrSges(2,3) * t632 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t652 - qJ(2) * t512 + t504 * t642 - t505 * t641;
t502 = mrSges(2,1) * g(3) + mrSges(2,3) * t633 + t652 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t663 + qJ(2) * t660 + t641 * t504 + t642 * t505;
t1 = [-m(1) * g(1) + t661; -m(1) * g(2) + t672; (-m(1) - m(2)) * g(3) + t663; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t672 - t646 * t502 + t650 * t503; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t661 + t650 * t502 + t646 * t503; pkin(1) * t512 + mrSges(2,1) * t632 - mrSges(2,2) * t633 + t645 * t506 + t649 * t513 + pkin(2) * t654 + pkin(7) * t659 + mrSges(3,1) * t603 - mrSges(3,2) * t604 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
