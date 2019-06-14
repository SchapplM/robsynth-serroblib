% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRP9
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-05-05 18:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRP9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP9_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP9_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP9_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP9_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP9_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP9_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP9_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:05:10
% EndTime: 2019-05-05 18:05:17
% DurationCPUTime: 4.69s
% Computational Cost: add. (49434->313), mult. (102611->371), div. (0->0), fcn. (64256->8), ass. (0->121)
t667 = Ifges(6,1) + Ifges(7,1);
t657 = Ifges(6,4) - Ifges(7,5);
t666 = -Ifges(6,5) - Ifges(7,4);
t665 = Ifges(6,2) + Ifges(7,3);
t654 = Ifges(6,6) - Ifges(7,6);
t664 = -Ifges(6,3) - Ifges(7,2);
t625 = sin(qJ(1));
t627 = cos(qJ(1));
t612 = -t627 * g(1) - t625 * g(2);
t663 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t612;
t662 = -pkin(1) - pkin(7);
t661 = cos(qJ(5));
t660 = mrSges(2,1) - mrSges(3,2);
t659 = -mrSges(6,3) - mrSges(7,2);
t658 = -Ifges(3,4) + Ifges(2,5);
t655 = (Ifges(3,5) - Ifges(2,6));
t611 = t625 * g(1) - t627 * g(2);
t629 = qJD(1) ^ 2;
t635 = -t629 * qJ(2) + qJDD(2) - t611;
t586 = t662 * qJDD(1) + t635;
t624 = sin(qJ(3));
t626 = cos(qJ(3));
t575 = -t626 * g(3) + t624 * t586;
t606 = (mrSges(4,1) * t624 + mrSges(4,2) * t626) * qJD(1);
t646 = qJD(1) * qJD(3);
t642 = t626 * t646;
t607 = t624 * qJDD(1) + t642;
t648 = qJD(1) * t626;
t610 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t648;
t580 = t662 * t629 - t663;
t643 = t624 * t646;
t608 = t626 * qJDD(1) - t643;
t556 = (-t608 + t643) * qJ(4) + (t607 + t642) * pkin(3) + t580;
t605 = (pkin(3) * t624 - qJ(4) * t626) * qJD(1);
t628 = qJD(3) ^ 2;
t647 = t624 * qJD(1);
t561 = -t628 * pkin(3) + qJDD(3) * qJ(4) - t605 * t647 + t575;
t621 = sin(pkin(9));
t622 = cos(pkin(9));
t602 = t621 * qJD(3) + t622 * t648;
t534 = -0.2e1 * qJD(4) * t602 + t622 * t556 - t621 * t561;
t584 = t621 * qJDD(3) + t622 * t608;
t601 = t622 * qJD(3) - t621 * t648;
t531 = (t601 * t647 - t584) * pkin(8) + (t601 * t602 + t607) * pkin(4) + t534;
t535 = 0.2e1 * qJD(4) * t601 + t621 * t556 + t622 * t561;
t583 = t622 * qJDD(3) - t621 * t608;
t585 = pkin(4) * t647 - t602 * pkin(8);
t600 = t601 ^ 2;
t533 = -t600 * pkin(4) + t583 * pkin(8) - t585 * t647 + t535;
t623 = sin(qJ(5));
t527 = t623 * t531 + t661 * t533;
t572 = t623 * t601 + t661 * t602;
t540 = t572 * qJD(5) - t661 * t583 + t623 * t584;
t614 = qJD(5) + t647;
t564 = t614 * mrSges(6,1) - t572 * mrSges(6,3);
t571 = -t661 * t601 + t623 * t602;
t604 = qJDD(5) + t607;
t551 = t571 * pkin(5) - t572 * qJ(6);
t613 = t614 ^ 2;
t524 = -t613 * pkin(5) + t604 * qJ(6) + 0.2e1 * qJD(6) * t614 - t571 * t551 + t527;
t565 = -t614 * mrSges(7,1) + t572 * mrSges(7,2);
t644 = m(7) * t524 + t604 * mrSges(7,3) + t614 * t565;
t552 = t571 * mrSges(7,1) - t572 * mrSges(7,3);
t649 = -t571 * mrSges(6,1) - t572 * mrSges(6,2) - t552;
t518 = m(6) * t527 - t604 * mrSges(6,2) + t659 * t540 - t614 * t564 + t649 * t571 + t644;
t526 = t661 * t531 - t623 * t533;
t541 = -t571 * qJD(5) + t623 * t583 + t661 * t584;
t562 = -t614 * mrSges(6,2) - t571 * mrSges(6,3);
t525 = -t604 * pkin(5) - t613 * qJ(6) + t572 * t551 + qJDD(6) - t526;
t563 = -t571 * mrSges(7,2) + t614 * mrSges(7,3);
t637 = -m(7) * t525 + t604 * mrSges(7,1) + t614 * t563;
t520 = m(6) * t526 + t604 * mrSges(6,1) + t659 * t541 + t614 * t562 + t649 * t572 + t637;
t515 = t623 * t518 + t661 * t520;
t573 = -t601 * mrSges(5,1) + t602 * mrSges(5,2);
t581 = -mrSges(5,2) * t647 + t601 * mrSges(5,3);
t511 = m(5) * t534 + t607 * mrSges(5,1) - t584 * mrSges(5,3) - t602 * t573 + t581 * t647 + t515;
t582 = mrSges(5,1) * t647 - t602 * mrSges(5,3);
t638 = t661 * t518 - t623 * t520;
t512 = m(5) * t535 - t607 * mrSges(5,2) + t583 * mrSges(5,3) + t601 * t573 - t582 * t647 + t638;
t639 = -t621 * t511 + t622 * t512;
t506 = m(4) * t575 - qJDD(3) * mrSges(4,2) - t607 * mrSges(4,3) - qJD(3) * t610 - t606 * t647 + t639;
t574 = t624 * g(3) + t626 * t586;
t609 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t647;
t558 = -qJDD(3) * pkin(3) - t628 * qJ(4) + t605 * t648 + qJDD(4) - t574;
t536 = -t583 * pkin(4) - t600 * pkin(8) + t602 * t585 + t558;
t529 = -0.2e1 * qJD(6) * t572 + (t571 * t614 - t541) * qJ(6) + (t572 * t614 + t540) * pkin(5) + t536;
t522 = m(7) * t529 + t540 * mrSges(7,1) - t541 * mrSges(7,3) + t571 * t563 - t572 * t565;
t631 = m(6) * t536 + t540 * mrSges(6,1) + t541 * mrSges(6,2) + t571 * t562 + t572 * t564 + t522;
t630 = -m(5) * t558 + t583 * mrSges(5,1) - t584 * mrSges(5,2) + t601 * t581 - t602 * t582 - t631;
t521 = m(4) * t574 + qJDD(3) * mrSges(4,1) - t608 * mrSges(4,3) + qJD(3) * t609 - t606 * t648 + t630;
t501 = t624 * t506 + t626 * t521;
t588 = -qJDD(1) * pkin(1) + t635;
t634 = -m(3) * t588 + (t629 * mrSges(3,3)) - t501;
t499 = m(2) * t611 - (t629 * mrSges(2,2)) + t660 * qJDD(1) + t634;
t587 = t629 * pkin(1) + t663;
t507 = t622 * t511 + t621 * t512;
t633 = -m(4) * t580 - t607 * mrSges(4,1) - t608 * mrSges(4,2) - t609 * t647 - t610 * t648 - t507;
t632 = -m(3) * t587 + (t629 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t633;
t504 = m(2) * t612 - (t629 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t632;
t653 = t627 * t499 + t625 * t504;
t652 = t665 * t571 - t657 * t572 - t654 * t614;
t651 = t654 * t571 + t666 * t572 + t664 * t614;
t650 = -t657 * t571 + t667 * t572 - t666 * t614;
t641 = -t625 * t499 + t627 * t504;
t640 = t626 * t506 - t624 * t521;
t595 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t626 - Ifges(4,4) * t624) * qJD(1);
t594 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t626 - Ifges(4,2) * t624) * qJD(1);
t593 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t626 - Ifges(4,6) * t624) * qJD(1);
t568 = Ifges(5,1) * t602 + Ifges(5,4) * t601 + Ifges(5,5) * t647;
t567 = Ifges(5,4) * t602 + Ifges(5,2) * t601 + Ifges(5,6) * t647;
t566 = Ifges(5,5) * t602 + Ifges(5,6) * t601 + Ifges(5,3) * t647;
t514 = mrSges(6,2) * t536 + mrSges(7,2) * t525 - mrSges(6,3) * t526 - mrSges(7,3) * t529 - qJ(6) * t522 - t657 * t540 + t667 * t541 + t651 * t571 - t604 * t666 + t652 * t614;
t513 = -mrSges(6,1) * t536 - mrSges(7,1) * t529 + mrSges(7,2) * t524 + mrSges(6,3) * t527 - pkin(5) * t522 - t665 * t540 + t657 * t541 + t651 * t572 + t654 * t604 + t650 * t614;
t500 = -m(3) * g(3) + t640;
t497 = mrSges(5,2) * t558 - mrSges(5,3) * t534 + Ifges(5,1) * t584 + Ifges(5,4) * t583 + Ifges(5,5) * t607 - pkin(8) * t515 - t623 * t513 + t661 * t514 + t601 * t566 - t567 * t647;
t496 = -mrSges(5,1) * t558 + mrSges(5,3) * t535 + Ifges(5,4) * t584 + Ifges(5,2) * t583 + Ifges(5,6) * t607 - pkin(4) * t631 + pkin(8) * t638 + t661 * t513 + t623 * t514 - t602 * t566 + t568 * t647;
t495 = -mrSges(6,1) * t526 + mrSges(6,2) * t527 + qJD(3) * t595 + t601 * t568 - t602 * t567 - mrSges(7,3) * t524 + mrSges(7,1) * t525 + (pkin(5) * mrSges(7,2) + t666) * t541 + Ifges(4,4) * t608 - pkin(5) * t637 - pkin(4) * t515 - mrSges(5,1) * t534 + mrSges(5,2) * t535 - pkin(3) * t507 - mrSges(4,1) * t580 - Ifges(5,6) * t583 - Ifges(5,5) * t584 + mrSges(4,3) * t575 + (qJ(6) * mrSges(7,2) + t654) * t540 - t593 * t648 + (qJ(6) * t552 - t650) * t571 + (pkin(5) * t552 + t652) * t572 - qJ(6) * t644 + Ifges(4,6) * qJDD(3) + t664 * t604 + (-Ifges(5,3) - Ifges(4,2)) * t607;
t494 = mrSges(4,2) * t580 - mrSges(4,3) * t574 + Ifges(4,1) * t608 - Ifges(4,4) * t607 + Ifges(4,5) * qJDD(3) - qJ(4) * t507 - qJD(3) * t594 - t621 * t496 + t622 * t497 - t593 * t647;
t493 = -qJ(2) * t500 - mrSges(2,3) * t611 + pkin(2) * t501 + mrSges(3,1) * t588 + pkin(3) * t630 + qJ(4) * t639 + t621 * t497 + t622 * t496 + Ifges(4,5) * t608 - Ifges(4,6) * t607 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t574 - mrSges(4,2) * t575 + (t655 * t629) + t658 * qJDD(1) + (t626 * t594 + t624 * t595) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t492 = -mrSges(3,1) * t587 + mrSges(2,3) * t612 - pkin(1) * t500 - pkin(2) * t633 - pkin(7) * t640 + t660 * g(3) - t655 * qJDD(1) - t624 * t494 - t626 * t495 + t658 * t629;
t1 = [-m(1) * g(1) + t641; -m(1) * g(2) + t653; (-m(1) - m(2) - m(3)) * g(3) + t640; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t653 - t625 * t492 + t627 * t493; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t641 + t627 * t492 + t625 * t493; pkin(1) * t634 + qJ(2) * t632 + t626 * t494 - t624 * t495 - pkin(7) * t501 + mrSges(2,1) * t611 - mrSges(2,2) * t612 + mrSges(3,2) * t588 - mrSges(3,3) * t587 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
