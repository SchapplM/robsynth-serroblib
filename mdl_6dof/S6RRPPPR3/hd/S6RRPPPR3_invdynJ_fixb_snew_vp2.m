% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 08:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:33:48
% EndTime: 2019-05-06 08:33:52
% DurationCPUTime: 2.70s
% Computational Cost: add. (14010->315), mult. (30607->377), div. (0->0), fcn. (16345->8), ass. (0->129)
t638 = Ifges(3,1) + Ifges(4,1) + Ifges(5,2);
t615 = Ifges(3,4) - Ifges(4,5) + Ifges(5,4);
t614 = Ifges(3,5) + Ifges(4,4) + Ifges(5,6);
t613 = Ifges(3,6) - Ifges(4,6) + Ifges(5,5);
t637 = Ifges(3,3) + Ifges(4,2) + Ifges(5,3);
t636 = Ifges(4,3) + Ifges(5,1) + Ifges(3,2);
t589 = cos(qJ(2));
t586 = sin(qJ(2));
t618 = qJD(1) * qJD(2);
t607 = t586 * t618;
t552 = qJDD(1) * t589 - t607;
t620 = qJD(1) * t586;
t557 = -qJD(2) * pkin(3) - qJ(4) * t620;
t617 = qJD(1) * qJD(4);
t579 = t589 ^ 2;
t592 = qJD(1) ^ 2;
t625 = t579 * t592;
t587 = sin(qJ(1));
t590 = cos(qJ(1));
t604 = -g(1) * t590 - g(2) * t587;
t533 = -pkin(1) * t592 + qJDD(1) * pkin(7) + t604;
t508 = -g(3) * t586 + t589 * t533;
t546 = (-pkin(2) * t589 - qJ(3) * t586) * qJD(1);
t619 = qJD(1) * t589;
t634 = qJDD(2) * qJ(3) + t546 * t619 + t508;
t635 = pkin(3) * t625 + t552 * qJ(4) - qJD(2) * t557 + 0.2e1 * t589 * t617 - t634;
t609 = t589 * t618;
t551 = qJDD(1) * t586 + t609;
t633 = -0.2e1 * t586 * t617 + (-t551 + t609) * qJ(4);
t632 = 2 * qJD(5);
t591 = qJD(2) ^ 2;
t631 = pkin(2) * t591;
t630 = -mrSges(4,1) + mrSges(5,2);
t629 = mrSges(4,2) - mrSges(5,3);
t628 = -pkin(2) - qJ(5);
t627 = pkin(3) + qJ(5);
t624 = t589 * t592;
t621 = t587 * g(1) - t590 * g(2);
t532 = -qJDD(1) * pkin(1) - t592 * pkin(7) - t621;
t602 = -t552 * pkin(2) + t532 + (-t551 - t609) * qJ(3);
t596 = -qJ(4) * t625 + qJDD(4) - t602 + ((2 * qJD(3)) + t557) * t620;
t475 = t596 + t627 * t552 + pkin(4) * t551 + (pkin(4) * t589 + t586 * t628) * t618;
t549 = (pkin(4) * t586 + qJ(5) * t589) * qJD(1);
t507 = -t589 * g(3) - t586 * t533;
t605 = t546 * t620 + qJDD(3) - t507;
t480 = (-pkin(4) - qJ(3)) * t591 + (-pkin(3) * t624 - qJD(1) * t549) * t586 + (-pkin(2) - t627) * qJDD(2) + t605 + t633;
t582 = sin(pkin(9));
t583 = cos(pkin(9));
t541 = qJD(2) * t582 + t583 * t619;
t469 = t583 * t475 - t480 * t582 + t541 * t632;
t513 = -qJDD(2) * t582 - t552 * t583;
t540 = -qJD(2) * t583 + t582 * t619;
t467 = (t540 * t620 - t513) * pkin(8) + (-t540 * t541 + t551) * pkin(5) + t469;
t470 = t582 * t475 + t583 * t480 + t540 * t632;
t512 = -qJDD(2) * t583 + t552 * t582;
t514 = pkin(5) * t620 + pkin(8) * t541;
t537 = t540 ^ 2;
t468 = -pkin(5) * t537 + pkin(8) * t512 - t514 * t620 + t470;
t585 = sin(qJ(6));
t588 = cos(qJ(6));
t465 = t467 * t588 - t468 * t585;
t504 = t540 * t588 + t541 * t585;
t488 = qJD(6) * t504 + t512 * t585 + t513 * t588;
t505 = t540 * t585 - t541 * t588;
t493 = -mrSges(7,1) * t504 + mrSges(7,2) * t505;
t566 = qJD(6) + t620;
t498 = -mrSges(7,2) * t566 + mrSges(7,3) * t504;
t544 = qJDD(6) + t551;
t461 = m(7) * t465 + mrSges(7,1) * t544 - mrSges(7,3) * t488 - t493 * t505 + t498 * t566;
t466 = t467 * t585 + t468 * t588;
t487 = -qJD(6) * t505 + t512 * t588 - t513 * t585;
t499 = mrSges(7,1) * t566 - mrSges(7,3) * t505;
t462 = m(7) * t466 - mrSges(7,2) * t544 + mrSges(7,3) * t487 + t493 * t504 - t499 * t566;
t455 = t588 * t461 + t585 * t462;
t506 = -mrSges(6,1) * t540 - mrSges(6,2) * t541;
t510 = -mrSges(6,2) * t620 + mrSges(6,3) * t540;
t453 = m(6) * t469 + mrSges(6,1) * t551 - mrSges(6,3) * t513 + t506 * t541 + t510 * t620 + t455;
t511 = mrSges(6,1) * t620 + mrSges(6,3) * t541;
t606 = -t461 * t585 + t588 * t462;
t454 = m(6) * t470 - mrSges(6,2) * t551 + mrSges(6,3) * t512 + t506 * t540 - t511 * t620 + t606;
t448 = t583 * t453 + t582 * t454;
t623 = -t582 * t453 + t583 * t454;
t558 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t620;
t560 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t620;
t622 = -t558 - t560;
t616 = qJD(3) * qJD(2);
t612 = -t637 * qJD(2) + (-t586 * t614 - t589 * t613) * qJD(1);
t611 = -t613 * qJD(2) + (-t586 * t615 - t636 * t589) * qJD(1);
t610 = t614 * qJD(2) + (t586 * t638 + t589 * t615) * qJD(1);
t481 = -pkin(2) * t607 + pkin(3) * t552 + t596;
t561 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t619;
t603 = m(5) * t481 + t551 * mrSges(5,1) - t561 * t619 + t448;
t571 = 0.2e1 * t616;
t479 = qJDD(2) * pkin(4) - t549 * t619 + t591 * t628 + qJDD(5) + t571 - t635;
t472 = -pkin(5) * t512 - pkin(8) * t537 - t514 * t541 + t479;
t601 = m(7) * t472 - t487 * mrSges(7,1) + t488 * mrSges(7,2) - t504 * t498 + t505 * t499;
t494 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t620 + t602;
t600 = m(4) * t494 - t603;
t497 = -qJDD(2) * pkin(2) - qJ(3) * t591 + t605;
t490 = Ifges(7,4) * t505 + Ifges(7,2) * t504 + Ifges(7,6) * t566;
t491 = Ifges(7,1) * t505 + Ifges(7,4) * t504 + Ifges(7,5) * t566;
t599 = mrSges(7,1) * t465 - mrSges(7,2) * t466 + Ifges(7,5) * t488 + Ifges(7,6) * t487 + Ifges(7,3) * t544 + t505 * t490 - t504 * t491;
t486 = (-t586 * t624 - qJDD(2)) * pkin(3) + t497 + t633;
t550 = (mrSges(5,1) * t586 - mrSges(5,2) * t589) * qJD(1);
t597 = m(5) * t486 + qJDD(2) * mrSges(5,2) + qJD(2) * t561 - t550 * t620 + t623;
t563 = mrSges(4,2) * t619 + qJD(2) * mrSges(4,3);
t595 = m(4) * t497 - qJDD(2) * mrSges(4,1) - qJD(2) * t563 + t597;
t463 = m(6) * t479 - t512 * mrSges(6,1) + t513 * mrSges(6,2) - t540 * t510 - t541 * t511 + t601;
t485 = -0.2e1 * t616 + t631 + t635;
t594 = -m(5) * t485 + qJDD(2) * mrSges(5,1) - t552 * mrSges(5,3) + qJD(2) * t558 + t463;
t496 = t571 - t631 + t634;
t547 = (-mrSges(4,1) * t589 - mrSges(4,3) * t586) * qJD(1);
t593 = m(4) * t496 + qJDD(2) * mrSges(4,3) + qJD(2) * t560 + t547 * t619 + t594;
t562 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t619;
t559 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t620;
t548 = (-t589 * mrSges(3,1) + t586 * mrSges(3,2)) * qJD(1);
t502 = -Ifges(6,1) * t541 + Ifges(6,4) * t540 + Ifges(6,5) * t620;
t501 = -Ifges(6,4) * t541 + Ifges(6,2) * t540 + Ifges(6,6) * t620;
t500 = -Ifges(6,5) * t541 + Ifges(6,6) * t540 + Ifges(6,3) * t620;
t489 = Ifges(7,5) * t505 + Ifges(7,6) * t504 + Ifges(7,3) * t566;
t457 = mrSges(7,2) * t472 - mrSges(7,3) * t465 + Ifges(7,1) * t488 + Ifges(7,4) * t487 + Ifges(7,5) * t544 + t489 * t504 - t490 * t566;
t456 = -mrSges(7,1) * t472 + mrSges(7,3) * t466 + Ifges(7,4) * t488 + Ifges(7,2) * t487 + Ifges(7,6) * t544 - t489 * t505 + t491 * t566;
t447 = -t551 * mrSges(5,3) + t597;
t446 = -mrSges(5,2) * t552 + t558 * t620 + t603;
t445 = t547 * t620 + t551 * t629 + t595;
t444 = -mrSges(4,3) * t551 + t630 * t552 + (-t563 * t589 + t622 * t586) * qJD(1) + t600;
t443 = mrSges(6,2) * t479 - mrSges(6,3) * t469 + Ifges(6,1) * t513 + Ifges(6,4) * t512 + Ifges(6,5) * t551 - pkin(8) * t455 - t456 * t585 + t457 * t588 + t500 * t540 - t501 * t620;
t442 = -mrSges(6,1) * t479 + mrSges(6,3) * t470 + Ifges(6,4) * t513 + Ifges(6,2) * t512 + Ifges(6,6) * t551 - pkin(5) * t601 + pkin(8) * t606 + t588 * t456 + t585 * t457 + t541 * t500 + t502 * t620;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t621 - mrSges(2,2) * t604 + t586 * (t599 - t541 * t501 + mrSges(3,2) * t532 - t540 * t502 + Ifges(6,6) * t512 + Ifges(6,5) * t513 - mrSges(3,3) * t507 + mrSges(4,2) * t497 - mrSges(5,3) * t486 - mrSges(4,3) * t494 + mrSges(5,1) * t481 + mrSges(6,1) * t469 - mrSges(6,2) * t470 + pkin(5) * t455 + pkin(4) * t448 - qJ(4) * t447 - qJ(3) * t444 - t612 * t619 + t615 * t552 + (Ifges(6,3) + t638) * t551 + t614 * qJDD(2) + t611 * qJD(2)) + t589 * (-qJ(4) * t594 + t582 * t442 - t583 * t443 - mrSges(3,1) * t532 + mrSges(3,3) * t508 + mrSges(4,2) * t496 + mrSges(5,3) * t485 - mrSges(4,1) * t494 - mrSges(5,2) * t481 + qJ(5) * t448 + pkin(3) * t446 - pkin(2) * t444 + t636 * t552 + t615 * t551 + t613 * qJDD(2) + t610 * qJD(2) + (qJ(4) * t550 * t589 + t586 * t612) * qJD(1)) + pkin(1) * (-m(3) * t532 + (-mrSges(3,2) + mrSges(4,3)) * t551 + (mrSges(3,1) - t630) * t552 + ((t562 + t563) * t589 + (-t559 - t622) * t586) * qJD(1) - t600) + pkin(7) * (t589 * (m(3) * t508 - qJDD(2) * mrSges(3,2) - qJD(2) * t559 + t593 + (mrSges(3,3) + mrSges(4,2)) * t552) + (t548 - t550) * t579 * qJD(1) + (-m(3) * t507 - qJDD(2) * mrSges(3,1) - qJD(2) * t562 + t595 - (-mrSges(3,3) - t629) * t551 + (t547 + t548) * t620) * t586); qJ(3) * t593 - t582 * t443 - t583 * t442 + mrSges(3,1) * t507 - mrSges(3,2) * t508 + mrSges(4,3) * t496 - mrSges(4,1) * t497 - mrSges(5,1) * t485 + mrSges(5,2) * t486 + pkin(4) * t463 - qJ(5) * t623 - pkin(3) * t447 - pkin(2) * t445 + (qJ(3) * mrSges(4,2) + t613) * t552 + t614 * t551 + t637 * qJDD(2) + (-t611 * t586 + (-qJ(3) * t550 - t610) * t589) * qJD(1); t445; t446; t463; t599;];
tauJ  = t1;
