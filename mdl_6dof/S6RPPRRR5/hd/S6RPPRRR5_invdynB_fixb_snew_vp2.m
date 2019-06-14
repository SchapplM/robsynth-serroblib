% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-05-05 15:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:50:12
% EndTime: 2019-05-05 15:50:15
% DurationCPUTime: 2.85s
% Computational Cost: add. (30341->290), mult. (58738->346), div. (0->0), fcn. (34473->8), ass. (0->108)
t622 = 2 * qJD(1);
t590 = sin(qJ(1));
t594 = cos(qJ(1));
t567 = -t594 * g(1) - t590 * g(2);
t621 = qJDD(1) * qJ(2) + (qJD(2) * t622) + t567;
t566 = t590 * g(1) - t594 * g(2);
t595 = qJD(1) ^ 2;
t548 = -qJDD(1) * pkin(1) - t595 * qJ(2) + qJDD(2) - t566;
t600 = qJDD(1) * qJ(3) + (qJD(3) * t622) - t548;
t588 = sin(qJ(5));
t589 = sin(qJ(4));
t592 = cos(qJ(5));
t593 = cos(qJ(4));
t553 = (t588 * t593 + t589 * t592) * qJD(1);
t620 = -m(3) - m(4);
t619 = mrSges(2,1) - mrSges(3,2);
t547 = t595 * pkin(1) - t621;
t544 = qJDD(3) + (-pkin(1) - qJ(3)) * t595 + t621;
t539 = -qJDD(1) * pkin(7) + t544;
t531 = t589 * g(3) + t593 * t539;
t615 = qJD(1) * qJD(4);
t610 = t589 * t615;
t562 = t593 * qJDD(1) - t610;
t514 = (-t562 - t610) * pkin(8) + (-t589 * t593 * t595 + qJDD(4)) * pkin(4) + t531;
t532 = -t593 * g(3) + t589 * t539;
t561 = -t589 * qJDD(1) - t593 * t615;
t616 = qJD(1) * t593;
t565 = qJD(4) * pkin(4) - pkin(8) * t616;
t582 = t589 ^ 2;
t515 = -t582 * t595 * pkin(4) + t561 * pkin(8) - qJD(4) * t565 + t532;
t504 = t588 * t514 + t592 * t515;
t554 = (-t588 * t589 + t592 * t593) * qJD(1);
t522 = -t554 * qJD(5) + t592 * t561 - t588 * t562;
t533 = t553 * mrSges(6,1) + t554 * mrSges(6,2);
t575 = qJD(4) + qJD(5);
t546 = t575 * mrSges(6,1) - t554 * mrSges(6,3);
t574 = qJDD(4) + qJDD(5);
t534 = t553 * pkin(5) - t554 * pkin(9);
t573 = t575 ^ 2;
t501 = -t573 * pkin(5) + t574 * pkin(9) - t553 * t534 + t504;
t517 = -t561 * pkin(4) + t565 * t616 + (-pkin(8) * t582 - pkin(7)) * t595 + t600;
t523 = -t553 * qJD(5) + t588 * t561 + t592 * t562;
t502 = (t553 * t575 - t523) * pkin(9) + (t554 * t575 - t522) * pkin(5) + t517;
t587 = sin(qJ(6));
t591 = cos(qJ(6));
t498 = -t587 * t501 + t591 * t502;
t540 = -t587 * t554 + t591 * t575;
t507 = t540 * qJD(6) + t591 * t523 + t587 * t574;
t541 = t591 * t554 + t587 * t575;
t518 = -t540 * mrSges(7,1) + t541 * mrSges(7,2);
t521 = qJDD(6) - t522;
t549 = qJD(6) + t553;
t524 = -t549 * mrSges(7,2) + t540 * mrSges(7,3);
t496 = m(7) * t498 + t521 * mrSges(7,1) - t507 * mrSges(7,3) - t541 * t518 + t549 * t524;
t499 = t591 * t501 + t587 * t502;
t506 = -t541 * qJD(6) - t587 * t523 + t591 * t574;
t525 = t549 * mrSges(7,1) - t541 * mrSges(7,3);
t497 = m(7) * t499 - t521 * mrSges(7,2) + t506 * mrSges(7,3) + t540 * t518 - t549 * t525;
t606 = -t587 * t496 + t591 * t497;
t487 = m(6) * t504 - t574 * mrSges(6,2) + t522 * mrSges(6,3) - t553 * t533 - t575 * t546 + t606;
t503 = t592 * t514 - t588 * t515;
t545 = -t575 * mrSges(6,2) - t553 * mrSges(6,3);
t500 = -t574 * pkin(5) - t573 * pkin(9) + t554 * t534 - t503;
t598 = -m(7) * t500 + t506 * mrSges(7,1) - t507 * mrSges(7,2) + t540 * t524 - t541 * t525;
t492 = m(6) * t503 + t574 * mrSges(6,1) - t523 * mrSges(6,3) - t554 * t533 + t575 * t545 + t598;
t481 = t588 * t487 + t592 * t492;
t560 = (mrSges(5,1) * t589 + mrSges(5,2) * t593) * qJD(1);
t617 = qJD(1) * t589;
t563 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t617;
t479 = m(5) * t531 + qJDD(4) * mrSges(5,1) - t562 * mrSges(5,3) + qJD(4) * t563 - t560 * t616 + t481;
t564 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t616;
t607 = t592 * t487 - t588 * t492;
t480 = m(5) * t532 - qJDD(4) * mrSges(5,2) + t561 * mrSges(5,3) - qJD(4) * t564 - t560 * t617 + t607;
t473 = t593 * t479 + t589 * t480;
t605 = -m(4) * t544 - qJDD(1) * mrSges(4,2) - t473;
t599 = -m(3) * t547 + (t595 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t605;
t471 = m(2) * t567 - qJDD(1) * mrSges(2,2) + ((-mrSges(2,1) - mrSges(4,3)) * t595) + t599;
t538 = -t595 * pkin(7) + t600;
t488 = t591 * t496 + t587 * t497;
t601 = m(6) * t517 - t522 * mrSges(6,1) + t523 * mrSges(6,2) + t553 * t545 + t554 * t546 + t488;
t597 = -m(5) * t538 + t561 * mrSges(5,1) - t562 * mrSges(5,2) - t563 * t617 - t564 * t616 - t601;
t484 = -m(4) * t600 - t595 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t597;
t596 = -m(3) * t548 + t595 * mrSges(3,3) - t484;
t483 = m(2) * t566 - t595 * mrSges(2,2) + t619 * qJDD(1) + t596;
t618 = t590 * t471 + t594 * t483;
t612 = (Ifges(2,5) - Ifges(3,4) + Ifges(4,5));
t611 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t609 = t594 * t471 - t590 * t483;
t608 = -t589 * t479 + t593 * t480;
t552 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t593 - Ifges(5,4) * t589) * qJD(1);
t551 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t593 - Ifges(5,2) * t589) * qJD(1);
t550 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t593 - Ifges(5,6) * t589) * qJD(1);
t528 = Ifges(6,1) * t554 - Ifges(6,4) * t553 + Ifges(6,5) * t575;
t527 = Ifges(6,4) * t554 - Ifges(6,2) * t553 + Ifges(6,6) * t575;
t526 = Ifges(6,5) * t554 - Ifges(6,6) * t553 + Ifges(6,3) * t575;
t510 = Ifges(7,1) * t541 + Ifges(7,4) * t540 + Ifges(7,5) * t549;
t509 = Ifges(7,4) * t541 + Ifges(7,2) * t540 + Ifges(7,6) * t549;
t508 = Ifges(7,5) * t541 + Ifges(7,6) * t540 + Ifges(7,3) * t549;
t490 = mrSges(7,2) * t500 - mrSges(7,3) * t498 + Ifges(7,1) * t507 + Ifges(7,4) * t506 + Ifges(7,5) * t521 + t540 * t508 - t549 * t509;
t489 = -mrSges(7,1) * t500 + mrSges(7,3) * t499 + Ifges(7,4) * t507 + Ifges(7,2) * t506 + Ifges(7,6) * t521 - t541 * t508 + t549 * t510;
t475 = -mrSges(6,1) * t517 - mrSges(7,1) * t498 + mrSges(7,2) * t499 + mrSges(6,3) * t504 + Ifges(6,4) * t523 - Ifges(7,5) * t507 + Ifges(6,2) * t522 + Ifges(6,6) * t574 - Ifges(7,6) * t506 - Ifges(7,3) * t521 - pkin(5) * t488 - t541 * t509 + t540 * t510 - t554 * t526 + t575 * t528;
t474 = mrSges(6,2) * t517 - mrSges(6,3) * t503 + Ifges(6,1) * t523 + Ifges(6,4) * t522 + Ifges(6,5) * t574 - pkin(9) * t488 - t587 * t489 + t591 * t490 - t553 * t526 - t575 * t527;
t472 = t620 * g(3) + t608;
t468 = mrSges(5,2) * t538 - mrSges(5,3) * t531 + Ifges(5,1) * t562 + Ifges(5,4) * t561 + Ifges(5,5) * qJDD(4) - pkin(8) * t481 - qJD(4) * t551 + t592 * t474 - t588 * t475 - t550 * t617;
t467 = -mrSges(5,1) * t538 + mrSges(5,3) * t532 + Ifges(5,4) * t562 + Ifges(5,2) * t561 + Ifges(5,6) * qJDD(4) - pkin(4) * t601 + pkin(8) * t607 + qJD(4) * t552 + t588 * t474 + t592 * t475 - t550 * t616;
t466 = (qJ(3) * m(4) + mrSges(4,3) + t619) * g(3) + t611 * qJDD(1) + ((-pkin(2) * mrSges(4,3) + t612) * t595) - pkin(2) * t605 + pkin(9) * t606 - qJ(3) * t608 + pkin(5) * t598 + (t593 * t551 + t589 * t552) * qJD(1) + Ifges(5,3) * qJDD(4) + t587 * t490 + t591 * t489 + mrSges(2,3) * t567 + Ifges(6,3) * t574 + t553 * t528 + t554 * t527 + Ifges(5,6) * t561 + Ifges(5,5) * t562 + mrSges(4,1) * t544 - mrSges(3,1) * t547 + mrSges(5,1) * t531 - mrSges(5,2) * t532 + Ifges(6,6) * t522 + Ifges(6,5) * t523 + mrSges(6,1) * t503 - mrSges(6,2) * t504 + pkin(4) * t481 + pkin(3) * t473 - pkin(1) * t472;
t465 = -qJ(2) * t472 - mrSges(2,3) * t566 + pkin(2) * t484 + mrSges(3,1) * t548 + t589 * t468 + t593 * t467 + pkin(3) * t597 + pkin(7) * t608 - mrSges(4,1) * t600 - t611 * t595 + t612 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t609; -m(1) * g(2) + t618; (-m(1) - m(2) + t620) * g(3) + t608; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t618 + t594 * t465 - t590 * t466; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t609 + t590 * t465 + t594 * t466; qJ(2) * (-(t595 * mrSges(4,3)) + t599) + pkin(1) * t596 + mrSges(2,1) * t566 - mrSges(2,2) * t567 - qJ(3) * t484 + mrSges(3,2) * t548 - mrSges(3,3) * t547 + t593 * t468 - t589 * t467 - pkin(7) * t473 + mrSges(4,2) * t544 + mrSges(4,3) * t600 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
