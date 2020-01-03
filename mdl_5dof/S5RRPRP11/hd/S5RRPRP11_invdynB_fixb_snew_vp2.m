% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRP11
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRP11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP11_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP11_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP11_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP11_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP11_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:12:34
% EndTime: 2019-12-31 20:12:37
% DurationCPUTime: 1.84s
% Computational Cost: add. (12754->271), mult. (25762->314), div. (0->0), fcn. (13338->6), ass. (0->106)
t621 = -2 * qJD(3);
t620 = Ifges(3,1) + Ifges(4,2);
t619 = Ifges(5,1) + Ifges(6,1);
t612 = Ifges(3,4) + Ifges(4,6);
t611 = Ifges(5,4) - Ifges(6,5);
t610 = Ifges(6,4) + Ifges(5,5);
t609 = Ifges(3,5) - Ifges(4,4);
t618 = Ifges(3,2) + Ifges(4,3);
t617 = Ifges(5,2) + Ifges(6,3);
t608 = Ifges(3,6) - Ifges(4,5);
t607 = Ifges(5,6) - Ifges(6,6);
t616 = (Ifges(3,3) + Ifges(4,1));
t615 = Ifges(5,3) + Ifges(6,2);
t574 = sin(qJ(1));
t577 = cos(qJ(1));
t560 = -t577 * g(1) - t574 * g(2);
t579 = qJD(1) ^ 2;
t535 = -t579 * pkin(1) + qJDD(1) * pkin(6) + t560;
t573 = sin(qJ(2));
t576 = cos(qJ(2));
t522 = -t573 * g(3) + t576 * t535;
t547 = (-pkin(2) * t576 - qJ(3) * t573) * qJD(1);
t578 = qJD(2) ^ 2;
t597 = qJD(1) * t576;
t498 = t578 * pkin(2) - qJDD(2) * qJ(3) + (qJD(2) * t621) - t547 * t597 - t522;
t614 = t579 * pkin(6);
t613 = -mrSges(5,3) - mrSges(6,2);
t521 = -t576 * g(3) - t573 * t535;
t548 = (mrSges(4,2) * t576 - mrSges(4,3) * t573) * qJD(1);
t549 = (-mrSges(3,1) * t576 + mrSges(3,2) * t573) * qJD(1);
t595 = qJD(1) * qJD(2);
t591 = t576 * t595;
t550 = t573 * qJDD(1) + t591;
t555 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t597;
t556 = -mrSges(4,1) * t597 - qJD(2) * mrSges(4,3);
t592 = t573 * t595;
t551 = t576 * qJDD(1) - t592;
t596 = t573 * qJD(1);
t558 = pkin(3) * t596 - (qJD(2) * pkin(7));
t571 = t576 ^ 2;
t559 = t574 * g(1) - t577 * g(2);
t587 = -qJDD(1) * pkin(1) - t559;
t583 = pkin(2) * t592 + t596 * t621 + (-t550 - t591) * qJ(3) + t587;
t492 = -t558 * t596 + (-pkin(3) * t571 - pkin(6)) * t579 + (-pkin(2) - pkin(7)) * t551 + t583;
t499 = -qJDD(2) * pkin(2) - t578 * qJ(3) + t547 * t596 + qJDD(3) - t521;
t496 = (-t573 * t576 * t579 - qJDD(2)) * pkin(7) + (t550 - t591) * pkin(3) + t499;
t572 = sin(qJ(4));
t575 = cos(qJ(4));
t490 = t575 * t492 + t572 * t496;
t546 = t575 * qJD(2) - t572 * t597;
t510 = t546 * qJD(4) + t572 * qJDD(2) + t575 * t551;
t563 = qJD(4) + t596;
t519 = t563 * mrSges(5,1) - t546 * mrSges(5,3);
t544 = qJDD(4) + t550;
t545 = t572 * qJD(2) + t575 * t597;
t514 = t545 * pkin(4) - t546 * qJ(5);
t561 = t563 ^ 2;
t485 = -t561 * pkin(4) + t544 * qJ(5) + 0.2e1 * qJD(5) * t563 - t545 * t514 + t490;
t520 = -t563 * mrSges(6,1) + t546 * mrSges(6,2);
t593 = m(6) * t485 + t544 * mrSges(6,3) + t563 * t520;
t515 = t545 * mrSges(6,1) - t546 * mrSges(6,3);
t601 = -t545 * mrSges(5,1) - t546 * mrSges(5,2) - t515;
t480 = m(5) * t490 - t544 * mrSges(5,2) + t613 * t510 - t563 * t519 + t601 * t545 + t593;
t489 = -t572 * t492 + t575 * t496;
t511 = -t545 * qJD(4) + t575 * qJDD(2) - t572 * t551;
t517 = -t563 * mrSges(5,2) - t545 * mrSges(5,3);
t486 = -t544 * pkin(4) - t561 * qJ(5) + t546 * t514 + qJDD(5) - t489;
t518 = -t545 * mrSges(6,2) + t563 * mrSges(6,3);
t588 = -m(6) * t486 + t544 * mrSges(6,1) + t563 * t518;
t482 = m(5) * t489 + t544 * mrSges(5,1) + t613 * t511 + t563 * t517 + t601 * t546 + t588;
t475 = t572 * t480 + t575 * t482;
t584 = -m(4) * t499 - t550 * mrSges(4,1) - t475;
t471 = m(3) * t521 - t550 * mrSges(3,3) + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t555 - t556) * qJD(2) + (-t548 - t549) * t596 + t584;
t554 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t596;
t557 = mrSges(4,1) * t596 + qJD(2) * mrSges(4,2);
t495 = -t571 * t579 * pkin(7) + t551 * pkin(3) + qJD(2) * t558 - t498;
t488 = -0.2e1 * qJD(5) * t546 + (t545 * t563 - t511) * qJ(5) + (t546 * t563 + t510) * pkin(4) + t495;
t483 = m(6) * t488 + t510 * mrSges(6,1) - t511 * mrSges(6,3) + t545 * t518 - t546 * t520;
t582 = m(5) * t495 + t510 * mrSges(5,1) + t511 * mrSges(5,2) + t545 * t517 + t546 * t519 + t483;
t581 = -m(4) * t498 + qJDD(2) * mrSges(4,3) + qJD(2) * t557 + t548 * t597 + t582;
t478 = t549 * t597 + m(3) * t522 + (mrSges(3,3) + mrSges(4,1)) * t551 + t581 - qJDD(2) * mrSges(3,2) - qJD(2) * t554;
t589 = -t573 * t471 + t576 * t478;
t466 = m(2) * t560 - t579 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t589;
t534 = t587 - t614;
t497 = -t551 * pkin(2) + t583 - t614;
t605 = t575 * t480 - t572 * t482;
t586 = -m(4) * t497 - t551 * mrSges(4,2) + t557 * t596 - t605;
t580 = -m(3) * t534 + t555 * t597 + t551 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t550 + (-t554 * t573 - t556 * t576) * qJD(1) + t586;
t469 = m(2) * t559 + qJDD(1) * mrSges(2,1) - t579 * mrSges(2,2) + t580;
t606 = t574 * t466 + t577 * t469;
t467 = t576 * t471 + t573 * t478;
t604 = t617 * t545 - t611 * t546 - t607 * t563;
t603 = t607 * t545 - t610 * t546 - t615 * t563;
t602 = -t611 * t545 + t619 * t546 + t610 * t563;
t600 = (t616 * qJD(2)) + (t609 * t573 + t608 * t576) * qJD(1);
t599 = -t608 * qJD(2) + (-t612 * t573 - t618 * t576) * qJD(1);
t598 = t609 * qJD(2) + (t620 * t573 + t612 * t576) * qJD(1);
t590 = t577 * t466 - t574 * t469;
t474 = mrSges(5,2) * t495 + mrSges(6,2) * t486 - mrSges(5,3) * t489 - mrSges(6,3) * t488 - qJ(5) * t483 - t611 * t510 + t619 * t511 + t610 * t544 + t603 * t545 + t604 * t563;
t473 = -mrSges(5,1) * t495 - mrSges(6,1) * t488 + mrSges(6,2) * t485 + mrSges(5,3) * t490 - pkin(4) * t483 - t617 * t510 + t611 * t511 + t607 * t544 + t603 * t546 + t602 * t563;
t472 = -t550 * mrSges(4,3) + t556 * t597 - t586;
t463 = -qJ(3) * t472 + pkin(3) * t475 - mrSges(3,3) * t521 + mrSges(3,2) * t534 + qJ(5) * t593 + pkin(4) * t588 + mrSges(6,3) * t485 - mrSges(6,1) * t486 + mrSges(5,1) * t489 - mrSges(5,2) * t490 - mrSges(4,3) * t497 + mrSges(4,1) * t499 + t612 * t551 + t620 * t550 + (-pkin(4) * t515 - t604) * t546 + (-qJ(5) * t515 + t602) * t545 + t615 * t544 + (-pkin(4) * mrSges(6,2) + t610) * t511 + (-qJ(5) * mrSges(6,2) - t607) * t510 + t609 * qJDD(2) + t599 * qJD(2) + t600 * t597;
t462 = -mrSges(3,1) * t534 - mrSges(4,1) * t498 + mrSges(4,2) * t497 + mrSges(3,3) * t522 - pkin(2) * t472 + pkin(3) * t582 - pkin(7) * t605 + t598 * qJD(2) + t608 * qJDD(2) - t575 * t473 - t572 * t474 + t612 * t550 + t618 * t551 - t600 * t596;
t461 = -pkin(1) * t467 + mrSges(2,3) * t560 - mrSges(3,1) * t521 + mrSges(3,2) * t522 - pkin(2) * (-qJD(2) * t556 + t584) - qJ(3) * t581 - mrSges(4,2) * t499 + mrSges(4,3) * t498 - t575 * t474 + t572 * t473 + pkin(7) * t475 + mrSges(2,1) * g(3) + t579 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-qJ(3) * mrSges(4,1) - t608) * t551 - t609 * t550 + (pkin(2) * mrSges(4,2) - t616) * qJDD(2) + (t598 * t576 + (pkin(2) * t548 + t599) * t573) * qJD(1);
t460 = -mrSges(2,2) * g(3) - mrSges(2,3) * t559 + Ifges(2,5) * qJDD(1) - t579 * Ifges(2,6) - pkin(6) * t467 - t573 * t462 + t576 * t463;
t1 = [-m(1) * g(1) + t590; -m(1) * g(2) + t606; (-m(1) - m(2)) * g(3) + t467; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t606 + t577 * t460 - t574 * t461; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t590 + t574 * t460 + t577 * t461; -mrSges(1,1) * g(2) + mrSges(2,1) * t559 + mrSges(1,2) * g(1) - mrSges(2,2) * t560 + Ifges(2,3) * qJDD(1) + pkin(1) * t580 + pkin(6) * t589 + t576 * t462 + t573 * t463;];
tauB = t1;
