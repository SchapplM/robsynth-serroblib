% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRRRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRRRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:26
% EndTime: 2019-12-05 18:58:29
% DurationCPUTime: 3.33s
% Computational Cost: add. (67969->233), mult. (69973->293), div. (0->0), fcn. (39543->10), ass. (0->103)
t579 = qJDD(1) + qJDD(2);
t572 = qJDD(3) + t579;
t585 = sin(qJ(4));
t590 = cos(qJ(4));
t581 = qJD(1) + qJD(2);
t573 = qJD(3) + t581;
t610 = qJD(4) * t573;
t555 = t585 * t572 + t590 * t610;
t588 = sin(qJ(1));
t593 = cos(qJ(1));
t567 = t593 * g(2) + t588 * g(3);
t563 = qJDD(1) * pkin(1) + t567;
t566 = t588 * g(2) - t593 * g(3);
t594 = qJD(1) ^ 2;
t564 = -t594 * pkin(1) + t566;
t587 = sin(qJ(2));
t592 = cos(qJ(2));
t543 = t592 * t563 - t587 * t564;
t540 = t579 * pkin(2) + t543;
t544 = t587 * t563 + t592 * t564;
t577 = t581 ^ 2;
t541 = -t577 * pkin(2) + t544;
t586 = sin(qJ(3));
t591 = cos(qJ(3));
t525 = t586 * t540 + t591 * t541;
t571 = t573 ^ 2;
t522 = -t571 * pkin(3) + t572 * pkin(8) + t525;
t611 = t585 * t522;
t614 = pkin(4) * t571;
t515 = qJDD(4) * pkin(4) - t555 * pkin(9) - t611 + (pkin(9) * t610 + t585 * t614 - g(1)) * t590;
t519 = -t585 * g(1) + t590 * t522;
t556 = t590 * t572 - t585 * t610;
t613 = t573 * t585;
t562 = qJD(4) * pkin(4) - pkin(9) * t613;
t583 = t590 ^ 2;
t516 = t556 * pkin(9) - qJD(4) * t562 - t583 * t614 + t519;
t584 = sin(qJ(5));
t589 = cos(qJ(5));
t513 = t589 * t515 - t584 * t516;
t550 = (-t584 * t585 + t589 * t590) * t573;
t531 = t550 * qJD(5) + t589 * t555 + t584 * t556;
t551 = (t584 * t590 + t585 * t589) * t573;
t536 = -t550 * mrSges(6,1) + t551 * mrSges(6,2);
t580 = qJD(4) + qJD(5);
t545 = -t580 * mrSges(6,2) + t550 * mrSges(6,3);
t578 = qJDD(4) + qJDD(5);
t510 = m(6) * t513 + t578 * mrSges(6,1) - t531 * mrSges(6,3) - t551 * t536 + t580 * t545;
t514 = t584 * t515 + t589 * t516;
t530 = -t551 * qJD(5) - t584 * t555 + t589 * t556;
t546 = t580 * mrSges(6,1) - t551 * mrSges(6,3);
t511 = m(6) * t514 - t578 * mrSges(6,2) + t530 * mrSges(6,3) + t550 * t536 - t580 * t546;
t501 = t589 * t510 + t584 * t511;
t518 = -t590 * g(1) - t611;
t548 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t585 + Ifges(5,2) * t590) * t573;
t549 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t585 + Ifges(5,4) * t590) * t573;
t533 = Ifges(6,4) * t551 + Ifges(6,2) * t550 + Ifges(6,6) * t580;
t534 = Ifges(6,1) * t551 + Ifges(6,4) * t550 + Ifges(6,5) * t580;
t599 = -mrSges(6,1) * t513 + mrSges(6,2) * t514 - Ifges(6,5) * t531 - Ifges(6,6) * t530 - Ifges(6,3) * t578 - t551 * t533 + t550 * t534;
t616 = mrSges(5,1) * t518 - mrSges(5,2) * t519 + Ifges(5,5) * t555 + Ifges(5,6) * t556 + Ifges(5,3) * qJDD(4) + pkin(4) * t501 + (t585 * t548 - t590 * t549) * t573 - t599;
t615 = -m(3) - m(4);
t612 = t573 * t590;
t554 = (-mrSges(5,1) * t590 + mrSges(5,2) * t585) * t573;
t561 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t612;
t499 = m(5) * t518 + qJDD(4) * mrSges(5,1) - t555 * mrSges(5,3) + qJD(4) * t561 - t554 * t613 + t501;
t560 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t613;
t605 = -t584 * t510 + t589 * t511;
t500 = m(5) * t519 - qJDD(4) * mrSges(5,2) + t556 * mrSges(5,3) - qJD(4) * t560 + t554 * t612 + t605;
t606 = -t585 * t499 + t590 * t500;
t493 = m(4) * t525 - t571 * mrSges(4,1) - t572 * mrSges(4,2) + t606;
t524 = t591 * t540 - t586 * t541;
t602 = -t572 * pkin(3) - t524;
t521 = -t571 * pkin(8) + t602;
t517 = t562 * t613 - t556 * pkin(4) + (-pkin(9) * t583 - pkin(8)) * t571 + t602;
t600 = m(6) * t517 - t530 * mrSges(6,1) + t531 * mrSges(6,2) - t550 * t545 + t551 * t546;
t596 = -m(5) * t521 + t556 * mrSges(5,1) - t555 * mrSges(5,2) - t560 * t613 + t561 * t612 - t600;
t505 = m(4) * t524 + t572 * mrSges(4,1) - t571 * mrSges(4,2) + t596;
t488 = t586 * t493 + t591 * t505;
t485 = m(3) * t543 + t579 * mrSges(3,1) - t577 * mrSges(3,2) + t488;
t607 = t591 * t493 - t586 * t505;
t486 = m(3) * t544 - t577 * mrSges(3,1) - t579 * mrSges(3,2) + t607;
t478 = t592 * t485 + t587 * t486;
t495 = t590 * t499 + t585 * t500;
t608 = -t587 * t485 + t592 * t486;
t475 = m(2) * t566 - t594 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t608;
t476 = m(2) * t567 + qJDD(1) * mrSges(2,1) - t594 * mrSges(2,2) + t478;
t609 = t593 * t475 - t588 * t476;
t604 = -t588 * t475 - t593 * t476;
t532 = Ifges(6,5) * t551 + Ifges(6,6) * t550 + Ifges(6,3) * t580;
t502 = -mrSges(6,1) * t517 + mrSges(6,3) * t514 + Ifges(6,4) * t531 + Ifges(6,2) * t530 + Ifges(6,6) * t578 - t551 * t532 + t580 * t534;
t503 = mrSges(6,2) * t517 - mrSges(6,3) * t513 + Ifges(6,1) * t531 + Ifges(6,4) * t530 + Ifges(6,5) * t578 + t550 * t532 - t580 * t533;
t547 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t585 + Ifges(5,6) * t590) * t573;
t481 = -mrSges(5,1) * t521 + mrSges(5,3) * t519 + Ifges(5,4) * t555 + Ifges(5,2) * t556 + Ifges(5,6) * qJDD(4) - pkin(4) * t600 + pkin(9) * t605 + qJD(4) * t549 + t589 * t502 + t584 * t503 - t547 * t613;
t490 = mrSges(5,2) * t521 - mrSges(5,3) * t518 + Ifges(5,1) * t555 + Ifges(5,4) * t556 + Ifges(5,5) * qJDD(4) - pkin(9) * t501 - qJD(4) * t548 - t584 * t502 + t589 * t503 + t547 * t612;
t601 = mrSges(4,1) * t524 - mrSges(4,2) * t525 + Ifges(4,3) * t572 + pkin(3) * t596 + pkin(8) * t606 + t590 * t481 + t585 * t490;
t598 = mrSges(3,1) * t543 - mrSges(3,2) * t544 + Ifges(3,3) * t579 + pkin(2) * t488 + t601;
t597 = mrSges(2,1) * t567 - mrSges(2,2) * t566 + Ifges(2,3) * qJDD(1) + pkin(1) * t478 + t598;
t479 = mrSges(4,1) * g(1) + mrSges(4,3) * t525 + t571 * Ifges(4,5) + Ifges(4,6) * t572 - pkin(3) * t495 - t616;
t473 = -mrSges(4,2) * g(1) - mrSges(4,3) * t524 + Ifges(4,5) * t572 - t571 * Ifges(4,6) - pkin(8) * t495 - t585 * t481 + t590 * t490;
t472 = -mrSges(3,2) * g(1) - mrSges(3,3) * t543 + Ifges(3,5) * t579 - t577 * Ifges(3,6) - pkin(7) * t488 + t591 * t473 - t586 * t479;
t471 = Ifges(3,6) * t579 + t577 * Ifges(3,5) + mrSges(3,1) * g(1) + mrSges(3,3) * t544 + t586 * t473 + t591 * t479 - pkin(2) * (-m(4) * g(1) + t495) + pkin(7) * t607;
t470 = -mrSges(2,2) * g(1) - mrSges(2,3) * t567 + Ifges(2,5) * qJDD(1) - t594 * Ifges(2,6) - pkin(6) * t478 - t587 * t471 + t592 * t472;
t469 = Ifges(2,6) * qJDD(1) + t594 * Ifges(2,5) + mrSges(2,3) * t566 + t587 * t472 + t592 * t471 - pkin(1) * t495 + pkin(6) * t608 + (-pkin(1) * t615 + mrSges(2,1)) * g(1);
t1 = [(-m(1) - m(2) + t615) * g(1) + t495; -m(1) * g(2) + t604; -m(1) * g(3) + t609; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t597; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t609 - t593 * t469 - t588 * t470; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t604 - t588 * t469 + t593 * t470; t597; t598; t601; t616; -t599;];
tauJB = t1;
