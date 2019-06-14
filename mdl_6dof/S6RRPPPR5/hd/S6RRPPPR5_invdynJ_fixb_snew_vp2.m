% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2019-05-06 08:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:52:00
% EndTime: 2019-05-06 08:52:04
% DurationCPUTime: 2.55s
% Computational Cost: add. (12902->316), mult. (28458->366), div. (0->0), fcn. (17557->8), ass. (0->128)
t544 = sin(pkin(9));
t546 = sin(qJ(2));
t581 = t546 * qJD(1);
t591 = cos(pkin(9));
t526 = t544 * qJD(2) + t591 * t581;
t602 = -0.2e1 * t526;
t601 = 2 * qJD(4);
t600 = Ifges(4,1) + Ifges(6,3) + Ifges(5,2);
t599 = Ifges(5,1) + Ifges(6,2) + Ifges(4,3);
t579 = -Ifges(4,5) + Ifges(6,6) + Ifges(5,4);
t578 = -Ifges(4,6) + Ifges(5,5) - Ifges(6,4);
t577 = -Ifges(5,6) + Ifges(6,5) - Ifges(4,4);
t598 = -Ifges(5,3) - Ifges(6,1) - Ifges(4,2);
t552 = qJD(1) ^ 2;
t547 = sin(qJ(1));
t550 = cos(qJ(1));
t566 = -g(1) * t550 - g(2) * t547;
t521 = -pkin(1) * t552 + qJDD(1) * pkin(7) + t566;
t549 = cos(qJ(2));
t491 = -t549 * g(3) - t546 * t521;
t529 = (-t549 * pkin(2) - t546 * qJ(3)) * qJD(1);
t551 = qJD(2) ^ 2;
t559 = qJDD(2) * pkin(2) + qJ(3) * t551 - t529 * t581 - qJDD(3) + t491;
t525 = -t591 * qJD(2) + t544 * t581;
t584 = qJD(1) * t549;
t572 = t525 * t584;
t597 = -qJ(4) * t572 - t559;
t487 = pkin(3) * t525 - qJ(4) * t526;
t580 = qJD(1) * qJD(2);
t538 = t546 * t580;
t532 = qJDD(1) * t549 - t538;
t568 = g(1) * t547 - t550 * g(2);
t520 = -qJDD(1) * pkin(1) - pkin(7) * t552 - t568;
t570 = t549 * t580;
t531 = qJDD(1) * t546 + t570;
t459 = (-t531 - t570) * qJ(3) + (-t532 + t538) * pkin(2) + t520;
t492 = -g(3) * t546 + t549 * t521;
t463 = -pkin(2) * t551 + qJDD(2) * qJ(3) + t529 * t584 + t492;
t586 = t544 * t459 + t591 * t463;
t588 = t549 ^ 2 * t552;
t596 = pkin(3) * t588 + t532 * qJ(4) + t525 * t487 + t584 * t601 - t586;
t446 = qJD(3) * t602 + t591 * t459 - t544 * t463;
t505 = mrSges(4,2) * t584 - mrSges(4,3) * t525;
t507 = -t591 * qJDD(2) + t531 * t544;
t508 = t544 * qJDD(2) + t591 * t531;
t571 = t526 * t584;
t554 = qJD(4) * t602 + (t507 - t571) * pkin(3) + t597;
t590 = qJ(4) * t508;
t445 = t554 - t590;
t502 = mrSges(5,1) * t525 + mrSges(5,3) * t584;
t500 = pkin(4) * t526 + qJ(5) * t584;
t522 = t525 ^ 2;
t594 = pkin(4) * t522;
t442 = t594 + t590 - 0.2e1 * qJD(5) * t525 + (-pkin(3) - qJ(5)) * t507 + (pkin(3) * t584 + t500 + t601) * t526 - t597;
t501 = -mrSges(6,2) * t526 + mrSges(6,3) * t584;
t504 = -mrSges(6,1) * t584 + mrSges(6,2) * t525;
t509 = -pkin(5) * t584 - pkin(8) * t525;
t523 = t526 ^ 2;
t595 = 0.2e1 * qJD(5);
t437 = -t594 + qJ(5) * t507 - t500 * t526 - pkin(8) * t523 + (-pkin(5) - qJ(4)) * t508 + (t595 + t509) * t525 + t554;
t545 = sin(qJ(6));
t548 = cos(qJ(6));
t485 = t525 * t548 + t526 * t545;
t450 = -qJD(6) * t485 - t507 * t545 + t508 * t548;
t484 = -t525 * t545 + t526 * t548;
t451 = qJD(6) * t484 + t507 * t548 + t508 * t545;
t535 = qJD(6) - t584;
t464 = -mrSges(7,2) * t535 + mrSges(7,3) * t484;
t465 = mrSges(7,1) * t535 - mrSges(7,3) * t485;
t563 = -m(7) * t437 + t450 * mrSges(7,1) - t451 * mrSges(7,2) + t484 * t464 - t485 * t465;
t558 = -m(6) * t442 - t508 * mrSges(6,1) + t507 * mrSges(6,3) - t526 * t501 + t525 * t504 - t563;
t556 = m(5) * t445 - t507 * mrSges(5,2) - t525 * t502 + t558;
t503 = mrSges(5,1) * t526 - mrSges(5,2) * t584;
t585 = mrSges(4,1) * t584 + mrSges(4,3) * t526 + t503;
t422 = -m(4) * t559 + t507 * mrSges(4,1) + t525 * t505 + (mrSges(4,2) - mrSges(5,3)) * t508 - t585 * t526 + t556;
t592 = -mrSges(5,2) + mrSges(6,3);
t589 = t525 * t526;
t444 = t532 * pkin(3) - qJ(4) * t588 + t526 * t487 + qJDD(4) - t446;
t438 = t584 * t595 + (t532 + t589) * qJ(5) + (t508 - t572) * pkin(4) + t444;
t434 = -t523 * pkin(5) + t508 * pkin(8) + t509 * t584 + t438;
t583 = qJD(3) * t525;
t515 = -0.2e1 * t583;
t440 = -pkin(4) * t507 - qJ(5) * t522 - t500 * t584 + qJDD(5) + t515 - t596;
t435 = (-t532 + t589) * pkin(5) + (-t507 - t571) * pkin(8) + t440;
t432 = -t434 * t545 + t435 * t548;
t457 = -mrSges(7,1) * t484 + mrSges(7,2) * t485;
t528 = qJDD(6) - t532;
t430 = m(7) * t432 + mrSges(7,1) * t528 - mrSges(7,3) * t451 - t457 * t485 + t464 * t535;
t433 = t434 * t548 + t435 * t545;
t431 = m(7) * t433 - mrSges(7,2) * t528 + mrSges(7,3) * t450 + t457 * t484 - t465 * t535;
t421 = t548 * t430 + t545 * t431;
t587 = -t545 * t430 + t548 * t431;
t576 = -t525 * t577 - t526 * t600 - t579 * t584;
t575 = -t525 * t578 + t526 * t579 + t584 * t599;
t574 = -t525 * t598 + t526 * t577 - t578 * t584;
t573 = t504 * t584;
t488 = mrSges(4,1) * t525 + mrSges(4,2) * t526;
t489 = -mrSges(5,2) * t525 - mrSges(5,3) * t526;
t486 = mrSges(6,1) * t526 - mrSges(6,3) * t525;
t562 = -m(6) * t438 + t508 * mrSges(6,2) + t526 * t486 - t587;
t557 = m(5) * t444 + t508 * mrSges(5,1) + t526 * t489 - t502 * t584 - t562;
t415 = m(4) * t446 - mrSges(4,3) * t508 - t488 * t526 + (-t504 - t505) * t584 + (-mrSges(4,1) - t592) * t532 - t557;
t447 = t515 + t586;
t443 = 0.2e1 * t583 + t596;
t564 = m(6) * t440 + t507 * mrSges(6,2) + t525 * t486 + t421;
t560 = -m(5) * t443 - t532 * mrSges(5,3) + t564;
t417 = m(4) * t447 + (-mrSges(6,1) + mrSges(4,2)) * t532 + (-t488 - t489) * t525 + (-mrSges(4,3) - mrSges(5,1)) * t507 + (-t501 - t585) * t584 + t560;
t567 = -t415 * t544 + t591 * t417;
t414 = t591 * t415 + t544 * t417;
t453 = Ifges(7,4) * t485 + Ifges(7,2) * t484 + Ifges(7,6) * t535;
t454 = Ifges(7,1) * t485 + Ifges(7,4) * t484 + Ifges(7,5) * t535;
t555 = mrSges(7,1) * t432 - mrSges(7,2) * t433 + Ifges(7,5) * t451 + Ifges(7,6) * t450 + Ifges(7,3) * t528 + t485 * t453 - t484 * t454;
t534 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t584;
t533 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t581;
t530 = (-t549 * mrSges(3,1) + t546 * mrSges(3,2)) * qJD(1);
t519 = Ifges(3,5) * qJD(2) + (t546 * Ifges(3,1) + t549 * Ifges(3,4)) * qJD(1);
t518 = Ifges(3,6) * qJD(2) + (t546 * Ifges(3,4) + Ifges(3,2) * t549) * qJD(1);
t517 = Ifges(3,3) * qJD(2) + (t546 * Ifges(3,5) + t549 * Ifges(3,6)) * qJD(1);
t452 = Ifges(7,5) * t485 + Ifges(7,6) * t484 + Ifges(7,3) * t535;
t425 = -t508 * mrSges(5,3) - t526 * t503 + t556;
t424 = mrSges(7,2) * t437 - mrSges(7,3) * t432 + Ifges(7,1) * t451 + Ifges(7,4) * t450 + Ifges(7,5) * t528 + t452 * t484 - t453 * t535;
t423 = -mrSges(7,1) * t437 + mrSges(7,3) * t433 + Ifges(7,4) * t451 + Ifges(7,2) * t450 + Ifges(7,6) * t528 - t452 * t485 + t454 * t535;
t420 = -mrSges(6,1) * t532 - t501 * t584 + t564;
t419 = mrSges(6,3) * t532 - t562 + t573;
t418 = t592 * t532 + t557 + t573;
t413 = t545 * t424 + t548 * t423 - mrSges(4,2) * t559 + pkin(5) * t563 + mrSges(5,1) * t444 - mrSges(5,3) * t445 - mrSges(4,3) * t446 + mrSges(6,1) * t442 - mrSges(6,2) * t438 + pkin(8) * t587 - qJ(4) * t425 + pkin(4) * t419 + t579 * t532 + t575 * t525 + t600 * t508 + t577 * t507 - t574 * t584;
t412 = t545 * t423 - t548 * t424 - qJ(5) * t558 + mrSges(4,1) * t559 + mrSges(5,2) * t445 + mrSges(4,3) * t447 - mrSges(6,2) * t440 + mrSges(6,3) * t442 - mrSges(5,1) * t443 - pkin(3) * t425 + pkin(4) * t420 + pkin(8) * t421 + t578 * t532 + t575 * t526 - t577 * t508 + t598 * t507 + t576 * t584;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t568 - mrSges(2,2) * t566 + t546 * (mrSges(3,2) * t520 - mrSges(3,3) * t491 + Ifges(3,1) * t531 + Ifges(3,4) * t532 + Ifges(3,5) * qJDD(2) - qJ(3) * t414 - qJD(2) * t518 - t544 * t412 + t591 * t413 + t517 * t584) + t549 * ((-t546 * t517 - qJ(4) * (-t501 - t503) * t549) * qJD(1) + t574 * t526 - qJ(4) * t560 + (mrSges(5,1) * qJ(4) - t578) * t507 + (qJ(4) * mrSges(6,1) + Ifges(3,2) + t599) * t532 + Ifges(3,6) * qJDD(2) + (qJ(4) * t489 + t576) * t525 + t579 * t508 - t555 + Ifges(3,4) * t531 + qJD(2) * t519 - mrSges(3,1) * t520 + mrSges(3,3) * t492 - mrSges(5,2) * t444 - mrSges(4,1) * t446 + mrSges(4,2) * t447 + mrSges(6,3) * t438 - mrSges(6,1) * t440 + mrSges(5,3) * t443 - pkin(2) * t414 + pkin(3) * t418 + qJ(5) * t419 - pkin(5) * t421) + pkin(1) * (-m(3) * t520 + t532 * mrSges(3,1) - t531 * mrSges(3,2) + (-t533 * t546 + t534 * t549) * qJD(1) - t414) + pkin(7) * (t549 * (m(3) * t492 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t532 - qJD(2) * t533 + t530 * t584 + t567) - t546 * (m(3) * t491 + qJDD(2) * mrSges(3,1) - t531 * mrSges(3,3) + qJD(2) * t534 - t530 * t581 - t422)); Ifges(3,5) * t531 + Ifges(3,6) * t532 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t491 - mrSges(3,2) * t492 + t544 * t413 + t591 * t412 - pkin(2) * t422 + qJ(3) * t567 + (t546 * t518 - t549 * t519) * qJD(1); t422; t418; t420; t555;];
tauJ  = t1;
