% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
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
% Datum: 2019-05-06 12:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:42:00
% EndTime: 2019-05-06 12:42:05
% DurationCPUTime: 2.71s
% Computational Cost: add. (7843->299), mult. (15799->328), div. (0->0), fcn. (8533->6), ass. (0->120)
t574 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t615 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t604 = Ifges(5,6) - Ifges(6,6) + Ifges(7,6);
t614 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t606 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t613 = Ifges(3,1) + Ifges(4,2);
t591 = Ifges(3,4) + Ifges(4,6);
t590 = Ifges(3,5) - Ifges(4,4);
t611 = Ifges(3,2) + Ifges(4,3);
t589 = Ifges(3,6) - Ifges(4,5);
t608 = Ifges(3,3) + Ifges(4,1);
t543 = sin(qJ(4));
t546 = cos(qJ(4));
t547 = cos(qJ(2));
t578 = qJD(1) * t547;
t513 = qJD(2) * t543 + t546 * t578;
t514 = qJD(2) * t546 - t543 * t578;
t544 = sin(qJ(2));
t577 = t544 * qJD(1);
t531 = qJD(4) + t577;
t607 = t615 * t513 + t574 * t514 + t604 * t531;
t576 = qJD(1) * qJD(2);
t569 = t544 * t576;
t519 = qJDD(1) * t547 - t569;
t475 = -qJD(4) * t513 + qJDD(2) * t546 - t519 * t543;
t486 = -mrSges(7,1) * t531 - mrSges(7,3) * t514;
t488 = -mrSges(6,1) * t531 + mrSges(6,2) * t514;
t601 = -(t486 + t488) * t514 - (mrSges(7,2) + mrSges(6,3)) * t475;
t605 = -Ifges(6,2) - Ifges(5,3) - Ifges(7,3);
t603 = -t574 * t513 + t614 * t514 + t606 * t531;
t587 = t513 * t531;
t602 = (-t475 + t587) * qJ(5);
t600 = t475 * mrSges(7,2) + t514 * t486;
t527 = pkin(3) * t577 - qJD(2) * pkin(8);
t541 = t547 ^ 2;
t550 = qJD(1) ^ 2;
t570 = t547 * t576;
t518 = qJDD(1) * t544 + t570;
t545 = sin(qJ(1));
t548 = cos(qJ(1));
t568 = g(1) * t545 - t548 * g(2);
t559 = -qJDD(1) * pkin(1) - t568;
t554 = pkin(2) * t569 - 0.2e1 * qJD(3) * t577 + (-t518 - t570) * qJ(3) + t559;
t438 = -t527 * t577 + (-pkin(3) * t541 - pkin(7)) * t550 + (-pkin(2) - pkin(8)) * t519 + t554;
t562 = -g(1) * t548 - t545 * g(2);
t503 = -pkin(1) * t550 + qJDD(1) * pkin(7) + t562;
t489 = -t547 * g(3) - t544 * t503;
t515 = (-t547 * pkin(2) - t544 * qJ(3)) * qJD(1);
t549 = qJD(2) ^ 2;
t448 = -qJDD(2) * pkin(2) - qJ(3) * t549 + t515 * t577 + qJDD(3) - t489;
t442 = (-t544 * t547 * t550 - qJDD(2)) * pkin(8) + (t518 - t570) * pkin(3) + t448;
t434 = -t543 * t438 + t442 * t546;
t478 = pkin(4) * t513 - qJ(5) * t514;
t512 = qJDD(4) + t518;
t528 = t531 ^ 2;
t430 = -pkin(4) * t512 - qJ(5) * t528 + t514 * t478 + qJDD(5) - t434;
t484 = -mrSges(6,2) * t513 + mrSges(6,3) * t531;
t599 = -m(6) * t430 + t512 * mrSges(6,1) + t531 * t484;
t598 = -0.2e1 * t514;
t597 = 2 * qJD(5);
t595 = pkin(7) * t550;
t594 = mrSges(3,1) - mrSges(4,2);
t592 = -mrSges(5,3) - mrSges(6,2);
t435 = t546 * t438 + t543 * t442;
t479 = mrSges(6,1) * t513 - mrSges(6,3) * t514;
t584 = -mrSges(5,1) * t513 - mrSges(5,2) * t514 - t479;
t490 = -t544 * g(3) + t547 * t503;
t582 = t608 * qJD(2) + (t544 * t590 + t547 * t589) * qJD(1);
t581 = t589 * qJD(2) + (t544 * t591 + t547 * t611) * qJD(1);
t580 = t590 * qJD(2) + (t544 * t613 + t547 * t591) * qJD(1);
t524 = -mrSges(4,1) * t578 - qJD(2) * mrSges(4,3);
t579 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t578 - t524;
t575 = qJD(3) * qJD(2);
t429 = -pkin(4) * t528 + t512 * qJ(5) - t513 * t478 + t531 * t597 + t435;
t474 = qJD(4) * t514 + qJDD(2) * t543 + t546 * t519;
t485 = -pkin(5) * t531 - qJ(6) * t514;
t511 = t513 ^ 2;
t425 = -pkin(5) * t511 + qJ(6) * t474 + 0.2e1 * qJD(6) * t513 + t485 * t531 + t429;
t480 = -mrSges(7,1) * t513 + mrSges(7,2) * t514;
t573 = m(7) * t425 + t474 * mrSges(7,3) + t513 * t480;
t572 = t604 * t513 - t606 * t514 + t605 * t531;
t482 = mrSges(7,2) * t531 + mrSges(7,3) * t513;
t483 = -mrSges(5,2) * t531 - mrSges(5,3) * t513;
t422 = qJD(6) * t598 + (-t475 - t587) * qJ(6) + (t513 * t514 - t512) * pkin(5) + t430;
t564 = -m(7) * t422 + t475 * mrSges(7,3) + t514 * t480;
t414 = m(5) * t434 + (t482 + t483) * t531 + t584 * t514 + (mrSges(5,1) + mrSges(7,1)) * t512 + t592 * t475 + t564 + t599;
t487 = mrSges(5,1) * t531 - mrSges(5,3) * t514;
t558 = m(6) * t429 + t512 * mrSges(6,3) + t531 * t488 + t573;
t415 = m(5) * t435 + (t486 - t487) * t531 + t584 * t513 + (-mrSges(5,2) + mrSges(7,2)) * t512 + t592 * t474 + t558;
t566 = -t414 * t543 + t546 * t415;
t535 = -0.2e1 * t575;
t561 = t549 * pkin(2) - qJDD(2) * qJ(3) - t515 * t578 - t490;
t555 = pkin(8) * t541 * t550 - pkin(3) * t519 - qJD(2) * t527 + t561;
t426 = -qJ(6) * t511 + qJDD(6) + t535 + (-pkin(4) - pkin(5)) * t474 - t602 + (-pkin(4) * t531 + t485 + t597) * t514 + t555;
t563 = m(7) * t426 - t474 * mrSges(7,1) - t513 * t482;
t412 = t546 * t414 + t543 * t415;
t443 = -pkin(2) * t519 + t554 - t595;
t560 = m(4) * t443 + t566;
t441 = 0.2e1 * t575 - t555;
t432 = qJD(5) * t598 + t602 + (t514 * t531 + t474) * pkin(4) + t441;
t557 = m(6) * t432 + t474 * mrSges(6,1) + t513 * t484 - t563;
t556 = m(4) * t448 + t518 * mrSges(4,1) + t412;
t420 = -mrSges(7,1) * t512 - t482 * t531 - t564;
t553 = -m(5) * t441 - t474 * mrSges(5,1) - t475 * mrSges(5,2) - t513 * t483 - t514 * t487 - t557;
t447 = t535 + t561;
t516 = (t547 * mrSges(4,2) - t544 * mrSges(4,3)) * qJD(1);
t525 = mrSges(4,1) * t577 + qJD(2) * mrSges(4,2);
t552 = -m(4) * t447 + qJDD(2) * mrSges(4,3) + qJD(2) * t525 + t516 * t578 - t553;
t417 = mrSges(6,2) * t475 + t479 * t514 + t420 - t599;
t551 = -mrSges(6,1) * t430 - mrSges(7,1) * t422 - mrSges(5,2) * t435 - pkin(5) * t420 - pkin(4) * t417 + qJ(5) * (t486 * t531 + t558) + mrSges(7,2) * t425 + mrSges(6,3) * t429 + mrSges(5,1) * t434 + t607 * t514 + (-qJ(5) * t479 + t603) * t513 + (mrSges(7,2) * qJ(5) - t605) * t512 + t606 * t475 + (-mrSges(6,2) * qJ(5) - t604) * t474;
t522 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t577;
t517 = (-t547 * mrSges(3,1) + t544 * mrSges(3,2)) * qJD(1);
t502 = t559 - t595;
t421 = t563 + t600;
t419 = t557 + t601;
t411 = qJDD(2) * mrSges(4,2) + qJD(2) * t524 + t516 * t577 + t556;
t410 = mrSges(4,2) * t519 - mrSges(4,3) * t518 + (t524 * t547 - t525 * t544) * qJD(1) + t560;
t409 = mrSges(5,2) * t441 + mrSges(6,2) * t430 + mrSges(7,2) * t426 - mrSges(5,3) * t434 - mrSges(6,3) * t432 - mrSges(7,3) * t422 - qJ(5) * t419 - qJ(6) * t420 - t574 * t474 + t614 * t475 + t606 * t512 + t572 * t513 - t607 * t531;
t408 = -mrSges(5,1) * t441 + mrSges(5,3) * t435 - mrSges(6,1) * t432 + mrSges(6,2) * t429 + mrSges(7,1) * t426 - mrSges(7,3) * t425 + pkin(5) * t421 - qJ(6) * t573 - pkin(4) * t419 + (-qJ(6) * t486 + t603) * t531 + t572 * t514 + (-qJ(6) * mrSges(7,2) + t604) * t512 + t574 * t475 + t615 * t474;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t568 - mrSges(2,2) * t562 + t544 * (mrSges(4,1) * t448 + mrSges(3,2) * t502 - mrSges(3,3) * t489 - mrSges(4,3) * t443 + pkin(3) * t412 - qJ(3) * t410 - t581 * qJD(2) + t590 * qJDD(2) + t613 * t518 + t591 * t519 + t582 * t578 + t551) + t547 * (-mrSges(3,1) * t502 + mrSges(3,3) * t490 - mrSges(4,1) * t447 + mrSges(4,2) * t443 - t543 * t409 - t546 * t408 - pkin(3) * (t553 - t601) - pkin(8) * t566 - pkin(2) * t410 + t611 * t519 + t591 * t518 + t589 * qJDD(2) + t580 * qJD(2) - t582 * t577) + pkin(1) * (-m(3) * t502 + t594 * t519 + (-mrSges(3,2) + mrSges(4,3)) * t518 + (t579 * t547 + (-t522 + t525) * t544) * qJD(1) - t560) + pkin(7) * (t547 * ((mrSges(3,3) + mrSges(4,1)) * t519 + t517 * t578 + t552 - qJD(2) * t522 + m(3) * t490 - qJDD(2) * mrSges(3,2) + t601) + (-m(3) * t489 + t518 * mrSges(3,3) - t594 * qJDD(2) - t579 * qJD(2) + (t516 + t517) * t577 + t556) * t544); mrSges(3,1) * t489 - mrSges(3,2) * t490 + mrSges(4,2) * t448 - mrSges(4,3) * t447 + t546 * t409 - t543 * t408 - pkin(8) * t412 - pkin(2) * t411 + qJ(3) * (-t475 * mrSges(6,3) - t514 * t488 + t552 - t600) + (qJ(3) * mrSges(4,1) + t589) * t519 + t590 * t518 + t608 * qJDD(2) + (t581 * t544 - t580 * t547) * qJD(1); t411; t551; t417; t421;];
tauJ  = t1;
