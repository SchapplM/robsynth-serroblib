% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 09:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRP12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP12_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP12_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:24:16
% EndTime: 2019-05-07 09:24:29
% DurationCPUTime: 4.85s
% Computational Cost: add. (40263->319), mult. (86114->386), div. (0->0), fcn. (65462->10), ass. (0->137)
t606 = Ifges(6,1) + Ifges(7,1);
t590 = Ifges(6,4) - Ifges(7,5);
t603 = Ifges(7,4) + Ifges(6,5);
t605 = Ifges(6,2) + Ifges(7,3);
t601 = Ifges(6,6) - Ifges(7,6);
t604 = Ifges(4,1) + Ifges(5,2);
t589 = -Ifges(4,5) + Ifges(5,4);
t602 = -Ifges(4,2) - Ifges(5,3);
t588 = Ifges(4,6) - Ifges(5,5);
t587 = -Ifges(5,6) - Ifges(4,4);
t600 = Ifges(4,3) + Ifges(5,1);
t599 = Ifges(6,3) + Ifges(7,2);
t546 = cos(pkin(6));
t542 = qJD(1) * t546 + qJD(2);
t548 = sin(qJ(3));
t549 = sin(qJ(2));
t545 = sin(pkin(6));
t574 = qJD(1) * t545;
t569 = t549 * t574;
t595 = cos(qJ(3));
t519 = -t542 * t595 + t548 * t569;
t551 = cos(qJ(2));
t573 = qJD(1) * t551;
t568 = t545 * t573;
t538 = -qJD(3) + t568;
t547 = sin(qJ(5));
t594 = cos(qJ(5));
t500 = -t519 * t594 - t538 * t547;
t501 = t547 * t519 - t538 * t594;
t520 = t548 * t542 + t569 * t595;
t517 = qJD(5) + t520;
t598 = t605 * t500 - t590 * t501 - t601 * t517;
t597 = -t590 * t500 + t606 * t501 + t603 * t517;
t571 = qJD(1) * qJD(2);
t533 = (-qJDD(1) * t551 + t549 * t571) * t545;
t531 = (-pkin(2) * t551 - pkin(9) * t549) * t574;
t540 = t542 ^ 2;
t541 = qJDD(1) * t546 + qJDD(2);
t553 = qJD(1) ^ 2;
t550 = sin(qJ(1));
t552 = cos(qJ(1));
t564 = -g(1) * t552 - g(2) * t550;
t593 = pkin(8) * t545;
t529 = -pkin(1) * t553 + qJDD(1) * t593 + t564;
t567 = t550 * g(1) - g(2) * t552;
t528 = qJDD(1) * pkin(1) + t553 * t593 + t567;
t584 = t528 * t546;
t575 = t551 * t529 + t549 * t584;
t460 = -pkin(2) * t540 + pkin(9) * t541 + (-g(3) * t549 + t531 * t573) * t545 + t575;
t532 = (qJDD(1) * t549 + t551 * t571) * t545;
t592 = g(3) * t546;
t461 = pkin(2) * t533 - pkin(9) * t532 - t592 + (-t528 + (pkin(2) * t549 - pkin(9) * t551) * t542 * qJD(1)) * t545;
t437 = t595 * t460 + t548 * t461;
t495 = pkin(3) * t519 - qJ(4) * t520;
t525 = qJDD(3) + t533;
t537 = t538 ^ 2;
t433 = pkin(3) * t537 - t525 * qJ(4) + 0.2e1 * qJD(4) * t538 + t519 * t495 - t437;
t504 = mrSges(5,1) * t519 + mrSges(5,3) * t538;
t436 = -t548 * t460 + t461 * t595;
t434 = -t525 * pkin(3) - t537 * qJ(4) + t520 * t495 + qJDD(4) - t436;
t492 = -t519 * qJD(3) + t532 * t595 + t548 * t541;
t585 = t519 * t538;
t428 = (t519 * t520 - t525) * pkin(10) + (t492 - t585) * pkin(4) + t434;
t491 = qJD(3) * t520 + t532 * t548 - t541 * t595;
t506 = pkin(4) * t520 + pkin(10) * t538;
t518 = t519 ^ 2;
t582 = t545 * t551;
t493 = -g(3) * t582 - t549 * t529 + t551 * t584;
t459 = -pkin(2) * t541 - pkin(9) * t540 + t531 * t569 - t493;
t557 = (-t492 - t585) * qJ(4) + t459 + (-t538 * pkin(3) - 0.2e1 * qJD(4)) * t520;
t432 = -pkin(4) * t518 - t506 * t520 + (pkin(3) + pkin(10)) * t491 + t557;
t424 = t547 * t428 + t594 * t432;
t446 = qJD(5) * t501 - t491 * t594 + t525 * t547;
t472 = mrSges(6,1) * t517 - mrSges(6,3) * t501;
t488 = qJDD(5) + t492;
t464 = pkin(5) * t500 - qJ(6) * t501;
t516 = t517 ^ 2;
t420 = -pkin(5) * t516 + qJ(6) * t488 + 0.2e1 * qJD(6) * t517 - t464 * t500 + t424;
t473 = -mrSges(7,1) * t517 + mrSges(7,2) * t501;
t570 = m(7) * t420 + t488 * mrSges(7,3) + t517 * t473;
t465 = mrSges(7,1) * t500 - mrSges(7,3) * t501;
t579 = -mrSges(6,1) * t500 - mrSges(6,2) * t501 - t465;
t591 = -mrSges(6,3) - mrSges(7,2);
t410 = m(6) * t424 - mrSges(6,2) * t488 + t446 * t591 - t472 * t517 + t500 * t579 + t570;
t423 = t428 * t594 - t547 * t432;
t447 = -t500 * qJD(5) + t547 * t491 + t525 * t594;
t471 = -mrSges(6,2) * t517 - mrSges(6,3) * t500;
t421 = -t488 * pkin(5) - t516 * qJ(6) + t501 * t464 + qJDD(6) - t423;
t470 = -mrSges(7,2) * t500 + mrSges(7,3) * t517;
t565 = -m(7) * t421 + t488 * mrSges(7,1) + t517 * t470;
t412 = m(6) * t423 + mrSges(6,1) * t488 + t447 * t591 + t471 * t517 + t501 * t579 + t565;
t405 = t547 * t410 + t412 * t594;
t497 = -mrSges(5,2) * t519 - mrSges(5,3) * t520;
t560 = m(5) * t434 + t492 * mrSges(5,1) + t520 * t497 + t405;
t402 = t525 * mrSges(5,2) - t538 * t504 + t560;
t431 = -pkin(4) * t491 - pkin(10) * t518 - t538 * t506 - t433;
t426 = -0.2e1 * qJD(6) * t501 + (t500 * t517 - t447) * qJ(6) + (t501 * t517 + t446) * pkin(5) + t431;
t417 = m(7) * t426 + t446 * mrSges(7,1) - mrSges(7,3) * t447 + t500 * t470 - t473 * t501;
t580 = t601 * t500 - t603 * t501 - t599 * t517;
t403 = -mrSges(6,1) * t431 - mrSges(7,1) * t426 + mrSges(7,2) * t420 + mrSges(6,3) * t424 - pkin(5) * t417 - t605 * t446 + t590 * t447 + t601 * t488 + t580 * t501 + t597 * t517;
t404 = mrSges(6,2) * t431 + mrSges(7,2) * t421 - mrSges(6,3) * t423 - mrSges(7,3) * t426 - qJ(6) * t417 - t590 * t446 + t606 * t447 + t603 * t488 + t580 * t500 + t598 * t517;
t505 = mrSges(5,1) * t520 - mrSges(5,2) * t538;
t559 = m(6) * t431 + mrSges(6,1) * t446 + t447 * mrSges(6,2) + t471 * t500 + t501 * t472 + t417;
t556 = -m(5) * t433 + t525 * mrSges(5,3) - t538 * t505 + t559;
t576 = t587 * t519 + t604 * t520 + t589 * t538;
t577 = t602 * t519 - t587 * t520 - t588 * t538;
t596 = -t491 * t588 - t492 * t589 + t519 * t576 + t520 * t577 + t600 * t525 + mrSges(4,1) * t436 - mrSges(4,2) * t437 + mrSges(5,2) * t434 - mrSges(5,3) * t433 - pkin(3) * t402 - pkin(10) * t405 + qJ(4) * (-mrSges(5,1) * t491 - t497 * t519 + t556) - t547 * t403 + t404 * t594;
t583 = t545 * t549;
t496 = mrSges(4,1) * t519 + mrSges(4,2) * t520;
t502 = mrSges(4,2) * t538 - mrSges(4,3) * t519;
t400 = m(4) * t436 - t492 * mrSges(4,3) - t520 * t496 + (-t502 + t504) * t538 + (mrSges(4,1) - mrSges(5,2)) * t525 - t560;
t503 = -mrSges(4,1) * t538 - mrSges(4,3) * t520;
t408 = t556 + (-t496 - t497) * t519 + (-mrSges(4,3) - mrSges(5,1)) * t491 + m(4) * t437 - mrSges(4,2) * t525 + t503 * t538;
t397 = t595 * t400 + t548 * t408;
t581 = t594 * t410 - t547 * t412;
t578 = t588 * t519 + t589 * t520 + t600 * t538;
t566 = -t400 * t548 + t595 * t408;
t435 = pkin(3) * t491 + t557;
t563 = -m(5) * t435 + t491 * mrSges(5,2) + t519 * t504 - t581;
t558 = -m(4) * t459 - t491 * mrSges(4,1) - t519 * t502 + (-t503 + t505) * t520 + (-mrSges(4,2) + mrSges(5,3)) * t492 + t563;
t416 = mrSges(7,2) * t447 + t465 * t501 - t565;
t555 = mrSges(6,1) * t423 - mrSges(7,1) * t421 - mrSges(6,2) * t424 + mrSges(7,3) * t420 - pkin(5) * t416 + qJ(6) * t570 - t598 * t501 + (-qJ(6) * t465 + t597) * t500 + t599 * t488 + t603 * t447 + (-qJ(6) * mrSges(7,2) - t601) * t446;
t530 = (-mrSges(3,1) * t551 + mrSges(3,2) * t549) * t574;
t527 = -mrSges(3,2) * t542 + mrSges(3,3) * t568;
t526 = mrSges(3,1) * t542 - mrSges(3,3) * t569;
t511 = -t528 * t545 - t592;
t510 = Ifges(3,5) * t542 + (Ifges(3,1) * t549 + Ifges(3,4) * t551) * t574;
t509 = Ifges(3,6) * t542 + (Ifges(3,4) * t549 + Ifges(3,2) * t551) * t574;
t508 = Ifges(3,3) * t542 + (Ifges(3,5) * t549 + Ifges(3,6) * t551) * t574;
t494 = -g(3) * t583 + t575;
t401 = -mrSges(5,3) * t492 - t505 * t520 - t563;
t398 = m(3) * t493 + mrSges(3,1) * t541 - mrSges(3,3) * t532 + t527 * t542 - t530 * t569 + t558;
t396 = m(3) * t494 - mrSges(3,2) * t541 - mrSges(3,3) * t533 - t526 * t542 + t530 * t568 + t566;
t395 = mrSges(5,1) * t434 + mrSges(4,2) * t459 - mrSges(4,3) * t436 - mrSges(5,3) * t435 + pkin(4) * t405 - qJ(4) * t401 + t587 * t491 + t604 * t492 + t578 * t519 - t589 * t525 + t577 * t538 + t555;
t394 = -mrSges(4,1) * t459 - mrSges(5,1) * t433 + mrSges(5,2) * t435 + mrSges(4,3) * t437 - pkin(3) * t401 + pkin(4) * t559 - pkin(10) * t581 - t594 * t403 - t547 * t404 + t602 * t491 - t587 * t492 + t578 * t520 + t588 * t525 - t576 * t538;
t393 = Ifges(3,5) * t532 - Ifges(3,6) * t533 + Ifges(3,3) * t541 + mrSges(3,1) * t493 - mrSges(3,2) * t494 + t548 * t395 + t595 * t394 + pkin(2) * t558 + pkin(9) * t566 + (t509 * t549 - t510 * t551) * t574;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t567 - mrSges(2,2) * t564 + (mrSges(3,2) * t511 - mrSges(3,3) * t493 + Ifges(3,1) * t532 - Ifges(3,4) * t533 + Ifges(3,5) * t541 - pkin(9) * t397 - t548 * t394 + t395 * t595 + t508 * t568 - t542 * t509) * t583 + (-mrSges(3,1) * t511 + mrSges(3,3) * t494 + Ifges(3,4) * t532 - Ifges(3,2) * t533 + Ifges(3,6) * t541 - pkin(2) * t397 - t508 * t569 + t542 * t510 - t596) * t582 + t546 * t393 + pkin(1) * ((t396 * t549 + t398 * t551) * t546 + (-m(3) * t511 - t533 * mrSges(3,1) - t532 * mrSges(3,2) + (-t526 * t549 + t527 * t551) * t574 - t397) * t545) + (t396 * t551 - t398 * t549) * t593; t393; t596; t402; t555; t416;];
tauJ  = t1;
