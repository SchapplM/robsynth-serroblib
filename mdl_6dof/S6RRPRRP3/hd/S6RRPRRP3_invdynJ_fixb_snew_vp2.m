% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:30:13
% EndTime: 2019-05-06 17:30:24
% DurationCPUTime: 5.92s
% Computational Cost: add. (55385->327), mult. (126347->403), div. (0->0), fcn. (90868->10), ass. (0->127)
t563 = Ifges(6,4) + Ifges(7,4);
t572 = Ifges(6,2) + Ifges(7,2);
t568 = Ifges(6,6) + Ifges(7,6);
t534 = sin(qJ(2));
t538 = cos(qJ(2));
t556 = qJD(1) * qJD(2);
t522 = qJDD(1) * t534 + t538 * t556;
t523 = qJDD(1) * t538 - t534 * t556;
t530 = sin(pkin(10));
t531 = cos(pkin(10));
t499 = t522 * t531 + t523 * t530;
t513 = (t530 * t538 + t531 * t534) * qJD(1);
t533 = sin(qJ(4));
t537 = cos(qJ(4));
t502 = qJD(2) * t533 + t513 * t537;
t472 = -qJD(4) * t502 + qJDD(2) * t537 - t499 * t533;
t501 = qJD(2) * t537 - t513 * t533;
t473 = qJD(4) * t501 + qJDD(2) * t533 + t499 * t537;
t532 = sin(qJ(5));
t536 = cos(qJ(5));
t478 = t501 * t536 - t502 * t532;
t437 = qJD(5) * t478 + t472 * t532 + t473 * t536;
t479 = t501 * t532 + t502 * t536;
t458 = -mrSges(7,1) * t478 + mrSges(7,2) * t479;
t541 = qJD(1) ^ 2;
t535 = sin(qJ(1));
t539 = cos(qJ(1));
t547 = -g(1) * t539 - g(2) * t535;
t519 = -pkin(1) * t541 + qJDD(1) * pkin(7) + t547;
t562 = t534 * t519;
t564 = pkin(2) * t541;
t480 = qJDD(2) * pkin(2) - t522 * qJ(3) - t562 + (qJ(3) * t556 + t534 * t564 - g(3)) * t538;
t504 = -g(3) * t534 + t538 * t519;
t559 = qJD(1) * t534;
t524 = qJD(2) * pkin(2) - qJ(3) * t559;
t529 = t538 ^ 2;
t481 = qJ(3) * t523 - qJD(2) * t524 - t529 * t564 + t504;
t558 = qJD(1) * t538;
t512 = -t530 * t559 + t531 * t558;
t457 = 0.2e1 * qJD(3) * t512 + t530 * t480 + t531 * t481;
t495 = -pkin(3) * t512 - pkin(8) * t513;
t540 = qJD(2) ^ 2;
t442 = -pkin(3) * t540 + qJDD(2) * pkin(8) + t495 * t512 + t457;
t553 = t535 * g(1) - t539 * g(2);
t546 = -qJDD(1) * pkin(1) - t553;
t484 = -t523 * pkin(2) + qJDD(3) + t524 * t559 + (-qJ(3) * t529 - pkin(7)) * t541 + t546;
t498 = -t530 * t522 + t523 * t531;
t446 = (-qJD(2) * t512 - t499) * pkin(8) + (qJD(2) * t513 - t498) * pkin(3) + t484;
t426 = -t533 * t442 + t537 * t446;
t497 = qJDD(4) - t498;
t511 = qJD(4) - t512;
t422 = (t501 * t511 - t473) * pkin(9) + (t501 * t502 + t497) * pkin(4) + t426;
t427 = t537 * t442 + t533 * t446;
t487 = pkin(4) * t511 - pkin(9) * t502;
t500 = t501 ^ 2;
t424 = -pkin(4) * t500 + pkin(9) * t472 - t487 * t511 + t427;
t416 = t536 * t422 - t532 * t424;
t496 = qJDD(5) + t497;
t507 = qJD(5) + t511;
t411 = -0.2e1 * qJD(6) * t479 + (t478 * t507 - t437) * qJ(6) + (t478 * t479 + t496) * pkin(5) + t416;
t461 = -mrSges(7,2) * t507 + mrSges(7,3) * t478;
t555 = m(7) * t411 + t496 * mrSges(7,1) + t507 * t461;
t408 = -t437 * mrSges(7,3) - t479 * t458 + t555;
t417 = t532 * t422 + t536 * t424;
t436 = -qJD(5) * t479 + t472 * t536 - t473 * t532;
t463 = pkin(5) * t507 - qJ(6) * t479;
t477 = t478 ^ 2;
t413 = -pkin(5) * t477 + qJ(6) * t436 + 0.2e1 * qJD(6) * t478 - t463 * t507 + t417;
t569 = Ifges(6,5) + Ifges(7,5);
t570 = Ifges(6,1) + Ifges(7,1);
t560 = -t563 * t478 - t570 * t479 - t569 * t507;
t566 = t572 * t478 + t563 * t479 + t568 * t507;
t567 = Ifges(6,3) + Ifges(7,3);
t571 = mrSges(6,1) * t416 + mrSges(7,1) * t411 - mrSges(6,2) * t417 - mrSges(7,2) * t413 + pkin(5) * t408 + t568 * t436 + t569 * t437 + t560 * t478 + t566 * t479 + t567 * t496;
t459 = -mrSges(6,1) * t478 + mrSges(6,2) * t479;
t462 = -mrSges(6,2) * t507 + mrSges(6,3) * t478;
t401 = m(6) * t416 + t496 * mrSges(6,1) + t507 * t462 + (-t458 - t459) * t479 + (-mrSges(6,3) - mrSges(7,3)) * t437 + t555;
t464 = mrSges(7,1) * t507 - mrSges(7,3) * t479;
t465 = mrSges(6,1) * t507 - mrSges(6,3) * t479;
t554 = m(7) * t413 + t436 * mrSges(7,3) + t478 * t458;
t404 = m(6) * t417 + t436 * mrSges(6,3) + t478 * t459 + (-t464 - t465) * t507 + (-mrSges(6,2) - mrSges(7,2)) * t496 + t554;
t399 = t536 * t401 + t532 * t404;
t467 = Ifges(5,4) * t502 + Ifges(5,2) * t501 + Ifges(5,6) * t511;
t468 = Ifges(5,1) * t502 + Ifges(5,4) * t501 + Ifges(5,5) * t511;
t565 = mrSges(5,1) * t426 - mrSges(5,2) * t427 + Ifges(5,5) * t473 + Ifges(5,6) * t472 + Ifges(5,3) * t497 + pkin(4) * t399 + t502 * t467 - t501 * t468 + t571;
t456 = -0.2e1 * qJD(3) * t513 + t480 * t531 - t530 * t481;
t491 = -mrSges(4,1) * t512 + mrSges(4,2) * t513;
t506 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t513;
t482 = -mrSges(5,1) * t501 + mrSges(5,2) * t502;
t485 = -mrSges(5,2) * t511 + mrSges(5,3) * t501;
t396 = m(5) * t426 + mrSges(5,1) * t497 - mrSges(5,3) * t473 - t482 * t502 + t485 * t511 + t399;
t486 = mrSges(5,1) * t511 - mrSges(5,3) * t502;
t549 = -t401 * t532 + t536 * t404;
t397 = m(5) * t427 - mrSges(5,2) * t497 + mrSges(5,3) * t472 + t482 * t501 - t486 * t511 + t549;
t550 = -t396 * t533 + t537 * t397;
t389 = m(4) * t457 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t498 - qJD(2) * t506 + t491 * t512 + t550;
t505 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t512;
t441 = -qJDD(2) * pkin(3) - pkin(8) * t540 + t513 * t495 - t456;
t425 = -pkin(4) * t472 - pkin(9) * t500 + t502 * t487 + t441;
t419 = -pkin(5) * t436 - qJ(6) * t477 + t463 * t479 + qJDD(6) + t425;
t414 = m(7) * t419 - t436 * mrSges(7,1) + t437 * mrSges(7,2) - t478 * t461 + t479 * t464;
t545 = m(6) * t425 - t436 * mrSges(6,1) + mrSges(6,2) * t437 - t478 * t462 + t465 * t479 + t414;
t543 = -m(5) * t441 + t472 * mrSges(5,1) - mrSges(5,2) * t473 + t501 * t485 - t486 * t502 - t545;
t406 = m(4) * t456 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t499 + qJD(2) * t505 - t491 * t513 + t543;
t386 = t530 * t389 + t531 * t406;
t391 = t537 * t396 + t533 * t397;
t561 = -t568 * t478 - t569 * t479 - t567 * t507;
t551 = t531 * t389 - t406 * t530;
t390 = m(4) * t484 - t498 * mrSges(4,1) + t499 * mrSges(4,2) - t512 * t505 + t513 * t506 + t391;
t526 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t558;
t525 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t559;
t521 = (-mrSges(3,1) * t538 + mrSges(3,2) * t534) * qJD(1);
t518 = -t541 * pkin(7) + t546;
t516 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t534 + Ifges(3,4) * t538) * qJD(1);
t515 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t534 + Ifges(3,2) * t538) * qJD(1);
t503 = -t538 * g(3) - t562;
t490 = Ifges(4,1) * t513 + Ifges(4,4) * t512 + Ifges(4,5) * qJD(2);
t489 = Ifges(4,4) * t513 + Ifges(4,2) * t512 + Ifges(4,6) * qJD(2);
t488 = Ifges(4,5) * t513 + Ifges(4,6) * t512 + Ifges(4,3) * qJD(2);
t466 = Ifges(5,5) * t502 + Ifges(5,6) * t501 + Ifges(5,3) * t511;
t398 = mrSges(6,2) * t425 + mrSges(7,2) * t419 - mrSges(6,3) * t416 - mrSges(7,3) * t411 - qJ(6) * t408 + t563 * t436 + t570 * t437 - t561 * t478 + t569 * t496 - t566 * t507;
t392 = -mrSges(6,1) * t425 + mrSges(6,3) * t417 - mrSges(7,1) * t419 + mrSges(7,3) * t413 - pkin(5) * t414 + qJ(6) * t554 + (-qJ(6) * t464 - t560) * t507 + (-mrSges(7,2) * qJ(6) + t568) * t496 + t561 * t479 + t563 * t437 + t572 * t436;
t385 = mrSges(5,2) * t441 - mrSges(5,3) * t426 + Ifges(5,1) * t473 + Ifges(5,4) * t472 + Ifges(5,5) * t497 - pkin(9) * t399 - t392 * t532 + t398 * t536 + t466 * t501 - t467 * t511;
t384 = -mrSges(5,1) * t441 + mrSges(5,3) * t427 + Ifges(5,4) * t473 + Ifges(5,2) * t472 + Ifges(5,6) * t497 - pkin(4) * t545 + pkin(9) * t549 + t536 * t392 + t532 * t398 - t502 * t466 + t511 * t468;
t383 = -mrSges(4,1) * t484 + mrSges(4,3) * t457 + Ifges(4,4) * t499 + Ifges(4,2) * t498 + Ifges(4,6) * qJDD(2) - pkin(3) * t391 + qJD(2) * t490 - t513 * t488 - t565;
t382 = mrSges(4,2) * t484 - mrSges(4,3) * t456 + Ifges(4,1) * t499 + Ifges(4,4) * t498 + Ifges(4,5) * qJDD(2) - pkin(8) * t391 - qJD(2) * t489 - t384 * t533 + t385 * t537 + t488 * t512;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t553 - mrSges(2,2) * t547 + t534 * (mrSges(3,2) * t518 - mrSges(3,3) * t503 + Ifges(3,1) * t522 + Ifges(3,4) * t523 + Ifges(3,5) * qJDD(2) - qJ(3) * t386 - qJD(2) * t515 + t531 * t382 - t530 * t383) + t538 * (-mrSges(3,1) * t518 + mrSges(3,3) * t504 + Ifges(3,4) * t522 + Ifges(3,2) * t523 + Ifges(3,6) * qJDD(2) - pkin(2) * t390 + qJ(3) * t551 + qJD(2) * t516 + t530 * t382 + t531 * t383) + pkin(1) * (-m(3) * t518 + t523 * mrSges(3,1) - t522 * mrSges(3,2) + (-t525 * t534 + t526 * t538) * qJD(1) - t390) + pkin(7) * (t538 * (m(3) * t504 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t523 - qJD(2) * t525 + t521 * t558 + t551) - t534 * (m(3) * t503 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t522 + qJD(2) * t526 - t521 * t559 + t386)); Ifges(3,5) * t522 + Ifges(3,6) * t523 + mrSges(3,1) * t503 - mrSges(3,2) * t504 + Ifges(4,5) * t499 + Ifges(4,6) * t498 + t513 * t489 - t512 * t490 + mrSges(4,1) * t456 - mrSges(4,2) * t457 + t533 * t385 + t537 * t384 + pkin(3) * t543 + pkin(8) * t550 + pkin(2) * t386 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t534 * t515 - t538 * t516) * qJD(1); t390; t565; t571; t414;];
tauJ  = t1;
