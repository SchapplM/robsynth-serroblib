% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-05-06 12:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:35:51
% EndTime: 2019-05-06 12:35:57
% DurationCPUTime: 3.41s
% Computational Cost: add. (21022->311), mult. (44114->366), div. (0->0), fcn. (26576->8), ass. (0->125)
t581 = Ifges(6,1) + Ifges(7,1);
t565 = Ifges(6,4) - Ifges(7,5);
t577 = Ifges(7,4) + Ifges(6,5);
t580 = Ifges(6,2) + Ifges(7,3);
t574 = Ifges(6,6) - Ifges(7,6);
t579 = -2 * qJD(3);
t578 = Ifges(3,1) + Ifges(4,2);
t566 = Ifges(3,4) + Ifges(4,6);
t564 = Ifges(3,5) - Ifges(4,4);
t576 = Ifges(3,2) + Ifges(4,3);
t575 = -Ifges(7,2) - Ifges(6,3);
t563 = Ifges(3,6) - Ifges(4,5);
t573 = Ifges(3,3) + Ifges(4,1);
t527 = sin(qJ(4));
t530 = cos(qJ(4));
t531 = cos(qJ(2));
t554 = qJD(1) * t531;
t501 = -qJD(2) * t527 - t530 * t554;
t502 = qJD(2) * t530 - t527 * t554;
t526 = sin(pkin(9));
t562 = cos(pkin(9));
t471 = -t562 * t501 + t502 * t526;
t472 = t526 * t501 + t562 * t502;
t528 = sin(qJ(2));
t555 = qJD(1) * t528;
t517 = qJD(4) + t555;
t572 = t580 * t471 - t565 * t472 - t574 * t517;
t571 = -t565 * t471 + t581 * t472 + t577 * t517;
t534 = qJD(1) ^ 2;
t529 = sin(qJ(1));
t532 = cos(qJ(1));
t544 = -g(1) * t532 - g(2) * t529;
t492 = -pkin(1) * t534 + qJDD(1) * pkin(7) + t544;
t478 = -g(3) * t528 + t531 * t492;
t503 = (-pkin(2) * t531 - qJ(3) * t528) * qJD(1);
t533 = qJD(2) ^ 2;
t453 = pkin(2) * t533 - qJDD(2) * qJ(3) + qJD(2) * t579 - t503 * t554 - t478;
t570 = -2 * qJD(5);
t569 = pkin(7) * t534;
t568 = mrSges(3,1) - mrSges(4,2);
t567 = -mrSges(6,3) - mrSges(7,2);
t553 = qJD(1) * qJD(2);
t549 = t528 * t553;
t507 = qJDD(1) * t531 - t549;
t514 = pkin(3) * t555 - qJD(2) * pkin(8);
t525 = t531 ^ 2;
t550 = t531 * t553;
t506 = qJDD(1) * t528 + t550;
t548 = g(1) * t529 - t532 * g(2);
t542 = -qJDD(1) * pkin(1) - t548;
t538 = pkin(2) * t549 + t555 * t579 + (-t506 - t550) * qJ(3) + t542;
t424 = -t514 * t555 + (-pkin(3) * t525 - pkin(7)) * t534 + (-pkin(2) - pkin(8)) * t507 + t538;
t477 = -t531 * g(3) - t528 * t492;
t454 = -qJDD(2) * pkin(2) - qJ(3) * t533 + t503 * t555 + qJDD(3) - t477;
t435 = (-t528 * t531 * t534 - qJDD(2)) * pkin(8) + (t506 - t550) * pkin(3) + t454;
t417 = -t424 * t527 + t530 * t435;
t470 = qJD(4) * t501 + qJDD(2) * t530 - t507 * t527;
t500 = qJDD(4) + t506;
t413 = (t501 * t517 - t470) * qJ(5) + (t501 * t502 + t500) * pkin(4) + t417;
t418 = t530 * t424 + t527 * t435;
t469 = -qJD(4) * t502 - qJDD(2) * t527 - t507 * t530;
t475 = pkin(4) * t517 - qJ(5) * t502;
t499 = t501 ^ 2;
t415 = -pkin(4) * t499 + qJ(5) * t469 - t475 * t517 + t418;
t409 = t526 * t413 + t562 * t415 + t471 * t570;
t442 = -t562 * t469 + t470 * t526;
t457 = mrSges(6,1) * t517 - mrSges(6,3) * t472;
t446 = pkin(5) * t471 - qJ(6) * t472;
t515 = t517 ^ 2;
t405 = -pkin(5) * t515 + qJ(6) * t500 + 0.2e1 * qJD(6) * t517 - t446 * t471 + t409;
t458 = -mrSges(7,1) * t517 + mrSges(7,2) * t472;
t551 = m(7) * t405 + t500 * mrSges(7,3) + t517 * t458;
t447 = mrSges(7,1) * t471 - mrSges(7,3) * t472;
t560 = -mrSges(6,1) * t471 - mrSges(6,2) * t472 - t447;
t395 = m(6) * t409 - mrSges(6,2) * t500 + t567 * t442 - t457 * t517 + t560 * t471 + t551;
t541 = t562 * t413 - t526 * t415;
t408 = t472 * t570 + t541;
t443 = t526 * t469 + t562 * t470;
t455 = -mrSges(6,2) * t517 - mrSges(6,3) * t471;
t406 = -t500 * pkin(5) - t515 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t446) * t472 - t541;
t456 = -mrSges(7,2) * t471 + mrSges(7,3) * t517;
t545 = -m(7) * t406 + t500 * mrSges(7,1) + t517 * t456;
t397 = m(6) * t408 + mrSges(6,1) * t500 + t567 * t443 + t455 * t517 + t560 * t472 + t545;
t390 = t526 * t395 + t562 * t397;
t561 = t471 * t574 - t472 * t577 + t517 * t575;
t559 = t573 * qJD(2) + (t528 * t564 + t531 * t563) * qJD(1);
t558 = t563 * qJD(2) + (t528 * t566 + t531 * t576) * qJD(1);
t557 = t564 * qJD(2) + (t528 * t578 + t531 * t566) * qJD(1);
t512 = -mrSges(4,1) * t554 - qJD(2) * mrSges(4,3);
t556 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t554 - t512;
t473 = -mrSges(5,1) * t501 + mrSges(5,2) * t502;
t474 = -mrSges(5,2) * t517 + mrSges(5,3) * t501;
t387 = m(5) * t417 + mrSges(5,1) * t500 - mrSges(5,3) * t470 - t473 * t502 + t474 * t517 + t390;
t476 = mrSges(5,1) * t517 - mrSges(5,3) * t502;
t546 = t562 * t395 - t397 * t526;
t388 = m(5) * t418 - mrSges(5,2) * t500 + mrSges(5,3) * t469 + t473 * t501 - t476 * t517 + t546;
t547 = -t387 * t527 + t530 * t388;
t385 = t530 * t387 + t527 * t388;
t449 = -pkin(2) * t507 + t538 - t569;
t543 = m(4) * t449 + t547;
t434 = -pkin(8) * t525 * t534 + pkin(3) * t507 + qJD(2) * t514 - t453;
t420 = -pkin(4) * t469 - qJ(5) * t499 + t502 * t475 + qJDD(5) + t434;
t411 = (t472 * t517 + t442) * pkin(5) + (t471 * t517 - t443) * qJ(6) - 0.2e1 * qJD(6) * t472 + t420;
t402 = m(7) * t411 + t442 * mrSges(7,1) - t443 * mrSges(7,3) + t471 * t456 - t472 * t458;
t539 = m(4) * t454 + t506 * mrSges(4,1) + t385;
t398 = m(6) * t420 + t442 * mrSges(6,1) + t443 * mrSges(6,2) + t471 * t455 + t472 * t457 + t402;
t537 = -m(5) * t434 + t469 * mrSges(5,1) - t470 * mrSges(5,2) + t501 * t474 - t502 * t476 - t398;
t504 = (mrSges(4,2) * t531 - mrSges(4,3) * t528) * qJD(1);
t513 = mrSges(4,1) * t555 + qJD(2) * mrSges(4,2);
t536 = -m(4) * t453 + qJDD(2) * mrSges(4,3) + qJD(2) * t513 + t504 * t554 - t537;
t401 = mrSges(7,2) * t443 + t447 * t472 - t545;
t461 = Ifges(5,4) * t502 + Ifges(5,2) * t501 + Ifges(5,6) * t517;
t462 = Ifges(5,1) * t502 + Ifges(5,4) * t501 + Ifges(5,5) * t517;
t535 = mrSges(5,1) * t417 + mrSges(6,1) * t408 - mrSges(7,1) * t406 - mrSges(5,2) * t418 - mrSges(6,2) * t409 + mrSges(7,3) * t405 + Ifges(5,5) * t470 + Ifges(5,6) * t469 + pkin(4) * t390 - pkin(5) * t401 + qJ(6) * t551 + t502 * t461 - t501 * t462 - t572 * t472 + (-qJ(6) * t447 + t571) * t471 + t577 * t443 + (-qJ(6) * mrSges(7,2) - t574) * t442 + (Ifges(5,3) - t575) * t500;
t510 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t555;
t505 = (-mrSges(3,1) * t531 + mrSges(3,2) * t528) * qJD(1);
t491 = t542 - t569;
t460 = Ifges(5,5) * t502 + Ifges(5,6) * t501 + Ifges(5,3) * t517;
t392 = mrSges(6,2) * t420 + mrSges(7,2) * t406 - mrSges(6,3) * t408 - mrSges(7,3) * t411 - qJ(6) * t402 - t565 * t442 + t581 * t443 + t561 * t471 + t577 * t500 + t572 * t517;
t391 = -mrSges(6,1) * t420 - mrSges(7,1) * t411 + mrSges(7,2) * t405 + mrSges(6,3) * t409 - pkin(5) * t402 - t580 * t442 + t565 * t443 + t561 * t472 + t574 * t500 + t571 * t517;
t384 = qJDD(2) * mrSges(4,2) + qJD(2) * t512 + t504 * t555 + t539;
t383 = mrSges(4,2) * t507 - mrSges(4,3) * t506 + (t512 * t531 - t513 * t528) * qJD(1) + t543;
t382 = mrSges(5,2) * t434 - mrSges(5,3) * t417 + Ifges(5,1) * t470 + Ifges(5,4) * t469 + Ifges(5,5) * t500 - qJ(5) * t390 - t526 * t391 + t562 * t392 + t501 * t460 - t517 * t461;
t381 = -mrSges(5,1) * t434 + mrSges(5,3) * t418 + Ifges(5,4) * t470 + Ifges(5,2) * t469 + Ifges(5,6) * t500 - pkin(4) * t398 + qJ(5) * t546 + t562 * t391 + t526 * t392 - t502 * t460 + t517 * t462;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t548 - mrSges(2,2) * t544 + t528 * (mrSges(4,1) * t454 + mrSges(3,2) * t491 - mrSges(3,3) * t477 - mrSges(4,3) * t449 + pkin(3) * t385 - qJ(3) * t383 - t558 * qJD(2) + t564 * qJDD(2) + t506 * t578 + t566 * t507 + t559 * t554 + t535) + t531 * (-mrSges(3,1) * t491 - mrSges(4,1) * t453 + mrSges(4,2) * t449 + mrSges(3,3) * t478 - pkin(2) * t383 - pkin(3) * t537 - pkin(8) * t547 + t557 * qJD(2) + t563 * qJDD(2) - t530 * t381 - t527 * t382 + t566 * t506 + t507 * t576 - t559 * t555) + pkin(1) * (-m(3) * t491 + t568 * t507 + (-mrSges(3,2) + mrSges(4,3)) * t506 + (t556 * t531 + (-t510 + t513) * t528) * qJD(1) - t543) + pkin(7) * (t531 * (t536 - qJDD(2) * mrSges(3,2) + (mrSges(3,3) + mrSges(4,1)) * t507 - qJD(2) * t510 + m(3) * t478 + t505 * t554) + (-m(3) * t477 + t506 * mrSges(3,3) - t568 * qJDD(2) - t556 * qJD(2) + (t504 + t505) * t555 + t539) * t528); mrSges(3,1) * t477 - mrSges(3,2) * t478 + mrSges(4,2) * t454 - mrSges(4,3) * t453 + t530 * t382 - t527 * t381 - pkin(8) * t385 - pkin(2) * t384 + qJ(3) * t536 + (mrSges(4,1) * qJ(3) + t563) * t507 + t564 * t506 + t573 * qJDD(2) + (t558 * t528 - t557 * t531) * qJD(1); t384; t535; t398; t401;];
tauJ  = t1;
