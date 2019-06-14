% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-05-05 04:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:00:55
% EndTime: 2019-05-05 04:00:58
% DurationCPUTime: 1.85s
% Computational Cost: add. (6397->256), mult. (12567->292), div. (0->0), fcn. (7463->10), ass. (0->113)
t544 = Ifges(6,4) + Ifges(7,4);
t568 = Ifges(6,2) + Ifges(7,2);
t562 = Ifges(6,6) + Ifges(7,6);
t567 = Ifges(6,1) + Ifges(7,1);
t563 = Ifges(6,5) + Ifges(7,5);
t502 = sin(qJ(3));
t505 = cos(qJ(3));
t545 = Ifges(4,4) + Ifges(5,6);
t566 = t502 * t545 + t505 * (Ifges(4,2) + Ifges(5,3));
t565 = t502 * (Ifges(4,1) + Ifges(5,2)) + t505 * t545;
t564 = -2 * qJD(4);
t543 = Ifges(4,5) - Ifges(5,4);
t542 = Ifges(4,6) - Ifges(5,5);
t561 = Ifges(6,3) + Ifges(7,3);
t501 = sin(qJ(5));
t504 = cos(qJ(5));
t529 = qJD(2) * t505;
t471 = -qJD(3) * t501 - t504 * t529;
t472 = qJD(3) * t504 - t501 * t529;
t530 = qJD(2) * t502;
t489 = qJD(5) + t530;
t560 = t568 * t471 + t544 * t472 + t562 * t489;
t528 = qJD(2) * qJD(3);
t522 = t502 * t528;
t477 = qJDD(2) * t505 - t522;
t438 = qJD(5) * t471 + qJDD(3) * t504 - t477 * t501;
t440 = -mrSges(7,1) * t471 + mrSges(7,2) * t472;
t497 = sin(pkin(10));
t499 = cos(pkin(10));
t480 = g(1) * t497 - g(2) * t499;
t496 = -g(3) + qJDD(1);
t498 = sin(pkin(6));
t500 = cos(pkin(6));
t448 = -t480 * t498 + t496 * t500;
t523 = t505 * t528;
t476 = qJDD(2) * t502 + t523;
t481 = -g(1) * t499 - g(2) * t497;
t503 = sin(qJ(2));
t506 = cos(qJ(2));
t558 = t480 * t500 + t496 * t498;
t421 = t506 * t481 + t503 * t558;
t508 = qJD(2) ^ 2;
t418 = -pkin(2) * t508 + qJDD(2) * pkin(8) + t421;
t415 = t502 * t418;
t473 = (-pkin(3) * t505 - qJ(4) * t502) * qJD(2);
t507 = qJD(3) ^ 2;
t519 = -qJ(4) * t507 + t473 * t530 + qJDD(4) + t415;
t548 = pkin(9) * t508;
t550 = -pkin(3) - pkin(9);
t405 = pkin(4) * t476 + t550 * qJDD(3) + (-pkin(4) * t528 - t502 * t548 - t448) * t505 + t519;
t487 = pkin(4) * t530 - qJD(3) * pkin(9);
t495 = t505 ^ 2;
t420 = -t503 * t481 + t506 * t558;
t512 = -qJDD(2) * pkin(2) - t420;
t509 = pkin(3) * t522 + t530 * t564 + (-t476 - t523) * qJ(4) + t512;
t407 = -t487 * t530 + (-pkin(4) * t495 - pkin(8)) * t508 + t550 * t477 + t509;
t397 = t504 * t405 - t407 * t501;
t468 = qJDD(5) + t476;
t392 = -0.2e1 * qJD(6) * t472 + (t471 * t489 - t438) * qJ(6) + (t471 * t472 + t468) * pkin(5) + t397;
t443 = -mrSges(7,2) * t489 + mrSges(7,3) * t471;
t526 = m(7) * t392 + t468 * mrSges(7,1) + t489 * t443;
t389 = -mrSges(7,3) * t438 - t440 * t472 + t526;
t398 = t501 * t405 + t504 * t407;
t437 = -qJD(5) * t472 - qJDD(3) * t501 - t477 * t504;
t445 = pkin(5) * t489 - qJ(6) * t472;
t467 = t471 ^ 2;
t395 = -pkin(5) * t467 + qJ(6) * t437 + 0.2e1 * qJD(6) * t471 - t445 * t489 + t398;
t534 = -t544 * t471 - t567 * t472 - t563 * t489;
t559 = mrSges(6,1) * t397 + mrSges(7,1) * t392 - mrSges(6,2) * t398 - mrSges(7,2) * t395 + pkin(5) * t389 + t437 * t562 + t438 * t563 + t468 * t561 + t471 * t534 + t472 * t560;
t412 = t505 * t418 + t502 * t448;
t408 = pkin(3) * t507 - qJDD(3) * qJ(4) + qJD(3) * t564 - t473 * t529 - t412;
t554 = -mrSges(7,1) * t437 - t443 * t471;
t549 = pkin(8) * t508;
t417 = t512 - t549;
t482 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t530;
t483 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t529;
t484 = -mrSges(5,1) * t529 - qJD(3) * mrSges(5,3);
t410 = -pkin(3) * t477 + t509 - t549;
t485 = mrSges(5,1) * t530 + qJD(3) * mrSges(5,2);
t441 = -mrSges(6,1) * t471 + mrSges(6,2) * t472;
t444 = -mrSges(6,2) * t489 + mrSges(6,3) * t471;
t385 = m(6) * t397 + mrSges(6,1) * t468 + t444 * t489 + (-t440 - t441) * t472 + (-mrSges(6,3) - mrSges(7,3)) * t438 + t526;
t446 = mrSges(7,1) * t489 - mrSges(7,3) * t472;
t447 = mrSges(6,1) * t489 - mrSges(6,3) * t472;
t525 = m(7) * t395 + t437 * mrSges(7,3) + t471 * t440;
t387 = m(6) * t398 + mrSges(6,3) * t437 + t441 * t471 + (-t446 - t447) * t489 + (-mrSges(6,2) - mrSges(7,2)) * t468 + t525;
t536 = -t501 * t385 + t504 * t387;
t517 = m(5) * t410 + t477 * mrSges(5,2) - t485 * t530 + t536;
t553 = -m(4) * t417 + t477 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t476 + t483 * t529 + (-t482 * t502 - t484 * t505) * qJD(2) - t517;
t552 = -(mrSges(6,1) + mrSges(7,1)) * t437 - (t443 + t444) * t471;
t551 = t502 * (t566 * qJD(2) + t542 * qJD(3)) - t505 * (t565 * qJD(2) + t543 * qJD(3));
t539 = t448 * t505;
t535 = -t471 * t562 - t472 * t563 - t489 * t561;
t404 = pkin(4) * t477 + qJD(3) * t487 - t495 * t548 - t408;
t400 = -pkin(5) * t437 - qJ(6) * t467 + t445 * t472 + qJDD(6) + t404;
t524 = m(7) * t400 + t438 * mrSges(7,2) + t472 * t446;
t411 = -t415 + t539;
t474 = (mrSges(5,2) * t505 - mrSges(5,3) * t502) * qJD(2);
t475 = (-mrSges(4,1) * t505 + mrSges(4,2) * t502) * qJD(2);
t381 = t385 * t504 + t387 * t501;
t409 = -qJDD(3) * pkin(3) + t519 - t539;
t514 = -m(5) * t409 - t476 * mrSges(5,1) - t381;
t377 = m(4) * t411 - mrSges(4,3) * t476 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t483 - t484) * qJD(3) + (-t474 - t475) * t530 + t514;
t518 = -m(6) * t404 - t438 * mrSges(6,2) - t472 * t447 - t524;
t511 = -m(5) * t408 + qJDD(3) * mrSges(5,3) + qJD(3) * t485 + t474 * t529 - t518;
t383 = t511 + m(4) * t412 - qJD(3) * t482 - qJDD(3) * mrSges(4,2) + (mrSges(4,3) + mrSges(5,1)) * t477 + t475 * t529 + t552;
t521 = -t377 * t502 + t505 * t383;
t393 = t524 + t554;
t380 = mrSges(6,2) * t404 + mrSges(7,2) * t400 - mrSges(6,3) * t397 - mrSges(7,3) * t392 - qJ(6) * t389 + t544 * t437 + t567 * t438 + t563 * t468 - t535 * t471 - t560 * t489;
t379 = qJDD(3) * mrSges(5,2) + qJD(3) * t484 + t474 * t530 - t514;
t378 = -mrSges(5,3) * t476 + t484 * t529 + t517;
t376 = -mrSges(6,1) * t404 + mrSges(6,3) * t398 - mrSges(7,1) * t400 + mrSges(7,3) * t395 - pkin(5) * t393 + qJ(6) * t525 + (-qJ(6) * t446 - t534) * t489 + t535 * t472 + (-mrSges(7,2) * qJ(6) + t562) * t468 + t544 * t438 + t568 * t437;
t1 = [m(2) * t496 + t500 * (m(3) * t448 + t377 * t505 + t383 * t502) + (t503 * (m(3) * t421 - mrSges(3,1) * t508 - qJDD(2) * mrSges(3,2) + t521) + t506 * (m(3) * t420 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t508 + t553)) * t498; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t420 - mrSges(3,2) * t421 + t502 * (mrSges(5,1) * t409 + mrSges(4,2) * t417 - mrSges(4,3) * t411 - mrSges(5,3) * t410 + pkin(4) * t381 - qJ(4) * t378 + t559) + t505 * (-mrSges(4,1) * t417 + mrSges(4,3) * t412 - mrSges(5,1) * t408 + mrSges(5,2) * t410 - t501 * t380 - t504 * t376 - pkin(4) * (t518 - t552) - pkin(9) * t536 - pkin(3) * t378) + pkin(8) * t521 + t566 * t477 + (t502 * t543 + t505 * t542) * qJDD(3) - t551 * qJD(3) + t565 * t476 + t553 * pkin(2); mrSges(4,1) * t411 - mrSges(4,2) * t412 + mrSges(5,2) * t409 - mrSges(5,3) * t408 + t504 * t380 - t501 * t376 - pkin(9) * t381 - pkin(3) * t379 + qJ(4) * (-mrSges(6,1) * t437 - t444 * t471 + t511 + t554) + (mrSges(5,1) * qJ(4) + t542) * t477 + t543 * t476 + (Ifges(4,3) + Ifges(5,1)) * qJDD(3) + t551 * qJD(2); t379; t559; t393;];
tauJ  = t1;
