% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:38
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:34:48
% EndTime: 2019-05-06 01:34:56
% DurationCPUTime: 5.25s
% Computational Cost: add. (49215->301), mult. (116307->371), div. (0->0), fcn. (88275->10), ass. (0->128)
t565 = Ifges(6,4) + Ifges(7,4);
t574 = Ifges(6,2) + Ifges(7,2);
t570 = Ifges(6,6) + Ifges(7,6);
t526 = sin(pkin(10));
t527 = cos(pkin(10));
t530 = sin(qJ(3));
t534 = cos(qJ(3));
t544 = t526 * t534 + t527 * t530;
t559 = qJD(1) * t527;
t560 = qJD(1) * t526;
t514 = -t530 * t560 + t534 * t559;
t558 = qJD(3) * t514;
t501 = qJDD(1) * t544 + t558;
t515 = t544 * qJD(1);
t529 = sin(qJ(4));
t533 = cos(qJ(4));
t506 = qJD(3) * t529 + t515 * t533;
t474 = -qJD(4) * t506 + qJDD(3) * t533 - t501 * t529;
t505 = qJD(3) * t533 - t515 * t529;
t475 = qJD(4) * t505 + qJDD(3) * t529 + t501 * t533;
t528 = sin(qJ(5));
t532 = cos(qJ(5));
t477 = t505 * t532 - t506 * t528;
t439 = qJD(5) * t477 + t474 * t528 + t475 * t532;
t478 = t505 * t528 + t506 * t532;
t458 = -mrSges(7,1) * t477 + mrSges(7,2) * t478;
t537 = qJD(1) ^ 2;
t531 = sin(qJ(1));
t535 = cos(qJ(1));
t547 = -g(1) * t535 - g(2) * t531;
t516 = -pkin(1) * t537 + qJDD(1) * qJ(2) + t547;
t557 = qJD(1) * qJD(2);
t552 = -g(3) * t527 - 0.2e1 * t526 * t557;
t566 = pkin(2) * t537;
t487 = (-pkin(7) * qJDD(1) + t527 * t566 - t516) * t526 + t552;
t503 = -g(3) * t526 + (t516 + 0.2e1 * t557) * t527;
t525 = t527 ^ 2;
t556 = qJDD(1) * t527;
t488 = pkin(7) * t556 - t525 * t566 + t503;
t462 = t530 * t487 + t534 * t488;
t498 = -pkin(3) * t514 - pkin(8) * t515;
t536 = qJD(3) ^ 2;
t445 = -pkin(3) * t536 + qJDD(3) * pkin(8) + t498 * t514 + t462;
t553 = g(1) * t531 - t535 * g(2);
t546 = qJDD(2) - t553;
t561 = -t526 ^ 2 - t525;
t499 = (-pkin(2) * t527 - pkin(1)) * qJDD(1) + (pkin(7) * t561 - qJ(2)) * t537 + t546;
t511 = t515 * qJD(3);
t500 = -t530 * t526 * qJDD(1) + t534 * t556 - t511;
t448 = (-t501 - t558) * pkin(8) + (-t500 + t511) * pkin(3) + t499;
t428 = -t445 * t529 + t533 * t448;
t497 = qJDD(4) - t500;
t512 = qJD(4) - t514;
t424 = (t505 * t512 - t475) * pkin(9) + (t505 * t506 + t497) * pkin(4) + t428;
t429 = t533 * t445 + t529 * t448;
t486 = pkin(4) * t512 - pkin(9) * t506;
t504 = t505 ^ 2;
t426 = -pkin(4) * t504 + pkin(9) * t474 - t486 * t512 + t429;
t418 = t532 * t424 - t426 * t528;
t495 = qJDD(5) + t497;
t510 = qJD(5) + t512;
t413 = -0.2e1 * qJD(6) * t478 + (t477 * t510 - t439) * qJ(6) + (t477 * t478 + t495) * pkin(5) + t418;
t463 = -mrSges(7,2) * t510 + mrSges(7,3) * t477;
t555 = m(7) * t413 + t495 * mrSges(7,1) + t510 * t463;
t410 = -mrSges(7,3) * t439 - t458 * t478 + t555;
t419 = t528 * t424 + t532 * t426;
t438 = -qJD(5) * t478 + t474 * t532 - t475 * t528;
t465 = pkin(5) * t510 - qJ(6) * t478;
t476 = t477 ^ 2;
t416 = -pkin(5) * t476 + qJ(6) * t438 + 0.2e1 * qJD(6) * t477 - t465 * t510 + t419;
t571 = Ifges(6,5) + Ifges(7,5);
t572 = Ifges(6,1) + Ifges(7,1);
t562 = -t565 * t477 - t478 * t572 - t571 * t510;
t568 = t574 * t477 + t565 * t478 + t570 * t510;
t569 = Ifges(6,3) + Ifges(7,3);
t573 = mrSges(6,1) * t418 + mrSges(7,1) * t413 - mrSges(6,2) * t419 - mrSges(7,2) * t416 + pkin(5) * t410 + t438 * t570 + t439 * t571 + t562 * t477 + t478 * t568 + t495 * t569;
t459 = -mrSges(6,1) * t477 + mrSges(6,2) * t478;
t464 = -mrSges(6,2) * t510 + mrSges(6,3) * t477;
t403 = m(6) * t418 + mrSges(6,1) * t495 + t464 * t510 + (-t458 - t459) * t478 + (-mrSges(6,3) - mrSges(7,3)) * t439 + t555;
t466 = mrSges(7,1) * t510 - mrSges(7,3) * t478;
t467 = mrSges(6,1) * t510 - mrSges(6,3) * t478;
t554 = m(7) * t416 + t438 * mrSges(7,3) + t477 * t458;
t406 = m(6) * t419 + mrSges(6,3) * t438 + t459 * t477 + (-t466 - t467) * t510 + (-mrSges(6,2) - mrSges(7,2)) * t495 + t554;
t401 = t532 * t403 + t528 * t406;
t469 = Ifges(5,4) * t506 + Ifges(5,2) * t505 + Ifges(5,6) * t512;
t470 = Ifges(5,1) * t506 + Ifges(5,4) * t505 + Ifges(5,5) * t512;
t567 = mrSges(5,1) * t428 - mrSges(5,2) * t429 + Ifges(5,5) * t475 + Ifges(5,6) * t474 + Ifges(5,3) * t497 + pkin(4) * t401 + t506 * t469 - t505 * t470 + t573;
t496 = -mrSges(4,1) * t514 + mrSges(4,2) * t515;
t508 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t515;
t479 = -mrSges(5,1) * t505 + mrSges(5,2) * t506;
t482 = -mrSges(5,2) * t512 + mrSges(5,3) * t505;
t398 = m(5) * t428 + mrSges(5,1) * t497 - mrSges(5,3) * t475 - t479 * t506 + t482 * t512 + t401;
t483 = mrSges(5,1) * t512 - mrSges(5,3) * t506;
t548 = -t403 * t528 + t532 * t406;
t399 = m(5) * t429 - mrSges(5,2) * t497 + mrSges(5,3) * t474 + t479 * t505 - t483 * t512 + t548;
t549 = -t398 * t529 + t533 * t399;
t392 = m(4) * t462 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t500 - qJD(3) * t508 + t496 * t514 + t549;
t461 = t487 * t534 - t530 * t488;
t507 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t514;
t444 = -qJDD(3) * pkin(3) - pkin(8) * t536 + t515 * t498 - t461;
t427 = -pkin(4) * t474 - pkin(9) * t504 + t506 * t486 + t444;
t421 = -pkin(5) * t438 - qJ(6) * t476 + t465 * t478 + qJDD(6) + t427;
t414 = m(7) * t421 - t438 * mrSges(7,1) + t439 * mrSges(7,2) - t477 * t463 + t478 * t466;
t541 = m(6) * t427 - t438 * mrSges(6,1) + mrSges(6,2) * t439 - t477 * t464 + t467 * t478 + t414;
t539 = -m(5) * t444 + t474 * mrSges(5,1) - mrSges(5,2) * t475 + t505 * t482 - t483 * t506 - t541;
t408 = m(4) * t461 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t501 + qJD(3) * t507 - t496 * t515 + t539;
t564 = t530 * t392 + t534 * t408;
t393 = t533 * t398 + t529 * t399;
t563 = -t477 * t570 - t478 * t571 - t510 * t569;
t550 = t534 * t392 - t530 * t408;
t545 = -t527 * mrSges(3,1) + t526 * mrSges(3,2);
t543 = mrSges(3,3) * qJDD(1) + t537 * t545;
t542 = m(4) * t499 - mrSges(4,1) * t500 + mrSges(4,2) * t501 - t507 * t514 + t508 * t515 + t393;
t518 = (Ifges(3,5) * t526 + Ifges(3,6) * t527) * qJD(1);
t513 = -qJDD(1) * pkin(1) - qJ(2) * t537 + t546;
t502 = -t516 * t526 + t552;
t491 = Ifges(4,1) * t515 + Ifges(4,4) * t514 + Ifges(4,5) * qJD(3);
t490 = Ifges(4,4) * t515 + Ifges(4,2) * t514 + Ifges(4,6) * qJD(3);
t489 = Ifges(4,5) * t515 + Ifges(4,6) * t514 + Ifges(4,3) * qJD(3);
t468 = Ifges(5,5) * t506 + Ifges(5,6) * t505 + Ifges(5,3) * t512;
t400 = mrSges(6,2) * t427 + mrSges(7,2) * t421 - mrSges(6,3) * t418 - mrSges(7,3) * t413 - qJ(6) * t410 + t565 * t438 + t439 * t572 - t563 * t477 + t571 * t495 - t568 * t510;
t394 = -mrSges(6,1) * t427 + mrSges(6,3) * t419 - mrSges(7,1) * t421 + mrSges(7,3) * t416 - pkin(5) * t414 + qJ(6) * t554 + (-qJ(6) * t466 - t562) * t510 + (-qJ(6) * mrSges(7,2) + t570) * t495 + t563 * t478 + t565 * t439 + t574 * t438;
t389 = mrSges(3,3) * t537 * t561 + m(3) * t513 + qJDD(1) * t545 + t542;
t388 = mrSges(5,2) * t444 - mrSges(5,3) * t428 + Ifges(5,1) * t475 + Ifges(5,4) * t474 + Ifges(5,5) * t497 - pkin(9) * t401 - t394 * t528 + t400 * t532 + t468 * t505 - t469 * t512;
t387 = -mrSges(5,1) * t444 + mrSges(5,3) * t429 + Ifges(5,4) * t475 + Ifges(5,2) * t474 + Ifges(5,6) * t497 - pkin(4) * t541 + pkin(9) * t548 + t532 * t394 + t528 * t400 - t506 * t468 + t512 * t470;
t386 = -mrSges(4,1) * t499 + mrSges(4,3) * t462 + Ifges(4,4) * t501 + Ifges(4,2) * t500 + Ifges(4,6) * qJDD(3) - pkin(3) * t393 + qJD(3) * t491 - t515 * t489 - t567;
t385 = mrSges(4,2) * t499 - mrSges(4,3) * t461 + Ifges(4,1) * t501 + Ifges(4,4) * t500 + Ifges(4,5) * qJDD(3) - pkin(8) * t393 - qJD(3) * t490 - t387 * t529 + t388 * t533 + t489 * t514;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t553 - mrSges(2,2) * t547 + t526 * (t518 * t559 + mrSges(3,2) * t513 - mrSges(3,3) * t502 + t534 * t385 - t530 * t386 - pkin(7) * t564 + (Ifges(3,1) * t526 + Ifges(3,4) * t527) * qJDD(1)) + t527 * (-t518 * t560 - mrSges(3,1) * t513 + mrSges(3,3) * t503 + t530 * t385 + t534 * t386 - pkin(2) * t542 + pkin(7) * t550 + (Ifges(3,4) * t526 + Ifges(3,2) * t527) * qJDD(1)) - pkin(1) * t389 + qJ(2) * ((m(3) * t503 + t527 * t543 + t550) * t527 + (-m(3) * t502 + t526 * t543 - t564) * t526); t389; mrSges(4,1) * t461 - mrSges(4,2) * t462 + Ifges(4,5) * t501 + Ifges(4,6) * t500 + Ifges(4,3) * qJDD(3) + pkin(3) * t539 + pkin(8) * t549 + t533 * t387 + t529 * t388 + t515 * t490 - t514 * t491; t567; t573; t414;];
tauJ  = t1;
