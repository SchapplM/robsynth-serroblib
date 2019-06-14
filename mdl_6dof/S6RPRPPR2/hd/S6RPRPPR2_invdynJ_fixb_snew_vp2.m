% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-05-05 16:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:33:01
% EndTime: 2019-05-05 16:33:03
% DurationCPUTime: 1.85s
% Computational Cost: add. (11775->275), mult. (25700->332), div. (0->0), fcn. (16125->10), ass. (0->115)
t515 = -2 * qJD(4);
t514 = Ifges(5,1) + Ifges(6,2);
t513 = -Ifges(6,1) - Ifges(5,3);
t510 = Ifges(5,4) + Ifges(6,6);
t509 = Ifges(5,5) - Ifges(6,4);
t512 = -Ifges(5,2) - Ifges(6,3);
t508 = Ifges(5,6) - Ifges(6,5);
t477 = sin(qJ(1));
t480 = cos(qJ(1));
t495 = t477 * g(1) - g(2) * t480;
t457 = qJDD(1) * pkin(1) + t495;
t482 = qJD(1) ^ 2;
t492 = -g(1) * t480 - g(2) * t477;
t459 = -pkin(1) * t482 + t492;
t473 = sin(pkin(9));
t474 = cos(pkin(9));
t428 = t473 * t457 + t474 * t459;
t414 = -pkin(2) * t482 + qJDD(1) * pkin(7) + t428;
t471 = -g(3) + qJDD(2);
t476 = sin(qJ(3));
t479 = cos(qJ(3));
t402 = -t414 * t476 + t479 * t471;
t497 = qJD(1) * qJD(3);
t496 = t479 * t497;
t460 = qJDD(1) * t476 + t496;
t387 = (-t460 + t496) * qJ(4) + (t476 * t479 * t482 + qJDD(3)) * pkin(3) + t402;
t403 = t479 * t414 + t476 * t471;
t461 = qJDD(1) * t479 - t476 * t497;
t501 = qJD(1) * t476;
t462 = qJD(3) * pkin(3) - qJ(4) * t501;
t470 = t479 ^ 2;
t388 = -pkin(3) * t470 * t482 + qJ(4) * t461 - qJD(3) * t462 + t403;
t472 = sin(pkin(10));
t507 = cos(pkin(10));
t446 = (t472 * t479 + t476 * t507) * qJD(1);
t381 = t387 * t507 - t472 * t388 + t446 * t515;
t511 = -2 * qJD(5);
t500 = qJD(1) * t479;
t445 = t472 * t501 - t500 * t507;
t419 = mrSges(5,1) * t445 + mrSges(5,2) * t446;
t430 = t460 * t507 + t472 * t461;
t434 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t445;
t436 = mrSges(6,1) * t445 - qJD(3) * mrSges(6,3);
t418 = pkin(4) * t445 - qJ(5) * t446;
t481 = qJD(3) ^ 2;
t378 = -qJDD(3) * pkin(4) - t481 * qJ(5) + t446 * t418 + qJDD(5) - t381;
t499 = qJD(3) * t445;
t373 = (t445 * t446 - qJDD(3)) * pkin(8) + (t430 + t499) * pkin(5) + t378;
t429 = t460 * t472 - t461 * t507;
t438 = pkin(5) * t446 - qJD(3) * pkin(8);
t444 = t445 ^ 2;
t427 = t457 * t474 - t473 * t459;
t490 = -qJDD(1) * pkin(2) - t427;
t389 = -pkin(3) * t461 + qJDD(4) + t462 * t501 + (-qJ(4) * t470 - pkin(7)) * t482 + t490;
t484 = (-t430 + t499) * qJ(5) + t389 + (qJD(3) * pkin(4) + t511) * t446;
t376 = -pkin(5) * t444 - t438 * t446 + (pkin(4) + pkin(8)) * t429 + t484;
t475 = sin(qJ(6));
t478 = cos(qJ(6));
t371 = t373 * t478 - t376 * t475;
t431 = -qJD(3) * t475 + t445 * t478;
t398 = qJD(6) * t431 + qJDD(3) * t478 + t429 * t475;
t432 = qJD(3) * t478 + t445 * t475;
t399 = -mrSges(7,1) * t431 + mrSges(7,2) * t432;
t443 = qJD(6) + t446;
t404 = -mrSges(7,2) * t443 + mrSges(7,3) * t431;
t426 = qJDD(6) + t430;
t368 = m(7) * t371 + mrSges(7,1) * t426 - mrSges(7,3) * t398 - t399 * t432 + t404 * t443;
t372 = t373 * t475 + t376 * t478;
t397 = -qJD(6) * t432 - qJDD(3) * t475 + t429 * t478;
t405 = mrSges(7,1) * t443 - mrSges(7,3) * t432;
t369 = m(7) * t372 - mrSges(7,2) * t426 + mrSges(7,3) * t397 + t399 * t431 - t405 * t443;
t360 = t368 * t478 + t369 * t475;
t420 = -mrSges(6,2) * t445 - mrSges(6,3) * t446;
t487 = -m(6) * t378 - t430 * mrSges(6,1) - t446 * t420 - t360;
t356 = m(5) * t381 - mrSges(5,3) * t430 - t419 * t446 + (mrSges(5,1) - mrSges(6,2)) * qJDD(3) + (t434 - t436) * qJD(3) + t487;
t441 = t445 * t515;
t505 = t472 * t387 + t507 * t388;
t382 = t441 + t505;
t435 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t446;
t489 = pkin(4) * t481 - qJDD(3) * qJ(5) - t505;
t377 = qJD(3) * t511 + ((2 * qJD(4)) + t418) * t445 + t489;
t437 = mrSges(6,1) * t446 + qJD(3) * mrSges(6,2);
t375 = -pkin(5) * t429 - pkin(8) * t444 - t418 * t445 + t441 + ((2 * qJD(5)) + t438) * qJD(3) - t489;
t488 = -m(7) * t375 + mrSges(7,1) * t397 - t398 * mrSges(7,2) + t404 * t431 - t432 * t405;
t485 = -m(6) * t377 + qJDD(3) * mrSges(6,3) + qJD(3) * t437 - t488;
t365 = m(5) * t382 - qJDD(3) * mrSges(5,2) - qJD(3) * t435 + (-t419 - t420) * t445 + (-mrSges(5,3) - mrSges(6,1)) * t429 + t485;
t354 = t507 * t356 + t472 * t365;
t506 = -t475 * t368 + t478 * t369;
t504 = t513 * qJD(3) + t508 * t445 - t509 * t446;
t503 = t508 * qJD(3) + t512 * t445 + t510 * t446;
t502 = t509 * qJD(3) - t510 * t445 + t514 * t446;
t458 = (-mrSges(4,1) * t479 + mrSges(4,2) * t476) * qJD(1);
t464 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t500;
t352 = m(4) * t402 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t460 + qJD(3) * t464 - t458 * t501 + t354;
t463 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t501;
t493 = -t356 * t472 + t507 * t365;
t353 = m(4) * t403 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t461 - qJD(3) * t463 + t458 * t500 + t493;
t494 = -t352 * t476 + t479 * t353;
t380 = pkin(4) * t429 + t484;
t359 = m(6) * t380 - t429 * mrSges(6,2) - t430 * mrSges(6,3) - t445 * t436 - t446 * t437 + t506;
t391 = Ifges(7,4) * t432 + Ifges(7,2) * t431 + Ifges(7,6) * t443;
t392 = Ifges(7,1) * t432 + Ifges(7,4) * t431 + Ifges(7,5) * t443;
t486 = mrSges(7,1) * t371 - mrSges(7,2) * t372 + Ifges(7,5) * t398 + Ifges(7,6) * t397 + Ifges(7,3) * t426 + t432 * t391 - t431 * t392;
t357 = m(5) * t389 + t429 * mrSges(5,1) + mrSges(5,2) * t430 + t445 * t434 + t435 * t446 + t359;
t413 = -pkin(7) * t482 + t490;
t483 = -m(4) * t413 + t461 * mrSges(4,1) - mrSges(4,2) * t460 - t463 * t501 + t464 * t500 - t357;
t452 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t476 + Ifges(4,4) * t479) * qJD(1);
t451 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t476 + Ifges(4,2) * t479) * qJD(1);
t390 = Ifges(7,5) * t432 + Ifges(7,6) * t431 + Ifges(7,3) * t443;
t362 = mrSges(7,2) * t375 - mrSges(7,3) * t371 + Ifges(7,1) * t398 + Ifges(7,4) * t397 + Ifges(7,5) * t426 + t390 * t431 - t391 * t443;
t361 = -mrSges(7,1) * t375 + mrSges(7,3) * t372 + Ifges(7,4) * t398 + Ifges(7,2) * t397 + Ifges(7,6) * t426 - t390 * t432 + t392 * t443;
t358 = qJDD(3) * mrSges(6,2) + qJD(3) * t436 - t487;
t350 = mrSges(6,1) * t378 + mrSges(5,2) * t389 - mrSges(5,3) * t381 - mrSges(6,3) * t380 + pkin(5) * t360 - qJ(5) * t359 - t503 * qJD(3) + t509 * qJDD(3) - t510 * t429 + t514 * t430 + t504 * t445 + t486;
t349 = -mrSges(5,1) * t389 - mrSges(6,1) * t377 + mrSges(6,2) * t380 + mrSges(5,3) * t382 - pkin(4) * t359 - pkin(5) * t488 - pkin(8) * t506 + t502 * qJD(3) + t508 * qJDD(3) - t478 * t361 - t475 * t362 + t512 * t429 + t510 * t430 + t504 * t446;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t495 - mrSges(2,2) * t492 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t427 - mrSges(3,2) * t428 + t476 * (mrSges(4,2) * t413 - mrSges(4,3) * t402 + Ifges(4,1) * t460 + Ifges(4,4) * t461 + Ifges(4,5) * qJDD(3) - qJ(4) * t354 - qJD(3) * t451 - t472 * t349 + t350 * t507) + t479 * (-mrSges(4,1) * t413 + mrSges(4,3) * t403 + Ifges(4,4) * t460 + Ifges(4,2) * t461 + Ifges(4,6) * qJDD(3) - pkin(3) * t357 + qJ(4) * t493 + qJD(3) * t452 + t349 * t507 + t472 * t350) + pkin(2) * t483 + pkin(7) * t494 + pkin(1) * (t473 * (m(3) * t428 - mrSges(3,1) * t482 - qJDD(1) * mrSges(3,2) + t494) + t474 * (m(3) * t427 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t482 + t483)); m(3) * t471 + t352 * t479 + t353 * t476; -t475 * t361 + t478 * t362 + qJ(5) * t485 + Ifges(4,5) * t460 + Ifges(4,6) * t461 + mrSges(4,1) * t402 - mrSges(4,2) * t403 + mrSges(5,1) * t381 - mrSges(5,2) * t382 - mrSges(6,3) * t377 + mrSges(6,2) * t378 - pkin(8) * t360 - pkin(4) * t358 + pkin(3) * t354 + t503 * t446 + (-qJ(5) * t420 + t502) * t445 + t509 * t430 + (-mrSges(6,1) * qJ(5) - t508) * t429 + (t451 * t476 - t452 * t479) * qJD(1) + (Ifges(4,3) - t513) * qJDD(3); t357; t358; t486;];
tauJ  = t1;
