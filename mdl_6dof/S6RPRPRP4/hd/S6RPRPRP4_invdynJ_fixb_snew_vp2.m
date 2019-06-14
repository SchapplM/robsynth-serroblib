% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-05-05 17:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:42:47
% EndTime: 2019-05-05 17:42:50
% DurationCPUTime: 1.78s
% Computational Cost: add. (6083->251), mult. (11849->291), div. (0->0), fcn. (6213->8), ass. (0->110)
t531 = Ifges(6,1) + Ifges(7,1);
t510 = Ifges(6,4) - Ifges(7,5);
t526 = Ifges(7,4) + Ifges(6,5);
t530 = Ifges(6,2) + Ifges(7,3);
t525 = Ifges(6,6) - Ifges(7,6);
t473 = sin(qJ(3));
t476 = cos(qJ(3));
t511 = Ifges(4,4) + Ifges(5,6);
t529 = t473 * t511 + t476 * (Ifges(4,2) + Ifges(5,3));
t528 = t473 * (Ifges(4,1) + Ifges(5,2)) + t476 * t511;
t527 = -2 * qJD(4);
t509 = Ifges(4,5) - Ifges(5,4);
t508 = Ifges(4,6) - Ifges(5,5);
t524 = Ifges(6,3) + Ifges(7,2);
t472 = sin(qJ(5));
t475 = cos(qJ(5));
t500 = qJD(1) * t476;
t442 = qJD(3) * t472 + t475 * t500;
t443 = qJD(3) * t475 - t472 * t500;
t501 = qJD(1) * t473;
t460 = qJD(5) + t501;
t523 = t530 * t442 - t510 * t443 - t525 * t460;
t522 = -t510 * t442 + t531 * t443 + t526 * t460;
t471 = cos(pkin(9));
t519 = pkin(1) * t471 + pkin(2);
t474 = sin(qJ(1));
t477 = cos(qJ(1));
t494 = t474 * g(1) - g(2) * t477;
t444 = qJDD(1) * pkin(1) + t494;
t479 = qJD(1) ^ 2;
t491 = -g(1) * t477 - g(2) * t474;
t448 = -pkin(1) * t479 + t491;
t470 = sin(pkin(9));
t409 = t470 * t444 + t471 * t448;
t391 = -pkin(2) * t479 + qJDD(1) * pkin(7) + t409;
t469 = -g(3) + qJDD(2);
t387 = t476 * t391 + t473 * t469;
t445 = (-t476 * pkin(3) - t473 * qJ(4)) * qJD(1);
t478 = qJD(3) ^ 2;
t382 = pkin(3) * t478 - qJDD(3) * qJ(4) + qJD(3) * t527 - t445 * t500 - t387;
t518 = t473 * (t529 * qJD(1) + t508 * qJD(3)) - t476 * (t528 * qJD(1) + t509 * qJD(3));
t517 = -pkin(3) - pkin(8);
t515 = pkin(7) * t479;
t514 = pkin(8) * t479;
t512 = -mrSges(6,3) - mrSges(7,2);
t507 = t469 * t476;
t499 = qJD(1) * qJD(3);
t495 = t473 * t499;
t450 = qJDD(1) * t476 - t495;
t457 = pkin(4) * t501 - qJD(3) * pkin(8);
t468 = t476 ^ 2;
t496 = t476 * t499;
t449 = qJDD(1) * t473 + t496;
t408 = t444 * t471 - t470 * t448;
t488 = -qJDD(1) * pkin(2) - t408;
t483 = pkin(3) * t495 + t501 * t527 + (-t449 - t496) * qJ(4) + t488;
t376 = -t457 * t501 + (-pkin(4) * t468 - pkin(7)) * t479 + t517 * t450 + t483;
t388 = t473 * t391;
t489 = -qJ(4) * t478 + t445 * t501 + qJDD(4) + t388;
t380 = pkin(4) * t449 + t517 * qJDD(3) + (-pkin(4) * t499 - t473 * t514 - t469) * t476 + t489;
t374 = t475 * t376 + t472 * t380;
t406 = qJD(5) * t443 + qJDD(3) * t472 + t475 * t450;
t417 = mrSges(6,1) * t460 - mrSges(6,3) * t443;
t441 = qJDD(5) + t449;
t412 = pkin(5) * t442 - qJ(6) * t443;
t458 = t460 ^ 2;
t368 = -pkin(5) * t458 + qJ(6) * t441 + 0.2e1 * qJD(6) * t460 - t412 * t442 + t374;
t418 = -mrSges(7,1) * t460 + mrSges(7,2) * t443;
t497 = m(7) * t368 + t441 * mrSges(7,3) + t460 * t418;
t413 = mrSges(7,1) * t442 - mrSges(7,3) * t443;
t504 = -mrSges(6,1) * t442 - mrSges(6,2) * t443 - t413;
t359 = m(6) * t374 - mrSges(6,2) * t441 + t406 * t512 - t417 * t460 + t442 * t504 + t497;
t373 = -t376 * t472 + t380 * t475;
t407 = -qJD(5) * t442 + qJDD(3) * t475 - t450 * t472;
t415 = -mrSges(6,2) * t460 - mrSges(6,3) * t442;
t369 = -pkin(5) * t441 - qJ(6) * t458 + t412 * t443 + qJDD(6) - t373;
t416 = -mrSges(7,2) * t442 + mrSges(7,3) * t460;
t492 = -m(7) * t369 + t441 * mrSges(7,1) + t460 * t416;
t361 = m(6) * t373 + mrSges(6,1) * t441 + t407 * t512 + t415 * t460 + t443 * t504 + t492;
t506 = t475 * t359 - t472 * t361;
t505 = t525 * t442 - t526 * t443 - t524 * t460;
t386 = -t388 + t507;
t446 = (mrSges(5,2) * t476 - mrSges(5,3) * t473) * qJD(1);
t447 = (-t476 * mrSges(4,1) + t473 * mrSges(4,2)) * qJD(1);
t454 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t500;
t455 = -mrSges(5,1) * t500 - qJD(3) * mrSges(5,3);
t355 = t359 * t472 + t361 * t475;
t383 = -qJDD(3) * pkin(3) + t489 - t507;
t485 = -m(5) * t383 - t449 * mrSges(5,1) - t355;
t350 = m(4) * t386 - mrSges(4,3) * t449 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t454 - t455) * qJD(3) + (-t446 - t447) * t501 + t485;
t453 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t501;
t456 = mrSges(5,1) * t501 + qJD(3) * mrSges(5,2);
t379 = pkin(4) * t450 + qJD(3) * t457 - t468 * t514 - t382;
t371 = -0.2e1 * qJD(6) * t443 + (t442 * t460 - t407) * qJ(6) + (t443 * t460 + t406) * pkin(5) + t379;
t365 = m(7) * t371 + t406 * mrSges(7,1) - mrSges(7,3) * t407 + t442 * t416 - t418 * t443;
t482 = m(6) * t379 + mrSges(6,1) * t406 + t407 * mrSges(6,2) + t415 * t442 + t443 * t417 + t365;
t481 = -m(5) * t382 + qJDD(3) * mrSges(5,3) + qJD(3) * t456 + t446 * t500 + t482;
t357 = -qJDD(3) * mrSges(4,2) + (mrSges(4,3) + mrSges(5,1)) * t450 + m(4) * t387 - qJD(3) * t453 + t481 + t447 * t500;
t493 = -t350 * t473 + t476 * t357;
t381 = -pkin(3) * t450 + t483 - t515;
t487 = m(5) * t381 + t450 * mrSges(5,2) - t456 * t501 + t506;
t390 = t488 - t515;
t484 = -m(4) * t390 + t450 * mrSges(4,1) + t454 * t500 - t487;
t364 = mrSges(7,2) * t407 + t413 * t443 - t492;
t480 = mrSges(6,1) * t373 - mrSges(7,1) * t369 - mrSges(6,2) * t374 + mrSges(7,3) * t368 - pkin(5) * t364 + qJ(6) * t497 - t523 * t443 + (-qJ(6) * t413 + t522) * t442 + t524 * t441 + t526 * t407 + (-mrSges(7,2) * qJ(6) - t525) * t406;
t354 = mrSges(6,2) * t379 + mrSges(7,2) * t369 - mrSges(6,3) * t373 - mrSges(7,3) * t371 - qJ(6) * t365 - t510 * t406 + t531 * t407 + t526 * t441 + t505 * t442 + t523 * t460;
t353 = -mrSges(6,1) * t379 - mrSges(7,1) * t371 + mrSges(7,2) * t368 + mrSges(6,3) * t374 - pkin(5) * t365 - t530 * t406 + t510 * t407 + t525 * t441 + t505 * t443 + t522 * t460;
t352 = qJDD(3) * mrSges(5,2) + qJD(3) * t455 + t446 * t501 - t485;
t351 = -mrSges(5,3) * t449 + t455 * t500 + t487;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t494 - mrSges(2,2) * t491 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t408 - mrSges(3,2) * t409 + t473 * (mrSges(5,1) * t383 + mrSges(4,2) * t390 - mrSges(4,3) * t386 - mrSges(5,3) * t381 + pkin(4) * t355 - qJ(4) * t351 + t480) + t476 * (-mrSges(4,1) * t390 - mrSges(5,1) * t382 + mrSges(5,2) * t381 + mrSges(4,3) * t387 - pkin(3) * t351 + pkin(4) * t482 - pkin(8) * t506 - t475 * t353 - t472 * t354) + pkin(2) * t484 + pkin(7) * t493 + pkin(1) * (t470 * (m(3) * t409 - mrSges(3,1) * t479 - qJDD(1) * mrSges(3,2) + t493) + t471 * (m(3) * t408 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t479 + t484)) + t529 * t450 + (t473 * t509 + t476 * t508) * qJDD(3) - t518 * qJD(3) + (t519 * (-mrSges(4,2) + mrSges(5,3)) + t528) * t449 + t519 * qJD(1) * (-t453 * t473 - t455 * t476); m(3) * t469 + t350 * t476 + t357 * t473; mrSges(4,1) * t386 - mrSges(4,2) * t387 + mrSges(5,2) * t383 - mrSges(5,3) * t382 + t475 * t354 - t472 * t353 - pkin(8) * t355 - pkin(3) * t352 + qJ(4) * t481 + (qJ(4) * mrSges(5,1) + t508) * t450 + t509 * t449 + (Ifges(4,3) + Ifges(5,1)) * qJDD(3) + t518 * qJD(1); t352; t480; t364;];
tauJ  = t1;
