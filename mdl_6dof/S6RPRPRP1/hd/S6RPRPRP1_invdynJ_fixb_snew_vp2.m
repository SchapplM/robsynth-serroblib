% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-05-05 17:31
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:28:45
% EndTime: 2019-05-05 17:28:48
% DurationCPUTime: 2.02s
% Computational Cost: add. (15426->270), mult. (33190->332), div. (0->0), fcn. (21590->10), ass. (0->111)
t503 = -2 * qJD(4);
t502 = Ifges(6,1) + Ifges(7,1);
t497 = Ifges(6,4) + Ifges(7,4);
t496 = Ifges(6,5) + Ifges(7,5);
t501 = Ifges(6,2) + Ifges(7,2);
t495 = Ifges(6,6) + Ifges(7,6);
t500 = Ifges(6,3) + Ifges(7,3);
t465 = sin(qJ(1));
t468 = cos(qJ(1));
t482 = g(1) * t465 - g(2) * t468;
t446 = qJDD(1) * pkin(1) + t482;
t470 = qJD(1) ^ 2;
t476 = -g(1) * t468 - g(2) * t465;
t448 = -pkin(1) * t470 + t476;
t460 = sin(pkin(9));
t462 = cos(pkin(9));
t424 = t446 * t460 + t448 * t462;
t416 = -pkin(2) * t470 + qJDD(1) * pkin(7) + t424;
t458 = -g(3) + qJDD(2);
t464 = sin(qJ(3));
t467 = cos(qJ(3));
t404 = -t464 * t416 + t458 * t467;
t486 = qJD(1) * qJD(3);
t483 = t467 * t486;
t449 = qJDD(1) * t464 + t483;
t383 = (-t449 + t483) * qJ(4) + (t464 * t467 * t470 + qJDD(3)) * pkin(3) + t404;
t405 = t416 * t467 + t458 * t464;
t450 = qJDD(1) * t467 - t464 * t486;
t489 = qJD(1) * t464;
t451 = qJD(3) * pkin(3) - qJ(4) * t489;
t457 = t467 ^ 2;
t384 = -pkin(3) * t457 * t470 + qJ(4) * t450 - qJD(3) * t451 + t405;
t459 = sin(pkin(10));
t461 = cos(pkin(10));
t436 = (t459 * t467 + t461 * t464) * qJD(1);
t375 = t383 * t461 - t384 * t459 + t436 * t503;
t435 = (t459 * t464 - t461 * t467) * qJD(1);
t426 = t449 * t461 + t450 * t459;
t463 = sin(qJ(5));
t466 = cos(qJ(5));
t428 = qJD(3) * t466 - t436 * t463;
t399 = qJD(5) * t428 + qJDD(3) * t463 + t426 * t466;
t429 = qJD(3) * t463 + t436 * t466;
t401 = -mrSges(7,1) * t428 + mrSges(7,2) * t429;
t376 = t383 * t459 + t384 * t461 + t435 * t503;
t419 = pkin(4) * t435 - pkin(8) * t436;
t469 = qJD(3) ^ 2;
t374 = -pkin(4) * t469 + qJDD(3) * pkin(8) - t419 * t435 + t376;
t423 = t462 * t446 - t448 * t460;
t474 = -qJDD(1) * pkin(2) - t423;
t386 = -t450 * pkin(3) + qJDD(4) + t451 * t489 + (-qJ(4) * t457 - pkin(7)) * t470 + t474;
t425 = -t449 * t459 + t450 * t461;
t379 = (qJD(3) * t435 - t426) * pkin(8) + (qJD(3) * t436 - t425) * pkin(4) + t386;
t369 = -t463 * t374 + t379 * t466;
t422 = qJDD(5) - t425;
t434 = qJD(5) + t435;
t366 = -0.2e1 * qJD(6) * t429 + (t428 * t434 - t399) * qJ(6) + (t428 * t429 + t422) * pkin(5) + t369;
t406 = -mrSges(7,2) * t434 + mrSges(7,3) * t428;
t485 = m(7) * t366 + mrSges(7,1) * t422 + t406 * t434;
t363 = -mrSges(7,3) * t399 - t429 * t401 + t485;
t370 = t374 * t466 + t379 * t463;
t398 = -qJD(5) * t429 + qJDD(3) * t466 - t426 * t463;
t408 = pkin(5) * t434 - qJ(6) * t429;
t427 = t428 ^ 2;
t368 = -pkin(5) * t427 + qJ(6) * t398 + 0.2e1 * qJD(6) * t428 - t408 * t434 + t370;
t491 = t428 * t497 + t429 * t502 + t434 * t496;
t492 = -t428 * t501 - t429 * t497 - t434 * t495;
t499 = mrSges(6,1) * t369 + mrSges(7,1) * t366 - mrSges(6,2) * t370 - mrSges(7,2) * t368 + pkin(5) * t363 + t398 * t495 + t399 * t496 + t422 * t500 - t428 * t491 - t429 * t492;
t498 = -mrSges(6,2) - mrSges(7,2);
t418 = mrSges(5,1) * t435 + mrSges(5,2) * t436;
t431 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t436;
t402 = -mrSges(6,1) * t428 + mrSges(6,2) * t429;
t407 = -mrSges(6,2) * t434 + mrSges(6,3) * t428;
t357 = m(6) * t369 + t422 * mrSges(6,1) + t434 * t407 + (-t401 - t402) * t429 + (-mrSges(6,3) - mrSges(7,3)) * t399 + t485;
t484 = m(7) * t368 + mrSges(7,3) * t398 + t401 * t428;
t409 = mrSges(7,1) * t434 - mrSges(7,3) * t429;
t490 = -mrSges(6,1) * t434 + mrSges(6,3) * t429 - t409;
t360 = m(6) * t370 + t398 * mrSges(6,3) + t428 * t402 + t422 * t498 + t434 * t490 + t484;
t479 = -t357 * t463 + t360 * t466;
t352 = m(5) * t376 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t425 - qJD(3) * t431 - t418 * t435 + t479;
t430 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t435;
t373 = -qJDD(3) * pkin(4) - pkin(8) * t469 + t419 * t436 - t375;
t371 = -pkin(5) * t398 - qJ(6) * t427 + t408 * t429 + qJDD(6) + t373;
t477 = -m(7) * t371 + mrSges(7,1) * t398 + t406 * t428;
t472 = -m(6) * t373 + mrSges(6,1) * t398 + t399 * t498 + t407 * t428 + t429 * t490 + t477;
t362 = m(5) * t375 + qJDD(3) * mrSges(5,1) - t426 * mrSges(5,3) + qJD(3) * t430 - t436 * t418 + t472;
t348 = t352 * t459 + t362 * t461;
t355 = t357 * t466 + t360 * t463;
t493 = -t428 * t495 - t429 * t496 - t434 * t500;
t488 = qJD(1) * t467;
t447 = (-mrSges(4,1) * t467 + mrSges(4,2) * t464) * qJD(1);
t453 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t488;
t346 = m(4) * t404 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t449 + qJD(3) * t453 - t447 * t489 + t348;
t452 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t489;
t480 = t352 * t461 - t362 * t459;
t347 = m(4) * t405 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t450 - qJD(3) * t452 + t447 * t488 + t480;
t481 = -t346 * t464 + t347 * t467;
t354 = m(5) * t386 - mrSges(5,1) * t425 + mrSges(5,2) * t426 + t430 * t435 + t431 * t436 + t355;
t415 = -t470 * pkin(7) + t474;
t471 = -m(4) * t415 + mrSges(4,1) * t450 - mrSges(4,2) * t449 - t452 * t489 + t453 * t488 - t354;
t442 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t464 + Ifges(4,4) * t467) * qJD(1);
t441 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t464 + Ifges(4,2) * t467) * qJD(1);
t414 = Ifges(5,1) * t436 - Ifges(5,4) * t435 + Ifges(5,5) * qJD(3);
t413 = Ifges(5,4) * t436 - Ifges(5,2) * t435 + Ifges(5,6) * qJD(3);
t412 = Ifges(5,5) * t436 - Ifges(5,6) * t435 + Ifges(5,3) * qJD(3);
t364 = t399 * mrSges(7,2) + t429 * t409 - t477;
t353 = mrSges(6,2) * t373 + mrSges(7,2) * t371 - mrSges(6,3) * t369 - mrSges(7,3) * t366 - qJ(6) * t363 + t398 * t497 + t399 * t502 + t422 * t496 - t428 * t493 + t434 * t492;
t349 = -mrSges(6,1) * t373 + mrSges(6,3) * t370 - mrSges(7,1) * t371 + mrSges(7,3) * t368 - pkin(5) * t364 + qJ(6) * t484 + (-qJ(6) * t409 + t491) * t434 + t493 * t429 + (-mrSges(7,2) * qJ(6) + t495) * t422 + t497 * t399 + t501 * t398;
t344 = -mrSges(5,1) * t386 + mrSges(5,3) * t376 + Ifges(5,4) * t426 + Ifges(5,2) * t425 + Ifges(5,6) * qJDD(3) - pkin(4) * t355 + qJD(3) * t414 - t436 * t412 - t499;
t343 = mrSges(5,2) * t386 - mrSges(5,3) * t375 + Ifges(5,1) * t426 + Ifges(5,4) * t425 + Ifges(5,5) * qJDD(3) - pkin(8) * t355 - qJD(3) * t413 - t349 * t463 + t353 * t466 - t412 * t435;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t482 - mrSges(2,2) * t476 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t423 - mrSges(3,2) * t424 + t464 * (mrSges(4,2) * t415 - mrSges(4,3) * t404 + Ifges(4,1) * t449 + Ifges(4,4) * t450 + Ifges(4,5) * qJDD(3) - qJ(4) * t348 - qJD(3) * t441 + t343 * t461 - t344 * t459) + t467 * (-mrSges(4,1) * t415 + mrSges(4,3) * t405 + Ifges(4,4) * t449 + Ifges(4,2) * t450 + Ifges(4,6) * qJDD(3) - pkin(3) * t354 + qJ(4) * t480 + qJD(3) * t442 + t459 * t343 + t461 * t344) + pkin(2) * t471 + pkin(7) * t481 + pkin(1) * (t460 * (m(3) * t424 - mrSges(3,1) * t470 - qJDD(1) * mrSges(3,2) + t481) + t462 * (m(3) * t423 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t470 + t471)); m(3) * t458 + t346 * t467 + t347 * t464; Ifges(4,5) * t449 + Ifges(4,6) * t450 + mrSges(4,1) * t404 - mrSges(4,2) * t405 + Ifges(5,5) * t426 + Ifges(5,6) * t425 + t436 * t413 + t435 * t414 + mrSges(5,1) * t375 - mrSges(5,2) * t376 + t463 * t353 + t466 * t349 + pkin(4) * t472 + pkin(8) * t479 + pkin(3) * t348 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t441 * t464 - t442 * t467) * qJD(1); t354; t499; t364;];
tauJ  = t1;
