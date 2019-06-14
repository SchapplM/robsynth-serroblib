% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRP3
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
% Datum: 2019-05-05 17:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:37:48
% EndTime: 2019-05-05 17:37:52
% DurationCPUTime: 2.40s
% Computational Cost: add. (17505->268), mult. (36013->328), div. (0->0), fcn. (23289->10), ass. (0->109)
t483 = Ifges(6,1) + Ifges(7,1);
t473 = Ifges(6,4) - Ifges(7,5);
t479 = -Ifges(6,5) - Ifges(7,4);
t482 = Ifges(6,2) + Ifges(7,3);
t471 = Ifges(6,6) - Ifges(7,6);
t442 = sin(pkin(10));
t444 = cos(pkin(10));
t447 = sin(qJ(3));
t466 = qJD(1) * t447;
t424 = qJD(3) * t444 - t442 * t466;
t425 = qJD(3) * t442 + t444 * t466;
t446 = sin(qJ(5));
t475 = cos(qJ(5));
t398 = -t475 * t424 + t425 * t446;
t449 = cos(qJ(3));
t464 = qJD(1) * qJD(3);
t462 = t449 * t464;
t432 = qJDD(1) * t447 + t462;
t407 = qJDD(3) * t444 - t432 * t442;
t408 = qJDD(3) * t442 + t432 * t444;
t362 = -t398 * qJD(5) + t446 * t407 + t475 * t408;
t399 = t446 * t424 + t475 * t425;
t378 = mrSges(7,1) * t398 - mrSges(7,3) * t399;
t448 = sin(qJ(1));
t450 = cos(qJ(1));
t461 = t448 * g(1) - g(2) * t450;
t428 = qJDD(1) * pkin(1) + t461;
t452 = qJD(1) ^ 2;
t456 = -g(1) * t450 - g(2) * t448;
t431 = -pkin(1) * t452 + t456;
t443 = sin(pkin(9));
t445 = cos(pkin(9));
t400 = t445 * t428 - t443 * t431;
t391 = -qJDD(1) * pkin(2) - t452 * pkin(7) - t400;
t439 = t447 * t464;
t433 = qJDD(1) * t449 - t439;
t374 = (-t432 - t462) * qJ(4) + (-t433 + t439) * pkin(3) + t391;
t401 = t443 * t428 + t445 * t431;
t392 = -pkin(2) * t452 + qJDD(1) * pkin(7) + t401;
t441 = -g(3) + qJDD(2);
t384 = t449 * t392 + t447 * t441;
t429 = (-pkin(3) * t449 - qJ(4) * t447) * qJD(1);
t451 = qJD(3) ^ 2;
t465 = qJD(1) * t449;
t380 = -pkin(3) * t451 + qJDD(3) * qJ(4) + t429 * t465 + t384;
t355 = -0.2e1 * qJD(4) * t425 + t444 * t374 - t442 * t380;
t352 = (-t424 * t465 - t408) * pkin(8) + (t424 * t425 - t433) * pkin(4) + t355;
t356 = 0.2e1 * qJD(4) * t424 + t442 * t374 + t444 * t380;
t409 = -pkin(4) * t465 - pkin(8) * t425;
t423 = t424 ^ 2;
t354 = -pkin(4) * t423 + pkin(8) * t407 + t409 * t465 + t356;
t347 = t475 * t352 - t446 * t354;
t377 = pkin(5) * t398 - qJ(6) * t399;
t427 = qJDD(5) - t433;
t437 = qJD(5) - t465;
t436 = t437 ^ 2;
t346 = -t427 * pkin(5) - t436 * qJ(6) + t399 * t377 + qJDD(6) - t347;
t388 = -mrSges(7,2) * t398 + mrSges(7,3) * t437;
t457 = -m(7) * t346 + t427 * mrSges(7,1) + t437 * t388;
t342 = t362 * mrSges(7,2) + t399 * t378 - t457;
t348 = t446 * t352 + t475 * t354;
t345 = -pkin(5) * t436 + qJ(6) * t427 + 0.2e1 * qJD(6) * t437 - t377 * t398 + t348;
t361 = qJD(5) * t399 - t475 * t407 + t408 * t446;
t387 = -mrSges(7,1) * t437 + mrSges(7,2) * t399;
t463 = m(7) * t345 + t427 * mrSges(7,3) + t437 * t387;
t469 = t482 * t398 - t473 * t399 - t471 * t437;
t476 = t473 * t398 - t483 * t399 + t479 * t437;
t478 = -Ifges(6,3) - Ifges(7,2);
t481 = t479 * t362 + t476 * t398 + t471 * t361 + t478 * t427 - mrSges(6,1) * t347 + mrSges(7,1) * t346 + mrSges(6,2) * t348 - mrSges(7,3) * t345 + pkin(5) * t342 - qJ(6) * (-t361 * mrSges(7,2) - t398 * t378 + t463) + t469 * t399;
t474 = -mrSges(6,3) - mrSges(7,2);
t386 = mrSges(6,1) * t437 - mrSges(6,3) * t399;
t467 = -mrSges(6,1) * t398 - mrSges(6,2) * t399 - t378;
t336 = m(6) * t348 - t427 * mrSges(6,2) + t474 * t361 - t437 * t386 + t467 * t398 + t463;
t385 = -mrSges(6,2) * t437 - mrSges(6,3) * t398;
t338 = m(6) * t347 + t427 * mrSges(6,1) + t474 * t362 + t437 * t385 + t467 * t399 + t457;
t331 = t446 * t336 + t475 * t338;
t470 = t471 * t398 + t479 * t399 + t478 * t437;
t430 = (-mrSges(4,1) * t449 + mrSges(4,2) * t447) * qJD(1);
t434 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t466;
t402 = -mrSges(5,1) * t424 + mrSges(5,2) * t425;
t405 = mrSges(5,2) * t465 + mrSges(5,3) * t424;
t329 = m(5) * t355 - mrSges(5,1) * t433 - mrSges(5,3) * t408 - t402 * t425 - t405 * t465 + t331;
t406 = -mrSges(5,1) * t465 - mrSges(5,3) * t425;
t458 = t475 * t336 - t338 * t446;
t330 = m(5) * t356 + mrSges(5,2) * t433 + mrSges(5,3) * t407 + t402 * t424 + t406 * t465 + t458;
t459 = -t329 * t442 + t444 * t330;
t326 = m(4) * t384 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t433 - qJD(3) * t434 + t430 * t465 + t459;
t383 = -t447 * t392 + t449 * t441;
t376 = -qJDD(3) * pkin(3) - t451 * qJ(4) + t429 * t466 + qJDD(4) - t383;
t357 = -t407 * pkin(4) - t423 * pkin(8) + t425 * t409 + t376;
t350 = -0.2e1 * qJD(6) * t399 + (t398 * t437 - t362) * qJ(6) + (t399 * t437 + t361) * pkin(5) + t357;
t343 = m(7) * t350 + t361 * mrSges(7,1) - t362 * mrSges(7,3) - t399 * t387 + t398 * t388;
t455 = m(6) * t357 + t361 * mrSges(6,1) + t362 * mrSges(6,2) + t398 * t385 + t399 * t386 + t343;
t340 = m(5) * t376 - t407 * mrSges(5,1) + t408 * mrSges(5,2) - t424 * t405 + t425 * t406 + t455;
t435 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t465;
t339 = m(4) * t383 + qJDD(3) * mrSges(4,1) - t432 * mrSges(4,3) + qJD(3) * t435 - t430 * t466 - t340;
t460 = t449 * t326 - t339 * t447;
t327 = t329 * t444 + t330 * t442;
t453 = -m(4) * t391 + t433 * mrSges(4,1) - mrSges(4,2) * t432 - t434 * t466 + t435 * t465 - t327;
t419 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t447 + Ifges(4,4) * t449) * qJD(1);
t418 = Ifges(4,6) * qJD(3) + (t447 * Ifges(4,4) + Ifges(4,2) * t449) * qJD(1);
t395 = Ifges(5,1) * t425 + Ifges(5,4) * t424 - Ifges(5,5) * t465;
t394 = Ifges(5,4) * t425 + Ifges(5,2) * t424 - Ifges(5,6) * t465;
t393 = Ifges(5,5) * t425 + Ifges(5,6) * t424 - Ifges(5,3) * t465;
t333 = mrSges(6,2) * t357 + mrSges(7,2) * t346 - mrSges(6,3) * t347 - mrSges(7,3) * t350 - qJ(6) * t343 - t473 * t361 + t483 * t362 + t470 * t398 - t479 * t427 + t469 * t437;
t332 = -mrSges(6,1) * t357 - mrSges(7,1) * t350 + mrSges(7,2) * t345 + mrSges(6,3) * t348 - pkin(5) * t343 - t482 * t361 + t473 * t362 + t470 * t399 + t471 * t427 - t476 * t437;
t324 = mrSges(5,2) * t376 - mrSges(5,3) * t355 + Ifges(5,1) * t408 + Ifges(5,4) * t407 - Ifges(5,5) * t433 - pkin(8) * t331 - t446 * t332 + t475 * t333 + t424 * t393 + t394 * t465;
t323 = -mrSges(5,1) * t376 + mrSges(5,3) * t356 + Ifges(5,4) * t408 + Ifges(5,2) * t407 - Ifges(5,6) * t433 - pkin(4) * t455 + pkin(8) * t458 + t475 * t332 + t446 * t333 - t425 * t393 - t395 * t465;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t461 - mrSges(2,2) * t456 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t400 - mrSges(3,2) * t401 + t447 * (mrSges(4,2) * t391 - mrSges(4,3) * t383 + Ifges(4,1) * t432 + Ifges(4,4) * t433 + Ifges(4,5) * qJDD(3) - qJ(4) * t327 - qJD(3) * t418 - t442 * t323 + t444 * t324) + t449 * (Ifges(4,6) * qJDD(3) + (Ifges(4,2) + Ifges(5,3)) * t433 + Ifges(4,4) * t432 + t424 * t395 - t425 * t394 - Ifges(5,6) * t407 - Ifges(5,5) * t408 + qJD(3) * t419 - mrSges(4,1) * t391 + mrSges(4,3) * t384 - mrSges(5,1) * t355 + mrSges(5,2) * t356 - pkin(4) * t331 - pkin(3) * t327 + t481) + pkin(2) * t453 + pkin(7) * t460 + pkin(1) * (t443 * (m(3) * t401 - mrSges(3,1) * t452 - qJDD(1) * mrSges(3,2) + t460) + t445 * (m(3) * t400 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t452 + t453)); m(3) * t441 + t326 * t447 + t339 * t449; Ifges(4,5) * t432 + Ifges(4,6) * t433 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t383 - mrSges(4,2) * t384 + t442 * t324 + t444 * t323 - pkin(3) * t340 + qJ(4) * t459 + (t418 * t447 - t419 * t449) * qJD(1); t340; -t481; t342;];
tauJ  = t1;
