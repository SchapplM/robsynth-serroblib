% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:35
% EndTime: 2019-12-05 17:44:38
% DurationCPUTime: 2.38s
% Computational Cost: add. (13893->231), mult. (38009->320), div. (0->0), fcn. (26039->10), ass. (0->112)
t433 = qJD(1) ^ 2;
t429 = sin(qJ(1));
t432 = cos(qJ(1));
t454 = t429 * g(2) - g(3) * t432;
t470 = -pkin(1) * t433 + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t454;
t450 = g(2) * t432 + g(3) * t429;
t438 = -qJ(2) * t433 + qJDD(2) - t450;
t424 = sin(pkin(8));
t426 = cos(pkin(8));
t447 = -t426 * pkin(2) - t424 * qJ(3);
t462 = qJD(1) * t424;
t469 = (-pkin(1) + t447) * qJDD(1) + t438 - 0.2e1 * qJD(3) * t462;
t387 = -t426 * g(1) - t470 * t424;
t425 = cos(pkin(9));
t468 = Ifges(4,4) * t425;
t423 = sin(pkin(9));
t467 = Ifges(4,2) * t423;
t466 = t423 * t424;
t465 = t424 * t425;
t388 = -g(1) * t424 + t470 * t426;
t406 = t447 * qJD(1);
t461 = qJD(1) * t426;
t377 = t406 * t461 + t388;
t422 = t424 ^ 2;
t445 = -t426 * pkin(3) - pkin(6) * t465;
t464 = t469 * t425;
t356 = t445 * qJDD(1) + (-t377 + (-pkin(3) * t422 * t425 + pkin(6) * t424 * t426) * t433) * t423 + t464;
t364 = t425 * t377 + t469 * t423;
t405 = t445 * qJD(1);
t459 = qJDD(1) * t424;
t455 = t423 * t459;
t457 = t423 ^ 2 * t422 * t433;
t357 = -pkin(3) * t457 - pkin(6) * t455 + t405 * t461 + t364;
t428 = sin(qJ(4));
t431 = cos(qJ(4));
t343 = t431 * t356 - t357 * t428;
t441 = t424 * (-t423 * t431 - t425 * t428);
t395 = qJD(1) * t441;
t440 = t424 * (-t423 * t428 + t425 * t431);
t381 = qJD(4) * t395 + qJDD(1) * t440;
t396 = qJD(1) * t440;
t458 = qJDD(1) * t426;
t415 = qJDD(4) - t458;
t416 = qJD(4) - t461;
t341 = (t395 * t416 - t381) * pkin(7) + (t395 * t396 + t415) * pkin(4) + t343;
t344 = t428 * t356 + t431 * t357;
t380 = -qJD(4) * t396 + qJDD(1) * t441;
t386 = pkin(4) * t416 - pkin(7) * t396;
t394 = t395 ^ 2;
t342 = -pkin(4) * t394 + pkin(7) * t380 - t386 * t416 + t344;
t427 = sin(qJ(5));
t430 = cos(qJ(5));
t339 = t341 * t430 - t342 * t427;
t374 = t395 * t430 - t396 * t427;
t352 = qJD(5) * t374 + t380 * t427 + t381 * t430;
t375 = t395 * t427 + t396 * t430;
t362 = -mrSges(6,1) * t374 + mrSges(6,2) * t375;
t414 = qJD(5) + t416;
t367 = -mrSges(6,2) * t414 + mrSges(6,3) * t374;
t412 = qJDD(5) + t415;
t336 = m(6) * t339 + mrSges(6,1) * t412 - mrSges(6,3) * t352 - t362 * t375 + t367 * t414;
t340 = t341 * t427 + t342 * t430;
t351 = -qJD(5) * t375 + t380 * t430 - t381 * t427;
t368 = mrSges(6,1) * t414 - mrSges(6,3) * t375;
t337 = m(6) * t340 - mrSges(6,2) * t412 + mrSges(6,3) * t351 + t362 * t374 - t368 * t414;
t330 = t430 * t336 + t427 * t337;
t378 = -mrSges(5,1) * t395 + mrSges(5,2) * t396;
t382 = -mrSges(5,2) * t416 + mrSges(5,3) * t395;
t328 = m(5) * t343 + mrSges(5,1) * t415 - mrSges(5,3) * t381 - t378 * t396 + t382 * t416 + t330;
t383 = mrSges(5,1) * t416 - mrSges(5,3) * t396;
t452 = -t336 * t427 + t430 * t337;
t329 = m(5) * t344 - mrSges(5,2) * t415 + mrSges(5,3) * t380 + t378 * t395 - t383 * t416 + t452;
t324 = t431 * t328 + t428 * t329;
t453 = -t328 * t428 + t431 * t329;
t449 = -t426 * mrSges(3,1) + t424 * mrSges(3,2);
t448 = t423 * mrSges(4,1) + t425 * mrSges(4,2);
t363 = -t377 * t423 + t464;
t399 = t448 * t462;
t442 = t426 * mrSges(4,2) - mrSges(4,3) * t466;
t402 = t442 * qJD(1);
t443 = -t426 * mrSges(4,1) - mrSges(4,3) * t465;
t322 = m(4) * t363 + t443 * qJDD(1) + (-t399 * t465 - t402 * t426) * qJD(1) + t324;
t403 = t443 * qJD(1);
t323 = m(4) * t364 + t442 * qJDD(1) + (-t399 * t466 + t403 * t426) * qJD(1) + t453;
t319 = t322 * t425 + t323 * t423;
t446 = t402 * t423 + t403 * t425;
t376 = t406 * t462 + qJDD(3) - t387;
t444 = -Ifges(4,5) * t425 + Ifges(4,6) * t423 + Ifges(3,4);
t366 = t425 * t405 * t462 + pkin(3) * t455 - pkin(6) * t457 + t376;
t346 = -pkin(4) * t380 - pkin(7) * t394 + t386 * t396 + t366;
t439 = m(6) * t346 - t351 * mrSges(6,1) + t352 * mrSges(6,2) - t374 * t367 + t375 * t368;
t359 = Ifges(6,4) * t375 + Ifges(6,2) * t374 + Ifges(6,6) * t414;
t360 = Ifges(6,1) * t375 + Ifges(6,4) * t374 + Ifges(6,5) * t414;
t437 = -mrSges(6,1) * t339 + mrSges(6,2) * t340 - Ifges(6,5) * t352 - Ifges(6,6) * t351 - Ifges(6,3) * t412 - t375 * t359 + t374 * t360;
t436 = m(5) * t366 - t380 * mrSges(5,1) + t381 * mrSges(5,2) - t395 * t382 + t396 * t383 + t439;
t435 = m(4) * t376 + t436;
t370 = Ifges(5,4) * t396 + Ifges(5,2) * t395 + Ifges(5,6) * t416;
t371 = Ifges(5,1) * t396 + Ifges(5,4) * t395 + Ifges(5,5) * t416;
t434 = mrSges(5,1) * t343 - mrSges(5,2) * t344 + Ifges(5,5) * t381 + Ifges(5,6) * t380 + Ifges(5,3) * t415 + pkin(4) * t330 + t396 * t370 - t395 * t371 - t437;
t408 = (Ifges(3,5) * t424 + Ifges(3,6) * t426) * qJD(1);
t407 = t449 * qJD(1);
t401 = -qJDD(1) * pkin(1) + t438;
t392 = (-Ifges(4,5) * t426 + (Ifges(4,1) * t425 - Ifges(4,4) * t423) * t424) * qJD(1);
t391 = (-Ifges(4,6) * t426 + (-t467 + t468) * t424) * qJD(1);
t369 = Ifges(5,5) * t396 + Ifges(5,6) * t395 + Ifges(5,3) * t416;
t358 = Ifges(6,5) * t375 + Ifges(6,6) * t374 + Ifges(6,3) * t414;
t332 = mrSges(6,2) * t346 - mrSges(6,3) * t339 + Ifges(6,1) * t352 + Ifges(6,4) * t351 + Ifges(6,5) * t412 + t358 * t374 - t359 * t414;
t331 = -mrSges(6,1) * t346 + mrSges(6,3) * t340 + Ifges(6,4) * t352 + Ifges(6,2) * t351 + Ifges(6,6) * t412 - t358 * t375 + t360 * t414;
t321 = mrSges(5,2) * t366 - mrSges(5,3) * t343 + Ifges(5,1) * t381 + Ifges(5,4) * t380 + Ifges(5,5) * t415 - pkin(7) * t330 - t331 * t427 + t332 * t430 + t369 * t395 - t370 * t416;
t320 = -mrSges(5,1) * t366 + mrSges(5,3) * t344 + Ifges(5,4) * t381 + Ifges(5,2) * t380 + Ifges(5,6) * t415 - pkin(4) * t439 + pkin(7) * t452 + t430 * t331 + t427 * t332 - t396 * t369 + t416 * t371;
t318 = m(3) * t401 + t449 * qJDD(1) + (-t426 ^ 2 - t422) * t433 * mrSges(3,3) + t319;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t450 - mrSges(2,2) * t454 + t424 * (t408 * t461 + mrSges(3,2) * t401 - mrSges(3,3) * t387 + t425 * (mrSges(4,2) * t376 - mrSges(4,3) * t363 - pkin(6) * t324 - t320 * t428 + t431 * t321 + t391 * t461) - t423 * (-mrSges(4,1) * t376 + mrSges(4,3) * t364 - pkin(3) * t436 + pkin(6) * t453 + t431 * t320 + t428 * t321 - t392 * t461) - qJ(3) * t319 + (t426 * t444 + (Ifges(4,1) * t425 ^ 2 + Ifges(3,1) + (t467 - 0.2e1 * t468) * t423) * t424) * qJDD(1)) + t426 * (-t434 + (t444 * qJDD(1) + (-t391 * t425 - t392 * t423 - t408) * qJD(1)) * t424 + (Ifges(3,2) + Ifges(4,3)) * t458 - mrSges(3,1) * t401 + mrSges(3,3) * t388 - mrSges(4,1) * t363 + mrSges(4,2) * t364 - pkin(3) * t324 - pkin(2) * t319) - pkin(1) * t318 + qJ(2) * ((m(3) * t388 - t423 * t322 + t425 * t323 + (qJDD(1) * mrSges(3,3) + qJD(1) * t407) * t426) * t426 + (t435 - m(3) * t387 + (t407 + t446) * t462 + (mrSges(3,3) + t448) * t459) * t424); t318; (qJD(1) * t446 + t448 * qJDD(1)) * t424 + t435; t434; -t437;];
tauJ = t1;
