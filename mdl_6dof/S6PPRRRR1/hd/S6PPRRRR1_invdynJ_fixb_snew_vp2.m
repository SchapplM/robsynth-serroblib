% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-05-04 20:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PPRRRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:39:41
% EndTime: 2019-05-04 20:39:43
% DurationCPUTime: 2.13s
% Computational Cost: add. (21070->229), mult. (38871->304), div. (0->0), fcn. (30541->16), ass. (0->108)
t424 = sin(pkin(12));
t428 = cos(pkin(12));
t413 = -g(1) * t428 - g(2) * t424;
t423 = sin(pkin(13));
t427 = cos(pkin(13));
t412 = g(1) * t424 - g(2) * t428;
t422 = -g(3) + qJDD(1);
t426 = sin(pkin(6));
t430 = cos(pkin(6));
t447 = t412 * t430 + t422 * t426;
t386 = -t413 * t423 + t447 * t427;
t399 = -t412 * t426 + t422 * t430 + qJDD(2);
t425 = sin(pkin(7));
t429 = cos(pkin(7));
t459 = t386 * t429 + t399 * t425;
t432 = sin(qJ(5));
t433 = sin(qJ(4));
t436 = cos(qJ(5));
t437 = cos(qJ(4));
t404 = (t433 * t432 - t437 * t436) * qJD(3);
t387 = t413 * t427 + t447 * t423;
t434 = sin(qJ(3));
t438 = cos(qJ(3));
t361 = -t434 * t387 + t459 * t438;
t362 = t438 * t387 + t459 * t434;
t439 = qJD(3) ^ 2;
t360 = -pkin(3) * t439 + qJDD(3) * pkin(9) + t362;
t372 = -t386 * t425 + t399 * t429;
t355 = -t360 * t433 + t437 * t372;
t454 = qJD(3) * qJD(4);
t453 = t437 * t454;
t410 = qJDD(3) * t433 + t453;
t353 = (-t410 + t453) * pkin(10) + (t433 * t437 * t439 + qJDD(4)) * pkin(4) + t355;
t356 = t437 * t360 + t433 * t372;
t411 = qJDD(3) * t437 - t433 * t454;
t455 = t433 * qJD(3);
t416 = qJD(4) * pkin(4) - pkin(10) * t455;
t421 = t437 ^ 2;
t354 = -pkin(4) * t421 * t439 + pkin(10) * t411 - qJD(4) * t416 + t356;
t349 = t432 * t353 + t436 * t354;
t405 = (t437 * t432 + t433 * t436) * qJD(3);
t379 = -qJD(5) * t405 - t410 * t432 + t411 * t436;
t392 = mrSges(6,1) * t404 + mrSges(6,2) * t405;
t420 = qJD(4) + qJD(5);
t398 = mrSges(6,1) * t420 - mrSges(6,3) * t405;
t419 = qJDD(4) + qJDD(5);
t393 = pkin(5) * t404 - pkin(11) * t405;
t418 = t420 ^ 2;
t346 = -pkin(5) * t418 + pkin(11) * t419 - t393 * t404 + t349;
t444 = -qJDD(3) * pkin(3) - t361;
t357 = -pkin(4) * t411 + t416 * t455 + (-pkin(10) * t421 - pkin(9)) * t439 + t444;
t380 = -qJD(5) * t404 + t410 * t436 + t411 * t432;
t350 = (t404 * t420 - t380) * pkin(11) + (t405 * t420 - t379) * pkin(5) + t357;
t431 = sin(qJ(6));
t435 = cos(qJ(6));
t343 = -t346 * t431 + t350 * t435;
t395 = -t405 * t431 + t420 * t435;
t365 = qJD(6) * t395 + t380 * t435 + t419 * t431;
t396 = t405 * t435 + t420 * t431;
t373 = -mrSges(7,1) * t395 + mrSges(7,2) * t396;
t378 = qJDD(6) - t379;
t400 = qJD(6) + t404;
t381 = -mrSges(7,2) * t400 + mrSges(7,3) * t395;
t340 = m(7) * t343 + mrSges(7,1) * t378 - mrSges(7,3) * t365 - t373 * t396 + t381 * t400;
t344 = t346 * t435 + t350 * t431;
t364 = -qJD(6) * t396 - t380 * t431 + t419 * t435;
t382 = mrSges(7,1) * t400 - mrSges(7,3) * t396;
t341 = m(7) * t344 - mrSges(7,2) * t378 + mrSges(7,3) * t364 + t373 * t395 - t382 * t400;
t450 = -t340 * t431 + t435 * t341;
t328 = m(6) * t349 - mrSges(6,2) * t419 + mrSges(6,3) * t379 - t392 * t404 - t398 * t420 + t450;
t348 = t353 * t436 - t354 * t432;
t397 = -mrSges(6,2) * t420 - mrSges(6,3) * t404;
t345 = -pkin(5) * t419 - pkin(11) * t418 + t393 * t405 - t348;
t445 = -m(7) * t345 + t364 * mrSges(7,1) - mrSges(7,2) * t365 + t395 * t381 - t382 * t396;
t336 = m(6) * t348 + mrSges(6,1) * t419 - mrSges(6,3) * t380 - t392 * t405 + t397 * t420 + t445;
t324 = t432 * t328 + t436 * t336;
t330 = t435 * t340 + t431 * t341;
t456 = qJD(3) * t437;
t409 = (-t437 * mrSges(5,1) + t433 * mrSges(5,2)) * qJD(3);
t415 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t456;
t322 = m(5) * t355 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t410 + qJD(4) * t415 - t409 * t455 + t324;
t414 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t455;
t451 = t436 * t328 - t336 * t432;
t323 = m(5) * t356 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t411 - qJD(4) * t414 + t409 * t456 + t451;
t452 = -t322 * t433 + t437 * t323;
t317 = m(4) * t362 - mrSges(4,1) * t439 - qJDD(3) * mrSges(4,2) + t452;
t359 = -pkin(9) * t439 + t444;
t443 = m(6) * t357 - t379 * mrSges(6,1) + mrSges(6,2) * t380 + t404 * t397 + t398 * t405 + t330;
t440 = -m(5) * t359 + t411 * mrSges(5,1) - mrSges(5,2) * t410 - t414 * t455 + t415 * t456 - t443;
t325 = m(4) * t361 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t439 + t440;
t449 = t317 * t434 + t325 * t438;
t366 = Ifges(7,5) * t396 + Ifges(7,6) * t395 + Ifges(7,3) * t400;
t368 = Ifges(7,1) * t396 + Ifges(7,4) * t395 + Ifges(7,5) * t400;
t333 = -mrSges(7,1) * t345 + mrSges(7,3) * t344 + Ifges(7,4) * t365 + Ifges(7,2) * t364 + Ifges(7,6) * t378 - t366 * t396 + t368 * t400;
t367 = Ifges(7,4) * t396 + Ifges(7,2) * t395 + Ifges(7,6) * t400;
t334 = mrSges(7,2) * t345 - mrSges(7,3) * t343 + Ifges(7,1) * t365 + Ifges(7,4) * t364 + Ifges(7,5) * t378 + t366 * t395 - t367 * t400;
t389 = Ifges(6,4) * t405 - Ifges(6,2) * t404 + Ifges(6,6) * t420;
t390 = Ifges(6,1) * t405 - Ifges(6,4) * t404 + Ifges(6,5) * t420;
t442 = mrSges(6,1) * t348 - mrSges(6,2) * t349 + Ifges(6,5) * t380 + Ifges(6,6) * t379 + Ifges(6,3) * t419 + pkin(5) * t445 + pkin(11) * t450 + t435 * t333 + t431 * t334 + t405 * t389 + t390 * t404;
t441 = mrSges(7,1) * t343 - mrSges(7,2) * t344 + Ifges(7,5) * t365 + Ifges(7,6) * t364 + Ifges(7,3) * t378 + t367 * t396 - t395 * t368;
t403 = Ifges(5,5) * qJD(4) + (t433 * Ifges(5,1) + t437 * Ifges(5,4)) * qJD(3);
t402 = Ifges(5,6) * qJD(4) + (t433 * Ifges(5,4) + t437 * Ifges(5,2)) * qJD(3);
t388 = Ifges(6,5) * t405 - Ifges(6,6) * t404 + Ifges(6,3) * t420;
t320 = -mrSges(6,1) * t357 + mrSges(6,3) * t349 + Ifges(6,4) * t380 + Ifges(6,2) * t379 + Ifges(6,6) * t419 - pkin(5) * t330 - t388 * t405 + t390 * t420 - t441;
t319 = mrSges(6,2) * t357 - mrSges(6,3) * t348 + Ifges(6,1) * t380 + Ifges(6,4) * t379 + Ifges(6,5) * t419 - pkin(11) * t330 - t333 * t431 + t334 * t435 - t388 * t404 - t389 * t420;
t318 = m(4) * t372 + t322 * t437 + t323 * t433;
t316 = m(3) * t399 + t318 * t429 + t449 * t425;
t1 = [m(2) * t422 + t430 * t316 + (t423 * (m(3) * t387 + t317 * t438 - t325 * t434) + t427 * (m(3) * t386 - t318 * t425 + t449 * t429)) * t426; t316; Ifges(4,3) * qJDD(3) + mrSges(4,1) * t361 - mrSges(4,2) * t362 + t433 * (mrSges(5,2) * t359 - mrSges(5,3) * t355 + Ifges(5,1) * t410 + Ifges(5,4) * t411 + Ifges(5,5) * qJDD(4) - pkin(10) * t324 - qJD(4) * t402 + t319 * t436 - t320 * t432) + t437 * (-mrSges(5,1) * t359 + mrSges(5,3) * t356 + Ifges(5,4) * t410 + Ifges(5,2) * t411 + Ifges(5,6) * qJDD(4) - pkin(4) * t443 + pkin(10) * t451 + qJD(4) * t403 + t432 * t319 + t436 * t320) + pkin(3) * t440 + pkin(9) * t452; mrSges(5,1) * t355 - mrSges(5,2) * t356 + Ifges(5,3) * qJDD(4) + t442 + Ifges(5,5) * t410 + Ifges(5,6) * t411 + pkin(4) * t324 + (t433 * t402 - t437 * t403) * qJD(3); t442; t441;];
tauJ  = t1;
