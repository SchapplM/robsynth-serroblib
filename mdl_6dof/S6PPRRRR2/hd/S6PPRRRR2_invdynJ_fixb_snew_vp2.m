% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PPRRRR2
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
% Datum: 2019-05-04 20:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PPRRRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:55:21
% EndTime: 2019-05-04 20:55:23
% DurationCPUTime: 2.20s
% Computational Cost: add. (22778->229), mult. (41533->299), div. (0->0), fcn. (32543->16), ass. (0->108)
t408 = sin(pkin(12));
t412 = cos(pkin(12));
t400 = -t412 * g(1) - t408 * g(2);
t407 = sin(pkin(13));
t411 = cos(pkin(13));
t399 = t408 * g(1) - t412 * g(2);
t406 = -g(3) + qJDD(1);
t410 = sin(pkin(6));
t414 = cos(pkin(6));
t430 = t399 * t414 + t406 * t410;
t366 = -t407 * t400 + t430 * t411;
t382 = -t410 * t399 + t414 * t406 + qJDD(2);
t409 = sin(pkin(7));
t413 = cos(pkin(7));
t442 = t366 * t413 + t382 * t409;
t367 = t411 * t400 + t430 * t407;
t418 = sin(qJ(3));
t422 = cos(qJ(3));
t346 = -t418 * t367 + t422 * t442;
t347 = t422 * t367 + t418 * t442;
t424 = qJD(3) ^ 2;
t345 = -t424 * pkin(3) + qJDD(3) * pkin(9) + t347;
t359 = -t409 * t366 + t413 * t382;
t417 = sin(qJ(4));
t421 = cos(qJ(4));
t338 = t421 * t345 + t417 * t359;
t396 = (-t421 * pkin(4) - t417 * pkin(10)) * qJD(3);
t423 = qJD(4) ^ 2;
t438 = t421 * qJD(3);
t336 = -t423 * pkin(4) + qJDD(4) * pkin(10) + t396 * t438 + t338;
t344 = -qJDD(3) * pkin(3) - t424 * pkin(9) - t346;
t437 = qJD(3) * qJD(4);
t436 = t421 * t437;
t397 = t417 * qJDD(3) + t436;
t405 = t417 * t437;
t398 = t421 * qJDD(3) - t405;
t341 = (-t397 - t436) * pkin(10) + (-t398 + t405) * pkin(4) + t344;
t416 = sin(qJ(5));
t420 = cos(qJ(5));
t331 = -t416 * t336 + t420 * t341;
t439 = t417 * qJD(3);
t393 = t420 * qJD(4) - t416 * t439;
t374 = t393 * qJD(5) + t416 * qJDD(4) + t420 * t397;
t392 = qJDD(5) - t398;
t394 = t416 * qJD(4) + t420 * t439;
t404 = qJD(5) - t438;
t329 = (t393 * t404 - t374) * pkin(11) + (t393 * t394 + t392) * pkin(5) + t331;
t332 = t420 * t336 + t416 * t341;
t373 = -t394 * qJD(5) + t420 * qJDD(4) - t416 * t397;
t381 = t404 * pkin(5) - t394 * pkin(11);
t391 = t393 ^ 2;
t330 = -t391 * pkin(5) + t373 * pkin(11) - t404 * t381 + t332;
t415 = sin(qJ(6));
t419 = cos(qJ(6));
t327 = t419 * t329 - t415 * t330;
t375 = t419 * t393 - t415 * t394;
t353 = t375 * qJD(6) + t415 * t373 + t419 * t374;
t376 = t415 * t393 + t419 * t394;
t360 = -t375 * mrSges(7,1) + t376 * mrSges(7,2);
t403 = qJD(6) + t404;
t362 = -t403 * mrSges(7,2) + t375 * mrSges(7,3);
t388 = qJDD(6) + t392;
t324 = m(7) * t327 + t388 * mrSges(7,1) - t353 * mrSges(7,3) - t376 * t360 + t403 * t362;
t328 = t415 * t329 + t419 * t330;
t352 = -t376 * qJD(6) + t419 * t373 - t415 * t374;
t363 = t403 * mrSges(7,1) - t376 * mrSges(7,3);
t325 = m(7) * t328 - t388 * mrSges(7,2) + t352 * mrSges(7,3) + t375 * t360 - t403 * t363;
t317 = t419 * t324 + t415 * t325;
t395 = (-t421 * mrSges(5,1) + t417 * mrSges(5,2)) * qJD(3);
t401 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t439;
t377 = -t393 * mrSges(6,1) + t394 * mrSges(6,2);
t379 = -t404 * mrSges(6,2) + t393 * mrSges(6,3);
t315 = m(6) * t331 + t392 * mrSges(6,1) - t374 * mrSges(6,3) - t394 * t377 + t404 * t379 + t317;
t380 = t404 * mrSges(6,1) - t394 * mrSges(6,3);
t433 = -t415 * t324 + t419 * t325;
t316 = m(6) * t332 - t392 * mrSges(6,2) + t373 * mrSges(6,3) + t393 * t377 - t404 * t380 + t433;
t434 = -t416 * t315 + t420 * t316;
t312 = m(5) * t338 - qJDD(4) * mrSges(5,2) + t398 * mrSges(5,3) - qJD(4) * t401 + t395 * t438 + t434;
t337 = -t417 * t345 + t421 * t359;
t402 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t438;
t335 = -qJDD(4) * pkin(4) - t423 * pkin(10) + t396 * t439 - t337;
t333 = -t373 * pkin(5) - t391 * pkin(11) + t394 * t381 + t335;
t429 = m(7) * t333 - t352 * mrSges(7,1) + t353 * mrSges(7,2) - t375 * t362 + t376 * t363;
t426 = -m(6) * t335 + t373 * mrSges(6,1) - t374 * mrSges(6,2) + t393 * t379 - t394 * t380 - t429;
t320 = m(5) * t337 + qJDD(4) * mrSges(5,1) - t397 * mrSges(5,3) + qJD(4) * t402 - t395 * t439 + t426;
t435 = t421 * t312 - t417 * t320;
t306 = m(4) * t347 - t424 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t435;
t313 = t420 * t315 + t416 * t316;
t427 = -m(5) * t344 + t398 * mrSges(5,1) - t397 * mrSges(5,2) - t401 * t439 + t402 * t438 - t313;
t310 = m(4) * t346 + qJDD(3) * mrSges(4,1) - t424 * mrSges(4,2) + t427;
t432 = t306 * t418 + t310 * t422;
t355 = Ifges(7,4) * t376 + Ifges(7,2) * t375 + Ifges(7,6) * t403;
t356 = Ifges(7,1) * t376 + Ifges(7,4) * t375 + Ifges(7,5) * t403;
t428 = -mrSges(7,1) * t327 + mrSges(7,2) * t328 - Ifges(7,5) * t353 - Ifges(7,6) * t352 - Ifges(7,3) * t388 - t376 * t355 + t375 * t356;
t369 = Ifges(6,4) * t394 + Ifges(6,2) * t393 + Ifges(6,6) * t404;
t370 = Ifges(6,1) * t394 + Ifges(6,4) * t393 + Ifges(6,5) * t404;
t425 = mrSges(6,1) * t331 - mrSges(6,2) * t332 + Ifges(6,5) * t374 + Ifges(6,6) * t373 + Ifges(6,3) * t392 + pkin(5) * t317 + t394 * t369 - t393 * t370 - t428;
t387 = Ifges(5,5) * qJD(4) + (t417 * Ifges(5,1) + t421 * Ifges(5,4)) * qJD(3);
t386 = Ifges(5,6) * qJD(4) + (t417 * Ifges(5,4) + t421 * Ifges(5,2)) * qJD(3);
t368 = Ifges(6,5) * t394 + Ifges(6,6) * t393 + Ifges(6,3) * t404;
t354 = Ifges(7,5) * t376 + Ifges(7,6) * t375 + Ifges(7,3) * t403;
t319 = mrSges(7,2) * t333 - mrSges(7,3) * t327 + Ifges(7,1) * t353 + Ifges(7,4) * t352 + Ifges(7,5) * t388 + t375 * t354 - t403 * t355;
t318 = -mrSges(7,1) * t333 + mrSges(7,3) * t328 + Ifges(7,4) * t353 + Ifges(7,2) * t352 + Ifges(7,6) * t388 - t376 * t354 + t403 * t356;
t309 = mrSges(6,2) * t335 - mrSges(6,3) * t331 + Ifges(6,1) * t374 + Ifges(6,4) * t373 + Ifges(6,5) * t392 - pkin(11) * t317 - t415 * t318 + t419 * t319 + t393 * t368 - t404 * t369;
t308 = -mrSges(6,1) * t335 + mrSges(6,3) * t332 + Ifges(6,4) * t374 + Ifges(6,2) * t373 + Ifges(6,6) * t392 - pkin(5) * t429 + pkin(11) * t433 + t419 * t318 + t415 * t319 - t394 * t368 + t404 * t370;
t307 = m(4) * t359 + t417 * t312 + t421 * t320;
t305 = m(3) * t382 + t413 * t307 + t432 * t409;
t1 = [m(2) * t406 + t414 * t305 + (t407 * (m(3) * t367 + t422 * t306 - t418 * t310) + t411 * (m(3) * t366 - t409 * t307 + t432 * t413)) * t410; t305; Ifges(4,3) * qJDD(3) + mrSges(4,1) * t346 - mrSges(4,2) * t347 + t417 * (mrSges(5,2) * t344 - mrSges(5,3) * t337 + Ifges(5,1) * t397 + Ifges(5,4) * t398 + Ifges(5,5) * qJDD(4) - pkin(10) * t313 - qJD(4) * t386 - t416 * t308 + t420 * t309) + t421 * (-mrSges(5,1) * t344 + mrSges(5,3) * t338 + Ifges(5,4) * t397 + Ifges(5,2) * t398 + Ifges(5,6) * qJDD(4) - pkin(4) * t313 + qJD(4) * t387 - t425) + pkin(3) * t427 + pkin(9) * t435; Ifges(5,5) * t397 + Ifges(5,6) * t398 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t337 - mrSges(5,2) * t338 + t416 * t309 + t420 * t308 + pkin(4) * t426 + pkin(10) * t434 + (t417 * t386 - t421 * t387) * qJD(3); t425; -t428;];
tauJ  = t1;
