% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 14:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:09:12
% EndTime: 2019-05-05 14:09:13
% DurationCPUTime: 1.32s
% Computational Cost: add. (9163->237), mult. (18638->298), div. (0->0), fcn. (11266->10), ass. (0->99)
t410 = sin(qJ(1));
t413 = cos(qJ(1));
t428 = t410 * g(1) - t413 * g(2);
t386 = qJDD(1) * pkin(1) + t428;
t415 = qJD(1) ^ 2;
t424 = -t413 * g(1) - t410 * g(2);
t388 = -t415 * pkin(1) + t424;
t405 = sin(pkin(9));
t407 = cos(pkin(9));
t365 = t405 * t386 + t407 * t388;
t425 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t365;
t404 = sin(pkin(10));
t406 = cos(pkin(10));
t409 = sin(qJ(4));
t412 = cos(qJ(4));
t374 = (t412 * t404 + t409 * t406) * qJD(1);
t436 = 2 * qJD(5);
t435 = -pkin(2) - pkin(7);
t434 = pkin(4) * t415;
t364 = t407 * t386 - t405 * t388;
t420 = -t415 * qJ(3) + qJDD(3) - t364;
t353 = t435 * qJDD(1) + t420;
t348 = t412 * t353;
t431 = qJD(1) * qJD(4);
t390 = t412 * qJDD(1) - t409 * t431;
t401 = -g(3) + qJDD(2);
t333 = qJDD(4) * pkin(4) - t390 * qJ(5) + t348 + (-qJ(5) * t431 - t412 * t434 - t401) * t409;
t345 = t409 * t353 + t412 * t401;
t389 = -t409 * qJDD(1) - t412 * t431;
t432 = qJD(1) * t412;
t392 = qJD(4) * pkin(4) - qJ(5) * t432;
t400 = t409 ^ 2;
t334 = t389 * qJ(5) - qJD(4) * t392 - t400 * t434 + t345;
t329 = t404 * t333 + t406 * t334 - t374 * t436;
t433 = qJD(1) * t409;
t375 = -t404 * t433 + t406 * t432;
t360 = t374 * mrSges(6,1) + t375 * mrSges(6,2);
t366 = t406 * t389 - t404 * t390;
t371 = qJD(4) * mrSges(6,1) - t375 * mrSges(6,3);
t361 = t374 * pkin(5) - t375 * pkin(8);
t414 = qJD(4) ^ 2;
t327 = -t414 * pkin(5) + qJDD(4) * pkin(8) - t374 * t361 + t329;
t336 = -t389 * pkin(4) + qJDD(5) + t392 * t432 + (-qJ(5) * t400 + t435) * t415 + t425;
t367 = t404 * t389 + t406 * t390;
t330 = (qJD(4) * t374 - t367) * pkin(8) + (qJD(4) * t375 - t366) * pkin(5) + t336;
t408 = sin(qJ(6));
t411 = cos(qJ(6));
t324 = -t408 * t327 + t411 * t330;
t368 = t411 * qJD(4) - t408 * t375;
t343 = t368 * qJD(6) + t408 * qJDD(4) + t411 * t367;
t369 = t408 * qJD(4) + t411 * t375;
t346 = -t368 * mrSges(7,1) + t369 * mrSges(7,2);
t373 = qJD(6) + t374;
t350 = -t373 * mrSges(7,2) + t368 * mrSges(7,3);
t363 = qJDD(6) - t366;
t322 = m(7) * t324 + t363 * mrSges(7,1) - t343 * mrSges(7,3) - t369 * t346 + t373 * t350;
t325 = t411 * t327 + t408 * t330;
t342 = -t369 * qJD(6) + t411 * qJDD(4) - t408 * t367;
t351 = t373 * mrSges(7,1) - t369 * mrSges(7,3);
t323 = m(7) * t325 - t363 * mrSges(7,2) + t342 * mrSges(7,3) + t368 * t346 - t373 * t351;
t426 = -t408 * t322 + t411 * t323;
t312 = m(6) * t329 - qJDD(4) * mrSges(6,2) + t366 * mrSges(6,3) - qJD(4) * t371 - t374 * t360 + t426;
t422 = -t406 * t333 + t404 * t334;
t328 = -0.2e1 * qJD(5) * t375 - t422;
t370 = -qJD(4) * mrSges(6,2) - t374 * mrSges(6,3);
t326 = -qJDD(4) * pkin(5) - t414 * pkin(8) + (t436 + t361) * t375 + t422;
t418 = -m(7) * t326 + t342 * mrSges(7,1) - t343 * mrSges(7,2) + t368 * t350 - t369 * t351;
t318 = m(6) * t328 + qJDD(4) * mrSges(6,1) - t367 * mrSges(6,3) + qJD(4) * t370 - t375 * t360 + t418;
t309 = t404 * t312 + t406 * t318;
t314 = t411 * t322 + t408 * t323;
t427 = t406 * t312 - t404 * t318;
t344 = -t409 * t401 + t348;
t387 = (t409 * mrSges(5,1) + t412 * mrSges(5,2)) * qJD(1);
t391 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t433;
t307 = m(5) * t344 + qJDD(4) * mrSges(5,1) - t390 * mrSges(5,3) + qJD(4) * t391 - t387 * t432 + t309;
t393 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t432;
t308 = m(5) * t345 - qJDD(4) * mrSges(5,2) + t389 * mrSges(5,3) - qJD(4) * t393 - t387 * t433 + t427;
t423 = t412 * t307 + t409 * t308;
t355 = -qJDD(1) * pkin(2) + t420;
t419 = m(4) * t355 - t415 * mrSges(4,3) + t423;
t313 = m(6) * t336 - t366 * mrSges(6,1) + t367 * mrSges(6,2) + t374 * t370 + t375 * t371 + t314;
t338 = Ifges(7,4) * t369 + Ifges(7,2) * t368 + Ifges(7,6) * t373;
t339 = Ifges(7,1) * t369 + Ifges(7,4) * t368 + Ifges(7,5) * t373;
t417 = mrSges(7,1) * t324 - mrSges(7,2) * t325 + Ifges(7,5) * t343 + Ifges(7,6) * t342 + Ifges(7,3) * t363 + t369 * t338 - t368 * t339;
t352 = t435 * t415 + t425;
t354 = t415 * pkin(2) - t425;
t416 = -m(4) * t354 + m(5) * t352 - t389 * mrSges(5,1) + t415 * mrSges(4,2) + t390 * mrSges(5,2) + qJDD(1) * mrSges(4,3) + t391 * t433 + t393 * t432 + t313;
t381 = (Ifges(5,5) * qJD(4)) + (t412 * Ifges(5,1) - t409 * Ifges(5,4)) * qJD(1);
t380 = (Ifges(5,6) * qJD(4)) + (t412 * Ifges(5,4) - t409 * Ifges(5,2)) * qJD(1);
t358 = Ifges(6,1) * t375 - Ifges(6,4) * t374 + (Ifges(6,5) * qJD(4));
t357 = Ifges(6,4) * t375 - Ifges(6,2) * t374 + (Ifges(6,6) * qJD(4));
t356 = Ifges(6,5) * t375 - Ifges(6,6) * t374 + (Ifges(6,3) * qJD(4));
t337 = Ifges(7,5) * t369 + Ifges(7,6) * t368 + Ifges(7,3) * t373;
t316 = mrSges(7,2) * t326 - mrSges(7,3) * t324 + Ifges(7,1) * t343 + Ifges(7,4) * t342 + Ifges(7,5) * t363 + t368 * t337 - t373 * t338;
t315 = -mrSges(7,1) * t326 + mrSges(7,3) * t325 + Ifges(7,4) * t343 + Ifges(7,2) * t342 + Ifges(7,6) * t363 - t369 * t337 + t373 * t339;
t306 = -mrSges(6,1) * t336 + mrSges(6,3) * t329 + Ifges(6,4) * t367 + Ifges(6,2) * t366 + Ifges(6,6) * qJDD(4) - pkin(5) * t314 + qJD(4) * t358 - t375 * t356 - t417;
t305 = mrSges(6,2) * t336 - mrSges(6,3) * t328 + Ifges(6,1) * t367 + Ifges(6,4) * t366 + Ifges(6,5) * qJDD(4) - pkin(8) * t314 - qJD(4) * t357 - t408 * t315 + t411 * t316 - t374 * t356;
t304 = qJDD(1) * mrSges(4,2) + t419;
t1 = [pkin(1) * (t405 * (m(3) * t365 - t415 * mrSges(3,1) + t416) + t407 * (m(3) * t364 - t415 * mrSges(3,2) - t419)) + mrSges(2,1) * t428 - mrSges(2,2) * t424 - pkin(2) * t304 + qJ(3) * t416 - pkin(7) * t423 + mrSges(3,1) * t364 - mrSges(3,2) * t365 + t412 * (mrSges(5,2) * t352 - mrSges(5,3) * t344 + Ifges(5,1) * t390 + Ifges(5,4) * t389 + Ifges(5,5) * qJDD(4) - qJ(5) * t309 - qJD(4) * t380 + t406 * t305 - t404 * t306) - t409 * (-mrSges(5,1) * t352 + mrSges(5,3) * t345 + Ifges(5,4) * t390 + Ifges(5,2) * t389 + Ifges(5,6) * qJDD(4) - pkin(4) * t313 + qJ(5) * t427 + qJD(4) * t381 + t404 * t305 + t406 * t306) + mrSges(4,2) * t355 - mrSges(4,3) * t354 + (pkin(1) * (-t405 * mrSges(3,2) + t407 * (mrSges(3,1) - mrSges(4,2))) + Ifges(2,3) + Ifges(3,3) + Ifges(4,1)) * qJDD(1); -t409 * t307 + t412 * t308 + (m(3) + m(4)) * t401; t304; Ifges(5,5) * t390 + Ifges(5,6) * t389 + mrSges(5,1) * t344 - mrSges(5,2) * t345 + Ifges(6,5) * t367 + Ifges(6,6) * t366 + t375 * t357 + t374 * t358 + mrSges(6,1) * t328 - mrSges(6,2) * t329 + t408 * t316 + t411 * t315 + pkin(5) * t418 + pkin(8) * t426 + pkin(4) * t309 + (Ifges(5,3) + Ifges(6,3)) * qJDD(4) + (t412 * t380 + t409 * t381) * qJD(1); t313; t417;];
tauJ  = t1;
