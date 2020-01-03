% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRP6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:53:08
% EndTime: 2019-12-31 21:53:11
% DurationCPUTime: 1.84s
% Computational Cost: add. (13797->254), mult. (27849->313), div. (0->0), fcn. (18759->8), ass. (0->103)
t443 = Ifges(5,1) + Ifges(6,1);
t437 = Ifges(5,4) + Ifges(6,4);
t436 = Ifges(5,5) + Ifges(6,5);
t442 = Ifges(5,2) + Ifges(6,2);
t435 = Ifges(5,6) + Ifges(6,6);
t441 = Ifges(5,3) + Ifges(6,3);
t405 = sin(qJ(3));
t406 = sin(qJ(2));
t409 = cos(qJ(3));
t410 = cos(qJ(2));
t387 = (t406 * t405 - t410 * t409) * qJD(1);
t426 = qJD(1) * qJD(2);
t393 = qJDD(1) * t406 + t410 * t426;
t394 = qJDD(1) * t410 - t406 * t426;
t363 = -qJD(3) * t387 + t393 * t409 + t394 * t405;
t388 = (t410 * t405 + t406 * t409) * qJD(1);
t402 = qJD(2) + qJD(3);
t404 = sin(qJ(4));
t408 = cos(qJ(4));
t376 = -t388 * t404 + t402 * t408;
t401 = qJDD(2) + qJDD(3);
t338 = qJD(4) * t376 + t363 * t408 + t401 * t404;
t377 = t388 * t408 + t402 * t404;
t352 = -mrSges(6,1) * t376 + mrSges(6,2) * t377;
t362 = -qJD(3) * t388 - t393 * t405 + t394 * t409;
t428 = qJD(1) * t406;
t397 = qJD(2) * pkin(2) - pkin(7) * t428;
t403 = t410 ^ 2;
t412 = qJD(1) ^ 2;
t407 = sin(qJ(1));
t411 = cos(qJ(1));
t423 = g(1) * t407 - t411 * g(2);
t417 = -qJDD(1) * pkin(1) - t423;
t364 = -pkin(2) * t394 + t397 * t428 + (-pkin(7) * t403 - pkin(6)) * t412 + t417;
t327 = (t387 * t402 - t363) * pkin(8) + (t388 * t402 - t362) * pkin(3) + t364;
t419 = -t411 * g(1) - t407 * g(2);
t390 = -pkin(1) * t412 + qJDD(1) * pkin(6) + t419;
t433 = t390 * t406;
t439 = pkin(2) * t412;
t355 = qJDD(2) * pkin(2) - pkin(7) * t393 - t433 + (pkin(7) * t426 + t406 * t439 - g(3)) * t410;
t379 = -g(3) * t406 + t410 * t390;
t356 = pkin(7) * t394 - qJD(2) * t397 - t403 * t439 + t379;
t333 = t405 * t355 + t409 * t356;
t374 = pkin(3) * t387 - pkin(8) * t388;
t400 = t402 ^ 2;
t330 = -pkin(3) * t400 + pkin(8) * t401 - t374 * t387 + t333;
t322 = t408 * t327 - t330 * t404;
t361 = qJDD(4) - t362;
t383 = qJD(4) + t387;
t319 = -0.2e1 * qJD(5) * t377 + (t376 * t383 - t338) * qJ(5) + (t376 * t377 + t361) * pkin(4) + t322;
t365 = -mrSges(6,2) * t383 + mrSges(6,3) * t376;
t425 = m(6) * t319 + t361 * mrSges(6,1) + t383 * t365;
t316 = -mrSges(6,3) * t338 - t352 * t377 + t425;
t323 = t404 * t327 + t408 * t330;
t337 = -qJD(4) * t377 - t363 * t404 + t401 * t408;
t367 = pkin(4) * t383 - qJ(5) * t377;
t375 = t376 ^ 2;
t321 = -pkin(4) * t375 + qJ(5) * t337 + 0.2e1 * qJD(5) * t376 - t367 * t383 + t323;
t430 = t437 * t376 + t443 * t377 + t436 * t383;
t431 = -t442 * t376 - t437 * t377 - t435 * t383;
t440 = mrSges(5,1) * t322 + mrSges(6,1) * t319 - mrSges(5,2) * t323 - mrSges(6,2) * t321 + pkin(4) * t316 + t337 * t435 + t338 * t436 + t441 * t361 - t376 * t430 - t377 * t431;
t438 = -mrSges(5,2) - mrSges(6,2);
t373 = mrSges(4,1) * t387 + mrSges(4,2) * t388;
t381 = mrSges(4,1) * t402 - mrSges(4,3) * t388;
t353 = -mrSges(5,1) * t376 + mrSges(5,2) * t377;
t366 = -mrSges(5,2) * t383 + mrSges(5,3) * t376;
t309 = m(5) * t322 + mrSges(5,1) * t361 + t366 * t383 + (-t352 - t353) * t377 + (-mrSges(5,3) - mrSges(6,3)) * t338 + t425;
t424 = m(6) * t321 + t337 * mrSges(6,3) + t376 * t352;
t368 = mrSges(6,1) * t383 - mrSges(6,3) * t377;
t429 = -mrSges(5,1) * t383 + mrSges(5,3) * t377 - t368;
t312 = m(5) * t323 + mrSges(5,3) * t337 + t353 * t376 + t361 * t438 + t383 * t429 + t424;
t421 = -t309 * t404 + t408 * t312;
t303 = m(4) * t333 - mrSges(4,2) * t401 + mrSges(4,3) * t362 - t373 * t387 - t381 * t402 + t421;
t332 = t355 * t409 - t405 * t356;
t380 = -mrSges(4,2) * t402 - mrSges(4,3) * t387;
t329 = -pkin(3) * t401 - pkin(8) * t400 + t388 * t374 - t332;
t324 = -pkin(4) * t337 - qJ(5) * t375 + t367 * t377 + qJDD(5) + t329;
t420 = -m(6) * t324 + t337 * mrSges(6,1) + t376 * t365;
t413 = -m(5) * t329 + t337 * mrSges(5,1) + t438 * t338 + t376 * t366 + t429 * t377 + t420;
t314 = m(4) * t332 + mrSges(4,1) * t401 - mrSges(4,3) * t363 - t373 * t388 + t380 * t402 + t413;
t298 = t405 * t303 + t409 * t314;
t307 = t408 * t309 + t404 * t312;
t432 = -t435 * t376 - t436 * t377 - t441 * t383;
t427 = qJD(1) * t410;
t422 = t409 * t303 - t314 * t405;
t317 = mrSges(6,2) * t338 + t368 * t377 - t420;
t300 = -mrSges(5,1) * t329 + mrSges(5,3) * t323 - mrSges(6,1) * t324 + mrSges(6,3) * t321 - pkin(4) * t317 + qJ(5) * t424 + (-qJ(5) * t368 + t430) * t383 + t432 * t377 + (-qJ(5) * mrSges(6,2) + t435) * t361 + t437 * t338 + t442 * t337;
t305 = mrSges(5,2) * t329 + mrSges(6,2) * t324 - mrSges(5,3) * t322 - mrSges(6,3) * t319 - qJ(5) * t316 + t437 * t337 + t443 * t338 + t436 * t361 - t432 * t376 + t431 * t383;
t371 = Ifges(4,4) * t388 - Ifges(4,2) * t387 + Ifges(4,6) * t402;
t372 = Ifges(4,1) * t388 - Ifges(4,4) * t387 + Ifges(4,5) * t402;
t415 = mrSges(4,1) * t332 - mrSges(4,2) * t333 + Ifges(4,5) * t363 + Ifges(4,6) * t362 + Ifges(4,3) * t401 + pkin(3) * t413 + pkin(8) * t421 + t408 * t300 + t404 * t305 + t388 * t371 + t372 * t387;
t414 = m(4) * t364 - mrSges(4,1) * t362 + mrSges(4,2) * t363 + t380 * t387 + t381 * t388 + t307;
t396 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t427;
t395 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t428;
t392 = (-t410 * mrSges(3,1) + t406 * mrSges(3,2)) * qJD(1);
t389 = -pkin(6) * t412 + t417;
t386 = Ifges(3,5) * qJD(2) + (t406 * Ifges(3,1) + t410 * Ifges(3,4)) * qJD(1);
t385 = Ifges(3,6) * qJD(2) + (t406 * Ifges(3,4) + t410 * Ifges(3,2)) * qJD(1);
t378 = -g(3) * t410 - t433;
t370 = Ifges(4,5) * t388 - Ifges(4,6) * t387 + Ifges(4,3) * t402;
t297 = -mrSges(4,1) * t364 + mrSges(4,3) * t333 + Ifges(4,4) * t363 + Ifges(4,2) * t362 + Ifges(4,6) * t401 - pkin(3) * t307 - t388 * t370 + t402 * t372 - t440;
t296 = mrSges(4,2) * t364 - mrSges(4,3) * t332 + Ifges(4,1) * t363 + Ifges(4,4) * t362 + Ifges(4,5) * t401 - pkin(8) * t307 - t300 * t404 + t305 * t408 - t370 * t387 - t371 * t402;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t423 - mrSges(2,2) * t419 + t406 * (mrSges(3,2) * t389 - mrSges(3,3) * t378 + Ifges(3,1) * t393 + Ifges(3,4) * t394 + Ifges(3,5) * qJDD(2) - pkin(7) * t298 - qJD(2) * t385 + t409 * t296 - t405 * t297) + t410 * (-mrSges(3,1) * t389 + mrSges(3,3) * t379 + Ifges(3,4) * t393 + Ifges(3,2) * t394 + Ifges(3,6) * qJDD(2) - pkin(2) * t414 + pkin(7) * t422 + qJD(2) * t386 + t405 * t296 + t409 * t297) + pkin(1) * (-m(3) * t389 + mrSges(3,1) * t394 - mrSges(3,2) * t393 + (-t395 * t406 + t396 * t410) * qJD(1) - t414) + pkin(6) * (t410 * (m(3) * t379 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t394 - qJD(2) * t395 + t392 * t427 + t422) - t406 * (m(3) * t378 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t393 + qJD(2) * t396 - t392 * t428 + t298)); (t406 * t385 - t410 * t386) * qJD(1) + mrSges(3,1) * t378 - mrSges(3,2) * t379 + Ifges(3,5) * t393 + Ifges(3,6) * t394 + pkin(2) * t298 + t415 + Ifges(3,3) * qJDD(2); t415; t440; t317;];
tauJ = t1;
