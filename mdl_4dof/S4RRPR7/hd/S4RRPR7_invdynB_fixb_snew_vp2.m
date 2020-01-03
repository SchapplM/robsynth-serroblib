% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRPR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR7_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR7_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR7_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:05
% EndTime: 2019-12-31 17:06:08
% DurationCPUTime: 1.83s
% Computational Cost: add. (16113->240), mult. (36298->310), div. (0->0), fcn. (22820->8), ass. (0->95)
t408 = sin(pkin(7));
t409 = cos(pkin(7));
t411 = sin(qJ(2));
t414 = cos(qJ(2));
t388 = (t408 * t411 - t409 * t414) * qJD(1);
t434 = 2 * qJD(3);
t417 = qJD(1) ^ 2;
t433 = pkin(2) * t417;
t412 = sin(qJ(1));
t415 = cos(qJ(1));
t405 = -t415 * g(1) - t412 * g(2);
t394 = -t417 * pkin(1) + qJDD(1) * pkin(5) + t405;
t432 = t411 * t394;
t428 = qJD(1) * qJD(2);
t399 = t411 * qJDD(1) + t414 * t428;
t363 = qJDD(2) * pkin(2) - t399 * qJ(3) - t432 + (qJ(3) * t428 + t411 * t433 - g(3)) * t414;
t382 = -t411 * g(3) + t414 * t394;
t400 = t414 * qJDD(1) - t411 * t428;
t430 = qJD(1) * t411;
t401 = qJD(2) * pkin(2) - qJ(3) * t430;
t407 = t414 ^ 2;
t364 = t400 * qJ(3) - qJD(2) * t401 - t407 * t433 + t382;
t353 = t408 * t363 + t409 * t364 - t388 * t434;
t389 = (t408 * t414 + t409 * t411) * qJD(1);
t373 = t388 * mrSges(4,1) + t389 * mrSges(4,2);
t377 = -t408 * t399 + t409 * t400;
t384 = qJD(2) * mrSges(4,1) - t389 * mrSges(4,3);
t374 = t388 * pkin(3) - t389 * pkin(6);
t416 = qJD(2) ^ 2;
t350 = -t416 * pkin(3) + qJDD(2) * pkin(6) - t388 * t374 + t353;
t404 = t412 * g(1) - t415 * g(2);
t421 = -qJDD(1) * pkin(1) - t404;
t366 = -t400 * pkin(2) + qJDD(3) + t401 * t430 + (-qJ(3) * t407 - pkin(5)) * t417 + t421;
t378 = t409 * t399 + t408 * t400;
t351 = (qJD(2) * t388 - t378) * pkin(6) + (qJD(2) * t389 - t377) * pkin(3) + t366;
t410 = sin(qJ(4));
t413 = cos(qJ(4));
t347 = -t410 * t350 + t413 * t351;
t379 = t413 * qJD(2) - t410 * t389;
t360 = t379 * qJD(4) + t410 * qJDD(2) + t413 * t378;
t380 = t410 * qJD(2) + t413 * t389;
t365 = -t379 * mrSges(5,1) + t380 * mrSges(5,2);
t387 = qJD(4) + t388;
t367 = -t387 * mrSges(5,2) + t379 * mrSges(5,3);
t376 = qJDD(4) - t377;
t345 = m(5) * t347 + t376 * mrSges(5,1) - t360 * mrSges(5,3) - t380 * t365 + t387 * t367;
t348 = t413 * t350 + t410 * t351;
t359 = -t380 * qJD(4) + t413 * qJDD(2) - t410 * t378;
t368 = t387 * mrSges(5,1) - t380 * mrSges(5,3);
t346 = m(5) * t348 - t376 * mrSges(5,2) + t359 * mrSges(5,3) + t379 * t365 - t387 * t368;
t424 = -t410 * t345 + t413 * t346;
t336 = m(4) * t353 - qJDD(2) * mrSges(4,2) + t377 * mrSges(4,3) - qJD(2) * t384 - t388 * t373 + t424;
t423 = -t409 * t363 + t408 * t364;
t352 = -0.2e1 * qJD(3) * t389 - t423;
t383 = -qJD(2) * mrSges(4,2) - t388 * mrSges(4,3);
t349 = -qJDD(2) * pkin(3) - t416 * pkin(6) + (t434 + t374) * t389 + t423;
t420 = -m(5) * t349 + t359 * mrSges(5,1) - t360 * mrSges(5,2) + t379 * t367 - t380 * t368;
t341 = m(4) * t352 + qJDD(2) * mrSges(4,1) - t378 * mrSges(4,3) + qJD(2) * t383 - t389 * t373 + t420;
t331 = t408 * t336 + t409 * t341;
t381 = -t414 * g(3) - t432;
t398 = (-mrSges(3,1) * t414 + mrSges(3,2) * t411) * qJD(1);
t429 = qJD(1) * t414;
t403 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t429;
t329 = m(3) * t381 + qJDD(2) * mrSges(3,1) - t399 * mrSges(3,3) + qJD(2) * t403 - t398 * t430 + t331;
t402 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t430;
t425 = t409 * t336 - t408 * t341;
t330 = m(3) * t382 - qJDD(2) * mrSges(3,2) + t400 * mrSges(3,3) - qJD(2) * t402 + t398 * t429 + t425;
t426 = -t411 * t329 + t414 * t330;
t322 = m(2) * t405 - t417 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t426;
t393 = -t417 * pkin(5) + t421;
t337 = t413 * t345 + t410 * t346;
t419 = m(4) * t366 - t377 * mrSges(4,1) + t378 * mrSges(4,2) + t388 * t383 + t389 * t384 + t337;
t418 = -m(3) * t393 + t400 * mrSges(3,1) - t399 * mrSges(3,2) - t402 * t430 + t403 * t429 - t419;
t333 = m(2) * t404 + qJDD(1) * mrSges(2,1) - t417 * mrSges(2,2) + t418;
t431 = t412 * t322 + t415 * t333;
t323 = t414 * t329 + t411 * t330;
t427 = t415 * t322 - t412 * t333;
t392 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t411 + Ifges(3,4) * t414) * qJD(1);
t391 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t411 + Ifges(3,2) * t414) * qJD(1);
t390 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t411 + Ifges(3,6) * t414) * qJD(1);
t371 = Ifges(4,1) * t389 - Ifges(4,4) * t388 + Ifges(4,5) * qJD(2);
t370 = Ifges(4,4) * t389 - Ifges(4,2) * t388 + Ifges(4,6) * qJD(2);
t369 = Ifges(4,5) * t389 - Ifges(4,6) * t388 + Ifges(4,3) * qJD(2);
t356 = Ifges(5,1) * t380 + Ifges(5,4) * t379 + Ifges(5,5) * t387;
t355 = Ifges(5,4) * t380 + Ifges(5,2) * t379 + Ifges(5,6) * t387;
t354 = Ifges(5,5) * t380 + Ifges(5,6) * t379 + Ifges(5,3) * t387;
t339 = mrSges(5,2) * t349 - mrSges(5,3) * t347 + Ifges(5,1) * t360 + Ifges(5,4) * t359 + Ifges(5,5) * t376 + t379 * t354 - t387 * t355;
t338 = -mrSges(5,1) * t349 + mrSges(5,3) * t348 + Ifges(5,4) * t360 + Ifges(5,2) * t359 + Ifges(5,6) * t376 - t380 * t354 + t387 * t356;
t325 = -mrSges(4,1) * t366 - mrSges(5,1) * t347 + mrSges(5,2) * t348 + mrSges(4,3) * t353 + Ifges(4,4) * t378 - Ifges(5,5) * t360 + Ifges(4,2) * t377 + Ifges(4,6) * qJDD(2) - Ifges(5,6) * t359 - Ifges(5,3) * t376 - pkin(3) * t337 + qJD(2) * t371 - t380 * t355 + t379 * t356 - t389 * t369;
t324 = mrSges(4,2) * t366 - mrSges(4,3) * t352 + Ifges(4,1) * t378 + Ifges(4,4) * t377 + Ifges(4,5) * qJDD(2) - pkin(6) * t337 - qJD(2) * t370 - t410 * t338 + t413 * t339 - t388 * t369;
t319 = mrSges(3,2) * t393 - mrSges(3,3) * t381 + Ifges(3,1) * t399 + Ifges(3,4) * t400 + Ifges(3,5) * qJDD(2) - qJ(3) * t331 - qJD(2) * t391 + t409 * t324 - t408 * t325 + t390 * t429;
t318 = -mrSges(3,1) * t393 + mrSges(3,3) * t382 + Ifges(3,4) * t399 + Ifges(3,2) * t400 + Ifges(3,6) * qJDD(2) - pkin(2) * t419 + qJ(3) * t425 + qJD(2) * t392 + t408 * t324 + t409 * t325 - t390 * t430;
t317 = -pkin(1) * t323 + mrSges(2,3) * t405 - pkin(2) * t331 - Ifges(3,5) * t399 - Ifges(3,6) * t400 - mrSges(3,1) * t381 + mrSges(3,2) * t382 - t410 * t339 - t413 * t338 - pkin(3) * t420 - pkin(6) * t424 - Ifges(4,5) * t378 - Ifges(4,6) * t377 - mrSges(4,1) * t352 + mrSges(4,2) * t353 + mrSges(2,1) * g(3) - t389 * t370 - t388 * t371 + t417 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t411 * t391 + t414 * t392) * qJD(1);
t316 = -mrSges(2,2) * g(3) - mrSges(2,3) * t404 + Ifges(2,5) * qJDD(1) - t417 * Ifges(2,6) - pkin(5) * t323 - t411 * t318 + t414 * t319;
t1 = [-m(1) * g(1) + t427; -m(1) * g(2) + t431; (-m(1) - m(2)) * g(3) + t323; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t431 + t415 * t316 - t412 * t317; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t427 + t412 * t316 + t415 * t317; -mrSges(1,1) * g(2) + mrSges(2,1) * t404 + mrSges(1,2) * g(1) - mrSges(2,2) * t405 + Ifges(2,3) * qJDD(1) + pkin(1) * t418 + pkin(5) * t426 + t414 * t318 + t411 * t319;];
tauB = t1;
