% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRPP3
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRPP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP3_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP3_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP3_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:28
% EndTime: 2019-12-31 16:57:30
% DurationCPUTime: 1.17s
% Computational Cost: add. (7875->216), mult. (17812->267), div. (0->0), fcn. (10206->6), ass. (0->85)
t445 = Ifges(4,1) + Ifges(5,1);
t439 = Ifges(4,4) - Ifges(5,5);
t438 = Ifges(4,5) + Ifges(5,4);
t444 = Ifges(4,2) + Ifges(5,3);
t443 = -Ifges(5,2) - Ifges(4,3);
t437 = Ifges(4,6) - Ifges(5,6);
t442 = -2 * qJD(3);
t417 = qJD(1) ^ 2;
t441 = pkin(2) * t417;
t440 = -mrSges(4,3) - mrSges(5,2);
t436 = cos(pkin(6));
t413 = sin(qJ(1));
t415 = cos(qJ(1));
t405 = -t415 * g(1) - t413 * g(2);
t394 = -t417 * pkin(1) + qJDD(1) * pkin(5) + t405;
t412 = sin(qJ(2));
t435 = t412 * t394;
t414 = cos(qJ(2));
t427 = qJD(1) * qJD(2);
t399 = t412 * qJDD(1) + t414 * t427;
t357 = qJDD(2) * pkin(2) - t399 * qJ(3) - t435 + (qJ(3) * t427 + t412 * t441 - g(3)) * t414;
t380 = -t412 * g(3) + t414 * t394;
t400 = t414 * qJDD(1) - t412 * t427;
t429 = qJD(1) * t412;
t401 = qJD(2) * pkin(2) - qJ(3) * t429;
t410 = t414 ^ 2;
t358 = t400 * qJ(3) - qJD(2) * t401 - t410 * t441 + t380;
t411 = sin(pkin(6));
t428 = qJD(1) * t414;
t387 = t411 * t429 - t436 * t428;
t354 = t411 * t357 + t436 * t358 + t387 * t442;
t375 = t411 * t399 - t436 * t400;
t388 = (t411 * t414 + t436 * t412) * qJD(1);
t382 = qJD(2) * mrSges(4,1) - t388 * mrSges(4,3);
t369 = t387 * pkin(3) - t388 * qJ(4);
t416 = qJD(2) ^ 2;
t349 = -t416 * pkin(3) + qJDD(2) * qJ(4) + 0.2e1 * qJD(4) * qJD(2) - t387 * t369 + t354;
t383 = -qJD(2) * mrSges(5,1) + t388 * mrSges(5,2);
t426 = m(5) * t349 + qJDD(2) * mrSges(5,3) + qJD(2) * t383;
t370 = t387 * mrSges(5,1) - t388 * mrSges(5,3);
t430 = -t387 * mrSges(4,1) - t388 * mrSges(4,2) - t370;
t345 = m(4) * t354 - qJDD(2) * mrSges(4,2) - qJD(2) * t382 + t440 * t375 + t430 * t387 + t426;
t420 = t436 * t357 - t411 * t358;
t353 = t388 * t442 + t420;
t376 = t436 * t399 + t411 * t400;
t381 = -qJD(2) * mrSges(4,2) - t387 * mrSges(4,3);
t350 = -qJDD(2) * pkin(3) - t416 * qJ(4) + qJDD(4) + ((2 * qJD(3)) + t369) * t388 - t420;
t384 = -t387 * mrSges(5,2) + qJD(2) * mrSges(5,3);
t422 = -m(5) * t350 + qJDD(2) * mrSges(5,1) + qJD(2) * t384;
t346 = m(4) * t353 + qJDD(2) * mrSges(4,1) + qJD(2) * t381 + t440 * t376 + t430 * t388 + t422;
t339 = t411 * t345 + t436 * t346;
t379 = -t414 * g(3) - t435;
t398 = (-mrSges(3,1) * t414 + mrSges(3,2) * t412) * qJD(1);
t403 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t428;
t337 = m(3) * t379 + qJDD(2) * mrSges(3,1) - t399 * mrSges(3,3) + qJD(2) * t403 - t398 * t429 + t339;
t402 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t429;
t423 = t436 * t345 - t411 * t346;
t338 = m(3) * t380 - qJDD(2) * mrSges(3,2) + t400 * mrSges(3,3) - qJD(2) * t402 + t398 * t428 + t423;
t424 = -t412 * t337 + t414 * t338;
t330 = m(2) * t405 - t417 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t424;
t404 = t413 * g(1) - t415 * g(2);
t421 = -qJDD(1) * pkin(1) - t404;
t393 = -t417 * pkin(5) + t421;
t359 = -t400 * pkin(2) + qJDD(3) + t401 * t429 + (-qJ(3) * t410 - pkin(5)) * t417 + t421;
t352 = -0.2e1 * qJD(4) * t388 + (qJD(2) * t387 - t376) * qJ(4) + (qJD(2) * t388 + t375) * pkin(3) + t359;
t347 = m(5) * t352 + t375 * mrSges(5,1) - t376 * mrSges(5,3) - t388 * t383 + t387 * t384;
t419 = m(4) * t359 + t375 * mrSges(4,1) + t376 * mrSges(4,2) + t387 * t381 + t388 * t382 + t347;
t418 = -m(3) * t393 + t400 * mrSges(3,1) - t399 * mrSges(3,2) - t402 * t429 + t403 * t428 - t419;
t341 = m(2) * t404 + qJDD(1) * mrSges(2,1) - t417 * mrSges(2,2) + t418;
t434 = t413 * t330 + t415 * t341;
t331 = t414 * t337 + t412 * t338;
t433 = -t437 * qJD(2) + t444 * t387 - t439 * t388;
t432 = t443 * qJD(2) + t437 * t387 - t438 * t388;
t431 = t438 * qJD(2) - t439 * t387 + t445 * t388;
t425 = t415 * t330 - t413 * t341;
t391 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t412 + Ifges(3,4) * t414) * qJD(1);
t390 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t412 + Ifges(3,2) * t414) * qJD(1);
t389 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t412 + Ifges(3,6) * t414) * qJD(1);
t333 = mrSges(4,2) * t359 + mrSges(5,2) * t350 - mrSges(4,3) * t353 - mrSges(5,3) * t352 - qJ(4) * t347 + t433 * qJD(2) + t438 * qJDD(2) - t439 * t375 + t445 * t376 + t432 * t387;
t332 = -mrSges(4,1) * t359 - mrSges(5,1) * t352 + mrSges(5,2) * t349 + mrSges(4,3) * t354 - pkin(3) * t347 + t431 * qJD(2) + t437 * qJDD(2) - t444 * t375 + t439 * t376 + t432 * t388;
t327 = mrSges(3,2) * t393 - mrSges(3,3) * t379 + Ifges(3,1) * t399 + Ifges(3,4) * t400 + Ifges(3,5) * qJDD(2) - qJ(3) * t339 - qJD(2) * t390 - t411 * t332 + t436 * t333 + t389 * t428;
t326 = -mrSges(3,1) * t393 + mrSges(3,3) * t380 + Ifges(3,4) * t399 + Ifges(3,2) * t400 + Ifges(3,6) * qJDD(2) - pkin(2) * t419 + qJ(3) * t423 + qJD(2) * t391 + t436 * t332 + t411 * t333 - t389 * t429;
t325 = Ifges(2,6) * qJDD(1) + t417 * Ifges(2,5) - qJ(4) * t426 - pkin(3) * t422 - Ifges(3,5) * t399 - Ifges(3,6) * t400 + mrSges(2,3) * t405 - mrSges(3,1) * t379 + mrSges(3,2) * t380 - mrSges(4,1) * t353 + mrSges(4,2) * t354 - mrSges(5,3) * t349 + mrSges(5,1) * t350 + mrSges(2,1) * g(3) - pkin(2) * t339 - pkin(1) * t331 + (pkin(3) * t370 + t433) * t388 + (qJ(4) * t370 - t431) * t387 + (pkin(3) * mrSges(5,2) - t438) * t376 + (qJ(4) * mrSges(5,2) + t437) * t375 + (-t412 * t390 + t414 * t391) * qJD(1) + (-Ifges(3,3) + t443) * qJDD(2);
t324 = -mrSges(2,2) * g(3) - mrSges(2,3) * t404 + Ifges(2,5) * qJDD(1) - t417 * Ifges(2,6) - pkin(5) * t331 - t412 * t326 + t414 * t327;
t1 = [-m(1) * g(1) + t425; -m(1) * g(2) + t434; (-m(1) - m(2)) * g(3) + t331; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t434 + t415 * t324 - t413 * t325; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t425 + t413 * t324 + t415 * t325; -mrSges(1,1) * g(2) + mrSges(2,1) * t404 + mrSges(1,2) * g(1) - mrSges(2,2) * t405 + Ifges(2,3) * qJDD(1) + pkin(1) * t418 + pkin(5) * t424 + t414 * t326 + t412 * t327;];
tauB = t1;
