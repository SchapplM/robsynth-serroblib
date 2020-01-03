% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR7_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR7_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR7_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:37
% EndTime: 2019-12-31 17:59:38
% DurationCPUTime: 1.35s
% Computational Cost: add. (12680->221), mult. (22747->268), div. (0->0), fcn. (11792->8), ass. (0->92)
t416 = sin(qJ(1));
t419 = cos(qJ(1));
t399 = t416 * g(1) - t419 * g(2);
t391 = qJDD(1) * pkin(1) + t399;
t400 = -t419 * g(1) - t416 * g(2);
t421 = qJD(1) ^ 2;
t393 = -t421 * pkin(1) + t400;
t412 = sin(pkin(8));
t413 = cos(pkin(8));
t375 = t412 * t391 + t413 * t393;
t445 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) - t375;
t444 = -pkin(2) - pkin(6);
t443 = mrSges(3,1) - mrSges(4,2);
t442 = -Ifges(4,4) + Ifges(3,5);
t441 = Ifges(4,5) - Ifges(3,6);
t409 = -g(3) + qJDD(2);
t415 = sin(qJ(4));
t440 = t415 * t409;
t374 = t413 * t391 - t412 * t393;
t426 = -t421 * qJ(3) + qJDD(3) - t374;
t364 = t444 * qJDD(1) + t426;
t418 = cos(qJ(4));
t360 = t415 * t364 + t418 * t409;
t392 = (mrSges(5,1) * t415 + mrSges(5,2) * t418) * qJD(1);
t436 = qJD(1) * qJD(4);
t432 = t418 * t436;
t395 = -t415 * qJDD(1) - t432;
t438 = qJD(1) * t418;
t398 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t438;
t363 = t444 * t421 - t445;
t433 = t415 * t436;
t396 = t418 * qJDD(1) - t433;
t356 = (-t396 + t433) * pkin(7) + (-t395 + t432) * pkin(4) + t363;
t394 = (pkin(4) * t415 - pkin(7) * t418) * qJD(1);
t420 = qJD(4) ^ 2;
t437 = t415 * qJD(1);
t358 = -t420 * pkin(4) + qJDD(4) * pkin(7) - t394 * t437 + t360;
t414 = sin(qJ(5));
t417 = cos(qJ(5));
t354 = t417 * t356 - t414 * t358;
t389 = t417 * qJD(4) - t414 * t438;
t373 = t389 * qJD(5) + t414 * qJDD(4) + t417 * t396;
t390 = t414 * qJD(4) + t417 * t438;
t376 = -t389 * mrSges(6,1) + t390 * mrSges(6,2);
t401 = qJD(5) + t437;
t377 = -t401 * mrSges(6,2) + t389 * mrSges(6,3);
t388 = qJDD(5) - t395;
t352 = m(6) * t354 + t388 * mrSges(6,1) - t373 * mrSges(6,3) - t390 * t376 + t401 * t377;
t355 = t414 * t356 + t417 * t358;
t372 = -t390 * qJD(5) + t417 * qJDD(4) - t414 * t396;
t378 = t401 * mrSges(6,1) - t390 * mrSges(6,3);
t353 = m(6) * t355 - t388 * mrSges(6,2) + t372 * mrSges(6,3) + t389 * t376 - t401 * t378;
t428 = -t414 * t352 + t417 * t353;
t344 = m(5) * t360 - qJDD(4) * mrSges(5,2) + t395 * mrSges(5,3) - qJD(4) * t398 - t392 * t437 + t428;
t359 = t418 * t364 - t440;
t397 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t437;
t357 = -qJDD(4) * pkin(4) - t420 * pkin(7) + t440 + (qJD(1) * t394 - t364) * t418;
t423 = -m(6) * t357 + t372 * mrSges(6,1) - t373 * mrSges(6,2) + t389 * t377 - t390 * t378;
t348 = m(5) * t359 + qJDD(4) * mrSges(5,1) - t396 * mrSges(5,3) + qJD(4) * t397 - t392 * t438 + t423;
t339 = t415 * t344 + t418 * t348;
t366 = -qJDD(1) * pkin(2) + t426;
t425 = -m(4) * t366 + t421 * mrSges(4,3) - t339;
t337 = m(3) * t374 - t421 * mrSges(3,2) + t443 * qJDD(1) + t425;
t365 = t421 * pkin(2) + t445;
t345 = t417 * t352 + t414 * t353;
t424 = -m(5) * t363 + t395 * mrSges(5,1) - t396 * mrSges(5,2) - t397 * t437 - t398 * t438 - t345;
t422 = -m(4) * t365 + t421 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t424;
t342 = m(3) * t375 - t421 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t422;
t333 = t413 * t337 + t412 * t342;
t331 = m(2) * t399 + qJDD(1) * mrSges(2,1) - t421 * mrSges(2,2) + t333;
t430 = -t412 * t337 + t413 * t342;
t332 = m(2) * t400 - t421 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t430;
t439 = t419 * t331 + t416 * t332;
t431 = -t416 * t331 + t419 * t332;
t429 = t418 * t344 - t415 * t348;
t338 = m(4) * t409 + t429;
t427 = m(3) * t409 + t338;
t384 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t418 - Ifges(5,4) * t415) * qJD(1);
t383 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t418 - Ifges(5,2) * t415) * qJD(1);
t382 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t418 - Ifges(5,6) * t415) * qJD(1);
t369 = Ifges(6,1) * t390 + Ifges(6,4) * t389 + Ifges(6,5) * t401;
t368 = Ifges(6,4) * t390 + Ifges(6,2) * t389 + Ifges(6,6) * t401;
t367 = Ifges(6,5) * t390 + Ifges(6,6) * t389 + Ifges(6,3) * t401;
t347 = mrSges(6,2) * t357 - mrSges(6,3) * t354 + Ifges(6,1) * t373 + Ifges(6,4) * t372 + Ifges(6,5) * t388 + t389 * t367 - t401 * t368;
t346 = -mrSges(6,1) * t357 + mrSges(6,3) * t355 + Ifges(6,4) * t373 + Ifges(6,2) * t372 + Ifges(6,6) * t388 - t390 * t367 + t401 * t369;
t335 = -mrSges(5,1) * t363 - mrSges(6,1) * t354 + mrSges(6,2) * t355 + mrSges(5,3) * t360 + Ifges(5,4) * t396 - Ifges(6,5) * t373 + Ifges(5,2) * t395 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t372 - Ifges(6,3) * t388 - pkin(4) * t345 + qJD(4) * t384 - t390 * t368 + t389 * t369 - t382 * t438;
t334 = mrSges(5,2) * t363 - mrSges(5,3) * t359 + Ifges(5,1) * t396 + Ifges(5,4) * t395 + Ifges(5,5) * qJDD(4) - pkin(7) * t345 - qJD(4) * t383 - t414 * t346 + t417 * t347 - t382 * t437;
t327 = -qJ(3) * t338 - mrSges(3,3) * t374 + pkin(3) * t339 + mrSges(4,1) * t366 + t414 * t347 + t417 * t346 + pkin(4) * t423 + pkin(7) * t428 + Ifges(5,5) * t396 + Ifges(5,6) * t395 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t359 - mrSges(5,2) * t360 + t441 * t421 + (mrSges(3,2) - mrSges(4,3)) * t409 + t442 * qJDD(1) + (t418 * t383 + t415 * t384) * qJD(1);
t326 = -mrSges(4,1) * t365 + mrSges(3,3) * t375 - pkin(2) * t338 - pkin(3) * t424 - pkin(6) * t429 - t441 * qJDD(1) - t415 * t334 - t418 * t335 - t443 * t409 + t442 * t421;
t325 = -mrSges(2,2) * g(3) - mrSges(2,3) * t399 + Ifges(2,5) * qJDD(1) - t421 * Ifges(2,6) - qJ(2) * t333 - t412 * t326 + t413 * t327;
t324 = mrSges(2,1) * g(3) + mrSges(2,3) * t400 + t421 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t427 + qJ(2) * t430 + t413 * t326 + t412 * t327;
t1 = [-m(1) * g(1) + t431; -m(1) * g(2) + t439; (-m(1) - m(2)) * g(3) + t427; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t439 - t416 * t324 + t419 * t325; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t431 + t419 * t324 + t416 * t325; pkin(1) * t333 + mrSges(2,1) * t399 - mrSges(2,2) * t400 + qJ(3) * t422 + pkin(2) * t425 + t418 * t334 - t415 * t335 - pkin(6) * t339 + mrSges(3,1) * t374 - mrSges(3,2) * t375 + mrSges(4,2) * t366 - mrSges(4,3) * t365 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(2) * mrSges(4,2) + Ifges(4,1) + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB = t1;
