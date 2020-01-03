% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRPR6
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRPR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR6_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR6_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR6_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:33
% EndTime: 2019-12-31 17:04:35
% DurationCPUTime: 2.06s
% Computational Cost: add. (21161->239), mult. (48930->310), div. (0->0), fcn. (31532->8), ass. (0->93)
t430 = qJD(1) ^ 2;
t444 = pkin(2) * t430;
t426 = sin(qJ(1));
t429 = cos(qJ(1));
t417 = -t429 * g(1) - t426 * g(2);
t406 = -t430 * pkin(1) + qJDD(1) * pkin(5) + t417;
t425 = sin(qJ(2));
t443 = t425 * t406;
t428 = cos(qJ(2));
t439 = qJD(1) * qJD(2);
t411 = t425 * qJDD(1) + t428 * t439;
t376 = qJDD(2) * pkin(2) - t411 * qJ(3) - t443 + (qJ(3) * t439 + t425 * t444 - g(3)) * t428;
t392 = -t425 * g(3) + t428 * t406;
t412 = t428 * qJDD(1) - t425 * t439;
t441 = qJD(1) * t425;
t413 = qJD(2) * pkin(2) - qJ(3) * t441;
t421 = t428 ^ 2;
t377 = t412 * qJ(3) - qJD(2) * t413 - t421 * t444 + t392;
t422 = sin(pkin(7));
t423 = cos(pkin(7));
t401 = (t422 * t428 + t423 * t425) * qJD(1);
t361 = -0.2e1 * qJD(3) * t401 + t423 * t376 - t422 * t377;
t390 = t423 * t411 + t422 * t412;
t400 = (-t422 * t425 + t423 * t428) * qJD(1);
t357 = (qJD(2) * t400 - t390) * pkin(6) + (t400 * t401 + qJDD(2)) * pkin(3) + t361;
t362 = 0.2e1 * qJD(3) * t400 + t422 * t376 + t423 * t377;
t389 = -t422 * t411 + t423 * t412;
t395 = qJD(2) * pkin(3) - t401 * pkin(6);
t399 = t400 ^ 2;
t358 = -t399 * pkin(3) + t389 * pkin(6) - qJD(2) * t395 + t362;
t424 = sin(qJ(4));
t427 = cos(qJ(4));
t355 = t427 * t357 - t424 * t358;
t384 = t427 * t400 - t424 * t401;
t366 = t384 * qJD(4) + t424 * t389 + t427 * t390;
t385 = t424 * t400 + t427 * t401;
t372 = -t384 * mrSges(5,1) + t385 * mrSges(5,2);
t420 = qJD(2) + qJD(4);
t379 = -t420 * mrSges(5,2) + t384 * mrSges(5,3);
t419 = qJDD(2) + qJDD(4);
t353 = m(5) * t355 + t419 * mrSges(5,1) - t366 * mrSges(5,3) - t385 * t372 + t420 * t379;
t356 = t424 * t357 + t427 * t358;
t365 = -t385 * qJD(4) + t427 * t389 - t424 * t390;
t380 = t420 * mrSges(5,1) - t385 * mrSges(5,3);
t354 = m(5) * t356 - t419 * mrSges(5,2) + t365 * mrSges(5,3) + t384 * t372 - t420 * t380;
t345 = t427 * t353 + t424 * t354;
t387 = -t400 * mrSges(4,1) + t401 * mrSges(4,2);
t393 = -qJD(2) * mrSges(4,2) + t400 * mrSges(4,3);
t343 = m(4) * t361 + qJDD(2) * mrSges(4,1) - t390 * mrSges(4,3) + qJD(2) * t393 - t401 * t387 + t345;
t394 = qJD(2) * mrSges(4,1) - t401 * mrSges(4,3);
t435 = -t424 * t353 + t427 * t354;
t344 = m(4) * t362 - qJDD(2) * mrSges(4,2) + t389 * mrSges(4,3) - qJD(2) * t394 + t400 * t387 + t435;
t339 = t423 * t343 + t422 * t344;
t391 = -t428 * g(3) - t443;
t410 = (-mrSges(3,1) * t428 + mrSges(3,2) * t425) * qJD(1);
t440 = qJD(1) * t428;
t415 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t440;
t337 = m(3) * t391 + qJDD(2) * mrSges(3,1) - t411 * mrSges(3,3) + qJD(2) * t415 - t410 * t441 + t339;
t414 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t441;
t436 = -t422 * t343 + t423 * t344;
t338 = m(3) * t392 - qJDD(2) * mrSges(3,2) + t412 * mrSges(3,3) - qJD(2) * t414 + t410 * t440 + t436;
t437 = -t425 * t337 + t428 * t338;
t330 = m(2) * t417 - t430 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t437;
t416 = t426 * g(1) - t429 * g(2);
t433 = -qJDD(1) * pkin(1) - t416;
t405 = -t430 * pkin(5) + t433;
t378 = -t412 * pkin(2) + qJDD(3) + t413 * t441 + (-qJ(3) * t421 - pkin(5)) * t430 + t433;
t360 = -t389 * pkin(3) - t399 * pkin(6) + t401 * t395 + t378;
t434 = m(5) * t360 - t365 * mrSges(5,1) + t366 * mrSges(5,2) - t384 * t379 + t385 * t380;
t432 = m(4) * t378 - t389 * mrSges(4,1) + t390 * mrSges(4,2) - t400 * t393 + t401 * t394 + t434;
t431 = -m(3) * t405 + t412 * mrSges(3,1) - t411 * mrSges(3,2) - t414 * t441 + t415 * t440 - t432;
t349 = m(2) * t416 + qJDD(1) * mrSges(2,1) - t430 * mrSges(2,2) + t431;
t442 = t426 * t330 + t429 * t349;
t331 = t428 * t337 + t425 * t338;
t438 = t429 * t330 - t426 * t349;
t404 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t425 + Ifges(3,4) * t428) * qJD(1);
t403 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t425 + Ifges(3,2) * t428) * qJD(1);
t402 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t425 + Ifges(3,6) * t428) * qJD(1);
t383 = Ifges(4,1) * t401 + Ifges(4,4) * t400 + Ifges(4,5) * qJD(2);
t382 = Ifges(4,4) * t401 + Ifges(4,2) * t400 + Ifges(4,6) * qJD(2);
t381 = Ifges(4,5) * t401 + Ifges(4,6) * t400 + Ifges(4,3) * qJD(2);
t369 = Ifges(5,1) * t385 + Ifges(5,4) * t384 + Ifges(5,5) * t420;
t368 = Ifges(5,4) * t385 + Ifges(5,2) * t384 + Ifges(5,6) * t420;
t367 = Ifges(5,5) * t385 + Ifges(5,6) * t384 + Ifges(5,3) * t420;
t347 = mrSges(5,2) * t360 - mrSges(5,3) * t355 + Ifges(5,1) * t366 + Ifges(5,4) * t365 + Ifges(5,5) * t419 + t384 * t367 - t420 * t368;
t346 = -mrSges(5,1) * t360 + mrSges(5,3) * t356 + Ifges(5,4) * t366 + Ifges(5,2) * t365 + Ifges(5,6) * t419 - t385 * t367 + t420 * t369;
t333 = mrSges(4,2) * t378 - mrSges(4,3) * t361 + Ifges(4,1) * t390 + Ifges(4,4) * t389 + Ifges(4,5) * qJDD(2) - pkin(6) * t345 - qJD(2) * t382 - t424 * t346 + t427 * t347 + t400 * t381;
t332 = -mrSges(4,1) * t378 + mrSges(4,3) * t362 + Ifges(4,4) * t390 + Ifges(4,2) * t389 + Ifges(4,6) * qJDD(2) - pkin(3) * t434 + pkin(6) * t435 + qJD(2) * t383 + t427 * t346 + t424 * t347 - t401 * t381;
t327 = mrSges(3,2) * t405 - mrSges(3,3) * t391 + Ifges(3,1) * t411 + Ifges(3,4) * t412 + Ifges(3,5) * qJDD(2) - qJ(3) * t339 - qJD(2) * t403 - t422 * t332 + t423 * t333 + t402 * t440;
t326 = -mrSges(3,1) * t405 + mrSges(3,3) * t392 + Ifges(3,4) * t411 + Ifges(3,2) * t412 + Ifges(3,6) * qJDD(2) - pkin(2) * t432 + qJ(3) * t436 + qJD(2) * t404 + t423 * t332 + t422 * t333 - t402 * t441;
t325 = Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t425 * t403 + t428 * t404) * qJD(1) + t430 * Ifges(2,5) - Ifges(5,3) * t419 - Ifges(3,5) * t411 - Ifges(3,6) * t412 + mrSges(2,3) * t417 + t400 * t383 - t401 * t382 - Ifges(4,6) * t389 - Ifges(4,5) * t390 - mrSges(3,1) * t391 + mrSges(3,2) * t392 + t384 * t369 - t385 * t368 - Ifges(5,6) * t365 - Ifges(5,5) * t366 - mrSges(4,1) * t361 + mrSges(4,2) * t362 - mrSges(5,1) * t355 + mrSges(5,2) * t356 - pkin(3) * t345 - pkin(2) * t339 - pkin(1) * t331;
t324 = -mrSges(2,2) * g(3) - mrSges(2,3) * t416 + Ifges(2,5) * qJDD(1) - t430 * Ifges(2,6) - pkin(5) * t331 - t425 * t326 + t428 * t327;
t1 = [-m(1) * g(1) + t438; -m(1) * g(2) + t442; (-m(1) - m(2)) * g(3) + t331; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t442 + t429 * t324 - t426 * t325; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t438 + t426 * t324 + t429 * t325; -mrSges(1,1) * g(2) + mrSges(2,1) * t416 + mrSges(1,2) * g(1) - mrSges(2,2) * t417 + Ifges(2,3) * qJDD(1) + pkin(1) * t431 + pkin(5) * t437 + t428 * t326 + t425 * t327;];
tauB = t1;
