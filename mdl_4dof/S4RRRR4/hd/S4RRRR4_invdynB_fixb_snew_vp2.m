% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRRR4
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRRR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR4_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR4_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR4_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:50
% EndTime: 2019-12-31 17:25:52
% DurationCPUTime: 1.94s
% Computational Cost: add. (18723->240), mult. (38232->309), div. (0->0), fcn. (24564->8), ass. (0->95)
t419 = sin(qJ(3));
t420 = sin(qJ(2));
t423 = cos(qJ(3));
t424 = cos(qJ(2));
t398 = (t419 * t420 - t423 * t424) * qJD(1);
t426 = qJD(1) ^ 2;
t441 = pkin(2) * t426;
t421 = sin(qJ(1));
t425 = cos(qJ(1));
t412 = -t425 * g(1) - t421 * g(2);
t401 = -t426 * pkin(1) + qJDD(1) * pkin(5) + t412;
t440 = t420 * t401;
t436 = qJD(1) * qJD(2);
t406 = t420 * qJDD(1) + t424 * t436;
t372 = qJDD(2) * pkin(2) - t406 * pkin(6) - t440 + (pkin(6) * t436 + t420 * t441 - g(3)) * t424;
t390 = -t420 * g(3) + t424 * t401;
t407 = t424 * qJDD(1) - t420 * t436;
t438 = qJD(1) * t420;
t410 = qJD(2) * pkin(2) - pkin(6) * t438;
t417 = t424 ^ 2;
t373 = t407 * pkin(6) - qJD(2) * t410 - t417 * t441 + t390;
t361 = t419 * t372 + t423 * t373;
t399 = (t419 * t424 + t420 * t423) * qJD(1);
t376 = -t399 * qJD(3) - t419 * t406 + t423 * t407;
t385 = t398 * mrSges(4,1) + t399 * mrSges(4,2);
t416 = qJD(2) + qJD(3);
t392 = t416 * mrSges(4,1) - t399 * mrSges(4,3);
t415 = qJDD(2) + qJDD(3);
t377 = -t398 * qJD(3) + t423 * t406 + t419 * t407;
t411 = t421 * g(1) - t425 * g(2);
t430 = -qJDD(1) * pkin(1) - t411;
t378 = -t407 * pkin(2) + t410 * t438 + (-pkin(6) * t417 - pkin(5)) * t426 + t430;
t357 = (t398 * t416 - t377) * pkin(7) + (t399 * t416 - t376) * pkin(3) + t378;
t386 = t398 * pkin(3) - t399 * pkin(7);
t414 = t416 ^ 2;
t359 = -t414 * pkin(3) + t415 * pkin(7) - t398 * t386 + t361;
t418 = sin(qJ(4));
t422 = cos(qJ(4));
t355 = t422 * t357 - t418 * t359;
t387 = -t418 * t399 + t422 * t416;
t364 = t387 * qJD(4) + t422 * t377 + t418 * t415;
t388 = t422 * t399 + t418 * t416;
t371 = -t387 * mrSges(5,1) + t388 * mrSges(5,2);
t375 = qJDD(4) - t376;
t394 = qJD(4) + t398;
t379 = -t394 * mrSges(5,2) + t387 * mrSges(5,3);
t353 = m(5) * t355 + t375 * mrSges(5,1) - t364 * mrSges(5,3) - t388 * t371 + t394 * t379;
t356 = t418 * t357 + t422 * t359;
t363 = -t388 * qJD(4) - t418 * t377 + t422 * t415;
t380 = t394 * mrSges(5,1) - t388 * mrSges(5,3);
t354 = m(5) * t356 - t375 * mrSges(5,2) + t363 * mrSges(5,3) + t387 * t371 - t394 * t380;
t432 = -t418 * t353 + t422 * t354;
t344 = m(4) * t361 - t415 * mrSges(4,2) + t376 * mrSges(4,3) - t398 * t385 - t416 * t392 + t432;
t360 = t423 * t372 - t419 * t373;
t391 = -t416 * mrSges(4,2) - t398 * mrSges(4,3);
t358 = -t415 * pkin(3) - t414 * pkin(7) + t399 * t386 - t360;
t429 = -m(5) * t358 + t363 * mrSges(5,1) - t364 * mrSges(5,2) + t387 * t379 - t388 * t380;
t349 = m(4) * t360 + t415 * mrSges(4,1) - t377 * mrSges(4,3) - t399 * t385 + t416 * t391 + t429;
t339 = t419 * t344 + t423 * t349;
t389 = -t424 * g(3) - t440;
t405 = (-mrSges(3,1) * t424 + mrSges(3,2) * t420) * qJD(1);
t437 = qJD(1) * t424;
t409 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t437;
t337 = m(3) * t389 + qJDD(2) * mrSges(3,1) - t406 * mrSges(3,3) + qJD(2) * t409 - t405 * t438 + t339;
t408 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t438;
t433 = t423 * t344 - t419 * t349;
t338 = m(3) * t390 - qJDD(2) * mrSges(3,2) + t407 * mrSges(3,3) - qJD(2) * t408 + t405 * t437 + t433;
t434 = -t420 * t337 + t424 * t338;
t330 = m(2) * t412 - t426 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t434;
t400 = -t426 * pkin(5) + t430;
t345 = t422 * t353 + t418 * t354;
t428 = m(4) * t378 - t376 * mrSges(4,1) + t377 * mrSges(4,2) + t398 * t391 + t399 * t392 + t345;
t427 = -m(3) * t400 + t407 * mrSges(3,1) - t406 * mrSges(3,2) - t408 * t438 + t409 * t437 - t428;
t341 = m(2) * t411 + qJDD(1) * mrSges(2,1) - t426 * mrSges(2,2) + t427;
t439 = t421 * t330 + t425 * t341;
t331 = t424 * t337 + t420 * t338;
t435 = t425 * t330 - t421 * t341;
t397 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t420 + Ifges(3,4) * t424) * qJD(1);
t396 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t420 + Ifges(3,2) * t424) * qJD(1);
t395 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t420 + Ifges(3,6) * t424) * qJD(1);
t383 = Ifges(4,1) * t399 - Ifges(4,4) * t398 + Ifges(4,5) * t416;
t382 = Ifges(4,4) * t399 - Ifges(4,2) * t398 + Ifges(4,6) * t416;
t381 = Ifges(4,5) * t399 - Ifges(4,6) * t398 + Ifges(4,3) * t416;
t367 = Ifges(5,1) * t388 + Ifges(5,4) * t387 + Ifges(5,5) * t394;
t366 = Ifges(5,4) * t388 + Ifges(5,2) * t387 + Ifges(5,6) * t394;
t365 = Ifges(5,5) * t388 + Ifges(5,6) * t387 + Ifges(5,3) * t394;
t347 = mrSges(5,2) * t358 - mrSges(5,3) * t355 + Ifges(5,1) * t364 + Ifges(5,4) * t363 + Ifges(5,5) * t375 + t387 * t365 - t394 * t366;
t346 = -mrSges(5,1) * t358 + mrSges(5,3) * t356 + Ifges(5,4) * t364 + Ifges(5,2) * t363 + Ifges(5,6) * t375 - t388 * t365 + t394 * t367;
t333 = -mrSges(4,1) * t378 - mrSges(5,1) * t355 + mrSges(5,2) * t356 + mrSges(4,3) * t361 + Ifges(4,4) * t377 - Ifges(5,5) * t364 + Ifges(4,2) * t376 + Ifges(4,6) * t415 - Ifges(5,6) * t363 - Ifges(5,3) * t375 - pkin(3) * t345 - t388 * t366 + t387 * t367 - t399 * t381 + t416 * t383;
t332 = mrSges(4,2) * t378 - mrSges(4,3) * t360 + Ifges(4,1) * t377 + Ifges(4,4) * t376 + Ifges(4,5) * t415 - pkin(7) * t345 - t418 * t346 + t422 * t347 - t398 * t381 - t416 * t382;
t327 = mrSges(3,2) * t400 - mrSges(3,3) * t389 + Ifges(3,1) * t406 + Ifges(3,4) * t407 + Ifges(3,5) * qJDD(2) - pkin(6) * t339 - qJD(2) * t396 + t423 * t332 - t419 * t333 + t395 * t437;
t326 = -mrSges(3,1) * t400 + mrSges(3,3) * t390 + Ifges(3,4) * t406 + Ifges(3,2) * t407 + Ifges(3,6) * qJDD(2) - pkin(2) * t428 + pkin(6) * t433 + qJD(2) * t397 + t419 * t332 + t423 * t333 - t395 * t438;
t325 = -pkin(1) * t331 + mrSges(2,3) * t412 - pkin(2) * t339 - Ifges(3,5) * t406 - Ifges(3,6) * t407 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t389 + mrSges(3,2) * t390 - Ifges(4,3) * t415 - mrSges(4,1) * t360 + mrSges(4,2) * t361 - t418 * t347 - t422 * t346 - pkin(3) * t429 - pkin(7) * t432 - Ifges(4,5) * t377 - Ifges(4,6) * t376 + mrSges(2,1) * g(3) - t399 * t382 - t398 * t383 + t426 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-t420 * t396 + t424 * t397) * qJD(1);
t324 = -mrSges(2,2) * g(3) - mrSges(2,3) * t411 + Ifges(2,5) * qJDD(1) - t426 * Ifges(2,6) - pkin(5) * t331 - t420 * t326 + t424 * t327;
t1 = [-m(1) * g(1) + t435; -m(1) * g(2) + t439; (-m(1) - m(2)) * g(3) + t331; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t439 + t425 * t324 - t421 * t325; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t435 + t421 * t324 + t425 * t325; -mrSges(1,1) * g(2) + mrSges(2,1) * t411 + mrSges(1,2) * g(1) - mrSges(2,2) * t412 + Ifges(2,3) * qJDD(1) + pkin(1) * t427 + pkin(5) * t434 + t424 * t326 + t420 * t327;];
tauB = t1;
