% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRR6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:47
% EndTime: 2019-12-05 17:09:51
% DurationCPUTime: 2.30s
% Computational Cost: add. (33587->219), mult. (42670->281), div. (0->0), fcn. (26262->10), ass. (0->93)
t450 = m(3) + m(4);
t449 = cos(pkin(9));
t423 = qJD(2) + qJD(3);
t428 = sin(qJ(4));
t448 = t423 * t428;
t432 = cos(qJ(4));
t447 = t423 * t432;
t426 = sin(pkin(9));
t415 = -t449 * g(1) - t426 * g(2);
t425 = -g(3) + qJDD(1);
t430 = sin(qJ(2));
t434 = cos(qJ(2));
t396 = -t430 * t415 + t434 * t425;
t394 = qJDD(2) * pkin(2) + t396;
t397 = t434 * t415 + t430 * t425;
t435 = qJD(2) ^ 2;
t395 = -t435 * pkin(2) + t397;
t429 = sin(qJ(3));
t433 = cos(qJ(3));
t383 = t429 * t394 + t433 * t395;
t419 = t423 ^ 2;
t421 = qJDD(2) + qJDD(3);
t378 = -t419 * pkin(3) + t421 * pkin(7) + t383;
t414 = t426 * g(1) - t449 * g(2);
t374 = -t428 * t378 - t432 * t414;
t445 = qJD(4) * t423;
t444 = t432 * t445;
t406 = t428 * t421 + t444;
t371 = (-t406 + t444) * pkin(8) + (t419 * t428 * t432 + qJDD(4)) * pkin(4) + t374;
t375 = t432 * t378 - t428 * t414;
t407 = t432 * t421 - t428 * t445;
t413 = qJD(4) * pkin(4) - pkin(8) * t448;
t424 = t432 ^ 2;
t372 = -t424 * t419 * pkin(4) + t407 * pkin(8) - qJD(4) * t413 + t375;
t427 = sin(qJ(5));
t431 = cos(qJ(5));
t369 = t431 * t371 - t427 * t372;
t401 = (-t427 * t428 + t431 * t432) * t423;
t381 = t401 * qJD(5) + t431 * t406 + t427 * t407;
t402 = (t427 * t432 + t428 * t431) * t423;
t388 = -t401 * mrSges(6,1) + t402 * mrSges(6,2);
t422 = qJD(4) + qJD(5);
t389 = -t422 * mrSges(6,2) + t401 * mrSges(6,3);
t420 = qJDD(4) + qJDD(5);
t367 = m(6) * t369 + t420 * mrSges(6,1) - t381 * mrSges(6,3) - t402 * t388 + t422 * t389;
t370 = t427 * t371 + t431 * t372;
t380 = -t402 * qJD(5) - t427 * t406 + t431 * t407;
t390 = t422 * mrSges(6,1) - t402 * mrSges(6,3);
t368 = m(6) * t370 - t420 * mrSges(6,2) + t380 * mrSges(6,3) + t401 * t388 - t422 * t390;
t359 = t431 * t367 + t427 * t368;
t405 = (-mrSges(5,1) * t432 + mrSges(5,2) * t428) * t423;
t412 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t447;
t357 = m(5) * t374 + qJDD(4) * mrSges(5,1) - t406 * mrSges(5,3) + qJD(4) * t412 - t405 * t448 + t359;
t411 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t448;
t439 = -t427 * t367 + t431 * t368;
t358 = m(5) * t375 - qJDD(4) * mrSges(5,2) + t407 * mrSges(5,3) - qJD(4) * t411 + t405 * t447 + t439;
t440 = -t428 * t357 + t432 * t358;
t350 = m(4) * t383 - t419 * mrSges(4,1) - t421 * mrSges(4,2) + t440;
t382 = t433 * t394 - t429 * t395;
t438 = -t421 * pkin(3) - t382;
t377 = -t419 * pkin(7) + t438;
t373 = t413 * t448 - t407 * pkin(4) + (-pkin(8) * t424 - pkin(7)) * t419 + t438;
t437 = m(6) * t373 - t380 * mrSges(6,1) + t381 * mrSges(6,2) - t401 * t389 + t402 * t390;
t436 = -m(5) * t377 + t407 * mrSges(5,1) - t406 * mrSges(5,2) - t411 * t448 + t412 * t447 - t437;
t363 = m(4) * t382 + t421 * mrSges(4,1) - t419 * mrSges(4,2) + t436;
t346 = t429 * t350 + t433 * t363;
t344 = m(3) * t396 + qJDD(2) * mrSges(3,1) - t435 * mrSges(3,2) + t346;
t441 = t433 * t350 - t429 * t363;
t345 = m(3) * t397 - t435 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t441;
t442 = -t430 * t344 + t434 * t345;
t337 = m(2) * t415 + t442;
t353 = t432 * t357 + t428 * t358;
t352 = (m(2) + t450) * t414 - t353;
t446 = t426 * t337 + t449 * t352;
t338 = t434 * t344 + t430 * t345;
t443 = t449 * t337 - t426 * t352;
t400 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t428 + Ifges(5,4) * t432) * t423;
t399 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t428 + Ifges(5,2) * t432) * t423;
t398 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t428 + Ifges(5,6) * t432) * t423;
t386 = Ifges(6,1) * t402 + Ifges(6,4) * t401 + Ifges(6,5) * t422;
t385 = Ifges(6,4) * t402 + Ifges(6,2) * t401 + Ifges(6,6) * t422;
t384 = Ifges(6,5) * t402 + Ifges(6,6) * t401 + Ifges(6,3) * t422;
t361 = mrSges(6,2) * t373 - mrSges(6,3) * t369 + Ifges(6,1) * t381 + Ifges(6,4) * t380 + Ifges(6,5) * t420 + t401 * t384 - t422 * t385;
t360 = -mrSges(6,1) * t373 + mrSges(6,3) * t370 + Ifges(6,4) * t381 + Ifges(6,2) * t380 + Ifges(6,6) * t420 - t402 * t384 + t422 * t386;
t347 = mrSges(5,2) * t377 - mrSges(5,3) * t374 + Ifges(5,1) * t406 + Ifges(5,4) * t407 + Ifges(5,5) * qJDD(4) - pkin(8) * t359 - qJD(4) * t399 - t427 * t360 + t431 * t361 + t398 * t447;
t340 = -mrSges(5,1) * t377 + mrSges(5,3) * t375 + Ifges(5,4) * t406 + Ifges(5,2) * t407 + Ifges(5,6) * qJDD(4) - pkin(4) * t437 + pkin(8) * t439 + qJD(4) * t400 + t431 * t360 + t427 * t361 - t398 * t448;
t339 = Ifges(4,6) * t421 + t419 * Ifges(4,5) + mrSges(4,1) * t414 + mrSges(4,3) * t383 - Ifges(5,5) * t406 - Ifges(5,6) * t407 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t374 + mrSges(5,2) * t375 - Ifges(6,5) * t381 - Ifges(6,6) * t380 - Ifges(6,3) * t420 - t402 * t385 + t401 * t386 - mrSges(6,1) * t369 + mrSges(6,2) * t370 - pkin(4) * t359 - pkin(3) * t353 + (-t428 * t399 + t432 * t400) * t423;
t334 = -mrSges(4,2) * t414 - mrSges(4,3) * t382 + Ifges(4,5) * t421 - t419 * Ifges(4,6) - pkin(7) * t353 - t428 * t340 + t432 * t347;
t333 = -mrSges(3,2) * t414 - mrSges(3,3) * t396 + Ifges(3,5) * qJDD(2) - t435 * Ifges(3,6) - pkin(6) * t346 + t433 * t334 - t429 * t339;
t332 = -mrSges(2,1) * t425 - mrSges(3,1) * t396 - mrSges(4,1) * t382 + mrSges(3,2) * t397 + mrSges(4,2) * t383 + mrSges(2,3) * t415 - Ifges(3,3) * qJDD(2) - Ifges(4,3) * t421 - pkin(1) * t338 - pkin(2) * t346 - pkin(3) * t436 - pkin(7) * t440 - t432 * t340 - t428 * t347;
t331 = Ifges(3,6) * qJDD(2) + t435 * Ifges(3,5) + mrSges(3,1) * t414 + mrSges(3,3) * t397 + t429 * t334 + t433 * t339 - pkin(2) * (-m(4) * t414 + t353) + pkin(6) * t441;
t330 = mrSges(2,2) * t425 - mrSges(2,3) * t414 - pkin(5) * t338 - t430 * t331 + t434 * t333;
t1 = [-m(1) * g(1) + t443; -m(1) * g(2) + t446; -m(1) * g(3) + m(2) * t425 + t338; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t446 + t449 * t330 - t426 * t332; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t443 + t426 * t330 + t449 * t332; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - mrSges(2,2) * t415 + t430 * t333 + t434 * t331 - pkin(1) * t353 + pkin(5) * t442 + (pkin(1) * t450 + mrSges(2,1)) * t414;];
tauB = t1;
