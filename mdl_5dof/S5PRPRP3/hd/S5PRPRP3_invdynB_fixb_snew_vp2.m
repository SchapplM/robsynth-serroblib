% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPRP3
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:50
% EndTime: 2019-12-05 15:32:54
% DurationCPUTime: 1.37s
% Computational Cost: add. (10331->198), mult. (17904->239), div. (0->0), fcn. (9426->8), ass. (0->86)
t446 = Ifges(5,1) + Ifges(6,1);
t439 = Ifges(5,4) + Ifges(6,4);
t438 = Ifges(5,5) + Ifges(6,5);
t445 = Ifges(5,2) + Ifges(6,2);
t444 = Ifges(5,6) + Ifges(6,6);
t443 = Ifges(5,3) + Ifges(6,3);
t416 = qJD(2) ^ 2;
t409 = sin(pkin(7));
t411 = cos(pkin(7));
t396 = -t411 * g(1) - t409 * g(2);
t407 = -g(3) + qJDD(1);
t413 = sin(qJ(2));
t415 = cos(qJ(2));
t372 = -t413 * t396 + t415 * t407;
t370 = qJDD(2) * pkin(2) + t372;
t373 = t415 * t396 + t413 * t407;
t371 = -t416 * pkin(2) + t373;
t408 = sin(pkin(8));
t410 = cos(pkin(8));
t365 = t410 * t370 - t408 * t371;
t418 = -qJDD(2) * pkin(3) - t365;
t363 = -t416 * pkin(6) + t418;
t412 = sin(qJ(4));
t414 = cos(qJ(4));
t429 = qJD(2) * qJD(4);
t425 = t414 * t429;
t391 = t412 * qJDD(2) + t425;
t392 = t414 * qJDD(2) - t412 * t429;
t430 = qJD(2) * t414;
t401 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t430;
t431 = qJD(2) * t412;
t397 = qJD(4) * pkin(4) - qJ(5) * t431;
t406 = t414 ^ 2;
t359 = t397 * t431 - t392 * pkin(4) + qJDD(5) + (-qJ(5) * t406 - pkin(6)) * t416 + t418;
t400 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t430;
t419 = m(6) * t359 - t392 * mrSges(6,1) - t400 * t430;
t440 = -mrSges(5,2) - mrSges(6,2);
t442 = -m(5) * t363 + t392 * mrSges(5,1) + t440 * t391 + t401 * t430 - t419;
t441 = pkin(4) * t416;
t366 = t408 * t370 + t410 * t371;
t364 = -t416 * pkin(3) + qJDD(2) * pkin(6) + t366;
t395 = t409 * g(1) - t411 * g(2);
t394 = qJDD(3) - t395;
t361 = t414 * t364 + t412 * t394;
t390 = (-mrSges(5,1) * t414 + mrSges(5,2) * t412) * qJD(2);
t428 = qJD(2) * qJD(5);
t358 = t392 * qJ(5) - qJD(4) * t397 - t406 * t441 + 0.2e1 * t414 * t428 + t361;
t389 = (-mrSges(6,1) * t414 + mrSges(6,2) * t412) * qJD(2);
t426 = m(6) * t358 + t392 * mrSges(6,3) + t389 * t430;
t398 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t431;
t432 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t431 - t398;
t353 = m(5) * t361 + t392 * mrSges(5,3) + t432 * qJD(4) + t440 * qJDD(4) + t390 * t430 + t426;
t351 = t414 * t353;
t382 = t414 * t394;
t360 = -t412 * t364 + t382;
t357 = qJDD(4) * pkin(4) + t382 + (-t391 + t425) * qJ(5) + (t414 * t441 - t364 - 0.2e1 * t428) * t412;
t427 = m(6) * t357 + qJDD(4) * mrSges(6,1) + qJD(4) * t400;
t352 = m(5) * t360 + qJDD(4) * mrSges(5,1) + qJD(4) * t401 + (-mrSges(5,3) - mrSges(6,3)) * t391 + (-t389 - t390) * t431 + t427;
t343 = m(4) * t366 - t416 * mrSges(4,1) - qJDD(2) * mrSges(4,2) - t412 * t352 + t351;
t421 = qJD(2) * t432;
t348 = m(4) * t365 + qJDD(2) * mrSges(4,1) - t416 * mrSges(4,2) + t412 * t421 + t442;
t338 = t408 * t343 + t410 * t348;
t336 = m(3) * t372 + qJDD(2) * mrSges(3,1) - t416 * mrSges(3,2) + t338;
t422 = t410 * t343 - t408 * t348;
t337 = m(3) * t373 - t416 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t422;
t423 = -t413 * t336 + t415 * t337;
t330 = m(2) * t396 + t423;
t346 = t414 * t352 + t412 * t353;
t420 = m(4) * t394 + t346;
t345 = (m(2) + m(3)) * t395 - t420;
t436 = t409 * t330 + t411 * t345;
t331 = t415 * t336 + t413 * t337;
t435 = t443 * qJD(4) + (t438 * t412 + t444 * t414) * qJD(2);
t434 = -t444 * qJD(4) + (-t439 * t412 - t445 * t414) * qJD(2);
t433 = t438 * qJD(4) + (t446 * t412 + t439 * t414) * qJD(2);
t424 = t411 * t330 - t409 * t345;
t354 = -t391 * mrSges(6,3) - t389 * t431 + t427;
t340 = mrSges(5,2) * t363 + mrSges(6,2) * t359 - mrSges(5,3) * t360 - mrSges(6,3) * t357 - qJ(5) * t354 + t434 * qJD(4) + t438 * qJDD(4) + t446 * t391 + t439 * t392 + t435 * t430;
t339 = -mrSges(5,1) * t363 + mrSges(5,3) * t361 - mrSges(6,1) * t359 + mrSges(6,3) * t358 - pkin(4) * t419 + qJ(5) * t426 + t445 * t392 + (-pkin(4) * mrSges(6,2) + t439) * t391 + (-qJ(5) * mrSges(6,2) + t444) * qJDD(4) + (-qJ(5) * t398 + t433) * qJD(4) + (-pkin(4) * t398 - t435) * t431;
t332 = -mrSges(4,1) * t394 - mrSges(5,1) * t360 - mrSges(6,1) * t357 + mrSges(5,2) * t361 + mrSges(6,2) * t358 + mrSges(4,3) * t366 + t416 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t346 - pkin(4) * t354 - t444 * t392 - t438 * t391 - t443 * qJDD(4) + (t434 * t412 + t433 * t414) * qJD(2);
t327 = mrSges(4,2) * t394 - mrSges(4,3) * t365 + Ifges(4,5) * qJDD(2) - t416 * Ifges(4,6) - pkin(6) * t346 - t412 * t339 + t414 * t340;
t326 = -mrSges(3,2) * t395 - mrSges(3,3) * t372 + Ifges(3,5) * qJDD(2) - t416 * Ifges(3,6) - qJ(3) * t338 + t410 * t327 - t408 * t332;
t325 = mrSges(3,1) * t395 + mrSges(3,3) * t373 + t416 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t420 + qJ(3) * t422 + t408 * t327 + t410 * t332;
t324 = -pkin(1) * t331 + mrSges(2,3) * t396 - pkin(2) * t338 - mrSges(3,1) * t372 + mrSges(3,2) * t373 - mrSges(4,1) * t365 + mrSges(4,2) * t366 - t414 * t339 - pkin(3) * t442 - pkin(6) * t351 - mrSges(2,1) * t407 + (-pkin(3) * t421 + pkin(6) * t352 - t340) * t412 + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2);
t323 = mrSges(2,2) * t407 - mrSges(2,3) * t395 - pkin(5) * t331 - t413 * t325 + t415 * t326;
t1 = [-m(1) * g(1) + t424; -m(1) * g(2) + t436; -m(1) * g(3) + m(2) * t407 + t331; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t436 + t411 * t323 - t409 * t324; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t424 + t409 * t323 + t411 * t324; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t395 - mrSges(2,2) * t396 + t413 * t326 + t415 * t325 + pkin(1) * (m(3) * t395 - t420) + pkin(5) * t423;];
tauB = t1;
