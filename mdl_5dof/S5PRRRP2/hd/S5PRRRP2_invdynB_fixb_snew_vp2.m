% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:47
% EndTime: 2019-12-05 16:41:50
% DurationCPUTime: 1.52s
% Computational Cost: add. (16556->193), mult. (22100->240), div. (0->0), fcn. (12352->8), ass. (0->83)
t437 = Ifges(5,1) + Ifges(6,1);
t433 = Ifges(5,4) - Ifges(6,5);
t432 = Ifges(5,5) + Ifges(6,4);
t436 = Ifges(5,2) + Ifges(6,3);
t431 = Ifges(5,6) - Ifges(6,6);
t435 = Ifges(5,3) + Ifges(6,2);
t434 = mrSges(5,3) + mrSges(6,2);
t400 = qJD(2) + qJD(3);
t406 = sin(qJ(4));
t430 = t400 * t406;
t409 = cos(qJ(4));
t429 = t400 * t409;
t403 = -g(3) + qJDD(1);
t428 = t409 * t403;
t404 = sin(pkin(8));
t405 = cos(pkin(8));
t392 = t404 * g(1) - t405 * g(2);
t393 = -t405 * g(1) - t404 * g(2);
t408 = sin(qJ(2));
t411 = cos(qJ(2));
t364 = t411 * t392 - t408 * t393;
t362 = qJDD(2) * pkin(2) + t364;
t365 = t408 * t392 + t411 * t393;
t413 = qJD(2) ^ 2;
t363 = -t413 * pkin(2) + t365;
t407 = sin(qJ(3));
t410 = cos(qJ(3));
t358 = t407 * t362 + t410 * t363;
t398 = t400 ^ 2;
t399 = qJDD(2) + qJDD(3);
t356 = -t398 * pkin(3) + t399 * pkin(7) + t358;
t353 = t409 * t356 + t406 * t403;
t380 = (-mrSges(5,1) * t409 + mrSges(5,2) * t406) * t400;
t423 = qJD(4) * t400;
t382 = t409 * t399 - t406 * t423;
t388 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t430;
t378 = (-pkin(4) * t409 - qJ(5) * t406) * t400;
t412 = qJD(4) ^ 2;
t350 = -t412 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t378 * t429 + t353;
t379 = (-mrSges(6,1) * t409 - mrSges(6,3) * t406) * t400;
t389 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t430;
t417 = m(6) * t350 + qJDD(4) * mrSges(6,3) + qJD(4) * t389 + t379 * t429;
t345 = m(5) * t353 - qJDD(4) * mrSges(5,2) - qJD(4) * t388 + t380 * t429 + t434 * t382 + t417;
t352 = -t406 * t356 + t428;
t381 = t406 * t399 + t409 * t423;
t390 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t429;
t351 = -qJDD(4) * pkin(4) - t412 * qJ(5) - t428 + qJDD(5) + (t378 * t400 + t356) * t406;
t391 = mrSges(6,2) * t429 + qJD(4) * mrSges(6,3);
t415 = -m(6) * t351 + qJDD(4) * mrSges(6,1) + qJD(4) * t391;
t346 = m(5) * t352 + qJDD(4) * mrSges(5,1) + qJD(4) * t390 + (-t379 - t380) * t430 - t434 * t381 + t415;
t418 = t409 * t345 - t406 * t346;
t338 = m(4) * t358 - t398 * mrSges(4,1) - t399 * mrSges(4,2) + t418;
t357 = t410 * t362 - t407 * t363;
t355 = -t399 * pkin(3) - t398 * pkin(7) - t357;
t348 = -t382 * pkin(4) - t381 * qJ(5) + (-0.2e1 * qJD(5) * t406 + (pkin(4) * t406 - qJ(5) * t409) * qJD(4)) * t400 + t355;
t347 = m(6) * t348 - t382 * mrSges(6,1) - t381 * mrSges(6,3) - t389 * t430 - t391 * t429;
t414 = -m(5) * t355 + t382 * mrSges(5,1) - t381 * mrSges(5,2) - t388 * t430 + t390 * t429 - t347;
t341 = m(4) * t357 + t399 * mrSges(4,1) - t398 * mrSges(4,2) + t414;
t333 = t407 * t338 + t410 * t341;
t331 = m(3) * t364 + qJDD(2) * mrSges(3,1) - t413 * mrSges(3,2) + t333;
t419 = t410 * t338 - t407 * t341;
t332 = m(3) * t365 - t413 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t419;
t325 = t411 * t331 + t408 * t332;
t323 = m(2) * t392 + t325;
t420 = -t408 * t331 + t411 * t332;
t324 = m(2) * t393 + t420;
t427 = t405 * t323 + t404 * t324;
t339 = t406 * t345 + t409 * t346;
t426 = (-t433 * t406 - t436 * t409) * t400 - t431 * qJD(4);
t425 = (t432 * t406 + t431 * t409) * t400 + t435 * qJD(4);
t424 = (t437 * t406 + t433 * t409) * t400 + t432 * qJD(4);
t422 = m(4) * t403 + t339;
t421 = -t404 * t323 + t405 * t324;
t416 = m(3) * t403 + t422;
t335 = mrSges(5,2) * t355 + mrSges(6,2) * t351 - mrSges(5,3) * t352 - mrSges(6,3) * t348 - qJ(5) * t347 + t426 * qJD(4) + t432 * qJDD(4) + t437 * t381 + t433 * t382 + t425 * t429;
t334 = -mrSges(5,1) * t355 - mrSges(6,1) * t348 + mrSges(6,2) * t350 + mrSges(5,3) * t353 - pkin(4) * t347 + t424 * qJD(4) + t431 * qJDD(4) + t433 * t381 + t436 * t382 - t425 * t430;
t327 = Ifges(4,6) * t399 + t398 * Ifges(4,5) - mrSges(4,1) * t403 + mrSges(4,3) * t358 - mrSges(5,1) * t352 + mrSges(5,2) * t353 + mrSges(6,1) * t351 - mrSges(6,3) * t350 - pkin(4) * t415 - qJ(5) * t417 - pkin(3) * t339 + (-qJ(5) * mrSges(6,2) - t431) * t382 + (pkin(4) * mrSges(6,2) - t432) * t381 - t435 * qJDD(4) + (t424 * t409 + (pkin(4) * t379 + t426) * t406) * t400;
t326 = mrSges(4,2) * t403 - mrSges(4,3) * t357 + Ifges(4,5) * t399 - t398 * Ifges(4,6) - pkin(7) * t339 - t406 * t334 + t409 * t335;
t319 = mrSges(3,2) * t403 - mrSges(3,3) * t364 + Ifges(3,5) * qJDD(2) - t413 * Ifges(3,6) - pkin(6) * t333 + t410 * t326 - t407 * t327;
t318 = -mrSges(3,1) * t403 + mrSges(3,3) * t365 + t413 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t422 + pkin(6) * t419 + t407 * t326 + t410 * t327;
t317 = mrSges(2,2) * t403 - mrSges(2,3) * t392 - pkin(5) * t325 - t408 * t318 + t411 * t319;
t316 = -mrSges(2,1) * t403 + mrSges(2,3) * t393 - pkin(1) * t416 + pkin(5) * t420 + t411 * t318 + t408 * t319;
t1 = [-m(1) * g(1) + t421; -m(1) * g(2) + t427; -m(1) * g(3) + m(2) * t403 + t416; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t427 - t404 * t316 + t405 * t317; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t421 + t405 * t316 + t404 * t317; pkin(1) * t325 - mrSges(2,2) * t393 + mrSges(2,1) * t392 + pkin(2) * t333 + mrSges(3,1) * t364 - mrSges(3,2) * t365 + t409 * t334 + pkin(3) * t414 + pkin(7) * t418 + mrSges(4,1) * t357 - mrSges(4,2) * t358 + t406 * t335 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(4,3) * t399 + Ifges(3,3) * qJDD(2);];
tauB = t1;
