% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:27
% EndTime: 2019-12-31 17:46:28
% DurationCPUTime: 1.24s
% Computational Cost: add. (13204->201), mult. (24984->243), div. (0->0), fcn. (11933->8), ass. (0->90)
t405 = qJD(1) ^ 2;
t435 = -pkin(1) - pkin(2);
t434 = mrSges(2,1) + mrSges(3,1);
t433 = Ifges(3,4) + Ifges(2,5);
t432 = Ifges(2,6) - Ifges(3,6);
t399 = cos(pkin(8));
t431 = mrSges(5,1) * t399;
t397 = sin(pkin(8));
t430 = mrSges(5,2) * t397;
t390 = t399 ^ 2;
t429 = t390 * t405;
t402 = sin(qJ(1));
t404 = cos(qJ(1));
t379 = -t404 * g(1) - t402 * g(2);
t411 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t379;
t372 = -t405 * pkin(1) + t411;
t367 = t435 * t405 + t411;
t378 = t402 * g(1) - t404 * g(2);
t410 = -t405 * qJ(2) + qJDD(2) - t378;
t371 = t435 * qJDD(1) + t410;
t398 = sin(pkin(7));
t400 = cos(pkin(7));
t355 = t400 * t367 + t398 * t371;
t353 = -t405 * pkin(3) - qJDD(1) * qJ(4) + t355;
t394 = g(3) + qJDD(3);
t425 = qJD(1) * qJD(4);
t427 = t399 * t394 + 0.2e1 * t397 * t425;
t346 = (pkin(4) * t399 * t405 + pkin(6) * qJDD(1) - t353) * t397 + t427;
t350 = t397 * t394 + (t353 - (2 * t425)) * t399;
t424 = qJDD(1) * t399;
t347 = -pkin(4) * t429 - pkin(6) * t424 + t350;
t401 = sin(qJ(5));
t403 = cos(qJ(5));
t344 = t403 * t346 - t401 * t347;
t415 = t397 * t401 - t399 * t403;
t374 = t415 * qJD(1);
t416 = -t397 * t403 - t399 * t401;
t375 = t416 * qJD(1);
t360 = -t374 * mrSges(6,1) + t375 * mrSges(6,2);
t363 = t374 * qJD(5) + t416 * qJDD(1);
t368 = -qJD(5) * mrSges(6,2) + t374 * mrSges(6,3);
t342 = m(6) * t344 + qJDD(5) * mrSges(6,1) - t363 * mrSges(6,3) + qJD(5) * t368 - t375 * t360;
t345 = t401 * t346 + t403 * t347;
t362 = -t375 * qJD(5) + t415 * qJDD(1);
t369 = qJD(5) * mrSges(6,1) - t375 * mrSges(6,3);
t343 = m(6) * t345 - qJDD(5) * mrSges(6,2) + t362 * mrSges(6,3) - qJD(5) * t369 + t374 * t360;
t335 = t403 * t342 + t401 * t343;
t349 = -t397 * t353 + t427;
t414 = mrSges(5,3) * qJDD(1) + t405 * (-t430 + t431);
t333 = m(5) * t349 + t414 * t397 + t335;
t420 = -t401 * t342 + t403 * t343;
t334 = m(5) * t350 - t414 * t399 + t420;
t421 = -t397 * t333 + t399 * t334;
t329 = m(4) * t355 - t405 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t421;
t354 = -t398 * t367 + t400 * t371;
t412 = qJDD(1) * pkin(3) + qJDD(4) - t354;
t352 = -t405 * qJ(4) + t412;
t389 = t397 ^ 2;
t348 = pkin(4) * t424 + (-qJ(4) + (-t389 - t390) * pkin(6)) * t405 + t412;
t408 = m(6) * t348 - t362 * mrSges(6,1) + t363 * mrSges(6,2) - t374 * t368 + t375 * t369;
t406 = -m(5) * t352 + qJDD(1) * t430 - t408 + (t389 * t405 + t429) * mrSges(5,3);
t338 = m(4) * t354 - t405 * mrSges(4,2) + (-mrSges(4,1) - t431) * qJDD(1) + t406;
t422 = t400 * t329 - t398 * t338;
t413 = m(3) * t372 + qJDD(1) * mrSges(3,3) + t422;
t323 = m(2) * t379 - qJDD(1) * mrSges(2,2) - t434 * t405 + t413;
t326 = t398 * t329 + t400 * t338;
t373 = -qJDD(1) * pkin(1) + t410;
t407 = -m(3) * t373 + qJDD(1) * mrSges(3,1) + t405 * mrSges(3,3) - t326;
t324 = m(2) * t378 + qJDD(1) * mrSges(2,1) - t405 * mrSges(2,2) + t407;
t428 = t402 * t323 + t404 * t324;
t417 = -Ifges(5,5) * t397 - Ifges(5,6) * t399;
t426 = t405 * t417;
t423 = t404 * t323 - t402 * t324;
t419 = -Ifges(5,1) * t397 - Ifges(5,4) * t399;
t418 = -Ifges(5,4) * t397 - Ifges(5,2) * t399;
t331 = t399 * t333 + t397 * t334;
t409 = -m(4) * t394 - t331;
t358 = Ifges(6,1) * t375 + Ifges(6,4) * t374 + Ifges(6,5) * qJD(5);
t357 = Ifges(6,4) * t375 + Ifges(6,2) * t374 + Ifges(6,6) * qJD(5);
t356 = Ifges(6,5) * t375 + Ifges(6,6) * t374 + Ifges(6,3) * qJD(5);
t337 = mrSges(6,2) * t348 - mrSges(6,3) * t344 + Ifges(6,1) * t363 + Ifges(6,4) * t362 + Ifges(6,5) * qJDD(5) - qJD(5) * t357 + t374 * t356;
t336 = -mrSges(6,1) * t348 + mrSges(6,3) * t345 + Ifges(6,4) * t363 + Ifges(6,2) * t362 + Ifges(6,6) * qJDD(5) + qJD(5) * t358 - t375 * t356;
t330 = -m(3) * g(3) + t409;
t327 = mrSges(5,2) * t352 - mrSges(5,3) * t349 - pkin(6) * t335 + t419 * qJDD(1) - t401 * t336 + t403 * t337 - t399 * t426;
t325 = -mrSges(5,1) * t352 + mrSges(5,3) * t350 - pkin(4) * t408 + pkin(6) * t420 + t418 * qJDD(1) + t403 * t336 + t401 * t337 + t397 * t426;
t319 = -mrSges(4,1) * t394 - mrSges(5,1) * t349 - mrSges(6,1) * t344 + mrSges(5,2) * t350 + mrSges(6,2) * t345 + mrSges(4,3) * t355 - Ifges(6,5) * t363 - Ifges(6,6) * t362 - Ifges(6,3) * qJDD(5) - pkin(3) * t331 - pkin(4) * t335 - t375 * t357 + t374 * t358 + (-Ifges(4,6) - t417) * qJDD(1) + (t397 * t418 - t399 * t419 + Ifges(4,5)) * t405;
t318 = mrSges(4,2) * t394 - mrSges(4,3) * t354 - Ifges(4,5) * qJDD(1) - t405 * Ifges(4,6) - qJ(4) * t331 - t397 * t325 + t399 * t327;
t317 = mrSges(3,2) * t373 - mrSges(2,3) * t378 - qJ(2) * t330 - qJ(3) * t326 + t400 * t318 - t398 * t319 - t432 * t405 + t433 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t316 = mrSges(3,2) * t372 + mrSges(2,3) * t379 - pkin(1) * t330 - pkin(2) * t409 + t434 * g(3) - qJ(3) * t422 + t432 * qJDD(1) - t398 * t318 - t400 * t319 + t433 * t405;
t1 = [-m(1) * g(1) + t423; -m(1) * g(2) + t428; (-m(1) - m(2) - m(3)) * g(3) + t409; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t428 - t402 * t316 + t404 * t317; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t423 + t404 * t316 + t402 * t317; pkin(1) * t407 + qJ(2) * (-t405 * mrSges(3,1) + t413) + mrSges(2,1) * t378 - mrSges(2,2) * t379 - pkin(2) * t326 - mrSges(3,1) * t373 + mrSges(3,3) * t372 - pkin(3) * t406 - qJ(4) * t421 - t397 * t327 - t399 * t325 - mrSges(4,1) * t354 + mrSges(4,2) * t355 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (pkin(3) * t431 + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1);];
tauB = t1;
