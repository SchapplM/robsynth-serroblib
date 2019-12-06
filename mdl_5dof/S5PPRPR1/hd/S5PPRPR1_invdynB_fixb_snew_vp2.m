% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPRPR1
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPRPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:02
% EndTime: 2019-12-05 15:01:04
% DurationCPUTime: 1.81s
% Computational Cost: add. (18528->184), mult. (35662->238), div. (0->0), fcn. (24486->10), ass. (0->88)
t392 = qJD(3) ^ 2;
t385 = cos(pkin(9));
t415 = pkin(4) * t385;
t382 = sin(pkin(9));
t414 = mrSges(5,2) * t382;
t380 = t385 ^ 2;
t413 = t380 * t392;
t384 = sin(pkin(7));
t387 = cos(pkin(7));
t373 = -t387 * g(1) - t384 * g(2);
t381 = -g(3) + qJDD(1);
t383 = sin(pkin(8));
t386 = cos(pkin(8));
t363 = -t383 * t373 + t386 * t381;
t364 = t386 * t373 + t383 * t381;
t389 = sin(qJ(3));
t391 = cos(qJ(3));
t349 = t389 * t363 + t391 * t364;
t347 = -t392 * pkin(3) + qJDD(3) * qJ(4) + t349;
t372 = t384 * g(1) - t387 * g(2);
t371 = qJDD(2) - t372;
t409 = qJD(3) * qJD(4);
t411 = t385 * t371 - 0.2e1 * t382 * t409;
t340 = (-pkin(6) * qJDD(3) + t392 * t415 - t347) * t382 + t411;
t343 = t382 * t371 + (t347 + 0.2e1 * t409) * t385;
t408 = qJDD(3) * t385;
t341 = -pkin(4) * t413 + pkin(6) * t408 + t343;
t388 = sin(qJ(5));
t390 = cos(qJ(5));
t338 = t390 * t340 - t388 * t341;
t397 = -t382 * t388 + t385 * t390;
t365 = t397 * qJD(3);
t398 = t382 * t390 + t385 * t388;
t366 = t398 * qJD(3);
t354 = -t365 * mrSges(6,1) + t366 * mrSges(6,2);
t357 = t365 * qJD(5) + qJDD(3) * t398;
t361 = -qJD(5) * mrSges(6,2) + t365 * mrSges(6,3);
t336 = m(6) * t338 + qJDD(5) * mrSges(6,1) - t357 * mrSges(6,3) + qJD(5) * t361 - t366 * t354;
t339 = t388 * t340 + t390 * t341;
t356 = -t366 * qJD(5) + qJDD(3) * t397;
t362 = qJD(5) * mrSges(6,1) - t366 * mrSges(6,3);
t337 = m(6) * t339 - qJDD(5) * mrSges(6,2) + t356 * mrSges(6,3) - qJD(5) * t362 + t365 * t354;
t328 = t390 * t336 + t388 * t337;
t342 = -t382 * t347 + t411;
t396 = mrSges(5,3) * qJDD(3) + t392 * (-mrSges(5,1) * t385 + t414);
t326 = m(5) * t342 - t382 * t396 + t328;
t403 = -t388 * t336 + t390 * t337;
t327 = m(5) * t343 + t385 * t396 + t403;
t404 = -t382 * t326 + t385 * t327;
t319 = m(4) * t349 - t392 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t404;
t348 = t391 * t363 - t389 * t364;
t399 = qJDD(4) - t348;
t346 = -qJDD(3) * pkin(3) - t392 * qJ(4) + t399;
t379 = t382 ^ 2;
t344 = (-pkin(3) - t415) * qJDD(3) + (-qJ(4) + (-t379 - t380) * pkin(6)) * t392 + t399;
t394 = m(6) * t344 - t356 * mrSges(6,1) + t357 * mrSges(6,2) - t365 * t361 + t366 * t362;
t393 = -m(5) * t346 + mrSges(5,1) * t408 - t394 + (t379 * t392 + t413) * mrSges(5,3);
t332 = m(4) * t348 - t392 * mrSges(4,2) + (mrSges(4,1) - t414) * qJDD(3) + t393;
t315 = t389 * t319 + t391 * t332;
t312 = m(3) * t363 + t315;
t405 = t391 * t319 - t389 * t332;
t313 = m(3) * t364 + t405;
t406 = -t383 * t312 + t386 * t313;
t306 = m(2) * t373 + t406;
t322 = t385 * t326 + t382 * t327;
t395 = (-m(3) - m(4)) * t371 - t322;
t321 = m(2) * t372 + t395;
t412 = t384 * t306 + t387 * t321;
t307 = t386 * t312 + t383 * t313;
t400 = Ifges(5,5) * t382 + Ifges(5,6) * t385;
t410 = t392 * t400;
t407 = t387 * t306 - t384 * t321;
t402 = Ifges(5,1) * t382 + Ifges(5,4) * t385;
t401 = Ifges(5,4) * t382 + Ifges(5,2) * t385;
t352 = Ifges(6,1) * t366 + Ifges(6,4) * t365 + Ifges(6,5) * qJD(5);
t351 = Ifges(6,4) * t366 + Ifges(6,2) * t365 + Ifges(6,6) * qJD(5);
t350 = Ifges(6,5) * t366 + Ifges(6,6) * t365 + Ifges(6,3) * qJD(5);
t330 = mrSges(6,2) * t344 - mrSges(6,3) * t338 + Ifges(6,1) * t357 + Ifges(6,4) * t356 + Ifges(6,5) * qJDD(5) - qJD(5) * t351 + t365 * t350;
t329 = -mrSges(6,1) * t344 + mrSges(6,3) * t339 + Ifges(6,4) * t357 + Ifges(6,2) * t356 + Ifges(6,6) * qJDD(5) + qJD(5) * t352 - t366 * t350;
t316 = mrSges(5,2) * t346 - mrSges(5,3) * t342 - pkin(6) * t328 + t402 * qJDD(3) - t388 * t329 + t390 * t330 + t385 * t410;
t314 = -mrSges(5,1) * t346 + mrSges(5,3) * t343 - pkin(4) * t394 + pkin(6) * t403 + t401 * qJDD(3) + t390 * t329 + t388 * t330 - t382 * t410;
t308 = -mrSges(4,1) * t371 - mrSges(5,1) * t342 - mrSges(6,1) * t338 + mrSges(5,2) * t343 + mrSges(6,2) * t339 + mrSges(4,3) * t349 - Ifges(6,5) * t357 - Ifges(6,6) * t356 - Ifges(6,3) * qJDD(5) - pkin(3) * t322 - pkin(4) * t328 - t366 * t351 + t365 * t352 + (Ifges(4,6) - t400) * qJDD(3) + (-t382 * t401 + t385 * t402 + Ifges(4,5)) * t392;
t303 = mrSges(4,2) * t371 - mrSges(4,3) * t348 + Ifges(4,5) * qJDD(3) - t392 * Ifges(4,6) - qJ(4) * t322 - t382 * t314 + t385 * t316;
t302 = mrSges(3,2) * t371 - mrSges(3,3) * t363 - pkin(5) * t315 + t391 * t303 - t389 * t308;
t301 = -mrSges(2,1) * t381 - pkin(1) * t307 + mrSges(2,3) * t373 - pkin(2) * t315 - mrSges(3,1) * t363 + mrSges(3,2) * t364 - qJ(4) * t404 - t382 * t316 - t385 * t314 - pkin(3) * (-qJDD(3) * t414 + t393) + mrSges(4,2) * t349 - mrSges(4,1) * t348 - Ifges(4,3) * qJDD(3);
t300 = -mrSges(3,1) * t371 + mrSges(3,3) * t364 + t389 * t303 + t391 * t308 - pkin(2) * (m(4) * t371 + t322) + pkin(5) * t405;
t299 = mrSges(2,2) * t381 - mrSges(2,3) * t372 - qJ(2) * t307 - t383 * t300 + t386 * t302;
t1 = [-m(1) * g(1) + t407; -m(1) * g(2) + t412; -m(1) * g(3) + m(2) * t381 + t307; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t412 + t387 * t299 - t384 * t301; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t407 + t384 * t299 + t387 * t301; -mrSges(1,1) * g(2) + mrSges(2,1) * t372 + mrSges(1,2) * g(1) - mrSges(2,2) * t373 + pkin(1) * t395 + qJ(2) * t406 + t386 * t300 + t383 * t302;];
tauB = t1;
