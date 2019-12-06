% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:30
% EndTime: 2019-12-05 17:53:32
% DurationCPUTime: 1.12s
% Computational Cost: add. (8911->218), mult. (19231->285), div. (0->0), fcn. (12248->10), ass. (0->90)
t375 = sin(qJ(1));
t378 = cos(qJ(1));
t392 = t378 * g(2) + t375 * g(3);
t351 = qJDD(1) * pkin(1) + t392;
t379 = qJD(1) ^ 2;
t387 = t375 * g(2) - t378 * g(3);
t353 = -t379 * pkin(1) + t387;
t370 = sin(pkin(8));
t372 = cos(pkin(8));
t331 = t370 * t351 + t372 * t353;
t326 = -t379 * pkin(2) + qJDD(1) * pkin(6) + t331;
t368 = -g(1) + qJDD(2);
t374 = sin(qJ(3));
t377 = cos(qJ(3));
t315 = -t374 * t326 + t377 * t368;
t389 = qJD(1) * qJD(3);
t388 = t377 * t389;
t354 = t374 * qJDD(1) + t388;
t312 = (-t354 + t388) * qJ(4) + (t374 * t377 * t379 + qJDD(3)) * pkin(3) + t315;
t316 = t377 * t326 + t374 * t368;
t355 = t377 * qJDD(1) - t374 * t389;
t390 = t374 * qJD(1);
t356 = qJD(3) * pkin(3) - qJ(4) * t390;
t367 = t377 ^ 2;
t313 = -t367 * t379 * pkin(3) + t355 * qJ(4) - qJD(3) * t356 + t316;
t369 = sin(pkin(9));
t371 = cos(pkin(9));
t341 = (t377 * t369 + t374 * t371) * qJD(1);
t292 = -0.2e1 * qJD(4) * t341 + t371 * t312 - t369 * t313;
t333 = t371 * t354 + t369 * t355;
t340 = (-t374 * t369 + t377 * t371) * qJD(1);
t290 = (qJD(3) * t340 - t333) * pkin(7) + (t340 * t341 + qJDD(3)) * pkin(4) + t292;
t293 = 0.2e1 * qJD(4) * t340 + t369 * t312 + t371 * t313;
t332 = -t369 * t354 + t371 * t355;
t336 = qJD(3) * pkin(4) - t341 * pkin(7);
t339 = t340 ^ 2;
t291 = -t339 * pkin(4) + t332 * pkin(7) - qJD(3) * t336 + t293;
t373 = sin(qJ(5));
t376 = cos(qJ(5));
t288 = t376 * t290 - t373 * t291;
t323 = t376 * t340 - t373 * t341;
t302 = t323 * qJD(5) + t373 * t332 + t376 * t333;
t324 = t373 * t340 + t376 * t341;
t311 = -t323 * mrSges(6,1) + t324 * mrSges(6,2);
t366 = qJD(3) + qJD(5);
t317 = -t366 * mrSges(6,2) + t323 * mrSges(6,3);
t365 = qJDD(3) + qJDD(5);
t284 = m(6) * t288 + t365 * mrSges(6,1) - t302 * mrSges(6,3) - t324 * t311 + t366 * t317;
t289 = t373 * t290 + t376 * t291;
t301 = -t324 * qJD(5) + t376 * t332 - t373 * t333;
t318 = t366 * mrSges(6,1) - t324 * mrSges(6,3);
t285 = m(6) * t289 - t365 * mrSges(6,2) + t301 * mrSges(6,3) + t323 * t311 - t366 * t318;
t278 = t376 * t284 + t373 * t285;
t328 = -t340 * mrSges(5,1) + t341 * mrSges(5,2);
t334 = -qJD(3) * mrSges(5,2) + t340 * mrSges(5,3);
t276 = m(5) * t292 + qJDD(3) * mrSges(5,1) - t333 * mrSges(5,3) + qJD(3) * t334 - t341 * t328 + t278;
t335 = qJD(3) * mrSges(5,1) - t341 * mrSges(5,3);
t384 = -t373 * t284 + t376 * t285;
t277 = m(5) * t293 - qJDD(3) * mrSges(5,2) + t332 * mrSges(5,3) - qJD(3) * t335 + t340 * t328 + t384;
t272 = t371 * t276 + t369 * t277;
t391 = qJD(1) * t377;
t352 = (-t377 * mrSges(4,1) + t374 * mrSges(4,2)) * qJD(1);
t358 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t391;
t270 = m(4) * t315 + qJDD(3) * mrSges(4,1) - t354 * mrSges(4,3) + qJD(3) * t358 - t352 * t390 + t272;
t357 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t390;
t385 = -t369 * t276 + t371 * t277;
t271 = m(4) * t316 - qJDD(3) * mrSges(4,2) + t355 * mrSges(4,3) - qJD(3) * t357 + t352 * t391 + t385;
t386 = -t374 * t270 + t377 * t271;
t330 = t372 * t351 - t370 * t353;
t382 = -qJDD(1) * pkin(2) - t330;
t314 = -t355 * pkin(3) + qJDD(4) + t356 * t390 + (-qJ(4) * t367 - pkin(6)) * t379 + t382;
t295 = -t332 * pkin(4) - t339 * pkin(7) + t341 * t336 + t314;
t383 = m(6) * t295 - t301 * mrSges(6,1) + t302 * mrSges(6,2) - t323 * t317 + t324 * t318;
t304 = Ifges(6,4) * t324 + Ifges(6,2) * t323 + Ifges(6,6) * t366;
t305 = Ifges(6,1) * t324 + Ifges(6,4) * t323 + Ifges(6,5) * t366;
t381 = mrSges(6,1) * t288 - mrSges(6,2) * t289 + Ifges(6,5) * t302 + Ifges(6,6) * t301 + Ifges(6,3) * t365 + t324 * t304 - t323 * t305;
t286 = m(5) * t314 - t332 * mrSges(5,1) + t333 * mrSges(5,2) - t340 * t334 + t341 * t335 + t383;
t325 = -t379 * pkin(6) + t382;
t380 = -m(4) * t325 + t355 * mrSges(4,1) - t354 * mrSges(4,2) - t357 * t390 + t358 * t391 - t286;
t347 = Ifges(4,5) * qJD(3) + (t374 * Ifges(4,1) + t377 * Ifges(4,4)) * qJD(1);
t346 = Ifges(4,6) * qJD(3) + (t374 * Ifges(4,4) + t377 * Ifges(4,2)) * qJD(1);
t322 = Ifges(5,1) * t341 + Ifges(5,4) * t340 + Ifges(5,5) * qJD(3);
t321 = Ifges(5,4) * t341 + Ifges(5,2) * t340 + Ifges(5,6) * qJD(3);
t320 = Ifges(5,5) * t341 + Ifges(5,6) * t340 + Ifges(5,3) * qJD(3);
t303 = Ifges(6,5) * t324 + Ifges(6,6) * t323 + Ifges(6,3) * t366;
t280 = mrSges(6,2) * t295 - mrSges(6,3) * t288 + Ifges(6,1) * t302 + Ifges(6,4) * t301 + Ifges(6,5) * t365 + t323 * t303 - t366 * t304;
t279 = -mrSges(6,1) * t295 + mrSges(6,3) * t289 + Ifges(6,4) * t302 + Ifges(6,2) * t301 + Ifges(6,6) * t365 - t324 * t303 + t366 * t305;
t268 = mrSges(5,2) * t314 - mrSges(5,3) * t292 + Ifges(5,1) * t333 + Ifges(5,4) * t332 + Ifges(5,5) * qJDD(3) - pkin(7) * t278 - qJD(3) * t321 - t373 * t279 + t376 * t280 + t340 * t320;
t267 = -mrSges(5,1) * t314 + mrSges(5,3) * t293 + Ifges(5,4) * t333 + Ifges(5,2) * t332 + Ifges(5,6) * qJDD(3) - pkin(4) * t383 + pkin(7) * t384 + qJD(3) * t322 + t376 * t279 + t373 * t280 - t341 * t320;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t392 - mrSges(2,2) * t387 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t330 - mrSges(3,2) * t331 + t374 * (mrSges(4,2) * t325 - mrSges(4,3) * t315 + Ifges(4,1) * t354 + Ifges(4,4) * t355 + Ifges(4,5) * qJDD(3) - qJ(4) * t272 - qJD(3) * t346 - t369 * t267 + t371 * t268) + t377 * (-mrSges(4,1) * t325 + mrSges(4,3) * t316 + Ifges(4,4) * t354 + Ifges(4,2) * t355 + Ifges(4,6) * qJDD(3) - pkin(3) * t286 + qJ(4) * t385 + qJD(3) * t347 + t371 * t267 + t369 * t268) + pkin(2) * t380 + pkin(6) * t386 + pkin(1) * (t370 * (m(3) * t331 - t379 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t386) + t372 * (m(3) * t330 + qJDD(1) * mrSges(3,1) - t379 * mrSges(3,2) + t380)); m(3) * t368 + t377 * t270 + t374 * t271; (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t374 * t346 - t377 * t347) * qJD(1) + t381 + Ifges(4,5) * t354 + Ifges(4,6) * t355 - t340 * t322 + t341 * t321 + Ifges(5,6) * t332 + Ifges(5,5) * t333 + mrSges(4,1) * t315 - mrSges(4,2) * t316 + mrSges(5,1) * t292 - mrSges(5,2) * t293 + pkin(4) * t278 + pkin(3) * t272; t286; t381;];
tauJ = t1;
