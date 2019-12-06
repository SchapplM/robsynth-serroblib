% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPRR2
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:44:56
% EndTime: 2019-12-05 15:44:58
% DurationCPUTime: 1.55s
% Computational Cost: add. (20937->174), mult. (27758->222), div. (0->0), fcn. (16760->10), ass. (0->77)
t342 = qJD(2) + qJD(4);
t348 = sin(qJ(5));
t366 = t342 * t348;
t351 = cos(qJ(5));
t365 = t342 * t351;
t345 = sin(pkin(8));
t347 = cos(pkin(8));
t336 = -t347 * g(1) - t345 * g(2);
t343 = -g(3) + qJDD(1);
t350 = sin(qJ(2));
t353 = cos(qJ(2));
t320 = -t350 * t336 + t353 * t343;
t318 = qJDD(2) * pkin(2) + t320;
t321 = t353 * t336 + t350 * t343;
t354 = qJD(2) ^ 2;
t319 = -t354 * pkin(2) + t321;
t344 = sin(pkin(9));
t346 = cos(pkin(9));
t313 = t346 * t318 - t344 * t319;
t311 = qJDD(2) * pkin(3) + t313;
t314 = t344 * t318 + t346 * t319;
t312 = -t354 * pkin(3) + t314;
t349 = sin(qJ(4));
t352 = cos(qJ(4));
t308 = t349 * t311 + t352 * t312;
t340 = t342 ^ 2;
t341 = qJDD(2) + qJDD(4);
t306 = -t340 * pkin(4) + t341 * pkin(7) + t308;
t335 = t345 * g(1) - t347 * g(2);
t334 = qJDD(3) - t335;
t303 = -t348 * t306 + t351 * t334;
t327 = (-mrSges(6,1) * t351 + mrSges(6,2) * t348) * t342;
t363 = qJD(5) * t342;
t328 = t348 * t341 + t351 * t363;
t333 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t365;
t301 = m(6) * t303 + qJDD(5) * mrSges(6,1) - t328 * mrSges(6,3) + qJD(5) * t333 - t327 * t366;
t304 = t351 * t306 + t348 * t334;
t329 = t351 * t341 - t348 * t363;
t332 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t366;
t302 = m(6) * t304 - qJDD(5) * mrSges(6,2) + t329 * mrSges(6,3) - qJD(5) * t332 + t327 * t365;
t357 = -t348 * t301 + t351 * t302;
t290 = m(5) * t308 - t340 * mrSges(5,1) - t341 * mrSges(5,2) + t357;
t307 = t352 * t311 - t349 * t312;
t305 = -t341 * pkin(4) - t340 * pkin(7) - t307;
t355 = -m(6) * t305 + t329 * mrSges(6,1) - t328 * mrSges(6,2) - t332 * t366 + t333 * t365;
t297 = m(5) * t307 + t341 * mrSges(5,1) - t340 * mrSges(5,2) + t355;
t287 = t349 * t290 + t352 * t297;
t285 = m(4) * t313 + qJDD(2) * mrSges(4,1) - t354 * mrSges(4,2) + t287;
t358 = t352 * t290 - t349 * t297;
t286 = m(4) * t314 - t354 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t358;
t279 = t346 * t285 + t344 * t286;
t277 = m(3) * t320 + qJDD(2) * mrSges(3,1) - t354 * mrSges(3,2) + t279;
t359 = -t344 * t285 + t346 * t286;
t278 = m(3) * t321 - t354 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t359;
t360 = -t350 * t277 + t353 * t278;
t270 = m(2) * t336 + t360;
t293 = t351 * t301 + t348 * t302;
t362 = m(5) * t334 + t293;
t356 = m(4) * t334 + t362;
t292 = (m(2) + m(3)) * t335 - t356;
t364 = t345 * t270 + t347 * t292;
t271 = t353 * t277 + t350 * t278;
t361 = t347 * t270 - t345 * t292;
t324 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t348 + Ifges(6,4) * t351) * t342;
t323 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t348 + Ifges(6,2) * t351) * t342;
t322 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t348 + Ifges(6,6) * t351) * t342;
t295 = mrSges(6,2) * t305 - mrSges(6,3) * t303 + Ifges(6,1) * t328 + Ifges(6,4) * t329 + Ifges(6,5) * qJDD(5) - qJD(5) * t323 + t322 * t365;
t294 = -mrSges(6,1) * t305 + mrSges(6,3) * t304 + Ifges(6,4) * t328 + Ifges(6,2) * t329 + Ifges(6,6) * qJDD(5) + qJD(5) * t324 - t322 * t366;
t284 = -mrSges(5,1) * t334 - mrSges(6,1) * t303 + mrSges(6,2) * t304 + mrSges(5,3) * t308 + t340 * Ifges(5,5) - Ifges(6,5) * t328 + Ifges(5,6) * t341 - Ifges(6,6) * t329 - Ifges(6,3) * qJDD(5) - pkin(4) * t293 + (-t323 * t348 + t324 * t351) * t342;
t280 = mrSges(5,2) * t334 - mrSges(5,3) * t307 + Ifges(5,5) * t341 - t340 * Ifges(5,6) - pkin(7) * t293 - t348 * t294 + t351 * t295;
t273 = mrSges(4,2) * t334 - mrSges(4,3) * t313 + Ifges(4,5) * qJDD(2) - t354 * Ifges(4,6) - pkin(6) * t287 + t352 * t280 - t349 * t284;
t272 = -mrSges(4,1) * t334 + mrSges(4,3) * t314 + t354 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t362 + pkin(6) * t358 + t349 * t280 + t352 * t284;
t267 = -pkin(1) * t271 - pkin(2) * t279 + mrSges(3,2) * t321 - mrSges(3,1) * t320 - pkin(3) * t287 - mrSges(4,1) * t313 + mrSges(4,2) * t314 - mrSges(5,1) * t307 + mrSges(5,2) * t308 - t348 * t295 - t351 * t294 - pkin(4) * t355 - pkin(7) * t357 + mrSges(2,3) * t336 - Ifges(5,3) * t341 - mrSges(2,1) * t343 + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2);
t266 = -mrSges(3,2) * t335 - mrSges(3,3) * t320 + Ifges(3,5) * qJDD(2) - t354 * Ifges(3,6) - qJ(3) * t279 - t344 * t272 + t346 * t273;
t265 = mrSges(3,1) * t335 + mrSges(3,3) * t321 + t354 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t356 + qJ(3) * t359 + t346 * t272 + t344 * t273;
t264 = mrSges(2,2) * t343 - mrSges(2,3) * t335 - pkin(5) * t271 - t350 * t265 + t353 * t266;
t1 = [-m(1) * g(1) + t361; -m(1) * g(2) + t364; -m(1) * g(3) + m(2) * t343 + t271; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t364 + t347 * t264 - t345 * t267; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t361 + t345 * t264 + t347 * t267; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t335 - mrSges(2,2) * t336 + t350 * t266 + t353 * t265 + pkin(1) * (m(3) * t335 - t356) + pkin(5) * t360;];
tauB = t1;
