% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR10_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR10_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR10_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:00
% EndTime: 2019-12-31 18:26:01
% DurationCPUTime: 1.03s
% Computational Cost: add. (14426->182), mult. (18963->217), div. (0->0), fcn. (6806->8), ass. (0->77)
t363 = -m(3) - m(4);
t362 = -pkin(1) - pkin(2);
t361 = -mrSges(2,1) - mrSges(3,1);
t360 = Ifges(3,4) + Ifges(2,5);
t359 = Ifges(2,6) - Ifges(3,6);
t330 = -qJD(1) + qJD(3);
t338 = sin(qJ(5));
t358 = t330 * t338;
t341 = cos(qJ(5));
t357 = t330 * t341;
t340 = sin(qJ(1));
t343 = cos(qJ(1));
t325 = -g(1) * t343 - g(2) * t340;
t344 = qJD(1) ^ 2;
t348 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t325;
t315 = -pkin(1) * t344 + t348;
t309 = t362 * t344 + t348;
t324 = g(1) * t340 - t343 * g(2);
t347 = -qJ(2) * t344 + qJDD(2) - t324;
t311 = t362 * qJDD(1) + t347;
t339 = sin(qJ(3));
t342 = cos(qJ(3));
t304 = -t309 * t339 + t342 * t311;
t329 = -qJDD(1) + qJDD(3);
t302 = pkin(3) * t329 + t304;
t305 = t342 * t309 + t339 * t311;
t328 = t330 ^ 2;
t303 = -pkin(3) * t328 + t305;
t336 = sin(pkin(8));
t337 = cos(pkin(8));
t299 = t336 * t302 + t337 * t303;
t297 = -pkin(4) * t328 + pkin(7) * t329 + t299;
t334 = g(3) + qJDD(4);
t294 = -t297 * t338 + t334 * t341;
t319 = (-mrSges(6,1) * t341 + mrSges(6,2) * t338) * t330;
t355 = qJD(5) * t330;
t320 = t329 * t338 + t341 * t355;
t323 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t357;
t292 = m(6) * t294 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t320 + qJD(5) * t323 - t319 * t358;
t295 = t297 * t341 + t334 * t338;
t321 = t329 * t341 - t338 * t355;
t322 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t358;
t293 = m(6) * t295 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t321 - qJD(5) * t322 + t319 * t357;
t351 = -t292 * t338 + t341 * t293;
t282 = m(5) * t299 - mrSges(5,1) * t328 - mrSges(5,2) * t329 + t351;
t298 = t302 * t337 - t303 * t336;
t296 = -pkin(4) * t329 - pkin(7) * t328 - t298;
t345 = -m(6) * t296 + t321 * mrSges(6,1) - mrSges(6,2) * t320 - t322 * t358 + t323 * t357;
t288 = m(5) * t298 + mrSges(5,1) * t329 - mrSges(5,2) * t328 + t345;
t279 = t336 * t282 + t337 * t288;
t276 = m(4) * t304 + mrSges(4,1) * t329 - mrSges(4,2) * t328 + t279;
t352 = t337 * t282 - t288 * t336;
t277 = m(4) * t305 - mrSges(4,1) * t328 - mrSges(4,2) * t329 + t352;
t353 = -t339 * t276 + t342 * t277;
t349 = m(3) * t315 + qJDD(1) * mrSges(3,3) + t353;
t271 = m(2) * t325 - qJDD(1) * mrSges(2,2) + t361 * t344 + t349;
t273 = t276 * t342 + t277 * t339;
t318 = -qJDD(1) * pkin(1) + t347;
t346 = -m(3) * t318 + qJDD(1) * mrSges(3,1) + t344 * mrSges(3,3) - t273;
t272 = m(2) * t324 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t344 + t346;
t356 = t340 * t271 + t343 * t272;
t284 = t341 * t292 + t338 * t293;
t354 = t343 * t271 - t272 * t340;
t350 = m(5) * t334 + t284;
t314 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t338 + Ifges(6,4) * t341) * t330;
t313 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t338 + Ifges(6,2) * t341) * t330;
t312 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t338 + Ifges(6,6) * t341) * t330;
t286 = mrSges(6,2) * t296 - mrSges(6,3) * t294 + Ifges(6,1) * t320 + Ifges(6,4) * t321 + Ifges(6,5) * qJDD(5) - qJD(5) * t313 + t312 * t357;
t285 = -mrSges(6,1) * t296 + mrSges(6,3) * t295 + Ifges(6,4) * t320 + Ifges(6,2) * t321 + Ifges(6,6) * qJDD(5) + qJD(5) * t314 - t312 * t358;
t283 = t363 * g(3) - t350;
t278 = -mrSges(5,1) * t334 - mrSges(6,1) * t294 + mrSges(6,2) * t295 + mrSges(5,3) * t299 + Ifges(5,5) * t328 - Ifges(6,5) * t320 + Ifges(5,6) * t329 - Ifges(6,6) * t321 - Ifges(6,3) * qJDD(5) - pkin(4) * t284 + (-t313 * t338 + t314 * t341) * t330;
t274 = mrSges(5,2) * t334 - mrSges(5,3) * t298 + Ifges(5,5) * t329 - Ifges(5,6) * t328 - pkin(7) * t284 - t285 * t338 + t286 * t341;
t267 = mrSges(4,2) * g(3) - mrSges(4,3) * t304 + Ifges(4,5) * t329 - Ifges(4,6) * t328 - qJ(4) * t279 + t274 * t337 - t278 * t336;
t266 = -mrSges(4,1) * g(3) + mrSges(4,3) * t305 + t328 * Ifges(4,5) + Ifges(4,6) * t329 - pkin(3) * t350 + qJ(4) * t352 + t336 * t274 + t337 * t278;
t265 = mrSges(3,2) * t318 - mrSges(2,3) * t324 - pkin(6) * t273 - qJ(2) * t283 - t266 * t339 + t267 * t342 - t359 * t344 + t360 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t264 = mrSges(2,3) * t325 + mrSges(3,2) * t315 - t339 * t267 - t342 * t266 + pkin(2) * t350 - pkin(6) * t353 - pkin(1) * t283 + t360 * t344 + t359 * qJDD(1) + (pkin(2) * m(4) - t361) * g(3);
t1 = [-m(1) * g(1) + t354; -m(1) * g(2) + t356; (-m(1) - m(2) + t363) * g(3) - t350; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t356 - t340 * t264 + t343 * t265; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t354 + t343 * t264 + t340 * t265; qJ(2) * (-mrSges(3,1) * t344 + t349) + pkin(1) * t346 + mrSges(2,1) * t324 - mrSges(2,2) * t325 - pkin(2) * t273 + mrSges(3,3) * t315 - mrSges(3,1) * t318 - pkin(3) * t279 - mrSges(4,1) * t304 + mrSges(4,2) * t305 - t338 * t286 - t341 * t285 - pkin(4) * t345 - pkin(7) * t351 - mrSges(5,1) * t298 + mrSges(5,2) * t299 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-Ifges(4,3) - Ifges(5,3)) * t329 + (Ifges(3,2) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
