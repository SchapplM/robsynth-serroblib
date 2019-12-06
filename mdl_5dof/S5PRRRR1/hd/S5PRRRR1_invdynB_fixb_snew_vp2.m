% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRR1
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 17:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_invdynB_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:03:00
% EndTime: 2019-12-05 17:03:03
% DurationCPUTime: 1.02s
% Computational Cost: add. (7743->222), mult. (14821->281), div. (0->0), fcn. (10535->8), ass. (0->81)
t334 = sin(qJ(4));
t335 = sin(qJ(3));
t338 = cos(qJ(4));
t339 = cos(qJ(3));
t319 = (t334 * t335 - t338 * t339) * qJD(2);
t354 = -m(2) - m(3);
t353 = mrSges(1,3) + mrSges(2,3);
t332 = -g(3) + qJDD(1);
t336 = sin(qJ(2));
t340 = cos(qJ(2));
t326 = -t340 * g(1) + t336 * t332;
t313 = t339 * g(2) - t335 * t326;
t322 = (-mrSges(4,1) * t339 + mrSges(4,2) * t335) * qJD(2);
t348 = qJD(2) * qJD(3);
t323 = t335 * qJDD(2) + t339 * t348;
t349 = qJD(2) * t339;
t328 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t349;
t350 = qJD(2) * t335;
t314 = t335 * g(2) + t339 * t326;
t341 = qJD(2) ^ 2;
t308 = (-t339 ^ 2 * t341 - qJD(3) ^ 2) * pkin(2) + t314;
t343 = (t335 * t339 * t341 + qJDD(3)) * pkin(2) + t313;
t295 = t338 * t308 + t334 * t343;
t347 = t335 * t348;
t324 = t339 * qJDD(2) - t347;
t325 = t336 * g(1) + t340 * t332;
t307 = (-t324 + t347) * pkin(2) - t325;
t333 = sin(qJ(5));
t337 = cos(qJ(5));
t287 = -t333 * t295 + t337 * t307;
t299 = -t319 * qJD(4) + t338 * t323 + t334 * t324;
t320 = (t334 * t339 + t335 * t338) * qJD(2);
t331 = qJD(3) + qJD(4);
t309 = -t333 * t320 + t337 * t331;
t330 = qJDD(3) + qJDD(4);
t290 = t309 * qJD(5) + t337 * t299 + t333 * t330;
t310 = t337 * t320 + t333 * t331;
t296 = -t309 * mrSges(6,1) + t310 * mrSges(6,2);
t298 = -t320 * qJD(4) - t334 * t323 + t338 * t324;
t297 = qJDD(5) - t298;
t315 = qJD(5) + t319;
t300 = -t315 * mrSges(6,2) + t309 * mrSges(6,3);
t285 = m(6) * t287 + t297 * mrSges(6,1) - t290 * mrSges(6,3) - t310 * t296 + t315 * t300;
t288 = t337 * t295 + t333 * t307;
t289 = -t310 * qJD(5) - t333 * t299 + t337 * t330;
t301 = t315 * mrSges(6,1) - t310 * mrSges(6,3);
t286 = m(6) * t288 - t297 * mrSges(6,2) + t289 * mrSges(6,3) + t309 * t296 - t315 * t301;
t305 = t319 * mrSges(5,1) + t320 * mrSges(5,2);
t312 = t331 * mrSges(5,1) - t320 * mrSges(5,3);
t279 = m(5) * t295 - t330 * mrSges(5,2) + t298 * mrSges(5,3) - t333 * t285 + t337 * t286 - t319 * t305 - t331 * t312;
t294 = t334 * t308 - t338 * t343;
t311 = -t331 * mrSges(5,2) - t319 * mrSges(5,3);
t284 = t330 * mrSges(5,1) + t289 * mrSges(6,1) - t290 * mrSges(6,2) - t299 * mrSges(5,3) + t309 * t300 - t310 * t301 - t320 * t305 + t331 * t311 + (-m(5) - m(6)) * t294;
t351 = t334 * t279 + t338 * t284;
t273 = m(4) * t313 + qJDD(3) * mrSges(4,1) - t323 * mrSges(4,3) + qJD(3) * t328 - t322 * t350 + t351;
t327 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t350;
t274 = m(4) * t314 - qJDD(3) * mrSges(4,2) + t324 * mrSges(4,3) - qJD(3) * t327 + t338 * t279 - t334 * t284 + t322 * t349;
t270 = m(3) * t326 - t341 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t335 * t273 + t339 * t274;
t342 = m(5) * t307 - t298 * mrSges(5,1) + t299 * mrSges(5,2) + t337 * t285 + t333 * t286 + t319 * t311 + t320 * t312;
t277 = qJDD(2) * mrSges(3,1) + t324 * mrSges(4,1) - t341 * mrSges(3,2) - t323 * mrSges(4,2) + (m(3) + m(4)) * t325 + (-t327 * t335 + t328 * t339) * qJD(2) - t342;
t352 = t336 * t270 + t340 * t277;
t346 = t340 * t270 - t336 * t277;
t345 = -t339 * t273 - t335 * t274;
t318 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t335 + Ifges(4,4) * t339) * qJD(2);
t317 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t335 + Ifges(4,2) * t339) * qJD(2);
t316 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t335 + Ifges(4,6) * t339) * qJD(2);
t304 = Ifges(5,1) * t320 - Ifges(5,4) * t319 + Ifges(5,5) * t331;
t303 = Ifges(5,4) * t320 - Ifges(5,2) * t319 + Ifges(5,6) * t331;
t302 = Ifges(5,5) * t320 - Ifges(5,6) * t319 + Ifges(5,3) * t331;
t293 = Ifges(6,1) * t310 + Ifges(6,4) * t309 + Ifges(6,5) * t315;
t292 = Ifges(6,4) * t310 + Ifges(6,2) * t309 + Ifges(6,6) * t315;
t291 = Ifges(6,5) * t310 + Ifges(6,6) * t309 + Ifges(6,3) * t315;
t282 = mrSges(6,2) * t294 - mrSges(6,3) * t287 + Ifges(6,1) * t290 + Ifges(6,4) * t289 + Ifges(6,5) * t297 + t309 * t291 - t315 * t292;
t281 = -mrSges(6,1) * t294 + mrSges(6,3) * t288 + Ifges(6,4) * t290 + Ifges(6,2) * t289 + Ifges(6,6) * t297 - t310 * t291 + t315 * t293;
t280 = -mrSges(5,1) * t307 - mrSges(6,1) * t287 + mrSges(6,2) * t288 + mrSges(5,3) * t295 + Ifges(5,4) * t299 - Ifges(6,5) * t290 + Ifges(5,2) * t298 + Ifges(5,6) * t330 - Ifges(6,6) * t289 - Ifges(6,3) * t297 - t310 * t292 + t309 * t293 - t320 * t302 + t331 * t304;
t275 = mrSges(5,2) * t307 + mrSges(5,3) * t294 + Ifges(5,1) * t299 + Ifges(5,4) * t298 + Ifges(5,5) * t330 - t333 * t281 + t337 * t282 - t319 * t302 - t331 * t303;
t272 = -mrSges(4,2) * t325 - mrSges(4,3) * t313 + Ifges(4,1) * t323 + Ifges(4,4) * t324 + Ifges(4,5) * qJDD(3) - qJD(3) * t317 + t338 * t275 - t334 * t280 + t316 * t349;
t271 = Ifges(3,6) * qJDD(2) + t341 * Ifges(3,5) - mrSges(3,1) * g(2) + mrSges(3,3) * t326 - Ifges(4,5) * t323 - Ifges(4,6) * t324 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t313 + mrSges(4,2) * t314 - Ifges(5,5) * t299 - Ifges(5,6) * t298 - Ifges(5,3) * t330 - t320 * t303 - t319 * t304 + mrSges(5,1) * t294 + mrSges(5,2) * t295 - t333 * t282 - t337 * t281 - pkin(2) * t351 + (-t335 * t317 + t339 * t318) * qJD(2);
t267 = mrSges(4,1) * t325 + mrSges(4,3) * t314 + Ifges(4,4) * t323 + Ifges(4,2) * t324 + Ifges(4,6) * qJDD(3) - pkin(2) * t342 + qJD(3) * t318 + t334 * t275 + t338 * t280 - t316 * t350;
t266 = mrSges(3,2) * g(2) - mrSges(3,3) * t325 + Ifges(3,5) * qJDD(2) - t341 * Ifges(3,6) - t335 * t267 + t339 * t272;
t1 = [(-m(1) - m(2)) * g(1) + t346; (-m(1) + t354) * g(2) + t345; -m(1) * g(3) + m(2) * t332 + t352; -mrSges(1,2) * g(3) + mrSges(2,2) * t332 + t340 * t266 - t336 * t271 - qJ(1) * t345 + (-qJ(1) * t354 + t353) * g(2); qJ(1) * t346 + mrSges(1,1) * g(3) - pkin(1) * t352 - mrSges(2,1) * t332 - mrSges(3,1) * t325 + mrSges(3,2) * t326 - t335 * t272 - t339 * t267 - Ifges(3,3) * qJDD(2) + (-qJ(1) * m(2) - t353) * g(1); t336 * t266 + t340 * t271 + pkin(1) * t345 + (-pkin(1) * m(3) - mrSges(1,1) - mrSges(2,1)) * g(2) + (mrSges(1,2) + mrSges(2,2)) * g(1);];
tauB = t1;
