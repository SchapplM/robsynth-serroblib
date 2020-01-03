% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:33
% EndTime: 2019-12-31 17:56:34
% DurationCPUTime: 1.00s
% Computational Cost: add. (11831->181), mult. (16473->217), div. (0->0), fcn. (6962->8), ass. (0->76)
t351 = -pkin(2) - pkin(3);
t350 = -mrSges(3,1) - mrSges(4,1);
t349 = Ifges(4,4) + Ifges(3,5);
t348 = Ifges(3,6) - Ifges(4,6);
t319 = -qJD(1) + qJD(4);
t327 = sin(qJ(5));
t347 = t319 * t327;
t330 = cos(qJ(5));
t346 = t319 * t330;
t329 = sin(qJ(1));
t332 = cos(qJ(1));
t311 = t329 * g(1) - t332 * g(2);
t307 = qJDD(1) * pkin(1) + t311;
t312 = -t332 * g(1) - t329 * g(2);
t333 = qJD(1) ^ 2;
t308 = -t333 * pkin(1) + t312;
t325 = sin(pkin(8));
t326 = cos(pkin(8));
t295 = t325 * t307 + t326 * t308;
t339 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t295;
t292 = -t333 * pkin(2) + t339;
t289 = t351 * t333 + t339;
t294 = t326 * t307 - t325 * t308;
t337 = -t333 * qJ(3) + qJDD(3) - t294;
t291 = t351 * qJDD(1) + t337;
t328 = sin(qJ(4));
t331 = cos(qJ(4));
t286 = t331 * t289 + t328 * t291;
t317 = t319 ^ 2;
t318 = -qJDD(1) + qJDD(4);
t284 = -t317 * pkin(4) + t318 * pkin(7) + t286;
t323 = g(3) - qJDD(2);
t281 = -t327 * t284 + t330 * t323;
t301 = (-mrSges(6,1) * t330 + mrSges(6,2) * t327) * t319;
t344 = qJD(5) * t319;
t302 = t327 * t318 + t330 * t344;
t310 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t346;
t279 = m(6) * t281 + qJDD(5) * mrSges(6,1) - t302 * mrSges(6,3) + qJD(5) * t310 - t301 * t347;
t282 = t330 * t284 + t327 * t323;
t303 = t330 * t318 - t327 * t344;
t309 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t347;
t280 = m(6) * t282 - qJDD(5) * mrSges(6,2) + t303 * mrSges(6,3) - qJD(5) * t309 + t301 * t346;
t340 = -t327 * t279 + t330 * t280;
t272 = m(5) * t286 - t317 * mrSges(5,1) - t318 * mrSges(5,2) + t340;
t285 = -t328 * t289 + t331 * t291;
t283 = -t318 * pkin(4) - t317 * pkin(7) - t285;
t334 = -m(6) * t283 + t303 * mrSges(6,1) - t302 * mrSges(6,2) - t309 * t347 + t310 * t346;
t277 = m(5) * t285 + t318 * mrSges(5,1) - t317 * mrSges(5,2) + t334;
t341 = t331 * t272 - t328 * t277;
t338 = m(4) * t292 + qJDD(1) * mrSges(4,3) + t341;
t267 = m(3) * t295 - qJDD(1) * mrSges(3,2) + t350 * t333 + t338;
t270 = t328 * t272 + t331 * t277;
t293 = -qJDD(1) * pkin(2) + t337;
t336 = -m(4) * t293 + qJDD(1) * mrSges(4,1) + t333 * mrSges(4,3) - t270;
t268 = m(3) * t294 + qJDD(1) * mrSges(3,1) - t333 * mrSges(3,2) + t336;
t262 = t325 * t267 + t326 * t268;
t260 = m(2) * t311 + qJDD(1) * mrSges(2,1) - t333 * mrSges(2,2) + t262;
t342 = t326 * t267 - t325 * t268;
t261 = m(2) * t312 - t333 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t342;
t345 = t332 * t260 + t329 * t261;
t343 = -t329 * t260 + t332 * t261;
t274 = t330 * t279 + t327 * t280;
t273 = -t274 + (-m(4) - m(5)) * t323;
t335 = -m(3) * t323 + t273;
t298 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t327 + Ifges(6,4) * t330) * t319;
t297 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t327 + Ifges(6,2) * t330) * t319;
t296 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t327 + Ifges(6,6) * t330) * t319;
t276 = mrSges(6,2) * t283 - mrSges(6,3) * t281 + Ifges(6,1) * t302 + Ifges(6,4) * t303 + Ifges(6,5) * qJDD(5) - qJD(5) * t297 + t296 * t346;
t275 = -mrSges(6,1) * t283 + mrSges(6,3) * t282 + Ifges(6,4) * t302 + Ifges(6,2) * t303 + Ifges(6,6) * qJDD(5) + qJD(5) * t298 - t296 * t347;
t269 = -mrSges(5,1) * t323 - mrSges(6,1) * t281 + mrSges(6,2) * t282 + mrSges(5,3) * t286 + t317 * Ifges(5,5) - Ifges(6,5) * t302 + Ifges(5,6) * t318 - Ifges(6,6) * t303 - Ifges(6,3) * qJDD(5) - pkin(4) * t274 + (-t297 * t327 + t298 * t330) * t319;
t263 = mrSges(5,2) * t323 - mrSges(5,3) * t285 + Ifges(5,5) * t318 - t317 * Ifges(5,6) - pkin(7) * t274 - t327 * t275 + t330 * t276;
t256 = mrSges(4,2) * t293 - mrSges(3,3) * t294 - pkin(6) * t270 - qJ(3) * t273 + t331 * t263 - t328 * t269 - t348 * t333 + (-mrSges(3,2) + mrSges(4,3)) * t323 + t349 * qJDD(1);
t255 = mrSges(3,3) * t295 + mrSges(4,2) * t292 - t328 * t263 - t331 * t269 + pkin(3) * t274 - pkin(6) * t341 - pkin(2) * t273 + t349 * t333 + (pkin(3) * m(5) - t350) * t323 + t348 * qJDD(1);
t254 = -mrSges(2,2) * g(3) - mrSges(2,3) * t311 + Ifges(2,5) * qJDD(1) - t333 * Ifges(2,6) - qJ(2) * t262 - t325 * t255 + t326 * t256;
t253 = mrSges(2,1) * g(3) + mrSges(2,3) * t312 + t333 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t335 + qJ(2) * t342 + t326 * t255 + t325 * t256;
t1 = [-m(1) * g(1) + t343; -m(1) * g(2) + t345; (-m(1) - m(2)) * g(3) + t335; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t345 - t329 * t253 + t332 * t254; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t343 + t332 * t253 + t329 * t254; pkin(1) * t262 + mrSges(2,1) * t311 - mrSges(2,2) * t312 + qJ(3) * (-t333 * mrSges(4,1) + t338) + pkin(2) * t336 - mrSges(3,2) * t295 + mrSges(3,1) * t294 - pkin(3) * t270 - mrSges(4,1) * t293 + mrSges(4,3) * t292 - t330 * t275 - pkin(4) * t334 - pkin(7) * t340 - mrSges(5,1) * t285 + mrSges(5,2) * t286 - t327 * t276 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - Ifges(5,3) * t318 + (Ifges(2,3) + Ifges(3,3) + Ifges(4,2)) * qJDD(1);];
tauB = t1;
