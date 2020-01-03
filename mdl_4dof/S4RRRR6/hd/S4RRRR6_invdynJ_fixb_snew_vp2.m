% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S4RRRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR6_invdynJ_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_invdynJ_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR6_invdynJ_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR6_invdynJ_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:38
% EndTime: 2019-12-31 17:29:40
% DurationCPUTime: 1.31s
% Computational Cost: add. (10313->214), mult. (22311->291), div. (0->0), fcn. (16547->10), ass. (0->97)
t308 = sin(pkin(4));
t312 = sin(qJ(2));
t316 = cos(qJ(2));
t328 = qJD(1) * qJD(2);
t300 = (-qJDD(1) * t316 + t312 * t328) * t308;
t336 = t308 * pkin(6);
t309 = cos(pkin(4));
t335 = t309 * g(3);
t318 = qJD(1) ^ 2;
t313 = sin(qJ(1));
t317 = cos(qJ(1));
t325 = t313 * g(1) - t317 * g(2);
t295 = qJDD(1) * pkin(1) + t318 * t336 + t325;
t334 = t295 * t309;
t333 = t308 * t312;
t332 = t308 * t316;
t330 = qJD(1) * t308;
t298 = (-t316 * pkin(2) - t312 * pkin(7)) * t330;
t305 = t309 * qJD(1) + qJD(2);
t303 = t305 ^ 2;
t304 = t309 * qJDD(1) + qJDD(2);
t329 = qJD(1) * t316;
t323 = -t317 * g(1) - t313 * g(2);
t296 = -t318 * pkin(1) + qJDD(1) * t336 + t323;
t331 = t316 * t296 + t312 * t334;
t260 = -t303 * pkin(2) + t304 * pkin(7) + (-g(3) * t312 + t298 * t329) * t308 + t331;
t299 = (qJDD(1) * t312 + t316 * t328) * t308;
t261 = t300 * pkin(2) - t299 * pkin(7) - t335 + (-t295 + (pkin(2) * t312 - pkin(7) * t316) * t305 * qJD(1)) * t308;
t311 = sin(qJ(3));
t315 = cos(qJ(3));
t249 = t315 * t260 + t311 * t261;
t327 = t312 * t330;
t288 = t315 * t305 - t311 * t327;
t289 = t311 * t305 + t315 * t327;
t276 = -t288 * pkin(3) - t289 * pkin(8);
t292 = qJDD(3) + t300;
t326 = t308 * t329;
t302 = qJD(3) - t326;
t301 = t302 ^ 2;
t246 = -t301 * pkin(3) + t292 * pkin(8) + t288 * t276 + t249;
t273 = -g(3) * t332 - t312 * t296 + t316 * t334;
t259 = -t304 * pkin(2) - t303 * pkin(7) + t298 * t327 - t273;
t271 = -t289 * qJD(3) - t311 * t299 + t315 * t304;
t272 = t288 * qJD(3) + t315 * t299 + t311 * t304;
t247 = (-t288 * t302 - t272) * pkin(8) + (t289 * t302 - t271) * pkin(3) + t259;
t310 = sin(qJ(4));
t314 = cos(qJ(4));
t243 = -t310 * t246 + t314 * t247;
t277 = -t310 * t289 + t314 * t302;
t252 = t277 * qJD(4) + t314 * t272 + t310 * t292;
t278 = t314 * t289 + t310 * t302;
t262 = -t277 * mrSges(5,1) + t278 * mrSges(5,2);
t287 = qJD(4) - t288;
t263 = -t287 * mrSges(5,2) + t277 * mrSges(5,3);
t269 = qJDD(4) - t271;
t240 = m(5) * t243 + t269 * mrSges(5,1) - t252 * mrSges(5,3) - t278 * t262 + t287 * t263;
t244 = t314 * t246 + t310 * t247;
t251 = -t278 * qJD(4) - t310 * t272 + t314 * t292;
t264 = t287 * mrSges(5,1) - t278 * mrSges(5,3);
t241 = m(5) * t244 - t269 * mrSges(5,2) + t251 * mrSges(5,3) + t277 * t262 - t287 * t264;
t234 = -t310 * t240 + t314 * t241;
t275 = -t288 * mrSges(4,1) + t289 * mrSges(4,2);
t280 = t302 * mrSges(4,1) - t289 * mrSges(4,3);
t232 = m(4) * t249 - t292 * mrSges(4,2) + t271 * mrSges(4,3) + t288 * t275 - t302 * t280 + t234;
t248 = -t311 * t260 + t315 * t261;
t245 = -t292 * pkin(3) - t301 * pkin(8) + t289 * t276 - t248;
t242 = -m(5) * t245 + t251 * mrSges(5,1) - t252 * mrSges(5,2) + t277 * t263 - t278 * t264;
t279 = -t302 * mrSges(4,2) + t288 * mrSges(4,3);
t238 = m(4) * t248 + t292 * mrSges(4,1) - t272 * mrSges(4,3) - t289 * t275 + t302 * t279 + t242;
t228 = t311 * t232 + t315 * t238;
t324 = t315 * t232 - t311 * t238;
t233 = t314 * t240 + t310 * t241;
t321 = -m(4) * t259 + t271 * mrSges(4,1) - t272 * mrSges(4,2) + t288 * t279 - t289 * t280 - t233;
t254 = Ifges(5,4) * t278 + Ifges(5,2) * t277 + Ifges(5,6) * t287;
t255 = Ifges(5,1) * t278 + Ifges(5,4) * t277 + Ifges(5,5) * t287;
t320 = mrSges(5,1) * t243 - mrSges(5,2) * t244 + Ifges(5,5) * t252 + Ifges(5,6) * t251 + Ifges(5,3) * t269 + t278 * t254 - t277 * t255;
t253 = Ifges(5,5) * t278 + Ifges(5,6) * t277 + Ifges(5,3) * t287;
t235 = -mrSges(5,1) * t245 + mrSges(5,3) * t244 + Ifges(5,4) * t252 + Ifges(5,2) * t251 + Ifges(5,6) * t269 - t278 * t253 + t287 * t255;
t236 = mrSges(5,2) * t245 - mrSges(5,3) * t243 + Ifges(5,1) * t252 + Ifges(5,4) * t251 + Ifges(5,5) * t269 + t277 * t253 - t287 * t254;
t266 = Ifges(4,4) * t289 + Ifges(4,2) * t288 + Ifges(4,6) * t302;
t267 = Ifges(4,1) * t289 + Ifges(4,4) * t288 + Ifges(4,5) * t302;
t319 = mrSges(4,1) * t248 - mrSges(4,2) * t249 + Ifges(4,5) * t272 + Ifges(4,6) * t271 + Ifges(4,3) * t292 + pkin(3) * t242 + pkin(8) * t234 + t314 * t235 + t310 * t236 + t289 * t266 - t288 * t267;
t297 = (-t316 * mrSges(3,1) + t312 * mrSges(3,2)) * t330;
t294 = -t305 * mrSges(3,2) + mrSges(3,3) * t326;
t293 = t305 * mrSges(3,1) - mrSges(3,3) * t327;
t284 = -t308 * t295 - t335;
t283 = Ifges(3,5) * t305 + (t312 * Ifges(3,1) + t316 * Ifges(3,4)) * t330;
t282 = Ifges(3,6) * t305 + (t312 * Ifges(3,4) + t316 * Ifges(3,2)) * t330;
t281 = Ifges(3,3) * t305 + (t312 * Ifges(3,5) + t316 * Ifges(3,6)) * t330;
t274 = -g(3) * t333 + t331;
t265 = Ifges(4,5) * t289 + Ifges(4,6) * t288 + Ifges(4,3) * t302;
t229 = m(3) * t273 + t304 * mrSges(3,1) - t299 * mrSges(3,3) + t305 * t294 - t297 * t327 + t321;
t227 = m(3) * t274 - t304 * mrSges(3,2) - t300 * mrSges(3,3) - t305 * t293 + t297 * t326 + t324;
t226 = -mrSges(4,1) * t259 + mrSges(4,3) * t249 + Ifges(4,4) * t272 + Ifges(4,2) * t271 + Ifges(4,6) * t292 - pkin(3) * t233 - t289 * t265 + t302 * t267 - t320;
t225 = mrSges(4,2) * t259 - mrSges(4,3) * t248 + Ifges(4,1) * t272 + Ifges(4,4) * t271 + Ifges(4,5) * t292 - pkin(8) * t233 - t310 * t235 + t314 * t236 + t288 * t265 - t302 * t266;
t224 = Ifges(3,5) * t299 - Ifges(3,6) * t300 + Ifges(3,3) * t304 + mrSges(3,1) * t273 - mrSges(3,2) * t274 + t311 * t225 + t315 * t226 + pkin(2) * t321 + pkin(7) * t324 + (t312 * t282 - t316 * t283) * t330;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t325 - mrSges(2,2) * t323 + (mrSges(3,2) * t284 - mrSges(3,3) * t273 + Ifges(3,1) * t299 - Ifges(3,4) * t300 + Ifges(3,5) * t304 - pkin(7) * t228 + t315 * t225 - t311 * t226 + t281 * t326 - t305 * t282) * t333 + (-mrSges(3,1) * t284 + mrSges(3,3) * t274 + Ifges(3,4) * t299 - Ifges(3,2) * t300 + Ifges(3,6) * t304 - pkin(2) * t228 - t281 * t327 + t305 * t283 - t319) * t332 + t309 * t224 + pkin(1) * ((t312 * t227 + t316 * t229) * t309 + (-m(3) * t284 - t300 * mrSges(3,1) - t299 * mrSges(3,2) + (-t293 * t312 + t294 * t316) * t330 - t228) * t308) + (t316 * t227 - t312 * t229) * t336; t224; t319; t320;];
tauJ = t1;
