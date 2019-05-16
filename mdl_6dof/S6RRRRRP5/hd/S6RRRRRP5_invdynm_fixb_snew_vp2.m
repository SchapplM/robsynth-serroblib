% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-08 05:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:56:58
% EndTime: 2019-05-08 04:58:10
% DurationCPUTime: 26.18s
% Computational Cost: add. (467995->378), mult. (960403->462), div. (0->0), fcn. (704413->10), ass. (0->143)
t337 = sin(qJ(3));
t342 = cos(qJ(3));
t338 = sin(qJ(2));
t369 = qJD(1) * t338;
t318 = qJD(2) * t337 + t342 * t369;
t343 = cos(qJ(2));
t367 = qJD(1) * qJD(2);
t364 = t343 * t367;
t321 = qJDD(1) * t338 + t364;
t286 = -qJD(3) * t318 + qJDD(2) * t342 - t321 * t337;
t317 = qJD(2) * t342 - t337 * t369;
t287 = qJD(3) * t317 + qJDD(2) * t337 + t321 * t342;
t336 = sin(qJ(4));
t341 = cos(qJ(4));
t290 = t317 * t336 + t318 * t341;
t252 = -qJD(4) * t290 + t286 * t341 - t287 * t336;
t289 = t317 * t341 - t318 * t336;
t253 = qJD(4) * t289 + t286 * t336 + t287 * t341;
t335 = sin(qJ(5));
t340 = cos(qJ(5));
t268 = t289 * t340 - t290 * t335;
t224 = qJD(5) * t268 + t252 * t335 + t253 * t340;
t269 = t289 * t335 + t290 * t340;
t245 = -mrSges(7,1) * t268 + mrSges(7,2) * t269;
t339 = sin(qJ(1));
t344 = cos(qJ(1));
t327 = t339 * g(1) - t344 * g(2);
t346 = qJD(1) ^ 2;
t310 = -qJDD(1) * pkin(1) - t346 * pkin(7) - t327;
t331 = t338 * t367;
t322 = qJDD(1) * t343 - t331;
t273 = (-t321 - t364) * pkin(8) + (-t322 + t331) * pkin(2) + t310;
t328 = -g(1) * t344 - g(2) * t339;
t311 = -pkin(1) * t346 + qJDD(1) * pkin(7) + t328;
t295 = -g(3) * t338 + t343 * t311;
t320 = (-pkin(2) * t343 - pkin(8) * t338) * qJD(1);
t345 = qJD(2) ^ 2;
t368 = qJD(1) * t343;
t276 = -pkin(2) * t345 + qJDD(2) * pkin(8) + t320 * t368 + t295;
t254 = t342 * t273 - t337 * t276;
t316 = qJDD(3) - t322;
t330 = qJD(3) - t368;
t231 = (t317 * t330 - t287) * pkin(9) + (t317 * t318 + t316) * pkin(3) + t254;
t255 = t337 * t273 + t342 * t276;
t296 = pkin(3) * t330 - pkin(9) * t318;
t315 = t317 ^ 2;
t233 = -pkin(3) * t315 + pkin(9) * t286 - t296 * t330 + t255;
t206 = t341 * t231 - t336 * t233;
t312 = qJDD(4) + t316;
t329 = qJD(4) + t330;
t202 = (t289 * t329 - t253) * pkin(10) + (t289 * t290 + t312) * pkin(4) + t206;
t207 = t336 * t231 + t341 * t233;
t279 = pkin(4) * t329 - pkin(10) * t290;
t288 = t289 ^ 2;
t204 = -pkin(4) * t288 + pkin(10) * t252 - t279 * t329 + t207;
t195 = t340 * t202 - t335 * t204;
t306 = qJDD(5) + t312;
t324 = qJD(5) + t329;
t190 = -0.2e1 * qJD(6) * t269 + (t268 * t324 - t224) * qJ(6) + (t268 * t269 + t306) * pkin(5) + t195;
t257 = -mrSges(7,2) * t324 + mrSges(7,3) * t268;
t366 = m(7) * t190 + t306 * mrSges(7,1) + t324 * t257;
t187 = -t224 * mrSges(7,3) - t269 * t245 + t366;
t196 = t335 * t202 + t340 * t204;
t223 = -qJD(5) * t269 + t252 * t340 - t253 * t335;
t238 = Ifges(6,4) * t269 + Ifges(6,2) * t268 + Ifges(6,6) * t324;
t239 = Ifges(7,1) * t269 + Ifges(7,4) * t268 + Ifges(7,5) * t324;
t240 = Ifges(6,1) * t269 + Ifges(6,4) * t268 + Ifges(6,5) * t324;
t259 = pkin(5) * t324 - qJ(6) * t269;
t267 = t268 ^ 2;
t193 = -pkin(5) * t267 + qJ(6) * t223 + 0.2e1 * qJD(6) * t268 - t259 * t324 + t196;
t237 = Ifges(7,4) * t269 + Ifges(7,2) * t268 + Ifges(7,6) * t324;
t356 = -mrSges(7,1) * t190 + mrSges(7,2) * t193 - Ifges(7,5) * t224 - Ifges(7,6) * t223 - Ifges(7,3) * t306 - t269 * t237;
t374 = mrSges(6,1) * t195 - mrSges(6,2) * t196 + Ifges(6,5) * t224 + Ifges(6,6) * t223 + Ifges(6,3) * t306 + pkin(5) * t187 + t269 * t238 - t356 - (t240 + t239) * t268;
t246 = -mrSges(6,1) * t268 + mrSges(6,2) * t269;
t258 = -mrSges(6,2) * t324 + mrSges(6,3) * t268;
t180 = m(6) * t195 + t306 * mrSges(6,1) + t324 * t258 + (-t245 - t246) * t269 + (-mrSges(6,3) - mrSges(7,3)) * t224 + t366;
t260 = mrSges(7,1) * t324 - mrSges(7,3) * t269;
t261 = mrSges(6,1) * t324 - mrSges(6,3) * t269;
t365 = m(7) * t193 + t223 * mrSges(7,3) + t268 * t245;
t183 = m(6) * t196 + t223 * mrSges(6,3) + t268 * t246 + (-t260 - t261) * t324 + (-mrSges(6,2) - mrSges(7,2)) * t306 + t365;
t178 = t340 * t180 + t335 * t183;
t263 = Ifges(5,4) * t290 + Ifges(5,2) * t289 + Ifges(5,6) * t329;
t264 = Ifges(5,1) * t290 + Ifges(5,4) * t289 + Ifges(5,5) * t329;
t373 = mrSges(5,1) * t206 - mrSges(5,2) * t207 + Ifges(5,5) * t253 + Ifges(5,6) * t252 + Ifges(5,3) * t312 + pkin(4) * t178 + t290 * t263 - t289 * t264 + t374;
t270 = -mrSges(5,1) * t289 + mrSges(5,2) * t290;
t277 = -mrSges(5,2) * t329 + mrSges(5,3) * t289;
t174 = m(5) * t206 + mrSges(5,1) * t312 - mrSges(5,3) * t253 - t270 * t290 + t277 * t329 + t178;
t278 = mrSges(5,1) * t329 - mrSges(5,3) * t290;
t360 = -t180 * t335 + t340 * t183;
t175 = m(5) * t207 - mrSges(5,2) * t312 + mrSges(5,3) * t252 + t270 * t289 - t278 * t329 + t360;
t169 = t341 * t174 + t336 * t175;
t281 = Ifges(4,4) * t318 + Ifges(4,2) * t317 + Ifges(4,6) * t330;
t282 = Ifges(4,1) * t318 + Ifges(4,4) * t317 + Ifges(4,5) * t330;
t372 = mrSges(4,1) * t254 - mrSges(4,2) * t255 + Ifges(4,5) * t287 + Ifges(4,6) * t286 + Ifges(4,3) * t316 + pkin(3) * t169 + t318 * t281 - t317 * t282 + t373;
t294 = -t343 * g(3) - t338 * t311;
t275 = -qJDD(2) * pkin(2) - pkin(8) * t345 + t320 * t369 - t294;
t247 = -pkin(3) * t286 - pkin(9) * t315 + t318 * t296 + t275;
t209 = -pkin(4) * t252 - pkin(10) * t288 + t290 * t279 + t247;
t235 = Ifges(7,5) * t269 + Ifges(7,6) * t268 + Ifges(7,3) * t324;
t236 = Ifges(6,5) * t269 + Ifges(6,6) * t268 + Ifges(6,3) * t324;
t199 = -pkin(5) * t223 - qJ(6) * t267 + t259 * t269 + qJDD(6) + t209;
t357 = -mrSges(7,1) * t199 + mrSges(7,3) * t193 + Ifges(7,4) * t224 + Ifges(7,2) * t223 + Ifges(7,6) * t306 + t324 * t239;
t359 = m(7) * t199 - t223 * mrSges(7,1) + t224 * mrSges(7,2) - t268 * t257 + t269 * t260;
t170 = Ifges(6,4) * t224 + Ifges(6,2) * t223 + Ifges(6,6) * t306 + t324 * t240 - mrSges(6,1) * t209 + mrSges(6,3) * t196 - pkin(5) * t359 + qJ(6) * (-t306 * mrSges(7,2) - t324 * t260 + t365) + (-t236 - t235) * t269 + t357;
t355 = mrSges(7,2) * t199 - mrSges(7,3) * t190 + Ifges(7,1) * t224 + Ifges(7,4) * t223 + Ifges(7,5) * t306 + t268 * t235;
t176 = mrSges(6,2) * t209 - mrSges(6,3) * t195 + Ifges(6,1) * t224 + Ifges(6,4) * t223 + Ifges(6,5) * t306 - qJ(6) * t187 + t268 * t236 + (-t237 - t238) * t324 + t355;
t262 = Ifges(5,5) * t290 + Ifges(5,6) * t289 + Ifges(5,3) * t329;
t353 = m(6) * t209 - t223 * mrSges(6,1) + t224 * mrSges(6,2) - t268 * t258 + t269 * t261 + t359;
t164 = -mrSges(5,1) * t247 + mrSges(5,3) * t207 + Ifges(5,4) * t253 + Ifges(5,2) * t252 + Ifges(5,6) * t312 - pkin(4) * t353 + pkin(10) * t360 + t340 * t170 + t335 * t176 - t290 * t262 + t329 * t264;
t165 = mrSges(5,2) * t247 - mrSges(5,3) * t206 + Ifges(5,1) * t253 + Ifges(5,4) * t252 + Ifges(5,5) * t312 - pkin(10) * t178 - t170 * t335 + t176 * t340 + t262 * t289 - t263 * t329;
t280 = Ifges(4,5) * t318 + Ifges(4,6) * t317 + Ifges(4,3) * t330;
t350 = m(5) * t247 - t252 * mrSges(5,1) + t253 * mrSges(5,2) - t289 * t277 + t290 * t278 + t353;
t361 = -t174 * t336 + t341 * t175;
t151 = -mrSges(4,1) * t275 + mrSges(4,3) * t255 + Ifges(4,4) * t287 + Ifges(4,2) * t286 + Ifges(4,6) * t316 - pkin(3) * t350 + pkin(9) * t361 + t341 * t164 + t336 * t165 - t318 * t280 + t330 * t282;
t152 = mrSges(4,2) * t275 - mrSges(4,3) * t254 + Ifges(4,1) * t287 + Ifges(4,4) * t286 + Ifges(4,5) * t316 - pkin(9) * t169 - t164 * t336 + t165 * t341 + t280 * t317 - t281 * t330;
t291 = -mrSges(4,1) * t317 + mrSges(4,2) * t318;
t292 = -mrSges(4,2) * t330 + mrSges(4,3) * t317;
t167 = m(4) * t254 + mrSges(4,1) * t316 - mrSges(4,3) * t287 - t291 * t318 + t292 * t330 + t169;
t293 = mrSges(4,1) * t330 - mrSges(4,3) * t318;
t168 = m(4) * t255 - mrSges(4,2) * t316 + mrSges(4,3) * t286 + t291 * t317 - t293 * t330 + t361;
t163 = -t167 * t337 + t342 * t168;
t185 = -m(4) * t275 + t286 * mrSges(4,1) - t287 * mrSges(4,2) + t317 * t292 - t318 * t293 - t350;
t308 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t338 + Ifges(3,2) * t343) * qJD(1);
t309 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t338 + Ifges(3,4) * t343) * qJD(1);
t371 = mrSges(3,1) * t294 - mrSges(3,2) * t295 + Ifges(3,5) * t321 + Ifges(3,6) * t322 + Ifges(3,3) * qJDD(2) + pkin(2) * t185 + pkin(8) * t163 + t342 * t151 + t337 * t152 + (t308 * t338 - t309 * t343) * qJD(1);
t319 = (-mrSges(3,1) * t343 + mrSges(3,2) * t338) * qJD(1);
t325 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t369;
t161 = m(3) * t295 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t322 - qJD(2) * t325 + t319 * t368 + t163;
t326 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t368;
t184 = m(3) * t294 + qJDD(2) * mrSges(3,1) - t321 * mrSges(3,3) + qJD(2) * t326 - t319 * t369 + t185;
t362 = t343 * t161 - t184 * t338;
t162 = t167 * t342 + t168 * t337;
t307 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t338 + Ifges(3,6) * t343) * qJD(1);
t150 = mrSges(3,2) * t310 - mrSges(3,3) * t294 + Ifges(3,1) * t321 + Ifges(3,4) * t322 + Ifges(3,5) * qJDD(2) - pkin(8) * t162 - qJD(2) * t308 - t151 * t337 + t152 * t342 + t307 * t368;
t154 = -mrSges(3,1) * t310 + mrSges(3,3) * t295 + Ifges(3,4) * t321 + Ifges(3,2) * t322 + Ifges(3,6) * qJDD(2) - pkin(2) * t162 + qJD(2) * t309 - t307 * t369 - t372;
t352 = -m(3) * t310 + t322 * mrSges(3,1) - mrSges(3,2) * t321 - t325 * t369 + t326 * t368 - t162;
t354 = mrSges(2,1) * t327 - mrSges(2,2) * t328 + Ifges(2,3) * qJDD(1) + pkin(1) * t352 + pkin(7) * t362 + t338 * t150 + t343 * t154;
t158 = m(2) * t327 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t346 + t352;
t157 = t161 * t338 + t184 * t343;
t155 = m(2) * t328 - mrSges(2,1) * t346 - qJDD(1) * mrSges(2,2) + t362;
t148 = mrSges(2,1) * g(3) + mrSges(2,3) * t328 + t346 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t157 - t371;
t147 = -mrSges(2,2) * g(3) - mrSges(2,3) * t327 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t346 - pkin(7) * t157 + t150 * t343 - t154 * t338;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t344 * t147 - t339 * t148 - pkin(6) * (t155 * t339 + t158 * t344), t147, t150, t152, t165, t176, -t237 * t324 + t355; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t339 * t147 + t344 * t148 + pkin(6) * (t155 * t344 - t158 * t339), t148, t154, t151, t164, t170, -t269 * t235 + t357; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t354, t354, t371, t372, t373, t374, -t268 * t239 - t356;];
m_new  = t1;
