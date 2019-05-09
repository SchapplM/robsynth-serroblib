% Calculate vector of cutting torques with Newton-Euler for
% S6PRRPPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% Datum: 2019-05-05 02:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRPPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:30:25
% EndTime: 2019-05-05 02:31:03
% DurationCPUTime: 28.92s
% Computational Cost: add. (509655->340), mult. (1134472->441), div. (0->0), fcn. (818654->14), ass. (0->143)
t353 = -2 * qJD(4);
t310 = sin(pkin(10));
t313 = cos(pkin(10));
t298 = g(1) * t310 - g(2) * t313;
t299 = -g(1) * t313 - g(2) * t310;
t307 = -g(3) + qJDD(1);
t320 = cos(qJ(2));
t314 = cos(pkin(6));
t317 = sin(qJ(2));
t345 = t314 * t317;
t311 = sin(pkin(6));
t346 = t311 * t317;
t258 = t298 * t345 + t299 * t320 + t307 * t346;
t322 = qJD(2) ^ 2;
t249 = -pkin(2) * t322 + qJDD(2) * pkin(8) + t258;
t278 = -t298 * t311 + t307 * t314;
t316 = sin(qJ(3));
t319 = cos(qJ(3));
t234 = -t316 * t249 + t278 * t319;
t341 = qJD(2) * qJD(3);
t340 = t319 * t341;
t295 = qJDD(2) * t316 + t340;
t229 = (-t295 + t340) * qJ(4) + (t316 * t319 * t322 + qJDD(3)) * pkin(3) + t234;
t235 = t249 * t319 + t278 * t316;
t296 = qJDD(2) * t319 - t316 * t341;
t344 = qJD(2) * t316;
t300 = qJD(3) * pkin(3) - qJ(4) * t344;
t306 = t319 ^ 2;
t230 = -pkin(3) * t306 * t322 + qJ(4) * t296 - qJD(3) * t300 + t235;
t309 = sin(pkin(11));
t343 = qJD(2) * t319;
t349 = cos(pkin(11));
t282 = t309 * t344 - t343 * t349;
t211 = t229 * t309 + t230 * t349 + t282 * t353;
t283 = (t309 * t319 + t316 * t349) * qJD(2);
t261 = mrSges(5,1) * t282 + mrSges(5,2) * t283;
t267 = t295 * t309 - t296 * t349;
t277 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t283;
t260 = pkin(4) * t282 - qJ(5) * t283;
t321 = qJD(3) ^ 2;
t208 = -pkin(4) * t321 + qJDD(3) * qJ(5) - t260 * t282 + t211;
t257 = -t317 * t299 + (t298 * t314 + t307 * t311) * t320;
t329 = -qJDD(2) * pkin(2) - t257;
t231 = -t296 * pkin(3) + qJDD(4) + t300 * t344 + (-qJ(4) * t306 - pkin(8)) * t322 + t329;
t268 = t295 * t349 + t296 * t309;
t215 = (qJD(3) * t282 - t268) * qJ(5) + (qJD(3) * t283 + t267) * pkin(4) + t231;
t308 = sin(pkin(12));
t312 = cos(pkin(12));
t273 = qJD(3) * t308 + t283 * t312;
t203 = -0.2e1 * qJD(5) * t273 - t308 * t208 + t215 * t312;
t253 = qJDD(3) * t308 + t268 * t312;
t272 = qJD(3) * t312 - t283 * t308;
t201 = (t272 * t282 - t253) * pkin(9) + (t272 * t273 + t267) * pkin(5) + t203;
t204 = 0.2e1 * qJD(5) * t272 + t208 * t312 + t215 * t308;
t250 = pkin(5) * t282 - pkin(9) * t273;
t252 = qJDD(3) * t312 - t268 * t308;
t271 = t272 ^ 2;
t202 = -pkin(5) * t271 + pkin(9) * t252 - t250 * t282 + t204;
t315 = sin(qJ(6));
t318 = cos(qJ(6));
t199 = t201 * t318 - t202 * t315;
t240 = t272 * t318 - t273 * t315;
t220 = qJD(6) * t240 + t252 * t315 + t253 * t318;
t241 = t272 * t315 + t273 * t318;
t225 = -mrSges(7,1) * t240 + mrSges(7,2) * t241;
t281 = qJD(6) + t282;
t232 = -mrSges(7,2) * t281 + mrSges(7,3) * t240;
t266 = qJDD(6) + t267;
t194 = m(7) * t199 + mrSges(7,1) * t266 - mrSges(7,3) * t220 - t225 * t241 + t232 * t281;
t200 = t201 * t315 + t202 * t318;
t219 = -qJD(6) * t241 + t252 * t318 - t253 * t315;
t233 = mrSges(7,1) * t281 - mrSges(7,3) * t241;
t195 = m(7) * t200 - mrSges(7,2) * t266 + mrSges(7,3) * t219 + t225 * t240 - t233 * t281;
t186 = t194 * t318 + t195 * t315;
t243 = -mrSges(6,1) * t272 + mrSges(6,2) * t273;
t246 = -mrSges(6,2) * t282 + mrSges(6,3) * t272;
t184 = m(6) * t203 + mrSges(6,1) * t267 - mrSges(6,3) * t253 - t243 * t273 + t246 * t282 + t186;
t247 = mrSges(6,1) * t282 - mrSges(6,3) * t273;
t336 = -t194 * t315 + t195 * t318;
t185 = m(6) * t204 - mrSges(6,2) * t267 + mrSges(6,3) * t252 + t243 * t272 - t247 * t282 + t336;
t337 = -t184 * t308 + t185 * t312;
t177 = m(5) * t211 - qJDD(3) * mrSges(5,2) - mrSges(5,3) * t267 - qJD(3) * t277 - t261 * t282 + t337;
t210 = t229 * t349 - t230 * t309 + t283 * t353;
t276 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t282;
t207 = -qJDD(3) * pkin(4) - qJ(5) * t321 + t260 * t283 + qJDD(5) - t210;
t205 = -pkin(5) * t252 - pkin(9) * t271 + t250 * t273 + t207;
t331 = m(7) * t205 - mrSges(7,1) * t219 + mrSges(7,2) * t220 - t232 * t240 + t233 * t241;
t326 = -m(6) * t207 + mrSges(6,1) * t252 - mrSges(6,2) * t253 + t246 * t272 - t247 * t273 - t331;
t190 = m(5) * t210 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t268 + qJD(3) * t276 - t261 * t283 + t326;
t168 = t177 * t309 + t190 * t349;
t294 = (-mrSges(4,1) * t319 + mrSges(4,2) * t316) * qJD(2);
t302 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t343;
t166 = m(4) * t234 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t295 + qJD(3) * t302 - t294 * t344 + t168;
t301 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t344;
t338 = t177 * t349 - t190 * t309;
t167 = m(4) * t235 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t296 - qJD(3) * t301 + t294 * t343 + t338;
t161 = t166 * t319 + t316 * t167;
t286 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t316 + Ifges(4,2) * t319) * qJD(2);
t287 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t316 + Ifges(4,4) * t319) * qJD(2);
t221 = Ifges(7,5) * t241 + Ifges(7,6) * t240 + Ifges(7,3) * t281;
t223 = Ifges(7,1) * t241 + Ifges(7,4) * t240 + Ifges(7,5) * t281;
t187 = -mrSges(7,1) * t205 + mrSges(7,3) * t200 + Ifges(7,4) * t220 + Ifges(7,2) * t219 + Ifges(7,6) * t266 - t221 * t241 + t223 * t281;
t222 = Ifges(7,4) * t241 + Ifges(7,2) * t240 + Ifges(7,6) * t281;
t188 = mrSges(7,2) * t205 - mrSges(7,3) * t199 + Ifges(7,1) * t220 + Ifges(7,4) * t219 + Ifges(7,5) * t266 + t221 * t240 - t222 * t281;
t236 = Ifges(6,5) * t273 + Ifges(6,6) * t272 + Ifges(6,3) * t282;
t238 = Ifges(6,1) * t273 + Ifges(6,4) * t272 + Ifges(6,5) * t282;
t170 = -mrSges(6,1) * t207 + mrSges(6,3) * t204 + Ifges(6,4) * t253 + Ifges(6,2) * t252 + Ifges(6,6) * t267 - pkin(5) * t331 + pkin(9) * t336 + t318 * t187 + t315 * t188 - t273 * t236 + t282 * t238;
t237 = Ifges(6,4) * t273 + Ifges(6,2) * t272 + Ifges(6,6) * t282;
t172 = mrSges(6,2) * t207 - mrSges(6,3) * t203 + Ifges(6,1) * t253 + Ifges(6,4) * t252 + Ifges(6,5) * t267 - pkin(9) * t186 - t187 * t315 + t188 * t318 + t236 * t272 - t237 * t282;
t255 = Ifges(5,4) * t283 - Ifges(5,2) * t282 + Ifges(5,6) * qJD(3);
t256 = Ifges(5,1) * t283 - Ifges(5,4) * t282 + Ifges(5,5) * qJD(3);
t327 = -mrSges(5,1) * t210 + mrSges(5,2) * t211 - Ifges(5,5) * t268 + Ifges(5,6) * t267 - Ifges(5,3) * qJDD(3) - pkin(4) * t326 - qJ(5) * t337 - t170 * t312 - t172 * t308 - t255 * t283 - t282 * t256;
t351 = mrSges(4,1) * t234 - mrSges(4,2) * t235 + Ifges(4,5) * t295 + Ifges(4,6) * t296 + Ifges(4,3) * qJDD(3) + pkin(3) * t168 + (t286 * t316 - t287 * t319) * qJD(2) - t327;
t147 = -mrSges(3,1) * t278 + mrSges(3,3) * t258 + t322 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t161 - t351;
t339 = -t166 * t316 + t167 * t319;
t159 = m(3) * t258 - mrSges(3,1) * t322 - qJDD(2) * mrSges(3,2) + t339;
t248 = -t322 * pkin(8) + t329;
t179 = t184 * t312 + t185 * t308;
t328 = m(5) * t231 + mrSges(5,1) * t267 + mrSges(5,2) * t268 + t276 * t282 + t277 * t283 + t179;
t325 = -m(4) * t248 + mrSges(4,1) * t296 - mrSges(4,2) * t295 - t301 * t344 + t302 * t343 - t328;
t174 = m(3) * t257 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t322 + t325;
t155 = t159 * t320 - t174 * t317;
t352 = pkin(7) * t155 + t147 * t320;
t347 = t174 * t320;
t160 = m(3) * t278 + t161;
t152 = t159 * t345 - t160 * t311 + t314 * t347;
t254 = Ifges(5,5) * t283 - Ifges(5,6) * t282 + Ifges(5,3) * qJD(3);
t156 = mrSges(5,2) * t231 - mrSges(5,3) * t210 + Ifges(5,1) * t268 - Ifges(5,4) * t267 + Ifges(5,5) * qJDD(3) - qJ(5) * t179 - qJD(3) * t255 - t170 * t308 + t172 * t312 - t254 * t282;
t330 = -mrSges(7,1) * t199 + mrSges(7,2) * t200 - Ifges(7,5) * t220 - Ifges(7,6) * t219 - Ifges(7,3) * t266 - t222 * t241 + t240 * t223;
t324 = -mrSges(6,1) * t203 + mrSges(6,2) * t204 - Ifges(6,5) * t253 - Ifges(6,6) * t252 - pkin(5) * t186 - t273 * t237 + t272 * t238 + t330;
t162 = Ifges(5,6) * qJDD(3) + t324 + (-Ifges(5,2) - Ifges(6,3)) * t267 - t283 * t254 + Ifges(5,4) * t268 + qJD(3) * t256 - mrSges(5,1) * t231 + mrSges(5,3) * t211 - pkin(4) * t179;
t285 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t316 + Ifges(4,6) * t319) * qJD(2);
t145 = -mrSges(4,1) * t248 + mrSges(4,3) * t235 + Ifges(4,4) * t295 + Ifges(4,2) * t296 + Ifges(4,6) * qJDD(3) - pkin(3) * t328 + qJ(4) * t338 + qJD(3) * t287 + t309 * t156 + t162 * t349 - t285 * t344;
t148 = mrSges(4,2) * t248 - mrSges(4,3) * t234 + Ifges(4,1) * t295 + Ifges(4,4) * t296 + Ifges(4,5) * qJDD(3) - qJ(4) * t168 - qJD(3) * t286 + t156 * t349 - t162 * t309 + t285 * t343;
t142 = mrSges(3,1) * t257 - mrSges(3,2) * t258 + Ifges(3,3) * qJDD(2) + pkin(2) * t325 + pkin(8) * t339 + t319 * t145 + t316 * t148;
t144 = mrSges(3,2) * t278 - mrSges(3,3) * t257 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t322 - pkin(8) * t161 - t145 * t316 + t148 * t319;
t332 = mrSges(2,1) * t298 - mrSges(2,2) * t299 + pkin(1) * t152 + t314 * t142 + t144 * t346 + t311 * t352;
t153 = m(2) * t299 + t155;
t151 = t314 * t160 + (t159 * t317 + t347) * t311;
t149 = m(2) * t298 + t152;
t140 = mrSges(2,2) * t307 - mrSges(2,3) * t298 + t320 * t144 - t317 * t147 + (-t151 * t311 - t152 * t314) * pkin(7);
t139 = -mrSges(2,1) * t307 + mrSges(2,3) * t299 - pkin(1) * t151 - t311 * t142 + (t144 * t317 + t352) * t314;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t313 * t140 - t310 * t139 - qJ(1) * (t149 * t313 + t153 * t310), t140, t144, t148, t156, t172, t188; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t310 * t140 + t313 * t139 + qJ(1) * (-t149 * t310 + t153 * t313), t139, t147, t145, t162, t170, t187; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t332, t332, t142, t351, -t327, Ifges(6,3) * t267 - t324, -t330;];
m_new  = t1;
