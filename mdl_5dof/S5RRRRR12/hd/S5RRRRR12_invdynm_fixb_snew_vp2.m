% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRR12
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRR12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR12_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR12_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR12_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR12_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR12_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:47:18
% EndTime: 2019-12-31 22:48:27
% DurationCPUTime: 48.62s
% Computational Cost: add. (839494->341), mult. (2086073->465), div. (0->0), fcn. (1723952->14), ass. (0->153)
t280 = cos(pkin(5));
t274 = qJD(1) * t280 + qJD(2);
t277 = sin(pkin(6));
t279 = cos(pkin(6));
t278 = sin(pkin(5));
t289 = cos(qJ(2));
t307 = qJD(1) * t289;
t303 = t278 * t307;
t258 = (t274 * t277 + t279 * t303) * pkin(9);
t284 = sin(qJ(2));
t309 = qJD(1) * t278;
t321 = pkin(9) * t277;
t262 = (-pkin(2) * t289 - t284 * t321) * t309;
t306 = qJD(1) * qJD(2);
t268 = (qJDD(1) * t284 + t289 * t306) * t278;
t273 = qJDD(1) * t280 + qJDD(2);
t285 = sin(qJ(1));
t290 = cos(qJ(1));
t271 = t285 * g(1) - g(2) * t290;
t291 = qJD(1) ^ 2;
t322 = pkin(8) * t278;
t265 = qJDD(1) * pkin(1) + t291 * t322 + t271;
t272 = -g(1) * t290 - g(2) * t285;
t266 = -pkin(1) * t291 + qJDD(1) * t322 + t272;
t311 = t280 * t289;
t301 = t265 * t311 - t284 * t266;
t308 = qJD(1) * t284;
t320 = pkin(9) * t279;
t219 = -t268 * t320 + t273 * pkin(2) + t274 * t258 + (-g(3) * t289 - t262 * t308) * t278 + t301;
t304 = t278 * t308;
t261 = pkin(2) * t274 - t304 * t320;
t269 = (qJDD(1) * t289 - t284 * t306) * t278;
t299 = t269 * t279 + t273 * t277;
t312 = t280 * t284;
t310 = t265 * t312 + t289 * t266;
t220 = -t274 * t261 + (-g(3) * t284 + t262 * t307) * t278 + t299 * pkin(9) + t310;
t319 = t280 * g(3);
t225 = -t268 * t321 - t269 * pkin(2) - t319 + (-t265 + (-t258 * t289 + t261 * t284) * qJD(1)) * t278;
t283 = sin(qJ(3));
t288 = cos(qJ(3));
t198 = -t283 * t220 + (t219 * t279 + t225 * t277) * t288;
t313 = t279 * t289;
t318 = t277 * t283;
t249 = t274 * t318 + (t283 * t313 + t284 * t288) * t309;
t235 = -t249 * qJD(3) - t283 * t268 + t288 * t299;
t317 = t277 * t288;
t248 = (-t283 * t284 + t288 * t313) * t309 + t274 * t317;
t316 = t278 * t284;
t315 = t278 * t289;
t314 = t279 * t283;
t199 = t219 * t314 + t288 * t220 + t225 * t318;
t238 = -pkin(3) * t248 - pkin(10) * t249;
t250 = -t269 * t277 + t273 * t279 + qJDD(3);
t259 = t274 * t279 - t277 * t303 + qJD(3);
t257 = t259 ^ 2;
t192 = -pkin(3) * t257 + pkin(10) * t250 + t238 * t248 + t199;
t204 = -t277 * t219 + t279 * t225;
t236 = t248 * qJD(3) + t288 * t268 + t283 * t299;
t194 = (-t248 * t259 - t236) * pkin(10) + (t249 * t259 - t235) * pkin(3) + t204;
t282 = sin(qJ(4));
t287 = cos(qJ(4));
t189 = t192 * t287 + t194 * t282;
t240 = -t249 * t282 + t259 * t287;
t241 = t249 * t287 + t259 * t282;
t222 = -pkin(4) * t240 - pkin(11) * t241;
t234 = qJDD(4) - t235;
t247 = qJD(4) - t248;
t246 = t247 ^ 2;
t186 = -pkin(4) * t246 + pkin(11) * t234 + t222 * t240 + t189;
t191 = -t250 * pkin(3) - t257 * pkin(10) + t249 * t238 - t198;
t208 = -qJD(4) * t241 - t236 * t282 + t250 * t287;
t209 = qJD(4) * t240 + t236 * t287 + t250 * t282;
t187 = (-t240 * t247 - t209) * pkin(11) + (t241 * t247 - t208) * pkin(4) + t191;
t281 = sin(qJ(5));
t286 = cos(qJ(5));
t183 = -t186 * t281 + t187 * t286;
t227 = -t241 * t281 + t247 * t286;
t197 = qJD(5) * t227 + t209 * t286 + t234 * t281;
t228 = t241 * t286 + t247 * t281;
t205 = -mrSges(6,1) * t227 + mrSges(6,2) * t228;
t207 = qJDD(5) - t208;
t239 = qJD(5) - t240;
t210 = -mrSges(6,2) * t239 + mrSges(6,3) * t227;
t180 = m(6) * t183 + mrSges(6,1) * t207 - mrSges(6,3) * t197 - t205 * t228 + t210 * t239;
t184 = t186 * t286 + t187 * t281;
t196 = -qJD(5) * t228 - t209 * t281 + t234 * t286;
t211 = mrSges(6,1) * t239 - mrSges(6,3) * t228;
t181 = m(6) * t184 - mrSges(6,2) * t207 + mrSges(6,3) * t196 + t205 * t227 - t211 * t239;
t174 = -t180 * t281 + t181 * t286;
t221 = -mrSges(5,1) * t240 + mrSges(5,2) * t241;
t230 = mrSges(5,1) * t247 - mrSges(5,3) * t241;
t172 = m(5) * t189 - mrSges(5,2) * t234 + mrSges(5,3) * t208 + t221 * t240 - t230 * t247 + t174;
t188 = -t192 * t282 + t194 * t287;
t185 = -pkin(4) * t234 - pkin(11) * t246 + t222 * t241 - t188;
t182 = -m(6) * t185 + mrSges(6,1) * t196 - mrSges(6,2) * t197 + t227 * t210 - t211 * t228;
t229 = -mrSges(5,2) * t247 + mrSges(5,3) * t240;
t178 = m(5) * t188 + mrSges(5,1) * t234 - mrSges(5,3) * t209 - t221 * t241 + t229 * t247 + t182;
t166 = t172 * t282 + t178 * t287;
t237 = -mrSges(4,1) * t248 + mrSges(4,2) * t249;
t243 = mrSges(4,1) * t259 - mrSges(4,3) * t249;
t302 = t172 * t287 - t178 * t282;
t163 = m(4) * t199 - mrSges(4,2) * t250 + mrSges(4,3) * t235 + t237 * t248 - t243 * t259 + t302;
t242 = -mrSges(4,2) * t259 + mrSges(4,3) * t248;
t165 = m(4) * t204 - mrSges(4,1) * t235 + mrSges(4,2) * t236 - t242 * t248 + t243 * t249 + t166;
t173 = t180 * t286 + t181 * t281;
t294 = -m(5) * t191 + t208 * mrSges(5,1) - mrSges(5,2) * t209 + t240 * t229 - t230 * t241 - t173;
t169 = m(4) * t198 + mrSges(4,1) * t250 - mrSges(4,3) * t236 - t237 * t249 + t242 * t259 + t294;
t152 = t163 * t318 + t165 * t279 + t169 * t317;
t153 = t169 * t279 * t288 + t163 * t314 - t165 * t277;
t244 = -g(3) * t315 + t301;
t264 = -mrSges(3,2) * t274 + mrSges(3,3) * t303;
t267 = (-mrSges(3,1) * t289 + mrSges(3,2) * t284) * t309;
t150 = m(3) * t244 + mrSges(3,1) * t273 - mrSges(3,3) * t268 + t264 * t274 - t267 * t304 + t153;
t157 = t163 * t288 - t169 * t283;
t245 = -g(3) * t316 + t310;
t263 = mrSges(3,1) * t274 - mrSges(3,3) * t304;
t156 = m(3) * t245 - mrSges(3,2) * t273 + mrSges(3,3) * t269 - t263 * t274 + t267 * t303 + t157;
t145 = -t150 * t284 + t156 * t289;
t254 = -t278 * t265 - t319;
t151 = m(3) * t254 - t269 * mrSges(3,1) + t268 * mrSges(3,2) + (t263 * t284 - t264 * t289) * t309 + t152;
t142 = t150 * t311 - t151 * t278 + t156 * t312;
t200 = Ifges(6,5) * t228 + Ifges(6,6) * t227 + Ifges(6,3) * t239;
t202 = Ifges(6,1) * t228 + Ifges(6,4) * t227 + Ifges(6,5) * t239;
t175 = -mrSges(6,1) * t185 + mrSges(6,3) * t184 + Ifges(6,4) * t197 + Ifges(6,2) * t196 + Ifges(6,6) * t207 - t200 * t228 + t202 * t239;
t201 = Ifges(6,4) * t228 + Ifges(6,2) * t227 + Ifges(6,6) * t239;
t176 = mrSges(6,2) * t185 - mrSges(6,3) * t183 + Ifges(6,1) * t197 + Ifges(6,4) * t196 + Ifges(6,5) * t207 + t200 * t227 - t201 * t239;
t212 = Ifges(5,5) * t241 + Ifges(5,6) * t240 + Ifges(5,3) * t247;
t213 = Ifges(5,4) * t241 + Ifges(5,2) * t240 + Ifges(5,6) * t247;
t158 = mrSges(5,2) * t191 - mrSges(5,3) * t188 + Ifges(5,1) * t209 + Ifges(5,4) * t208 + Ifges(5,5) * t234 - pkin(11) * t173 - t175 * t281 + t176 * t286 + t212 * t240 - t213 * t247;
t214 = Ifges(5,1) * t241 + Ifges(5,4) * t240 + Ifges(5,5) * t247;
t293 = mrSges(6,1) * t183 - mrSges(6,2) * t184 + Ifges(6,5) * t197 + Ifges(6,6) * t196 + Ifges(6,3) * t207 + t201 * t228 - t202 * t227;
t159 = -mrSges(5,1) * t191 + mrSges(5,3) * t189 + Ifges(5,4) * t209 + Ifges(5,2) * t208 + Ifges(5,6) * t234 - pkin(4) * t173 - t212 * t241 + t214 * t247 - t293;
t231 = Ifges(4,5) * t249 + Ifges(4,6) * t248 + Ifges(4,3) * t259;
t232 = Ifges(4,4) * t249 + Ifges(4,2) * t248 + Ifges(4,6) * t259;
t147 = mrSges(4,2) * t204 - mrSges(4,3) * t198 + Ifges(4,1) * t236 + Ifges(4,4) * t235 + Ifges(4,5) * t250 - pkin(10) * t166 + t158 * t287 - t159 * t282 + t231 * t248 - t232 * t259;
t233 = Ifges(4,1) * t249 + Ifges(4,4) * t248 + Ifges(4,5) * t259;
t292 = mrSges(5,1) * t188 - mrSges(5,2) * t189 + Ifges(5,5) * t209 + Ifges(5,6) * t208 + Ifges(5,3) * t234 + pkin(4) * t182 + pkin(11) * t174 + t175 * t286 + t176 * t281 + t213 * t241 - t214 * t240;
t148 = -mrSges(4,1) * t204 + mrSges(4,3) * t199 + Ifges(4,4) * t236 + Ifges(4,2) * t235 + Ifges(4,6) * t250 - pkin(3) * t166 - t231 * t249 + t233 * t259 - t292;
t296 = pkin(9) * t157 + t147 * t283 + t148 * t288;
t146 = mrSges(4,1) * t198 - mrSges(4,2) * t199 + Ifges(4,5) * t236 + Ifges(4,6) * t235 + Ifges(4,3) * t250 + pkin(3) * t294 + pkin(10) * t302 + t282 * t158 + t287 * t159 + t249 * t232 - t248 * t233;
t252 = Ifges(3,6) * t274 + (Ifges(3,4) * t284 + Ifges(3,2) * t289) * t309;
t253 = Ifges(3,5) * t274 + (Ifges(3,1) * t284 + Ifges(3,4) * t289) * t309;
t134 = mrSges(3,1) * t244 - mrSges(3,2) * t245 + Ifges(3,5) * t268 + Ifges(3,6) * t269 + Ifges(3,3) * t273 + pkin(2) * t153 + t279 * t146 + (t252 * t284 - t253 * t289) * t309 + t296 * t277;
t251 = Ifges(3,3) * t274 + (Ifges(3,5) * t284 + Ifges(3,6) * t289) * t309;
t136 = -mrSges(3,1) * t254 + mrSges(3,3) * t245 + Ifges(3,4) * t268 + Ifges(3,2) * t269 + Ifges(3,6) * t273 - pkin(2) * t152 - t277 * t146 - t251 * t304 + t274 * t253 + t279 * t296;
t138 = t251 * t303 + mrSges(3,2) * t254 - mrSges(3,3) * t244 + Ifges(3,1) * t268 + Ifges(3,4) * t269 + Ifges(3,5) * t273 + t288 * t147 - t283 * t148 - t274 * t252 + (-t152 * t277 - t153 * t279) * pkin(9);
t295 = mrSges(2,1) * t271 - mrSges(2,2) * t272 + Ifges(2,3) * qJDD(1) + pkin(1) * t142 + t134 * t280 + t136 * t315 + t138 * t316 + t145 * t322;
t143 = m(2) * t272 - mrSges(2,1) * t291 - qJDD(1) * mrSges(2,2) + t145;
t141 = t280 * t151 + (t150 * t289 + t156 * t284) * t278;
t139 = m(2) * t271 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t291 + t142;
t132 = -mrSges(2,2) * g(3) - mrSges(2,3) * t271 + Ifges(2,5) * qJDD(1) - t291 * Ifges(2,6) - t284 * t136 + t289 * t138 + (-t141 * t278 - t142 * t280) * pkin(8);
t131 = mrSges(2,1) * g(3) + mrSges(2,3) * t272 + t291 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t141 - t278 * t134 + (pkin(8) * t145 + t136 * t289 + t138 * t284) * t280;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t290 * t132 - t285 * t131 - pkin(7) * (t139 * t290 + t143 * t285), t132, t138, t147, t158, t176; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t285 * t132 + t290 * t131 + pkin(7) * (-t139 * t285 + t143 * t290), t131, t136, t148, t159, t175; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t295, t295, t134, t146, t292, t293;];
m_new = t1;
