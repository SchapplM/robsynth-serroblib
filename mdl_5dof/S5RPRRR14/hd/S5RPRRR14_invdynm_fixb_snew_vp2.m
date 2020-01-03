% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRR14
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
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRR14_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR14_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR14_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR14_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:17:26
% EndTime: 2019-12-31 19:18:11
% DurationCPUTime: 41.08s
% Computational Cost: add. (622978->327), mult. (1941365->460), div. (0->0), fcn. (1620898->14), ass. (0->156)
t275 = sin(pkin(11));
t277 = sin(pkin(5));
t278 = cos(pkin(11));
t280 = cos(pkin(5));
t283 = sin(qJ(3));
t279 = cos(pkin(6));
t287 = cos(qJ(3));
t317 = t279 * t287;
t276 = sin(pkin(6));
t322 = t276 * t287;
t293 = t277 * (-t275 * t283 + t278 * t317) + t280 * t322;
t248 = t293 * qJD(1);
t318 = t279 * t283;
t323 = t276 * t283;
t295 = t280 * t323 + (t275 * t287 + t278 * t318) * t277;
t249 = t295 * qJD(1);
t237 = -t249 * qJD(3) + qJDD(1) * t293;
t320 = t277 * t279;
t260 = (t276 * t280 + t278 * t320) * qJD(1) * pkin(8);
t284 = sin(qJ(1));
t288 = cos(qJ(1));
t272 = -g(1) * t288 - g(2) * t284;
t289 = qJD(1) ^ 2;
t326 = qJ(2) * t277;
t264 = -pkin(1) * t289 + qJDD(1) * t326 + t272;
t329 = pkin(8) * t275;
t305 = -pkin(2) * t278 - t276 * t329;
t316 = qJD(1) * t277;
t327 = pkin(8) * qJDD(1);
t300 = qJD(1) * t305 * t316 + t279 * t327;
t271 = t284 * g(1) - g(2) * t288;
t263 = qJDD(1) * pkin(1) + t289 * t326 + t271;
t312 = qJD(2) * t316;
t319 = t278 * t280;
t321 = t277 * t278;
t306 = -g(3) * t321 + t263 * t319 - 0.2e1 * t275 * t312;
t219 = (pkin(2) * qJDD(1) + qJD(1) * t260) * t280 + (-t277 * t300 - t264) * t275 + t306;
t265 = (pkin(2) * t280 - t320 * t329) * qJD(1);
t324 = t275 * t280;
t313 = t263 * t324 + (t264 + 0.2e1 * t312) * t278;
t220 = (-qJD(1) * t265 + t276 * t327) * t280 + (-g(3) * t275 + t278 * t300) * t277 + t313;
t311 = -t280 * g(3) + qJDD(2);
t228 = (-t263 + t305 * qJDD(1) + (-t260 * t278 + t265 * t275) * qJD(1)) * t277 + t311;
t195 = -t283 * t220 + (t219 * t279 + t228 * t276) * t287;
t328 = Ifges(3,3) * t280;
t325 = t275 * t277;
t196 = t219 * t318 + t220 * t287 + t228 * t323;
t236 = -pkin(3) * t248 - pkin(9) * t249;
t301 = -t276 * t321 + t279 * t280;
t261 = qJD(1) * t301 + qJD(3);
t257 = t261 ^ 2;
t258 = qJDD(1) * t301 + qJDD(3);
t192 = -pkin(3) * t257 + pkin(9) * t258 + t236 * t248 + t196;
t204 = -t276 * t219 + t279 * t228;
t238 = t248 * qJD(3) + qJDD(1) * t295;
t194 = (-t248 * t261 - t238) * pkin(9) + (t249 * t261 - t237) * pkin(3) + t204;
t282 = sin(qJ(4));
t286 = cos(qJ(4));
t189 = t192 * t286 + t194 * t282;
t242 = -t249 * t282 + t261 * t286;
t243 = t249 * t286 + t261 * t282;
t222 = -pkin(4) * t242 - pkin(10) * t243;
t234 = qJDD(4) - t237;
t247 = qJD(4) - t248;
t246 = t247 ^ 2;
t186 = -pkin(4) * t246 + pkin(10) * t234 + t222 * t242 + t189;
t191 = -t258 * pkin(3) - t257 * pkin(9) + t249 * t236 - t195;
t214 = -qJD(4) * t243 - t238 * t282 + t258 * t286;
t215 = qJD(4) * t242 + t238 * t286 + t258 * t282;
t187 = (-t242 * t247 - t215) * pkin(10) + (t243 * t247 - t214) * pkin(4) + t191;
t281 = sin(qJ(5));
t285 = cos(qJ(5));
t183 = -t186 * t281 + t187 * t285;
t226 = -t243 * t281 + t247 * t285;
t199 = qJD(5) * t226 + t215 * t285 + t234 * t281;
t227 = t243 * t285 + t247 * t281;
t205 = -mrSges(6,1) * t226 + mrSges(6,2) * t227;
t239 = qJD(5) - t242;
t206 = -mrSges(6,2) * t239 + mrSges(6,3) * t226;
t212 = qJDD(5) - t214;
t180 = m(6) * t183 + mrSges(6,1) * t212 - mrSges(6,3) * t199 - t205 * t227 + t206 * t239;
t184 = t186 * t285 + t187 * t281;
t198 = -qJD(5) * t227 - t215 * t281 + t234 * t285;
t207 = mrSges(6,1) * t239 - mrSges(6,3) * t227;
t181 = m(6) * t184 - mrSges(6,2) * t212 + mrSges(6,3) * t198 + t205 * t226 - t207 * t239;
t174 = -t180 * t281 + t181 * t285;
t221 = -mrSges(5,1) * t242 + mrSges(5,2) * t243;
t230 = mrSges(5,1) * t247 - mrSges(5,3) * t243;
t172 = m(5) * t189 - mrSges(5,2) * t234 + mrSges(5,3) * t214 + t221 * t242 - t230 * t247 + t174;
t188 = -t192 * t282 + t194 * t286;
t185 = -pkin(4) * t234 - pkin(10) * t246 + t222 * t243 - t188;
t182 = -m(6) * t185 + mrSges(6,1) * t198 - mrSges(6,2) * t199 + t206 * t226 - t207 * t227;
t229 = -mrSges(5,2) * t247 + mrSges(5,3) * t242;
t178 = m(5) * t188 + mrSges(5,1) * t234 - mrSges(5,3) * t215 - t221 * t243 + t229 * t247 + t182;
t166 = t172 * t282 + t178 * t286;
t235 = -mrSges(4,1) * t248 + mrSges(4,2) * t249;
t245 = mrSges(4,1) * t261 - mrSges(4,3) * t249;
t310 = t172 * t286 - t178 * t282;
t163 = m(4) * t196 - mrSges(4,2) * t258 + mrSges(4,3) * t237 + t235 * t248 - t245 * t261 + t310;
t244 = -mrSges(4,2) * t261 + mrSges(4,3) * t248;
t165 = m(4) * t204 - mrSges(4,1) * t237 + mrSges(4,2) * t238 - t244 * t248 + t245 * t249 + t166;
t173 = t180 * t285 + t181 * t281;
t292 = -m(5) * t191 + mrSges(5,1) * t214 - mrSges(5,2) * t215 + t229 * t242 - t230 * t243 - t173;
t169 = m(4) * t195 + mrSges(4,1) * t258 - mrSges(4,3) * t238 - t235 * t249 + t244 * t261 + t292;
t152 = t163 * t323 + t165 * t279 + t169 * t322;
t153 = t163 * t318 - t276 * t165 + t169 * t317;
t240 = -t275 * t264 + t306;
t309 = -mrSges(3,1) * t278 + mrSges(3,2) * t275;
t262 = t309 * t316;
t303 = -mrSges(3,2) * t280 + mrSges(3,3) * t321;
t267 = t303 * qJD(1);
t304 = mrSges(3,1) * t280 - mrSges(3,3) * t325;
t150 = m(3) * t240 + t304 * qJDD(1) + (-t262 * t325 + t267 * t280) * qJD(1) + t153;
t157 = t163 * t287 - t283 * t169;
t241 = -g(3) * t325 + t313;
t266 = t304 * qJD(1);
t156 = m(3) * t241 + t303 * qJDD(1) + (t262 * t321 - t266 * t280) * qJD(1) + t157;
t145 = -t150 * t275 + t156 * t278;
t250 = -t277 * t263 + t311;
t151 = m(3) * t250 + (t309 * qJDD(1) + (t266 * t275 - t267 * t278) * qJD(1)) * t277 + t152;
t142 = t150 * t319 - t151 * t277 + t156 * t324;
t308 = Ifges(3,5) * t275 + Ifges(3,6) * t278;
t200 = Ifges(6,5) * t227 + Ifges(6,6) * t226 + Ifges(6,3) * t239;
t202 = Ifges(6,1) * t227 + Ifges(6,4) * t226 + Ifges(6,5) * t239;
t175 = -mrSges(6,1) * t185 + mrSges(6,3) * t184 + Ifges(6,4) * t199 + Ifges(6,2) * t198 + Ifges(6,6) * t212 - t200 * t227 + t202 * t239;
t201 = Ifges(6,4) * t227 + Ifges(6,2) * t226 + Ifges(6,6) * t239;
t176 = mrSges(6,2) * t185 - mrSges(6,3) * t183 + Ifges(6,1) * t199 + Ifges(6,4) * t198 + Ifges(6,5) * t212 + t200 * t226 - t201 * t239;
t208 = Ifges(5,5) * t243 + Ifges(5,6) * t242 + Ifges(5,3) * t247;
t209 = Ifges(5,4) * t243 + Ifges(5,2) * t242 + Ifges(5,6) * t247;
t158 = mrSges(5,2) * t191 - mrSges(5,3) * t188 + Ifges(5,1) * t215 + Ifges(5,4) * t214 + Ifges(5,5) * t234 - pkin(10) * t173 - t175 * t281 + t176 * t285 + t208 * t242 - t209 * t247;
t210 = Ifges(5,1) * t243 + Ifges(5,4) * t242 + Ifges(5,5) * t247;
t291 = mrSges(6,1) * t183 - mrSges(6,2) * t184 + Ifges(6,5) * t199 + Ifges(6,6) * t198 + Ifges(6,3) * t212 + t201 * t227 - t202 * t226;
t159 = -mrSges(5,1) * t191 + mrSges(5,3) * t189 + Ifges(5,4) * t215 + Ifges(5,2) * t214 + Ifges(5,6) * t234 - pkin(4) * t173 - t208 * t243 + t210 * t247 - t291;
t231 = Ifges(4,5) * t249 + Ifges(4,6) * t248 + Ifges(4,3) * t261;
t232 = Ifges(4,4) * t249 + Ifges(4,2) * t248 + Ifges(4,6) * t261;
t147 = mrSges(4,2) * t204 - mrSges(4,3) * t195 + Ifges(4,1) * t238 + Ifges(4,4) * t237 + Ifges(4,5) * t258 - pkin(9) * t166 + t158 * t286 - t159 * t282 + t231 * t248 - t232 * t261;
t233 = Ifges(4,1) * t249 + Ifges(4,4) * t248 + Ifges(4,5) * t261;
t290 = mrSges(5,1) * t188 - mrSges(5,2) * t189 + Ifges(5,5) * t215 + Ifges(5,6) * t214 + Ifges(5,3) * t234 + pkin(4) * t182 + pkin(10) * t174 + t175 * t285 + t176 * t281 + t209 * t243 - t210 * t242;
t148 = -mrSges(4,1) * t204 + mrSges(4,3) * t196 + Ifges(4,4) * t238 + Ifges(4,2) * t237 + Ifges(4,6) * t258 - pkin(3) * t166 - t231 * t249 + t233 * t261 - t290;
t299 = pkin(8) * t157 + t147 * t283 + t148 * t287;
t298 = Ifges(3,5) * t280 + (Ifges(3,1) * t275 + Ifges(3,4) * t278) * t277;
t297 = Ifges(3,6) * t280 + (Ifges(3,4) * t275 + Ifges(3,2) * t278) * t277;
t146 = mrSges(4,1) * t195 - mrSges(4,2) * t196 + Ifges(4,5) * t238 + Ifges(4,6) * t237 + Ifges(4,3) * t258 + pkin(3) * t292 + pkin(9) * t310 + t282 * t158 + t286 * t159 + t249 * t232 - t248 * t233;
t254 = t297 * qJD(1);
t255 = t298 * qJD(1);
t134 = qJDD(1) * t328 + mrSges(3,1) * t240 - mrSges(3,2) * t241 + pkin(2) * t153 + t279 * t146 + t299 * t276 + (t308 * qJDD(1) + (t254 * t275 - t255 * t278) * qJD(1)) * t277;
t253 = (t277 * t308 + t328) * qJD(1);
t136 = -mrSges(3,1) * t250 + mrSges(3,3) * t241 - pkin(2) * t152 - t276 * t146 + (-t253 * t325 + t255 * t280) * qJD(1) + t299 * t279 + t297 * qJDD(1);
t138 = mrSges(3,2) * t250 - mrSges(3,3) * t240 + t287 * t147 - t283 * t148 + (t253 * t321 - t254 * t280) * qJD(1) + (-t152 * t276 - t153 * t279) * pkin(8) + t298 * qJDD(1);
t296 = mrSges(2,1) * t271 - mrSges(2,2) * t272 + Ifges(2,3) * qJDD(1) + pkin(1) * t142 + t134 * t280 + t136 * t321 + t138 * t325 + t145 * t326;
t143 = m(2) * t272 - mrSges(2,1) * t289 - qJDD(1) * mrSges(2,2) + t145;
t141 = t280 * t151 + (t150 * t278 + t156 * t275) * t277;
t139 = m(2) * t271 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t289 + t142;
t132 = -mrSges(2,2) * g(3) - mrSges(2,3) * t271 + Ifges(2,5) * qJDD(1) - t289 * Ifges(2,6) - t275 * t136 + t278 * t138 + (-t141 * t277 - t142 * t280) * qJ(2);
t131 = mrSges(2,1) * g(3) + mrSges(2,3) * t272 + t289 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t141 - t277 * t134 + (qJ(2) * t145 + t136 * t278 + t138 * t275) * t280;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t288 * t132 - t284 * t131 - pkin(7) * (t139 * t288 + t143 * t284), t132, t138, t147, t158, t176; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t284 * t132 + t288 * t131 + pkin(7) * (-t139 * t284 + t143 * t288), t131, t136, t148, t159, t175; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t296, t296, t134, t146, t290, t291;];
m_new = t1;
