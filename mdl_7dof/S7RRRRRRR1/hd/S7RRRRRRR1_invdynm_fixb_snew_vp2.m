% Calculate vector of cutting torques with Newton-Euler for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% qJDD [7x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% m_mdh [8x1]
%   mass of all robot links (including the base)
% mrSges [8x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [8x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x8]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-09 00:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S7RRRRRRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(7,1),zeros(3,1),zeros(4,1),zeros(8,1),zeros(8,3),zeros(8,6)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_invdynm_fixb_snew_vp2: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_invdynm_fixb_snew_vp2: qJD has to be [7x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [7 1]), ...
  'S7RRRRRRR1_invdynm_fixb_snew_vp2: qJDD has to be [7x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S7RRRRRRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_invdynm_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [8 1]), ...
  'S7RRRRRRR1_invdynm_fixb_snew_vp2: m has to be [8x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [8,3]), ...
  'S7RRRRRRR1_invdynm_fixb_snew_vp2: mrSges has to be [8x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [8 6]), ...
  'S7RRRRRRR1_invdynm_fixb_snew_vp2: Ifges has to be [8x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 21:40:55
% EndTime: 2019-05-08 21:42:23
% DurationCPUTime: 30.93s
% Computational Cost: add. (569027->397), mult. (1031315->502), div. (0->0), fcn. (854424->14), ass. (0->150)
t295 = sin(qJ(2));
t302 = cos(qJ(2));
t318 = qJD(1) * qJD(2);
t317 = t302 * t318;
t279 = qJDD(1) * t295 + t317;
t296 = sin(qJ(1));
t303 = cos(qJ(1));
t284 = t296 * g(1) - g(2) * t303;
t259 = (t279 + t317) * pkin(2) - t284;
t285 = -g(1) * t303 - g(2) * t296;
t268 = -t295 * g(3) + t285 * t302;
t304 = qJD(1) ^ 2;
t260 = (t295 * t302 * t304 - qJDD(2)) * pkin(2) + t268;
t294 = sin(qJ(3));
t301 = cos(qJ(3));
t240 = -t259 * t294 + t260 * t301;
t267 = -t302 * g(3) - t295 * t285;
t264 = (-t295 ^ 2 * t304 - qJD(2) ^ 2) * pkin(2) + t267;
t293 = sin(qJ(4));
t300 = cos(qJ(4));
t222 = t240 * t300 - t293 * t264;
t320 = qJD(1) * t295;
t276 = -qJD(2) * t294 + t301 * t320;
t255 = -qJD(3) * t276 - qJDD(2) * t301 - t279 * t294;
t254 = qJDD(4) + t255;
t319 = qJD(1) * t302;
t286 = qJD(3) + t319;
t262 = -t276 * t293 - t286 * t300;
t263 = t276 * t300 - t286 * t293;
t203 = (-t262 * t263 + t254) * pkin(3) + t222;
t292 = sin(qJ(5));
t299 = cos(qJ(5));
t275 = -qJD(2) * t301 - t294 * t320;
t256 = qJD(3) * t275 - qJDD(2) * t294 + t279 * t301;
t280 = qJDD(1) * t302 - t295 * t318;
t274 = qJDD(3) + t280;
t229 = qJD(4) * t262 + t256 * t300 - t274 * t293;
t239 = -t301 * t259 - t294 * t260;
t273 = qJD(4) + t275;
t310 = (-t262 * t273 - t229) * pkin(3) + t239;
t188 = t299 * t203 + t292 * t310;
t221 = -t293 * t240 - t300 * t264;
t210 = (-t263 ^ 2 - t273 ^ 2) * pkin(3) - t221;
t291 = sin(qJ(6));
t298 = cos(qJ(6));
t184 = t188 * t298 + t210 * t291;
t245 = t263 * t299 + t273 * t292;
t205 = -qJD(5) * t245 - t229 * t292 + t254 * t299;
t204 = qJDD(6) - t205;
t261 = qJD(5) - t262;
t230 = -t245 * t291 + t261 * t298;
t231 = t245 * t298 + t261 * t291;
t176 = (t230 * t231 - t204) * pkin(4) + t184;
t187 = t292 * t203 - t299 * t310;
t244 = -t263 * t292 + t273 * t299;
t206 = qJD(5) * t244 + t229 * t299 + t254 * t292;
t228 = -qJD(4) * t263 - t256 * t293 - t274 * t300;
t226 = qJDD(5) - t228;
t195 = qJD(6) * t230 + t206 * t298 + t226 * t291;
t243 = qJD(6) - t244;
t177 = (t230 * t243 + t195) * pkin(4) + t187;
t290 = sin(qJ(7));
t297 = cos(qJ(7));
t175 = t176 * t297 - t177 * t290;
t212 = t231 * t297 - t243 * t290;
t180 = -qJD(7) * t212 - t195 * t290 - t204 * t297;
t211 = -t231 * t290 - t243 * t297;
t181 = qJD(7) * t211 + t195 * t297 - t204 * t290;
t183 = -t291 * t188 + t210 * t298;
t182 = (-t231 ^ 2 - t243 ^ 2) * pkin(4) + t183;
t227 = qJD(7) + t230;
t189 = Ifges(8,5) * t212 + Ifges(8,6) * t211 + Ifges(8,3) * t227;
t191 = Ifges(8,1) * t212 + Ifges(8,4) * t211 + Ifges(8,5) * t227;
t194 = -qJD(6) * t231 - t206 * t291 + t226 * t298;
t193 = qJDD(7) + t194;
t167 = -mrSges(8,1) * t182 + mrSges(8,3) * t175 + Ifges(8,4) * t181 + Ifges(8,2) * t180 + Ifges(8,6) * t193 - t189 * t212 + t191 * t227;
t174 = -t176 * t290 - t177 * t297;
t190 = Ifges(8,4) * t212 + Ifges(8,2) * t211 + Ifges(8,6) * t227;
t168 = mrSges(8,2) * t182 - mrSges(8,3) * t174 + Ifges(8,1) * t181 + Ifges(8,4) * t180 + Ifges(8,5) * t193 + t189 * t211 - t190 * t227;
t199 = Ifges(7,5) * t231 + Ifges(7,6) * t230 + Ifges(7,3) * t243;
t200 = Ifges(7,4) * t231 + Ifges(7,2) * t230 + Ifges(7,6) * t243;
t196 = -mrSges(8,1) * t211 + mrSges(8,2) * t212;
t197 = -mrSges(8,2) * t227 + mrSges(8,3) * t211;
t171 = m(8) * t174 + mrSges(8,1) * t193 - mrSges(8,3) * t181 - t196 * t212 + t197 * t227;
t198 = mrSges(8,1) * t227 - mrSges(8,3) * t212;
t172 = m(8) * t175 - mrSges(8,2) * t193 + mrSges(8,3) * t180 + t196 * t211 - t198 * t227;
t314 = -t297 * t171 - t290 * t172;
t157 = mrSges(7,2) * t187 - mrSges(7,3) * t183 + Ifges(7,1) * t195 + Ifges(7,4) * t194 + Ifges(7,5) * t204 + pkin(4) * t314 - t290 * t167 + t297 * t168 + t230 * t199 - t243 * t200;
t201 = Ifges(7,1) * t231 + Ifges(7,4) * t230 + Ifges(7,5) * t243;
t311 = mrSges(8,1) * t174 - mrSges(8,2) * t175 + Ifges(8,5) * t181 + Ifges(8,6) * t180 + Ifges(8,3) * t193 + t190 * t212 - t191 * t211;
t164 = -mrSges(7,1) * t187 + mrSges(7,3) * t184 + Ifges(7,4) * t195 + Ifges(7,2) * t194 + Ifges(7,6) * t204 - t199 * t231 + t201 * t243 + t311;
t215 = Ifges(6,5) * t245 + Ifges(6,6) * t244 + Ifges(6,3) * t261;
t216 = Ifges(6,4) * t245 + Ifges(6,2) * t244 + Ifges(6,6) * t261;
t151 = mrSges(6,2) * t210 + mrSges(6,3) * t187 + Ifges(6,1) * t206 + Ifges(6,4) * t205 + Ifges(6,5) * t226 + t157 * t298 - t164 * t291 + t215 * t244 - t216 * t261;
t217 = Ifges(6,1) * t245 + Ifges(6,4) * t244 + Ifges(6,5) * t261;
t166 = -t171 * t290 + t172 * t297;
t305 = mrSges(7,1) * t183 - mrSges(7,2) * t184 + Ifges(7,5) * t195 + Ifges(7,6) * t194 + Ifges(7,3) * t204 - pkin(4) * t166 - t167 * t297 - t168 * t290 + t200 * t231 - t201 * t230;
t156 = -mrSges(6,1) * t210 + mrSges(6,3) * t188 + Ifges(6,4) * t206 + Ifges(6,2) * t205 + Ifges(6,6) * t226 - t215 * t245 + t217 * t261 - t305;
t234 = Ifges(5,5) * t263 + Ifges(5,6) * t262 + Ifges(5,3) * t273;
t235 = Ifges(5,4) * t263 + Ifges(5,2) * t262 + Ifges(5,6) * t273;
t207 = -mrSges(7,1) * t230 + mrSges(7,2) * t231;
t214 = mrSges(7,1) * t243 - mrSges(7,3) * t231;
t165 = m(7) * t184 - mrSges(7,2) * t204 + mrSges(7,3) * t194 + t207 * t230 - t214 * t243 + t166;
t213 = -mrSges(7,2) * t243 + mrSges(7,3) * t230;
t169 = m(7) * t183 + m(8) * t182 + mrSges(7,1) * t204 - mrSges(8,1) * t180 + mrSges(8,2) * t181 - mrSges(7,3) * t195 - t197 * t211 + t198 * t212 - t207 * t231 + t213 * t243;
t219 = -mrSges(6,1) * t244 + mrSges(6,2) * t245;
t233 = mrSges(6,1) * t261 - mrSges(6,3) * t245;
t161 = m(6) * t188 - mrSges(6,2) * t226 + mrSges(6,3) * t205 + t165 * t298 - t169 * t291 + t219 * t244 - t233 * t261;
t232 = -mrSges(6,2) * t261 + mrSges(6,3) * t244;
t163 = t226 * mrSges(6,1) + t194 * mrSges(7,1) - t195 * mrSges(7,2) - t206 * mrSges(6,3) + t230 * t213 - t231 * t214 - t245 * t219 + t261 * t232 + (-m(6) - m(7)) * t187 - t314;
t321 = t161 * t292 + t163 * t299;
t144 = mrSges(5,2) * t239 - mrSges(5,3) * t221 + Ifges(5,1) * t229 + Ifges(5,4) * t228 + Ifges(5,5) * t254 - pkin(3) * t321 + t151 * t299 - t156 * t292 + t234 * t262 - t235 * t273;
t236 = Ifges(5,1) * t263 + Ifges(5,4) * t262 + Ifges(5,5) * t273;
t306 = mrSges(6,1) * t187 + mrSges(6,2) * t188 - Ifges(6,5) * t206 - Ifges(6,6) * t205 - Ifges(6,3) * t226 - t157 * t291 - t164 * t298 - t216 * t245 + t217 * t244;
t149 = -mrSges(5,1) * t239 + mrSges(5,3) * t222 + Ifges(5,4) * t229 + Ifges(5,2) * t228 + Ifges(5,6) * t254 - t234 * t263 + t236 * t273 + t306;
t249 = Ifges(4,5) * t276 + Ifges(4,6) * t275 + Ifges(4,3) * t286;
t250 = Ifges(4,4) * t276 + Ifges(4,2) * t275 + Ifges(4,6) * t286;
t141 = mrSges(4,2) * t264 - mrSges(4,3) * t239 + Ifges(4,1) * t256 + Ifges(4,4) * t255 + Ifges(4,5) * t274 + t144 * t300 - t149 * t293 + t249 * t275 - t250 * t286;
t251 = Ifges(4,1) * t276 + Ifges(4,4) * t275 + Ifges(4,5) * t286;
t316 = t161 * t299 - t163 * t292;
t309 = mrSges(5,1) * t221 - mrSges(5,2) * t222 + Ifges(5,5) * t229 + Ifges(5,6) * t228 + Ifges(5,3) * t254 + pkin(3) * t316 + t151 * t292 + t156 * t299 + t235 * t263 - t236 * t262;
t143 = -mrSges(4,1) * t264 + mrSges(4,3) * t240 + Ifges(4,4) * t256 + Ifges(4,2) * t255 + Ifges(4,6) * t274 - t249 * t276 + t251 * t286 + t309;
t241 = -mrSges(5,1) * t262 + mrSges(5,2) * t263;
t247 = mrSges(5,1) * t273 - mrSges(5,3) * t263;
t153 = m(5) * t222 - mrSges(5,2) * t254 + mrSges(5,3) * t228 + t241 * t262 - t247 * t273 + t316;
t246 = -mrSges(5,2) * t273 + mrSges(5,3) * t262;
t158 = m(5) * t221 - m(6) * t210 + mrSges(5,1) * t254 + mrSges(6,1) * t205 - mrSges(6,2) * t206 - mrSges(5,3) * t229 - t165 * t291 - t169 * t298 + t232 * t244 - t233 * t245 - t241 * t263 + t246 * t273;
t258 = -mrSges(4,1) * t275 + mrSges(4,2) * t276;
t266 = mrSges(4,1) * t286 - mrSges(4,3) * t276;
t148 = m(4) * t240 - mrSges(4,2) * t274 + mrSges(4,3) * t255 + t153 * t300 - t158 * t293 + t258 * t275 - t266 * t286;
t265 = -mrSges(4,2) * t286 + mrSges(4,3) * t275;
t152 = t274 * mrSges(4,1) - t228 * mrSges(5,1) + t229 * mrSges(5,2) - t256 * mrSges(4,3) - t262 * t246 + t263 * t247 - t276 * t258 + t286 * t265 + (m(4) + m(5)) * t239 + t321;
t146 = t148 * t301 - t152 * t294;
t270 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t295 + Ifges(3,2) * t302) * qJD(1);
t271 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t295 + Ifges(3,4) * t302) * qJD(1);
t322 = mrSges(3,1) * t267 - mrSges(3,2) * t268 + Ifges(3,5) * t279 + Ifges(3,6) * t280 + Ifges(3,3) * qJDD(2) - pkin(2) * t146 - t294 * t141 - t301 * t143 + (t270 * t295 - t271 * t302) * qJD(1);
t315 = -t294 * t148 - t301 * t152;
t269 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t295 + Ifges(3,6) * t302) * qJD(1);
t138 = -mrSges(3,2) * t284 - mrSges(3,3) * t267 + Ifges(3,1) * t279 + Ifges(3,4) * t280 + Ifges(3,5) * qJDD(2) + pkin(2) * t315 - qJD(2) * t270 + t301 * t141 - t294 * t143 + t269 * t319;
t308 = mrSges(4,1) * t239 - mrSges(4,2) * t240 + Ifges(4,5) * t256 + Ifges(4,6) * t255 + Ifges(4,3) * t274 - t144 * t293 - t149 * t300 + t250 * t276 - t251 * t275;
t140 = mrSges(3,1) * t284 + mrSges(3,3) * t268 + Ifges(3,4) * t279 + Ifges(3,2) * t280 + Ifges(3,6) * qJDD(2) + qJD(2) * t271 - t269 * t320 + t308;
t312 = mrSges(2,1) * t284 - mrSges(2,2) * t285 + Ifges(2,3) * qJDD(1) + t138 * t295 + t140 * t302;
t283 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t319;
t282 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t320;
t278 = (-mrSges(3,1) * t302 + mrSges(3,2) * t295) * qJD(1);
t145 = qJDD(1) * mrSges(2,1) + t280 * mrSges(3,1) - t304 * mrSges(2,2) - t279 * mrSges(3,2) + (m(2) + m(3)) * t284 + (-t282 * t295 + t283 * t302) * qJD(1) - t315;
t142 = m(2) * t285 - qJDD(1) * mrSges(2,2) - t304 * mrSges(2,1) + t302 * (m(3) * t268 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t280 - qJD(2) * t282 + t278 * t319 + t146) - t295 * (m(3) * t267 + m(4) * t264 + qJDD(2) * mrSges(3,1) - mrSges(4,1) * t255 + mrSges(4,2) * t256 - mrSges(3,3) * t279 + qJD(2) * t283 - t153 * t293 - t158 * t300 - t265 * t275 + t266 * t276 - t278 * t320);
t136 = mrSges(2,1) * g(3) + mrSges(2,3) * t285 + t304 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - t322;
t135 = -mrSges(2,2) * g(3) - mrSges(2,3) * t284 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t304 + t138 * t302 - t140 * t295;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t303 * t135 - t296 * t136 - pkin(1) * (t142 * t296 + t145 * t303), t135, t138, t141, t144, t151, t157, t168; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t296 * t135 + t303 * t136 + pkin(1) * (t142 * t303 - t145 * t296), t136, t140, t143, t149, t156, t164, t167; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t312, t312, t322, t308, t309, -t306, t305, t311;];
m_new  = t1;
