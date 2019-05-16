% Calculate vector of cutting torques with Newton-Euler for
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-05-04 20:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PPRRRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:39:10
% EndTime: 2019-05-04 20:39:41
% DurationCPUTime: 30.81s
% Computational Cost: add. (625080->292), mult. (1153039->387), div. (0->0), fcn. (905697->16), ass. (0->136)
t269 = sin(pkin(12));
t273 = cos(pkin(12));
t258 = -t273 * g(1) - t269 * g(2);
t268 = sin(pkin(13));
t272 = cos(pkin(13));
t257 = t269 * g(1) - t273 * g(2);
t267 = -g(3) + qJDD(1);
t271 = sin(pkin(6));
t275 = cos(pkin(6));
t295 = t257 * t275 + t267 * t271;
t230 = -t268 * t258 + t295 * t272;
t231 = t272 * t258 + t295 * t268;
t243 = -t271 * t257 + t275 * t267 + qJDD(2);
t283 = cos(qJ(3));
t274 = cos(pkin(7));
t279 = sin(qJ(3));
t305 = t274 * t279;
t270 = sin(pkin(7));
t306 = t270 * t279;
t206 = t230 * t305 + t283 * t231 + t243 * t306;
t284 = qJD(3) ^ 2;
t204 = -t284 * pkin(3) + qJDD(3) * pkin(9) + t206;
t216 = -t270 * t230 + t274 * t243;
t278 = sin(qJ(4));
t282 = cos(qJ(4));
t199 = -t278 * t204 + t282 * t216;
t302 = qJD(3) * qJD(4);
t301 = t282 * t302;
t254 = t278 * qJDD(3) + t301;
t197 = (-t254 + t301) * pkin(10) + (t278 * t282 * t284 + qJDD(4)) * pkin(4) + t199;
t200 = t282 * t204 + t278 * t216;
t255 = t282 * qJDD(3) - t278 * t302;
t304 = qJD(3) * t278;
t261 = qJD(4) * pkin(4) - pkin(10) * t304;
t266 = t282 ^ 2;
t198 = -t266 * t284 * pkin(4) + t255 * pkin(10) - qJD(4) * t261 + t200;
t277 = sin(qJ(5));
t281 = cos(qJ(5));
t193 = t277 * t197 + t281 * t198;
t249 = (t277 * t282 + t278 * t281) * qJD(3);
t223 = -t249 * qJD(5) - t277 * t254 + t281 * t255;
t248 = (t277 * t278 - t281 * t282) * qJD(3);
t236 = t248 * mrSges(6,1) + t249 * mrSges(6,2);
t265 = qJD(4) + qJD(5);
t242 = t265 * mrSges(6,1) - t249 * mrSges(6,3);
t264 = qJDD(4) + qJDD(5);
t237 = t248 * pkin(5) - t249 * pkin(11);
t263 = t265 ^ 2;
t190 = -t263 * pkin(5) + t264 * pkin(11) - t248 * t237 + t193;
t205 = -t279 * t231 + (t230 * t274 + t243 * t270) * t283;
t290 = -qJDD(3) * pkin(3) - t205;
t201 = -t255 * pkin(4) + t261 * t304 + (-pkin(10) * t266 - pkin(9)) * t284 + t290;
t224 = -t248 * qJD(5) + t281 * t254 + t277 * t255;
t194 = (t248 * t265 - t224) * pkin(11) + (t249 * t265 - t223) * pkin(5) + t201;
t276 = sin(qJ(6));
t280 = cos(qJ(6));
t187 = -t276 * t190 + t280 * t194;
t239 = -t276 * t249 + t280 * t265;
t209 = t239 * qJD(6) + t280 * t224 + t276 * t264;
t240 = t280 * t249 + t276 * t265;
t217 = -t239 * mrSges(7,1) + t240 * mrSges(7,2);
t222 = qJDD(6) - t223;
t244 = qJD(6) + t248;
t225 = -t244 * mrSges(7,2) + t239 * mrSges(7,3);
t183 = m(7) * t187 + t222 * mrSges(7,1) - t209 * mrSges(7,3) - t240 * t217 + t244 * t225;
t188 = t280 * t190 + t276 * t194;
t208 = -t240 * qJD(6) - t276 * t224 + t280 * t264;
t226 = t244 * mrSges(7,1) - t240 * mrSges(7,3);
t184 = m(7) * t188 - t222 * mrSges(7,2) + t208 * mrSges(7,3) + t239 * t217 - t244 * t226;
t298 = -t276 * t183 + t280 * t184;
t170 = m(6) * t193 - t264 * mrSges(6,2) + t223 * mrSges(6,3) - t248 * t236 - t265 * t242 + t298;
t192 = t281 * t197 - t277 * t198;
t241 = -t265 * mrSges(6,2) - t248 * mrSges(6,3);
t189 = -t264 * pkin(5) - t263 * pkin(11) + t249 * t237 - t192;
t291 = -m(7) * t189 + t208 * mrSges(7,1) - t209 * mrSges(7,2) + t239 * t225 - t240 * t226;
t179 = m(6) * t192 + t264 * mrSges(6,1) - t224 * mrSges(6,3) - t249 * t236 + t265 * t241 + t291;
t164 = t277 * t170 + t281 * t179;
t253 = (-mrSges(5,1) * t282 + mrSges(5,2) * t278) * qJD(3);
t303 = qJD(3) * t282;
t260 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t303;
t162 = m(5) * t199 + qJDD(4) * mrSges(5,1) - t254 * mrSges(5,3) + qJD(4) * t260 - t253 * t304 + t164;
t259 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t304;
t299 = t281 * t170 - t277 * t179;
t163 = m(5) * t200 - qJDD(4) * mrSges(5,2) + t255 * mrSges(5,3) - qJD(4) * t259 + t253 * t303 + t299;
t300 = -t278 * t162 + t282 * t163;
t153 = m(4) * t206 - t284 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t300;
t156 = t282 * t162 + t278 * t163;
t155 = m(4) * t216 + t156;
t203 = -t284 * pkin(9) + t290;
t172 = t280 * t183 + t276 * t184;
t289 = m(6) * t201 - t223 * mrSges(6,1) + t224 * mrSges(6,2) + t248 * t241 + t249 * t242 + t172;
t286 = -m(5) * t203 + t255 * mrSges(5,1) - t254 * mrSges(5,2) - t259 * t304 + t260 * t303 - t289;
t167 = m(4) * t205 + qJDD(3) * mrSges(4,1) - t284 * mrSges(4,2) + t286;
t307 = t167 * t283;
t143 = t153 * t305 - t270 * t155 + t274 * t307;
t140 = m(3) * t230 + t143;
t149 = t283 * t153 - t279 * t167;
t148 = m(3) * t231 + t149;
t315 = t140 * t272 + t148 * t268;
t210 = Ifges(7,5) * t240 + Ifges(7,6) * t239 + Ifges(7,3) * t244;
t212 = Ifges(7,1) * t240 + Ifges(7,4) * t239 + Ifges(7,5) * t244;
t176 = -mrSges(7,1) * t189 + mrSges(7,3) * t188 + Ifges(7,4) * t209 + Ifges(7,2) * t208 + Ifges(7,6) * t222 - t240 * t210 + t244 * t212;
t211 = Ifges(7,4) * t240 + Ifges(7,2) * t239 + Ifges(7,6) * t244;
t177 = mrSges(7,2) * t189 - mrSges(7,3) * t187 + Ifges(7,1) * t209 + Ifges(7,4) * t208 + Ifges(7,5) * t222 + t239 * t210 - t244 * t211;
t232 = Ifges(6,5) * t249 - Ifges(6,6) * t248 + Ifges(6,3) * t265;
t233 = Ifges(6,4) * t249 - Ifges(6,2) * t248 + Ifges(6,6) * t265;
t157 = mrSges(6,2) * t201 - mrSges(6,3) * t192 + Ifges(6,1) * t224 + Ifges(6,4) * t223 + Ifges(6,5) * t264 - pkin(11) * t172 - t276 * t176 + t280 * t177 - t248 * t232 - t265 * t233;
t234 = Ifges(6,1) * t249 - Ifges(6,4) * t248 + Ifges(6,5) * t265;
t287 = mrSges(7,1) * t187 - mrSges(7,2) * t188 + Ifges(7,5) * t209 + Ifges(7,6) * t208 + Ifges(7,3) * t222 + t240 * t211 - t239 * t212;
t158 = -mrSges(6,1) * t201 + mrSges(6,3) * t193 + Ifges(6,4) * t224 + Ifges(6,2) * t223 + Ifges(6,6) * t264 - pkin(5) * t172 - t249 * t232 + t265 * t234 - t287;
t245 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t278 + Ifges(5,6) * t282) * qJD(3);
t247 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t278 + Ifges(5,4) * t282) * qJD(3);
t144 = -mrSges(5,1) * t203 + mrSges(5,3) * t200 + Ifges(5,4) * t254 + Ifges(5,2) * t255 + Ifges(5,6) * qJDD(4) - pkin(4) * t289 + pkin(10) * t299 + qJD(4) * t247 + t277 * t157 + t281 * t158 - t245 * t304;
t246 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t278 + Ifges(5,2) * t282) * qJD(3);
t145 = mrSges(5,2) * t203 - mrSges(5,3) * t199 + Ifges(5,1) * t254 + Ifges(5,4) * t255 + Ifges(5,5) * qJDD(4) - pkin(10) * t164 - qJD(4) * t246 + t281 * t157 - t277 * t158 + t245 * t303;
t133 = mrSges(4,1) * t205 - mrSges(4,2) * t206 + Ifges(4,3) * qJDD(3) + pkin(3) * t286 + pkin(9) * t300 + t282 * t144 + t278 * t145;
t142 = t153 * t306 + t274 * t155 + t270 * t307;
t134 = mrSges(4,2) * t216 - mrSges(4,3) * t205 + Ifges(4,5) * qJDD(3) - t284 * Ifges(4,6) - pkin(9) * t156 - t278 * t144 + t282 * t145;
t288 = -mrSges(6,1) * t192 + mrSges(6,2) * t193 - Ifges(6,5) * t224 - Ifges(6,6) * t223 - Ifges(6,3) * t264 - pkin(5) * t291 - pkin(11) * t298 - t280 * t176 - t276 * t177 - t249 * t233 - t248 * t234;
t313 = mrSges(5,1) * t199 - mrSges(5,2) * t200 + Ifges(5,5) * t254 + Ifges(5,6) * t255 + Ifges(5,3) * qJDD(4) + pkin(4) * t164 + (t278 * t246 - t282 * t247) * qJD(3) - t288;
t138 = -mrSges(4,1) * t216 + mrSges(4,3) * t206 + t284 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t156 - t313;
t293 = pkin(8) * t149 + t134 * t279 + t138 * t283;
t126 = -mrSges(3,1) * t243 + mrSges(3,3) * t231 - pkin(2) * t142 - t270 * t133 + t293 * t274;
t128 = mrSges(3,2) * t243 - mrSges(3,3) * t230 + t283 * t134 - t279 * t138 + (-t142 * t270 - t143 * t274) * pkin(8);
t137 = -t268 * t140 + t272 * t148;
t314 = qJ(2) * t137 + t126 * t272 + t128 * t268;
t141 = m(3) * t243 + t142;
t132 = -t271 * t141 + t315 * t275;
t124 = mrSges(3,1) * t230 - mrSges(3,2) * t231 + pkin(2) * t143 + t274 * t133 + t293 * t270;
t292 = mrSges(2,1) * t257 - mrSges(2,2) * t258 + pkin(1) * t132 + t275 * t124 + t314 * t271;
t135 = m(2) * t258 + t137;
t131 = t275 * t141 + t315 * t271;
t129 = m(2) * t257 + t132;
t122 = mrSges(2,2) * t267 - mrSges(2,3) * t257 - t268 * t126 + t272 * t128 + (-t131 * t271 - t132 * t275) * qJ(2);
t121 = -mrSges(2,1) * t267 + mrSges(2,3) * t258 - pkin(1) * t131 - t271 * t124 + t314 * t275;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t273 * t122 - t269 * t121 - qJ(1) * (t273 * t129 + t269 * t135), t122, t128, t134, t145, t157, t177; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t269 * t122 + t273 * t121 + qJ(1) * (-t269 * t129 + t273 * t135), t121, t126, t138, t144, t158, t176; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t292, t292, t124, t133, t313, -t288, t287;];
m_new  = t1;
