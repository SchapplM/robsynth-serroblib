% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-05-04 22:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:20:21
% EndTime: 2019-05-04 22:20:44
% DurationCPUTime: 17.06s
% Computational Cost: add. (328966->297), mult. (640120->384), div. (0->0), fcn. (451997->14), ass. (0->131)
t271 = sin(pkin(10));
t275 = cos(pkin(10));
t258 = t271 * g(1) - t275 * g(2);
t259 = -t275 * g(1) - t271 * g(2);
t268 = -g(3) + qJDD(1);
t279 = sin(qJ(2));
t276 = cos(pkin(6));
t282 = cos(qJ(2));
t301 = t276 * t282;
t272 = sin(pkin(6));
t303 = t272 * t282;
t218 = t258 * t301 - t279 * t259 + t268 * t303;
t216 = qJDD(2) * pkin(2) + t218;
t302 = t276 * t279;
t304 = t272 * t279;
t219 = t258 * t302 + t282 * t259 + t268 * t304;
t284 = qJD(2) ^ 2;
t217 = -t284 * pkin(2) + t219;
t270 = sin(pkin(11));
t274 = cos(pkin(11));
t200 = t270 * t216 + t274 * t217;
t197 = -t284 * pkin(3) + qJDD(2) * pkin(8) + t200;
t236 = -t272 * t258 + t276 * t268;
t230 = qJDD(3) + t236;
t278 = sin(qJ(4));
t281 = cos(qJ(4));
t192 = t281 * t197 + t278 * t230;
t253 = (-pkin(4) * t281 - qJ(5) * t278) * qJD(2);
t283 = qJD(4) ^ 2;
t299 = t281 * qJD(2);
t187 = -t283 * pkin(4) + qJDD(4) * qJ(5) + t253 * t299 + t192;
t199 = t274 * t216 - t270 * t217;
t196 = -qJDD(2) * pkin(3) - t284 * pkin(8) - t199;
t298 = qJD(2) * qJD(4);
t296 = t281 * t298;
t255 = t278 * qJDD(2) + t296;
t266 = t278 * t298;
t256 = t281 * qJDD(2) - t266;
t190 = (-t255 - t296) * qJ(5) + (-t256 + t266) * pkin(4) + t196;
t269 = sin(pkin(12));
t273 = cos(pkin(12));
t300 = qJD(2) * t278;
t249 = t269 * qJD(4) + t273 * t300;
t182 = -0.2e1 * qJD(5) * t249 - t269 * t187 + t273 * t190;
t234 = t269 * qJDD(4) + t273 * t255;
t248 = t273 * qJD(4) - t269 * t300;
t180 = (-t248 * t299 - t234) * pkin(9) + (t248 * t249 - t256) * pkin(5) + t182;
t183 = 0.2e1 * qJD(5) * t248 + t273 * t187 + t269 * t190;
t233 = t273 * qJDD(4) - t269 * t255;
t235 = -pkin(5) * t299 - t249 * pkin(9);
t247 = t248 ^ 2;
t181 = -t247 * pkin(5) + t233 * pkin(9) + t235 * t299 + t183;
t277 = sin(qJ(6));
t280 = cos(qJ(6));
t179 = t277 * t180 + t280 * t181;
t191 = -t278 * t197 + t281 * t230;
t186 = -qJDD(4) * pkin(4) - t283 * qJ(5) + t253 * t300 + qJDD(5) - t191;
t184 = -t233 * pkin(5) - t247 * pkin(9) + t249 * t235 + t186;
t225 = t277 * t248 + t280 * t249;
t204 = -t225 * qJD(6) + t280 * t233 - t277 * t234;
t224 = t280 * t248 - t277 * t249;
t205 = t224 * qJD(6) + t277 * t233 + t280 * t234;
t264 = qJD(6) - t299;
t206 = Ifges(7,5) * t225 + Ifges(7,6) * t224 + Ifges(7,3) * t264;
t208 = Ifges(7,1) * t225 + Ifges(7,4) * t224 + Ifges(7,5) * t264;
t251 = qJDD(6) - t256;
t167 = -mrSges(7,1) * t184 + mrSges(7,3) * t179 + Ifges(7,4) * t205 + Ifges(7,2) * t204 + Ifges(7,6) * t251 - t225 * t206 + t264 * t208;
t178 = t280 * t180 - t277 * t181;
t207 = Ifges(7,4) * t225 + Ifges(7,2) * t224 + Ifges(7,6) * t264;
t168 = mrSges(7,2) * t184 - mrSges(7,3) * t178 + Ifges(7,1) * t205 + Ifges(7,4) * t204 + Ifges(7,5) * t251 + t224 * t206 - t264 * t207;
t220 = Ifges(6,5) * t249 + Ifges(6,6) * t248 - Ifges(6,3) * t299;
t222 = Ifges(6,1) * t249 + Ifges(6,4) * t248 - Ifges(6,5) * t299;
t214 = -t264 * mrSges(7,2) + t224 * mrSges(7,3);
t215 = t264 * mrSges(7,1) - t225 * mrSges(7,3);
t289 = m(7) * t184 - t204 * mrSges(7,1) + t205 * mrSges(7,2) - t224 * t214 + t225 * t215;
t210 = -t224 * mrSges(7,1) + t225 * mrSges(7,2);
t172 = m(7) * t178 + t251 * mrSges(7,1) - t205 * mrSges(7,3) - t225 * t210 + t264 * t214;
t173 = m(7) * t179 - t251 * mrSges(7,2) + t204 * mrSges(7,3) + t224 * t210 - t264 * t215;
t293 = -t277 * t172 + t280 * t173;
t153 = -mrSges(6,1) * t186 + mrSges(6,3) * t183 + Ifges(6,4) * t234 + Ifges(6,2) * t233 - Ifges(6,6) * t256 - pkin(5) * t289 + pkin(9) * t293 + t280 * t167 + t277 * t168 - t249 * t220 - t222 * t299;
t166 = t280 * t172 + t277 * t173;
t221 = Ifges(6,4) * t249 + Ifges(6,2) * t248 - Ifges(6,6) * t299;
t154 = mrSges(6,2) * t186 - mrSges(6,3) * t182 + Ifges(6,1) * t234 + Ifges(6,4) * t233 - Ifges(6,5) * t256 - pkin(9) * t166 - t277 * t167 + t280 * t168 + t248 * t220 + t221 * t299;
t226 = -t248 * mrSges(6,1) + t249 * mrSges(6,2);
t231 = mrSges(6,2) * t299 + t248 * mrSges(6,3);
t164 = m(6) * t182 - t256 * mrSges(6,1) - t234 * mrSges(6,3) - t249 * t226 - t231 * t299 + t166;
t232 = -mrSges(6,1) * t299 - t249 * mrSges(6,3);
t165 = m(6) * t183 + t256 * mrSges(6,2) + t233 * mrSges(6,3) + t248 * t226 + t232 * t299 + t293;
t162 = -t269 * t164 + t273 * t165;
t176 = -m(6) * t186 + t233 * mrSges(6,1) - t234 * mrSges(6,2) + t248 * t231 - t249 * t232 - t289;
t243 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t278 + Ifges(5,2) * t281) * qJD(2);
t244 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t278 + Ifges(5,4) * t281) * qJD(2);
t306 = mrSges(5,1) * t191 - mrSges(5,2) * t192 + Ifges(5,5) * t255 + Ifges(5,6) * t256 + Ifges(5,3) * qJDD(4) + pkin(4) * t176 + qJ(5) * t162 + t273 * t153 + t269 * t154 + (t278 * t243 - t281 * t244) * qJD(2);
t254 = (-mrSges(5,1) * t281 + mrSges(5,2) * t278) * qJD(2);
t260 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t300;
t159 = m(5) * t192 - qJDD(4) * mrSges(5,2) + t256 * mrSges(5,3) - qJD(4) * t260 + t254 * t299 + t162;
t261 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t299;
t175 = m(5) * t191 + qJDD(4) * mrSges(5,1) - t255 * mrSges(5,3) + qJD(4) * t261 - t254 * t300 + t176;
t294 = t281 * t159 - t278 * t175;
t149 = m(4) * t200 - t284 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t294;
t161 = t273 * t164 + t269 * t165;
t287 = -m(5) * t196 + t256 * mrSges(5,1) - t255 * mrSges(5,2) - t260 * t300 + t261 * t299 - t161;
t156 = m(4) * t199 + qJDD(2) * mrSges(4,1) - t284 * mrSges(4,2) + t287;
t144 = t270 * t149 + t274 * t156;
t142 = m(3) * t218 + qJDD(2) * mrSges(3,1) - t284 * mrSges(3,2) + t144;
t295 = t274 * t149 - t270 * t156;
t143 = m(3) * t219 - t284 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t295;
t136 = -t279 * t142 + t282 * t143;
t305 = pkin(7) * t136;
t152 = t278 * t159 + t281 * t175;
t297 = m(4) * t230 + t152;
t150 = m(3) * t236 + t297;
t132 = t142 * t301 + t143 * t302 - t272 * t150;
t242 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t278 + Ifges(5,6) * t281) * qJD(2);
t138 = mrSges(5,2) * t196 - mrSges(5,3) * t191 + Ifges(5,1) * t255 + Ifges(5,4) * t256 + Ifges(5,5) * qJDD(4) - qJ(5) * t161 - qJD(4) * t243 - t269 * t153 + t273 * t154 + t242 * t299;
t288 = -mrSges(7,1) * t178 + mrSges(7,2) * t179 - Ifges(7,5) * t205 - Ifges(7,6) * t204 - Ifges(7,3) * t251 - t225 * t207 + t224 * t208;
t285 = -mrSges(6,1) * t182 + mrSges(6,2) * t183 - Ifges(6,5) * t234 - Ifges(6,6) * t233 - pkin(5) * t166 - t249 * t221 + t248 * t222 + t288;
t146 = -t242 * t300 + Ifges(5,6) * qJDD(4) + (Ifges(5,2) + Ifges(6,3)) * t256 + Ifges(5,4) * t255 + qJD(4) * t244 + mrSges(5,3) * t192 - mrSges(5,1) * t196 - pkin(4) * t161 + t285;
t128 = mrSges(4,2) * t230 - mrSges(4,3) * t199 + Ifges(4,5) * qJDD(2) - t284 * Ifges(4,6) - pkin(8) * t152 + t281 * t138 - t278 * t146;
t133 = -mrSges(4,1) * t230 + mrSges(4,3) * t200 + t284 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t152 - t306;
t123 = -mrSges(3,1) * t236 + mrSges(3,3) * t219 + t284 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t297 + qJ(3) * t295 + t270 * t128 + t274 * t133;
t125 = mrSges(3,2) * t236 - mrSges(3,3) * t218 + Ifges(3,5) * qJDD(2) - t284 * Ifges(3,6) - qJ(3) * t144 + t274 * t128 - t270 * t133;
t290 = mrSges(4,1) * t199 - mrSges(4,2) * t200 + Ifges(4,3) * qJDD(2) + pkin(3) * t287 + pkin(8) * t294 + t278 * t138 + t281 * t146;
t127 = mrSges(3,1) * t218 - mrSges(3,2) * t219 + Ifges(3,3) * qJDD(2) + pkin(2) * t144 + t290;
t291 = mrSges(2,1) * t258 - mrSges(2,2) * t259 + pkin(1) * t132 + t123 * t303 + t125 * t304 + t276 * t127 + t272 * t305;
t134 = m(2) * t259 + t136;
t131 = t276 * t150 + (t142 * t282 + t143 * t279) * t272;
t129 = m(2) * t258 + t132;
t121 = mrSges(2,2) * t268 - mrSges(2,3) * t258 - t279 * t123 + t282 * t125 + (-t131 * t272 - t132 * t276) * pkin(7);
t120 = -mrSges(2,1) * t268 + mrSges(2,3) * t259 - pkin(1) * t131 - t272 * t127 + (t123 * t282 + t125 * t279 + t305) * t276;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t275 * t121 - t271 * t120 - qJ(1) * (t275 * t129 + t271 * t134), t121, t125, t128, t138, t154, t168; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t271 * t121 + t275 * t120 + qJ(1) * (-t271 * t129 + t275 * t134), t120, t123, t133, t146, t153, t167; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t291, t291, t127, t290, t306, -Ifges(6,3) * t256 - t285, -t288;];
m_new  = t1;
