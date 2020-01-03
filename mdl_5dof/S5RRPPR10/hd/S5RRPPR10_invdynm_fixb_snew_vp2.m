% Calculate vector of cutting torques with Newton-Euler for
% S5RRPPR10
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPPR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR10_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR10_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR10_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:43:07
% EndTime: 2019-12-31 19:43:16
% DurationCPUTime: 4.46s
% Computational Cost: add. (47453->310), mult. (103719->373), div. (0->0), fcn. (64023->8), ass. (0->116)
t303 = -2 * qJD(3);
t266 = sin(qJ(1));
t269 = cos(qJ(1));
t252 = t266 * g(1) - t269 * g(2);
t271 = qJD(1) ^ 2;
t231 = -qJDD(1) * pkin(1) - t271 * pkin(6) - t252;
t265 = sin(qJ(2));
t268 = cos(qJ(2));
t289 = qJD(1) * qJD(2);
t287 = t268 * t289;
t247 = t265 * qJDD(1) + t287;
t256 = t265 * t289;
t248 = t268 * qJDD(1) - t256;
t185 = (-t247 - t287) * qJ(3) + (-t248 + t256) * pkin(2) + t231;
t253 = -t269 * g(1) - t266 * g(2);
t232 = -t271 * pkin(1) + qJDD(1) * pkin(6) + t253;
t209 = -t265 * g(3) + t268 * t232;
t245 = (-pkin(2) * t268 - qJ(3) * t265) * qJD(1);
t270 = qJD(2) ^ 2;
t290 = t268 * qJD(1);
t190 = -t270 * pkin(2) + qJDD(2) * qJ(3) + t245 * t290 + t209;
t263 = sin(pkin(8));
t292 = qJD(1) * t265;
t297 = cos(pkin(8));
t238 = t263 * qJD(2) + t297 * t292;
t168 = t297 * t185 - t263 * t190 + t238 * t303;
t221 = t263 * qJDD(2) + t297 * t247;
t208 = -t268 * g(3) - t265 * t232;
t278 = qJDD(2) * pkin(2) + t270 * qJ(3) - t245 * t292 - qJDD(3) + t208;
t237 = -t297 * qJD(2) + t263 * t292;
t288 = t237 * t290;
t302 = -(t221 + t288) * qJ(4) - t278;
t169 = t263 * t185 + t297 * t190 + t237 * t303;
t200 = Ifges(4,1) * t238 - Ifges(4,4) * t237 - Ifges(4,5) * t290;
t206 = t237 * mrSges(5,1) - t238 * mrSges(5,3);
t216 = -t237 * mrSges(5,2) - mrSges(5,3) * t290;
t219 = mrSges(5,1) * t290 + t238 * mrSges(5,2);
t220 = -t297 * qJDD(2) + t263 * t247;
t205 = t237 * pkin(3) - t238 * qJ(4);
t296 = t268 ^ 2 * t271;
t166 = t248 * pkin(3) - qJ(4) * t296 + t238 * t205 + qJDD(4) - t168;
t158 = (-t221 + t288) * pkin(7) + (t237 * t238 + t248) * pkin(4) + t166;
t299 = -2 * qJD(4);
t164 = -pkin(3) * t296 - t248 * qJ(4) - t237 * t205 + t290 * t299 + t169;
t222 = pkin(4) * t290 - t238 * pkin(7);
t235 = t237 ^ 2;
t159 = -t235 * pkin(4) + t220 * pkin(7) - t222 * t290 + t164;
t264 = sin(qJ(5));
t267 = cos(qJ(5));
t156 = t267 * t158 - t264 * t159;
t203 = t267 * t237 - t264 * t238;
t177 = t203 * qJD(5) + t264 * t220 + t267 * t221;
t204 = t264 * t237 + t267 * t238;
t183 = -t203 * mrSges(6,1) + t204 * mrSges(6,2);
t254 = qJD(5) + t290;
t191 = -t254 * mrSges(6,2) + t203 * mrSges(6,3);
t244 = qJDD(5) + t248;
t151 = m(6) * t156 + t244 * mrSges(6,1) - t177 * mrSges(6,3) - t204 * t183 + t254 * t191;
t157 = t264 * t158 + t267 * t159;
t176 = -t204 * qJD(5) + t267 * t220 - t264 * t221;
t192 = t254 * mrSges(6,1) - t204 * mrSges(6,3);
t152 = m(6) * t157 - t244 * mrSges(6,2) + t176 * mrSges(6,3) + t203 * t183 - t254 * t192;
t142 = t267 * t151 + t264 * t152;
t199 = Ifges(5,1) * t238 - Ifges(5,4) * t290 + Ifges(5,5) * t237;
t179 = Ifges(6,4) * t204 + Ifges(6,2) * t203 + Ifges(6,6) * t254;
t180 = Ifges(6,1) * t204 + Ifges(6,4) * t203 + Ifges(6,5) * t254;
t281 = mrSges(6,1) * t156 - mrSges(6,2) * t157 + Ifges(6,5) * t177 + Ifges(6,6) * t176 + Ifges(6,3) * t244 + t204 * t179 - t203 * t180;
t275 = mrSges(5,1) * t166 - mrSges(5,3) * t164 - Ifges(5,4) * t221 + Ifges(5,2) * t248 - Ifges(5,6) * t220 + pkin(4) * t142 - t237 * t199 + t281;
t280 = -m(5) * t166 - t248 * mrSges(5,1) - t142;
t143 = -t264 * t151 + t267 * t152;
t282 = m(5) * t164 - t248 * mrSges(5,3) + t143;
t195 = Ifges(5,5) * t238 - Ifges(5,6) * t290 + Ifges(5,3) * t237;
t294 = Ifges(4,4) * t238 - Ifges(4,2) * t237 - Ifges(4,6) * t290 - t195;
t301 = -t294 * t238 - mrSges(4,1) * t168 + mrSges(4,2) * t169 - Ifges(4,5) * t221 + Ifges(4,6) * t220 - pkin(3) * (-t221 * mrSges(5,2) - t238 * t206 - t216 * t290 + t280) - qJ(4) * (-t220 * mrSges(5,2) - t237 * t206 - t219 * t290 + t282) - t237 * t200 + t275;
t161 = -t235 * pkin(7) + (-pkin(3) - pkin(4)) * t220 + (pkin(3) * t290 + (2 * qJD(4)) + t222) * t238 - t302;
t153 = -m(6) * t161 + t176 * mrSges(6,1) - t177 * mrSges(6,2) + t203 * t191 - t204 * t192;
t167 = t238 * t299 + (-t238 * t290 + t220) * pkin(3) + t302;
t149 = m(5) * t167 + t220 * mrSges(5,1) - t221 * mrSges(5,3) + t237 * t216 - t238 * t219 + t153;
t178 = Ifges(6,5) * t204 + Ifges(6,6) * t203 + Ifges(6,3) * t254;
t145 = -mrSges(6,1) * t161 + mrSges(6,3) * t157 + Ifges(6,4) * t177 + Ifges(6,2) * t176 + Ifges(6,6) * t244 - t204 * t178 + t254 * t180;
t146 = mrSges(6,2) * t161 - mrSges(6,3) * t156 + Ifges(6,1) * t177 + Ifges(6,4) * t176 + Ifges(6,5) * t244 + t203 * t178 - t254 * t179;
t276 = -mrSges(5,1) * t167 + mrSges(5,2) * t164 - pkin(4) * t153 - pkin(7) * t143 - t267 * t145 - t264 * t146;
t197 = Ifges(5,4) * t238 - Ifges(5,2) * t290 + Ifges(5,6) * t237;
t295 = -Ifges(4,5) * t238 + Ifges(4,6) * t237 + Ifges(4,3) * t290 - t197;
t126 = mrSges(4,1) * t278 + mrSges(4,3) * t169 - pkin(3) * t149 + (-Ifges(4,6) + Ifges(5,6)) * t248 + t295 * t238 + (Ifges(4,4) - Ifges(5,5)) * t221 + (-Ifges(4,2) - Ifges(5,3)) * t220 + (-t199 - t200) * t290 + t276;
t277 = mrSges(5,2) * t166 - mrSges(5,3) * t167 + Ifges(5,1) * t221 - Ifges(5,4) * t248 + Ifges(5,5) * t220 - pkin(7) * t142 - t264 * t145 + t267 * t146;
t127 = -mrSges(4,2) * t278 - mrSges(4,3) * t168 + Ifges(4,1) * t221 - Ifges(4,4) * t220 - Ifges(4,5) * t248 - qJ(4) * t149 + t295 * t237 + t294 * t290 + t277;
t218 = -mrSges(4,1) * t290 - t238 * mrSges(4,3);
t293 = -t237 * mrSges(4,1) - t238 * mrSges(4,2) - t206;
t298 = -mrSges(4,3) - mrSges(5,2);
t138 = m(4) * t169 + t248 * mrSges(4,2) + t293 * t237 + t298 * t220 + (t218 - t219) * t290 + t282;
t217 = mrSges(4,2) * t290 - t237 * mrSges(4,3);
t139 = m(4) * t168 - t248 * mrSges(4,1) + t293 * t238 + t298 * t221 + (-t216 - t217) * t290 + t280;
t136 = t297 * t138 - t263 * t139;
t148 = m(4) * t278 - t220 * mrSges(4,1) - t221 * mrSges(4,2) - t237 * t217 - t238 * t218 - t149;
t229 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t265 + Ifges(3,2) * t268) * qJD(1);
t230 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t265 + Ifges(3,4) * t268) * qJD(1);
t300 = mrSges(3,1) * t208 - mrSges(3,2) * t209 + Ifges(3,5) * t247 + Ifges(3,6) * t248 + Ifges(3,3) * qJDD(2) + pkin(2) * t148 + qJ(3) * t136 + (t229 * t265 - t230 * t268) * qJD(1) + t297 * t126 + t263 * t127;
t246 = (-mrSges(3,1) * t268 + mrSges(3,2) * t265) * qJD(1);
t250 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t292;
t134 = m(3) * t209 - qJDD(2) * mrSges(3,2) + t248 * mrSges(3,3) - qJD(2) * t250 + t246 * t290 + t136;
t251 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t290;
t147 = m(3) * t208 + qJDD(2) * mrSges(3,1) - t247 * mrSges(3,3) + qJD(2) * t251 - t246 * t292 + t148;
t286 = t268 * t134 - t265 * t147;
t135 = t263 * t138 + t297 * t139;
t228 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t265 + Ifges(3,6) * t268) * qJD(1);
t123 = mrSges(3,2) * t231 - mrSges(3,3) * t208 + Ifges(3,1) * t247 + Ifges(3,4) * t248 + Ifges(3,5) * qJDD(2) - qJ(3) * t135 - qJD(2) * t229 - t263 * t126 + t297 * t127 + t228 * t290;
t125 = Ifges(3,4) * t247 - t228 * t292 + qJD(2) * t230 - mrSges(3,1) * t231 + mrSges(3,3) * t209 - pkin(2) * t135 + Ifges(3,6) * qJDD(2) + (Ifges(3,2) + Ifges(4,3)) * t248 + t301;
t274 = -m(3) * t231 + t248 * mrSges(3,1) - t247 * mrSges(3,2) - t250 * t292 + t251 * t290 - t135;
t279 = mrSges(2,1) * t252 - mrSges(2,2) * t253 + Ifges(2,3) * qJDD(1) + pkin(1) * t274 + pkin(6) * t286 + t265 * t123 + t268 * t125;
t131 = m(2) * t252 + qJDD(1) * mrSges(2,1) - t271 * mrSges(2,2) + t274;
t130 = t265 * t134 + t268 * t147;
t128 = m(2) * t253 - t271 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t286;
t121 = mrSges(2,1) * g(3) + mrSges(2,3) * t253 + t271 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t130 - t300;
t120 = -mrSges(2,2) * g(3) - mrSges(2,3) * t252 + Ifges(2,5) * qJDD(1) - t271 * Ifges(2,6) - pkin(6) * t130 + t268 * t123 - t265 * t125;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t269 * t120 - t266 * t121 - pkin(5) * (t266 * t128 + t269 * t131), t120, t123, t127, -t195 * t290 - t237 * t197 + t277, t146; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t266 * t120 + t269 * t121 + pkin(5) * (t269 * t128 - t266 * t131), t121, t125, t126, -t238 * t195 - t275, t145; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t279, t279, t300, -Ifges(4,3) * t248 - t301, Ifges(5,5) * t221 - Ifges(5,6) * t248 + Ifges(5,3) * t220 + t238 * t197 + t199 * t290 - t276, t281;];
m_new = t1;
