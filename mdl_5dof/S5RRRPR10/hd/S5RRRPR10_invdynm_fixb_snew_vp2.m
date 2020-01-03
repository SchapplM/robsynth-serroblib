% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPR10
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR10_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR10_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR10_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR10_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:26:57
% EndTime: 2019-12-31 21:27:32
% DurationCPUTime: 18.51s
% Computational Cost: add. (320557->324), mult. (709882->428), div. (0->0), fcn. (555197->12), ass. (0->138)
t270 = sin(pkin(5));
t275 = sin(qJ(2));
t279 = cos(qJ(2));
t296 = qJD(1) * qJD(2);
t257 = (-qJDD(1) * t279 + t275 * t296) * t270;
t306 = 2 * qJD(4);
t305 = pkin(7) * t270;
t272 = cos(pkin(5));
t304 = t272 * g(3);
t303 = t270 * t275;
t302 = t270 * t279;
t301 = t272 * t275;
t300 = t272 * t279;
t298 = qJD(1) * t270;
t255 = (-pkin(2) * t279 - pkin(8) * t275) * t298;
t265 = qJD(1) * t272 + qJD(2);
t263 = t265 ^ 2;
t264 = qJDD(1) * t272 + qJDD(2);
t297 = qJD(1) * t279;
t276 = sin(qJ(1));
t280 = cos(qJ(1));
t261 = t276 * g(1) - g(2) * t280;
t281 = qJD(1) ^ 2;
t252 = qJDD(1) * pkin(1) + t281 * t305 + t261;
t262 = -g(1) * t280 - g(2) * t276;
t253 = -pkin(1) * t281 + qJDD(1) * t305 + t262;
t299 = t252 * t301 + t279 * t253;
t213 = -t263 * pkin(2) + t264 * pkin(8) + (-g(3) * t275 + t255 * t297) * t270 + t299;
t256 = (qJDD(1) * t275 + t279 * t296) * t270;
t214 = t257 * pkin(2) - t256 * pkin(8) - t304 + (-t252 + (pkin(2) * t275 - pkin(8) * t279) * t265 * qJD(1)) * t270;
t274 = sin(qJ(3));
t278 = cos(qJ(3));
t187 = -t274 * t213 + t278 * t214;
t295 = t275 * t298;
t244 = t265 * t278 - t274 * t295;
t226 = qJD(3) * t244 + t256 * t278 + t264 * t274;
t245 = t265 * t274 + t278 * t295;
t249 = qJDD(3) + t257;
t294 = t270 * t297;
t260 = qJD(3) - t294;
t180 = (t244 * t260 - t226) * qJ(4) + (t244 * t245 + t249) * pkin(3) + t187;
t188 = t278 * t213 + t274 * t214;
t225 = -qJD(3) * t245 - t256 * t274 + t264 * t278;
t235 = pkin(3) * t260 - qJ(4) * t245;
t243 = t244 ^ 2;
t182 = -pkin(3) * t243 + qJ(4) * t225 - t235 * t260 + t188;
t269 = sin(pkin(10));
t271 = cos(pkin(10));
t231 = t244 * t271 - t245 * t269;
t177 = t180 * t269 + t271 * t182 + t231 * t306;
t201 = t225 * t271 - t226 * t269;
t232 = t244 * t269 + t245 * t271;
t207 = -mrSges(5,1) * t231 + mrSges(5,2) * t232;
t218 = mrSges(5,1) * t260 - mrSges(5,3) * t232;
t208 = -pkin(4) * t231 - pkin(9) * t232;
t259 = t260 ^ 2;
t174 = -pkin(4) * t259 + pkin(9) * t249 + t208 * t231 + t177;
t227 = -g(3) * t302 + t252 * t300 - t275 * t253;
t212 = -t264 * pkin(2) - t263 * pkin(8) + t255 * t295 - t227;
t183 = -t225 * pkin(3) - t243 * qJ(4) + t245 * t235 + qJDD(4) + t212;
t202 = t225 * t269 + t226 * t271;
t178 = (-t231 * t260 - t202) * pkin(9) + (t232 * t260 - t201) * pkin(4) + t183;
t273 = sin(qJ(5));
t277 = cos(qJ(5));
t171 = -t174 * t273 + t178 * t277;
t215 = -t232 * t273 + t260 * t277;
t186 = qJD(5) * t215 + t202 * t277 + t249 * t273;
t216 = t232 * t277 + t260 * t273;
t194 = -mrSges(6,1) * t215 + mrSges(6,2) * t216;
t230 = qJD(5) - t231;
t195 = -mrSges(6,2) * t230 + mrSges(6,3) * t215;
t200 = qJDD(5) - t201;
t167 = m(6) * t171 + mrSges(6,1) * t200 - mrSges(6,3) * t186 - t194 * t216 + t195 * t230;
t172 = t174 * t277 + t178 * t273;
t185 = -qJD(5) * t216 - t202 * t273 + t249 * t277;
t196 = mrSges(6,1) * t230 - mrSges(6,3) * t216;
t168 = m(6) * t172 - mrSges(6,2) * t200 + mrSges(6,3) * t185 + t194 * t215 - t196 * t230;
t291 = -t167 * t273 + t168 * t277;
t154 = m(5) * t177 - mrSges(5,2) * t249 + mrSges(5,3) * t201 + t207 * t231 - t218 * t260 + t291;
t290 = -t271 * t180 + t269 * t182;
t176 = -0.2e1 * qJD(4) * t232 - t290;
t217 = -mrSges(5,2) * t260 + mrSges(5,3) * t231;
t173 = -t249 * pkin(4) - t259 * pkin(9) + (t306 + t208) * t232 + t290;
t288 = -m(6) * t173 + t185 * mrSges(6,1) - mrSges(6,2) * t186 + t215 * t195 - t196 * t216;
t163 = m(5) * t176 + mrSges(5,1) * t249 - mrSges(5,3) * t202 - t207 * t232 + t217 * t260 + t288;
t149 = t154 * t269 + t163 * t271;
t233 = -mrSges(4,1) * t244 + mrSges(4,2) * t245;
t234 = -mrSges(4,2) * t260 + mrSges(4,3) * t244;
t147 = m(4) * t187 + mrSges(4,1) * t249 - mrSges(4,3) * t226 - t233 * t245 + t234 * t260 + t149;
t236 = mrSges(4,1) * t260 - mrSges(4,3) * t245;
t292 = t154 * t271 - t163 * t269;
t148 = m(4) * t188 - mrSges(4,2) * t249 + mrSges(4,3) * t225 + t233 * t244 - t236 * t260 + t292;
t141 = t147 * t278 + t148 * t274;
t156 = t167 * t277 + t168 * t273;
t228 = -g(3) * t303 + t299;
t250 = mrSges(3,1) * t265 - mrSges(3,3) * t295;
t254 = (-mrSges(3,1) * t279 + mrSges(3,2) * t275) * t298;
t293 = -t147 * t274 + t148 * t278;
t139 = m(3) * t228 - mrSges(3,2) * t264 - mrSges(3,3) * t257 - t250 * t265 + t254 * t294 + t293;
t251 = -mrSges(3,2) * t265 + mrSges(3,3) * t294;
t286 = m(5) * t183 - t201 * mrSges(5,1) + mrSges(5,2) * t202 - t231 * t217 + t218 * t232 + t156;
t283 = -m(4) * t212 + t225 * mrSges(4,1) - mrSges(4,2) * t226 + t244 * t234 - t236 * t245 - t286;
t151 = m(3) * t227 + mrSges(3,1) * t264 - mrSges(3,3) * t256 + t251 * t265 - t254 * t295 + t283;
t136 = t139 * t279 - t151 * t275;
t240 = -t270 * t252 - t304;
t140 = m(3) * t240 + t257 * mrSges(3,1) + t256 * mrSges(3,2) + (t250 * t275 - t251 * t279) * t298 + t141;
t131 = t139 * t301 - t140 * t270 + t151 * t300;
t189 = Ifges(6,5) * t216 + Ifges(6,6) * t215 + Ifges(6,3) * t230;
t191 = Ifges(6,1) * t216 + Ifges(6,4) * t215 + Ifges(6,5) * t230;
t160 = -mrSges(6,1) * t173 + mrSges(6,3) * t172 + Ifges(6,4) * t186 + Ifges(6,2) * t185 + Ifges(6,6) * t200 - t189 * t216 + t191 * t230;
t190 = Ifges(6,4) * t216 + Ifges(6,2) * t215 + Ifges(6,6) * t230;
t161 = mrSges(6,2) * t173 - mrSges(6,3) * t171 + Ifges(6,1) * t186 + Ifges(6,4) * t185 + Ifges(6,5) * t200 + t189 * t215 - t190 * t230;
t203 = Ifges(5,5) * t232 + Ifges(5,6) * t231 + Ifges(5,3) * t260;
t204 = Ifges(5,4) * t232 + Ifges(5,2) * t231 + Ifges(5,6) * t260;
t142 = mrSges(5,2) * t183 - mrSges(5,3) * t176 + Ifges(5,1) * t202 + Ifges(5,4) * t201 + Ifges(5,5) * t249 - pkin(9) * t156 - t160 * t273 + t161 * t277 + t203 * t231 - t204 * t260;
t205 = Ifges(5,1) * t232 + Ifges(5,4) * t231 + Ifges(5,5) * t260;
t284 = mrSges(6,1) * t171 - mrSges(6,2) * t172 + Ifges(6,5) * t186 + Ifges(6,6) * t185 + Ifges(6,3) * t200 + t190 * t216 - t191 * t215;
t143 = -mrSges(5,1) * t183 + mrSges(5,3) * t177 + Ifges(5,4) * t202 + Ifges(5,2) * t201 + Ifges(5,6) * t249 - pkin(4) * t156 - t203 * t232 + t205 * t260 - t284;
t219 = Ifges(4,5) * t245 + Ifges(4,6) * t244 + Ifges(4,3) * t260;
t221 = Ifges(4,1) * t245 + Ifges(4,4) * t244 + Ifges(4,5) * t260;
t132 = -mrSges(4,1) * t212 + mrSges(4,3) * t188 + Ifges(4,4) * t226 + Ifges(4,2) * t225 + Ifges(4,6) * t249 - pkin(3) * t286 + qJ(4) * t292 + t269 * t142 + t271 * t143 - t245 * t219 + t260 * t221;
t220 = Ifges(4,4) * t245 + Ifges(4,2) * t244 + Ifges(4,6) * t260;
t133 = mrSges(4,2) * t212 - mrSges(4,3) * t187 + Ifges(4,1) * t226 + Ifges(4,4) * t225 + Ifges(4,5) * t249 - qJ(4) * t149 + t142 * t271 - t143 * t269 + t219 * t244 - t220 * t260;
t238 = Ifges(3,6) * t265 + (Ifges(3,4) * t275 + Ifges(3,2) * t279) * t298;
t239 = Ifges(3,5) * t265 + (Ifges(3,1) * t275 + Ifges(3,4) * t279) * t298;
t123 = Ifges(3,5) * t256 - Ifges(3,6) * t257 + Ifges(3,3) * t264 + mrSges(3,1) * t227 - mrSges(3,2) * t228 + t274 * t133 + t278 * t132 + pkin(2) * t283 + pkin(8) * t293 + (t238 * t275 - t239 * t279) * t298;
t237 = Ifges(3,3) * t265 + (Ifges(3,5) * t275 + Ifges(3,6) * t279) * t298;
t125 = mrSges(3,2) * t240 - mrSges(3,3) * t227 + Ifges(3,1) * t256 - Ifges(3,4) * t257 + Ifges(3,5) * t264 - pkin(8) * t141 - t132 * t274 + t133 * t278 + t237 * t294 - t238 * t265;
t285 = -mrSges(5,1) * t176 + mrSges(5,2) * t177 - Ifges(5,5) * t202 - Ifges(5,6) * t201 - Ifges(5,3) * t249 - pkin(4) * t288 - pkin(9) * t291 - t160 * t277 - t161 * t273 - t232 * t204 + t205 * t231;
t282 = mrSges(4,1) * t187 - mrSges(4,2) * t188 + Ifges(4,5) * t226 + Ifges(4,6) * t225 + Ifges(4,3) * t249 + pkin(3) * t149 + t220 * t245 - t221 * t244 - t285;
t127 = -mrSges(3,1) * t240 + mrSges(3,3) * t228 + Ifges(3,4) * t256 - Ifges(3,2) * t257 + Ifges(3,6) * t264 - pkin(2) * t141 - t237 * t295 + t239 * t265 - t282;
t287 = mrSges(2,1) * t261 - mrSges(2,2) * t262 + Ifges(2,3) * qJDD(1) + pkin(1) * t131 + t123 * t272 + t125 * t303 + t127 * t302 + t136 * t305;
t134 = m(2) * t262 - mrSges(2,1) * t281 - qJDD(1) * mrSges(2,2) + t136;
t130 = t272 * t140 + (t139 * t275 + t151 * t279) * t270;
t128 = m(2) * t261 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t281 + t131;
t121 = -mrSges(2,2) * g(3) - mrSges(2,3) * t261 + Ifges(2,5) * qJDD(1) - t281 * Ifges(2,6) + t279 * t125 - t275 * t127 + (-t130 * t270 - t131 * t272) * pkin(7);
t120 = mrSges(2,1) * g(3) + mrSges(2,3) * t262 + t281 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t130 - t270 * t123 + (pkin(7) * t136 + t125 * t275 + t127 * t279) * t272;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t280 * t121 - t276 * t120 - pkin(6) * (t128 * t280 + t134 * t276), t121, t125, t133, t142, t161; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t276 * t121 + t280 * t120 + pkin(6) * (-t128 * t276 + t134 * t280), t120, t127, t132, t143, t160; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t287, t287, t123, t282, -t285, t284;];
m_new = t1;
