% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRR11
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
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRR11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR11_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR11_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR11_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR11_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:39:09
% EndTime: 2019-12-31 22:39:46
% DurationCPUTime: 20.17s
% Computational Cost: add. (353356->324), mult. (755505->426), div. (0->0), fcn. (592079->12), ass. (0->138)
t269 = sin(pkin(5));
t274 = sin(qJ(2));
t279 = cos(qJ(2));
t293 = qJD(1) * qJD(2);
t256 = (-qJDD(1) * t279 + t274 * t293) * t269;
t302 = pkin(7) * t269;
t270 = cos(pkin(5));
t301 = t270 * g(3);
t300 = t269 * t274;
t299 = t269 * t279;
t298 = t270 * t274;
t297 = t270 * t279;
t295 = qJD(1) * t269;
t254 = (-pkin(2) * t279 - pkin(8) * t274) * t295;
t265 = t270 * qJD(1) + qJD(2);
t263 = t265 ^ 2;
t264 = t270 * qJDD(1) + qJDD(2);
t294 = qJD(1) * t279;
t275 = sin(qJ(1));
t280 = cos(qJ(1));
t261 = t275 * g(1) - t280 * g(2);
t281 = qJD(1) ^ 2;
t251 = qJDD(1) * pkin(1) + t281 * t302 + t261;
t262 = -t280 * g(1) - t275 * g(2);
t252 = -t281 * pkin(1) + qJDD(1) * t302 + t262;
t296 = t251 * t298 + t279 * t252;
t206 = -t263 * pkin(2) + t264 * pkin(8) + (-g(3) * t274 + t254 * t294) * t269 + t296;
t255 = (qJDD(1) * t274 + t279 * t293) * t269;
t207 = t256 * pkin(2) - t255 * pkin(8) - t301 + (-t251 + (pkin(2) * t274 - pkin(8) * t279) * t265 * qJD(1)) * t269;
t273 = sin(qJ(3));
t278 = cos(qJ(3));
t187 = t278 * t206 + t273 * t207;
t292 = t274 * t295;
t243 = t278 * t265 - t273 * t292;
t244 = t273 * t265 + t278 * t292;
t228 = -t243 * pkin(3) - t244 * pkin(9);
t248 = qJDD(3) + t256;
t291 = t269 * t294;
t260 = qJD(3) - t291;
t258 = t260 ^ 2;
t181 = -t258 * pkin(3) + t248 * pkin(9) + t243 * t228 + t187;
t225 = -g(3) * t299 + t251 * t297 - t274 * t252;
t205 = -t264 * pkin(2) - t263 * pkin(8) + t254 * t292 - t225;
t223 = -t244 * qJD(3) - t273 * t255 + t278 * t264;
t224 = t243 * qJD(3) + t278 * t255 + t273 * t264;
t185 = (-t243 * t260 - t224) * pkin(9) + (t244 * t260 - t223) * pkin(3) + t205;
t272 = sin(qJ(4));
t277 = cos(qJ(4));
t171 = -t272 * t181 + t277 * t185;
t230 = -t272 * t244 + t277 * t260;
t195 = t230 * qJD(4) + t277 * t224 + t272 * t248;
t221 = qJDD(4) - t223;
t231 = t277 * t244 + t272 * t260;
t242 = qJD(4) - t243;
t169 = (t230 * t242 - t195) * pkin(10) + (t230 * t231 + t221) * pkin(4) + t171;
t172 = t277 * t181 + t272 * t185;
t194 = -t231 * qJD(4) - t272 * t224 + t277 * t248;
t214 = t242 * pkin(4) - t231 * pkin(10);
t229 = t230 ^ 2;
t170 = -t229 * pkin(4) + t194 * pkin(10) - t242 * t214 + t172;
t271 = sin(qJ(5));
t276 = cos(qJ(5));
t167 = t276 * t169 - t271 * t170;
t208 = t276 * t230 - t271 * t231;
t178 = t208 * qJD(5) + t271 * t194 + t276 * t195;
t209 = t271 * t230 + t276 * t231;
t192 = -t208 * mrSges(6,1) + t209 * mrSges(6,2);
t240 = qJD(5) + t242;
t196 = -t240 * mrSges(6,2) + t208 * mrSges(6,3);
t216 = qJDD(5) + t221;
t163 = m(6) * t167 + t216 * mrSges(6,1) - t178 * mrSges(6,3) - t209 * t192 + t240 * t196;
t168 = t271 * t169 + t276 * t170;
t177 = -t209 * qJD(5) + t276 * t194 - t271 * t195;
t197 = t240 * mrSges(6,1) - t209 * mrSges(6,3);
t164 = m(6) * t168 - t216 * mrSges(6,2) + t177 * mrSges(6,3) + t208 * t192 - t240 * t197;
t155 = t276 * t163 + t271 * t164;
t210 = -t230 * mrSges(5,1) + t231 * mrSges(5,2);
t212 = -t242 * mrSges(5,2) + t230 * mrSges(5,3);
t153 = m(5) * t171 + t221 * mrSges(5,1) - t195 * mrSges(5,3) - t231 * t210 + t242 * t212 + t155;
t213 = t242 * mrSges(5,1) - t231 * mrSges(5,3);
t289 = -t271 * t163 + t276 * t164;
t154 = m(5) * t172 - t221 * mrSges(5,2) + t194 * mrSges(5,3) + t230 * t210 - t242 * t213 + t289;
t151 = -t272 * t153 + t277 * t154;
t227 = -t243 * mrSges(4,1) + t244 * mrSges(4,2);
t233 = t260 * mrSges(4,1) - t244 * mrSges(4,3);
t149 = m(4) * t187 - t248 * mrSges(4,2) + t223 * mrSges(4,3) + t243 * t227 - t260 * t233 + t151;
t186 = -t273 * t206 + t278 * t207;
t180 = -t248 * pkin(3) - t258 * pkin(9) + t244 * t228 - t186;
t173 = -t194 * pkin(4) - t229 * pkin(10) + t231 * t214 + t180;
t287 = m(6) * t173 - t177 * mrSges(6,1) + t178 * mrSges(6,2) - t208 * t196 + t209 * t197;
t165 = -m(5) * t180 + t194 * mrSges(5,1) - t195 * mrSges(5,2) + t230 * t212 - t231 * t213 - t287;
t232 = -t260 * mrSges(4,2) + t243 * mrSges(4,3);
t159 = m(4) * t186 + t248 * mrSges(4,1) - t224 * mrSges(4,3) - t244 * t227 + t260 * t232 + t165;
t142 = t273 * t149 + t278 * t159;
t226 = -g(3) * t300 + t296;
t249 = t265 * mrSges(3,1) - mrSges(3,3) * t292;
t253 = (-mrSges(3,1) * t279 + mrSges(3,2) * t274) * t295;
t290 = t278 * t149 - t273 * t159;
t140 = m(3) * t226 - t264 * mrSges(3,2) - t256 * mrSges(3,3) - t265 * t249 + t253 * t291 + t290;
t250 = -t265 * mrSges(3,2) + mrSges(3,3) * t291;
t150 = t277 * t153 + t272 * t154;
t284 = -m(4) * t205 + t223 * mrSges(4,1) - t224 * mrSges(4,2) + t243 * t232 - t244 * t233 - t150;
t146 = m(3) * t225 + t264 * mrSges(3,1) - t255 * mrSges(3,3) + t265 * t250 - t253 * t292 + t284;
t136 = t279 * t140 - t274 * t146;
t237 = -t269 * t251 - t301;
t141 = m(3) * t237 + t256 * mrSges(3,1) + t255 * mrSges(3,2) + (t249 * t274 - t250 * t279) * t295 + t142;
t132 = t140 * t298 - t269 * t141 + t146 * t297;
t188 = Ifges(6,5) * t209 + Ifges(6,6) * t208 + Ifges(6,3) * t240;
t190 = Ifges(6,1) * t209 + Ifges(6,4) * t208 + Ifges(6,5) * t240;
t156 = -mrSges(6,1) * t173 + mrSges(6,3) * t168 + Ifges(6,4) * t178 + Ifges(6,2) * t177 + Ifges(6,6) * t216 - t209 * t188 + t240 * t190;
t189 = Ifges(6,4) * t209 + Ifges(6,2) * t208 + Ifges(6,6) * t240;
t157 = mrSges(6,2) * t173 - mrSges(6,3) * t167 + Ifges(6,1) * t178 + Ifges(6,4) * t177 + Ifges(6,5) * t216 + t208 * t188 - t240 * t189;
t198 = Ifges(5,5) * t231 + Ifges(5,6) * t230 + Ifges(5,3) * t242;
t200 = Ifges(5,1) * t231 + Ifges(5,4) * t230 + Ifges(5,5) * t242;
t143 = -mrSges(5,1) * t180 + mrSges(5,3) * t172 + Ifges(5,4) * t195 + Ifges(5,2) * t194 + Ifges(5,6) * t221 - pkin(4) * t287 + pkin(10) * t289 + t276 * t156 + t271 * t157 - t231 * t198 + t242 * t200;
t199 = Ifges(5,4) * t231 + Ifges(5,2) * t230 + Ifges(5,6) * t242;
t144 = mrSges(5,2) * t180 - mrSges(5,3) * t171 + Ifges(5,1) * t195 + Ifges(5,4) * t194 + Ifges(5,5) * t221 - pkin(10) * t155 - t271 * t156 + t276 * t157 + t230 * t198 - t242 * t199;
t217 = Ifges(4,5) * t244 + Ifges(4,6) * t243 + Ifges(4,3) * t260;
t218 = Ifges(4,4) * t244 + Ifges(4,2) * t243 + Ifges(4,6) * t260;
t133 = mrSges(4,2) * t205 - mrSges(4,3) * t186 + Ifges(4,1) * t224 + Ifges(4,4) * t223 + Ifges(4,5) * t248 - pkin(9) * t150 - t272 * t143 + t277 * t144 + t243 * t217 - t260 * t218;
t219 = Ifges(4,1) * t244 + Ifges(4,4) * t243 + Ifges(4,5) * t260;
t285 = -mrSges(6,1) * t167 + mrSges(6,2) * t168 - Ifges(6,5) * t178 - Ifges(6,6) * t177 - Ifges(6,3) * t216 - t209 * t189 + t208 * t190;
t282 = mrSges(5,1) * t171 - mrSges(5,2) * t172 + Ifges(5,5) * t195 + Ifges(5,6) * t194 + Ifges(5,3) * t221 + pkin(4) * t155 + t231 * t199 - t230 * t200 - t285;
t137 = -mrSges(4,1) * t205 + mrSges(4,3) * t187 + Ifges(4,4) * t224 + Ifges(4,2) * t223 + Ifges(4,6) * t248 - pkin(3) * t150 - t244 * t217 + t260 * t219 - t282;
t235 = Ifges(3,6) * t265 + (Ifges(3,4) * t274 + Ifges(3,2) * t279) * t295;
t236 = Ifges(3,5) * t265 + (Ifges(3,1) * t274 + Ifges(3,4) * t279) * t295;
t124 = Ifges(3,5) * t255 - Ifges(3,6) * t256 + Ifges(3,3) * t264 + mrSges(3,1) * t225 - mrSges(3,2) * t226 + t273 * t133 + t278 * t137 + pkin(2) * t284 + pkin(8) * t290 + (t235 * t274 - t236 * t279) * t295;
t234 = Ifges(3,3) * t265 + (Ifges(3,5) * t274 + Ifges(3,6) * t279) * t295;
t126 = mrSges(3,2) * t237 - mrSges(3,3) * t225 + Ifges(3,1) * t255 - Ifges(3,4) * t256 + Ifges(3,5) * t264 - pkin(8) * t142 + t278 * t133 - t273 * t137 + t234 * t291 - t265 * t235;
t283 = mrSges(4,1) * t186 - mrSges(4,2) * t187 + Ifges(4,5) * t224 + Ifges(4,6) * t223 + Ifges(4,3) * t248 + pkin(3) * t165 + pkin(9) * t151 + t277 * t143 + t272 * t144 + t244 * t218 - t243 * t219;
t128 = -mrSges(3,1) * t237 + mrSges(3,3) * t226 + Ifges(3,4) * t255 - Ifges(3,2) * t256 + Ifges(3,6) * t264 - pkin(2) * t142 - t234 * t292 + t265 * t236 - t283;
t286 = mrSges(2,1) * t261 - mrSges(2,2) * t262 + Ifges(2,3) * qJDD(1) + pkin(1) * t132 + t270 * t124 + t126 * t300 + t128 * t299 + t136 * t302;
t134 = m(2) * t262 - t281 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t136;
t131 = t270 * t141 + (t140 * t274 + t146 * t279) * t269;
t129 = m(2) * t261 + qJDD(1) * mrSges(2,1) - t281 * mrSges(2,2) + t132;
t122 = -mrSges(2,2) * g(3) - mrSges(2,3) * t261 + Ifges(2,5) * qJDD(1) - t281 * Ifges(2,6) + t279 * t126 - t274 * t128 + (-t131 * t269 - t132 * t270) * pkin(7);
t121 = mrSges(2,1) * g(3) + mrSges(2,3) * t262 + t281 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t131 - t269 * t124 + (pkin(7) * t136 + t126 * t274 + t128 * t279) * t270;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t280 * t122 - t275 * t121 - pkin(6) * (t280 * t129 + t275 * t134), t122, t126, t133, t144, t157; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t275 * t122 + t280 * t121 + pkin(6) * (-t275 * t129 + t280 * t134), t121, t128, t137, t143, t156; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t286, t286, t124, t283, t282, -t285;];
m_new = t1;
