% Calculate vector of cutting torques with Newton-Euler for
% S6PPRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-05-04 20:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PPRRRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:27:04
% EndTime: 2019-05-04 20:27:23
% DurationCPUTime: 14.93s
% Computational Cost: add. (293863->286), mult. (527935->364), div. (0->0), fcn. (402678->14), ass. (0->126)
t260 = sin(pkin(11));
t264 = cos(pkin(11));
t254 = -g(1) * t264 - g(2) * t260;
t259 = sin(pkin(12));
t263 = cos(pkin(12));
t253 = g(1) * t260 - g(2) * t264;
t258 = -g(3) + qJDD(1);
t262 = sin(pkin(6));
t266 = cos(pkin(6));
t283 = t253 * t266 + t258 * t262;
t203 = -t259 * t254 + t263 * t283;
t204 = t263 * t254 + t259 * t283;
t233 = -t253 * t262 + t258 * t266 + qJDD(2);
t272 = cos(qJ(3));
t265 = cos(pkin(7));
t269 = sin(qJ(3));
t297 = t265 * t269;
t261 = sin(pkin(7));
t298 = t261 * t269;
t193 = t203 * t297 + t272 * t204 + t233 * t298;
t274 = qJD(3) ^ 2;
t191 = -pkin(3) * t274 + qJDD(3) * pkin(9) + t193;
t195 = -t203 * t261 + t233 * t265;
t268 = sin(qJ(4));
t271 = cos(qJ(4));
t184 = t271 * t191 + t268 * t195;
t249 = (-pkin(4) * t271 - pkin(10) * t268) * qJD(3);
t273 = qJD(4) ^ 2;
t293 = qJD(3) * t271;
t182 = -pkin(4) * t273 + qJDD(4) * pkin(10) + t249 * t293 + t184;
t192 = -t269 * t204 + (t203 * t265 + t233 * t261) * t272;
t190 = -qJDD(3) * pkin(3) - t274 * pkin(9) - t192;
t292 = qJD(3) * qJD(4);
t288 = t271 * t292;
t250 = qJDD(3) * t268 + t288;
t289 = t268 * t292;
t251 = qJDD(3) * t271 - t289;
t187 = (-t250 - t288) * pkin(10) + (-t251 + t289) * pkin(4) + t190;
t267 = sin(qJ(5));
t270 = cos(qJ(5));
t176 = -t267 * t182 + t270 * t187;
t294 = qJD(3) * t268;
t246 = qJD(4) * t270 - t267 * t294;
t223 = qJD(5) * t246 + qJDD(4) * t267 + t250 * t270;
t247 = qJD(4) * t267 + t270 * t294;
t225 = -mrSges(7,1) * t246 + mrSges(7,2) * t247;
t226 = -mrSges(6,1) * t246 + mrSges(6,2) * t247;
t257 = qJD(5) - t293;
t229 = -mrSges(6,2) * t257 + mrSges(6,3) * t246;
t245 = qJDD(5) - t251;
t172 = -0.2e1 * qJD(6) * t247 + (t246 * t257 - t223) * qJ(6) + (t246 * t247 + t245) * pkin(5) + t176;
t228 = -mrSges(7,2) * t257 + mrSges(7,3) * t246;
t291 = m(7) * t172 + t245 * mrSges(7,1) + t257 * t228;
t163 = m(6) * t176 + t245 * mrSges(6,1) + t257 * t229 + (-t225 - t226) * t247 + (-mrSges(6,3) - mrSges(7,3)) * t223 + t291;
t177 = t270 * t182 + t267 * t187;
t222 = -qJD(5) * t247 + qJDD(4) * t270 - t250 * t267;
t230 = pkin(5) * t257 - qJ(6) * t247;
t244 = t246 ^ 2;
t175 = -pkin(5) * t244 + qJ(6) * t222 + 0.2e1 * qJD(6) * t246 - t230 * t257 + t177;
t290 = m(7) * t175 + t222 * mrSges(7,3) + t246 * t225;
t231 = mrSges(7,1) * t257 - mrSges(7,3) * t247;
t295 = -mrSges(6,1) * t257 + mrSges(6,3) * t247 - t231;
t305 = -mrSges(6,2) - mrSges(7,2);
t165 = m(6) * t177 + t222 * mrSges(6,3) + t246 * t226 + t245 * t305 + t295 * t257 + t290;
t162 = -t163 * t267 + t270 * t165;
t248 = (-mrSges(5,1) * t271 + mrSges(5,2) * t268) * qJD(3);
t255 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t294;
t159 = m(5) * t184 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t251 - qJD(4) * t255 + t248 * t293 + t162;
t183 = -t268 * t191 + t195 * t271;
t181 = -qJDD(4) * pkin(4) - pkin(10) * t273 + t249 * t294 - t183;
t179 = -pkin(5) * t222 - qJ(6) * t244 + t230 * t247 + qJDD(6) + t181;
t286 = -m(7) * t179 + t222 * mrSges(7,1) + t246 * t228;
t168 = -m(6) * t181 + t222 * mrSges(6,1) + t223 * t305 + t246 * t229 + t295 * t247 + t286;
t256 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t293;
t167 = m(5) * t183 + qJDD(4) * mrSges(5,1) - t250 * mrSges(5,3) + qJD(4) * t256 - t248 * t294 + t168;
t287 = t271 * t159 - t167 * t268;
t149 = m(4) * t193 - mrSges(4,1) * t274 - qJDD(3) * mrSges(4,2) + t287;
t152 = t268 * t159 + t271 * t167;
t151 = m(4) * t195 + t152;
t161 = t163 * t270 + t165 * t267;
t277 = -m(5) * t190 + t251 * mrSges(5,1) - mrSges(5,2) * t250 - t255 * t294 + t256 * t293 - t161;
t156 = m(4) * t192 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t274 + t277;
t299 = t156 * t272;
t139 = t149 * t297 - t151 * t261 + t265 * t299;
t136 = m(3) * t203 + t139;
t144 = t272 * t149 - t156 * t269;
t143 = m(3) * t204 + t144;
t309 = t136 * t263 + t143 * t259;
t205 = Ifges(7,5) * t247 + Ifges(7,6) * t246 + Ifges(7,3) * t257;
t206 = Ifges(6,5) * t247 + Ifges(6,6) * t246 + Ifges(6,3) * t257;
t210 = Ifges(6,1) * t247 + Ifges(6,4) * t246 + Ifges(6,5) * t257;
t209 = Ifges(7,1) * t247 + Ifges(7,4) * t246 + Ifges(7,5) * t257;
t281 = -mrSges(7,1) * t179 + mrSges(7,3) * t175 + Ifges(7,4) * t223 + Ifges(7,2) * t222 + Ifges(7,6) * t245 + t257 * t209;
t153 = Ifges(6,4) * t223 + Ifges(6,2) * t222 + Ifges(6,6) * t245 + t257 * t210 - mrSges(6,1) * t181 + mrSges(6,3) * t177 - pkin(5) * (t223 * mrSges(7,2) - t286) + qJ(6) * (-t245 * mrSges(7,2) - t257 * t231 + t290) + (-pkin(5) * t231 - t205 - t206) * t247 + t281;
t169 = -t223 * mrSges(7,3) - t247 * t225 + t291;
t207 = Ifges(7,4) * t247 + Ifges(7,2) * t246 + Ifges(7,6) * t257;
t208 = Ifges(6,4) * t247 + Ifges(6,2) * t246 + Ifges(6,6) * t257;
t279 = mrSges(7,2) * t179 - mrSges(7,3) * t172 + Ifges(7,1) * t223 + Ifges(7,4) * t222 + Ifges(7,5) * t245 + t246 * t205;
t160 = mrSges(6,2) * t181 - mrSges(6,3) * t176 + Ifges(6,1) * t223 + Ifges(6,4) * t222 + Ifges(6,5) * t245 - qJ(6) * t169 + t246 * t206 + (-t207 - t208) * t257 + t279;
t235 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t268 + Ifges(5,6) * t271) * qJD(3);
t236 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t268 + Ifges(5,2) * t271) * qJD(3);
t140 = mrSges(5,2) * t190 - mrSges(5,3) * t183 + Ifges(5,1) * t250 + Ifges(5,4) * t251 + Ifges(5,5) * qJDD(4) - pkin(10) * t161 - qJD(4) * t236 - t153 * t267 + t160 * t270 + t235 * t293;
t237 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t268 + Ifges(5,4) * t271) * qJD(3);
t280 = -mrSges(7,1) * t172 + mrSges(7,2) * t175 - Ifges(7,5) * t223 - Ifges(7,6) * t222 - Ifges(7,3) * t245 - t247 * t207;
t307 = mrSges(6,1) * t176 - mrSges(6,2) * t177 + Ifges(6,5) * t223 + Ifges(6,6) * t222 + Ifges(6,3) * t245 + pkin(5) * t169 + t247 * t208 - (t210 + t209) * t246 - t280;
t145 = -mrSges(5,1) * t190 + mrSges(5,3) * t184 + Ifges(5,4) * t250 + Ifges(5,2) * t251 + Ifges(5,6) * qJDD(4) - pkin(4) * t161 + qJD(4) * t237 - t235 * t294 - t307;
t129 = mrSges(4,1) * t192 - mrSges(4,2) * t193 + Ifges(4,3) * qJDD(3) + pkin(3) * t277 + pkin(9) * t287 + t268 * t140 + t271 * t145;
t138 = t149 * t298 + t265 * t151 + t261 * t299;
t130 = mrSges(4,2) * t195 - mrSges(4,3) * t192 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t274 - pkin(9) * t152 + t140 * t271 - t145 * t268;
t306 = mrSges(5,1) * t183 - mrSges(5,2) * t184 + Ifges(5,5) * t250 + Ifges(5,6) * t251 + Ifges(5,3) * qJDD(4) + pkin(4) * t168 + pkin(10) * t162 + t270 * t153 + t267 * t160 + (t236 * t268 - t237 * t271) * qJD(3);
t134 = -mrSges(4,1) * t195 + mrSges(4,3) * t193 + t274 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t152 - t306;
t282 = pkin(8) * t144 + t130 * t269 + t134 * t272;
t122 = -mrSges(3,1) * t233 + mrSges(3,3) * t204 - pkin(2) * t138 - t261 * t129 + t265 * t282;
t124 = mrSges(3,2) * t233 - mrSges(3,3) * t203 + t272 * t130 - t269 * t134 + (-t138 * t261 - t139 * t265) * pkin(8);
t133 = -t136 * t259 + t263 * t143;
t308 = qJ(2) * t133 + t122 * t263 + t124 * t259;
t137 = m(3) * t233 + t138;
t128 = -t137 * t262 + t309 * t266;
t120 = mrSges(3,1) * t203 - mrSges(3,2) * t204 + pkin(2) * t139 + t265 * t129 + t261 * t282;
t278 = mrSges(2,1) * t253 - mrSges(2,2) * t254 + pkin(1) * t128 + t266 * t120 + t308 * t262;
t131 = m(2) * t254 + t133;
t127 = t266 * t137 + t309 * t262;
t125 = m(2) * t253 + t128;
t118 = mrSges(2,2) * t258 - mrSges(2,3) * t253 - t259 * t122 + t263 * t124 + (-t127 * t262 - t128 * t266) * qJ(2);
t117 = -mrSges(2,1) * t258 + mrSges(2,3) * t254 - pkin(1) * t127 - t262 * t120 + t308 * t266;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t264 * t118 - t260 * t117 - qJ(1) * (t125 * t264 + t131 * t260), t118, t124, t130, t140, t160, -t207 * t257 + t279; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t260 * t118 + t264 * t117 + qJ(1) * (-t125 * t260 + t131 * t264), t117, t122, t134, t145, t153, -t247 * t205 + t281; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t278, t278, t120, t129, t306, t307, -t246 * t209 - t280;];
m_new  = t1;
