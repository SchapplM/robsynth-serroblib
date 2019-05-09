% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRRP1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-05-04 23:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:29:14
% EndTime: 2019-05-04 23:29:27
% DurationCPUTime: 8.06s
% Computational Cost: add. (154584->292), mult. (282974->364), div. (0->0), fcn. (192222->12), ass. (0->123)
t274 = sin(qJ(5));
t277 = cos(qJ(5));
t275 = sin(qJ(4));
t301 = qJD(2) * t275;
t251 = qJD(4) * t277 - t274 * t301;
t278 = cos(qJ(4));
t299 = qJD(2) * qJD(4);
t294 = t278 * t299;
t255 = qJDD(2) * t275 + t294;
t223 = qJD(5) * t251 + qJDD(4) * t274 + t255 * t277;
t252 = qJD(4) * t274 + t277 * t301;
t225 = -mrSges(7,1) * t251 + mrSges(7,2) * t252;
t269 = sin(pkin(10));
t272 = cos(pkin(10));
t258 = g(1) * t269 - g(2) * t272;
t259 = -g(1) * t272 - g(2) * t269;
t267 = -g(3) + qJDD(1);
t276 = sin(qJ(2));
t273 = cos(pkin(6));
t279 = cos(qJ(2));
t304 = t273 * t279;
t270 = sin(pkin(6));
t306 = t270 * t279;
t203 = t258 * t304 - t259 * t276 + t267 * t306;
t201 = qJDD(2) * pkin(2) + t203;
t305 = t273 * t276;
t307 = t270 * t276;
t204 = t258 * t305 + t279 * t259 + t267 * t307;
t281 = qJD(2) ^ 2;
t202 = -pkin(2) * t281 + t204;
t268 = sin(pkin(11));
t271 = cos(pkin(11));
t193 = t268 * t201 + t271 * t202;
t190 = -pkin(3) * t281 + qJDD(2) * pkin(8) + t193;
t235 = -t258 * t270 + t273 * t267;
t234 = qJDD(3) + t235;
t186 = t278 * t190 + t275 * t234;
t254 = (-pkin(4) * t278 - pkin(9) * t275) * qJD(2);
t280 = qJD(4) ^ 2;
t300 = qJD(2) * t278;
t181 = -pkin(4) * t280 + qJDD(4) * pkin(9) + t254 * t300 + t186;
t192 = t271 * t201 - t268 * t202;
t189 = -qJDD(2) * pkin(3) - t281 * pkin(8) - t192;
t295 = t275 * t299;
t256 = qJDD(2) * t278 - t295;
t184 = (-t255 - t294) * pkin(9) + (-t256 + t295) * pkin(4) + t189;
t175 = -t274 * t181 + t277 * t184;
t249 = qJDD(5) - t256;
t264 = qJD(5) - t300;
t171 = -0.2e1 * qJD(6) * t252 + (t251 * t264 - t223) * qJ(6) + (t251 * t252 + t249) * pkin(5) + t175;
t229 = -mrSges(7,2) * t264 + mrSges(7,3) * t251;
t298 = m(7) * t171 + t249 * mrSges(7,1) + t264 * t229;
t168 = -t223 * mrSges(7,3) - t252 * t225 + t298;
t176 = t277 * t181 + t274 * t184;
t208 = Ifges(6,4) * t252 + Ifges(6,2) * t251 + Ifges(6,6) * t264;
t209 = Ifges(7,1) * t252 + Ifges(7,4) * t251 + Ifges(7,5) * t264;
t210 = Ifges(6,1) * t252 + Ifges(6,4) * t251 + Ifges(6,5) * t264;
t222 = -qJD(5) * t252 + qJDD(4) * t277 - t255 * t274;
t231 = pkin(5) * t264 - qJ(6) * t252;
t248 = t251 ^ 2;
t174 = -pkin(5) * t248 + qJ(6) * t222 + 0.2e1 * qJD(6) * t251 - t231 * t264 + t176;
t207 = Ifges(7,4) * t252 + Ifges(7,2) * t251 + Ifges(7,6) * t264;
t288 = -mrSges(7,1) * t171 + mrSges(7,2) * t174 - Ifges(7,5) * t223 - Ifges(7,6) * t222 - Ifges(7,3) * t249 - t252 * t207;
t311 = mrSges(6,1) * t175 - mrSges(6,2) * t176 + Ifges(6,5) * t223 + Ifges(6,6) * t222 + Ifges(6,3) * t249 + pkin(5) * t168 + t252 * t208 - (t210 + t209) * t251 - t288;
t185 = -t275 * t190 + t234 * t278;
t180 = -qJDD(4) * pkin(4) - pkin(9) * t280 + t254 * t301 - t185;
t205 = Ifges(7,5) * t252 + Ifges(7,6) * t251 + Ifges(7,3) * t264;
t206 = Ifges(6,5) * t252 + Ifges(6,6) * t251 + Ifges(6,3) * t264;
t232 = mrSges(7,1) * t264 - mrSges(7,3) * t252;
t178 = -pkin(5) * t222 - qJ(6) * t248 + t231 * t252 + qJDD(6) + t180;
t289 = -mrSges(7,1) * t178 + mrSges(7,3) * t174 + Ifges(7,4) * t223 + Ifges(7,2) * t222 + Ifges(7,6) * t249 + t264 * t209;
t291 = -m(7) * t178 + t222 * mrSges(7,1) + t251 * t229;
t297 = m(7) * t174 + t222 * mrSges(7,3) + t251 * t225;
t152 = Ifges(6,4) * t223 + Ifges(6,2) * t222 + Ifges(6,6) * t249 + t264 * t210 - mrSges(6,1) * t180 + mrSges(6,3) * t176 - pkin(5) * (t223 * mrSges(7,2) - t291) + qJ(6) * (-t249 * mrSges(7,2) - t264 * t232 + t297) + (-pkin(5) * t232 - t205 - t206) * t252 + t289;
t287 = mrSges(7,2) * t178 - mrSges(7,3) * t171 + Ifges(7,1) * t223 + Ifges(7,4) * t222 + Ifges(7,5) * t249 + t251 * t205;
t159 = mrSges(6,2) * t180 - mrSges(6,3) * t175 + Ifges(6,1) * t223 + Ifges(6,4) * t222 + Ifges(6,5) * t249 - qJ(6) * t168 + t251 * t206 + (-t207 - t208) * t264 + t287;
t226 = -mrSges(6,1) * t251 + mrSges(6,2) * t252;
t230 = -mrSges(6,2) * t264 + mrSges(6,3) * t251;
t162 = m(6) * t175 + t249 * mrSges(6,1) + t264 * t230 + (-t225 - t226) * t252 + (-mrSges(6,3) - mrSges(7,3)) * t223 + t298;
t302 = -mrSges(6,1) * t264 + mrSges(6,3) * t252 - t232;
t308 = -mrSges(6,2) - mrSges(7,2);
t164 = m(6) * t176 + t222 * mrSges(6,3) + t251 * t226 + t308 * t249 + t302 * t264 + t297;
t161 = -t162 * t274 + t277 * t164;
t167 = -m(6) * t180 + t222 * mrSges(6,1) + t308 * t223 + t251 * t230 + t302 * t252 + t291;
t240 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t275 + Ifges(5,2) * t278) * qJD(2);
t241 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t275 + Ifges(5,4) * t278) * qJD(2);
t310 = mrSges(5,1) * t185 - mrSges(5,2) * t186 + Ifges(5,5) * t255 + Ifges(5,6) * t256 + Ifges(5,3) * qJDD(4) + pkin(4) * t167 + pkin(9) * t161 + t277 * t152 + t274 * t159 + (t240 * t275 - t241 * t278) * qJD(2);
t253 = (-mrSges(5,1) * t278 + mrSges(5,2) * t275) * qJD(2);
t260 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t301;
t157 = m(5) * t186 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t256 - qJD(4) * t260 + t253 * t300 + t161;
t261 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t300;
t166 = m(5) * t185 + qJDD(4) * mrSges(5,1) - t255 * mrSges(5,3) + qJD(4) * t261 - t253 * t301 + t167;
t292 = t278 * t157 - t166 * t275;
t148 = m(4) * t193 - mrSges(4,1) * t281 - qJDD(2) * mrSges(4,2) + t292;
t160 = t162 * t277 + t164 * t274;
t284 = -m(5) * t189 + t256 * mrSges(5,1) - mrSges(5,2) * t255 - t260 * t301 + t261 * t300 - t160;
t154 = m(4) * t192 + qJDD(2) * mrSges(4,1) - mrSges(4,2) * t281 + t284;
t143 = t268 * t148 + t271 * t154;
t141 = m(3) * t203 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t281 + t143;
t293 = t271 * t148 - t154 * t268;
t142 = m(3) * t204 - mrSges(3,1) * t281 - qJDD(2) * mrSges(3,2) + t293;
t134 = -t141 * t276 + t279 * t142;
t309 = pkin(7) * t134;
t151 = t275 * t157 + t278 * t166;
t296 = m(4) * t234 + t151;
t149 = m(3) * t235 + t296;
t131 = t141 * t304 + t142 * t305 - t149 * t270;
t239 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t275 + Ifges(5,6) * t278) * qJD(2);
t137 = mrSges(5,2) * t189 - mrSges(5,3) * t185 + Ifges(5,1) * t255 + Ifges(5,4) * t256 + Ifges(5,5) * qJDD(4) - pkin(9) * t160 - qJD(4) * t240 - t152 * t274 + t159 * t277 + t239 * t300;
t145 = -mrSges(5,1) * t189 + mrSges(5,3) * t186 + Ifges(5,4) * t255 + Ifges(5,2) * t256 + Ifges(5,6) * qJDD(4) - pkin(4) * t160 + qJD(4) * t241 - t239 * t301 - t311;
t127 = mrSges(4,2) * t234 - mrSges(4,3) * t192 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t281 - pkin(8) * t151 + t137 * t278 - t145 * t275;
t135 = -mrSges(4,1) * t234 + mrSges(4,3) * t193 + t281 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t151 - t310;
t122 = -mrSges(3,1) * t235 + mrSges(3,3) * t204 + t281 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t296 + qJ(3) * t293 + t268 * t127 + t271 * t135;
t124 = mrSges(3,2) * t235 - mrSges(3,3) * t203 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t281 - qJ(3) * t143 + t127 * t271 - t135 * t268;
t285 = mrSges(4,1) * t192 - mrSges(4,2) * t193 + Ifges(4,3) * qJDD(2) + pkin(3) * t284 + pkin(8) * t292 + t275 * t137 + t278 * t145;
t126 = mrSges(3,1) * t203 - mrSges(3,2) * t204 + Ifges(3,3) * qJDD(2) + pkin(2) * t143 + t285;
t286 = mrSges(2,1) * t258 - mrSges(2,2) * t259 + pkin(1) * t131 + t122 * t306 + t124 * t307 + t273 * t126 + t270 * t309;
t132 = m(2) * t259 + t134;
t130 = t273 * t149 + (t141 * t279 + t142 * t276) * t270;
t128 = m(2) * t258 + t131;
t120 = mrSges(2,2) * t267 - mrSges(2,3) * t258 - t276 * t122 + t279 * t124 + (-t130 * t270 - t131 * t273) * pkin(7);
t119 = -mrSges(2,1) * t267 + mrSges(2,3) * t259 - pkin(1) * t130 - t270 * t126 + (t122 * t279 + t124 * t276 + t309) * t273;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t272 * t120 - t269 * t119 - qJ(1) * (t128 * t272 + t132 * t269), t120, t124, t127, t137, t159, -t207 * t264 + t287; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t269 * t120 + t272 * t119 + qJ(1) * (-t128 * t269 + t132 * t272), t119, t122, t135, t145, t152, -t252 * t205 + t289; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t286, t286, t126, t285, t310, t311, -t251 * t209 - t288;];
m_new  = t1;
