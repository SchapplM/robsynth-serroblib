% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-05-04 23:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRPR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:18:58
% EndTime: 2019-05-04 23:19:07
% DurationCPUTime: 4.48s
% Computational Cost: add. (58179->300), mult. (109417->362), div. (0->0), fcn. (61563->10), ass. (0->126)
t267 = sin(pkin(10));
t269 = cos(pkin(10));
t243 = g(1) * t267 - g(2) * t269;
t244 = -g(1) * t269 - g(2) * t267;
t264 = -g(3) + qJDD(1);
t268 = sin(pkin(6));
t270 = cos(pkin(6));
t273 = sin(qJ(2));
t276 = cos(qJ(2));
t193 = -t273 * t244 + (t243 * t270 + t264 * t268) * t276;
t278 = qJD(2) ^ 2;
t288 = -t278 * qJ(3) + qJDD(3) - t193;
t323 = -pkin(2) - pkin(8);
t189 = t323 * qJDD(2) + t288;
t210 = -t243 * t268 + t264 * t270;
t272 = sin(qJ(4));
t275 = cos(qJ(4));
t183 = t275 * t189 - t272 * t210;
t308 = qJD(2) * qJD(4);
t255 = t272 * t308;
t240 = qJDD(2) * t275 - t255;
t310 = qJD(2) * t272;
t245 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t310;
t247 = mrSges(6,1) * t310 - qJD(4) * mrSges(6,3);
t236 = (pkin(4) * t272 - qJ(5) * t275) * qJD(2);
t277 = qJD(4) ^ 2;
t309 = qJD(2) * t275;
t179 = -qJDD(4) * pkin(4) - t277 * qJ(5) + t236 * t309 + qJDD(5) - t183;
t175 = (t272 * t275 * t278 - qJDD(4)) * pkin(9) + (t240 + t255) * pkin(5) + t179;
t305 = t275 * t308;
t239 = qJDD(2) * t272 + t305;
t250 = pkin(5) * t309 - qJD(4) * pkin(9);
t263 = t272 ^ 2;
t311 = t270 * t273;
t312 = t268 * t273;
t194 = t243 * t311 + t276 * t244 + t264 * t312;
t301 = qJDD(2) * qJ(3) + 0.2e1 * qJD(3) * qJD(2) + t194;
t324 = -2 * qJD(5);
t290 = pkin(4) * t305 + t309 * t324 + t301 + (-t240 + t255) * qJ(5);
t176 = -t250 * t309 + (pkin(4) + pkin(9)) * t239 + (-pkin(5) * t263 + t323) * t278 + t290;
t271 = sin(qJ(6));
t274 = cos(qJ(6));
t171 = t175 * t274 - t176 * t271;
t234 = -qJD(4) * t271 + t274 * t310;
t203 = qJD(6) * t234 + qJDD(4) * t274 + t239 * t271;
t235 = qJD(4) * t274 + t271 * t310;
t204 = -mrSges(7,1) * t234 + mrSges(7,2) * t235;
t253 = qJD(6) + t309;
t207 = -mrSges(7,2) * t253 + mrSges(7,3) * t234;
t231 = qJDD(6) + t240;
t167 = m(7) * t171 + mrSges(7,1) * t231 - mrSges(7,3) * t203 - t204 * t235 + t207 * t253;
t172 = t175 * t271 + t176 * t274;
t202 = -qJD(6) * t235 - qJDD(4) * t271 + t239 * t274;
t208 = mrSges(7,1) * t253 - mrSges(7,3) * t235;
t168 = m(7) * t172 - mrSges(7,2) * t231 + mrSges(7,3) * t202 + t204 * t234 - t208 * t253;
t157 = t274 * t167 + t271 * t168;
t294 = -m(6) * t179 - t240 * mrSges(6,1) - t157;
t237 = (-mrSges(6,2) * t272 - mrSges(6,3) * t275) * qJD(2);
t303 = qJD(2) * (-t237 - (mrSges(5,1) * t272 + mrSges(5,2) * t275) * qJD(2));
t319 = mrSges(5,1) - mrSges(6,2);
t152 = m(5) * t183 - t240 * mrSges(5,3) + t319 * qJDD(4) + (t245 - t247) * qJD(4) + t275 * t303 + t294;
t184 = t272 * t189 + t275 * t210;
t246 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t309;
t289 = -t277 * pkin(4) + qJDD(4) * qJ(5) - t236 * t310 + t184;
t174 = -t263 * t278 * pkin(9) - t239 * pkin(5) + ((2 * qJD(5)) + t250) * qJD(4) + t289;
t169 = -m(7) * t174 + t202 * mrSges(7,1) - t203 * mrSges(7,2) + t234 * t207 - t235 * t208;
t177 = qJD(4) * t324 - t289;
t248 = mrSges(6,1) * t309 + qJD(4) * mrSges(6,2);
t286 = -m(6) * t177 + qJDD(4) * mrSges(6,3) + qJD(4) * t248 - t169;
t163 = m(5) * t184 - qJDD(4) * mrSges(5,2) - qJD(4) * t246 + (-mrSges(5,3) - mrSges(6,1)) * t239 + t272 * t303 + t286;
t147 = t275 * t152 + t272 * t163;
t192 = -qJDD(2) * pkin(2) + t288;
t216 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t275 - Ifges(5,2) * t272) * qJD(2);
t217 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t275 - Ifges(5,4) * t272) * qJD(2);
t195 = Ifges(7,5) * t235 + Ifges(7,6) * t234 + Ifges(7,3) * t253;
t197 = Ifges(7,1) * t235 + Ifges(7,4) * t234 + Ifges(7,5) * t253;
t160 = -mrSges(7,1) * t174 + mrSges(7,3) * t172 + Ifges(7,4) * t203 + Ifges(7,2) * t202 + Ifges(7,6) * t231 - t195 * t235 + t197 * t253;
t196 = Ifges(7,4) * t235 + Ifges(7,2) * t234 + Ifges(7,6) * t253;
t161 = mrSges(7,2) * t174 - mrSges(7,3) * t171 + Ifges(7,1) * t203 + Ifges(7,4) * t202 + Ifges(7,5) * t231 + t195 * t234 - t196 * t253;
t218 = Ifges(6,5) * qJD(4) + (-Ifges(6,6) * t275 + Ifges(6,3) * t272) * qJD(2);
t219 = Ifges(6,4) * qJD(4) + (-Ifges(6,2) * t275 + Ifges(6,6) * t272) * qJD(2);
t328 = mrSges(6,2) * t179 - mrSges(6,3) * t177 + Ifges(6,1) * qJDD(4) - Ifges(6,4) * t240 + Ifges(6,5) * t239 - pkin(9) * t157 - t271 * t160 + t274 * t161 - (t218 * t275 + t219 * t272) * qJD(2);
t279 = -mrSges(5,2) * t184 + pkin(4) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t247 - t237 * t309 + t294) + qJ(5) * (-t239 * mrSges(6,1) - t237 * t310 + t286) + mrSges(5,1) * t183 + t217 * t310 + t216 * t309 - Ifges(5,6) * t239 + Ifges(5,5) * t240 + Ifges(5,3) * qJDD(4) + t328;
t330 = mrSges(4,1) * t192 + pkin(3) * t147 + t279;
t306 = t323 * t278;
t187 = t306 + t301;
t158 = -t271 * t167 + t274 * t168;
t181 = t239 * pkin(4) + t290 + t306;
t284 = m(6) * t181 - t240 * mrSges(6,3) - (t247 * t272 + t248 * t275) * qJD(2) + t158;
t329 = -m(5) * t187 - t240 * mrSges(5,2) - t319 * t239 - t245 * t310 - t246 * t309 - t284;
t148 = -t152 * t272 + t275 * t163;
t146 = m(4) * t210 + t148;
t154 = -t239 * mrSges(6,2) + t284;
t283 = -mrSges(6,1) * t177 + mrSges(6,2) * t181 - pkin(5) * t169 - pkin(9) * t158 - t274 * t160 - t271 * t161;
t220 = Ifges(6,1) * qJD(4) + (-Ifges(6,4) * t275 + Ifges(6,5) * t272) * qJD(2);
t304 = qJD(2) * (-Ifges(5,3) * qJD(4) - (Ifges(5,5) * t275 - Ifges(5,6) * t272) * qJD(2) - t220);
t318 = Ifges(5,4) + Ifges(6,6);
t137 = -mrSges(5,1) * t187 + mrSges(5,3) * t184 - pkin(4) * t154 + t318 * t240 + (-Ifges(5,2) - Ifges(6,3)) * t239 + (Ifges(5,6) - Ifges(6,5)) * qJDD(4) + (t217 - t219) * qJD(4) + t275 * t304 + t283;
t292 = mrSges(7,1) * t171 - mrSges(7,2) * t172 + Ifges(7,5) * t203 + Ifges(7,6) * t202 + Ifges(7,3) * t231 + t235 * t196 - t234 * t197;
t282 = mrSges(6,1) * t179 - mrSges(6,3) * t181 + pkin(5) * t157 + t292;
t139 = -mrSges(5,3) * t183 + mrSges(5,2) * t187 - qJ(5) * t154 - t318 * t239 + (-t216 + t218) * qJD(4) + (Ifges(5,5) - Ifges(6,4)) * qJDD(4) + (Ifges(5,1) + Ifges(6,2)) * t240 + t272 * t304 + t282;
t190 = t278 * pkin(2) - t301;
t287 = -mrSges(4,1) * t190 - pkin(3) * t329 - pkin(8) * t148 - t275 * t137 - t272 * t139;
t316 = -Ifges(3,6) + Ifges(4,5);
t317 = Ifges(3,5) - Ifges(4,4);
t320 = mrSges(3,1) - mrSges(4,2);
t130 = mrSges(3,3) * t194 - pkin(2) * t146 - t316 * qJDD(2) - t320 * t210 + t317 * t278 + t287;
t295 = -m(4) * t192 + t278 * mrSges(4,3) - t147;
t144 = m(3) * t193 - t278 * mrSges(3,2) + t320 * qJDD(2) + t295;
t280 = -m(4) * t190 + t278 * mrSges(4,2) + qJDD(2) * mrSges(4,3) - t329;
t151 = m(3) * t194 - t278 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t280;
t142 = -t144 * t273 + t276 * t151;
t327 = pkin(7) * t142 + t130 * t276;
t313 = t144 * t276;
t145 = m(3) * t210 + t146;
t136 = -t145 * t268 + t151 * t311 + t270 * t313;
t291 = mrSges(4,2) * t192 - mrSges(4,3) * t190 + Ifges(4,1) * qJDD(2) - pkin(8) * t147 - t272 * t137 + t275 * t139;
t128 = Ifges(3,3) * qJDD(2) + mrSges(3,1) * t193 - mrSges(3,2) * t194 + pkin(2) * (-qJDD(2) * mrSges(4,2) + t295) + qJ(3) * t280 + t291;
t132 = -mrSges(3,3) * t193 - qJ(3) * t146 + t317 * qJDD(2) + (mrSges(3,2) - mrSges(4,3)) * t210 + t316 * t278 + t330;
t293 = mrSges(2,1) * t243 - mrSges(2,2) * t244 + pkin(1) * t136 + t270 * t128 + t132 * t312 + t327 * t268;
t140 = m(2) * t244 + t142;
t135 = t270 * t145 + (t151 * t273 + t313) * t268;
t133 = m(2) * t243 + t136;
t126 = mrSges(2,2) * t264 - mrSges(2,3) * t243 - t273 * t130 + t276 * t132 + (-t135 * t268 - t136 * t270) * pkin(7);
t125 = -mrSges(2,1) * t264 + mrSges(2,3) * t244 - pkin(1) * t135 - t268 * t128 + (t132 * t273 + t327) * t270;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t269 * t126 - t267 * t125 - qJ(1) * (t133 * t269 + t140 * t267), t126, t132, t291, t139, t328, t161; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t267 * t126 + t269 * t125 + qJ(1) * (-t133 * t267 + t140 * t269), t125, t130, mrSges(4,3) * t210 + Ifges(4,4) * qJDD(2) - t278 * Ifges(4,5) - t330, t137, Ifges(6,4) * qJDD(4) - Ifges(6,2) * t240 + Ifges(6,6) * t239 - qJD(4) * t218 + t220 * t310 - t282, t160; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t293, t293, t128, -mrSges(4,2) * t210 + t278 * Ifges(4,4) + Ifges(4,5) * qJDD(2) - t287, t279, Ifges(6,5) * qJDD(4) - Ifges(6,6) * t240 + Ifges(6,3) * t239 + qJD(4) * t219 + t220 * t309 - t283, t292;];
m_new  = t1;
