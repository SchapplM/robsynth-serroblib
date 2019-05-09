% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-05-05 14:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRRP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:57:28
% EndTime: 2019-05-05 14:57:31
% DurationCPUTime: 2.59s
% Computational Cost: add. (32511->294), mult. (60197->333), div. (0->0), fcn. (30818->6), ass. (0->106)
t273 = sin(qJ(1));
t276 = cos(qJ(1));
t247 = t273 * g(1) - t276 * g(2);
t279 = qJD(1) ^ 2;
t221 = -qJDD(1) * pkin(1) - t279 * qJ(2) + qJDD(2) - t247;
t210 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t221;
t248 = -t276 * g(1) - t273 * g(2);
t319 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t248;
t271 = sin(qJ(5));
t274 = cos(qJ(5));
t275 = cos(qJ(4));
t309 = qJD(1) * t275;
t238 = qJD(4) * t274 - t271 * t309;
t272 = sin(qJ(4));
t308 = qJD(1) * qJD(4);
t303 = t272 * t308;
t243 = qJDD(1) * t275 - t303;
t197 = qJD(5) * t238 + qJDD(4) * t271 + t243 * t274;
t239 = qJD(4) * t271 + t274 * t309;
t202 = -mrSges(7,1) * t238 + mrSges(7,2) * t239;
t205 = -t279 * pkin(7) - t210;
t302 = t275 * t308;
t242 = -qJDD(1) * t272 - t302;
t169 = (-t243 + t303) * pkin(8) + (-t242 + t302) * pkin(4) + t205;
t211 = qJDD(3) + (-pkin(1) - qJ(3)) * t279 + t319;
t206 = -qJDD(1) * pkin(7) + t211;
t199 = -g(3) * t275 + t272 * t206;
t241 = (pkin(4) * t272 - pkin(8) * t275) * qJD(1);
t278 = qJD(4) ^ 2;
t310 = qJD(1) * t272;
t172 = -pkin(4) * t278 + qJDD(4) * pkin(8) - t241 * t310 + t199;
t165 = t274 * t169 - t271 * t172;
t237 = qJDD(5) - t242;
t249 = qJD(5) + t310;
t159 = -0.2e1 * qJD(6) * t239 + (t238 * t249 - t197) * qJ(6) + (t238 * t239 + t237) * pkin(5) + t165;
t212 = -mrSges(7,2) * t249 + mrSges(7,3) * t238;
t305 = m(7) * t159 + t237 * mrSges(7,1) + t249 * t212;
t156 = -t197 * mrSges(7,3) - t239 * t202 + t305;
t166 = t271 * t169 + t274 * t172;
t180 = Ifges(6,4) * t239 + Ifges(6,2) * t238 + Ifges(6,6) * t249;
t181 = Ifges(7,1) * t239 + Ifges(7,4) * t238 + Ifges(7,5) * t249;
t182 = Ifges(6,1) * t239 + Ifges(6,4) * t238 + Ifges(6,5) * t249;
t196 = -qJD(5) * t239 + qJDD(4) * t274 - t243 * t271;
t214 = pkin(5) * t249 - qJ(6) * t239;
t236 = t238 ^ 2;
t162 = -pkin(5) * t236 + qJ(6) * t196 + 0.2e1 * qJD(6) * t238 - t214 * t249 + t166;
t179 = Ifges(7,4) * t239 + Ifges(7,2) * t238 + Ifges(7,6) * t249;
t293 = -mrSges(7,1) * t159 + mrSges(7,2) * t162 - Ifges(7,5) * t197 - Ifges(7,6) * t196 - Ifges(7,3) * t237 - t239 * t179;
t318 = mrSges(6,1) * t165 - mrSges(6,2) * t166 + Ifges(6,5) * t197 + Ifges(6,6) * t196 + Ifges(6,3) * t237 + pkin(5) * t156 + t239 * t180 - (t182 + t181) * t238 - t293;
t317 = mrSges(3,2) - mrSges(4,3);
t316 = -mrSges(6,2) - mrSges(7,2);
t315 = Ifges(3,4) - Ifges(4,5);
t314 = Ifges(2,6) - Ifges(3,5);
t313 = mrSges(4,3) * t279;
t240 = (mrSges(5,1) * t272 + mrSges(5,2) * t275) * qJD(1);
t246 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t309;
t203 = -mrSges(6,1) * t238 + mrSges(6,2) * t239;
t213 = -mrSges(6,2) * t249 + mrSges(6,3) * t238;
t148 = m(6) * t165 + t237 * mrSges(6,1) + t249 * t213 + (-t202 - t203) * t239 + (-mrSges(6,3) - mrSges(7,3)) * t197 + t305;
t304 = m(7) * t162 + t196 * mrSges(7,3) + t238 * t202;
t215 = mrSges(7,1) * t249 - mrSges(7,3) * t239;
t311 = -mrSges(6,1) * t249 + mrSges(6,3) * t239 - t215;
t151 = m(6) * t166 + t196 * mrSges(6,3) + t238 * t203 + t316 * t237 + t311 * t249 + t304;
t300 = -t148 * t271 + t274 * t151;
t142 = m(5) * t199 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t242 - qJD(4) * t246 - t240 * t310 + t300;
t198 = g(3) * t272 + t206 * t275;
t245 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t310;
t171 = -qJDD(4) * pkin(4) - pkin(8) * t278 + t241 * t309 - t198;
t164 = -pkin(5) * t196 - qJ(6) * t236 + t214 * t239 + qJDD(6) + t171;
t298 = -m(7) * t164 + t196 * mrSges(7,1) + t238 * t212;
t285 = -m(6) * t171 + t196 * mrSges(6,1) + t316 * t197 + t238 * t213 + t311 * t239 + t298;
t153 = m(5) * t198 + qJDD(4) * mrSges(5,1) - t243 * mrSges(5,3) + qJD(4) * t245 - t240 * t309 + t285;
t130 = t272 * t142 + t275 * t153;
t145 = t274 * t148 + t271 * t151;
t301 = t275 * t142 - t272 * t153;
t299 = m(4) * t211 + qJDD(1) * mrSges(4,2) + t130;
t295 = -m(5) * t205 + t242 * mrSges(5,1) - t243 * mrSges(5,2) - t245 * t310 - t246 * t309 - t145;
t294 = -mrSges(7,1) * t164 + mrSges(7,3) * t162 + Ifges(7,4) * t197 + Ifges(7,2) * t196 + Ifges(7,6) * t237 + t249 * t181;
t177 = Ifges(7,5) * t239 + Ifges(7,6) * t238 + Ifges(7,3) * t249;
t292 = mrSges(7,2) * t164 - mrSges(7,3) * t159 + Ifges(7,1) * t197 + Ifges(7,4) * t196 + Ifges(7,5) * t237 + t238 * t177;
t178 = Ifges(6,5) * t239 + Ifges(6,6) * t238 + Ifges(6,3) * t249;
t133 = Ifges(6,4) * t197 + Ifges(6,2) * t196 + Ifges(6,6) * t237 + t249 * t182 - mrSges(6,1) * t171 + mrSges(6,3) * t166 - pkin(5) * (t197 * mrSges(7,2) - t298) + qJ(6) * (-t237 * mrSges(7,2) - t249 * t215 + t304) + (-pkin(5) * t215 - t177 - t178) * t239 + t294;
t143 = mrSges(6,2) * t171 - mrSges(6,3) * t165 + Ifges(6,1) * t197 + Ifges(6,4) * t196 + Ifges(6,5) * t237 - qJ(6) * t156 + t238 * t178 + (-t179 - t180) * t249 + t292;
t223 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t275 - Ifges(5,6) * t272) * qJD(1);
t224 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t275 - Ifges(5,2) * t272) * qJD(1);
t121 = mrSges(5,2) * t205 - mrSges(5,3) * t198 + Ifges(5,1) * t243 + Ifges(5,4) * t242 + Ifges(5,5) * qJDD(4) - pkin(8) * t145 - qJD(4) * t224 - t133 * t271 + t143 * t274 - t223 * t310;
t225 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t275 - Ifges(5,4) * t272) * qJD(1);
t123 = -mrSges(5,1) * t205 + mrSges(5,3) * t199 + Ifges(5,4) * t243 + Ifges(5,2) * t242 + Ifges(5,6) * qJDD(4) - pkin(4) * t145 + qJD(4) * t225 - t223 * t309 - t318;
t291 = mrSges(4,1) * t210 + mrSges(4,2) * g(3) + t279 * Ifges(4,4) + Ifges(4,5) * qJDD(1) + pkin(3) * t295 + pkin(7) * t301 + t272 * t121 + t275 * t123;
t219 = pkin(1) * t279 - t319;
t290 = -m(3) * t219 + t279 * mrSges(3,2) + qJDD(1) * mrSges(3,3) + t299;
t289 = mrSges(4,2) * t211 - mrSges(4,3) * t210 + Ifges(4,1) * qJDD(1) - pkin(7) * t130 + t275 * t121 - t123 * t272;
t288 = mrSges(5,1) * t198 - mrSges(5,2) * t199 + Ifges(5,5) * t243 + Ifges(5,6) * t242 + Ifges(5,3) * qJDD(4) + pkin(4) * t285 + pkin(8) * t300 + t274 * t133 + t271 * t143 + t224 * t309 + t225 * t310;
t136 = m(4) * t210 - t279 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t295;
t287 = mrSges(3,1) * t221 + pkin(2) * t136 + t291;
t286 = -m(3) * t221 + t279 * mrSges(3,3) - t136;
t284 = mrSges(4,1) * t211 - Ifges(4,4) * qJDD(1) + pkin(3) * t130 + t288;
t283 = mrSges(3,2) * t221 - mrSges(3,3) * t219 + Ifges(3,1) * qJDD(1) - qJ(3) * t136 + t289;
t282 = -mrSges(2,2) * t248 + qJ(2) * (t290 - t313) + pkin(1) * (-qJDD(1) * mrSges(3,2) + t286) + mrSges(2,1) * t247 + Ifges(2,3) * qJDD(1) + t283;
t281 = mrSges(3,1) * t219 + pkin(2) * (-t299 + t313) + qJ(3) * (-m(4) * g(3) + t301) - t284;
t134 = t286 + (mrSges(2,1) - mrSges(3,2)) * qJDD(1) - t279 * mrSges(2,2) + m(2) * t247;
t127 = (-m(3) - m(4)) * g(3) + t301;
t124 = m(2) * t248 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t279 + t290;
t118 = -t281 + (Ifges(2,5) - t315) * t279 + t314 * qJDD(1) + (mrSges(2,1) - t317) * g(3) + mrSges(2,3) * t248 - pkin(1) * t127;
t117 = t287 - t314 * t279 + (-Ifges(3,4) + Ifges(2,5)) * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t247 - qJ(2) * t127;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t276 * t117 - t273 * t118 - pkin(6) * (t124 * t273 + t134 * t276), t117, t283, t289, t121, t143, -t179 * t249 + t292; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t273 * t117 + t276 * t118 + pkin(6) * (t124 * t276 - t134 * t273), t118, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t279 * Ifges(3,5) - t287, -mrSges(4,3) * g(3) - t279 * Ifges(4,5) - t284, t123, t133, -t239 * t177 + t294; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t282, t282, Ifges(3,5) * qJDD(1) + t317 * g(3) + t315 * t279 + t281, t291, t288, t318, -t238 * t181 - t293;];
m_new  = t1;
