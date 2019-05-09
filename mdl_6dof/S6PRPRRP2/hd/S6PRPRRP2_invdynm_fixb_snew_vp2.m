% Calculate vector of cutting torques with Newton-Euler for
% S6PRPRRP2
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
% Datum: 2019-05-04 23:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRPRRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:35:04
% EndTime: 2019-05-04 23:35:18
% DurationCPUTime: 7.77s
% Computational Cost: add. (150951->293), mult. (274422->364), div. (0->0), fcn. (185948->12), ass. (0->122)
t267 = sin(pkin(10));
t270 = cos(pkin(10));
t253 = g(1) * t267 - g(2) * t270;
t254 = -g(1) * t270 - g(2) * t267;
t265 = -g(3) + qJDD(1);
t274 = sin(qJ(2));
t271 = cos(pkin(6));
t276 = cos(qJ(2));
t301 = t271 * t276;
t268 = sin(pkin(6));
t303 = t268 * t276;
t201 = t253 * t301 - t254 * t274 + t265 * t303;
t199 = qJDD(2) * pkin(2) + t201;
t302 = t271 * t274;
t304 = t268 * t274;
t202 = t253 * t302 + t276 * t254 + t265 * t304;
t278 = qJD(2) ^ 2;
t200 = -pkin(2) * t278 + t202;
t266 = sin(pkin(11));
t269 = cos(pkin(11));
t193 = t266 * t199 + t269 * t200;
t190 = -pkin(3) * t278 + qJDD(2) * pkin(8) + t193;
t231 = -t253 * t268 + t271 * t265;
t230 = qJDD(3) + t231;
t273 = sin(qJ(4));
t275 = cos(qJ(4));
t186 = t275 * t190 + t273 * t230;
t249 = (-pkin(4) * t275 - pkin(9) * t273) * qJD(2);
t277 = qJD(4) ^ 2;
t296 = qJD(2) * t275;
t182 = -pkin(4) * t277 + qJDD(4) * pkin(9) + t249 * t296 + t186;
t192 = t269 * t199 - t266 * t200;
t189 = -qJDD(2) * pkin(3) - t278 * pkin(8) - t192;
t295 = qJD(2) * qJD(4);
t291 = t275 * t295;
t250 = qJDD(2) * t273 + t291;
t292 = t273 * t295;
t251 = qJDD(2) * t275 - t292;
t184 = (-t250 - t291) * pkin(9) + (-t251 + t292) * pkin(4) + t189;
t272 = sin(qJ(5));
t307 = cos(qJ(5));
t177 = -t272 * t182 + t184 * t307;
t178 = t307 * t182 + t272 * t184;
t297 = qJD(2) * t273;
t246 = -qJD(4) * t307 + t272 * t297;
t247 = t272 * qJD(4) + t297 * t307;
t260 = qJD(5) - t296;
t203 = Ifges(7,5) * t247 + Ifges(7,6) * t260 + Ifges(7,3) * t246;
t206 = Ifges(6,4) * t247 - Ifges(6,2) * t246 + Ifges(6,6) * t260;
t208 = Ifges(6,1) * t247 - Ifges(6,4) * t246 + Ifges(6,5) * t260;
t217 = t247 * qJD(5) - qJDD(4) * t307 + t272 * t250;
t218 = -t246 * qJD(5) + t272 * qJDD(4) + t250 * t307;
t222 = mrSges(7,1) * t246 - mrSges(7,3) * t247;
t244 = qJDD(5) - t251;
t221 = pkin(5) * t246 - qJ(6) * t247;
t259 = t260 ^ 2;
t174 = -pkin(5) * t259 + qJ(6) * t244 + 0.2e1 * qJD(6) * t260 - t221 * t246 + t178;
t176 = -t244 * pkin(5) - t259 * qJ(6) + t247 * t221 + qJDD(6) - t177;
t207 = Ifges(7,1) * t247 + Ifges(7,4) * t260 + Ifges(7,5) * t246;
t285 = mrSges(7,1) * t176 - mrSges(7,3) * t174 - Ifges(7,4) * t218 - Ifges(7,2) * t244 - Ifges(7,6) * t217 - t246 * t207;
t229 = -mrSges(7,2) * t246 + mrSges(7,3) * t260;
t288 = -m(7) * t176 + t244 * mrSges(7,1) + t260 * t229;
t228 = -mrSges(7,1) * t260 + mrSges(7,2) * t247;
t294 = m(7) * t174 + t244 * mrSges(7,3) + t260 * t228;
t309 = -(-t206 + t203) * t247 + mrSges(6,1) * t177 - mrSges(6,2) * t178 + Ifges(6,5) * t218 - Ifges(6,6) * t217 + Ifges(6,3) * t244 + pkin(5) * (-t218 * mrSges(7,2) - t247 * t222 + t288) + qJ(6) * (-t217 * mrSges(7,2) - t246 * t222 + t294) + t246 * t208 - t285;
t185 = -t273 * t190 + t275 * t230;
t181 = -qJDD(4) * pkin(4) - t277 * pkin(9) + t249 * t297 - t185;
t179 = -0.2e1 * qJD(6) * t247 + (t246 * t260 - t218) * qJ(6) + (t247 * t260 + t217) * pkin(5) + t181;
t171 = m(7) * t179 + mrSges(7,1) * t217 - t218 * mrSges(7,3) - t247 * t228 + t229 * t246;
t287 = -mrSges(7,1) * t179 + mrSges(7,2) * t174;
t205 = Ifges(7,4) * t247 + Ifges(7,2) * t260 + Ifges(7,6) * t246;
t300 = -Ifges(6,5) * t247 + Ifges(6,6) * t246 - Ifges(6,3) * t260 - t205;
t159 = -mrSges(6,1) * t181 + mrSges(6,3) * t178 - pkin(5) * t171 + (t207 + t208) * t260 + t300 * t247 + (Ifges(6,6) - Ifges(7,6)) * t244 + (Ifges(6,4) - Ifges(7,5)) * t218 + (-Ifges(6,2) - Ifges(7,3)) * t217 + t287;
t284 = mrSges(7,2) * t176 - mrSges(7,3) * t179 + Ifges(7,1) * t218 + Ifges(7,4) * t244 + Ifges(7,5) * t217 + t260 * t203;
t160 = mrSges(6,2) * t181 - mrSges(6,3) * t177 + Ifges(6,1) * t218 - Ifges(6,4) * t217 + Ifges(6,5) * t244 - qJ(6) * t171 - t260 * t206 + t246 * t300 + t284;
t227 = mrSges(6,1) * t260 - mrSges(6,3) * t247;
t298 = -mrSges(6,1) * t246 - mrSges(6,2) * t247 - t222;
t305 = -mrSges(6,3) - mrSges(7,2);
t164 = m(6) * t178 - t244 * mrSges(6,2) + t217 * t305 - t260 * t227 + t246 * t298 + t294;
t226 = -mrSges(6,2) * t260 - mrSges(6,3) * t246;
t165 = m(6) * t177 + t244 * mrSges(6,1) + t218 * t305 + t260 * t226 + t247 * t298 + t288;
t162 = t307 * t164 - t165 * t272;
t168 = -m(6) * t181 - t217 * mrSges(6,1) - mrSges(6,2) * t218 - t246 * t226 - t227 * t247 - t171;
t236 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t273 + Ifges(5,2) * t275) * qJD(2);
t237 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t273 + Ifges(5,4) * t275) * qJD(2);
t308 = mrSges(5,1) * t185 - mrSges(5,2) * t186 + Ifges(5,5) * t250 + Ifges(5,6) * t251 + Ifges(5,3) * qJDD(4) + pkin(4) * t168 + pkin(9) * t162 + (t236 * t273 - t237 * t275) * qJD(2) + t159 * t307 + t272 * t160;
t248 = (-mrSges(5,1) * t275 + mrSges(5,2) * t273) * qJD(2);
t255 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t297;
t157 = m(5) * t186 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t251 - qJD(4) * t255 + t248 * t296 + t162;
t256 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t296;
t167 = m(5) * t185 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t250 + qJD(4) * t256 - t248 * t297 + t168;
t289 = t275 * t157 - t167 * t273;
t149 = m(4) * t193 - mrSges(4,1) * t278 - qJDD(2) * mrSges(4,2) + t289;
t161 = t272 * t164 + t165 * t307;
t281 = -m(5) * t189 + t251 * mrSges(5,1) - t250 * mrSges(5,2) - t255 * t297 + t256 * t296 - t161;
t154 = m(4) * t192 + qJDD(2) * mrSges(4,1) - t278 * mrSges(4,2) + t281;
t144 = t266 * t149 + t269 * t154;
t142 = m(3) * t201 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t278 + t144;
t290 = t269 * t149 - t154 * t266;
t143 = m(3) * t202 - mrSges(3,1) * t278 - qJDD(2) * mrSges(3,2) + t290;
t135 = -t142 * t274 + t276 * t143;
t306 = pkin(7) * t135;
t152 = t273 * t157 + t275 * t167;
t293 = m(4) * t230 + t152;
t150 = m(3) * t231 + t293;
t132 = t142 * t301 + t143 * t302 - t150 * t268;
t235 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t273 + Ifges(5,6) * t275) * qJD(2);
t138 = mrSges(5,2) * t189 - mrSges(5,3) * t185 + Ifges(5,1) * t250 + Ifges(5,4) * t251 + Ifges(5,5) * qJDD(4) - pkin(9) * t161 - qJD(4) * t236 - t272 * t159 + t160 * t307 + t235 * t296;
t146 = -mrSges(5,1) * t189 + mrSges(5,3) * t186 + Ifges(5,4) * t250 + Ifges(5,2) * t251 + Ifges(5,6) * qJDD(4) - pkin(4) * t161 + qJD(4) * t237 - t235 * t297 - t309;
t128 = mrSges(4,2) * t230 - mrSges(4,3) * t192 + Ifges(4,5) * qJDD(2) - Ifges(4,6) * t278 - pkin(8) * t152 + t138 * t275 - t146 * t273;
t136 = -mrSges(4,1) * t230 + mrSges(4,3) * t193 + t278 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t152 - t308;
t123 = -mrSges(3,1) * t231 + mrSges(3,3) * t202 + t278 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t293 + qJ(3) * t290 + t266 * t128 + t269 * t136;
t125 = mrSges(3,2) * t231 - mrSges(3,3) * t201 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t278 - qJ(3) * t144 + t128 * t269 - t136 * t266;
t282 = mrSges(4,1) * t192 - mrSges(4,2) * t193 + Ifges(4,3) * qJDD(2) + pkin(3) * t281 + pkin(8) * t289 + t273 * t138 + t275 * t146;
t127 = mrSges(3,1) * t201 - mrSges(3,2) * t202 + Ifges(3,3) * qJDD(2) + pkin(2) * t144 + t282;
t283 = mrSges(2,1) * t253 - mrSges(2,2) * t254 + pkin(1) * t132 + t123 * t303 + t125 * t304 + t271 * t127 + t268 * t306;
t133 = m(2) * t254 + t135;
t131 = t271 * t150 + (t142 * t276 + t143 * t274) * t268;
t129 = m(2) * t253 + t132;
t121 = mrSges(2,2) * t265 - mrSges(2,3) * t253 - t274 * t123 + t276 * t125 + (-t131 * t268 - t132 * t271) * pkin(7);
t120 = -mrSges(2,1) * t265 + mrSges(2,3) * t254 - pkin(1) * t131 - t268 * t127 + (t123 * t276 + t125 * t274 + t306) * t271;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t270 * t121 - t267 * t120 - qJ(1) * (t129 * t270 + t133 * t267), t121, t125, t128, t138, t160, -t205 * t246 + t284; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t267 * t121 + t270 * t120 + qJ(1) * (-t129 * t267 + t133 * t270), t120, t123, t136, t146, t159, -t247 * t203 - t285; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t283, t283, t127, t282, t308, t309, Ifges(7,5) * t218 + Ifges(7,6) * t244 + Ifges(7,3) * t217 + t247 * t205 - t260 * t207 - t287;];
m_new  = t1;
