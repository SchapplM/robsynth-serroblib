% Calculate vector of cutting torques with Newton-Euler for
% S5RPPPR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:31:06
% EndTime: 2019-12-05 17:31:17
% DurationCPUTime: 6.93s
% Computational Cost: add. (64644->281), mult. (181142->386), div. (0->0), fcn. (122199->10), ass. (0->133)
t243 = sin(pkin(7));
t246 = cos(pkin(7));
t248 = sin(qJ(1));
t250 = cos(qJ(1));
t226 = t248 * g(2) - t250 * g(3);
t251 = qJD(1) ^ 2;
t298 = -t251 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t226;
t195 = -t243 * g(1) + t298 * t246;
t268 = -pkin(2) * t246 - qJ(3) * t243;
t218 = t268 * qJD(1);
t284 = qJD(1) * t246;
t186 = t218 * t284 + t195;
t242 = sin(pkin(8));
t245 = cos(pkin(8));
t227 = t250 * g(2) + t248 * g(3);
t259 = -t251 * qJ(2) + qJDD(2) - t227;
t285 = qJD(1) * t243;
t297 = (-pkin(1) + t268) * qJDD(1) + t259 - 0.2e1 * qJD(3) * t285;
t166 = -t242 * t186 + t297 * t245;
t194 = -t246 * g(1) - t298 * t243;
t241 = sin(pkin(9));
t244 = cos(pkin(9));
t288 = t243 * t245;
t261 = t241 * t288 + t244 * t246;
t207 = t261 * qJD(1);
t205 = t261 * qJDD(1);
t167 = t245 * t186 + t297 * t242;
t211 = (pkin(3) * t242 - qJ(4) * t245) * t285;
t279 = t242 * t285;
t281 = qJDD(1) * t246;
t290 = t246 ^ 2 * t251;
t164 = -pkin(3) * t290 - qJ(4) * t281 - t211 * t279 + t167;
t185 = t218 * t285 + qJDD(3) - t194;
t287 = t246 * t251;
t172 = ((-qJDD(1) * t245 - t242 * t287) * qJ(4) + (qJDD(1) * t242 - t245 * t287) * pkin(3)) * t243 + t185;
t295 = 2 * qJD(4);
t160 = t244 * t164 + t241 * t172 - t207 * t295;
t278 = t245 * t285;
t208 = -t241 * t284 + t244 * t278;
t188 = t207 * pkin(4) - t208 * pkin(6);
t282 = qJDD(1) * t243;
t275 = t242 * t282;
t291 = t243 ^ 2 * t251;
t280 = t242 ^ 2 * t291;
t158 = -pkin(4) * t280 + pkin(6) * t275 - t207 * t188 + t160;
t163 = pkin(3) * t281 - qJ(4) * t290 + t211 * t278 + qJDD(4) - t166;
t206 = (-t241 * t246 + t244 * t288) * qJDD(1);
t161 = (t207 * t279 - t206) * pkin(6) + (t208 * t279 + t205) * pkin(4) + t163;
t247 = sin(qJ(5));
t249 = cos(qJ(5));
t155 = -t247 * t158 + t249 * t161;
t189 = -t247 * t208 + t249 * t279;
t190 = t249 * t208 + t247 * t279;
t174 = -t189 * mrSges(6,1) + t190 * mrSges(6,2);
t176 = t189 * qJD(5) + t249 * t206 + t247 * t275;
t204 = qJD(5) + t207;
t177 = -t204 * mrSges(6,2) + t189 * mrSges(6,3);
t203 = qJDD(5) + t205;
t152 = m(6) * t155 + t203 * mrSges(6,1) - t176 * mrSges(6,3) - t190 * t174 + t204 * t177;
t156 = t249 * t158 + t247 * t161;
t175 = -t190 * qJD(5) - t247 * t206 + t249 * t275;
t178 = t204 * mrSges(6,1) - t190 * mrSges(6,3);
t153 = m(6) * t156 - t203 * mrSges(6,2) + t175 * mrSges(6,3) + t189 * t174 - t204 * t178;
t146 = t249 * t152 + t247 * t153;
t267 = t241 * t164 - t244 * t172;
t157 = -pkin(4) * t275 - pkin(6) * t280 + (t295 + t188) * t208 + t267;
t168 = Ifges(6,5) * t190 + Ifges(6,6) * t189 + Ifges(6,3) * t204;
t170 = Ifges(6,1) * t190 + Ifges(6,4) * t189 + Ifges(6,5) * t204;
t148 = -mrSges(6,1) * t157 + mrSges(6,3) * t156 + Ifges(6,4) * t176 + Ifges(6,2) * t175 + Ifges(6,6) * t203 - t190 * t168 + t204 * t170;
t169 = Ifges(6,4) * t190 + Ifges(6,2) * t189 + Ifges(6,6) * t204;
t149 = mrSges(6,2) * t157 - mrSges(6,3) * t155 + Ifges(6,1) * t176 + Ifges(6,4) * t175 + Ifges(6,5) * t203 + t189 * t168 - t204 * t169;
t159 = -0.2e1 * qJD(4) * t208 - t267;
t179 = Ifges(5,5) * t208 - Ifges(5,6) * t207 + Ifges(5,3) * t279;
t180 = Ifges(5,4) * t208 - Ifges(5,2) * t207 + Ifges(5,6) * t279;
t289 = t242 * t243;
t134 = mrSges(5,2) * t163 - mrSges(5,3) * t159 + Ifges(5,1) * t206 - Ifges(5,4) * t205 - pkin(6) * t146 - t247 * t148 + t249 * t149 - t207 * t179 + (Ifges(5,5) * qJDD(1) - qJD(1) * t180) * t289;
t181 = Ifges(5,1) * t208 - Ifges(5,4) * t207 + Ifges(5,5) * t279;
t253 = mrSges(6,1) * t155 - mrSges(6,2) * t156 + Ifges(6,5) * t176 + Ifges(6,6) * t175 + Ifges(6,3) * t203 + t190 * t169 - t189 * t170;
t135 = -mrSges(5,1) * t163 + mrSges(5,3) * t160 + Ifges(5,4) * t206 - Ifges(5,2) * t205 - pkin(4) * t146 - t208 * t179 + (Ifges(5,6) * qJDD(1) + qJD(1) * t181) * t289 - t253;
t147 = -t247 * t152 + t249 * t153;
t187 = t207 * mrSges(5,1) + t208 * mrSges(5,2);
t193 = mrSges(5,1) * t279 - t208 * mrSges(5,3);
t144 = m(5) * t160 - t205 * mrSges(5,3) - t207 * t187 + (-mrSges(5,2) * qJDD(1) - qJD(1) * t193) * t289 + t147;
t154 = -m(6) * t157 + t175 * mrSges(6,1) - t176 * mrSges(6,2) + t189 * t177 - t190 * t178;
t192 = -mrSges(5,2) * t279 - t207 * mrSges(5,3);
t150 = m(5) * t159 - t206 * mrSges(5,3) - t208 * t187 + (mrSges(5,1) * qJDD(1) + qJD(1) * t192) * t289 + t154;
t140 = t241 * t144 + t244 * t150;
t269 = Ifges(4,5) * t245 - Ifges(4,6) * t242;
t198 = (-Ifges(4,3) * t246 + t269 * t243) * qJD(1);
t292 = Ifges(4,6) * t246;
t293 = Ifges(4,4) * t245;
t199 = (-t292 + (-Ifges(4,2) * t242 + t293) * t243) * qJD(1);
t257 = -Ifges(4,5) * t246 + (Ifges(4,1) * t245 - Ifges(4,4) * t242) * t243;
t123 = mrSges(4,2) * t185 - mrSges(4,3) * t166 - qJ(4) * t140 + t244 * t134 - t241 * t135 + (-t198 * t289 + t199 * t246) * qJD(1) + t257 * qJDD(1);
t200 = t257 * qJD(1);
t252 = -mrSges(5,1) * t159 + mrSges(5,2) * t160 - Ifges(5,5) * t206 + Ifges(5,6) * t205 - pkin(4) * t154 - pkin(6) * t147 - t249 * t148 - t247 * t149 - t208 * t180 - t207 * t181;
t124 = -pkin(3) * t140 - mrSges(4,1) * t185 + t252 + mrSges(4,3) * t167 + (-t292 + (t293 + (-Ifges(4,2) - Ifges(5,3)) * t242) * t243) * qJDD(1) + (-t198 * t288 - t246 * t200) * qJD(1);
t141 = t244 * t144 - t241 * t150;
t272 = mrSges(4,1) * t242 + mrSges(4,2) * t245;
t212 = t272 * t285;
t263 = -mrSges(4,1) * t246 - mrSges(4,3) * t288;
t216 = t263 * qJD(1);
t262 = mrSges(4,2) * t246 - mrSges(4,3) * t289;
t138 = m(4) * t167 + t262 * qJDD(1) + (-t212 * t289 + t216 * t246) * qJD(1) + t141;
t145 = -m(5) * t163 - t205 * mrSges(5,1) - t206 * mrSges(5,2) - t207 * t192 - t208 * t193 - t146;
t215 = t262 * qJD(1);
t142 = m(4) * t166 + t263 * qJDD(1) + (-t212 * t288 - t215 * t246) * qJD(1) + t145;
t133 = t245 * t138 - t242 * t142;
t260 = -m(4) * t185 - t140;
t265 = -t215 * t242 - t216 * t245;
t271 = Ifges(3,1) * t243 + Ifges(3,4) * t246;
t296 = -((Ifges(3,4) * t243 + Ifges(3,2) * t246) * t285 - t271 * t284) * qJD(1) - mrSges(3,1) * t194 + mrSges(3,2) * t195 - pkin(2) * ((t265 * qJD(1) - t272 * qJDD(1)) * t243 + t260) - qJ(3) * t133 - t242 * t123 - t245 * t124;
t294 = mrSges(3,2) * t243;
t219 = (-mrSges(3,1) * t246 + t294) * qJD(1);
t130 = m(3) * t195 + (qJDD(1) * mrSges(3,3) + qJD(1) * t219) * t246 + t133;
t136 = m(3) * t194 + ((-mrSges(3,3) - t272) * qJDD(1) + (-t219 + t265) * qJD(1)) * t243 + t260;
t274 = t246 * t130 - t243 * t136;
t270 = Ifges(3,5) * t243 + Ifges(3,6) * t246;
t132 = t242 * t138 + t245 * t142;
t266 = t199 * t245 + t200 * t242;
t214 = -qJDD(1) * pkin(1) + t259;
t220 = t270 * qJD(1);
t120 = mrSges(3,2) * t214 - mrSges(3,3) * t194 - qJ(3) * t132 + t271 * qJDD(1) + t245 * t123 - t242 * t124 + t220 * t284;
t254 = mrSges(4,1) * t166 - mrSges(4,2) * t167 + pkin(3) * t145 + qJ(4) * t141 + t241 * t134 + t244 * t135;
t122 = -mrSges(3,1) * t214 + mrSges(3,3) * t195 - pkin(2) * t132 + (Ifges(3,2) + Ifges(4,3)) * t281 + ((Ifges(3,4) - t269) * qJDD(1) + (-t220 - t266) * qJD(1)) * t243 - t254;
t256 = -m(3) * t214 + mrSges(3,1) * t281 - t132 + (t290 + t291) * mrSges(3,3);
t258 = -mrSges(2,2) * t226 + qJ(2) * t274 + t243 * t120 + t246 * t122 + pkin(1) * (-mrSges(3,2) * t282 + t256) + mrSges(2,1) * t227 + Ifges(2,3) * qJDD(1);
t128 = m(2) * t227 - t251 * mrSges(2,2) + (mrSges(2,1) - t294) * qJDD(1) + t256;
t127 = t243 * t130 + t246 * t136;
t125 = m(2) * t226 - t251 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t274;
t118 = mrSges(2,1) * g(1) + mrSges(2,3) * t226 + t251 * Ifges(2,5) - pkin(1) * t127 + (Ifges(2,6) - t270) * qJDD(1) + t296;
t117 = -mrSges(2,2) * g(1) - mrSges(2,3) * t227 + Ifges(2,5) * qJDD(1) - t251 * Ifges(2,6) - qJ(2) * t127 + t246 * t120 - t243 * t122;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t258, t117, t120, t123, t134, t149; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t248 * t117 - t250 * t118 - pkin(5) * (t250 * t125 - t248 * t128), t118, t122, t124, t135, t148; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t250 * t117 - t248 * t118 + pkin(5) * (-t248 * t125 - t250 * t128), t258, t270 * qJDD(1) - t296, -Ifges(4,3) * t281 + (t266 * qJD(1) + t269 * qJDD(1)) * t243 + t254, Ifges(5,3) * t275 - t252, t253;];
m_new = t1;
