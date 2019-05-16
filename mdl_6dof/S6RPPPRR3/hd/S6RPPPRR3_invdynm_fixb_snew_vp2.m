% Calculate vector of cutting torques with Newton-Euler for
% S6RPPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Datum: 2019-05-05 13:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPPRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR3_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR3_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR3_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:40:47
% EndTime: 2019-05-05 13:40:57
% DurationCPUTime: 7.26s
% Computational Cost: add. (107989->277), mult. (218020->338), div. (0->0), fcn. (124946->10), ass. (0->124)
t260 = qJD(1) ^ 2;
t248 = sin(pkin(10));
t250 = cos(pkin(10));
t253 = sin(qJ(5));
t256 = cos(qJ(5));
t278 = t248 * t253 - t250 * t256;
t215 = t278 * qJD(1);
t254 = sin(qJ(1));
t257 = cos(qJ(1));
t223 = -t257 * g(1) - t254 * g(2);
t274 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t223;
t296 = -pkin(1) - pkin(2);
t205 = t296 * t260 + t274;
t222 = t254 * g(1) - t257 * g(2);
t273 = -t260 * qJ(2) + qJDD(2) - t222;
t210 = t296 * qJDD(1) + t273;
t249 = sin(pkin(9));
t251 = cos(pkin(9));
t185 = t251 * t205 + t249 * t210;
t180 = -t260 * pkin(3) - qJDD(1) * qJ(4) + t185;
t244 = g(3) + qJDD(3);
t287 = qJD(1) * qJD(4);
t291 = t250 * t244 + 0.2e1 * t248 * t287;
t166 = (pkin(4) * t250 * t260 + pkin(7) * qJDD(1) - t180) * t248 + t291;
t170 = t248 * t244 + (t180 - 0.2e1 * t287) * t250;
t286 = qJDD(1) * t250;
t235 = t250 ^ 2;
t292 = t235 * t260;
t167 = -pkin(4) * t292 - pkin(7) * t286 + t170;
t163 = t253 * t166 + t256 * t167;
t279 = -t248 * t256 - t250 * t253;
t216 = t279 * qJD(1);
t192 = -t215 * mrSges(6,1) + t216 * mrSges(6,2);
t288 = t216 * qJD(5);
t198 = t278 * qJDD(1) - t288;
t207 = qJD(5) * mrSges(6,1) - t216 * mrSges(6,3);
t197 = -t215 * pkin(5) - t216 * pkin(8);
t259 = qJD(5) ^ 2;
t159 = -t259 * pkin(5) + qJDD(5) * pkin(8) + t215 * t197 + t163;
t234 = t248 ^ 2;
t184 = -t249 * t205 + t251 * t210;
t275 = qJDD(1) * pkin(3) + qJDD(4) - t184;
t168 = pkin(4) * t286 + (-qJ(4) + (-t234 - t235) * pkin(7)) * t260 + t275;
t289 = t215 * qJD(5);
t199 = t279 * qJDD(1) + t289;
t160 = (-t199 - t289) * pkin(8) + (-t198 + t288) * pkin(5) + t168;
t252 = sin(qJ(6));
t255 = cos(qJ(6));
t156 = -t252 * t159 + t255 * t160;
t202 = t255 * qJD(5) - t252 * t216;
t177 = t202 * qJD(6) + t252 * qJDD(5) + t255 * t199;
t203 = t252 * qJD(5) + t255 * t216;
t182 = -t202 * mrSges(7,1) + t203 * mrSges(7,2);
t213 = qJD(6) - t215;
t186 = -t213 * mrSges(7,2) + t202 * mrSges(7,3);
t196 = qJDD(6) - t198;
t152 = m(7) * t156 + t196 * mrSges(7,1) - t177 * mrSges(7,3) - t203 * t182 + t213 * t186;
t157 = t255 * t159 + t252 * t160;
t176 = -t203 * qJD(6) + t255 * qJDD(5) - t252 * t199;
t187 = t213 * mrSges(7,1) - t203 * mrSges(7,3);
t153 = m(7) * t157 - t196 * mrSges(7,2) + t176 * mrSges(7,3) + t202 * t182 - t213 * t187;
t284 = -t252 * t152 + t255 * t153;
t139 = m(6) * t163 - qJDD(5) * mrSges(6,2) + t198 * mrSges(6,3) - qJD(5) * t207 + t215 * t192 + t284;
t162 = t256 * t166 - t253 * t167;
t206 = -qJD(5) * mrSges(6,2) + t215 * mrSges(6,3);
t158 = -qJDD(5) * pkin(5) - t259 * pkin(8) + t216 * t197 - t162;
t272 = -m(7) * t158 + t176 * mrSges(7,1) - t177 * mrSges(7,2) + t202 * t186 - t203 * t187;
t148 = m(6) * t162 + qJDD(5) * mrSges(6,1) - t199 * mrSges(6,3) + qJD(5) * t206 - t216 * t192 + t272;
t134 = t253 * t139 + t256 * t148;
t169 = -t248 * t180 + t291;
t171 = Ifges(7,5) * t203 + Ifges(7,6) * t202 + Ifges(7,3) * t213;
t173 = Ifges(7,1) * t203 + Ifges(7,4) * t202 + Ifges(7,5) * t213;
t145 = -mrSges(7,1) * t158 + mrSges(7,3) * t157 + Ifges(7,4) * t177 + Ifges(7,2) * t176 + Ifges(7,6) * t196 - t203 * t171 + t213 * t173;
t172 = Ifges(7,4) * t203 + Ifges(7,2) * t202 + Ifges(7,6) * t213;
t146 = mrSges(7,2) * t158 - mrSges(7,3) * t156 + Ifges(7,1) * t177 + Ifges(7,4) * t176 + Ifges(7,5) * t196 + t202 * t171 - t213 * t172;
t189 = Ifges(6,4) * t216 + Ifges(6,2) * t215 + Ifges(6,6) * qJD(5);
t190 = Ifges(6,1) * t216 + Ifges(6,4) * t215 + Ifges(6,5) * qJD(5);
t267 = -mrSges(6,1) * t162 + mrSges(6,2) * t163 - Ifges(6,5) * t199 - Ifges(6,6) * t198 - Ifges(6,3) * qJDD(5) - pkin(5) * t272 - pkin(8) * t284 - t255 * t145 - t252 * t146 - t216 * t189 + t215 * t190;
t282 = -Ifges(5,4) * t248 - Ifges(5,2) * t250;
t283 = -Ifges(5,1) * t248 - Ifges(5,4) * t250;
t297 = -mrSges(5,1) * t169 + mrSges(5,2) * t170 - pkin(4) * t134 + (t248 * t282 - t250 * t283) * t260 + t267;
t295 = mrSges(2,1) + mrSges(3,1);
t294 = mrSges(5,1) * t250;
t293 = mrSges(5,2) * t248;
t141 = t255 * t152 + t252 * t153;
t281 = -Ifges(5,5) * t248 - Ifges(5,6) * t250;
t290 = t260 * t281;
t277 = mrSges(5,3) * qJDD(1) + t260 * (-t293 + t294);
t132 = m(5) * t169 + t277 * t248 + t134;
t285 = t256 * t139 - t253 * t148;
t133 = m(5) * t170 - t277 * t250 + t285;
t128 = -t248 * t132 + t250 * t133;
t124 = m(4) * t185 - t260 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t128;
t179 = -t260 * qJ(4) + t275;
t269 = m(6) * t168 - t198 * mrSges(6,1) + t199 * mrSges(6,2) - t215 * t206 + t216 * t207 + t141;
t265 = -m(5) * t179 + qJDD(1) * t293 - t269 + (t234 * t260 + t292) * mrSges(5,3);
t135 = t265 - t260 * mrSges(4,2) + m(4) * t184 + (-mrSges(4,1) - t294) * qJDD(1);
t122 = t251 * t124 - t249 * t135;
t121 = t249 * t124 + t251 * t135;
t127 = t250 * t132 + t248 * t133;
t211 = -t260 * pkin(1) + t274;
t276 = m(3) * t211 + qJDD(1) * mrSges(3,3) + t122;
t126 = -m(4) * t244 - t127;
t214 = -qJDD(1) * pkin(1) + t273;
t271 = -m(3) * t214 + qJDD(1) * mrSges(3,1) + t260 * mrSges(3,3) - t121;
t188 = Ifges(6,5) * t216 + Ifges(6,6) * t215 + Ifges(6,3) * qJD(5);
t129 = mrSges(6,2) * t168 - mrSges(6,3) * t162 + Ifges(6,1) * t199 + Ifges(6,4) * t198 + Ifges(6,5) * qJDD(5) - pkin(8) * t141 - qJD(5) * t189 - t252 * t145 + t255 * t146 + t215 * t188;
t264 = mrSges(7,1) * t156 - mrSges(7,2) * t157 + Ifges(7,5) * t177 + Ifges(7,6) * t176 + Ifges(7,3) * t196 + t203 * t172 - t202 * t173;
t130 = -mrSges(6,1) * t168 + mrSges(6,3) * t163 + Ifges(6,4) * t199 + Ifges(6,2) * t198 + Ifges(6,6) * qJDD(5) - pkin(5) * t141 + qJD(5) * t190 - t216 * t188 - t264;
t115 = -mrSges(5,1) * t179 + mrSges(5,3) * t170 - pkin(4) * t269 + pkin(7) * t285 + t282 * qJDD(1) + t253 * t129 + t256 * t130 + t248 * t290;
t116 = mrSges(5,2) * t179 - mrSges(5,3) * t169 - pkin(7) * t134 + t283 * qJDD(1) + t256 * t129 - t253 * t130 - t250 * t290;
t113 = mrSges(4,2) * t244 - mrSges(4,3) * t184 - Ifges(4,5) * qJDD(1) - t260 * Ifges(4,6) - qJ(4) * t127 - t248 * t115 + t250 * t116;
t114 = (-Ifges(4,6) - t281) * qJDD(1) + t260 * Ifges(4,5) - mrSges(4,1) * t244 + mrSges(4,3) * t185 - pkin(3) * t127 + t297;
t270 = mrSges(3,2) * t214 + mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) + t260 * Ifges(3,6) - qJ(3) * t121 + t251 * t113 - t249 * t114;
t268 = mrSges(3,2) * t211 - pkin(2) * t126 - qJ(3) * t122 - t249 * t113 - t251 * t114;
t266 = mrSges(4,1) * t184 + pkin(3) * (-mrSges(5,1) * t286 + t265) + qJ(4) * t128 + t250 * t115 + t248 * t116 - mrSges(4,2) * t185 - Ifges(4,3) * qJDD(1);
t263 = -mrSges(3,1) * t214 + mrSges(3,3) * t211 + Ifges(3,2) * qJDD(1) - pkin(2) * t121 - t266;
t261 = -mrSges(2,2) * t223 + qJ(2) * (-t260 * mrSges(3,1) + t276) + pkin(1) * t271 + mrSges(2,1) * t222 + Ifges(2,3) * qJDD(1) + t263;
t125 = -m(3) * g(3) + t126;
t118 = m(2) * t222 + qJDD(1) * mrSges(2,1) - t260 * mrSges(2,2) + t271;
t117 = m(2) * t223 - qJDD(1) * mrSges(2,2) - t295 * t260 + t276;
t111 = -mrSges(2,2) * g(3) - mrSges(2,3) * t222 + Ifges(2,5) * qJDD(1) - t260 * Ifges(2,6) - qJ(2) * t125 + t270;
t110 = mrSges(2,3) * t223 - pkin(1) * t125 + (Ifges(3,4) + Ifges(2,5)) * t260 + (Ifges(2,6) - Ifges(3,6)) * qJDD(1) + t295 * g(3) + t268;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t257 * t111 - t254 * t110 - pkin(6) * (t254 * t117 + t257 * t118), t111, t270, t113, t116, t129, t146; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t254 * t111 + t257 * t110 + pkin(6) * (t257 * t117 - t254 * t118), t110, t263, t114, t115, t130, t145; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t261, t261, -mrSges(3,1) * g(3) - t260 * Ifges(3,4) + Ifges(3,6) * qJDD(1) - t268, t266, t281 * qJDD(1) - t297, -t267, t264;];
m_new  = t1;
