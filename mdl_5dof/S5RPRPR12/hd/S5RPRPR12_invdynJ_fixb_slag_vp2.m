% Calculate vector of inverse dynamics joint torques for
% S5RPRPR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR12_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR12_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:27
% EndTime: 2019-12-31 18:29:48
% DurationCPUTime: 12.72s
% Computational Cost: add. (5476->510), mult. (13378->681), div. (0->0), fcn. (10060->14), ass. (0->222)
t178 = sin(pkin(9));
t180 = cos(pkin(9));
t184 = sin(qJ(5));
t187 = cos(qJ(5));
t150 = t178 * t187 + t180 * t184;
t144 = t150 * qJD(5);
t179 = sin(pkin(8));
t185 = sin(qJ(3));
t181 = cos(pkin(8));
t261 = cos(qJ(3));
t221 = t261 * t181;
t197 = -t185 * t179 + t221;
t141 = t197 * qJD(1);
t90 = t150 * t141;
t299 = t90 - t144;
t198 = t178 * t184 - t180 * t187;
t143 = t198 * qJD(5);
t91 = t198 * t141;
t298 = t91 - t143;
t139 = qJD(5) - t141;
t151 = t179 * t261 + t185 * t181;
t142 = t151 * qJD(1);
t122 = qJD(3) * t178 + t142 * t180;
t212 = t180 * qJD(3) - t142 * t178;
t71 = t122 * t187 + t184 * t212;
t259 = Ifges(6,4) * t71;
t312 = -t122 * t184 + t187 * t212;
t27 = Ifges(6,2) * t312 + Ifges(6,6) * t139 + t259;
t282 = t27 / 0.2e1;
t68 = Ifges(6,4) * t312;
t28 = Ifges(6,1) * t71 + Ifges(6,5) * t139 + t68;
t281 = t28 / 0.2e1;
t230 = qJD(1) * qJD(2);
t162 = qJDD(1) * qJ(2) + t230;
t232 = t179 ^ 2 + t181 ^ 2;
t255 = pkin(6) + qJ(2);
t159 = t255 * t181;
t153 = qJD(1) * t159;
t140 = t185 * t153;
t157 = t255 * t179;
t152 = qJD(1) * t157;
t113 = -t152 * t261 - t140;
t104 = -qJD(3) * pkin(3) + qJD(4) - t113;
t167 = pkin(2) * t181 + pkin(1);
t155 = -qJD(1) * t167 + qJD(2);
t251 = Ifges(5,4) * t180;
t202 = -Ifges(5,2) * t178 + t251;
t252 = Ifges(5,4) * t178;
t203 = Ifges(5,1) * t180 - t252;
t204 = mrSges(5,1) * t178 + mrSges(5,2) * t180;
t262 = t180 / 0.2e1;
t297 = Ifges(4,5) * qJD(3);
t305 = t212 / 0.2e1;
t306 = t122 / 0.2e1;
t318 = t155 * mrSges(4,2) + t202 * t305 + t203 * t306 + (t122 * Ifges(5,1) + Ifges(5,4) * t212 - Ifges(5,5) * t141) * t262 - t178 * (t122 * Ifges(5,4) + Ifges(5,2) * t212 - t141 * Ifges(5,6)) / 0.2e1 - t113 * mrSges(4,3) + t104 * t204 + t297 / 0.2e1;
t114 = -t185 * t152 + t153 * t261;
t106 = qJD(3) * qJ(4) + t114;
t81 = -pkin(3) * t141 - qJ(4) * t142 + t155;
t45 = -t106 * t178 + t180 * t81;
t29 = -pkin(4) * t141 - pkin(7) * t122 + t45;
t46 = t180 * t106 + t178 * t81;
t35 = pkin(7) * t212 + t46;
t10 = t184 * t29 + t187 * t35;
t296 = Ifges(4,6) * qJD(3);
t9 = -t184 * t35 + t187 * t29;
t317 = -t155 * mrSges(4,1) - t45 * mrSges(5,1) - t9 * mrSges(6,1) + t46 * mrSges(5,2) + t10 * mrSges(6,2) + t114 * mrSges(4,3) + t296 / 0.2e1;
t301 = t212 * Ifges(5,6);
t302 = t122 * Ifges(5,5);
t315 = t71 * Ifges(6,5) + Ifges(6,6) * t312 - t141 * Ifges(5,3) + t139 * Ifges(6,3) + t301 + t302;
t300 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t212 - mrSges(5,2) * t122 - mrSges(4,3) * t142;
t186 = sin(qJ(1));
t188 = cos(qJ(1));
t314 = g(1) * t188 + g(2) * t186;
t218 = m(3) * qJ(2) + mrSges(3,3);
t311 = -m(5) * t255 + mrSges(2,2) - mrSges(4,3) - t204 - t218;
t177 = pkin(8) + qJ(3);
t171 = sin(t177);
t173 = cos(t177);
t207 = mrSges(4,1) * t173 - mrSges(4,2) * t171;
t208 = -mrSges(3,1) * t181 + mrSges(3,2) * t179;
t310 = m(3) * pkin(1) + t171 * mrSges(6,3) + mrSges(2,1) + t207 - t208;
t145 = t197 * qJD(3);
t111 = qJD(1) * t145 + qJDD(1) * t151;
t146 = t151 * qJD(3);
t228 = qJDD(1) * t179;
t112 = qJD(1) * t146 - qJDD(1) * t221 + t185 * t228;
t154 = -qJDD(1) * t167 + qJDD(2);
t42 = pkin(3) * t112 - qJ(4) * t111 - qJD(4) * t142 + t154;
t213 = pkin(6) * qJDD(1) + t162;
t133 = t213 * t179;
t134 = t213 * t181;
t215 = qJD(3) * t261;
t222 = -t185 * t133 + t261 * t134 - t152 * t215;
t51 = qJDD(3) * qJ(4) + (qJD(4) - t140) * qJD(3) + t222;
t19 = t178 * t42 + t180 * t51;
t92 = qJDD(3) * t180 - t111 * t178;
t11 = pkin(7) * t92 + t19;
t18 = -t178 * t51 + t180 * t42;
t93 = qJDD(3) * t178 + t111 * t180;
t8 = pkin(4) * t112 - pkin(7) * t93 + t18;
t1 = qJD(5) * t9 + t11 * t187 + t184 * t8;
t2 = -qJD(5) * t10 - t11 * t184 + t187 * t8;
t309 = t2 * mrSges(6,1) - t1 * mrSges(6,2);
t23 = qJD(5) * t312 + t184 * t92 + t187 * t93;
t284 = t23 / 0.2e1;
t24 = -qJD(5) * t71 - t184 * t93 + t187 * t92;
t283 = t24 / 0.2e1;
t275 = t92 / 0.2e1;
t274 = t93 / 0.2e1;
t109 = qJDD(5) + t112;
t272 = t109 / 0.2e1;
t271 = t112 / 0.2e1;
t254 = pkin(7) + qJ(4);
t156 = t254 * t178;
t158 = t254 * t180;
t118 = -t156 * t184 + t158 * t187;
t241 = t141 * t180;
t107 = pkin(3) * t142 - qJ(4) * t141;
t59 = t180 * t107 - t113 * t178;
t36 = pkin(4) * t142 - pkin(7) * t241 + t59;
t242 = t141 * t178;
t60 = t178 * t107 + t180 * t113;
t43 = -pkin(7) * t242 + t60;
t304 = -qJD(4) * t150 - qJD(5) * t118 + t184 * t43 - t187 * t36;
t116 = -t156 * t187 - t158 * t184;
t303 = -qJD(4) * t198 + qJD(5) * t116 - t184 * t36 - t187 * t43;
t294 = t171 * t314;
t205 = -t180 * mrSges(5,1) + t178 * mrSges(5,2);
t195 = m(5) * pkin(3) - t205;
t291 = t195 * t173;
t290 = -t261 * t157 - t185 * t159;
t87 = mrSges(5,2) * t141 + mrSges(5,3) * t212;
t88 = -mrSges(5,1) * t141 - mrSges(5,3) * t122;
t289 = -t178 * t88 + t180 * t87;
t288 = -t178 * t18 + t180 * t19;
t286 = Ifges(6,4) * t284 + Ifges(6,2) * t283 + Ifges(6,6) * t272;
t285 = Ifges(6,1) * t284 + Ifges(6,4) * t283 + Ifges(6,5) * t272;
t280 = Ifges(5,1) * t274 + Ifges(5,4) * t275 + Ifges(5,5) * t271;
t279 = -t312 / 0.2e1;
t278 = t312 / 0.2e1;
t277 = -t71 / 0.2e1;
t276 = t71 / 0.2e1;
t273 = m(5) + m(6);
t270 = -t139 / 0.2e1;
t269 = t139 / 0.2e1;
t268 = t141 / 0.2e1;
t267 = -t141 / 0.2e1;
t265 = t142 / 0.2e1;
t258 = pkin(4) * t178;
t77 = pkin(3) * t146 - qJ(4) * t145 - qJD(4) * t151;
t85 = t197 * qJD(2) + qJD(3) * t290;
t40 = t178 * t77 + t180 * t85;
t253 = Ifges(4,4) * t142;
t110 = -pkin(3) * t197 - qJ(4) * t151 - t167;
t119 = -t185 * t157 + t159 * t261;
t65 = t178 * t110 + t180 * t119;
t240 = t145 * t178;
t239 = t145 * t180;
t238 = t151 * t178;
t237 = t151 * t180;
t236 = t173 * t186;
t235 = t173 * t188;
t231 = qJD(3) * t185;
t227 = qJDD(1) * t181;
t226 = Ifges(6,5) * t23 + Ifges(6,6) * t24 + Ifges(6,3) * t109;
t217 = m(5) * qJ(4) + mrSges(5,3);
t216 = m(6) * t254 + mrSges(6,3);
t44 = -t92 * mrSges(5,1) + t93 * mrSges(5,2);
t7 = -t24 * mrSges(6,1) + t23 * mrSges(6,2);
t39 = -t178 * t85 + t180 * t77;
t214 = t112 * mrSges(4,1) + t111 * mrSges(4,2);
t64 = t180 * t110 - t119 * t178;
t209 = -mrSges(3,1) * t227 + mrSges(3,2) * t228;
t201 = Ifges(5,5) * t180 - Ifges(5,6) * t178;
t200 = t178 * t45 - t180 * t46;
t41 = -pkin(4) * t197 - pkin(7) * t237 + t64;
t47 = -pkin(7) * t238 + t65;
t14 = -t184 * t47 + t187 * t41;
t15 = t184 * t41 + t187 * t47;
t166 = pkin(4) * t180 + pkin(3);
t199 = t166 * t173 + t171 * t254;
t58 = -t133 * t261 - t185 * t134 + t152 * t231 - t153 * t215;
t176 = pkin(9) + qJ(5);
t170 = sin(t176);
t172 = cos(t176);
t193 = m(6) * t166 + mrSges(6,1) * t172 - mrSges(6,2) * t170;
t52 = -qJDD(3) * pkin(3) + qJDD(4) - t58;
t192 = t171 * t217 + t291;
t86 = qJD(2) * t151 + qJD(3) * t119;
t168 = -qJDD(1) * pkin(1) + qJDD(2);
t160 = t188 * t167;
t135 = Ifges(4,4) * t141;
t128 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t141;
t127 = t170 * t186 + t172 * t235;
t126 = -t170 * t235 + t172 * t186;
t125 = t170 * t188 - t172 * t236;
t124 = t170 * t236 + t172 * t188;
t102 = t142 * Ifges(4,1) + t135 + t297;
t101 = t141 * Ifges(4,2) + t253 + t296;
t98 = t198 * t151;
t97 = t150 * t151;
t89 = pkin(4) * t238 - t290;
t76 = pkin(4) * t242 + t114;
t67 = -pkin(4) * t212 + t104;
t66 = pkin(4) * t240 + t86;
t57 = -t153 * t231 + t222;
t56 = mrSges(6,1) * t139 - mrSges(6,3) * t71;
t55 = -mrSges(6,2) * t139 + mrSges(6,3) * t312;
t54 = mrSges(5,1) * t112 - mrSges(5,3) * t93;
t53 = -mrSges(5,2) * t112 + mrSges(5,3) * t92;
t50 = t143 * t151 - t145 * t150;
t49 = -t144 * t151 - t145 * t198;
t34 = -mrSges(6,1) * t312 + mrSges(6,2) * t71;
t33 = -pkin(7) * t240 + t40;
t31 = t93 * Ifges(5,4) + t92 * Ifges(5,2) + t112 * Ifges(5,6);
t30 = -t92 * pkin(4) + t52;
t25 = pkin(4) * t146 - pkin(7) * t239 + t39;
t17 = -mrSges(6,2) * t109 + mrSges(6,3) * t24;
t16 = mrSges(6,1) * t109 - mrSges(6,3) * t23;
t4 = -qJD(5) * t15 - t184 * t33 + t187 * t25;
t3 = qJD(5) * t14 + t184 * t25 + t187 * t33;
t5 = [-(Ifges(5,5) * t93 + Ifges(5,6) * t92 + Ifges(5,3) * t112 + t226) * t197 / 0.2e1 + m(5) * (t18 * t64 + t19 * t65 + t39 * t45 + t40 * t46) + m(4) * (t114 * t85 + t119 * t57 - t154 * t167) + (-Ifges(6,1) * t98 - Ifges(6,4) * t97) * t284 + m(6) * (t1 * t15 + t10 * t3 + t14 * t2 + t30 * t89 + t4 * t9 + t66 * t67) + (Ifges(6,1) * t49 + Ifges(6,4) * t50) * t276 + (-t125 * mrSges(6,1) - t124 * mrSges(6,2) + (-m(4) * t255 - m(6) * (t255 + t258) + t311) * t188 + (m(4) * t167 - m(6) * (-t167 - t199) - m(5) * (-qJ(4) * t171 - t167) + t171 * mrSges(5,3) + t291 + t310) * t186) * g(1) + (-m(5) * t160 - t127 * mrSges(6,1) - t126 * mrSges(6,2) + (-m(4) - m(6)) * (t186 * t255 + t160) + (-m(6) * t258 + t311) * t186 + (-m(6) * t199 - t192 - t310) * t188) * g(2) + (-t154 * mrSges(4,1) - t18 * mrSges(5,1) + t19 * mrSges(5,2) + t57 * mrSges(4,3) + Ifges(4,4) * t111 - Ifges(5,5) * t274 - Ifges(6,5) * t284 - Ifges(4,2) * t112 + Ifges(4,6) * qJDD(3) - Ifges(5,6) * t275 - Ifges(6,6) * t283 - Ifges(5,3) * t271 - Ifges(6,3) * t272 - t309) * t197 + (-t1 * t97 + t10 * t50 + t2 * t98 - t49 * t9) * mrSges(6,3) + (Ifges(6,5) * t49 + Ifges(6,6) * t50) * t269 + (Ifges(3,4) * t179 + Ifges(3,2) * t181) * t227 + (Ifges(3,1) * t179 + Ifges(3,4) * t181) * t228 + (-m(4) * t113 + m(5) * t104 - t300) * t86 + t30 * (mrSges(6,1) * t97 - mrSges(6,2) * t98) - t167 * t214 - pkin(1) * t209 + t237 * t280 + t49 * t281 + t50 * t282 - t98 * t285 - t97 * t286 + (t154 * mrSges(4,2) - t58 * mrSges(4,3) + Ifges(4,1) * t111 - Ifges(4,4) * t112 + Ifges(4,5) * qJDD(3) + t201 * t271 + t202 * t275 + t203 * t274 + t204 * t52) * t151 + (-Ifges(6,5) * t98 - Ifges(6,6) * t97) * t272 - (-m(4) * t58 + m(5) * t52 - qJDD(3) * mrSges(4,1) + mrSges(4,3) * t111 + t44) * t290 + (-t18 * t237 - t19 * t238 - t239 * t45 - t240 * t46) * mrSges(5,3) + 0.2e1 * t232 * t162 * mrSges(3,3) + t168 * t208 + (-Ifges(6,4) * t98 - Ifges(6,2) * t97) * t283 + (Ifges(6,4) * t49 + Ifges(6,2) * t50) * t278 - t31 * t238 / 0.2e1 + m(3) * (-pkin(1) * t168 + (t162 + t230) * qJ(2) * t232) + t14 * t16 + t15 * t17 + t3 * t55 + t4 * t56 + t64 * t54 + t65 * t53 + t66 * t34 + t67 * (-mrSges(6,1) * t50 + mrSges(6,2) * t49) + t40 * t87 + t39 * t88 + t89 * t7 + Ifges(2,3) * qJDD(1) + t119 * (-qJDD(3) * mrSges(4,2) - mrSges(4,3) * t112) + t85 * t128 + (t315 / 0.2e1 - Ifges(4,4) * t265 + Ifges(5,3) * t267 - Ifges(4,2) * t268 + Ifges(6,3) * t269 + Ifges(6,5) * t276 + Ifges(6,6) * t278 + Ifges(5,6) * t305 + Ifges(5,5) * t306 - t101 / 0.2e1 - t317) * t146 + (Ifges(4,1) * t265 + t201 * t267 + Ifges(4,4) * t268 + t102 / 0.2e1 + t318) * t145; (-t34 + t300) * t142 - (t128 + t289) * t141 + m(3) * t168 + t214 + t299 * t56 + t298 * t55 + t209 - t198 * t16 + t150 * t17 + t178 * t53 + t180 * t54 + (-g(1) * t186 + g(2) * t188) * (m(3) + m(4) + t273) - t218 * t232 * qJD(1) ^ 2 + (t1 * t150 + t10 * t298 - t142 * t67 - t198 * t2 + t299 * t9) * m(6) + (-t104 * t142 + t141 * t200 + t19 * t178 + t18 * t180) * m(5) + (t113 * t142 - t114 * t141 + t154) * m(4); (Ifges(6,5) * t150 - Ifges(6,6) * t198) * t272 + (Ifges(6,4) * t150 - Ifges(6,2) * t198) * t283 + (Ifges(6,1) * t150 - Ifges(6,4) * t198) * t284 + t30 * (mrSges(6,1) * t198 + mrSges(6,2) * t150) - t198 * t286 + (-t1 * t198 + t10 * t299 - t150 * t2 - t298 * t9) * mrSges(6,3) + t298 * t281 + t299 * t282 + (t314 * (mrSges(4,2) - t216 - t217) - t193 * g(3)) * t173 + t303 * t55 + (t1 * t118 + t10 * t303 + t116 * t2 - t166 * t30 + t304 * t9 - t67 * t76) * m(6) + t304 * t56 + (t314 * (mrSges(4,1) + t193 + t195) - t216 * g(3)) * t171 + t31 * t262 + t101 * t265 + (Ifges(5,5) * t178 + Ifges(5,6) * t180) * t271 + (Ifges(5,1) * t178 + t251) * t274 + (Ifges(5,2) * t180 + t252) * t275 + (-Ifges(6,5) * t143 - Ifges(6,6) * t144) * t269 + (-Ifges(6,1) * t143 - Ifges(6,4) * t144) * t276 + (-Ifges(6,4) * t143 - Ifges(6,2) * t144) * t278 + (-Ifges(6,1) * t91 - Ifges(6,4) * t90) * t277 + (-Ifges(6,4) * t91 - Ifges(6,2) * t90) * t279 + (-Ifges(6,5) * t91 - Ifges(6,6) * t90) * t270 + (t102 + t135) * t267 + t178 * t280 + t150 * t285 - (Ifges(4,1) * t141 - t253 + t315) * t142 / 0.2e1 + (t241 * t45 + t242 * t46 + t288) * mrSges(5,3) + (-pkin(3) * t52 + qJ(4) * t288 - t200 * qJD(4) - t104 * t114 - t45 * t59 - t46 * t60) * m(5) + t289 * qJD(4) + t52 * t205 + (-t192 - t207) * g(3) + (-t178 * t54 + t180 * t53) * qJ(4) + (-mrSges(6,1) * t299 + mrSges(6,2) * t298) * t67 + t300 * t114 + Ifges(4,3) * qJDD(3) - pkin(3) * t44 - t57 * mrSges(4,2) + t58 * mrSges(4,1) - t76 * t34 - t60 * t87 - t59 * t88 + Ifges(4,5) * t111 - Ifges(4,6) * t112 + t116 * t16 + t118 * t17 - t113 * t128 + (-Ifges(4,2) * t267 + Ifges(5,3) * t268 + Ifges(6,3) * t270 + Ifges(6,5) * t277 + Ifges(6,6) * t279 - t301 / 0.2e1 - t302 / 0.2e1 + t317) * t142 - (-t201 * t268 + t318) * t141 - t166 * t7; t273 * t173 * g(3) + t122 * t88 - t212 * t87 - t312 * t55 + t71 * t56 + t44 + t7 + (-t10 * t312 + t71 * t9 - t294 + t30) * m(6) + (t122 * t45 - t212 * t46 - t294 + t52) * m(5); -t67 * (mrSges(6,1) * t71 + mrSges(6,2) * t312) + (Ifges(6,1) * t312 - t259) * t277 + t27 * t276 + (Ifges(6,5) * t312 - Ifges(6,6) * t71) * t270 - t9 * t55 + t10 * t56 - g(1) * (mrSges(6,1) * t126 - mrSges(6,2) * t127) - g(2) * (-mrSges(6,1) * t124 + mrSges(6,2) * t125) - g(3) * (-mrSges(6,1) * t170 - mrSges(6,2) * t172) * t171 + (t10 * t71 + t312 * t9) * mrSges(6,3) + t226 + (-Ifges(6,2) * t71 + t28 + t68) * t279 + t309;];
tau = t5;
