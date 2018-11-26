% Calculate vector of centrifugal and coriolis load on the joints for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:06
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6PRPRRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:06:09
% EndTime: 2018-11-23 15:06:13
% DurationCPUTime: 4.72s
% Computational Cost: add. (5495->429), mult. (12238->595), div. (0->0), fcn. (8410->10), ass. (0->217)
t169 = sin(qJ(5));
t173 = cos(qJ(5));
t174 = cos(qJ(4));
t228 = qJD(4) * t174;
t238 = t169 * t174;
t170 = sin(qJ(4));
t229 = qJD(4) * t170;
t294 = -qJD(5) * t170 - t229;
t102 = -qJD(5) * t238 - t169 * t228 + t173 * t294;
t165 = qJD(4) + qJD(5);
t237 = t173 * t174;
t209 = t165 * t237;
t103 = t169 * t294 + t209;
t163 = pkin(4) * t228;
t149 = qJD(3) + t163;
t175 = cos(qJ(2));
t166 = sin(pkin(6));
t235 = qJD(1) * t166;
t217 = t175 * t235;
t311 = pkin(5) * t103 - pkin(10) * t102 + t149 - t217;
t176 = -pkin(2) - pkin(8);
t270 = pkin(9) - t176;
t143 = t270 * t170;
t144 = t270 * t174;
t106 = -t143 * t169 + t173 * t144;
t134 = qJD(4) * t144;
t139 = t170 * t173 + t238;
t210 = t270 * t229;
t171 = sin(qJ(2));
t221 = t171 * t235;
t301 = -qJD(5) * t106 - t173 * t134 - t139 * t221 + t169 * t210;
t310 = t165 * Ifges(6,6) / 0.2e1;
t168 = sin(qJ(6));
t172 = cos(qJ(6));
t199 = qJD(3) - t217;
t127 = qJD(2) * t176 + t199;
t167 = cos(pkin(6));
t234 = qJD(1) * t174;
t215 = t167 * t234;
t291 = t170 * (pkin(9) * qJD(2) - t127) - t215;
t251 = t173 * t291;
t120 = t174 * t127;
t240 = t167 * t170;
t220 = qJD(1) * t240;
t231 = qJD(2) * t174;
t185 = -pkin(9) * t231 - t220;
t98 = t120 + t185;
t91 = qJD(4) * pkin(4) + t98;
t42 = t169 * t91 - t251;
t36 = pkin(10) * t165 + t42;
t141 = qJD(2) * qJ(3) + t221;
t233 = qJD(2) * t170;
t121 = pkin(4) * t233 + t141;
t130 = t139 * qJD(2);
t216 = t169 * t233;
t131 = t173 * t231 - t216;
t67 = pkin(5) * t130 - pkin(10) * t131 + t121;
t17 = -t168 * t36 + t172 * t67;
t135 = (qJD(3) + t217) * qJD(2);
t118 = qJD(2) * t163 + t135;
t92 = t165 * t130;
t93 = qJD(2) * t209 - t165 * t216;
t29 = pkin(5) * t93 + pkin(10) * t92 + t118;
t232 = qJD(2) * t171;
t219 = t166 * t232;
t212 = qJD(1) * t219;
t295 = t291 * qJD(4) + t174 * t212;
t253 = t169 * t291;
t41 = t173 * t91 + t253;
t236 = t127 * t228 + t170 * t212;
t68 = qJD(4) * t185 + t236;
t8 = qJD(5) * t41 + t169 * t295 + t173 * t68;
t2 = qJD(6) * t17 + t168 * t29 + t172 * t8;
t18 = t168 * t67 + t172 * t36;
t3 = -qJD(6) * t18 - t168 * t8 + t172 * t29;
t309 = -t168 * t3 + t172 * t2;
t256 = t165 * Ifges(6,5);
t308 = t121 * mrSges(6,2) + t256 / 0.2e1;
t257 = t131 * Ifges(6,4);
t307 = t310 + t257 / 0.2e1 - t130 * Ifges(6,2) / 0.2e1;
t107 = -t143 * t173 - t144 * t169;
t138 = t169 * t170 - t237;
t158 = t170 * pkin(4) + qJ(3);
t97 = pkin(5) * t139 + pkin(10) * t138 + t158;
t51 = -t107 * t168 + t172 * t97;
t305 = qJD(6) * t51 + t311 * t168 + t301 * t172;
t52 = t107 * t172 + t168 * t97;
t304 = -qJD(6) * t52 - t301 * t168 + t311 * t172;
t113 = -t131 * t168 + t165 * t172;
t55 = qJD(6) * t113 - t172 * t92;
t114 = t131 * t172 + t165 * t168;
t56 = -qJD(6) * t114 + t168 * t92;
t16 = -mrSges(7,1) * t56 + mrSges(7,2) * t55;
t9 = t42 * qJD(5) + t169 * t68 - t173 * t295;
t302 = m(7) * t9 + t16;
t300 = -qJD(5) * t107 + t134 * t169 + t138 * t221 + t173 * t210;
t197 = t168 * t18 + t17 * t172;
t299 = -qJD(6) * t197 + t309;
t206 = mrSges(7,1) * t168 + mrSges(7,2) * t172;
t35 = -pkin(5) * t165 - t41;
t187 = t35 * t206;
t201 = Ifges(7,5) * t172 - Ifges(7,6) * t168;
t263 = Ifges(7,4) * t172;
t203 = -Ifges(7,2) * t168 + t263;
t264 = Ifges(7,4) * t168;
t205 = Ifges(7,1) * t172 - t264;
t277 = t172 / 0.2e1;
t123 = qJD(6) + t130;
t282 = t123 / 0.2e1;
t284 = t114 / 0.2e1;
t286 = t113 / 0.2e1;
t265 = Ifges(7,4) * t114;
t46 = Ifges(7,2) * t113 + Ifges(7,6) * t123 + t265;
t296 = -t46 / 0.2e1;
t108 = Ifges(7,4) * t113;
t47 = t114 * Ifges(7,1) + t123 * Ifges(7,5) + t108;
t298 = t168 * t296 + t201 * t282 + t203 * t286 + t205 * t284 + t47 * t277 + t187;
t297 = -t121 * mrSges(6,1) - t17 * mrSges(7,1) + t18 * mrSges(7,2) + t307;
t196 = -t168 * t17 + t172 * t18;
t293 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t55 + Ifges(7,6) * t56;
t292 = m(5) * t176 - mrSges(5,3);
t241 = t166 * t175;
t125 = -t174 * t241 - t240;
t188 = -t167 * t174 + t170 * t241;
t192 = t173 * t125 + t169 * t188;
t288 = t192 * t9;
t287 = -t113 / 0.2e1;
t285 = -t114 / 0.2e1;
t283 = -t123 / 0.2e1;
t281 = t130 / 0.2e1;
t280 = -t131 / 0.2e1;
t279 = t168 / 0.2e1;
t278 = t170 / 0.2e1;
t276 = -t174 / 0.2e1;
t275 = m(6) * t121;
t274 = t106 * t9;
t273 = t138 * t9;
t269 = mrSges(6,3) * t130;
t268 = mrSges(6,3) * t131;
t267 = Ifges(5,4) * t170;
t266 = Ifges(5,4) * t174;
t122 = Ifges(6,4) * t130;
t262 = t113 * Ifges(7,6);
t261 = t114 * Ifges(7,5);
t260 = t123 * Ifges(7,3);
t258 = t131 * Ifges(6,1);
t73 = -t127 * t229 + (-qJD(4) * t167 + t219) * t234;
t250 = t174 * t73;
t249 = mrSges(6,1) * t165 + mrSges(7,1) * t113 - mrSges(7,2) * t114 - t268;
t140 = (mrSges(5,1) * t170 + mrSges(5,2) * t174) * qJD(2);
t95 = mrSges(6,1) * t130 + mrSges(6,2) * t131;
t248 = t140 + t95;
t247 = Ifges(5,5) * qJD(4);
t246 = Ifges(5,6) * qJD(4);
t245 = t130 * t168;
t244 = t130 * t172;
t243 = t141 * t175;
t242 = t166 * t171;
t230 = qJD(2) * t175;
t226 = qJD(6) * t168;
t225 = qJD(6) * t172;
t224 = qJD(2) * qJD(4);
t218 = t166 * t230;
t214 = -t95 - t275;
t96 = pkin(5) * t131 + pkin(10) * t130;
t208 = mrSges(5,1) * t174 - mrSges(5,2) * t170;
t207 = -mrSges(7,1) * t172 + mrSges(7,2) * t168;
t204 = Ifges(7,1) * t168 + t263;
t202 = Ifges(7,2) * t172 + t264;
t200 = Ifges(7,5) * t168 + Ifges(7,6) * t172;
t198 = -t102 * t41 - t103 * t42;
t23 = mrSges(7,1) * t93 - mrSges(7,3) * t55;
t24 = -mrSges(7,2) * t93 + mrSges(7,3) * t56;
t195 = -t168 * t23 + t172 * t24;
t69 = -mrSges(7,2) * t123 + mrSges(7,3) * t113;
t70 = mrSges(7,1) * t123 - mrSges(7,3) * t114;
t194 = -t168 * t69 - t172 * t70;
t109 = t120 - t220;
t110 = t127 * t170 + t215;
t193 = -t109 * t170 + t110 * t174;
t84 = t125 * t169 - t173 * t188;
t190 = qJ(3) * t135 + qJD(3) * t141;
t115 = -mrSges(6,2) * t165 - t269;
t189 = -t168 * t70 + t172 * t69 + t115;
t65 = -t168 * t84 + t172 * t242;
t66 = t168 * t242 + t172 * t84;
t186 = t135 * t171 + t141 * t230;
t179 = m(7) * (-t17 * t225 - t18 * t226 + t309) - t69 * t226 - t70 * t225 + t195;
t12 = t55 * Ifges(7,4) + t56 * Ifges(7,2) + t93 * Ifges(7,6);
t13 = t55 * Ifges(7,1) + t56 * Ifges(7,4) + t93 * Ifges(7,5);
t45 = t260 + t261 + t262;
t81 = -t122 + t256 + t258;
t178 = t55 * t204 / 0.2e1 + t56 * t202 / 0.2e1 + t47 * t244 / 0.2e1 + t12 * t277 + t13 * t279 + t245 * t296 - t41 * t269 - t8 * mrSges(6,2) - Ifges(6,5) * t92 + (t200 / 0.2e1 - Ifges(6,6)) * t93 + (t207 - mrSges(6,1)) * t9 + (t81 - t122) * t281 + (-t257 + t45) * t280 + (-Ifges(6,1) * t280 - t201 * t283 - t203 * t287 - t205 * t285 + t187 + t308) * t130 + (Ifges(7,5) * t285 - Ifges(6,2) * t281 + Ifges(7,6) * t287 + Ifges(7,3) * t283 + t297 + t310) * t131 + (-t17 * t244 - t18 * t245 + t299) * mrSges(7,3) + t298 * qJD(6);
t177 = qJD(2) ^ 2;
t146 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t231;
t145 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t233;
t137 = -qJD(2) * pkin(2) + t199;
t132 = t208 * t224;
t129 = t247 + (t174 * Ifges(5,1) - t267) * qJD(2);
t128 = t246 + (-Ifges(5,2) * t170 + t266) * qJD(2);
t105 = qJD(4) * t188 + t174 * t219;
t104 = qJD(4) * t125 + t170 * t219;
t87 = Ifges(7,3) * t93;
t79 = pkin(4) * t231 + t96;
t72 = -qJD(4) * t220 + t236;
t49 = t173 * t98 + t253;
t48 = t169 * t98 - t251;
t44 = mrSges(6,1) * t93 - mrSges(6,2) * t92;
t26 = qJD(5) * t84 + t104 * t169 - t173 * t105;
t25 = qJD(5) * t192 + t104 * t173 + t105 * t169;
t22 = t168 * t96 + t172 * t41;
t21 = -t168 * t41 + t172 * t96;
t20 = t168 * t79 + t172 * t49;
t19 = -t168 * t49 + t172 * t79;
t15 = -qJD(6) * t66 - t168 * t25 + t172 * t218;
t14 = qJD(6) * t65 + t168 * t218 + t172 * t25;
t1 = [t104 * t145 + t105 * t146 + t25 * t115 + t14 * t69 + t15 * t70 - t192 * t16 + t65 * t23 + t66 * t24 - t249 * t26 + (t192 * t92 - t84 * t93) * mrSges(6,3) + (t125 * t170 + t174 * t188) * mrSges(5,3) * t224 + m(7) * (t14 * t18 + t15 * t17 + t2 * t66 + t26 * t35 + t3 * t65 - t288) + m(6) * (t25 * t42 - t26 * t41 + t8 * t84 - t288) + m(5) * (t104 * t110 + t105 * t109 + t125 * t73 - t188 * t72) + (((-mrSges(3,2) + mrSges(4,3)) * t177 + t248 * qJD(2)) * t175 + (t132 + t44 + (-mrSges(3,1) + mrSges(4,2)) * t177) * t171 + m(6) * (t118 * t171 + t121 * t230) + m(5) * t186 + (t137 * t232 - t175 * t212 + t186) * m(4)) * t166; t304 * t70 + (t17 * t304 + t18 * t305 + t2 * t52 + t3 * t51 - t300 * t35 + t274) * m(7) + t305 * t69 + t300 * t249 + (t107 * t8 + t118 * t158 + t121 * t149 + t300 * t41 + t301 * t42 + t274) * m(6) + t301 * t115 + (t135 * mrSges(5,1) + (-t141 * mrSges(5,2) - t129 / 0.2e1 - t176 * t146 - t247 / 0.2e1 - t292 * t109 + (0.3e1 / 0.2e1 * t267 + (-0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(5,2)) * t174) * qJD(2)) * qJD(4) + t292 * t72 + (-m(5) * t110 - t145) * t221) * t170 + (qJD(2) * t199 + t135) * mrSges(4,3) + (t81 / 0.2e1 - t122 / 0.2e1 + t258 / 0.2e1 - t197 * mrSges(7,3) + t298 + t308) * t102 + m(5) * (t190 + (t110 * t228 + t250) * t176) + m(4) * t190 + 0.2e1 * (-m(5) * (t109 * t171 * t174 + t243) / 0.2e1 - t175 * t275 / 0.2e1 - (pkin(2) * t232 + t137 * t171 + t243) * m(4) / 0.2e1) * t235 + (t262 / 0.2e1 + t261 / 0.2e1 + t260 / 0.2e1 + t45 / 0.2e1 - t297 - t307) * t103 + (-t8 * mrSges(6,3) + t87 / 0.2e1 + Ifges(6,4) * t92 + t118 * mrSges(6,1) + (Ifges(7,3) / 0.2e1 + Ifges(6,2)) * t93 + t293) * t139 + (-t106 * t92 - t107 * t93 + t198) * mrSges(6,3) + (-t55 * t205 / 0.2e1 - t56 * t203 / 0.2e1 + t12 * t279 - t172 * t13 / 0.2e1 + Ifges(6,1) * t92 - t118 * mrSges(6,2) + (-mrSges(6,3) - t206) * t9 + (t168 * t2 + t172 * t3) * mrSges(7,3) + (mrSges(7,3) * t196 + t200 * t282 + t202 * t286 + t204 * t284 + t207 * t35 + t277 * t46 + t279 * t47) * qJD(6) + (-t201 / 0.2e1 + Ifges(6,4)) * t93) * t138 + (-t146 * t221 + t135 * mrSges(5,2) - t73 * mrSges(5,3) + (t141 * mrSges(5,1) - t128 / 0.2e1 - t110 * mrSges(5,3) + t176 * t145 - 0.3e1 / 0.2e1 * Ifges(5,4) * t231 - t246 / 0.2e1) * qJD(4)) * t174 - t248 * t217 + t51 * t23 + t52 * t24 + t106 * t16 + qJ(3) * t132 + qJD(3) * t140 + t149 * t95 + t158 * t44; -t177 * mrSges(4,3) + (-t92 * mrSges(6,3) + t16) * t138 + t249 * t102 + (t145 * t174 - t146 * t170) * qJD(4) + t189 * t103 + (-t93 * mrSges(6,3) + qJD(6) * t194 + t195) * t139 + m(6) * (t139 * t8 - t198 + t273) + m(7) * (-t102 * t35 + t196 * t103 + t139 * t299 + t273) + m(5) * (qJD(4) * t193 + t170 * t72 + t250) + (-m(5) * t141 - m(7) * t197 - t140 + (-t141 + t221) * m(4) + t194 + t214) * qJD(2); t179 * (pkin(4) * t169 + pkin(10)) + (-t141 * t208 + t129 * t278 + t174 * t128 / 0.2e1 + ((-Ifges(5,2) * t174 - t267) * t278 + (-Ifges(5,1) * t170 - t266) * t276) * qJD(2) + t214 * t174 * pkin(4) + t193 * mrSges(5,3) + (-Ifges(5,5) * t170 / 0.2e1 + Ifges(5,6) * t276) * qJD(4)) * qJD(2) + (m(6) * (t169 * t8 - t173 * t9) + (-t169 * t93 + t173 * t92) * mrSges(6,3) + ((-m(6) * t41 + m(7) * t35 - t249) * t169 + (m(6) * t42 + m(7) * t196 + t189) * t173) * qJD(5)) * pkin(4) - m(6) * (-t41 * t48 + t42 * t49) + t178 - m(7) * (t17 * t19 + t18 * t20 + t35 * t48) + t249 * t48 + t42 * t268 - t20 * t69 - t19 * t70 - t72 * mrSges(5,2) + t73 * mrSges(5,1) - t49 * t115 - t109 * t145 + t110 * t146 + t302 * (-pkin(4) * t173 - pkin(5)); t179 * pkin(10) - m(7) * (t17 * t21 + t18 * t22 + t35 * t42) + t178 + (t249 + t268) * t42 - t22 * t69 - t21 * t70 - t41 * t115 - t302 * pkin(5); t87 - t35 * (mrSges(7,1) * t114 + mrSges(7,2) * t113) + (Ifges(7,1) * t113 - t265) * t285 + t46 * t284 + (Ifges(7,5) * t113 - Ifges(7,6) * t114) * t283 - t17 * t69 + t18 * t70 + (t113 * t17 + t114 * t18) * mrSges(7,3) + (-Ifges(7,2) * t114 + t108 + t47) * t287 + t293;];
tauc  = t1(:);
