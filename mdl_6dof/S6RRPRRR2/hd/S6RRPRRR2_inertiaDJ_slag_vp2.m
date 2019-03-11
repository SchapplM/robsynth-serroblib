% Calculate time derivative of joint inertia matrix for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:16:44
% EndTime: 2019-03-09 13:16:53
% DurationCPUTime: 4.04s
% Computational Cost: add. (12549->445), mult. (26713->653), div. (0->0), fcn. (27858->10), ass. (0->197)
t192 = sin(qJ(5));
t196 = cos(qJ(5));
t171 = -mrSges(6,1) * t196 + mrSges(6,2) * t192;
t284 = -mrSges(5,1) + t171;
t235 = qJD(5) * t196;
t191 = sin(qJ(6));
t195 = cos(qJ(6));
t208 = t191 * t192 - t195 * t196;
t277 = qJD(5) + qJD(6);
t123 = t277 * t208;
t165 = t191 * t196 + t192 * t195;
t124 = t277 * t165;
t238 = -Ifges(7,5) * t123 - Ifges(7,6) * t124;
t283 = Ifges(6,5) * t235 + t238;
t189 = sin(pkin(11));
t190 = cos(pkin(11));
t194 = sin(qJ(2));
t198 = cos(qJ(2));
t159 = -t189 * t194 + t190 * t198;
t160 = t189 * t198 + t190 * t194;
t193 = sin(qJ(4));
t197 = cos(qJ(4));
t119 = t159 * t193 + t160 * t197;
t153 = t160 * qJD(2);
t154 = t159 * qJD(2);
t209 = t197 * t159 - t160 * t193;
t89 = t209 * qJD(4) - t153 * t193 + t154 * t197;
t205 = t119 * t235 + t192 * t89;
t246 = t196 * t89;
t90 = t119 * qJD(4) + t197 * t153 + t154 * t193;
t282 = Ifges(6,5) * t246 + Ifges(6,3) * t90;
t181 = pkin(2) * t190 + pkin(3);
t262 = pkin(2) * t189;
t149 = t193 * t181 + t197 * t262;
t146 = pkin(9) + t149;
t258 = -pkin(10) - t146;
t218 = qJD(5) * t258;
t148 = t181 * t197 - t193 * t262;
t140 = t148 * qJD(4);
t240 = t140 * t196;
t105 = t192 * t218 + t240;
t241 = t140 * t192;
t106 = t196 * t218 - t241;
t130 = t258 * t192;
t186 = t196 * pkin(10);
t131 = t146 * t196 + t186;
t98 = t130 * t195 - t131 * t191;
t53 = t98 * qJD(6) + t105 * t195 + t106 * t191;
t99 = t130 * t191 + t131 * t195;
t54 = -t99 * qJD(6) - t105 * t191 + t106 * t195;
t281 = t54 * mrSges(7,1) - t53 * mrSges(7,2);
t265 = -pkin(10) - pkin(9);
t177 = t265 * t192;
t178 = pkin(9) * t196 + t186;
t133 = t177 * t195 - t178 * t191;
t226 = qJD(5) * t265;
t169 = t192 * t226;
t170 = t196 * t226;
t103 = t133 * qJD(6) + t169 * t195 + t170 * t191;
t134 = t177 * t191 + t178 * t195;
t104 = -t134 * qJD(6) - t169 * t191 + t170 * t195;
t280 = t104 * mrSges(7,1) - t103 * mrSges(7,2);
t257 = -qJ(3) - pkin(7);
t172 = t257 * t194;
t173 = t257 * t198;
t128 = t190 * t172 + t173 * t189;
t115 = -pkin(8) * t160 + t128;
t129 = t189 * t172 - t190 * t173;
t116 = pkin(8) * t159 + t129;
t74 = t115 * t193 + t116 * t197;
t71 = t196 * t74;
t183 = -pkin(2) * t198 - pkin(1);
t137 = -t159 * pkin(3) + t183;
t72 = -pkin(4) * t209 - t119 * pkin(9) + t137;
t48 = t192 * t72 + t71;
t279 = qJD(5) * t48;
t79 = t208 * t119;
t278 = t197 * t115 - t116 * t193;
t213 = mrSges(6,1) * t192 + mrSges(6,2) * t196;
t166 = t213 * qJD(5);
t276 = 2 * m(5);
t275 = 2 * m(6);
t274 = 2 * m(7);
t273 = -2 * mrSges(5,3);
t219 = qJD(2) * t257;
t150 = qJD(3) * t198 + t194 * t219;
t151 = -t194 * qJD(3) + t198 * t219;
t114 = t190 * t150 + t189 * t151;
t102 = -pkin(8) * t153 + t114;
t113 = -t150 * t189 + t151 * t190;
t202 = -pkin(8) * t154 + t113;
t41 = t74 * qJD(4) + t102 * t193 - t197 * t202;
t272 = 0.2e1 * t41;
t271 = -0.2e1 * t278;
t93 = t124 * mrSges(7,1) - t123 * mrSges(7,2);
t270 = 0.2e1 * t93;
t185 = qJD(2) * t194 * pkin(2);
t135 = pkin(3) * t153 + t185;
t269 = 0.2e1 * t135;
t268 = 0.2e1 * t183;
t267 = m(4) * pkin(2);
t261 = pkin(5) * t196;
t260 = t41 * t278;
t256 = Ifges(6,4) * t192;
t255 = Ifges(6,4) * t196;
t254 = Ifges(6,6) * t192;
t253 = pkin(5) * qJD(6);
t251 = t124 * mrSges(7,3);
t236 = qJD(5) * t192;
t40 = qJD(4) * t278 + t197 * t102 + t193 * t202;
t51 = pkin(4) * t90 - pkin(9) * t89 + t135;
t13 = t192 * t51 + t196 * t40 + t72 * t235 - t236 * t74;
t250 = t13 * t196;
t220 = -t192 * t40 + t196 * t51;
t14 = t220 - t279;
t249 = t14 * t192;
t141 = t149 * qJD(4);
t248 = t141 * t278;
t244 = t119 * t192;
t243 = t119 * t196;
t125 = mrSges(7,1) * t208 + mrSges(7,2) * t165;
t228 = pkin(5) * t236;
t132 = t141 + t228;
t242 = t132 * t125;
t237 = t192 ^ 2 + t196 ^ 2;
t234 = qJD(6) * t191;
t233 = qJD(6) * t195;
t232 = 0.2e1 * mrSges(7,3);
t231 = 0.2e1 * t198;
t26 = -t124 * t119 - t208 * t89;
t27 = -t165 * t89 + t277 * t79;
t230 = Ifges(7,5) * t26 + Ifges(7,6) * t27 + Ifges(7,3) * t90;
t229 = mrSges(7,3) * t253;
t227 = t195 * t123 * mrSges(7,3);
t225 = t119 * t236;
t223 = t90 * mrSges(5,1) + t89 * mrSges(5,2);
t222 = -t236 / 0.2e1;
t221 = -(2 * Ifges(5,4)) - t254;
t47 = -t192 * t74 + t196 * t72;
t217 = t284 * t141;
t216 = t237 * mrSges(6,3);
t215 = t153 * mrSges(4,1) + t154 * mrSges(4,2);
t214 = t237 * t140;
t145 = -pkin(4) - t148;
t212 = Ifges(6,1) * t196 - t256;
t211 = -Ifges(6,2) * t192 + t255;
t210 = Ifges(6,5) * t192 + Ifges(6,6) * t196;
t32 = -pkin(5) * t209 - pkin(10) * t243 + t47;
t37 = -pkin(10) * t244 + t48;
t16 = -t191 * t37 + t195 * t32;
t17 = t191 * t32 + t195 * t37;
t10 = -t205 * pkin(10) + t13;
t5 = -pkin(10) * t246 + pkin(5) * t90 + (-t71 + (pkin(10) * t119 - t72) * t192) * qJD(5) + t220;
t3 = t16 * qJD(6) + t10 * t195 + t191 * t5;
t4 = -t17 * qJD(6) - t10 * t191 + t195 * t5;
t207 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + t230;
t206 = -t195 * t208 * t229 + (-pkin(5) * t251 + t165 * t229) * t191 + t283;
t204 = t225 - t246;
t126 = Ifges(7,4) * t165 - Ifges(7,2) * t208;
t127 = Ifges(7,1) * t165 - Ifges(7,4) * t208;
t167 = t211 * qJD(5);
t168 = t212 * qJD(5);
t175 = Ifges(6,1) * t192 + t255;
t94 = -Ifges(7,4) * t123 - Ifges(7,2) * t124;
t95 = -Ifges(7,1) * t123 - Ifges(7,4) * t124;
t203 = -t123 * t127 - t124 * t126 + t165 * t95 + t196 * t167 + t192 * t168 + t175 * t235 - t208 * t94;
t201 = -t249 + (-t192 * t48 - t196 * t47) * qJD(5);
t174 = Ifges(6,2) * t196 + t256;
t200 = -t174 * t236 + t203;
t23 = t205 * pkin(5) + t41;
t30 = -t204 * Ifges(6,4) - t205 * Ifges(6,2) + t90 * Ifges(6,6);
t31 = -t204 * Ifges(6,1) - t205 * Ifges(6,4) + t90 * Ifges(6,5);
t78 = t165 * t119;
t42 = -Ifges(7,4) * t79 - Ifges(7,2) * t78 - Ifges(7,6) * t209;
t43 = -Ifges(7,1) * t79 - Ifges(7,4) * t78 - Ifges(7,5) * t209;
t58 = pkin(5) * t244 - t278;
t63 = -Ifges(6,6) * t209 + t211 * t119;
t64 = -Ifges(6,5) * t209 + t212 * t119;
t8 = Ifges(7,4) * t26 + Ifges(7,2) * t27 + Ifges(7,6) * t90;
t9 = Ifges(7,1) * t26 + Ifges(7,4) * t27 + Ifges(7,5) * t90;
t199 = t284 * t41 + t63 * t222 - t205 * t174 / 0.2e1 - (-Ifges(6,6) * t236 + t283) * t209 / 0.2e1 + mrSges(6,3) * t250 - t17 * t251 - t278 * t166 + (Ifges(7,5) * t165 - Ifges(7,6) * t208 + t210) * t90 / 0.2e1 - t208 * t8 / 0.2e1 + (t123 * t16 - t165 * t4 - t208 * t3) * mrSges(7,3) + t196 * t30 / 0.2e1 + t192 * t31 / 0.2e1 + (t119 * t222 + t246 / 0.2e1) * t175 + t64 * t235 / 0.2e1 - t40 * mrSges(5,2) + t168 * t243 / 0.2e1 - t167 * t244 / 0.2e1 + Ifges(5,5) * t89 - Ifges(5,6) * t90 + t58 * t93 - t78 * t94 / 0.2e1 - t79 * t95 / 0.2e1 - t123 * t43 / 0.2e1 - t124 * t42 / 0.2e1 + t23 * t125 + t27 * t126 / 0.2e1 + t26 * t127 / 0.2e1 + t165 * t9 / 0.2e1;
t182 = -pkin(4) - t261;
t158 = (-mrSges(7,1) * t191 - mrSges(7,2) * t195) * t253;
t136 = t145 - t261;
t92 = -mrSges(6,1) * t209 - mrSges(6,3) * t243;
t91 = mrSges(6,2) * t209 - mrSges(6,3) * t244;
t83 = t213 * t119;
t62 = -mrSges(7,1) * t209 + mrSges(7,3) * t79;
t61 = mrSges(7,2) * t209 - mrSges(7,3) * t78;
t55 = mrSges(7,1) * t78 - mrSges(7,2) * t79;
t45 = -mrSges(6,2) * t90 - t205 * mrSges(6,3);
t44 = mrSges(6,1) * t90 + t204 * mrSges(6,3);
t33 = t205 * mrSges(6,1) - t204 * mrSges(6,2);
t19 = -mrSges(7,2) * t90 + mrSges(7,3) * t27;
t18 = mrSges(7,1) * t90 - mrSges(7,3) * t26;
t11 = -mrSges(7,1) * t27 + mrSges(7,2) * t26;
t1 = [0.2e1 * m(4) * (t113 * t128 + t114 * t129) + t215 * t268 + t33 * t271 + t83 * t272 + (t16 * t4 + t17 * t3 + t23 * t58) * t274 - 0.2e1 * t159 * Ifges(4,2) * t153 + 0.2e1 * t160 * t154 * Ifges(4,1) + (t13 * t48 + t14 * t47 - t260) * t275 + (t135 * t137 + t40 * t74 - t260) * t276 - (mrSges(5,1) * t269 + t221 * t89 + t40 * t273 + t230 + t282) * t209 + (mrSges(5,2) * t269 + mrSges(5,3) * t272 + 0.2e1 * Ifges(5,1) * t89 - t192 * t30 + t196 * t31 + (-t192 * t64 - t196 * t63 + t209 * t210) * qJD(5)) * t119 + (-Ifges(7,5) * t79 - Ifges(7,6) * t78 + t74 * t273 + (Ifges(6,5) * t196 + t221) * t119 - ((2 * Ifges(5,2)) + Ifges(6,3) + Ifges(7,3)) * t209) * t90 + 0.2e1 * t137 * t223 + 0.2e1 * t16 * t18 + 0.2e1 * t17 * t19 + 0.2e1 * (-t153 * t160 + t154 * t159) * Ifges(4,4) + 0.2e1 * (-t113 * t160 + t114 * t159 - t128 * t154 - t129 * t153) * mrSges(4,3) + t27 * t42 + t26 * t43 + 0.2e1 * t47 * t44 + 0.2e1 * t48 * t45 + 0.2e1 * t23 * t55 + 0.2e1 * t58 * t11 + 0.2e1 * t3 * t61 + 0.2e1 * t4 * t62 - t78 * t8 - t79 * t9 + 0.2e1 * t13 * t91 + 0.2e1 * t14 * t92 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t198) * t231 + (t267 * t268 + 0.2e1 * pkin(2) * (-mrSges(4,1) * t159 + mrSges(4,2) * t160) - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t194 + (-Ifges(3,2) + Ifges(3,1)) * t231) * t194) * qJD(2) + (mrSges(5,3) * t271 - t192 * t63 + t196 * t64) * t89; (t113 * t190 + t114 * t189) * t267 + (-t153 * t189 - t154 * t190) * pkin(2) * mrSges(4,3) + t199 + (-t140 * t92 + (-qJD(5) * t91 - t44) * t146 + (-t14 - t279) * mrSges(6,3)) * t192 + (Ifges(3,5) * t198 - Ifges(3,6) * t194 + (-mrSges(3,1) * t198 + mrSges(3,2) * t194) * pkin(7)) * qJD(2) + (t140 * t91 + t146 * t45 + (-mrSges(6,3) * t47 - t146 * t92) * qJD(5)) * t196 + m(5) * (t140 * t74 - t148 * t41 + t149 * t40 - t248) + m(7) * (t132 * t58 + t136 * t23 + t16 * t54 + t17 * t53 + t3 * t99 + t4 * t98) + (t119 * t141 + t140 * t209 - t148 * t89 - t149 * t90) * mrSges(5,3) + t53 * t61 + t54 * t62 + t98 * t18 + t99 * t19 + t113 * mrSges(4,1) - t114 * mrSges(4,2) + t132 * t55 + t136 * t11 + t141 * t83 + t145 * t33 - Ifges(4,6) * t153 + Ifges(4,5) * t154 + ((t201 + t250) * t146 + t145 * t41 + t48 * t240 - t47 * t241 - t248) * m(6); (t123 * t98 - t124 * t99 - t165 * t54 - t208 * t53) * t232 + (t141 * t145 + t146 * t214) * t275 + (t132 * t136 + t53 * t99 + t54 * t98) * t274 - t141 * t148 * t276 + t200 + 0.2e1 * t242 + t136 * t270 + 0.2e1 * t145 * t166 + 0.2e1 * t217 + (t149 * t276 - 0.2e1 * mrSges(5,2) + 0.2e1 * t216) * t140; m(4) * t185 - t123 * t61 - t124 * t62 - t208 * t18 + t165 * t19 + t192 * t45 + t196 * t44 + (-t192 * t92 + t196 * t91) * qJD(5) + m(7) * (-t123 * t17 - t124 * t16 + t165 * t3 - t208 * t4) + m(6) * (t13 * t192 + t14 * t196 + (-t192 * t47 + t196 * t48) * qJD(5)) + m(5) * t135 + t215 + t223; m(7) * (-t123 * t99 - t124 * t98 + t165 * t53 - t208 * t54); (-t123 * t165 + t124 * t208) * t274; (m(6) * (-t235 * t47 - t236 * t48 - t249 + t250) + t196 * t45 - t192 * t44 - t92 * t235 - t91 * t236) * pkin(9) + t201 * mrSges(6,3) + t199 + m(7) * (t103 * t17 + t104 * t16 + t133 * t4 + t134 * t3 + t182 * t23 + t58 * t228) + t55 * t228 + t182 * t11 + t103 * t61 + t104 * t62 + t133 * t18 + t134 * t19 + (-m(6) * t41 - t33) * pkin(4); t242 + (t182 + t136) * t93 + (-pkin(4) + t145) * t166 + t217 + (pkin(5) * t125 - t174) * t236 + (-mrSges(5,2) + t216) * t140 + m(7) * (t103 * t99 + t104 * t98 + t132 * t182 + t133 * t54 + t134 * t53 + t136 * t228) + m(6) * (-pkin(4) * t141 + pkin(9) * t214) + ((-t104 - t54) * t165 - (t103 + t53) * t208 - (t134 + t99) * t124 - (-t133 - t98) * t123) * mrSges(7,3) + t203; m(7) * (t103 * t165 - t104 * t208 - t123 * t134 - t124 * t133); 0.2e1 * t125 * t228 + t182 * t270 - 0.2e1 * pkin(4) * t166 + (t103 * t134 + t104 * t133 + t182 * t228) * t274 + (-t103 * t208 - t104 * t165 + t123 * t133 - t124 * t134) * t232 + t200; -Ifges(6,5) * t225 + t14 * mrSges(6,1) - t13 * mrSges(6,2) - t205 * Ifges(6,6) + (m(7) * (-t16 * t234 + t17 * t233 + t191 * t3 + t195 * t4) + t61 * t233 + t191 * t19 - t62 * t234 + t195 * t18) * pkin(5) + t207 + t282; -t213 * t140 + (t146 * t171 - t254) * qJD(5) + (m(7) * (t191 * t53 + t195 * t54 + t233 * t99 - t234 * t98) + t227) * pkin(5) + t206 + t281; -t166 + m(7) * (-t123 * t191 - t124 * t195 + (t165 * t195 + t191 * t208) * qJD(6)) * pkin(5) - t93; (pkin(9) * t171 - t254) * qJD(5) + (m(7) * (t103 * t191 + t104 * t195 - t133 * t234 + t134 * t233) + t227) * pkin(5) + t206 + t280; 0.2e1 * t158; t207; t238 + t281; -t93; t238 + t280; t158; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
