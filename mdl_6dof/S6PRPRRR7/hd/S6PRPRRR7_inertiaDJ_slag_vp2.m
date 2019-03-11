% Calculate time derivative of joint inertia matrix for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_inertiaDJ_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:50:45
% EndTime: 2019-03-08 20:50:58
% DurationCPUTime: 5.62s
% Computational Cost: add. (10456->546), mult. (31889->838), div. (0->0), fcn. (34466->16), ass. (0->237)
t171 = sin(qJ(5));
t175 = cos(qJ(5));
t162 = sin(pkin(14));
t164 = sin(pkin(7));
t166 = cos(pkin(14));
t168 = cos(pkin(7));
t176 = cos(qJ(4));
t167 = cos(pkin(8));
t172 = sin(qJ(4));
t226 = t167 * t172;
t163 = sin(pkin(8));
t231 = t163 * t172;
t104 = t168 * t231 + (t162 * t176 + t166 * t226) * t164;
t225 = t167 * t176;
t230 = t163 * t176;
t284 = t164 * (-t162 * t172 + t166 * t225) + t168 * t230;
t259 = pkin(2) * t168;
t157 = t166 * t259;
t233 = t162 * t164;
t106 = pkin(3) * t168 + t157 + (-pkin(10) * t167 - qJ(3)) * t233;
t116 = (-pkin(10) * t162 * t163 - pkin(3) * t166 - pkin(2)) * t164;
t77 = -t106 * t163 + t167 * t116;
t53 = -pkin(4) * t284 - pkin(11) * t104 + t77;
t229 = t164 * t166;
t124 = -t163 * t229 + t167 * t168;
t128 = qJ(3) * t229 + t162 * t259;
t100 = (t163 * t168 + t167 * t229) * pkin(10) + t128;
t93 = t176 * t100;
t59 = t106 * t226 + t116 * t231 + t93;
t56 = pkin(11) * t124 + t59;
t251 = t171 * t53 + t175 * t56;
t280 = qJD(5) * t251;
t221 = qJD(3) * t164;
t190 = t106 * t167 + t116 * t163;
t58 = -t172 * t100 + t190 * t176;
t49 = (-t162 * t226 + t166 * t176) * t221 + t58 * qJD(4);
t207 = t162 * t221;
t199 = t163 * t207;
t97 = t284 * qJD(4);
t98 = t104 * qJD(4);
t72 = pkin(4) * t98 - pkin(11) * t97 + t199;
t13 = -t171 * t49 + t175 * t72 - t280;
t11 = -pkin(5) * t98 - t13;
t170 = sin(qJ(6));
t174 = cos(qJ(6));
t80 = t104 * t171 - t175 * t124;
t63 = -qJD(5) * t80 + t175 * t97;
t81 = t104 * t175 + t124 * t171;
t67 = -t170 * t81 - t174 * t284;
t32 = qJD(6) * t67 + t170 * t98 + t174 * t63;
t68 = -t170 * t284 + t174 * t81;
t33 = -qJD(6) * t68 - t170 * t63 + t174 * t98;
t14 = -mrSges(7,1) * t33 + mrSges(7,2) * t32;
t209 = -m(7) * t11 - t14;
t46 = mrSges(6,1) * t98 - mrSges(6,3) * t63;
t286 = -t209 - t46;
t215 = qJD(6) * t171;
t217 = qJD(5) * t175;
t183 = -t170 * t215 + t174 * t217;
t24 = -t171 * t56 + t175 * t53;
t20 = pkin(5) * t284 - t24;
t37 = -mrSges(7,1) * t67 + mrSges(7,2) * t68;
t70 = -mrSges(6,1) * t284 - mrSges(6,3) * t81;
t252 = t37 - t70;
t285 = -m(6) * t24 + m(7) * t20 + t252;
t283 = t170 / 0.2e1;
t260 = t174 / 0.2e1;
t250 = -mrSges(5,1) * t124 + mrSges(6,1) * t80 + mrSges(6,2) * t81 + mrSges(5,3) * t104;
t238 = -mrSges(6,1) * t175 + mrSges(6,2) * t171 - mrSges(5,1);
t269 = m(6) * pkin(4);
t281 = t238 - t269;
t279 = m(7) * pkin(12) + mrSges(7,3);
t146 = -mrSges(7,1) * t174 + mrSges(7,2) * t170;
t278 = -m(7) * pkin(5) - mrSges(6,1) + t146;
t276 = 0.2e1 * m(7);
t275 = 0.2e1 * pkin(11);
t274 = -0.2e1 * mrSges(5,3);
t273 = -2 * Ifges(5,4);
t272 = t164 ^ 2;
t271 = m(6) / 0.2e1;
t270 = m(7) / 0.2e1;
t268 = t33 / 0.2e1;
t245 = Ifges(7,4) * t170;
t149 = Ifges(7,2) * t174 + t245;
t244 = Ifges(7,4) * t174;
t195 = -Ifges(7,2) * t170 + t244;
t88 = -t149 * t215 + (Ifges(7,6) * t171 + t175 * t195) * qJD(5);
t267 = t88 / 0.2e1;
t151 = Ifges(7,1) * t170 + t244;
t196 = Ifges(7,1) * t174 - t245;
t89 = -t151 * t215 + (Ifges(7,5) * t171 + t175 * t196) * qJD(5);
t266 = t89 / 0.2e1;
t123 = -Ifges(7,5) * t175 + t171 * t196;
t265 = t123 / 0.2e1;
t264 = Ifges(7,5) * t283 + Ifges(7,6) * t260;
t263 = t151 / 0.2e1;
t262 = -t170 / 0.2e1;
t261 = -t174 / 0.2e1;
t258 = pkin(11) * t175;
t177 = cos(qJ(2));
t173 = sin(qJ(2));
t224 = t168 * t173;
t186 = t162 * t177 + t166 * t224;
t165 = sin(pkin(6));
t222 = qJD(2) * t165;
t227 = t164 * t173;
t178 = (t163 * t186 + t167 * t227) * t222;
t113 = (-t162 * t224 + t166 * t177) * t222;
t223 = t168 * t177;
t169 = cos(pkin(6));
t228 = t164 * t169;
t101 = t166 * t228 + (-t162 * t173 + t166 * t223) * t165;
t125 = -t164 * t165 * t177 + t169 * t168;
t191 = t101 * t167 + t125 * t163;
t102 = t165 * t173 * t166 + (t165 * t223 + t228) * t162;
t237 = t102 * t172;
t39 = t113 * t176 + (t176 * t191 - t237) * qJD(4) + (t163 * t227 - t167 * t186) * t172 * t222;
t66 = t102 * t176 + t172 * t191;
t79 = -t101 * t163 + t125 * t167;
t41 = t171 * t79 + t175 * t66;
t15 = qJD(5) * t41 + t171 * t39 - t175 * t178;
t40 = t171 * t66 - t79 * t175;
t257 = t15 * t40;
t50 = (t162 * t225 + t166 * t172) * t221 + (t172 * t190 + t93) * qJD(4);
t65 = -t101 * t225 - t125 * t230 + t237;
t256 = t50 * t65;
t181 = t186 * t222;
t208 = t173 * t222;
t200 = t164 * t208;
t38 = qJD(4) * t66 + t113 * t172 + t181 * t225 - t200 * t230;
t255 = t65 * t38;
t254 = t98 * Ifges(6,5);
t253 = t98 * Ifges(6,6);
t249 = Ifges(5,5) * t97 - Ifges(5,6) * t98;
t248 = mrSges(7,3) * t171;
t247 = Ifges(6,4) * t171;
t246 = Ifges(6,4) * t175;
t243 = Ifges(7,6) * t170;
t242 = t284 * Ifges(6,6);
t241 = t15 * t171;
t218 = qJD(5) * t171;
t16 = t171 * t178 + t39 * t175 + t217 * t79 - t218 * t66;
t240 = t16 * t175;
t239 = t176 * t50;
t129 = -t175 * t167 + t171 * t231;
t219 = qJD(4) * t176;
t205 = t163 * t219;
t108 = -qJD(5) * t129 + t175 * t205;
t236 = t108 * t175;
t130 = t167 * t171 + t175 * t231;
t109 = qJD(5) * t130 + t171 * t205;
t235 = t109 * t129;
t234 = t109 * t171;
t220 = qJD(4) * t172;
t216 = qJD(6) * t170;
t214 = qJD(6) * t174;
t213 = qJD(6) * t175;
t64 = qJD(5) * t81 + t171 * t97;
t5 = Ifges(7,5) * t32 + Ifges(7,6) * t33 + Ifges(7,3) * t64;
t212 = Ifges(6,5) * t63 - Ifges(6,6) * t64 + Ifges(6,3) * t98;
t206 = t163 * t220;
t145 = -pkin(5) * t175 - pkin(12) * t171 - pkin(4);
t118 = t145 * t174 - t170 * t258;
t144 = (pkin(5) * t171 - pkin(12) * t175) * qJD(5);
t84 = t145 * t214 + t144 * t170 + (-t170 * t213 - t174 * t218) * pkin(11);
t202 = -qJD(6) * t118 + t84;
t119 = t145 * t170 + t174 * t258;
t85 = -t145 * t216 + t144 * t174 + (t170 * t218 - t174 * t213) * pkin(11);
t201 = -qJD(6) * t119 - t85;
t12 = t171 * t72 + t175 * t49 + t53 * t217 - t218 * t56;
t10 = pkin(12) * t98 + t12;
t19 = pkin(5) * t64 - pkin(12) * t63 + t50;
t21 = -pkin(12) * t284 + t251;
t55 = -pkin(4) * t124 - t58;
t31 = pkin(5) * t80 - pkin(12) * t81 + t55;
t8 = -t170 * t21 + t174 * t31;
t1 = qJD(6) * t8 + t10 * t174 + t170 * t19;
t9 = t170 * t31 + t174 * t21;
t2 = -qJD(6) * t9 - t10 * t170 + t174 * t19;
t198 = t1 * t174 - t2 * t170;
t197 = mrSges(7,1) * t170 + mrSges(7,2) * t174;
t194 = t109 * t40 + t129 * t15;
t23 = t170 * t65 + t174 * t41;
t22 = -t170 * t41 + t174 * t65;
t27 = Ifges(7,4) * t68 + Ifges(7,2) * t67 + Ifges(7,6) * t80;
t28 = Ifges(7,1) * t68 + Ifges(7,4) * t67 + Ifges(7,5) * t80;
t189 = t260 * t28 + t262 * t27;
t188 = t217 * t40 + t241;
t110 = -t130 * t170 - t174 * t230;
t187 = -t130 * t174 + t170 * t230;
t185 = t129 * t217 + t234;
t182 = t170 * t217 + t171 * t214;
t87 = t183 * Ifges(7,5) - Ifges(7,6) * t182 + Ifges(7,3) * t218;
t161 = Ifges(6,5) * t217;
t160 = Ifges(7,5) * t214;
t152 = Ifges(6,1) * t171 + t246;
t150 = Ifges(6,2) * t175 + t247;
t143 = -mrSges(7,1) * t175 - t174 * t248;
t142 = mrSges(7,2) * t175 - t170 * t248;
t140 = (Ifges(6,1) * t175 - t247) * qJD(5);
t139 = t196 * qJD(6);
t138 = (-Ifges(6,2) * t171 + t246) * qJD(5);
t137 = t195 * qJD(6);
t136 = -Ifges(7,6) * t216 + t160;
t135 = (mrSges(6,1) * t171 + mrSges(6,2) * t175) * qJD(5);
t134 = t197 * qJD(6);
t133 = -mrSges(4,2) * t168 + mrSges(4,3) * t229;
t132 = mrSges(4,1) * t168 - mrSges(4,3) * t233;
t131 = t197 * t171;
t127 = -qJ(3) * t233 + t157;
t122 = -Ifges(7,6) * t175 + t171 * t195;
t121 = -Ifges(7,3) * t175 + (Ifges(7,5) * t174 - t243) * t171;
t115 = -mrSges(7,2) * t218 - mrSges(7,3) * t182;
t114 = mrSges(7,1) * t218 - mrSges(7,3) * t183;
t99 = mrSges(7,1) * t182 + mrSges(7,2) * t183;
t82 = -mrSges(5,2) * t124 + mrSges(5,3) * t284;
t76 = -mrSges(5,1) * t284 + mrSges(5,2) * t104;
t75 = qJD(6) * t187 - t108 * t170 + t174 * t206;
t74 = qJD(6) * t110 + t108 * t174 + t170 * t206;
t73 = mrSges(5,1) * t98 + mrSges(5,2) * t97;
t69 = mrSges(6,2) * t284 - mrSges(6,3) * t80;
t47 = -mrSges(6,2) * t98 - mrSges(6,3) * t64;
t45 = Ifges(6,1) * t81 - Ifges(6,4) * t80 - Ifges(6,5) * t284;
t44 = Ifges(6,4) * t81 - Ifges(6,2) * t80 - t242;
t43 = mrSges(7,1) * t80 - mrSges(7,3) * t68;
t42 = -mrSges(7,2) * t80 + mrSges(7,3) * t67;
t36 = mrSges(6,1) * t64 + mrSges(6,2) * t63;
t35 = Ifges(6,1) * t63 - Ifges(6,4) * t64 + t254;
t34 = Ifges(6,4) * t63 - Ifges(6,2) * t64 + t253;
t26 = Ifges(7,5) * t68 + Ifges(7,6) * t67 + Ifges(7,3) * t80;
t18 = -mrSges(7,2) * t64 + mrSges(7,3) * t33;
t17 = mrSges(7,1) * t64 - mrSges(7,3) * t32;
t7 = Ifges(7,1) * t32 + Ifges(7,4) * t33 + Ifges(7,5) * t64;
t6 = Ifges(7,4) * t32 + Ifges(7,2) * t33 + Ifges(7,6) * t64;
t4 = qJD(6) * t22 + t16 * t174 + t170 * t38;
t3 = -qJD(6) * t23 - t16 * t170 + t174 * t38;
t25 = [0.2e1 * m(7) * (t22 * t3 + t23 * t4 + t257) + 0.2e1 * m(6) * (t16 * t41 + t255 + t257) + 0.2e1 * m(5) * (t178 * t79 + t66 * t39 + t255) + 0.2e1 * m(4) * (t102 * t113 + (-t101 * t186 + t125 * t227) * t222); t113 * t133 + m(7) * (t1 * t23 + t2 * t22 + t3 * t8 + t4 * t9) + t76 * t178 + t39 * t82 + t79 * t73 + t16 * t69 + t65 * t36 + t4 * t42 + t3 * t43 + t41 * t47 + t22 * t17 + t23 * t18 - t177 * mrSges(3,2) * t222 - t132 * t181 + m(6) * (t12 * t41 + t16 * t251 + t256) + m(5) * (t178 * t77 + t199 * t79 + t59 * t39 + t49 * t66 + t256) + m(4) * (t128 * t113 + (-t101 * t162 + t102 * t166) * t221 + (-pkin(2) * t173 * t272 - t127 * t186) * t222) + (-m(6) * t13 + t286) * t40 + (-m(5) * t58 + m(6) * t55 + t250) * t38 + (t272 * (-mrSges(4,1) * t166 + mrSges(4,2) * t162) - mrSges(3,1)) * t208 + t285 * t15 + (t65 * t97 - t66 * t98) * mrSges(5,3); 0.2e1 * t250 * t50 + 0.2e1 * t49 * t82 + t81 * t35 + 0.2e1 * t77 * t73 + t67 * t6 + t68 * t7 + 0.2e1 * t12 * t69 + 0.2e1 * t13 * t70 + t63 * t45 + 0.2e1 * t55 * t36 + 0.2e1 * t2 * t43 + 0.2e1 * t24 * t46 + 0.2e1 * t1 * t42 + 0.2e1 * t11 * t37 + t32 * t28 + t33 * t27 + 0.2e1 * t20 * t14 + 0.2e1 * t8 * t17 + 0.2e1 * t9 * t18 + (t1 * t9 + t11 * t20 + t2 * t8) * t276 + (t5 - t34) * t80 + (t26 - t44) * t64 + 0.2e1 * (t133 * t166 + (t163 * t76 - t132) * t162 + m(4) * (-t127 * t162 + t128 * t166)) * t221 + (t59 * t274 + t104 * t273 + Ifges(6,5) * t81 - Ifges(5,6) * t124 - Ifges(6,6) * t80 - ((2 * Ifges(5,2)) + Ifges(6,3)) * t284) * t98 + (0.2e1 * Ifges(5,1) * t104 + Ifges(5,5) * t124 - t273 * t284 + t274 * t58) * t97 - t284 * t212 + 0.2e1 * t251 * t47 + 0.2e1 * m(6) * (t12 * t251 + t13 * t24 + t50 * t55) + 0.2e1 * m(5) * (t199 * t77 + t49 * t59 - t50 * t58) + t124 * t249; m(7) * (t110 * t3 - t187 * t4 + t22 * t75 + t23 * t74 + t194) + m(6) * (t108 * t41 + t130 * t16 + (-t176 * t38 + t220 * t65) * t163 + t194) + m(5) * (t167 * t178 + t205 * t66 + t206 * t65 - t230 * t38 + t231 * t39) + m(4) * t200; t108 * t69 + t110 * t17 - t187 * t18 + t130 * t47 + t167 * t73 + t74 * t42 + t75 * t43 + (t14 - t46) * t129 + t252 * t109 + m(7) * (-t1 * t187 + t109 * t20 + t11 * t129 + t110 * t2 + t74 * t9 + t75 * t8) + m(6) * (t108 * t251 - t109 * t24 + t12 * t130 - t129 * t13) + (-t176 * t36 + (-t172 * t98 - t176 * t97) * mrSges(5,3) + (t172 * t250 + t176 * t82) * qJD(4) + m(6) * (t220 * t55 - t239) + m(5) * (t167 * t207 + t172 * t49 + t219 * t59 - t220 * t58 - t239)) * t163; 0.2e1 * m(6) * (-t163 ^ 2 * t172 * t219 + t108 * t130 + t235) + 0.2e1 * m(7) * (t110 * t75 - t187 * t74 + t235); -t39 * mrSges(5,2) + t22 * t114 + t23 * t115 + t15 * t131 + t65 * t135 + t4 * t142 + t3 * t143 + t40 * t99 + m(7) * (t118 * t3 + t119 * t4 + t22 * t85 + t23 * t84) + (t188 * t270 + (-t218 * t41 + t188 + t240) * t271) * t275 + (t241 + t240 + (-t171 * t41 + t175 * t40) * qJD(5)) * mrSges(6,3) + t281 * t38; m(7) * (t1 * t119 + t118 * t2 + t8 * t85 + t84 * t9) + t81 * t140 / 0.2e1 + t1 * t142 + t2 * t143 + t63 * t152 / 0.2e1 + t11 * t131 + t55 * t135 + t8 * t114 + t9 * t115 + t118 * t17 + t119 * t18 + t20 * t99 + t84 * t42 + t85 * t43 - t49 * mrSges(5,2) - pkin(4) * t36 + (-t138 / 0.2e1 + t87 / 0.2e1) * t80 + (-t150 / 0.2e1 + t121 / 0.2e1) * t64 + t68 * t266 + t67 * t267 + t122 * t268 + t32 * t265 + (t7 * t260 + t6 * t262 - t13 * mrSges(6,3) + t35 / 0.2e1 + t254 / 0.2e1 + (t261 * t27 + t262 * t28) * qJD(6) + (-t44 / 0.2e1 + t26 / 0.2e1 - t251 * mrSges(6,3) + t242 / 0.2e1) * qJD(5) + (-qJD(5) * t69 + m(6) * (-t13 - t280) + t286) * pkin(11)) * t171 + t249 + t281 * t50 - t284 * t161 / 0.2e1 + (t12 * mrSges(6,3) + t34 / 0.2e1 - t5 / 0.2e1 + t253 / 0.2e1 + (t45 / 0.2e1 - t24 * mrSges(6,3) + t189) * qJD(5) + (m(6) * t12 + t285 * qJD(5) + t47) * pkin(11)) * t175; t109 * t131 + t110 * t114 - t187 * t115 + t129 * t99 + t74 * t142 + t75 * t143 + (-t176 * t135 + (-mrSges(5,2) * t176 + t172 * t238) * qJD(4)) * t163 - t206 * t269 + m(7) * (t110 * t85 + t118 * t75 + t119 * t74 - t187 * t84) + ((-t130 * t218 + t185 + t236) * t271 + t185 * t270) * t275 + (t236 + t234 + (t129 * t175 - t130 * t171) * qJD(5)) * mrSges(6,3); 0.2e1 * t84 * t142 + 0.2e1 * t119 * t115 + 0.2e1 * t85 * t143 + 0.2e1 * t118 * t114 + (t118 * t85 + t119 * t84) * t276 - 0.2e1 * pkin(4) * t135 + (t138 - t87 + (-t170 * t122 + t174 * t123 + t131 * t275 + t152) * qJD(5)) * t175 + (t99 * t275 - t170 * t88 + t174 * t89 + t140 + (-t122 * t174 - t123 * t170) * qJD(6) + (pkin(11) ^ 2 * t175 * t276 + t121 - t150) * qJD(5)) * t171; -t16 * mrSges(6,2) + t40 * t134 + t279 * (-t3 * t170 + t4 * t174 + (-t170 * t23 - t174 * t22) * qJD(6)) + t278 * t15; t67 * t137 / 0.2e1 + t68 * t139 / 0.2e1 + t11 * t146 + t64 * t264 + t149 * t268 + t32 * t263 + t20 * t134 + t80 * t136 / 0.2e1 + t13 * mrSges(6,1) - t12 * mrSges(6,2) + t7 * t283 + t6 * t260 + t189 * qJD(6) + t209 * pkin(5) + ((-t170 * t9 - t174 * t8) * qJD(6) + t198) * mrSges(7,3) + (m(7) * (-t214 * t8 - t216 * t9 + t198) + t174 * t18 - t170 * t17 - t43 * t214 - t42 * t216) * pkin(12) + t212; -t108 * mrSges(6,2) + t129 * t134 + t279 * (-t75 * t170 + t74 * t174 + (-t110 * t174 + t170 * t187) * qJD(6)) + t278 * t109; -pkin(5) * t99 + t161 + (-t136 / 0.2e1 + t278 * qJD(5) * pkin(11)) * t175 + (qJD(6) * t265 + t217 * t263 + t267 + t202 * mrSges(7,3) + (m(7) * t202 - qJD(6) * t143 + t115) * pkin(12)) * t174 + (-qJD(6) * t122 / 0.2e1 - t149 * t217 / 0.2e1 + t266 + t201 * mrSges(7,3) + (m(7) * t201 - qJD(6) * t142 - t114) * pkin(12)) * t170 + (t139 * t260 + t137 * t262 + pkin(11) * t134 + (t149 * t261 + t151 * t262) * qJD(6) + (pkin(11) * mrSges(6,2) - Ifges(6,6) + t264) * qJD(5)) * t171; -0.2e1 * pkin(5) * t134 + t137 * t174 + t139 * t170 + (-t149 * t170 + t151 * t174) * qJD(6); mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t5; mrSges(7,1) * t75 - mrSges(7,2) * t74; mrSges(7,1) * t85 - mrSges(7,2) * t84 + t87; t160 + (pkin(12) * t146 - t243) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
