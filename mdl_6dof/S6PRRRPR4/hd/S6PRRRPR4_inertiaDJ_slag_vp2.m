% Calculate time derivative of joint inertia matrix for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:15:44
% EndTime: 2019-03-08 23:15:53
% DurationCPUTime: 4.01s
% Computational Cost: add. (5505->501), mult. (14224->756), div. (0->0), fcn. (13272->12), ass. (0->208)
t250 = -Ifges(5,3) - Ifges(6,3);
t183 = sin(qJ(4));
t184 = sin(qJ(3));
t187 = cos(qJ(4));
t216 = qJD(4) * t187;
t188 = cos(qJ(3));
t219 = qJD(3) * t188;
t191 = t183 * t219 + t184 * t216;
t178 = sin(pkin(12));
t180 = cos(pkin(12));
t143 = t178 * t187 + t180 * t183;
t196 = t178 * t183 - t180 * t187;
t217 = qJD(4) * t184;
t79 = -t143 * t219 + t196 * t217;
t134 = t143 * qJD(4);
t80 = -t134 * t184 - t196 * t219;
t42 = -t79 * mrSges(6,1) + t80 * mrSges(6,2);
t182 = sin(qJ(6));
t186 = cos(qJ(6));
t125 = t143 * t184;
t126 = t196 * t184;
t74 = -t125 * t186 + t126 * t182;
t31 = qJD(6) * t74 + t182 * t79 + t186 * t80;
t75 = -t125 * t182 - t126 * t186;
t32 = -qJD(6) * t75 - t182 * t80 + t186 * t79;
t9 = -t32 * mrSges(7,1) + t31 * mrSges(7,2);
t249 = t42 + t9;
t157 = -pkin(3) * t188 - pkin(9) * t184 - pkin(2);
t224 = t187 * t188;
t168 = pkin(8) * t224;
t215 = qJD(5) * t187;
t155 = (pkin(3) * t184 - pkin(9) * t188) * qJD(3);
t220 = qJD(3) * t184;
t238 = pkin(8) * t183;
t222 = t187 * t155 + t220 * t238;
t48 = -t184 * t215 + (pkin(4) * t184 - qJ(5) * t224) * qJD(3) + (-t168 + (qJ(5) * t184 - t157) * t183) * qJD(4) + t222;
t223 = t183 * t155 + t157 * t216;
t225 = t184 * t187;
t57 = (-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t225 + (-qJD(5) * t184 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t188) * t183 + t223;
t17 = -t178 * t57 + t180 * t48;
t13 = pkin(5) * t220 - pkin(10) * t80 + t17;
t19 = t178 * t48 + t180 * t57;
t14 = pkin(10) * t79 + t19;
t116 = t183 * t157 + t168;
t226 = t183 * t184;
t106 = -qJ(5) * t226 + t116;
t145 = t187 * t157;
t96 = -qJ(5) * t225 + t145 + (-pkin(4) - t238) * t188;
t58 = -t106 * t178 + t180 * t96;
t38 = -pkin(5) * t188 + pkin(10) * t126 + t58;
t59 = t180 * t106 + t178 * t96;
t41 = -pkin(10) * t125 + t59;
t15 = -t182 * t41 + t186 * t38;
t2 = qJD(6) * t15 + t13 * t182 + t14 * t186;
t16 = t182 * t38 + t186 * t41;
t3 = -qJD(6) * t16 + t13 * t186 - t14 * t182;
t248 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t135 = t196 * qJD(4);
t94 = -t143 * t182 - t186 * t196;
t49 = qJD(6) * t94 - t134 * t182 - t135 * t186;
t95 = t143 * t186 - t182 * t196;
t50 = -qJD(6) * t95 - t134 * t186 + t135 * t182;
t21 = -t50 * mrSges(7,1) + t49 * mrSges(7,2);
t88 = t134 * mrSges(6,1) - t135 * mrSges(6,2);
t247 = t21 + t88;
t246 = 2 * m(5);
t245 = 2 * m(6);
t244 = 2 * m(7);
t243 = 0.2e1 * pkin(8);
t235 = Ifges(5,4) * t183;
t161 = Ifges(5,2) * t187 + t235;
t241 = -t161 / 0.2e1;
t240 = -t183 / 0.2e1;
t239 = pkin(4) * t178;
t236 = -qJ(5) - pkin(9);
t234 = Ifges(5,4) * t187;
t233 = Ifges(5,6) * t183;
t232 = t188 * Ifges(5,6);
t159 = -mrSges(5,1) * t187 + mrSges(5,2) * t183;
t231 = -mrSges(4,1) + t159;
t181 = cos(pkin(6));
t179 = sin(pkin(6));
t185 = sin(qJ(2));
t228 = t179 * t185;
t137 = t181 * t184 + t188 * t228;
t189 = cos(qJ(2));
t221 = qJD(2) * t189;
t210 = t179 * t221;
t102 = qJD(3) * t137 + t184 * t210;
t230 = t102 * t184;
t136 = -t181 * t188 + t184 * t228;
t103 = -qJD(3) * t136 + t188 * t210;
t229 = t103 * t188;
t73 = t136 * t102;
t227 = t179 * t189;
t204 = qJD(4) * t236;
t132 = t183 * t204 + t215;
t133 = -qJD(5) * t183 + t187 * t204;
t82 = t180 * t132 + t178 * t133;
t158 = t236 * t183;
t160 = t236 * t187;
t108 = t178 * t158 - t180 * t160;
t156 = pkin(4) * t226 + t184 * pkin(8);
t218 = qJD(4) * t183;
t214 = -Ifges(7,5) * t31 - Ifges(7,6) * t32 - Ifges(7,3) * t220;
t213 = pkin(4) * t218;
t176 = pkin(8) * t219;
t195 = -t137 * t187 + t183 * t227;
t211 = qJD(2) * t228;
t52 = qJD(4) * t195 - t103 * t183 + t187 * t211;
t104 = -t137 * t183 - t187 * t227;
t53 = qJD(4) * t104 + t103 * t187 + t183 * t211;
t18 = -t178 * t53 + t180 * t52;
t20 = t178 * t52 + t180 * t53;
t61 = t104 * t180 + t178 * t195;
t62 = t104 * t178 - t180 * t195;
t29 = -t182 * t62 + t186 * t61;
t5 = qJD(6) * t29 + t18 * t182 + t186 * t20;
t30 = t182 * t61 + t186 * t62;
t6 = -qJD(6) * t30 + t18 * t186 - t182 * t20;
t212 = t6 * mrSges(7,1) - t5 * mrSges(7,2);
t114 = pkin(4) * t191 + t176;
t171 = -pkin(4) * t187 - pkin(3);
t208 = t183 * t217;
t206 = t187 * t219;
t205 = (2 * Ifges(4,4)) + t233;
t115 = -t188 * t238 + t145;
t69 = (-t187 * t220 - t188 * t218) * pkin(8) + t223;
t203 = -qJD(4) * t115 + t69;
t81 = -t132 * t178 + t180 * t133;
t107 = t180 * t158 + t160 * t178;
t202 = mrSges(5,1) * t183 + mrSges(5,2) * t187;
t170 = pkin(4) * t180 + pkin(5);
t130 = t170 * t186 - t182 * t239;
t121 = t130 * qJD(6);
t131 = t170 * t182 + t186 * t239;
t122 = t131 * qJD(6);
t201 = -t122 * mrSges(7,1) - t121 * mrSges(7,2);
t200 = Ifges(5,1) * t187 - t235;
t162 = Ifges(5,1) * t183 + t234;
t199 = -Ifges(5,2) * t183 + t234;
t198 = Ifges(5,5) * t183 + Ifges(5,6) * t187;
t83 = -pkin(10) * t143 + t107;
t84 = -pkin(10) * t196 + t108;
t39 = -t182 * t84 + t186 * t83;
t40 = t182 * t83 + t186 * t84;
t63 = pkin(10) * t135 + t81;
t64 = -pkin(10) * t134 + t82;
t11 = qJD(6) * t39 + t182 * t63 + t186 * t64;
t12 = -qJD(6) * t40 - t182 * t64 + t186 * t63;
t46 = Ifges(7,6) * t50;
t47 = Ifges(7,5) * t49;
t197 = t12 * mrSges(7,1) - t11 * mrSges(7,2) + t46 + t47;
t194 = t136 * t219 + t230;
t193 = -Ifges(5,5) * t206 - Ifges(6,5) * t80 - Ifges(6,6) * t79 + t220 * t250 + t214;
t192 = t206 - t208;
t190 = -t52 * t183 + t53 * t187 + (-t104 * t187 + t183 * t195) * qJD(4);
t175 = Ifges(5,5) * t216;
t154 = -mrSges(5,1) * t188 - mrSges(5,3) * t225;
t153 = mrSges(5,2) * t188 - mrSges(5,3) * t226;
t152 = t200 * qJD(4);
t151 = t199 * qJD(4);
t150 = (mrSges(4,1) * t184 + mrSges(4,2) * t188) * qJD(3);
t149 = t202 * qJD(4);
t140 = t202 * t184;
t129 = Ifges(6,5) * t135;
t128 = Ifges(6,6) * t134;
t124 = -Ifges(5,5) * t188 + t184 * t200;
t123 = t184 * t199 - t232;
t117 = pkin(5) * t196 + t171;
t113 = -mrSges(5,2) * t220 - mrSges(5,3) * t191;
t112 = mrSges(5,1) * t220 - mrSges(5,3) * t192;
t111 = pkin(5) * t134 + t213;
t110 = -mrSges(6,1) * t188 + mrSges(6,3) * t126;
t109 = mrSges(6,2) * t188 - mrSges(6,3) * t125;
t101 = Ifges(6,1) * t143 - Ifges(6,4) * t196;
t100 = Ifges(6,4) * t143 - Ifges(6,2) * t196;
t99 = mrSges(6,1) * t196 + mrSges(6,2) * t143;
t97 = pkin(5) * t125 + t156;
t93 = mrSges(5,1) * t191 + mrSges(5,2) * t192;
t90 = -Ifges(6,1) * t135 - Ifges(6,4) * t134;
t89 = -Ifges(6,4) * t135 - Ifges(6,2) * t134;
t87 = -t162 * t217 + (Ifges(5,5) * t184 + t188 * t200) * qJD(3);
t86 = -t161 * t217 + (Ifges(5,6) * t184 + t188 * t199) * qJD(3);
t85 = mrSges(6,1) * t125 - mrSges(6,2) * t126;
t72 = -Ifges(6,1) * t126 - Ifges(6,4) * t125 - Ifges(6,5) * t188;
t71 = -Ifges(6,4) * t126 - Ifges(6,2) * t125 - Ifges(6,6) * t188;
t70 = -qJD(4) * t116 + t222;
t68 = mrSges(6,1) * t220 - mrSges(6,3) * t80;
t67 = -mrSges(6,2) * t220 + mrSges(6,3) * t79;
t66 = -mrSges(7,1) * t188 - mrSges(7,3) * t75;
t65 = mrSges(7,2) * t188 + mrSges(7,3) * t74;
t60 = -pkin(5) * t79 + t114;
t56 = Ifges(7,1) * t95 + Ifges(7,4) * t94;
t55 = Ifges(7,4) * t95 + Ifges(7,2) * t94;
t54 = -mrSges(7,1) * t94 + mrSges(7,2) * t95;
t37 = -mrSges(7,1) * t74 + mrSges(7,2) * t75;
t36 = Ifges(6,1) * t80 + Ifges(6,4) * t79 + Ifges(6,5) * t220;
t35 = Ifges(6,4) * t80 + Ifges(6,2) * t79 + Ifges(6,6) * t220;
t34 = Ifges(7,1) * t75 + Ifges(7,4) * t74 - Ifges(7,5) * t188;
t33 = Ifges(7,4) * t75 + Ifges(7,2) * t74 - Ifges(7,6) * t188;
t25 = -mrSges(7,2) * t220 + mrSges(7,3) * t32;
t24 = mrSges(7,1) * t220 - mrSges(7,3) * t31;
t23 = Ifges(7,1) * t49 + Ifges(7,4) * t50;
t22 = Ifges(7,4) * t49 + Ifges(7,2) * t50;
t8 = Ifges(7,1) * t31 + Ifges(7,4) * t32 + Ifges(7,5) * t220;
t7 = Ifges(7,4) * t31 + Ifges(7,2) * t32 + Ifges(7,6) * t220;
t1 = [0.2e1 * m(7) * (t29 * t6 + t30 * t5 + t73) + 0.2e1 * m(6) * (t18 * t61 + t20 * t62 + t73) + 0.2e1 * m(5) * (t104 * t52 - t195 * t53 + t73) + 0.2e1 * m(4) * (-t179 ^ 2 * t185 * t221 + t137 * t103 + t73); t104 * t112 - t195 * t113 + t20 * t109 + t18 * t110 + t53 * t153 + t52 * t154 + t29 * t24 + t30 * t25 + t5 * t65 + t6 * t66 + t61 * t68 + t62 * t67 + (t93 + t249) * t136 + (t140 + t85 + t37) * t102 + (-t189 * t150 + (-t189 * mrSges(3,2) + (-mrSges(4,1) * t188 + mrSges(4,2) * t184 - mrSges(3,1)) * t185) * qJD(2)) * t179 + (t230 + t229 + (t136 * t188 - t137 * t184) * qJD(3)) * mrSges(4,3) + m(5) * (t104 * t70 + t115 * t52 + t116 * t53 - t195 * t69) - m(4) * pkin(2) * t211 + m(6) * (t102 * t156 + t114 * t136 + t17 * t61 + t18 * t58 + t19 * t62 + t20 * t59) + m(7) * (t102 * t97 + t136 * t60 + t15 * t6 + t16 * t5 + t2 * t30 + t29 * t3) + (m(5) * t194 / 0.2e1 + m(4) * (-t137 * t220 + t194 + t229) / 0.2e1) * t243; (t15 * t3 + t16 * t2 + t60 * t97) * t244 + (t114 * t156 + t17 * t58 + t19 * t59) * t245 + (t115 * t70 + t116 * t69) * t246 - 0.2e1 * pkin(2) * t150 + 0.2e1 * t69 * t153 + 0.2e1 * t70 * t154 + 0.2e1 * t156 * t42 - t125 * t35 - t126 * t36 + 0.2e1 * t114 * t85 + 0.2e1 * t115 * t112 + 0.2e1 * t116 * t113 + 0.2e1 * t19 * t109 + 0.2e1 * t17 * t110 + 0.2e1 * t97 * t9 + t80 * t72 + t74 * t7 + t75 * t8 + t79 * t71 + 0.2e1 * t2 * t65 + 0.2e1 * t3 * t66 + 0.2e1 * t59 * t67 + 0.2e1 * t58 * t68 + 0.2e1 * t60 * t37 + t32 * t33 + t31 * t34 + 0.2e1 * t15 * t24 + 0.2e1 * t16 * t25 + ((-t183 * t123 + t187 * t124 + t140 * t243 + t188 * t205) * qJD(3) + t193) * t188 + (t93 * t243 - t183 * t86 + t187 * t87 + (-t187 * t123 - t183 * t124 + t188 * t198) * qJD(4) + (-Ifges(6,5) * t126 - Ifges(6,6) * t125 + Ifges(7,5) * t75 + Ifges(7,6) * t74 + (Ifges(5,5) * t187 - t205) * t184 + (pkin(8) ^ 2 * t246 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(7,3) + t250) * t188) * qJD(3)) * t184; -t103 * mrSges(4,2) + (t149 + t247) * t136 + (t54 + t99 + t231) * t102 + m(7) * (t102 * t117 + t11 * t30 + t111 * t136 + t12 * t29 + t39 * t6 + t40 * t5) + m(6) * (t102 * t171 + t107 * t18 + t108 * t20 + t136 * t213 + t61 * t81 + t62 * t82) + (-t29 * t49 + t30 * t50 + t5 * t94 - t6 * t95) * mrSges(7,3) + (-t134 * t62 + t135 * t61 - t143 * t18 - t196 * t20) * mrSges(6,3) + t190 * mrSges(5,3) + (-pkin(3) * t102 + pkin(9) * t190) * m(5); -t196 * t35 / 0.2e1 + (-t134 * t59 + t135 * t58 - t143 * t17 - t19 * t196) * mrSges(6,3) + (t187 * t152 / 0.2e1 + t151 * t240 - Ifges(4,6) * qJD(3) + (t162 * t240 + t187 * t241) * qJD(4) + (mrSges(4,2) * qJD(3) + t149) * pkin(8) + (Ifges(6,5) * t143 + Ifges(7,5) * t95 - Ifges(6,6) * t196 + Ifges(7,6) * t94 + t198) * qJD(3) / 0.2e1) * t184 + (t129 / 0.2e1 + t128 / 0.2e1 - t47 / 0.2e1 - t46 / 0.2e1 - t175 / 0.2e1 + (pkin(8) * t231 + Ifges(4,5)) * qJD(3)) * t188 + t171 * t42 + t143 * t36 / 0.2e1 + t156 * t88 - t134 * t71 / 0.2e1 - t135 * t72 / 0.2e1 - t125 * t89 / 0.2e1 - t126 * t90 / 0.2e1 + t114 * t99 + t117 * t9 + t107 * t68 + t108 * t67 + t82 * t109 + t81 * t110 + t111 * t37 - pkin(3) * t93 + t94 * t7 / 0.2e1 + t95 * t8 / 0.2e1 + t97 * t21 + t79 * t100 / 0.2e1 + t80 * t101 / 0.2e1 + t74 * t22 / 0.2e1 + t75 * t23 / 0.2e1 + t11 * t65 + t12 * t66 + t60 * t54 + t49 * t34 / 0.2e1 + t50 * t33 / 0.2e1 + t32 * t55 / 0.2e1 + t31 * t56 / 0.2e1 + t40 * t25 + m(7) * (t11 * t16 + t111 * t97 + t117 * t60 + t12 * t15 + t2 * t40 + t3 * t39) + t39 * t24 + (t87 / 0.2e1 - pkin(9) * t112 - t70 * mrSges(5,3) + t219 * t241 + (t232 / 0.2e1 - pkin(9) * t153 - t116 * mrSges(5,3) + pkin(4) * t85 - t123 / 0.2e1) * qJD(4)) * t183 + (t86 / 0.2e1 + qJD(4) * t124 / 0.2e1 + t162 * t219 / 0.2e1 + t203 * mrSges(5,3) + (m(5) * t203 - qJD(4) * t154 + t113) * pkin(9)) * t187 + m(6) * (t107 * t17 + t108 * t19 + t114 * t171 + t156 * t213 + t58 * t81 + t59 * t82) + m(5) * (-pkin(3) * t176 + (-t116 * t218 - t70 * t183) * pkin(9)) + (-t15 * t49 + t16 * t50 + t2 * t94 - t3 * t95) * mrSges(7,3); -0.2e1 * pkin(3) * t149 - t134 * t100 - t135 * t101 + 0.2e1 * t111 * t54 + 0.2e1 * t117 * t21 - t196 * t89 + t143 * t90 + t187 * t151 + t183 * t152 + 0.2e1 * t171 * t88 + t94 * t22 + t95 * t23 + t49 * t56 + t50 * t55 + (t187 * t162 + (0.2e1 * pkin(4) * t99 - t161) * t183) * qJD(4) + (t107 * t81 + t108 * t82 + t171 * t213) * t245 + (t11 * t40 + t111 * t117 + t12 * t39) * t244 + 0.2e1 * (t11 * t94 - t12 * t95 - t39 * t49 + t40 * t50) * mrSges(7,3) + 0.2e1 * (t107 * t135 - t108 * t134 - t143 * t81 - t196 * t82) * mrSges(6,3); t52 * mrSges(5,1) + t18 * mrSges(6,1) - t53 * mrSges(5,2) - t20 * mrSges(6,2) + m(7) * (t121 * t30 - t122 * t29 + t130 * t6 + t131 * t5) + m(6) * (t178 * t20 + t18 * t180) * pkin(4) + t212; (m(6) * (t17 * t180 + t178 * t19) + t180 * t68 + t178 * t67) * pkin(4) + t130 * t24 + t131 * t25 + t121 * t65 - t122 * t66 - t69 * mrSges(5,2) + t70 * mrSges(5,1) + t17 * mrSges(6,1) - t19 * mrSges(6,2) - Ifges(5,5) * t208 - t191 * Ifges(5,6) - t193 + m(7) * (t121 * t16 - t122 * t15 + t130 * t3 + t131 * t2) + t248; m(7) * (t11 * t131 + t12 * t130 + t121 * t40 - t122 * t39) + t81 * mrSges(6,1) - t82 * mrSges(6,2) - t129 - t128 + t175 + (pkin(9) * t159 - t233) * qJD(4) + (m(6) * (t178 * t82 + t180 * t81) + (-t134 * t178 + t135 * t180) * mrSges(6,3)) * pkin(4) + (t121 * t94 + t122 * t95 - t130 * t49 + t131 * t50) * mrSges(7,3) + t197; 0.2e1 * m(7) * (t121 * t131 - t122 * t130) + 0.2e1 * t201; 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t102; m(6) * t114 + m(7) * t60 + t249; m(6) * t213 + m(7) * t111 + t247; 0; 0; t212; -t214 + t248; t197; t201; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
