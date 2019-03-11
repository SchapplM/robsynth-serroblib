% Calculate time derivative of joint inertia matrix for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:00:17
% EndTime: 2019-03-09 05:00:24
% DurationCPUTime: 3.56s
% Computational Cost: add. (5244->465), mult. (12262->685), div. (0->0), fcn. (10927->10), ass. (0->192)
t233 = -Ifges(5,3) - Ifges(6,3);
t161 = sin(pkin(11));
t162 = cos(pkin(11));
t165 = sin(qJ(4));
t168 = cos(qJ(4));
t130 = t161 * t168 + t162 * t165;
t174 = t161 * t165 - t162 * t168;
t166 = sin(qJ(3));
t195 = qJD(4) * t166;
t169 = cos(qJ(3));
t197 = qJD(3) * t169;
t69 = -t130 * t197 + t174 * t195;
t122 = t130 * qJD(4);
t70 = -t122 * t166 - t174 * t197;
t39 = -t69 * mrSges(6,1) + t70 * mrSges(6,2);
t164 = sin(qJ(6));
t167 = cos(qJ(6));
t110 = t130 * t166;
t111 = t174 * t166;
t64 = -t110 * t167 + t111 * t164;
t24 = qJD(6) * t64 + t164 * t69 + t167 * t70;
t65 = -t110 * t164 - t111 * t167;
t25 = -qJD(6) * t65 - t164 * t70 + t167 * t69;
t6 = -t25 * mrSges(7,1) + t24 * mrSges(7,2);
t173 = -t39 - t6;
t194 = qJD(4) * t168;
t170 = t165 * t197 + t166 * t194;
t186 = t165 * t195;
t187 = t168 * t197;
t171 = -t186 + t187;
t85 = t170 * mrSges(5,1) + t171 * mrSges(5,2);
t232 = -t85 + t173;
t190 = -cos(pkin(10)) * pkin(1) - pkin(2);
t231 = 0.2e1 * t190;
t127 = -pkin(3) * t169 - t166 * pkin(8) + t190;
t152 = sin(pkin(10)) * pkin(1) + pkin(7);
t202 = t168 * t169;
t139 = t152 * t202;
t193 = qJD(5) * t168;
t215 = pkin(8) * t169;
t217 = pkin(3) * t166;
t140 = (-t215 + t217) * qJD(3);
t198 = qJD(3) * t166;
t189 = t152 * t198;
t200 = t168 * t140 + t165 * t189;
t35 = -t166 * t193 + (pkin(4) * t166 - qJ(5) * t202) * qJD(3) + (-t139 + (qJ(5) * t166 - t127) * t165) * qJD(4) + t200;
t201 = t127 * t194 + t165 * t140;
t203 = t166 * t168;
t42 = (-qJ(5) * qJD(4) - qJD(3) * t152) * t203 + (-qJD(5) * t166 + (-qJ(5) * qJD(3) - qJD(4) * t152) * t169) * t165 + t201;
t14 = -t161 * t42 + t162 * t35;
t7 = pkin(5) * t198 - pkin(9) * t70 + t14;
t15 = t161 * t35 + t162 * t42;
t8 = pkin(9) * t69 + t15;
t113 = t168 * t127;
t71 = -qJ(5) * t203 + t113 + (-t152 * t165 - pkin(4)) * t169;
t204 = t165 * t166;
t92 = t165 * t127 + t139;
t80 = -qJ(5) * t204 + t92;
t40 = -t161 * t80 + t162 * t71;
t26 = -pkin(5) * t169 + t111 * pkin(9) + t40;
t41 = t161 * t71 + t162 * t80;
t27 = -pkin(9) * t110 + t41;
t9 = -t164 * t27 + t167 * t26;
t2 = qJD(6) * t9 + t164 * t7 + t167 * t8;
t10 = t164 * t26 + t167 * t27;
t3 = -qJD(6) * t10 - t164 * t8 + t167 * t7;
t230 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t123 = t174 * qJD(4);
t86 = -t130 * t164 - t167 * t174;
t46 = qJD(6) * t86 - t122 * t164 - t123 * t167;
t87 = t130 * t167 - t164 * t174;
t47 = -qJD(6) * t87 - t122 * t167 + t123 * t164;
t16 = -t47 * mrSges(7,1) + t46 * mrSges(7,2);
t82 = t122 * mrSges(6,1) - t123 * mrSges(6,2);
t229 = -t16 - t82;
t199 = t165 ^ 2 + t168 ^ 2;
t228 = 2 * m(5);
t227 = 2 * m(6);
t226 = 2 * m(7);
t225 = 0.2e1 * t152;
t224 = t86 / 0.2e1;
t223 = t87 / 0.2e1;
t221 = -t174 / 0.2e1;
t220 = t130 / 0.2e1;
t213 = Ifges(5,4) * t165;
t145 = Ifges(5,2) * t168 + t213;
t219 = -t145 / 0.2e1;
t218 = -t165 / 0.2e1;
t216 = pkin(4) * t161;
t214 = -qJ(5) - pkin(8);
t212 = Ifges(5,4) * t168;
t211 = Ifges(5,5) * t165;
t210 = Ifges(5,6) * t165;
t209 = Ifges(5,6) * t168;
t208 = t169 * Ifges(5,6);
t53 = -t92 * qJD(4) + t200;
t207 = t53 * t165;
t143 = -mrSges(5,1) * t168 + mrSges(5,2) * t165;
t206 = -mrSges(4,1) + t143;
t205 = t152 * t169;
t181 = qJD(4) * t214;
t119 = t165 * t181 + t193;
t120 = -qJD(5) * t165 + t168 * t181;
t73 = t162 * t119 + t161 * t120;
t142 = t214 * t165;
t144 = t214 * t168;
t95 = t161 * t142 - t162 * t144;
t121 = pkin(4) * t204 + t166 * t152;
t196 = qJD(4) * t165;
t192 = -Ifges(7,5) * t24 - Ifges(7,6) * t25 - Ifges(7,3) * t198;
t191 = pkin(4) * t196;
t141 = t152 * t197;
t93 = pkin(4) * t170 + t141;
t154 = -pkin(4) * t168 - pkin(3);
t188 = t166 * t197;
t183 = t169 * t196;
t182 = (2 * Ifges(4,4)) + t210;
t52 = (-t168 * t198 - t183) * t152 + t201;
t91 = -t165 * t205 + t113;
t180 = -t91 * qJD(4) + t52;
t72 = -t119 * t161 + t162 * t120;
t94 = t162 * t142 + t144 * t161;
t179 = mrSges(5,1) * t165 + mrSges(5,2) * t168;
t153 = pkin(4) * t162 + pkin(5);
t117 = t153 * t167 - t164 * t216;
t105 = t117 * qJD(6);
t118 = t153 * t164 + t167 * t216;
t106 = t118 * qJD(6);
t178 = -t106 * mrSges(7,1) - t105 * mrSges(7,2);
t177 = Ifges(5,1) * t168 - t213;
t146 = Ifges(5,1) * t165 + t212;
t176 = -Ifges(5,2) * t165 + t212;
t74 = -pkin(9) * t130 + t94;
t75 = -pkin(9) * t174 + t95;
t36 = -t164 * t75 + t167 * t74;
t37 = t164 * t74 + t167 * t75;
t54 = pkin(9) * t123 + t72;
t55 = -pkin(9) * t122 + t73;
t12 = qJD(6) * t36 + t164 * t54 + t167 * t55;
t13 = -qJD(6) * t37 - t164 * t55 + t167 * t54;
t44 = Ifges(7,6) * t47;
t45 = Ifges(7,5) * t46;
t175 = t13 * mrSges(7,1) - t12 * mrSges(7,2) + t44 + t45;
t172 = -Ifges(5,5) * t187 - Ifges(6,5) * t70 - Ifges(6,6) * t69 + t198 * t233 + t192;
t158 = Ifges(5,5) * t194;
t138 = -mrSges(5,1) * t169 - mrSges(5,3) * t203;
t137 = mrSges(5,2) * t169 - mrSges(5,3) * t204;
t136 = t177 * qJD(4);
t135 = t176 * qJD(4);
t134 = t179 * qJD(4);
t126 = t179 * t166;
t116 = Ifges(6,5) * t123;
t115 = Ifges(6,6) * t122;
t109 = -Ifges(5,5) * t169 + t166 * t177;
t108 = t166 * t176 - t208;
t101 = pkin(5) * t174 + t154;
t100 = -mrSges(5,2) * t198 - mrSges(5,3) * t170;
t99 = mrSges(5,1) * t198 - mrSges(5,3) * t171;
t98 = pkin(5) * t122 + t191;
t97 = -mrSges(6,1) * t169 + t111 * mrSges(6,3);
t96 = mrSges(6,2) * t169 - t110 * mrSges(6,3);
t90 = Ifges(6,1) * t130 - Ifges(6,4) * t174;
t89 = Ifges(6,4) * t130 - Ifges(6,2) * t174;
t88 = mrSges(6,1) * t174 + mrSges(6,2) * t130;
t84 = -Ifges(6,1) * t123 - Ifges(6,4) * t122;
t83 = -Ifges(6,4) * t123 - Ifges(6,2) * t122;
t81 = pkin(5) * t110 + t121;
t79 = -t146 * t195 + (Ifges(5,5) * t166 + t169 * t177) * qJD(3);
t78 = -t145 * t195 + (Ifges(5,6) * t166 + t169 * t176) * qJD(3);
t77 = mrSges(6,1) * t110 - mrSges(6,2) * t111;
t61 = -Ifges(6,1) * t111 - Ifges(6,4) * t110 - Ifges(6,5) * t169;
t60 = -Ifges(6,4) * t111 - Ifges(6,2) * t110 - Ifges(6,6) * t169;
t59 = mrSges(6,1) * t198 - mrSges(6,3) * t70;
t58 = -mrSges(6,2) * t198 + mrSges(6,3) * t69;
t57 = -mrSges(7,1) * t169 - t65 * mrSges(7,3);
t56 = mrSges(7,2) * t169 + t64 * mrSges(7,3);
t51 = -pkin(5) * t69 + t93;
t50 = Ifges(7,1) * t87 + Ifges(7,4) * t86;
t49 = Ifges(7,4) * t87 + Ifges(7,2) * t86;
t48 = -mrSges(7,1) * t86 + mrSges(7,2) * t87;
t34 = -mrSges(7,1) * t64 + mrSges(7,2) * t65;
t33 = Ifges(6,1) * t70 + Ifges(6,4) * t69 + Ifges(6,5) * t198;
t32 = Ifges(6,4) * t70 + Ifges(6,2) * t69 + Ifges(6,6) * t198;
t29 = Ifges(7,1) * t65 + Ifges(7,4) * t64 - Ifges(7,5) * t169;
t28 = Ifges(7,4) * t65 + Ifges(7,2) * t64 - Ifges(7,6) * t169;
t20 = -mrSges(7,2) * t198 + mrSges(7,3) * t25;
t19 = mrSges(7,1) * t198 - mrSges(7,3) * t24;
t18 = Ifges(7,1) * t46 + Ifges(7,4) * t47;
t17 = Ifges(7,4) * t46 + Ifges(7,2) * t47;
t5 = Ifges(7,1) * t24 + Ifges(7,4) * t25 + Ifges(7,5) * t198;
t4 = Ifges(7,4) * t24 + Ifges(7,2) * t25 + Ifges(7,6) * t198;
t1 = [(t85 * t225 - t165 * t78 + t168 * t79 + (-t168 * t108 - t165 * t109 - t169 * (-t209 - t211)) * qJD(4) + (mrSges(4,1) * t231 - Ifges(6,5) * t111 - Ifges(6,6) * t110 + Ifges(7,5) * t65 + Ifges(7,6) * t64 + (Ifges(5,5) * t168 - t182) * t166 + (t152 ^ 2 * t228 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(7,3) + t233) * t169) * qJD(3)) * t166 + (t10 * t2 + t3 * t9 + t51 * t81) * t226 + (t121 * t93 + t14 * t40 + t15 * t41) * t227 + (t92 * t52 + t91 * t53) * t228 + ((mrSges(4,2) * t231 - t165 * t108 + t168 * t109 + t126 * t225 + t169 * t182) * qJD(3) + t172) * t169 + 0.2e1 * t52 * t137 + 0.2e1 * t53 * t138 - t111 * t33 + 0.2e1 * t121 * t39 - t110 * t32 + 0.2e1 * t93 * t77 + 0.2e1 * t15 * t96 + 0.2e1 * t14 * t97 + 0.2e1 * t91 * t99 + 0.2e1 * t92 * t100 + 0.2e1 * t81 * t6 + t69 * t60 + t70 * t61 + t64 * t4 + t65 * t5 + 0.2e1 * t2 * t56 + 0.2e1 * t3 * t57 + 0.2e1 * t41 * t58 + 0.2e1 * t40 * t59 + 0.2e1 * t51 * t34 + t25 * t28 + t24 * t29 + 0.2e1 * t9 * t19 + 0.2e1 * t10 * t20; -t110 * t59 - t111 * t58 + t64 * t19 + t65 * t20 + t24 * t56 + t25 * t57 + t69 * t97 + t70 * t96 + ((t137 * t168 - t138 * t165) * qJD(3) + t232) * t169 + (t168 * t100 - t165 * t99 + (-t137 * t165 - t138 * t168) * qJD(4) + (t126 + t34 + t77) * qJD(3)) * t166 + m(6) * (-t14 * t110 - t15 * t111 + t121 * t198 - t169 * t93 + t40 * t69 + t41 * t70) + m(7) * (t10 * t24 - t169 * t51 + t198 * t81 + t2 * t65 + t9 * t25 + t3 * t64) + m(5) * ((-t165 * t91 + t168 * t92 - t205) * t197 + (t189 - t207 + t168 * t52 + (-t165 * t92 - t168 * t91) * qJD(4)) * t166); 0.2e1 * m(6) * (-t110 * t69 - t111 * t70 - t188) + 0.2e1 * m(7) * (t65 * t24 + t64 * t25 - t188) + 0.2e1 * m(5) * (-0.1e1 + t199) * t188; (-t122 * t41 + t123 * t40 - t130 * t14 - t15 * t174) * mrSges(6,3) + (-t158 / 0.2e1 + t116 / 0.2e1 + t115 / 0.2e1 - t45 / 0.2e1 - t44 / 0.2e1 + (t152 * t206 + Ifges(4,5)) * qJD(3)) * t169 + (t10 * t47 + t2 * t86 - t3 * t87 - t46 * t9) * mrSges(7,3) + t5 * t223 + t4 * t224 + t154 * t39 - t122 * t60 / 0.2e1 - t123 * t61 / 0.2e1 - t111 * t84 / 0.2e1 + t121 * t82 - t110 * t83 / 0.2e1 + t69 * t89 / 0.2e1 + t70 * t90 / 0.2e1 + t93 * t88 + t94 * t59 + t95 * t58 + t73 * t96 + t72 * t97 + t98 * t34 + t101 * t6 + t81 * t16 - pkin(3) * t85 + t64 * t17 / 0.2e1 + t65 * t18 / 0.2e1 + t12 * t56 + t13 * t57 + t46 * t29 / 0.2e1 + t47 * t28 / 0.2e1 + t25 * t49 / 0.2e1 + t24 * t50 / 0.2e1 + t51 * t48 + t37 * t20 + t36 * t19 + t33 * t220 + t32 * t221 + m(7) * (t10 * t12 + t101 * t51 + t13 * t9 + t2 * t37 + t3 * t36 + t81 * t98) + m(5) * (-pkin(3) * t141 + (-t196 * t92 - t207) * pkin(8)) + m(6) * (t121 * t191 + t14 * t94 + t15 * t95 + t154 * t93 + t40 * t72 + t41 * t73) + (t146 * t197 / 0.2e1 + t78 / 0.2e1 + qJD(4) * t109 / 0.2e1 + t180 * mrSges(5,3) + (m(5) * t180 - qJD(4) * t138 + t100) * pkin(8)) * t168 + (t197 * t219 + t79 / 0.2e1 - pkin(8) * t99 - t53 * mrSges(5,3) + (-pkin(8) * t137 + pkin(4) * t77 - t92 * mrSges(5,3) + t208 / 0.2e1 - t108 / 0.2e1) * qJD(4)) * t165 + (t168 * t136 / 0.2e1 + t135 * t218 + t152 * t134 + (t146 * t218 + t168 * t219) * qJD(4) + (t152 * mrSges(4,2) + t211 / 0.2e1 + t209 / 0.2e1 + Ifges(6,5) * t220 + Ifges(6,6) * t221 + Ifges(7,5) * t223 + Ifges(7,6) * t224 - Ifges(4,6)) * qJD(3)) * t166; (-t134 + t229) * t169 + m(6) * (-pkin(4) * t183 - t72 * t110 - t73 * t111 + t94 * t69 + t95 * t70) + m(7) * (t12 * t65 + t13 * t64 - t169 * t98 + t37 * t24 + t36 * t25) + (t24 * t86 - t25 * t87 - t46 * t64 + t47 * t65) * mrSges(7,3) + (-t110 * t123 + t111 * t122 - t130 * t69 - t174 * t70) * mrSges(6,3) + ((mrSges(5,3) * t199 - mrSges(4,2)) * t169 + m(5) * (t199 * t215 - t217) + (m(6) * t154 + m(7) * t101 + t206 + t48 + t88) * t166) * qJD(3); -0.2e1 * pkin(3) * t134 + 0.2e1 * t101 * t16 - t122 * t89 - t123 * t90 - t174 * t83 + t130 * t84 + t168 * t135 + t165 * t136 + 0.2e1 * t154 * t82 + t86 * t17 + t87 * t18 + t46 * t50 + t47 * t49 + 0.2e1 * t98 * t48 + (t168 * t146 + (0.2e1 * pkin(4) * t88 - t145) * t165) * qJD(4) + (t154 * t191 + t72 * t94 + t73 * t95) * t227 + (t101 * t98 + t12 * t37 + t13 * t36) * t226 + 0.2e1 * (t12 * t86 - t13 * t87 - t36 * t46 + t37 * t47) * mrSges(7,3) + 0.2e1 * (-t122 * t95 + t123 * t94 - t130 * t72 - t174 * t73) * mrSges(6,3); -t170 * Ifges(5,6) - Ifges(5,5) * t186 - t172 + (t161 * t58 + m(6) * (t14 * t162 + t15 * t161) + t162 * t59) * pkin(4) + t117 * t19 + t118 * t20 + t105 * t56 - t106 * t57 - t52 * mrSges(5,2) + t53 * mrSges(5,1) + t14 * mrSges(6,1) - t15 * mrSges(6,2) + m(7) * (t10 * t105 - t106 * t9 + t117 * t3 + t118 * t2) + t230; m(7) * (t105 * t65 - t106 * t64 + t117 * t25 + t118 * t24) + m(6) * (t161 * t70 + t162 * t69) * pkin(4) + t232; m(7) * (t105 * t37 - t106 * t36 + t117 * t13 + t118 * t12) - t73 * mrSges(6,2) + t72 * mrSges(6,1) - t116 - t115 + t158 + (pkin(8) * t143 - t210) * qJD(4) + (m(6) * (t161 * t73 + t162 * t72) + (-t122 * t161 + t123 * t162) * mrSges(6,3)) * pkin(4) + (t105 * t86 + t106 * t87 - t117 * t46 + t118 * t47) * mrSges(7,3) + t175; 0.2e1 * m(7) * (t105 * t118 - t106 * t117) + 0.2e1 * t178; m(6) * t93 + m(7) * t51 - t173; (m(6) + m(7)) * t198; m(6) * t191 + m(7) * t98 - t229; 0; 0; -t192 + t230; -t6; t175; t178; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
