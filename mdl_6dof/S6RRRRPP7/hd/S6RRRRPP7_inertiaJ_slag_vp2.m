% Calculate joint inertia matrix for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:17:12
% EndTime: 2019-03-09 21:17:18
% DurationCPUTime: 2.07s
% Computational Cost: add. (3568->454), mult. (8089->627), div. (0->0), fcn. (8762->10), ass. (0->157)
t219 = 2 * pkin(9);
t174 = cos(pkin(6));
t176 = sin(qJ(3));
t179 = cos(qJ(3));
t172 = sin(pkin(6));
t177 = sin(qJ(2));
t205 = t172 * t177;
t117 = -t174 * t179 + t176 * t205;
t171 = sin(pkin(11));
t173 = cos(pkin(11));
t118 = t174 * t176 + t179 * t205;
t175 = sin(qJ(4));
t178 = cos(qJ(4));
t180 = cos(qJ(2));
t204 = t172 * t180;
t79 = -t118 * t175 - t178 * t204;
t80 = t118 * t178 - t175 * t204;
t39 = t171 * t80 - t173 * t79;
t40 = t171 * t79 + t173 * t80;
t10 = Ifges(7,4) * t40 + Ifges(7,2) * t117 + Ifges(7,6) * t39;
t25 = Ifges(5,5) * t80 + Ifges(5,6) * t79 + Ifges(5,3) * t117;
t9 = Ifges(6,5) * t40 - Ifges(6,6) * t39 + Ifges(6,3) * t117;
t218 = t10 + t25 + t9;
t27 = Ifges(5,1) * t80 + Ifges(5,4) * t79 + Ifges(5,5) * t117;
t217 = t27 / 0.2e1;
t208 = Ifges(5,4) * t178;
t108 = -Ifges(5,6) * t179 + (-Ifges(5,2) * t175 + t208) * t176;
t216 = t108 / 0.2e1;
t209 = Ifges(5,4) * t175;
t109 = -Ifges(5,5) * t179 + (Ifges(5,1) * t178 - t209) * t176;
t215 = t109 / 0.2e1;
t143 = Ifges(5,1) * t175 + t208;
t214 = t143 / 0.2e1;
t147 = pkin(8) * t205;
t213 = pkin(1) * t180;
t104 = t147 + (-pkin(2) - t213) * t174;
t45 = pkin(3) * t117 - pkin(10) * t118 + t104;
t120 = t174 * t177 * pkin(1) + pkin(8) * t204;
t105 = pkin(9) * t174 + t120;
t106 = (-pkin(2) * t180 - pkin(9) * t177 - pkin(1)) * t172;
t52 = t179 * t105 + t176 * t106;
t47 = -pkin(10) * t204 + t52;
t19 = t175 * t45 + t178 * t47;
t15 = qJ(5) * t79 + t19;
t18 = -t175 * t47 + t178 * t45;
t7 = pkin(4) * t117 - qJ(5) * t80 + t18;
t4 = t173 * t15 + t171 * t7;
t212 = pkin(9) * t179;
t166 = t176 * pkin(9);
t210 = -qJ(5) - pkin(10);
t136 = -pkin(3) * t179 - pkin(10) * t176 - pkin(2);
t131 = t178 * t136;
t202 = t176 * t178;
t67 = -qJ(5) * t202 + t131 + (-pkin(9) * t175 - pkin(4)) * t179;
t203 = t175 * t176;
t96 = t175 * t136 + t178 * t212;
t81 = -qJ(5) * t203 + t96;
t31 = t171 * t67 + t173 * t81;
t119 = t174 * t213 - t147;
t207 = t119 * mrSges(3,1);
t206 = t120 * mrSges(3,2);
t128 = t171 * t175 - t173 * t178;
t111 = t128 * t176;
t91 = t179 * mrSges(7,1) - t111 * mrSges(7,2);
t129 = t171 * t178 + t173 * t175;
t110 = t129 * t176;
t201 = -Ifges(7,4) * t111 + Ifges(7,6) * t110;
t200 = -Ifges(6,5) * t111 - Ifges(6,6) * t110;
t199 = -Ifges(4,5) * t118 + Ifges(4,6) * t117;
t72 = Ifges(6,5) * t129 - Ifges(6,6) * t128;
t73 = Ifges(7,4) * t129 + Ifges(7,6) * t128;
t140 = Ifges(5,5) * t175 + Ifges(5,6) * t178;
t198 = Ifges(4,5) * t176 + Ifges(4,6) * t179;
t135 = pkin(4) * t203 + t166;
t197 = t175 ^ 2 + t178 ^ 2;
t139 = t210 * t178;
t187 = t210 * t175;
t83 = -t139 * t171 - t173 * t187;
t85 = -t173 * t139 + t171 * t187;
t196 = t83 ^ 2 + t85 ^ 2;
t195 = -Ifges(6,3) - Ifges(5,3) - Ifges(7,2);
t11 = Ifges(6,4) * t40 - Ifges(6,2) * t39 + Ifges(6,6) * t117;
t8 = Ifges(7,5) * t40 + Ifges(7,6) * t117 + Ifges(7,3) * t39;
t194 = t8 / 0.2e1 - t11 / 0.2e1;
t12 = Ifges(7,1) * t40 + Ifges(7,4) * t117 + Ifges(7,5) * t39;
t13 = Ifges(6,1) * t40 - Ifges(6,4) * t39 + Ifges(6,5) * t117;
t193 = t12 / 0.2e1 + t13 / 0.2e1;
t53 = -Ifges(7,5) * t111 - Ifges(7,6) * t179 + Ifges(7,3) * t110;
t56 = -Ifges(6,4) * t111 - Ifges(6,2) * t110 - Ifges(6,6) * t179;
t192 = t53 / 0.2e1 - t56 / 0.2e1;
t57 = -Ifges(7,1) * t111 - Ifges(7,4) * t179 + Ifges(7,5) * t110;
t58 = -Ifges(6,1) * t111 - Ifges(6,4) * t110 - Ifges(6,5) * t179;
t191 = t57 / 0.2e1 + t58 / 0.2e1;
t71 = Ifges(7,5) * t129 + Ifges(7,3) * t128;
t74 = Ifges(6,4) * t129 - Ifges(6,2) * t128;
t190 = t71 / 0.2e1 - t74 / 0.2e1;
t75 = Ifges(7,1) * t129 + Ifges(7,5) * t128;
t76 = Ifges(6,1) * t129 - Ifges(6,4) * t128;
t189 = t75 / 0.2e1 + t76 / 0.2e1;
t188 = Ifges(3,5) * t205 + Ifges(3,6) * t204 + Ifges(3,3) * t174;
t158 = -pkin(4) * t178 - pkin(3);
t17 = t39 * mrSges(6,1) + t40 * mrSges(6,2);
t16 = t39 * mrSges(7,1) - t40 * mrSges(7,3);
t62 = t110 * mrSges(6,1) - t111 * mrSges(6,2);
t23 = -t117 * mrSges(7,1) + t40 * mrSges(7,2);
t61 = t110 * mrSges(7,1) + t111 * mrSges(7,3);
t70 = t128 * mrSges(6,1) + t129 * mrSges(6,2);
t69 = t128 * mrSges(7,1) - t129 * mrSges(7,3);
t51 = -t176 * t105 + t106 * t179;
t186 = t72 / 0.2e1 + t73 / 0.2e1 + t140 / 0.2e1;
t185 = Ifges(5,5) * t202 - Ifges(5,6) * t203;
t46 = pkin(3) * t204 - t51;
t3 = -t15 * t171 + t173 * t7;
t184 = mrSges(5,1) * t175 + mrSges(5,2) * t178;
t30 = -t171 * t81 + t173 * t67;
t24 = -pkin(4) * t79 + t46;
t182 = pkin(9) ^ 2;
t170 = t179 ^ 2;
t168 = t176 ^ 2;
t165 = t168 * t182;
t157 = -pkin(4) * t173 - pkin(5);
t154 = pkin(4) * t171 + qJ(6);
t144 = Ifges(4,1) * t176 + Ifges(4,4) * t179;
t142 = Ifges(4,4) * t176 + Ifges(4,2) * t179;
t141 = Ifges(5,2) * t178 + t209;
t138 = -mrSges(4,1) * t179 + mrSges(4,2) * t176;
t137 = -mrSges(5,1) * t178 + mrSges(5,2) * t175;
t134 = -mrSges(5,1) * t179 - mrSges(5,3) * t202;
t133 = mrSges(5,2) * t179 - mrSges(5,3) * t203;
t127 = t184 * t176;
t107 = -Ifges(5,3) * t179 + t185;
t95 = -t175 * t212 + t131;
t90 = -mrSges(6,1) * t179 + mrSges(6,3) * t111;
t89 = mrSges(6,2) * t179 - mrSges(6,3) * t110;
t88 = -mrSges(7,2) * t110 - mrSges(7,3) * t179;
t87 = -mrSges(4,1) * t204 - mrSges(4,3) * t118;
t86 = mrSges(4,2) * t204 - mrSges(4,3) * t117;
t64 = pkin(5) * t128 - qJ(6) * t129 + t158;
t63 = mrSges(4,1) * t117 + mrSges(4,2) * t118;
t60 = Ifges(4,1) * t118 - Ifges(4,4) * t117 - Ifges(4,5) * t204;
t59 = Ifges(4,4) * t118 - Ifges(4,2) * t117 - Ifges(4,6) * t204;
t55 = -Ifges(7,2) * t179 + t201;
t54 = -Ifges(6,3) * t179 + t200;
t50 = mrSges(5,1) * t117 - mrSges(5,3) * t80;
t49 = -mrSges(5,2) * t117 + mrSges(5,3) * t79;
t48 = pkin(5) * t110 + qJ(6) * t111 + t135;
t41 = -mrSges(5,1) * t79 + mrSges(5,2) * t80;
t29 = pkin(5) * t179 - t30;
t28 = -qJ(6) * t179 + t31;
t26 = Ifges(5,4) * t80 + Ifges(5,2) * t79 + Ifges(5,6) * t117;
t22 = mrSges(6,1) * t117 - mrSges(6,3) * t40;
t21 = -mrSges(6,2) * t117 - mrSges(6,3) * t39;
t20 = -mrSges(7,2) * t39 + mrSges(7,3) * t117;
t5 = pkin(5) * t39 - qJ(6) * t40 + t24;
t2 = -pkin(5) * t117 - t3;
t1 = qJ(6) * t117 + t4;
t6 = [m(3) * (pkin(1) ^ 2 * t172 ^ 2 + t119 ^ 2 + t120 ^ 2) + m(4) * (t104 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2 + t46 ^ 2) + m(6) * (t24 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + 0.2e1 * t46 * t41 + 0.2e1 * t19 * t49 + 0.2e1 * t18 * t50 + 0.2e1 * t1 * t20 + 0.2e1 * t4 * t21 + 0.2e1 * t3 * t22 + 0.2e1 * t2 * t23 + 0.2e1 * t24 * t17 + 0.2e1 * t5 * t16 + (t12 + t13) * t40 + (t8 - t11) * t39 + (-t59 + t218) * t117 + (t188 - 0.2e1 * t206 + 0.2e1 * t207) * t174 + ((-0.2e1 * t119 * mrSges(3,3) + Ifges(3,5) * t174 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t177) * t172) * t177 + (0.2e1 * t120 * mrSges(3,3) + Ifges(3,6) * t174 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t177 + (Ifges(4,3) + Ifges(3,2)) * t180) * t172 + t199) * t180) * t172 + t79 * t26 + t80 * t27 + 0.2e1 * t52 * t86 + 0.2e1 * t51 * t87 + 0.2e1 * t104 * t63 + Ifges(2,3) + t118 * t60; t5 * t61 + t24 * t62 - pkin(2) * t63 + t48 * t16 + t29 * t23 + t30 * t22 + t31 * t21 + t28 * t20 + (t54 / 0.2e1 + t55 / 0.2e1 + t107 / 0.2e1 - t142 / 0.2e1) * t117 + (pkin(9) * t86 + t52 * mrSges(4,3) + t59 / 0.2e1 - t25 / 0.2e1 - t9 / 0.2e1 - t10 / 0.2e1) * t179 + m(6) * (t135 * t24 + t3 * t30 + t31 * t4) + m(7) * (t1 * t28 + t2 * t29 + t48 * t5) + t188 + m(4) * (-pkin(2) * t104 + (-t51 * t176 + t52 * t179) * pkin(9)) + t207 - t206 + m(5) * (t46 * t166 + t18 * t95 + t19 * t96) - t198 * t204 / 0.2e1 + (-t175 * t26 / 0.2e1 - t51 * mrSges(4,3) + t178 * t217 + t60 / 0.2e1 + (-t87 + t41) * pkin(9)) * t176 + t1 * t88 + t4 * t89 + t3 * t90 + t2 * t91 + t95 * t50 + t96 * t49 + t46 * t127 + t19 * t133 + t18 * t134 + t135 * t17 + t104 * t138 + t118 * t144 / 0.2e1 + t80 * t215 + t79 * t216 + t191 * t40 + t192 * t39 - t193 * t111 + t194 * t110; -0.2e1 * pkin(2) * t138 + 0.2e1 * t96 * t133 + 0.2e1 * t95 * t134 + 0.2e1 * t135 * t62 + 0.2e1 * t28 * t88 + 0.2e1 * t29 * t91 + 0.2e1 * t30 * t90 + 0.2e1 * t31 * t89 + 0.2e1 * t48 * t61 + Ifges(3,3) - (t57 + t58) * t111 + (t53 - t56) * t110 + (t168 + t170) * mrSges(4,3) * t219 + (-t54 - t55 - t107 + t142) * t179 + (-t108 * t175 + t109 * t178 + t127 * t219 + t144) * t176 + m(5) * (t95 ^ 2 + t96 ^ 2 + t165) + m(6) * (t135 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(7) * (t28 ^ 2 + t29 ^ 2 + t48 ^ 2) + m(4) * (pkin(2) ^ 2 + t170 * t182 + t165); t64 * t16 + t51 * mrSges(4,1) - t52 * mrSges(4,2) - pkin(3) * t41 + (t20 + t21) * t85 + (t23 - t22) * t83 + (pkin(10) * t49 + t19 * mrSges(5,3) + t26 / 0.2e1) * t178 + (-t18 * mrSges(5,3) - pkin(10) * t50 + t217) * t175 - t199 + m(5) * (-pkin(3) * t46 + (-t18 * t175 + t19 * t178) * pkin(10)) - Ifges(4,3) * t204 + t5 * t69 + t24 * t70 + t46 * t137 + t79 * t141 / 0.2e1 + t80 * t214 + t158 * t17 + t186 * t117 + t189 * t40 + t190 * t39 + (t2 * mrSges(7,2) - t3 * mrSges(6,3) + t193) * t129 + (-t1 * mrSges(7,2) - t4 * mrSges(6,3) + t194) * t128 + m(6) * (t158 * t24 - t3 * t83 + t4 * t85) + m(7) * (t1 * t85 + t2 * t83 + t5 * t64); -pkin(3) * t127 + t135 * t70 + t158 * t62 + t48 * t69 + t64 * t61 + (t88 + t89) * t85 + (-t90 + t91) * t83 - t189 * t111 + t190 * t110 + (-mrSges(4,1) + t137) * t166 + (t96 * mrSges(5,3) + pkin(10) * t133 + t176 * t214 + t216) * t178 + (-pkin(10) * t134 - t95 * mrSges(5,3) - t176 * t141 / 0.2e1 + t215) * t175 + m(6) * (t135 * t158 - t30 * t83 + t31 * t85) + m(7) * (t28 * t85 + t29 * t83 + t48 * t64) + m(5) * (-pkin(3) * t166 + (-t175 * t95 + t178 * t96) * pkin(10)) + (-pkin(9) * mrSges(4,2) - t186) * t179 + (t29 * mrSges(7,2) - t30 * mrSges(6,3) + t191) * t129 + (-t28 * mrSges(7,2) - t31 * mrSges(6,3) + t192) * t128 + t198; -0.2e1 * pkin(3) * t137 + t178 * t141 + t175 * t143 + 0.2e1 * t158 * t70 + 0.2e1 * t64 * t69 + Ifges(4,3) + m(7) * (t64 ^ 2 + t196) + m(6) * (t158 ^ 2 + t196) + m(5) * (t197 * pkin(10) ^ 2 + pkin(3) ^ 2) + (t75 + t76) * t129 + (t71 - t74) * t128 + 0.2e1 * t197 * pkin(10) * mrSges(5,3) + 0.2e1 * (-t128 * t85 + t129 * t83) * (mrSges(7,2) + mrSges(6,3)); -t4 * mrSges(6,2) - t2 * mrSges(7,1) + t157 * t23 + t1 * mrSges(7,3) + t154 * t20 + m(7) * (t1 * t154 + t157 * t2) + t3 * mrSges(6,1) - t19 * mrSges(5,2) + t18 * mrSges(5,1) + (m(6) * (t171 * t4 + t173 * t3) + t171 * t21 + t173 * t22) * pkin(4) + t218; m(7) * (t154 * t28 + t157 * t29) + t28 * mrSges(7,3) + t154 * t88 - t31 * mrSges(6,2) + t30 * mrSges(6,1) - t29 * mrSges(7,1) + t157 * t91 - t96 * mrSges(5,2) + t95 * mrSges(5,1) + t195 * t179 + (m(6) * (t171 * t31 + t173 * t30) + t171 * t89 + t173 * t90) * pkin(4) + t185 + t200 + t201; m(7) * (t154 * t85 + t157 * t83) - t85 * mrSges(6,2) - t83 * mrSges(6,1) - t83 * mrSges(7,1) + t85 * mrSges(7,3) - t184 * pkin(10) + (-t128 * t154 + t129 * t157) * mrSges(7,2) + (m(6) * (t171 * t85 - t173 * t83) + (-t128 * t171 - t129 * t173) * mrSges(6,3)) * pkin(4) + t140 + t73 + t72; -0.2e1 * t157 * mrSges(7,1) + 0.2e1 * t154 * mrSges(7,3) + m(7) * (t154 ^ 2 + t157 ^ 2) - t195 + (0.2e1 * mrSges(6,1) * t173 - 0.2e1 * mrSges(6,2) * t171 + m(6) * (t171 ^ 2 + t173 ^ 2) * pkin(4)) * pkin(4); m(6) * t24 + m(7) * t5 + t16 + t17; m(6) * t135 + m(7) * t48 + t61 + t62; m(6) * t158 + m(7) * t64 + t69 + t70; 0; m(6) + m(7); m(7) * t2 + t23; m(7) * t29 + t91; m(7) * t83 + t129 * mrSges(7,2); m(7) * t157 - mrSges(7,1); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
