% Calculate joint inertia matrix for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:34
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP10_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP10_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP10_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:33:53
% EndTime: 2018-11-23 18:33:55
% DurationCPUTime: 2.07s
% Computational Cost: add. (3888->472), mult. (8716->648), div. (0->0), fcn. (9422->10), ass. (0->163)
t225 = 2 * pkin(9);
t172 = cos(pkin(6));
t175 = sin(qJ(3));
t179 = cos(qJ(3));
t171 = sin(pkin(6));
t176 = sin(qJ(2));
t209 = t171 * t176;
t119 = -t172 * t179 + t175 * t209;
t173 = sin(qJ(5));
t177 = cos(qJ(5));
t120 = t172 * t175 + t179 * t209;
t174 = sin(qJ(4));
t178 = cos(qJ(4));
t180 = cos(qJ(2));
t208 = t171 * t180;
t81 = -t120 * t174 - t178 * t208;
t82 = t120 * t178 - t174 * t208;
t41 = t173 * t82 - t177 * t81;
t42 = t173 * t81 + t177 * t82;
t11 = Ifges(6,5) * t42 - Ifges(6,6) * t41 + Ifges(6,3) * t119;
t12 = Ifges(7,4) * t42 + Ifges(7,2) * t119 + Ifges(7,6) * t41;
t224 = t11 + t12;
t223 = 2 * mrSges(7,3);
t29 = Ifges(5,1) * t82 + Ifges(5,4) * t81 + Ifges(5,5) * t119;
t222 = t29 / 0.2e1;
t221 = -pkin(11) - pkin(10);
t212 = Ifges(5,4) * t178;
t110 = -Ifges(5,6) * t179 + (-Ifges(5,2) * t174 + t212) * t175;
t220 = t110 / 0.2e1;
t213 = Ifges(5,4) * t174;
t111 = -Ifges(5,5) * t179 + (Ifges(5,1) * t178 - t213) * t175;
t219 = t111 / 0.2e1;
t142 = Ifges(5,1) * t174 + t212;
t218 = t142 / 0.2e1;
t147 = pkin(8) * t209;
t217 = pkin(1) * t180;
t101 = t147 + (-pkin(2) - t217) * t172;
t47 = pkin(3) * t119 - pkin(10) * t120 + t101;
t122 = t172 * t176 * pkin(1) + pkin(8) * t208;
t102 = pkin(9) * t172 + t122;
t103 = (-pkin(2) * t180 - pkin(9) * t176 - pkin(1)) * t171;
t54 = t179 * t102 + t175 * t103;
t49 = -pkin(10) * t208 + t54;
t21 = t174 * t47 + t178 * t49;
t17 = pkin(11) * t81 + t21;
t20 = -t174 * t49 + t178 * t47;
t9 = pkin(4) * t119 - pkin(11) * t82 + t20;
t6 = t177 * t17 + t173 * t9;
t216 = pkin(9) * t179;
t166 = t175 * pkin(9);
t214 = -Ifges(6,3) - Ifges(7,2);
t136 = -pkin(3) * t179 - pkin(10) * t175 - pkin(2);
t129 = t178 * t136;
t206 = t175 * t178;
t69 = -pkin(11) * t206 + t129 + (-pkin(9) * t174 - pkin(4)) * t179;
t100 = t174 * t136 + t178 * t216;
t207 = t174 * t175;
t83 = -pkin(11) * t207 + t100;
t35 = t173 * t69 + t177 * t83;
t121 = t172 * t217 - t147;
t211 = t121 * mrSges(3,1);
t210 = t122 * mrSges(3,2);
t130 = t173 * t174 - t177 * t178;
t113 = t130 * t175;
t95 = t179 * mrSges(7,1) - t113 * mrSges(7,2);
t131 = t173 * t178 + t174 * t177;
t112 = t131 * t175;
t205 = -Ifges(7,4) * t113 + Ifges(7,6) * t112;
t204 = -Ifges(6,5) * t113 - Ifges(6,6) * t112;
t203 = -Ifges(4,5) * t120 + Ifges(4,6) * t119;
t76 = Ifges(6,5) * t131 - Ifges(6,6) * t130;
t77 = Ifges(7,4) * t131 + Ifges(7,6) * t130;
t139 = Ifges(5,5) * t174 + Ifges(5,6) * t178;
t202 = Ifges(4,5) * t175 + Ifges(4,6) * t179;
t135 = pkin(4) * t207 + t166;
t201 = t174 ^ 2 + t178 ^ 2;
t144 = t221 * t178;
t192 = t221 * t174;
t89 = -t144 * t173 - t177 * t192;
t91 = -t177 * t144 + t173 * t192;
t200 = t89 ^ 2 + t91 ^ 2;
t199 = -Ifges(5,3) + t214;
t27 = Ifges(5,5) * t82 + Ifges(5,6) * t81 + Ifges(5,3) * t119;
t10 = Ifges(7,5) * t42 + Ifges(7,6) * t119 + Ifges(7,3) * t41;
t13 = Ifges(6,4) * t42 - Ifges(6,2) * t41 + Ifges(6,6) * t119;
t198 = -t13 / 0.2e1 + t10 / 0.2e1;
t14 = Ifges(7,1) * t42 + Ifges(7,4) * t119 + Ifges(7,5) * t41;
t15 = Ifges(6,1) * t42 - Ifges(6,4) * t41 + Ifges(6,5) * t119;
t197 = t14 / 0.2e1 + t15 / 0.2e1;
t55 = -Ifges(7,5) * t113 - Ifges(7,6) * t179 + Ifges(7,3) * t112;
t58 = -Ifges(6,4) * t113 - Ifges(6,2) * t112 - Ifges(6,6) * t179;
t196 = t55 / 0.2e1 - t58 / 0.2e1;
t59 = -Ifges(7,1) * t113 - Ifges(7,4) * t179 + Ifges(7,5) * t112;
t60 = -Ifges(6,1) * t113 - Ifges(6,4) * t112 - Ifges(6,5) * t179;
t195 = t60 / 0.2e1 + t59 / 0.2e1;
t75 = Ifges(7,5) * t131 + Ifges(7,3) * t130;
t78 = Ifges(6,4) * t131 - Ifges(6,2) * t130;
t194 = -t78 / 0.2e1 + t75 / 0.2e1;
t79 = Ifges(7,1) * t131 + Ifges(7,5) * t130;
t80 = Ifges(6,1) * t131 - Ifges(6,4) * t130;
t193 = t79 / 0.2e1 + t80 / 0.2e1;
t191 = Ifges(3,5) * t209 + Ifges(3,6) * t208 + Ifges(3,3) * t172;
t158 = -pkin(4) * t178 - pkin(3);
t25 = -t119 * mrSges(7,1) + t42 * mrSges(7,2);
t53 = -t175 * t102 + t103 * t179;
t190 = t76 / 0.2e1 + t77 / 0.2e1 + t139 / 0.2e1;
t189 = Ifges(5,5) * t206 - Ifges(5,6) * t207;
t48 = pkin(3) * t208 - t53;
t5 = -t17 * t173 + t177 * t9;
t188 = mrSges(5,1) * t174 + mrSges(5,2) * t178;
t34 = -t173 * t83 + t177 * t69;
t187 = (mrSges(6,1) * t177 - mrSges(6,2) * t173) * pkin(4);
t26 = -pkin(4) * t81 + t48;
t31 = -qJ(6) * t179 + t35;
t32 = pkin(5) * t179 - t34;
t186 = t34 * mrSges(6,1) - t32 * mrSges(7,1) - t35 * mrSges(6,2) + t31 * mrSges(7,3) + t204 + t205;
t185 = t76 + t77 + (-mrSges(6,2) + mrSges(7,3)) * t91 + (-mrSges(6,1) - mrSges(7,1)) * t89;
t2 = qJ(6) * t119 + t6;
t3 = -pkin(5) * t119 - t5;
t184 = t5 * mrSges(6,1) - t3 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3) + t224;
t182 = pkin(9) ^ 2;
t170 = t179 ^ 2;
t168 = t175 ^ 2;
t165 = t168 * t182;
t157 = -pkin(4) * t177 - pkin(5);
t155 = pkin(4) * t173 + qJ(6);
t143 = Ifges(4,1) * t175 + Ifges(4,4) * t179;
t141 = Ifges(4,4) * t175 + Ifges(4,2) * t179;
t140 = Ifges(5,2) * t178 + t213;
t138 = -mrSges(4,1) * t179 + mrSges(4,2) * t175;
t137 = -mrSges(5,1) * t178 + mrSges(5,2) * t174;
t134 = -mrSges(5,1) * t179 - mrSges(5,3) * t206;
t133 = mrSges(5,2) * t179 - mrSges(5,3) * t207;
t123 = t188 * t175;
t109 = -Ifges(5,3) * t179 + t189;
t99 = -t174 * t216 + t129;
t94 = -mrSges(6,1) * t179 + mrSges(6,3) * t113;
t93 = mrSges(6,2) * t179 - mrSges(6,3) * t112;
t92 = -mrSges(7,2) * t112 - mrSges(7,3) * t179;
t88 = -mrSges(4,1) * t208 - mrSges(4,3) * t120;
t87 = mrSges(4,2) * t208 - mrSges(4,3) * t119;
t74 = mrSges(6,1) * t130 + mrSges(6,2) * t131;
t73 = mrSges(7,1) * t130 - mrSges(7,3) * t131;
t67 = pkin(5) * t130 - qJ(6) * t131 + t158;
t65 = mrSges(4,1) * t119 + mrSges(4,2) * t120;
t64 = mrSges(6,1) * t112 - mrSges(6,2) * t113;
t63 = mrSges(7,1) * t112 + mrSges(7,3) * t113;
t62 = Ifges(4,1) * t120 - Ifges(4,4) * t119 - Ifges(4,5) * t208;
t61 = Ifges(4,4) * t120 - Ifges(4,2) * t119 - Ifges(4,6) * t208;
t57 = -Ifges(7,2) * t179 + t205;
t56 = -Ifges(6,3) * t179 + t204;
t52 = mrSges(5,1) * t119 - mrSges(5,3) * t82;
t51 = -mrSges(5,2) * t119 + mrSges(5,3) * t81;
t50 = pkin(5) * t112 + qJ(6) * t113 + t135;
t43 = -mrSges(5,1) * t81 + mrSges(5,2) * t82;
t28 = Ifges(5,4) * t82 + Ifges(5,2) * t81 + Ifges(5,6) * t119;
t24 = mrSges(6,1) * t119 - mrSges(6,3) * t42;
t23 = -mrSges(6,2) * t119 - mrSges(6,3) * t41;
t22 = -mrSges(7,2) * t41 + mrSges(7,3) * t119;
t19 = mrSges(6,1) * t41 + mrSges(6,2) * t42;
t18 = mrSges(7,1) * t41 - mrSges(7,3) * t42;
t7 = pkin(5) * t41 - qJ(6) * t42 + t26;
t1 = [0.2e1 * t101 * t65 + 0.2e1 * t54 * t87 + 0.2e1 * t53 * t88 + t81 * t28 + t82 * t29 + 0.2e1 * t48 * t43 + 0.2e1 * t21 * t51 + 0.2e1 * t20 * t52 + 0.2e1 * t2 * t22 + 0.2e1 * t6 * t23 + 0.2e1 * t5 * t24 + 0.2e1 * t3 * t25 + 0.2e1 * t26 * t19 + 0.2e1 * t7 * t18 + ((-0.2e1 * t121 * mrSges(3,3) + Ifges(3,5) * t172 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t176) * t171) * t176 + (0.2e1 * t122 * mrSges(3,3) + Ifges(3,6) * t172 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t176 + (Ifges(4,3) + Ifges(3,2)) * t180) * t171 + t203) * t180) * t171 + m(6) * (t26 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t7 ^ 2) + m(5) * (t20 ^ 2 + t21 ^ 2 + t48 ^ 2) + m(4) * (t101 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(3) * (pkin(1) ^ 2 * t171 ^ 2 + t121 ^ 2 + t122 ^ 2) + (t14 + t15) * t42 + (-t13 + t10) * t41 + (t191 - 0.2e1 * t210 + 0.2e1 * t211) * t172 + t120 * t62 + Ifges(2,3) + (t27 - t61 + t224) * t119; t99 * t52 + t100 * t51 + t2 * t92 + t6 * t93 + t5 * t94 + t3 * t95 + t7 * t63 + t26 * t64 - pkin(2) * t65 + t50 * t18 + t32 * t25 + t34 * t24 + t35 * t23 + t31 * t22 + t191 + m(4) * (-pkin(2) * t101 + (-t53 * t175 + t54 * t179) * pkin(9)) + t82 * t219 + t81 * t220 + (t109 / 0.2e1 + t56 / 0.2e1 + t57 / 0.2e1 - t141 / 0.2e1) * t119 - t202 * t208 / 0.2e1 - t197 * t113 + t198 * t112 + t195 * t42 + t196 * t41 - t210 + t211 + (t54 * mrSges(4,3) + pkin(9) * t87 + t61 / 0.2e1 - t27 / 0.2e1 - t11 / 0.2e1 - t12 / 0.2e1) * t179 + t48 * t123 + t21 * t133 + t20 * t134 + t135 * t19 + t101 * t138 + t120 * t143 / 0.2e1 + m(5) * (t100 * t21 + t166 * t48 + t20 * t99) + (-t53 * mrSges(4,3) - t174 * t28 / 0.2e1 + t178 * t222 + t62 / 0.2e1 + (t43 - t88) * pkin(9)) * t175 + m(7) * (t2 * t31 + t3 * t32 + t50 * t7) + m(6) * (t135 * t26 + t34 * t5 + t35 * t6); -0.2e1 * pkin(2) * t138 + 0.2e1 * t100 * t133 + 0.2e1 * t99 * t134 + 0.2e1 * t135 * t64 + 0.2e1 * t31 * t92 + 0.2e1 * t32 * t95 + 0.2e1 * t34 * t94 + 0.2e1 * t35 * t93 + 0.2e1 * t50 * t63 + Ifges(3,3) - (t59 + t60) * t113 + (-t58 + t55) * t112 + (t168 + t170) * mrSges(4,3) * t225 + (-t109 - t56 - t57 + t141) * t179 + (-t110 * t174 + t111 * t178 + t123 * t225 + t143) * t175 + m(7) * (t31 ^ 2 + t32 ^ 2 + t50 ^ 2) + m(6) * (t135 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(5) * (t100 ^ 2 + t99 ^ 2 + t165) + m(4) * (pkin(2) ^ 2 + t170 * t182 + t165); t7 * t73 + t26 * t74 + t67 * t18 + t53 * mrSges(4,1) - t54 * mrSges(4,2) - pkin(3) * t43 + m(5) * (-pkin(3) * t48 + (-t20 * t174 + t21 * t178) * pkin(10)) + (t22 + t23) * t91 + (t25 - t24) * t89 + (t21 * mrSges(5,3) + pkin(10) * t51 + t28 / 0.2e1) * t178 + (-t20 * mrSges(5,3) - pkin(10) * t52 + t222) * t174 + m(7) * (t2 * t91 + t3 * t89 + t67 * t7) + m(6) * (t158 * t26 - t5 * t89 + t6 * t91) + (t3 * mrSges(7,2) - t5 * mrSges(6,3) + t197) * t131 + (-t2 * mrSges(7,2) - t6 * mrSges(6,3) + t198) * t130 + t193 * t42 + t194 * t41 + t190 * t119 - t203 + t48 * t137 + t81 * t140 / 0.2e1 + t82 * t218 + t158 * t19 - Ifges(4,3) * t208; -pkin(3) * t123 + t135 * t74 + t158 * t64 + t50 * t73 + t67 * t63 + (t92 + t93) * t91 + (-t94 + t95) * t89 - t193 * t113 + t194 * t112 + (-mrSges(4,1) + t137) * t166 + (t100 * mrSges(5,3) + pkin(10) * t133 + t175 * t218 + t220) * t178 + (-t99 * mrSges(5,3) - pkin(10) * t134 - t175 * t140 / 0.2e1 + t219) * t174 + m(7) * (t31 * t91 + t32 * t89 + t50 * t67) + m(6) * (t135 * t158 - t34 * t89 + t35 * t91) + m(5) * (-pkin(3) * t166 + (t100 * t178 - t99 * t174) * pkin(10)) + (-pkin(9) * mrSges(4,2) - t190) * t179 + (t32 * mrSges(7,2) - t34 * mrSges(6,3) + t195) * t131 + (-t31 * mrSges(7,2) - t35 * mrSges(6,3) + t196) * t130 + t202; -0.2e1 * pkin(3) * t137 + t178 * t140 + t174 * t142 + 0.2e1 * t158 * t74 + 0.2e1 * t67 * t73 + Ifges(4,3) + m(7) * (t67 ^ 2 + t200) + m(6) * (t158 ^ 2 + t200) + m(5) * (pkin(10) ^ 2 * t201 + pkin(3) ^ 2) + (t79 + t80) * t131 + (t75 - t78) * t130 + 0.2e1 * t201 * pkin(10) * mrSges(5,3) + 0.2e1 * (-t130 * t91 + t131 * t89) * (mrSges(7,2) + mrSges(6,3)); t20 * mrSges(5,1) - t21 * mrSges(5,2) + m(7) * (t155 * t2 + t157 * t3) + t184 + t155 * t22 + t157 * t25 + (m(6) * (t173 * t6 + t177 * t5) + t173 * t23 + t177 * t24) * pkin(4) + t27; t99 * mrSges(5,1) - t100 * mrSges(5,2) + m(7) * (t155 * t31 + t157 * t32) + t199 * t179 + t186 + t155 * t92 + t157 * t95 + (m(6) * (t173 * t35 + t177 * t34) + t173 * t93 + t177 * t94) * pkin(4) + t189; m(7) * (t155 * t91 + t157 * t89) - t188 * pkin(10) + (-t130 * t155 + t131 * t157) * mrSges(7,2) + (m(6) * (t173 * t91 - t177 * t89) + (-t130 * t173 - t131 * t177) * mrSges(6,3)) * pkin(4) + t185 + t139; -0.2e1 * t157 * mrSges(7,1) + t155 * t223 + 0.2e1 * t187 + m(7) * (t155 ^ 2 + t157 ^ 2) + m(6) * (t173 ^ 2 + t177 ^ 2) * pkin(4) ^ 2 - t199; qJ(6) * t22 + m(7) * (-pkin(5) * t3 + qJ(6) * t2) - pkin(5) * t25 + t184; m(7) * (-pkin(5) * t32 + qJ(6) * t31) - pkin(5) * t95 + qJ(6) * t92 + t214 * t179 + t186; m(7) * (-pkin(5) * t89 + qJ(6) * t91) + (-pkin(5) * t131 - qJ(6) * t130) * mrSges(7,2) + t185; m(7) * (-pkin(5) * t157 + qJ(6) * t155) + t187 + (t155 + qJ(6)) * mrSges(7,3) + (-t157 + pkin(5)) * mrSges(7,1) - t214; 0.2e1 * pkin(5) * mrSges(7,1) + qJ(6) * t223 + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) - t214; m(7) * t3 + t25; m(7) * t32 + t95; m(7) * t89 + t131 * mrSges(7,2); m(7) * t157 - mrSges(7,1); -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
