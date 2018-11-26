% Calculate joint inertia matrix for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2018-11-23 17:30
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR13_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR13_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:29:55
% EndTime: 2018-11-23 17:29:57
% DurationCPUTime: 1.88s
% Computational Cost: add. (2810->406), mult. (6116->575), div. (0->0), fcn. (6379->10), ass. (0->152)
t195 = Ifges(4,1) + Ifges(3,3);
t148 = (-pkin(2) - pkin(9));
t194 = -2 * t148;
t141 = sin(qJ(5));
t145 = cos(qJ(5));
t138 = sin(pkin(6));
t143 = sin(qJ(2));
t175 = t138 * t143;
t139 = cos(pkin(6));
t142 = sin(qJ(4));
t146 = cos(qJ(4));
t147 = cos(qJ(2));
t174 = t138 * t147;
t86 = t139 * t146 - t142 * t174;
t55 = -t141 * t86 + t145 * t175;
t56 = t141 * t175 + t145 * t86;
t85 = t139 * t142 + t146 * t174;
t17 = Ifges(6,1) * t56 + Ifges(6,4) * t55 + Ifges(6,5) * t85;
t193 = t17 / 0.2e1;
t140 = sin(qJ(6));
t144 = cos(qJ(6));
t30 = -t140 * t56 + t144 * t55;
t192 = t30 / 0.2e1;
t31 = t140 * t55 + t144 * t56;
t191 = t31 / 0.2e1;
t177 = Ifges(6,4) * t145;
t75 = Ifges(6,6) * t142 + (-Ifges(6,2) * t141 + t177) * t146;
t190 = t75 / 0.2e1;
t178 = Ifges(6,4) * t141;
t76 = Ifges(6,5) * t142 + (Ifges(6,1) * t145 - t178) * t146;
t189 = t76 / 0.2e1;
t98 = t140 * t145 + t141 * t144;
t78 = t98 * t146;
t188 = -t78 / 0.2e1;
t97 = -t140 * t141 + t144 * t145;
t80 = t97 * t146;
t187 = t80 / 0.2e1;
t186 = t97 / 0.2e1;
t185 = t98 / 0.2e1;
t184 = -pkin(11) - pkin(10);
t107 = Ifges(6,1) * t141 + t177;
t183 = t107 / 0.2e1;
t182 = pkin(1) * t147;
t116 = pkin(8) * t175;
t87 = t139 * t182 - t116;
t181 = t87 * mrSges(3,1);
t88 = t139 * t143 * pkin(1) + pkin(8) * t174;
t180 = t88 * mrSges(3,2);
t165 = -pkin(2) - t182;
t45 = pkin(3) * t175 + t116 + (-pkin(9) + t165) * t139;
t163 = -qJ(3) * t143 - pkin(1);
t59 = (t147 * t148 + t163) * t138;
t27 = t142 * t45 + t146 * t59;
t20 = pkin(10) * t175 + t27;
t68 = -t139 * qJ(3) - t88;
t58 = pkin(3) * t174 - t68;
t32 = pkin(4) * t85 - pkin(10) * t86 + t58;
t10 = t141 * t32 + t145 * t20;
t33 = -mrSges(6,1) * t55 + mrSges(6,2) * t56;
t61 = mrSges(5,1) * t175 - mrSges(5,3) * t86;
t179 = t61 - t33;
t52 = Ifges(7,5) * t98 + Ifges(7,6) * t97;
t102 = -mrSges(6,1) * t145 + mrSges(6,2) * t141;
t176 = mrSges(5,1) - t102;
t101 = pkin(4) * t142 - pkin(10) * t146 + qJ(3);
t172 = t142 * t148;
t67 = t141 * t101 + t145 * t172;
t173 = t141 * t146;
t171 = t145 * t146;
t170 = t146 * t148;
t96 = mrSges(4,1) * t175 + t139 * mrSges(4,2);
t104 = Ifges(6,5) * t141 + Ifges(6,6) * t145;
t169 = t141 ^ 2 + t145 ^ 2;
t134 = t142 ^ 2;
t136 = t146 ^ 2;
t168 = t134 + t136;
t6 = Ifges(7,5) * t31 + Ifges(7,6) * t30 + Ifges(7,3) * t85;
t15 = Ifges(6,5) * t56 + Ifges(6,6) * t55 + Ifges(6,3) * t85;
t36 = Ifges(7,5) * t80 - Ifges(7,6) * t78 + Ifges(7,3) * t142;
t167 = Ifges(5,5) * t86 - Ifges(5,6) * t85 + Ifges(5,3) * t175;
t166 = t104 / 0.2e1 + t52 / 0.2e1;
t77 = t98 * t142;
t79 = t97 * t142;
t164 = -t77 * mrSges(7,1) - t79 * mrSges(7,2);
t9 = -t141 * t20 + t145 * t32;
t26 = -t142 * t59 + t146 * t45;
t162 = t169 * mrSges(6,3);
t161 = t168 * mrSges(5,3);
t160 = Ifges(3,5) * t175 + Ifges(3,6) * t174 + t139 * t195;
t159 = t10 * t145 - t141 * t9;
t158 = mrSges(6,1) * t141 + mrSges(6,2) * t145;
t93 = t145 * t101;
t66 = -t141 * t172 + t93;
t157 = -t141 * t66 + t145 * t67;
t156 = t142 * t27 + t146 * t26;
t109 = t184 * t141;
t110 = t184 * t145;
t62 = t109 * t144 + t110 * t140;
t63 = t109 * t140 - t110 * t144;
t155 = t62 * mrSges(7,1) - t63 * mrSges(7,2) + t52;
t74 = Ifges(6,5) * t171 - Ifges(6,6) * t173 + Ifges(6,3) * t142;
t4 = pkin(5) * t85 - pkin(11) * t56 + t9;
t5 = pkin(11) * t55 + t10;
t2 = -t140 * t5 + t144 * t4;
t3 = t140 * t4 + t144 * t5;
t154 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t6;
t44 = -pkin(11) * t171 + t93 + (-t141 * t148 + pkin(5)) * t142;
t50 = -pkin(11) * t173 + t67;
t22 = -t140 * t50 + t144 * t44;
t23 = t140 * t44 + t144 * t50;
t153 = t22 * mrSges(7,1) - t23 * mrSges(7,2) + t36;
t152 = (mrSges(7,1) * t144 - mrSges(7,2) * t140) * pkin(5);
t19 = -pkin(4) * t175 - t26;
t149 = qJ(3) ^ 2;
t137 = t148 ^ 2;
t132 = Ifges(5,5) * t146;
t123 = t136 * t148;
t122 = t136 * t137;
t121 = -pkin(5) * t145 - pkin(4);
t108 = Ifges(5,1) * t146 - Ifges(5,4) * t142;
t106 = Ifges(5,4) * t146 - Ifges(5,2) * t142;
t105 = Ifges(6,2) * t145 + t178;
t103 = mrSges(5,1) * t142 + mrSges(5,2) * t146;
t100 = mrSges(6,1) * t142 - mrSges(6,3) * t171;
t99 = -mrSges(6,2) * t142 - mrSges(6,3) * t173;
t95 = -mrSges(4,1) * t174 - mrSges(4,3) * t139;
t94 = (pkin(5) * t141 - t148) * t146;
t89 = t158 * t146;
t70 = t139 * t165 + t116;
t69 = (-pkin(2) * t147 + t163) * t138;
t65 = mrSges(7,1) * t142 - mrSges(7,3) * t80;
t64 = -mrSges(7,2) * t142 - mrSges(7,3) * t78;
t60 = -mrSges(5,2) * t175 - mrSges(5,3) * t85;
t54 = Ifges(7,1) * t98 + Ifges(7,4) * t97;
t53 = Ifges(7,4) * t98 + Ifges(7,2) * t97;
t51 = -mrSges(7,1) * t97 + mrSges(7,2) * t98;
t42 = mrSges(5,1) * t85 + mrSges(5,2) * t86;
t41 = mrSges(7,1) * t78 + mrSges(7,2) * t80;
t40 = Ifges(5,1) * t86 - Ifges(5,4) * t85 + Ifges(5,5) * t175;
t39 = Ifges(5,4) * t86 - Ifges(5,2) * t85 + Ifges(5,6) * t175;
t38 = Ifges(7,1) * t80 - Ifges(7,4) * t78 + Ifges(7,5) * t142;
t37 = Ifges(7,4) * t80 - Ifges(7,2) * t78 + Ifges(7,6) * t142;
t35 = mrSges(6,1) * t85 - mrSges(6,3) * t56;
t34 = -mrSges(6,2) * t85 + mrSges(6,3) * t55;
t16 = Ifges(6,4) * t56 + Ifges(6,2) * t55 + Ifges(6,6) * t85;
t14 = mrSges(7,1) * t85 - mrSges(7,3) * t31;
t13 = -mrSges(7,2) * t85 + mrSges(7,3) * t30;
t12 = -pkin(5) * t55 + t19;
t11 = -mrSges(7,1) * t30 + mrSges(7,2) * t31;
t8 = Ifges(7,1) * t31 + Ifges(7,4) * t30 + Ifges(7,5) * t85;
t7 = Ifges(7,4) * t31 + Ifges(7,2) * t30 + Ifges(7,6) * t85;
t1 = [(0.2e1 * t69 * (mrSges(4,2) * t147 - mrSges(4,3) * t143) + t143 * t167 + (m(3) * pkin(1) ^ 2 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2)) * t147) * t147 + (-0.2e1 * pkin(1) * mrSges(3,2) + (Ifges(4,2) + Ifges(3,1)) * t143 + 0.2e1 * (Ifges(3,4) + Ifges(4,6)) * t147) * t143) * t138 + 0.2e1 * (-t143 * t87 + t147 * t88) * mrSges(3,3) + ((-(2 * Ifges(4,5)) + Ifges(3,6)) * t147 + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t143) * t139) * t138 + (t160 - 0.2e1 * t180 + 0.2e1 * t181) * t139 + m(3) * (t87 ^ 2 + t88 ^ 2) + t86 * t40 + 0.2e1 * t68 * t95 + 0.2e1 * t70 * t96 + t56 * t17 + 0.2e1 * t58 * t42 + 0.2e1 * t27 * t60 + 0.2e1 * t26 * t61 + t55 * t16 + 0.2e1 * t9 * t35 + t30 * t7 + t31 * t8 + 0.2e1 * t19 * t33 + 0.2e1 * t10 * t34 + 0.2e1 * t12 * t11 + 0.2e1 * t3 * t13 + 0.2e1 * t2 * t14 + m(4) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(5) * (t26 ^ 2 + t27 ^ 2 + t58 ^ 2) + m(6) * (t10 ^ 2 + t19 ^ 2 + t9 ^ 2) + m(7) * (t12 ^ 2 + t2 ^ 2 + t3 ^ 2) + (t15 + t6 - t39) * t85 + Ifges(2,3); (t40 / 0.2e1 + t145 * t193 - t141 * t16 / 0.2e1 - t26 * mrSges(5,3) + t179 * t148) * t146 - t180 + t181 + t160 + m(6) * (t10 * t67 - t170 * t19 + t66 * t9) + m(5) * (qJ(3) * t58 + t148 * t156) + t10 * t99 + t9 * t100 + t58 * t103 + t86 * t108 / 0.2e1 + t19 * t89 + t94 * t11 - pkin(2) * t96 + t67 * t34 - t68 * mrSges(4,3) + t70 * mrSges(4,2) + t3 * t64 + t2 * t65 + t66 * t35 + t12 * t41 + t22 * t14 + t23 * t13 + (-t39 / 0.2e1 + t15 / 0.2e1 + t6 / 0.2e1 + t148 * t60 - t27 * mrSges(5,3)) * t142 + (-t95 + t42) * qJ(3) + m(7) * (t12 * t94 + t2 * t22 + t23 * t3) + m(4) * (-pkin(2) * t70 - qJ(3) * t68) + (-t106 / 0.2e1 + t74 / 0.2e1 + t36 / 0.2e1) * t85 + t8 * t187 + t7 * t188 + t56 * t189 + t55 * t190 + t38 * t191 + t37 * t192 + (-Ifges(4,5) * t147 + t143 * (-Ifges(5,6) * t142 + t132) / 0.2e1 - Ifges(4,4) * t143) * t138; -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t66 * t100 + 0.2e1 * t22 * t65 + 0.2e1 * t23 * t64 - t78 * t37 + t80 * t38 + 0.2e1 * t94 * t41 + 0.2e1 * t67 * t99 + (-t106 + t74 + t36) * t142 + (-t141 * t75 + t145 * t76 + t89 * t194 + t108) * t146 + m(5) * (t134 * t137 + t122 + t149) + m(4) * (pkin(2) ^ 2 + t149) + m(6) * (t66 ^ 2 + t67 ^ 2 + t122) + m(7) * (t22 ^ 2 + t23 ^ 2 + t94 ^ 2) + 0.2e1 * (t103 + mrSges(4,3)) * qJ(3) + t161 * t194 + t195; t79 * t13 - t77 * t14 + (-t11 + t179) * t146 + (-t141 * t35 + t145 * t34 + t60) * t142 + m(7) * (-t12 * t146 - t2 * t77 + t3 * t79) + m(6) * (t142 * t159 - t146 * t19) + m(5) * t156 + m(4) * t70 + t96; -m(4) * pkin(2) + t79 * t64 - t77 * t65 + mrSges(4,2) + (-t41 - t89) * t146 + (-t141 * t100 + t145 * t99) * t142 - t161 + m(7) * (-t146 * t94 - t22 * t77 + t23 * t79) + m(6) * (t142 * t157 + t123) + m(5) * (t134 * t148 + t123); m(4) + m(5) * t168 + m(6) * (t134 * t169 + t136) + m(7) * (t77 ^ 2 + t79 ^ 2 + t136); t167 + t166 * t85 + m(6) * (-pkin(4) * t19 + pkin(10) * t159) + t121 * t11 + t8 * t185 + t19 * t102 + t55 * t105 / 0.2e1 + t56 * t183 + t7 * t186 + t62 * t14 + t63 * t13 + t12 * t51 + t53 * t192 + t54 * t191 - t27 * mrSges(5,2) - pkin(4) * t33 + t26 * mrSges(5,1) + (t16 / 0.2e1 + pkin(10) * t34 + t10 * mrSges(6,3)) * t145 + (-t9 * mrSges(6,3) - pkin(10) * t35 + t193) * t141 + m(7) * (t12 * t121 + t2 * t62 + t3 * t63) + (-t2 * t98 + t3 * t97) * mrSges(7,3); t132 + t121 * t41 - pkin(4) * t89 + t94 * t51 + t37 * t186 + t38 * t185 + t54 * t187 + t53 * t188 + t63 * t64 + t62 * t65 + t176 * t170 + (-t22 * t98 + t23 * t97) * mrSges(7,3) + (t67 * mrSges(6,3) + pkin(10) * t99 + t146 * t183 + t190) * t145 + (t189 - t146 * t105 / 0.2e1 - pkin(10) * t100 - t66 * mrSges(6,3)) * t141 + m(6) * (pkin(4) * t170 + pkin(10) * t157) + m(7) * (t121 * t94 + t22 * t62 + t23 * t63) + (-t148 * mrSges(5,2) - Ifges(5,6) + t166) * t142; (t77 * t98 + t79 * t97) * mrSges(7,3) + (-t51 + t176) * t146 + (-mrSges(5,2) + t162) * t142 + m(6) * (pkin(10) * t142 * t169 + pkin(4) * t146) + m(7) * (-t121 * t146 - t62 * t77 + t63 * t79); -0.2e1 * pkin(4) * t102 + t145 * t105 + t141 * t107 + 0.2e1 * t121 * t51 + t97 * t53 + t98 * t54 + Ifges(5,3) + m(7) * (t121 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(6) * (pkin(10) ^ 2 * t169 + pkin(4) ^ 2) + 0.2e1 * (-t62 * t98 + t63 * t97) * mrSges(7,3) + 0.2e1 * pkin(10) * t162; t9 * mrSges(6,1) - t10 * mrSges(6,2) + (t144 * t14 + m(7) * (t140 * t3 + t144 * t2) + t140 * t13) * pkin(5) + t154 + t15; t66 * mrSges(6,1) - t67 * mrSges(6,2) + (m(7) * (t140 * t23 + t144 * t22) + t140 * t64 + t144 * t65) * pkin(5) + t153 + t74; -t158 * t142 + m(7) * (t140 * t79 - t144 * t77) * pkin(5) + t164; -t158 * pkin(10) + (m(7) * (t140 * t63 + t144 * t62) + (t140 * t97 - t144 * t98) * mrSges(7,3)) * pkin(5) + t155 + t104; Ifges(6,3) + Ifges(7,3) + m(7) * (t140 ^ 2 + t144 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t152; t154; t153; t164; t155; Ifges(7,3) + t152; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
