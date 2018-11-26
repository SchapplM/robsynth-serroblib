% Calculate joint inertia matrix for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2018-11-23 17:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR13_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR13_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR13_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:09:13
% EndTime: 2018-11-23 17:09:14
% DurationCPUTime: 1.66s
% Computational Cost: add. (2503->387), mult. (5481->545), div. (0->0), fcn. (5717->10), ass. (0->147)
t186 = Ifges(4,1) + Ifges(3,3);
t144 = (-pkin(2) - pkin(9));
t185 = -2 * t144;
t134 = sin(pkin(11));
t136 = cos(pkin(11));
t135 = sin(pkin(6));
t140 = sin(qJ(2));
t163 = t135 * t140;
t137 = cos(pkin(6));
t139 = sin(qJ(4));
t142 = cos(qJ(4));
t143 = cos(qJ(2));
t162 = t135 * t143;
t82 = t137 * t142 - t139 * t162;
t53 = -t134 * t82 + t136 * t163;
t54 = t134 * t163 + t136 * t82;
t81 = t137 * t139 + t142 * t162;
t16 = Ifges(6,1) * t54 + Ifges(6,4) * t53 + Ifges(6,5) * t81;
t184 = t16 / 0.2e1;
t138 = sin(qJ(6));
t141 = cos(qJ(6));
t25 = -t138 * t54 + t141 * t53;
t183 = t25 / 0.2e1;
t26 = t138 * t53 + t141 * t54;
t182 = t26 / 0.2e1;
t166 = Ifges(6,4) * t136;
t71 = Ifges(6,6) * t139 + (-Ifges(6,2) * t134 + t166) * t142;
t181 = t71 / 0.2e1;
t167 = Ifges(6,4) * t134;
t72 = Ifges(6,5) * t139 + (Ifges(6,1) * t136 - t167) * t142;
t180 = t72 / 0.2e1;
t93 = t134 * t141 + t136 * t138;
t75 = t93 * t142;
t179 = -t75 / 0.2e1;
t92 = -t134 * t138 + t136 * t141;
t77 = t92 * t142;
t178 = t77 / 0.2e1;
t177 = t92 / 0.2e1;
t176 = t93 / 0.2e1;
t175 = m(6) + m(7);
t104 = Ifges(6,1) * t134 + t166;
t174 = t104 / 0.2e1;
t173 = pkin(1) * t143;
t115 = pkin(8) * t163;
t84 = t137 * t173 - t115;
t172 = t84 * mrSges(3,1);
t85 = t137 * t140 * pkin(1) + pkin(8) * t162;
t171 = t85 * mrSges(3,2);
t170 = pkin(10) + qJ(5);
t154 = -pkin(2) - t173;
t44 = pkin(3) * t163 + t115 + (-pkin(9) + t154) * t137;
t153 = -qJ(3) * t140 - pkin(1);
t58 = (t144 * t143 + t153) * t135;
t30 = t139 * t44 + t142 * t58;
t20 = qJ(5) * t163 + t30;
t68 = -t137 * qJ(3) - t85;
t57 = pkin(3) * t162 - t68;
t32 = pkin(4) * t81 - qJ(5) * t82 + t57;
t9 = t134 * t32 + t136 * t20;
t40 = t75 * mrSges(7,1) + t77 * mrSges(7,2);
t161 = t136 * t142;
t164 = t134 * t142;
t83 = mrSges(6,1) * t164 + mrSges(6,2) * t161;
t169 = -t40 - t83;
t31 = -t53 * mrSges(6,1) + t54 * mrSges(6,2);
t60 = mrSges(5,1) * t163 - mrSges(5,3) * t82;
t168 = t60 - t31;
t48 = Ifges(7,5) * t93 + Ifges(7,6) * t92;
t100 = -t136 * mrSges(6,1) + t134 * mrSges(6,2);
t165 = mrSges(5,1) - t100;
t160 = t139 * t144;
t98 = pkin(4) * t139 - qJ(5) * t142 + qJ(3);
t64 = t134 * t98 + t136 * t160;
t159 = t142 * t144;
t95 = mrSges(4,1) * t163 + t137 * mrSges(4,2);
t158 = t134 ^ 2 + t136 ^ 2;
t131 = t139 ^ 2;
t132 = t142 ^ 2;
t157 = t131 + t132;
t5 = Ifges(7,5) * t26 + Ifges(7,6) * t25 + Ifges(7,3) * t81;
t35 = Ifges(7,5) * t77 - Ifges(7,6) * t75 + Ifges(7,3) * t139;
t156 = Ifges(5,5) * t82 - Ifges(5,6) * t81 + Ifges(5,3) * t163;
t155 = Ifges(6,5) * t134 / 0.2e1 + Ifges(6,6) * t136 / 0.2e1 + t48 / 0.2e1;
t10 = -t25 * mrSges(7,1) + t26 * mrSges(7,2);
t47 = -t92 * mrSges(7,1) + t93 * mrSges(7,2);
t8 = -t134 * t20 + t136 * t32;
t29 = -t139 * t58 + t142 * t44;
t152 = t157 * mrSges(5,3);
t151 = qJ(5) * t158;
t150 = Ifges(3,5) * t163 + Ifges(3,6) * t162 + t137 * t186;
t149 = -t8 * t134 + t9 * t136;
t90 = t136 * t98;
t63 = -t134 * t160 + t90;
t148 = -t63 * t134 + t64 * t136;
t147 = t30 * t139 + t29 * t142;
t21 = -pkin(4) * t163 - t29;
t146 = qJ(3) ^ 2;
t133 = t144 ^ 2;
t128 = Ifges(5,5) * t142;
t121 = t132 * t144;
t120 = t132 * t133;
t119 = -pkin(5) * t136 - pkin(4);
t107 = Ifges(5,1) * t142 - Ifges(5,4) * t139;
t106 = Ifges(5,4) * t142 - Ifges(5,2) * t139;
t105 = mrSges(5,1) * t139 + mrSges(5,2) * t142;
t103 = Ifges(6,2) * t136 + t167;
t101 = t170 * t136;
t99 = t170 * t134;
t97 = mrSges(6,1) * t139 - mrSges(6,3) * t161;
t96 = -mrSges(6,2) * t139 - mrSges(6,3) * t164;
t94 = -mrSges(4,1) * t162 - mrSges(4,3) * t137;
t91 = (pkin(5) * t134 - t144) * t142;
t76 = t92 * t139;
t74 = t93 * t139;
t73 = t154 * t137 + t115;
t70 = Ifges(6,3) * t139 + (Ifges(6,5) * t136 - Ifges(6,6) * t134) * t142;
t69 = (-pkin(2) * t143 + t153) * t135;
t62 = mrSges(7,1) * t139 - mrSges(7,3) * t77;
t61 = -mrSges(7,2) * t139 - mrSges(7,3) * t75;
t59 = -mrSges(5,2) * t163 - mrSges(5,3) * t81;
t56 = t101 * t141 - t138 * t99;
t55 = -t101 * t138 - t141 * t99;
t50 = Ifges(7,1) * t93 + Ifges(7,4) * t92;
t49 = Ifges(7,4) * t93 + Ifges(7,2) * t92;
t46 = -pkin(10) * t164 + t64;
t42 = -pkin(10) * t161 + t90 + (-t134 * t144 + pkin(5)) * t139;
t41 = mrSges(5,1) * t81 + mrSges(5,2) * t82;
t39 = Ifges(5,1) * t82 - Ifges(5,4) * t81 + Ifges(5,5) * t163;
t38 = Ifges(5,4) * t82 - Ifges(5,2) * t81 + Ifges(5,6) * t163;
t37 = Ifges(7,1) * t77 - Ifges(7,4) * t75 + Ifges(7,5) * t139;
t36 = Ifges(7,4) * t77 - Ifges(7,2) * t75 + Ifges(7,6) * t139;
t34 = mrSges(6,1) * t81 - mrSges(6,3) * t54;
t33 = -mrSges(6,2) * t81 + mrSges(6,3) * t53;
t19 = t138 * t42 + t141 * t46;
t18 = -t138 * t46 + t141 * t42;
t15 = Ifges(6,4) * t54 + Ifges(6,2) * t53 + Ifges(6,6) * t81;
t14 = Ifges(6,5) * t54 + Ifges(6,6) * t53 + Ifges(6,3) * t81;
t13 = mrSges(7,1) * t81 - mrSges(7,3) * t26;
t12 = -mrSges(7,2) * t81 + mrSges(7,3) * t25;
t11 = -pkin(5) * t53 + t21;
t7 = Ifges(7,1) * t26 + Ifges(7,4) * t25 + Ifges(7,5) * t81;
t6 = Ifges(7,4) * t26 + Ifges(7,2) * t25 + Ifges(7,6) * t81;
t4 = pkin(10) * t53 + t9;
t3 = pkin(5) * t81 - pkin(10) * t54 + t8;
t2 = t138 * t3 + t141 * t4;
t1 = -t138 * t4 + t141 * t3;
t17 = [0.2e1 * t68 * t94 + 0.2e1 * t73 * t95 + t82 * t39 + t53 * t15 + t54 * t16 + 0.2e1 * t57 * t41 + 0.2e1 * t30 * t59 + 0.2e1 * t29 * t60 + 0.2e1 * t21 * t31 + 0.2e1 * t9 * t33 + 0.2e1 * t8 * t34 + t25 * t6 + t26 * t7 + 0.2e1 * t2 * t12 + 0.2e1 * t1 * t13 + 0.2e1 * t11 * t10 + m(3) * (t84 ^ 2 + t85 ^ 2) + (0.2e1 * t69 * (mrSges(4,2) * t143 - mrSges(4,3) * t140) + t140 * t156 + (m(3) * pkin(1) ^ 2 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2)) * t143) * t143 + (-0.2e1 * pkin(1) * mrSges(3,2) + (Ifges(4,2) + Ifges(3,1)) * t140 + 0.2e1 * (Ifges(3,4) + Ifges(4,6)) * t143) * t140) * t135 + 0.2e1 * (-t84 * t140 + t85 * t143) * mrSges(3,3) + ((-(2 * Ifges(4,5)) + Ifges(3,6)) * t143 + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t140) * t137) * t135 + m(4) * (t68 ^ 2 + t69 ^ 2 + t73 ^ 2) + m(5) * (t29 ^ 2 + t30 ^ 2 + t57 ^ 2) + m(6) * (t21 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(7) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2) + (-t38 + t14 + t5) * t81 + (t150 - 0.2e1 * t171 + 0.2e1 * t172) * t137 + Ifges(2,3); (-t94 + t41) * qJ(3) + m(7) * (t1 * t18 + t11 * t91 + t19 * t2) + m(4) * (-pkin(2) * t73 - qJ(3) * t68) + (-t106 / 0.2e1 + t70 / 0.2e1 + t35 / 0.2e1) * t81 + (-Ifges(4,5) * t143 + t140 * (-Ifges(5,6) * t139 + t128) / 0.2e1 - Ifges(4,4) * t140) * t135 - pkin(2) * t95 + t9 * t96 + t8 * t97 + t57 * t105 + t82 * t107 / 0.2e1 + t21 * t83 + t91 * t10 + t1 * t62 + t63 * t34 + t64 * t33 - t68 * mrSges(4,3) + t73 * mrSges(4,2) + t2 * t61 + t11 * t40 + t18 * t13 + t19 * t12 + t150 + t7 * t178 + t6 * t179 + t54 * t180 + t53 * t181 + t37 * t182 + t36 * t183 + (t144 * t59 - t38 / 0.2e1 + t14 / 0.2e1 + t5 / 0.2e1 - t30 * mrSges(5,3)) * t139 + m(6) * (-t21 * t159 + t63 * t8 + t64 * t9) + m(5) * (qJ(3) * t57 + t147 * t144) + (t136 * t184 - t134 * t15 / 0.2e1 - t29 * mrSges(5,3) + t39 / 0.2e1 + t168 * t144) * t142 - t171 + t172; -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t18 * t62 + 0.2e1 * t19 * t61 - t75 * t36 + t77 * t37 + 0.2e1 * t91 * t40 + 0.2e1 * t63 * t97 + 0.2e1 * t64 * t96 + (-t106 + t70 + t35) * t139 + (-t134 * t71 + t136 * t72 + t83 * t185 + t107) * t142 + m(5) * (t131 * t133 + t120 + t146) + m(4) * (pkin(2) ^ 2 + t146) + m(6) * (t63 ^ 2 + t64 ^ 2 + t120) + m(7) * (t18 ^ 2 + t19 ^ 2 + t91 ^ 2) + 0.2e1 * (t105 + mrSges(4,3)) * qJ(3) + t152 * t185 + t186; t76 * t12 - t74 * t13 + (-t10 + t168) * t142 + (-t134 * t34 + t136 * t33 + t59) * t139 + m(7) * (-t1 * t74 - t11 * t142 + t2 * t76) + m(6) * (t149 * t139 - t142 * t21) + m(5) * t147 + m(4) * t73 + t95; -m(4) * pkin(2) + t76 * t61 - t74 * t62 + mrSges(4,2) + t169 * t142 + (-t134 * t97 + t136 * t96) * t139 - t152 + m(7) * (-t142 * t91 - t18 * t74 + t19 * t76) + m(6) * (t148 * t139 + t121) + m(5) * (t131 * t144 + t121); m(4) + m(5) * t157 + m(6) * (t158 * t131 + t132) + m(7) * (t74 ^ 2 + t76 ^ 2 + t132); (t15 / 0.2e1 + qJ(5) * t33 + t9 * mrSges(6,3)) * t136 + (-t8 * mrSges(6,3) - qJ(5) * t34 + t184) * t134 + m(7) * (t1 * t55 + t11 * t119 + t2 * t56) + t156 + (-t1 * t93 + t2 * t92) * mrSges(7,3) + t119 * t10 + t7 * t176 + t21 * t100 + t53 * t103 / 0.2e1 + t54 * t174 + t6 * t177 + t55 * t13 + t56 * t12 + t11 * t47 + t49 * t183 + t50 * t182 + t29 * mrSges(5,1) - t30 * mrSges(5,2) - pkin(4) * t31 + t155 * t81 + m(6) * (-pkin(4) * t21 + t149 * qJ(5)); t119 * t40 + t128 + t37 * t176 - pkin(4) * t83 + t91 * t47 + t36 * t177 + t49 * t179 + t50 * t178 + t55 * t62 + t56 * t61 + t165 * t159 + (-t18 * t93 + t19 * t92) * mrSges(7,3) + (t64 * mrSges(6,3) + qJ(5) * t96 + t142 * t174 + t181) * t136 + (-t142 * t103 / 0.2e1 + t180 - qJ(5) * t97 - t63 * mrSges(6,3)) * t134 + m(6) * (pkin(4) * t159 + t148 * qJ(5)) + m(7) * (t119 * t91 + t18 * t55 + t19 * t56) + (-t144 * mrSges(5,2) - Ifges(5,6) + t155) * t139; (t74 * t93 + t76 * t92) * mrSges(7,3) + (-t47 + t165) * t142 + (t158 * mrSges(6,3) - mrSges(5,2)) * t139 + m(6) * (pkin(4) * t142 + t139 * t151) + m(7) * (-t119 * t142 - t55 * t74 + t56 * t76); -0.2e1 * pkin(4) * t100 + t136 * t103 + t134 * t104 + 0.2e1 * t119 * t47 + t92 * t49 + t93 * t50 + Ifges(5,3) + m(7) * (t119 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(6) * (t158 * qJ(5) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t55 * t93 + t56 * t92) * mrSges(7,3) + 0.2e1 * mrSges(6,3) * t151; m(6) * t21 + m(7) * t11 + t10 + t31; -m(6) * t159 + m(7) * t91 - t169; -t175 * t142; -m(6) * pkin(4) + m(7) * t119 + t100 + t47; t175; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t5; mrSges(7,1) * t18 - mrSges(7,2) * t19 + t35; -mrSges(7,1) * t74 - mrSges(7,2) * t76; mrSges(7,1) * t55 - t56 * mrSges(7,2) + t48; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t17(1) t17(2) t17(4) t17(7) t17(11) t17(16); t17(2) t17(3) t17(5) t17(8) t17(12) t17(17); t17(4) t17(5) t17(6) t17(9) t17(13) t17(18); t17(7) t17(8) t17(9) t17(10) t17(14) t17(19); t17(11) t17(12) t17(13) t17(14) t17(15) t17(20); t17(16) t17(17) t17(18) t17(19) t17(20) t17(21);];
Mq  = res;
