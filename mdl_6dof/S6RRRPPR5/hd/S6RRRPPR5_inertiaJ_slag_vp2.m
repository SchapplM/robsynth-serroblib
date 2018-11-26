% Calculate joint inertia matrix for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2018-11-23 17:36
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:35:30
% EndTime: 2018-11-23 17:35:32
% DurationCPUTime: 1.75s
% Computational Cost: add. (4401->405), mult. (9788->596), div. (0->0), fcn. (11196->12), ass. (0->145)
t158 = sin(pkin(6));
t167 = cos(qJ(2));
t178 = t158 * t167;
t157 = sin(pkin(11));
t160 = cos(pkin(11));
t163 = sin(qJ(3));
t166 = cos(qJ(3));
t123 = t157 * t163 - t160 * t166;
t125 = t157 * t166 + t160 * t163;
t201 = Ifges(4,5) * t163 + Ifges(5,5) * t125 + Ifges(4,6) * t166 - Ifges(5,6) * t123;
t186 = -qJ(4) - pkin(9);
t134 = t186 * t166;
t171 = t186 * t163;
t95 = -t134 * t157 - t160 * t171;
t200 = t95 ^ 2;
t199 = 0.2e1 * t95;
t156 = sin(pkin(12));
t159 = cos(pkin(12));
t161 = cos(pkin(6));
t164 = sin(qJ(2));
t179 = t158 * t164;
t111 = t161 * t166 - t163 * t179;
t112 = t161 * t163 + t166 * t179;
t76 = t111 * t157 + t112 * t160;
t54 = -t156 * t76 - t159 * t178;
t55 = -t156 * t178 + t159 * t76;
t75 = -t160 * t111 + t112 * t157;
t19 = Ifges(6,1) * t55 + Ifges(6,4) * t54 + Ifges(6,5) * t75;
t198 = t19 / 0.2e1;
t184 = Ifges(6,4) * t159;
t57 = Ifges(6,6) * t123 + (-Ifges(6,2) * t156 + t184) * t125;
t197 = t57 / 0.2e1;
t185 = Ifges(6,4) * t156;
t58 = Ifges(6,5) * t123 + (Ifges(6,1) * t159 - t185) * t125;
t196 = t58 / 0.2e1;
t162 = sin(qJ(6));
t165 = cos(qJ(6));
t124 = -t156 * t162 + t159 * t165;
t126 = t156 * t165 + t159 * t162;
t90 = Ifges(7,4) * t126 + Ifges(7,2) * t124;
t195 = t90 / 0.2e1;
t92 = Ifges(7,1) * t126 + Ifges(7,4) * t124;
t194 = t92 / 0.2e1;
t193 = t124 / 0.2e1;
t192 = t126 / 0.2e1;
t132 = Ifges(6,1) * t156 + t184;
t191 = t132 / 0.2e1;
t189 = pkin(1) * t167;
t188 = Ifges(4,3) + Ifges(5,3);
t144 = pkin(3) * t157 + qJ(5);
t187 = pkin(10) + t144;
t116 = t161 * t164 * pkin(1) + pkin(8) * t178;
t104 = pkin(9) * t161 + t116;
t105 = (-pkin(2) * t167 - pkin(9) * t164 - pkin(1)) * t158;
t61 = -t104 * t163 + t166 * t105;
t45 = -pkin(3) * t178 - qJ(4) * t112 + t61;
t62 = t166 * t104 + t163 * t105;
t50 = qJ(4) * t111 + t62;
t23 = t157 * t45 + t160 * t50;
t20 = -qJ(5) * t178 + t23;
t139 = pkin(8) * t179;
t103 = t139 + (-pkin(2) - t189) * t161;
t81 = -pkin(3) * t111 + t103;
t26 = pkin(4) * t75 - qJ(5) * t76 + t81;
t6 = t156 * t26 + t159 * t20;
t147 = -pkin(3) * t166 - pkin(2);
t84 = pkin(4) * t123 - qJ(5) * t125 + t147;
t97 = -t160 * t134 + t157 * t171;
t47 = t156 * t84 + t159 * t97;
t115 = t161 * t189 - t139;
t183 = t115 * mrSges(3,1);
t182 = t116 * mrSges(3,2);
t181 = t125 * t156;
t180 = t125 * t159;
t82 = mrSges(6,1) * t181 + mrSges(6,2) * t180;
t89 = Ifges(7,5) * t126 + Ifges(7,6) * t124;
t175 = t156 ^ 2 + t159 ^ 2;
t174 = t163 ^ 2 + t166 ^ 2;
t30 = -t162 * t55 + t165 * t54;
t31 = t162 * t54 + t165 * t55;
t7 = Ifges(7,5) * t31 + Ifges(7,6) * t30 + Ifges(7,3) * t75;
t73 = t126 * t125;
t74 = t124 * t125;
t34 = Ifges(7,5) * t74 - Ifges(7,6) * t73 + Ifges(7,3) * t123;
t173 = t89 / 0.2e1 + Ifges(6,5) * t156 / 0.2e1 + Ifges(6,6) * t159 / 0.2e1;
t172 = Ifges(3,5) * t179 + Ifges(3,6) * t178 + Ifges(3,3) * t161;
t146 = -pkin(3) * t160 - pkin(4);
t43 = t75 * mrSges(5,1) + t76 * mrSges(5,2);
t32 = -t54 * mrSges(6,1) + t55 * mrSges(6,2);
t10 = -t30 * mrSges(7,1) + t31 * mrSges(7,2);
t42 = t73 * mrSges(7,1) + t74 * mrSges(7,2);
t5 = -t156 * t20 + t159 * t26;
t46 = -t156 * t97 + t159 * t84;
t22 = -t157 * t50 + t160 * t45;
t87 = t123 * mrSges(5,1) + t125 * mrSges(5,2);
t129 = -t159 * mrSges(6,1) + t156 * mrSges(6,2);
t170 = -Ifges(4,5) * t112 - Ifges(5,5) * t76 - Ifges(4,6) * t111 + Ifges(5,6) * t75;
t21 = pkin(4) * t178 - t22;
t88 = -t124 * mrSges(7,1) + t126 * mrSges(7,2);
t136 = Ifges(4,1) * t163 + Ifges(4,4) * t166;
t135 = Ifges(4,4) * t163 + Ifges(4,2) * t166;
t133 = -mrSges(4,1) * t166 + mrSges(4,2) * t163;
t131 = Ifges(6,2) * t159 + t185;
t128 = -pkin(5) * t159 + t146;
t114 = t187 * t159;
t113 = t187 * t156;
t99 = -mrSges(4,1) * t178 - mrSges(4,3) * t112;
t98 = mrSges(4,2) * t178 + mrSges(4,3) * t111;
t93 = Ifges(5,1) * t125 - Ifges(5,4) * t123;
t91 = Ifges(5,4) * t125 - Ifges(5,2) * t123;
t86 = mrSges(6,1) * t123 - mrSges(6,3) * t180;
t85 = -mrSges(6,2) * t123 - mrSges(6,3) * t181;
t83 = -mrSges(4,1) * t111 + mrSges(4,2) * t112;
t80 = -t113 * t162 + t114 * t165;
t79 = -t113 * t165 - t114 * t162;
t65 = Ifges(4,1) * t112 + Ifges(4,4) * t111 - Ifges(4,5) * t178;
t64 = Ifges(4,4) * t112 + Ifges(4,2) * t111 - Ifges(4,6) * t178;
t63 = pkin(5) * t181 + t95;
t60 = -mrSges(5,1) * t178 - mrSges(5,3) * t76;
t59 = mrSges(5,2) * t178 - mrSges(5,3) * t75;
t56 = Ifges(6,3) * t123 + (Ifges(6,5) * t159 - Ifges(6,6) * t156) * t125;
t52 = mrSges(7,1) * t123 - mrSges(7,3) * t74;
t51 = -mrSges(7,2) * t123 - mrSges(7,3) * t73;
t41 = -pkin(10) * t181 + t47;
t40 = Ifges(5,1) * t76 - Ifges(5,4) * t75 - Ifges(5,5) * t178;
t39 = Ifges(5,4) * t76 - Ifges(5,2) * t75 - Ifges(5,6) * t178;
t38 = mrSges(6,1) * t75 - mrSges(6,3) * t55;
t37 = -mrSges(6,2) * t75 + mrSges(6,3) * t54;
t36 = Ifges(7,1) * t74 - Ifges(7,4) * t73 + Ifges(7,5) * t123;
t35 = Ifges(7,4) * t74 - Ifges(7,2) * t73 + Ifges(7,6) * t123;
t33 = pkin(5) * t123 - pkin(10) * t180 + t46;
t18 = Ifges(6,4) * t55 + Ifges(6,2) * t54 + Ifges(6,6) * t75;
t17 = Ifges(6,5) * t55 + Ifges(6,6) * t54 + Ifges(6,3) * t75;
t15 = mrSges(7,1) * t75 - mrSges(7,3) * t31;
t14 = -mrSges(7,2) * t75 + mrSges(7,3) * t30;
t13 = t162 * t33 + t165 * t41;
t12 = -t162 * t41 + t165 * t33;
t11 = -pkin(5) * t54 + t21;
t9 = Ifges(7,1) * t31 + Ifges(7,4) * t30 + Ifges(7,5) * t75;
t8 = Ifges(7,4) * t31 + Ifges(7,2) * t30 + Ifges(7,6) * t75;
t4 = pkin(10) * t54 + t6;
t3 = pkin(5) * t75 - pkin(10) * t55 + t5;
t2 = t162 * t3 + t165 * t4;
t1 = -t162 * t4 + t165 * t3;
t16 = [m(6) * (t21 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2 + t81 ^ 2) + m(4) * (t103 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(3) * (pkin(1) ^ 2 * t158 ^ 2 + t115 ^ 2 + t116 ^ 2) + m(7) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2) + ((-0.2e1 * t115 * mrSges(3,3) + Ifges(3,5) * t161 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t164) * t158) * t164 + (0.2e1 * t116 * mrSges(3,3) + Ifges(3,6) * t161 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t164 + (Ifges(3,2) + t188) * t167) * t158 + t170) * t167) * t158 + (t7 + t17 - t39) * t75 + (t172 - 0.2e1 * t182 + 0.2e1 * t183) * t161 + 0.2e1 * t11 * t10 + Ifges(2,3) + 0.2e1 * t2 * t14 + 0.2e1 * t1 * t15 + t30 * t8 + t31 * t9 + 0.2e1 * t21 * t32 + 0.2e1 * t6 * t37 + 0.2e1 * t5 * t38 + t54 * t18 + t55 * t19 + 0.2e1 * t23 * t59 + 0.2e1 * t22 * t60 + t76 * t40 + 0.2e1 * t81 * t43 + 0.2e1 * t62 * t98 + 0.2e1 * t61 * t99 + 0.2e1 * t103 * t83 + t111 * t64 + t112 * t65; (t34 / 0.2e1 + t56 / 0.2e1 - t91 / 0.2e1) * t75 + (t62 * mrSges(4,3) + pkin(9) * t98 + t64 / 0.2e1) * t166 + (-t61 * mrSges(4,3) - pkin(9) * t99 + t65 / 0.2e1) * t163 + (-t22 * mrSges(5,3) - t156 * t18 / 0.2e1 + t159 * t198 + t40 / 0.2e1) * t125 + m(7) * (t1 * t12 + t11 * t63 + t13 * t2) + m(6) * (t21 * t95 + t46 * t5 + t47 * t6) + m(5) * (t147 * t81 - t22 * t95 + t23 * t97) + (-t23 * mrSges(5,3) + t7 / 0.2e1 + t17 / 0.2e1 - t39 / 0.2e1) * t123 + t172 - t201 * t178 / 0.2e1 + t55 * t196 + t54 * t197 - t182 + t183 + m(4) * (-pkin(2) * t103 + (-t61 * t163 + t62 * t166) * pkin(9)) + t13 * t14 + t12 * t15 + t30 * t35 / 0.2e1 + t31 * t36 / 0.2e1 + t11 * t42 + t46 * t38 + t47 * t37 + t2 * t51 + t1 * t52 + t63 * t10 - t73 * t8 / 0.2e1 + t74 * t9 / 0.2e1 + t21 * t82 - pkin(2) * t83 + t6 * t85 + t5 * t86 + t81 * t87 + t76 * t93 / 0.2e1 + t97 * t59 + t103 * t133 + t111 * t135 / 0.2e1 + t112 * t136 / 0.2e1 + t147 * t43 + (t32 - t60) * t95; -0.2e1 * pkin(2) * t133 + 0.2e1 * t12 * t52 + 0.2e1 * t13 * t51 + t166 * t135 + t163 * t136 + 0.2e1 * t147 * t87 - t73 * t35 + t74 * t36 + 0.2e1 * t63 * t42 + 0.2e1 * t46 * t86 + 0.2e1 * t47 * t85 + t82 * t199 + Ifges(3,3) + 0.2e1 * t174 * pkin(9) * mrSges(4,3) + (mrSges(5,3) * t199 - t156 * t57 + t159 * t58 + t93) * t125 + (-0.2e1 * mrSges(5,3) * t97 + t34 + t56 - t91) * t123 + m(7) * (t12 ^ 2 + t13 ^ 2 + t63 ^ 2) + m(6) * (t46 ^ 2 + t47 ^ 2 + t200) + m(5) * (t147 ^ 2 + t97 ^ 2 + t200) + m(4) * (pkin(9) ^ 2 * t174 + pkin(2) ^ 2); (-t5 * mrSges(6,3) - t144 * t38 + t198) * t156 - t188 * t178 + (-t1 * t126 + t2 * t124) * mrSges(7,3) + t8 * t193 + t31 * t194 + t30 * t195 + t55 * t191 + t9 * t192 - t170 + (t157 * t59 + t160 * t60 + m(5) * (t157 * t23 + t160 * t22)) * pkin(3) + m(6) * (t146 * t21 + (-t5 * t156 + t6 * t159) * t144) + m(7) * (t1 * t79 + t11 * t128 + t2 * t80) + t22 * mrSges(5,1) - t23 * mrSges(5,2) + t61 * mrSges(4,1) - t62 * mrSges(4,2) + t79 * t15 + t80 * t14 + t11 * t88 + t128 * t10 + t21 * t129 + t54 * t131 / 0.2e1 + t146 * t32 + t173 * t75 + (t6 * mrSges(6,3) + t144 * t37 + t18 / 0.2e1) * t159; m(7) * (t12 * t79 + t128 * t63 + t13 * t80) + (-mrSges(4,1) * t163 - mrSges(4,2) * t166) * pkin(9) + (-t12 * t126 + t124 * t13) * mrSges(7,3) + (m(5) * (t157 * t97 - t160 * t95) + (-t123 * t157 - t125 * t160) * mrSges(5,3)) * pkin(3) + m(6) * (t146 * t95 + (-t156 * t46 + t159 * t47) * t144) + (-mrSges(5,1) + t129) * t95 + (t47 * mrSges(6,3) + t125 * t191 + t144 * t85 + t197) * t159 + (-t46 * mrSges(6,3) - t125 * t131 / 0.2e1 - t144 * t86 + t196) * t156 + t79 * t52 + t80 * t51 + t63 * t88 - t73 * t195 + t74 * t194 - t97 * mrSges(5,2) + t35 * t193 + t36 * t192 + t128 * t42 + t146 * t82 + t173 * t123 + t201; t124 * t90 + t126 * t92 + 0.2e1 * t128 * t88 + 0.2e1 * t146 * t129 + t159 * t131 + t156 * t132 + m(7) * (t128 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(6) * (t144 ^ 2 * t175 + t146 ^ 2) + m(5) * (t157 ^ 2 + t160 ^ 2) * pkin(3) ^ 2 + t188 + 0.2e1 * (mrSges(5,1) * t160 - mrSges(5,2) * t157) * pkin(3) + 0.2e1 * (t124 * t80 - t126 * t79) * mrSges(7,3) + 0.2e1 * t175 * t144 * mrSges(6,3); t124 * t15 + t126 * t14 + t156 * t37 + t159 * t38 + m(7) * (t1 * t124 + t126 * t2) + m(6) * (t156 * t6 + t159 * t5) + m(5) * t81 + t43; t124 * t52 + t126 * t51 + t156 * t85 + t159 * t86 + m(7) * (t12 * t124 + t126 * t13) + m(6) * (t156 * t47 + t159 * t46) + m(5) * t147 + t87; m(7) * (t124 * t79 + t126 * t80); m(5) + m(6) * t175 + m(7) * (t124 ^ 2 + t126 ^ 2); m(6) * t21 + m(7) * t11 + t10 + t32; m(6) * t95 + m(7) * t63 + t42 + t82; m(6) * t146 + m(7) * t128 + t129 + t88; 0; m(6) + m(7); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t12 - mrSges(7,2) * t13 + t34; mrSges(7,1) * t79 - t80 * mrSges(7,2) + t89; -t88; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
