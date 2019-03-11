% Calculate joint inertia matrix for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR10_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR10_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:18:24
% EndTime: 2019-03-09 14:18:28
% DurationCPUTime: 1.82s
% Computational Cost: add. (4797->410), mult. (10671->599), div. (0->0), fcn. (12271->12), ass. (0->149)
t157 = sin(pkin(6));
t166 = cos(qJ(2));
t182 = t157 * t166;
t158 = cos(pkin(12));
t191 = pkin(9) + qJ(3);
t130 = t191 * t158;
t162 = sin(qJ(4));
t156 = sin(pkin(12));
t176 = t191 * t156;
t193 = cos(qJ(4));
t95 = t130 * t162 + t193 * t176;
t205 = t95 ^ 2;
t204 = 0.2e1 * t95;
t161 = sin(qJ(5));
t165 = cos(qJ(5));
t159 = cos(pkin(6));
t163 = sin(qJ(2));
t183 = t157 * t163;
t112 = -t156 * t183 + t158 * t159;
t113 = t156 * t159 + t158 * t183;
t76 = t162 * t112 + t113 * t193;
t56 = -t161 * t76 - t165 * t182;
t57 = -t161 * t182 + t165 * t76;
t75 = -t112 * t193 + t113 * t162;
t23 = Ifges(6,1) * t57 + Ifges(6,4) * t56 + Ifges(6,5) * t75;
t203 = t23 / 0.2e1;
t125 = t156 * t162 - t158 * t193;
t126 = t156 * t193 + t162 * t158;
t188 = Ifges(6,4) * t165;
t59 = Ifges(6,6) * t125 + (-Ifges(6,2) * t161 + t188) * t126;
t202 = t59 / 0.2e1;
t189 = Ifges(6,4) * t161;
t60 = Ifges(6,5) * t125 + (Ifges(6,1) * t165 - t189) * t126;
t201 = t60 / 0.2e1;
t160 = sin(qJ(6));
t164 = cos(qJ(6));
t127 = -t160 * t161 + t164 * t165;
t128 = t160 * t165 + t161 * t164;
t93 = Ifges(7,4) * t128 + Ifges(7,2) * t127;
t200 = t93 / 0.2e1;
t94 = Ifges(7,1) * t128 + Ifges(7,4) * t127;
t199 = t94 / 0.2e1;
t198 = -pkin(11) - pkin(10);
t197 = t127 / 0.2e1;
t196 = t128 / 0.2e1;
t136 = Ifges(6,1) * t161 + t188;
t195 = t136 / 0.2e1;
t192 = pkin(1) * t166;
t115 = t159 * t163 * pkin(1) + pkin(8) * t182;
t106 = qJ(3) * t159 + t115;
t107 = (-pkin(2) * t166 - qJ(3) * t163 - pkin(1)) * t157;
t63 = -t106 * t156 + t158 * t107;
t46 = -pkin(3) * t182 - pkin(9) * t113 + t63;
t64 = t158 * t106 + t156 * t107;
t51 = pkin(9) * t112 + t64;
t25 = t162 * t46 + t193 * t51;
t20 = -pkin(10) * t182 + t25;
t141 = pkin(8) * t183;
t109 = t141 + (-pkin(2) - t192) * t159;
t81 = -pkin(3) * t112 + t109;
t28 = pkin(4) * t75 - pkin(10) * t76 + t81;
t7 = t161 * t28 + t165 * t20;
t190 = -Ifges(5,5) * t76 + Ifges(5,6) * t75;
t146 = -pkin(3) * t158 - pkin(2);
t84 = pkin(4) * t125 - pkin(10) * t126 + t146;
t97 = t130 * t193 - t162 * t176;
t48 = t161 * t84 + t165 * t97;
t114 = t159 * t192 - t141;
t187 = t114 * mrSges(3,1);
t186 = t115 * mrSges(3,2);
t185 = t126 * t161;
t184 = t126 * t165;
t181 = Ifges(5,5) * t126 - Ifges(5,6) * t125;
t92 = Ifges(7,5) * t128 + Ifges(7,6) * t127;
t134 = Ifges(6,5) * t161 + Ifges(6,6) * t165;
t180 = t156 ^ 2 + t158 ^ 2;
t179 = t161 ^ 2 + t165 ^ 2;
t31 = -t160 * t57 + t164 * t56;
t32 = t160 * t56 + t164 * t57;
t8 = Ifges(7,5) * t32 + Ifges(7,6) * t31 + Ifges(7,3) * t75;
t21 = Ifges(6,5) * t57 + Ifges(6,6) * t56 + Ifges(6,3) * t75;
t77 = t128 * t126;
t78 = t127 * t126;
t37 = Ifges(7,5) * t78 - Ifges(7,6) * t77 + Ifges(7,3) * t125;
t178 = t92 / 0.2e1 + t134 / 0.2e1;
t177 = Ifges(3,5) * t183 + Ifges(3,6) * t182 + Ifges(3,3) * t159;
t43 = t75 * mrSges(5,1) + t76 * mrSges(5,2);
t6 = -t161 * t20 + t165 * t28;
t47 = -t161 * t97 + t165 * t84;
t82 = -t112 * mrSges(4,1) + t113 * mrSges(4,2);
t129 = -t158 * mrSges(4,1) + t156 * mrSges(4,2);
t88 = t125 * mrSges(5,1) + t126 * mrSges(5,2);
t91 = -t127 * mrSges(7,1) + t128 * mrSges(7,2);
t24 = -t162 * t51 + t193 * t46;
t133 = -t165 * mrSges(6,1) + t161 * mrSges(6,2);
t175 = mrSges(6,1) * t161 + mrSges(6,2) * t165;
t174 = -t63 * t156 + t64 * t158;
t137 = t198 * t161;
t138 = t198 * t165;
t101 = t137 * t164 + t138 * t160;
t102 = t137 * t160 - t138 * t164;
t173 = t101 * mrSges(7,1) - t102 * mrSges(7,2) + t92;
t19 = pkin(4) * t182 - t24;
t58 = Ifges(6,5) * t184 - Ifges(6,6) * t185 + Ifges(6,3) * t125;
t4 = pkin(5) * t75 - pkin(11) * t57 + t6;
t5 = pkin(11) * t56 + t7;
t2 = -t160 * t5 + t164 * t4;
t3 = t160 * t4 + t164 * t5;
t172 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t8;
t34 = pkin(5) * t125 - pkin(11) * t184 + t47;
t42 = -pkin(11) * t185 + t48;
t14 = -t160 * t42 + t164 * t34;
t15 = t160 * t34 + t164 * t42;
t171 = t14 * mrSges(7,1) - t15 * mrSges(7,2) + t37;
t170 = (mrSges(7,1) * t164 - mrSges(7,2) * t160) * pkin(5);
t147 = -pkin(5) * t165 - pkin(4);
t135 = Ifges(6,2) * t165 + t189;
t132 = Ifges(4,1) * t156 + Ifges(4,4) * t158;
t131 = Ifges(4,4) * t156 + Ifges(4,2) * t158;
t100 = -mrSges(4,1) * t182 - mrSges(4,3) * t113;
t99 = mrSges(4,2) * t182 + mrSges(4,3) * t112;
t90 = Ifges(5,1) * t126 - Ifges(5,4) * t125;
t89 = Ifges(5,4) * t126 - Ifges(5,2) * t125;
t86 = mrSges(6,1) * t125 - mrSges(6,3) * t184;
t85 = -mrSges(6,2) * t125 - mrSges(6,3) * t185;
t83 = t175 * t126;
t67 = Ifges(4,1) * t113 + Ifges(4,4) * t112 - Ifges(4,5) * t182;
t66 = Ifges(4,4) * t113 + Ifges(4,2) * t112 - Ifges(4,6) * t182;
t65 = pkin(5) * t185 + t95;
t62 = -mrSges(5,1) * t182 - mrSges(5,3) * t76;
t61 = mrSges(5,2) * t182 - mrSges(5,3) * t75;
t53 = mrSges(7,1) * t125 - mrSges(7,3) * t78;
t52 = -mrSges(7,2) * t125 - mrSges(7,3) * t77;
t44 = mrSges(7,1) * t77 + mrSges(7,2) * t78;
t41 = Ifges(5,1) * t76 - Ifges(5,4) * t75 - Ifges(5,5) * t182;
t40 = Ifges(5,4) * t76 - Ifges(5,2) * t75 - Ifges(5,6) * t182;
t39 = Ifges(7,1) * t78 - Ifges(7,4) * t77 + Ifges(7,5) * t125;
t38 = Ifges(7,4) * t78 - Ifges(7,2) * t77 + Ifges(7,6) * t125;
t36 = mrSges(6,1) * t75 - mrSges(6,3) * t57;
t35 = -mrSges(6,2) * t75 + mrSges(6,3) * t56;
t33 = -mrSges(6,1) * t56 + mrSges(6,2) * t57;
t22 = Ifges(6,4) * t57 + Ifges(6,2) * t56 + Ifges(6,6) * t75;
t17 = mrSges(7,1) * t75 - mrSges(7,3) * t32;
t16 = -mrSges(7,2) * t75 + mrSges(7,3) * t31;
t12 = -t56 * pkin(5) + t19;
t11 = -mrSges(7,1) * t31 + mrSges(7,2) * t32;
t10 = Ifges(7,1) * t32 + Ifges(7,4) * t31 + Ifges(7,5) * t75;
t9 = Ifges(7,4) * t32 + Ifges(7,2) * t31 + Ifges(7,6) * t75;
t1 = [m(7) * (t12 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t19 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(5) * (t24 ^ 2 + t25 ^ 2 + t81 ^ 2) + m(4) * (t109 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(3) * (pkin(1) ^ 2 * t157 ^ 2 + t114 ^ 2 + t115 ^ 2) + (t8 + t21 - t40) * t75 + ((-0.2e1 * t114 * mrSges(3,3) + Ifges(3,5) * t159 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t163) * t157) * t163 + (0.2e1 * t115 * mrSges(3,3) - Ifges(4,5) * t113 + Ifges(3,6) * t159 - Ifges(4,6) * t112 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t163 + (Ifges(5,3) + Ifges(4,3) + Ifges(3,2)) * t166) * t157 + t190) * t166) * t157 + (t177 - 0.2e1 * t186 + 0.2e1 * t187) * t159 + Ifges(2,3) + 0.2e1 * t12 * t11 + 0.2e1 * t3 * t16 + 0.2e1 * t2 * t17 + t31 * t9 + t32 * t10 + 0.2e1 * t19 * t33 + 0.2e1 * t7 * t35 + 0.2e1 * t6 * t36 + t56 * t22 + t57 * t23 + 0.2e1 * t25 * t61 + 0.2e1 * t24 * t62 + t76 * t41 + 0.2e1 * t81 * t43 + 0.2e1 * t64 * t99 + 0.2e1 * t63 * t100 + 0.2e1 * t109 * t82 + t112 * t66 + t113 * t67; m(4) * (-pkin(2) * t109 + qJ(3) * t174) + t174 * mrSges(4,3) + m(7) * (t12 * t65 + t14 * t2 + t15 * t3) + m(6) * (t19 * t95 + t47 * t6 + t48 * t7) + m(5) * (t146 * t81 - t24 * t95 + t25 * t97) + (-t24 * mrSges(5,3) - t161 * t22 / 0.2e1 + t165 * t203 + t41 / 0.2e1) * t126 + (t37 / 0.2e1 + t58 / 0.2e1 - t89 / 0.2e1) * t75 + (t33 - t62) * t95 - (Ifges(4,5) * t156 + Ifges(4,6) * t158 + t181) * t182 / 0.2e1 + t177 + t57 * t201 + t56 * t202 - t186 + t187 + (-t156 * t100 + t158 * t99) * qJ(3) + t15 * t16 + t14 * t17 + t31 * t38 / 0.2e1 + t32 * t39 / 0.2e1 + t12 * t44 + t47 * t36 + t48 * t35 + t3 * t52 + t2 * t53 + t65 * t11 - t77 * t9 / 0.2e1 + t78 * t10 / 0.2e1 - pkin(2) * t82 + t19 * t83 + t7 * t85 + t6 * t86 + t81 * t88 + t76 * t90 / 0.2e1 + t97 * t61 + t109 * t129 + t112 * t131 / 0.2e1 + t113 * t132 / 0.2e1 + t146 * t43 + t156 * t67 / 0.2e1 + t158 * t66 / 0.2e1 + (-t25 * mrSges(5,3) + t8 / 0.2e1 + t21 / 0.2e1 - t40 / 0.2e1) * t125; -0.2e1 * pkin(2) * t129 + t158 * t131 + t156 * t132 + 0.2e1 * t14 * t53 + 0.2e1 * t146 * t88 + 0.2e1 * t15 * t52 - t77 * t38 + t78 * t39 + 0.2e1 * t65 * t44 + 0.2e1 * t47 * t86 + 0.2e1 * t48 * t85 + t83 * t204 + Ifges(3,3) + 0.2e1 * t180 * qJ(3) * mrSges(4,3) + (mrSges(5,3) * t204 - t161 * t59 + t165 * t60 + t90) * t126 + (-0.2e1 * mrSges(5,3) * t97 + t37 + t58 - t89) * t125 + m(7) * (t14 ^ 2 + t15 ^ 2 + t65 ^ 2) + m(6) * (t47 ^ 2 + t48 ^ 2 + t205) + m(5) * (t146 ^ 2 + t97 ^ 2 + t205) + m(4) * (qJ(3) ^ 2 * t180 + pkin(2) ^ 2); t127 * t17 + t128 * t16 + t161 * t35 + t165 * t36 + m(7) * (t127 * t2 + t128 * t3) + m(6) * (t161 * t7 + t165 * t6) + m(5) * t81 + m(4) * t109 + t82 + t43; -m(4) * pkin(2) + t127 * t53 + t128 * t52 + t161 * t85 + t165 * t86 + m(7) * (t127 * t14 + t128 * t15) + m(6) * (t161 * t48 + t165 * t47) + m(5) * t146 + t88 + t129; m(4) + m(5) + m(6) * t179 + m(7) * (t127 ^ 2 + t128 ^ 2); t178 * t75 - Ifges(5,3) * t182 + (pkin(10) * t35 + t7 * mrSges(6,3) + t22 / 0.2e1) * t165 + (-t6 * mrSges(6,3) - pkin(10) * t36 + t203) * t161 + m(7) * (t101 * t2 + t102 * t3 + t12 * t147) + (t127 * t3 - t128 * t2) * mrSges(7,3) - t190 + m(6) * (-pkin(4) * t19 + (-t161 * t6 + t165 * t7) * pkin(10)) + t24 * mrSges(5,1) - t25 * mrSges(5,2) - pkin(4) * t33 + t12 * t91 + t31 * t200 + t32 * t199 + t101 * t17 + t102 * t16 + t9 * t197 + t10 * t196 + t19 * t133 + t56 * t135 / 0.2e1 + t57 * t195 + t147 * t11; -pkin(4) * t83 + t65 * t91 - t77 * t200 + t78 * t199 - t97 * mrSges(5,2) + t101 * t53 + t102 * t52 + t38 * t197 + t39 * t196 + t147 * t44 + (-mrSges(5,1) + t133) * t95 + t178 * t125 + (t127 * t15 - t128 * t14) * mrSges(7,3) + (t48 * mrSges(6,3) + pkin(10) * t85 + t126 * t195 + t202) * t165 + (-t47 * mrSges(6,3) - pkin(10) * t86 - t126 * t135 / 0.2e1 + t201) * t161 + m(7) * (t101 * t14 + t102 * t15 + t147 * t65) + m(6) * (-pkin(4) * t95 + (-t161 * t47 + t165 * t48) * pkin(10)) + t181; m(7) * (t101 * t127 + t102 * t128); -0.2e1 * pkin(4) * t133 + t127 * t93 + t128 * t94 + t165 * t135 + t161 * t136 + 0.2e1 * t147 * t91 + Ifges(5,3) + m(7) * (t101 ^ 2 + t102 ^ 2 + t147 ^ 2) + m(6) * (pkin(10) ^ 2 * t179 + pkin(4) ^ 2) + 0.2e1 * (-t101 * t128 + t102 * t127) * mrSges(7,3) + 0.2e1 * t179 * pkin(10) * mrSges(6,3); t6 * mrSges(6,1) - t7 * mrSges(6,2) + (m(7) * (t160 * t3 + t164 * t2) + t160 * t16 + t164 * t17) * pkin(5) + t172 + t21; t47 * mrSges(6,1) - t48 * mrSges(6,2) + (m(7) * (t14 * t164 + t15 * t160) + t160 * t52 + t164 * t53) * pkin(5) + t171 + t58; m(7) * (t127 * t164 + t128 * t160) * pkin(5) - t133 - t91; -t175 * pkin(10) + (m(7) * (t101 * t164 + t102 * t160) + (t127 * t160 - t128 * t164) * mrSges(7,3)) * pkin(5) + t173 + t134; Ifges(6,3) + Ifges(7,3) + m(7) * (t160 ^ 2 + t164 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t170; t172; t171; -t91; t173; Ifges(7,3) + t170; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
