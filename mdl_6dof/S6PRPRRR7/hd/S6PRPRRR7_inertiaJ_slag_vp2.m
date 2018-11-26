% Calculate joint inertia matrix for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:07
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_inertiaJ_slag_vp2: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:07:09
% EndTime: 2018-11-23 15:07:10
% DurationCPUTime: 1.44s
% Computational Cost: add. (3578->339), mult. (9838->521), div. (0->0), fcn. (11378->16), ass. (0->138)
t177 = 2 * pkin(11);
t176 = m(7) * pkin(12) + mrSges(7,3);
t117 = sin(pkin(8));
t121 = cos(pkin(8));
t126 = sin(qJ(4));
t130 = cos(qJ(4));
t122 = cos(pkin(7));
t118 = sin(pkin(7));
t120 = cos(pkin(14));
t152 = t118 * t120;
t142 = t121 * t152;
t116 = sin(pkin(14));
t166 = pkin(2) * t122;
t79 = qJ(3) * t152 + t116 * t166;
t53 = (t117 * t122 + t142) * pkin(10) + t79;
t101 = t120 * t166;
t155 = t116 * t118;
t61 = t122 * pkin(3) + t101 + (-pkin(10) * t121 - qJ(3)) * t155;
t66 = (-pkin(10) * t116 * t117 - pkin(3) * t120 - pkin(2)) * t118;
t27 = -t126 * t53 + (t117 * t66 + t121 * t61) * t130;
t124 = sin(qJ(6));
t128 = cos(qJ(6));
t90 = -t128 * mrSges(7,1) + t124 * mrSges(7,2);
t175 = -m(7) * pkin(5) - mrSges(6,1) + t90;
t125 = sin(qJ(5));
t129 = cos(qJ(5));
t119 = sin(pkin(6));
t127 = sin(qJ(2));
t131 = cos(qJ(2));
t149 = t122 * t131;
t123 = cos(pkin(6));
t151 = t118 * t123;
t57 = t120 * t151 + (-t116 * t127 + t120 * t149) * t119;
t158 = t121 * t57;
t58 = t120 * t119 * t127 + (t119 * t149 + t151) * t116;
t76 = -t118 * t119 * t131 + t123 * t122;
t31 = t130 * t58 + (t117 * t76 + t158) * t126;
t44 = -t117 * t57 + t121 * t76;
t14 = t125 * t31 - t129 * t44;
t174 = t14 ^ 2;
t153 = t117 * t130;
t29 = t126 * t58 - t130 * t158 - t76 * t153;
t173 = t29 ^ 2;
t154 = t117 * t126;
t80 = -t129 * t121 + t125 * t154;
t172 = t80 ^ 2;
t150 = t121 * t126;
t60 = t122 * t154 + (t116 * t130 + t120 * t150) * t118;
t75 = -t117 * t152 + t121 * t122;
t46 = t125 * t75 + t129 * t60;
t59 = -t122 * t153 + t126 * t155 - t130 * t142;
t34 = -t124 * t46 + t128 * t59;
t171 = t34 / 0.2e1;
t35 = t124 * t59 + t128 * t46;
t170 = t35 / 0.2e1;
t169 = -t124 / 0.2e1;
t168 = t124 / 0.2e1;
t167 = t128 / 0.2e1;
t165 = pkin(11) * t129;
t164 = t80 * t14;
t91 = -t129 * mrSges(6,1) + t125 * mrSges(6,2);
t163 = mrSges(5,1) - t91;
t13 = -t34 * mrSges(7,1) + t35 * mrSges(7,2);
t37 = t59 * mrSges(6,1) - t46 * mrSges(6,3);
t162 = t13 - t37;
t39 = -t117 * t61 + t121 * t66;
t22 = t59 * pkin(4) - t60 * pkin(11) + t39;
t28 = t130 * t53 + t61 * t150 + t66 * t154;
t25 = t75 * pkin(11) + t28;
t8 = t125 * t22 + t129 * t25;
t45 = t125 * t60 - t129 * t75;
t26 = t45 * mrSges(6,1) + t46 * mrSges(6,2);
t48 = t75 * mrSges(5,1) - t60 * mrSges(5,3);
t161 = t26 - t48;
t160 = Ifges(7,4) * t124;
t159 = Ifges(7,4) * t128;
t157 = t14 * t125;
t156 = t80 * t125;
t148 = t124 * t125;
t147 = t125 * t128;
t92 = Ifges(7,5) * t124 + Ifges(7,6) * t128;
t146 = Ifges(6,5) * t125 + Ifges(6,6) * t129;
t145 = t124 ^ 2 + t128 ^ 2;
t9 = Ifges(7,5) * t35 + Ifges(7,6) * t34 + Ifges(7,3) * t45;
t144 = Ifges(6,5) * t46 - Ifges(6,6) * t45 + Ifges(6,3) * t59;
t143 = Ifges(5,5) * t60 - Ifges(5,6) * t59 + Ifges(5,3) * t75;
t24 = -t75 * pkin(4) - t27;
t12 = t45 * pkin(5) - t46 * pkin(12) + t24;
t4 = t59 * pkin(12) + t8;
t1 = t128 * t12 - t124 * t4;
t2 = t124 * t12 + t128 * t4;
t141 = -t1 * t124 + t2 * t128;
t139 = mrSges(7,1) * t124 + mrSges(7,2) * t128;
t89 = -t129 * pkin(5) - t125 * pkin(12) - pkin(4);
t68 = -t124 * t165 + t128 * t89;
t69 = t124 * t89 + t128 * t165;
t136 = -t68 * t124 + t69 * t128;
t16 = t125 * t44 + t129 * t31;
t135 = t16 * t129 + t157;
t7 = -t125 * t25 + t129 * t22;
t82 = t125 * t121 + t129 * t154;
t134 = t82 * t129 + t156;
t71 = Ifges(7,5) * t147 - Ifges(7,6) * t148 - Ifges(7,3) * t129;
t133 = pkin(11) ^ 2;
t115 = t129 ^ 2;
t113 = t125 ^ 2;
t110 = t117 ^ 2;
t109 = t113 * t133;
t103 = t110 * t130 ^ 2;
t99 = mrSges(4,2) * t155;
t96 = Ifges(6,1) * t125 + Ifges(6,4) * t129;
t95 = Ifges(7,1) * t124 + t159;
t94 = Ifges(6,4) * t125 + Ifges(6,2) * t129;
t93 = Ifges(7,2) * t128 + t160;
t88 = -t129 * mrSges(7,1) - mrSges(7,3) * t147;
t87 = t129 * mrSges(7,2) - mrSges(7,3) * t148;
t85 = -t122 * mrSges(4,2) + mrSges(4,3) * t152;
t84 = t122 * mrSges(4,1) - mrSges(4,3) * t155;
t83 = t139 * t125;
t78 = -mrSges(4,1) * t152 + t99;
t77 = -qJ(3) * t155 + t101;
t73 = -Ifges(7,5) * t129 + (Ifges(7,1) * t128 - t160) * t125;
t72 = -Ifges(7,6) * t129 + (-Ifges(7,2) * t124 + t159) * t125;
t64 = -t124 * t153 + t128 * t82;
t63 = -t124 * t82 - t128 * t153;
t47 = -t75 * mrSges(5,2) - t59 * mrSges(5,3);
t38 = t59 * mrSges(5,1) + t60 * mrSges(5,2);
t36 = -t59 * mrSges(6,2) - t45 * mrSges(6,3);
t20 = Ifges(6,1) * t46 - Ifges(6,4) * t45 + Ifges(6,5) * t59;
t19 = Ifges(6,4) * t46 - Ifges(6,2) * t45 + Ifges(6,6) * t59;
t18 = t45 * mrSges(7,1) - t35 * mrSges(7,3);
t17 = -t45 * mrSges(7,2) + t34 * mrSges(7,3);
t11 = Ifges(7,1) * t35 + Ifges(7,4) * t34 + Ifges(7,5) * t45;
t10 = Ifges(7,4) * t35 + Ifges(7,2) * t34 + Ifges(7,6) * t45;
t6 = t124 * t29 + t128 * t16;
t5 = -t124 * t16 + t128 * t29;
t3 = -t59 * pkin(5) - t7;
t15 = [m(2) + m(7) * (t5 ^ 2 + t6 ^ 2 + t174) + m(6) * (t16 ^ 2 + t173 + t174) + m(5) * (t31 ^ 2 + t44 ^ 2 + t173) + m(4) * (t57 ^ 2 + t58 ^ 2 + t76 ^ 2) + m(3) * (t123 ^ 2 + (t127 ^ 2 + t131 ^ 2) * t119 ^ 2); t16 * t36 + t6 * t17 + t5 * t18 + t31 * t47 + t44 * t38 + t57 * t84 + t58 * t85 + t76 * t78 + t161 * t29 + t162 * t14 + (t131 * mrSges(3,1) - t127 * mrSges(3,2)) * t119 + m(7) * (t1 * t5 + t3 * t14 + t2 * t6) + m(6) * (-t7 * t14 + t8 * t16 + t24 * t29) + m(5) * (-t27 * t29 + t28 * t31 + t39 * t44) + m(4) * (-t118 * pkin(2) * t76 + t77 * t57 + t79 * t58); m(5) * (t27 ^ 2 + t28 ^ 2 + t39 ^ 2) + m(6) * (t24 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + (t9 - t19) * t45 + Ifges(4,3) * t122 ^ 2 + t60 * (Ifges(5,1) * t60 + Ifges(5,5) * t75) + m(4) * (t77 ^ 2 + t79 ^ 2) + t75 * t143 + (-0.2e1 * Ifges(5,4) * t60 + Ifges(5,2) * t59 - Ifges(5,6) * t75 + t144) * t59 + 0.2e1 * t77 * t84 + 0.2e1 * t79 * t85 + Ifges(3,3) + 0.2e1 * t3 * t13 + 0.2e1 * t2 * t17 + 0.2e1 * t1 * t18 + 0.2e1 * t24 * t26 + t34 * t10 + t35 * t11 + 0.2e1 * t8 * t36 + 0.2e1 * t7 * t37 + 0.2e1 * t39 * t38 + t46 * t20 + 0.2e1 * t28 * t47 + 0.2e1 * t27 * t48 + (-0.2e1 * pkin(2) * t78 + 0.2e1 * t122 * (Ifges(4,5) * t116 + Ifges(4,6) * t120) + (m(4) * pkin(2) ^ 2 + Ifges(4,2) * t120 ^ 2 + (Ifges(4,1) * t116 + 0.2e1 * Ifges(4,4) * t120) * t116) * t118) * t118; m(7) * (t63 * t5 + t64 * t6 + t164) + m(6) * (-t29 * t153 + t82 * t16 + t164) + m(5) * (t121 * t44 + (t126 * t31 - t130 * t29) * t117) + m(4) * t76; t121 * t38 + t64 * t17 + t63 * t18 + t82 * t36 + t99 + t162 * t80 + (-m(4) * pkin(2) - t120 * mrSges(4,1)) * t118 + (t126 * t47 - t161 * t130) * t117 + m(7) * (t63 * t1 + t64 * t2 + t80 * t3) + m(6) * (-t24 * t153 - t80 * t7 + t82 * t8) + m(5) * (t121 * t39 + (t126 * t28 + t130 * t27) * t117); m(4) + m(7) * (t63 ^ 2 + t64 ^ 2 + t172) + m(6) * (t82 ^ 2 + t103 + t172) + m(5) * (t110 * t126 ^ 2 + t121 ^ 2 + t103); -t31 * mrSges(5,2) + t14 * t83 + t5 * t88 + t6 * t87 - t163 * t29 + t135 * mrSges(6,3) + m(7) * (pkin(11) * t157 + t68 * t5 + t69 * t6) + m(6) * (-pkin(4) * t29 + t135 * pkin(11)); m(6) * (-pkin(4) * t24 + (-t7 * t125 + t8 * t129) * pkin(11)) + t143 + (t20 / 0.2e1 - t7 * mrSges(6,3) + t10 * t169 + t11 * t167 + t162 * pkin(11)) * t125 + (-t9 / 0.2e1 + t19 / 0.2e1 + pkin(11) * t36 + t8 * mrSges(6,3)) * t129 + m(7) * (t125 * pkin(11) * t3 + t68 * t1 + t69 * t2) + (-t94 / 0.2e1 + t71 / 0.2e1) * t45 + t3 * t83 + t2 * t87 + t1 * t88 + t24 * t91 + t46 * t96 / 0.2e1 + t59 * t146 / 0.2e1 - pkin(4) * t26 + t27 * mrSges(5,1) - t28 * mrSges(5,2) + t68 * t18 + t69 * t17 + t72 * t171 + t73 * t170; t63 * t88 + t64 * t87 + t80 * t83 + t134 * mrSges(6,3) + (-t126 * mrSges(5,2) + t163 * t130) * t117 + m(7) * (pkin(11) * t156 + t68 * t63 + t69 * t64) + m(6) * (pkin(4) * t153 + t134 * pkin(11)); -0.2e1 * pkin(4) * t91 + 0.2e1 * t68 * t88 + 0.2e1 * t69 * t87 + Ifges(5,3) + (-t71 + t94) * t129 + (t113 + t115) * mrSges(6,3) * t177 + m(7) * (t68 ^ 2 + t69 ^ 2 + t109) + m(6) * (pkin(4) ^ 2 + t115 * t133 + t109) + (-t124 * t72 + t128 * t73 + t83 * t177 + t96) * t125; -t16 * mrSges(6,2) + t176 * (-t5 * t124 + t6 * t128) + t175 * t14; t95 * t170 + t93 * t171 + t45 * t92 / 0.2e1 + t11 * t168 + t3 * t90 + t10 * t167 - t8 * mrSges(6,2) + t7 * mrSges(6,1) + t141 * mrSges(7,3) + t144 + (-m(7) * t3 - t13) * pkin(5) + (m(7) * t141 - t124 * t18 + t128 * t17) * pkin(12); -t82 * mrSges(6,2) + t176 * (-t63 * t124 + t64 * t128) + t175 * t80; t73 * t168 + t72 * t167 - pkin(5) * t83 + (m(7) * t136 - t124 * t88 + t128 * t87) * pkin(12) + (t175 * pkin(11) + t95 * t167 + t93 * t169) * t125 + (-t92 / 0.2e1 - pkin(11) * mrSges(6,2)) * t129 + t136 * mrSges(7,3) + t146; Ifges(6,3) - 0.2e1 * pkin(5) * t90 + t124 * t95 + t128 * t93 + m(7) * (t145 * pkin(12) ^ 2 + pkin(5) ^ 2) + 0.2e1 * t145 * pkin(12) * mrSges(7,3); t5 * mrSges(7,1) - t6 * mrSges(7,2); t1 * mrSges(7,1) - t2 * mrSges(7,2) + t9; t63 * mrSges(7,1) - t64 * mrSges(7,2); t68 * mrSges(7,1) - t69 * mrSges(7,2) + t71; -t139 * pkin(12) + t92; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
