% Calculate joint inertia matrix for
% S5RRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR11_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR11_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR11_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR11_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR11_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:38:50
% EndTime: 2019-12-31 22:38:54
% DurationCPUTime: 1.27s
% Computational Cost: add. (2071->326), mult. (4725->481), div. (0->0), fcn. (5009->10), ass. (0->134)
t172 = 2 * pkin(8);
t126 = sin(qJ(4));
t130 = cos(qJ(4));
t123 = sin(pkin(5));
t132 = cos(qJ(2));
t148 = t123 * t132;
t124 = cos(pkin(5));
t127 = sin(qJ(3));
t131 = cos(qJ(3));
t128 = sin(qJ(2));
t149 = t123 * t128;
t80 = t124 * t127 + t131 * t149;
t50 = -t126 * t80 - t130 * t148;
t51 = -t126 * t148 + t130 * t80;
t79 = -t124 * t131 + t127 * t149;
t17 = Ifges(5,1) * t51 + Ifges(5,4) * t50 + Ifges(5,5) * t79;
t171 = t17 / 0.2e1;
t125 = sin(qJ(5));
t129 = cos(qJ(5));
t23 = -t125 * t51 + t129 * t50;
t170 = t23 / 0.2e1;
t24 = t125 * t50 + t129 * t51;
t169 = t24 / 0.2e1;
t151 = Ifges(5,4) * t130;
t71 = -Ifges(5,6) * t131 + (-Ifges(5,2) * t126 + t151) * t127;
t168 = t71 / 0.2e1;
t152 = Ifges(5,4) * t126;
t72 = -Ifges(5,5) * t131 + (Ifges(5,1) * t130 - t152) * t127;
t167 = t72 / 0.2e1;
t89 = t125 * t130 + t126 * t129;
t73 = t89 * t127;
t166 = -t73 / 0.2e1;
t88 = -t125 * t126 + t129 * t130;
t74 = t88 * t127;
t165 = t74 / 0.2e1;
t164 = t88 / 0.2e1;
t163 = t89 / 0.2e1;
t99 = Ifges(5,1) * t126 + t151;
t162 = t99 / 0.2e1;
t161 = -pkin(10) - pkin(9);
t160 = pkin(1) * t132;
t159 = pkin(8) * t127;
t158 = pkin(8) * t131;
t105 = pkin(7) * t149;
t81 = t124 * t160 - t105;
t157 = t81 * mrSges(3,1);
t82 = t124 * t128 * pkin(1) + pkin(7) * t148;
t156 = t82 * mrSges(3,2);
t155 = -Ifges(6,3) - Ifges(5,3);
t65 = t105 + (-pkin(2) - t160) * t124;
t29 = pkin(3) * t79 - pkin(9) * t80 + t65;
t66 = pkin(8) * t124 + t82;
t67 = (-pkin(2) * t132 - pkin(8) * t128 - pkin(1)) * t123;
t35 = t127 * t67 + t131 * t66;
t31 = -pkin(9) * t148 + t35;
t11 = t126 * t29 + t130 * t31;
t154 = Ifges(6,5) * t74 - Ifges(6,6) * t73;
t153 = -Ifges(4,5) * t80 + Ifges(4,6) * t79;
t47 = Ifges(6,5) * t89 + Ifges(6,6) * t88;
t150 = Ifges(6,3) * t131;
t93 = -pkin(3) * t131 - pkin(9) * t127 - pkin(2);
t64 = t126 * t93 + t130 * t158;
t147 = t126 * t127;
t146 = t127 * t130;
t96 = Ifges(5,5) * t126 + Ifges(5,6) * t130;
t145 = Ifges(4,5) * t127 + Ifges(4,6) * t131;
t144 = t126 ^ 2 + t130 ^ 2;
t5 = Ifges(6,5) * t24 + Ifges(6,6) * t23 + Ifges(6,3) * t79;
t15 = Ifges(5,5) * t51 + Ifges(5,6) * t50 + Ifges(5,3) * t79;
t143 = t96 / 0.2e1 + t47 / 0.2e1;
t142 = Ifges(3,5) * t149 + Ifges(3,6) * t148 + Ifges(3,3) * t124;
t10 = -t126 * t31 + t130 * t29;
t34 = -t127 * t66 + t131 * t67;
t30 = pkin(3) * t148 - t34;
t141 = Ifges(5,5) * t146 - Ifges(5,6) * t147;
t140 = mrSges(5,1) * t126 + mrSges(5,2) * t130;
t87 = t130 * t93;
t43 = -pkin(10) * t146 + t87 + (-pkin(8) * t126 - pkin(4)) * t131;
t52 = -pkin(10) * t147 + t64;
t19 = -t125 * t52 + t129 * t43;
t20 = t125 * t43 + t129 * t52;
t139 = t19 * mrSges(6,1) - t20 * mrSges(6,2) + t154;
t101 = t161 * t126;
t102 = t161 * t130;
t56 = t101 * t129 + t102 * t125;
t57 = t101 * t125 - t102 * t129;
t138 = t56 * mrSges(6,1) - t57 * mrSges(6,2) + t47;
t4 = pkin(4) * t79 - pkin(10) * t51 + t10;
t8 = pkin(10) * t50 + t11;
t2 = -t125 * t8 + t129 * t4;
t3 = t125 * t4 + t129 * t8;
t137 = t2 * mrSges(6,1) - t3 * mrSges(6,2) + t5;
t136 = (mrSges(6,1) * t129 - mrSges(6,2) * t125) * pkin(4);
t134 = pkin(8) ^ 2;
t122 = t131 ^ 2;
t120 = t127 ^ 2;
t118 = t120 * t134;
t112 = -pkin(4) * t130 - pkin(3);
t100 = Ifges(4,1) * t127 + Ifges(4,4) * t131;
t98 = Ifges(4,4) * t127 + Ifges(4,2) * t131;
t97 = Ifges(5,2) * t130 + t152;
t95 = -mrSges(4,1) * t131 + mrSges(4,2) * t127;
t94 = -mrSges(5,1) * t130 + mrSges(5,2) * t126;
t92 = (pkin(4) * t126 + pkin(8)) * t127;
t91 = -mrSges(5,1) * t131 - mrSges(5,3) * t146;
t90 = mrSges(5,2) * t131 - mrSges(5,3) * t147;
t83 = t140 * t127;
t70 = -Ifges(5,3) * t131 + t141;
t63 = -t126 * t158 + t87;
t59 = -mrSges(6,1) * t131 - mrSges(6,3) * t74;
t58 = mrSges(6,2) * t131 - mrSges(6,3) * t73;
t55 = -mrSges(4,1) * t148 - mrSges(4,3) * t80;
t54 = mrSges(4,2) * t148 - mrSges(4,3) * t79;
t49 = Ifges(6,1) * t89 + Ifges(6,4) * t88;
t48 = Ifges(6,4) * t89 + Ifges(6,2) * t88;
t46 = -mrSges(6,1) * t88 + mrSges(6,2) * t89;
t42 = mrSges(4,1) * t79 + mrSges(4,2) * t80;
t41 = mrSges(6,1) * t73 + mrSges(6,2) * t74;
t40 = Ifges(4,1) * t80 - Ifges(4,4) * t79 - Ifges(4,5) * t148;
t39 = Ifges(4,4) * t80 - Ifges(4,2) * t79 - Ifges(4,6) * t148;
t38 = Ifges(6,1) * t74 - Ifges(6,4) * t73 - Ifges(6,5) * t131;
t37 = Ifges(6,4) * t74 - Ifges(6,2) * t73 - Ifges(6,6) * t131;
t36 = -t150 + t154;
t33 = mrSges(5,1) * t79 - mrSges(5,3) * t51;
t32 = -mrSges(5,2) * t79 + mrSges(5,3) * t50;
t25 = -mrSges(5,1) * t50 + mrSges(5,2) * t51;
t16 = Ifges(5,4) * t51 + Ifges(5,2) * t50 + Ifges(5,6) * t79;
t14 = -pkin(4) * t50 + t30;
t13 = mrSges(6,1) * t79 - mrSges(6,3) * t24;
t12 = -mrSges(6,2) * t79 + mrSges(6,3) * t23;
t9 = -mrSges(6,1) * t23 + mrSges(6,2) * t24;
t7 = Ifges(6,1) * t24 + Ifges(6,4) * t23 + Ifges(6,5) * t79;
t6 = Ifges(6,4) * t24 + Ifges(6,2) * t23 + Ifges(6,6) * t79;
t1 = [0.2e1 * t10 * t33 + 0.2e1 * t11 * t32 + 0.2e1 * t3 * t12 + 0.2e1 * t2 * t13 + 0.2e1 * t14 * t9 + t50 * t16 + t51 * t17 + t23 * t6 + t24 * t7 + 0.2e1 * t30 * t25 + 0.2e1 * t34 * t55 + 0.2e1 * t35 * t54 + t80 * t40 + 0.2e1 * t65 * t42 + Ifges(2,3) + (t15 + t5 - t39) * t79 + (t142 - 0.2e1 * t156 + 0.2e1 * t157) * t124 + ((-0.2e1 * t81 * mrSges(3,3) + Ifges(3,5) * t124 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t128) * t123) * t128 + (0.2e1 * t82 * mrSges(3,3) + Ifges(3,6) * t124 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t128 + (Ifges(3,2) + Ifges(4,3)) * t132) * t123 + t153) * t132) * t123 + m(3) * (pkin(1) ^ 2 * t123 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(4) * (t34 ^ 2 + t35 ^ 2 + t65 ^ 2) + m(5) * (t10 ^ 2 + t11 ^ 2 + t30 ^ 2) + m(6) * (t14 ^ 2 + t2 ^ 2 + t3 ^ 2); (t39 / 0.2e1 - t15 / 0.2e1 - t5 / 0.2e1 + pkin(8) * t54 + t35 * mrSges(4,3)) * t131 + m(6) * (t14 * t92 + t19 * t2 + t20 * t3) + (-t98 / 0.2e1 + t70 / 0.2e1 + t36 / 0.2e1) * t79 + t80 * t100 / 0.2e1 + t30 * t83 + t11 * t90 + t10 * t91 + t92 * t9 + t65 * t95 + t63 * t33 + t64 * t32 + t3 * t58 + t2 * t59 + t14 * t41 - pkin(2) * t42 - t145 * t148 / 0.2e1 + t142 + t19 * t13 + t20 * t12 + t7 * t165 + t6 * t166 + t51 * t167 + t50 * t168 + t38 * t169 + t37 * t170 + m(4) * (-pkin(2) * t65 + (-t34 * t127 + t35 * t131) * pkin(8)) + m(5) * (t10 * t63 + t11 * t64 + t30 * t159) + (t40 / 0.2e1 + t130 * t171 - t126 * t16 / 0.2e1 - t34 * mrSges(4,3) + (-t55 + t25) * pkin(8)) * t127 - t156 + t157; -0.2e1 * pkin(2) * t95 + 0.2e1 * t19 * t59 + 0.2e1 * t20 * t58 - t73 * t37 + t74 * t38 + 0.2e1 * t92 * t41 + 0.2e1 * t63 * t91 + 0.2e1 * t64 * t90 + Ifges(3,3) + (t120 + t122) * mrSges(4,3) * t172 + (-t36 - t70 + t98) * t131 + m(6) * (t19 ^ 2 + t20 ^ 2 + t92 ^ 2) + m(5) * (t63 ^ 2 + t64 ^ 2 + t118) + m(4) * (pkin(2) ^ 2 + t122 * t134 + t118) + (-t126 * t71 + t130 * t72 + t83 * t172 + t100) * t127; t50 * t97 / 0.2e1 + t51 * t162 + t112 * t9 + t6 * t164 + t7 * t163 + t30 * t94 + t48 * t170 + t49 * t169 + t56 * t13 + t57 * t12 + t14 * t46 + t143 * t79 - pkin(3) * t25 + t34 * mrSges(4,1) - t35 * mrSges(4,2) + m(5) * (-pkin(3) * t30 + (-t10 * t126 + t11 * t130) * pkin(9)) - Ifges(4,3) * t148 + m(6) * (t112 * t14 + t2 * t56 + t3 * t57) - t153 + (-t10 * mrSges(5,3) - pkin(9) * t33 + t171) * t126 + (t16 / 0.2e1 + pkin(9) * t32 + t11 * mrSges(5,3)) * t130 + (-t2 * t89 + t3 * t88) * mrSges(6,3); t112 * t41 - pkin(3) * t83 + t37 * t164 + t38 * t163 + t92 * t46 + t48 * t166 + t49 * t165 + t57 * t58 + t56 * t59 + (-mrSges(4,1) + t94) * t159 + (-t19 * t89 + t20 * t88) * mrSges(6,3) + (t64 * mrSges(5,3) + pkin(9) * t90 + t127 * t162 + t168) * t130 + (t167 - t127 * t97 / 0.2e1 - pkin(9) * t91 - t63 * mrSges(5,3)) * t126 + m(5) * (-pkin(3) * t159 + (-t126 * t63 + t130 * t64) * pkin(9)) + m(6) * (t112 * t92 + t19 * t56 + t20 * t57) + (-pkin(8) * mrSges(4,2) - t143) * t131 + t145; -0.2e1 * pkin(3) * t94 + 0.2e1 * t112 * t46 + t126 * t99 + t130 * t97 + t88 * t48 + t89 * t49 + Ifges(4,3) + m(6) * (t112 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(5) * (t144 * pkin(9) ^ 2 + pkin(3) ^ 2) + 0.2e1 * (-t56 * t89 + t57 * t88) * mrSges(6,3) + 0.2e1 * t144 * pkin(9) * mrSges(5,3); t10 * mrSges(5,1) - t11 * mrSges(5,2) + (t129 * t13 + m(6) * (t125 * t3 + t129 * t2) + t125 * t12) * pkin(4) + t137 + t15; t63 * mrSges(5,1) - t64 * mrSges(5,2) + t155 * t131 + (t125 * t58 + m(6) * (t125 * t20 + t129 * t19) + t129 * t59) * pkin(4) + t139 + t141; -t140 * pkin(9) + (m(6) * (t125 * t57 + t129 * t56) + (t125 * t88 - t129 * t89) * mrSges(6,3)) * pkin(4) + t138 + t96; m(6) * (t125 ^ 2 + t129 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t136 - t155; t137; t139 - t150; t138; Ifges(6,3) + t136; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
