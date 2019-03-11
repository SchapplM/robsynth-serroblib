% Calculate joint inertia matrix for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR10_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR10_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR10_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:04:27
% EndTime: 2019-03-09 11:04:31
% DurationCPUTime: 1.71s
% Computational Cost: add. (2668->344), mult. (5894->472), div. (0->0), fcn. (6522->10), ass. (0->131)
t175 = Ifges(6,4) - Ifges(5,5);
t174 = Ifges(6,5) - Ifges(5,6);
t125 = sin(pkin(6));
t132 = cos(qJ(2));
t148 = t125 * t132;
t124 = sin(pkin(11));
t126 = cos(pkin(11));
t129 = sin(qJ(4));
t164 = cos(qJ(4));
t97 = t124 * t129 - t164 * t126;
t98 = t164 * t124 + t129 * t126;
t173 = t174 * t97 - t175 * t98;
t128 = sin(qJ(6));
t131 = cos(qJ(6));
t127 = cos(pkin(6));
t130 = sin(qJ(2));
t149 = t125 * t130;
t85 = -t124 * t149 + t126 * t127;
t86 = t124 * t127 + t126 * t149;
t55 = t129 * t86 - t164 * t85;
t32 = t128 * t148 + t131 * t55;
t172 = t32 / 0.2e1;
t152 = Ifges(7,4) * t131;
t37 = Ifges(7,5) * t98 + (Ifges(7,1) * t128 + t152) * t97;
t171 = t37 / 0.2e1;
t170 = pkin(4) + pkin(10);
t118 = Ifges(7,5) * t131;
t169 = -Ifges(7,6) * t128 / 0.2e1 + t118 / 0.2e1;
t168 = t128 / 0.2e1;
t167 = t131 / 0.2e1;
t146 = t128 ^ 2 + t131 ^ 2;
t99 = m(7) * t146;
t165 = m(6) + t99;
t163 = pkin(1) * t132;
t110 = pkin(8) * t149;
t87 = t127 * t163 - t110;
t162 = t87 * mrSges(3,1);
t88 = t127 * t130 * pkin(1) + pkin(8) * t148;
t161 = t88 * mrSges(3,2);
t159 = -mrSges(6,2) + mrSges(5,1);
t158 = Ifges(6,1) + Ifges(5,3);
t157 = pkin(9) + qJ(3);
t78 = qJ(3) * t127 + t88;
t79 = (-pkin(2) * t132 - qJ(3) * t130 - pkin(1)) * t125;
t42 = -t124 * t78 + t126 * t79;
t26 = -pkin(3) * t148 - pkin(9) * t86 + t42;
t43 = t124 * t79 + t126 * t78;
t29 = pkin(9) * t85 + t43;
t12 = t129 * t26 + t164 * t29;
t154 = mrSges(7,1) * t131;
t153 = Ifges(7,4) * t128;
t151 = t128 * t97;
t150 = t131 * t97;
t147 = t124 ^ 2 + t126 ^ 2;
t33 = t128 * t55 - t131 * t148;
t56 = t129 * t85 + t164 * t86;
t8 = Ifges(7,5) * t33 + Ifges(7,6) * t32 + Ifges(7,3) * t56;
t35 = Ifges(7,5) * t151 + Ifges(7,6) * t150 + Ifges(7,3) * t98;
t101 = t157 * t126;
t143 = t157 * t124;
t69 = t101 * t129 + t164 * t143;
t71 = t164 * t101 - t129 * t143;
t145 = t69 ^ 2 + t71 ^ 2;
t144 = Ifges(3,5) * t149 + Ifges(3,6) * t148 + Ifges(3,3) * t127;
t115 = -pkin(3) * t126 - pkin(2);
t58 = -t85 * mrSges(4,1) + t86 * mrSges(4,2);
t142 = t146 * mrSges(7,3);
t141 = -t174 * t55 + t175 * t56;
t100 = -t126 * mrSges(4,1) + t124 * mrSges(4,2);
t11 = -t129 * t29 + t164 * t26;
t7 = pkin(4) * t148 - t11;
t3 = t56 * pkin(5) + pkin(10) * t148 + t7;
t81 = t110 + (-pkin(2) - t163) * t127;
t57 = -pkin(3) * t85 + t81;
t136 = -qJ(5) * t56 + t57;
t5 = t170 * t55 + t136;
t1 = -t128 * t5 + t131 * t3;
t2 = t128 * t3 + t131 * t5;
t140 = t1 * t131 + t128 * t2;
t139 = -mrSges(7,2) * t128 + t154;
t137 = -qJ(5) * t98 + t115;
t34 = t170 * t97 + t137;
t46 = pkin(5) * t98 + t69;
t15 = -t128 * t34 + t131 * t46;
t16 = t128 * t46 + t131 * t34;
t138 = t128 * t16 + t131 * t15;
t6 = qJ(5) * t148 - t12;
t39 = t56 * mrSges(6,1) - mrSges(6,2) * t148;
t134 = qJ(5) ^ 2;
t107 = Ifges(7,1) * t131 - t153;
t106 = -Ifges(7,2) * t128 + t152;
t104 = mrSges(7,1) * t128 + mrSges(7,2) * t131;
t103 = Ifges(4,1) * t124 + Ifges(4,4) * t126;
t102 = Ifges(4,4) * t124 + Ifges(4,2) * t126;
t90 = t98 * mrSges(6,3);
t89 = t98 * mrSges(5,2);
t74 = -mrSges(4,1) * t148 - mrSges(4,3) * t86;
t73 = mrSges(4,2) * t148 + mrSges(4,3) * t85;
t68 = Ifges(5,1) * t98 - Ifges(5,4) * t97;
t67 = Ifges(5,4) * t98 - Ifges(5,2) * t97;
t66 = -Ifges(6,2) * t98 + Ifges(6,6) * t97;
t65 = -Ifges(6,6) * t98 + Ifges(6,3) * t97;
t64 = -t97 * mrSges(6,2) - t90;
t63 = t97 * mrSges(5,1) + t89;
t62 = -mrSges(7,2) * t98 + mrSges(7,3) * t150;
t61 = mrSges(7,1) * t98 - mrSges(7,3) * t151;
t60 = t139 * t97;
t59 = pkin(4) * t97 + t137;
t49 = t56 * mrSges(6,3);
t48 = t56 * mrSges(5,2);
t47 = -t97 * pkin(5) + t71;
t45 = Ifges(4,1) * t86 + Ifges(4,4) * t85 - Ifges(4,5) * t148;
t44 = Ifges(4,4) * t86 + Ifges(4,2) * t85 - Ifges(4,6) * t148;
t41 = -mrSges(5,1) * t148 - mrSges(5,3) * t56;
t40 = mrSges(5,2) * t148 - mrSges(5,3) * t55;
t38 = mrSges(6,1) * t55 + mrSges(6,3) * t148;
t36 = Ifges(7,6) * t98 + (Ifges(7,2) * t131 + t153) * t97;
t24 = -t55 * mrSges(6,2) - t49;
t23 = t55 * mrSges(5,1) + t48;
t22 = Ifges(5,1) * t56 - Ifges(5,4) * t55 - Ifges(5,5) * t148;
t21 = Ifges(5,4) * t56 - Ifges(5,2) * t55 - Ifges(5,6) * t148;
t20 = -Ifges(6,4) * t148 - Ifges(6,2) * t56 + Ifges(6,6) * t55;
t19 = -Ifges(6,5) * t148 - Ifges(6,6) * t56 + Ifges(6,3) * t55;
t18 = mrSges(7,1) * t56 - mrSges(7,3) * t33;
t17 = -mrSges(7,2) * t56 + mrSges(7,3) * t32;
t14 = -mrSges(7,1) * t32 + mrSges(7,2) * t33;
t13 = pkin(4) * t55 + t136;
t10 = Ifges(7,1) * t33 + Ifges(7,4) * t32 + Ifges(7,5) * t56;
t9 = Ifges(7,4) * t33 + Ifges(7,2) * t32 + Ifges(7,6) * t56;
t4 = -pkin(5) * t55 - t6;
t25 = [t86 * t45 + 0.2e1 * t43 * t73 + 0.2e1 * t42 * t74 + 0.2e1 * t81 * t58 + t85 * t44 + 0.2e1 * t57 * t23 + t32 * t9 + t33 * t10 + 0.2e1 * t6 * t38 + 0.2e1 * t7 * t39 + 0.2e1 * t12 * t40 + 0.2e1 * t11 * t41 + 0.2e1 * t2 * t17 + 0.2e1 * t1 * t18 + 0.2e1 * t13 * t24 + 0.2e1 * t4 * t14 + m(3) * (pkin(1) ^ 2 * t125 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(4) * (t42 ^ 2 + t43 ^ 2 + t81 ^ 2) + m(5) * (t11 ^ 2 + t12 ^ 2 + t57 ^ 2) + m(6) * (t13 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) + (t22 + t8 - t20) * t56 + (t19 - t21) * t55 + Ifges(2,3) + ((-0.2e1 * t87 * mrSges(3,3) + Ifges(3,5) * t127 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t130) * t125) * t130 + (0.2e1 * t88 * mrSges(3,3) - Ifges(4,5) * t86 + Ifges(3,6) * t127 - Ifges(4,6) * t85 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t130 + (Ifges(3,2) + Ifges(4,3) + t158) * t132) * t125 + t141) * t132) * t125 + (t144 - 0.2e1 * t161 + 0.2e1 * t162) * t127; m(5) * (-t11 * t69 + t115 * t57 + t12 * t71) + m(6) * (t13 * t59 - t6 * t71 + t69 * t7) + m(7) * (t1 * t15 + t16 * t2 + t4 * t47) + (t22 / 0.2e1 + t8 / 0.2e1 - t20 / 0.2e1 - t11 * mrSges(5,3) + t7 * mrSges(6,1)) * t98 + t162 - t161 + m(4) * (-pkin(2) * t81 + (-t124 * t42 + t126 * t43) * qJ(3)) + t115 * t23 + t81 * t100 + t85 * t102 / 0.2e1 + t86 * t103 / 0.2e1 + t59 * t24 - t4 * t60 + t1 * t61 + t2 * t62 + t57 * t63 + t13 * t64 - pkin(2) * t58 + t47 * t14 + t16 * t17 + t15 * t18 + (-t66 / 0.2e1 + t68 / 0.2e1 + t35 / 0.2e1) * t56 + (t65 / 0.2e1 - t67 / 0.2e1) * t55 + (t44 / 0.2e1 + qJ(3) * t73 + t43 * mrSges(4,3)) * t126 + (t45 / 0.2e1 - qJ(3) * t74 - t42 * mrSges(4,3)) * t124 + (t40 - t38) * t71 + (t39 - t41) * t69 + t144 + t33 * t171 + t36 * t172 + (t19 / 0.2e1 - t21 / 0.2e1 + t9 * t167 + t10 * t168 + t6 * mrSges(6,1) - t12 * mrSges(5,3)) * t97 - (Ifges(4,5) * t124 + Ifges(4,6) * t126 + t173) * t148 / 0.2e1; -0.2e1 * pkin(2) * t100 + t126 * t102 + t124 * t103 + 0.2e1 * t115 * t63 + 0.2e1 * t15 * t61 + 0.2e1 * t16 * t62 - 0.2e1 * t47 * t60 + 0.2e1 * t59 * t64 + Ifges(3,3) + (t35 - t66 + t68) * t98 + m(4) * (t147 * qJ(3) ^ 2 + pkin(2) ^ 2) + m(5) * (t115 ^ 2 + t145) + m(6) * (t59 ^ 2 + t145) + m(7) * (t15 ^ 2 + t16 ^ 2 + t47 ^ 2) + (t128 * t37 + t131 * t36 + t65 - t67) * t97 + 0.2e1 * t147 * qJ(3) * mrSges(4,3) + 0.2e1 * (t69 * t98 - t71 * t97) * (mrSges(6,1) + mrSges(5,3)); -t128 * t18 + t131 * t17 + t48 - t49 + t159 * t55 + m(7) * (-t1 * t128 + t131 * t2) + m(6) * t13 + m(5) * t57 + m(4) * t81 + t58; -m(4) * pkin(2) - t128 * t61 + t131 * t62 + t89 - t90 + t159 * t97 + m(7) * (-t128 * t15 + t131 * t16) + m(6) * t59 + m(5) * t115 + t100; m(4) + m(5) + t165; t4 * t104 + t56 * t169 + t106 * t172 + t33 * t107 / 0.2e1 - pkin(4) * t39 + t11 * mrSges(5,1) - t12 * mrSges(5,2) - t6 * mrSges(6,3) + t7 * mrSges(6,2) - t158 * t148 + (-t38 + t14) * qJ(5) + (t10 / 0.2e1 - t170 * t18 - t1 * mrSges(7,3)) * t131 + (-t9 / 0.2e1 - t170 * t17 - t2 * mrSges(7,3)) * t128 + m(7) * (qJ(5) * t4 - t140 * t170) + m(6) * (-pkin(4) * t7 - qJ(5) * t6) - t141; -qJ(5) * t60 + t47 * t104 + (-pkin(4) * mrSges(6,1) + t169) * t98 + (-mrSges(5,2) + mrSges(6,3)) * t71 - t159 * t69 + (-t15 * mrSges(7,3) - t170 * t61 + t171) * t131 + (-t36 / 0.2e1 - t170 * t62 - t16 * mrSges(7,3)) * t128 + m(7) * (qJ(5) * t47 - t138 * t170) + m(6) * (-pkin(4) * t69 + qJ(5) * t71) + (-qJ(5) * mrSges(6,1) + t106 * t167 + t107 * t168) * t97 + t173; 0; -0.2e1 * pkin(4) * mrSges(6,2) - t128 * t106 + t131 * t107 + m(7) * (t146 * t170 ^ 2 + t134) + m(6) * (pkin(4) ^ 2 + t134) + t158 + 0.2e1 * (t104 + mrSges(6,3)) * qJ(5) + 0.2e1 * t170 * t142; m(6) * t7 + m(7) * t140 + t128 * t17 + t131 * t18 + t39; m(6) * t69 + m(7) * t138 + t98 * mrSges(6,1) + t128 * t62 + t131 * t61; 0; -m(6) * pkin(4) - t170 * t99 + mrSges(6,2) - t142; t165; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t8; mrSges(7,1) * t15 - mrSges(7,2) * t16 + t35; -t104; -t170 * t154 + t118 + (mrSges(7,2) * t170 - Ifges(7,6)) * t128; t139; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
