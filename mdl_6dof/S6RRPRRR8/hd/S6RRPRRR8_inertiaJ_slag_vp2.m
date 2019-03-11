% Calculate joint inertia matrix for
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:02:23
% EndTime: 2019-03-09 14:02:26
% DurationCPUTime: 1.28s
% Computational Cost: add. (3561->321), mult. (7254->462), div. (0->0), fcn. (8061->10), ass. (0->131)
t169 = 2 * pkin(7);
t133 = sin(qJ(5));
t137 = cos(qJ(5));
t130 = sin(pkin(11));
t131 = cos(pkin(11));
t134 = sin(qJ(4));
t138 = cos(qJ(4));
t106 = t130 * t138 + t131 * t134;
t135 = sin(qJ(2));
t91 = t106 * t135;
t105 = -t130 * t134 + t131 * t138;
t92 = t105 * t135;
t60 = -t133 * t92 - t137 * t91;
t61 = -t133 * t91 + t137 * t92;
t168 = -Ifges(6,5) * t61 - Ifges(6,6) * t60;
t167 = -Ifges(5,5) * t92 + Ifges(5,6) * t91;
t166 = -t130 / 0.2e1;
t165 = t131 / 0.2e1;
t164 = pkin(4) * t133;
t139 = cos(qJ(2));
t163 = pkin(7) * t139;
t125 = t135 * pkin(7);
t121 = pkin(4) * t137 + pkin(5);
t132 = sin(qJ(6));
t136 = cos(qJ(6));
t95 = t121 * t132 + t136 * t164;
t162 = t95 * mrSges(7,2);
t161 = -Ifges(7,3) - Ifges(6,3);
t160 = pkin(8) + qJ(3);
t29 = -t132 * t61 + t136 * t60;
t30 = t132 * t60 + t136 * t61;
t159 = -Ifges(7,5) * t30 - Ifges(7,6) * t29;
t110 = -pkin(2) * t139 - qJ(3) * t135 - pkin(1);
t101 = t131 * t110;
t154 = t131 * t135;
t72 = -pkin(8) * t154 + t101 + (-pkin(7) * t130 - pkin(3)) * t139;
t155 = t130 * t135;
t84 = t130 * t110 + t131 * t163;
t78 = -pkin(8) * t155 + t84;
t47 = -t134 * t78 + t138 * t72;
t32 = -pkin(4) * t139 - pkin(9) * t92 + t47;
t48 = t134 * t72 + t138 * t78;
t37 = -pkin(9) * t91 + t48;
t14 = t133 * t32 + t137 * t37;
t111 = t160 * t130;
t113 = t160 * t131;
t79 = -t138 * t111 - t113 * t134;
t62 = -pkin(9) * t106 + t79;
t80 = -t134 * t111 + t138 * t113;
t63 = pkin(9) * t105 + t80;
t36 = t133 * t62 + t137 * t63;
t158 = Ifges(4,4) * t130;
t157 = Ifges(4,4) * t131;
t156 = t132 * mrSges(7,2);
t96 = mrSges(4,1) * t155 + mrSges(4,2) * t154;
t109 = pkin(3) * t155 + t125;
t153 = t130 ^ 2 + t131 ^ 2;
t152 = pkin(5) * t156;
t151 = -Ifges(5,3) + t161;
t120 = -pkin(3) * t131 - pkin(2);
t64 = t91 * mrSges(5,1) + t92 * mrSges(5,2);
t33 = -t60 * mrSges(6,1) + t61 * mrSges(6,2);
t70 = t105 * t137 - t106 * t133;
t71 = t105 * t133 + t106 * t137;
t44 = -t70 * mrSges(6,1) + t71 * mrSges(6,2);
t11 = -t29 * mrSges(7,1) + t30 * mrSges(7,2);
t42 = -t132 * t71 + t136 * t70;
t43 = t132 * t70 + t136 * t71;
t15 = -t42 * mrSges(7,1) + t43 * mrSges(7,2);
t75 = -t105 * mrSges(5,1) + t106 * mrSges(5,2);
t13 = -t133 * t37 + t137 * t32;
t35 = -t133 * t63 + t137 * t62;
t112 = -t131 * mrSges(4,1) + t130 * mrSges(4,2);
t94 = t121 * t136 - t132 * t164;
t93 = t94 * mrSges(7,1);
t150 = Ifges(7,3) + t93 - t162;
t74 = pkin(4) * t91 + t109;
t7 = -pkin(5) * t139 - pkin(10) * t61 + t13;
t8 = pkin(10) * t60 + t14;
t2 = -t132 * t8 + t136 * t7;
t3 = t132 * t7 + t136 * t8;
t149 = t2 * mrSges(7,1) - t3 * mrSges(7,2) - t159;
t40 = Ifges(7,6) * t42;
t41 = Ifges(7,5) * t43;
t18 = -pkin(10) * t71 + t35;
t19 = pkin(10) * t70 + t36;
t5 = -t132 * t19 + t136 * t18;
t6 = t132 * t18 + t136 * t19;
t148 = t5 * mrSges(7,1) - t6 * mrSges(7,2) + t40 + t41;
t83 = -t130 * t163 + t101;
t147 = -t130 * t83 + t131 * t84;
t85 = -pkin(4) * t105 + t120;
t146 = (mrSges(6,1) * t137 - mrSges(6,2) * t133) * pkin(4);
t145 = t13 * mrSges(6,1) - t14 * mrSges(6,2) + t149 - t168;
t68 = Ifges(6,6) * t70;
t69 = Ifges(6,5) * t71;
t144 = t35 * mrSges(6,1) - t36 * mrSges(6,2) + t148 + t68 + t69;
t141 = pkin(7) ^ 2;
t129 = t139 ^ 2;
t128 = t135 ^ 2;
t124 = t128 * t141;
t122 = t136 * pkin(5) * mrSges(7,1);
t115 = Ifges(4,1) * t130 + t157;
t114 = Ifges(4,2) * t131 + t158;
t108 = -mrSges(4,1) * t139 - mrSges(4,3) * t154;
t107 = mrSges(4,2) * t139 - mrSges(4,3) * t155;
t99 = Ifges(5,5) * t106;
t98 = Ifges(5,6) * t105;
t90 = -Ifges(4,5) * t139 + (Ifges(4,1) * t131 - t158) * t135;
t89 = -Ifges(4,6) * t139 + (-Ifges(4,2) * t130 + t157) * t135;
t82 = -mrSges(5,1) * t139 - mrSges(5,3) * t92;
t81 = mrSges(5,2) * t139 - mrSges(5,3) * t91;
t77 = Ifges(5,1) * t106 + Ifges(5,4) * t105;
t76 = Ifges(5,4) * t106 + Ifges(5,2) * t105;
t59 = Ifges(5,1) * t92 - Ifges(5,4) * t91 - Ifges(5,5) * t139;
t58 = Ifges(5,4) * t92 - Ifges(5,2) * t91 - Ifges(5,6) * t139;
t51 = -mrSges(6,1) * t139 - mrSges(6,3) * t61;
t50 = mrSges(6,2) * t139 + mrSges(6,3) * t60;
t49 = -pkin(5) * t70 + t85;
t46 = Ifges(6,1) * t71 + Ifges(6,4) * t70;
t45 = Ifges(6,4) * t71 + Ifges(6,2) * t70;
t38 = -pkin(5) * t60 + t74;
t26 = Ifges(6,1) * t61 + Ifges(6,4) * t60 - Ifges(6,5) * t139;
t25 = Ifges(6,4) * t61 + Ifges(6,2) * t60 - Ifges(6,6) * t139;
t21 = -mrSges(7,1) * t139 - mrSges(7,3) * t30;
t20 = mrSges(7,2) * t139 + mrSges(7,3) * t29;
t17 = Ifges(7,1) * t43 + Ifges(7,4) * t42;
t16 = Ifges(7,4) * t43 + Ifges(7,2) * t42;
t10 = Ifges(7,1) * t30 + Ifges(7,4) * t29 - Ifges(7,5) * t139;
t9 = Ifges(7,4) * t30 + Ifges(7,2) * t29 - Ifges(7,6) * t139;
t1 = [t29 * t9 + t30 * t10 + 0.2e1 * t3 * t20 + 0.2e1 * t2 * t21 + m(6) * (t13 ^ 2 + t14 ^ 2 + t74 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t38 ^ 2) + m(5) * (t109 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(4) * (t83 ^ 2 + t84 ^ 2 + t124) + m(3) * (pkin(1) ^ 2 + t129 * t141 + t124) + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t135 - t130 * t89 + t131 * t90 + t169 * t96) * t135 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + Ifges(4,3) - t151) * t139 + (-Ifges(4,5) * t131 + Ifges(4,6) * t130 + (2 * Ifges(3,4))) * t135 + t159 + t167 + t168) * t139 + 0.2e1 * t38 * t11 + Ifges(2,3) + (t128 + t129) * mrSges(3,3) * t169 + 0.2e1 * t14 * t50 + 0.2e1 * t13 * t51 + t60 * t25 + t61 * t26 + 0.2e1 * t74 * t33 + 0.2e1 * t48 * t81 + 0.2e1 * t47 * t82 - t91 * t58 + t92 * t59 + 0.2e1 * t84 * t107 + 0.2e1 * t83 * t108 + 0.2e1 * t109 * t64; t29 * t16 / 0.2e1 + t30 * t17 / 0.2e1 + t6 * t20 + t5 * t21 + (t114 * t166 + t115 * t165 + Ifges(3,5) + (t112 - mrSges(3,1)) * pkin(7)) * t135 + (-pkin(7) * mrSges(3,2) - t41 / 0.2e1 - t40 / 0.2e1 - t69 / 0.2e1 - t68 / 0.2e1 - t99 / 0.2e1 - t98 / 0.2e1 + Ifges(4,5) * t166 - Ifges(4,6) * t131 / 0.2e1 + Ifges(3,6)) * t139 + m(4) * (-pkin(2) * t125 + qJ(3) * t147) + (-t2 * t43 + t3 * t42) * mrSges(7,3) + (-t13 * t71 + t14 * t70) * mrSges(6,3) + (t105 * t48 - t106 * t47) * mrSges(5,3) + t147 * mrSges(4,3) + t89 * t165 + (t131 * t107 - t130 * t108) * qJ(3) + m(7) * (t2 * t5 + t3 * t6 + t38 * t49) + m(6) * (t13 * t35 + t14 * t36 + t74 * t85) + m(5) * (t109 * t120 + t47 * t79 + t48 * t80) + t38 * t15 + t42 * t9 / 0.2e1 + t43 * t10 / 0.2e1 + t49 * t11 + t36 * t50 + t35 * t51 + t60 * t45 / 0.2e1 + t61 * t46 / 0.2e1 + t70 * t25 / 0.2e1 + t71 * t26 / 0.2e1 + t74 * t44 + t80 * t81 + t79 * t82 + t85 * t33 - t91 * t76 / 0.2e1 + t92 * t77 / 0.2e1 - pkin(2) * t96 + t105 * t58 / 0.2e1 + t106 * t59 / 0.2e1 + t109 * t75 + t120 * t64 + t130 * t90 / 0.2e1; -0.2e1 * pkin(2) * t112 + t105 * t76 + t106 * t77 + t131 * t114 + t130 * t115 + 0.2e1 * t120 * t75 + 0.2e1 * t49 * t15 + t42 * t16 + t43 * t17 + 0.2e1 * t85 * t44 + t70 * t45 + t71 * t46 + Ifges(3,3) + m(7) * (t49 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(6) * (t35 ^ 2 + t36 ^ 2 + t85 ^ 2) + m(5) * (t120 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(4) * (qJ(3) ^ 2 * t153 + pkin(2) ^ 2) + 0.2e1 * (t42 * t6 - t43 * t5) * mrSges(7,3) + 0.2e1 * (-t35 * t71 + t36 * t70) * mrSges(6,3) + 0.2e1 * (t105 * t80 - t106 * t79) * mrSges(5,3) + 0.2e1 * t153 * qJ(3) * mrSges(4,3); m(4) * t125 + m(5) * t109 + m(6) * t74 + m(7) * t38 + t11 + t33 + t64 + t96; -m(4) * pkin(2) + m(5) * t120 + m(6) * t85 + m(7) * t49 + t112 + t15 + t44 + t75; m(4) + m(5) + m(6) + m(7); t145 + (t133 * t50 + t137 * t51 + m(6) * (t13 * t137 + t133 * t14)) * pkin(4) + t47 * mrSges(5,1) - t48 * mrSges(5,2) + t94 * t21 + t95 * t20 + t151 * t139 + m(7) * (t2 * t94 + t3 * t95) - t167; m(7) * (t5 * t94 + t6 * t95) + t79 * mrSges(5,1) - t80 * mrSges(5,2) + t99 + t98 + (t42 * t95 - t43 * t94) * mrSges(7,3) + (m(6) * (t133 * t36 + t137 * t35) + (t133 * t70 - t137 * t71) * mrSges(6,3)) * pkin(4) + t144; 0; -0.2e1 * t162 + 0.2e1 * t93 + 0.2e1 * t146 + m(7) * (t94 ^ 2 + t95 ^ 2) + m(6) * (t133 ^ 2 + t137 ^ 2) * pkin(4) ^ 2 - t151; t161 * t139 + (t136 * t21 + m(7) * (t132 * t3 + t136 * t2) + t132 * t20) * pkin(5) + t145; (m(7) * (t132 * t6 + t136 * t5) + (t132 * t42 - t136 * t43) * mrSges(7,3)) * pkin(5) + t144; 0; Ifges(6,3) + t122 + t146 + (m(7) * (t132 * t95 + t136 * t94) - t156) * pkin(5) + t150; -0.2e1 * t152 + 0.2e1 * t122 + m(7) * (t132 ^ 2 + t136 ^ 2) * pkin(5) ^ 2 - t161; -Ifges(7,3) * t139 + t149; t148; 0; t150; Ifges(7,3) + t122 - t152; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
