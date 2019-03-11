% Calculate joint inertia matrix for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:22:54
% EndTime: 2019-03-09 09:22:57
% DurationCPUTime: 1.20s
% Computational Cost: add. (1588->300), mult. (3147->417), div. (0->0), fcn. (3154->8), ass. (0->113)
t109 = sin(pkin(10));
t110 = cos(pkin(10));
t143 = t109 ^ 2 + t110 ^ 2;
t142 = 2 * pkin(7);
t116 = cos(qJ(2));
t113 = sin(qJ(2));
t112 = sin(qJ(5));
t115 = cos(qJ(5));
t73 = t109 * t115 - t110 * t112;
t60 = t73 * t113;
t122 = t109 * t112 + t110 * t115;
t61 = t122 * t113;
t141 = Ifges(6,5) * t61 + Ifges(6,6) * t60 + Ifges(6,3) * t116;
t140 = pkin(3) + pkin(4);
t139 = pkin(7) * t116;
t138 = -pkin(8) + qJ(3);
t104 = t116 * pkin(3);
t81 = -pkin(2) * t116 - qJ(3) * t113 - pkin(1);
t95 = t109 * t139;
t32 = pkin(4) * t116 + t104 + t95 + (-pkin(8) * t113 - t81) * t110;
t131 = t109 * t113;
t52 = t109 * t81 + t110 * t139;
t47 = -qJ(4) * t116 + t52;
t38 = pkin(8) * t131 + t47;
t13 = t112 * t32 + t115 * t38;
t82 = t138 * t109;
t85 = t138 * t110;
t44 = t112 * t82 + t115 * t85;
t130 = t110 * t113;
t64 = mrSges(4,1) * t131 + mrSges(4,2) * t130;
t137 = t143 * qJ(3) ^ 2;
t136 = Ifges(4,4) * t109;
t135 = Ifges(4,4) * t110;
t134 = Ifges(5,5) * t109;
t133 = Ifges(5,5) * t110;
t132 = t109 * mrSges(5,3);
t78 = t116 * mrSges(5,1) + mrSges(5,2) * t130;
t111 = sin(qJ(6));
t114 = cos(qJ(6));
t25 = -t111 * t61 + t114 * t60;
t26 = t111 * t60 + t114 * t61;
t129 = Ifges(7,5) * t26 + Ifges(7,6) * t25 + Ifges(7,3) * t116;
t74 = -t111 * t112 + t114 * t115;
t75 = t111 * t115 + t112 * t114;
t128 = t74 * mrSges(7,1) - t75 * mrSges(7,2);
t127 = qJ(4) * t109 + pkin(2);
t51 = t110 * t81 - t95;
t12 = -t112 * t38 + t115 * t32;
t43 = -t112 * t85 + t115 * t82;
t29 = -t60 * mrSges(6,1) + t61 * mrSges(6,2);
t39 = mrSges(6,1) * t122 + t73 * mrSges(6,2);
t9 = -t25 * mrSges(7,1) + t26 * mrSges(7,2);
t36 = -t111 * t73 - t114 * t122;
t37 = -t111 * t122 + t114 * t73;
t14 = -t36 * mrSges(7,1) + t37 * mrSges(7,2);
t27 = -pkin(9) * t73 + t43;
t28 = -pkin(9) * t122 + t44;
t10 = -t111 * t28 + t114 * t27;
t11 = t111 * t27 + t114 * t28;
t34 = Ifges(7,6) * t36;
t35 = Ifges(7,5) * t37;
t125 = t10 * mrSges(7,1) - t11 * mrSges(7,2) + t34 + t35;
t48 = t104 - t51;
t124 = t48 * t109 + t47 * t110;
t123 = -t51 * t109 + t52 * t110;
t4 = pkin(5) * t116 - pkin(9) * t61 + t12;
t5 = pkin(9) * t60 + t13;
t2 = -t111 * t5 + t114 * t4;
t3 = t111 * t4 + t114 * t5;
t121 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t129;
t120 = (mrSges(7,1) * t114 - mrSges(7,2) * t111) * pkin(5);
t62 = t140 * t110 + t127;
t90 = qJ(4) * t130;
t45 = -t90 - (-t140 * t109 - pkin(7)) * t113;
t118 = pkin(7) ^ 2;
t108 = t116 ^ 2;
t107 = t113 ^ 2;
t103 = t107 * t118;
t99 = t109 * mrSges(4,2);
t91 = mrSges(5,1) * t131;
t89 = Ifges(4,1) * t109 + t135;
t88 = Ifges(5,1) * t109 - t133;
t87 = Ifges(4,2) * t110 + t136;
t86 = -Ifges(5,3) * t110 + t134;
t84 = -t110 * mrSges(4,1) + t99;
t83 = -t110 * mrSges(5,1) - t132;
t80 = -pkin(3) * t110 - t127;
t79 = -mrSges(5,2) * t131 - mrSges(5,3) * t116;
t77 = -mrSges(4,1) * t116 - mrSges(4,3) * t130;
t76 = mrSges(4,2) * t116 - mrSges(4,3) * t131;
t66 = Ifges(6,5) * t73;
t65 = Ifges(6,6) * t122;
t63 = -mrSges(5,3) * t130 + t91;
t59 = -t90 + (pkin(3) * t109 + pkin(7)) * t113;
t58 = -Ifges(4,5) * t116 + (Ifges(4,1) * t110 - t136) * t113;
t57 = -Ifges(5,4) * t116 + (Ifges(5,1) * t110 + t134) * t113;
t56 = -Ifges(4,6) * t116 + (-Ifges(4,2) * t109 + t135) * t113;
t55 = -Ifges(5,6) * t116 + (Ifges(5,3) * t109 + t133) * t113;
t50 = mrSges(6,1) * t116 - mrSges(6,3) * t61;
t49 = -mrSges(6,2) * t116 + mrSges(6,3) * t60;
t42 = pkin(5) * t122 + t62;
t41 = Ifges(6,1) * t73 - Ifges(6,4) * t122;
t40 = Ifges(6,4) * t73 - Ifges(6,2) * t122;
t24 = Ifges(6,1) * t61 + Ifges(6,4) * t60 + Ifges(6,5) * t116;
t23 = Ifges(6,4) * t61 + Ifges(6,2) * t60 + Ifges(6,6) * t116;
t19 = pkin(5) * t60 + t45;
t18 = mrSges(7,1) * t116 - mrSges(7,3) * t26;
t17 = -mrSges(7,2) * t116 + mrSges(7,3) * t25;
t16 = Ifges(7,1) * t37 + Ifges(7,4) * t36;
t15 = Ifges(7,4) * t37 + Ifges(7,2) * t36;
t7 = Ifges(7,1) * t26 + Ifges(7,4) * t25 + Ifges(7,5) * t116;
t6 = Ifges(7,4) * t26 + Ifges(7,2) * t25 + Ifges(7,6) * t116;
t1 = [0.2e1 * t12 * t50 + 0.2e1 * t13 * t49 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 - 0.2e1 * t19 * t9 + t60 * t23 + t61 * t24 + t25 * t6 + t26 * t7 - 0.2e1 * t45 * t29 + 0.2e1 * t47 * t79 + 0.2e1 * t48 * t78 + 0.2e1 * t51 * t77 + 0.2e1 * t52 * t76 + 0.2e1 * t59 * t63 + Ifges(2,3) + (t107 + t108) * mrSges(3,3) * t142 + m(3) * (pkin(1) ^ 2 + t108 * t118 + t103) + m(4) * (t51 ^ 2 + t52 ^ 2 + t103) + m(5) * (t47 ^ 2 + t48 ^ 2 + t59 ^ 2) + m(6) * (t12 ^ 2 + t13 ^ 2 + t45 ^ 2) + m(7) * (t19 ^ 2 + t2 ^ 2 + t3 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(5,2) + Ifges(3,2)) * t116 + t129 + t141) * t116 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t113 + t64 * t142 + (t57 + t58) * t110 + (t55 - t56) * t109 + ((2 * Ifges(3,4)) + (-Ifges(5,4) - Ifges(4,5)) * t110 + (Ifges(4,6) - Ifges(5,6)) * t109) * t116) * t113; (-t55 / 0.2e1 + t56 / 0.2e1) * t110 + (t57 / 0.2e1 + t58 / 0.2e1) * t109 + (-t2 * t37 + t3 * t36) * mrSges(7,3) - t122 * t23 / 0.2e1 + (-t12 * t73 - t122 * t13) * mrSges(6,3) + (-pkin(7) * mrSges(3,2) + t66 / 0.2e1 - t65 / 0.2e1 + t35 / 0.2e1 + t34 / 0.2e1 + Ifges(3,6) + (-Ifges(4,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t110 + (-Ifges(4,5) / 0.2e1 - Ifges(5,4) / 0.2e1) * t109) * t116 + t80 * t63 + t59 * t83 - pkin(2) * t64 + t73 * t24 / 0.2e1 + t60 * t40 / 0.2e1 + t61 * t41 / 0.2e1 + t62 * t29 + t37 * t7 / 0.2e1 + t42 * t9 - t45 * t39 + t44 * t49 + t43 * t50 + t25 * t15 / 0.2e1 + t26 * t16 / 0.2e1 + t36 * t6 / 0.2e1 + t11 * t17 + t10 * t18 - t19 * t14 + (Ifges(3,5) + (t88 / 0.2e1 + t89 / 0.2e1) * t110 + (t86 / 0.2e1 - t87 / 0.2e1) * t109 + (t84 - mrSges(3,1)) * pkin(7)) * t113 + ((t76 + t79) * t110 + (-t77 + t78) * t109) * qJ(3) + m(5) * (t124 * qJ(3) + t59 * t80) + t123 * mrSges(4,3) + m(4) * (-pkin(2) * pkin(7) * t113 + t123 * qJ(3)) + t124 * mrSges(5,2) + m(6) * (t12 * t43 + t13 * t44 - t45 * t62) + m(7) * (t10 * t2 + t11 * t3 - t19 * t42); -0.2e1 * pkin(2) * t84 + 0.2e1 * t42 * t14 + t36 * t15 + t37 * t16 + 0.2e1 * t62 * t39 - t122 * t40 + t73 * t41 + 0.2e1 * t80 * t83 + Ifges(3,3) + (-t86 + t87) * t110 + (t88 + t89) * t109 + m(5) * (t80 ^ 2 + t137) + m(4) * (pkin(2) ^ 2 + t137) + m(6) * (t43 ^ 2 + t44 ^ 2 + t62 ^ 2) + m(7) * (t10 ^ 2 + t11 ^ 2 + t42 ^ 2) + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * qJ(3) * t143 + 0.2e1 * (-t10 * t37 + t11 * t36) * mrSges(7,3) + 0.2e1 * (-t122 * t44 - t43 * t73) * mrSges(6,3); t91 + (m(4) * pkin(7) - mrSges(5,3) * t110) * t113 + m(5) * t59 + m(6) * t45 + m(7) * t19 - t9 - t29 + t64; -m(4) * pkin(2) - t132 + t99 + (-mrSges(5,1) - mrSges(4,1)) * t110 + m(5) * t80 - m(6) * t62 - m(7) * t42 - t14 - t39; m(4) + m(5) + m(6) + m(7); t112 * t49 + t115 * t50 + t75 * t17 + t74 * t18 + m(7) * (t2 * t74 + t3 * t75) + m(6) * (t112 * t13 + t115 * t12) + m(5) * t48 + t78; (m(5) * qJ(3) + mrSges(5,2)) * t109 + (t36 * t75 - t37 * t74) * mrSges(7,3) + (-t112 * t122 - t115 * t73) * mrSges(6,3) + m(7) * (t10 * t74 + t11 * t75) + m(6) * (t112 * t44 + t115 * t43); 0; m(5) + m(6) * (t112 ^ 2 + t115 ^ 2) + m(7) * (t74 ^ 2 + t75 ^ 2); t12 * mrSges(6,1) - t13 * mrSges(6,2) + (m(7) * (t111 * t3 + t114 * t2) + t111 * t17 + t114 * t18) * pkin(5) + t121 + t141; t43 * mrSges(6,1) - t44 * mrSges(6,2) - t65 + t66 + (m(7) * (t10 * t114 + t11 * t111) + (t111 * t36 - t114 * t37) * mrSges(7,3)) * pkin(5) + t125; 0; t115 * mrSges(6,1) - t112 * mrSges(6,2) + m(7) * (t111 * t75 + t114 * t74) * pkin(5) + t128; Ifges(6,3) + Ifges(7,3) + m(7) * (t111 ^ 2 + t114 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t120; t121; t125; 0; t128; Ifges(7,3) + t120; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
