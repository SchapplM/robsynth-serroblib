% Calculate joint inertia matrix for
% S6RPRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2018-11-23 16:17
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:17:24
% EndTime: 2018-11-23 16:17:25
% DurationCPUTime: 0.92s
% Computational Cost: add. (2706->244), mult. (5131->350), div. (0->0), fcn. (6053->10), ass. (0->96)
t107 = sin(qJ(4));
t110 = cos(qJ(4));
t108 = sin(qJ(3));
t111 = cos(qJ(3));
t103 = sin(pkin(10));
t137 = pkin(7) + qJ(2);
t85 = t137 * t103;
t105 = cos(pkin(10));
t87 = t137 * t105;
t63 = -t108 * t87 - t111 * t85;
t81 = t103 * t111 + t105 * t108;
t116 = -pkin(8) * t81 + t63;
t65 = -t108 * t85 + t111 * t87;
t79 = -t103 * t108 + t105 * t111;
t48 = pkin(8) * t79 + t65;
t30 = t107 * t48 - t110 * t116;
t144 = t30 ^ 2;
t143 = 0.2e1 * t30;
t102 = sin(pkin(11));
t104 = cos(pkin(11));
t106 = sin(qJ(6));
t109 = cos(qJ(6));
t78 = -t102 * t106 + t104 * t109;
t80 = t102 * t109 + t104 * t106;
t59 = -t78 * mrSges(7,1) + t80 * mrSges(7,2);
t142 = 0.2e1 * t59;
t93 = -pkin(2) * t105 - pkin(1);
t66 = -pkin(3) * t79 + t93;
t141 = 0.2e1 * t66;
t140 = 0.2e1 * t79;
t101 = t105 ^ 2;
t138 = pkin(3) * t110;
t57 = t107 * t81 - t110 * t79;
t58 = t107 * t79 + t110 * t81;
t29 = pkin(4) * t57 - qJ(5) * t58 + t66;
t32 = t107 * t116 + t110 * t48;
t12 = t102 * t29 + t104 * t32;
t130 = t104 * t58;
t132 = t102 * t58;
t38 = mrSges(6,1) * t132 + mrSges(6,2) * t130;
t136 = Ifges(7,5) * t80 + Ifges(7,6) * t78;
t135 = Ifges(6,4) * t102;
t134 = Ifges(6,4) * t104;
t11 = -t102 * t32 + t104 * t29;
t133 = t102 * t11;
t131 = t104 * t12;
t91 = pkin(3) * t107 + qJ(5);
t129 = t104 * t91;
t128 = t103 ^ 2 + t101;
t127 = t102 ^ 2 + t104 ^ 2;
t126 = 2 * mrSges(7,3);
t36 = t80 * t58;
t37 = t78 * t58;
t125 = Ifges(7,5) * t37 - Ifges(7,6) * t36 + Ifges(7,3) * t57;
t92 = -pkin(5) * t104 - pkin(4);
t124 = -t79 * mrSges(4,1) + t81 * mrSges(4,2);
t13 = t36 * mrSges(7,1) + t37 * mrSges(7,2);
t123 = -t105 * mrSges(3,1) + t103 * mrSges(3,2);
t84 = -t104 * mrSges(6,1) + t102 * mrSges(6,2);
t122 = t127 * qJ(5);
t60 = Ifges(7,4) * t80 + Ifges(7,2) * t78;
t61 = Ifges(7,1) * t80 + Ifges(7,4) * t78;
t88 = Ifges(6,2) * t104 + t135;
t89 = Ifges(6,1) * t102 + t134;
t121 = t102 * t89 + t104 * t88 + t78 * t60 + t80 * t61 + Ifges(5,3);
t120 = 0.2e1 * t127 * mrSges(6,3);
t119 = t131 - t133;
t118 = (mrSges(5,1) * t110 - mrSges(5,2) * t107) * pkin(3);
t117 = t59 + t84;
t15 = pkin(5) * t132 + t30;
t4 = pkin(5) * t57 - pkin(9) * t130 + t11;
t5 = -pkin(9) * t132 + t12;
t2 = -t106 * t5 + t109 * t4;
t22 = t57 * Ifges(6,6) + (-Ifges(6,2) * t102 + t134) * t58;
t23 = t57 * Ifges(6,5) + (Ifges(6,1) * t104 - t135) * t58;
t3 = t106 * t4 + t109 * t5;
t8 = Ifges(7,4) * t37 - Ifges(7,2) * t36 + Ifges(7,6) * t57;
t9 = Ifges(7,1) * t37 - Ifges(7,4) * t36 + Ifges(7,5) * t57;
t115 = -t32 * mrSges(5,2) + mrSges(6,3) * t131 + t15 * t59 + t102 * t23 / 0.2e1 + t104 * t22 / 0.2e1 - t36 * t60 / 0.2e1 + t37 * t61 / 0.2e1 - t88 * t132 / 0.2e1 + t89 * t130 / 0.2e1 - Ifges(5,6) * t57 + Ifges(5,5) * t58 + t78 * t8 / 0.2e1 + t80 * t9 / 0.2e1 + (t84 - mrSges(5,1)) * t30 + (Ifges(6,5) * t102 + Ifges(6,6) * t104 + t136) * t57 / 0.2e1 + (-t2 * t80 + t3 * t78) * mrSges(7,3);
t97 = t104 * pkin(9);
t94 = -pkin(4) - t138;
t86 = qJ(5) * t104 + t97;
t83 = (-pkin(9) - qJ(5)) * t102;
t82 = t92 - t138;
t69 = t97 + t129;
t68 = (-pkin(9) - t91) * t102;
t64 = t106 * t83 + t109 * t86;
t62 = -t106 * t86 + t109 * t83;
t53 = t58 * mrSges(5,2);
t52 = t106 * t68 + t109 * t69;
t51 = -t106 * t69 + t109 * t68;
t40 = mrSges(6,1) * t57 - mrSges(6,3) * t130;
t39 = -mrSges(6,2) * t57 - mrSges(6,3) * t132;
t17 = mrSges(7,1) * t57 - mrSges(7,3) * t37;
t16 = -mrSges(7,2) * t57 - mrSges(7,3) * t36;
t1 = [m(4) * (t63 ^ 2 + t65 ^ 2 + t93 ^ 2) + m(7) * (t15 ^ 2 + t2 ^ 2 + t3 ^ 2) + 0.2e1 * t11 * t40 - t36 * t8 + t37 * t9 + 0.2e1 * t12 * t39 + 0.2e1 * t15 * t13 + 0.2e1 * t3 * t16 + 0.2e1 * t2 * t17 + Ifges(3,2) * t101 + (-0.2e1 * t63 * mrSges(4,3) + Ifges(4,1) * t81 + Ifges(4,4) * t140) * t81 + (mrSges(5,1) * t141 - 0.2e1 * t32 * mrSges(5,3) + (Ifges(6,3) + Ifges(5,2)) * t57 + t125) * t57 + (mrSges(5,3) * t143 + Ifges(5,1) * t58 - t102 * t22 + t104 * t23 + (Ifges(6,5) * t104 - Ifges(6,6) * t102 - (2 * Ifges(5,4))) * t57) * t58 + m(5) * (t32 ^ 2 + t66 ^ 2 + t144) + m(6) * (t11 ^ 2 + t12 ^ 2 + t144) + m(3) * (qJ(2) ^ 2 * t128 + pkin(1) ^ 2) - 0.2e1 * pkin(1) * t123 + 0.2e1 * t93 * t124 + 0.2e1 * t128 * qJ(2) * mrSges(3,3) + Ifges(2,3) + t65 * mrSges(4,3) * t140 + t53 * t141 + t38 * t143 + Ifges(4,2) * t79 ^ 2 + (Ifges(3,1) * t103 + 0.2e1 * Ifges(3,4) * t105) * t103; -m(3) * pkin(1) + t57 * mrSges(5,1) + t102 * t39 + t104 * t40 + t80 * t16 + t78 * t17 + t53 + m(7) * (t2 * t78 + t3 * t80) + m(6) * (t102 * t12 + t104 * t11) + m(5) * t66 + m(4) * t93 + t123 + t124; m(3) + m(4) + m(5) + m(6) * t127 + m(7) * (t78 ^ 2 + t80 ^ 2); t82 * t13 + t94 * t38 + Ifges(4,6) * t79 + Ifges(4,5) * t81 + t63 * mrSges(4,1) - t65 * mrSges(4,2) + t51 * t17 + t52 * t16 + t39 * t129 + t115 + (m(5) * (t107 * t32 - t110 * t30) + (-t107 * t57 - t110 * t58) * mrSges(5,3)) * pkin(3) + m(6) * (t119 * t91 + t30 * t94) + (-t11 * mrSges(6,3) - t91 * t40) * t102 + m(7) * (t15 * t82 + t2 * t51 + t3 * t52); m(7) * (t51 * t78 + t52 * t80); t82 * t142 + 0.2e1 * t94 * t84 + Ifges(4,3) + 0.2e1 * t118 + (-t51 * t80 + t52 * t78) * t126 + t91 * t120 + m(7) * (t51 ^ 2 + t52 ^ 2 + t82 ^ 2) + m(6) * (t127 * t91 ^ 2 + t94 ^ 2) + m(5) * (t107 ^ 2 + t110 ^ 2) * pkin(3) ^ 2 + t121; t92 * t13 + t62 * t17 + t64 * t16 - pkin(4) * t38 - mrSges(6,3) * t133 + t115 + m(6) * (-pkin(4) * t30 + qJ(5) * t119) + m(7) * (t15 * t92 + t2 * t62 + t3 * t64) + (-t102 * t40 + t104 * t39) * qJ(5); m(7) * (t62 * t78 + t64 * t80); (t94 - pkin(4)) * t84 + (t82 + t92) * t59 + t118 + m(7) * (t51 * t62 + t52 * t64 + t82 * t92) + m(6) * (-pkin(4) * t94 + t122 * t91) + ((-t51 - t62) * t80 + (t52 + t64) * t78) * mrSges(7,3) + (t127 * t91 + t122) * mrSges(6,3) + t121; -0.2e1 * pkin(4) * t84 + t92 * t142 + (-t62 * t80 + t64 * t78) * t126 + qJ(5) * t120 + m(7) * (t62 ^ 2 + t64 ^ 2 + t92 ^ 2) + m(6) * (qJ(5) ^ 2 * t127 + pkin(4) ^ 2) + t121; m(6) * t30 + m(7) * t15 + t13 + t38; 0; m(6) * t94 + m(7) * t82 + t117; -m(6) * pkin(4) + m(7) * t92 + t117; m(6) + m(7); mrSges(7,1) * t2 - mrSges(7,2) * t3 + t125; -t59; mrSges(7,1) * t51 - mrSges(7,2) * t52 + t136; mrSges(7,1) * t62 - t64 * mrSges(7,2) + t136; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
