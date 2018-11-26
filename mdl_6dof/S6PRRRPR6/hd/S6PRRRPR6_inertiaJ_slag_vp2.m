% Calculate joint inertia matrix for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2018-11-23 15:25
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:25:19
% EndTime: 2018-11-23 15:25:20
% DurationCPUTime: 1.10s
% Computational Cost: add. (1032->297), mult. (2193->402), div. (0->0), fcn. (2089->10), ass. (0->117)
t91 = sin(qJ(4));
t95 = cos(qJ(4));
t145 = t91 ^ 2 + t95 ^ 2;
t88 = sin(pkin(6));
t97 = cos(qJ(2));
t123 = t88 * t97;
t93 = sin(qJ(2));
t124 = t88 * t93;
t89 = cos(pkin(6));
t92 = sin(qJ(3));
t96 = cos(qJ(3));
t40 = t124 * t96 + t89 * t92;
t15 = t123 * t95 + t40 * t91;
t17 = -t123 * t91 + t40 * t95;
t144 = t15 * t91 + t17 * t95;
t143 = 2 * pkin(8);
t142 = mrSges(5,3) + mrSges(6,2);
t121 = t92 * t95;
t122 = t91 * t92;
t141 = -Ifges(6,6) * t122 + (-Ifges(6,4) - Ifges(5,5)) * t121;
t38 = t124 * t92 - t89 * t96;
t140 = t38 ^ 2;
t139 = pkin(4) + pkin(5);
t138 = pkin(9) - pkin(10);
t137 = pkin(8) * t96;
t136 = Ifges(5,4) * t91;
t135 = Ifges(5,4) * t95;
t134 = Ifges(6,5) * t91;
t133 = Ifges(6,5) * t95;
t132 = Ifges(5,6) * t96;
t131 = Ifges(6,6) * t96;
t128 = t38 * t92;
t127 = t40 * t96;
t90 = sin(qJ(6));
t94 = cos(qJ(6));
t53 = -qJ(5) * t90 - t139 * t94;
t126 = t53 * mrSges(7,1);
t54 = qJ(5) * t94 - t139 * t90;
t125 = t54 * mrSges(7,2);
t120 = Ifges(6,2) + Ifges(5,3);
t49 = mrSges(5,2) * t96 - mrSges(5,3) * t122;
t52 = -mrSges(6,2) * t122 - mrSges(6,3) * t96;
t119 = t49 + t52;
t50 = -mrSges(5,1) * t96 - mrSges(5,3) * t121;
t51 = t96 * mrSges(6,1) + mrSges(6,2) * t121;
t118 = -t50 + t51;
t58 = -mrSges(5,1) * t95 + mrSges(5,2) * t91;
t117 = t58 - mrSges(4,1);
t56 = -pkin(3) * t96 - pkin(9) * t92 - pkin(2);
t26 = t95 * t137 + t91 * t56;
t116 = t145 * pkin(9) ^ 2;
t115 = qJ(5) * t95;
t48 = -t90 * t95 + t91 * t94;
t34 = t48 * t92;
t103 = t90 * t91 + t94 * t95;
t35 = t103 * t92;
t114 = Ifges(7,5) * t35 + Ifges(7,6) * t34 + Ifges(7,3) * t96;
t113 = qJ(5) * t91 + pkin(3);
t70 = t91 * t137;
t25 = t56 * t95 - t70;
t111 = t144 * pkin(9);
t21 = -qJ(5) * t96 + t26;
t3 = t15 * t94 - t17 * t90;
t4 = t15 * t90 + t17 * t94;
t109 = t3 * mrSges(7,1) - t4 * mrSges(7,2);
t108 = t91 * mrSges(5,1) + t95 * mrSges(5,2);
t107 = t91 * mrSges(6,1) - t95 * mrSges(6,3);
t106 = t94 * mrSges(7,1) - t90 * mrSges(7,2);
t105 = -pkin(4) * t91 + t115;
t64 = t138 * t91;
t65 = t138 * t95;
t19 = t64 * t94 - t65 * t90;
t20 = t64 * t90 + t65 * t94;
t44 = Ifges(7,6) * t103;
t45 = Ifges(7,5) * t48;
t102 = t19 * mrSges(7,1) - t20 * mrSges(7,2) - t44 + t45;
t10 = pkin(10) * t122 + t21;
t82 = t96 * pkin(4);
t8 = pkin(5) * t96 + t70 + t82 + (-pkin(10) * t92 - t56) * t95;
t1 = -t10 * t90 + t8 * t94;
t2 = t10 * t94 + t8 * t90;
t101 = t1 * mrSges(7,1) - t2 * mrSges(7,2) + t114;
t100 = pkin(8) ^ 2;
t87 = t96 ^ 2;
t85 = t92 ^ 2;
t83 = t88 ^ 2;
t80 = t85 * t100;
t78 = Ifges(6,4) * t91;
t77 = Ifges(5,5) * t91;
t76 = Ifges(5,6) * t95;
t72 = t83 * t97 ^ 2;
t63 = Ifges(5,1) * t91 + t135;
t62 = Ifges(6,1) * t91 - t133;
t61 = Ifges(5,2) * t95 + t136;
t60 = -Ifges(6,3) * t95 + t134;
t59 = -mrSges(4,1) * t96 + mrSges(4,2) * t92;
t57 = -mrSges(6,1) * t95 - mrSges(6,3) * t91;
t55 = -pkin(4) * t95 - t113;
t43 = t139 * t95 + t113;
t42 = t108 * t92;
t41 = t107 * t92;
t33 = (pkin(8) - t105) * t92;
t32 = -Ifges(5,5) * t96 + (Ifges(5,1) * t95 - t136) * t92;
t31 = -Ifges(6,4) * t96 + (Ifges(6,1) * t95 + t134) * t92;
t30 = -t132 + (-Ifges(5,2) * t91 + t135) * t92;
t29 = -t131 + (Ifges(6,3) * t91 + t133) * t92;
t24 = mrSges(7,1) * t96 - mrSges(7,3) * t35;
t23 = -mrSges(7,2) * t96 + mrSges(7,3) * t34;
t22 = -t25 + t82;
t18 = (-t139 * t91 - pkin(8) + t115) * t92;
t14 = Ifges(7,1) * t48 - Ifges(7,4) * t103;
t13 = Ifges(7,4) * t48 - Ifges(7,2) * t103;
t12 = mrSges(7,1) * t103 + mrSges(7,2) * t48;
t7 = -mrSges(7,1) * t34 + mrSges(7,2) * t35;
t6 = Ifges(7,1) * t35 + Ifges(7,4) * t34 + Ifges(7,5) * t96;
t5 = Ifges(7,4) * t35 + Ifges(7,2) * t34 + Ifges(7,6) * t96;
t9 = [m(2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t140) + m(4) * (t40 ^ 2 + t140 + t72) + m(3) * (t83 * t93 ^ 2 + t89 ^ 2 + t72) + (m(6) + m(5)) * (t15 ^ 2 + t17 ^ 2 + t140); mrSges(4,3) * t127 + t4 * t23 + t3 * t24 + t119 * t17 + t118 * t15 + (-t93 * mrSges(3,2) + (mrSges(3,1) - t59) * t97) * t88 + (t92 * mrSges(4,3) + t41 + t42 - t7) * t38 + m(7) * (t1 * t3 - t18 * t38 + t2 * t4) + m(6) * (t15 * t22 + t17 * t21 + t33 * t38) + m(5) * (pkin(8) * t128 - t15 * t25 + t17 * t26) + m(4) * (pkin(2) * t123 + (t127 + t128) * pkin(8)); -0.2e1 * pkin(2) * t59 + 0.2e1 * t1 * t24 + 0.2e1 * t18 * t7 + 0.2e1 * t2 * t23 + 0.2e1 * t21 * t52 + 0.2e1 * t22 * t51 + 0.2e1 * t25 * t50 + 0.2e1 * t26 * t49 + 0.2e1 * t33 * t41 + t34 * t5 + t35 * t6 + Ifges(3,3) + (t85 + t87) * mrSges(4,3) * t143 + m(4) * (pkin(2) ^ 2 + t100 * t87 + t80) + m(5) * (t25 ^ 2 + t26 ^ 2 + t80) + m(6) * (t21 ^ 2 + t22 ^ 2 + t33 ^ 2) + m(7) * (t1 ^ 2 + t18 ^ 2 + t2 ^ 2) + ((Ifges(4,2) + t120) * t96 + t114 + t141) * t96 + (Ifges(4,1) * t92 + 0.2e1 * Ifges(4,4) * t96 + t42 * t143 + (t31 + t32) * t95 + (t29 - t30 + t132) * t91) * t92; -t40 * mrSges(4,2) + (-t103 * t4 - t3 * t48) * mrSges(7,3) + (-t12 + t57 + t117) * t38 + m(7) * (t19 * t3 + t20 * t4 - t38 * t43) + m(6) * (t38 * t55 + t111) + m(5) * (-pkin(3) * t38 + t111) + t142 * t144; t55 * t41 - t103 * t5 / 0.2e1 + t48 * t6 / 0.2e1 + t34 * t13 / 0.2e1 + t35 * t14 / 0.2e1 - pkin(3) * t42 + t43 * t7 + t20 * t23 + t19 * t24 + t18 * t12 + (-t1 * t48 - t103 * t2) * mrSges(7,3) + m(7) * (t1 * t19 + t18 * t43 + t2 * t20) + (-t78 / 0.2e1 - t77 / 0.2e1 - t76 / 0.2e1 + t45 / 0.2e1 - t44 / 0.2e1 + Ifges(4,6) - pkin(8) * mrSges(4,2)) * t96 + (-t29 / 0.2e1 + t30 / 0.2e1 + t131 / 0.2e1 + t21 * mrSges(6,2) + t26 * mrSges(5,3)) * t95 + (t31 / 0.2e1 + t32 / 0.2e1 + t22 * mrSges(6,2) - t25 * mrSges(5,3)) * t91 + (t119 * t95 + t118 * t91 + m(5) * (-t25 * t91 + t26 * t95) + m(6) * (t21 * t95 + t22 * t91)) * pkin(9) + (Ifges(4,5) + (t62 / 0.2e1 + t63 / 0.2e1) * t95 + (t60 / 0.2e1 - t61 / 0.2e1) * t91 + (-m(5) * pkin(3) + t117) * pkin(8)) * t92 + (m(6) * t55 + t57) * t33; -0.2e1 * pkin(3) * t58 + 0.2e1 * t43 * t12 - t103 * t13 + t48 * t14 + 0.2e1 * t55 * t57 + Ifges(4,3) + (-t60 + t61) * t95 + (t63 + t62) * t91 + 0.2e1 * (-t103 * t20 - t19 * t48) * mrSges(7,3) + m(7) * (t19 ^ 2 + t20 ^ 2 + t43 ^ 2) + m(6) * (t55 ^ 2 + t116) + m(5) * (pkin(3) ^ 2 + t116) + 0.2e1 * t142 * pkin(9) * t145; (-mrSges(5,2) + mrSges(6,3)) * t17 + (-mrSges(5,1) - mrSges(6,1)) * t15 + m(7) * (t3 * t53 + t4 * t54) + m(6) * (-pkin(4) * t15 + qJ(5) * t17) - t109; -t101 + m(7) * (t1 * t53 + t2 * t54) + m(6) * (-pkin(4) * t22 + qJ(5) * t21) - t120 * t96 + t54 * t23 - pkin(4) * t51 + qJ(5) * t52 + t53 * t24 - Ifges(5,6) * t122 + t21 * mrSges(6,3) - t22 * mrSges(6,1) + t25 * mrSges(5,1) - t26 * mrSges(5,2) - t141; m(7) * (t19 * t53 + t20 * t54) - Ifges(6,6) * t95 + t77 + t76 + t78 + (-t103 * t54 - t48 * t53) * mrSges(7,3) + t105 * mrSges(6,2) + (m(6) * t105 - t107 - t108) * pkin(9) - t102; 0.2e1 * pkin(4) * mrSges(6,1) - 0.2e1 * t126 + 0.2e1 * t125 + 0.2e1 * qJ(5) * mrSges(6,3) + Ifges(7,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + m(7) * (t53 ^ 2 + t54 ^ 2) + t120; m(7) * (t3 * t94 + t4 * t90) + m(6) * t15; t90 * t23 + t94 * t24 + m(7) * (t1 * t94 + t2 * t90) + m(6) * t22 + t51; m(7) * (t19 * t94 + t20 * t90) + (m(6) * pkin(9) + mrSges(6,2)) * t91 + (-t103 * t90 - t48 * t94) * mrSges(7,3); -mrSges(6,1) - m(6) * pkin(4) + m(7) * (t53 * t94 + t54 * t90) - t106; m(6) + m(7) * (t90 ^ 2 + t94 ^ 2); t109; t101; t102; -Ifges(7,3) - t125 + t126; t106; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t9(1) t9(2) t9(4) t9(7) t9(11) t9(16); t9(2) t9(3) t9(5) t9(8) t9(12) t9(17); t9(4) t9(5) t9(6) t9(9) t9(13) t9(18); t9(7) t9(8) t9(9) t9(10) t9(14) t9(19); t9(11) t9(12) t9(13) t9(14) t9(15) t9(20); t9(16) t9(17) t9(18) t9(19) t9(20) t9(21);];
Mq  = res;
