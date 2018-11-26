% Calculate joint inertia matrix for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:24
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:24:33
% EndTime: 2018-11-23 16:24:34
% DurationCPUTime: 0.99s
% Computational Cost: add. (1288->272), mult. (2522->366), div. (0->0), fcn. (2391->8), ass. (0->107)
t92 = sin(qJ(5));
t93 = sin(qJ(4));
t95 = cos(qJ(5));
t96 = cos(qJ(4));
t62 = t92 * t96 + t93 * t95;
t94 = sin(qJ(3));
t49 = t62 * t94;
t61 = -t92 * t93 + t95 * t96;
t50 = t61 * t94;
t19 = mrSges(7,1) * t49 + mrSges(7,2) * t50;
t20 = mrSges(6,1) * t49 + t50 * mrSges(6,2);
t103 = -t20 - t19;
t106 = mrSges(5,1) * t93 + mrSges(5,2) * t96;
t54 = t106 * t94;
t135 = -t54 + t103;
t90 = sin(pkin(10));
t79 = pkin(1) * t90 + pkin(7);
t134 = 0.2e1 * t79;
t133 = (t95 * mrSges(6,1) + (-mrSges(6,2) - mrSges(7,2)) * t92) * pkin(4);
t132 = 0.2e1 * mrSges(7,1);
t131 = m(7) * pkin(5);
t130 = -pkin(9) - pkin(8);
t129 = m(7) * t92;
t97 = cos(qJ(3));
t128 = pkin(3) * t97;
t127 = pkin(4) * t95;
t126 = pkin(5) * t49;
t125 = pkin(8) * t94;
t119 = t94 * t96;
t91 = cos(pkin(10));
t80 = -pkin(1) * t91 - pkin(2);
t55 = -t125 + t80 - t128;
t52 = t96 * t55;
t14 = -pkin(9) * t119 + t52 + (-t79 * t93 - pkin(4)) * t97;
t120 = t93 * t94;
t121 = t79 * t97;
t24 = t121 * t96 + t55 * t93;
t21 = -pkin(9) * t120 + t24;
t6 = t14 * t92 + t21 * t95;
t124 = Ifges(5,4) * t93;
t123 = Ifges(5,4) * t96;
t122 = t61 * t92;
t73 = t94 * t79;
t117 = Ifges(6,3) + Ifges(7,3);
t34 = mrSges(7,2) * t97 - mrSges(7,3) * t49;
t35 = mrSges(6,2) * t97 - mrSges(6,3) * t49;
t116 = t34 + t35;
t74 = t130 * t93;
t75 = t130 * t96;
t33 = t74 * t92 - t75 * t95;
t69 = -mrSges(5,1) * t96 + mrSges(5,2) * t93;
t115 = t69 - mrSges(4,1);
t53 = pkin(4) * t120 + t73;
t114 = t93 ^ 2 + t96 ^ 2;
t87 = t94 ^ 2;
t89 = t97 ^ 2;
t113 = t87 + t89;
t112 = Ifges(5,3) + t117;
t5 = t14 * t95 - t21 * t92;
t2 = -pkin(5) * t97 - qJ(6) * t50 + t5;
t36 = -mrSges(7,1) * t97 - mrSges(7,3) * t50;
t111 = m(7) * t2 + t36;
t82 = -pkin(4) * t96 - pkin(3);
t109 = t114 * mrSges(5,3);
t25 = -mrSges(7,1) * t61 + mrSges(7,2) * t62;
t32 = t74 * t95 + t75 * t92;
t108 = (-Ifges(6,5) - Ifges(7,5)) * t50 + (Ifges(6,6) + Ifges(7,6)) * t49;
t17 = -qJ(6) * t62 + t32;
t107 = m(7) * t17 - t62 * mrSges(7,3);
t23 = -t121 * t93 + t52;
t105 = -t23 * t93 + t24 * t96;
t3 = -qJ(6) * t49 + t6;
t102 = mrSges(6,1) * t5 + mrSges(7,1) * t2 - t6 * mrSges(6,2) - t3 * mrSges(7,2) - t108;
t18 = qJ(6) * t61 + t33;
t57 = Ifges(7,6) * t61;
t58 = Ifges(6,6) * t61;
t59 = Ifges(7,5) * t62;
t60 = Ifges(6,5) * t62;
t101 = mrSges(6,1) * t32 + mrSges(7,1) * t17 - t33 * mrSges(6,2) - t18 * mrSges(7,2) + t57 + t58 + t59 + t60;
t99 = pkin(4) ^ 2;
t85 = t92 ^ 2 * t99;
t84 = Ifges(5,5) * t93;
t83 = Ifges(5,6) * t96;
t81 = pkin(5) + t127;
t78 = t79 ^ 2;
t76 = Ifges(5,5) * t119;
t72 = t87 * t78;
t71 = Ifges(5,1) * t93 + t123;
t70 = Ifges(5,2) * t96 + t124;
t67 = -mrSges(5,1) * t97 - mrSges(5,3) * t119;
t66 = mrSges(5,2) * t97 - mrSges(5,3) * t120;
t48 = -Ifges(5,5) * t97 + (Ifges(5,1) * t96 - t124) * t94;
t47 = -Ifges(5,6) * t97 + (-Ifges(5,2) * t93 + t123) * t94;
t39 = -pkin(5) * t61 + t82;
t38 = t92 * pkin(4) * t50;
t37 = -mrSges(6,1) * t97 - mrSges(6,3) * t50;
t30 = Ifges(6,1) * t62 + Ifges(6,4) * t61;
t29 = Ifges(7,1) * t62 + Ifges(7,4) * t61;
t28 = Ifges(6,4) * t62 + Ifges(6,2) * t61;
t27 = Ifges(7,4) * t62 + Ifges(7,2) * t61;
t26 = -mrSges(6,1) * t61 + mrSges(6,2) * t62;
t22 = t53 + t126;
t12 = Ifges(6,1) * t50 - Ifges(6,4) * t49 - Ifges(6,5) * t97;
t11 = Ifges(7,1) * t50 - Ifges(7,4) * t49 - Ifges(7,5) * t97;
t10 = Ifges(6,4) * t50 - Ifges(6,2) * t49 - Ifges(6,6) * t97;
t9 = Ifges(7,4) * t50 - Ifges(7,2) * t49 - Ifges(7,6) * t97;
t1 = [0.2e1 * t22 * t19 + 0.2e1 * t2 * t36 + 0.2e1 * t53 * t20 + 0.2e1 * t23 * t67 + 0.2e1 * t24 * t66 + 0.2e1 * t3 * t34 + 0.2e1 * t6 * t35 + 0.2e1 * t5 * t37 + Ifges(2,3) + Ifges(3,3) + (t11 + t12) * t50 - (t9 + t10) * t49 + (0.2e1 * mrSges(4,2) * t80 + Ifges(4,1) * t94 + t54 * t134 - t47 * t93 + t48 * t96) * t94 + (-0.2e1 * t80 * mrSges(4,1) - t76 + (Ifges(4,2) + t112) * t97 + (Ifges(5,6) * t93 + (2 * Ifges(4,4))) * t94 + t108) * t97 + m(4) * (t78 * t89 + t80 ^ 2 + t72) + m(5) * (t23 ^ 2 + t24 ^ 2 + t72) + m(6) * (t5 ^ 2 + t53 ^ 2 + t6 ^ 2) + m(7) * (t2 ^ 2 + t22 ^ 2 + t3 ^ 2) + m(3) * (t90 ^ 2 + t91 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (mrSges(3,1) * t91 - mrSges(3,2) * t90) * pkin(1) + t113 * mrSges(4,3) * t134; t116 * t50 - (t36 + t37) * t49 + t135 * t97 + m(6) * (-t49 * t5 + t50 * t6 - t53 * t97) + m(7) * (-t2 * t49 - t22 * t97 + t3 * t50) + (-t93 * t67 + t96 * t66 + m(5) * (t105 - t121)) * t94; m(3) + m(5) * (t114 * t87 + t89) + m(4) * t113 + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t49 ^ 2 + t50 ^ 2 + t89); -pkin(3) * t54 + t17 * t36 + t18 * t34 + t39 * t19 + t82 * t20 + t22 * t25 + t53 * t26 + t32 * t37 + t33 * t35 + (-t84 / 0.2e1 - t83 / 0.2e1 - t59 / 0.2e1 - t57 / 0.2e1 - t60 / 0.2e1 - t58 / 0.2e1 + Ifges(4,6) - t79 * mrSges(4,2)) * t97 + (t29 / 0.2e1 + t30 / 0.2e1) * t50 - (t27 / 0.2e1 + t28 / 0.2e1) * t49 + (t47 / 0.2e1 + pkin(8) * t66 + t24 * mrSges(5,3)) * t96 + (t48 / 0.2e1 - pkin(8) * t67 - t23 * mrSges(5,3)) * t93 + (Ifges(4,5) - t93 * t70 / 0.2e1 + t96 * t71 / 0.2e1 + t115 * t79) * t94 + m(5) * (-pkin(3) * t73 + pkin(8) * t105) + m(6) * (t32 * t5 + t33 * t6 + t53 * t82) + m(7) * (t17 * t2 + t18 * t3 + t22 * t39) + (t11 / 0.2e1 + t12 / 0.2e1 - t5 * mrSges(6,3) - t2 * mrSges(7,3)) * t62 + (t9 / 0.2e1 + t10 / 0.2e1 + t3 * mrSges(7,3) + t6 * mrSges(6,3)) * t61; (-mrSges(4,2) + t109) * t94 + (-t25 - t26 - t115) * t97 + m(6) * (-t32 * t49 + t33 * t50 - t82 * t97) + m(7) * (-t17 * t49 + t18 * t50 - t39 * t97) + m(5) * (t114 * t125 + t128) + (mrSges(7,3) + mrSges(6,3)) * (t49 * t62 + t50 * t61); -0.2e1 * pkin(3) * t69 + 0.2e1 * t39 * t25 + 0.2e1 * t82 * t26 + t96 * t70 + t93 * t71 + Ifges(4,3) + 0.2e1 * pkin(8) * t109 + m(6) * (t32 ^ 2 + t33 ^ 2 + t82 ^ 2) + m(7) * (t17 ^ 2 + t18 ^ 2 + t39 ^ 2) + m(5) * (pkin(8) ^ 2 * t114 + pkin(3) ^ 2) + (-0.2e1 * mrSges(6,3) * t32 - 0.2e1 * mrSges(7,3) * t17 + t29 + t30) * t62 + (0.2e1 * mrSges(6,3) * t33 + 0.2e1 * mrSges(7,3) * t18 + t27 + t28) * t61; -Ifges(5,6) * t120 + t23 * mrSges(5,1) - t24 * mrSges(5,2) + t76 + t111 * t81 + (t95 * t37 + t116 * t92 + t3 * t129 + m(6) * (t5 * t95 + t6 * t92)) * pkin(4) - t112 * t97 + t102; m(6) * (-t127 * t49 + t38) + m(7) * (-t49 * t81 + t38) + t135; t83 + t84 + t107 * t81 - t106 * pkin(8) + (mrSges(7,3) * t122 + (-t62 * t95 + t122) * mrSges(6,3) + m(6) * (t32 * t95 + t33 * t92) + t18 * t129) * pkin(4) + t101; t81 * t132 + m(6) * (t95 ^ 2 * t99 + t85) + m(7) * (t81 ^ 2 + t85) + 0.2e1 * t133 + t112; pkin(5) * t111 - t117 * t97 + t102; -m(7) * t126 + t103; pkin(5) * t107 + t101; t81 * t131 + (pkin(5) + t81) * mrSges(7,1) + t133 + t117; (t132 + t131) * pkin(5) + t117; m(7) * t22 + t19; -m(7) * t97; m(7) * t39 + t25; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
