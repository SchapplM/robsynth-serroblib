% Calculate joint inertia matrix for
% S6PRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2018-11-23 15:12
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:12:20
% EndTime: 2018-11-23 15:12:21
% DurationCPUTime: 0.87s
% Computational Cost: add. (1083->256), mult. (2403->353), div. (0->0), fcn. (2423->10), ass. (0->106)
t140 = 2 * pkin(8);
t139 = m(6) + m(7);
t138 = mrSges(7,2) + mrSges(6,3);
t130 = cos(qJ(5));
t94 = sin(pkin(11));
t96 = cos(pkin(11));
t98 = sin(qJ(5));
t105 = t130 * t96 - t98 * t94;
t67 = t130 * t94 + t98 * t96;
t22 = -mrSges(7,1) * t105 - t67 * mrSges(7,3);
t23 = -mrSges(6,1) * t105 + t67 * mrSges(6,2);
t137 = t22 + t23;
t99 = sin(qJ(3));
t50 = t67 * t99;
t51 = t105 * t99;
t15 = t50 * mrSges(7,1) - t51 * mrSges(7,3);
t16 = t50 * mrSges(6,1) + t51 * mrSges(6,2);
t124 = t96 * t99;
t125 = t94 * t99;
t56 = mrSges(5,1) * t125 + mrSges(5,2) * t124;
t136 = t15 + t16 + t56;
t135 = -m(7) * pkin(5) - mrSges(7,1);
t134 = -mrSges(6,1) + t135;
t133 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t101 = cos(qJ(3));
t100 = sin(qJ(2));
t95 = sin(pkin(6));
t116 = t100 * t95;
t97 = cos(pkin(6));
t53 = -t97 * t101 + t116 * t99;
t52 = t53 ^ 2;
t132 = -t94 / 0.2e1;
t131 = t96 / 0.2e1;
t88 = t99 * pkin(8);
t129 = Ifges(5,4) * t94;
t128 = Ifges(5,4) * t96;
t127 = pkin(8) * t101;
t126 = t53 * t99;
t122 = -Ifges(6,3) - Ifges(7,2);
t121 = pkin(9) + qJ(4);
t71 = -pkin(3) * t101 - qJ(4) * t99 - pkin(2);
t64 = t96 * t71;
t19 = -pkin(9) * t124 + t64 + (-pkin(8) * t94 - pkin(4)) * t101;
t40 = t96 * t127 + t94 * t71;
t30 = -pkin(9) * t125 + t40;
t5 = t130 * t30 + t98 * t19;
t35 = -mrSges(7,2) * t50 - mrSges(7,3) * t101;
t36 = mrSges(6,2) * t101 - mrSges(6,3) * t50;
t120 = t35 + t36;
t37 = -mrSges(6,1) * t101 - mrSges(6,3) * t51;
t38 = t101 * mrSges(7,1) + t51 * mrSges(7,2);
t119 = -t37 + t38;
t72 = -t96 * mrSges(5,1) + t94 * mrSges(5,2);
t118 = t72 - mrSges(4,1);
t70 = pkin(4) * t125 + t88;
t117 = t94 ^ 2 + t96 ^ 2;
t55 = t101 * t116 + t97 * t99;
t115 = t101 * t55;
t102 = cos(qJ(2));
t114 = t102 * t95;
t110 = t121 * t94;
t73 = t121 * t96;
t32 = t110 * t130 + t73 * t98;
t34 = -t110 * t98 + t130 * t73;
t113 = t32 ^ 2 + t34 ^ 2;
t84 = -pkin(4) * t96 - pkin(3);
t28 = -t114 * t96 - t55 * t94;
t29 = -t114 * t94 + t55 * t96;
t7 = -t130 * t28 + t29 * t98;
t9 = t130 * t29 + t98 * t28;
t112 = t32 * t7 + t34 * t9;
t109 = (-Ifges(7,4) - Ifges(6,5)) * t51 + (Ifges(6,6) - Ifges(7,6)) * t50;
t107 = -t28 * t94 + t29 * t96;
t39 = -t127 * t94 + t64;
t106 = -t39 * t94 + t40 * t96;
t4 = t130 * t19 - t98 * t30;
t104 = pkin(8) ^ 2;
t93 = t101 ^ 2;
t92 = t99 ^ 2;
t90 = t95 ^ 2;
t87 = t92 * t104;
t81 = t90 * t102 ^ 2;
t76 = -mrSges(4,1) * t101 + mrSges(4,2) * t99;
t75 = Ifges(5,1) * t94 + t128;
t74 = Ifges(5,2) * t96 + t129;
t69 = -mrSges(5,1) * t101 - mrSges(5,3) * t124;
t68 = mrSges(5,2) * t101 - mrSges(5,3) * t125;
t62 = Ifges(7,4) * t67;
t61 = Ifges(6,5) * t67;
t60 = Ifges(6,6) * t105;
t59 = Ifges(7,6) * t105;
t49 = -Ifges(5,5) * t101 + (Ifges(5,1) * t96 - t129) * t99;
t48 = -Ifges(5,6) * t101 + (-Ifges(5,2) * t94 + t128) * t99;
t27 = Ifges(6,1) * t67 + Ifges(6,4) * t105;
t26 = Ifges(7,1) * t67 - Ifges(7,5) * t105;
t25 = Ifges(6,4) * t67 + Ifges(6,2) * t105;
t24 = Ifges(7,5) * t67 - Ifges(7,3) * t105;
t18 = -pkin(5) * t105 - qJ(6) * t67 + t84;
t14 = Ifges(6,1) * t51 - Ifges(6,4) * t50 - Ifges(6,5) * t101;
t13 = Ifges(7,1) * t51 - Ifges(7,4) * t101 + Ifges(7,5) * t50;
t12 = Ifges(6,4) * t51 - Ifges(6,2) * t50 - Ifges(6,6) * t101;
t11 = Ifges(7,5) * t51 - Ifges(7,6) * t101 + Ifges(7,3) * t50;
t10 = pkin(5) * t50 - qJ(6) * t51 + t70;
t3 = t101 * pkin(5) - t4;
t2 = -qJ(6) * t101 + t5;
t1 = [m(2) + m(5) * (t28 ^ 2 + t29 ^ 2 + t52) + m(4) * (t55 ^ 2 + t52 + t81) + m(3) * (t100 ^ 2 * t90 + t97 ^ 2 + t81) + t139 * (t7 ^ 2 + t9 ^ 2 + t52); mrSges(4,3) * t115 + t28 * t69 + t29 * t68 + t120 * t9 + t119 * t7 + (-t100 * mrSges(3,2) + (mrSges(3,1) - t76) * t102) * t95 + (t99 * mrSges(4,3) + t136) * t53 + m(6) * (-t4 * t7 + t5 * t9 + t53 * t70) + m(7) * (t10 * t53 + t2 * t9 + t3 * t7) + m(5) * (pkin(8) * t126 + t28 * t39 + t29 * t40) + m(4) * (pkin(2) * t114 + (t115 + t126) * pkin(8)); -0.2e1 * pkin(2) * t76 + 0.2e1 * t10 * t15 + 0.2e1 * t70 * t16 + 0.2e1 * t2 * t35 + 0.2e1 * t3 * t38 + 0.2e1 * t5 * t36 + 0.2e1 * t4 * t37 + 0.2e1 * t39 * t69 + 0.2e1 * t40 * t68 + Ifges(3,3) + (t13 + t14) * t51 + (t11 - t12) * t50 + (t92 + t93) * mrSges(4,3) * t140 + (Ifges(4,1) * t99 + t56 * t140 - t94 * t48 + t96 * t49) * t99 + m(4) * (pkin(2) ^ 2 + t104 * t93 + t87) + m(5) * (t39 ^ 2 + t40 ^ 2 + t87) + m(6) * (t4 ^ 2 + t5 ^ 2 + t70 ^ 2) + m(7) * (t10 ^ 2 + t2 ^ 2 + t3 ^ 2) + ((Ifges(5,3) + Ifges(4,2) - t122) * t101 + (-Ifges(5,5) * t96 + Ifges(5,6) * t94 + (2 * Ifges(4,4))) * t99 + t109) * t101; -t55 * mrSges(4,2) + t107 * mrSges(5,3) + (t118 + t137) * t53 + m(6) * (t53 * t84 + t112) + m(7) * (t18 * t53 + t112) + m(5) * (-pkin(3) * t53 + qJ(4) * t107) + t138 * (t105 * t9 + t67 * t7); t94 * t49 / 0.2e1 + t48 * t131 + t84 * t16 + t70 * t23 - pkin(3) * t56 + t18 * t15 + t10 * t22 + (t26 / 0.2e1 + t27 / 0.2e1) * t51 + (t24 / 0.2e1 - t25 / 0.2e1) * t50 + t120 * t34 + t119 * t32 + (-pkin(8) * mrSges(4,2) + Ifges(5,5) * t132 - Ifges(5,6) * t96 / 0.2e1 - t61 / 0.2e1 - t60 / 0.2e1 - t62 / 0.2e1 + t59 / 0.2e1 + Ifges(4,6)) * t101 + (t96 * t68 - t94 * t69) * qJ(4) + t106 * mrSges(5,3) + (pkin(8) * t118 + t131 * t75 + t132 * t74 + Ifges(4,5)) * t99 + m(5) * (-pkin(3) * t88 + qJ(4) * t106) + m(6) * (-t32 * t4 + t34 * t5 + t70 * t84) + m(7) * (t10 * t18 + t2 * t34 + t3 * t32) + (t3 * mrSges(7,2) - t4 * mrSges(6,3) + t13 / 0.2e1 + t14 / 0.2e1) * t67 - (-t2 * mrSges(7,2) - t5 * mrSges(6,3) + t11 / 0.2e1 - t12 / 0.2e1) * t105; -0.2e1 * pkin(3) * t72 + 0.2e1 * t18 * t22 + 0.2e1 * t84 * t23 + t96 * t74 + t94 * t75 + Ifges(4,3) + m(7) * (t18 ^ 2 + t113) + m(6) * (t84 ^ 2 + t113) + m(5) * (qJ(4) ^ 2 * t117 + pkin(3) ^ 2) + (t26 + t27) * t67 - (t24 - t25) * t105 + 0.2e1 * t117 * qJ(4) * mrSges(5,3) + 0.2e1 * (t105 * t34 + t32 * t67) * t138; 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1 + m(5) / 0.2e1) * t53; m(5) * t88 + m(6) * t70 + m(7) * t10 + t136; -m(5) * pkin(3) + m(6) * t84 + m(7) * t18 + t137 + t72; m(5) + t139; t133 * t9 + t134 * t7; -t5 * mrSges(6,2) + t4 * mrSges(6,1) - pkin(5) * t38 + t2 * mrSges(7,3) + qJ(6) * t35 + m(7) * (-pkin(5) * t3 + qJ(6) * t2) - t3 * mrSges(7,1) + t122 * t101 - t109; -t59 + t61 + t60 + t62 + (-pkin(5) * t67 + qJ(6) * t105) * mrSges(7,2) + t133 * t34 + t134 * t32; 0; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) - t122; m(7) * t7; m(7) * t3 + t38; m(7) * t32 + t67 * mrSges(7,2); 0; t135; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
