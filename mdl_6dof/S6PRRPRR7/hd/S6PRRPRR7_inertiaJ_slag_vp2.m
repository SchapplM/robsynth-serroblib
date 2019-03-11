% Calculate joint inertia matrix for
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:30:17
% EndTime: 2019-03-08 22:30:19
% DurationCPUTime: 0.86s
% Computational Cost: add. (979->247), mult. (1975->340), div. (0->0), fcn. (1900->10), ass. (0->102)
t132 = pkin(4) + pkin(8);
t82 = sin(qJ(3));
t86 = cos(qJ(3));
t131 = t82 ^ 2 + t86 ^ 2;
t118 = mrSges(5,1) + mrSges(4,3);
t129 = -m(5) * pkin(3) + mrSges(5,2);
t78 = sin(pkin(6));
t83 = sin(qJ(2));
t122 = t78 * t83;
t79 = cos(pkin(6));
t35 = t86 * t122 + t79 * t82;
t31 = t35 ^ 2;
t128 = m(7) * pkin(5);
t81 = sin(qJ(5));
t127 = -t81 / 0.2e1;
t88 = -pkin(3) - pkin(9);
t126 = -pkin(10) + t88;
t125 = Ifges(6,4) * t81;
t85 = cos(qJ(5));
t124 = Ifges(6,4) * t85;
t33 = t82 * t122 - t79 * t86;
t123 = t33 * t82;
t87 = cos(qJ(2));
t121 = t78 * t87;
t120 = t85 * mrSges(6,1);
t119 = t85 * t86;
t106 = -qJ(4) * t82 - pkin(2);
t38 = t88 * t86 + t106;
t58 = t132 * t82;
t13 = t85 * t38 + t81 * t58;
t55 = mrSges(6,1) * t81 + mrSges(6,2) * t85;
t117 = t55 + mrSges(5,3);
t116 = t131 * pkin(8) ^ 2;
t59 = t132 * t86;
t115 = t81 ^ 2 + t85 ^ 2;
t114 = qJ(4) * t35;
t80 = sin(qJ(6));
t84 = cos(qJ(6));
t43 = -t80 * t85 - t84 * t81;
t96 = t80 * t81 - t84 * t85;
t112 = t43 ^ 2 + t96 ^ 2;
t28 = t96 * t86;
t29 = t43 * t86;
t111 = Ifges(7,5) * t29 + Ifges(7,6) * t28 + Ifges(7,3) * t82;
t18 = t81 * t121 + t33 * t85;
t19 = -t85 * t121 + t33 * t81;
t5 = t18 * t84 - t19 * t80;
t6 = t18 * t80 + t19 * t84;
t110 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t109 = m(6) * t115;
t108 = t115 * mrSges(6,3);
t107 = -mrSges(7,1) * t96 + t43 * mrSges(7,2);
t11 = -pkin(10) * t119 + t13;
t47 = t85 * t58;
t9 = pkin(5) * t82 + t47 + (pkin(10) * t86 - t38) * t81;
t2 = -t11 * t80 + t84 * t9;
t3 = t11 * t84 + t80 * t9;
t104 = -t2 * t96 - t3 * t43;
t103 = t43 * t6 + t5 * t96;
t102 = -t81 * mrSges(6,2) + t120;
t101 = -Ifges(6,5) * t81 - Ifges(6,6) * t85;
t12 = -t38 * t81 + t47;
t100 = t12 * t85 + t13 * t81;
t99 = t18 * t85 + t19 * t81;
t50 = t126 * t81;
t51 = t126 * t85;
t20 = -t50 * t80 + t51 * t84;
t21 = t50 * t84 + t51 * t80;
t98 = -t20 * t96 - t21 * t43;
t97 = t43 * t80 + t84 * t96;
t40 = Ifges(7,6) * t43;
t41 = Ifges(7,5) * t96;
t95 = t20 * mrSges(7,1) - t21 * mrSges(7,2) + t40 - t41;
t94 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t111;
t93 = (mrSges(7,1) * t84 - mrSges(7,2) * t80) * pkin(5);
t92 = (t35 * t86 + t123) * pkin(8);
t89 = qJ(4) ^ 2;
t72 = t78 ^ 2;
t67 = Ifges(6,5) * t85;
t66 = Ifges(6,3) * t82;
t62 = pkin(5) * t81 + qJ(4);
t60 = t72 * t87 ^ 2;
t57 = Ifges(6,1) * t85 - t125;
t56 = -Ifges(6,2) * t81 + t124;
t54 = -mrSges(4,1) * t86 + mrSges(4,2) * t82;
t53 = mrSges(5,2) * t86 - mrSges(5,3) * t82;
t52 = -pkin(3) * t86 + t106;
t49 = -mrSges(6,2) * t82 - mrSges(6,3) * t119;
t48 = mrSges(6,3) * t81 * t86 + mrSges(6,1) * t82;
t37 = t102 * t86;
t32 = pkin(5) * t119 + t59;
t27 = Ifges(6,5) * t82 + (-Ifges(6,1) * t81 - t124) * t86;
t26 = Ifges(6,6) * t82 + (-Ifges(6,2) * t85 - t125) * t86;
t23 = mrSges(7,1) * t82 - mrSges(7,3) * t29;
t22 = -mrSges(7,2) * t82 + mrSges(7,3) * t28;
t17 = -Ifges(7,1) * t96 + Ifges(7,4) * t43;
t16 = -Ifges(7,4) * t96 + Ifges(7,2) * t43;
t15 = -mrSges(7,1) * t43 - mrSges(7,2) * t96;
t10 = -mrSges(7,1) * t28 + mrSges(7,2) * t29;
t8 = Ifges(7,1) * t29 + Ifges(7,4) * t28 + Ifges(7,5) * t82;
t7 = Ifges(7,4) * t29 + Ifges(7,2) * t28 + Ifges(7,6) * t82;
t1 = [m(2) + m(7) * (t5 ^ 2 + t6 ^ 2 + t31) + m(6) * (t18 ^ 2 + t19 ^ 2 + t31) + m(3) * (t72 * t83 ^ 2 + t79 ^ 2 + t60) + (m(5) + m(4)) * (t33 ^ 2 + t31 + t60); t18 * t48 + t19 * t49 + t6 * t22 + t5 * t23 + t118 * t123 + (-t83 * mrSges(3,2) + (mrSges(3,1) - t53 - t54) * t87) * t78 + (t118 * t86 + t10 + t37) * t35 + m(7) * (t2 * t5 + t3 * t6 + t32 * t35) + m(6) * (t12 * t18 + t13 * t19 + t35 * t59) + m(5) * (-t52 * t121 + t92) + m(4) * (pkin(2) * t121 + t92); -0.2e1 * pkin(2) * t54 + 0.2e1 * t32 * t10 + 0.2e1 * t12 * t48 + 0.2e1 * t13 * t49 + 0.2e1 * t2 * t23 + 0.2e1 * t3 * t22 + t28 * t7 + t29 * t8 + 0.2e1 * t59 * t37 + 0.2e1 * t52 * t53 + Ifges(3,3) + (t66 + (Ifges(5,2) + Ifges(4,1)) * t82 + t111) * t82 + m(5) * (t52 ^ 2 + t116) + m(4) * (pkin(2) ^ 2 + t116) + m(6) * (t12 ^ 2 + t13 ^ 2 + t59 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t32 ^ 2) + (-t85 * t26 - t81 * t27 + (Ifges(4,2) + Ifges(5,3)) * t86 + ((2 * Ifges(4,4)) + (2 * Ifges(5,6)) + t101) * t82) * t86 + 0.2e1 * t118 * pkin(8) * t131; (-mrSges(4,1) + mrSges(5,2)) * t33 + t103 * mrSges(7,3) - t99 * mrSges(6,3) + (-mrSges(4,2) + t15 + t117) * t35 + m(7) * (t20 * t5 + t21 * t6 + t35 * t62) + m(6) * (t99 * t88 + t114) + m(5) * (-pkin(3) * t33 + t114); -t96 * t8 / 0.2e1 + t59 * t55 + t62 * t10 + qJ(4) * t37 + t43 * t7 / 0.2e1 + t28 * t16 / 0.2e1 + t29 * t17 / 0.2e1 + t32 * t15 + t21 * t22 + t20 * t23 - t104 * mrSges(7,3) + (t88 * t48 - t12 * mrSges(6,3) + t27 / 0.2e1) * t85 + (t88 * t49 - t13 * mrSges(6,3) - t26 / 0.2e1) * t81 + m(6) * (qJ(4) * t59 + t100 * t88) + m(7) * (t2 * t20 + t21 * t3 + t32 * t62) + (-pkin(3) * mrSges(5,1) + Ifges(6,6) * t127 + t67 / 0.2e1 - t41 / 0.2e1 + t40 / 0.2e1 + Ifges(4,5) - Ifges(5,4)) * t82 + (-t85 * t56 / 0.2e1 + t57 * t127 + qJ(4) * mrSges(5,1) + Ifges(4,6) - Ifges(5,5)) * t86 + ((m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t86 + (-mrSges(4,1) + t129) * t82) * pkin(8); -0.2e1 * pkin(3) * mrSges(5,2) + 0.2e1 * t62 * t15 + t43 * t16 - t96 * t17 - t81 * t56 + t85 * t57 + Ifges(5,1) + Ifges(4,3) + m(7) * (t20 ^ 2 + t21 ^ 2 + t62 ^ 2) + m(6) * (t115 * t88 ^ 2 + t89) + m(5) * (pkin(3) ^ 2 + t89) - 0.2e1 * t98 * mrSges(7,3) + 0.2e1 * t117 * qJ(4) - 0.2e1 * t88 * t108; m(5) * t33 + m(6) * t99 - m(7) * t103; -t43 * t22 - t96 * t23 + t85 * t48 + t81 * t49 + (m(5) * pkin(8) + mrSges(5,1)) * t82 + m(7) * t104 + m(6) * t100; m(7) * t98 - t112 * mrSges(7,3) + t88 * t109 - t108 + t129; m(7) * t112 + m(5) + t109; t18 * mrSges(6,1) - t19 * mrSges(6,2) + (t5 * t84 + t6 * t80) * t128 + t110; t12 * mrSges(6,1) - t13 * mrSges(6,2) + t66 + t101 * t86 + (t84 * t23 + m(7) * (t2 * t84 + t3 * t80) + t80 * t22) * pkin(5) + t94; t88 * t120 + t67 + (-mrSges(6,2) * t88 - Ifges(6,6)) * t81 + (m(7) * (t20 * t84 + t21 * t80) + t97 * mrSges(7,3)) * pkin(5) + t95; -t97 * t128 + t102 + t107; Ifges(6,3) + Ifges(7,3) + m(7) * (t80 ^ 2 + t84 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t93; t110; t94; t95; t107; Ifges(7,3) + t93; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
