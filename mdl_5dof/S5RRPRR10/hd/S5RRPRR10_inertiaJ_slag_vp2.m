% Calculate joint inertia matrix for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR10_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR10_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR10_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:23:50
% EndTime: 2019-12-31 20:23:53
% DurationCPUTime: 0.92s
% Computational Cost: add. (1492->240), mult. (3611->356), div. (0->0), fcn. (3795->10), ass. (0->103)
t134 = Ifges(3,3) + Ifges(4,3);
t89 = sin(pkin(10));
t77 = pkin(2) * t89 + pkin(8);
t133 = 0.2e1 * t77;
t90 = sin(pkin(5));
t98 = cos(qJ(2));
t116 = t90 * t98;
t95 = sin(qJ(2));
t117 = t90 * t95;
t91 = cos(pkin(10));
t48 = -t91 * t116 + t89 * t117;
t49 = (t89 * t98 + t91 * t95) * t90;
t132 = Ifges(4,5) * t49 - Ifges(4,6) * t48;
t92 = cos(pkin(5));
t94 = sin(qJ(4));
t97 = cos(qJ(4));
t37 = t49 * t97 + t92 * t94;
t93 = sin(qJ(5));
t96 = cos(qJ(5));
t25 = -t37 * t93 + t48 * t96;
t36 = t49 * t94 - t92 * t97;
t12 = -mrSges(6,2) * t36 + mrSges(6,3) * t25;
t26 = t37 * t96 + t48 * t93;
t13 = mrSges(6,1) * t36 - mrSges(6,3) * t26;
t131 = t96 * t12 - t93 * t13;
t57 = (-pkin(2) * t98 - pkin(1)) * t90;
t130 = 0.2e1 * t57;
t129 = t25 / 0.2e1;
t128 = t26 / 0.2e1;
t127 = -t93 / 0.2e1;
t126 = t93 / 0.2e1;
t125 = t96 / 0.2e1;
t124 = pkin(1) * t92;
t123 = pkin(4) * t97;
t122 = Ifges(6,4) * t93;
t121 = Ifges(6,4) * t96;
t73 = t98 * t124;
t53 = -pkin(7) * t117 + t73;
t120 = t53 * mrSges(3,1);
t54 = pkin(7) * t116 + t95 * t124;
t119 = t54 * mrSges(3,2);
t118 = t77 * t97;
t114 = t93 * t94;
t113 = t94 * t96;
t38 = pkin(2) * t92 + t73 + (-pkin(7) - qJ(3)) * t117;
t43 = qJ(3) * t116 + t54;
t23 = t89 * t38 + t91 * t43;
t19 = pkin(8) * t92 + t23;
t27 = pkin(3) * t48 - pkin(8) * t49 + t57;
t10 = t97 * t19 + t94 * t27;
t11 = -mrSges(6,1) * t25 + mrSges(6,2) * t26;
t29 = mrSges(5,1) * t48 - mrSges(5,3) * t37;
t111 = -t29 + t11;
t62 = Ifges(6,5) * t93 + Ifges(6,6) * t96;
t110 = Ifges(5,5) * t94 + Ifges(5,6) * t97;
t109 = t93 ^ 2 + t96 ^ 2;
t86 = t94 ^ 2;
t88 = t97 ^ 2;
t108 = t86 + t88;
t3 = Ifges(6,5) * t26 + Ifges(6,6) * t25 + Ifges(6,3) * t36;
t107 = Ifges(5,5) * t37 - Ifges(5,6) * t36 + Ifges(5,3) * t48;
t78 = -pkin(2) * t91 - pkin(3);
t106 = t109 * t94;
t22 = t38 * t91 - t89 * t43;
t7 = pkin(9) * t48 + t10;
t18 = -pkin(3) * t92 - t22;
t8 = pkin(4) * t36 - pkin(9) * t37 + t18;
t1 = -t7 * t93 + t8 * t96;
t2 = t7 * t96 + t8 * t93;
t105 = -t1 * t93 + t2 * t96;
t61 = -t97 * mrSges(5,1) + t94 * mrSges(5,2);
t104 = mrSges(6,1) * t93 + mrSges(6,2) * t96;
t9 = -t19 * t94 + t27 * t97;
t56 = -pkin(9) * t94 - t123 + t78;
t34 = -t93 * t118 + t56 * t96;
t35 = t96 * t118 + t56 * t93;
t103 = -t34 * t93 + t35 * t96;
t58 = mrSges(6,2) * t97 - mrSges(6,3) * t114;
t59 = -mrSges(6,1) * t97 - mrSges(6,3) * t113;
t102 = t96 * t58 - t93 * t59;
t101 = Ifges(3,5) * t117 + Ifges(3,6) * t116 + t134 * t92 + t132;
t50 = Ifges(6,5) * t113 - Ifges(6,6) * t114 - Ifges(6,3) * t97;
t75 = t77 ^ 2;
t67 = t86 * t75;
t66 = Ifges(5,1) * t94 + Ifges(5,4) * t97;
t65 = Ifges(6,1) * t93 + t121;
t64 = Ifges(5,4) * t94 + Ifges(5,2) * t97;
t63 = Ifges(6,2) * t96 + t122;
t60 = -mrSges(6,1) * t96 + mrSges(6,2) * t93;
t55 = t104 * t94;
t52 = -Ifges(6,5) * t97 + (Ifges(6,1) * t96 - t122) * t94;
t51 = -Ifges(6,6) * t97 + (-Ifges(6,2) * t93 + t121) * t94;
t44 = t49 * mrSges(4,2);
t40 = mrSges(4,1) * t92 - mrSges(4,3) * t49;
t39 = -mrSges(4,2) * t92 - mrSges(4,3) * t48;
t28 = -mrSges(5,2) * t48 - mrSges(5,3) * t36;
t17 = mrSges(5,1) * t36 + mrSges(5,2) * t37;
t15 = Ifges(5,1) * t37 - Ifges(5,4) * t36 + Ifges(5,5) * t48;
t14 = Ifges(5,4) * t37 - Ifges(5,2) * t36 + Ifges(5,6) * t48;
t6 = -pkin(4) * t48 - t9;
t5 = Ifges(6,1) * t26 + Ifges(6,4) * t25 + Ifges(6,5) * t36;
t4 = Ifges(6,4) * t26 + Ifges(6,2) * t25 + Ifges(6,6) * t36;
t16 = [Ifges(4,1) * t49 ^ 2 + 0.2e1 * t1 * t13 + 0.2e1 * t10 * t28 + 0.2e1 * t6 * t11 + 0.2e1 * t2 * t12 + t37 * t15 + 0.2e1 * t18 * t17 + 0.2e1 * t22 * t40 + 0.2e1 * t23 * t39 + t25 * t4 + t26 * t5 + 0.2e1 * t9 * t29 + t44 * t130 + Ifges(2,3) + (mrSges(4,1) * t130 - 0.2e1 * Ifges(4,4) * t49 + Ifges(4,2) * t48 + t107) * t48 + (t3 - t14) * t36 + (t101 - 0.2e1 * t119 + 0.2e1 * t120 + t132) * t92 + m(3) * (t53 ^ 2 + t54 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2 + t57 ^ 2) + m(5) * (t10 ^ 2 + t18 ^ 2 + t9 ^ 2) + m(6) * (t1 ^ 2 + t2 ^ 2 + t6 ^ 2) + ((Ifges(3,5) * t95 + Ifges(3,6) * t98) * t92 + 0.2e1 * (-t53 * t95 + t54 * t98) * mrSges(3,3) + (-0.2e1 * pkin(1) * (-mrSges(3,1) * t98 + mrSges(3,2) * t95) + t98 * (Ifges(3,4) * t95 + Ifges(3,2) * t98) + t95 * (Ifges(3,1) * t95 + Ifges(3,4) * t98) + m(3) * pkin(1) ^ 2) * t90) * t90; (t14 / 0.2e1 - t3 / 0.2e1 + t77 * t28 + t10 * mrSges(5,3)) * t97 + m(6) * (t6 * t77 * t94 + t1 * t34 + t2 * t35) + (-t64 / 0.2e1 + t50 / 0.2e1) * t36 + t37 * t66 / 0.2e1 + t78 * t17 + t48 * t110 / 0.2e1 + t51 * t129 + t52 * t128 + t120 - t119 + t6 * t55 + t2 * t58 + t1 * t59 + t18 * t61 + t34 * t13 + t35 * t12 - t23 * mrSges(4,2) + t22 * mrSges(4,1) + t101 + (t15 / 0.2e1 + t4 * t127 - t9 * mrSges(5,3) + t5 * t125 + t111 * t77) * t94 + m(5) * (t18 * t78 + (t10 * t97 - t9 * t94) * t77) + (t89 * t39 + m(4) * (t22 * t91 + t23 * t89) + t91 * t40) * pkin(2); 0.2e1 * t34 * t59 + 0.2e1 * t35 * t58 + 0.2e1 * t78 * t61 + (-t50 + t64) * t97 + m(6) * (t34 ^ 2 + t35 ^ 2 + t67) + m(5) * (t75 * t88 + t78 ^ 2 + t67) + m(4) * (t89 ^ 2 + t91 ^ 2) * pkin(2) ^ 2 + (t55 * t133 - t51 * t93 + t52 * t96 + t66) * t94 + 0.2e1 * (mrSges(4,1) * t91 - mrSges(4,2) * t89) * pkin(2) + t108 * mrSges(5,3) * t133 + t134; t48 * mrSges(4,1) + t44 - t111 * t97 + (t28 + t131) * t94 + m(6) * (t105 * t94 - t6 * t97) + m(5) * (t10 * t94 + t9 * t97) + m(4) * t57; -t97 * t55 + (m(6) * (t103 - t118) + t102) * t94; m(4) + m(5) * t108 + m(6) * (t109 * t86 + t88); t6 * t60 + t5 * t126 + t4 * t125 + t65 * t128 + t63 * t129 + t36 * t62 / 0.2e1 + t9 * mrSges(5,1) - t10 * mrSges(5,2) + t105 * mrSges(6,3) + t107 + (-m(6) * t6 - t11) * pkin(4) + (m(6) * t105 + t131) * pkin(9); t52 * t126 + t51 * t125 - pkin(4) * t55 + (m(6) * t103 + t102) * pkin(9) + (t63 * t127 + t65 * t125 + (-m(6) * pkin(4) - mrSges(5,1) + t60) * t77) * t94 + (-t77 * mrSges(5,2) - t62 / 0.2e1) * t97 + t103 * mrSges(6,3) + t110; -t97 * t60 + m(6) * (pkin(9) * t106 + t123) + mrSges(6,3) * t106 - t61; Ifges(5,3) + t96 * t63 - 0.2e1 * pkin(4) * t60 + t93 * t65 + m(6) * (t109 * pkin(9) ^ 2 + pkin(4) ^ 2) + 0.2e1 * t109 * pkin(9) * mrSges(6,3); mrSges(6,1) * t1 - mrSges(6,2) * t2 + t3; mrSges(6,1) * t34 - mrSges(6,2) * t35 + t50; -t55; -t104 * pkin(9) + t62; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t16(1), t16(2), t16(4), t16(7), t16(11); t16(2), t16(3), t16(5), t16(8), t16(12); t16(4), t16(5), t16(6), t16(9), t16(13); t16(7), t16(8), t16(9), t16(10), t16(14); t16(11), t16(12), t16(13), t16(14), t16(15);];
Mq = res;
