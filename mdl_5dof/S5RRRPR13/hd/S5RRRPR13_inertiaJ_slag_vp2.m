% Calculate joint inertia matrix for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR13_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR13_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR13_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:19
% EndTime: 2019-12-31 21:43:23
% DurationCPUTime: 1.01s
% Computational Cost: add. (1084->264), mult. (2423->367), div. (0->0), fcn. (2358->8), ass. (0->101)
t90 = sin(qJ(3));
t93 = cos(qJ(3));
t131 = t90 ^ 2 + t93 ^ 2;
t89 = sin(qJ(5));
t126 = -t89 / 0.2e1;
t130 = -m(5) * pkin(3) + mrSges(5,2);
t87 = sin(pkin(5));
t94 = cos(qJ(2));
t115 = t87 * t94;
t91 = sin(qJ(2));
t116 = t87 * t91;
t88 = cos(pkin(5));
t47 = t90 * t116 - t88 * t93;
t92 = cos(qJ(5));
t27 = t89 * t115 + t47 * t92;
t129 = t27 / 0.2e1;
t119 = Ifges(6,4) * t92;
t41 = Ifges(6,5) * t90 + (-Ifges(6,1) * t89 - t119) * t93;
t128 = t41 / 0.2e1;
t78 = Ifges(6,5) * t92;
t127 = Ifges(6,6) * t126 + t78 / 0.2e1;
t125 = -t92 / 0.2e1;
t124 = pkin(3) + pkin(9);
t123 = pkin(4) + pkin(8);
t122 = pkin(1) * t94;
t121 = mrSges(6,3) * t93;
t120 = Ifges(6,4) * t89;
t70 = pkin(7) * t116;
t49 = t88 * t122 - t70;
t118 = t49 * mrSges(3,1);
t50 = t88 * t91 * pkin(1) + pkin(7) * t115;
t117 = t50 * mrSges(3,2);
t54 = -t90 * mrSges(6,2) - t92 * t121;
t114 = t89 * t54;
t113 = t92 * t124;
t112 = Ifges(5,1) + Ifges(4,3);
t37 = t88 * pkin(8) + t50;
t38 = (-pkin(2) * t94 - pkin(8) * t91 - pkin(1)) * t87;
t16 = t93 * t37 + t90 * t38;
t111 = Ifges(4,5) * t90 + Ifges(4,6) * t93;
t110 = t131 * pkin(8) ^ 2;
t109 = t89 ^ 2 + t92 ^ 2;
t28 = -t92 * t115 + t47 * t89;
t48 = t93 * t116 + t88 * t90;
t3 = Ifges(6,5) * t28 + Ifges(6,6) * t27 + Ifges(6,3) * t48;
t108 = Ifges(3,5) * t116 + Ifges(3,6) * t115 + Ifges(3,3) * t88;
t107 = t115 / 0.2e1;
t106 = m(6) * t109;
t105 = t109 * mrSges(6,3);
t104 = -t90 * qJ(4) - pkin(2);
t15 = -t90 * t37 + t93 * t38;
t103 = (Ifges(5,4) - Ifges(4,5)) * t48 + (-Ifges(5,5) + Ifges(4,6)) * t47;
t12 = pkin(3) * t115 - t15;
t6 = pkin(4) * t48 + pkin(9) * t115 + t12;
t36 = t70 + (-pkin(2) - t122) * t88;
t98 = -t48 * qJ(4) + t36;
t7 = t124 * t47 + t98;
t1 = t6 * t92 - t7 * t89;
t2 = t6 * t89 + t7 * t92;
t101 = t1 * t92 + t2 * t89;
t100 = t92 * mrSges(6,1) - t89 * mrSges(6,2);
t52 = -t124 * t93 + t104;
t66 = t123 * t90;
t23 = -t89 * t52 + t92 * t66;
t24 = t92 * t52 + t89 * t66;
t99 = t23 * t92 + t24 * t89;
t11 = qJ(4) * t115 - t16;
t30 = t48 * mrSges(5,1) - mrSges(5,2) * t115;
t39 = Ifges(6,3) * t90 + (-Ifges(6,5) * t89 - Ifges(6,6) * t92) * t93;
t96 = qJ(4) ^ 2;
t67 = t123 * t93;
t65 = Ifges(4,1) * t90 + Ifges(4,4) * t93;
t64 = Ifges(6,1) * t92 - t120;
t63 = Ifges(4,4) * t90 + Ifges(4,2) * t93;
t62 = -Ifges(6,2) * t89 + t119;
t60 = -Ifges(5,2) * t90 - Ifges(5,6) * t93;
t59 = -Ifges(5,6) * t90 - Ifges(5,3) * t93;
t58 = t89 * mrSges(6,1) + t92 * mrSges(6,2);
t57 = -t93 * mrSges(4,1) + t90 * mrSges(4,2);
t56 = t93 * mrSges(5,2) - t90 * mrSges(5,3);
t55 = -t93 * pkin(3) + t104;
t53 = t90 * mrSges(6,1) + t89 * t121;
t51 = t100 * t93;
t40 = Ifges(6,6) * t90 + (-Ifges(6,2) * t92 - t120) * t93;
t32 = -mrSges(4,1) * t115 - t48 * mrSges(4,3);
t31 = mrSges(4,2) * t115 - t47 * mrSges(4,3);
t29 = t47 * mrSges(5,1) + mrSges(5,3) * t115;
t22 = -t47 * mrSges(5,2) - t48 * mrSges(5,3);
t21 = t47 * mrSges(4,1) + t48 * mrSges(4,2);
t20 = Ifges(4,1) * t48 - Ifges(4,4) * t47 - Ifges(4,5) * t115;
t19 = Ifges(4,4) * t48 - Ifges(4,2) * t47 - Ifges(4,6) * t115;
t18 = -Ifges(5,4) * t115 - Ifges(5,2) * t48 + Ifges(5,6) * t47;
t17 = -Ifges(5,5) * t115 - Ifges(5,6) * t48 + Ifges(5,3) * t47;
t14 = t48 * mrSges(6,1) - t28 * mrSges(6,3);
t13 = -t48 * mrSges(6,2) + t27 * mrSges(6,3);
t10 = t47 * pkin(3) + t98;
t9 = -t27 * mrSges(6,1) + t28 * mrSges(6,2);
t8 = -t47 * pkin(4) - t11;
t5 = Ifges(6,1) * t28 + Ifges(6,4) * t27 + Ifges(6,5) * t48;
t4 = Ifges(6,4) * t28 + Ifges(6,2) * t27 + Ifges(6,6) * t48;
t25 = [0.2e1 * t1 * t14 + 0.2e1 * t10 * t22 + 0.2e1 * t11 * t29 + 0.2e1 * t12 * t30 + 0.2e1 * t2 * t13 + 0.2e1 * t15 * t32 + 0.2e1 * t16 * t31 + 0.2e1 * t36 * t21 + t27 * t4 + t28 * t5 + 0.2e1 * t8 * t9 + Ifges(2,3) + (t17 - t19) * t47 + (t108 - 0.2e1 * t117 + 0.2e1 * t118) * t88 + (t20 + t3 - t18) * t48 + ((-0.2e1 * t49 * mrSges(3,3) + Ifges(3,5) * t88 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t91) * t87) * t91 + (0.2e1 * t50 * mrSges(3,3) + Ifges(3,6) * t88 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t91 + (Ifges(3,2) + t112) * t94) * t87 + t103) * t94) * t87 + m(6) * (t1 ^ 2 + t2 ^ 2 + t8 ^ 2) + m(5) * (t10 ^ 2 + t11 ^ 2 + t12 ^ 2) + m(4) * (t15 ^ 2 + t16 ^ 2 + t36 ^ 2) + m(3) * (pkin(1) ^ 2 * t87 ^ 2 + t49 ^ 2 + t50 ^ 2); m(4) * (-pkin(2) * t36 + (-t15 * t90 + t16 * t93) * pkin(8)) + (-t17 / 0.2e1 + t19 / 0.2e1 - t11 * mrSges(5,1) + t16 * mrSges(4,3) + t5 * t126 + t4 * t125 + Ifges(5,5) * t107 + (-t29 + t31) * pkin(8)) * t93 + (-t18 / 0.2e1 + t20 / 0.2e1 + t3 / 0.2e1 + t12 * mrSges(5,1) - t15 * mrSges(4,3) + Ifges(5,4) * t107 + (t30 - t32) * pkin(8)) * t90 + m(6) * (t1 * t23 + t2 * t24 + t67 * t8) + (-t60 / 0.2e1 + t65 / 0.2e1 + t39 / 0.2e1) * t48 + (t59 / 0.2e1 - t63 / 0.2e1) * t47 + t108 + t1 * t53 + t2 * t54 + t55 * t22 + t10 * t56 + t36 * t57 + t67 * t9 + t118 - t117 + t8 * t51 + t40 * t129 + t28 * t128 + t24 * t13 - pkin(2) * t21 + t23 * t14 + m(5) * (t55 * t10 + (-t11 * t93 + t12 * t90) * pkin(8)) - t111 * t115 / 0.2e1; -0.2e1 * pkin(2) * t57 + 0.2e1 * t23 * t53 + 0.2e1 * t24 * t54 + 0.2e1 * t67 * t51 + 0.2e1 * t55 * t56 + Ifges(3,3) + (t39 + t65 - t60) * t90 + m(6) * (t23 ^ 2 + t24 ^ 2 + t67 ^ 2) + m(5) * (t55 ^ 2 + t110) + m(4) * (pkin(2) ^ 2 + t110) + (-t40 * t92 - t41 * t89 - t59 + t63) * t93 + 0.2e1 * (mrSges(5,1) + mrSges(4,3)) * pkin(8) * t131; t8 * t58 + t48 * t127 + t62 * t129 + t28 * t64 / 0.2e1 - pkin(3) * t30 + t12 * mrSges(5,2) + t15 * mrSges(4,1) - t16 * mrSges(4,2) - t11 * mrSges(5,3) - t112 * t115 + (-t29 + t9) * qJ(4) + (t5 / 0.2e1 - t1 * mrSges(6,3) - t124 * t14) * t92 + (-t4 / 0.2e1 - t2 * mrSges(6,3) - t124 * t13) * t89 + m(5) * (-pkin(3) * t12 - qJ(4) * t11) + m(6) * (qJ(4) * t8 - t101 * t124) - t103; -t124 * t114 - t53 * t113 + t67 * t58 + t92 * t128 + t40 * t126 + qJ(4) * t51 + m(6) * (qJ(4) * t67 - t124 * t99) - t99 * mrSges(6,3) + (-pkin(3) * mrSges(5,1) - Ifges(5,4) + t127) * t90 + (qJ(4) * mrSges(5,1) + t62 * t125 + t64 * t126 - Ifges(5,5)) * t93 + ((m(5) * qJ(4) - mrSges(4,2) + mrSges(5,3)) * t93 + (-mrSges(4,1) + t130) * t90) * pkin(8) + t111; -0.2e1 * pkin(3) * mrSges(5,2) - t89 * t62 + t92 * t64 + m(6) * (t109 * t124 ^ 2 + t96) + m(5) * (pkin(3) ^ 2 + t96) + t112 + 0.2e1 * (t58 + mrSges(5,3)) * qJ(4) + 0.2e1 * t124 * t105; m(5) * t12 + m(6) * t101 + t89 * t13 + t92 * t14 + t30; m(6) * t99 + t114 + t92 * t53 + (m(5) * pkin(8) + mrSges(5,1)) * t90; -t106 * t124 - t105 + t130; m(5) + t106; mrSges(6,1) * t1 - mrSges(6,2) * t2 + t3; t23 * mrSges(6,1) - t24 * mrSges(6,2) + t39; -mrSges(6,1) * t113 + t78 + (mrSges(6,2) * t124 - Ifges(6,6)) * t89; t100; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t25(1), t25(2), t25(4), t25(7), t25(11); t25(2), t25(3), t25(5), t25(8), t25(12); t25(4), t25(5), t25(6), t25(9), t25(13); t25(7), t25(8), t25(9), t25(10), t25(14); t25(11), t25(12), t25(13), t25(14), t25(15);];
Mq = res;
