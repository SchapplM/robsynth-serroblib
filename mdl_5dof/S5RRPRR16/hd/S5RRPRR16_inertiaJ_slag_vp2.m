% Calculate joint inertia matrix for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR16_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR16_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:45
% EndTime: 2019-12-31 20:44:47
% DurationCPUTime: 0.90s
% Computational Cost: add. (983->239), mult. (2178->338), div. (0->0), fcn. (2070->8), ass. (0->100)
t127 = Ifges(4,1) + Ifges(3,3);
t88 = (-pkin(2) - pkin(8));
t126 = -2 * t88;
t80 = sin(pkin(5));
t84 = sin(qJ(2));
t111 = t80 * t84;
t87 = cos(qJ(2));
t110 = t80 * t87;
t81 = cos(pkin(5));
t83 = sin(qJ(4));
t86 = cos(qJ(4));
t40 = -t83 * t110 + t81 * t86;
t82 = sin(qJ(5));
t85 = cos(qJ(5));
t22 = t85 * t111 - t40 * t82;
t39 = t86 * t110 + t81 * t83;
t12 = -t39 * mrSges(6,2) + t22 * mrSges(6,3);
t23 = t82 * t111 + t40 * t85;
t13 = t39 * mrSges(6,1) - t23 * mrSges(6,3);
t125 = t85 * t12 - t82 * t13;
t49 = -t85 * mrSges(6,1) + t82 * mrSges(6,2);
t124 = m(6) * pkin(4) + mrSges(5,1) - t49;
t123 = t22 / 0.2e1;
t122 = t23 / 0.2e1;
t51 = Ifges(6,5) * t82 + Ifges(6,6) * t85;
t121 = t51 / 0.2e1;
t120 = -t82 / 0.2e1;
t119 = t82 / 0.2e1;
t118 = t85 / 0.2e1;
t117 = pkin(1) * t87;
t100 = -pkin(2) - t117;
t60 = pkin(7) * t111;
t18 = pkin(3) * t111 + t60 + (-pkin(8) + t100) * t81;
t97 = -qJ(3) * t84 - pkin(1);
t25 = (t88 * t87 + t97) * t80;
t8 = t86 * t18 - t83 * t25;
t6 = -pkin(4) * t111 - t8;
t116 = t86 * t6;
t115 = Ifges(6,4) * t82;
t114 = Ifges(6,4) * t85;
t41 = t81 * t117 - t60;
t113 = t41 * mrSges(3,1);
t42 = t81 * t84 * pkin(1) + pkin(7) * t110;
t112 = t42 * mrSges(3,2);
t108 = t82 * t86;
t107 = t83 * t88;
t105 = t85 * t86;
t11 = -t22 * mrSges(6,1) + t23 * mrSges(6,2);
t27 = mrSges(5,1) * t111 - t40 * mrSges(5,3);
t104 = -t11 + t27;
t9 = t83 * t18 + t86 * t25;
t45 = mrSges(4,1) * t111 + t81 * mrSges(4,2);
t103 = t82 ^ 2 + t85 ^ 2;
t76 = t83 ^ 2;
t78 = t86 ^ 2;
t102 = t76 + t78;
t3 = Ifges(6,5) * t23 + Ifges(6,6) * t22 + Ifges(6,3) * t39;
t101 = Ifges(5,5) * t40 - Ifges(5,6) * t39 + Ifges(5,3) * t111;
t30 = -t81 * qJ(3) - t42;
t98 = t102 * mrSges(5,3);
t96 = Ifges(3,5) * t111 + Ifges(3,6) * t110 + t127 * t81;
t24 = pkin(3) * t110 - t30;
t10 = t39 * pkin(4) - t40 * pkin(9) + t24;
t7 = pkin(9) * t111 + t9;
t1 = t85 * t10 - t82 * t7;
t2 = t82 * t10 + t85 * t7;
t95 = -t1 * t82 + t2 * t85;
t94 = t86 * t8 + t83 * t9;
t93 = mrSges(6,1) * t82 + mrSges(6,2) * t85;
t48 = t83 * pkin(4) - t86 * pkin(9) + qJ(3);
t28 = -t82 * t107 + t85 * t48;
t29 = t85 * t107 + t82 * t48;
t92 = -t28 * t82 + t29 * t85;
t46 = -t83 * mrSges(6,2) - mrSges(6,3) * t108;
t47 = t83 * mrSges(6,1) - mrSges(6,3) * t105;
t91 = t85 * t46 - t82 * t47;
t33 = Ifges(6,5) * t105 - Ifges(6,6) * t108 + Ifges(6,3) * t83;
t89 = qJ(3) ^ 2;
t79 = t88 ^ 2;
t74 = Ifges(5,5) * t86;
t66 = t78 * t88;
t65 = t78 * t79;
t55 = Ifges(5,1) * t86 - Ifges(5,4) * t83;
t54 = Ifges(6,1) * t82 + t114;
t53 = Ifges(5,4) * t86 - Ifges(5,2) * t83;
t52 = Ifges(6,2) * t85 + t115;
t50 = t83 * mrSges(5,1) + t86 * mrSges(5,2);
t44 = -mrSges(4,1) * t110 - t81 * mrSges(4,3);
t43 = t93 * t86;
t35 = Ifges(6,5) * t83 + (Ifges(6,1) * t85 - t115) * t86;
t34 = Ifges(6,6) * t83 + (-Ifges(6,2) * t82 + t114) * t86;
t32 = t100 * t81 + t60;
t31 = (-pkin(2) * t87 + t97) * t80;
t26 = -mrSges(5,2) * t111 - t39 * mrSges(5,3);
t16 = t39 * mrSges(5,1) + t40 * mrSges(5,2);
t15 = Ifges(5,1) * t40 - Ifges(5,4) * t39 + Ifges(5,5) * t111;
t14 = Ifges(5,4) * t40 - Ifges(5,2) * t39 + Ifges(5,6) * t111;
t5 = Ifges(6,1) * t23 + Ifges(6,4) * t22 + Ifges(6,5) * t39;
t4 = Ifges(6,4) * t23 + Ifges(6,2) * t22 + Ifges(6,6) * t39;
t17 = [0.2e1 * t1 * t13 + 0.2e1 * t6 * t11 + 0.2e1 * t2 * t12 + t40 * t15 + 0.2e1 * t24 * t16 + t22 * t4 + t23 * t5 + 0.2e1 * t9 * t26 + 0.2e1 * t8 * t27 + 0.2e1 * t30 * t44 + 0.2e1 * t32 * t45 + Ifges(2,3) + (t3 - t14) * t39 + (-0.2e1 * t112 + t96 + 0.2e1 * t113) * t81 + m(3) * (t41 ^ 2 + t42 ^ 2) + m(6) * (t1 ^ 2 + t2 ^ 2 + t6 ^ 2) + m(5) * (t24 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(4) * (t30 ^ 2 + t31 ^ 2 + t32 ^ 2) + (0.2e1 * t31 * (mrSges(4,2) * t87 - mrSges(4,3) * t84) + t84 * t101 + (m(3) * pkin(1) ^ 2 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + Ifges(4,3)) * t87) * t87 + (-0.2e1 * pkin(1) * mrSges(3,2) + (Ifges(4,2) + Ifges(3,1)) * t84 + 0.2e1 * (Ifges(3,4) + Ifges(4,6)) * t87) * t84) * t80 + 0.2e1 * (-t41 * t84 + t42 * t87) * mrSges(3,3) + ((-(2 * Ifges(4,5)) + Ifges(3,6)) * t87 + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t84) * t81) * t80; (t15 / 0.2e1 + t5 * t118 + t4 * t120 - t8 * mrSges(5,3) + t104 * t88) * t86 + (-t14 / 0.2e1 + t3 / 0.2e1 + t88 * t26 - t9 * mrSges(5,3)) * t83 + (-t44 + t16) * qJ(3) + m(6) * (t28 * t1 - t88 * t116 + t29 * t2) + m(4) * (-pkin(2) * t32 - qJ(3) * t30) + (-Ifges(4,5) * t87 + t84 * (-Ifges(5,6) * t83 + t74) / 0.2e1 - Ifges(4,4) * t84) * t80 + (-t53 / 0.2e1 + t33 / 0.2e1) * t39 + m(5) * (qJ(3) * t24 + t94 * t88) + t96 + t40 * t55 / 0.2e1 + t113 - t112 + t6 * t43 - pkin(2) * t45 + t2 * t46 + t1 * t47 + t24 * t50 - t30 * mrSges(4,3) + t32 * mrSges(4,2) + t34 * t123 + t35 * t122 + t28 * t13 + t29 * t12; -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t28 * t47 + 0.2e1 * t29 * t46 + (t33 - t53) * t83 + m(6) * (t28 ^ 2 + t29 ^ 2 + t65) + m(5) * (t76 * t79 + t65 + t89) + m(4) * (pkin(2) ^ 2 + t89) + (t43 * t126 - t82 * t34 + t85 * t35 + t55) * t86 + 0.2e1 * (t50 + mrSges(4,3)) * qJ(3) + t98 * t126 + t127; t104 * t86 + (t26 + t125) * t83 + m(6) * (t95 * t83 - t116) + m(5) * t94 + m(4) * t32 + t45; -m(4) * pkin(2) - t86 * t43 + mrSges(4,2) + t91 * t83 - t98 + m(6) * (t92 * t83 + t66) + m(5) * (t76 * t88 + t66); m(4) + m(5) * t102 + m(6) * (t103 * t76 + t78); t8 * mrSges(5,1) - t9 * mrSges(5,2) + t95 * mrSges(6,3) + t4 * t118 + t5 * t119 + t39 * t121 + t54 * t122 + t52 * t123 + t6 * t49 + t101 + (-m(6) * t6 - t11) * pkin(4) + (m(6) * t95 + t125) * pkin(9); t35 * t119 + t34 * t118 - pkin(4) * t43 + t74 + (m(6) * t92 + t91) * pkin(9) + (t54 * t118 + t52 * t120 + t124 * t88) * t86 + t92 * mrSges(6,3) + (-t88 * mrSges(5,2) - Ifges(5,6) + t121) * t83; t124 * t86 + (-mrSges(5,2) + (m(6) * pkin(9) + mrSges(6,3)) * t103) * t83; Ifges(5,3) - 0.2e1 * pkin(4) * t49 + t82 * t54 + t85 * t52 + m(6) * (t103 * pkin(9) ^ 2 + pkin(4) ^ 2) + 0.2e1 * t103 * pkin(9) * mrSges(6,3); t1 * mrSges(6,1) - t2 * mrSges(6,2) + t3; t28 * mrSges(6,1) - t29 * mrSges(6,2) + t33; -t93 * t83; -t93 * pkin(9) + t51; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t17(1), t17(2), t17(4), t17(7), t17(11); t17(2), t17(3), t17(5), t17(8), t17(12); t17(4), t17(5), t17(6), t17(9), t17(13); t17(7), t17(8), t17(9), t17(10), t17(14); t17(11), t17(12), t17(13), t17(14), t17(15);];
Mq = res;
