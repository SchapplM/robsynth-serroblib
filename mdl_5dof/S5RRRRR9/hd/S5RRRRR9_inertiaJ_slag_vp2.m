% Calculate joint inertia matrix for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR9_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR9_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:45
% EndTime: 2019-12-31 22:27:47
% DurationCPUTime: 0.85s
% Computational Cost: add. (1476->246), mult. (2979->367), div. (0->0), fcn. (3019->8), ass. (0->105)
t126 = 2 * pkin(6);
t93 = sin(qJ(4));
t94 = sin(qJ(3));
t97 = cos(qJ(4));
t98 = cos(qJ(3));
t65 = t93 * t98 + t94 * t97;
t95 = sin(qJ(2));
t55 = t65 * t95;
t64 = -t93 * t94 + t97 * t98;
t56 = t64 * t95;
t125 = -Ifges(5,5) * t56 + Ifges(5,6) * t55;
t124 = -pkin(8) - pkin(7);
t123 = pkin(3) * t93;
t99 = cos(qJ(2));
t122 = pkin(6) * t99;
t87 = t95 * pkin(6);
t121 = Ifges(4,4) * t94;
t120 = Ifges(4,4) * t98;
t81 = pkin(3) * t97 + pkin(4);
t92 = sin(qJ(5));
t96 = cos(qJ(5));
t58 = t96 * t123 + t81 * t92;
t119 = t58 * mrSges(6,2);
t118 = t92 * mrSges(6,2);
t117 = t94 * t95;
t116 = t95 * t98;
t115 = -Ifges(6,3) - Ifges(5,3);
t24 = -t55 * t96 - t56 * t92;
t25 = -t55 * t92 + t56 * t96;
t114 = -Ifges(6,5) * t25 - Ifges(6,6) * t24;
t72 = -pkin(2) * t99 - pkin(7) * t95 - pkin(1);
t63 = t98 * t72;
t35 = -pkin(8) * t116 + t63 + (-pkin(6) * t94 - pkin(3)) * t99;
t48 = t98 * t122 + t94 * t72;
t41 = -pkin(8) * t117 + t48;
t17 = t93 * t35 + t97 * t41;
t76 = t124 * t94;
t77 = t124 * t98;
t44 = t93 * t76 - t97 * t77;
t71 = pkin(3) * t117 + t87;
t113 = t94 ^ 2 + t98 ^ 2;
t112 = pkin(4) * t118;
t111 = -Ifges(4,3) + t115;
t82 = -pkin(3) * t98 - pkin(2);
t16 = t97 * t35 - t41 * t93;
t43 = t97 * t76 + t77 * t93;
t57 = -t92 * t123 + t81 * t96;
t54 = t57 * mrSges(6,1);
t110 = Ifges(6,3) + t54 - t119;
t109 = mrSges(4,1) * t94 + mrSges(4,2) * t98;
t11 = -pkin(9) * t55 + t17;
t6 = -pkin(4) * t99 - pkin(9) * t56 + t16;
t2 = -t11 * t92 + t6 * t96;
t3 = t11 * t96 + t6 * t92;
t108 = t2 * mrSges(6,1) - t3 * mrSges(6,2) - t114;
t27 = -pkin(9) * t65 + t43;
t28 = pkin(9) * t64 + t44;
t10 = t27 * t92 + t28 * t96;
t33 = t64 * t96 - t65 * t92;
t31 = Ifges(6,6) * t33;
t34 = t64 * t92 + t65 * t96;
t32 = Ifges(6,5) * t34;
t9 = t27 * t96 - t28 * t92;
t107 = t9 * mrSges(6,1) - t10 * mrSges(6,2) + t31 + t32;
t106 = (mrSges(5,1) * t97 - mrSges(5,2) * t93) * pkin(3);
t105 = t16 * mrSges(5,1) - t17 * mrSges(5,2) + t108 - t125;
t60 = Ifges(5,6) * t64;
t61 = Ifges(5,5) * t65;
t104 = t43 * mrSges(5,1) - t44 * mrSges(5,2) + t107 + t60 + t61;
t101 = pkin(6) ^ 2;
t91 = t99 ^ 2;
t89 = t95 ^ 2;
t86 = t89 * t101;
t85 = Ifges(4,5) * t94;
t84 = Ifges(4,6) * t98;
t83 = t96 * pkin(4) * mrSges(6,1);
t78 = Ifges(4,5) * t116;
t75 = Ifges(4,1) * t94 + t120;
t74 = Ifges(4,2) * t98 + t121;
t73 = -mrSges(4,1) * t98 + mrSges(4,2) * t94;
t70 = -mrSges(4,1) * t99 - mrSges(4,3) * t116;
t69 = mrSges(4,2) * t99 - mrSges(4,3) * t117;
t59 = t109 * t95;
t53 = -Ifges(4,5) * t99 + (Ifges(4,1) * t98 - t121) * t95;
t52 = -Ifges(4,6) * t99 + (-Ifges(4,2) * t94 + t120) * t95;
t49 = -pkin(4) * t64 + t82;
t47 = -t94 * t122 + t63;
t46 = -mrSges(5,1) * t99 - mrSges(5,3) * t56;
t45 = mrSges(5,2) * t99 - mrSges(5,3) * t55;
t40 = Ifges(5,1) * t65 + Ifges(5,4) * t64;
t39 = Ifges(5,4) * t65 + Ifges(5,2) * t64;
t38 = -mrSges(5,1) * t64 + mrSges(5,2) * t65;
t36 = pkin(4) * t55 + t71;
t26 = mrSges(5,1) * t55 + mrSges(5,2) * t56;
t23 = Ifges(5,1) * t56 - Ifges(5,4) * t55 - Ifges(5,5) * t99;
t22 = Ifges(5,4) * t56 - Ifges(5,2) * t55 - Ifges(5,6) * t99;
t19 = -mrSges(6,1) * t99 - mrSges(6,3) * t25;
t18 = mrSges(6,2) * t99 + mrSges(6,3) * t24;
t14 = Ifges(6,1) * t34 + Ifges(6,4) * t33;
t13 = Ifges(6,4) * t34 + Ifges(6,2) * t33;
t12 = -mrSges(6,1) * t33 + mrSges(6,2) * t34;
t7 = -mrSges(6,1) * t24 + mrSges(6,2) * t25;
t5 = Ifges(6,1) * t25 + Ifges(6,4) * t24 - Ifges(6,5) * t99;
t4 = Ifges(6,4) * t25 + Ifges(6,2) * t24 - Ifges(6,6) * t99;
t1 = [0.2e1 * t16 * t46 + 0.2e1 * t17 * t45 + 0.2e1 * t3 * t18 + 0.2e1 * t2 * t19 - t55 * t22 + t56 * t23 + t24 * t4 + t25 * t5 + 0.2e1 * t71 * t26 + 0.2e1 * t36 * t7 + 0.2e1 * t47 * t70 + 0.2e1 * t48 * t69 + Ifges(2,3) + (t89 + t91) * mrSges(3,3) * t126 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t95 + t59 * t126 - t52 * t94 + t53 * t98) * t95 + m(3) * (pkin(1) ^ 2 + t101 * t91 + t86) + m(4) * (t47 ^ 2 + t48 ^ 2 + t86) + m(5) * (t16 ^ 2 + t17 ^ 2 + t71 ^ 2) + m(6) * (t2 ^ 2 + t3 ^ 2 + t36 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) - t78 + (Ifges(3,2) - t111) * t99 + (Ifges(4,6) * t94 + (2 * Ifges(3,4))) * t95 + t114 + t125) * t99; (pkin(7) * t69 + t48 * mrSges(4,3) + t52 / 0.2e1) * t98 + (-t47 * mrSges(4,3) - pkin(7) * t70 + t53 / 0.2e1) * t94 + m(5) * (t16 * t43 + t17 * t44 + t71 * t82) + m(6) * (t10 * t3 + t2 * t9 + t36 * t49) + (-pkin(6) * mrSges(3,2) - t85 / 0.2e1 - t84 / 0.2e1 - t61 / 0.2e1 - t60 / 0.2e1 - t32 / 0.2e1 - t31 / 0.2e1 + Ifges(3,6)) * t99 + (-t2 * t34 + t3 * t33) * mrSges(6,3) + (-t16 * t65 + t17 * t64) * mrSges(5,3) + m(4) * (-pkin(2) * t87 + (-t47 * t94 + t48 * t98) * pkin(7)) + t82 * t26 + t64 * t22 / 0.2e1 + t65 * t23 / 0.2e1 + t71 * t38 - t55 * t39 / 0.2e1 + t56 * t40 / 0.2e1 - pkin(2) * t59 + t44 * t45 + t43 * t46 + t49 * t7 + t25 * t14 / 0.2e1 + t33 * t4 / 0.2e1 + t34 * t5 / 0.2e1 + t36 * t12 + t10 * t18 + t9 * t19 + t24 * t13 / 0.2e1 + (t98 * t75 / 0.2e1 - t94 * t74 / 0.2e1 + Ifges(3,5) + (t73 - mrSges(3,1)) * pkin(6)) * t95; -0.2e1 * pkin(2) * t73 + 0.2e1 * t49 * t12 + t33 * t13 + t34 * t14 + 0.2e1 * t82 * t38 + t64 * t39 + t65 * t40 + t98 * t74 + t94 * t75 + Ifges(3,3) + m(6) * (t10 ^ 2 + t49 ^ 2 + t9 ^ 2) + m(5) * (t43 ^ 2 + t44 ^ 2 + t82 ^ 2) + m(4) * (t113 * pkin(7) ^ 2 + pkin(2) ^ 2) + 0.2e1 * (t10 * t33 - t34 * t9) * mrSges(6,3) + 0.2e1 * (-t43 * t65 + t44 * t64) * mrSges(5,3) + 0.2e1 * t113 * pkin(7) * mrSges(4,3); m(6) * (t2 * t57 + t3 * t58) + t111 * t99 + t78 - Ifges(4,6) * t117 + t57 * t19 + t58 * t18 + t47 * mrSges(4,1) - t48 * mrSges(4,2) + t105 + (t97 * t46 + t93 * t45 + m(5) * (t16 * t97 + t17 * t93)) * pkin(3); m(6) * (t10 * t58 + t57 * t9) + t85 + t84 - t109 * pkin(7) + (t33 * t58 - t34 * t57) * mrSges(6,3) + (m(5) * (t43 * t97 + t44 * t93) + (t64 * t93 - t65 * t97) * mrSges(5,3)) * pkin(3) + t104; -0.2e1 * t119 + 0.2e1 * t54 + 0.2e1 * t106 + m(6) * (t57 ^ 2 + t58 ^ 2) + m(5) * (t93 ^ 2 + t97 ^ 2) * pkin(3) ^ 2 - t111; t115 * t99 + (m(6) * (t2 * t96 + t3 * t92) + t92 * t18 + t96 * t19) * pkin(4) + t105; (m(6) * (t10 * t92 + t9 * t96) + (t33 * t92 - t34 * t96) * mrSges(6,3)) * pkin(4) + t104; Ifges(5,3) + t83 + t106 + (m(6) * (t57 * t96 + t58 * t92) - t118) * pkin(4) + t110; -0.2e1 * t112 + 0.2e1 * t83 + m(6) * (t92 ^ 2 + t96 ^ 2) * pkin(4) ^ 2 - t115; -Ifges(6,3) * t99 + t108; t107; t110; Ifges(6,3) + t83 - t112; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
