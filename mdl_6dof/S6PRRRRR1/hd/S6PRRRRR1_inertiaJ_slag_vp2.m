% Calculate joint inertia matrix for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:37:38
% EndTime: 2019-03-09 00:37:40
% DurationCPUTime: 0.99s
% Computational Cost: add. (1979->243), mult. (4072->366), div. (0->0), fcn. (4590->12), ass. (0->116)
t93 = sin(qJ(6));
t98 = cos(qJ(6));
t69 = -mrSges(7,1) * t98 + mrSges(7,2) * t93;
t156 = t69 - mrSges(6,1);
t132 = t93 * mrSges(7,3);
t100 = cos(qJ(4));
t101 = cos(qJ(3));
t95 = sin(qJ(4));
t96 = sin(qJ(3));
t63 = t100 * t101 - t95 * t96;
t64 = t100 * t96 + t101 * t95;
t94 = sin(qJ(5));
t99 = cos(qJ(5));
t39 = -t99 * t63 + t64 * t94;
t40 = t63 * t94 + t64 * t99;
t23 = -mrSges(7,2) * t39 - t132 * t40;
t135 = t40 * t98;
t24 = mrSges(7,1) * t39 - mrSges(7,3) * t135;
t155 = t98 * t23 - t93 * t24;
t144 = pkin(5) * t69;
t87 = t93 ^ 2;
t142 = mrSges(7,3) * t87;
t81 = pkin(11) * t142;
t89 = t98 ^ 2;
t141 = mrSges(7,3) * t89;
t82 = pkin(11) * t141;
t154 = t81 + t82 - t144;
t147 = -pkin(9) - pkin(8);
t120 = t147 * t101;
t121 = t147 * t96;
t44 = t100 * t121 + t120 * t95;
t109 = -t64 * pkin(10) + t44;
t45 = -t100 * t120 + t95 * t121;
t32 = pkin(10) * t63 + t45;
t16 = -t109 * t99 + t32 * t94;
t115 = mrSges(7,1) * t93 + mrSges(7,2) * t98;
t22 = t115 * t40;
t153 = m(7) * t16 + t22;
t145 = pkin(4) * t99;
t78 = -pkin(5) - t145;
t58 = t78 * t69;
t77 = pkin(4) * t94 + pkin(11);
t67 = t77 * t142;
t68 = t77 * t141;
t83 = mrSges(6,1) * t145;
t152 = t58 + t67 + t68 + t83;
t80 = -pkin(3) * t101 - pkin(2);
t48 = -pkin(4) * t63 + t80;
t14 = pkin(5) * t39 - pkin(11) * t40 + t48;
t18 = t109 * t94 + t99 * t32;
t3 = t14 * t93 + t18 * t98;
t143 = t3 * t98;
t2 = t14 * t98 - t18 * t93;
t117 = -t2 * t93 + t143;
t151 = m(7) * t117 + t155;
t150 = t16 ^ 2;
t91 = sin(pkin(6));
t97 = sin(qJ(2));
t133 = t91 * t97;
t92 = cos(pkin(6));
t55 = t101 * t92 - t133 * t96;
t56 = t101 * t133 + t92 * t96;
t33 = t100 * t55 - t56 * t95;
t34 = t100 * t56 + t55 * t95;
t19 = -t99 * t33 + t34 * t94;
t149 = t19 ^ 2;
t88 = t96 ^ 2;
t148 = 0.2e1 * t16;
t146 = pkin(3) * t95;
t140 = Ifges(7,4) * t93;
t139 = Ifges(7,4) * t98;
t102 = cos(qJ(2));
t124 = t102 * t91;
t21 = t33 * t94 + t34 * t99;
t10 = -t124 * t93 + t21 * t98;
t138 = t10 * t98;
t137 = t16 * t19;
t136 = t40 * t93;
t79 = pkin(3) * t100 + pkin(4);
t54 = t99 * t146 + t94 * t79;
t134 = t54 * mrSges(6,2);
t130 = t94 * mrSges(6,2);
t128 = Ifges(7,5) * t135 + Ifges(7,3) * t39;
t127 = Ifges(7,5) * t93 + Ifges(7,6) * t98;
t126 = t87 + t89;
t125 = t101 ^ 2 + t88;
t123 = pkin(4) * t130;
t71 = Ifges(7,2) * t98 + t140;
t72 = Ifges(7,1) * t93 + t139;
t122 = t98 * t71 + t93 * t72 + Ifges(6,3);
t119 = t126 * t77;
t118 = Ifges(5,3) + t122;
t9 = -t124 * t98 - t21 * t93;
t116 = -t9 * t93 + t138;
t114 = t101 * t56 - t55 * t96;
t53 = -t146 * t94 + t79 * t99;
t113 = (mrSges(5,1) * t100 - mrSges(5,2) * t95) * pkin(3);
t112 = -t21 * mrSges(6,2) + mrSges(7,3) * t138 - t132 * t9 + t156 * t19;
t51 = -pkin(5) - t53;
t42 = t51 * t69;
t52 = pkin(11) + t54;
t46 = t52 * t142;
t47 = t52 * t141;
t49 = t53 * mrSges(6,1);
t111 = t122 + t42 + t46 + t47 + t49 - t134;
t110 = t33 * mrSges(5,1) - t34 * mrSges(5,2) + t112;
t6 = Ifges(7,6) * t39 + (-Ifges(7,2) * t93 + t139) * t40;
t7 = Ifges(7,5) * t39 + (Ifges(7,1) * t98 - t140) * t40;
t108 = -t18 * mrSges(6,2) + mrSges(7,3) * t143 - t132 * t2 - t71 * t136 / 0.2e1 + t72 * t135 / 0.2e1 + Ifges(6,5) * t40 + t93 * t7 / 0.2e1 + t98 * t6 / 0.2e1 + (t127 / 0.2e1 - Ifges(6,6)) * t39 + t156 * t16;
t107 = t44 * mrSges(5,1) - t45 * mrSges(5,2) + Ifges(5,5) * t64 + Ifges(5,6) * t63 + t108;
t86 = t91 ^ 2;
t76 = t86 * t102 ^ 2;
t70 = -mrSges(4,1) * t101 + mrSges(4,2) * t96;
t41 = -mrSges(5,1) * t63 + mrSges(5,2) * t64;
t25 = mrSges(6,1) * t39 + mrSges(6,2) * t40;
t1 = [m(2) + m(7) * (t10 ^ 2 + t9 ^ 2 + t149) + m(6) * (t21 ^ 2 + t149 + t76) + m(5) * (t33 ^ 2 + t34 ^ 2 + t76) + m(4) * (t55 ^ 2 + t56 ^ 2 + t76) + m(3) * (t86 * t97 ^ 2 + t92 ^ 2 + t76); t10 * t23 + t19 * t22 + t9 * t24 + (t19 * t40 - t21 * t39) * mrSges(6,3) + (-t33 * t64 + t34 * t63) * mrSges(5,3) + t114 * mrSges(4,3) + (-t97 * mrSges(3,2) + (mrSges(3,1) - t25 - t41 - t70) * t102) * t91 + m(7) * (t10 * t3 + t2 * t9 + t137) + m(6) * (-t124 * t48 + t18 * t21 + t137) + m(5) * (-t124 * t80 + t33 * t44 + t34 * t45) + m(4) * (pkin(2) * t124 + pkin(8) * t114); Ifges(4,1) * t88 - 0.2e1 * pkin(2) * t70 + t22 * t148 + 0.2e1 * t2 * t24 + 0.2e1 * t3 * t23 + 0.2e1 * t48 * t25 + 0.2e1 * t80 * t41 + Ifges(3,3) + (0.2e1 * Ifges(4,4) * t96 + Ifges(4,2) * t101) * t101 + (-0.2e1 * t44 * mrSges(5,3) + Ifges(5,1) * t64) * t64 + 0.2e1 * t125 * pkin(8) * mrSges(4,3) + (0.2e1 * t45 * mrSges(5,3) + 0.2e1 * Ifges(5,4) * t64 + Ifges(5,2) * t63) * t63 + (-0.2e1 * t18 * mrSges(6,3) + Ifges(6,2) * t39 + t128) * t39 + (mrSges(6,3) * t148 + Ifges(6,1) * t40 - t93 * t6 + t98 * t7 + (-Ifges(7,6) * t93 - (2 * Ifges(6,4))) * t39) * t40 + m(4) * (pkin(8) ^ 2 * t125 + pkin(2) ^ 2) + m(5) * (t44 ^ 2 + t45 ^ 2 + t80 ^ 2) + m(6) * (t18 ^ 2 + t48 ^ 2 + t150) + m(7) * (t2 ^ 2 + t3 ^ 2 + t150); t55 * mrSges(4,1) - t56 * mrSges(4,2) + m(7) * (t116 * t52 + t19 * t51) + m(6) * (-t19 * t53 + t21 * t54) + m(5) * (t100 * t33 + t34 * t95) * pkin(3) + t110; t107 + (m(5) * (t100 * t44 + t45 * t95) + (-t100 * t64 + t95 * t63) * mrSges(5,3)) * pkin(3) + m(7) * (t117 * t52 + t16 * t51) + t155 * t52 + m(6) * (-t16 * t53 + t18 * t54) + Ifges(4,6) * t101 + Ifges(4,5) * t96 + t51 * t22 + (-t54 * t39 - t53 * t40) * mrSges(6,3) + (-t96 * mrSges(4,1) - t101 * mrSges(4,2)) * pkin(8); -0.2e1 * t134 + Ifges(4,3) + 0.2e1 * t42 + 0.2e1 * t46 + 0.2e1 * t47 + 0.2e1 * t49 + 0.2e1 * t113 + m(7) * (t126 * t52 ^ 2 + t51 ^ 2) + m(6) * (t53 ^ 2 + t54 ^ 2) + m(5) * (t100 ^ 2 + t95 ^ 2) * pkin(3) ^ 2 + t118; m(7) * (t116 * t77 + t78 * t19) + m(6) * (-t19 * t99 + t21 * t94) * pkin(4) + t110; t107 + (m(6) * (-t16 * t99 + t18 * t94) + (-t39 * t94 - t40 * t99) * mrSges(6,3)) * pkin(4) + t153 * t78 + t151 * t77; t111 + m(7) * (t119 * t52 + t51 * t78) + (m(6) * (t53 * t99 + t54 * t94) - t130) * pkin(4) + t113 + Ifges(5,3) + t152; -0.2e1 * t123 + 0.2e1 * t58 + 0.2e1 * t67 + 0.2e1 * t68 + 0.2e1 * t83 + m(7) * (t126 * t77 ^ 2 + t78 ^ 2) + m(6) * (t94 ^ 2 + t99 ^ 2) * pkin(4) ^ 2 + t118; m(7) * (-pkin(5) * t19 + pkin(11) * t116) + t112; -t153 * pkin(5) + t151 * pkin(11) + t108; m(7) * (pkin(11) * t126 * t52 - pkin(5) * t51) + t111 + t154; m(7) * (-pkin(5) * t78 + pkin(11) * t119) - t123 + t122 + t152 + t154; -0.2e1 * t144 + m(7) * (pkin(11) ^ 2 * t126 + pkin(5) ^ 2) + 0.2e1 * t81 + 0.2e1 * t82 + t122; mrSges(7,1) * t9 - mrSges(7,2) * t10; mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t136 + t128; -t115 * t52 + t127; -t115 * t77 + t127; -pkin(11) * t115 + t127; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
