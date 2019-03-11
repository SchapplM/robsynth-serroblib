% Calculate joint inertia matrix for
% S6RRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:31:52
% EndTime: 2019-03-09 15:31:56
% DurationCPUTime: 1.29s
% Computational Cost: add. (1862->322), mult. (3675->442), div. (0->0), fcn. (3717->8), ass. (0->113)
t148 = 2 * pkin(7);
t119 = sin(qJ(2));
t115 = sin(pkin(10));
t116 = cos(pkin(10));
t118 = sin(qJ(3));
t121 = cos(qJ(3));
t84 = t115 * t121 + t116 * t118;
t74 = t84 * t119;
t83 = t115 * t118 - t116 * t121;
t75 = t83 * t119;
t92 = (pkin(3) * t118 + pkin(7)) * t119;
t19 = t74 * pkin(4) + t75 * qJ(5) + t92;
t137 = t119 * t121;
t147 = -Ifges(4,5) * t137 + (Ifges(6,4) + Ifges(5,5)) * t75 + (Ifges(5,6) - Ifges(6,6)) * t74;
t122 = cos(qJ(2));
t146 = pkin(7) * t122;
t101 = pkin(3) * t115 + qJ(5);
t117 = sin(qJ(6));
t120 = cos(qJ(6));
t103 = -pkin(3) * t116 - pkin(4);
t98 = -pkin(5) + t103;
t61 = -t101 * t117 + t120 * t98;
t145 = t61 * mrSges(7,1);
t62 = t101 * t120 + t117 * t98;
t144 = t62 * mrSges(7,2);
t142 = -qJ(4) - pkin(8);
t93 = -pkin(2) * t122 - pkin(8) * t119 - pkin(1);
t86 = t121 * t93;
t41 = -qJ(4) * t137 + t86 + (-pkin(7) * t118 - pkin(3)) * t122;
t138 = t118 * t119;
t60 = t118 * t93 + t121 * t146;
t50 = -qJ(4) * t138 + t60;
t18 = t115 * t41 + t116 * t50;
t131 = t142 * t118;
t95 = t142 * t121;
t54 = t115 * t131 - t116 * t95;
t140 = Ifges(4,4) * t118;
t139 = Ifges(4,4) * t121;
t58 = t122 * mrSges(6,1) - t75 * mrSges(6,2);
t136 = t118 ^ 2 + t121 ^ 2;
t52 = -t115 * t95 - t116 * t131;
t135 = t52 ^ 2 + t54 ^ 2;
t134 = Ifges(5,3) + Ifges(4,3) + Ifges(6,2);
t29 = t117 * t75 + t120 * t74;
t30 = t117 * t74 - t120 * t75;
t133 = Ifges(7,5) * t30 + Ifges(7,6) * t29 + Ifges(7,3) * t122;
t104 = -pkin(3) * t121 - pkin(2);
t34 = t74 * mrSges(5,1) - t75 * mrSges(5,2);
t45 = t83 * mrSges(5,1) + t84 * mrSges(5,2);
t33 = t74 * mrSges(6,1) + t75 * mrSges(6,3);
t44 = t83 * mrSges(6,1) - t84 * mrSges(6,3);
t17 = -t115 * t50 + t116 * t41;
t10 = -qJ(5) * t122 + t18;
t13 = t122 * pkin(4) - t17;
t6 = -t29 * mrSges(7,1) + t30 * mrSges(7,2);
t39 = -t117 * t84 + t120 * t83;
t40 = t117 * t83 + t120 * t84;
t14 = -t39 * mrSges(7,1) + t40 * mrSges(7,2);
t130 = mrSges(4,1) * t118 + mrSges(4,2) * t121;
t129 = t120 * mrSges(7,1) - t117 * mrSges(7,2);
t128 = qJ(5) * t84 - t104;
t37 = Ifges(7,6) * t39;
t38 = Ifges(7,5) * t40;
t31 = -pkin(9) * t84 + t52;
t32 = pkin(9) * t83 + t54;
t8 = -t117 * t32 + t120 * t31;
t9 = t117 * t31 + t120 * t32;
t127 = t8 * mrSges(7,1) - t9 * mrSges(7,2) + t37 + t38;
t5 = pkin(5) * t122 + pkin(9) * t75 + t13;
t7 = pkin(9) * t74 + t10;
t1 = -t117 * t7 + t120 * t5;
t2 = t117 * t5 + t120 * t7;
t126 = t1 * mrSges(7,1) - t2 * mrSges(7,2) + t133;
t124 = pkin(7) ^ 2;
t114 = t122 ^ 2;
t112 = t119 ^ 2;
t109 = t112 * t124;
t108 = Ifges(4,5) * t118;
t107 = Ifges(4,6) * t121;
t97 = Ifges(4,1) * t118 + t139;
t96 = Ifges(4,2) * t121 + t140;
t94 = -mrSges(4,1) * t121 + mrSges(4,2) * t118;
t91 = -mrSges(4,1) * t122 - mrSges(4,3) * t137;
t90 = mrSges(4,2) * t122 - mrSges(4,3) * t138;
t82 = t130 * t119;
t81 = Ifges(6,4) * t84;
t80 = Ifges(5,5) * t84;
t79 = Ifges(5,6) * t83;
t78 = Ifges(6,6) * t83;
t73 = -Ifges(4,5) * t122 + (Ifges(4,1) * t121 - t140) * t119;
t72 = -Ifges(4,6) * t122 + (-Ifges(4,2) * t118 + t139) * t119;
t59 = -t118 * t146 + t86;
t57 = -mrSges(5,1) * t122 + mrSges(5,3) * t75;
t56 = mrSges(5,2) * t122 - mrSges(5,3) * t74;
t55 = -mrSges(6,2) * t74 - mrSges(6,3) * t122;
t49 = Ifges(5,1) * t84 - Ifges(5,4) * t83;
t48 = Ifges(6,1) * t84 + Ifges(6,5) * t83;
t47 = Ifges(5,4) * t84 - Ifges(5,2) * t83;
t46 = Ifges(6,5) * t84 + Ifges(6,3) * t83;
t35 = pkin(4) * t83 - t128;
t28 = -Ifges(5,1) * t75 - Ifges(5,4) * t74 - Ifges(5,5) * t122;
t27 = -Ifges(6,1) * t75 - Ifges(6,4) * t122 + Ifges(6,5) * t74;
t26 = -Ifges(5,4) * t75 - Ifges(5,2) * t74 - Ifges(5,6) * t122;
t25 = -Ifges(6,5) * t75 - Ifges(6,6) * t122 + Ifges(6,3) * t74;
t22 = (-pkin(4) - pkin(5)) * t83 + t128;
t21 = mrSges(7,1) * t122 - mrSges(7,3) * t30;
t20 = -mrSges(7,2) * t122 + mrSges(7,3) * t29;
t16 = Ifges(7,1) * t40 + Ifges(7,4) * t39;
t15 = Ifges(7,4) * t40 + Ifges(7,2) * t39;
t11 = pkin(5) * t74 + t19;
t4 = Ifges(7,1) * t30 + Ifges(7,4) * t29 + Ifges(7,5) * t122;
t3 = Ifges(7,4) * t30 + Ifges(7,2) * t29 + Ifges(7,6) * t122;
t12 = [0.2e1 * t1 * t21 + 0.2e1 * t10 * t55 - 0.2e1 * t11 * t6 + 0.2e1 * t13 * t58 + 0.2e1 * t17 * t57 + 0.2e1 * t18 * t56 + 0.2e1 * t19 * t33 + 0.2e1 * t2 * t20 + t29 * t3 + t30 * t4 + 0.2e1 * t92 * t34 + 0.2e1 * t59 * t91 + 0.2e1 * t60 * t90 + Ifges(2,3) - (t27 + t28) * t75 + (t25 - t26) * t74 + (t112 + t114) * mrSges(3,3) * t148 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t119 - t118 * t72 + t121 * t73 + t82 * t148) * t119 + m(3) * (pkin(1) ^ 2 + t114 * t124 + t109) + m(4) * (t59 ^ 2 + t60 ^ 2 + t109) + m(5) * (t17 ^ 2 + t18 ^ 2 + t92 ^ 2) + m(7) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2) + m(6) * (t10 ^ 2 + t13 ^ 2 + t19 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + t134) * t122 + (Ifges(4,6) * t118 + (2 * Ifges(3,4))) * t119 + t133 + t147) * t122; m(5) * (t104 * t92 - t17 * t52 + t18 * t54) + m(6) * (t10 * t54 + t13 * t52 + t19 * t35) + m(7) * (t1 * t8 - t11 * t22 + t2 * t9) + (-t108 / 0.2e1 - t107 / 0.2e1 - t80 / 0.2e1 + t79 / 0.2e1 - t81 / 0.2e1 - t78 / 0.2e1 + t38 / 0.2e1 + t37 / 0.2e1 + Ifges(3,6) - pkin(7) * mrSges(3,2)) * t122 + m(4) * (-pkin(2) * pkin(7) * t119 + (-t59 * t118 + t60 * t121) * pkin(8)) + (t27 / 0.2e1 + t28 / 0.2e1 + t13 * mrSges(6,2) - t17 * mrSges(5,3)) * t84 + (t25 / 0.2e1 - t26 / 0.2e1 - t10 * mrSges(6,2) - t18 * mrSges(5,3)) * t83 - (t48 / 0.2e1 + t49 / 0.2e1) * t75 + (t46 / 0.2e1 - t47 / 0.2e1) * t74 + (Ifges(3,5) + t121 * t97 / 0.2e1 - t118 * t96 / 0.2e1 + (t94 - mrSges(3,1)) * pkin(7)) * t119 + (-t1 * t40 + t2 * t39) * mrSges(7,3) + t104 * t34 + t92 * t45 - pkin(2) * t82 + t39 * t3 / 0.2e1 + t40 * t4 / 0.2e1 + t19 * t44 + t29 * t15 / 0.2e1 + t30 * t16 / 0.2e1 + t35 * t33 + t9 * t20 + t8 * t21 + t22 * t6 - t11 * t14 + (t72 / 0.2e1 + pkin(8) * t90 + t60 * mrSges(4,3)) * t121 + (t73 / 0.2e1 - pkin(8) * t91 - t59 * mrSges(4,3)) * t118 + (t55 + t56) * t54 + (-t57 + t58) * t52; -0.2e1 * pkin(2) * t94 + 0.2e1 * t104 * t45 + t118 * t97 + t121 * t96 + 0.2e1 * t22 * t14 + t39 * t15 + t40 * t16 + 0.2e1 * t35 * t44 + Ifges(3,3) + (t48 + t49) * t84 + (t46 - t47) * t83 + m(4) * (t136 * pkin(8) ^ 2 + pkin(2) ^ 2) + m(5) * (t104 ^ 2 + t135) + m(6) * (t35 ^ 2 + t135) + m(7) * (t22 ^ 2 + t8 ^ 2 + t9 ^ 2) + 0.2e1 * (t39 * t9 - t40 * t8) * mrSges(7,3) + 0.2e1 * t136 * pkin(8) * mrSges(4,3) + 0.2e1 * (t52 * t84 - t54 * t83) * (mrSges(6,2) + mrSges(5,3)); (t115 * t56 + m(5) * (t115 * t18 + t116 * t17) + t116 * t57) * pkin(3) + m(6) * (t10 * t101 + t103 * t13) + m(7) * (t1 * t61 + t2 * t62) - t147 + t101 * t55 + t103 * t58 + t59 * mrSges(4,1) - t60 * mrSges(4,2) + t61 * t21 + t62 * t20 + t10 * mrSges(6,3) - t13 * mrSges(6,1) + t17 * mrSges(5,1) - t18 * mrSges(5,2) - t126 - t134 * t122 - Ifges(4,6) * t138; t107 + t108 + t78 - t79 + t80 + t81 + (-mrSges(5,2) + mrSges(6,3)) * t54 + (-mrSges(6,1) - mrSges(5,1)) * t52 - t130 * pkin(8) + (t39 * t62 - t40 * t61) * mrSges(7,3) + (-t101 * t83 + t103 * t84) * mrSges(6,2) + m(6) * (t101 * t54 + t103 * t52) + m(7) * (t61 * t8 + t62 * t9) + (m(5) * (t115 * t54 - t116 * t52) + (-t115 * t83 - t116 * t84) * mrSges(5,3)) * pkin(3) - t127; -0.2e1 * t103 * mrSges(6,1) - 0.2e1 * t145 + 0.2e1 * t144 + 0.2e1 * t101 * mrSges(6,3) + Ifges(7,3) + m(7) * (t61 ^ 2 + t62 ^ 2) + m(6) * (t101 ^ 2 + t103 ^ 2) + t134 + (0.2e1 * mrSges(5,1) * t116 - 0.2e1 * mrSges(5,2) * t115 + m(5) * (t115 ^ 2 + t116 ^ 2) * pkin(3)) * pkin(3); m(5) * t92 + m(6) * t19 + m(7) * t11 + t33 + t34 - t6; m(5) * t104 + m(6) * t35 - m(7) * t22 - t14 + t44 + t45; 0; m(5) + m(6) + m(7); t117 * t20 + t120 * t21 + m(7) * (t1 * t120 + t117 * t2) + m(6) * t13 + t58; t84 * mrSges(6,2) + (t117 * t39 - t120 * t40) * mrSges(7,3) + m(7) * (t117 * t9 + t120 * t8) + m(6) * t52; -mrSges(6,1) + m(7) * (t117 * t62 + t120 * t61) + m(6) * t103 - t129; 0; m(6) + m(7) * (t117 ^ 2 + t120 ^ 2); t126; t127; -Ifges(7,3) - t144 + t145; 0; t129; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
