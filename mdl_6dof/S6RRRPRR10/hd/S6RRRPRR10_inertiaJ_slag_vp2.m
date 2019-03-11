% Calculate joint inertia matrix for
% S6RRRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR10_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR10_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:15:40
% EndTime: 2019-03-09 19:15:45
% DurationCPUTime: 1.45s
% Computational Cost: add. (2060->351), mult. (3854->490), div. (0->0), fcn. (3816->8), ass. (0->132)
t121 = sin(qJ(3));
t125 = cos(qJ(3));
t165 = t121 ^ 2 + t125 ^ 2;
t164 = 2 * pkin(7);
t122 = sin(qJ(2));
t147 = t122 * t125;
t148 = t121 * t122;
t163 = -Ifges(5,6) * t148 + (-Ifges(5,4) - Ifges(4,5)) * t147;
t126 = cos(qJ(2));
t120 = sin(qJ(5));
t124 = cos(qJ(5));
t73 = -t120 * t125 + t121 * t124;
t61 = t73 * t122;
t136 = t120 * t121 + t124 * t125;
t62 = t136 * t122;
t162 = Ifges(6,5) * t62 + Ifges(6,6) * t61 + Ifges(6,3) * t126;
t127 = -pkin(3) - pkin(4);
t161 = pkin(8) - pkin(9);
t160 = pkin(3) * t121;
t159 = pkin(7) * t126;
t119 = sin(qJ(6));
t123 = cos(qJ(6));
t82 = -qJ(4) * t120 + t124 * t127;
t81 = -pkin(5) + t82;
t83 = qJ(4) * t124 + t120 * t127;
t40 = t119 * t81 + t123 * t83;
t158 = t40 * mrSges(7,2);
t157 = t82 * mrSges(6,1);
t156 = t83 * mrSges(6,2);
t155 = Ifges(5,2) + Ifges(4,3);
t154 = Ifges(6,3) + Ifges(7,3);
t114 = t126 * pkin(3);
t85 = -pkin(2) * t126 - pkin(8) * t122 - pkin(1);
t99 = t121 * t159;
t31 = pkin(4) * t126 + t114 + t99 + (-pkin(9) * t122 - t85) * t125;
t53 = t121 * t85 + t125 * t159;
t48 = -qJ(4) * t126 + t53;
t37 = pkin(9) * t148 + t48;
t13 = t120 * t31 + t124 * t37;
t92 = t161 * t121;
t93 = t161 * t125;
t47 = t120 * t92 + t124 * t93;
t153 = Ifges(4,4) * t121;
t152 = Ifges(4,4) * t125;
t151 = Ifges(5,5) * t121;
t150 = Ifges(5,5) * t125;
t149 = t126 * Ifges(5,6);
t79 = t126 * mrSges(5,1) + mrSges(5,2) * t147;
t146 = t165 * pkin(8) ^ 2;
t84 = -t125 * pkin(3) - t121 * qJ(4) - pkin(2);
t24 = -t119 * t62 + t123 * t61;
t25 = t119 * t61 + t123 * t62;
t145 = Ifges(7,5) * t25 + Ifges(7,6) * t24 + Ifges(7,3) * t126;
t70 = -t119 * t120 + t123 * t124;
t72 = t119 * t124 + t120 * t123;
t144 = t70 * mrSges(7,1) - t72 * mrSges(7,2);
t12 = -t120 * t37 + t124 * t31;
t46 = -t120 * t93 + t124 * t92;
t52 = t125 * t85 - t99;
t65 = t125 * pkin(4) - t84;
t39 = -t119 * t83 + t123 * t81;
t38 = t39 * mrSges(7,1);
t143 = -Ifges(7,3) + t38 - t158;
t27 = -pkin(10) * t73 + t46;
t28 = -pkin(10) * t136 + t47;
t10 = -t119 * t28 + t123 * t27;
t11 = t119 * t27 + t123 * t28;
t35 = -t119 * t73 - t123 * t136;
t33 = Ifges(7,6) * t35;
t36 = -t119 * t136 + t123 * t73;
t34 = Ifges(7,5) * t36;
t141 = -t10 * mrSges(7,1) + t11 * mrSges(7,2) - t33 - t34;
t140 = t121 * mrSges(4,1) + t125 * mrSges(4,2);
t139 = t121 * mrSges(5,1) - t125 * mrSges(5,3);
t138 = t123 * mrSges(7,1) - t119 * mrSges(7,2);
t137 = qJ(4) * t125 - t160;
t4 = pkin(5) * t126 - pkin(10) * t62 + t12;
t5 = pkin(10) * t61 + t13;
t2 = -t119 * t5 + t123 * t4;
t3 = t119 * t4 + t123 * t5;
t135 = -t2 * mrSges(7,1) + t3 * mrSges(7,2) - t145;
t134 = t138 * pkin(5);
t94 = qJ(4) * t147;
t45 = t94 + (t127 * t121 - pkin(7)) * t122;
t133 = t124 * mrSges(6,1) - t120 * mrSges(6,2) + t144;
t67 = Ifges(6,6) * t136;
t68 = Ifges(6,5) * t73;
t132 = t46 * mrSges(6,1) - t47 * mrSges(6,2) - t141 - t67 + t68;
t131 = t12 * mrSges(6,1) - t13 * mrSges(6,2) - t135 + t162;
t129 = pkin(7) ^ 2;
t118 = t126 ^ 2;
t116 = t122 ^ 2;
t110 = t116 * t129;
t108 = Ifges(5,4) * t121;
t107 = Ifges(4,5) * t121;
t106 = Ifges(4,6) * t125;
t91 = Ifges(4,1) * t121 + t152;
t90 = Ifges(5,1) * t121 - t150;
t89 = Ifges(4,2) * t125 + t153;
t88 = -Ifges(5,3) * t125 + t151;
t87 = -mrSges(4,1) * t125 + mrSges(4,2) * t121;
t86 = -mrSges(5,1) * t125 - mrSges(5,3) * t121;
t80 = -mrSges(5,2) * t148 - mrSges(5,3) * t126;
t78 = -mrSges(4,1) * t126 - mrSges(4,3) * t147;
t77 = mrSges(4,2) * t126 - mrSges(4,3) * t148;
t64 = t140 * t122;
t63 = t139 * t122;
t60 = -t94 + (pkin(7) + t160) * t122;
t59 = -Ifges(4,5) * t126 + (Ifges(4,1) * t125 - t153) * t122;
t58 = -Ifges(5,4) * t126 + (Ifges(5,1) * t125 + t151) * t122;
t57 = -Ifges(4,6) * t126 + (-Ifges(4,2) * t121 + t152) * t122;
t56 = -t149 + (Ifges(5,3) * t121 + t150) * t122;
t51 = mrSges(6,1) * t126 - mrSges(6,3) * t62;
t50 = -mrSges(6,2) * t126 + mrSges(6,3) * t61;
t49 = t114 - t52;
t44 = pkin(5) * t136 + t65;
t43 = Ifges(6,1) * t73 - Ifges(6,4) * t136;
t42 = Ifges(6,4) * t73 - Ifges(6,2) * t136;
t41 = mrSges(6,1) * t136 + mrSges(6,2) * t73;
t26 = -mrSges(6,1) * t61 + mrSges(6,2) * t62;
t23 = Ifges(6,1) * t62 + Ifges(6,4) * t61 + Ifges(6,5) * t126;
t22 = Ifges(6,4) * t62 + Ifges(6,2) * t61 + Ifges(6,6) * t126;
t19 = -pkin(5) * t61 + t45;
t18 = mrSges(7,1) * t126 - mrSges(7,3) * t25;
t17 = -mrSges(7,2) * t126 + mrSges(7,3) * t24;
t16 = Ifges(7,1) * t36 + Ifges(7,4) * t35;
t15 = Ifges(7,4) * t36 + Ifges(7,2) * t35;
t14 = -mrSges(7,1) * t35 + mrSges(7,2) * t36;
t8 = -mrSges(7,1) * t24 + mrSges(7,2) * t25;
t7 = Ifges(7,1) * t25 + Ifges(7,4) * t24 + Ifges(7,5) * t126;
t6 = Ifges(7,4) * t25 + Ifges(7,2) * t24 + Ifges(7,6) * t126;
t1 = [0.2e1 * t12 * t51 + 0.2e1 * t13 * t50 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 + 0.2e1 * t19 * t8 + t61 * t22 + t62 * t23 + t24 * t6 + t25 * t7 + 0.2e1 * t45 * t26 + 0.2e1 * t48 * t80 + 0.2e1 * t49 * t79 + 0.2e1 * t52 * t78 + 0.2e1 * t53 * t77 + 0.2e1 * t60 * t63 + Ifges(2,3) + (t116 + t118) * mrSges(3,3) * t164 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t122 + t64 * t164 + (t58 + t59) * t125 + (t56 - t57) * t121) * t122 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + t155) * t126 + (Ifges(4,6) * t121 + (2 * Ifges(3,4))) * t122 + t145 + t162 + t163) * t126 + m(3) * (pkin(1) ^ 2 + t118 * t129 + t110) + m(4) * (t52 ^ 2 + t53 ^ 2 + t110) + m(5) * (t48 ^ 2 + t49 ^ 2 + t60 ^ 2) + m(6) * (t12 ^ 2 + t13 ^ 2 + t45 ^ 2) + m(7) * (t19 ^ 2 + t2 ^ 2 + t3 ^ 2); m(6) * (t12 * t46 + t13 * t47 + t45 * t65) + m(7) * (t10 * t2 + t11 * t3 + t19 * t44) + (-t12 * t73 - t13 * t136) * mrSges(6,3) - t136 * t22 / 0.2e1 + (-t108 / 0.2e1 - t107 / 0.2e1 - t106 / 0.2e1 + t68 / 0.2e1 - t67 / 0.2e1 + t34 / 0.2e1 + t33 / 0.2e1 + Ifges(3,6)) * t126 + m(4) * (-pkin(2) * pkin(7) * t122 + (-t52 * t121 + t53 * t125) * pkin(8)) + m(5) * (t60 * t84 + (t49 * t121 + t48 * t125) * pkin(8)) + (-t2 * t36 + t3 * t35) * mrSges(7,3) + Ifges(3,5) * t122 + t84 * t63 + t60 * t86 + (-t56 / 0.2e1 + t57 / 0.2e1 + t149 / 0.2e1 + t53 * mrSges(4,3) + t48 * mrSges(5,2) + (t91 / 0.2e1 + t90 / 0.2e1) * t122 + (t77 + t80) * pkin(8)) * t125 + (t58 / 0.2e1 + t59 / 0.2e1 + t49 * mrSges(5,2) - t52 * mrSges(4,3) + (t88 / 0.2e1 - t89 / 0.2e1) * t122 + (-t78 + t79) * pkin(8)) * t121 + (-t126 * mrSges(3,2) + (-mrSges(3,1) + t87) * t122) * pkin(7) + t73 * t23 / 0.2e1 + t61 * t42 / 0.2e1 + t62 * t43 / 0.2e1 - pkin(2) * t64 + t65 * t26 + t45 * t41 + t47 * t50 + t46 * t51 + t35 * t6 / 0.2e1 + t36 * t7 / 0.2e1 + t44 * t8 + t24 * t15 / 0.2e1 + t25 * t16 / 0.2e1 + t11 * t17 + t10 * t18 + t19 * t14; -0.2e1 * pkin(2) * t87 + 0.2e1 * t44 * t14 + t35 * t15 + t36 * t16 + 0.2e1 * t65 * t41 - t136 * t42 + t73 * t43 + 0.2e1 * t84 * t86 + Ifges(3,3) + (-t88 + t89) * t125 + (t90 + t91) * t121 + m(5) * (t84 ^ 2 + t146) + m(4) * (pkin(2) ^ 2 + t146) + m(6) * (t46 ^ 2 + t47 ^ 2 + t65 ^ 2) + m(7) * (t10 ^ 2 + t11 ^ 2 + t44 ^ 2) + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * pkin(8) * t165 + 0.2e1 * (-t10 * t36 + t11 * t35) * mrSges(7,3) + 0.2e1 * (-t136 * t47 - t46 * t73) * mrSges(6,3); m(6) * (t12 * t82 + t13 * t83) + m(5) * (-pkin(3) * t49 + qJ(4) * t48) + m(7) * (t2 * t39 + t3 * t40) - t155 * t126 - pkin(3) * t79 + qJ(4) * t80 + t82 * t51 + t83 * t50 - t131 - Ifges(4,6) * t148 + t48 * mrSges(5,3) - t49 * mrSges(5,1) + t52 * mrSges(4,1) - t53 * mrSges(4,2) + t39 * t18 + t40 * t17 - t163; m(6) * (t46 * t82 + t47 * t83) + m(7) * (t10 * t39 + t11 * t40) + (t35 * t40 - t36 * t39) * mrSges(7,3) + (-t136 * t83 - t73 * t82) * mrSges(6,3) + t137 * mrSges(5,2) + (m(5) * t137 - t139 - t140) * pkin(8) + t106 + t107 + t108 - Ifges(5,6) * t125 - t132; 0.2e1 * pkin(3) * mrSges(5,1) - 0.2e1 * t157 + 0.2e1 * t156 + 0.2e1 * t158 + 0.2e1 * qJ(4) * mrSges(5,3) - 0.2e1 * t38 + m(7) * (t39 ^ 2 + t40 ^ 2) + m(6) * (t82 ^ 2 + t83 ^ 2) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2) + t154 + t155; t120 * t50 + t124 * t51 + t72 * t17 + t70 * t18 + m(7) * (t2 * t70 + t3 * t72) + m(6) * (t12 * t124 + t120 * t13) + m(5) * t49 + t79; (m(5) * pkin(8) + mrSges(5,2)) * t121 + (t35 * t72 - t36 * t70) * mrSges(7,3) + (-t120 * t136 - t124 * t73) * mrSges(6,3) + m(7) * (t10 * t70 + t11 * t72) + m(6) * (t120 * t47 + t124 * t46); -m(5) * pkin(3) - mrSges(5,1) + m(7) * (t39 * t70 + t40 * t72) + m(6) * (t120 * t83 + t124 * t82) - t133; m(5) + m(6) * (t120 ^ 2 + t124 ^ 2) + m(7) * (t70 ^ 2 + t72 ^ 2); (m(7) * (t119 * t3 + t123 * t2) + t119 * t17 + t123 * t18) * pkin(5) + t131; (m(7) * (t10 * t123 + t11 * t119) + (t119 * t35 - t123 * t36) * mrSges(7,3)) * pkin(5) + t132; t157 - t156 - Ifges(6,3) + (m(7) * (t119 * t40 + t123 * t39) - t138) * pkin(5) + t143; m(7) * (t119 * t72 + t123 * t70) * pkin(5) + t133; m(7) * (t119 ^ 2 + t123 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t134 + t154; -t135; -t141; t143; t144; Ifges(7,3) + t134; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
