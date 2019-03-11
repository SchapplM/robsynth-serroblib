% Calculate joint inertia matrix for
% S6RRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:38:58
% EndTime: 2019-03-09 22:39:01
% DurationCPUTime: 1.29s
% Computational Cost: add. (2210->358), mult. (4225->489), div. (0->0), fcn. (4230->8), ass. (0->128)
t162 = 2 * pkin(7);
t119 = sin(qJ(6));
t123 = cos(qJ(6));
t126 = cos(qJ(2));
t120 = sin(qJ(4));
t124 = cos(qJ(4));
t121 = sin(qJ(3));
t122 = sin(qJ(2));
t125 = cos(qJ(3));
t147 = t122 * t125;
t94 = -pkin(2) * t126 - pkin(8) * t122 - pkin(1);
t83 = t125 * t94;
t44 = -pkin(9) * t147 + t83 + (-pkin(7) * t121 - pkin(3)) * t126;
t148 = t121 * t122;
t157 = pkin(7) * t126;
t65 = t121 * t94 + t125 * t157;
t53 = -pkin(9) * t148 + t65;
t20 = -t120 * t53 + t124 * t44;
t14 = t126 * pkin(4) - t20;
t84 = t120 * t121 - t124 * t125;
t76 = t84 * t122;
t4 = pkin(5) * t126 + pkin(10) * t76 + t14;
t21 = t120 * t44 + t124 * t53;
t13 = -qJ(5) * t126 + t21;
t85 = t120 * t125 + t121 * t124;
t75 = t85 * t122;
t7 = pkin(10) * t75 + t13;
t2 = -t119 * t7 + t123 * t4;
t3 = t119 * t4 + t123 * t7;
t161 = t2 * mrSges(7,1) - t3 * mrSges(7,2);
t160 = 2 * mrSges(6,3);
t127 = -pkin(4) - pkin(5);
t159 = -pkin(9) - pkin(8);
t113 = t122 * pkin(7);
t105 = -pkin(3) * t124 - pkin(4);
t102 = -pkin(5) + t105;
t103 = pkin(3) * t120 + qJ(5);
t66 = t102 * t123 - t103 * t119;
t156 = t66 * mrSges(7,1);
t67 = t102 * t119 + t103 * t123;
t155 = t67 * mrSges(7,2);
t92 = -qJ(5) * t119 + t123 * t127;
t154 = t92 * mrSges(7,1);
t93 = qJ(5) * t123 + t119 * t127;
t153 = t93 * mrSges(7,2);
t151 = -Ifges(5,3) - Ifges(6,2);
t141 = t159 * t121;
t98 = t159 * t125;
t59 = t120 * t141 - t124 * t98;
t150 = Ifges(4,4) * t121;
t149 = Ifges(4,4) * t125;
t63 = t126 * mrSges(6,1) - t76 * mrSges(6,2);
t91 = pkin(3) * t148 + t113;
t146 = t121 ^ 2 + t125 ^ 2;
t57 = -t120 * t98 - t124 * t141;
t145 = t57 ^ 2 + t59 ^ 2;
t144 = Ifges(7,3) - t151;
t143 = -Ifges(4,3) + t151;
t32 = t119 * t76 + t123 * t75;
t33 = t119 * t75 - t123 * t76;
t142 = Ifges(7,5) * t33 + Ifges(7,6) * t32 + Ifges(7,3) * t126;
t106 = -pkin(3) * t125 - pkin(2);
t140 = t123 * mrSges(7,1) - t119 * mrSges(7,2);
t139 = -mrSges(6,1) - t140;
t138 = -qJ(5) * t76 - t91;
t36 = -pkin(10) * t85 + t57;
t37 = pkin(10) * t84 + t59;
t10 = -t119 * t37 + t123 * t36;
t11 = t119 * t36 + t123 * t37;
t42 = -t119 * t85 + t123 * t84;
t40 = Ifges(7,6) * t42;
t43 = t119 * t84 + t123 * t85;
t41 = Ifges(7,5) * t43;
t137 = t10 * mrSges(7,1) - t11 * mrSges(7,2) + t40 + t41;
t136 = mrSges(4,1) * t121 + mrSges(4,2) * t125;
t135 = qJ(5) * t85 - t106;
t134 = (mrSges(5,1) * t124 - mrSges(5,2) * t120) * pkin(3);
t133 = t142 + (Ifges(6,4) + Ifges(5,5)) * t76 + (Ifges(5,6) - Ifges(6,6)) * t75;
t78 = Ifges(6,6) * t84;
t79 = Ifges(5,6) * t84;
t80 = Ifges(5,5) * t85;
t81 = Ifges(6,4) * t85;
t132 = -t137 + t78 - t79 + t80 + t81 + (-mrSges(5,2) + mrSges(6,3)) * t59 + (-mrSges(5,1) - mrSges(6,1)) * t57;
t131 = t20 * mrSges(5,1) - t14 * mrSges(6,1) - t21 * mrSges(5,2) + t13 * mrSges(6,3) - t133 - t161;
t129 = pkin(7) ^ 2;
t118 = t126 ^ 2;
t116 = t122 ^ 2;
t112 = t116 * t129;
t111 = Ifges(4,5) * t121;
t110 = Ifges(4,6) * t125;
t99 = Ifges(4,5) * t147;
t97 = Ifges(4,1) * t121 + t149;
t96 = Ifges(4,2) * t125 + t150;
t95 = -mrSges(4,1) * t125 + mrSges(4,2) * t121;
t90 = -mrSges(4,1) * t126 - mrSges(4,3) * t147;
t89 = mrSges(4,2) * t126 - mrSges(4,3) * t148;
t77 = t136 * t122;
t74 = -Ifges(4,5) * t126 + (Ifges(4,1) * t125 - t150) * t122;
t73 = -Ifges(4,6) * t126 + (-Ifges(4,2) * t121 + t149) * t122;
t64 = -t121 * t157 + t83;
t62 = -mrSges(5,1) * t126 + mrSges(5,3) * t76;
t61 = mrSges(5,2) * t126 - mrSges(5,3) * t75;
t60 = -mrSges(6,2) * t75 - mrSges(6,3) * t126;
t52 = Ifges(5,1) * t85 - Ifges(5,4) * t84;
t51 = Ifges(6,1) * t85 + Ifges(6,5) * t84;
t50 = Ifges(5,4) * t85 - Ifges(5,2) * t84;
t49 = Ifges(6,5) * t85 + Ifges(6,3) * t84;
t48 = mrSges(5,1) * t84 + mrSges(5,2) * t85;
t47 = mrSges(6,1) * t84 - mrSges(6,3) * t85;
t39 = pkin(4) * t84 - t135;
t35 = mrSges(5,1) * t75 - mrSges(5,2) * t76;
t34 = mrSges(6,1) * t75 + mrSges(6,3) * t76;
t31 = -Ifges(5,1) * t76 - Ifges(5,4) * t75 - Ifges(5,5) * t126;
t30 = -Ifges(6,1) * t76 - Ifges(6,4) * t126 + Ifges(6,5) * t75;
t29 = -Ifges(5,4) * t76 - Ifges(5,2) * t75 - Ifges(5,6) * t126;
t28 = -Ifges(6,5) * t76 - Ifges(6,6) * t126 + Ifges(6,3) * t75;
t25 = t127 * t84 + t135;
t24 = mrSges(7,1) * t126 - mrSges(7,3) * t33;
t23 = -mrSges(7,2) * t126 + mrSges(7,3) * t32;
t22 = pkin(4) * t75 - t138;
t18 = Ifges(7,1) * t43 + Ifges(7,4) * t42;
t17 = Ifges(7,4) * t43 + Ifges(7,2) * t42;
t16 = -mrSges(7,1) * t42 + mrSges(7,2) * t43;
t15 = t127 * t75 + t138;
t8 = -mrSges(7,1) * t32 + mrSges(7,2) * t33;
t6 = Ifges(7,1) * t33 + Ifges(7,4) * t32 + Ifges(7,5) * t126;
t5 = Ifges(7,4) * t33 + Ifges(7,2) * t32 + Ifges(7,6) * t126;
t1 = [0.2e1 * t13 * t60 + 0.2e1 * t14 * t63 + 0.2e1 * t15 * t8 + 0.2e1 * t2 * t24 + 0.2e1 * t20 * t62 + 0.2e1 * t21 * t61 + 0.2e1 * t22 * t34 + 0.2e1 * t3 * t23 + t32 * t5 + t33 * t6 + 0.2e1 * t91 * t35 + 0.2e1 * t64 * t90 + 0.2e1 * t65 * t89 + Ifges(2,3) - (t30 + t31) * t76 + (-t29 + t28) * t75 + (t116 + t118) * mrSges(3,3) * t162 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t122 - t121 * t73 + t125 * t74 + t77 * t162) * t122 + m(3) * (pkin(1) ^ 2 + t118 * t129 + t112) + m(4) * (t64 ^ 2 + t65 ^ 2 + t112) + m(5) * (t20 ^ 2 + t21 ^ 2 + t91 ^ 2) + m(6) * (t13 ^ 2 + t14 ^ 2 + t22 ^ 2) + m(7) * (t15 ^ 2 + t2 ^ 2 + t3 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) - t99 + (Ifges(3,2) - t143) * t126 + (Ifges(4,6) * t121 + (2 * Ifges(3,4))) * t122 + t133) * t126; (-t2 * t43 + t3 * t42) * mrSges(7,3) + (pkin(8) * t89 + t65 * mrSges(4,3) + t73 / 0.2e1) * t125 + (t74 / 0.2e1 - pkin(8) * t90 - t64 * mrSges(4,3)) * t121 + (t60 + t61) * t59 + t106 * t35 + t91 * t48 - pkin(2) * t77 + t39 * t34 + t42 * t5 / 0.2e1 + t43 * t6 / 0.2e1 + t22 * t47 + t32 * t17 / 0.2e1 + t33 * t18 / 0.2e1 + (-t62 + t63) * t57 + (-pkin(7) * mrSges(3,2) - t111 / 0.2e1 - t110 / 0.2e1 - t80 / 0.2e1 + t79 / 0.2e1 - t81 / 0.2e1 - t78 / 0.2e1 + t41 / 0.2e1 + t40 / 0.2e1 + Ifges(3,6)) * t126 + m(5) * (t106 * t91 - t20 * t57 + t21 * t59) + m(6) * (t13 * t59 + t14 * t57 + t22 * t39) + m(7) * (t10 * t2 + t11 * t3 + t15 * t25) + (t30 / 0.2e1 + t31 / 0.2e1 + t14 * mrSges(6,2) - t20 * mrSges(5,3)) * t85 + (t28 / 0.2e1 - t29 / 0.2e1 - t21 * mrSges(5,3) - t13 * mrSges(6,2)) * t84 - (t51 / 0.2e1 + t52 / 0.2e1) * t76 + (-t50 / 0.2e1 + t49 / 0.2e1) * t75 + t11 * t23 + t10 * t24 + t25 * t8 + t15 * t16 + (t125 * t97 / 0.2e1 + Ifges(3,5) - t121 * t96 / 0.2e1 + (t95 - mrSges(3,1)) * pkin(7)) * t122 + m(4) * (-pkin(2) * t113 + (-t64 * t121 + t65 * t125) * pkin(8)); -0.2e1 * pkin(2) * t95 + 0.2e1 * t106 * t48 + t121 * t97 + t125 * t96 + 0.2e1 * t25 * t16 + t42 * t17 + t43 * t18 + 0.2e1 * t39 * t47 + Ifges(3,3) + (t51 + t52) * t85 + (t49 - t50) * t84 + m(4) * (pkin(8) ^ 2 * t146 + pkin(2) ^ 2) + m(5) * (t106 ^ 2 + t145) + m(6) * (t39 ^ 2 + t145) + m(7) * (t10 ^ 2 + t11 ^ 2 + t25 ^ 2) + 0.2e1 * (-t10 * t43 + t11 * t42) * mrSges(7,3) + 0.2e1 * t146 * pkin(8) * mrSges(4,3) + 0.2e1 * (t57 * t85 - t59 * t84) * (mrSges(6,2) + mrSges(5,3)); t99 + t131 + t105 * t63 + t103 * t60 - t65 * mrSges(4,2) + t66 * t24 + t67 * t23 + t64 * mrSges(4,1) + (t124 * t62 + m(5) * (t120 * t21 + t124 * t20) + t120 * t61) * pkin(3) - Ifges(4,6) * t148 + t143 * t126 + m(6) * (t103 * t13 + t105 * t14) + m(7) * (t2 * t66 + t3 * t67); m(6) * (t103 * t59 + t105 * t57) + m(7) * (t10 * t66 + t11 * t67) - t136 * pkin(8) + (t42 * t67 - t43 * t66) * mrSges(7,3) + (-t103 * t84 + t105 * t85) * mrSges(6,2) + t110 + t111 + t132 + (m(5) * (t120 * t59 - t124 * t57) + (-t120 * t84 - t124 * t85) * mrSges(5,3)) * pkin(3); -0.2e1 * t105 * mrSges(6,1) - 0.2e1 * t156 + 0.2e1 * t155 + t103 * t160 + Ifges(7,3) + 0.2e1 * t134 + m(7) * (t66 ^ 2 + t67 ^ 2) + m(6) * (t103 ^ 2 + t105 ^ 2) + m(5) * (t120 ^ 2 + t124 ^ 2) * pkin(3) ^ 2 - t143; m(7) * (t2 * t92 + t3 * t93) + m(6) * (-pkin(4) * t14 + qJ(5) * t13) + t151 * t126 + t131 + t92 * t24 + t93 * t23 + qJ(5) * t60 - pkin(4) * t63; m(7) * (t10 * t92 + t11 * t93) + m(6) * (-pkin(4) * t57 + qJ(5) * t59) + (t42 * t93 - t43 * t92) * mrSges(7,3) + (-pkin(4) * t85 - qJ(5) * t84) * mrSges(6,2) + t132; t134 + (qJ(5) + t103) * mrSges(6,3) + (t93 + t67) * mrSges(7,2) + (-t92 - t66) * mrSges(7,1) + (pkin(4) - t105) * mrSges(6,1) + m(7) * (t66 * t92 + t67 * t93) + m(6) * (-pkin(4) * t105 + qJ(5) * t103) + t144; 0.2e1 * pkin(4) * mrSges(6,1) - 0.2e1 * t154 + 0.2e1 * t153 + qJ(5) * t160 + m(7) * (t92 ^ 2 + t93 ^ 2) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + t144; t119 * t23 + t123 * t24 + m(7) * (t119 * t3 + t123 * t2) + m(6) * t14 + t63; t85 * mrSges(6,2) + (t119 * t42 - t123 * t43) * mrSges(7,3) + m(7) * (t10 * t123 + t11 * t119) + m(6) * t57; m(7) * (t119 * t67 + t123 * t66) + m(6) * t105 + t139; m(7) * (t119 * t93 + t123 * t92) - m(6) * pkin(4) + t139; m(6) + m(7) * (t119 ^ 2 + t123 ^ 2); t142 + t161; t137; -Ifges(7,3) - t155 + t156; -Ifges(7,3) - t153 + t154; t140; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
