% Calculate joint inertia matrix for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:02:23
% EndTime: 2019-03-08 22:02:26
% DurationCPUTime: 1.30s
% Computational Cost: add. (2234->314), mult. (5722->466), div. (0->0), fcn. (6374->14), ass. (0->125)
t161 = Ifges(4,3) + Ifges(5,3);
t108 = sin(pkin(13));
t95 = pkin(3) * t108 + pkin(10);
t160 = 0.2e1 * t95;
t114 = sin(qJ(6));
t118 = cos(qJ(6));
t112 = cos(pkin(7));
t115 = sin(qJ(5));
t119 = cos(qJ(5));
t109 = sin(pkin(7));
t111 = cos(pkin(13));
t116 = sin(qJ(3));
t120 = cos(qJ(3));
t62 = (t108 * t120 + t111 * t116) * t109;
t49 = t112 * t115 + t119 * t62;
t139 = t109 * t120;
t140 = t109 * t116;
t61 = t108 * t140 - t111 * t139;
t33 = -t114 * t49 + t118 * t61;
t48 = -t119 * t112 + t115 * t62;
t17 = -mrSges(7,2) * t48 + mrSges(7,3) * t33;
t34 = t114 * t61 + t118 * t49;
t18 = mrSges(7,1) * t48 - mrSges(7,3) * t34;
t159 = -t114 * t18 + t118 * t17;
t78 = -mrSges(7,1) * t118 + mrSges(7,2) * t114;
t158 = -m(7) * pkin(5) - mrSges(6,1) + t78;
t110 = sin(pkin(6));
t113 = cos(pkin(6));
t117 = sin(qJ(2));
t121 = cos(qJ(2));
t138 = t112 * t121;
t43 = t113 * t139 + (-t116 * t117 + t120 * t138) * t110;
t44 = t113 * t140 + (t116 * t138 + t117 * t120) * t110;
t23 = t108 * t43 + t111 * t44;
t67 = -t109 * t110 * t121 + t112 * t113;
t14 = t115 * t23 - t67 * t119;
t157 = t14 ^ 2;
t21 = t108 * t44 - t111 * t43;
t156 = t21 ^ 2;
t155 = t33 / 0.2e1;
t154 = t34 / 0.2e1;
t153 = -t114 / 0.2e1;
t152 = t114 / 0.2e1;
t151 = t118 / 0.2e1;
t59 = Ifges(5,5) * t62;
t150 = pkin(2) * t112;
t149 = pkin(5) * t119;
t13 = -mrSges(7,1) * t33 + mrSges(7,2) * t34;
t37 = mrSges(6,1) * t61 - mrSges(6,3) * t49;
t148 = t13 - t37;
t91 = t120 * t150;
t50 = pkin(3) * t112 + t91 + (-pkin(9) - qJ(4)) * t140;
t70 = pkin(9) * t139 + t116 * t150;
t55 = qJ(4) * t139 + t70;
t31 = t108 * t50 + t111 * t55;
t27 = pkin(10) * t112 + t31;
t75 = (-pkin(3) * t120 - pkin(2)) * t109;
t35 = pkin(4) * t61 - pkin(10) * t62 + t75;
t12 = t115 * t35 + t119 * t27;
t147 = Ifges(7,4) * t114;
t146 = Ifges(7,4) * t118;
t144 = t115 * t95;
t142 = t119 * t14;
t141 = t119 * t95;
t80 = Ifges(7,5) * t114 + Ifges(7,6) * t118;
t137 = t114 * t115;
t136 = t115 * t118;
t135 = Ifges(6,5) * t115 + Ifges(6,6) * t119;
t134 = t114 ^ 2 + t118 ^ 2;
t105 = t115 ^ 2;
t107 = t119 ^ 2;
t133 = t105 + t107;
t5 = Ifges(7,5) * t34 + Ifges(7,6) * t33 + Ifges(7,3) * t48;
t132 = Ifges(6,5) * t49 - Ifges(6,6) * t48 + Ifges(6,3) * t61;
t96 = -pkin(3) * t111 - pkin(4);
t38 = t61 * mrSges(5,1) + t62 * mrSges(5,2);
t30 = -t108 * t55 + t111 * t50;
t131 = t134 * t115;
t26 = -pkin(4) * t112 - t30;
t10 = pkin(5) * t48 - pkin(11) * t49 + t26;
t9 = pkin(11) * t61 + t12;
t1 = t10 * t118 - t114 * t9;
t2 = t10 * t114 + t118 * t9;
t130 = -t1 * t114 + t118 * t2;
t16 = t115 * t67 + t119 * t23;
t3 = -t114 * t16 + t118 * t21;
t4 = t114 * t21 + t118 * t16;
t129 = -t114 * t3 + t118 * t4;
t79 = -t119 * mrSges(6,1) + t115 * mrSges(6,2);
t128 = mrSges(7,1) * t114 + mrSges(7,2) * t118;
t72 = -pkin(11) * t115 - t149 + t96;
t46 = -t114 * t141 + t118 * t72;
t47 = t114 * t72 + t118 * t141;
t127 = -t114 * t46 + t118 * t47;
t76 = mrSges(7,2) * t119 - mrSges(7,3) * t137;
t77 = -mrSges(7,1) * t119 - mrSges(7,3) * t136;
t126 = -t114 * t77 + t118 * t76;
t125 = t115 * t14 + t119 * t16;
t11 = -t115 * t27 + t119 * t35;
t124 = Ifges(4,5) * t140 + Ifges(4,6) * t139 - Ifges(5,6) * t61 + t112 * t161 + t59;
t63 = Ifges(7,5) * t136 - Ifges(7,6) * t137 - Ifges(7,3) * t119;
t93 = t95 ^ 2;
t85 = t105 * t93;
t84 = Ifges(6,1) * t115 + Ifges(6,4) * t119;
t83 = Ifges(7,1) * t114 + t146;
t82 = Ifges(6,4) * t115 + Ifges(6,2) * t119;
t81 = Ifges(7,2) * t118 + t147;
t74 = -mrSges(4,2) * t112 + mrSges(4,3) * t139;
t73 = mrSges(4,1) * t112 - mrSges(4,3) * t140;
t71 = t128 * t115;
t69 = -pkin(9) * t140 + t91;
t68 = (-mrSges(4,1) * t120 + mrSges(4,2) * t116) * t109;
t66 = t67 ^ 2;
t65 = -Ifges(7,5) * t119 + (Ifges(7,1) * t118 - t147) * t115;
t64 = -Ifges(7,6) * t119 + (-Ifges(7,2) * t114 + t146) * t115;
t52 = mrSges(5,1) * t112 - mrSges(5,3) * t62;
t51 = -mrSges(5,2) * t112 - mrSges(5,3) * t61;
t36 = -mrSges(6,2) * t61 - mrSges(6,3) * t48;
t25 = mrSges(6,1) * t48 + mrSges(6,2) * t49;
t20 = Ifges(6,1) * t49 - Ifges(6,4) * t48 + Ifges(6,5) * t61;
t19 = Ifges(6,4) * t49 - Ifges(6,2) * t48 + Ifges(6,6) * t61;
t8 = -pkin(5) * t61 - t11;
t7 = Ifges(7,1) * t34 + Ifges(7,4) * t33 + Ifges(7,5) * t48;
t6 = Ifges(7,4) * t34 + Ifges(7,2) * t33 + Ifges(7,6) * t48;
t15 = [m(2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t157) + m(6) * (t16 ^ 2 + t156 + t157) + m(5) * (t23 ^ 2 + t156 + t66) + m(4) * (t43 ^ 2 + t44 ^ 2 + t66) + m(3) * (t113 ^ 2 + (t117 ^ 2 + t121 ^ 2) * t110 ^ 2); t16 * t36 + t4 * t17 + t3 * t18 + t23 * t51 + t43 * t73 + t44 * t74 + (t38 + t68) * t67 + (t25 - t52) * t21 + t148 * t14 + (mrSges(3,1) * t121 - mrSges(3,2) * t117) * t110 + m(7) * (t1 * t3 + t14 * t8 + t2 * t4) + m(6) * (-t11 * t14 + t12 * t16 + t21 * t26) + m(5) * (-t21 * t30 + t23 * t31 + t67 * t75) + m(4) * (-pkin(2) * t109 * t67 + t43 * t69 + t44 * t70); Ifges(5,1) * t62 ^ 2 + 0.2e1 * t1 * t18 + 0.2e1 * t11 * t37 + 0.2e1 * t12 * t36 + 0.2e1 * t8 * t13 + 0.2e1 * t2 * t17 + t49 * t20 + 0.2e1 * t26 * t25 + 0.2e1 * t30 * t52 + 0.2e1 * t31 * t51 + t33 * t6 + t34 * t7 + 0.2e1 * t75 * t38 + 0.2e1 * t69 * t73 + 0.2e1 * t70 * t74 + Ifges(3,3) + (t5 - t19) * t48 + (t124 + t59) * t112 + (-0.2e1 * Ifges(5,4) * t62 + Ifges(5,2) * t61 - Ifges(5,6) * t112 + t132) * t61 + m(4) * (t69 ^ 2 + t70 ^ 2) + m(5) * (t30 ^ 2 + t31 ^ 2 + t75 ^ 2) + m(6) * (t11 ^ 2 + t12 ^ 2 + t26 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t8 ^ 2) + (-0.2e1 * pkin(2) * t68 + (t116 * Ifges(4,5) + t120 * Ifges(4,6)) * t112 + (Ifges(4,2) * t120 ^ 2 + (Ifges(4,1) * t116 + 0.2e1 * Ifges(4,4) * t120) * t116 + m(4) * pkin(2) ^ 2) * t109) * t109; t43 * mrSges(4,1) - t44 * mrSges(4,2) - t23 * mrSges(5,2) + t14 * t71 + t3 * t77 + t4 * t76 + (-mrSges(5,1) + t79) * t21 + t125 * mrSges(6,3) + m(7) * (t14 * t144 + t3 * t46 + t4 * t47) + m(6) * (t125 * t95 + t21 * t96) + m(5) * (t108 * t23 - t111 * t21) * pkin(3); m(7) * (t1 * t46 + t8 * t144 + t2 * t47) + (-t82 / 0.2e1 + t63 / 0.2e1) * t48 + t96 * t25 + t61 * t135 / 0.2e1 + t124 + t2 * t76 + t1 * t77 + t26 * t79 + t49 * t84 / 0.2e1 + t64 * t155 + t65 * t154 + t69 * mrSges(4,1) - t70 * mrSges(4,2) + t8 * t71 + t46 * t18 + t47 * t17 - t31 * mrSges(5,2) + t30 * mrSges(5,1) + (t108 * t51 + m(5) * (t108 * t31 + t111 * t30) + t111 * t52) * pkin(3) + m(6) * (t26 * t96 + (-t11 * t115 + t12 * t119) * t95) + (t20 / 0.2e1 + t7 * t151 + t6 * t153 - t11 * mrSges(6,3) + t148 * t95) * t115 + (t19 / 0.2e1 - t5 / 0.2e1 + t95 * t36 + t12 * mrSges(6,3)) * t119; 0.2e1 * t46 * t77 + 0.2e1 * t47 * t76 + 0.2e1 * t96 * t79 + (-t63 + t82) * t119 + m(7) * (t46 ^ 2 + t47 ^ 2 + t85) + m(6) * (t107 * t93 + t96 ^ 2 + t85) + m(5) * (t108 ^ 2 + t111 ^ 2) * pkin(3) ^ 2 + (-t114 * t64 + t118 * t65 + t71 * t160 + t84) * t115 + 0.2e1 * (t111 * mrSges(5,1) - t108 * mrSges(5,2)) * pkin(3) + t133 * mrSges(6,3) * t160 + t161; m(7) * (t129 * t115 - t142) + m(6) * (t115 * t16 - t142) + m(5) * t67; -t148 * t119 + (t36 + t159) * t115 + m(7) * (t130 * t115 - t119 * t8) + m(6) * (t11 * t119 + t115 * t12) + m(5) * t75 + t38; -t119 * t71 + (m(7) * (t127 - t141) + t126) * t115; m(5) + m(6) * t133 + m(7) * (t134 * t105 + t107); -t16 * mrSges(6,2) + (m(7) * pkin(11) + mrSges(7,3)) * t129 + t158 * t14; t8 * t78 + t6 * t151 + t7 * t152 + t83 * t154 + t81 * t155 + t48 * t80 / 0.2e1 - t12 * mrSges(6,2) + t11 * mrSges(6,1) + t130 * mrSges(7,3) + t132 + (-m(7) * t8 - t13) * pkin(5) + (m(7) * t130 + t159) * pkin(11); t65 * t152 + t64 * t151 - pkin(5) * t71 + (m(7) * t127 + t126) * pkin(11) + (t83 * t151 + t81 * t153 + t158 * t95) * t115 + (-t80 / 0.2e1 - t95 * mrSges(6,2)) * t119 + t127 * mrSges(7,3) + t135; -t119 * t78 + m(7) * (pkin(11) * t131 + t149) + mrSges(7,3) * t131 - t79; Ifges(6,3) + t118 * t81 - 0.2e1 * pkin(5) * t78 + t114 * t83 + m(7) * (t134 * pkin(11) ^ 2 + pkin(5) ^ 2) + 0.2e1 * t134 * pkin(11) * mrSges(7,3); mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t5; mrSges(7,1) * t46 - mrSges(7,2) * t47 + t63; -t71; -t128 * pkin(11) + t80; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
