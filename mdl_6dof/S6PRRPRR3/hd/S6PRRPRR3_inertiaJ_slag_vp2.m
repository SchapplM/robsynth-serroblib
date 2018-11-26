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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:16
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:15:55
% EndTime: 2018-11-23 15:15:56
% DurationCPUTime: 1.30s
% Computational Cost: add. (2234->314), mult. (5722->468), div. (0->0), fcn. (6374->14), ass. (0->123)
t161 = Ifges(4,3) + Ifges(5,3);
t116 = sin(qJ(6));
t120 = cos(qJ(6));
t114 = cos(pkin(7));
t117 = sin(qJ(5));
t121 = cos(qJ(5));
t110 = sin(pkin(13));
t111 = sin(pkin(7));
t113 = cos(pkin(13));
t118 = sin(qJ(3));
t122 = cos(qJ(3));
t62 = (t110 * t122 + t113 * t118) * t111;
t49 = t114 * t117 + t121 * t62;
t140 = t111 * t122;
t141 = t111 * t118;
t61 = t110 * t141 - t113 * t140;
t33 = -t116 * t49 + t120 * t61;
t48 = -t121 * t114 + t117 * t62;
t17 = -mrSges(7,2) * t48 + mrSges(7,3) * t33;
t34 = t116 * t61 + t120 * t49;
t18 = mrSges(7,1) * t48 - mrSges(7,3) * t34;
t160 = -t116 * t18 + t120 * t17;
t78 = -mrSges(7,1) * t120 + mrSges(7,2) * t116;
t159 = -m(7) * pkin(5) - mrSges(6,1) + t78;
t112 = sin(pkin(6));
t115 = cos(pkin(6));
t119 = sin(qJ(2));
t123 = cos(qJ(2));
t139 = t114 * t123;
t43 = t115 * t140 + (-t118 * t119 + t122 * t139) * t112;
t44 = t115 * t141 + (t118 * t139 + t119 * t122) * t112;
t23 = t110 * t43 + t113 * t44;
t67 = -t111 * t112 * t123 + t114 * t115;
t14 = t117 * t23 - t121 * t67;
t158 = t14 ^ 2;
t21 = t110 * t44 - t113 * t43;
t157 = t21 ^ 2;
t156 = t33 / 0.2e1;
t155 = t34 / 0.2e1;
t154 = -t116 / 0.2e1;
t153 = t116 / 0.2e1;
t152 = t120 / 0.2e1;
t59 = Ifges(5,5) * t62;
t151 = pkin(2) * t114;
t150 = pkin(5) * t121;
t13 = -mrSges(7,1) * t33 + mrSges(7,2) * t34;
t37 = mrSges(6,1) * t61 - mrSges(6,3) * t49;
t149 = t13 - t37;
t91 = t122 * t151;
t50 = pkin(3) * t114 + t91 + (-pkin(9) - qJ(4)) * t141;
t70 = pkin(9) * t140 + t118 * t151;
t55 = qJ(4) * t140 + t70;
t31 = t110 * t50 + t113 * t55;
t27 = pkin(10) * t114 + t31;
t75 = (-pkin(3) * t122 - pkin(2)) * t111;
t35 = pkin(4) * t61 - pkin(10) * t62 + t75;
t12 = t117 * t35 + t121 * t27;
t148 = Ifges(7,4) * t116;
t147 = Ifges(7,4) * t120;
t97 = pkin(3) * t110 + pkin(10);
t145 = t117 * t97;
t143 = t121 * t14;
t142 = t121 * t97;
t138 = t116 * t117;
t137 = t117 * t120;
t80 = Ifges(7,5) * t116 + Ifges(7,6) * t120;
t136 = Ifges(6,5) * t117 + Ifges(6,6) * t121;
t135 = t116 ^ 2 + t120 ^ 2;
t107 = t117 ^ 2;
t109 = t121 ^ 2;
t134 = t107 + t109;
t5 = Ifges(7,5) * t34 + Ifges(7,6) * t33 + Ifges(7,3) * t48;
t133 = Ifges(6,5) * t49 - Ifges(6,6) * t48 + Ifges(6,3) * t61;
t98 = -pkin(3) * t113 - pkin(4);
t38 = t61 * mrSges(5,1) + t62 * mrSges(5,2);
t30 = -t110 * t55 + t113 * t50;
t132 = t135 * t117;
t26 = -pkin(4) * t114 - t30;
t10 = pkin(5) * t48 - pkin(11) * t49 + t26;
t9 = pkin(11) * t61 + t12;
t1 = t10 * t120 - t116 * t9;
t2 = t10 * t116 + t120 * t9;
t131 = -t1 * t116 + t120 * t2;
t16 = t117 * t67 + t121 * t23;
t3 = -t116 * t16 + t120 * t21;
t4 = t116 * t21 + t120 * t16;
t130 = -t116 * t3 + t120 * t4;
t79 = -t121 * mrSges(6,1) + t117 * mrSges(6,2);
t72 = -pkin(11) * t117 - t150 + t98;
t46 = -t116 * t142 + t120 * t72;
t47 = t116 * t72 + t120 * t142;
t129 = -t116 * t46 + t120 * t47;
t76 = mrSges(7,2) * t121 - mrSges(7,3) * t138;
t77 = -mrSges(7,1) * t121 - mrSges(7,3) * t137;
t128 = -t116 * t77 + t120 * t76;
t127 = t117 * t14 + t121 * t16;
t11 = -t117 * t27 + t121 * t35;
t126 = Ifges(4,5) * t141 + Ifges(4,6) * t140 - Ifges(5,6) * t61 + t114 * t161 + t59;
t63 = Ifges(7,5) * t137 - Ifges(7,6) * t138 - Ifges(7,3) * t121;
t95 = t97 ^ 2;
t85 = t107 * t95;
t84 = Ifges(6,1) * t117 + Ifges(6,4) * t121;
t83 = Ifges(7,1) * t116 + t147;
t82 = Ifges(6,4) * t117 + Ifges(6,2) * t121;
t81 = Ifges(7,2) * t120 + t148;
t74 = -mrSges(4,2) * t114 + mrSges(4,3) * t140;
t73 = mrSges(4,1) * t114 - mrSges(4,3) * t141;
t71 = -mrSges(7,1) * t138 - mrSges(7,2) * t137;
t69 = -pkin(9) * t141 + t91;
t68 = (-mrSges(4,1) * t122 + mrSges(4,2) * t118) * t111;
t66 = t67 ^ 2;
t65 = -Ifges(7,5) * t121 + (Ifges(7,1) * t120 - t148) * t117;
t64 = -Ifges(7,6) * t121 + (-Ifges(7,2) * t116 + t147) * t117;
t52 = mrSges(5,1) * t114 - mrSges(5,3) * t62;
t51 = -mrSges(5,2) * t114 - mrSges(5,3) * t61;
t36 = -mrSges(6,2) * t61 - mrSges(6,3) * t48;
t25 = mrSges(6,1) * t48 + mrSges(6,2) * t49;
t20 = Ifges(6,1) * t49 - Ifges(6,4) * t48 + Ifges(6,5) * t61;
t19 = Ifges(6,4) * t49 - Ifges(6,2) * t48 + Ifges(6,6) * t61;
t8 = -pkin(5) * t61 - t11;
t7 = Ifges(7,1) * t34 + Ifges(7,4) * t33 + Ifges(7,5) * t48;
t6 = Ifges(7,4) * t34 + Ifges(7,2) * t33 + Ifges(7,6) * t48;
t15 = [m(2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t158) + m(6) * (t16 ^ 2 + t157 + t158) + m(5) * (t23 ^ 2 + t157 + t66) + m(4) * (t43 ^ 2 + t44 ^ 2 + t66) + m(3) * (t115 ^ 2 + (t119 ^ 2 + t123 ^ 2) * t112 ^ 2); t16 * t36 + t4 * t17 + t3 * t18 + t23 * t51 + t43 * t73 + t44 * t74 + (t38 + t68) * t67 + (t25 - t52) * t21 + t149 * t14 + (t123 * mrSges(3,1) - t119 * mrSges(3,2)) * t112 + m(7) * (t1 * t3 + t14 * t8 + t2 * t4) + m(6) * (-t11 * t14 + t12 * t16 + t21 * t26) + m(5) * (-t21 * t30 + t23 * t31 + t67 * t75) + m(4) * (-pkin(2) * t111 * t67 + t43 * t69 + t44 * t70); Ifges(5,1) * t62 ^ 2 + 0.2e1 * t1 * t18 + 0.2e1 * t11 * t37 + 0.2e1 * t12 * t36 + 0.2e1 * t8 * t13 + 0.2e1 * t2 * t17 + t49 * t20 + 0.2e1 * t26 * t25 + 0.2e1 * t30 * t52 + 0.2e1 * t31 * t51 + t33 * t6 + t34 * t7 + 0.2e1 * t75 * t38 + 0.2e1 * t69 * t73 + 0.2e1 * t70 * t74 + Ifges(3,3) + (t5 - t19) * t48 + (t126 + t59) * t114 + (-0.2e1 * Ifges(5,4) * t62 + Ifges(5,2) * t61 - Ifges(5,6) * t114 + t133) * t61 + m(4) * (t69 ^ 2 + t70 ^ 2) + m(5) * (t30 ^ 2 + t31 ^ 2 + t75 ^ 2) + m(6) * (t11 ^ 2 + t12 ^ 2 + t26 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t8 ^ 2) + (-0.2e1 * pkin(2) * t68 + (t118 * Ifges(4,5) + t122 * Ifges(4,6)) * t114 + (Ifges(4,2) * t122 ^ 2 + (Ifges(4,1) * t118 + 0.2e1 * Ifges(4,4) * t122) * t118 + m(4) * pkin(2) ^ 2) * t111) * t111; t43 * mrSges(4,1) - t44 * mrSges(4,2) - t23 * mrSges(5,2) - t14 * t71 + t3 * t77 + t4 * t76 + (-mrSges(5,1) + t79) * t21 + t127 * mrSges(6,3) + m(7) * (t14 * t145 + t3 * t46 + t4 * t47) + m(6) * (t127 * t97 + t21 * t98) + m(5) * (t110 * t23 - t113 * t21) * pkin(3); m(7) * (t1 * t46 + t8 * t145 + t2 * t47) + (-t82 / 0.2e1 + t63 / 0.2e1) * t48 + (m(5) * (t110 * t31 + t113 * t30) + t113 * t52 + t110 * t51) * pkin(3) + t49 * t84 / 0.2e1 + t98 * t25 + t61 * t136 / 0.2e1 + t69 * mrSges(4,1) - t70 * mrSges(4,2) - t8 * t71 + t2 * t76 + t1 * t77 + t26 * t79 + t64 * t156 + t65 * t155 + t126 + t46 * t18 + t47 * t17 + t30 * mrSges(5,1) - t31 * mrSges(5,2) + m(6) * (t26 * t98 + (-t11 * t117 + t12 * t121) * t97) + (t20 / 0.2e1 + t7 * t152 + t6 * t154 - t11 * mrSges(6,3) + t149 * t97) * t117 + (t19 / 0.2e1 - t5 / 0.2e1 + t97 * t36 + t12 * mrSges(6,3)) * t121; 0.2e1 * t46 * t77 + 0.2e1 * t47 * t76 + 0.2e1 * t98 * t79 + (-t63 + t82) * t121 + m(7) * (t46 ^ 2 + t47 ^ 2 + t85) + m(6) * (t109 * t95 + t98 ^ 2 + t85) + m(5) * (t110 ^ 2 + t113 ^ 2) * pkin(3) ^ 2 + (-t116 * t64 + t120 * t65 - 0.2e1 * t71 * t97 + t84) * t117 + 0.2e1 * (mrSges(5,1) * t113 - mrSges(5,2) * t110) * pkin(3) + 0.2e1 * t134 * t97 * mrSges(6,3) + t161; m(7) * (t130 * t117 - t143) + m(6) * (t117 * t16 - t143) + m(5) * t67; -t149 * t121 + (t36 + t160) * t117 + m(7) * (t131 * t117 - t121 * t8) + m(6) * (t11 * t121 + t117 * t12) + m(5) * t75 + t38; t121 * t71 + (m(7) * (t129 - t142) + t128) * t117; m(5) + m(7) * (t135 * t107 + t109) + m(6) * t134; -t16 * mrSges(6,2) + (m(7) * pkin(11) + mrSges(7,3)) * t130 + t159 * t14; t8 * t78 + t81 * t156 + t48 * t80 / 0.2e1 + t83 * t155 + t6 * t152 + t7 * t153 - t12 * mrSges(6,2) + t11 * mrSges(6,1) + t131 * mrSges(7,3) + t133 + (-m(7) * t8 - t13) * pkin(5) + (m(7) * t131 + t160) * pkin(11); t64 * t152 + t65 * t153 + pkin(5) * t71 + (m(7) * t129 + t128) * pkin(11) + (t83 * t152 + t81 * t154 + t159 * t97) * t117 + (-t80 / 0.2e1 - t97 * mrSges(6,2)) * t121 + t129 * mrSges(7,3) + t136; -t121 * t78 + m(7) * (pkin(11) * t132 + t150) + mrSges(7,3) * t132 - t79; Ifges(6,3) + m(7) * (t135 * pkin(11) ^ 2 + pkin(5) ^ 2) + t116 * t83 + t120 * t81 - 0.2e1 * pkin(5) * t78 + 0.2e1 * t135 * pkin(11) * mrSges(7,3); mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t5; mrSges(7,1) * t46 - mrSges(7,2) * t47 + t63; t71; (-mrSges(7,1) * t116 - mrSges(7,2) * t120) * pkin(11) + t80; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
