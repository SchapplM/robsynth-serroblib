% Calculate joint inertia matrix for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2018-11-23 16:38
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:38:11
% EndTime: 2018-11-23 16:38:12
% DurationCPUTime: 1.02s
% Computational Cost: add. (1934->309), mult. (3751->445), div. (0->0), fcn. (3822->8), ass. (0->118)
t105 = sin(qJ(3));
t108 = cos(qJ(4));
t109 = cos(qJ(3));
t130 = t108 * t109;
t145 = Ifges(5,5) * t130 + Ifges(5,3) * t105;
t103 = sin(qJ(5));
t104 = sin(qJ(4));
t107 = cos(qJ(5));
t72 = t103 * t108 + t104 * t107;
t60 = t72 * t109;
t71 = -t103 * t104 + t107 * t108;
t62 = t71 * t109;
t144 = Ifges(6,5) * t62 - Ifges(6,6) * t60 + Ifges(6,3) * t105;
t143 = 2 * qJ(2);
t142 = -pkin(9) - pkin(8);
t141 = pkin(4) * t103;
t102 = sin(qJ(6));
t106 = cos(qJ(6));
t87 = pkin(4) * t107 + pkin(5);
t64 = t102 * t87 + t106 * t141;
t140 = t64 * mrSges(7,2);
t139 = Ifges(6,3) + Ifges(7,3);
t110 = -pkin(1) - pkin(7);
t78 = pkin(3) * t105 - pkin(8) * t109 + qJ(2);
t69 = t108 * t78;
t36 = -pkin(9) * t130 + t69 + (-t104 * t110 + pkin(4)) * t105;
t132 = t104 * t109;
t131 = t105 * t110;
t51 = t104 * t78 + t108 * t131;
t41 = -pkin(9) * t132 + t51;
t14 = t103 * t36 + t107 * t41;
t82 = t142 * t104;
t83 = t142 * t108;
t47 = t103 * t82 - t107 * t83;
t79 = -mrSges(5,1) * t108 + mrSges(5,2) * t104;
t138 = -t79 + mrSges(4,1);
t137 = t104 ^ 2 + t108 ^ 2;
t136 = Ifges(5,4) * t104;
t135 = Ifges(5,4) * t108;
t134 = t102 * mrSges(7,2);
t100 = t109 ^ 2;
t98 = t105 ^ 2;
t133 = -t98 - t100;
t129 = t109 * t110;
t128 = pkin(5) * t134;
t26 = -t102 * t62 - t106 * t60;
t28 = -t102 * t60 + t106 * t62;
t127 = Ifges(7,5) * t28 + Ifges(7,6) * t26 + Ifges(7,3) * t105;
t88 = -pkin(4) * t108 - pkin(3);
t126 = t137 * mrSges(5,3);
t59 = t72 * t105;
t61 = t71 * t105;
t25 = -t102 * t61 - t106 * t59;
t27 = -t102 * t59 + t106 * t61;
t125 = mrSges(7,1) * t25 - t27 * mrSges(7,2);
t124 = t133 * mrSges(4,3);
t13 = -t103 * t41 + t107 * t36;
t46 = t103 * t83 + t107 * t82;
t70 = pkin(4) * t132 - t129;
t63 = -t102 * t141 + t106 * t87;
t58 = t63 * mrSges(7,1);
t123 = Ifges(7,3) + t58 - t140;
t30 = -pkin(10) * t72 + t46;
t31 = pkin(10) * t71 + t47;
t10 = -t102 * t31 + t106 * t30;
t11 = t102 * t30 + t106 * t31;
t39 = -t102 * t72 + t106 * t71;
t34 = Ifges(7,6) * t39;
t40 = t102 * t71 + t106 * t72;
t35 = Ifges(7,5) * t40;
t122 = mrSges(7,1) * t10 - t11 * mrSges(7,2) + t34 + t35;
t121 = mrSges(5,1) * t104 + mrSges(5,2) * t108;
t50 = -t104 * t131 + t69;
t120 = -t104 * t50 + t108 * t51;
t6 = pkin(5) * t105 - pkin(10) * t62 + t13;
t8 = -pkin(10) * t60 + t14;
t2 = -t102 * t8 + t106 * t6;
t3 = t102 * t6 + t106 * t8;
t119 = mrSges(7,1) * t2 - t3 * mrSges(7,2) + t127;
t118 = (mrSges(6,1) * t107 - mrSges(6,2) * t103) * pkin(4);
t117 = -mrSges(6,1) * t59 - t61 * mrSges(6,2) + t125;
t66 = Ifges(6,6) * t71;
t67 = Ifges(6,5) * t72;
t116 = mrSges(6,1) * t46 - t47 * mrSges(6,2) + t122 + t66 + t67;
t115 = mrSges(6,1) * t13 - t14 * mrSges(6,2) + t119 + t144;
t111 = qJ(2) ^ 2;
t101 = t110 ^ 2;
t96 = Ifges(5,5) * t104;
t95 = Ifges(5,6) * t108;
t91 = t106 * pkin(5) * mrSges(7,1);
t90 = t100 * t110;
t89 = t100 * t101;
t81 = Ifges(5,1) * t104 + t135;
t80 = Ifges(5,2) * t108 + t136;
t77 = mrSges(5,1) * t105 - mrSges(5,3) * t130;
t76 = -mrSges(5,2) * t105 - mrSges(5,3) * t132;
t65 = t121 * t109;
t57 = Ifges(5,5) * t105 + (Ifges(5,1) * t108 - t136) * t109;
t56 = Ifges(5,6) * t105 + (-Ifges(5,2) * t104 + t135) * t109;
t52 = -pkin(5) * t71 + t88;
t49 = mrSges(6,1) * t105 - mrSges(6,3) * t62;
t48 = -mrSges(6,2) * t105 - mrSges(6,3) * t60;
t44 = Ifges(6,1) * t72 + Ifges(6,4) * t71;
t43 = Ifges(6,4) * t72 + Ifges(6,2) * t71;
t42 = -mrSges(6,1) * t71 + mrSges(6,2) * t72;
t38 = pkin(5) * t60 + t70;
t29 = mrSges(6,1) * t60 + mrSges(6,2) * t62;
t24 = Ifges(6,1) * t62 - Ifges(6,4) * t60 + Ifges(6,5) * t105;
t23 = Ifges(6,4) * t62 - Ifges(6,2) * t60 + Ifges(6,6) * t105;
t19 = mrSges(7,1) * t105 - mrSges(7,3) * t28;
t18 = -mrSges(7,2) * t105 + mrSges(7,3) * t26;
t17 = Ifges(7,1) * t40 + Ifges(7,4) * t39;
t16 = Ifges(7,4) * t40 + Ifges(7,2) * t39;
t15 = -mrSges(7,1) * t39 + mrSges(7,2) * t40;
t7 = -mrSges(7,1) * t26 + mrSges(7,2) * t28;
t5 = Ifges(7,1) * t28 + Ifges(7,4) * t26 + Ifges(7,5) * t105;
t4 = Ifges(7,4) * t28 + Ifges(7,2) * t26 + Ifges(7,6) * t105;
t1 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t143) + 0.2e1 * t13 * t49 + 0.2e1 * t14 * t48 + 0.2e1 * t3 * t18 + 0.2e1 * t2 * t19 - t60 * t23 + t62 * t24 + t26 * t4 + t28 * t5 + 0.2e1 * t70 * t29 + 0.2e1 * t38 * t7 + 0.2e1 * t50 * t77 + 0.2e1 * t51 * t76 + Ifges(3,1) + Ifges(2,3) + 0.2e1 * t110 * t124 + ((mrSges(4,2) * t143) + Ifges(4,1) * t109 - t104 * t56 + t108 * t57 - 0.2e1 * t110 * t65) * t109 + (mrSges(4,1) * t143 + Ifges(4,2) * t105 + (-Ifges(5,6) * t104 - (2 * Ifges(4,4))) * t109 + t127 + t144 + t145) * t105 + m(4) * (t101 * t98 + t111 + t89) + (m(3) * (pkin(1) ^ 2 + t111)) + m(5) * (t50 ^ 2 + t51 ^ 2 + t89) + m(6) * (t13 ^ 2 + t14 ^ 2 + t70 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t38 ^ 2); -(m(3) * pkin(1)) + t27 * t18 + t25 * t19 + t61 * t48 - t59 * t49 + mrSges(3,2) + (-t104 * t77 + t108 * t76) * t105 + t124 + (-t29 - t65 - t7) * t109 + m(7) * (-t109 * t38 + t2 * t25 + t27 * t3) + m(6) * (-t109 * t70 - t13 * t59 + t14 * t61) + m(5) * (t105 * t120 + t90) + m(4) * (t110 * t98 + t90); m(3) + m(7) * (t25 ^ 2 + t27 ^ 2 + t100) + m(6) * (t59 ^ 2 + t61 ^ 2 + t100) + m(5) * (t137 * t98 + t100) - m(4) * t133; (t56 / 0.2e1 + pkin(8) * t76 + t51 * mrSges(5,3)) * t108 + (t57 / 0.2e1 - pkin(8) * t77 - t50 * mrSges(5,3)) * t104 + t88 * t29 + t70 * t42 + t71 * t23 / 0.2e1 + t72 * t24 / 0.2e1 - t60 * t43 / 0.2e1 + t62 * t44 / 0.2e1 - pkin(3) * t65 + t46 * t49 + t52 * t7 + t38 * t15 + t39 * t4 / 0.2e1 + t40 * t5 / 0.2e1 + t47 * t48 + t26 * t16 / 0.2e1 + t28 * t17 / 0.2e1 + t11 * t18 + t10 * t19 + (Ifges(4,5) + t108 * t81 / 0.2e1 - t104 * t80 / 0.2e1 + t138 * t110) * t109 + m(5) * (pkin(3) * t129 + pkin(8) * t120) + m(6) * (t13 * t46 + t14 * t47 + t70 * t88) + m(7) * (t10 * t2 + t11 * t3 + t38 * t52) + (-t13 * t72 + t14 * t71) * mrSges(6,3) + (-t2 * t40 + t3 * t39) * mrSges(7,3) + (t96 / 0.2e1 + t95 / 0.2e1 + t67 / 0.2e1 + t66 / 0.2e1 + t35 / 0.2e1 + t34 / 0.2e1 - Ifges(4,6) - t110 * mrSges(4,2)) * t105; (-t25 * t40 + t27 * t39) * mrSges(7,3) + (t59 * t72 + t61 * t71) * mrSges(6,3) + (-mrSges(4,2) + t126) * t105 + (-t15 - t42 + t138) * t109 + m(7) * (t10 * t25 - t109 * t52 + t11 * t27) + m(6) * (-t109 * t88 - t46 * t59 + t47 * t61) + m(5) * (pkin(8) * t105 * t137 + pkin(3) * t109); -0.2e1 * pkin(3) * t79 + t104 * t81 + t108 * t80 + 0.2e1 * t52 * t15 + t39 * t16 + t40 * t17 + 0.2e1 * t88 * t42 + t71 * t43 + t72 * t44 + Ifges(4,3) + m(7) * (t10 ^ 2 + t11 ^ 2 + t52 ^ 2) + m(6) * (t46 ^ 2 + t47 ^ 2 + t88 ^ 2) + m(5) * (pkin(8) ^ 2 * t137 + pkin(3) ^ 2) + 0.2e1 * (-t10 * t40 + t11 * t39) * mrSges(7,3) + 0.2e1 * (-t46 * t72 + t47 * t71) * mrSges(6,3) + 0.2e1 * pkin(8) * t126; t115 + m(7) * (t2 * t63 + t3 * t64) + (m(6) * (t103 * t14 + t107 * t13) + t103 * t48 + t107 * t49) * pkin(4) - Ifges(5,6) * t132 + t63 * t19 + t64 * t18 + t50 * mrSges(5,1) - t51 * mrSges(5,2) + t145; -t121 * t105 + m(7) * (t25 * t63 + t27 * t64) + m(6) * (t103 * t61 - t107 * t59) * pkin(4) + t117; m(7) * (t10 * t63 + t11 * t64) + t96 + t95 - t121 * pkin(8) + (t39 * t64 - t40 * t63) * mrSges(7,3) + (m(6) * (t103 * t47 + t107 * t46) + (t103 * t71 - t107 * t72) * mrSges(6,3)) * pkin(4) + t116; -0.2e1 * t140 + Ifges(5,3) + 0.2e1 * t58 + 0.2e1 * t118 + m(7) * (t63 ^ 2 + t64 ^ 2) + m(6) * (t103 ^ 2 + t107 ^ 2) * pkin(4) ^ 2 + t139; (m(7) * (t102 * t3 + t106 * t2) + t102 * t18 + t106 * t19) * pkin(5) + t115; m(7) * (t102 * t27 + t106 * t25) * pkin(5) + t117; (m(7) * (t10 * t106 + t102 * t11) + (t102 * t39 - t106 * t40) * mrSges(7,3)) * pkin(5) + t116; Ifges(6,3) + t91 + t118 + (m(7) * (t102 * t64 + t106 * t63) - t134) * pkin(5) + t123; -0.2e1 * t128 + 0.2e1 * t91 + m(7) * (t102 ^ 2 + t106 ^ 2) * pkin(5) ^ 2 + t139; t119; t125; t122; t123; Ifges(7,3) + t91 - t128; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
