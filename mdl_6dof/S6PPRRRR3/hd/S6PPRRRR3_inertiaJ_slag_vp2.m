% Calculate joint inertia matrix for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2018-11-23 14:53
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PPRRRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_inertiaJ_slag_vp2: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:52:57
% EndTime: 2018-11-23 14:52:58
% DurationCPUTime: 1.07s
% Computational Cost: add. (1957->313), mult. (5333->475), div. (0->0), fcn. (6227->16), ass. (0->126)
t160 = 2 * pkin(11);
t159 = m(7) * pkin(12) + mrSges(7,3);
t108 = sin(qJ(6));
t112 = cos(qJ(6));
t73 = -mrSges(7,1) * t112 + mrSges(7,2) * t108;
t158 = -m(7) * pkin(5) - mrSges(6,1) + t73;
t101 = sin(pkin(8));
t105 = cos(pkin(8));
t110 = sin(qJ(4));
t114 = cos(qJ(4));
t100 = sin(pkin(14));
t103 = sin(pkin(6));
t107 = cos(pkin(6));
t111 = sin(qJ(3));
t115 = cos(qJ(3));
t104 = cos(pkin(14));
t106 = cos(pkin(7));
t130 = t104 * t106;
t102 = sin(pkin(7));
t131 = t102 * t115;
t33 = t107 * t131 + (-t100 * t111 + t115 * t130) * t103;
t132 = t102 * t111;
t34 = t107 * t132 + (t100 * t115 + t111 * t130) * t103;
t59 = -t102 * t103 * t104 + t106 * t107;
t10 = t114 * t34 + (t101 * t59 + t105 * t33) * t110;
t109 = sin(qJ(5));
t113 = cos(qJ(5));
t18 = -t101 * t33 + t105 * t59;
t3 = t10 * t109 - t113 * t18;
t157 = t3 ^ 2;
t129 = t105 * t114;
t133 = t101 * t114;
t8 = t110 * t34 - t33 * t129 - t59 * t133;
t156 = t8 ^ 2;
t128 = t105 * t115;
t134 = t101 * t110;
t37 = t106 * t134 + (t110 * t128 + t111 * t114) * t102;
t60 = -t101 * t131 + t105 * t106;
t19 = t109 * t37 - t113 * t60;
t155 = t19 ^ 2;
t35 = -t114 * t102 * t128 - t106 * t133 + t110 * t132;
t154 = t35 ^ 2;
t62 = t105 * t109 + t113 * t134;
t40 = -t108 * t62 - t112 * t133;
t153 = t40 / 0.2e1;
t41 = -t108 * t133 + t112 * t62;
t152 = t41 / 0.2e1;
t151 = t19 * t3;
t150 = t35 * t8;
t149 = -t108 / 0.2e1;
t148 = t108 / 0.2e1;
t147 = t112 / 0.2e1;
t146 = pkin(3) * t101;
t145 = pkin(11) * t113;
t144 = t109 * t3;
t74 = -mrSges(6,1) * t113 + mrSges(6,2) * t109;
t143 = -mrSges(5,1) + t74;
t16 = -mrSges(7,1) * t40 + mrSges(7,2) * t41;
t43 = -mrSges(6,1) * t133 - mrSges(6,3) * t62;
t142 = t16 - t43;
t61 = -t113 * t105 + t109 * t134;
t32 = mrSges(6,1) * t61 + mrSges(6,2) * t62;
t67 = mrSges(5,1) * t105 - mrSges(5,3) * t134;
t141 = t32 - t67;
t65 = t105 * t110 * pkin(3) + pkin(10) * t133;
t51 = pkin(11) * t105 + t65;
t52 = (-pkin(4) * t114 - pkin(11) * t110 - pkin(3)) * t101;
t28 = t109 * t52 + t113 * t51;
t140 = -Ifges(6,5) * t62 + Ifges(6,6) * t61;
t75 = Ifges(7,5) * t108 + Ifges(7,6) * t112;
t139 = Ifges(6,5) * t109 + Ifges(6,6) * t113;
t138 = t108 ^ 2 + t112 ^ 2;
t137 = Ifges(7,4) * t108;
t136 = Ifges(7,4) * t112;
t135 = t109 * t19;
t127 = t108 * t109;
t126 = t109 * t112;
t13 = Ifges(7,5) * t41 + Ifges(7,6) * t40 + Ifges(7,3) * t61;
t125 = Ifges(5,5) * t134 + Ifges(5,6) * t133 + Ifges(5,3) * t105;
t83 = pkin(10) * t134;
t50 = t83 + (-pkin(3) * t114 - pkin(4)) * t105;
t22 = pkin(5) * t61 - pkin(12) * t62 + t50;
t24 = -pkin(12) * t133 + t28;
t6 = -t108 * t24 + t112 * t22;
t7 = t108 * t22 + t112 * t24;
t123 = -t108 * t6 + t112 * t7;
t5 = t10 * t113 + t109 * t18;
t122 = t113 * t5 + t144;
t121 = mrSges(7,1) * t108 + mrSges(7,2) * t112;
t72 = -pkin(5) * t113 - pkin(12) * t109 - pkin(4);
t47 = -t108 * t145 + t112 * t72;
t48 = t108 * t72 + t112 * t145;
t119 = -t108 * t47 + t112 * t48;
t21 = t109 * t60 + t113 * t37;
t118 = t113 * t21 + t135;
t27 = -t109 * t51 + t113 * t52;
t53 = Ifges(7,5) * t126 - Ifges(7,6) * t127 - Ifges(7,3) * t113;
t117 = pkin(11) ^ 2;
t99 = t113 ^ 2;
t97 = t109 ^ 2;
t93 = t97 * t117;
t79 = Ifges(6,1) * t109 + Ifges(6,4) * t113;
t78 = Ifges(7,1) * t108 + t136;
t77 = Ifges(6,4) * t109 + Ifges(6,2) * t113;
t76 = Ifges(7,2) * t112 + t137;
t70 = -mrSges(7,1) * t113 - mrSges(7,3) * t126;
t69 = mrSges(7,2) * t113 - mrSges(7,3) * t127;
t68 = -mrSges(5,2) * t105 + mrSges(5,3) * t133;
t66 = t121 * t109;
t64 = pkin(3) * t129 - t83;
t63 = (-mrSges(5,1) * t114 + mrSges(5,2) * t110) * t101;
t55 = -Ifges(7,5) * t113 + (Ifges(7,1) * t112 - t137) * t109;
t54 = -Ifges(7,6) * t113 + (-Ifges(7,2) * t108 + t136) * t109;
t42 = mrSges(6,2) * t133 - mrSges(6,3) * t61;
t30 = Ifges(6,1) * t62 - Ifges(6,4) * t61 - Ifges(6,5) * t133;
t29 = Ifges(6,4) * t62 - Ifges(6,2) * t61 - Ifges(6,6) * t133;
t26 = mrSges(7,1) * t61 - mrSges(7,3) * t41;
t25 = -mrSges(7,2) * t61 + mrSges(7,3) * t40;
t23 = pkin(5) * t133 - t27;
t15 = Ifges(7,1) * t41 + Ifges(7,4) * t40 + Ifges(7,5) * t61;
t14 = Ifges(7,4) * t41 + Ifges(7,2) * t40 + Ifges(7,6) * t61;
t12 = t108 * t35 + t112 * t21;
t11 = -t108 * t21 + t112 * t35;
t2 = t108 * t8 + t112 * t5;
t1 = -t108 * t5 + t112 * t8;
t4 = [m(2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t157) + m(6) * (t5 ^ 2 + t156 + t157) + m(5) * (t10 ^ 2 + t18 ^ 2 + t156) + m(4) * (t33 ^ 2 + t34 ^ 2 + t59 ^ 2) + m(3) * (t107 ^ 2 + (t100 ^ 2 + t104 ^ 2) * t103 ^ 2); m(3) * t107 + m(7) * (t1 * t11 + t12 * t2 + t151) + m(6) * (t21 * t5 + t150 + t151) + m(5) * (t10 * t37 + t18 * t60 + t150) + m(4) * (t106 * t59 + (t111 * t34 + t115 * t33) * t102); m(3) + m(7) * (t11 ^ 2 + t12 ^ 2 + t155) + m(6) * (t21 ^ 2 + t154 + t155) + m(5) * (t37 ^ 2 + t60 ^ 2 + t154) + m(4) * (t106 ^ 2 + (t111 ^ 2 + t115 ^ 2) * t102 ^ 2); t33 * mrSges(4,1) - t34 * mrSges(4,2) + t1 * t26 + t10 * t68 + t18 * t63 + t2 * t25 + t5 * t42 + t141 * t8 + t142 * t3 + m(7) * (t1 * t6 + t2 * t7 + t23 * t3) + m(6) * (-t27 * t3 + t28 * t5 + t50 * t8) + m(5) * (t10 * t65 - t18 * t146 - t64 * t8); t11 * t26 + t12 * t25 + t21 * t42 + t37 * t68 + t60 * t63 + t141 * t35 + t142 * t19 + (mrSges(4,1) * t115 - mrSges(4,2) * t111) * t102 + m(7) * (t11 * t6 + t12 * t7 + t19 * t23) + m(6) * (-t19 * t27 + t21 * t28 + t35 * t50) + m(5) * (-t60 * t146 - t35 * t64 + t37 * t65); Ifges(4,3) + 0.2e1 * t6 * t26 + 0.2e1 * t7 * t25 + t41 * t15 + t40 * t14 + 0.2e1 * t23 * t16 + 0.2e1 * t50 * t32 + t62 * t30 + 0.2e1 * t27 * t43 + 0.2e1 * t28 * t42 + t105 * t125 + 0.2e1 * t64 * t67 + 0.2e1 * t65 * t68 + (t13 - t29) * t61 + (-0.2e1 * pkin(3) * t63 + (Ifges(5,1) * t134 + Ifges(5,5) * t105) * t110 + (Ifges(5,6) * t105 + (0.2e1 * Ifges(5,4) * t110 + (Ifges(6,3) + Ifges(5,2)) * t114) * t101 + t140) * t114) * t101 + m(7) * (t23 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(6) * (t27 ^ 2 + t28 ^ 2 + t50 ^ 2) + m(5) * (pkin(3) ^ 2 * t101 ^ 2 + t64 ^ 2 + t65 ^ 2); -t10 * mrSges(5,2) + t1 * t70 + t2 * t69 + t3 * t66 + t143 * t8 + t122 * mrSges(6,3) + m(7) * (pkin(11) * t144 + t1 * t47 + t2 * t48) + m(6) * (-pkin(4) * t8 + t122 * pkin(11)); -t37 * mrSges(5,2) + t11 * t70 + t12 * t69 + t19 * t66 + t143 * t35 + t118 * mrSges(6,3) + m(7) * (pkin(11) * t135 + t11 * t47 + t12 * t48) + m(6) * (-pkin(4) * t35 + t118 * pkin(11)); m(6) * (-pkin(4) * t50 + (-t27 * t109 + t28 * t113) * pkin(11)) - t139 * t133 / 0.2e1 + t23 * t66 + t7 * t69 + t6 * t70 + t50 * t74 + t62 * t79 / 0.2e1 + t64 * mrSges(5,1) - t65 * mrSges(5,2) + t125 + t47 * t26 + t48 * t25 + t54 * t153 + t55 * t152 + (t30 / 0.2e1 + t15 * t147 + t14 * t149 - t27 * mrSges(6,3) + t142 * pkin(11)) * t109 - pkin(4) * t32 + (t29 / 0.2e1 - t13 / 0.2e1 + pkin(11) * t42 + t28 * mrSges(6,3)) * t113 + (-t77 / 0.2e1 + t53 / 0.2e1) * t61 + m(7) * (pkin(11) * t109 * t23 + t47 * t6 + t48 * t7); -0.2e1 * pkin(4) * t74 + 0.2e1 * t47 * t70 + 0.2e1 * t48 * t69 + Ifges(5,3) + (-t53 + t77) * t113 + (t97 + t99) * mrSges(6,3) * t160 + m(7) * (t47 ^ 2 + t48 ^ 2 + t93) + m(6) * (pkin(4) ^ 2 + t117 * t99 + t93) + (-t108 * t54 + t112 * t55 + t160 * t66 + t79) * t109; -t5 * mrSges(6,2) + t159 * (-t1 * t108 + t112 * t2) + t158 * t3; -t21 * mrSges(6,2) + t159 * (-t108 * t11 + t112 * t12) + t158 * t19; t61 * t75 / 0.2e1 + t78 * t152 + t76 * t153 + t23 * t73 + t15 * t148 + t14 * t147 - Ifges(6,3) * t133 + t27 * mrSges(6,1) - t28 * mrSges(6,2) + t123 * mrSges(7,3) - t140 + (-m(7) * t23 - t16) * pkin(5) + (m(7) * t123 - t108 * t26 + t112 * t25) * pkin(12); t54 * t147 + t55 * t148 - pkin(5) * t66 + (m(7) * t119 - t108 * t70 + t112 * t69) * pkin(12) + (pkin(11) * t158 + t78 * t147 + t76 * t149) * t109 + (-t75 / 0.2e1 - pkin(11) * mrSges(6,2)) * t113 + t119 * mrSges(7,3) + t139; Ifges(6,3) - 0.2e1 * pkin(5) * t73 + t108 * t78 + t112 * t76 + m(7) * (t138 * pkin(12) ^ 2 + pkin(5) ^ 2) + 0.2e1 * t138 * pkin(12) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2; mrSges(7,1) * t11 - mrSges(7,2) * t12; mrSges(7,1) * t6 - mrSges(7,2) * t7 + t13; mrSges(7,1) * t47 - mrSges(7,2) * t48 + t53; -t121 * pkin(12) + t75; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t4(1) t4(2) t4(4) t4(7) t4(11) t4(16); t4(2) t4(3) t4(5) t4(8) t4(12) t4(17); t4(4) t4(5) t4(6) t4(9) t4(13) t4(18); t4(7) t4(8) t4(9) t4(10) t4(14) t4(19); t4(11) t4(12) t4(13) t4(14) t4(15) t4(20); t4(16) t4(17) t4(18) t4(19) t4(20) t4(21);];
Mq  = res;
