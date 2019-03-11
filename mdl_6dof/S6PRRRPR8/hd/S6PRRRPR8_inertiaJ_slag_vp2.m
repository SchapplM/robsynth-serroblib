% Calculate joint inertia matrix for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:48:17
% EndTime: 2019-03-08 23:48:20
% DurationCPUTime: 1.34s
% Computational Cost: add. (1585->336), mult. (3811->464), div. (0->0), fcn. (4000->12), ass. (0->124)
t111 = sin(qJ(4));
t115 = cos(qJ(4));
t163 = t111 ^ 2 + t115 ^ 2;
t110 = sin(qJ(6));
t155 = -t110 / 0.2e1;
t152 = mrSges(6,1) + mrSges(5,3);
t161 = -m(6) * pkin(4) + mrSges(6,2);
t107 = sin(pkin(6));
t109 = cos(pkin(6));
t112 = sin(qJ(3));
t113 = sin(qJ(2));
t116 = cos(qJ(3));
t108 = cos(pkin(7));
t117 = cos(qJ(2));
t137 = t108 * t117;
t106 = sin(pkin(7));
t139 = t106 * t112;
t35 = t109 * t139 + (t112 * t137 + t113 * t116) * t107;
t59 = -t106 * t107 * t117 + t108 * t109;
t16 = t111 * t59 + t115 * t35;
t12 = t16 ^ 2;
t138 = t106 * t116;
t33 = -t109 * t138 + (t112 * t113 - t116 * t137) * t107;
t30 = t33 ^ 2;
t114 = cos(qJ(6));
t60 = -t108 * t115 + t111 * t139;
t38 = t110 * t138 + t114 * t60;
t160 = t38 / 0.2e1;
t143 = Ifges(7,4) * t114;
t53 = Ifges(7,5) * t111 + (-Ifges(7,1) * t110 - t143) * t115;
t159 = t53 / 0.2e1;
t96 = Ifges(7,5) * t114;
t158 = Ifges(7,6) * t155 + t96 / 0.2e1;
t157 = pkin(4) + pkin(11);
t156 = pkin(5) + pkin(10);
t154 = -t114 / 0.2e1;
t153 = pkin(2) * t116;
t75 = mrSges(7,1) * t110 + mrSges(7,2) * t114;
t151 = mrSges(6,3) + t75;
t150 = Ifges(6,1) + Ifges(5,3);
t40 = mrSges(6,1) * t60 + mrSges(6,3) * t138;
t42 = mrSges(5,2) * t138 - mrSges(5,3) * t60;
t149 = -t40 + t42;
t61 = t108 * t111 + t115 * t139;
t41 = t61 * mrSges(6,1) - mrSges(6,2) * t138;
t43 = -mrSges(5,1) * t138 - mrSges(5,3) * t61;
t148 = t41 - t43;
t64 = pkin(2) * t108 * t112 + pkin(9) * t138;
t49 = pkin(10) * t108 + t64;
t50 = (-pkin(3) * t116 - pkin(10) * t112 - pkin(2)) * t106;
t23 = t111 * t50 + t115 * t49;
t147 = Ifges(5,5) * t111 + Ifges(5,6) * t115;
t146 = t163 * pkin(10) ^ 2;
t145 = mrSges(7,3) * t115;
t144 = Ifges(7,4) * t110;
t142 = qJ(5) * t16;
t70 = -mrSges(7,2) * t111 - t114 * t145;
t141 = t110 * t70;
t14 = t111 * t35 - t115 * t59;
t140 = t111 * t14;
t136 = t114 * t157;
t135 = -t110 ^ 2 - t114 ^ 2;
t39 = t110 * t60 - t114 * t138;
t5 = Ifges(7,5) * t39 + Ifges(7,6) * t38 + Ifges(7,3) * t61;
t133 = Ifges(4,5) * t139 + Ifges(4,6) * t138 + Ifges(4,3) * t108;
t132 = t138 / 0.2e1;
t131 = m(7) * t135;
t130 = -qJ(5) * t111 - pkin(3);
t22 = -t111 * t49 + t115 * t50;
t129 = t135 * mrSges(7,3);
t128 = (Ifges(6,4) - Ifges(5,5)) * t61 + (-Ifges(6,5) + Ifges(5,6)) * t60;
t19 = pkin(4) * t138 - t22;
t8 = pkin(5) * t61 + pkin(11) * t138 + t19;
t88 = pkin(9) * t139;
t48 = t88 + (-pkin(3) - t153) * t108;
t121 = -qJ(5) * t61 + t48;
t9 = t157 * t60 + t121;
t1 = -t110 * t9 + t114 * t8;
t2 = t110 * t8 + t114 * t9;
t126 = t1 * t114 + t110 * t2;
t3 = -t110 * t33 + t114 * t14;
t4 = t110 * t14 + t114 * t33;
t125 = t110 * t4 + t114 * t3;
t124 = mrSges(7,1) * t114 - mrSges(7,2) * t110;
t66 = -t115 * t157 + t130;
t84 = t156 * t111;
t31 = -t110 * t66 + t114 * t84;
t32 = t110 * t84 + t114 * t66;
t123 = t110 * t32 + t114 * t31;
t18 = qJ(5) * t138 - t23;
t122 = (t115 * t16 + t140) * pkin(10);
t51 = Ifges(7,3) * t111 + (-Ifges(7,5) * t110 - Ifges(7,6) * t114) * t115;
t119 = qJ(5) ^ 2;
t85 = t156 * t115;
t82 = Ifges(5,1) * t111 + Ifges(5,4) * t115;
t81 = Ifges(7,1) * t114 - t144;
t80 = Ifges(5,4) * t111 + Ifges(5,2) * t115;
t79 = -Ifges(7,2) * t110 + t143;
t77 = -Ifges(6,2) * t111 - Ifges(6,6) * t115;
t76 = -Ifges(6,6) * t111 - Ifges(6,3) * t115;
t74 = -mrSges(5,1) * t115 + mrSges(5,2) * t111;
t73 = mrSges(6,2) * t115 - mrSges(6,3) * t111;
t72 = -pkin(4) * t115 + t130;
t69 = mrSges(7,1) * t111 + t110 * t145;
t68 = -mrSges(4,2) * t108 + mrSges(4,3) * t138;
t67 = mrSges(4,1) * t108 - mrSges(4,3) * t139;
t65 = t124 * t115;
t63 = t108 * t153 - t88;
t62 = (-mrSges(4,1) * t116 + mrSges(4,2) * t112) * t106;
t52 = Ifges(7,6) * t111 + (-Ifges(7,2) * t114 - t144) * t115;
t29 = -mrSges(6,2) * t60 - mrSges(6,3) * t61;
t28 = mrSges(5,1) * t60 + mrSges(5,2) * t61;
t27 = Ifges(5,1) * t61 - Ifges(5,4) * t60 - Ifges(5,5) * t138;
t26 = Ifges(5,4) * t61 - Ifges(5,2) * t60 - Ifges(5,6) * t138;
t25 = -Ifges(6,4) * t138 - Ifges(6,2) * t61 + Ifges(6,6) * t60;
t24 = -Ifges(6,5) * t138 - Ifges(6,6) * t61 + Ifges(6,3) * t60;
t21 = mrSges(7,1) * t61 - mrSges(7,3) * t39;
t20 = -mrSges(7,2) * t61 + mrSges(7,3) * t38;
t13 = pkin(4) * t60 + t121;
t11 = -mrSges(7,1) * t38 + mrSges(7,2) * t39;
t10 = -pkin(5) * t60 - t18;
t7 = Ifges(7,1) * t39 + Ifges(7,4) * t38 + Ifges(7,5) * t61;
t6 = Ifges(7,4) * t39 + Ifges(7,2) * t38 + Ifges(7,6) * t61;
t15 = [m(2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t12) + m(4) * (t35 ^ 2 + t59 ^ 2 + t30) + m(3) * (t109 ^ 2 + (t113 ^ 2 + t117 ^ 2) * t107 ^ 2) + (m(6) + m(5)) * (t14 ^ 2 + t12 + t30); t4 * t20 + t3 * t21 + t35 * t68 + t59 * t62 + t148 * t14 + (mrSges(3,1) * t117 - mrSges(3,2) * t113) * t107 + (t28 + t29 - t67) * t33 + (t11 + t149) * t16 + m(6) * (t13 * t33 + t14 * t19 - t16 * t18) + m(7) * (t1 * t3 + t10 * t16 + t2 * t4) + m(5) * (-t14 * t22 + t16 * t23 + t33 * t48) + m(4) * (-pkin(2) * t106 * t59 - t33 * t63 + t35 * t64); t108 * t133 + 0.2e1 * t63 * t67 + 0.2e1 * t64 * t68 + 0.2e1 * t22 * t43 + 0.2e1 * t48 * t28 + t38 * t6 + t39 * t7 + 0.2e1 * t18 * t40 + 0.2e1 * t19 * t41 + 0.2e1 * t23 * t42 + 0.2e1 * t13 * t29 + 0.2e1 * t2 * t20 + 0.2e1 * t1 * t21 + 0.2e1 * t10 * t11 + Ifges(3,3) + (t24 - t26) * t60 + (-t25 + t27 + t5) * t61 + (-0.2e1 * pkin(2) * t62 + (Ifges(4,1) * t139 + Ifges(4,5) * t108) * t112 + (Ifges(4,6) * t108 + (0.2e1 * Ifges(4,4) * t112 + (Ifges(4,2) + t150) * t116) * t106 + t128) * t116) * t106 + m(4) * (pkin(2) ^ 2 * t106 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2 + t48 ^ 2) + m(6) * (t13 ^ 2 + t18 ^ 2 + t19 ^ 2) + m(7) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2); -t35 * mrSges(4,2) + t3 * t69 + t4 * t70 + t152 * t140 + (-mrSges(4,1) + t73 + t74) * t33 + (t115 * t152 + t65) * t16 + m(6) * (t33 * t72 + t122) + m(7) * (t16 * t85 + t3 * t31 + t32 * t4) + m(5) * (-pkin(3) * t33 + t122); m(6) * (t13 * t72 + (t111 * t19 - t115 * t18) * pkin(10)) + m(5) * (-pkin(3) * t48 + (-t111 * t22 + t115 * t23) * pkin(10)) + m(7) * (t1 * t31 + t10 * t85 + t2 * t32) + (-t77 / 0.2e1 + t82 / 0.2e1 + t51 / 0.2e1) * t61 + (t76 / 0.2e1 - t80 / 0.2e1) * t60 + t85 * t11 + t63 * mrSges(4,1) - t64 * mrSges(4,2) + t10 * t65 + t1 * t69 + t2 * t70 + t72 * t29 + t13 * t73 + t48 * t74 + t52 * t160 + t39 * t159 - pkin(3) * t28 + t31 * t21 + t32 * t20 + (-t24 / 0.2e1 + t26 / 0.2e1 + Ifges(6,5) * t132 + t6 * t154 + t7 * t155 - t18 * mrSges(6,1) + t23 * mrSges(5,3) + t149 * pkin(10)) * t115 + (-t25 / 0.2e1 + t27 / 0.2e1 + t5 / 0.2e1 + Ifges(6,4) * t132 + t19 * mrSges(6,1) - t22 * mrSges(5,3) + t148 * pkin(10)) * t111 - t147 * t138 / 0.2e1 + t133; -0.2e1 * pkin(3) * t74 + 0.2e1 * t31 * t69 + 0.2e1 * t32 * t70 + 0.2e1 * t85 * t65 + 0.2e1 * t72 * t73 + Ifges(4,3) + (t51 - t77 + t82) * t111 + m(7) * (t31 ^ 2 + t32 ^ 2 + t85 ^ 2) + m(6) * (t72 ^ 2 + t146) + m(5) * (pkin(3) ^ 2 + t146) + (-t110 * t53 - t114 * t52 - t76 + t80) * t115 + 0.2e1 * t152 * pkin(10) * t163; (-mrSges(5,1) + mrSges(6,2)) * t14 - t125 * mrSges(7,3) + (-mrSges(5,2) + t151) * t16 + m(6) * (-pkin(4) * t14 + t142) + m(7) * (-t125 * t157 + t142); t10 * t75 + t61 * t158 + t79 * t160 + t39 * t81 / 0.2e1 - pkin(4) * t41 + t22 * mrSges(5,1) - t23 * mrSges(5,2) - t18 * mrSges(6,3) + t19 * mrSges(6,2) - t150 * t138 + (-t40 + t11) * qJ(5) + (t7 / 0.2e1 - t157 * t21 - t1 * mrSges(7,3)) * t114 + (-t6 / 0.2e1 - t157 * t20 - t2 * mrSges(7,3)) * t110 + m(7) * (qJ(5) * t10 - t126 * t157) + m(6) * (-pkin(4) * t19 - qJ(5) * t18) - t128; -t157 * t141 - t69 * t136 + m(7) * (qJ(5) * t85 - t123 * t157) + t114 * t159 + t52 * t155 + t85 * t75 + qJ(5) * t65 - t123 * mrSges(7,3) + (-pkin(4) * mrSges(6,1) - Ifges(6,4) + t158) * t111 + (qJ(5) * mrSges(6,1) + t79 * t154 + t81 * t155 - Ifges(6,5)) * t115 + ((m(6) * qJ(5) - mrSges(5,2) + mrSges(6,3)) * t115 + (-mrSges(5,1) + t161) * t111) * pkin(10) + t147; -0.2e1 * pkin(4) * mrSges(6,2) - t110 * t79 + t114 * t81 + m(6) * (pkin(4) ^ 2 + t119) + m(7) * (-t135 * t157 ^ 2 + t119) + t150 + 0.2e1 * t151 * qJ(5) - 0.2e1 * t157 * t129; m(6) * t14 + m(7) * t125; m(6) * t19 + m(7) * t126 + t110 * t20 + t114 * t21 + t41; m(7) * t123 + t141 + t114 * t69 + (m(6) * pkin(10) + mrSges(6,1)) * t111; t131 * t157 + t129 + t161; m(6) - t131; mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t5; mrSges(7,1) * t31 - mrSges(7,2) * t32 + t51; -mrSges(7,1) * t136 + t96 + (mrSges(7,2) * t157 - Ifges(7,6)) * t110; t124; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
