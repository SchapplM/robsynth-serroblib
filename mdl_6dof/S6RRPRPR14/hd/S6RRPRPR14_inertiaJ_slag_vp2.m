% Calculate joint inertia matrix for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR14_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR14_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR14_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:32:36
% EndTime: 2019-03-09 11:32:40
% DurationCPUTime: 1.36s
% Computational Cost: add. (1530->326), mult. (3246->433), div. (0->0), fcn. (3098->8), ass. (0->127)
t108 = sin(qJ(6));
t111 = cos(qJ(6));
t136 = t108 ^ 2 + t111 ^ 2;
t132 = m(7) * t136;
t167 = m(6) + t132;
t166 = Ifges(4,1) + Ifges(3,3);
t165 = Ifges(6,1) + Ifges(5,3);
t109 = sin(qJ(4));
t101 = t109 ^ 2;
t112 = cos(qJ(4));
t103 = t112 ^ 2;
t135 = t101 + t103;
t164 = (mrSges(6,1) + mrSges(5,3)) * t135;
t163 = -t108 / 0.2e1;
t94 = qJ(5) * t109;
t161 = m(6) * (pkin(4) * t112 + t94);
t160 = t136 * mrSges(7,3) - mrSges(6,2);
t106 = sin(pkin(6));
t110 = sin(qJ(2));
t141 = t106 * t110;
t107 = cos(pkin(6));
t113 = cos(qJ(2));
t140 = t106 * t113;
t50 = t107 * t109 + t112 * t140;
t31 = -t108 * t141 + t111 * t50;
t159 = t31 / 0.2e1;
t32 = t108 * t50 + t111 * t141;
t158 = t32 / 0.2e1;
t96 = Ifges(7,5) * t111;
t157 = Ifges(7,6) * t163 + t96 / 0.2e1;
t115 = -pkin(2) - pkin(9);
t156 = pkin(4) + pkin(10);
t155 = t108 / 0.2e1;
t154 = t111 / 0.2e1;
t153 = pkin(1) * t113;
t82 = pkin(8) * t141;
t52 = t107 * t153 - t82;
t151 = t52 * mrSges(3,1);
t53 = t107 * t110 * pkin(1) + pkin(8) * t140;
t150 = t53 * mrSges(3,2);
t149 = pkin(5) - t115;
t134 = -pkin(2) - t153;
t26 = pkin(3) * t141 + t82 + (-pkin(9) + t134) * t107;
t131 = -qJ(3) * t110 - pkin(1);
t34 = (t115 * t113 + t131) * t106;
t12 = t109 * t26 + t112 * t34;
t35 = mrSges(6,1) * t50 - mrSges(6,3) * t141;
t37 = -mrSges(5,2) * t141 - mrSges(5,3) * t50;
t148 = -t35 + t37;
t51 = t107 * t112 - t109 * t140;
t36 = t51 * mrSges(6,1) + mrSges(6,2) * t141;
t38 = mrSges(5,1) * t141 - mrSges(5,3) * t51;
t147 = -t36 + t38;
t64 = mrSges(7,1) * t108 + mrSges(7,2) * t111;
t146 = t64 + mrSges(6,3);
t145 = t135 * t115 ^ 2;
t57 = mrSges(4,1) * t141 + t107 * mrSges(4,2);
t144 = Ifges(7,4) * t108;
t143 = Ifges(7,4) * t111;
t138 = t109 * t111;
t60 = -mrSges(7,2) * t112 + mrSges(7,3) * t138;
t142 = t108 * t60;
t139 = t108 * t109;
t137 = t111 * t156;
t6 = Ifges(7,5) * t32 + Ifges(7,6) * t31 + Ifges(7,3) * t51;
t39 = -t107 * qJ(3) - t53;
t42 = Ifges(7,5) * t139 + Ifges(7,6) * t138 + Ifges(7,3) * t112;
t133 = -t141 / 0.2e1;
t11 = -t109 * t34 + t112 * t26;
t130 = Ifges(3,5) * t141 + Ifges(3,6) * t140 + t166 * t107;
t33 = pkin(3) * t140 - t39;
t129 = t136 * t156;
t61 = t109 * pkin(4) - qJ(5) * t112 + qJ(3);
t3 = pkin(5) * t51 - t156 * t141 - t11;
t119 = -qJ(5) * t51 + t33;
t5 = t156 * t50 + t119;
t1 = -t108 * t5 + t111 * t3;
t2 = t108 * t3 + t111 * t5;
t127 = t1 * t111 + t108 * t2;
t10 = -pkin(4) * t141 - t11;
t9 = -qJ(5) * t141 - t12;
t126 = -t10 * t112 - t9 * t109;
t125 = mrSges(7,1) * t111 - t108 * mrSges(7,2);
t15 = -mrSges(7,2) * t51 + mrSges(7,3) * t31;
t16 = mrSges(7,1) * t51 - mrSges(7,3) * t32;
t124 = t108 * t15 + t111 * t16;
t54 = pkin(10) * t109 + t61;
t63 = t149 * t112;
t24 = -t108 * t54 + t111 * t63;
t25 = t108 * t63 + t111 * t54;
t123 = t108 * t25 + t111 * t24;
t59 = mrSges(7,1) * t112 - mrSges(7,3) * t139;
t122 = -t111 * t59 - t142;
t121 = t12 * t109 + t11 * t112;
t120 = (-Ifges(6,4) + Ifges(5,5)) * t51 + (Ifges(6,5) - Ifges(5,6)) * t50 + t165 * t141;
t118 = 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t135;
t117 = qJ(3) ^ 2;
t116 = qJ(5) ^ 2;
t98 = Ifges(5,5) * t112;
t97 = Ifges(6,5) * t109;
t73 = Ifges(5,1) * t112 - Ifges(5,4) * t109;
t72 = Ifges(7,1) * t111 - t144;
t71 = Ifges(5,4) * t112 - Ifges(5,2) * t109;
t70 = -Ifges(7,2) * t108 + t143;
t68 = -Ifges(6,2) * t112 + Ifges(6,6) * t109;
t67 = -Ifges(6,6) * t112 + Ifges(6,3) * t109;
t66 = mrSges(5,1) * t109 + mrSges(5,2) * t112;
t65 = -mrSges(6,2) * t109 - mrSges(6,3) * t112;
t62 = t149 * t109;
t56 = -mrSges(4,1) * t140 - mrSges(4,3) * t107;
t55 = t125 * t109;
t44 = Ifges(7,5) * t112 + (Ifges(7,1) * t108 + t143) * t109;
t43 = Ifges(7,6) * t112 + (Ifges(7,2) * t111 + t144) * t109;
t41 = t134 * t107 + t82;
t40 = (-pkin(2) * t113 + t131) * t106;
t22 = -mrSges(6,2) * t50 - mrSges(6,3) * t51;
t21 = mrSges(5,1) * t50 + mrSges(5,2) * t51;
t20 = Ifges(5,1) * t51 - Ifges(5,4) * t50 + Ifges(5,5) * t141;
t19 = Ifges(5,4) * t51 - Ifges(5,2) * t50 + Ifges(5,6) * t141;
t18 = Ifges(6,4) * t141 - Ifges(6,2) * t51 + Ifges(6,6) * t50;
t17 = Ifges(6,5) * t141 - Ifges(6,6) * t51 + Ifges(6,3) * t50;
t14 = -mrSges(7,1) * t31 + mrSges(7,2) * t32;
t13 = pkin(4) * t50 + t119;
t8 = Ifges(7,1) * t32 + Ifges(7,4) * t31 + Ifges(7,5) * t51;
t7 = Ifges(7,4) * t32 + Ifges(7,2) * t31 + Ifges(7,6) * t51;
t4 = -pkin(5) * t50 - t9;
t23 = [0.2e1 * t1 * t16 + 0.2e1 * t10 * t36 + 0.2e1 * t11 * t38 + 0.2e1 * t12 * t37 + 0.2e1 * t13 * t22 + 0.2e1 * t4 * t14 + 0.2e1 * t2 * t15 + 0.2e1 * t33 * t21 + t31 * t7 + t32 * t8 + 0.2e1 * t9 * t35 + 0.2e1 * t39 * t56 + 0.2e1 * t41 * t57 + Ifges(2,3) + (t17 - t19) * t50 + (t20 + t6 - t18) * t51 + (t130 - 0.2e1 * t150 + 0.2e1 * t151) * t107 + m(3) * (t52 ^ 2 + t53 ^ 2) + m(5) * (t11 ^ 2 + t12 ^ 2 + t33 ^ 2) + m(4) * (t39 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t4 ^ 2) + m(6) * (t10 ^ 2 + t13 ^ 2 + t9 ^ 2) + (m(3) * pkin(1) ^ 2 * t106 + (0.2e1 * t40 * mrSges(4,2) + 0.2e1 * t53 * mrSges(3,3) + (-(2 * Ifges(4,5)) + Ifges(3,6)) * t107 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + Ifges(4,3)) * t113) * t106) * t113 + (-0.2e1 * t52 * mrSges(3,3) - 0.2e1 * t40 * mrSges(4,3) + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t107 + (-0.2e1 * pkin(1) * mrSges(3,2) + (Ifges(4,2) + Ifges(3,1)) * t110 + 0.2e1 * (Ifges(3,4) + Ifges(4,6)) * t113) * t106 + t120) * t110) * t106; m(7) * (t1 * t24 + t2 * t25 - t4 * t62) + (-Ifges(4,5) * t113 + (-Ifges(4,4) + t98 / 0.2e1 + t97 / 0.2e1) * t110) * t106 + t44 * t158 + t43 * t159 + (-t56 + t21) * qJ(3) + m(4) * (-pkin(2) * t41 - qJ(3) * t39) + (-t68 / 0.2e1 + t73 / 0.2e1 + t42 / 0.2e1) * t51 + (-t71 / 0.2e1 + t67 / 0.2e1) * t50 - t150 + t151 + (t17 / 0.2e1 - t19 / 0.2e1 + t7 * t154 + t8 * t155 + t9 * mrSges(6,1) - t12 * mrSges(5,3) + Ifges(5,6) * t133 + t148 * t115) * t109 + (-t18 / 0.2e1 + t20 / 0.2e1 + t6 / 0.2e1 + Ifges(6,4) * t133 + t10 * mrSges(6,1) - t11 * mrSges(5,3) + t147 * t115) * t112 + m(6) * (t126 * t115 + t13 * t61) + m(5) * (qJ(3) * t33 + t121 * t115) - t4 * t55 - pkin(2) * t57 + t1 * t59 + t2 * t60 + t61 * t22 - t62 * t14 + t13 * t65 + t33 * t66 - t39 * mrSges(4,3) + t41 * mrSges(4,2) + t24 * t16 + t25 * t15 + t130; -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t24 * t59 + 0.2e1 * t25 * t60 + 0.2e1 * t62 * t55 + 0.2e1 * t61 * t65 + 0.2e1 * (t66 + mrSges(4,3)) * qJ(3) + (t42 + t73 - t68) * t112 + (t108 * t44 + t111 * t43 + t67 - t71) * t109 + m(7) * (t24 ^ 2 + t25 ^ 2 + t62 ^ 2) + m(6) * (t61 ^ 2 + t145) + m(5) * (t117 + t145) + m(4) * (pkin(2) ^ 2 + t117) - 0.2e1 * t115 * t164 + t166; (t14 + t148) * t109 + (-t124 + t147) * t112 + m(7) * (t109 * t4 - t127 * t112) + m(6) * t126 + m(5) * t121 + m(4) * t41 + t57; -m(4) * pkin(2) - t109 * t55 + mrSges(4,2) + t122 * t112 + m(7) * (-t109 * t62 - t123 * t112) + t115 * t118 - t164; m(4) + m(7) * (t136 * t103 + t101) + t118; t120 + (t8 / 0.2e1 - t156 * t16 - t1 * mrSges(7,3)) * t111 + (-t7 / 0.2e1 - t156 * t15 - t2 * mrSges(7,3)) * t108 + (-t35 + t14) * qJ(5) + m(6) * (-pkin(4) * t10 - qJ(5) * t9) + t51 * t157 + t70 * t159 + t72 * t158 + t4 * t64 - pkin(4) * t36 - t9 * mrSges(6,3) + t10 * mrSges(6,2) + t11 * mrSges(5,1) - t12 * mrSges(5,2) + m(7) * (qJ(5) * t4 - t127 * t156); -t59 * t137 - t156 * t142 + t98 + t97 + m(7) * (-qJ(5) * t62 - t123 * t156) - t62 * t64 + t44 * t154 + t43 * t163 - qJ(5) * t55 - t123 * mrSges(7,3) + (-pkin(4) * mrSges(6,1) - Ifges(6,4) + t157) * t112 + (-qJ(5) * mrSges(6,1) + t70 * t154 + t72 * t155 - Ifges(5,6)) * t109 + (t161 + (mrSges(5,1) - mrSges(6,2)) * t112 + (-mrSges(5,2) + mrSges(6,3)) * t109) * t115; (mrSges(5,1) + t160) * t112 + t161 + m(7) * (t112 * t129 + t94) + (-mrSges(5,2) + t146) * t109; -0.2e1 * pkin(4) * mrSges(6,2) - t108 * t70 + t111 * t72 + m(7) * (t136 * t156 ^ 2 + t116) + m(6) * (pkin(4) ^ 2 + t116) + 0.2e1 * mrSges(7,3) * t129 + 0.2e1 * t146 * qJ(5) + t165; m(6) * t10 + m(7) * t127 + t124 + t36; m(7) * t123 + (-m(6) * t115 + mrSges(6,1)) * t112 - t122; -t167 * t112; -m(6) * pkin(4) - t132 * t156 - t160; t167; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t6; mrSges(7,1) * t24 - mrSges(7,2) * t25 + t42; -t125 * t112; -mrSges(7,1) * t137 + t96 + (mrSges(7,2) * t156 - Ifges(7,6)) * t108; t125; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t23(1) t23(2) t23(4) t23(7) t23(11) t23(16); t23(2) t23(3) t23(5) t23(8) t23(12) t23(17); t23(4) t23(5) t23(6) t23(9) t23(13) t23(18); t23(7) t23(8) t23(9) t23(10) t23(14) t23(19); t23(11) t23(12) t23(13) t23(14) t23(15) t23(20); t23(16) t23(17) t23(18) t23(19) t23(20) t23(21);];
Mq  = res;
