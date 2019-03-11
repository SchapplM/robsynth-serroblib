% Calculate time derivative of joint inertia matrix for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:20:20
% EndTime: 2019-03-09 02:20:24
% DurationCPUTime: 2.19s
% Computational Cost: add. (4530->338), mult. (9740->507), div. (0->0), fcn. (9475->10), ass. (0->139)
t109 = sin(qJ(5));
t138 = qJD(5) * t109;
t105 = sin(pkin(11));
t106 = cos(pkin(11));
t110 = sin(qJ(4));
t113 = cos(qJ(4));
t84 = t105 * t113 + t110 * t106;
t132 = t84 * t138;
t112 = cos(qJ(5));
t83 = t105 * t110 - t113 * t106;
t78 = t83 * qJD(4);
t140 = t112 * t78;
t114 = t132 + t140;
t137 = qJD(5) * t112;
t115 = -t109 * t78 + t137 * t84;
t30 = t115 * mrSges(6,1) - t114 * mrSges(6,2);
t108 = sin(qJ(6));
t111 = cos(qJ(6));
t118 = t108 * t109 - t111 * t112;
t163 = qJD(5) + qJD(6);
t86 = t108 * t112 + t109 * t111;
t66 = t163 * t86;
t18 = t118 * t78 - t66 * t84;
t55 = t118 * t84;
t19 = t163 * t55 + t86 * t78;
t6 = -t19 * mrSges(7,1) + t18 * mrSges(7,2);
t167 = t30 + t6;
t65 = t163 * t118;
t31 = t66 * mrSges(7,1) - t65 * mrSges(7,2);
t124 = mrSges(6,1) * t109 + mrSges(6,2) * t112;
t87 = t124 * qJD(5);
t166 = t31 + t87;
t79 = t84 * qJD(4);
t165 = -Ifges(6,5) * t140 + Ifges(6,3) * t79;
t97 = sin(pkin(10)) * pkin(1) + qJ(3);
t152 = pkin(7) + t97;
t80 = t152 * t105;
t81 = t152 * t106;
t164 = -t110 * t81 - t113 * t80;
t162 = 2 * m(7);
t161 = -2 * mrSges(5,3);
t57 = -t110 * t80 + t113 * t81;
t43 = qJD(3) * t84 + qJD(4) * t57;
t160 = 0.2e1 * t43;
t159 = -0.2e1 * t164;
t158 = m(7) * pkin(5);
t156 = -t84 / 0.2e1;
t145 = Ifges(6,4) * t109;
t93 = Ifges(6,2) * t112 + t145;
t155 = -t93 / 0.2e1;
t154 = -pkin(9) - pkin(8);
t153 = pkin(4) * t79;
t151 = t43 * t164;
t150 = t79 * Ifges(6,5);
t149 = t79 * Ifges(6,6);
t148 = t83 * Ifges(6,6);
t67 = t83 * t79;
t50 = t112 * t57;
t117 = -cos(pkin(10)) * pkin(1) - pkin(3) * t106 - pkin(2);
t51 = pkin(4) * t83 - pkin(8) * t84 + t117;
t28 = t109 * t51 + t50;
t147 = -Ifges(7,5) * t65 - Ifges(7,6) * t66;
t92 = -mrSges(6,1) * t112 + mrSges(6,2) * t109;
t146 = t92 - mrSges(5,1);
t144 = Ifges(6,4) * t112;
t143 = Ifges(6,6) * t109;
t142 = t109 * t84;
t139 = t112 * t84;
t136 = qJD(6) * t108;
t135 = qJD(6) * t111;
t134 = Ifges(7,5) * t18 + Ifges(7,6) * t19 + Ifges(7,3) * t79;
t133 = pkin(5) * t138;
t131 = qJD(5) * t154;
t75 = t78 * mrSges(5,2);
t130 = t79 * mrSges(5,1) - t75;
t129 = -(2 * Ifges(5,4)) - t143;
t42 = -t83 * qJD(3) + t164 * qJD(4);
t59 = pkin(8) * t78 + t153;
t128 = -t109 * t42 + t112 * t59;
t27 = -t109 * t57 + t112 * t51;
t127 = (t109 ^ 2 + t112 ^ 2) * t78;
t125 = -t164 * t79 + t43 * t83;
t123 = Ifges(6,1) * t112 - t145;
t122 = -Ifges(6,2) * t109 + t144;
t20 = pkin(5) * t83 - pkin(9) * t139 + t27;
t21 = -pkin(9) * t142 + t28;
t8 = -t108 * t21 + t111 * t20;
t9 = t108 * t20 + t111 * t21;
t95 = t154 * t109;
t96 = t154 * t112;
t71 = t108 * t96 + t111 * t95;
t72 = t108 * t95 - t111 * t96;
t121 = -t109 * t27 + t112 * t28;
t60 = -mrSges(6,2) * t83 - mrSges(6,3) * t142;
t61 = mrSges(6,1) * t83 - mrSges(6,3) * t139;
t120 = -t109 * t61 + t112 * t60;
t90 = t109 * t131;
t91 = t112 * t131;
t36 = qJD(6) * t71 + t108 * t91 + t111 * t90;
t37 = -qJD(6) * t72 - t108 * t90 + t111 * t91;
t119 = t37 * mrSges(7,1) - t36 * mrSges(7,2) + t147;
t11 = t109 * t59 + t112 * t42 + t51 * t137 - t138 * t57;
t10 = -pkin(9) * t115 + t11;
t7 = pkin(9) * t140 + pkin(5) * t79 + (-t50 + (pkin(9) * t84 - t51) * t109) * qJD(5) + t128;
t2 = qJD(6) * t8 + t10 * t111 + t108 * t7;
t3 = -qJD(6) * t9 - t10 * t108 + t111 * t7;
t116 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t134;
t100 = Ifges(6,5) * t137;
t99 = -pkin(5) * t112 - pkin(4);
t94 = Ifges(6,1) * t109 + t144;
t89 = t123 * qJD(5);
t88 = t122 * qJD(5);
t82 = (-mrSges(7,1) * t108 - mrSges(7,2) * t111) * qJD(6) * pkin(5);
t70 = Ifges(7,1) * t86 - Ifges(7,4) * t118;
t69 = Ifges(7,4) * t86 - Ifges(7,2) * t118;
t68 = mrSges(7,1) * t118 + mrSges(7,2) * t86;
t58 = t124 * t84;
t54 = t86 * t84;
t46 = t83 * Ifges(6,5) + t123 * t84;
t45 = t122 * t84 + t148;
t44 = pkin(5) * t142 - t164;
t41 = mrSges(7,1) * t83 + mrSges(7,3) * t55;
t40 = -mrSges(7,2) * t83 - mrSges(7,3) * t54;
t39 = -mrSges(6,2) * t79 - mrSges(6,3) * t115;
t38 = mrSges(6,1) * t79 + mrSges(6,3) * t114;
t33 = -Ifges(7,1) * t65 - Ifges(7,4) * t66;
t32 = -Ifges(7,4) * t65 - Ifges(7,2) * t66;
t29 = mrSges(7,1) * t54 - mrSges(7,2) * t55;
t26 = pkin(5) * t115 + t43;
t25 = -Ifges(7,1) * t55 - Ifges(7,4) * t54 + Ifges(7,5) * t83;
t24 = -Ifges(7,4) * t55 - Ifges(7,2) * t54 + Ifges(7,6) * t83;
t23 = -Ifges(6,1) * t114 - Ifges(6,4) * t115 + t150;
t22 = -Ifges(6,4) * t114 - Ifges(6,2) * t115 + t149;
t14 = -mrSges(7,2) * t79 + mrSges(7,3) * t19;
t13 = mrSges(7,1) * t79 - mrSges(7,3) * t18;
t12 = -qJD(5) * t28 + t128;
t5 = Ifges(7,1) * t18 + Ifges(7,4) * t19 + t79 * Ifges(7,5);
t4 = Ifges(7,4) * t18 + Ifges(7,2) * t19 + t79 * Ifges(7,6);
t1 = [0.2e1 * (m(4) * t97 + mrSges(4,3)) * qJD(3) * (t105 ^ 2 + t106 ^ 2) + (t57 * t161 + (Ifges(6,5) * t112 + t129) * t84 + ((2 * Ifges(5,2)) + Ifges(6,3) + Ifges(7,3)) * t83 - Ifges(7,5) * t55 - Ifges(7,6) * t54) * t79 + (mrSges(5,3) * t160 - 0.2e1 * Ifges(5,1) * t78 - t109 * t22 + t112 * t23 + (-t112 * t45 - t109 * t46 + t83 * (-Ifges(6,5) * t109 - Ifges(6,6) * t112)) * qJD(5)) * t84 - (mrSges(5,3) * t159 - t109 * t45 + t112 * t46) * t78 + 0.2e1 * t117 * t130 + (-t129 * t78 + t42 * t161 + t134 + t165) * t83 + t30 * t159 + t58 * t160 + (t2 * t9 + t26 * t44 + t3 * t8) * t162 + 0.2e1 * m(5) * (t42 * t57 - t151) + 0.2e1 * m(6) * (t11 * t28 + t12 * t27 - t151) + 0.2e1 * t11 * t60 + 0.2e1 * t12 * t61 - t54 * t4 - t55 * t5 + 0.2e1 * t28 * t39 + 0.2e1 * t2 * t40 + 0.2e1 * t3 * t41 + 0.2e1 * t44 * t6 + 0.2e1 * t27 * t38 + t19 * t24 + t18 * t25 + 0.2e1 * t26 * t29 + 0.2e1 * t8 * t13 + 0.2e1 * t9 * t14; -t54 * t13 - t55 * t14 + t18 * t40 + t19 * t41 + t167 * t83 + (t29 + t58) * t79 - t120 * t78 + (-t109 * t38 + t112 * t39 + (-t109 * t60 - t112 * t61) * qJD(5)) * t84 + m(7) * (t18 * t9 + t19 * t8 - t2 * t55 + t26 * t83 - t3 * t54 + t44 * t79) + m(6) * (-t121 * t78 + (-t109 * t12 + t11 * t112 + (-t109 * t28 - t112 * t27) * qJD(5)) * t84 + t125) + m(5) * (t42 * t84 - t57 * t78 + t125); 0.2e1 * m(7) * (-t18 * t55 - t19 * t54 + t67) + 0.2e1 * m(6) * (-t127 * t84 + t67) + 0.2e1 * m(5) * (-t78 * t84 + t67); t109 * t39 + t112 * t38 - t118 * t13 + t86 * t14 - t65 * t40 - t66 * t41 + t120 * qJD(5) + m(7) * (-t118 * t3 + t2 * t86 - t65 * t9 - t66 * t8) + m(6) * (qJD(5) * t121 + t109 * t11 + t112 * t12) + t130; m(7) * (-t118 * t19 + t18 * t86 + t54 * t66 + t55 * t65); (t118 * t66 - t65 * t86) * t162; (t149 / 0.2e1 + t22 / 0.2e1 + t11 * mrSges(6,3) + t84 * t89 / 0.2e1 - t78 * t94 / 0.2e1 + (t84 * t155 - t27 * mrSges(6,3) + t46 / 0.2e1) * qJD(5) + (-qJD(5) * t61 + m(6) * (-qJD(5) * t27 + t11) + t39) * pkin(8)) * t112 + (t150 / 0.2e1 + t23 / 0.2e1 - t12 * mrSges(6,3) + t88 * t156 - t78 * t155 + (t94 * t156 - t28 * mrSges(6,3) + pkin(5) * t29 - t148 / 0.2e1 - t45 / 0.2e1 + t44 * t158) * qJD(5) + (-m(6) * t12 - t38 + (-m(6) * t28 - t60) * qJD(5)) * pkin(8)) * t109 - t164 * t87 + (-m(6) * pkin(4) + t146) * t43 + (-t118 * t2 - t3 * t86 + t65 * t8 - t66 * t9) * mrSges(7,3) + (t147 + t100) * t83 / 0.2e1 + t79 * (Ifges(7,5) * t86 - Ifges(7,6) * t118) / 0.2e1 - t118 * t4 / 0.2e1 + m(7) * (t2 * t72 + t26 * t99 + t3 * t71 + t36 * t9 + t37 * t8) + t99 * t6 + t86 * t5 / 0.2e1 - Ifges(5,5) * t78 - Ifges(5,6) * t79 + t26 * t68 + t19 * t69 / 0.2e1 + t18 * t70 / 0.2e1 + t71 * t13 + t72 * t14 - t65 * t25 / 0.2e1 - t66 * t24 / 0.2e1 - t54 * t32 / 0.2e1 - t55 * t33 / 0.2e1 + t36 * t40 + t37 * t41 - t42 * mrSges(5,2) + t44 * t31 - pkin(4) * t30; t75 + t166 * t83 - mrSges(6,3) * t127 + (t68 + t146) * t79 + m(7) * (t133 * t83 + t18 * t72 + t19 * t71 - t36 * t55 - t37 * t54 + t79 * t99) + m(6) * (-pkin(8) * t127 - t153) + (-t118 * t18 - t19 * t86 - t54 * t65 + t55 * t66) * mrSges(7,3); m(7) * (-t118 * t37 + t36 * t86 - t65 * t72 - t66 * t71); (t133 * t99 + t36 * t72 + t37 * t71) * t162 - t66 * t69 - t118 * t32 + 0.2e1 * t68 * t133 + 0.2e1 * t99 * t31 - t65 * t70 + t86 * t33 - t93 * t138 - 0.2e1 * pkin(4) * t87 + t109 * t89 + (qJD(5) * t94 + t88) * t112 + 0.2e1 * (-t118 * t36 - t37 * t86 + t65 * t71 - t66 * t72) * mrSges(7,3); -Ifges(6,5) * t132 + t12 * mrSges(6,1) - t11 * mrSges(6,2) - t115 * Ifges(6,6) + (-t41 * t136 + t111 * t13 + m(7) * (t108 * t2 + t111 * t3 + t135 * t9 - t136 * t8) + t40 * t135 + t108 * t14) * pkin(5) + t116 + t165; (t108 * t18 + t111 * t19 + (t108 * t54 - t111 * t55) * qJD(6)) * t158 - t167; (-t108 * t65 - t111 * t66 + (t108 * t118 + t111 * t86) * qJD(6)) * t158 - t166; t100 + (pkin(8) * t92 - t143) * qJD(5) + (m(7) * (t108 * t36 + t111 * t37 + (-t108 * t71 + t111 * t72) * qJD(6)) + (-t108 * t66 + t111 * t65 + (t108 * t86 - t111 * t118) * qJD(6)) * mrSges(7,3)) * pkin(5) + t119; 0.2e1 * t82; t116; -t6; -t31; t119; t82; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
