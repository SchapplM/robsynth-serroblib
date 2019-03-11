% Calculate time derivative of joint inertia matrix for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRRPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:45:27
% EndTime: 2019-03-08 18:45:32
% DurationCPUTime: 2.41s
% Computational Cost: add. (2693->351), mult. (8032->556), div. (0->0), fcn. (8195->14), ass. (0->154)
t111 = sin(pkin(13));
t115 = cos(pkin(13));
t184 = mrSges(6,1) * t111 + mrSges(6,2) * t115;
t112 = sin(pkin(12));
t114 = sin(pkin(6));
t118 = cos(pkin(6));
t121 = sin(qJ(3));
t124 = cos(qJ(3));
t116 = cos(pkin(12));
t117 = cos(pkin(7));
t147 = t116 * t117;
t113 = sin(pkin(7));
t150 = t113 * t124;
t183 = (-t112 * t121 + t124 * t147) * t114 + t118 * t150;
t176 = m(6) / 0.2e1;
t182 = m(7) + 0.2e1 * t176;
t120 = sin(qJ(4));
t123 = cos(qJ(4));
t125 = -t113 * t114 * t116 + t117 * t118;
t151 = t113 * t121;
t51 = t118 * t151 + (t112 * t124 + t121 * t147) * t114;
t26 = t125 * t120 + t51 * t123;
t47 = t183 * qJD(3);
t13 = t26 * qJD(4) + t47 * t120;
t25 = t51 * t120 - t125 * t123;
t144 = qJD(3) * t124;
t133 = t113 * t144;
t82 = t117 * t120 + t123 * t151;
t60 = t82 * qJD(4) + t120 * t133;
t81 = -t123 * t117 + t120 * t151;
t169 = t81 * t13 + t60 * t25;
t119 = sin(qJ(6));
t122 = cos(qJ(6));
t129 = t111 * t119 - t115 * t122;
t79 = t129 * qJD(6);
t180 = 0.2e1 * m(6);
t179 = 0.2e1 * m(7);
t178 = 2 * pkin(9);
t110 = t115 ^ 2;
t177 = m(5) / 0.2e1;
t175 = m(5) * pkin(3);
t174 = -t129 / 0.2e1;
t90 = t111 * t122 + t115 * t119;
t173 = t90 / 0.2e1;
t172 = t115 / 0.2e1;
t9 = t25 * t13;
t48 = t51 * qJD(3);
t171 = t48 * t183;
t33 = t81 * t60;
t170 = pkin(10) + qJ(5);
t142 = qJD(4) * t123;
t80 = t90 * qJD(6);
t40 = -t120 * t80 - t129 * t142;
t41 = t120 * t79 - t90 * t142;
t19 = -t41 * mrSges(7,1) + t40 * mrSges(7,2);
t75 = t184 * t142;
t168 = t19 + t75;
t73 = t90 * t120;
t74 = t129 * t120;
t42 = mrSges(7,1) * t73 - mrSges(7,2) * t74;
t83 = t184 * t120;
t167 = t42 + t83;
t166 = -Ifges(7,5) * t79 - Ifges(7,6) * t80;
t165 = -mrSges(6,1) * t115 + mrSges(6,2) * t111 - mrSges(5,1);
t162 = Ifges(6,4) * t111;
t161 = Ifges(6,4) * t115;
t160 = Ifges(6,2) * t111;
t159 = t13 * t120;
t14 = -t25 * qJD(4) + t47 * t123;
t158 = t14 * t123;
t59 = -t81 * qJD(4) + t123 * t133;
t157 = t59 * t123;
t156 = t60 * t120;
t155 = -t123 * mrSges(5,1) + t120 * mrSges(5,2) - mrSges(4,1);
t143 = qJD(4) * t120;
t138 = pkin(9) * t143;
t78 = -qJD(5) * t120 + (pkin(4) * t120 - qJ(5) * t123) * qJD(4);
t55 = t111 * t138 + t115 * t78;
t148 = t115 * t123;
t96 = -pkin(4) * t123 - qJ(5) * t120 - pkin(3);
t69 = pkin(9) * t148 + t111 * t96;
t154 = t111 * t120;
t153 = t111 * t123;
t149 = t115 * t120;
t145 = qJD(3) * t121;
t52 = mrSges(7,1) * t129 + mrSges(7,2) * t90;
t140 = t52 + t165;
t139 = -Ifges(7,5) * t40 - Ifges(7,6) * t41 - Ifges(7,3) * t143;
t135 = pkin(5) * t111 + pkin(9);
t134 = t113 * t145;
t132 = -Ifges(6,5) * t115 + Ifges(6,6) * t111;
t10 = -t111 * t14 + t115 * t48;
t11 = t111 * t48 + t115 * t14;
t131 = -t10 * t111 + t11 * t115;
t35 = -t111 * t59 + t115 * t134;
t36 = t111 * t134 + t115 * t59;
t130 = -t111 * t35 + t115 * t36;
t15 = -t111 * t26 - t115 * t183;
t16 = -t111 * t183 + t115 * t26;
t3 = -t119 * t16 + t122 * t15;
t4 = t119 * t15 + t122 * t16;
t88 = t115 * t96;
t49 = -pkin(10) * t149 + t88 + (-pkin(9) * t111 - pkin(5)) * t123;
t61 = -pkin(10) * t154 + t69;
t21 = -t119 * t61 + t122 * t49;
t22 = t119 * t49 + t122 * t61;
t57 = -t111 * t82 - t115 * t150;
t58 = -t111 * t150 + t115 * t82;
t23 = -t119 * t58 + t122 * t57;
t24 = t119 * t57 + t122 * t58;
t97 = t170 * t111;
t99 = t170 * t115;
t62 = -t119 * t99 - t122 * t97;
t63 = -t119 * t97 + t122 * t99;
t128 = -t124 * t48 - t145 * t183;
t127 = t25 * t142 + t159;
t126 = t81 * t142 + t156;
t107 = -pkin(5) * t115 - pkin(4);
t95 = t135 * t120;
t93 = (mrSges(5,1) * t120 + mrSges(5,2) * t123) * qJD(4);
t92 = -mrSges(6,1) * t123 - mrSges(6,3) * t149;
t91 = mrSges(6,2) * t123 - mrSges(6,3) * t154;
t86 = t135 * t142;
t85 = (mrSges(6,1) * t120 - mrSges(6,3) * t148) * qJD(4);
t84 = (-mrSges(6,2) * t120 - mrSges(6,3) * t153) * qJD(4);
t70 = t111 * t78;
t68 = -pkin(9) * t153 + t88;
t67 = (t120 * Ifges(6,5) + (Ifges(6,1) * t115 - t162) * t123) * qJD(4);
t66 = (t120 * Ifges(6,6) + (-t160 + t161) * t123) * qJD(4);
t65 = -mrSges(7,1) * t123 + mrSges(7,3) * t74;
t64 = mrSges(7,2) * t123 - mrSges(7,3) * t73;
t56 = -t115 * t138 + t70;
t54 = Ifges(7,1) * t90 - Ifges(7,4) * t129;
t53 = Ifges(7,4) * t90 - Ifges(7,2) * t129;
t46 = -Ifges(7,1) * t79 - Ifges(7,4) * t80;
t45 = -Ifges(7,4) * t79 - Ifges(7,2) * t80;
t44 = mrSges(7,1) * t80 - mrSges(7,2) * t79;
t43 = t70 + (-pkin(9) * t149 - pkin(10) * t153) * qJD(4);
t34 = (pkin(5) * t120 - pkin(10) * t148) * qJD(4) + t55;
t32 = -Ifges(7,1) * t74 - Ifges(7,4) * t73 - Ifges(7,5) * t123;
t31 = -Ifges(7,4) * t74 - Ifges(7,2) * t73 - Ifges(7,6) * t123;
t30 = -t90 * qJD(5) - t63 * qJD(6);
t29 = -t129 * qJD(5) + t62 * qJD(6);
t28 = -mrSges(7,2) * t143 + mrSges(7,3) * t41;
t27 = mrSges(7,1) * t143 - mrSges(7,3) * t40;
t18 = Ifges(7,1) * t40 + Ifges(7,4) * t41 + Ifges(7,5) * t143;
t17 = Ifges(7,4) * t40 + Ifges(7,2) * t41 + Ifges(7,6) * t143;
t8 = -t22 * qJD(6) - t119 * t43 + t122 * t34;
t7 = t21 * qJD(6) + t119 * t34 + t122 * t43;
t6 = -t24 * qJD(6) - t119 * t36 + t122 * t35;
t5 = t23 * qJD(6) + t119 * t35 + t122 * t36;
t2 = -t4 * qJD(6) + t10 * t122 - t11 * t119;
t1 = t3 * qJD(6) + t10 * t119 + t11 * t122;
t12 = [0.2e1 * m(7) * (t1 * t4 + t2 * t3 + t9) + 0.2e1 * m(6) * (t10 * t15 + t11 * t16 + t9) + 0.2e1 * m(5) * (t14 * t26 - t171 + t9) + 0.2e1 * m(4) * (t47 * t51 - t171); m(7) * (t1 * t24 + t2 * t23 + t3 * t6 + t4 * t5 + t169) + m(6) * (t10 * t57 + t11 * t58 + t15 * t35 + t16 * t36 + t169) + m(5) * (t82 * t14 + t59 * t26 + t169) + 0.2e1 * (t128 * t177 + m(4) * (t121 * t47 + t51 * t144 + t128) / 0.2e1) * t113; 0.2e1 * m(7) * (t23 * t6 + t24 * t5 + t33) + 0.2e1 * m(5) * (-t113 ^ 2 * t121 * t144 + t82 * t59 + t33) + 0.2e1 * m(6) * (t35 * t57 + t36 * t58 + t33); -t47 * mrSges(4,2) + t1 * t64 + t10 * t92 + t11 * t91 + t15 * t85 + t16 * t84 + t2 * t65 + t3 * t27 + t4 * t28 - t183 * t93 + t168 * t25 + t167 * t13 + m(7) * (t1 * t22 + t13 * t95 + t2 * t21 + t25 * t86 + t3 * t8 + t4 * t7) + m(6) * (t10 * t68 + t11 * t69 + t15 * t55 + t16 * t56) + (t127 * t176 + (-t143 * t26 + t127 + t158) * t177) * t178 + (t159 + t158 + (-t120 * t26 + t123 * t25) * qJD(4)) * mrSges(5,3) + (t155 - t175) * t48; t23 * t27 + t24 * t28 + t35 * t92 + t36 * t91 + t5 * t64 + t57 * t85 + t58 * t84 + t6 * t65 + t168 * t81 + t167 * t60 + (-t124 * t93 + (-mrSges(4,2) * t124 + t121 * t155) * qJD(3)) * t113 + m(7) * (t21 * t6 + t22 * t5 + t23 * t8 + t24 * t7 + t60 * t95 + t81 * t86) - t134 * t175 + m(6) * (t35 * t68 + t36 * t69 + t55 * t57 + t56 * t58) + ((-t143 * t82 + t126 + t157) * t177 + t126 * t176) * t178 + (t156 + t157 + (-t120 * t82 + t123 * t81) * qJD(4)) * mrSges(5,3); -0.2e1 * pkin(3) * t93 - t73 * t17 - t74 * t18 + 0.2e1 * t95 * t19 + 0.2e1 * t21 * t27 + 0.2e1 * t22 * t28 + t41 * t31 + t40 * t32 + 0.2e1 * t86 * t42 + 0.2e1 * t55 * t92 + 0.2e1 * t56 * t91 + 0.2e1 * t7 * t64 + 0.2e1 * t8 * t65 + 0.2e1 * t68 * t85 + 0.2e1 * t69 * t84 + (t21 * t8 + t22 * t7 + t86 * t95) * t179 + (t55 * t68 + t56 * t69) * t180 + (t75 * t178 - t111 * t66 + t115 * t67 + (-Ifges(7,5) * t74 - Ifges(7,6) * t73 + (-(2 * Ifges(5,4)) - t132) * t120) * qJD(4)) * t120 + ((t83 * t178 + 0.2e1 * (Ifges(5,4) + t132) * t123 + (Ifges(6,1) * t110 - (2 * Ifges(5,2)) + (2 * Ifges(5,1)) + (pkin(9) ^ 2) * t180 - Ifges(7,3) - (2 * Ifges(6,3)) + (t160 - 0.2e1 * t161) * t111) * t120) * qJD(4) + t139) * t123; -t14 * mrSges(5,2) + t25 * t44 + t131 * mrSges(6,3) + t140 * t13 + m(7) * (t1 * t63 + t107 * t13 + t2 * t62 + t29 * t4 + t3 * t30) + m(6) * (-pkin(4) * t13 + (-t111 * t15 + t115 * t16) * qJD(5) + t131 * qJ(5)) + (-t1 * t129 - t2 * t90 + t3 * t79 - t4 * t80) * mrSges(7,3); -t59 * mrSges(5,2) + t81 * t44 + t130 * mrSges(6,3) + t140 * t60 + m(7) * (t107 * t60 + t23 * t30 + t24 * t29 + t5 * t63 + t6 * t62) + m(6) * (-pkin(4) * t60 + (-t111 * t57 + t115 * t58) * qJD(5) + t130 * qJ(5)) + (-t129 * t5 + t23 * t79 - t24 * t80 - t6 * t90) * mrSges(7,3); m(7) * (t107 * t86 + t30 * t21 + t29 * t22 + t62 * t8 + t63 * t7) + t41 * t53 / 0.2e1 + t40 * t54 / 0.2e1 + t62 * t27 + t63 * t28 + t29 * t64 + t30 * t65 - t73 * t45 / 0.2e1 - t74 * t46 / 0.2e1 - pkin(4) * t75 - t79 * t32 / 0.2e1 - t80 * t31 / 0.2e1 + t86 * t52 + t17 * t174 + t18 * t173 + t95 * t44 + t107 * t19 - t123 * t166 / 0.2e1 + (m(6) * (qJ(5) * t56 + qJD(5) * t69) + t56 * mrSges(6,3) + qJ(5) * t84 + qJD(5) * t91 + t66 / 0.2e1) * t115 + (m(6) * (-qJ(5) * t55 - qJD(5) * t68) - t55 * mrSges(6,3) - qJ(5) * t85 - qJD(5) * t92 + t67 / 0.2e1) * t111 + (-t129 * t7 + t21 * t79 - t22 * t80 - t8 * t90) * mrSges(7,3) + ((pkin(9) * mrSges(5,2) + Ifges(7,5) * t173 + Ifges(7,6) * t174 + Ifges(6,5) * t111 / 0.2e1 + Ifges(6,6) * t172 - Ifges(5,6)) * t120 + (-t111 * (Ifges(6,2) * t115 + t162) / 0.2e1 + (Ifges(6,1) * t111 + t161) * t172 + Ifges(5,5) + (-m(6) * pkin(4) + t165) * pkin(9)) * t123) * qJD(4); (t29 * t63 + t30 * t62) * t179 + 0.2e1 * t107 * t44 - t80 * t53 - t129 * t45 - t79 * t54 + t90 * t46 + 0.2e1 * (-t129 * t29 - t30 * t90 + t62 * t79 - t63 * t80) * mrSges(7,3) + (qJ(5) * t180 + 0.2e1 * mrSges(6,3)) * qJD(5) * (t111 ^ 2 + t110); t13 * t182; t60 * t182; m(6) * pkin(9) * t142 + m(7) * t86 + t168; t44; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1; mrSges(7,1) * t6 - mrSges(7,2) * t5; mrSges(7,1) * t8 - mrSges(7,2) * t7 - t139; mrSges(7,1) * t30 - mrSges(7,2) * t29 + t166; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
