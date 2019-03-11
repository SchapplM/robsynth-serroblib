% Calculate time derivative of joint inertia matrix for
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPPRRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_inertiaDJ_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPPRRR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:14
% EndTime: 2019-03-08 18:39:19
% DurationCPUTime: 1.94s
% Computational Cost: add. (3097->291), mult. (9698->471), div. (0->0), fcn. (11131->16), ass. (0->145)
t77 = cos(pkin(8));
t85 = cos(qJ(4));
t125 = t77 * t85;
t72 = sin(pkin(8));
t128 = t72 * t85;
t70 = sin(pkin(14));
t73 = sin(pkin(7));
t75 = cos(pkin(14));
t78 = cos(pkin(7));
t82 = sin(qJ(4));
t163 = t78 * t128 + (t75 * t125 - t70 * t82) * t73;
t81 = sin(qJ(5));
t84 = cos(qJ(5));
t124 = -mrSges(6,1) * t84 + mrSges(6,2) * t81 - mrSges(5,1);
t152 = m(6) * pkin(4);
t161 = t124 - t152;
t160 = m(7) * pkin(11) + mrSges(7,3);
t80 = sin(qJ(6));
t83 = cos(qJ(6));
t61 = -mrSges(7,1) * t83 + mrSges(7,2) * t80;
t159 = -m(7) * pkin(5) - mrSges(6,1) + t61;
t158 = 0.2e1 * m(7);
t157 = 2 * pkin(10);
t156 = 0.2e1 * t72;
t155 = m(5) / 0.2e1;
t154 = m(6) / 0.2e1;
t153 = m(7) / 0.2e1;
t76 = cos(pkin(13));
t126 = t76 * t78;
t79 = cos(pkin(6));
t127 = t73 * t79;
t71 = sin(pkin(13));
t74 = sin(pkin(6));
t33 = t71 * t74 * t75 + (t126 * t74 + t127) * t70;
t32 = t75 * t127 + (t126 * t75 - t70 * t71) * t74;
t48 = -t73 * t74 * t76 + t78 * t79;
t99 = t32 * t77 + t48 * t72;
t14 = t33 * t85 + t82 * t99;
t98 = -t32 * t72 + t48 * t77;
t10 = t14 * t84 + t81 * t98;
t135 = t33 * t82;
t11 = (t85 * t99 - t135) * qJD(4);
t3 = qJD(5) * t10 + t11 * t81;
t9 = t14 * t81 - t84 * t98;
t151 = t3 * t9;
t144 = Ifges(7,4) * t80;
t63 = Ifges(7,2) * t83 + t144;
t150 = -t63 / 0.2e1;
t149 = -t80 / 0.2e1;
t148 = pkin(10) * t84;
t147 = t3 * t81;
t4 = -qJD(5) * t9 + t11 * t84;
t146 = t4 * t84;
t145 = mrSges(7,3) * t81;
t143 = Ifges(7,4) * t83;
t142 = Ifges(7,6) * t80;
t141 = Ifges(7,6) * t84;
t12 = t14 * qJD(4);
t13 = -t125 * t32 - t128 * t48 + t135;
t140 = t12 * t13;
t129 = t72 * t82;
t35 = t78 * t129 + (t75 * t77 * t82 + t70 * t85) * t73;
t95 = -t72 * t73 * t75 + t77 * t78;
t22 = t81 * t35 - t84 * t95;
t29 = t163 * qJD(4);
t15 = -qJD(5) * t22 + t84 * t29;
t139 = t15 * t84;
t23 = t84 * t35 + t81 * t95;
t16 = qJD(5) * t23 + t81 * t29;
t138 = t16 * t22;
t137 = t16 * t81;
t30 = t35 * qJD(4);
t136 = t30 * t163;
t121 = qJD(4) * t85;
t112 = t72 * t121;
t49 = t129 * t81 - t84 * t77;
t36 = -qJD(5) * t49 + t112 * t84;
t134 = t36 * t84;
t50 = t129 * t84 + t77 * t81;
t37 = qJD(5) * t50 + t112 * t81;
t133 = t37 * t81;
t132 = t49 * t37;
t60 = -pkin(5) * t84 - pkin(11) * t81 - pkin(4);
t131 = t60 * t80;
t119 = qJD(5) * t84;
t110 = t83 * t119;
t120 = qJD(5) * t81;
t123 = -Ifges(7,5) * t110 - Ifges(7,3) * t120;
t122 = qJD(4) * t82;
t118 = qJD(6) * t81;
t117 = qJD(6) * t83;
t116 = qJD(6) * t84;
t113 = t72 * t122;
t111 = t80 * t118;
t109 = (2 * Ifges(6,4)) + t142;
t59 = (pkin(5) * t81 - pkin(11) * t84) * qJD(5);
t24 = t60 * t117 + t59 * t80 + (-t116 * t80 - t120 * t83) * pkin(10);
t43 = -t148 * t80 + t60 * t83;
t108 = -t43 * qJD(6) + t24;
t107 = t16 * t9 + t22 * t3;
t106 = t49 * t3 + t37 * t9;
t105 = mrSges(7,1) * t80 + mrSges(7,2) * t83;
t104 = Ifges(7,1) * t83 - t144;
t64 = Ifges(7,1) * t80 + t143;
t103 = -Ifges(7,2) * t80 + t143;
t102 = Ifges(7,5) * t80 + Ifges(7,6) * t83;
t6 = t10 * t83 + t13 * t80;
t5 = -t10 * t80 + t13 * t83;
t101 = -t12 * t163 + t13 * t30;
t100 = t49 * t16 + t37 * t22;
t18 = -t163 * t80 + t23 * t83;
t17 = -t163 * t83 - t23 * t80;
t97 = t119 * t9 + t147;
t38 = -t128 * t83 - t80 * t50;
t96 = t128 * t80 - t83 * t50;
t94 = -t12 * t85 + t122 * t13;
t93 = -t122 * t163 - t30 * t85;
t92 = t119 * t22 + t137;
t91 = t119 * t49 + t133;
t90 = t110 - t111;
t89 = t117 * t81 + t119 * t80;
t69 = Ifges(7,5) * t117;
t58 = -mrSges(7,1) * t84 - t145 * t83;
t57 = mrSges(7,2) * t84 - t145 * t80;
t55 = t104 * qJD(6);
t54 = t103 * qJD(6);
t53 = (mrSges(6,1) * t81 + mrSges(6,2) * t84) * qJD(5);
t52 = t105 * qJD(6);
t51 = t105 * t81;
t47 = -t84 * Ifges(7,5) + t104 * t81;
t46 = t103 * t81 - t141;
t44 = t148 * t83 + t131;
t41 = -mrSges(7,2) * t120 - mrSges(7,3) * t89;
t40 = mrSges(7,1) * t120 - mrSges(7,3) * t90;
t31 = mrSges(7,1) * t89 + mrSges(7,2) * t90;
t27 = -t64 * t118 + (Ifges(7,5) * t81 + t104 * t84) * qJD(5);
t26 = -t63 * t118 + (Ifges(7,6) * t81 + t103 * t84) * qJD(5);
t25 = -qJD(6) * t131 + t59 * t83 + (-t116 * t83 + t120 * t80) * pkin(10);
t20 = qJD(6) * t96 + t113 * t83 - t80 * t36;
t19 = qJD(6) * t38 + t113 * t80 + t83 * t36;
t8 = -qJD(6) * t18 - t15 * t80 + t30 * t83;
t7 = qJD(6) * t17 + t15 * t83 + t30 * t80;
t2 = qJD(6) * t5 + t12 * t80 + t4 * t83;
t1 = -qJD(6) * t6 + t12 * t83 - t4 * t80;
t21 = [0.2e1 * m(7) * (t1 * t5 + t2 * t6 + t151) + 0.2e1 * m(6) * (t10 * t4 + t140 + t151) + 0.2e1 * m(5) * (t11 * t14 + t140); m(7) * (t1 * t17 + t18 * t2 + t5 * t8 + t6 * t7 + t107) + m(6) * (t10 * t15 + t23 * t4 + t101 + t107) + m(5) * (t11 * t35 + t14 * t29 + t101); 0.2e1 * m(7) * (t17 * t8 + t18 * t7 + t138) + 0.2e1 * m(6) * (t15 * t23 - t136 + t138) + 0.2e1 * m(5) * (t29 * t35 - t136); m(7) * (t1 * t38 + t19 * t6 - t2 * t96 + t20 * t5 + t106) + m(6) * (t36 * t10 + t50 * t4 + t106) + (t94 * t154 + (t11 * t82 + t121 * t14 + t94) * t155) * t156; m(7) * (t17 * t20 + t18 * t19 + t38 * t8 - t7 * t96 + t100) + m(6) * (t50 * t15 + t36 * t23 + t100) + (t93 * t154 + (t121 * t35 + t29 * t82 + t93) * t155) * t156; 0.2e1 * m(7) * (-t19 * t96 + t20 * t38 + t132) + 0.2e1 * m(6) * (-t121 * t72 ^ 2 * t82 + t50 * t36 + t132); -t11 * mrSges(5,2) + t1 * t58 + t13 * t53 + t2 * t57 + t3 * t51 + t9 * t31 + t5 * t40 + t6 * t41 + m(7) * (t1 * t43 + t2 * t44 + t24 * t6 + t25 * t5) + (t97 * t153 + (-t10 * t120 + t146 + t97) * t154) * t157 + (t147 + t146 + (-t10 * t81 + t84 * t9) * qJD(5)) * mrSges(6,3) + t161 * t12; -t29 * mrSges(5,2) + t16 * t51 + t17 * t40 + t18 * t41 + t22 * t31 - t163 * t53 + t7 * t57 + t8 * t58 + m(7) * (t17 * t25 + t18 * t24 + t43 * t8 + t44 * t7) + (t92 * t153 + (-t120 * t23 + t139 + t92) * t154) * t157 + (t139 + t137 + (t22 * t84 - t23 * t81) * qJD(5)) * mrSges(6,3) + t161 * t30; t19 * t57 + t20 * t58 + t49 * t31 + t37 * t51 + t38 * t40 - t96 * t41 + (-t85 * t53 + (-mrSges(5,2) * t85 + t124 * t82) * qJD(4)) * t72 + m(7) * (t19 * t44 + t20 * t43 - t24 * t96 + t25 * t38) - t113 * t152 + (t91 * t153 + (-t120 * t50 + t134 + t91) * t154) * t157 + (t134 + t133 + (t49 * t84 - t50 * t81) * qJD(5)) * mrSges(6,3); 0.2e1 * t24 * t57 + 0.2e1 * t44 * t41 + 0.2e1 * t25 * t58 + 0.2e1 * t43 * t40 + (t24 * t44 + t25 * t43) * t158 - 0.2e1 * pkin(4) * t53 + ((t109 * t84 + t157 * t51 - t80 * t46 + t83 * t47) * qJD(5) + t123) * t84 + (t31 * t157 - t80 * t26 + t83 * t27 + (t102 * t84 - t83 * t46 - t80 * t47) * qJD(6) + ((Ifges(7,5) * t83 - t109) * t81 + ((pkin(10) ^ 2) * t158 + (2 * Ifges(6,1)) - (2 * Ifges(6,2)) - Ifges(7,3)) * t84) * qJD(5)) * t81; -t4 * mrSges(6,2) + t9 * t52 + t160 * (-t1 * t80 + t2 * t83 + (-t5 * t83 - t6 * t80) * qJD(6)) + t159 * t3; -t15 * mrSges(6,2) + t22 * t52 + t160 * (t7 * t83 - t8 * t80 + (-t17 * t83 - t18 * t80) * qJD(6)) + t159 * t16; -t36 * mrSges(6,2) + t49 * t52 + t160 * (t19 * t83 - t20 * t80 + (-t38 * t83 + t80 * t96) * qJD(6)) + t159 * t37; -pkin(5) * t31 + (-t69 / 0.2e1 + (pkin(10) * t159 + Ifges(6,5)) * qJD(5)) * t84 + (t119 * t150 - t25 * mrSges(7,3) + t27 / 0.2e1 + (-t44 * mrSges(7,3) - t46 / 0.2e1 + t141 / 0.2e1) * qJD(6) + (-qJD(6) * t57 - t40 + m(7) * (-t44 * qJD(6) - t25)) * pkin(11)) * t80 + (t64 * t119 / 0.2e1 + qJD(6) * t47 / 0.2e1 + t26 / 0.2e1 + t108 * mrSges(7,3) + (m(7) * t108 - qJD(6) * t58 + t41) * pkin(11)) * t83 + (t54 * t149 + t83 * t55 / 0.2e1 + (t149 * t64 + t150 * t83) * qJD(6) + pkin(10) * t52 + (-Ifges(6,6) + t102 / 0.2e1 + pkin(10) * mrSges(6,2)) * qJD(5)) * t81; -0.2e1 * pkin(5) * t52 + t54 * t83 + t55 * t80 + (-t63 * t80 + t64 * t83) * qJD(6); mrSges(7,1) * t1 - mrSges(7,2) * t2; mrSges(7,1) * t8 - mrSges(7,2) * t7; mrSges(7,1) * t20 - mrSges(7,2) * t19; mrSges(7,1) * t25 - mrSges(7,2) * t24 - Ifges(7,5) * t111 - Ifges(7,6) * t89 - t123; t69 + (pkin(11) * t61 - t142) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
