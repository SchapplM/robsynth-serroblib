% Calculate time derivative of joint inertia matrix for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:29:37
% EndTime: 2019-03-08 19:29:41
% DurationCPUTime: 1.76s
% Computational Cost: add. (2085->319), mult. (5677->501), div. (0->0), fcn. (5366->12), ass. (0->134)
t94 = sin(pkin(12));
t97 = cos(pkin(12));
t151 = mrSges(6,1) * t94 + mrSges(6,2) * t97;
t93 = t97 ^ 2;
t130 = t94 ^ 2 + t93;
t102 = cos(qJ(6));
t99 = sin(qJ(6));
t76 = t102 * t97 - t94 * t99;
t69 = t76 * qJD(6);
t149 = 2 * m(6);
t148 = 2 * m(7);
t95 = sin(pkin(11));
t88 = pkin(2) * t95 + pkin(8);
t147 = 0.2e1 * t88;
t146 = m(6) / 0.2e1;
t145 = t76 / 0.2e1;
t77 = t102 * t94 + t97 * t99;
t144 = t77 / 0.2e1;
t143 = t97 / 0.2e1;
t140 = Ifges(6,4) * t94;
t139 = Ifges(6,4) * t97;
t100 = sin(qJ(4));
t138 = pkin(4) * t100;
t103 = cos(qJ(4));
t122 = cos(pkin(6));
t101 = sin(qJ(2));
t104 = cos(qJ(2));
t96 = sin(pkin(6));
t98 = cos(pkin(11));
t59 = (t101 * t98 + t104 * t95) * t96;
t105 = -t59 * t100 + t103 * t122;
t46 = t100 * t122 + t59 * t103;
t107 = t101 * t95 - t104 * t98;
t120 = qJD(2) * t96;
t54 = t107 * t120;
t15 = qJD(4) * t46 - t54 * t100;
t12 = t105 * t15;
t53 = qJD(2) * t59;
t58 = t107 * t96;
t137 = t53 * t58;
t136 = t94 * Ifges(6,2);
t135 = pkin(9) + qJ(5);
t118 = qJD(4) * t103;
t70 = t77 * qJD(6);
t30 = -t100 * t70 + t118 * t76;
t31 = -t100 * t69 - t118 * t77;
t11 = -t31 * mrSges(7,1) + t30 * mrSges(7,2);
t64 = t151 * t118;
t134 = t11 + t64;
t60 = t77 * t100;
t61 = t76 * t100;
t32 = mrSges(7,1) * t60 + mrSges(7,2) * t61;
t71 = t151 * t100;
t133 = t32 + t71;
t132 = Ifges(7,5) * t69 - Ifges(7,6) * t70;
t119 = qJD(4) * t100;
t113 = t88 * t119;
t117 = qJD(5) * t100;
t121 = qJ(5) * t103;
t68 = -t117 + (-t121 + t138) * qJD(4);
t37 = t94 * t113 + t97 * t68;
t123 = t97 * t103;
t89 = -pkin(2) * t98 - pkin(3);
t74 = -pkin(4) * t103 - qJ(5) * t100 + t89;
t41 = t88 * t123 + t94 * t74;
t131 = -mrSges(6,1) * t97 + mrSges(6,2) * t94 - mrSges(5,1);
t129 = t100 * t94;
t128 = t100 * t97;
t127 = t103 * t15;
t126 = t15 * t100;
t16 = qJD(4) * t105 - t54 * t103;
t125 = t16 * t103;
t124 = t94 * t103;
t42 = -mrSges(7,1) * t76 + mrSges(7,2) * t77;
t116 = t42 + t131;
t115 = -Ifges(7,5) * t30 - Ifges(7,6) * t31 - Ifges(7,3) * t119;
t114 = m(6) * t130;
t112 = pkin(5) * t94 + t88;
t111 = t130 * mrSges(6,3);
t110 = -Ifges(6,5) * t97 + Ifges(6,6) * t94;
t13 = -t16 * t94 + t53 * t97;
t14 = t16 * t97 + t53 * t94;
t109 = -t13 * t94 + t14 * t97;
t17 = -t46 * t94 + t58 * t97;
t18 = t46 * t97 + t58 * t94;
t108 = -t17 * t94 + t18 * t97;
t6 = t102 * t18 + t17 * t99;
t5 = t102 * t17 - t18 * t99;
t63 = t97 * t74;
t29 = -pkin(9) * t128 + t63 + (-t88 * t94 - pkin(5)) * t103;
t33 = -pkin(9) * t129 + t41;
t10 = t102 * t33 + t29 * t99;
t9 = t102 * t29 - t33 * t99;
t82 = t135 * t94;
t84 = t135 * t97;
t48 = t102 * t84 - t82 * t99;
t47 = -t102 * t82 - t84 * t99;
t106 = -t105 * t118 + t126;
t90 = -pkin(5) * t97 - pkin(4);
t81 = (t100 * mrSges(5,1) + t103 * mrSges(5,2)) * qJD(4);
t79 = -mrSges(6,1) * t103 - mrSges(6,3) * t128;
t78 = mrSges(6,2) * t103 - mrSges(6,3) * t129;
t73 = (mrSges(6,1) * t100 - mrSges(6,3) * t123) * qJD(4);
t72 = (-mrSges(6,2) * t100 - mrSges(6,3) * t124) * qJD(4);
t67 = t112 * t100;
t57 = t112 * t118;
t55 = t94 * t68;
t52 = (t100 * Ifges(6,5) + (t97 * Ifges(6,1) - t140) * t103) * qJD(4);
t51 = (t100 * Ifges(6,6) + (-t136 + t139) * t103) * qJD(4);
t50 = -mrSges(7,1) * t103 - mrSges(7,3) * t61;
t49 = mrSges(7,2) * t103 - mrSges(7,3) * t60;
t44 = Ifges(7,1) * t77 + Ifges(7,4) * t76;
t43 = Ifges(7,4) * t77 + Ifges(7,2) * t76;
t40 = -t124 * t88 + t63;
t39 = t105 * t119;
t38 = -t113 * t97 + t55;
t36 = Ifges(7,1) * t69 - Ifges(7,4) * t70;
t35 = Ifges(7,4) * t69 - Ifges(7,2) * t70;
t34 = mrSges(7,1) * t70 + mrSges(7,2) * t69;
t26 = t55 + (-pkin(9) * t124 - t128 * t88) * qJD(4);
t25 = Ifges(7,1) * t61 - Ifges(7,4) * t60 - Ifges(7,5) * t103;
t24 = Ifges(7,4) * t61 - Ifges(7,2) * t60 - Ifges(7,6) * t103;
t23 = -qJD(5) * t77 - qJD(6) * t48;
t22 = qJD(5) * t76 + qJD(6) * t47;
t21 = (pkin(5) * t100 - pkin(9) * t123) * qJD(4) + t37;
t20 = -mrSges(7,2) * t119 + mrSges(7,3) * t31;
t19 = mrSges(7,1) * t119 - mrSges(7,3) * t30;
t8 = Ifges(7,1) * t30 + Ifges(7,4) * t31 + Ifges(7,5) * t119;
t7 = Ifges(7,4) * t30 + Ifges(7,2) * t31 + Ifges(7,6) * t119;
t4 = -qJD(6) * t10 + t102 * t21 - t26 * t99;
t3 = qJD(6) * t9 + t102 * t26 + t21 * t99;
t2 = -qJD(6) * t6 + t102 * t13 - t14 * t99;
t1 = qJD(6) * t5 + t102 * t14 + t13 * t99;
t27 = [0.2e1 * m(7) * (t1 * t6 + t2 * t5 - t12) + 0.2e1 * m(5) * (t16 * t46 - t12 + t137) + 0.2e1 * m(6) * (t13 * t17 + t14 * t18 - t12) + 0.2e1 * m(4) * (-t54 * t59 + t137); t54 * mrSges(4,2) + t1 * t49 + t13 * t79 + t14 * t78 + t17 * t73 + t18 * t72 + t5 * t19 + t2 * t50 + t6 * t20 + t58 * t81 + (-t103 * mrSges(5,1) + t100 * mrSges(5,2) - mrSges(4,1)) * t53 - t134 * t105 + t133 * t15 + (-mrSges(3,1) * t101 - mrSges(3,2) * t104) * t120 + (t126 + t125 + (-t100 * t46 - t103 * t105) * qJD(4)) * mrSges(5,3) + m(7) * (t1 * t10 - t105 * t57 + t15 * t67 + t2 * t9 + t3 * t6 + t4 * t5) + m(6) * (t13 * t40 + t14 * t41 + t17 * t37 + t18 * t38) + m(5) * t89 * t53 + (t106 * t146 + m(5) * (-t46 * t119 + t106 + t125) / 0.2e1) * t147 + m(4) * (-t53 * t98 - t54 * t95) * pkin(2); 0.2e1 * t10 * t20 + 0.2e1 * t67 * t11 + 0.2e1 * t9 * t19 + t31 * t24 + t30 * t25 + 0.2e1 * t3 * t49 + 0.2e1 * t57 * t32 + 0.2e1 * t37 * t79 + 0.2e1 * t38 * t78 + 0.2e1 * t4 * t50 + 0.2e1 * t40 * t73 + 0.2e1 * t41 * t72 - t60 * t7 + t61 * t8 + 0.2e1 * t89 * t81 + (t10 * t3 + t4 * t9 + t57 * t67) * t148 + (t37 * t40 + t38 * t41) * t149 + (-t94 * t51 + t97 * t52 + t64 * t147 + (Ifges(7,5) * t61 - Ifges(7,6) * t60 + (-(2 * Ifges(5,4)) - t110) * t100) * qJD(4)) * t100 + ((t71 * t147 + 0.2e1 * (Ifges(5,4) + t110) * t103 + (Ifges(6,1) * t93 - Ifges(7,3) + (2 * Ifges(5,1)) - (2 * Ifges(5,2)) + t88 ^ 2 * t149 - (2 * Ifges(6,3)) + (t136 - 0.2e1 * t139) * t94) * t100) * qJD(4) + t115) * t103; m(7) * (t1 * t61 - t2 * t60 + t30 * t6 + t31 * t5 - t127 - t39) + m(5) * (t100 * t16 - t127 + (-t100 * t105 + t103 * t46) * qJD(4)) + m(6) * (-t39 + t109 * t100 + (qJD(4) * t108 - t15) * t103); -t60 * t19 + t61 * t20 + t30 * t49 + t31 * t50 - t134 * t103 + (t97 * t72 - t94 * t73) * t100 + ((t78 * t97 - t79 * t94) * t103 + t133 * t100) * qJD(4) + m(7) * (t10 * t30 - t103 * t57 + t119 * t67 + t3 * t61 + t31 * t9 - t4 * t60) + m(6) * ((-t37 * t94 + t38 * t97) * t100 + (t100 ^ 2 * t88 + (-t103 * t88 - t40 * t94 + t41 * t97) * t103) * qJD(4)); (t30 * t61 - t31 * t60) * t148 + 0.4e1 * ((-0.1e1 + t130) * t146 - m(7) / 0.2e1) * t100 * t118; -t16 * mrSges(5,2) - t105 * t34 + t109 * mrSges(6,3) + t116 * t15 + m(7) * (t1 * t48 + t15 * t90 + t2 * t47 + t22 * t6 + t23 * t5) + m(6) * (-pkin(4) * t15 + qJ(5) * t109 + qJD(5) * t108) + (t1 * t76 - t2 * t77 - t5 * t69 - t6 * t70) * mrSges(7,3); m(7) * (t10 * t22 + t23 * t9 + t3 * t48 + t4 * t47 + t57 * t90) + t31 * t43 / 0.2e1 + t30 * t44 / 0.2e1 + t47 * t19 + t48 * t20 + t22 * t49 + t23 * t50 + t57 * t42 - t60 * t35 / 0.2e1 + t61 * t36 / 0.2e1 - pkin(4) * t64 + t67 * t34 + t69 * t25 / 0.2e1 - t70 * t24 / 0.2e1 + t7 * t145 + t8 * t144 + t90 * t11 - t103 * t132 / 0.2e1 + (m(6) * (qJ(5) * t38 + qJD(5) * t41) + qJ(5) * t72 + t38 * mrSges(6,3) + qJD(5) * t78 + t51 / 0.2e1) * t97 + (m(6) * (-qJ(5) * t37 - qJD(5) * t40) - qJ(5) * t73 - t37 * mrSges(6,3) - qJD(5) * t79 + t52 / 0.2e1) * t94 + (-t10 * t70 + t3 * t76 - t4 * t77 - t69 * t9) * mrSges(7,3) + ((t88 * mrSges(5,2) + Ifges(7,5) * t144 + Ifges(7,6) * t145 + Ifges(6,5) * t94 / 0.2e1 + Ifges(6,6) * t143 - Ifges(5,6)) * t100 + ((Ifges(6,1) * t94 + t139) * t143 - t94 * (Ifges(6,2) * t97 + t140) / 0.2e1 + Ifges(5,5) + (-m(6) * pkin(4) + t131) * t88) * t103) * qJD(4); -t103 * t34 + m(7) * (t22 * t61 - t23 * t60 + t30 * t48 + t31 * t47) + t114 * t117 + ((-mrSges(5,2) + t111) * t103 + m(6) * (t130 * t121 - t138) + (m(7) * t90 + t116) * t100) * qJD(4) + (t30 * t76 - t31 * t77 + t60 * t69 - t61 * t70) * mrSges(7,3); (t22 * t48 + t23 * t47) * t148 - t70 * t43 + t76 * t35 + t69 * t44 + t77 * t36 + 0.2e1 * t90 * t34 + 0.2e1 * (t22 * t76 - t23 * t77 - t47 * t69 - t48 * t70) * mrSges(7,3) + 0.2e1 * (qJ(5) * t114 + t111) * qJD(5); 0.2e1 * (m(7) / 0.2e1 + t146) * t15; m(6) * t118 * t88 + m(7) * t57 + t134; (m(6) + m(7)) * t119; t34; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1; mrSges(7,1) * t4 - mrSges(7,2) * t3 - t115; -t11; mrSges(7,1) * t23 - t22 * mrSges(7,2) + t132; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t27(1) t27(2) t27(4) t27(7) t27(11) t27(16); t27(2) t27(3) t27(5) t27(8) t27(12) t27(17); t27(4) t27(5) t27(6) t27(9) t27(13) t27(18); t27(7) t27(8) t27(9) t27(10) t27(14) t27(19); t27(11) t27(12) t27(13) t27(14) t27(15) t27(20); t27(16) t27(17) t27(18) t27(19) t27(20) t27(21);];
Mq  = res;
