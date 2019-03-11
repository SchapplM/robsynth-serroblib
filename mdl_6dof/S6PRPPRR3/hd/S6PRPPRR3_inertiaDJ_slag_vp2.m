% Calculate time derivative of joint inertia matrix for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPPRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:30
% EndTime: 2019-03-08 19:21:33
% DurationCPUTime: 1.49s
% Computational Cost: add. (1263->257), mult. (3190->408), div. (0->0), fcn. (2825->10), ass. (0->123)
t66 = cos(qJ(6));
t102 = qJD(6) * t66;
t67 = cos(qJ(5));
t105 = qJD(5) * t67;
t63 = sin(qJ(6));
t64 = sin(qJ(5));
t73 = t64 * t102 + t63 * t105;
t111 = cos(pkin(6));
t60 = sin(pkin(11));
t68 = cos(qJ(2));
t117 = t60 * t68;
t61 = sin(pkin(6));
t62 = cos(pkin(11));
t65 = sin(qJ(2));
t25 = (t62 * t65 - t117) * t61;
t75 = -t111 * t67 - t25 * t64;
t107 = t75 * qJD(5);
t142 = m(4) * qJ(3);
t112 = t63 ^ 2 + t66 ^ 2;
t141 = m(7) * pkin(9) + mrSges(7,3);
t24 = (t60 * t65 + t62 * t68) * t61;
t21 = qJD(2) * t24;
t6 = t21 * t67 + t107;
t126 = t6 * t67;
t88 = t111 * t64;
t5 = -qJD(5) * t88 + t25 * t105 + t21 * t64;
t127 = t5 * t64;
t15 = t25 * t67 - t88;
t140 = (-t15 * t64 - t67 * t75) * qJD(5) + t126 + t127;
t44 = -mrSges(7,1) * t66 + mrSges(7,2) * t63;
t85 = -m(7) * pkin(5) - mrSges(6,1) + t44;
t115 = t66 * t67;
t69 = -pkin(2) - pkin(3);
t40 = -t60 * qJ(3) + t62 * t69;
t32 = pkin(4) - t40;
t22 = t67 * pkin(5) + t64 * pkin(9) + t32;
t41 = t62 * qJ(3) + t60 * t69;
t33 = -pkin(8) + t41;
t10 = t33 * t115 + t22 * t63;
t106 = qJD(5) * t64;
t103 = qJD(6) * t64;
t93 = t63 * t103;
t74 = t66 * t105 - t93;
t18 = -mrSges(7,1) * t106 + t74 * mrSges(7,3);
t108 = qJD(3) * t62;
t97 = t67 * t108;
t72 = qJD(6) * t22 - t33 * t106 + t97;
t109 = qJD(3) * t60;
t129 = pkin(9) * t67;
t130 = pkin(5) * t64;
t84 = -qJD(6) * t33 * t67 + t109 + (t129 - t130) * qJD(5);
t4 = -t72 * t63 + t84 * t66;
t139 = m(7) * (-t10 * qJD(6) - t4) - t18;
t138 = 2 * m(6);
t137 = 0.2e1 * m(7);
t136 = 0.2e1 * t33;
t135 = m(6) / 0.2e1;
t134 = m(7) / 0.2e1;
t124 = Ifges(7,4) * t63;
t46 = Ifges(7,2) * t66 + t124;
t133 = t46 / 0.2e1;
t132 = t63 / 0.2e1;
t131 = -t66 / 0.2e1;
t128 = t75 * t5;
t125 = mrSges(7,3) * t64;
t123 = Ifges(7,4) * t66;
t122 = Ifges(7,6) * t63;
t110 = qJD(2) * t61;
t99 = t65 * t110;
t20 = t110 * t117 - t62 * t99;
t121 = t20 * t24;
t120 = t20 * t62;
t57 = t64 ^ 2;
t119 = t33 * t57;
t59 = t67 ^ 2;
t118 = t33 * t59;
t116 = t63 * t67;
t114 = t67 * Ifges(7,6);
t113 = mrSges(6,1) * t67 - mrSges(6,2) * t64 + mrSges(5,1);
t104 = qJD(6) * t63;
t100 = Ifges(7,5) * t93 + t73 * Ifges(7,6);
t98 = t60 * t108;
t96 = t60 * t106;
t94 = t64 * t105;
t90 = -Ifges(7,5) * t66 + (2 * Ifges(6,4));
t3 = t84 * t63 + t72 * t66;
t9 = -t33 * t116 + t22 * t66;
t89 = -t9 * qJD(6) + t3;
t87 = 0.2e1 * t94;
t86 = t75 * t64 * t108;
t8 = t15 * t66 + t24 * t63;
t1 = -t8 * qJD(6) + t20 * t66 - t6 * t63;
t7 = -t15 * t63 + t24 * t66;
t2 = t7 * qJD(6) + t20 * t63 + t6 * t66;
t83 = -t1 * t63 + t2 * t66;
t81 = -t64 * mrSges(6,1) - t67 * mrSges(6,2);
t80 = mrSges(7,1) * t63 + mrSges(7,2) * t66;
t79 = Ifges(7,1) * t66 - t124;
t47 = Ifges(7,1) * t63 + t123;
t78 = -Ifges(7,2) * t63 + t123;
t30 = t60 * t115 - t62 * t63;
t29 = -t60 * t116 - t62 * t66;
t76 = -t105 * t75 + t127;
t16 = t29 * qJD(6) - t66 * t96;
t17 = -t30 * qJD(6) + t63 * t96;
t70 = t16 * t66 - t17 * t63 + (-t29 * t66 - t30 * t63) * qJD(6);
t53 = Ifges(7,5) * t102;
t42 = t57 * t98;
t39 = mrSges(7,1) * t67 + t66 * t125;
t38 = -mrSges(7,2) * t67 + t63 * t125;
t37 = t79 * qJD(6);
t36 = t78 * qJD(6);
t35 = t81 * qJD(5);
t34 = t80 * qJD(6);
t31 = t80 * t64;
t27 = Ifges(7,5) * t67 - t79 * t64;
t26 = -t78 * t64 + t114;
t23 = t108 * t119;
t19 = mrSges(7,2) * t106 + t73 * mrSges(7,3);
t13 = -t73 * mrSges(7,1) - t74 * mrSges(7,2);
t12 = t47 * t103 + (-Ifges(7,5) * t64 - t79 * t67) * qJD(5);
t11 = t46 * t103 + (-Ifges(7,6) * t64 - t78 * t67) * qJD(5);
t14 = [0.2e1 * m(7) * (t1 * t7 + t2 * t8 - t128) + 0.2e1 * m(6) * (t15 * t6 + t121 - t128) + 0.2e1 * m(5) * (t21 * t25 + t121); t21 * mrSges(5,2) + t1 * t39 - t75 * t13 + t7 * t18 + t8 * t19 + t2 * t38 + t24 * t35 - t5 * t31 + t113 * t20 + m(7) * (t1 * t9 + t10 * t2 + t3 * t8 + t4 * t7 - t86) + m(6) * (t24 * t109 + t15 * t97 + t20 * t32 - t86) + m(5) * (-t20 * t40 + t21 * t41 + (t24 * t60 + t25 * t62) * qJD(3)) + (t76 * t134 + (-t15 * t106 + t126 + t76) * t135) * t136 - t140 * mrSges(6,3) + ((-mrSges(3,2) + mrSges(4,3) + t142) * qJD(2) * t68 + (m(4) * qJD(3) + (-m(4) * pkin(2) - mrSges(3,1) - mrSges(4,1)) * qJD(2)) * t65) * t61; 0.2e1 * t10 * t19 + 0.2e1 * t9 * t18 + 0.2e1 * t3 * t38 + 0.2e1 * t32 * t35 + 0.2e1 * t4 * t39 + (t10 * t3 + t4 * t9 + t23) * t137 + t23 * t138 + ((-t31 * t136 + t63 * t26 - t66 * t27 + t90 * t67) * qJD(5) + t100) * t67 + (0.2e1 * t142 + 0.2e1 * mrSges(4,3) + (t62 * t118 + t32 * t60) * t138 + 0.2e1 * m(5) * (-t40 * t60 + t41 * t62) + 0.2e1 * t113 * t60 + 0.2e1 * (mrSges(5,2) + (-t57 - t59) * mrSges(6,3)) * t62) * qJD(3) + (-0.2e1 * t31 * t108 + t63 * t11 - t66 * t12 + t13 * t136 + (t66 * t26 + t63 * t27) * qJD(6) + ((-t90 - t122) * t64 + 0.2e1 * (m(7) * t33 ^ 2 + Ifges(6,1) - Ifges(6,2) - Ifges(7,3)) * t67) * qJD(5)) * t64; m(4) * t99 + m(7) * (t1 * t29 + t16 * t8 + t17 * t7 + t2 * t30 + t76 * t60) + m(6) * (t140 * t60 - t120) + m(5) * (t21 * t60 - t120); t16 * t38 + t17 * t39 + t29 * t18 + t30 * t19 - t62 * t35 + (-t31 * t105 + t64 * t13) * t60 + m(7) * (t33 * t60 * t87 + t10 * t16 + t17 * t9 + t29 * t4 + t3 * t30 + t42) + m(6) * (t42 + (t59 - 0.1e1) * t98); (t60 ^ 2 * t94 + t16 * t30 + t17 * t29) * t137; 0.2e1 * ((-t5 + (-t63 * t7 + t66 * t8) * qJD(5)) * t134 + (qJD(5) * t15 - t5) * t135) * t67 + 0.2e1 * ((-t7 * t102 - t8 * t104 - t107 + t83) * t134 + (t6 - t107) * t135) * t64; -t67 * t13 + (m(7) * (t10 * t115 - t9 * t116 - t118 + t119) + t38 * t115 - t39 * t116) * qJD(5) + (m(7) * (-t9 * t102 + t3 * t66 - t97) - t38 * t104 + t66 * t19 - t39 * t102 - qJD(5) * t31 + t139 * t63) * t64; m(7) * (t70 * t64 + ((-t29 * t63 + t30 * t66) * t67 + (t57 - t59) * t60) * qJD(5)); m(7) * (-0.1e1 + t112) * t87; -t6 * mrSges(6,2) - t75 * t34 + t141 * ((-t63 * t8 - t66 * t7) * qJD(6) + t83) + t85 * t5; -pkin(5) * t13 + (t53 / 0.2e1 - mrSges(6,2) * t108 + (t85 * t33 - Ifges(6,5)) * qJD(5)) * t67 + (-t4 * mrSges(7,3) + t12 / 0.2e1 + t105 * t133 + (-t26 / 0.2e1 - t114 / 0.2e1 - t10 * mrSges(7,3)) * qJD(6) + (-qJD(6) * t38 + t139) * pkin(9)) * t63 + (qJD(6) * t27 / 0.2e1 + t11 / 0.2e1 - t47 * t105 / 0.2e1 + t89 * mrSges(7,3) + (m(7) * t89 - qJD(6) * t39 + t19) * pkin(9)) * t66 + (t37 * t131 + t36 * t132 + t33 * t34 + (t47 * t132 + t66 * t133) * qJD(6) + (-Ifges(7,5) * t63 / 0.2e1 + Ifges(7,6) * t131 + Ifges(6,6) + t33 * mrSges(6,2)) * qJD(5) + t85 * t108) * t64; t141 * t70 + ((qJD(5) * mrSges(6,2) + t34) * t64 + t85 * t105) * t60; -t67 * t34 + (t64 * t44 + m(7) * (t112 * t129 - t130) + t112 * t67 * mrSges(7,3) + t81) * qJD(5); -0.2e1 * pkin(5) * t34 + t36 * t66 + t37 * t63 + (-t46 * t63 + t47 * t66) * qJD(6); mrSges(7,1) * t1 - mrSges(7,2) * t2; mrSges(7,1) * t4 - mrSges(7,2) * t3 + (-Ifges(7,5) * t115 - Ifges(7,3) * t64) * qJD(5) + t100; mrSges(7,1) * t17 - mrSges(7,2) * t16; t13; t53 + (t44 * pkin(9) - t122) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
