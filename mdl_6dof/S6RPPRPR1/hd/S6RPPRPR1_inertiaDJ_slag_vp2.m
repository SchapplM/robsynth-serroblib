% Calculate time derivative of joint inertia matrix for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:39:09
% EndTime: 2019-03-09 01:39:13
% DurationCPUTime: 1.49s
% Computational Cost: add. (2854->263), mult. (6344->406), div. (0->0), fcn. (6148->10), ass. (0->108)
t82 = sin(pkin(9)) * pkin(1) + qJ(3);
t126 = pkin(7) + t82;
t90 = sin(pkin(10));
t73 = t126 * t90;
t92 = cos(pkin(10));
t74 = t126 * t92;
t95 = sin(qJ(4));
t97 = cos(qJ(4));
t131 = -t97 * t73 - t74 * t95;
t89 = sin(pkin(11));
t91 = cos(pkin(11));
t94 = sin(qJ(6));
t96 = cos(qJ(6));
t99 = t89 * t94 - t91 * t96;
t69 = t99 * qJD(6);
t130 = 2 * m(6);
t129 = 2 * m(7);
t87 = t91 ^ 2;
t49 = -t95 * t73 + t74 * t97;
t78 = t90 * t97 + t95 * t92;
t32 = t78 * qJD(3) + t49 * qJD(4);
t128 = 0.2e1 * t32;
t72 = t78 * qJD(4);
t127 = pkin(4) * t72;
t76 = t90 * t95 - t97 * t92;
t71 = t76 * qJD(4);
t114 = t91 * t71;
t115 = t89 * t71;
t41 = -mrSges(6,1) * t115 - mrSges(6,2) * t114;
t77 = t89 * t96 + t91 * t94;
t70 = t77 * qJD(6);
t20 = -t78 * t70 + t99 * t71;
t21 = t78 * t69 + t77 * t71;
t7 = -t21 * mrSges(7,1) + t20 * mrSges(7,2);
t125 = t41 + t7;
t124 = Ifges(6,4) * t89;
t123 = Ifges(6,4) * t91;
t122 = Ifges(6,5) * t72;
t121 = Ifges(6,6) * t72;
t120 = t32 * t131;
t59 = t76 * t72;
t118 = t78 * t89;
t117 = t78 * t91;
t116 = t89 * Ifges(6,2);
t113 = -mrSges(6,1) * t91 + mrSges(6,2) * t89 - mrSges(5,1);
t112 = pkin(8) + qJ(5);
t31 = -t76 * qJD(3) + t131 * qJD(4);
t36 = qJ(5) * t71 - qJD(5) * t78 + t127;
t12 = t91 * t31 + t89 * t36;
t98 = -cos(pkin(9)) * pkin(1) - pkin(3) * t92 - pkin(2);
t43 = pkin(4) * t76 - qJ(5) * t78 + t98;
t23 = t89 * t43 + t91 * t49;
t111 = -Ifges(7,5) * t69 - Ifges(7,6) * t70;
t110 = t89 ^ 2 + t87;
t109 = Ifges(7,5) * t20 + Ifges(7,6) * t21 + Ifges(7,3) * t72;
t108 = mrSges(6,3) * t110;
t107 = t110 * t78;
t65 = t71 * mrSges(5,2);
t106 = t72 * mrSges(5,1) - t65;
t11 = -t31 * t89 + t91 * t36;
t22 = t91 * t43 - t49 * t89;
t105 = t110 * qJ(5);
t103 = Ifges(6,5) * t91 - Ifges(6,6) * t89;
t10 = pkin(5) * t76 - pkin(8) * t117 + t22;
t15 = -pkin(8) * t118 + t23;
t3 = t10 * t96 - t15 * t94;
t4 = t10 * t94 + t15 * t96;
t102 = -t11 * t89 + t12 * t91;
t101 = -t22 * t89 + t23 * t91;
t100 = -t131 * t72 + t32 * t76;
t79 = t112 * t89;
t81 = t112 * t91;
t60 = -t79 * t96 - t81 * t94;
t61 = -t79 * t94 + t81 * t96;
t84 = -pkin(5) * t91 - pkin(4);
t58 = Ifges(7,1) * t77 - Ifges(7,4) * t99;
t57 = Ifges(7,4) * t77 - Ifges(7,2) * t99;
t56 = mrSges(7,1) * t99 + mrSges(7,2) * t77;
t55 = mrSges(6,1) * t76 - mrSges(6,3) * t117;
t54 = -mrSges(6,2) * t76 - mrSges(6,3) * t118;
t53 = -Ifges(7,1) * t69 - Ifges(7,4) * t70;
t52 = -Ifges(7,4) * t69 - Ifges(7,2) * t70;
t51 = mrSges(7,1) * t70 - mrSges(7,2) * t69;
t50 = (mrSges(6,1) * t89 + mrSges(6,2) * t91) * t78;
t47 = mrSges(6,1) * t72 + mrSges(6,3) * t114;
t46 = -mrSges(6,2) * t72 + mrSges(6,3) * t115;
t45 = t99 * t78;
t44 = t77 * t78;
t38 = -t77 * qJD(5) - t61 * qJD(6);
t37 = -t99 * qJD(5) + t60 * qJD(6);
t33 = pkin(5) * t118 - t131;
t30 = mrSges(7,1) * t76 + mrSges(7,3) * t45;
t29 = -mrSges(7,2) * t76 - mrSges(7,3) * t44;
t27 = t122 - (t91 * Ifges(6,1) - t124) * t71;
t26 = t121 - (-t116 + t123) * t71;
t25 = -pkin(5) * t115 + t32;
t24 = mrSges(7,1) * t44 - mrSges(7,2) * t45;
t17 = -Ifges(7,1) * t45 - Ifges(7,4) * t44 + Ifges(7,5) * t76;
t16 = -Ifges(7,4) * t45 - Ifges(7,2) * t44 + Ifges(7,6) * t76;
t14 = -mrSges(7,2) * t72 + mrSges(7,3) * t21;
t13 = mrSges(7,1) * t72 - mrSges(7,3) * t20;
t9 = pkin(8) * t115 + t12;
t8 = pkin(5) * t72 + pkin(8) * t114 + t11;
t6 = Ifges(7,1) * t20 + Ifges(7,4) * t21 + Ifges(7,5) * t72;
t5 = Ifges(7,4) * t20 + Ifges(7,2) * t21 + Ifges(7,6) * t72;
t2 = -t4 * qJD(6) + t8 * t96 - t9 * t94;
t1 = t3 * qJD(6) + t8 * t94 + t9 * t96;
t18 = [0.2e1 * (m(4) * t82 + mrSges(4,3)) * qJD(3) * (t90 ^ 2 + t92 ^ 2) + 0.2e1 * t98 * t106 + 0.2e1 * m(5) * (t31 * t49 - t120) + (t11 * t22 + t12 * t23 - t120) * t130 + (mrSges(5,3) * t128 - t89 * t26 + t91 * t27 + (-(2 * Ifges(5,4)) + t103) * t72 - (Ifges(6,1) * t87 + (2 * Ifges(5,1)) + (t116 - 0.2e1 * t123) * t89) * t71) * t78 + t72 * (-Ifges(7,5) * t45 - Ifges(7,6) * t44) + (t1 * t4 + t2 * t3 + t25 * t33) * t129 + t50 * t128 + 0.2e1 * (t131 * t71 - t49 * t72) * mrSges(5,3) - 0.2e1 * t131 * t41 + (-0.2e1 * t31 * mrSges(5,3) + ((2 * Ifges(5,2)) + (2 * Ifges(6,3)) + Ifges(7,3)) * t72 + t109 - 0.2e1 * (-Ifges(5,4) + t103) * t71) * t76 + 0.2e1 * t3 * t13 + 0.2e1 * t4 * t14 + t20 * t17 + t21 * t16 + 0.2e1 * t25 * t24 + 0.2e1 * t1 * t29 + 0.2e1 * t2 * t30 + 0.2e1 * t33 * t7 - t44 * t5 - t45 * t6 + 0.2e1 * t23 * t46 + 0.2e1 * t22 * t47 + 0.2e1 * t12 * t54 + 0.2e1 * t11 * t55; -t44 * t13 - t45 * t14 + t20 * t29 + t21 * t30 + (t91 * t46 - t89 * t47) * t78 + t125 * t76 + (t24 + t50) * t72 - (t91 * t54 - t89 * t55) * t71 + m(7) * (-t1 * t45 - t2 * t44 + t20 * t4 + t21 * t3 + t25 * t76 + t33 * t72) + m(6) * (-t101 * t71 + t102 * t78 + t100) + m(5) * (t31 * t78 - t49 * t71 + t100); 0.2e1 * m(7) * (-t20 * t45 - t21 * t44 + t59) + 0.2e1 * m(6) * (-t71 * t107 + t59) + 0.2e1 * m(5) * (-t71 * t78 + t59); -t99 * t13 + t77 * t14 - t69 * t29 - t70 * t30 + t89 * t46 + t91 * t47 + m(7) * (t1 * t77 - t2 * t99 - t3 * t70 - t4 * t69) + m(6) * (t11 * t91 + t12 * t89) + t106; m(7) * (t20 * t77 - t21 * t99 + t44 * t70 + t45 * t69); (-t69 * t77 + t70 * t99) * t129; m(7) * (t1 * t61 + t2 * t60 + t25 * t84 + t3 * t38 + t37 * t4) + (qJ(5) * t46 + qJD(5) * t54 + t12 * mrSges(6,3) + t26 / 0.2e1 + t121 / 0.2e1) * t91 + (-qJ(5) * t47 - qJD(5) * t55 - t11 * mrSges(6,3) + t27 / 0.2e1 + t122 / 0.2e1) * t89 - (t91 * (Ifges(6,1) * t89 + t123) / 0.2e1 - t89 * (Ifges(6,2) * t91 + t124) / 0.2e1 + Ifges(5,5)) * t71 + t76 * t111 / 0.2e1 + t113 * t32 + m(6) * (-pkin(4) * t32 + t102 * qJ(5) + t101 * qJD(5)) + (-t1 * t99 - t2 * t77 + t3 * t69 - t4 * t70) * mrSges(7,3) + t72 * (Ifges(7,5) * t77 - Ifges(7,6) * t99) / 0.2e1 - t99 * t5 / 0.2e1 - t31 * mrSges(5,2) + t37 * t29 + t38 * t30 - pkin(4) * t41 + t33 * t51 - t44 * t52 / 0.2e1 - t45 * t53 / 0.2e1 + t25 * t56 + t21 * t57 / 0.2e1 + t20 * t58 / 0.2e1 + t60 * t13 + t61 * t14 - t69 * t17 / 0.2e1 - t70 * t16 / 0.2e1 - Ifges(5,6) * t72 + t77 * t6 / 0.2e1 + t84 * t7; t76 * t51 + t65 - t71 * t108 + (t56 + t113) * t72 + m(7) * (t20 * t61 + t21 * t60 - t37 * t45 - t38 * t44 + t72 * t84) + m(6) * (qJD(5) * t107 - t71 * t105 - t127) + (-t20 * t99 - t21 * t77 - t44 * t69 + t45 * t70) * mrSges(7,3); m(7) * (t37 * t77 - t38 * t99 - t60 * t70 - t61 * t69); 0.2e1 * t84 * t51 + (t37 * t61 + t38 * t60) * t129 - t69 * t58 + t77 * t53 - t70 * t57 - t99 * t52 + 0.2e1 * (-t37 * t99 - t38 * t77 + t60 * t69 - t61 * t70) * mrSges(7,3) + (t105 * t130 + 0.2e1 * t108) * qJD(5); m(6) * t32 + m(7) * t25 + t125; 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t72; 0; t51; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t109; -t7; -t51; mrSges(7,1) * t38 - t37 * mrSges(7,2) + t111; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t18(1) t18(2) t18(4) t18(7) t18(11) t18(16); t18(2) t18(3) t18(5) t18(8) t18(12) t18(17); t18(4) t18(5) t18(6) t18(9) t18(13) t18(18); t18(7) t18(8) t18(9) t18(10) t18(14) t18(19); t18(11) t18(12) t18(13) t18(14) t18(15) t18(20); t18(16) t18(17) t18(18) t18(19) t18(20) t18(21);];
Mq  = res;
