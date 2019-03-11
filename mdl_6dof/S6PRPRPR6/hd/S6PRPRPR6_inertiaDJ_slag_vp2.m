% Calculate time derivative of joint inertia matrix for
% S6PRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:46:44
% EndTime: 2019-03-08 19:46:47
% DurationCPUTime: 1.53s
% Computational Cost: add. (1743->315), mult. (4407->496), div. (0->0), fcn. (3881->10), ass. (0->135)
t105 = (-pkin(2) - pkin(8));
t150 = 2 * t105;
t100 = sin(qJ(4));
t103 = cos(qJ(4));
t156 = t100 * mrSges(5,1) + t103 * mrSges(5,2) + mrSges(4,3);
t97 = cos(pkin(11));
t94 = t97 ^ 2;
t95 = sin(pkin(11));
t135 = t95 ^ 2 + t94;
t101 = sin(qJ(2));
t96 = sin(pkin(6));
t133 = t101 * t96;
t119 = qJD(2) * t133;
t104 = cos(qJ(2));
t129 = t104 * t96;
t121 = t100 * t129;
t124 = qJD(4) * t103;
t98 = cos(pkin(6));
t44 = -qJD(4) * t121 - t103 * t119 + t124 * t98;
t132 = t103 * t44;
t125 = qJD(4) * t100;
t67 = t100 * t98 + t103 * t129;
t55 = t67 * t125;
t155 = t55 - t132;
t102 = cos(qJ(6));
t99 = sin(qJ(6));
t75 = t102 * t97 - t95 * t99;
t154 = t75 * qJD(6);
t153 = qJD(5) * m(6) * t135;
t43 = -qJD(4) * t67 + t100 * t119;
t68 = t103 * t98 - t121;
t152 = (t100 * t67 + t103 * t68) * qJD(4) + t100 * t43 - t132;
t151 = 2 * m(7);
t149 = m(6) / 0.2e1;
t148 = t75 / 0.2e1;
t76 = t102 * t95 + t97 * t99;
t147 = t76 / 0.2e1;
t146 = t95 / 0.2e1;
t113 = mrSges(6,1) * t95 + mrSges(6,2) * t97;
t61 = t113 * t125;
t66 = t76 * qJD(6);
t26 = -t103 * t66 - t125 * t75;
t28 = -t103 * t154 + t125 * t76;
t7 = -t28 * mrSges(7,1) + t26 * mrSges(7,2);
t145 = -t61 + t7;
t144 = Ifges(6,4) * t95;
t143 = Ifges(6,4) * t97;
t142 = Ifges(6,2) * t95;
t141 = pkin(4) * t100;
t19 = t67 * t44;
t140 = pkin(9) + qJ(5);
t57 = t76 * t103;
t59 = t75 * t103;
t30 = mrSges(7,1) * t57 + mrSges(7,2) * t59;
t69 = t113 * t103;
t139 = t30 + t69;
t138 = Ifges(7,5) * t154 - Ifges(7,6) * t66;
t137 = -mrSges(6,1) * t97 + mrSges(6,2) * t95 - mrSges(5,1);
t116 = t105 * t124;
t60 = -qJD(5) * t103 + qJD(3) + (pkin(4) * t103 + qJ(5) * t100) * qJD(4);
t36 = t97 * t116 + t95 * t60;
t126 = qJD(2) * t104;
t118 = t96 * t126;
t136 = qJ(3) * t118 + qJD(3) * t133;
t127 = t100 * t105;
t128 = qJ(5) * t103;
t80 = qJ(3) - t128 + t141;
t52 = t97 * t127 + t95 * t80;
t134 = t100 * t97;
t131 = t103 * t95;
t130 = t103 * t97;
t38 = -mrSges(7,1) * t75 + mrSges(7,2) * t76;
t123 = t38 + t137;
t122 = Ifges(7,5) * t26 + Ifges(7,6) * t28 + Ifges(7,3) * t124;
t117 = t135 * mrSges(6,3);
t115 = -t105 * t95 + pkin(5);
t114 = pkin(5) * t95 - t105;
t112 = -Ifges(6,5) * t97 + Ifges(6,6) * t95;
t20 = t118 * t97 - t43 * t95;
t21 = t118 * t95 + t43 * t97;
t111 = -t20 * t95 + t21 * t97;
t41 = t133 * t97 - t68 * t95;
t42 = t133 * t95 + t68 * t97;
t110 = -t41 * t95 + t42 * t97;
t73 = t97 * t80;
t34 = -pkin(9) * t130 + t100 * t115 + t73;
t37 = -pkin(9) * t131 + t52;
t9 = t102 * t37 + t34 * t99;
t8 = t102 * t34 - t37 * t99;
t11 = t102 * t42 + t41 * t99;
t10 = t102 * t41 - t42 * t99;
t83 = t140 * t95;
t85 = t140 * t97;
t46 = t102 * t85 - t83 * t99;
t45 = -t102 * t83 - t85 * t99;
t106 = m(5) * t152;
t91 = -pkin(5) * t97 - pkin(4);
t79 = (mrSges(5,1) * t103 - mrSges(5,2) * t100) * qJD(4);
t78 = mrSges(6,1) * t100 - mrSges(6,3) * t130;
t77 = -mrSges(6,2) * t100 - mrSges(6,3) * t131;
t74 = t114 * t103;
t71 = (mrSges(6,1) * t103 + mrSges(6,3) * t134) * qJD(4);
t70 = (mrSges(6,3) * t100 * t95 - mrSges(6,2) * t103) * qJD(4);
t64 = t114 * t125;
t58 = t75 * t100;
t56 = t76 * t100;
t54 = t97 * t60;
t51 = -t127 * t95 + t73;
t50 = (Ifges(6,5) * t103 + (-Ifges(6,1) * t97 + t144) * t100) * qJD(4);
t49 = (Ifges(6,6) * t103 + (t142 - t143) * t100) * qJD(4);
t48 = mrSges(7,1) * t100 - mrSges(7,3) * t59;
t47 = -mrSges(7,2) * t100 - mrSges(7,3) * t57;
t40 = Ifges(7,1) * t76 + Ifges(7,4) * t75;
t39 = Ifges(7,4) * t76 + Ifges(7,2) * t75;
t35 = -t116 * t95 + t54;
t33 = Ifges(7,1) * t154 - Ifges(7,4) * t66;
t32 = Ifges(7,4) * t154 - Ifges(7,2) * t66;
t31 = mrSges(7,1) * t66 + mrSges(7,2) * t154;
t29 = pkin(9) * t125 * t95 + t36;
t27 = -qJD(4) * t57 - t100 * t154;
t25 = qJD(4) * t59 - t100 * t66;
t18 = Ifges(7,1) * t59 - Ifges(7,4) * t57 + Ifges(7,5) * t100;
t17 = Ifges(7,4) * t59 - Ifges(7,2) * t57 + Ifges(7,6) * t100;
t16 = -qJD(5) * t76 - qJD(6) * t46;
t15 = qJD(5) * t75 + qJD(6) * t45;
t14 = t54 + (pkin(9) * t134 + t103 * t115) * qJD(4);
t13 = -mrSges(7,2) * t124 + mrSges(7,3) * t28;
t12 = mrSges(7,1) * t124 - mrSges(7,3) * t26;
t6 = Ifges(7,1) * t26 + Ifges(7,4) * t28 + Ifges(7,5) * t124;
t5 = Ifges(7,4) * t26 + Ifges(7,2) * t28 + Ifges(7,6) * t124;
t4 = -qJD(6) * t11 + t102 * t20 - t21 * t99;
t3 = qJD(6) * t10 + t102 * t21 + t20 * t99;
t2 = -qJD(6) * t9 + t102 * t14 - t29 * t99;
t1 = qJD(6) * t8 + t102 * t29 + t14 * t99;
t22 = [0.2e1 * m(6) * (t20 * t41 + t21 * t42 + t19) + 0.2e1 * m(7) * (t10 * t4 + t11 * t3 + t19) + 0.2e1 * m(5) * (t101 * t126 * t96 ^ 2 + t43 * t68 + t19); t10 * t12 + t11 * t13 + t20 * t78 + t21 * t77 + t3 * t47 + t4 * t48 + t41 * t71 + t42 * t70 + t145 * t67 + t139 * t44 - t152 * mrSges(5,3) + (t101 * t79 + ((-mrSges(3,1) + mrSges(4,2)) * t101 + (-mrSges(3,2) + t156) * t104) * qJD(2)) * t96 + m(5) * t136 + m(6) * (t51 * t20 + t52 * t21 + t35 * t41 + t36 * t42) + m(4) * (-pkin(2) * t119 + t136) + m(7) * (t1 * t11 + t10 * t2 + t3 * t9 + t4 * t8 + t44 * t74 - t64 * t67) + (t106 / 0.2e1 + t155 * t149) * t150; 0.2e1 * qJ(3) * t79 + 0.2e1 * t1 * t47 + 0.2e1 * t8 * t12 + 0.2e1 * t9 * t13 + t28 * t17 + t26 * t18 + 0.2e1 * t2 * t48 - 0.2e1 * t64 * t30 + 0.2e1 * t35 * t78 + 0.2e1 * t36 * t77 - t57 * t5 + 0.2e1 * t51 * t71 + 0.2e1 * t52 * t70 + t59 * t6 + 0.2e1 * t74 * t7 + (t1 * t9 + t2 * t8 - t64 * t74) * t151 + 0.2e1 * m(6) * (t51 * t35 + t52 * t36) + 0.2e1 * ((m(4) + m(5)) * qJ(3) + t156) * qJD(3) + (t61 * t150 - t95 * t49 + t97 * t50 + (Ifges(7,5) * t59 - Ifges(7,6) * t57 + (-(2 * Ifges(5,4)) - t112) * t103) * qJD(4)) * t103 + ((t69 * t150 + 0.2e1 * (Ifges(5,4) + t112) * t100 + ((2 * Ifges(6,3)) - (2 * Ifges(5,1)) + (2 * Ifges(5,2)) + Ifges(7,3) - 0.2e1 * m(6) * (t105 ^ 2) - Ifges(6,1) * t94 + (-t142 + 0.2e1 * t143) * t95) * t103) * qJD(4) + t122) * t100; m(4) * t119 + m(6) * (t55 + t111 * t100 + (qJD(4) * t110 - t44) * t103) + m(7) * (t10 * t27 + t11 * t25 + t3 * t58 - t4 * t56 + t155) + t106; -t56 * t12 + t58 * t13 + t25 * t47 + t27 * t48 - t145 * t103 + (t97 * t70 - t95 * t71) * t100 + ((t77 * t97 - t78 * t95) * t103 + t139 * t100) * qJD(4) + m(7) * (t1 * t58 + t103 * t64 + t125 * t74 - t2 * t56 + t25 * t9 + t27 * t8) + m(6) * ((-t35 * t95 + t36 * t97) * t100 + (-t51 * t95 + t52 * t97 - 0.2e1 * t127) * t124); (t25 * t58 - t27 * t56) * t151 + 0.4e1 * ((-0.1e1 + t135) * t149 - m(7) / 0.2e1) * t100 * t124; -t43 * mrSges(5,2) + t67 * t31 + t111 * mrSges(6,3) + t123 * t44 + m(6) * (-pkin(4) * t44 + qJ(5) * t111 + qJD(5) * t110) + m(7) * (t10 * t16 + t11 * t15 + t3 * t46 + t4 * t45 + t44 * t91) + (-t10 * t154 - t11 * t66 + t3 * t75 - t4 * t76) * mrSges(7,3); t100 * t138 / 0.2e1 + t91 * t7 + t74 * t31 + t5 * t148 + t6 * t147 - t57 * t32 / 0.2e1 + t59 * t33 / 0.2e1 + pkin(4) * t61 - t64 * t38 + t154 * t18 / 0.2e1 - t66 * t17 / 0.2e1 + t45 * t12 + t46 * t13 + t15 * t47 + t16 * t48 + t28 * t39 / 0.2e1 + t26 * t40 / 0.2e1 + m(7) * (t1 * t46 + t15 * t9 + t16 * t8 + t2 * t45 - t64 * t91) + (t49 / 0.2e1 + qJ(5) * t70 + qJD(5) * t77 + t36 * mrSges(6,3) + m(6) * (qJ(5) * t36 + qJD(5) * t52)) * t97 + (t50 / 0.2e1 - qJ(5) * t71 - qJD(5) * t78 - t35 * mrSges(6,3) + m(6) * (-qJ(5) * t35 - qJD(5) * t51)) * t95 + (t1 * t75 - t154 * t8 - t2 * t76 - t66 * t9) * mrSges(7,3) + ((Ifges(6,5) * t146 + Ifges(6,6) * t97 / 0.2e1 + Ifges(7,5) * t147 + Ifges(7,6) * t148 - Ifges(5,6) - t105 * mrSges(5,2)) * t103 + (-Ifges(5,5) + (Ifges(6,2) * t97 + t144) * t146 - t97 * (Ifges(6,1) * t95 + t143) / 0.2e1 + (-m(6) * pkin(4) + t137) * t105) * t100) * qJD(4); -t103 * t31 + m(7) * (t15 * t58 - t16 * t56 + t25 * t46 + t27 * t45) + t100 * t153 + ((-mrSges(5,2) + t117) * t103 + m(6) * (t128 * t135 - t141) + (m(7) * t91 + t123) * t100) * qJD(4) + (t154 * t56 + t25 * t75 - t27 * t76 - t58 * t66) * mrSges(7,3); t154 * t40 + t76 * t33 - t66 * t39 + t75 * t32 + 0.2e1 * t91 * t31 + (t15 * t46 + t16 * t45) * t151 + 0.2e1 * qJ(5) * t153 + 0.2e1 * t117 * qJD(5) + 0.2e1 * (t15 * t75 - t154 * t45 - t16 * t76 - t46 * t66) * mrSges(7,3); 0.2e1 * (t149 + m(7) / 0.2e1) * t44; -m(7) * t64 + (m(6) * t105 - t113) * t125 + t7; (m(6) + m(7)) * t125; t31; 0; mrSges(7,1) * t4 - mrSges(7,2) * t3; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t122; mrSges(7,1) * t27 - mrSges(7,2) * t25; mrSges(7,1) * t16 - t15 * mrSges(7,2) + t138; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t22(1) t22(2) t22(4) t22(7) t22(11) t22(16); t22(2) t22(3) t22(5) t22(8) t22(12) t22(17); t22(4) t22(5) t22(6) t22(9) t22(13) t22(18); t22(7) t22(8) t22(9) t22(10) t22(14) t22(19); t22(11) t22(12) t22(13) t22(14) t22(15) t22(20); t22(16) t22(17) t22(18) t22(19) t22(20) t22(21);];
Mq  = res;
