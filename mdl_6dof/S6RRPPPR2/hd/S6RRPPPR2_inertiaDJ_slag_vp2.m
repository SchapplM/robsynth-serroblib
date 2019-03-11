% Calculate time derivative of joint inertia matrix for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:09:59
% EndTime: 2019-03-09 08:10:03
% DurationCPUTime: 1.74s
% Computational Cost: add. (3129->292), mult. (6813->429), div. (0->0), fcn. (6397->8), ass. (0->133)
t162 = 2 * mrSges(5,1) + 2 * mrSges(4,3);
t105 = sin(pkin(10));
t107 = cos(pkin(10));
t109 = sin(qJ(6));
t111 = cos(qJ(6));
t118 = t111 * t105 + t109 * t107;
t160 = -t105 * t109 + t107 * t111;
t82 = t160 * qJD(6);
t83 = t118 * qJD(6);
t161 = t118 * t82 - t160 * t83;
t104 = t107 ^ 2;
t125 = (t105 ^ 2 + t104) * qJD(5);
t159 = 2 * m(6);
t158 = 2 * m(7);
t110 = sin(qJ(2));
t131 = qJD(2) * t110;
t102 = pkin(2) * t131;
t106 = sin(pkin(9));
t108 = cos(pkin(9));
t112 = cos(qJ(2));
t133 = t108 * t112;
t81 = qJD(2) * t133 - t106 * t131;
t86 = t106 * t112 + t108 * t110;
t115 = -qJ(4) * t81 - qJD(4) * t86 + t102;
t80 = t86 * qJD(2);
t38 = pkin(3) * t80 + t115;
t157 = -0.2e1 * t38;
t101 = -pkin(2) * t112 - pkin(1);
t116 = -t86 * qJ(4) + t101;
t84 = t106 * t110 - t133;
t60 = t84 * pkin(3) + t116;
t156 = -0.2e1 * t60;
t155 = 0.2e1 * t101;
t154 = m(4) * pkin(2);
t153 = -t118 / 0.2e1;
t152 = t160 / 0.2e1;
t151 = t107 / 0.2e1;
t148 = pkin(2) * t108;
t100 = -pkin(3) - t148;
t95 = -qJ(5) + t100;
t150 = -pkin(8) + t95;
t149 = pkin(2) * t106;
t145 = mrSges(5,2) - mrSges(4,1);
t144 = pkin(3) + qJ(5);
t143 = -qJ(3) - pkin(7);
t21 = qJD(5) * t84 + t144 * t80 + t115;
t126 = qJD(2) * t143;
t76 = qJD(3) * t112 + t110 * t126;
t77 = -t110 * qJD(3) + t112 * t126;
t46 = t106 * t76 - t108 * t77;
t34 = pkin(4) * t81 + t46;
t11 = t105 * t34 + t107 * t21;
t42 = t144 * t84 + t116;
t91 = t143 * t110;
t92 = t143 * t112;
t67 = -t106 * t92 - t108 * t91;
t48 = pkin(4) * t86 + t67;
t17 = t105 * t48 + t107 * t42;
t142 = -Ifges(7,5) * t83 - Ifges(7,6) * t82;
t141 = Ifges(6,1) * t105;
t140 = Ifges(6,4) * t105;
t139 = Ifges(6,4) * t107;
t138 = t105 * t80;
t137 = t107 * t80;
t136 = t107 * t84;
t130 = 2 * mrSges(7,3);
t129 = 0.2e1 * t112;
t50 = t160 * t84;
t25 = qJD(6) * t50 + t118 * t80;
t26 = t160 * t80 - t84 * t83;
t128 = Ifges(7,5) * t25 + Ifges(7,6) * t26 + Ifges(7,3) * t81;
t127 = -pkin(5) * t107 - pkin(4);
t9 = -t26 * mrSges(7,1) + t25 * mrSges(7,2);
t96 = qJ(4) + t149;
t45 = -mrSges(6,1) * t137 + mrSges(6,2) * t138;
t57 = mrSges(7,1) * t82 - mrSges(7,2) * t83;
t47 = t106 * t77 + t108 * t76;
t68 = t106 * t91 - t108 * t92;
t122 = t46 * t67 + t47 * t68;
t121 = Ifges(6,5) * t105 + Ifges(6,6) * t107;
t30 = t107 * t34;
t10 = -t105 * t21 + t30;
t120 = t10 * t107 + t105 * t11;
t44 = t107 * t48;
t12 = pkin(5) * t86 + t44 + (-pkin(8) * t84 - t42) * t105;
t13 = pkin(8) * t136 + t17;
t3 = -t109 * t13 + t111 * t12;
t4 = t109 * t12 + t111 * t13;
t78 = t150 * t105;
t79 = t150 * t107;
t53 = t109 * t79 + t111 * t78;
t52 = -t109 * t78 + t111 * t79;
t119 = t161 * t158;
t7 = pkin(5) * t81 + t30 + (-pkin(8) * t80 - t21) * t105;
t8 = pkin(8) * t137 + t11;
t1 = t3 * qJD(6) + t109 * t7 + t111 * t8;
t2 = -t4 * qJD(6) - t109 * t8 + t111 * t7;
t114 = -t1 * t118 - t160 * t2 + t3 * t83 - t4 * t82;
t36 = -t118 * qJD(5) + t52 * qJD(6);
t37 = -qJD(5) * t160 - t53 * qJD(6);
t113 = -t118 * t36 - t160 * t37 + t52 * t83 - t53 * t82;
t90 = mrSges(6,1) * t105 + mrSges(6,2) * t107;
t89 = pkin(5) * t105 + t96;
t72 = t81 * mrSges(5,3);
t71 = t81 * mrSges(4,2);
t65 = Ifges(7,1) * t160 - Ifges(7,4) * t118;
t64 = Ifges(7,4) * t160 - Ifges(7,2) * t118;
t63 = mrSges(7,1) * t118 + mrSges(7,2) * t160;
t62 = -mrSges(6,2) * t86 + mrSges(6,3) * t136;
t61 = -mrSges(6,3) * t105 * t84 + mrSges(6,1) * t86;
t59 = -Ifges(7,1) * t83 - Ifges(7,4) * t82;
t58 = -Ifges(7,4) * t83 - Ifges(7,2) * t82;
t56 = (-mrSges(6,1) * t107 + mrSges(6,2) * t105) * t84;
t55 = -mrSges(6,2) * t81 + mrSges(6,3) * t137;
t54 = mrSges(6,1) * t81 - mrSges(6,3) * t138;
t51 = t118 * t84;
t49 = -pkin(4) * t84 + t68;
t40 = mrSges(7,1) * t86 - mrSges(7,3) * t51;
t39 = -mrSges(7,2) * t86 + mrSges(7,3) * t50;
t35 = -pkin(4) * t80 + t47;
t33 = t127 * t84 + t68;
t32 = t81 * Ifges(6,5) + (t139 + t141) * t80;
t31 = t81 * Ifges(6,6) + (Ifges(6,2) * t107 + t140) * t80;
t28 = t127 * t80 + t47;
t27 = -mrSges(7,1) * t50 + mrSges(7,2) * t51;
t19 = Ifges(7,1) * t51 + Ifges(7,4) * t50 + Ifges(7,5) * t86;
t18 = Ifges(7,4) * t51 + Ifges(7,2) * t50 + Ifges(7,6) * t86;
t16 = -t105 * t42 + t44;
t15 = -mrSges(7,2) * t81 + mrSges(7,3) * t26;
t14 = mrSges(7,1) * t81 - mrSges(7,3) * t25;
t6 = Ifges(7,1) * t25 + Ifges(7,4) * t26 + t81 * Ifges(7,5);
t5 = Ifges(7,4) * t25 + Ifges(7,2) * t26 + t81 * Ifges(7,6);
t20 = [((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t112) * t129 + (0.2e1 * pkin(2) * (mrSges(4,1) * t84 + mrSges(4,2) * t86) + t154 * t155 - 0.2e1 * pkin(1) * mrSges(3,1) - 0.2e1 * Ifges(3,4) * t110 + (Ifges(3,1) - Ifges(3,2)) * t129) * t110) * qJD(2) + 0.2e1 * m(5) * (t38 * t60 + t122) + 0.2e1 * m(4) * t122 + (mrSges(5,3) * t157 + t162 * t46 + t128) * t86 + (mrSges(4,1) * t155 + mrSges(5,2) * t156 - t68 * t162 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6) + t121) * t86 + (Ifges(6,2) * t104 + (2 * Ifges(4,2)) + (2 * Ifges(5,3)) + (0.2e1 * t139 + t141) * t105) * t84) * t80 + (mrSges(5,2) * t157 + t105 * t32 + t107 * t31 - t162 * t47) * t84 + t71 * t155 + t72 * t156 + (t1 * t4 + t2 * t3 + t28 * t33) * t158 + (t10 * t16 + t11 * t17 + t35 * t49) * t159 + (Ifges(7,5) * t51 + Ifges(7,6) * t50 + t67 * t162 + ((2 * Ifges(4,1)) + (2 * Ifges(5,2)) + (2 * Ifges(6,3)) + Ifges(7,3)) * t86 + (-0.2e1 * Ifges(4,4) - 0.2e1 * Ifges(5,6) + t121) * t84) * t81 + 0.2e1 * t3 * t14 + 0.2e1 * t4 * t15 + t25 * t19 + t26 * t18 + 0.2e1 * t28 * t27 + 0.2e1 * t33 * t9 + 0.2e1 * t1 * t39 + 0.2e1 * t2 * t40 + 0.2e1 * t49 * t45 + t50 * t5 + t51 * t6 + 0.2e1 * t16 * t54 + 0.2e1 * t17 * t55 + 0.2e1 * t35 * t56 + 0.2e1 * t10 * t61 + 0.2e1 * t11 * t62; (-t84 * mrSges(5,1) + t27 + t56) * qJD(4) + m(7) * (qJD(4) * t33 + t1 * t53 + t2 * t52 + t28 * t89 + t37 * t3 + t36 * t4) + m(5) * (qJD(4) * t68 + t100 * t46 + t47 * t96) + (-qJD(5) * t61 - t10 * mrSges(6,3) + t95 * t54 + t32 / 0.2e1) * t107 + (-qJD(5) * t62 - t11 * mrSges(6,3) + t95 * t55 - t31 / 0.2e1) * t105 + (-mrSges(4,3) * t148 + t100 * mrSges(5,1) + Ifges(4,5) - Ifges(5,4) + Ifges(7,5) * t152 + Ifges(7,6) * t153 + Ifges(6,5) * t151 - Ifges(6,6) * t105 / 0.2e1) * t81 + (-mrSges(4,3) * t149 - t96 * mrSges(5,1) + t105 * (Ifges(6,1) * t107 - t140) / 0.2e1 + (-Ifges(6,2) * t105 + t139) * t151 + Ifges(5,5) - Ifges(4,6)) * t80 + t86 * t142 / 0.2e1 + t145 * t46 + m(6) * (qJD(4) * t49 + t35 * t96 + t120 * t95 + (-t105 * t17 - t107 * t16) * qJD(5)) + t6 * t152 + t5 * t153 + (t106 * t47 - t108 * t46) * t154 + t114 * mrSges(7,3) + (Ifges(3,5) * t112 - Ifges(3,6) * t110 + (-mrSges(3,1) * t112 + mrSges(3,2) * t110) * pkin(7)) * qJD(2) + (-mrSges(4,2) + mrSges(5,3)) * t47 + t36 * t39 + t37 * t40 + t52 * t14 + t53 * t15 + t33 * t57 + t50 * t58 / 0.2e1 + t51 * t59 / 0.2e1 + t28 * t63 + t26 * t64 / 0.2e1 + t25 * t65 / 0.2e1 - t82 * t18 / 0.2e1 - t83 * t19 / 0.2e1 + t89 * t9 + t35 * t90 + t96 * t45; 0.2e1 * t89 * t57 - t118 * t58 + t160 * t59 - t82 * t64 - t83 * t65 + (qJD(4) * t89 + t36 * t53 + t37 * t52) * t158 + (qJD(4) * t96 - t95 * t125) * t159 + t113 * t130 + 0.2e1 * mrSges(6,3) * t125 + 0.2e1 * (m(5) * t96 + mrSges(5,3) + t63 + t90) * qJD(4); m(4) * t102 - t105 * t54 + t107 * t55 - t118 * t14 + t160 * t15 - t83 * t39 - t82 * t40 + t71 - t72 - t145 * t80 + m(7) * (t1 * t160 - t118 * t2 - t3 * t82 - t4 * t83) + m(6) * (-t10 * t105 + t107 * t11) + m(5) * t38; m(7) * (-t118 * t37 + t160 * t36 - t52 * t82 - t53 * t83); t119; m(5) * t46 + m(6) * t120 - m(7) * t114 + t81 * mrSges(5,1) + t105 * t55 + t107 * t54 + t118 * t15 + t14 * t160 + t82 * t39 - t83 * t40; -m(6) * t125 - m(7) * t113 - t161 * t130; 0; t119; m(6) * t35 + m(7) * t28 + t45 + t9; (m(6) + m(7)) * qJD(4) + t57; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t128; mrSges(7,1) * t37 - mrSges(7,2) * t36 + t142; -t57; -mrSges(7,1) * t83 - t82 * mrSges(7,2); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;
