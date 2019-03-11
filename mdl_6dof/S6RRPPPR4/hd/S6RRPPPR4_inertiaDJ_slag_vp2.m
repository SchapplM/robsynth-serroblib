% Calculate time derivative of joint inertia matrix for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:17:19
% EndTime: 2019-03-09 08:17:23
% DurationCPUTime: 1.65s
% Computational Cost: add. (1752->320), mult. (3818->475), div. (0->0), fcn. (2904->6), ass. (0->138)
t107 = cos(pkin(9));
t94 = -t107 * qJD(5) + qJD(3);
t165 = 0.2e1 * t94;
t164 = Ifges(5,1) + Ifges(6,1);
t112 = cos(qJ(2));
t106 = sin(pkin(9));
t109 = sin(qJ(6));
t111 = cos(qJ(6));
t73 = t106 * t111 - t107 * t109;
t53 = t73 * t112;
t105 = t107 ^ 2;
t125 = (t106 ^ 2 + t105) * qJD(4);
t142 = qJ(5) * t106;
t154 = pkin(4) + pkin(5);
t163 = -t107 * t154 - t142;
t162 = 2 * m(5);
t161 = 2 * m(6);
t160 = 2 * m(7);
t159 = -2 * pkin(1);
t110 = sin(qJ(2));
t134 = qJD(2) * t110;
t126 = pkin(2) * t134 - t110 * qJD(3);
t133 = qJD(2) * t112;
t60 = -qJ(3) * t133 + t126;
t158 = 0.2e1 * t60;
t128 = -t110 * qJ(3) - pkin(1);
t85 = -pkin(2) * t112 + t128;
t157 = -0.2e1 * t85;
t156 = m(6) / 0.2e1;
t155 = pkin(3) + pkin(7);
t102 = t112 * pkin(7);
t108 = -pkin(2) - qJ(4);
t152 = pkin(8) + t108;
t117 = t109 * t106 + t111 * t107;
t25 = qJD(6) * t53 - t117 * t134;
t52 = t117 * t112;
t26 = qJD(6) * t52 + t134 * t73;
t151 = -Ifges(7,5) * t26 - Ifges(7,6) * t25;
t44 = -qJD(4) * t112 + (-qJ(3) * t112 + qJ(4) * t110) * qJD(2) + t126;
t80 = t155 * t133;
t18 = t106 * t80 + t107 * t44;
t71 = t108 * t112 + t128;
t88 = t155 * t110;
t35 = t106 * t88 + t107 * t71;
t61 = t73 * qJD(6);
t62 = t117 * qJD(6);
t150 = Ifges(7,5) * t61 - Ifges(7,6) * t62;
t141 = t106 * t110;
t67 = (mrSges(5,1) * t112 - mrSges(5,3) * t141) * qJD(2);
t131 = t106 * t134;
t68 = -mrSges(6,1) * t133 + mrSges(6,2) * t131;
t149 = t67 - t68;
t137 = t107 * t110;
t69 = (-mrSges(5,2) * t112 + mrSges(5,3) * t137) * qJD(2);
t70 = (mrSges(6,2) * t137 + mrSges(6,3) * t112) * qJD(2);
t148 = t69 + t70;
t147 = mrSges(6,3) * t106;
t146 = Ifges(5,4) * t106;
t145 = Ifges(5,4) * t107;
t144 = Ifges(6,5) * t106;
t143 = Ifges(6,5) * t107;
t139 = t106 * t112;
t136 = t107 * t112;
t89 = t112 * pkin(3) + t102;
t132 = 2 * mrSges(7,3);
t29 = t110 * qJ(5) + t35;
t130 = t107 * t134;
t129 = -m(4) * qJ(3) - mrSges(4,3);
t17 = -t106 * t44 + t107 * t80;
t34 = -t106 * t71 + t107 * t88;
t10 = qJ(5) * t133 + t110 * qJD(5) + t18;
t127 = qJ(5) * t107 - qJ(3);
t7 = -t25 * mrSges(7,1) + t26 * mrSges(7,2);
t30 = t62 * mrSges(7,1) + t61 * mrSges(7,2);
t123 = t117 * t61 - t62 * t73;
t122 = t108 * t125;
t121 = mrSges(6,1) * t107 + t147;
t120 = pkin(4) * t107 + t142;
t11 = -pkin(4) * t133 - t17;
t119 = t10 * t106 - t107 * t11;
t118 = t106 * t18 + t107 * t17;
t16 = pkin(8) * t139 - t110 * t154 - t34;
t22 = pkin(8) * t136 + t29;
t3 = -t109 * t22 + t111 * t16;
t4 = t109 * t16 + t111 * t22;
t81 = t152 * t106;
t82 = t152 * t107;
t40 = -t109 * t82 + t111 * t81;
t39 = -t109 * t81 - t111 * t82;
t8 = (-pkin(8) * t141 - t112 * t154) * qJD(2) - t17;
t9 = -pkin(8) * t130 + t10;
t1 = qJD(6) * t3 + t109 * t8 + t111 * t9;
t2 = -qJD(6) * t4 - t109 * t9 + t111 * t8;
t116 = t1 * t73 - t117 * t2 - t3 * t61 - t4 * t62;
t14 = -qJD(4) * t73 + qJD(6) * t39;
t15 = qJD(4) * t117 - qJD(6) * t40;
t115 = -t117 * t15 + t14 * t73 - t39 * t61 - t40 * t62;
t114 = (-Ifges(5,6) + Ifges(6,6)) * t107 + (-Ifges(6,4) - Ifges(5,5)) * t106;
t113 = -t109 * t62 - t111 * t61 + (t109 * t117 + t111 * t73) * qJD(6);
t92 = qJD(5) * t139;
t91 = mrSges(5,2) * t131;
t87 = mrSges(5,1) * t106 + mrSges(5,2) * t107;
t86 = mrSges(6,1) * t106 - mrSges(6,3) * t107;
t83 = pkin(4) * t106 - t127;
t79 = t155 * t134;
t78 = -mrSges(6,2) * t136 + t110 * mrSges(6,3);
t77 = -t110 * mrSges(5,2) - mrSges(5,3) * t136;
t76 = -t110 * mrSges(6,1) - mrSges(6,2) * t139;
t75 = t110 * mrSges(5,1) + mrSges(5,3) * t139;
t66 = (mrSges(5,1) * t107 - mrSges(5,2) * t106) * t112;
t65 = t121 * t112;
t63 = -t106 * t154 + t127;
t57 = -mrSges(5,1) * t130 + t91;
t56 = t121 * t134;
t51 = t112 * t120 + t89;
t50 = (t112 * Ifges(5,5) + (t106 * Ifges(5,1) + t145) * t110) * qJD(2);
t49 = (t112 * Ifges(6,4) + (t106 * Ifges(6,1) - t143) * t110) * qJD(2);
t48 = (t112 * Ifges(5,6) + (t107 * Ifges(5,2) + t146) * t110) * qJD(2);
t47 = (t112 * Ifges(6,6) + (-t107 * Ifges(6,3) + t144) * t110) * qJD(2);
t46 = -mrSges(7,1) * t110 + mrSges(7,3) * t53;
t45 = mrSges(7,2) * t110 + mrSges(7,3) * t52;
t43 = Ifges(7,1) * t117 + Ifges(7,4) * t73;
t42 = Ifges(7,4) * t117 + Ifges(7,2) * t73;
t41 = -mrSges(7,1) * t73 + mrSges(7,2) * t117;
t38 = t163 * t112 - t89;
t33 = -t110 * pkin(4) - t34;
t32 = Ifges(7,1) * t61 - Ifges(7,4) * t62;
t31 = Ifges(7,4) * t61 - Ifges(7,2) * t62;
t28 = t92 + (-t120 - t155) * t134;
t27 = -mrSges(7,1) * t52 - mrSges(7,2) * t53;
t21 = -Ifges(7,1) * t53 + Ifges(7,4) * t52 - Ifges(7,5) * t110;
t20 = -Ifges(7,4) * t53 + Ifges(7,2) * t52 - Ifges(7,6) * t110;
t19 = -t92 + (t155 - t163) * t134;
t13 = -mrSges(7,1) * t133 - t26 * mrSges(7,3);
t12 = mrSges(7,2) * t133 + t25 * mrSges(7,3);
t6 = Ifges(7,1) * t26 + Ifges(7,4) * t25 - Ifges(7,5) * t133;
t5 = Ifges(7,4) * t26 + Ifges(7,2) * t25 - Ifges(7,6) * t133;
t23 = [m(4) * t85 * t158 + 0.2e1 * t89 * t57 + 0.2e1 * t17 * t75 + 0.2e1 * t11 * t76 + 0.2e1 * t18 * t77 + 0.2e1 * t10 * t78 - 0.2e1 * t79 * t66 + 0.2e1 * t28 * t65 + 0.2e1 * t34 * t67 + 0.2e1 * t33 * t68 + 0.2e1 * t35 * t69 + 0.2e1 * t29 * t70 + t52 * t5 - t53 * t6 - 0.2e1 * t51 * t56 + 0.2e1 * t1 * t45 + 0.2e1 * t2 * t46 + t26 * t21 + 0.2e1 * t19 * t27 + 0.2e1 * t38 * t7 + t25 * t20 + 0.2e1 * t4 * t12 + 0.2e1 * t3 * t13 + (t1 * t4 + t19 * t38 + t2 * t3) * t160 + (t10 * t29 + t11 * t33 + t28 * t51) * t161 + (t17 * t34 + t18 * t35 - t79 * t89) * t162 + (((mrSges(3,2) * t159) + mrSges(4,3) * t157 + Ifges(7,5) * t53 - Ifges(7,6) * t52 + ((2 * Ifges(3,4)) + (2 * Ifges(4,6)) + t114) * t112) * t112 + (mrSges(4,2) * t157 + (mrSges(3,1) * t159) + 0.2e1 * (-Ifges(3,4) - Ifges(4,6) - t114) * t110 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) + (2 * Ifges(6,2)) - (2 * Ifges(4,3)) + (2 * Ifges(5,3)) + (2 * Ifges(7,3)) + (-Ifges(6,3) - Ifges(5,2)) * t105 + (-t164 * t106 + 0.2e1 * (-Ifges(5,4) + Ifges(6,5)) * t107) * t106) * t112) * t110) * qJD(2) + (mrSges(4,2) * t158 + (t47 - t48) * t107 + (-t49 - t50) * t106) * t112 + (-0.2e1 * t60 * mrSges(4,3) + t151) * t110; t117 * t6 / 0.2e1 + ((-Ifges(7,5) * t117 / 0.2e1 - Ifges(7,6) * t73 / 0.2e1 - Ifges(4,4) + Ifges(3,5) - pkin(2) * mrSges(4,1) + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t107 + (Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * t106 + (-m(4) * pkin(2) - mrSges(3,1) + mrSges(4,2)) * pkin(7)) * t112 + (Ifges(4,5) - Ifges(3,6) - t107 * (Ifges(6,3) * t106 + t143) / 0.2e1 + t107 * (-Ifges(5,2) * t106 + t145) / 0.2e1 - qJ(3) * mrSges(4,1) + (mrSges(3,2) + t129) * pkin(7) + (t164 * t107 + t144 - t146) * t106 / 0.2e1) * t110) * qJD(2) + t73 * t5 / 0.2e1 - t83 * t56 + t28 * t86 - t79 * t87 - t62 * t20 / 0.2e1 + t63 * t7 + t52 * t31 / 0.2e1 - t53 * t32 / 0.2e1 + qJ(3) * t57 + t61 * t21 / 0.2e1 + t25 * t42 / 0.2e1 + t26 * t43 / 0.2e1 + t14 * t45 + t15 * t46 + t38 * t30 + t39 * t13 + t40 * t12 + t19 * t41 + m(7) * (t1 * t40 + t14 * t4 + t15 * t3 + t19 * t63 + t2 * t39 - t38 * t94) + (m(4) * t102 + m(5) * t89 + t112 * mrSges(4,1) + t66) * qJD(3) + m(5) * (-qJ(3) * t79 + t118 * t108 + (-t106 * t35 - t107 * t34) * qJD(4)) + (-t18 * mrSges(5,3) - t10 * mrSges(6,2) + t47 / 0.2e1 - t48 / 0.2e1 + t148 * t108 + (-t77 - t78) * qJD(4)) * t106 + (-t17 * mrSges(5,3) + t11 * mrSges(6,2) + t49 / 0.2e1 + t50 / 0.2e1 + t149 * t108 + (-t75 + t76) * qJD(4)) * t107 - t110 * t150 / 0.2e1 + m(6) * (t28 * t83 + t51 * t94 + t119 * t108 + (-t106 * t29 + t107 * t33) * qJD(4)) + t116 * mrSges(7,3) + (t65 - t27) * t94; 0.2e1 * t63 * t30 + t73 * t31 + t117 * t32 - t62 * t42 + t61 * t43 + (t14 * t40 + t15 * t39 - t63 * t94) * t160 + (t83 * t94 - t122) * t161 + (qJ(3) * qJD(3) - t122) * t162 + t115 * t132 + 0.2e1 * (mrSges(6,2) + mrSges(5,3)) * t125 + (-t41 + t86) * t165 + 0.2e1 * (-t129 + t87) * qJD(3); t73 * t12 - t117 * t13 - t62 * t45 - t61 * t46 + t149 * t107 + t148 * t106 + (m(4) * pkin(7) + mrSges(4,1)) * t133 + m(7) * t116 + m(6) * t119 + m(5) * t118; t123 * t132 + m(7) * t115 - 0.2e1 * (t156 + m(5) / 0.2e1) * t125; t123 * t160; t91 + (-t147 + (-mrSges(5,1) - mrSges(6,1)) * t107) * t134 - m(5) * t79 + m(6) * t28 - m(7) * t19 - t7; m(5) * qJD(3) + (m(7) / 0.2e1 + t156) * t165 - t30; 0; 0; t109 * t12 + t111 * t13 + (-t109 * t46 + t111 * t45) * qJD(6) + m(7) * (t1 * t109 + t111 * t2 + (-t109 * t3 + t111 * t4) * qJD(6)) + m(6) * t11 + t68; m(7) * (t109 * t14 + t111 * t15 + (-t109 * t39 + t111 * t40) * qJD(6)) + m(6) * t107 * qJD(4) + t113 * mrSges(7,3); m(7) * t113; 0; 0; t2 * mrSges(7,1) - t1 * mrSges(7,2) - Ifges(7,3) * t133 - t151; mrSges(7,1) * t15 - mrSges(7,2) * t14 + t150; -mrSges(7,1) * t61 + mrSges(7,2) * t62; 0; (-mrSges(7,1) * t109 - mrSges(7,2) * t111) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t23(1) t23(2) t23(4) t23(7) t23(11) t23(16); t23(2) t23(3) t23(5) t23(8) t23(12) t23(17); t23(4) t23(5) t23(6) t23(9) t23(13) t23(18); t23(7) t23(8) t23(9) t23(10) t23(14) t23(19); t23(11) t23(12) t23(13) t23(14) t23(15) t23(20); t23(16) t23(17) t23(18) t23(19) t23(20) t23(21);];
Mq  = res;
