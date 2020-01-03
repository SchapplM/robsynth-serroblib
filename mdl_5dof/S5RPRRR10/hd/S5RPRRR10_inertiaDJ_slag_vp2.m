% Calculate time derivative of joint inertia matrix for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR10_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR10_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR10_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:09:45
% EndTime: 2019-12-31 19:09:50
% DurationCPUTime: 1.94s
% Computational Cost: add. (3154->278), mult. (7131->430), div. (0->0), fcn. (6831->8), ass. (0->125)
t104 = cos(qJ(4));
t102 = sin(qJ(3));
t105 = cos(qJ(3));
t98 = sin(pkin(9));
t99 = cos(pkin(9));
t78 = t102 * t98 - t105 * t99;
t74 = t78 * qJD(3);
t128 = t104 * t74;
t79 = t102 * t99 + t105 * t98;
t75 = t79 * qJD(3);
t151 = -Ifges(5,5) * t128 + Ifges(5,3) * t75;
t100 = sin(qJ(5));
t101 = sin(qJ(4));
t103 = cos(qJ(5));
t109 = t100 * t101 - t103 * t104;
t50 = t109 * t79;
t135 = pkin(6) + qJ(2);
t87 = t135 * t98;
t88 = t135 * t99;
t150 = -t102 * t88 - t105 * t87;
t113 = mrSges(5,1) * t101 + mrSges(5,2) * t104;
t82 = t113 * qJD(4);
t149 = qJD(4) + qJD(5);
t148 = 2 * m(6);
t147 = -2 * mrSges(4,3);
t145 = -0.2e1 * t150;
t144 = m(6) * pkin(4);
t142 = -t79 / 0.2e1;
t133 = Ifges(5,4) * t101;
t89 = Ifges(5,2) * t104 + t133;
t141 = -t89 / 0.2e1;
t140 = -pkin(8) - pkin(7);
t139 = Ifges(5,5) * t75;
t138 = Ifges(5,6) * t75;
t67 = -t102 * t87 + t105 * t88;
t44 = t79 * qJD(2) + t67 * qJD(3);
t137 = t44 * t150;
t136 = t78 * Ifges(5,6);
t61 = t149 * t109;
t81 = t100 * t104 + t101 * t103;
t62 = t149 * t81;
t134 = -Ifges(6,5) * t61 - Ifges(6,6) * t62;
t120 = -pkin(2) * t99 - pkin(1);
t54 = pkin(3) * t78 - pkin(7) * t79 + t120;
t60 = t104 * t67;
t29 = t101 * t54 + t60;
t132 = Ifges(5,4) * t104;
t131 = Ifges(5,6) * t101;
t130 = t101 * t79;
t127 = t104 * t79;
t126 = qJD(4) * t101;
t125 = qJD(4) * t104;
t124 = qJD(5) * t100;
t123 = qJD(5) * t103;
t17 = t109 * t74 - t62 * t79;
t18 = t149 * t50 + t81 * t74;
t122 = Ifges(6,5) * t17 + Ifges(6,6) * t18 + Ifges(6,3) * t75;
t121 = pkin(4) * t126;
t119 = t79 * t126;
t118 = qJD(4) * t140;
t117 = t75 * mrSges(4,1) - t74 * mrSges(4,2);
t30 = t62 * mrSges(6,1) - t61 * mrSges(6,2);
t116 = -(2 * Ifges(4,4)) - t131;
t43 = -t78 * qJD(2) + qJD(3) * t150;
t53 = pkin(3) * t75 + pkin(7) * t74;
t115 = -t101 * t43 + t104 * t53;
t28 = -t101 * t67 + t104 * t54;
t112 = Ifges(5,1) * t104 - t133;
t111 = -Ifges(5,2) * t101 + t132;
t19 = pkin(4) * t78 - pkin(8) * t127 + t28;
t24 = -pkin(8) * t130 + t29;
t9 = -t100 * t24 + t103 * t19;
t10 = t100 * t19 + t103 * t24;
t91 = t140 * t101;
t92 = t140 * t104;
t68 = t100 * t92 + t103 * t91;
t69 = t100 * t91 - t103 * t92;
t85 = t101 * t118;
t86 = t104 * t118;
t34 = qJD(5) * t68 + t100 * t86 + t103 * t85;
t35 = -qJD(5) * t69 - t100 * t85 + t103 * t86;
t110 = t35 * mrSges(6,1) - t34 * mrSges(6,2) + t134;
t7 = pkin(8) * t128 + pkin(4) * t75 + (-t60 + (pkin(8) * t79 - t54) * t101) * qJD(4) + t115;
t107 = -t101 * t74 + t79 * t125;
t11 = t101 * t53 + t104 * t43 + t54 * t125 - t67 * t126;
t8 = -pkin(8) * t107 + t11;
t2 = qJD(5) * t9 + t100 * t7 + t103 * t8;
t3 = -qJD(5) * t10 - t100 * t8 + t103 * t7;
t108 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + t122;
t106 = t119 + t128;
t95 = Ifges(5,5) * t125;
t94 = -pkin(4) * t104 - pkin(3);
t90 = Ifges(5,1) * t101 + t132;
t84 = t112 * qJD(4);
t83 = t111 * qJD(4);
t76 = (-mrSges(6,1) * t100 - mrSges(6,2) * t103) * qJD(5) * pkin(4);
t65 = Ifges(6,1) * t81 - Ifges(6,4) * t109;
t64 = Ifges(6,4) * t81 - Ifges(6,2) * t109;
t63 = mrSges(6,1) * t109 + mrSges(6,2) * t81;
t56 = mrSges(5,1) * t78 - mrSges(5,3) * t127;
t55 = -mrSges(5,2) * t78 - mrSges(5,3) * t130;
t49 = t81 * t79;
t45 = pkin(4) * t130 - t150;
t41 = Ifges(5,5) * t78 + t112 * t79;
t40 = t111 * t79 + t136;
t39 = mrSges(6,1) * t78 + mrSges(6,3) * t50;
t38 = -mrSges(6,2) * t78 - mrSges(6,3) * t49;
t37 = -mrSges(5,2) * t75 - t107 * mrSges(5,3);
t36 = mrSges(5,1) * t75 + mrSges(5,3) * t106;
t32 = -Ifges(6,1) * t61 - Ifges(6,4) * t62;
t31 = -Ifges(6,4) * t61 - Ifges(6,2) * t62;
t27 = mrSges(5,1) * t107 - mrSges(5,2) * t106;
t26 = mrSges(6,1) * t49 - mrSges(6,2) * t50;
t25 = pkin(4) * t107 + t44;
t23 = -Ifges(6,1) * t50 - Ifges(6,4) * t49 + Ifges(6,5) * t78;
t22 = -Ifges(6,4) * t50 - Ifges(6,2) * t49 + Ifges(6,6) * t78;
t21 = -Ifges(5,1) * t106 - Ifges(5,4) * t107 + t139;
t20 = -Ifges(5,4) * t106 - Ifges(5,2) * t107 + t138;
t14 = -mrSges(6,2) * t75 + mrSges(6,3) * t18;
t13 = mrSges(6,1) * t75 - mrSges(6,3) * t17;
t12 = -qJD(4) * t29 + t115;
t6 = -mrSges(6,1) * t18 + mrSges(6,2) * t17;
t5 = Ifges(6,1) * t17 + Ifges(6,4) * t18 + Ifges(6,5) * t75;
t4 = Ifges(6,4) * t17 + Ifges(6,2) * t18 + Ifges(6,6) * t75;
t1 = [t18 * t22 + t17 * t23 + 0.2e1 * t9 * t13 + 0.2e1 * t10 * t14 - (mrSges(4,3) * t145 - t101 * t40 + t104 * t41) * t74 + 0.2e1 * m(5) * (t11 * t29 + t12 * t28 - t137) + 0.2e1 * m(4) * (t43 * t67 - t137) + 0.2e1 * t120 * t117 + t27 * t145 + (t10 * t2 + t25 * t45 + t3 * t9) * t148 + (-t116 * t74 + t43 * t147 + t122 + t151) * t78 + (-t101 * t20 + t104 * t21 - 0.2e1 * Ifges(4,1) * t74 + (-t101 * t41 - t104 * t40 + t78 * (-Ifges(5,5) * t101 - Ifges(5,6) * t104)) * qJD(4) + 0.2e1 * (mrSges(4,3) + t113) * t44) * t79 + 0.2e1 * t25 * t26 + 0.2e1 * t28 * t36 + 0.2e1 * t29 * t37 + 0.2e1 * t2 * t38 + 0.2e1 * t3 * t39 + 0.2e1 * t45 * t6 - t49 * t4 - t50 * t5 + 0.2e1 * t11 * t55 + 0.2e1 * t12 * t56 + (((2 * Ifges(4,2)) + Ifges(5,3) + Ifges(6,3)) * t78 - Ifges(6,5) * t50 - Ifges(6,6) * t49 + (Ifges(5,5) * t104 + t116) * t79 + t67 * t147) * t75 + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t98 ^ 2 + t99 ^ 2) * qJD(2); t101 * t37 + t104 * t36 - t109 * t13 + t81 * t14 - t61 * t38 - t62 * t39 + (-t101 * t56 + t104 * t55) * qJD(4) + m(6) * (-t10 * t61 - t109 * t3 + t2 * t81 - t62 * t9) + m(5) * (t101 * t11 + t104 * t12 + (-t101 * t28 + t104 * t29) * qJD(4)) + t117; (t109 * t62 - t61 * t81) * t148; (-t10 * t62 - t109 * t2 - t3 * t81 + t9 * t61) * mrSges(6,3) + (t134 + t95) * t78 / 0.2e1 + t75 * (Ifges(6,5) * t81 - Ifges(6,6) * t109) / 0.2e1 - t109 * t4 / 0.2e1 + (-m(5) * t44 - t27) * pkin(3) - t150 * t82 + (-t74 * t90 / 0.2e1 + t79 * t84 / 0.2e1 + t11 * mrSges(5,3) + t20 / 0.2e1 - t44 * mrSges(5,1) + t138 / 0.2e1 + (t79 * t141 - t28 * mrSges(5,3) + t41 / 0.2e1) * qJD(4) + (m(5) * (-t28 * qJD(4) + t11) + t37 - qJD(4) * t56) * pkin(7)) * t104 + (-t74 * t141 + t83 * t142 - t12 * mrSges(5,3) + t21 / 0.2e1 + t44 * mrSges(5,2) + t139 / 0.2e1 + (-m(5) * t12 - t36) * pkin(7) + (-t40 / 0.2e1 + t90 * t142 - t29 * mrSges(5,3) + pkin(4) * t26 - t136 / 0.2e1 + t45 * t144 + (-m(5) * t29 - t55) * pkin(7)) * qJD(4)) * t101 + t34 * t38 + t35 * t39 - t43 * mrSges(4,2) - t44 * mrSges(4,1) + t45 * t30 - t49 * t31 / 0.2e1 - t50 * t32 / 0.2e1 - t61 * t23 / 0.2e1 - t62 * t22 / 0.2e1 + t25 * t63 + t18 * t64 / 0.2e1 + t17 * t65 / 0.2e1 + t68 * t13 + t69 * t14 - Ifges(4,5) * t74 - Ifges(4,6) * t75 + t81 * t5 / 0.2e1 + t94 * t6 + m(6) * (t10 * t34 + t2 * t69 + t25 * t94 + t3 * t68 + t35 * t9); m(6) * (-t109 * t35 + t34 * t81 - t61 * t69 - t62 * t68); (t121 * t94 + t34 * t69 + t35 * t68) * t148 - t61 * t65 + t81 * t32 - t62 * t64 - t109 * t31 + 0.2e1 * t63 * t121 + 0.2e1 * t94 * t30 - t89 * t126 - 0.2e1 * pkin(3) * t82 + t101 * t84 + (qJD(4) * t90 + t83) * t104 + 0.2e1 * (-t109 * t34 - t35 * t81 + t61 * t68 - t62 * t69) * mrSges(6,3); -Ifges(5,5) * t119 + t12 * mrSges(5,1) - t11 * mrSges(5,2) - t107 * Ifges(5,6) + (t38 * t123 + t100 * t14 - t39 * t124 + t103 * t13 + m(6) * (t10 * t123 + t100 * t2 + t103 * t3 - t124 * t9)) * pkin(4) + t108 + t151; -t82 + (-t100 * t61 - t103 * t62 + (t100 * t109 + t103 * t81) * qJD(5)) * t144 - t30; t95 + (-t131 + (-t104 * mrSges(5,1) + t101 * mrSges(5,2)) * pkin(7)) * qJD(4) + (m(6) * (t100 * t34 + t103 * t35 + (-t100 * t68 + t103 * t69) * qJD(5)) + (-t100 * t62 + t103 * t61 + (t100 * t81 - t103 * t109) * qJD(5)) * mrSges(6,3)) * pkin(4) + t110; 0.2e1 * t76; t108; -t30; t110; t76; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
