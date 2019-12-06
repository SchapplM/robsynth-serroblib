% Calculate time derivative of joint inertia matrix for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP8_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP8_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:16
% EndTime: 2019-12-05 16:58:20
% DurationCPUTime: 1.57s
% Computational Cost: add. (963->271), mult. (2719->397), div. (0->0), fcn. (2098->8), ass. (0->126)
t125 = Ifges(6,4) + Ifges(5,5);
t151 = -Ifges(6,2) - Ifges(5,3);
t82 = cos(qJ(3));
t116 = qJD(3) * t82;
t78 = sin(qJ(4));
t106 = t78 * t116;
t81 = cos(qJ(4));
t113 = qJD(4) * t81;
t79 = sin(qJ(3));
t87 = t79 * t113 + t106;
t150 = -Ifges(5,6) * t81 - t125 * t78;
t137 = Ifges(6,5) * t78;
t97 = Ifges(6,1) * t81 + t137;
t139 = Ifges(5,4) * t78;
t98 = Ifges(5,1) * t81 - t139;
t149 = (t97 + t98) * qJD(4);
t112 = qJD(4) * t82;
t115 = qJD(4) * t78;
t117 = qJD(3) * t79;
t53 = (pkin(3) * t79 - pkin(8) * t82) * qJD(3);
t57 = -pkin(3) * t82 - pkin(8) * t79 - pkin(2);
t9 = (-t81 * t112 + t78 * t117) * pkin(7) - t57 * t115 + t53 * t81;
t148 = m(6) * qJ(5) + mrSges(6,3);
t111 = qJD(5) * t81;
t94 = pkin(4) * t81 + qJ(5) * t78;
t147 = t94 * qJD(4) - t111;
t146 = 2 * m(5);
t145 = 0.2e1 * pkin(7);
t144 = m(5) / 0.2e1;
t142 = pkin(7) * t82;
t76 = sin(pkin(5));
t80 = sin(qJ(2));
t130 = t76 * t80;
t109 = qJD(2) * t130;
t83 = cos(qJ(2));
t129 = t76 * t83;
t110 = t78 * t129;
t118 = qJD(2) * t83;
t108 = t76 * t118;
t77 = cos(pkin(5));
t35 = t79 * t130 - t77 * t82;
t18 = -t35 * qJD(3) + t82 * t108;
t36 = t82 * t130 + t77 * t79;
t3 = -qJD(4) * t110 - t81 * t109 + t36 * t113 + t18 * t78;
t141 = t3 * t78;
t19 = t81 * t129 + t36 * t78;
t4 = -t19 * qJD(4) + t78 * t109 + t18 * t81;
t140 = t4 * t81;
t138 = Ifges(5,4) * t81;
t136 = Ifges(6,5) * t81;
t134 = Ifges(5,6) * t82;
t17 = t36 * qJD(3) + t79 * t108;
t133 = t17 * t79;
t132 = t18 * t82;
t10 = t35 * t17;
t131 = t57 * t81;
t128 = t78 * t79;
t127 = t79 * t81;
t105 = t81 * t116;
t114 = qJD(4) * t79;
t88 = -t78 * t114 + t105;
t23 = mrSges(5,1) * t117 - t88 * mrSges(5,3);
t24 = mrSges(6,2) * t105 + (-mrSges(6,1) * qJD(3) - mrSges(6,2) * t115) * t79;
t124 = -t23 + t24;
t25 = -mrSges(5,2) * t117 - t87 * mrSges(5,3);
t26 = -t87 * mrSges(6,2) + mrSges(6,3) * t117;
t123 = t25 + t26;
t31 = -t82 * Ifges(6,4) + t97 * t79;
t32 = -t82 * Ifges(5,5) + t98 * t79;
t122 = t31 + t32;
t49 = mrSges(5,2) * t82 - mrSges(5,3) * t128;
t52 = -mrSges(6,2) * t128 - mrSges(6,3) * t82;
t121 = t49 + t52;
t50 = -mrSges(5,1) * t82 - mrSges(5,3) * t127;
t51 = mrSges(6,1) * t82 + mrSges(6,2) * t127;
t120 = -t50 + t51;
t59 = -t81 * mrSges(5,1) + t78 * mrSges(5,2);
t119 = t59 - mrSges(4,1);
t28 = t81 * t142 + t78 * t57;
t60 = -Ifges(6,3) * t81 + t137;
t61 = Ifges(5,2) * t81 + t139;
t104 = t60 / 0.2e1 - t61 / 0.2e1;
t62 = Ifges(6,1) * t78 - t136;
t63 = Ifges(5,1) * t78 + t138;
t103 = t62 / 0.2e1 + t63 / 0.2e1;
t95 = Ifges(6,3) * t78 + t136;
t29 = -Ifges(6,6) * t82 + t95 * t79;
t96 = -Ifges(5,2) * t78 + t138;
t30 = t96 * t79 - t134;
t101 = t29 - t30 + t134;
t100 = mrSges(5,1) * t78 + mrSges(5,2) * t81;
t58 = -t81 * mrSges(6,1) - t78 * mrSges(6,3);
t99 = mrSges(6,1) * t78 - mrSges(6,3) * t81;
t93 = pkin(4) * t78 - qJ(5) * t81;
t92 = -t87 * Ifges(6,6) - t125 * t105 + t151 * t117;
t91 = pkin(7) + t93;
t90 = t35 * t116 + t133;
t8 = t78 * t53 + t57 * t113 + (-t78 * t112 - t81 * t117) * pkin(7);
t75 = Ifges(6,4) * t113;
t74 = Ifges(5,5) * t113;
t72 = Ifges(6,6) * t115;
t54 = -pkin(3) - t94;
t46 = t96 * qJD(4);
t45 = t95 * qJD(4);
t44 = (mrSges(4,1) * t79 + mrSges(4,2) * t82) * qJD(3);
t43 = t100 * qJD(4);
t42 = t99 * qJD(4);
t39 = t100 * t79;
t38 = t99 * t79;
t34 = t93 * qJD(4) - qJD(5) * t78;
t33 = t91 * t79;
t27 = -t78 * t142 + t131;
t22 = -t131 + (pkin(7) * t78 + pkin(4)) * t82;
t21 = -qJ(5) * t82 + t28;
t20 = t36 * t81 - t110;
t16 = t87 * mrSges(5,1) + t88 * mrSges(5,2);
t15 = t87 * mrSges(6,1) - t88 * mrSges(6,3);
t14 = -t63 * t114 + (Ifges(5,5) * t79 + t98 * t82) * qJD(3);
t13 = -t62 * t114 + (Ifges(6,4) * t79 + t97 * t82) * qJD(3);
t12 = -t61 * t114 + (Ifges(5,6) * t79 + t96 * t82) * qJD(3);
t11 = -t60 * t114 + (Ifges(6,6) * t79 + t95 * t82) * qJD(3);
t7 = t91 * t116 + t147 * t79;
t6 = -pkin(4) * t117 - t9;
t5 = qJ(5) * t117 - qJD(5) * t82 + t8;
t2 = pkin(8) * t140;
t1 = [0.2e1 * m(4) * (-t76 ^ 2 * t80 * t118 + t36 * t18 + t10) + 0.2e1 * (m(5) + m(6)) * (t19 * t3 + t20 * t4 + t10); t121 * t4 + (t16 + t15) * t35 + t120 * t3 + t123 * t20 + t124 * t19 + (t39 + t38) * t17 + (-t83 * t44 + (-t83 * mrSges(3,2) + (-mrSges(4,1) * t82 + mrSges(4,2) * t79 - mrSges(3,1)) * t80) * qJD(2)) * t76 + m(5) * (-t19 * t9 + t20 * t8 - t27 * t3 + t28 * t4) + m(6) * (t17 * t33 + t19 * t6 + t20 * t5 + t21 * t4 + t22 * t3 + t35 * t7) - m(4) * pkin(2) * t109 + (t90 * t144 + m(4) * (-t36 * t117 + t132 + t90) / 0.2e1) * t145 + (t133 + t132 + (t35 * t82 - t36 * t79) * qJD(3)) * mrSges(4,3); -0.2e1 * pkin(2) * t44 + 0.2e1 * t33 * t15 + 0.2e1 * t21 * t26 + 0.2e1 * t22 * t24 + 0.2e1 * t27 * t23 + 0.2e1 * t28 * t25 + 0.2e1 * t7 * t38 + 0.2e1 * t8 * t49 + 0.2e1 * t5 * t52 + 0.2e1 * t9 * t50 + 0.2e1 * t6 * t51 + (t27 * t9 + t28 * t8) * t146 + 0.2e1 * m(6) * (t21 * t5 + t22 * t6 + t33 * t7) + ((0.2e1 * Ifges(4,4) * t82 + t101 * t78 + t122 * t81 + t39 * t145) * qJD(3) + t92) * t82 + (t16 * t145 + (t13 + t14) * t81 + (t11 - t12) * t78 + (t101 * t81 + (t125 * t82 - t122) * t78) * qJD(4) + ((-0.2e1 * Ifges(4,4) + t125 * t81 + (-Ifges(5,6) + Ifges(6,6)) * t78) * t79 + (pkin(7) ^ 2 * t146 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) + t151) * t82) * qJD(3)) * t79; -t18 * mrSges(4,2) + (t42 + t43) * t35 + (t58 + t119) * t17 + m(5) * (-pkin(3) * t17 + t2) + m(6) * (t17 * t54 + t34 * t35 + t2) + 0.2e1 * (t144 + m(6) / 0.2e1) * (t19 * t113 - t20 * t115 + t141) * pkin(8) + (mrSges(5,3) + mrSges(6,2)) * (t141 + t140 + (t19 * t81 - t20 * t78) * qJD(4)); t54 * t15 + t7 * t58 + t34 * t38 + t33 * t42 - pkin(3) * t16 + m(6) * (t33 * t34 + t54 * t7) + (-t74 / 0.2e1 - t75 / 0.2e1 - t72 / 0.2e1 + (Ifges(4,5) + (-m(5) * pkin(3) + t119) * pkin(7)) * qJD(3)) * t82 + (-t9 * mrSges(5,3) + t6 * mrSges(6,2) + t13 / 0.2e1 + t14 / 0.2e1 + t104 * t116 + (t29 / 0.2e1 - t30 / 0.2e1 + t134 / 0.2e1 - t21 * mrSges(6,2) - t28 * mrSges(5,3)) * qJD(4) + (-t121 * qJD(4) + m(6) * (-t21 * qJD(4) + t6) + m(5) * (-t28 * qJD(4) - t9) + t124) * pkin(8)) * t78 + (t8 * mrSges(5,3) + t5 * mrSges(6,2) - t11 / 0.2e1 + t12 / 0.2e1 + t103 * t116 + (t31 / 0.2e1 + t32 / 0.2e1 + t22 * mrSges(6,2) - t27 * mrSges(5,3)) * qJD(4) + (t120 * qJD(4) + m(6) * (t22 * qJD(4) + t5) + m(5) * (-t27 * qJD(4) + t8) + t123) * pkin(8)) * t81 + (-Ifges(4,6) * qJD(3) + (qJD(3) * mrSges(4,2) + t43) * pkin(7) + t104 * t113 + (t45 / 0.2e1 - t46 / 0.2e1 - t103 * qJD(4)) * t78 + t149 * t81 / 0.2e1 + (-Ifges(6,6) * t81 - t150) * qJD(3) / 0.2e1) * t79; -0.2e1 * pkin(3) * t43 + 0.2e1 * t42 * t54 + (-t45 + t46) * t81 + t149 * t78 + 0.2e1 * (m(6) * t54 + t58) * t34 + ((t62 + t63) * t81 + (t60 - t61) * t78) * qJD(4); m(6) * qJD(5) * t20 + (-mrSges(5,2) + t148) * t4 + (-m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1)) * t3; -Ifges(5,6) * t106 - pkin(4) * t24 + m(6) * (-pkin(4) * t6 + qJ(5) * t5 + qJD(5) * t21) + qJD(5) * t52 + qJ(5) * t26 + t5 * mrSges(6,3) - t8 * mrSges(5,2) + t9 * mrSges(5,1) - t6 * mrSges(6,1) + t150 * t114 - t92; -Ifges(5,6) * t115 + t72 + t74 + t75 - t147 * mrSges(6,2) + (m(6) * t111 + (-m(6) * t94 + t58 + t59) * qJD(4)) * pkin(8); 0.2e1 * t148 * qJD(5); m(6) * t3; m(6) * t6 + t24; (m(6) * pkin(8) + mrSges(6,2)) * t113; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
