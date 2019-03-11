% Calculate time derivative of joint inertia matrix for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
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
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:13:54
% EndTime: 2019-03-09 08:13:57
% DurationCPUTime: 1.27s
% Computational Cost: add. (1800->283), mult. (3652->409), div. (0->0), fcn. (2817->6), ass. (0->116)
t97 = sin(pkin(9));
t98 = cos(pkin(9));
t149 = -mrSges(6,1) * t98 + mrSges(6,2) * t97 - mrSges(5,1);
t100 = sin(qJ(6));
t138 = cos(qJ(6));
t106 = -t100 * t97 + t138 * t98;
t112 = qJD(6) * t138;
t117 = qJD(6) * t100;
t60 = -t98 * t112 + t97 * t117;
t61 = -t97 * t112 - t98 * t117;
t67 = t100 * t98 + t138 * t97;
t148 = -t106 * t61 + t60 * t67;
t31 = mrSges(7,1) * t61 + t60 * mrSges(7,2);
t102 = cos(qJ(2));
t51 = t106 * t102;
t96 = t98 ^ 2;
t111 = (t97 ^ 2 + t96) * qJD(5);
t147 = 2 * m(6);
t146 = 2 * m(7);
t145 = -2 * pkin(1);
t101 = sin(qJ(2));
t73 = -t102 * pkin(2) - t101 * qJ(3) - pkin(1);
t144 = 0.2e1 * t73;
t133 = pkin(7) * t102;
t91 = t102 * qJ(4);
t76 = -t91 + t133;
t143 = 0.2e1 * t76;
t142 = -t106 / 0.2e1;
t141 = -t67 / 0.2e1;
t140 = -t97 / 0.2e1;
t103 = -pkin(2) - pkin(3);
t94 = -qJ(5) + t103;
t139 = pkin(8) - t94;
t136 = mrSges(6,2) * t98;
t135 = Ifges(6,4) * t97;
t134 = Ifges(6,4) * t98;
t119 = qJD(2) * t101;
t110 = qJ(4) * t119 - qJD(4) * t102;
t56 = pkin(7) * t119 - t110;
t132 = t56 * t76;
t130 = t97 * Ifges(6,2);
t129 = mrSges(5,3) - mrSges(4,2);
t128 = pkin(7) - qJ(4);
t99 = qJ(3) + pkin(4);
t118 = qJD(2) * t102;
t126 = qJ(3) * t118 + t101 * qJD(3);
t30 = qJD(5) * t102 + (pkin(4) * t102 + t94 * t101) * qJD(2) + t126;
t58 = -qJD(4) * t101 + t128 * t118;
t11 = t97 * t30 + t98 * t58;
t127 = Ifges(7,5) * t60 - Ifges(7,6) * t61;
t65 = t102 * pkin(3) - t73;
t47 = pkin(4) * t101 + qJ(5) * t102 + t65;
t75 = t128 * t101;
t29 = t97 * t47 + t98 * t75;
t113 = t97 * t119;
t52 = mrSges(6,1) * t113 + t119 * t136;
t125 = mrSges(5,1) * t118 + mrSges(5,2) * t119;
t123 = t102 * t97;
t122 = t102 * t98;
t121 = t98 * t101;
t120 = qJD(3) * t76;
t116 = 2 * mrSges(7,3);
t50 = t67 * t102;
t25 = qJD(6) * t50 + t106 * t119;
t26 = qJD(6) * t51 - t67 * t119;
t115 = Ifges(7,5) * t25 + Ifges(7,6) * t26 + Ifges(7,3) * t118;
t114 = pkin(5) * t97 - pkin(7);
t7 = -t26 * mrSges(7,1) + t25 * mrSges(7,2);
t10 = t98 * t30 - t58 * t97;
t28 = t98 * t47 - t75 * t97;
t109 = -Ifges(6,5) * t98 + Ifges(6,6) * t97;
t108 = -t10 * t97 + t11 * t98;
t107 = t148 * t146;
t12 = pkin(5) * t101 + pkin(8) * t122 + t28;
t17 = pkin(8) * t123 + t29;
t3 = -t100 * t17 + t138 * t12;
t69 = t139 * t97;
t70 = t139 * t98;
t34 = t100 * t70 + t138 * t69;
t4 = t100 * t12 + t138 * t17;
t35 = t100 * t69 - t138 * t70;
t8 = (pkin(5) * t102 - pkin(8) * t121) * qJD(2) + t10;
t9 = -pkin(8) * t113 + t11;
t1 = t3 * qJD(6) + t100 * t8 + t138 * t9;
t2 = -t4 * qJD(6) - t100 * t9 + t138 * t8;
t105 = -t1 * t106 + t2 * t67 - t3 * t60 - t4 * t61;
t13 = -qJD(5) * t106 + qJD(6) * t34;
t14 = qJD(5) * t67 - qJD(6) * t35;
t104 = -t106 * t13 + t14 * t67 - t34 * t60 - t35 * t61;
t82 = pkin(5) * t98 + t99;
t72 = mrSges(6,1) * t101 + mrSges(6,3) * t122;
t71 = -mrSges(6,2) * t101 + mrSges(6,3) * t123;
t64 = (mrSges(6,1) * t102 - mrSges(6,3) * t121) * qJD(2);
t63 = (-mrSges(6,3) * t101 * t97 - mrSges(6,2) * t102) * qJD(2);
t62 = (-mrSges(6,1) * t97 - t136) * t102;
t59 = -t114 * t102 - t91;
t57 = pkin(2) * t119 - t126;
t48 = t103 * t119 + t126;
t46 = (t102 * Ifges(6,5) + (t98 * Ifges(6,1) - t135) * t101) * qJD(2);
t45 = (t102 * Ifges(6,6) + (-t130 + t134) * t101) * qJD(2);
t44 = mrSges(7,1) * t101 + mrSges(7,3) * t51;
t43 = -mrSges(7,2) * t101 + mrSges(7,3) * t50;
t42 = t114 * t119 + t110;
t38 = -Ifges(7,1) * t67 - Ifges(7,4) * t106;
t37 = -Ifges(7,4) * t67 - Ifges(7,2) * t106;
t36 = mrSges(7,1) * t106 - mrSges(7,2) * t67;
t33 = Ifges(7,1) * t60 - Ifges(7,4) * t61;
t32 = Ifges(7,4) * t60 - Ifges(7,2) * t61;
t27 = -mrSges(7,1) * t50 - mrSges(7,2) * t51;
t19 = -Ifges(7,1) * t51 + Ifges(7,4) * t50 + Ifges(7,5) * t101;
t18 = -Ifges(7,4) * t51 + Ifges(7,2) * t50 + Ifges(7,6) * t101;
t16 = -mrSges(7,2) * t118 + mrSges(7,3) * t26;
t15 = mrSges(7,1) * t118 - mrSges(7,3) * t25;
t6 = Ifges(7,1) * t25 + Ifges(7,4) * t26 + Ifges(7,5) * t118;
t5 = Ifges(7,4) * t25 + Ifges(7,2) * t26 + Ifges(7,6) * t118;
t20 = [0.2e1 * t57 * (-mrSges(4,1) * t102 - mrSges(4,3) * t101) + 0.2e1 * t48 * (mrSges(5,1) * t101 - mrSges(5,2) * t102) + 0.2e1 * t11 * t71 + 0.2e1 * t10 * t72 - 0.2e1 * t56 * t62 + 0.2e1 * t29 * t63 + 0.2e1 * t28 * t64 + 0.2e1 * t59 * t7 + 0.2e1 * t2 * t44 + t50 * t5 - t51 * t6 + 0.2e1 * t42 * t27 + 0.2e1 * t1 * t43 + t25 * t19 + t26 * t18 + 0.2e1 * t3 * t15 + 0.2e1 * t4 * t16 + 0.2e1 * t65 * t125 - t46 * t122 + t101 * t115 + (((mrSges(3,2) * t145) - 0.2e1 * t73 * mrSges(4,3) - 0.2e1 * t75 * mrSges(5,3) - Ifges(7,5) * t51 + Ifges(7,6) * t50 + ((2 * Ifges(3,4)) + (2 * Ifges(5,4)) - (2 * Ifges(4,5)) + t109) * t102) * t102 + (mrSges(4,1) * t144 + (mrSges(3,1) * t145) + mrSges(5,3) * t143 + 0.2e1 * (-Ifges(3,4) - Ifges(5,4) + Ifges(4,5) - t109) * t101 + (-t96 * Ifges(6,1) + (2 * Ifges(3,1)) + (2 * Ifges(4,1)) - (2 * Ifges(5,1)) - (2 * Ifges(3,2)) + (2 * Ifges(5,2)) - (2 * Ifges(4,3)) + (2 * Ifges(6,3)) + Ifges(7,3) + (-t130 + 0.2e1 * t134) * t97) * t102) * t101) * qJD(2) + m(4) * t57 * t144 + 0.2e1 * (-t101 * t58 + t102 * t56) * mrSges(5,3) + (t1 * t4 + t2 * t3 + t42 * t59) * t146 + (t10 * t28 + t11 * t29 - t132) * t147 + t45 * t123 + t52 * t143 + 0.2e1 * m(5) * (t48 * t65 + t58 * t75 - t132); t99 * t52 + m(7) * (t1 * t35 + t13 * t4 + t14 * t3 + t2 * t34 + t42 * t82) + t82 * t7 + t58 * mrSges(5,2) + t59 * t31 + t60 * t19 / 0.2e1 - t61 * t18 / 0.2e1 + t14 * t44 + t50 * t32 / 0.2e1 - t51 * t33 / 0.2e1 + t34 * t15 + t35 * t16 + t26 * t37 / 0.2e1 + t25 * t38 / 0.2e1 + t42 * t36 + t13 * t43 + t101 * t127 / 0.2e1 + m(6) * (t120 - t56 * t99 + t108 * t94 + (t28 * t97 - t29 * t98) * qJD(5)) + m(5) * (-qJ(3) * t56 + t103 * t58 + t120) + t105 * mrSges(7,3) + (-t45 / 0.2e1 + t94 * t63 - qJD(5) * t71 - t11 * mrSges(6,3)) * t98 + (-t46 / 0.2e1 - t94 * t64 + qJD(5) * t72 + t10 * mrSges(6,3)) * t97 + (m(4) * t133 + m(7) * t59 - t129 * t102 + t27 + t62) * qJD(3) + t149 * t56 + t6 * t141 + t5 * t142 + ((Ifges(6,5) * t140 - Ifges(6,6) * t98 / 0.2e1 + Ifges(7,5) * t141 + Ifges(7,6) * t142 + Ifges(5,6) + Ifges(4,4) + Ifges(3,5) - pkin(2) * mrSges(4,2) - t103 * mrSges(5,3) + (-m(4) * pkin(2) - mrSges(3,1) - mrSges(4,1)) * pkin(7)) * t102 + (-Ifges(5,5) + Ifges(4,6) - Ifges(3,6) + (-Ifges(6,2) * t98 - t135) * t140 + t98 * (-Ifges(6,1) * t97 - t134) / 0.2e1 + t129 * qJ(3) + (-m(4) * qJ(3) + mrSges(3,2) - mrSges(4,3)) * pkin(7)) * t101) * qJD(2); 0.2e1 * t82 * t31 - t106 * t32 - t67 * t33 - t61 * t37 + t60 * t38 + (t13 * t35 + t14 * t34) * t146 + t104 * t116 + (-t94 * t147 + 0.2e1 * mrSges(6,3)) * t111 + (0.2e1 * mrSges(4,3) + 0.2e1 * t36 - 0.2e1 * t149 + t82 * t146 + t99 * t147 + 0.2e1 * (m(4) + m(5)) * qJ(3)) * qJD(3); -t67 * t15 + t106 * t16 + t61 * t43 + t60 * t44 + t98 * t63 - t97 * t64 + (m(4) * pkin(7) - t129) * t118 - m(7) * t105 + m(6) * t108 + m(5) * t58; -m(6) * t111 - m(7) * t104 + t148 * t116; -t107; t106 * t15 + t67 * t16 - t60 * t43 + t61 * t44 + t97 * t63 + t98 * t64 + m(7) * (t1 * t67 + t106 * t2 + t3 * t61 - t4 * t60) + m(6) * (t10 * t98 + t11 * t97) + m(5) * t48 + t125; m(7) * (t106 * t14 + t13 * t67 + t34 * t61 - t35 * t60); 0; -t107; -m(6) * t56 + m(7) * t42 + t52 + t7; (m(6) + m(7)) * qJD(3) + t31; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t115; mrSges(7,1) * t14 - mrSges(7,2) * t13 + t127; mrSges(7,1) * t60 - mrSges(7,2) * t61; t31; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;
