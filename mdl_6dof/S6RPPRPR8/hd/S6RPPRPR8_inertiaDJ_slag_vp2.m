% Calculate time derivative of joint inertia matrix for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
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
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR8_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR8_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:55:20
% EndTime: 2019-03-09 01:55:23
% DurationCPUTime: 1.12s
% Computational Cost: add. (1669->212), mult. (3348->307), div. (0->0), fcn. (3020->6), ass. (0->98)
t125 = m(6) + m(5);
t110 = sin(qJ(4));
t58 = sin(pkin(9));
t59 = cos(pkin(9));
t63 = cos(qJ(4));
t39 = -t110 * t58 + t59 * t63;
t87 = qJD(4) * t110;
t94 = qJD(4) * t63;
t36 = -t58 * t87 + t59 * t94;
t61 = sin(qJ(6));
t103 = t61 * t36;
t38 = t110 * t59 + t58 * t63;
t62 = cos(qJ(6));
t92 = qJD(6) * t62;
t69 = t38 * t92 + t103;
t99 = (mrSges(6,1) + mrSges(5,3));
t84 = 2 * t99;
t119 = 2 * qJD(2);
t86 = (t58 ^ 2 + t59 ^ 2) * qJD(3);
t60 = -pkin(1) - qJ(3);
t111 = -pkin(7) + t60;
t43 = t111 * t58;
t44 = t111 * t59;
t26 = t110 * t44 + t43 * t63;
t18 = qJD(3) * t39 + qJD(4) * t26;
t35 = t38 * qJD(4);
t10 = -pkin(5) * t35 + t18;
t114 = pkin(4) + pkin(8);
t50 = pkin(3) * t58 + qJ(2);
t80 = -qJ(5) * t39 + t50;
t14 = t114 * t38 + t80;
t25 = t110 * t43 - t44 * t63;
t19 = pkin(5) * t39 + t25;
t3 = -t14 * t61 + t19 * t62;
t67 = qJ(5) * t35 - qJD(5) * t39 + qJD(2);
t8 = t114 * t36 + t67;
t1 = qJD(6) * t3 + t10 * t61 + t62 * t8;
t4 = t14 * t62 + t19 * t61;
t2 = -qJD(6) * t4 + t10 * t62 - t61 * t8;
t123 = t1 * t61 + t2 * t62;
t17 = qJD(3) * t38 + t43 * t87 - t44 * t94;
t101 = t62 * t36;
t93 = qJD(6) * t61;
t90 = t38 * t93;
t68 = t90 - t101;
t11 = mrSges(7,2) * t35 - mrSges(7,3) * t68;
t12 = -mrSges(7,1) * t35 - mrSges(7,3) * t69;
t109 = mrSges(7,3) * t38;
t23 = mrSges(7,1) * t39 - t109 * t61;
t24 = -mrSges(7,2) * t39 + t109 * t62;
t73 = t61 * t23 - t62 * t24;
t122 = -t73 * qJD(6) + t61 * t11 + t62 * t12;
t81 = t3 * t61 - t4 * t62;
t121 = qJD(6) * t81 - t123;
t13 = pkin(4) * t36 + t67;
t120 = -0.2e1 * t13;
t118 = -t61 / 0.2e1;
t117 = t61 / 0.2e1;
t116 = -t62 / 0.2e1;
t115 = t62 / 0.2e1;
t108 = Ifges(7,4) * t61;
t107 = Ifges(7,4) * t62;
t106 = t35 * t39;
t27 = t38 * t36;
t105 = t39 * Ifges(7,6);
t47 = Ifges(7,1) * t62 - t108;
t102 = t61 * t47;
t46 = -Ifges(7,2) * t61 + t107;
t100 = t62 * t46;
t98 = mrSges(6,2) - mrSges(5,1);
t97 = -mrSges(6,3) + mrSges(5,2);
t95 = t61 ^ 2 + t62 ^ 2;
t91 = 0.2e1 * t36;
t85 = Ifges(7,5) * t69 + Ifges(7,6) * t101 - Ifges(7,3) * t35;
t82 = t3 * t62 + t4 * t61;
t45 = mrSges(7,1) * t61 + mrSges(7,2) * t62;
t79 = Ifges(7,1) * t61 + t107;
t78 = Ifges(7,2) * t62 + t108;
t77 = -Ifges(7,5) * t61 - Ifges(7,6) * t62;
t75 = -t17 * t26 + t18 * t25;
t74 = t23 * t62 + t24 * t61;
t72 = -t25 * t35 - t26 * t36;
t71 = qJ(5) * t36 + qJD(5) * t38;
t42 = t79 * qJD(6);
t41 = t78 * qJD(6);
t40 = -mrSges(7,1) * t92 + mrSges(7,2) * t93;
t32 = t35 * mrSges(5,2);
t31 = t35 * mrSges(6,3);
t22 = (-mrSges(7,1) * t62 + mrSges(7,2) * t61) * t38;
t21 = pkin(4) * t38 + t80;
t20 = -pkin(5) * t38 + t26;
t16 = t39 * Ifges(7,5) + t38 * t79;
t15 = t38 * t78 + t105;
t9 = -pkin(5) * t36 - t17;
t7 = mrSges(7,1) * t68 + mrSges(7,2) * t69;
t6 = Ifges(7,1) * t69 - Ifges(7,4) * t68 - t35 * Ifges(7,5);
t5 = Ifges(7,4) * t69 - Ifges(7,2) * t68 - t35 * Ifges(7,6);
t28 = [t15 * t101 + t16 * t103 + 0.2e1 * t50 * (t36 * mrSges(5,1) - t32) + 0.2e1 * t21 * (-t36 * mrSges(6,2) + t31) + 0.2e1 * t20 * t7 + 0.2e1 * t9 * t22 + 0.2e1 * t2 * t23 + 0.2e1 * t1 * t24 + 0.2e1 * t4 * t11 + 0.2e1 * t3 * t12 + 0.2e1 * m(4) * (qJ(2) * qJD(2) - t60 * t86) + 0.2e1 * m(5) * (qJD(2) * t50 + t75) + 0.2e1 * m(6) * (t13 * t21 + t75) + 0.2e1 * m(7) * (t1 * t4 + t2 * t3 + t20 * t9) + (mrSges(5,2) * t119 + mrSges(6,3) * t120 + (-Ifges(5,4) - Ifges(6,6)) * t91 + t18 * t84 + (-(2 * Ifges(5,1)) - (2 * Ifges(6,2)) - Ifges(7,3)) * t35 + t85) * t39 + (mrSges(5,1) * t119 + mrSges(6,2) * t120 + t62 * t5 + t61 * t6 + (Ifges(5,2) + Ifges(6,3)) * t91 + t17 * t84 + (0.2e1 * Ifges(5,4) + 0.2e1 * Ifges(6,6) + t77) * t35 + (t62 * t16 + (-t15 - t105) * t61) * qJD(6)) * t38 + t72 * t84 + (m(3) * qJ(2) + mrSges(4,1) * t58 + mrSges(4,2) * t59 + mrSges(3,3)) * t119 + 0.2e1 * mrSges(4,3) * t86; t38 * t7 + t74 * t35 + (-0.2e1 * t38 * t99 + t22) * t36 + (t35 * t84 - t122) * t39 + m(7) * (t121 * t39 + t20 * t36 + t82 * t35 + t38 * t9) - m(4) * t86 + t125 * (-t17 * t38 - t18 * t39 - t72); 0.2e1 * m(7) * (-t106 * t95 + t27) + 0.2e1 * t125 * (t27 - t106); t62 * t11 - t61 * t12 + t31 - t32 - t98 * t36 - t74 * qJD(6) + (m(5) + m(4)) * qJD(2) + m(7) * (-qJD(6) * t82 + t1 * t62 - t2 * t61) + m(6) * t13; 0; 0; qJ(5) * t7 + qJD(5) * t22 - t20 * t40 + t9 * t45 + t98 * t18 + t97 * t17 + (t6 / 0.2e1 - t2 * mrSges(7,3) - t114 * t12) * t62 + (-t5 / 0.2e1 - t114 * t11 - t1 * mrSges(7,3)) * t61 + m(7) * (qJ(5) * t9 + qJD(5) * t20 - t114 * t123) + m(6) * (-pkin(4) * t18 - qJ(5) * t17 + qJD(5) * t26) + (-qJD(5) * mrSges(6,1) - t41 * t115 - t42 * t117) * t38 + (mrSges(6,1) * pkin(4) + Ifges(7,5) * t116 + Ifges(7,6) * t117 + Ifges(6,4) - Ifges(5,5)) * t35 + (t100 / 0.2e1 + t102 / 0.2e1 - qJ(5) * mrSges(6,1) + Ifges(6,5) - Ifges(5,6)) * t36 + (t39 * t77 / 0.2e1 + t15 * t116 + t16 * t118 + (t115 * t47 + t118 * t46) * t38 + t81 * mrSges(7,3) - (-m(7) * t81 - t73) * t114) * qJD(6); -t38 * t40 + (t45 - t97) * t36 + (-mrSges(7,3) * t95 + t98) * t35 + m(6) * (-pkin(4) * t35 + t71) + m(7) * (-t114 * t35 * t95 + t71); 0; -0.2e1 * qJ(5) * t40 + t41 * t61 - t42 * t62 + (-t100 - t102) * qJD(6) + 0.2e1 * (mrSges(6,3) + t45 + (m(6) + m(7)) * qJ(5)) * qJD(5); m(6) * t18 - m(7) * t121 - t35 * mrSges(6,1) + t122; 0.2e1 * (m(6) / 0.2e1 + m(7) * t95 / 0.2e1) * t35; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 - Ifges(7,6) * t90 + t85; (-t35 * t61 + t39 * t92) * mrSges(7,2) + (t35 * t62 + t39 * t93) * mrSges(7,1); t40; ((mrSges(7,2) * t114 - Ifges(7,6)) * t62 + (mrSges(7,1) * t114 - Ifges(7,5)) * t61) * qJD(6); -t45 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t28(1) t28(2) t28(4) t28(7) t28(11) t28(16); t28(2) t28(3) t28(5) t28(8) t28(12) t28(17); t28(4) t28(5) t28(6) t28(9) t28(13) t28(18); t28(7) t28(8) t28(9) t28(10) t28(14) t28(19); t28(11) t28(12) t28(13) t28(14) t28(15) t28(20); t28(16) t28(17) t28(18) t28(19) t28(20) t28(21);];
Mq  = res;
