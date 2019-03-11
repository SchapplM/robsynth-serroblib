% Calculate time derivative of joint inertia matrix for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:08:12
% EndTime: 2019-03-08 21:08:15
% DurationCPUTime: 1.25s
% Computational Cost: add. (1031->249), mult. (2512->362), div. (0->0), fcn. (1897->8), ass. (0->112)
t141 = m(5) + m(6);
t73 = cos(qJ(3));
t105 = qJD(6) * t73;
t70 = sin(qJ(3));
t110 = qJD(3) * t70;
t69 = sin(qJ(6));
t72 = cos(qJ(6));
t80 = t69 * t105 + t72 * t110;
t131 = -t72 / 0.2e1;
t116 = pkin(8) - qJ(5);
t74 = cos(qJ(2));
t111 = qJD(2) * t74;
t66 = sin(pkin(6));
t101 = t66 * t111;
t71 = sin(qJ(2));
t123 = t66 * t71;
t104 = t70 * t123;
t67 = cos(pkin(6));
t15 = -qJD(3) * t104 + (qJD(3) * t67 + t101) * t73;
t29 = t123 * t73 + t67 * t70;
t24 = t29 * qJD(4);
t139 = t15 * qJ(4) + t24;
t122 = t66 * t74;
t28 = -t67 * t73 + t104;
t16 = t122 * t72 - t28 * t69;
t17 = t122 * t69 + t28 * t72;
t112 = qJD(2) * t71;
t102 = t66 * t112;
t14 = qJD(3) * t29 + t101 * t70;
t3 = -qJD(6) * t17 - t102 * t72 - t14 * t69;
t4 = qJD(6) * t16 - t102 * t69 + t14 * t72;
t138 = qJD(6) * (t16 * t72 + t17 * t69) + t69 * t3 - t72 * t4;
t109 = qJD(3) * t73;
t114 = qJ(4) * t109 + t70 * qJD(4);
t75 = -pkin(3) - pkin(4);
t65 = -pkin(9) + t75;
t12 = (pkin(5) * t73 + t65 * t70) * qJD(3) + t114;
t27 = -qJD(5) * t70 + t109 * t116;
t40 = -t73 * pkin(3) - t70 * qJ(4) - pkin(2);
t31 = t73 * pkin(4) - t40;
t20 = pkin(5) * t70 + pkin(9) * t73 + t31;
t43 = t116 * t70;
t9 = t20 * t72 - t43 * t69;
t1 = qJD(6) * t9 + t12 * t69 + t27 * t72;
t10 = t20 * t69 + t43 * t72;
t2 = -qJD(6) * t10 + t12 * t72 - t27 * t69;
t137 = -t1 * t72 + t2 * t69 + (t10 * t69 + t72 * t9) * qJD(6);
t136 = -2 * mrSges(6,3);
t130 = pkin(8) * t73;
t45 = -t73 * qJ(5) + t130;
t135 = 0.2e1 * t45;
t68 = qJ(4) + pkin(5);
t134 = 0.2e1 * t68;
t133 = m(5) / 0.2e1;
t132 = t69 / 0.2e1;
t129 = mrSges(7,3) * t73;
t128 = Ifges(7,4) * t69;
t127 = Ifges(7,4) * t72;
t126 = Ifges(7,5) * t72;
t125 = t14 * t70;
t25 = qJD(5) * t73 + t110 * t116;
t124 = t25 * t45;
t6 = t29 * t15;
t85 = Ifges(7,2) * t72 + t128;
t121 = t69 * t85;
t87 = Ifges(7,1) * t69 + t127;
t120 = t72 * t87;
t119 = -mrSges(4,1) - mrSges(5,1);
t118 = mrSges(5,3) - mrSges(4,2);
t117 = mrSges(6,3) - mrSges(5,2);
t42 = mrSges(7,1) * t72 - mrSges(7,2) * t69;
t115 = -t42 - mrSges(6,1);
t33 = mrSges(6,1) * t109 + mrSges(6,2) * t110;
t108 = qJD(4) * t45;
t107 = qJD(6) * t69;
t106 = qJD(6) * t72;
t103 = mrSges(4,3) - t117;
t100 = t69 * t110;
t98 = t72 * t105;
t96 = t80 * Ifges(7,5) + Ifges(7,6) * t98 + Ifges(7,3) * t109;
t94 = m(5) * pkin(8) - t117;
t93 = t103 * t73;
t89 = -mrSges(7,1) * t69 - mrSges(7,2) * t72;
t88 = Ifges(7,1) * t72 - t128;
t86 = -Ifges(7,2) * t69 + t127;
t84 = t15 * t45 - t25 * t29;
t82 = -Ifges(7,6) * t69 - (2 * Ifges(4,4)) - (2 * Ifges(6,4)) + (2 * Ifges(5,5));
t79 = -t98 + t100;
t77 = m(7) * t138;
t18 = mrSges(7,1) * t109 - mrSges(7,3) * t80;
t19 = -mrSges(7,2) * t109 - mrSges(7,3) * t79;
t76 = -m(7) * t137 - t69 * t18 + t72 * t19;
t58 = Ifges(7,6) * t107;
t46 = mrSges(6,1) * t70 - mrSges(6,2) * t73;
t44 = -t73 * mrSges(5,1) - t70 * mrSges(5,3);
t39 = mrSges(7,1) * t70 + t72 * t129;
t38 = -mrSges(7,2) * t70 + t69 * t129;
t37 = t88 * qJD(6);
t36 = t86 * qJD(6);
t35 = (mrSges(4,1) * t70 + mrSges(4,2) * t73) * qJD(3);
t34 = (mrSges(5,1) * t70 - mrSges(5,3) * t73) * qJD(3);
t32 = t89 * qJD(6);
t30 = t89 * t73;
t26 = pkin(3) * t110 - t114;
t23 = Ifges(7,5) * t70 - t73 * t88;
t22 = Ifges(7,6) * t70 - t73 * t86;
t21 = t75 * t110 + t114;
t13 = t15 * t130;
t11 = mrSges(7,1) * t79 + mrSges(7,2) * t80;
t8 = t87 * t105 + (Ifges(7,5) * t73 + t70 * t88) * qJD(3);
t7 = t85 * t105 + (Ifges(7,6) * t73 + t70 * t86) * qJD(3);
t5 = [0.2e1 * m(7) * (t16 * t3 + t17 * t4 + t6) + 0.2e1 * (m(4) + t141) * (-t66 ^ 2 * t71 * t111 + t28 * t14 + t6); t29 * t11 + t16 * t18 + t17 * t19 + t3 * t39 + t4 * t38 + t103 * t125 + (t30 + t93) * t15 + (-t103 * t29 * t70 + t28 * t93) * qJD(3) + ((t33 - t34 - t35) * t74 + (-t74 * mrSges(3,2) + (-mrSges(4,1) * t73 + t70 * mrSges(4,2) - mrSges(3,1) + t44 - t46) * t71) * qJD(2)) * t66 + m(6) * (t14 * t43 + t27 * t28 + (-t112 * t31 + t21 * t74) * t66 + t84) + m(4) * (-pkin(2) * t102 + t13) + m(5) * (t102 * t40 - t122 * t26 + t13) + m(7) * (t1 * t17 + t10 * t4 + t16 * t2 + t3 * t9 + t84) + 0.2e1 * (m(4) / 0.2e1 + t133) * (t109 * t28 - t110 * t29 + t125) * pkin(8); -0.2e1 * pkin(2) * t35 + 0.2e1 * t1 * t38 + 0.2e1 * t10 * t19 + t11 * t135 + 0.2e1 * t9 * t18 + 0.2e1 * t2 * t39 + 0.2e1 * t21 * t46 - 0.2e1 * t25 * t30 + 0.2e1 * t31 * t33 + 0.2e1 * t40 * t34 + 0.2e1 * (m(5) * t40 + t44) * t26 + 0.2e1 * m(7) * (t1 * t10 + t2 * t9 - t124) + 0.2e1 * m(6) * (t21 * t31 + t27 * t43 - t124) + (t27 * t136 + (mrSges(6,3) * t135 - t69 * t22 + t72 * t23 + t70 * t82) * qJD(3) + t96) * t70 + (0.2e1 * t25 * mrSges(6,3) + t69 * t7 - t72 * t8 + (t72 * t22 + t69 * t23) * qJD(6) + (t43 * t136 + (-t82 - t126) * t73 + ((2 * Ifges(4,1)) + (2 * Ifges(5,1)) - (2 * Ifges(6,1)) - (2 * Ifges(4,2)) + (2 * Ifges(6,2)) - (2 * Ifges(5,3)) + Ifges(7,3)) * t70) * qJD(3)) * t73; t29 * t32 + (mrSges(6,2) + t119) * t14 + (-t115 + t118) * t15 + m(6) * (t14 * t75 + t139) + m(7) * (t15 * t68 + t24) + m(5) * (-pkin(3) * t14 + t139) - t65 * t77 + t138 * mrSges(7,3); t70 * (-Ifges(7,5) * t106 + t58) / 0.2e1 + t7 * t131 + t68 * t11 - t23 * t106 / 0.2e1 + qJD(4) * t30 + t45 * t32 + t27 * mrSges(6,2) + (-t8 / 0.2e1 + qJD(6) * t22 / 0.2e1) * t69 + t115 * t25 + m(6) * (-qJ(4) * t25 + t27 * t75 + t108) + m(7) * (-t25 * t68 + t108) + t137 * mrSges(7,3) + (-t106 * t39 - t107 * t38 + t76) * t65 + (-t37 * t131 - t36 * t132 + (t131 * t85 - t132 * t87) * qJD(6) + t94 * qJD(4)) * t73 + ((-pkin(3) * mrSges(5,2) - t75 * mrSges(6,3) - Ifges(7,5) * t69 / 0.2e1 + Ifges(7,6) * t131 + Ifges(6,6) + Ifges(5,4) + Ifges(4,5) + (-m(5) * pkin(3) + t119) * pkin(8)) * t73 + (-t120 / 0.2e1 + t121 / 0.2e1 - Ifges(6,5) + Ifges(5,6) - Ifges(4,6) + t117 * qJ(4) + (-m(5) * qJ(4) - t118) * pkin(8)) * t70) * qJD(3); t32 * t134 + t36 * t72 + t37 * t69 + (t120 - t121) * qJD(6) + (m(7) * t134 + 0.2e1 * t141 * qJ(4) + 0.2e1 * mrSges(6,1) + 0.2e1 * mrSges(5,3) + 0.2e1 * t42) * qJD(4); -t77 + 0.2e1 * (m(6) / 0.2e1 + t133) * t14; (-t69 * t38 - t72 * t39) * qJD(6) + m(6) * t27 + t94 * t109 + t76; 0; 0; -m(6) * t102 + m(7) * (t3 * t72 + t4 * t69 + (-t16 * t69 + t17 * t72) * qJD(6)); t72 * t18 + t69 * t19 + (t72 * t38 - t69 * t39) * qJD(6) + m(7) * (t1 * t69 + t2 * t72 + (t10 * t72 - t69 * t9) * qJD(6)) + m(6) * t21 + t33; 0; 0; 0; mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t2 - mrSges(7,2) * t1 - Ifges(7,6) * t100 + t96; t58 + (-t42 * t65 - t126) * qJD(6); -t42 * qJD(6); t32; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t5(1) t5(2) t5(4) t5(7) t5(11) t5(16); t5(2) t5(3) t5(5) t5(8) t5(12) t5(17); t5(4) t5(5) t5(6) t5(9) t5(13) t5(18); t5(7) t5(8) t5(9) t5(10) t5(14) t5(19); t5(11) t5(12) t5(13) t5(14) t5(15) t5(20); t5(16) t5(17) t5(18) t5(19) t5(20) t5(21);];
Mq  = res;
