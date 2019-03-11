% Calculate time derivative of joint inertia matrix for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:34:32
% EndTime: 2019-03-08 19:34:35
% DurationCPUTime: 1.20s
% Computational Cost: add. (1169->235), mult. (3110->353), div. (0->0), fcn. (2720->10), ass. (0->115)
t66 = cos(qJ(4));
t101 = qJD(6) * t66;
t65 = cos(qJ(6));
t105 = qJD(4) * t65;
t62 = sin(qJ(6));
t63 = sin(qJ(4));
t74 = t62 * t101 + t63 * t105;
t45 = mrSges(7,1) * t62 + mrSges(7,2) * t65;
t134 = mrSges(6,3) + t45;
t44 = t66 * mrSges(6,2) - t63 * mrSges(6,3);
t133 = -mrSges(5,1) * t66 + mrSges(5,2) * t63 + t44;
t106 = qJD(4) * t63;
t96 = t62 * t106;
t97 = t65 * t101;
t132 = -t96 + t97;
t109 = cos(pkin(6));
t59 = sin(pkin(11));
t60 = sin(pkin(6));
t61 = cos(pkin(11));
t64 = sin(qJ(2));
t67 = cos(qJ(2));
t76 = (t59 * t67 + t61 * t64) * t60;
t75 = t63 * t76;
t18 = -t109 * t66 + t75;
t79 = t59 * t64 - t61 * t67;
t26 = t79 * t60;
t10 = t18 * t62 + t26 * t65;
t24 = qJD(2) * t76;
t19 = t109 * t63 + t66 * t76;
t108 = qJD(2) * t60;
t25 = t79 * t108;
t7 = qJD(4) * t19 - t25 * t63;
t1 = -qJD(6) * t10 - t24 * t62 + t65 * t7;
t9 = t18 * t65 - t26 * t62;
t2 = qJD(6) * t9 + t24 * t65 + t62 * t7;
t131 = t1 * t65 + t2 * t62 + (t10 * t65 - t62 * t9) * qJD(6);
t130 = 2 * m(7);
t129 = -t62 / 0.2e1;
t128 = -t65 / 0.2e1;
t127 = t65 / 0.2e1;
t68 = -pkin(4) - pkin(9);
t126 = pkin(4) * t66;
t8 = -qJD(4) * t75 + (qJD(4) * t109 - t25) * t66;
t5 = t19 * t8;
t6 = t63 * t8;
t125 = t66 * t7;
t53 = pkin(2) * t59 + pkin(8);
t124 = pkin(5) + t53;
t104 = qJD(4) * t66;
t123 = t19 * t104 + t6;
t122 = mrSges(7,3) * t66;
t121 = Ifges(7,4) * t62;
t120 = Ifges(7,4) * t65;
t119 = Ifges(7,5) * t63;
t11 = t26 * t24;
t83 = Ifges(7,1) * t62 + t120;
t29 = -t66 * t83 + t119;
t118 = t62 * t29;
t47 = Ifges(7,1) * t65 - t121;
t117 = t62 * t47;
t116 = t62 * t68;
t82 = Ifges(7,2) * t65 + t121;
t28 = t63 * Ifges(7,6) - t66 * t82;
t115 = t65 * t28;
t46 = -Ifges(7,2) * t62 + t120;
t114 = t65 * t46;
t113 = t65 * t68;
t112 = -mrSges(5,1) + mrSges(6,2);
t111 = t62 ^ 2 + t65 ^ 2;
t110 = qJ(5) * t63;
t36 = t124 * t66;
t31 = qJD(4) * t36;
t107 = qJD(4) * t62;
t103 = qJD(6) * t62;
t102 = qJD(6) * t65;
t100 = t63 * qJD(5);
t99 = -mrSges(5,2) + t134;
t54 = -pkin(2) * t61 - pkin(3);
t94 = m(7) * t111;
t93 = m(6) * t53 + mrSges(6,1);
t91 = Ifges(7,5) * t96 + t74 * Ifges(7,6) + Ifges(7,3) * t104;
t90 = pkin(4) * t106 - t100;
t77 = t54 - t110;
t27 = t66 * t68 + t77;
t35 = t124 * t63;
t12 = -t27 * t62 + t35 * t65;
t23 = (pkin(9) * t63 - qJ(5) * t66) * qJD(4) + t90;
t3 = qJD(6) * t12 + t23 * t65 + t31 * t62;
t13 = t27 * t65 + t35 * t62;
t4 = -qJD(6) * t13 - t23 * t62 + t31 * t65;
t88 = t3 * t62 + t4 * t65;
t87 = t7 * t63 + t8 * t66;
t84 = mrSges(7,1) * t65 - mrSges(7,2) * t62;
t81 = -Ifges(7,5) * t62 - Ifges(7,6) * t65;
t80 = t12 * t62 - t13 * t65;
t78 = t8 * qJ(5) + t19 * qJD(5);
t21 = -mrSges(7,2) * t104 + mrSges(7,3) * t74;
t22 = mrSges(7,1) * t104 + t132 * mrSges(7,3);
t42 = mrSges(7,1) * t63 + t122 * t62;
t43 = -mrSges(7,2) * t63 - t122 * t65;
t70 = -t102 * t43 + t103 * t42 - t62 * t21 - t65 * t22;
t69 = m(7) * t131;
t41 = t83 * qJD(6);
t40 = t82 * qJD(6);
t39 = (t63 * mrSges(5,1) + t66 * mrSges(5,2)) * qJD(4);
t38 = (-t63 * mrSges(6,2) - t66 * mrSges(6,3)) * qJD(4);
t37 = t84 * qJD(6);
t34 = t84 * t66;
t33 = t77 - t126;
t32 = qJ(5) * t104 - t90;
t30 = t124 * t106;
t16 = t74 * mrSges(7,1) + t132 * mrSges(7,2);
t15 = -t47 * t101 + (Ifges(7,5) * t66 + t63 * t83) * qJD(4);
t14 = -t46 * t101 + (Ifges(7,6) * t66 + t63 * t82) * qJD(4);
t17 = [0.2e1 * m(7) * (t1 * t9 + t10 * t2 + t5) + 0.2e1 * m(4) * (-t25 * t76 + t11) + 0.2e1 * (m(6) + m(5)) * (t18 * t7 + t11 + t5); t25 * mrSges(4,2) + t1 * t42 + t10 * t21 - t19 * t16 + t2 * t43 + t9 * t22 + t8 * t34 + (t38 + t39) * t26 + (-mrSges(3,1) * t64 - mrSges(3,2) * t67) * t108 + (-mrSges(4,1) + t133) * t24 + m(6) * (t24 * t33 - t26 * t32) + m(5) * t24 * t54 + m(7) * (t1 * t12 + t10 * t3 + t13 * t2 - t19 * t30 + t36 * t8 + t4 * t9) + m(4) * (-t24 * t61 - t25 * t59) * pkin(2) + 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * (t104 * t18 - t106 * t19 + t87) * t53 + (mrSges(5,3) + mrSges(6,1)) * ((t18 * t66 - t19 * t63) * qJD(4) + t87); (t12 * t4 + t13 * t3 - t30 * t36) * t130 + 0.2e1 * t54 * t39 + 0.2e1 * t33 * t38 + 0.2e1 * t4 * t42 + 0.2e1 * t3 * t43 - 0.2e1 * t30 * t34 - 0.2e1 * t36 * t16 + 0.2e1 * t13 * t21 + 0.2e1 * t12 * t22 + 0.2e1 * (-m(6) * t33 - t44) * t32 + ((t115 + t118 + 0.2e1 * (-Ifges(5,4) - Ifges(6,6)) * t63) * qJD(4) + t91) * t63 + (-t65 * t14 - t62 * t15 + (t62 * t28 + (-t29 - t119) * t65) * qJD(6) + ((0.2e1 * Ifges(5,4) + 0.2e1 * Ifges(6,6) + t81) * t66 + ((2 * Ifges(5,1)) - (2 * Ifges(5,2)) + (2 * Ifges(6,2)) - (2 * Ifges(6,3)) + Ifges(7,3)) * t63) * qJD(4)) * t66; m(6) * (t6 - t125 + (t18 * t63 + t19 * t66) * qJD(4)) + m(7) * ((t10 * t62 + t65 * t9) * t106 - t131 * t66 + t123) + m(5) * (t106 * t18 + t123 - t125); (-t16 + t42 * t105 + m(7) * (t105 * t12 + t107 * t13 - t30) + t43 * t107) * t63 + (qJD(4) * t34 + m(7) * (-t102 * t13 + t103 * t12 + t31 - t88) + t70) * t66; (0.1e1 - t111) * t63 * t104 * t130; t19 * t37 + t112 * t7 + t99 * t8 + m(6) * (-pkin(4) * t7 + t78) + m(7) * t78 + t68 * t69 - t131 * mrSges(7,3); m(7) * (-qJ(5) * t30 + qJD(5) * t36 + t113 * t4 + t116 * t3) + t22 * t113 + t21 * t116 + t15 * t127 + t14 * t129 - t30 * t45 + qJD(5) * t34 + t36 * t37 - qJ(5) * t16 - t88 * mrSges(7,3) + (t93 * qJD(5) - t40 * t128 - t41 * t129) * t66 + (-t115 / 0.2e1 - t118 / 0.2e1 + t63 * t81 / 0.2e1 + (t47 * t128 + t62 * t46 / 0.2e1) * t66 + t80 * mrSges(7,3) + (-m(7) * t80 - t62 * t42 + t65 * t43) * t68) * qJD(6) + ((-pkin(4) * mrSges(6,1) + Ifges(7,5) * t127 + Ifges(7,6) * t129 - Ifges(6,4) + Ifges(5,5)) * t66 + (Ifges(6,5) - Ifges(5,6) + t114 / 0.2e1 + t117 / 0.2e1 - qJ(5) * mrSges(6,1)) * t63 + (m(6) * (-t110 - t126) + t133) * t53) * qJD(4); t63 * t37 + m(7) * t100 + m(6) * t32 + ((m(7) * qJ(5) + t99) * t66 + (-mrSges(7,3) * t111 + t68 * t94 + t112) * t63) * qJD(4); 0.2e1 * qJ(5) * t37 + t40 * t62 - t41 * t65 + (-t114 - t117) * qJD(6) + 0.2e1 * ((m(6) + m(7)) * qJ(5) + t134) * qJD(5); m(6) * t7 + t69; m(7) * (-qJD(6) * t80 + t88) + t93 * t104 - t70; (m(6) + t94) * t106; 0; 0; mrSges(7,1) * t1 - mrSges(7,2) * t2; mrSges(7,1) * t4 - mrSges(7,2) * t3 - Ifges(7,5) * t97 + t91; t16; ((-mrSges(7,2) * t68 - Ifges(7,6)) * t65 + (-mrSges(7,1) * t68 - Ifges(7,5)) * t62) * qJD(6); -t45 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t17(1) t17(2) t17(4) t17(7) t17(11) t17(16); t17(2) t17(3) t17(5) t17(8) t17(12) t17(17); t17(4) t17(5) t17(6) t17(9) t17(13) t17(18); t17(7) t17(8) t17(9) t17(10) t17(14) t17(19); t17(11) t17(12) t17(13) t17(14) t17(15) t17(20); t17(16) t17(17) t17(18) t17(19) t17(20) t17(21);];
Mq  = res;
