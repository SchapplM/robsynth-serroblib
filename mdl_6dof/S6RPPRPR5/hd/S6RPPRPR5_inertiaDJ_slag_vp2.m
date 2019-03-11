% Calculate time derivative of joint inertia matrix for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:48:32
% EndTime: 2019-03-09 01:48:35
% DurationCPUTime: 1.48s
% Computational Cost: add. (1436->271), mult. (3162->423), div. (0->0), fcn. (2554->6), ass. (0->110)
t87 = sin(qJ(4));
t102 = t87 * qJD(2);
t89 = cos(qJ(4));
t104 = qJD(4) * t89;
t84 = qJ(2) - pkin(7);
t132 = t84 * t104 + t102;
t131 = 0.2e1 * t84;
t130 = 2 * qJD(3);
t83 = cos(pkin(9));
t79 = t83 ^ 2;
t82 = sin(pkin(9));
t108 = t82 ^ 2 + t79;
t129 = qJD(5) * m(6) * t108;
t86 = sin(qJ(6));
t88 = cos(qJ(6));
t61 = t82 * t86 - t88 * t83;
t128 = qJD(6) * t61;
t62 = t82 * t88 + t83 * t86;
t127 = t62 * qJD(6);
t126 = 2 * m(7);
t125 = -t61 / 0.2e1;
t124 = t62 / 0.2e1;
t123 = t82 / 0.2e1;
t122 = pkin(8) * t89;
t121 = t87 * pkin(4);
t120 = mrSges(5,2) * t87;
t119 = Ifges(6,4) * t82;
t118 = Ifges(6,4) * t83;
t117 = Ifges(6,2) * t82;
t115 = t83 * t87;
t114 = t84 * t87;
t112 = t89 * mrSges(6,3);
t85 = pkin(1) + qJ(3);
t111 = pkin(8) + qJ(5);
t110 = -Ifges(7,5) * t128 - Ifges(7,6) * t127;
t109 = -mrSges(6,1) * t83 + mrSges(6,2) * t82 - mrSges(5,1);
t107 = qJ(5) * t89;
t65 = -t107 + t85 + t121;
t39 = t83 * t114 + t82 * t65;
t80 = t87 ^ 2;
t106 = qJD(2) * t80;
t105 = qJD(4) * t87;
t81 = t89 ^ 2;
t77 = t81 * qJD(2);
t103 = t85 * qJD(3);
t101 = t89 * qJD(2);
t20 = t61 * t105 - t127 * t89;
t22 = t62 * t105 + t128 * t89;
t100 = Ifges(7,5) * t20 + Ifges(7,6) * t22 + Ifges(7,3) * t104;
t49 = -qJD(5) * t89 + qJD(3) + (pkin(4) * t89 + qJ(5) * t87) * qJD(4);
t25 = t132 * t83 + t82 * t49;
t97 = -t82 * t84 + pkin(5);
t96 = pkin(5) * t82 - t84;
t95 = t108 * mrSges(6,3);
t5 = -t22 * mrSges(7,1) + t20 * mrSges(7,2);
t94 = m(6) * pkin(4) - t109;
t93 = mrSges(6,1) * t82 + mrSges(6,2) * t83;
t27 = mrSges(7,1) * t127 - mrSges(7,2) * t128;
t92 = -Ifges(6,5) * t83 + Ifges(6,6) * t82;
t59 = t83 * t65;
t26 = -t83 * t122 + t97 * t87 + t59;
t30 = -t82 * t122 + t39;
t6 = t26 * t88 - t30 * t86;
t7 = t26 * t86 + t30 * t88;
t67 = t111 * t82;
t69 = t111 * t83;
t34 = -t67 * t88 - t69 * t86;
t35 = -t67 * t86 + t69 * t88;
t46 = t62 * t89;
t48 = t61 * t89;
t75 = -pkin(5) * t83 - pkin(4);
t73 = t84 * t77;
t64 = t87 * mrSges(6,1) - t83 * t112;
t63 = -t87 * mrSges(6,2) - t82 * t112;
t60 = t96 * t89;
t57 = (mrSges(6,1) * t89 + mrSges(6,3) * t115) * qJD(4);
t56 = (mrSges(6,3) * t82 * t87 - mrSges(6,2) * t89) * qJD(4);
t55 = t93 * t89;
t50 = t93 * t105;
t47 = t61 * t87;
t45 = t62 * t87;
t44 = -t96 * t105 - t101;
t43 = t83 * t49;
t41 = (Ifges(6,5) * t89 + (-Ifges(6,1) * t83 + t119) * t87) * qJD(4);
t40 = (Ifges(6,6) * t89 + (t117 - t118) * t87) * qJD(4);
t38 = -t82 * t114 + t59;
t37 = mrSges(7,1) * t87 + mrSges(7,3) * t48;
t36 = -mrSges(7,2) * t87 - mrSges(7,3) * t46;
t33 = Ifges(7,1) * t62 - Ifges(7,4) * t61;
t32 = Ifges(7,4) * t62 - Ifges(7,2) * t61;
t31 = mrSges(7,1) * t61 + mrSges(7,2) * t62;
t29 = -Ifges(7,1) * t128 - Ifges(7,4) * t127;
t28 = -Ifges(7,4) * t128 - Ifges(7,2) * t127;
t24 = -t132 * t82 + t43;
t23 = mrSges(7,1) * t46 - mrSges(7,2) * t48;
t21 = -qJD(4) * t46 + t128 * t87;
t19 = -qJD(4) * t48 - t127 * t87;
t15 = -Ifges(7,1) * t48 - Ifges(7,4) * t46 + Ifges(7,5) * t87;
t14 = -Ifges(7,4) * t48 - Ifges(7,2) * t46 + Ifges(7,6) * t87;
t13 = -t62 * qJD(5) - t35 * qJD(6);
t12 = -t61 * qJD(5) + t34 * qJD(6);
t11 = pkin(8) * t82 * t105 + t25;
t10 = -mrSges(7,2) * t104 + t22 * mrSges(7,3);
t9 = mrSges(7,1) * t104 - t20 * mrSges(7,3);
t8 = -t82 * t102 + t43 + (pkin(8) * t115 + t97 * t89) * qJD(4);
t4 = Ifges(7,1) * t20 + Ifges(7,4) * t22 + Ifges(7,5) * t104;
t3 = Ifges(7,4) * t20 + Ifges(7,2) * t22 + Ifges(7,6) * t104;
t2 = -t7 * qJD(6) - t11 * t86 + t8 * t88;
t1 = t6 * qJD(6) + t11 * t88 + t8 * t86;
t16 = [t87 * t100 + 0.2e1 * t25 * t63 + 0.2e1 * t24 * t64 + 0.2e1 * t39 * t56 + 0.2e1 * t38 * t57 + 0.2e1 * t60 * t5 - t46 * t3 - t48 * t4 + 0.2e1 * t1 * t36 + 0.2e1 * t2 * t37 + 0.2e1 * t44 * t23 + t22 * t14 + t20 * t15 + 0.2e1 * t6 * t9 + 0.2e1 * t7 * t10 + (t87 * mrSges(5,1) + mrSges(4,3)) * t130 + (mrSges(5,2) * t130 + t50 * t131 - t82 * t40 + t83 * t41) * t89 + 0.2e1 * (m(3) * qJ(2) - t89 * t55 + mrSges(4,2) + mrSges(3,3) + (-t80 - t81) * mrSges(5,3)) * qJD(2) + 0.2e1 * m(5) * (t84 * t106 + t103 + t73) + 0.2e1 * m(4) * (qJ(2) * qJD(2) + t103) + (t1 * t7 + t2 * t6 + t44 * t60) * t126 + 0.2e1 * m(6) * (t38 * t24 + t39 * t25 + t73) + ((0.2e1 * t85 * mrSges(5,1) - Ifges(7,5) * t48 - Ifges(7,6) * t46 + (-(2 * Ifges(5,4)) - t92) * t89) * t89 + (t55 * t131 - 0.2e1 * t85 * mrSges(5,2) + 0.2e1 * (Ifges(5,4) + t92) * t87 + (-Ifges(6,1) * t79 + (2 * Ifges(6,3)) - (2 * Ifges(5,1)) + (2 * Ifges(5,2)) - 0.2e1 * m(6) * t84 ^ 2 + Ifges(7,3) + (-t117 + 0.2e1 * t118) * t82) * t89) * t87) * qJD(4); -t62 * t10 + t128 * t36 + t127 * t37 - t82 * t56 - t83 * t57 + t61 * t9 + (-mrSges(5,1) * t89 + t120) * qJD(4) + (-m(5) - m(4)) * qJD(3) + m(7) * (-t1 * t62 + t127 * t6 + t128 * t7 + t2 * t61) + m(6) * (-t24 * t83 - t25 * t82); (t127 * t61 - t128 * t62) * t126; m(4) * qJD(2) - t47 * t10 + t19 * t36 + t21 * t37 - t45 * t9 + (-t5 + t50) * t89 + (t83 * t56 - t82 * t57) * t87 + ((t83 * t63 - t82 * t64) * t89 + (t23 + t55) * t87) * qJD(4) + m(7) * (-t47 * t1 + t60 * t105 + t19 * t7 - t45 * t2 + t21 * t6 - t44 * t89) + m(6) * (t77 + (-t24 * t82 + t25 * t83) * t87 + (-t38 * t82 + t39 * t83 - 0.2e1 * t114) * t104) + m(5) * (t77 + t106); m(7) * (-t127 * t45 - t128 * t47 - t19 * t62 + t21 * t61); (-t47 * t19 - t45 * t21) * t126 + 0.4e1 * (m(6) * (-0.1e1 + t108) / 0.2e1 - m(7) / 0.2e1) * t87 * t104; t87 * t110 / 0.2e1 + t75 * t5 + t4 * t124 - t128 * t15 / 0.2e1 - t127 * t14 / 0.2e1 + t60 * t27 + t3 * t125 - t46 * t28 / 0.2e1 - t48 * t29 / 0.2e1 + pkin(4) * t50 + t22 * t32 / 0.2e1 + t20 * t33 / 0.2e1 + t34 * t9 + t35 * t10 + t12 * t36 + t13 * t37 + t44 * t31 + m(7) * (t1 * t35 + t12 * t7 + t13 * t6 + t2 * t34 + t44 * t75) + (t94 * t89 - t120) * qJD(2) + (t40 / 0.2e1 + qJD(5) * t63 + qJ(5) * t56 + t25 * mrSges(6,3) + m(6) * (qJ(5) * t25 + qJD(5) * t39)) * t83 + (t41 / 0.2e1 - qJD(5) * t64 - qJ(5) * t57 - t24 * mrSges(6,3) + m(6) * (-qJ(5) * t24 - qJD(5) * t38)) * t82 + (-t1 * t61 - t127 * t7 + t128 * t6 - t2 * t62) * mrSges(7,3) + ((-t84 * mrSges(5,2) + Ifges(6,5) * t123 + Ifges(6,6) * t83 / 0.2e1 + Ifges(7,5) * t124 + Ifges(7,6) * t125 - Ifges(5,6)) * t89 + (-t83 * (Ifges(6,1) * t82 + t118) / 0.2e1 + (Ifges(6,2) * t83 + t119) * t123 - Ifges(5,5) - t94 * t84) * t87) * qJD(4); m(7) * (-t12 * t62 + t127 * t34 + t128 * t35 + t13 * t61); -t89 * t27 + m(7) * (-t12 * t47 - t13 * t45 + t19 * t35 + t21 * t34) + t87 * t129 + ((-mrSges(5,2) + t95) * t89 + m(6) * (t108 * t107 - t121) + (m(7) * t75 + t109 + t31) * t87) * qJD(4) + (t127 * t47 - t128 * t45 - t19 * t61 - t21 * t62) * mrSges(7,3); (t12 * t35 + t13 * t34) * t126 - t128 * t33 + t62 * t29 - t127 * t32 - t61 * t28 + 0.2e1 * t75 * t27 + 0.2e1 * qJ(5) * t129 + 0.2e1 * t95 * qJD(5) + 0.2e1 * (-t12 * t61 - t127 * t35 + t128 * t34 - t13 * t62) * mrSges(7,3); m(7) * t44 - m(6) * t101 + (m(6) * t84 - t93) * t105 + t5; 0; (m(6) + m(7)) * t105; t27; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t100; t27; mrSges(7,1) * t21 - mrSges(7,2) * t19; mrSges(7,1) * t13 - mrSges(7,2) * t12 + t110; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
