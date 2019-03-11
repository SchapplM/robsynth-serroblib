% Calculate time derivative of joint inertia matrix for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPPRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:14:01
% EndTime: 2019-03-08 19:14:03
% DurationCPUTime: 1.33s
% Computational Cost: add. (1950->233), mult. (5012->362), div. (0->0), fcn. (5067->12), ass. (0->113)
t65 = sin(pkin(11));
t131 = pkin(2) * t65;
t57 = qJ(4) + t131;
t64 = sin(pkin(12));
t67 = cos(pkin(12));
t133 = (m(5) * t57 + mrSges(5,3)) * (t64 ^ 2 + t67 ^ 2);
t68 = cos(pkin(11));
t130 = pkin(2) * t68;
t73 = cos(qJ(6));
t100 = qJD(6) * t73;
t122 = cos(qJ(5));
t71 = sin(qJ(5));
t76 = t122 * t67 - t64 * t71;
t45 = t76 * qJD(5);
t49 = t122 * t64 + t71 * t67;
t70 = sin(qJ(6));
t79 = -t49 * t100 - t45 * t70;
t114 = t45 * t73;
t101 = qJD(6) * t70;
t98 = t49 * t101;
t78 = t98 - t114;
t54 = -mrSges(7,1) * t73 + mrSges(7,2) * t70;
t129 = -m(7) * pkin(5) - mrSges(6,1) + t54;
t90 = t122 * qJD(4);
t123 = pkin(8) + t57;
t47 = t123 * t67;
t94 = t122 * t47;
t95 = t71 * t123;
t99 = t71 * qJD(4);
t19 = qJD(5) * t94 + t67 * t99 + (-qJD(5) * t95 + t90) * t64;
t128 = 0.2e1 * t19;
t127 = -t49 / 0.2e1;
t120 = Ifges(7,4) * t70;
t55 = Ifges(7,2) * t73 + t120;
t126 = -t55 / 0.2e1;
t46 = t49 * qJD(5);
t125 = pkin(5) * t46;
t66 = sin(pkin(6));
t72 = sin(qJ(2));
t74 = cos(qJ(2));
t40 = (t65 * t74 + t68 * t72) * t66;
t69 = cos(pkin(6));
t32 = -t40 * t64 + t67 * t69;
t33 = t40 * t67 + t64 * t69;
t15 = t122 * t33 + t71 * t32;
t105 = qJD(2) * t66;
t80 = t65 * t72 - t68 * t74;
t35 = t80 * t105;
t6 = qJD(5) * t15 - t35 * t49;
t77 = t122 * t32 - t33 * t71;
t124 = t77 * t6;
t121 = mrSges(7,3) * t49;
t119 = Ifges(7,4) * t73;
t118 = Ifges(7,6) * t70;
t108 = t71 * t47;
t89 = t123 * t122;
t24 = t64 * t89 + t108;
t117 = t19 * t24;
t34 = qJD(2) * t40;
t39 = t80 * t66;
t22 = t39 * t34;
t113 = t46 * mrSges(6,3);
t112 = t46 * Ifges(7,5);
t111 = t46 * Ifges(7,6);
t110 = t46 * t76;
t109 = t76 * Ifges(7,6);
t107 = Ifges(7,5) * t114 + Ifges(7,3) * t46;
t96 = -pkin(3) - t130;
t53 = -pkin(4) * t67 + t96;
t23 = -pkin(5) * t76 - pkin(9) * t49 + t53;
t25 = -t64 * t95 + t94;
t11 = t23 * t73 - t25 * t70;
t104 = qJD(6) * t11;
t12 = t23 * t70 + t25 * t73;
t103 = qJD(6) * t12;
t102 = qJD(6) * t49;
t93 = t45 * (t70 ^ 2 + t73 ^ 2);
t26 = t46 * mrSges(6,1) + t45 * mrSges(6,2);
t91 = -(2 * Ifges(6,4)) - t118;
t10 = t15 * t73 + t39 * t70;
t9 = -t15 * t70 + t39 * t73;
t88 = t10 * t73 - t70 * t9;
t87 = -t19 * t77 + t24 * t6;
t86 = -t46 * t77 - t6 * t76;
t85 = mrSges(7,1) * t70 + mrSges(7,2) * t73;
t84 = Ifges(7,1) * t73 - t120;
t83 = -Ifges(7,2) * t70 + t119;
t82 = -t19 * t76 + t24 * t46;
t81 = -t32 * t64 + t33 * t67;
t5 = qJD(5) * t77 - t35 * t76;
t1 = qJD(6) * t9 + t34 * t70 + t5 * t73;
t2 = -qJD(6) * t10 + t34 * t73 - t5 * t70;
t75 = t1 * t73 - t2 * t70 + (-t10 * t70 - t73 * t9) * qJD(6);
t59 = Ifges(7,5) * t100;
t56 = Ifges(7,1) * t70 + t119;
t52 = t84 * qJD(6);
t51 = t83 * qJD(6);
t50 = t85 * qJD(6);
t31 = -mrSges(7,1) * t76 - t121 * t73;
t30 = mrSges(7,2) * t76 - t121 * t70;
t28 = -pkin(9) * t45 + t125;
t27 = t85 * t49;
t21 = -Ifges(7,5) * t76 + t49 * t84;
t20 = t49 * t83 - t109;
t18 = t67 * t90 - qJD(5) * t108 + (-qJD(5) * t89 - t99) * t64;
t17 = -mrSges(7,2) * t46 + mrSges(7,3) * t79;
t16 = mrSges(7,1) * t46 + mrSges(7,3) * t78;
t13 = t79 * mrSges(7,1) + t78 * mrSges(7,2);
t8 = -Ifges(7,1) * t78 + Ifges(7,4) * t79 + t112;
t7 = -Ifges(7,4) * t78 + Ifges(7,2) * t79 + t111;
t4 = -t18 * t70 + t28 * t73 - t103;
t3 = t18 * t73 + t28 * t70 + t104;
t14 = [0.2e1 * m(7) * (t1 * t10 + t2 * t9 - t124) + 0.2e1 * m(6) * (t15 * t5 - t124 + t22) + 0.2e1 * m(5) * (-t35 * t81 + t22) + 0.2e1 * m(4) * (-t35 * t40 + t22); m(7) * (t1 * t12 + t10 * t3 + t11 * t2 + t4 * t9 + t87) + t1 * t30 + t10 * t17 + t2 * t31 + t9 * t16 + t6 * t27 + t77 * t13 + m(6) * (t15 * t18 + t25 * t5 + t87) - t15 * t113 + t39 * t26 + m(5) * t81 * qJD(4) + (-m(4) * t130 + m(5) * t96 + m(6) * t53 - mrSges(5,1) * t67 - mrSges(6,1) * t76 + mrSges(5,2) * t64 + mrSges(6,2) * t49 - mrSges(4,1)) * t34 + (-mrSges(3,1) * t72 - mrSges(3,2) * t74) * t105 + (-t45 * t77 + t49 * t6 + t5 * t76) * mrSges(6,3) + (-m(4) * t131 + mrSges(4,2) - t133) * t35; -0.2e1 * t25 * t113 + 0.2e1 * t11 * t16 + 0.2e1 * t12 * t17 - 0.2e1 * t24 * t13 + t27 * t128 + 0.2e1 * t53 * t26 + 0.2e1 * t3 * t30 + 0.2e1 * t4 * t31 + 0.2e1 * m(7) * (t11 * t4 + t12 * t3 + t117) + 0.2e1 * m(6) * (t18 * t25 + t117) + (0.2e1 * t24 * mrSges(6,3) - t70 * t20 + t73 * t21) * t45 + 0.2e1 * t133 * qJD(4) - (-0.2e1 * t18 * mrSges(6,3) + ((2 * Ifges(6,2)) + Ifges(7,3)) * t46 + t91 * t45 + t107) * t76 + (mrSges(6,3) * t128 + 0.2e1 * Ifges(6,1) * t45 - t70 * t7 + t73 * t8 + (Ifges(7,5) * t73 + t91) * t46 + (-t76 * (-Ifges(7,5) * t70 - Ifges(7,6) * t73) - t70 * t21 - t73 * t20) * qJD(6)) * t49; m(7) * (t45 * t88 + t49 * t75 + t86) + m(6) * (t15 * t45 + t49 * t5 + t86); t76 * t13 + t46 * t27 + m(7) * t82 + m(6) * (t18 * t49 + t25 * t45 + t82) + (m(7) * (-t102 * t11 + t12 * t45 + t3 * t49) + t45 * t30 + t49 * t17 - t31 * t102) * t73 + (m(7) * (-t102 * t12 - t11 * t45 - t4 * t49) - t30 * t102 - t45 * t31 - t49 * t16) * t70; 0.2e1 * m(6) * (t45 * t49 - t110) + 0.2e1 * m(7) * (t49 * t93 - t110); m(7) * (qJD(6) * t88 + t1 * t70 + t2 * t73) + 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * t34; m(7) * (t3 * t70 + t4 * t73 + (-t11 * t70 + t12 * t73) * qJD(6)) + t30 * t100 + t70 * t17 - t31 * t101 + t73 * t16 + t26; 0; 0; -t5 * mrSges(6,2) - t77 * t50 + (m(7) * pkin(9) + mrSges(7,3)) * t75 + t129 * t6; -t76 * t59 / 0.2e1 + t24 * t50 + Ifges(6,5) * t45 - Ifges(6,6) * t46 - t18 * mrSges(6,2) + pkin(5) * t13 + t129 * t19 + (t112 / 0.2e1 + t8 / 0.2e1 + t51 * t127 + t45 * t126 - t4 * mrSges(7,3) + (t109 / 0.2e1 - t20 / 0.2e1 + t56 * t127 - t12 * mrSges(7,3)) * qJD(6) + (m(7) * (-t4 - t103) - qJD(6) * t30 - t16) * pkin(9)) * t70 + (t111 / 0.2e1 + t7 / 0.2e1 + t49 * t52 / 0.2e1 + t45 * t56 / 0.2e1 + t3 * mrSges(7,3) + (t21 / 0.2e1 + t49 * t126 - t11 * mrSges(7,3)) * qJD(6) + (m(7) * (t3 - t104) + t17 - qJD(6) * t31) * pkin(9)) * t73; t46 * t54 - t76 * t50 + m(7) * (pkin(9) * t93 - t125) + mrSges(7,3) * t93 - t26; 0; -0.2e1 * pkin(5) * t50 + t51 * t73 + t52 * t70 + (-t70 * t55 + t73 * t56) * qJD(6); mrSges(7,1) * t2 - mrSges(7,2) * t1; mrSges(7,1) * t4 - mrSges(7,2) * t3 - Ifges(7,5) * t98 + Ifges(7,6) * t79 + t107; t13; -t50; t59 + (pkin(9) * t54 - t118) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
