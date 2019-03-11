% Calculate time derivative of joint inertia matrix for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:52:54
% EndTime: 2019-03-09 01:52:58
% DurationCPUTime: 1.54s
% Computational Cost: add. (2857->284), mult. (6072->439), div. (0->0), fcn. (5856->8), ass. (0->118)
t100 = sin(pkin(9));
t102 = cos(pkin(9));
t137 = sin(qJ(4));
t138 = cos(qJ(4));
t79 = t137 * t100 - t138 * t102;
t147 = 2 * mrSges(5,3);
t80 = t138 * t100 + t137 * t102;
t146 = 0.2e1 * t80;
t103 = -pkin(1) - qJ(3);
t125 = -pkin(7) + t103;
t82 = t125 * t100;
t83 = t125 * t102;
t145 = -t137 * t82 + t138 * t83;
t101 = cos(pkin(10));
t104 = sin(qJ(6));
t105 = cos(qJ(6));
t99 = sin(pkin(10));
t78 = t101 * t105 - t104 * t99;
t73 = t78 * qJD(6);
t107 = (t100 ^ 2 + t102 ^ 2) * qJD(3);
t144 = 2 * m(6);
t143 = 2 * m(7);
t95 = t99 ^ 2;
t142 = -2 * mrSges(5,3);
t141 = -0.2e1 * t145;
t91 = t100 * pkin(3) + qJ(2);
t140 = 0.2e1 * t91;
t109 = qJD(4) * t137;
t110 = qJD(4) * t138;
t75 = -t100 * t110 - t102 * t109;
t118 = t101 * t75;
t129 = t75 * t99;
t42 = mrSges(6,1) * t129 + mrSges(6,2) * t118;
t81 = t104 * t101 + t105 * t99;
t46 = t81 * t79;
t22 = qJD(6) * t46 + t75 * t78;
t24 = t73 * t79 - t75 * t81;
t7 = -t24 * mrSges(7,1) + t22 * mrSges(7,2);
t139 = t42 + t7;
t136 = Ifges(6,4) * t99;
t76 = -t100 * t109 + t102 * t110;
t135 = Ifges(6,5) * t76;
t134 = Ifges(6,6) * t76;
t133 = Ifges(6,6) * t99;
t37 = -t80 * qJD(3) + t145 * qJD(4);
t132 = t37 * t80;
t60 = t137 * t83 + t138 * t82;
t38 = -t79 * qJD(3) + t60 * qJD(4);
t131 = t38 * t145;
t130 = t38 * t79;
t64 = t79 * t75;
t128 = t79 * t99;
t127 = t99 * mrSges(6,3);
t126 = -mrSges(6,1) * t101 + mrSges(6,2) * t99 - mrSges(5,1);
t124 = pkin(8) + qJ(5);
t35 = pkin(4) * t76 - qJ(5) * t75 + qJD(5) * t79 + qJD(2);
t11 = t101 * t37 + t99 * t35;
t52 = pkin(4) * t80 + qJ(5) * t79 + t91;
t27 = t101 * t60 + t99 * t52;
t74 = t81 * qJD(6);
t123 = Ifges(7,5) * t73 - Ifges(7,6) * t74;
t122 = t101 ^ 2 + t95;
t120 = Ifges(6,4) * t101;
t119 = t101 * Ifges(6,1);
t117 = t101 * t79;
t116 = Ifges(7,5) * t22 + Ifges(7,6) * t24 + Ifges(7,3) * t76;
t114 = t122 * mrSges(6,3);
t113 = t122 * t80;
t10 = t101 * t35 - t37 * t99;
t26 = t101 * t52 - t60 * t99;
t108 = t122 * qJ(5);
t106 = t145 * t75 + t130;
t12 = pkin(5) * t80 + pkin(8) * t117 + t26;
t17 = pkin(8) * t128 + t27;
t3 = -t104 * t17 + t105 * t12;
t4 = t104 * t12 + t105 * t17;
t84 = t124 * t99;
t86 = t124 * t101;
t65 = -t104 * t86 - t105 * t84;
t66 = -t104 * t84 + t105 * t86;
t93 = -pkin(5) * t101 - pkin(4);
t69 = t75 * mrSges(5,2);
t63 = Ifges(7,1) * t81 + Ifges(7,4) * t78;
t62 = Ifges(7,4) * t81 + Ifges(7,2) * t78;
t61 = -mrSges(7,1) * t78 + mrSges(7,2) * t81;
t58 = mrSges(6,1) * t80 + mrSges(6,3) * t117;
t57 = -mrSges(6,2) * t80 + t79 * t127;
t55 = Ifges(7,1) * t73 - Ifges(7,4) * t74;
t54 = Ifges(7,4) * t73 - Ifges(7,2) * t74;
t53 = mrSges(7,1) * t74 + mrSges(7,2) * t73;
t51 = (-mrSges(6,1) * t99 - mrSges(6,2) * t101) * t79;
t50 = mrSges(6,1) * t76 - mrSges(6,3) * t118;
t49 = -mrSges(6,2) * t76 - t75 * t127;
t48 = t78 * t79;
t47 = t78 * t80;
t45 = t81 * t80;
t41 = -qJD(5) * t81 - qJD(6) * t66;
t40 = qJD(5) * t78 + qJD(6) * t65;
t39 = -pkin(5) * t128 - t145;
t32 = mrSges(7,1) * t80 + mrSges(7,3) * t48;
t31 = -mrSges(7,2) * t80 + mrSges(7,3) * t46;
t30 = t135 + (t119 - t136) * t75;
t29 = t134 + (-t99 * Ifges(6,2) + t120) * t75;
t28 = pkin(5) * t129 + t38;
t25 = -mrSges(7,1) * t46 - mrSges(7,2) * t48;
t23 = -t80 * t73 - t76 * t81;
t21 = -t74 * t80 + t76 * t78;
t16 = -Ifges(7,1) * t48 + Ifges(7,4) * t46 + Ifges(7,5) * t80;
t15 = -Ifges(7,4) * t48 + Ifges(7,2) * t46 + Ifges(7,6) * t80;
t14 = -mrSges(7,2) * t76 + mrSges(7,3) * t24;
t13 = mrSges(7,1) * t76 - mrSges(7,3) * t22;
t9 = -pkin(8) * t129 + t11;
t8 = pkin(5) * t76 - pkin(8) * t118 + t10;
t6 = Ifges(7,1) * t22 + Ifges(7,4) * t24 + Ifges(7,5) * t76;
t5 = Ifges(7,4) * t22 + Ifges(7,2) * t24 + Ifges(7,6) * t76;
t2 = -qJD(6) * t4 - t104 * t9 + t105 * t8;
t1 = t3 * qJD(6) + t104 * t8 + t105 * t9;
t18 = [t80 * t116 + 0.2e1 * m(5) * (qJD(2) * t91 + t37 * t60 - t131) + (t10 * t26 + t11 * t27 - t131) * t144 + 0.2e1 * m(4) * (qJ(2) * qJD(2) - t103 * t107) + (mrSges(5,1) * t140 + t60 * t142 - Ifges(7,5) * t48 + Ifges(7,6) * t46 + ((2 * Ifges(5,2)) + (2 * Ifges(6,3)) + Ifges(7,3)) * t80 + (-Ifges(6,5) * t101 + (2 * Ifges(5,4)) + t133) * t79) * t76 + (mrSges(5,3) * t141 + (-t95 * Ifges(6,2) - (2 * Ifges(5,1))) * t79 + (Ifges(6,5) * t146 + (-t119 + 0.2e1 * t136) * t79) * t101 + (-Ifges(5,4) - t133) * t146) * t75 + 0.2e1 * (m(3) * qJ(2) + t100 * mrSges(4,1) + t80 * mrSges(5,1) + t102 * mrSges(4,2) - t79 * mrSges(5,2) + mrSges(3,3)) * qJD(2) + (-t130 - t132) * t147 + (-t101 * t30 + t99 * t29) * t79 + 0.2e1 * mrSges(4,3) * t107 + 0.2e1 * t3 * t13 + 0.2e1 * t4 * t14 + t22 * t16 + t24 * t15 + t69 * t140 + t42 * t141 + (t1 * t4 + t2 * t3 + t28 * t39) * t143 + 0.2e1 * t28 * t25 + 0.2e1 * t1 * t31 + 0.2e1 * t2 * t32 + 0.2e1 * t39 * t7 + t46 * t5 - t48 * t6 + 0.2e1 * t27 * t49 + 0.2e1 * t26 * t50 + 0.2e1 * t38 * t51 + 0.2e1 * t11 * t57 + 0.2e1 * t10 * t58; -t45 * t13 + t47 * t14 + t21 * t31 + t23 * t32 + (t101 * t49 - t99 * t50) * t80 + t139 * t79 + (t101 * t57 + t80 * t142 - t99 * t58) * t76 + (t79 * t147 - t25 - t51) * t75 + m(7) * (t1 * t47 - t2 * t45 + t21 * t4 + t23 * t3 + t28 * t79 - t39 * t75) + m(6) * ((-t10 * t80 - t26 * t76) * t99 + (t11 * t80 + t27 * t76) * t101 + t106) + m(5) * (t60 * t76 + t106 + t132) - m(4) * t107; 0.2e1 * m(7) * (t21 * t47 - t23 * t45 - t64) + 0.2e1 * m(6) * (t76 * t113 - t64) + 0.2e1 * m(5) * (t76 * t80 - t64); t76 * mrSges(5,1) + t101 * t50 + t78 * t13 + t81 * t14 + t73 * t31 - t74 * t32 + t99 * t49 + t69 + (m(5) + m(4)) * qJD(2) + m(7) * (t1 * t81 + t2 * t78 - t3 * t74 + t4 * t73) + m(6) * (t10 * t101 + t11 * t99); m(7) * (t21 * t81 + t23 * t78 + t45 * t74 + t47 * t73); (t73 * t81 - t74 * t78) * t143; (t1 * t78 - t2 * t81 - t3 * t73 - t4 * t74) * mrSges(7,3) + m(6) * (-pkin(4) * t38 + (t101 * t27 - t26 * t99) * qJD(5) + (-t10 * t99 + t101 * t11) * qJ(5)) + m(7) * (t1 * t66 + t2 * t65 + t28 * t93 + t3 * t41 + t4 * t40) + t80 * t123 / 0.2e1 + t126 * t38 + (qJ(5) * t49 + qJD(5) * t57 + t11 * mrSges(6,3) + t29 / 0.2e1 + t134 / 0.2e1) * t101 + (-t10 * mrSges(6,3) - qJ(5) * t50 - qJD(5) * t58 + t30 / 0.2e1 + t135 / 0.2e1) * t99 + (-t99 * (Ifges(6,2) * t101 + t136) / 0.2e1 + t101 * (Ifges(6,1) * t99 + t120) / 0.2e1 + Ifges(5,5)) * t75 - t37 * mrSges(5,2) + t40 * t31 + t41 * t32 - pkin(4) * t42 + t39 * t53 + t46 * t54 / 0.2e1 - t48 * t55 / 0.2e1 + t28 * t61 + t24 * t62 / 0.2e1 + t22 * t63 / 0.2e1 + t65 * t13 + t66 * t14 + t73 * t16 / 0.2e1 - t74 * t15 / 0.2e1 - Ifges(5,6) * t76 + t78 * t5 / 0.2e1 + t81 * t6 / 0.2e1 + t76 * (Ifges(7,5) * t81 + Ifges(7,6) * t78) / 0.2e1 + t93 * t7; t79 * t53 + (-mrSges(5,2) + t114) * t76 + (-t61 - t126) * t75 + m(7) * (t21 * t66 + t23 * t65 + t40 * t47 - t41 * t45 - t75 * t93) + m(6) * (pkin(4) * t75 + qJD(5) * t113 + t76 * t108) + (t21 * t78 - t23 * t81 + t45 * t73 - t47 * t74) * mrSges(7,3); m(7) * (t40 * t81 + t41 * t78 - t65 * t74 + t66 * t73); t73 * t63 + t81 * t55 - t74 * t62 + t78 * t54 + 0.2e1 * t93 * t53 + (t40 * t66 + t41 * t65) * t143 + 0.2e1 * (t40 * t78 - t41 * t81 - t65 * t73 - t66 * t74) * mrSges(7,3) + (t108 * t144 + 0.2e1 * t114) * qJD(5); m(6) * t38 + m(7) * t28 + t139; 0.2e1 * (-m(7) / 0.2e1 - m(6) / 0.2e1) * t75; 0; t53; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t116; mrSges(7,1) * t23 - mrSges(7,2) * t21; -t53; mrSges(7,1) * t41 - mrSges(7,2) * t40 + t123; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t18(1) t18(2) t18(4) t18(7) t18(11) t18(16); t18(2) t18(3) t18(5) t18(8) t18(12) t18(17); t18(4) t18(5) t18(6) t18(9) t18(13) t18(18); t18(7) t18(8) t18(9) t18(10) t18(14) t18(19); t18(11) t18(12) t18(13) t18(14) t18(15) t18(20); t18(16) t18(17) t18(18) t18(19) t18(20) t18(21);];
Mq  = res;
