% Calculate time derivative of joint inertia matrix for
% S6RPRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
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
% Datum: 2019-03-09 03:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:11:14
% EndTime: 2019-03-09 03:11:17
% DurationCPUTime: 1.86s
% Computational Cost: add. (1288->291), mult. (2750->406), div. (0->0), fcn. (1839->6), ass. (0->129)
t68 = sin(pkin(9)) * pkin(1) + pkin(7);
t146 = pkin(4) + t68;
t159 = Ifges(7,4) + Ifges(6,5);
t158 = Ifges(7,2) + Ifges(6,3);
t77 = sin(qJ(5));
t79 = cos(qJ(5));
t126 = t77 ^ 2 + t79 ^ 2;
t123 = qJD(3) * t79;
t78 = sin(qJ(3));
t115 = t78 * t123;
t80 = cos(qJ(3));
t118 = qJD(5) * t80;
t85 = t77 * t118 + t115;
t122 = qJD(3) * t80;
t124 = qJD(3) * t78;
t112 = t77 * t124;
t84 = -t118 * t79 + t112;
t20 = mrSges(6,1) * t122 - mrSges(6,3) * t84;
t119 = qJD(5) * t79;
t21 = mrSges(7,2) * t112 + (-qJD(3) * mrSges(7,1) - mrSges(7,2) * t119) * t80;
t131 = t20 - t21;
t19 = -mrSges(6,2) * t122 + mrSges(6,3) * t85;
t22 = mrSges(7,2) * t85 + mrSges(7,3) * t122;
t132 = t19 + t22;
t137 = t79 * t80;
t51 = -t78 * mrSges(6,2) - mrSges(6,3) * t137;
t52 = -mrSges(7,2) * t137 + t78 * mrSges(7,3);
t128 = t51 + t52;
t140 = t77 * t80;
t49 = t78 * mrSges(6,1) + mrSges(6,3) * t140;
t50 = -t78 * mrSges(7,1) - mrSges(7,2) * t140;
t129 = t49 - t50;
t83 = t128 * t79 - t129 * t77;
t157 = t83 * qJD(5) + t131 * t79 + t132 * t77;
t156 = 2 * qJ(4);
t108 = -cos(pkin(9)) * pkin(1) - pkin(2);
t155 = 0.2e1 * t108;
t81 = -pkin(3) - pkin(8);
t86 = -qJ(4) * t78 + t108;
t26 = t80 * t81 + t86;
t41 = t146 * t78;
t154 = t79 * t26 + t77 * t41;
t152 = m(7) * qJ(6) + mrSges(7,3);
t121 = qJD(4) * t78;
t105 = pkin(3) * t124 - t121;
t23 = (pkin(8) * t78 - qJ(4) * t80) * qJD(3) + t105;
t35 = t146 * t122;
t4 = -qJD(5) * t154 - t23 * t77 + t35 * t79;
t151 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t149 = 2 * mrSges(5,3);
t148 = 0.2e1 * pkin(3) * t80 - 0.2e1 * t86;
t36 = qJ(4) * t122 - t105;
t147 = m(5) * t36;
t145 = Ifges(6,4) * t77;
t144 = Ifges(6,4) * t79;
t143 = Ifges(7,5) * t77;
t142 = Ifges(7,5) * t79;
t141 = Ifges(7,6) * t78;
t139 = t77 * t81;
t138 = t78 * Ifges(6,6);
t136 = t79 * t81;
t135 = mrSges(4,2) - mrSges(5,3);
t134 = mrSges(5,2) - mrSges(4,1);
t95 = Ifges(7,1) * t77 - t142;
t29 = t78 * Ifges(7,4) - t80 * t95;
t96 = Ifges(6,1) * t77 + t144;
t30 = t78 * Ifges(6,5) - t80 * t96;
t130 = t29 + t30;
t127 = t126 * t81 * t124;
t42 = t146 * t80;
t125 = qJD(3) * t77;
t120 = qJD(5) * t77;
t117 = qJD(6) * t77;
t116 = qJ(6) * qJD(3);
t111 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t58 = Ifges(7,3) * t77 + t142;
t59 = -Ifges(6,2) * t77 + t144;
t110 = -t58 / 0.2e1 + t59 / 0.2e1;
t60 = Ifges(7,1) * t79 + t143;
t61 = Ifges(6,1) * t79 - t145;
t109 = t60 / 0.2e1 + t61 / 0.2e1;
t107 = m(5) * t68 + mrSges(5,1);
t91 = -pkin(5) * t77 + qJ(6) * t79;
t53 = qJ(4) - t91;
t56 = t77 * mrSges(7,1) - t79 * mrSges(7,3);
t106 = m(7) * t53 + t56;
t93 = -Ifges(7,3) * t79 + t143;
t27 = -t80 * t93 + t141;
t94 = Ifges(6,2) * t79 + t145;
t28 = -t80 * t94 + t138;
t104 = -t27 + t28 - t141;
t3 = t41 * t119 - t120 * t26 + t79 * t23 + t77 * t35;
t1 = qJD(6) * t78 + t116 * t80 + t3;
t2 = -pkin(5) * t122 - t4;
t103 = t1 * t77 - t2 * t79;
t102 = t3 * t77 + t4 * t79;
t6 = qJ(6) * t78 + t154;
t8 = -t26 * t77 + t41 * t79;
t7 = -pkin(5) * t78 - t8;
t101 = t6 * t79 + t7 * t77;
t100 = t154 * t79 - t77 * t8;
t99 = (m(7) * t81 - mrSges(7,2)) * t77;
t98 = mrSges(6,1) * t79 - mrSges(6,2) * t77;
t57 = t77 * mrSges(6,1) + t79 * mrSges(6,2);
t97 = mrSges(7,1) * t79 + mrSges(7,3) * t77;
t92 = pkin(5) * t79 + qJ(6) * t77;
t88 = (m(7) / 0.2e1 + m(6) / 0.2e1) * t124;
t87 = Ifges(6,6) * t85 + t112 * t159 + t122 * t158;
t82 = m(7) * t91 - t56 - t57;
t70 = Ifges(7,6) * t119;
t48 = t96 * qJD(5);
t47 = t95 * qJD(5);
t46 = t94 * qJD(5);
t45 = t93 * qJD(5);
t44 = t98 * qJD(5);
t43 = t97 * qJD(5);
t40 = t98 * t80;
t39 = t97 * t80;
t34 = t146 * t124;
t32 = qJD(5) * t92 - qJD(6) * t79 + qJD(4);
t16 = t80 * t92 + t42;
t15 = -mrSges(6,1) * t85 + mrSges(6,2) * t84;
t14 = -mrSges(7,1) * t85 - mrSges(7,3) * t84;
t13 = -t61 * t118 + (Ifges(6,5) * t80 + t78 * t96) * qJD(3);
t12 = -t60 * t118 + (Ifges(7,4) * t80 + t78 * t95) * qJD(3);
t11 = -t59 * t118 + (Ifges(6,6) * t80 + t78 * t94) * qJD(3);
t10 = -t58 * t118 + (Ifges(7,6) * t80 + t78 * t93) * qJD(3);
t5 = (qJD(5) * t91 + t117) * t80 + (-t92 - t146) * t124;
t9 = [t147 * t148 + 0.2e1 * t1 * t52 + 0.2e1 * t16 * t14 + 0.2e1 * t42 * t15 + 0.2e1 * t154 * t19 + 0.2e1 * t2 * t50 + 0.2e1 * t8 * t20 + 0.2e1 * t7 * t21 + 0.2e1 * t6 * t22 + 0.2e1 * t3 * t51 - 0.2e1 * t34 * t40 + 0.2e1 * t5 * t39 + 0.2e1 * t4 * t49 + 0.2e1 * m(6) * (t154 * t3 - t34 * t42 + t4 * t8) + 0.2e1 * m(7) * (t1 * t6 + t16 * t5 + t2 * t7) + (t36 * t149 + (mrSges(4,1) * t155 + mrSges(5,2) * t148 + 0.2e1 * (-Ifges(4,4) - Ifges(5,6)) * t78 + t104 * t79 + t130 * t77) * qJD(3) + t87) * t78 + (-0.2e1 * t36 * mrSges(5,2) + (t10 - t11) * t79 + (-t12 - t13) * t77 + (t104 * t77 + (-t159 * t78 - t130) * t79) * qJD(5) + (mrSges(4,2) * t155 + mrSges(5,3) * t148 + (0.2e1 * Ifges(4,4) + 0.2e1 * Ifges(5,6) + (-Ifges(6,6) + Ifges(7,6)) * t79 - t159 * t77) * t80 + (-(2 * Ifges(5,3)) + (2 * Ifges(5,2)) + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) + t158) * t78) * qJD(3)) * t80; (t14 + t15 + (t128 * t77 + t129 * t79) * qJD(3) + m(7) * (-t123 * t7 + t125 * t6 + t5) + m(6) * (t123 * t8 + t125 * t154 - t34)) * t78 + ((t39 + t40) * qJD(3) + m(7) * (qJD(3) * t16 - t119 * t6 - t120 * t7 - t103) + m(6) * (qJD(3) * t42 - t119 * t154 + t120 * t8 - t102) - t157) * t80; 0.4e1 * (0.1e1 - t126) * t80 * t88; t78 * t70 / 0.2e1 + t5 * t56 - t34 * t57 + t32 * t39 + qJD(4) * t40 + t16 * t43 + t42 * t44 + t53 * t14 + qJ(4) * t15 + (t12 / 0.2e1 + t13 / 0.2e1 + t2 * mrSges(7,2) - t4 * mrSges(6,3) + t131 * t81) * t79 + (t10 / 0.2e1 - t11 / 0.2e1 - t1 * mrSges(7,2) - t3 * mrSges(6,3) + t132 * t81) * t77 + m(6) * (-qJ(4) * t34 + qJD(4) * t42 + t136 * t4 + t139 * t3) + m(7) * (t1 * t139 - t136 * t2 + t32 * t16 + t53 * t5) + ((-t45 / 0.2e1 + t46 / 0.2e1) * t79 + (t47 / 0.2e1 + t48 / 0.2e1) * t77 + t107 * qJD(4)) * t80 + ((-t138 / 0.2e1 + t27 / 0.2e1 - t28 / 0.2e1 - t6 * mrSges(7,2) - t154 * mrSges(6,3) - t109 * t80) * t79 + (-t29 / 0.2e1 - t30 / 0.2e1 - t7 * mrSges(7,2) + t8 * mrSges(6,3) + t110 * t80 - t111 * t78) * t77 + (m(6) * t100 + m(7) * t101 + t83) * t81) * qJD(5) + ((-pkin(3) * mrSges(5,1) - Ifges(5,4) + Ifges(4,5) + t111 * t79 + (-Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1) * t77 + (-m(5) * pkin(3) + t134) * t68) * t80 + (-qJ(4) * mrSges(5,1) + Ifges(5,5) - Ifges(4,6) + t110 * t79 + t109 * t77 + (-m(5) * qJ(4) + t135) * t68) * t78) * qJD(3); (t43 + t44) * t78 + m(7) * (t32 * t78 + t127) + t147 + m(6) * (t121 + t127) + ((m(6) * qJ(4) + t106 - t135 + t57) * t80 + (t134 + (-mrSges(7,2) - mrSges(6,3)) * t126) * t78) * qJD(3); t44 * t156 + 0.2e1 * t43 * t53 + (-t47 - t48) * t79 + (-t45 + t46) * t77 + ((t58 - t59) * t79 + (-t60 - t61) * t77) * qJD(5) + 0.2e1 * t106 * t32 + (t149 + 0.2e1 * t57 + (m(5) + m(6)) * t156) * qJD(4); t107 * t122 + m(7) * (qJD(5) * t101 + t103) + m(6) * (qJD(5) * t100 + t102) + t157; m(5) * t124 + 0.2e1 * t126 * t88; 0; 0; -Ifges(7,6) * t115 - pkin(5) * t21 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t6) + qJD(6) * t52 + qJ(6) * t22 + t1 * mrSges(7,3) + t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + (-Ifges(7,6) * t77 - t159 * t79) * t118 + t87; m(7) * (-qJD(6) * t80 + t116 * t78) * t77 + (t151 * t79 + (mrSges(7,3) - mrSges(6,2)) * t77) * t124 + ((mrSges(6,2) - t152) * t79 + t151 * t77) * t118; t70 + qJD(6) * t99 + ((-qJ(6) * mrSges(7,2) - Ifges(6,6)) * t79 + (mrSges(7,2) * pkin(5) - t159) * t77 + t82 * t81) * qJD(5); m(7) * t117 + qJD(5) * t82; 0.2e1 * t152 * qJD(6); m(7) * t2 + t21; -t85 * m(7); qJD(5) * t99; m(7) * t120; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t9(1) t9(2) t9(4) t9(7) t9(11) t9(16); t9(2) t9(3) t9(5) t9(8) t9(12) t9(17); t9(4) t9(5) t9(6) t9(9) t9(13) t9(18); t9(7) t9(8) t9(9) t9(10) t9(14) t9(19); t9(11) t9(12) t9(13) t9(14) t9(15) t9(20); t9(16) t9(17) t9(18) t9(19) t9(20) t9(21);];
Mq  = res;
