% Calculate time derivative of joint inertia matrix for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:01:22
% EndTime: 2019-03-09 03:01:27
% DurationCPUTime: 1.94s
% Computational Cost: add. (2404->319), mult. (5208->451), div. (0->0), fcn. (4571->8), ass. (0->143)
t148 = Ifges(6,5) + Ifges(7,5);
t169 = Ifges(6,3) + Ifges(7,3);
t94 = cos(qJ(5));
t133 = qJD(5) * t94;
t89 = sin(pkin(10));
t90 = cos(pkin(10));
t93 = sin(qJ(3));
t95 = cos(qJ(3));
t64 = t89 * t95 + t90 * t93;
t126 = t64 * t133;
t63 = t89 * t93 - t90 * t95;
t59 = t63 * qJD(3);
t92 = sin(qJ(5));
t151 = t92 * t59;
t101 = t126 - t151;
t81 = sin(pkin(9)) * pkin(1) + pkin(7);
t137 = qJ(4) + t81;
t115 = t137 * t93;
t62 = t137 * t95;
t38 = -t115 * t89 + t90 * t62;
t35 = t94 * t38;
t121 = -cos(pkin(9)) * pkin(1) - pkin(2);
t72 = -pkin(3) * t95 + t121;
t36 = t63 * pkin(4) - t64 * pkin(8) + t72;
t14 = t92 * t36 + t35;
t119 = t59 * (t92 ^ 2 + t94 ^ 2);
t80 = pkin(3) * t89 + pkin(8);
t120 = m(6) * t80 + mrSges(6,3);
t168 = 2 * m(5);
t167 = 2 * m(7);
t166 = -2 * mrSges(5,3);
t165 = -2 * mrSges(7,3);
t112 = qJD(3) * t137;
t48 = qJD(4) * t95 - t112 * t93;
t96 = -t93 * qJD(4) - t112 * t95;
t23 = t48 * t89 - t90 * t96;
t164 = 0.2e1 * t23;
t37 = t90 * t115 + t62 * t89;
t163 = 0.2e1 * t37;
t162 = m(7) * pkin(5);
t160 = mrSges(6,2) * t94;
t159 = Ifges(6,4) * t92;
t158 = Ifges(6,4) * t94;
t157 = Ifges(7,4) * t92;
t156 = Ifges(7,4) * t94;
t155 = t23 * t37;
t58 = t64 * qJD(3);
t154 = t58 * t63;
t153 = t64 * t92;
t152 = t64 * t94;
t149 = t94 * t59;
t147 = -Ifges(6,6) - Ifges(7,6);
t134 = qJD(5) * t92;
t127 = t64 * t134;
t100 = t127 + t149;
t17 = mrSges(7,1) * t58 + mrSges(7,3) * t100;
t18 = mrSges(6,1) * t58 + mrSges(6,3) * t100;
t146 = t17 + t18;
t19 = -mrSges(7,2) * t58 - mrSges(7,3) * t101;
t20 = -mrSges(6,2) * t58 - mrSges(6,3) * t101;
t145 = t19 + t20;
t104 = -Ifges(7,2) * t92 + t156;
t25 = t63 * Ifges(7,6) + t104 * t64;
t105 = -Ifges(6,2) * t92 + t158;
t26 = t63 * Ifges(6,6) + t105 * t64;
t144 = -t25 - t26;
t106 = Ifges(7,1) * t94 - t157;
t27 = Ifges(7,5) * t63 + t106 * t64;
t107 = Ifges(6,1) * t94 - t159;
t28 = Ifges(6,5) * t63 + t107 * t64;
t143 = t27 + t28;
t41 = -mrSges(7,2) * t63 - mrSges(7,3) * t153;
t42 = -mrSges(6,2) * t63 - mrSges(6,3) * t153;
t142 = t41 + t42;
t43 = mrSges(7,1) * t63 - mrSges(7,3) * t152;
t44 = mrSges(6,1) * t63 - mrSges(6,3) * t152;
t141 = -t43 - t44;
t140 = -mrSges(6,1) * t94 + mrSges(6,2) * t92 - mrSges(5,1);
t65 = mrSges(7,1) * t134 + mrSges(7,2) * t133;
t138 = qJ(6) * t64;
t136 = qJ(6) + t80;
t135 = qJD(5) * t64;
t24 = t90 * t48 + t89 * t96;
t129 = pkin(3) * qJD(3) * t93;
t34 = pkin(4) * t58 + pkin(8) * t59 + t129;
t132 = t36 * t133 + t94 * t24 + t92 * t34;
t131 = -t101 * mrSges(7,1) + mrSges(7,2) * t149;
t128 = pkin(5) * t134;
t125 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t75 = Ifges(7,2) * t94 + t157;
t76 = Ifges(6,2) * t94 + t159;
t124 = -t75 / 0.2e1 - t76 / 0.2e1;
t77 = Ifges(7,1) * t92 + t156;
t78 = Ifges(6,1) * t92 + t158;
t123 = t77 / 0.2e1 + t78 / 0.2e1;
t122 = -mrSges(6,1) - t162;
t82 = -pkin(3) * t90 - pkin(4);
t118 = t147 * t92;
t55 = t59 * mrSges(5,2);
t117 = t58 * mrSges(5,1) - t55;
t116 = -t24 * t92 + t94 * t34;
t13 = t94 * t36 - t38 * t92;
t114 = 0.2e1 * t129;
t113 = -t148 * t149 + t169 * t58;
t111 = qJD(5) * t136;
t110 = -(2 * Ifges(5,4)) + t118;
t109 = mrSges(4,1) * t93 + mrSges(4,2) * t95;
t108 = mrSges(6,1) * t92 + t160;
t103 = t23 * t63 + t37 * t58;
t102 = qJ(6) * t59 - qJD(6) * t64;
t15 = -mrSges(7,2) * t127 - t131;
t99 = t122 * t92 - t160;
t98 = t147 * t94 - t148 * t92;
t97 = t141 * t92 + t142 * t94;
t86 = Ifges(6,5) * t133;
t85 = Ifges(7,5) * t133;
t73 = -mrSges(7,1) * t94 + mrSges(7,2) * t92;
t71 = -pkin(5) * t94 + t82;
t70 = t107 * qJD(5);
t69 = t106 * qJD(5);
t68 = t105 * qJD(5);
t67 = t104 * qJD(5);
t66 = t108 * qJD(5);
t61 = t136 * t94;
t60 = t136 * t92;
t47 = -qJD(6) * t92 - t111 * t94;
t46 = qJD(6) * t94 - t111 * t92;
t40 = t108 * t64;
t39 = (mrSges(7,1) * t92 + mrSges(7,2) * t94) * t64;
t22 = pkin(5) * t153 + t37;
t16 = mrSges(6,1) * t101 - mrSges(6,2) * t100;
t12 = pkin(5) * t101 + t23;
t10 = -Ifges(6,1) * t100 - Ifges(6,4) * t101 + t58 * Ifges(6,5);
t9 = -Ifges(7,1) * t100 - Ifges(7,4) * t101 + t58 * Ifges(7,5);
t8 = -Ifges(6,4) * t100 - Ifges(6,2) * t101 + t58 * Ifges(6,6);
t7 = -Ifges(7,4) * t100 - Ifges(7,2) * t101 + t58 * Ifges(7,6);
t6 = -t138 * t92 + t14;
t5 = pkin(5) * t63 - t138 * t94 + t13;
t4 = -t14 * qJD(5) + t116;
t3 = -t134 * t38 + t132;
t2 = -qJ(6) * t126 + (-qJD(5) * t38 + t102) * t92 + t132;
t1 = pkin(5) * t58 + t102 * t94 + (-t35 + (-t36 + t138) * t92) * qJD(5) + t116;
t11 = [0.2e1 * t72 * t117 + 0.2e1 * t12 * t39 + t40 * t164 + 0.2e1 * t2 * t41 + 0.2e1 * t3 * t42 + 0.2e1 * t1 * t43 + 0.2e1 * t4 * t44 + t16 * t163 + 0.2e1 * t5 * t17 + 0.2e1 * t13 * t18 + 0.2e1 * t6 * t19 + 0.2e1 * t14 * t20 + 0.2e1 * t22 * t15 + t38 * t58 * t166 + (t1 * t5 + t12 * t22 + t2 * t6) * t167 + 0.2e1 * m(6) * (t13 * t4 + t14 * t3 + t155) + (t129 * t72 + t24 * t38 + t155) * t168 - (mrSges(5,3) * t163 + t143 * t94 + t144 * t92) * t59 + (mrSges(5,1) * t114 + t24 * t166 - t110 * t59 + ((2 * Ifges(5,2)) + t169) * t58 + t113) * t63 + (mrSges(5,2) * t114 + mrSges(5,3) * t164 - 0.2e1 * Ifges(5,1) * t59 + (t10 + t9) * t94 + (-t7 - t8) * t92 + (t148 * t94 + t110) * t58 + (-t143 * t92 + t144 * t94 + t63 * t98) * qJD(5)) * t64 + 0.2e1 * (t109 * t121 + (-t93 ^ 2 + t95 ^ 2) * Ifges(4,4) + (Ifges(4,1) - Ifges(4,2)) * t93 * t95) * qJD(3); (t15 + t16) * t63 + (t39 + t40) * t58 - t97 * t59 + m(7) * (t12 * t63 - t149 * t6 + t151 * t5 + t22 * t58) + m(5) * (-t38 * t59 + t103) + m(6) * (t13 * t151 - t14 * t149 + t103) + (t145 * t94 - t146 * t92 + (t141 * t94 - t142 * t92) * qJD(5) + m(7) * (-t1 * t92 - t133 * t5 - t134 * t6 + t2 * t94) + m(5) * t24 + m(6) * (-t13 * t133 - t134 * t14 + t3 * t94 - t4 * t92)) * t64; (-t59 * t64 + t154) * t168 + 0.4e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * (-t119 * t64 + t154); t71 * t15 + t12 * t73 + t82 * t16 + t22 * t65 + t37 * t66 - Ifges(5,6) * t58 - Ifges(5,5) * t59 - t60 * t17 + t61 * t19 + t46 * t41 + t47 * t43 - t24 * mrSges(5,2) + m(7) * (-t1 * t60 + t12 * t71 + t2 * t61 + t46 * t6 + t47 * t5) + (t86 / 0.2e1 + t85 / 0.2e1) * t63 + (m(6) * t82 + t140) * t23 + (Ifges(4,5) * t95 - Ifges(4,6) * t93 + (-mrSges(4,1) * t95 + mrSges(4,2) * t93) * t81) * qJD(3) + (m(5) * (-t23 * t90 + t24 * t89) + (-t89 * t58 + t90 * t59) * mrSges(5,3)) * pkin(3) + (t9 / 0.2e1 + t10 / 0.2e1 - t1 * mrSges(7,3) - t80 * t18 + (-t67 / 0.2e1 - t68 / 0.2e1) * t64 - t124 * t59 + (Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1) * t58 - t120 * t4 + (-t80 * t42 + pkin(5) * t39 - t6 * mrSges(7,3) - t25 / 0.2e1 - t26 / 0.2e1 - t123 * t64 - t125 * t63 + t22 * t162 - t120 * t14) * qJD(5)) * t92 + (t7 / 0.2e1 + t8 / 0.2e1 + t2 * mrSges(7,3) + t80 * t20 + (t69 / 0.2e1 + t70 / 0.2e1) * t64 - t123 * t59 + t125 * t58 + t120 * t3 + (-t80 * t44 - t5 * mrSges(7,3) + t27 / 0.2e1 + t28 / 0.2e1 + t124 * t64 - t120 * t13) * qJD(5)) * t94; t55 + (t65 + t66) * t63 - t109 * qJD(3) + (t73 + t140) * t58 + m(6) * (-t119 * t80 + t58 * t82) + m(5) * (-t58 * t90 - t59 * t89) * pkin(3) - (mrSges(6,3) + mrSges(7,3)) * t119 + (t128 * t63 - t149 * t61 - t151 * t60 + t58 * t71 + (t46 * t94 - t47 * t92 + (t60 * t94 - t61 * t92) * qJD(5)) * t64) * m(7); 0.2e1 * t71 * t65 + (t46 * t61 - t47 * t60) * t167 + 0.2e1 * t82 * t66 + (t47 * t165 + t69 + t70 + (t61 * t165 - t75 - t76 + 0.2e1 * (m(7) * t71 + t73) * pkin(5)) * qJD(5)) * t92 + (0.2e1 * t46 * mrSges(7,3) + t67 + t68 + (-t165 * t60 + t77 + t78) * qJD(5)) * t94; m(5) * t129 + t146 * t94 + t145 * t92 + t97 * qJD(5) + m(7) * (t1 * t94 + t2 * t92 + (-t5 * t92 + t6 * t94) * qJD(5)) + m(6) * (t3 * t92 + t4 * t94 + (-t13 * t92 + t14 * t94) * qJD(5)) + t117; 0; m(7) * (t46 * t92 + t47 * t94 + (t60 * t92 + t61 * t94) * qJD(5)); 0; mrSges(6,1) * t4 + mrSges(7,1) * t1 - mrSges(6,2) * t3 - mrSges(7,2) * t2 - t59 * t118 + (m(7) * t1 + t17) * pkin(5) + t98 * t135 + t113; -t99 * t59 + (t122 * t94 + (mrSges(6,2) + mrSges(7,2)) * t92) * t135 + t131; -mrSges(7,2) * t46 + t85 + t86 + (mrSges(7,1) + t162) * t47 + ((-mrSges(6,1) * t80 - (mrSges(7,3) * pkin(5))) * t94 + (mrSges(6,2) * t80 + t147) * t92) * qJD(5); qJD(5) * t99 - t65; 0; m(7) * t12 + t15; m(7) * t58; m(7) * t128 + t65; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t11(1) t11(2) t11(4) t11(7) t11(11) t11(16); t11(2) t11(3) t11(5) t11(8) t11(12) t11(17); t11(4) t11(5) t11(6) t11(9) t11(13) t11(18); t11(7) t11(8) t11(9) t11(10) t11(14) t11(19); t11(11) t11(12) t11(13) t11(14) t11(15) t11(20); t11(16) t11(17) t11(18) t11(19) t11(20) t11(21);];
Mq  = res;
