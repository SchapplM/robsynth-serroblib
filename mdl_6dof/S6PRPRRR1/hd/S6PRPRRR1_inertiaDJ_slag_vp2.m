% Calculate time derivative of joint inertia matrix for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:22:37
% EndTime: 2019-03-08 20:22:43
% DurationCPUTime: 2.71s
% Computational Cost: add. (3273->303), mult. (8196->471), div. (0->0), fcn. (8022->12), ass. (0->139)
t102 = sin(qJ(6));
t106 = cos(qJ(6));
t165 = t102 ^ 2 + t106 ^ 2;
t184 = qJD(4) + qJD(5);
t103 = sin(qJ(5));
t104 = sin(qJ(4));
t107 = cos(qJ(5));
t108 = cos(qJ(4));
t80 = t103 * t104 - t107 * t108;
t58 = t184 * t80;
t188 = t165 * t58;
t185 = t165 * t107;
t87 = -mrSges(7,1) * t106 + mrSges(7,2) * t102;
t187 = t87 - mrSges(6,1);
t141 = qJD(6) * t106;
t81 = t103 * t108 + t104 * t107;
t119 = -t102 * t58 + t81 * t141;
t83 = (mrSges(5,1) * t104 + mrSges(5,2) * t108) * qJD(4);
t183 = 2 * m(7);
t182 = 2 * pkin(4);
t181 = -2 * mrSges(6,3);
t98 = sin(pkin(12));
t90 = pkin(2) * t98 + pkin(8);
t174 = pkin(9) + t90;
t133 = qJD(4) * t174;
t129 = t108 * t133;
t135 = t174 * t104;
t78 = t174 * t108;
t46 = -t103 * t135 + t107 * t78;
t75 = t104 * t133;
t25 = t46 * qJD(5) - t103 * t75 + t107 * t129;
t180 = 0.2e1 * t25;
t45 = t103 * t78 + t107 * t135;
t179 = 0.2e1 * t45;
t178 = m(6) / 0.2e1;
t177 = m(7) / 0.2e1;
t59 = t184 * t81;
t175 = pkin(5) * t59;
t101 = cos(pkin(6));
t100 = cos(pkin(12));
t105 = sin(qJ(2));
t109 = cos(qJ(2));
t99 = sin(pkin(6));
t73 = (t100 * t105 + t109 * t98) * t99;
t60 = t101 * t108 - t104 * t73;
t61 = t101 * t104 + t108 * t73;
t29 = t103 * t60 + t107 * t61;
t120 = t100 * t109 - t105 * t98;
t146 = qJD(2) * t99;
t70 = t120 * t146;
t32 = -t61 * qJD(4) - t104 * t70;
t33 = t60 * qJD(4) + t108 * t70;
t10 = t29 * qJD(5) + t103 * t33 - t107 * t32;
t28 = t103 * t61 - t107 * t60;
t173 = t10 * t28;
t72 = t120 * t99;
t17 = -t102 * t72 + t106 * t29;
t69 = qJD(2) * t73;
t9 = -t28 * qJD(5) + t103 * t32 + t107 * t33;
t3 = -t17 * qJD(6) - t102 * t9 + t106 * t69;
t172 = t102 * t3;
t171 = t45 * t25;
t91 = -pkin(2) * t100 - pkin(3);
t86 = -pkin(4) * t108 + t91;
t44 = t80 * pkin(5) - t81 * pkin(10) + t86;
t22 = t102 * t44 + t106 * t46;
t24 = -t45 * qJD(5) - t103 * t129 - t107 * t75;
t139 = pkin(4) * qJD(4) * t104;
t27 = t58 * pkin(10) + t139 + t175;
t6 = -qJD(6) * t22 - t102 * t24 + t106 * t27;
t168 = t6 * t102;
t41 = t72 * t69;
t167 = t80 * t59;
t152 = t106 * t58;
t166 = -Ifges(7,5) * t152 + Ifges(7,3) * t59;
t164 = mrSges(7,3) * t106;
t163 = Ifges(7,4) * t102;
t162 = Ifges(7,4) * t106;
t161 = Ifges(7,6) * t102;
t160 = pkin(4) * qJD(5);
t158 = t102 * t81;
t157 = t103 * mrSges(6,1);
t156 = t103 * t28;
t155 = t103 * t45;
t154 = t103 * t80;
t153 = t103 * t87;
t151 = t106 * t81;
t93 = pkin(4) * t103 + pkin(10);
t150 = t106 * t93;
t149 = t107 * mrSges(6,2);
t145 = qJD(6) * t81;
t144 = t102 * t107;
t143 = t106 * t107;
t142 = qJD(6) * t102;
t140 = 0.2e1 * t104;
t138 = t81 * t142;
t118 = t138 + t152;
t15 = t119 * mrSges(7,1) - t118 * mrSges(7,2);
t136 = m(7) * t25 + t15;
t30 = t59 * mrSges(6,1) - t58 * mrSges(6,2);
t131 = -(2 * Ifges(6,4)) - t161;
t130 = mrSges(7,3) * t185;
t128 = t45 * t10 + t25 * t28;
t127 = t80 * t10 + t59 * t28;
t126 = t80 * t25 + t59 * t45;
t124 = mrSges(7,1) * t102 + mrSges(7,2) * t106;
t123 = Ifges(7,1) * t106 - t163;
t122 = -Ifges(7,2) * t102 + t162;
t121 = Ifges(7,5) * t102 + Ifges(7,6) * t106;
t16 = -t102 * t29 - t106 * t72;
t21 = -t102 * t46 + t106 * t44;
t84 = t122 * qJD(6);
t85 = t123 * qJD(6);
t88 = Ifges(7,2) * t106 + t163;
t89 = Ifges(7,1) * t102 + t162;
t117 = t102 * t85 + t106 * t84 + t89 * t141 - t88 * t142;
t82 = t124 * qJD(6);
t116 = -mrSges(7,3) * t188 + t59 * t87 + t80 * t82 - t30;
t115 = -t172 + (-t102 * t17 - t106 * t16) * qJD(6);
t2 = t16 * qJD(6) + t102 * t69 + t106 * t9;
t114 = t106 * t2 + t115;
t113 = -t104 * t32 + t108 * t33 + (-t104 * t61 - t108 * t60) * qJD(4);
t112 = -t9 * mrSges(6,2) + t115 * mrSges(7,3) + t10 * t187 + t2 * t164 + t28 * t82;
t19 = t59 * mrSges(7,1) + t118 * mrSges(7,3);
t20 = -t59 * mrSges(7,2) - t119 * mrSges(7,3);
t5 = qJD(6) * t21 + t102 * t27 + t106 * t24;
t52 = -t80 * mrSges(7,2) - mrSges(7,3) * t158;
t53 = t80 * mrSges(7,1) - mrSges(7,3) * t151;
t111 = -t53 * t141 - t52 * t142 - t102 * t19 + t106 * t20 + m(7) * (t106 * t5 - t21 * t141 - t22 * t142 - t168);
t13 = -t118 * Ifges(7,4) - t119 * Ifges(7,2) + Ifges(7,6) * t59;
t14 = -t118 * Ifges(7,1) - t119 * Ifges(7,4) + Ifges(7,5) * t59;
t36 = Ifges(7,6) * t80 + t122 * t81;
t37 = Ifges(7,5) * t80 + t123 * t81;
t95 = Ifges(7,5) * t141;
t110 = t102 * t14 / 0.2e1 + t106 * t13 / 0.2e1 + t37 * t141 / 0.2e1 + t45 * t82 + t5 * t164 - t89 * t152 / 0.2e1 - Ifges(6,5) * t58 - t84 * t158 / 0.2e1 + t85 * t151 / 0.2e1 + t80 * (-Ifges(7,6) * t142 + t95) / 0.2e1 + (-t168 + (-t102 * t22 - t106 * t21) * qJD(6)) * mrSges(7,3) - t24 * mrSges(6,2) + (t121 / 0.2e1 - Ifges(6,6)) * t59 + t187 * t25 - t119 * t88 / 0.2e1 - (t81 * t89 + t36) * t142 / 0.2e1;
t94 = -pkin(4) * t107 - pkin(5);
t64 = mrSges(6,1) * t80 + mrSges(6,2) * t81;
t47 = t124 * t81;
t1 = [0.2e1 * m(7) * (t16 * t3 + t17 * t2 + t173) + 0.2e1 * m(6) * (t29 * t9 + t173 - t41) + 0.2e1 * m(5) * (t32 * t60 + t33 * t61 - t41) + 0.2e1 * m(4) * (t70 * t73 - t41); -t70 * mrSges(4,2) + t10 * t47 + t28 * t15 + t16 * t19 + t17 * t20 + t2 * t52 + t3 * t53 - (t30 + t83) * t72 + (-mrSges(3,1) * t105 - mrSges(3,2) * t109) * t146 + (-t108 * mrSges(5,1) + t104 * mrSges(5,2) - mrSges(4,1) + t64) * t69 + (t10 * t81 - t28 * t58 - t29 * t59 - t9 * t80) * mrSges(6,3) + t113 * mrSges(5,3) + m(7) * (t16 * t6 + t17 * t5 + t2 * t22 + t21 * t3 + t128) + m(6) * (-t72 * t139 + t24 * t29 + t46 * t9 + t69 * t86 + t128) + m(4) * (-t100 * t69 + t70 * t98) * pkin(2) + (t113 * t90 + t69 * t91) * m(5); t46 * t59 * t181 + t15 * t179 + 0.2e1 * t21 * t19 + 0.2e1 * t22 * t20 + t47 * t180 + 0.2e1 * t86 * t30 + 0.2e1 * t5 * t52 + 0.2e1 * t6 * t53 + 0.2e1 * t91 * t83 + (t21 * t6 + t22 * t5 + t171) * t183 + 0.2e1 * m(6) * (t86 * t139 + t46 * t24 + t171) - (mrSges(6,3) * t179 - t102 * t36 + t106 * t37) * t58 + (t24 * t181 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t59 - t131 * t58 + t166) * t80 + (mrSges(6,3) * t180 - 0.2e1 * Ifges(6,1) * t58 - t102 * t13 + t106 * t14 + (Ifges(7,5) * t106 + t131) * t59 + (-t102 * t37 - t106 * t36 - t80 * t121) * qJD(6)) * t81 + ((-Ifges(5,4) * t104 + pkin(4) * t64) * t140 + (0.2e1 * Ifges(5,4) * t108 + (Ifges(5,1) - Ifges(5,2)) * t140) * t108) * qJD(4); m(7) * (-(-t102 * t16 + t106 * t17) * t58 + t114 * t81 + t127) + m(6) * (-t29 * t58 + t81 * t9 + t127) + m(5) * (t104 * t33 + t108 * t32 + (-t104 * t60 + t108 * t61) * qJD(4)); t80 * t15 + t59 * t47 + m(7) * t126 + m(6) * (t24 * t81 - t46 * t58 + t126) + (m(7) * (-t21 * t145 - t22 * t58 + t5 * t81) - t53 * t145 - t58 * t52 + t81 * t20) * t106 + (m(7) * (-t22 * t145 + t21 * t58 - t6 * t81) + t58 * t53 - t81 * t19 - t52 * t145) * t102; 0.2e1 * m(6) * (-t81 * t58 + t167) + 0.2e1 * m(7) * (-t188 * t81 + t167); -t33 * mrSges(5,2) + t32 * mrSges(5,1) + ((-t10 * t107 + t103 * t9) * t178 + ((t17 * t143 - t16 * t144 + t156) * t177 + (t107 * t29 + t156) * t178) * qJD(5)) * t182 + t112 + (t10 * t94 + t2 * t150 + (-t141 * t16 - t142 * t17 - t172) * t93) * m(7); t110 + ((-mrSges(5,1) * t90 + Ifges(5,5)) * t108 + (mrSges(5,2) * t90 - Ifges(5,6)) * t104) * qJD(4) + t136 * t94 + t111 * t93 + (m(6) * (t103 * t24 - t107 * t25) + (-t103 * t59 + t107 * t58) * mrSges(6,3) + ((t81 * mrSges(6,3) + t47) * t103 + (-t80 * mrSges(6,3) - t102 * t53 + t106 * t52) * t107 + m(6) * (t107 * t46 + t155) + m(7) * (t22 * t143 - t21 * t144 + t155)) * qJD(5)) * pkin(4); m(7) * (-t188 * t93 + t59 * t94) - t83 + ((-t103 * t58 - t107 * t59) * t178 + ((t107 * t81 + t154) * t178 + (t185 * t81 + t154) * t177) * qJD(5)) * t182 + t116; 0.2e1 * t94 * t82 + (-0.2e1 * t149 - 0.2e1 * t157 + 0.2e1 * t153 + (t103 * t94 + t185 * t93) * t183 + 0.2e1 * t130) * t160 + t117; m(7) * (-pkin(5) * t10 + t114 * pkin(10)) + t112; -t136 * pkin(5) + t111 * pkin(10) + t110; m(7) * (-pkin(10) * t188 - t175) + t116; (t94 - pkin(5)) * t82 + (-t149 - t157 + t153 + m(7) * (-pkin(5) * t103 + pkin(10) * t185) + t130) * t160 + t117; -0.2e1 * pkin(5) * t82 + t117; mrSges(7,1) * t3 - mrSges(7,2) * t2; t6 * mrSges(7,1) - t5 * mrSges(7,2) - Ifges(7,5) * t138 - t119 * Ifges(7,6) + t166; -t15; t95 - t124 * t107 * t160 + (-mrSges(7,1) * t150 + (mrSges(7,2) * t93 - Ifges(7,6)) * t102) * qJD(6); t95 + (t87 * pkin(10) - t161) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
