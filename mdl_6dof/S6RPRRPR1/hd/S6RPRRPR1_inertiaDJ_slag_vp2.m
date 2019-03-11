% Calculate time derivative of joint inertia matrix for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:34
% EndTime: 2019-03-09 04:57:38
% DurationCPUTime: 2.41s
% Computational Cost: add. (5149->291), mult. (10833->443), div. (0->0), fcn. (10312->10), ass. (0->138)
t106 = sin(qJ(6));
t109 = cos(qJ(6));
t150 = t106 ^ 2 + t109 ^ 2;
t168 = mrSges(7,1) * t109;
t92 = mrSges(7,2) * t106 - t168;
t198 = -mrSges(6,1) + t92;
t111 = cos(qJ(3));
t136 = sin(pkin(10)) * pkin(1) + pkin(7);
t127 = pkin(8) + t136;
t121 = qJD(3) * t127;
t122 = t111 * t127;
t202 = -qJD(4) * t122 - t111 * t121;
t108 = sin(qJ(3));
t123 = t108 * t127;
t201 = qJD(4) * t123 + t108 * t121;
t148 = qJD(6) * t109;
t104 = sin(pkin(11));
t155 = cos(pkin(11));
t195 = qJD(3) + qJD(4);
t107 = sin(qJ(4));
t110 = cos(qJ(4));
t86 = -t107 * t108 + t110 * t111;
t68 = t195 * t86;
t87 = t107 * t111 + t110 * t108;
t69 = t195 * t87;
t49 = -t104 * t69 + t155 * t68;
t64 = t104 * t86 + t155 * t87;
t200 = t106 * t49 + t64 * t148;
t199 = t150 * mrSges(7,3);
t60 = -t107 * t122 - t110 * t123;
t114 = -t87 * qJ(5) + t60;
t61 = -t107 * t123 + t110 * t122;
t51 = qJ(5) * t86 + t61;
t26 = t104 * t114 + t155 * t51;
t63 = t104 * t87 - t155 * t86;
t143 = -cos(pkin(10)) * pkin(1) - pkin(2);
t91 = -pkin(3) * t111 + t143;
t70 = -t86 * pkin(4) + t91;
t27 = t63 * pkin(5) - t64 * pkin(9) + t70;
t14 = t106 * t27 + t109 * t26;
t151 = t14 * qJD(6);
t32 = t107 * t201 + t110 * t202;
t112 = -t68 * qJ(5) - t87 * qJD(5) + t32;
t31 = t107 * t202 - t110 * t201;
t19 = -t69 * qJ(5) + t86 * qJD(5) + t31;
t11 = t104 * t112 + t155 * t19;
t48 = t104 * t68 + t155 * t69;
t59 = qJD(3) * t108 * pkin(3) + pkin(4) * t69;
t17 = pkin(5) * t48 - pkin(9) * t49 + t59;
t3 = -t106 * t11 + t109 * t17 - t151;
t197 = -t3 - t151;
t152 = t104 * t107;
t164 = pkin(3) * qJD(4);
t77 = (t155 * t110 - t152) * t164;
t196 = t150 * t77;
t149 = qJD(6) * t106;
t145 = t64 * t149;
t158 = t109 * t49;
t125 = t145 - t158;
t194 = (-mrSges(5,1) * t107 - mrSges(5,2) * t110) * t164;
t193 = 2 * m(6);
t192 = 2 * m(7);
t10 = t104 * t19 - t155 * t112;
t191 = 0.2e1 * t10;
t190 = 0.2e1 * t59;
t189 = 0.2e1 * t91;
t188 = m(5) * pkin(3);
t187 = m(6) * pkin(4);
t185 = pkin(4) * t104;
t25 = t104 * t51 - t155 * t114;
t184 = t10 * t25;
t183 = t106 * t3;
t13 = -t106 * t26 + t109 * t27;
t154 = qJD(6) * t13;
t2 = t106 * t17 + t109 * t11 + t154;
t182 = t109 * t2;
t137 = t155 * t107;
t76 = (t104 * t110 + t137) * t164;
t181 = t25 * t76;
t180 = t48 * mrSges(6,3);
t179 = t48 * t63;
t178 = t49 * t64;
t177 = t63 * t76;
t176 = t64 * t77;
t175 = t68 * t87;
t174 = t69 * t86;
t99 = pkin(3) * t110 + pkin(4);
t80 = -pkin(3) * t152 + t155 * t99;
t74 = -pkin(5) - t80;
t133 = mrSges(7,1) * t106 + mrSges(7,2) * t109;
t88 = t133 * qJD(6);
t173 = t74 * t88;
t142 = t155 * pkin(4);
t98 = -t142 - pkin(5);
t172 = t98 * t88;
t171 = Ifges(7,5) * t158 + Ifges(7,3) * t48;
t170 = t48 * mrSges(6,1) + t49 * mrSges(6,2);
t81 = pkin(3) * t137 + t104 * t99;
t166 = Ifges(7,4) * t106;
t165 = Ifges(7,4) * t109;
t15 = mrSges(7,1) * t48 + t125 * mrSges(7,3);
t161 = t106 * t15;
t159 = t106 * t64;
t157 = t109 * t64;
t97 = pkin(9) + t185;
t156 = t109 * t97;
t75 = pkin(9) + t81;
t153 = qJD(6) * t75;
t147 = 0.2e1 * t111;
t146 = mrSges(7,3) * t154;
t141 = t69 * mrSges(5,1) + t68 * mrSges(5,2);
t139 = -Ifges(7,6) * t106 - (2 * Ifges(6,4));
t138 = t150 * t97;
t135 = t10 * t63 + t25 * t48;
t132 = Ifges(7,1) * t109 - t166;
t131 = -Ifges(7,2) * t106 + t165;
t130 = Ifges(7,5) * t106 + Ifges(7,6) * t109;
t129 = -t106 * t13 + t109 * t14;
t36 = -mrSges(7,2) * t63 - mrSges(7,3) * t159;
t37 = mrSges(7,1) * t63 - mrSges(7,3) * t157;
t128 = -t106 * t37 + t109 * t36;
t89 = t131 * qJD(6);
t90 = t132 * qJD(6);
t93 = Ifges(7,2) * t109 + t166;
t94 = Ifges(7,1) * t106 + t165;
t124 = t106 * t90 + t109 * t89 + t94 * t148 - t93 * t149;
t116 = t199 * t49 + t48 * t92 + t63 * t88 - t141 - t170;
t115 = -t183 + t182 + (-t106 * t14 - t109 * t13) * qJD(6);
t100 = Ifges(7,5) * t148;
t23 = t63 * Ifges(7,6) + t131 * t64;
t24 = t63 * Ifges(7,5) + t132 * t64;
t6 = -t125 * Ifges(7,4) - Ifges(7,2) * t200 + t48 * Ifges(7,6);
t7 = -t125 * Ifges(7,1) - Ifges(7,4) * t200 + t48 * Ifges(7,5);
t113 = -t31 * mrSges(5,2) - t11 * mrSges(6,2) + mrSges(7,3) * t182 + t25 * t88 + t24 * t148 / 0.2e1 + t94 * t158 / 0.2e1 + t32 * mrSges(5,1) + t106 * t7 / 0.2e1 + Ifges(6,5) * t49 + t109 * t6 / 0.2e1 - t89 * t159 / 0.2e1 + t90 * t157 / 0.2e1 + t63 * (-Ifges(7,6) * t149 + t100) / 0.2e1 - Ifges(5,6) * t69 + Ifges(5,5) * t68 + (t130 / 0.2e1 - Ifges(6,6)) * t48 - t200 * t93 / 0.2e1 - (t64 * t94 + t23) * t149 / 0.2e1 + t198 * t10;
t35 = t133 * t64;
t16 = -mrSges(7,2) * t48 - mrSges(7,3) * t200;
t12 = -mrSges(7,1) * t200 + mrSges(7,2) * t125;
t1 = [0.2e1 * Ifges(5,1) * t175 - 0.2e1 * Ifges(5,2) * t174 + t141 * t189 + 0.2e1 * t70 * t170 + 0.2e1 * t3 * t37 + t35 * t191 + 0.2e1 * t2 * t36 - 0.2e1 * t25 * t12 + 0.2e1 * t13 * t15 + 0.2e1 * t14 * t16 - 0.2e1 * t26 * t180 + (0.2e1 * t25 * mrSges(6,3) - t106 * t23 + t109 * t24) * t49 + (t13 * t3 + t14 * t2 + t184) * t192 + (t11 * t26 + t59 * t70 + t184) * t193 + 0.2e1 * m(5) * (t31 * t61 + t32 * t60) + ((t143 * mrSges(4,2) + Ifges(4,4) * t111) * t147 + (0.2e1 * pkin(3) * (-mrSges(5,1) * t86 + mrSges(5,2) * t87) + 0.2e1 * t143 * mrSges(4,1) + t188 * t189 - 0.2e1 * Ifges(4,4) * t108 + (Ifges(4,1) - Ifges(4,2)) * t147) * t108) * qJD(3) + (mrSges(6,1) * t190 - 0.2e1 * mrSges(6,3) * t11 + t139 * t49 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t48 + t171) * t63 + (mrSges(6,2) * t190 + mrSges(6,3) * t191 + 0.2e1 * Ifges(6,1) * t49 - t106 * t6 + t109 * t7 + (Ifges(7,5) * t109 + t139) * t48 + (-t106 * t24 - t109 * t23 - t63 * t130) * qJD(6)) * t64 + 0.2e1 * (t68 * t86 - t69 * t87) * Ifges(5,4) + 0.2e1 * (t31 * t86 - t32 * t87 - t60 * t68 - t61 * t69) * mrSges(5,3); -t63 * t12 + t48 * t35 + t128 * t49 + (-t161 + t109 * t16 + (-t106 * t36 - t109 * t37) * qJD(6)) * t64 + m(7) * (t115 * t64 + t129 * t49 + t135) + m(5) * (t31 * t87 + t32 * t86 - t60 * t69 + t61 * t68) + m(6) * (t11 * t64 + t26 * t49 + t135); 0.2e1 * m(7) * (t150 * t178 + t179) + 0.2e1 * m(5) * (-t174 + t175) + 0.2e1 * m(6) * (t178 + t179); t113 + ((-t136 * mrSges(4,1) + Ifges(4,5)) * t111 + (t136 * mrSges(4,2) - Ifges(4,6)) * t108) * qJD(3) + m(7) * (t10 * t74 + t181) + (m(5) * (t107 * t31 + t110 * t32 + (-t107 * t60 + t110 * t61) * qJD(4)) + (-t107 * t69 - t110 * t68 + (t107 * t87 + t110 * t86) * qJD(4)) * mrSges(5,3)) * pkin(3) + (-t36 * t153 + t197 * mrSges(7,3) + (-m(7) * t13 - t37) * t77 + (m(7) * t197 - t15) * t75) * t106 + (-t146 - t37 * t153 + t75 * t16 + t77 * t36 + m(7) * (-t13 * t153 + t14 * t77 + t2 * t75)) * t109 - t74 * t12 + t76 * t35 + m(6) * (-t10 * t80 + t11 * t81 + t26 * t77 + t181) + (-t81 * t48 - t80 * t49 - t77 * t63 + t76 * t64) * mrSges(6,3); (-mrSges(4,1) * t108 - mrSges(4,2) * t111) * qJD(3) + m(7) * (t48 * t74 + t177 + t150 * (t49 * t75 + t176)) + m(6) * (-t48 * t80 + t49 * t81 + t176 + t177) + (t107 * t68 - t110 * t69 + (-t107 * t86 + t110 * t87) * qJD(4)) * t188 + t116; 0.2e1 * t173 - 0.2e1 * t77 * mrSges(6,2) + (-t76 * t80 + t77 * t81) * t193 + (t196 * t75 + t74 * t76) * t192 + t124 + 0.2e1 * t198 * t76 + 0.2e1 * t194 + 0.2e1 * mrSges(7,3) * t196; (-t155 * t10 + t104 * t11) * t187 + t113 - t109 * t146 - t180 * t185 - t49 * mrSges(6,3) * t142 + t16 * t156 + (m(7) * t10 - t12) * t98 + (-t14 * t149 - t183) * mrSges(7,3) + (m(7) * t115 - t37 * t148 - t36 * t149 - t161) * t97; m(7) * (t49 * t138 + t98 * t48) + (t104 * t49 - t155 * t48) * t187 + t116; t124 + t172 + t173 + t194 + (m(7) * t98 - t155 * t187 + t198) * t76 + (m(7) * t138 + t104 * t187 - mrSges(6,2) + t199) * t77; t124 + 0.2e1 * t172; t106 * t16 + t109 * t15 + t128 * qJD(6) + m(7) * (t129 * qJD(6) + t106 * t2 + t109 * t3) + m(6) * t59 + t170; 0; 0; 0; 0; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t145 - Ifges(7,6) * t200 + t171; t12; t100 - t133 * t77 + (-t75 * t168 + (mrSges(7,2) * t75 - Ifges(7,6)) * t106) * qJD(6); t100 + (-mrSges(7,1) * t156 + (mrSges(7,2) * t97 - Ifges(7,6)) * t106) * qJD(6); -t88; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
