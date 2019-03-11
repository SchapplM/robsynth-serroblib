% Calculate time derivative of joint inertia matrix for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:12:07
% EndTime: 2019-03-09 10:12:13
% DurationCPUTime: 2.75s
% Computational Cost: add. (5267->289), mult. (11308->419), div. (0->0), fcn. (11383->8), ass. (0->135)
t204 = -mrSges(5,1) + mrSges(6,2);
t105 = sin(qJ(6));
t107 = cos(qJ(6));
t160 = qJD(6) * t107;
t106 = sin(qJ(4));
t166 = sin(pkin(10));
t167 = cos(pkin(10));
t186 = sin(qJ(2));
t188 = cos(qJ(2));
t126 = t166 * t186 - t167 * t188;
t120 = qJD(2) * t126;
t85 = t166 * t188 + t167 * t186;
t121 = qJD(2) * t85;
t122 = t106 * t126;
t187 = cos(qJ(4));
t142 = qJD(4) * t187;
t47 = -qJD(4) * t122 - t106 * t120 + t121 * t187 + t142 * t85;
t117 = t187 * t126;
t64 = t106 * t85 + t117;
t131 = t105 * t47 + t64 * t160;
t176 = Ifges(7,4) * t107;
t96 = -Ifges(7,2) * t105 + t176;
t177 = Ifges(7,4) * t105;
t97 = Ifges(7,1) * t107 - t177;
t205 = t105 * t97 + t107 * t96;
t161 = qJD(6) * t105;
t199 = t107 * t47 - t64 * t161;
t203 = mrSges(5,2) - mrSges(6,3);
t202 = Ifges(3,1) - Ifges(3,2);
t201 = -qJ(3) - pkin(7);
t94 = mrSges(7,1) * t105 + mrSges(7,2) * t107;
t200 = t94 + mrSges(6,3);
t162 = qJD(4) * t106;
t46 = t106 * t121 + t162 * t85 - (-qJD(2) - qJD(4)) * t117;
t65 = t187 * t85 - t122;
t143 = qJD(2) * t186;
t102 = pkin(2) * t143;
t68 = pkin(3) * t121 + t102;
t19 = t47 * pkin(4) + t46 * qJ(5) - t65 * qJD(5) + t68;
t10 = t47 * pkin(9) + t19;
t100 = -pkin(2) * t188 - pkin(1);
t70 = pkin(3) * t126 + t100;
t111 = -t65 * qJ(5) + t70;
t189 = pkin(4) + pkin(9);
t23 = t189 * t64 + t111;
t93 = t201 * t186;
t95 = t201 * t188;
t66 = t166 * t95 + t167 * t93;
t129 = -t85 * pkin(8) + t66;
t59 = t187 * t129;
t67 = t166 * t93 - t167 * t95;
t62 = -pkin(8) * t126 + t67;
t31 = t106 * t62 - t59;
t24 = pkin(5) * t65 + t31;
t12 = -t105 * t23 + t107 * t24;
t115 = qJD(2) * t95 - qJD(3) * t186;
t116 = qJD(2) * t93 + qJD(3) * t188;
t60 = t115 * t167 - t116 * t166;
t109 = pkin(8) * t120 + t60;
t61 = t115 * t166 + t116 * t167;
t110 = -pkin(8) * t121 + t61;
t127 = t106 * t129;
t18 = qJD(4) * t127 + t106 * t110 - t187 * t109 + t142 * t62;
t9 = -t46 * pkin(5) + t18;
t1 = qJD(6) * t12 + t10 * t107 + t105 * t9;
t184 = t1 * t105;
t13 = t105 * t24 + t107 * t23;
t165 = qJD(6) * t13;
t2 = -t10 * t105 + t107 * t9 - t165;
t198 = t2 * t107 + t184;
t20 = mrSges(7,2) * t46 + mrSges(7,3) * t199;
t21 = -mrSges(7,1) * t46 - mrSges(7,3) * t131;
t112 = m(7) * ((-t105 * t12 + t107 * t13) * qJD(6) + t198) + t107 * t21 + t105 * t20;
t173 = t105 * t64;
t48 = mrSges(7,1) * t65 - mrSges(7,3) * t173;
t169 = t107 * t64;
t49 = -mrSges(7,2) * t65 + mrSges(7,3) * t169;
t197 = t49 * t160 - t48 * t161 + t112;
t195 = 2 * m(5);
t194 = 2 * m(6);
t193 = 0.2e1 * m(7);
t89 = -mrSges(7,1) * t160 + mrSges(7,2) * t161;
t191 = -0.2e1 * t89;
t185 = Ifges(4,4) * t85;
t148 = pkin(2) * t166;
t147 = t167 * pkin(2);
t99 = t147 + pkin(3);
t82 = t106 * t99 + t148 * t187;
t72 = t82 * qJD(4);
t182 = t31 * t72;
t181 = t46 * mrSges(6,1);
t180 = t65 * t72;
t98 = t106 * t148;
t71 = -qJD(4) * t98 + t99 * t142;
t69 = -qJD(5) - t71;
t77 = qJ(5) + t82;
t179 = t69 * t77;
t178 = t71 * mrSges(5,2);
t164 = t105 ^ 2 + t107 ^ 2;
t149 = t169 / 0.2e1;
t146 = t47 * mrSges(5,1) - t46 * mrSges(5,2);
t145 = -t47 * mrSges(6,2) + t46 * mrSges(6,3);
t144 = qJD(2) * t188;
t140 = t164 * mrSges(7,3);
t139 = t164 * t72;
t81 = t187 * t99 - t98;
t17 = -qJD(4) * t59 - t106 * t109 - t187 * t110 + t162 * t62;
t32 = t187 * t62 + t127;
t138 = -t17 * t32 + t18 * t31;
t137 = mrSges(7,1) * t107 - mrSges(7,2) * t105;
t136 = Ifges(7,1) * t105 + t176;
t135 = Ifges(7,2) * t107 + t177;
t134 = Ifges(7,5) * t105 + Ifges(7,6) * t107;
t133 = t105 * t13 + t107 * t12;
t132 = -qJ(5) * t69 + qJD(5) * t77;
t78 = -pkin(4) - t81;
t128 = t131 * Ifges(7,5) + t199 * Ifges(7,6) - Ifges(7,3) * t46;
t90 = t135 * qJD(6);
t91 = t136 * qJD(6);
t125 = -t205 * qJD(6) + t105 * t90 - t107 * t91;
t119 = mrSges(4,3) * t121;
t118 = mrSges(4,3) * t120;
t114 = mrSges(4,1) * t121 - mrSges(4,2) * t120;
t25 = -t64 * pkin(5) + t32;
t28 = Ifges(7,6) * t65 + t135 * t64;
t29 = Ifges(7,5) * t65 + t136 * t64;
t5 = Ifges(7,4) * t131 + Ifges(7,2) * t199 - t46 * Ifges(7,6);
t6 = Ifges(7,1) * t131 + Ifges(7,4) * t199 - t46 * Ifges(7,5);
t8 = -pkin(5) * t47 - t17;
t113 = t12 * mrSges(7,3) * t161 - t28 * t160 / 0.2e1 - t105 * t5 / 0.2e1 + t107 * t6 / 0.2e1 - t91 * t173 / 0.2e1 - t90 * t149 + t8 * t94 - t25 * t89 + (-Ifges(7,5) * t107 / 0.2e1 + Ifges(7,6) * t105 / 0.2e1 + Ifges(6,4) - Ifges(5,5)) * t46 + t204 * t18 - (t64 * t96 + t29) * t161 / 0.2e1 + (t97 * t149 - t65 * t134 / 0.2e1) * qJD(6) + (Ifges(6,5) - Ifges(5,6) + t205 / 0.2e1) * t47;
t76 = -pkin(9) + t78;
t38 = t137 * t64;
t30 = t64 * pkin(4) + t111;
t14 = -mrSges(7,1) * t199 + mrSges(7,2) * t131;
t3 = [0.2e1 * ((-Ifges(5,4) - Ifges(6,6)) * t65 + (Ifges(6,3) + Ifges(5,2)) * t64) * t47 + 0.2e1 * t66 * t118 - 0.2e1 * t67 * t119 + 0.2e1 * t19 * (-mrSges(6,2) * t64 - mrSges(6,3) * t65) + 0.2e1 * t68 * (mrSges(5,1) * t64 + mrSges(5,2) * t65) + 0.2e1 * t2 * t48 + 0.2e1 * t1 * t49 - 0.2e1 * t8 * t38 + 0.2e1 * t25 * t14 + 0.2e1 * t13 * t20 + 0.2e1 * t12 * t21 + (-0.2e1 * Ifges(3,4) * t186 + t188 * t202) * t143 + (0.2e1 * Ifges(3,4) * t188 + t186 * t202) * t144 + t199 * t28 + t131 * t29 - (-0.2e1 * Ifges(4,4) * t126 + (Ifges(4,1) - Ifges(4,2)) * t85) * t120 - 0.2e1 * (t126 * t61 + t60 * t85) * mrSges(4,3) + ((-(2 * Ifges(5,1)) - (2 * Ifges(6,2)) - Ifges(7,3)) * t65 + (0.2e1 * Ifges(5,4) + 0.2e1 * Ifges(6,6) - t134) * t64) * t46 + 0.2e1 * (mrSges(5,3) + mrSges(6,1)) * (t17 * t64 + t18 * t65 - t31 * t46 - t32 * t47) + (t85 * (-Ifges(4,1) * t126 - t185) - 0.2e1 * pkin(1) * (mrSges(3,1) * t186 + mrSges(3,2) * t188)) * qJD(2) + 0.2e1 * t100 * t114 + 0.2e1 * (mrSges(4,1) * t126 + t85 * mrSges(4,2)) * t102 + 0.2e1 * m(4) * (t100 * t102 + t66 * t60 + t67 * t61) + t5 * t169 + t6 * t173 + (t1 * t13 + t12 * t2 + t25 * t8) * t193 + (t19 * t30 + t138) * t194 + (t68 * t70 + t138) * t195 + 0.2e1 * t30 * t145 + 0.2e1 * t70 * t146 - (-Ifges(4,2) * t126 + t185) * t121 + t65 * t128; t113 + Ifges(3,5) * t144 - t119 * t148 - Ifges(4,5) * t120 - Ifges(4,6) * t121 + t77 * t14 + t69 * t38 + t60 * mrSges(4,1) - t61 * mrSges(4,2) + m(4) * (t166 * t61 + t167 * t60) * pkin(2) + m(7) * (-t25 * t69 + t77 * t8) - Ifges(3,6) * t143 + t118 * t147 - t78 * t181 + m(6) * (t18 * t78 - t32 * t69 + t182) + m(5) * (-t18 * t81 + t32 * t71 + t182) + (m(7) * t133 + t105 * t49 + t107 * t48) * t72 + (-m(5) * t82 - m(6) * t77 + t203) * t17 + (-mrSges(3,1) * t144 + mrSges(3,2) * t143) * pkin(7) + (-t13 * t160 - t198) * mrSges(7,3) + (-t47 * t77 + t64 * t69 + t180) * mrSges(6,1) + t197 * t76 + (t46 * t81 - t47 * t82 - t64 * t71 + t180) * mrSges(5,3); -0.2e1 * t178 + t77 * t191 + (t139 * t76 - t179) * t193 - t179 * t194 + t71 * t82 * t195 + t125 - 0.2e1 * t200 * t69 + (t194 * t78 - t195 * t81 - 0.2e1 * t140 + 0.2e1 * t204) * t72; -t49 * t161 + t107 * t20 + m(7) * (-qJD(6) * t133 + t1 * t107 - t105 * t2) - t48 * t160 - t105 * t21 + m(6) * t19 + m(5) * t68 + m(4) * t102 + t114 + t145 + t146; 0; 0; t113 - qJD(5) * t38 + qJ(5) * t14 + (-t184 + (-t2 - t165) * t107) * mrSges(7,3) + t203 * t17 + m(7) * (qJ(5) * t8 + qJD(5) * t25) + m(6) * (-pkin(4) * t18 - qJ(5) * t17 + qJD(5) * t32) + (pkin(4) * t46 - qJ(5) * t47 - qJD(5) * t64) * mrSges(6,1) - t197 * t189; -t178 + (-t77 - qJ(5)) * t89 + (-t140 + t204) * t72 + m(7) * (-t139 * t189 + t132) + m(6) * (-pkin(4) * t72 + t132) + t125 + t200 * (-t69 + qJD(5)); 0; qJ(5) * t191 + 0.2e1 * ((m(6) + m(7)) * qJ(5) + t200) * qJD(5) + t125; -t181 + (-t105 * t48 + t107 * t49) * qJD(6) + m(6) * t18 + t112; 0.2e1 * (m(7) * t164 / 0.2e1 + m(6) / 0.2e1) * t72; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t128; t137 * t72 + ((-mrSges(7,2) * t76 - Ifges(7,6)) * t107 + (-mrSges(7,1) * t76 - Ifges(7,5)) * t105) * qJD(6); t89; ((mrSges(7,2) * t189 - Ifges(7,6)) * t107 + (mrSges(7,1) * t189 - Ifges(7,5)) * t105) * qJD(6); -t94 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
