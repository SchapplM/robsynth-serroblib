% Calculate time derivative of joint inertia matrix for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
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
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP5_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP5_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP5_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP5_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:41:32
% EndTime: 2019-03-09 08:41:36
% DurationCPUTime: 2.41s
% Computational Cost: add. (2522->355), mult. (5548->498), div. (0->0), fcn. (4491->6), ass. (0->146)
t206 = Ifges(7,4) + Ifges(6,5);
t205 = Ifges(6,6) - Ifges(7,6);
t204 = -2 * mrSges(6,3) - 2 * mrSges(7,2);
t203 = m(6) + m(7);
t202 = mrSges(6,1) + mrSges(7,1);
t132 = sin(pkin(9));
t133 = cos(pkin(9));
t186 = sin(qJ(5));
t187 = cos(qJ(5));
t101 = t186 * t132 - t187 * t133;
t102 = t132 * t187 + t133 * t186;
t157 = qJD(5) * t186;
t158 = qJD(5) * t187;
t92 = -t132 * t158 - t133 * t157;
t93 = -t132 * t157 + t133 * t158;
t200 = -t101 * t92 + t102 * t93;
t199 = -t205 * t93 + t206 * t92;
t131 = t133 ^ 2;
t153 = (t132 ^ 2 + t131) * qJD(4);
t134 = -pkin(2) - qJ(4);
t185 = -pkin(8) + t134;
t109 = t185 * t132;
t151 = t185 * t187;
t155 = t186 * qJD(4);
t156 = t187 * qJD(4);
t28 = -t132 * t156 - t109 * t157 + (qJD(5) * t151 - t155) * t133;
t150 = t185 * t186;
t29 = t109 * t158 - t132 * t155 + (qJD(5) * t150 + t156) * t133;
t63 = t109 * t186 - t133 * t151;
t64 = t109 * t187 + t133 * t150;
t198 = t101 * t29 + t102 * t28 - t63 * t92 + t93 * t64;
t197 = m(7) * qJ(6) + mrSges(7,3);
t136 = cos(qJ(2));
t135 = sin(qJ(2));
t170 = t135 * qJ(3);
t100 = t134 * t136 - pkin(1) - t170;
t190 = pkin(3) + pkin(7);
t113 = t190 * t135;
t104 = t133 * t113;
t46 = t135 * pkin(4) + t104 + (pkin(8) * t136 - t100) * t132;
t171 = t133 * t136;
t61 = t133 * t100 + t132 * t113;
t52 = -pkin(8) * t171 + t61;
t182 = t186 * t46 + t187 * t52;
t173 = t132 * t135;
t167 = qJD(2) * t136;
t108 = t190 * t167;
t168 = qJD(2) * t135;
t152 = pkin(2) * t168 - t135 * qJD(3);
t72 = -qJD(4) * t136 + (-qJ(3) * t136 + qJ(4) * t135) * qJD(2) + t152;
t31 = t133 * t108 - t132 * t72;
t19 = (pkin(4) * t136 - pkin(8) * t173) * qJD(2) + t31;
t160 = t133 * t168;
t32 = t132 * t108 + t133 * t72;
t22 = pkin(8) * t160 + t32;
t4 = -qJD(5) * t182 - t186 * t22 + t187 * t19;
t196 = 2 * m(5);
t195 = 2 * m(6);
t194 = 2 * m(7);
t193 = -0.2e1 * pkin(1);
t146 = -pkin(2) * t136 - t170;
t110 = -pkin(1) + t146;
t191 = -0.2e1 * t110;
t189 = t133 / 0.2e1;
t128 = t136 * pkin(7);
t161 = t132 * t168;
t166 = qJD(5) * t136;
t49 = t102 * t166 + t160 * t187 - t161 * t186;
t24 = t49 * mrSges(7,2) + mrSges(7,3) * t167;
t27 = -mrSges(6,2) * t167 + t49 * mrSges(6,3);
t184 = t24 + t27;
t48 = t101 * t166 + t102 * t168;
t25 = mrSges(6,1) * t167 - t48 * mrSges(6,3);
t26 = -mrSges(7,1) * t167 + t48 * mrSges(7,2);
t183 = -t25 + t26;
t80 = t101 * t136;
t73 = mrSges(7,2) * t80 + mrSges(7,3) * t135;
t74 = -mrSges(6,2) * t135 + mrSges(6,3) * t80;
t181 = t73 + t74;
t81 = t102 * t136;
t75 = mrSges(6,1) * t135 + mrSges(6,3) * t81;
t76 = -mrSges(7,1) * t135 - mrSges(7,2) * t81;
t180 = t75 - t76;
t177 = Ifges(5,4) * t132;
t176 = Ifges(5,4) * t133;
t174 = t132 * Ifges(5,1);
t172 = t132 * t136;
t122 = t132 * pkin(4) + qJ(3);
t114 = t136 * pkin(3) + t128;
t165 = Ifges(6,5) * t48 + Ifges(6,6) * t49 + Ifges(6,3) * t167;
t164 = Ifges(7,4) * t48 + Ifges(7,2) * t167 - Ifges(7,6) * t49;
t91 = pkin(4) * t171 + t114;
t15 = -t49 * mrSges(6,1) + t48 * mrSges(6,2);
t55 = t93 * mrSges(6,1) + t92 * mrSges(6,2);
t14 = -t49 * mrSges(7,1) - t48 * mrSges(7,3);
t54 = t93 * mrSges(7,1) - t92 * mrSges(7,3);
t159 = t64 * t28 + t29 * t63;
t148 = t136 * mrSges(4,2) - t135 * mrSges(4,3);
t147 = -Ifges(5,5) * t132 - Ifges(5,6) * t133;
t145 = t132 * t32 + t133 * t31;
t83 = -mrSges(5,1) * t160 + mrSges(5,2) * t161;
t79 = (-pkin(4) * t133 - t190) * t168;
t12 = -t186 * t52 + t187 * t46;
t3 = -t157 * t52 + t46 * t158 + t186 * t19 + t187 * t22;
t140 = pkin(5) * t92 + qJ(6) * t93 + qJD(6) * t102;
t1 = qJ(6) * t167 + t135 * qJD(6) + t3;
t2 = -pkin(5) * t167 - t4;
t6 = qJ(6) * t135 + t182;
t7 = -t135 * pkin(5) - t12;
t139 = t1 * t102 + t101 * t2 + t6 * t93 - t7 * t92;
t138 = t101 * t4 - t102 * t3 - t12 * t92 - t182 * t93;
t112 = mrSges(5,1) * t132 + mrSges(5,2) * t133;
t107 = t190 * t168;
t106 = -t135 * mrSges(5,2) - mrSges(5,3) * t171;
t105 = t135 * mrSges(5,1) + mrSges(5,3) * t172;
t98 = (mrSges(5,3) * t133 * t135 - mrSges(5,2) * t136) * qJD(2);
t97 = (mrSges(5,1) * t136 - mrSges(5,3) * t173) * qJD(2);
t96 = (mrSges(5,1) * t133 - mrSges(5,2) * t132) * t136;
t78 = (t136 * Ifges(5,5) + (t174 + t176) * t135) * qJD(2);
t77 = (t136 * Ifges(5,6) + (t133 * Ifges(5,2) + t177) * t135) * qJD(2);
t70 = -Ifges(6,1) * t101 - Ifges(6,4) * t102;
t69 = -Ifges(7,1) * t101 + Ifges(7,5) * t102;
t68 = -Ifges(6,4) * t101 - Ifges(6,2) * t102;
t67 = -Ifges(7,5) * t101 + Ifges(7,3) * t102;
t66 = mrSges(6,1) * t102 - mrSges(6,2) * t101;
t65 = mrSges(7,1) * t102 + mrSges(7,3) * t101;
t60 = -t132 * t100 + t104;
t59 = Ifges(6,1) * t92 - Ifges(6,4) * t93;
t58 = Ifges(7,1) * t92 + Ifges(7,5) * t93;
t57 = Ifges(6,4) * t92 - Ifges(6,2) * t93;
t56 = Ifges(7,5) * t92 + Ifges(7,3) * t93;
t53 = pkin(5) * t102 + qJ(6) * t101 + t122;
t51 = -mrSges(6,1) * t80 - mrSges(6,2) * t81;
t50 = -mrSges(7,1) * t80 + mrSges(7,3) * t81;
t38 = -Ifges(6,1) * t81 + Ifges(6,4) * t80 + Ifges(6,5) * t135;
t37 = -Ifges(7,1) * t81 + Ifges(7,4) * t135 - Ifges(7,5) * t80;
t36 = -Ifges(6,4) * t81 + Ifges(6,2) * t80 + Ifges(6,6) * t135;
t35 = -Ifges(7,5) * t81 + Ifges(7,6) * t135 - Ifges(7,3) * t80;
t23 = pkin(5) * t93 - qJ(6) * t92 + qJD(6) * t101 + qJD(3);
t21 = -pkin(5) * t80 + qJ(6) * t81 + t91;
t11 = Ifges(6,1) * t48 + Ifges(6,4) * t49 + Ifges(6,5) * t167;
t10 = Ifges(7,1) * t48 + Ifges(7,4) * t167 - Ifges(7,5) * t49;
t9 = Ifges(6,4) * t48 + Ifges(6,2) * t49 + Ifges(6,6) * t167;
t8 = Ifges(7,5) * t48 + Ifges(7,6) * t167 - Ifges(7,3) * t49;
t5 = -pkin(5) * t49 - qJ(6) * t48 + qJD(6) * t81 + t79;
t13 = [(t36 - t35) * t49 - (t10 + t11) * t81 + (t9 - t8) * t80 + (t37 + t38) * t48 + (t12 * t4 + t182 * t3 + t79 * t91) * t195 + 0.2e1 * t182 * t27 + (-t107 * t114 + t31 * t60 + t32 * t61) * t196 + ((mrSges(3,2) * t193 + mrSges(4,3) * t191 + ((2 * Ifges(3,4)) + (2 * Ifges(4,6)) + t147) * t136 - t206 * t81 + t205 * t80) * t136 + (mrSges(3,1) * t193 + mrSges(4,2) * t191 + 0.2e1 * (-Ifges(3,4) - Ifges(4,6) - t147) * t135 + (-t131 * Ifges(5,2) + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) + Ifges(7,2) - (2 * Ifges(4,3)) + (2 * Ifges(5,3)) + Ifges(6,3) + (-t174 - 0.2e1 * t176) * t132) * t136) * t135) * qJD(2) + (t1 * t6 + t2 * t7 + t21 * t5) * t194 - t77 * t171 - t78 * t172 + t135 * t164 + t135 * t165 + 0.2e1 * (m(4) * t110 + t148) * (-qJ(3) * t167 + t152) + 0.2e1 * t21 * t14 + 0.2e1 * t6 * t24 + 0.2e1 * t12 * t25 + 0.2e1 * t7 * t26 + 0.2e1 * t5 * t50 + 0.2e1 * t1 * t73 + 0.2e1 * t3 * t74 + 0.2e1 * t4 * t75 + 0.2e1 * t2 * t76 + 0.2e1 * t79 * t51 + 0.2e1 * t91 * t15 + 0.2e1 * t60 * t97 + 0.2e1 * t61 * t98 + 0.2e1 * t31 * t105 + 0.2e1 * t32 * t106 - 0.2e1 * t107 * t96 + 0.2e1 * t114 * t83; (-qJD(4) * t105 - t31 * mrSges(5,3) + t134 * t97 + t78 / 0.2e1) * t133 + (t134 * t98 - qJD(4) * t106 - t32 * mrSges(5,3) - t77 / 0.2e1) * t132 + (t8 / 0.2e1 - t9 / 0.2e1) * t102 + (-t10 / 0.2e1 - t11 / 0.2e1) * t101 + m(6) * (-t12 * t29 + t122 * t79 + t182 * t28 + t3 * t64 - t4 * t63) + m(5) * (-qJ(3) * t107 + t145 * t134 + (-t132 * t61 - t133 * t60) * qJD(4)) + m(7) * (t1 * t64 + t2 * t63 + t21 * t23 + t28 * t6 + t29 * t7 + t5 * t53) + t199 * t135 / 0.2e1 - (t58 / 0.2e1 + t59 / 0.2e1) * t81 + (t57 / 0.2e1 - t56 / 0.2e1) * t80 + (-t67 / 0.2e1 + t68 / 0.2e1) * t49 + (t69 / 0.2e1 + t70 / 0.2e1) * t48 + ((Ifges(5,5) * t189 - Ifges(5,6) * t132 / 0.2e1 - Ifges(4,4) + Ifges(3,5) - pkin(2) * mrSges(4,1) + (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t102 + (-Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1) * t101) * t136 + (Ifges(4,5) - Ifges(3,6) + t132 * (Ifges(5,1) * t133 - t177) / 0.2e1 + (-Ifges(5,2) * t132 + t176) * t189 - qJ(3) * mrSges(4,1)) * t135 + (m(4) * t146 - t136 * mrSges(3,1) + t135 * mrSges(3,2) + t148) * pkin(7)) * qJD(2) - t180 * t29 + t181 * t28 + t183 * t63 + t184 * t64 + t138 * mrSges(6,3) - t139 * mrSges(7,2) + (m(4) * t128 + m(5) * t114 + m(6) * t91 + t136 * mrSges(4,1) + t51 + t96) * qJD(3) + t23 * t50 + t53 * t14 + t21 * t54 + t5 * t65 + t79 * t66 + qJ(3) * t83 + t91 * t55 + (t37 / 0.2e1 + t38 / 0.2e1) * t92 + (t35 / 0.2e1 - t36 / 0.2e1) * t93 - t107 * t112 + t122 * t15; 0.2e1 * t122 * t55 + 0.2e1 * t23 * t65 + 0.2e1 * t53 * t54 + (t67 - t68) * t93 + (t69 + t70) * t92 + (t56 - t57) * t102 + (-t58 - t59) * t101 + (qJ(3) * qJD(3) - t134 * t153) * t196 + (qJD(3) * t122 + t159) * t195 + (t23 * t53 + t159) * t194 + t198 * t204 + 0.2e1 * mrSges(5,3) * t153 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3) + t112 + t66) * qJD(3); t132 * t98 + t133 * t97 + t181 * t93 + t180 * t92 + t184 * t102 + t183 * t101 + (m(4) * pkin(7) + mrSges(4,1)) * t167 + m(7) * t139 - m(6) * t138 + m(5) * t145; -m(5) * t153 + t203 * t198 + t200 * t204; 0.2e1 * t203 * t200; -m(5) * t107 + m(6) * t79 + m(7) * t5 + t14 + t15 + t83; m(7) * t23 + (m(6) + m(5)) * qJD(3) + t54 + t55; 0; 0; -t2 * mrSges(7,1) - pkin(5) * t26 + t1 * mrSges(7,3) + qJD(6) * t73 + qJ(6) * t24 + m(7) * (-pkin(5) * t2 + qJ(6) * t1 + qJD(6) * t6) - t3 * mrSges(6,2) + t4 * mrSges(6,1) + t164 + t165; m(7) * qJD(6) * t64 - t140 * mrSges(7,2) + (-m(7) * pkin(5) - t202) * t29 + (-mrSges(6,2) + t197) * t28 + t199; m(7) * t140 + (-mrSges(6,2) + mrSges(7,3)) * t93 + t202 * t92; 0; 0.2e1 * t197 * qJD(6); m(7) * t2 + t26; m(7) * t29 + t92 * mrSges(7,2); -m(7) * t92; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t13(1) t13(2) t13(4) t13(7) t13(11) t13(16); t13(2) t13(3) t13(5) t13(8) t13(12) t13(17); t13(4) t13(5) t13(6) t13(9) t13(13) t13(18); t13(7) t13(8) t13(9) t13(10) t13(14) t13(19); t13(11) t13(12) t13(13) t13(14) t13(15) t13(20); t13(16) t13(17) t13(18) t13(19) t13(20) t13(21);];
Mq  = res;
