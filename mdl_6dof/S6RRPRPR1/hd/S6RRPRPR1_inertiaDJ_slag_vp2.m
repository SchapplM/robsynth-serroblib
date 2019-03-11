% Calculate time derivative of joint inertia matrix for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:07:56
% EndTime: 2019-03-09 10:08:01
% DurationCPUTime: 2.74s
% Computational Cost: add. (8803->358), mult. (19120->533), div. (0->0), fcn. (19648->10), ass. (0->154)
t159 = sin(pkin(11));
t161 = cos(pkin(11));
t218 = -mrSges(6,1) * t161 + mrSges(6,2) * t159 - mrSges(5,1);
t162 = cos(pkin(10));
t152 = pkin(2) * t162 + pkin(3);
t164 = sin(qJ(4));
t167 = cos(qJ(4));
t160 = sin(pkin(10));
t205 = pkin(2) * t160;
t127 = t152 * t167 - t164 * t205;
t165 = sin(qJ(2));
t201 = -qJ(3) - pkin(7);
t147 = t201 * t165;
t168 = cos(qJ(2));
t148 = t201 * t168;
t111 = t162 * t147 + t148 * t160;
t138 = t160 * t168 + t162 * t165;
t97 = -pkin(8) * t138 + t111;
t112 = t160 * t147 - t162 * t148;
t136 = -t160 * t165 + t162 * t168;
t98 = pkin(8) * t136 + t112;
t217 = -t164 * t98 + t167 * t97;
t128 = t164 * t152 + t167 * t205;
t163 = sin(qJ(6));
t166 = cos(qJ(6));
t171 = t159 * t163 - t161 * t166;
t134 = t171 * qJD(6);
t216 = 2 * m(5);
t215 = 2 * m(6);
t214 = 2 * m(7);
t213 = -2 * mrSges(5,3);
t133 = t136 * qJD(2);
t182 = qJD(2) * t201;
t129 = qJD(3) * t168 + t165 * t182;
t130 = -t165 * qJD(3) + t168 * t182;
t95 = -t129 * t160 + t130 * t162;
t170 = -pkin(8) * t133 + t95;
t66 = t164 * t97 + t167 * t98;
t132 = t138 * qJD(2);
t96 = t162 * t129 + t160 * t130;
t87 = -pkin(8) * t132 + t96;
t40 = qJD(4) * t66 + t164 * t87 - t167 * t170;
t212 = 0.2e1 * t40;
t211 = -0.2e1 * t217;
t139 = t159 * t166 + t161 * t163;
t135 = t139 * qJD(6);
t99 = t135 * mrSges(7,1) - t134 * mrSges(7,2);
t210 = 0.2e1 * t99;
t158 = t161 ^ 2;
t155 = qJD(2) * t165 * pkin(2);
t113 = pkin(3) * t132 + t155;
t209 = 0.2e1 * t113;
t154 = -pkin(2) * t168 - pkin(1);
t208 = 0.2e1 * t154;
t207 = m(4) * pkin(2);
t204 = pkin(5) * t161;
t39 = qJD(4) * t217 + t164 * t170 + t167 * t87;
t103 = t136 * t164 + t138 * t167;
t172 = t167 * t136 - t138 * t164;
t79 = qJD(4) * t172 - t132 * t164 + t133 * t167;
t80 = qJD(4) * t103 + t167 * t132 + t133 * t164;
t43 = pkin(4) * t80 - qJ(5) * t79 - qJD(5) * t103 + t113;
t15 = -t159 * t39 + t161 * t43;
t203 = t15 * mrSges(6,3);
t202 = t40 * t217;
t16 = t159 * t43 + t161 * t39;
t116 = -t136 * pkin(3) + t154;
t64 = -pkin(4) * t172 - t103 * qJ(5) + t116;
t47 = t159 * t64 + t161 * t66;
t193 = t161 * t79;
t195 = t159 * t79;
t50 = mrSges(6,1) * t195 + mrSges(6,2) * t193;
t200 = Ifges(6,4) * t159;
t199 = Ifges(6,4) * t161;
t198 = Ifges(6,2) * t159;
t119 = t127 * qJD(4);
t197 = t119 * mrSges(5,2);
t120 = t128 * qJD(4);
t196 = t120 * t217;
t194 = t16 * t161;
t191 = t103 * t159;
t190 = t103 * t161;
t188 = -Ifges(7,5) * t134 - Ifges(7,6) * t135;
t187 = t159 ^ 2 + t158;
t186 = 2 * mrSges(7,3);
t185 = 0.2e1 * t168;
t27 = -t103 * t135 - t171 * t79;
t28 = t103 * t134 - t139 * t79;
t184 = Ifges(7,5) * t27 + Ifges(7,6) * t28 + Ifges(7,3) * t80;
t183 = t80 * mrSges(5,1) + t79 * mrSges(5,2);
t13 = -t28 * mrSges(7,1) + t27 * mrSges(7,2);
t46 = -t159 * t66 + t161 * t64;
t100 = -Ifges(7,4) * t134 - Ifges(7,2) * t135;
t101 = -Ifges(7,1) * t134 - Ifges(7,4) * t135;
t105 = Ifges(7,4) * t139 - Ifges(7,2) * t171;
t106 = Ifges(7,1) * t139 - Ifges(7,4) * t171;
t181 = -t100 * t171 + t139 * t101 - t135 * t105 - t134 * t106;
t180 = t132 * mrSges(4,1) + t133 * mrSges(4,2);
t115 = qJD(5) + t119;
t179 = t187 * t115;
t178 = t187 * qJD(5);
t123 = -pkin(4) - t127;
t104 = mrSges(7,1) * t171 + mrSges(7,2) * t139;
t177 = (t104 + t218) * t120;
t176 = Ifges(6,5) * t161 - Ifges(6,6) * t159;
t175 = -t15 * t159 + t194;
t174 = -t159 * t46 + t161 * t47;
t23 = -pkin(5) * t172 - pkin(9) * t190 + t46;
t34 = -pkin(9) * t191 + t47;
t11 = -t163 * t34 + t166 * t23;
t12 = t163 * t23 + t166 * t34;
t173 = 0.2e1 * t187 * mrSges(6,3);
t122 = qJ(5) + t128;
t109 = (-pkin(9) - t122) * t159;
t156 = t161 * pkin(9);
t110 = t122 * t161 + t156;
t85 = t109 * t166 - t110 * t163;
t86 = t109 * t163 + t110 * t166;
t144 = (-pkin(9) - qJ(5)) * t159;
t146 = qJ(5) * t161 + t156;
t107 = t144 * t166 - t146 * t163;
t108 = t144 * t163 + t146 * t166;
t4 = pkin(5) * t80 - pkin(9) * t193 + t15;
t9 = -pkin(9) * t195 + t16;
t2 = qJD(6) * t11 + t163 * t4 + t166 * t9;
t20 = pkin(5) * t195 + t40;
t3 = -qJD(6) * t12 - t163 * t9 + t166 * t4;
t31 = t80 * Ifges(6,6) + (-t198 + t199) * t79;
t32 = t80 * Ifges(6,5) + (Ifges(6,1) * t161 - t200) * t79;
t70 = t139 * t103;
t71 = t171 * t103;
t44 = -Ifges(7,4) * t71 - Ifges(7,2) * t70 - Ifges(7,6) * t172;
t45 = -Ifges(7,1) * t71 - Ifges(7,4) * t70 - Ifges(7,5) * t172;
t53 = pkin(5) * t191 - t217;
t7 = Ifges(7,4) * t27 + Ifges(7,2) * t28 + Ifges(7,6) * t80;
t8 = Ifges(7,1) * t27 + Ifges(7,4) * t28 + Ifges(7,5) * t80;
t169 = -t172 * t188 / 0.2e1 + (Ifges(6,1) * t159 + t199) * t193 / 0.2e1 - (Ifges(6,2) * t161 + t200) * t195 / 0.2e1 + mrSges(6,3) * t194 - t39 * mrSges(5,2) + Ifges(5,5) * t79 - Ifges(5,6) * t80 + t53 * t99 - t70 * t100 / 0.2e1 - t71 * t101 / 0.2e1 + t20 * t104 + t28 * t105 / 0.2e1 + t27 * t106 / 0.2e1 - t134 * t45 / 0.2e1 - t135 * t44 / 0.2e1 - t171 * t7 / 0.2e1 + t139 * t8 / 0.2e1 + t159 * t32 / 0.2e1 + t161 * t31 / 0.2e1 + t218 * t40 + (Ifges(6,5) * t159 + Ifges(7,5) * t139 + Ifges(6,6) * t161 - Ifges(7,6) * t171) * t80 / 0.2e1 + (t11 * t134 - t12 * t135 - t139 * t3 - t171 * t2) * mrSges(7,3);
t153 = -pkin(4) - t204;
t114 = t123 - t204;
t93 = -qJD(5) * t139 - qJD(6) * t108;
t92 = -qJD(5) * t171 + qJD(6) * t107;
t82 = -mrSges(6,1) * t172 - mrSges(6,3) * t190;
t81 = mrSges(6,2) * t172 - mrSges(6,3) * t191;
t74 = (mrSges(6,1) * t159 + mrSges(6,2) * t161) * t103;
t59 = -mrSges(7,1) * t172 + mrSges(7,3) * t71;
t58 = mrSges(7,2) * t172 - mrSges(7,3) * t70;
t55 = -qJD(6) * t86 - t115 * t139;
t54 = qJD(6) * t85 - t115 * t171;
t52 = mrSges(6,1) * t80 - mrSges(6,3) * t193;
t51 = -mrSges(6,2) * t80 - mrSges(6,3) * t195;
t48 = mrSges(7,1) * t70 - mrSges(7,2) * t71;
t19 = -mrSges(7,2) * t80 + mrSges(7,3) * t28;
t18 = mrSges(7,1) * t80 - mrSges(7,3) * t27;
t1 = [0.2e1 * t116 * t183 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t168) * t185 + (t207 * t208 - 0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * pkin(2) * (-mrSges(4,1) * t136 + mrSges(4,2) * t138) - 0.2e1 * Ifges(3,4) * t165 + (-Ifges(3,2) + Ifges(3,1)) * t185) * t165) * qJD(2) - 0.2e1 * t136 * Ifges(4,2) * t132 + 0.2e1 * t138 * t133 * Ifges(4,1) + (-mrSges(5,1) * t209 - t39 * t213 - t184 + (-(2 * Ifges(5,2)) - (2 * Ifges(6,3)) - Ifges(7,3)) * t80 - 0.2e1 * (-Ifges(5,4) + t176) * t79) * t172 + (t15 * t46 + t16 * t47 - t202) * t215 + (t113 * t116 + t39 * t66 - t202) * t216 + t180 * t208 + t50 * t211 + t74 * t212 + (t11 * t3 + t12 * t2 + t20 * t53) * t214 + (mrSges(5,2) * t209 + mrSges(5,3) * t212 - t159 * t31 + t161 * t32) * t103 + 0.2e1 * (-t132 * t138 + t133 * t136) * Ifges(4,4) + 0.2e1 * (-t111 * t133 - t112 * t132 + t136 * t96 - t138 * t95) * mrSges(4,3) + (mrSges(5,3) * t211 + (Ifges(6,1) * t158 + (2 * Ifges(5,1)) + (t198 - 0.2e1 * t199) * t159) * t103) * t79 + 0.2e1 * m(4) * (t111 * t95 + t112 * t96) + 0.2e1 * t11 * t18 + 0.2e1 * t12 * t19 + t28 * t44 + t27 * t45 + 0.2e1 * t20 * t48 + 0.2e1 * t47 * t51 + 0.2e1 * t46 * t52 + 0.2e1 * t53 * t13 + 0.2e1 * t2 * t58 + 0.2e1 * t3 * t59 - t70 * t7 - t71 * t8 + 0.2e1 * t16 * t81 + 0.2e1 * t15 * t82 + ((-0.2e1 * Ifges(5,4) + t176) * t103 - Ifges(7,5) * t71 - Ifges(7,6) * t70 + t66 * t213) * t80; m(6) * (t115 * t174 + t122 * t175 + t123 * t40 - t196) + (t103 * t120 + t119 * t172 - t127 * t79 - t128 * t80) * mrSges(5,3) + (t48 + t74) * t120 + t169 + (Ifges(3,5) * t168 - Ifges(3,6) * t165 + (-mrSges(3,1) * t168 + mrSges(3,2) * t165) * pkin(7)) * qJD(2) + (t115 * t81 + t122 * t51) * t161 + (-t115 * t82 - t122 * t52 - t203) * t159 + (t160 * t96 + t162 * t95) * t207 + m(7) * (t11 * t55 + t114 * t20 + t12 * t54 + t120 * t53 + t2 * t86 + t3 * t85) + m(5) * (t119 * t66 - t127 * t40 + t128 * t39 - t196) + (-t132 * t160 - t133 * t162) * pkin(2) * mrSges(4,3) + t54 * t58 + t55 * t59 + t85 * t18 + t86 * t19 + t95 * mrSges(4,1) - t96 * mrSges(4,2) + t114 * t13 + t123 * t50 - Ifges(4,6) * t132 + Ifges(4,5) * t133; -0.2e1 * t197 + t114 * t210 + t115 * t173 + 0.2e1 * t177 + (t114 * t120 + t54 * t86 + t55 * t85) * t214 + (t120 * t123 + t122 * t179) * t215 + (t119 * t128 - t120 * t127) * t216 + (t85 * t134 - t86 * t135 - t55 * t139 - t171 * t54) * t186 + t181; m(4) * t155 - t134 * t58 - t135 * t59 - t171 * t18 + t139 * t19 + t159 * t51 + t161 * t52 + m(7) * (-t11 * t135 - t12 * t134 + t139 * t2 - t171 * t3) + m(6) * (t15 * t161 + t159 * t16) + m(5) * t113 + t180 + t183; m(7) * (-t134 * t86 - t135 * t85 + t139 * t54 - t171 * t55); (-t134 * t139 + t135 * t171) * t214; t169 + m(6) * (-pkin(4) * t40 + qJ(5) * t175 + qJD(5) * t174) + (qJ(5) * t51 + qJD(5) * t81) * t161 + (-qJ(5) * t52 - qJD(5) * t82 - t203) * t159 + m(7) * (t107 * t3 + t108 * t2 + t11 * t93 + t12 * t92 + t153 * t20) - pkin(4) * t50 + t92 * t58 + t93 * t59 + t107 * t18 + t108 * t19 + t153 * t13; -t197 + (t114 + t153) * t99 + t177 + m(7) * (t107 * t55 + t108 * t54 + t120 * t153 + t85 * t93 + t86 * t92) + m(6) * (-pkin(4) * t120 + qJ(5) * t179 + t122 * t178) + (t179 + t178) * mrSges(6,3) + ((-t55 - t93) * t139 - (t54 + t92) * t171 - (t108 + t86) * t135 - (-t107 - t85) * t134) * mrSges(7,3) + t181; m(7) * (-t107 * t135 - t108 * t134 + t139 * t92 - t171 * t93); t153 * t210 + (t107 * t93 + t108 * t92) * t214 + (qJ(5) * t187 * t215 + t173) * qJD(5) + (t107 * t134 - t108 * t135 - t93 * t139 - t171 * t92) * t186 + t181; m(6) * t40 + m(7) * t20 + t13 + t50; 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t120 + t99; 0; t99; 0; mrSges(7,1) * t3 - mrSges(7,2) * t2 + t184; mrSges(7,1) * t55 - mrSges(7,2) * t54 + t188; -t99; mrSges(7,1) * t93 - mrSges(7,2) * t92 + t188; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
