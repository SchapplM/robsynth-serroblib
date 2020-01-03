% Calculate time derivative of joint inertia matrix for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR8_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR8_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:25
% EndTime: 2019-12-31 22:24:34
% DurationCPUTime: 3.03s
% Computational Cost: add. (5180->359), mult. (11638->535), div. (0->0), fcn. (10619->8), ass. (0->166)
t159 = sin(qJ(3));
t160 = sin(qJ(2));
t163 = cos(qJ(3));
t164 = cos(qJ(2));
t131 = t159 * t164 + t163 * t160;
t158 = sin(qJ(4));
t197 = qJD(4) * t158;
t184 = t131 * t197;
t129 = t159 * t160 - t163 * t164;
t236 = qJD(2) + qJD(3);
t102 = t236 * t129;
t162 = cos(qJ(4));
t205 = t102 * t162;
t168 = t184 + t205;
t196 = qJD(4) * t162;
t206 = t102 * t158;
t169 = t131 * t196 - t206;
t35 = mrSges(5,1) * t169 - mrSges(5,2) * t168;
t228 = -pkin(7) - pkin(6);
t145 = t228 * t160;
t147 = t228 * t164;
t116 = t145 * t159 - t147 * t163;
t186 = qJD(2) * t228;
t139 = t160 * t186;
t178 = t164 * t186;
t69 = qJD(3) * t116 + t139 * t159 - t163 * t178;
t244 = m(5) * t69 + t35;
t237 = (t158 ^ 2 + t162 ^ 2) * t163;
t157 = sin(qJ(5));
t161 = cos(qJ(5));
t172 = t157 * t158 - t161 * t162;
t235 = qJD(4) + qJD(5);
t100 = t235 * t172;
t130 = t157 * t162 + t158 * t161;
t101 = t235 * t130;
t218 = -Ifges(6,5) * t100 - Ifges(6,6) * t101;
t243 = Ifges(5,5) * t196 + t218;
t103 = t236 * t131;
t53 = pkin(2) * qJD(2) * t160 + pkin(3) * t103 + pkin(8) * t102;
t238 = t163 * t145 + t147 * t159;
t68 = t238 * qJD(3) + t163 * t139 + t159 * t178;
t152 = -pkin(2) * t164 - pkin(1);
t91 = t129 * pkin(3) - t131 * pkin(8) + t152;
t14 = -t116 * t197 + t158 * t53 + t162 * t68 + t91 * t196;
t180 = -t158 * t68 + t162 * t53;
t109 = t162 * t116;
t56 = t158 * t91 + t109;
t15 = -t56 * qJD(4) + t180;
t242 = t14 * t162 - t15 * t158;
t241 = -Ifges(5,5) * t205 + Ifges(5,3) * t103;
t149 = pkin(2) * t159 + pkin(8);
t219 = -pkin(9) - t149;
t179 = qJD(4) * t219;
t214 = pkin(2) * qJD(3);
t189 = t163 * t214;
t117 = t158 * t179 + t162 * t189;
t118 = -t158 * t189 + t162 * t179;
t126 = t219 * t158;
t154 = t162 * pkin(9);
t127 = t149 * t162 + t154;
t86 = t126 * t161 - t127 * t157;
t40 = qJD(5) * t86 + t117 * t161 + t118 * t157;
t87 = t126 * t157 + t127 * t161;
t41 = -qJD(5) * t87 - t117 * t157 + t118 * t161;
t240 = t41 * mrSges(6,1) - t40 * mrSges(6,2);
t227 = -pkin(9) - pkin(8);
t144 = t227 * t158;
t146 = pkin(8) * t162 + t154;
t113 = t144 * t161 - t146 * t157;
t185 = qJD(4) * t227;
t136 = t158 * t185;
t137 = t162 * t185;
t66 = qJD(5) * t113 + t136 * t161 + t137 * t157;
t115 = t144 * t157 + t146 * t161;
t67 = -qJD(5) * t115 - t136 * t157 + t137 * t161;
t239 = t67 * mrSges(6,1) - t66 * mrSges(6,2);
t83 = t172 * t131;
t234 = 0.2e1 * m(5);
t233 = 2 * m(6);
t57 = mrSges(6,1) * t101 - mrSges(6,2) * t100;
t232 = 0.2e1 * t57;
t231 = 0.2e1 * t69;
t106 = mrSges(6,1) * t172 + mrSges(6,2) * t130;
t230 = 0.2e1 * t106;
t229 = 0.2e1 * t152;
t222 = pkin(2) * t163;
t217 = Ifges(5,4) * t158;
t216 = Ifges(5,4) * t162;
t215 = Ifges(5,6) * t158;
t213 = pkin(4) * qJD(5);
t212 = t101 * mrSges(6,3);
t211 = t238 * t69;
t208 = t159 * mrSges(4,1);
t207 = t163 * mrSges(4,2);
t203 = t131 * t158;
t202 = t131 * t162;
t141 = -mrSges(5,1) * t162 + mrSges(5,2) * t158;
t198 = t159 * t141;
t195 = qJD(5) * t157;
t194 = qJD(5) * t161;
t193 = 0.2e1 * mrSges(6,3);
t192 = 0.2e1 * t164;
t24 = -t101 * t131 + t102 * t172;
t25 = t130 * t102 + t235 * t83;
t191 = Ifges(6,5) * t24 + Ifges(6,6) * t25 + Ifges(6,3) * t103;
t190 = mrSges(6,3) * t213;
t188 = pkin(4) * t197;
t187 = t161 * t100 * mrSges(6,3);
t151 = -pkin(4) * t162 - pkin(3);
t182 = -t197 / 0.2e1;
t181 = -(2 * Ifges(4,4)) - t215;
t55 = -t116 * t158 + t162 * t91;
t177 = mrSges(5,3) * t237;
t176 = mrSges(5,1) * t158 + mrSges(5,2) * t162;
t175 = Ifges(5,1) * t162 - t217;
t174 = -Ifges(5,2) * t158 + t216;
t173 = Ifges(5,5) * t158 + Ifges(5,6) * t162;
t38 = pkin(4) * t129 - pkin(9) * t202 + t55;
t46 = -pkin(9) * t203 + t56;
t16 = -t157 * t46 + t161 * t38;
t17 = t157 * t38 + t161 * t46;
t10 = pkin(9) * t205 + pkin(4) * t103 + (-t109 + (pkin(9) * t131 - t91) * t158) * qJD(4) + t180;
t11 = -pkin(9) * t169 + t14;
t3 = qJD(5) * t16 + t10 * t157 + t11 * t161;
t4 = -qJD(5) * t17 + t10 * t161 - t11 * t157;
t171 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + t191;
t170 = -t161 * t172 * t190 + (-pkin(4) * t212 + t130 * t190) * t157 + t243;
t107 = Ifges(6,4) * t130 - Ifges(6,2) * t172;
t108 = Ifges(6,1) * t130 - Ifges(6,4) * t172;
t134 = t174 * qJD(4);
t135 = t175 * qJD(4);
t142 = Ifges(5,2) * t162 + t217;
t143 = Ifges(5,1) * t158 + t216;
t58 = -Ifges(6,4) * t100 - Ifges(6,2) * t101;
t59 = -Ifges(6,1) * t100 - Ifges(6,4) * t101;
t167 = -t100 * t108 - t101 * t107 + t130 * t59 + t162 * t134 + t158 * t135 - t142 * t197 + t143 * t196 - t172 * t58;
t44 = mrSges(5,1) * t103 + mrSges(5,3) * t168;
t45 = -mrSges(5,2) * t103 - mrSges(5,3) * t169;
t92 = -mrSges(5,2) * t129 - mrSges(5,3) * t203;
t93 = mrSges(5,1) * t129 - mrSges(5,3) * t202;
t166 = -t93 * t196 - t92 * t197 - t158 * t44 + t162 * t45 + m(5) * (-t196 * t55 - t197 * t56 + t242);
t133 = t176 * qJD(4);
t29 = -Ifges(5,4) * t168 - Ifges(5,2) * t169 + Ifges(5,6) * t103;
t30 = -Ifges(5,1) * t168 - Ifges(5,4) * t169 + Ifges(5,5) * t103;
t33 = pkin(4) * t169 + t69;
t82 = t130 * t131;
t42 = -Ifges(6,4) * t83 - Ifges(6,2) * t82 + Ifges(6,6) * t129;
t43 = -Ifges(6,1) * t83 - Ifges(6,4) * t82 + Ifges(6,5) * t129;
t7 = Ifges(6,4) * t24 + Ifges(6,2) * t25 + Ifges(6,6) * t103;
t74 = Ifges(5,6) * t129 + t131 * t174;
t75 = Ifges(5,5) * t129 + t131 * t175;
t78 = pkin(4) * t203 - t238;
t8 = Ifges(6,1) * t24 + Ifges(6,4) * t25 + Ifges(6,5) * t103;
t165 = (t131 * t182 - t205 / 0.2e1) * t143 + (Ifges(6,5) * t130 - Ifges(6,6) * t172 + t173) * t103 / 0.2e1 - t172 * t7 / 0.2e1 + (t16 * t100 - t4 * t130 - t172 * t3) * mrSges(6,3) - t238 * t133 + t135 * t202 / 0.2e1 - t134 * t203 / 0.2e1 + t75 * t196 / 0.2e1 - t17 * t212 + (t141 - mrSges(4,1)) * t69 + t162 * t29 / 0.2e1 + t158 * t30 / 0.2e1 + t130 * t8 / 0.2e1 + t33 * t106 + t25 * t107 / 0.2e1 + t24 * t108 / 0.2e1 - t100 * t43 / 0.2e1 - t101 * t42 / 0.2e1 - Ifges(4,5) * t102 - Ifges(4,6) * t103 + t78 * t57 - t82 * t58 / 0.2e1 - t83 * t59 / 0.2e1 - t68 * mrSges(4,2) + t74 * t182 - t169 * t142 / 0.2e1 + (-Ifges(5,6) * t197 + t243) * t129 / 0.2e1 + ((-t158 * t56 - t162 * t55) * qJD(4) + t242) * mrSges(5,3);
t150 = -pkin(3) - t222;
t140 = t151 - t222;
t138 = t159 * t214 + t188;
t125 = (-mrSges(6,1) * t157 - mrSges(6,2) * t161) * t213;
t89 = t176 * t131;
t73 = mrSges(6,1) * t129 + mrSges(6,3) * t83;
t72 = -mrSges(6,2) * t129 - mrSges(6,3) * t82;
t47 = mrSges(6,1) * t82 - mrSges(6,2) * t83;
t19 = -mrSges(6,2) * t103 + mrSges(6,3) * t25;
t18 = mrSges(6,1) * t103 - mrSges(6,3) * t24;
t9 = -mrSges(6,1) * t25 + mrSges(6,2) * t24;
t1 = [0.2e1 * m(4) * (t116 * t68 - t211) + (t14 * t56 + t15 * t55 - t211) * t234 + (-0.2e1 * t68 * mrSges(4,3) - t181 * t102 + ((2 * Ifges(4,2)) + Ifges(5,3) + Ifges(6,3)) * t103 + t191 + t241) * t129 + (mrSges(4,3) * t231 - 0.2e1 * Ifges(4,1) * t102 - t158 * t29 + t162 * t30 + (Ifges(5,5) * t162 + t181) * t103 + (-t129 * t173 - t158 * t75 - t162 * t74) * qJD(4)) * t131 + (mrSges(4,1) * t103 - mrSges(4,2) * t102) * t229 - 0.2e1 * t238 * t35 + 0.2e1 * (t102 * t238 - t103 * t116) * mrSges(4,3) + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t164) * t192 + (0.2e1 * pkin(2) * (mrSges(4,1) * t129 + mrSges(4,2) * t131) - 0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t229 - 0.2e1 * Ifges(3,4) * t160 + (Ifges(3,1) - Ifges(3,2)) * t192) * t160) * qJD(2) + t74 * t206 + t89 * t231 + (t16 * t4 + t17 * t3 + t33 * t78) * t233 - t75 * t205 + t103 * (-Ifges(6,5) * t83 - Ifges(6,6) * t82) + 0.2e1 * t14 * t92 + 0.2e1 * t15 * t93 + 0.2e1 * t78 * t9 - t82 * t7 - t83 * t8 + 0.2e1 * t3 * t72 + 0.2e1 * t4 * t73 + 0.2e1 * t55 * t44 + 0.2e1 * t56 * t45 + 0.2e1 * t33 * t47 + t25 * t42 + t24 * t43 + 0.2e1 * t16 * t18 + 0.2e1 * t17 * t19; t165 + t166 * t149 + (m(4) * (t159 * t68 - t163 * t69) + (t163 * t102 - t159 * t103) * mrSges(4,3) + ((-t129 * mrSges(4,3) - t158 * t93 + t162 * t92 + m(4) * t116 + m(5) * (-t158 * t55 + t162 * t56)) * t163 + (t131 * mrSges(4,3) + t89 - (m(4) + m(5)) * t238) * t159) * qJD(3)) * pkin(2) + t138 * t47 + t140 * t9 + t86 * t18 + t87 * t19 + t40 * t72 + t41 * t73 + (Ifges(3,5) * t164 - Ifges(3,6) * t160 + (-mrSges(3,1) * t164 + mrSges(3,2) * t160) * pkin(6)) * qJD(2) + m(6) * (t138 * t78 + t140 * t33 + t16 * t41 + t17 * t40 + t3 * t87 + t4 * t86) + t244 * t150; t138 * t230 + t140 * t232 + (t138 * t140 + t40 * t87 + t41 * t86) * t233 + 0.2e1 * t150 * t133 + (t86 * t100 - t87 * t101 - t41 * t130 - t172 * t40) * t193 + (-0.2e1 * t207 - 0.2e1 * t208 + 0.2e1 * t198 + (t237 * t149 + t150 * t159) * t234 + 0.2e1 * t177) * t214 + t167; t165 + t166 * pkin(8) + m(6) * (t113 * t4 + t115 * t3 + t151 * t33 + t16 * t67 + t17 * t66 + t188 * t78) + t47 * t188 + t151 * t9 + t113 * t18 + t115 * t19 + t66 * t72 + t67 * t73 - t244 * pkin(3); m(6) * (t113 * t41 + t115 * t40 + t138 * t151 + t140 * t188 + t66 * t87 + t67 * t86) + (t151 + t140) * t57 + (t150 - pkin(3)) * t133 + (t138 + t188) * t106 + (-t207 - t208 + t198 + m(5) * (-pkin(3) * t159 + t237 * pkin(8)) + t177) * t214 + ((-t41 - t67) * t130 - (t40 + t66) * t172 - (t115 + t87) * t101 - (-t113 - t86) * t100) * mrSges(6,3) + t167; t188 * t230 + t151 * t232 - 0.2e1 * pkin(3) * t133 + (t113 * t67 + t115 * t66 + t151 * t188) * t233 + (t113 * t100 - t115 * t101 - t67 * t130 - t172 * t66) * t193 + t167; -Ifges(5,5) * t184 + t15 * mrSges(5,1) - t14 * mrSges(5,2) - t169 * Ifges(5,6) + (m(6) * (t157 * t3 - t16 * t195 + t161 * t4 + t17 * t194) + t72 * t194 + t157 * t19 - t73 * t195 + t161 * t18) * pkin(4) + t171 + t241; -t176 * t189 + (t141 * t149 - t215) * qJD(4) + (t187 + m(6) * (t157 * t40 + t161 * t41 + t194 * t87 - t195 * t86)) * pkin(4) + t170 + t240; (pkin(8) * t141 - t215) * qJD(4) + (t187 + m(6) * (-t113 * t195 + t115 * t194 + t157 * t66 + t161 * t67)) * pkin(4) + t170 + t239; 0.2e1 * t125; t171; t218 + t240; t218 + t239; t125; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
