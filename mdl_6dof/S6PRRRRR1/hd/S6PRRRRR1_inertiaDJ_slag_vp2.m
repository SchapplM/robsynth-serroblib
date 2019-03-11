% Calculate time derivative of joint inertia matrix for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:37:38
% EndTime: 2019-03-09 00:37:49
% DurationCPUTime: 4.61s
% Computational Cost: add. (7840->391), mult. (18566->600), div. (0->0), fcn. (18722->12), ass. (0->186)
t238 = -qJD(4) - qJD(3);
t135 = sin(qJ(6));
t131 = t135 ^ 2;
t140 = cos(qJ(6));
t132 = t140 ^ 2;
t245 = t131 + t132;
t137 = sin(qJ(4));
t142 = cos(qJ(4));
t138 = sin(qJ(3));
t229 = -pkin(9) - pkin(8);
t182 = t229 * t138;
t143 = cos(qJ(3));
t240 = t229 * t143;
t90 = t137 * t240 + t142 * t182;
t216 = mrSges(7,1) * t140;
t116 = mrSges(7,2) * t135 - t216;
t244 = t116 - mrSges(6,1);
t185 = qJD(6) * t140;
t136 = sin(qJ(5));
t141 = cos(qJ(5));
t108 = -t137 * t138 + t142 * t143;
t109 = t137 * t143 + t138 * t142;
t162 = t108 * t141 - t109 * t136;
t85 = t238 * t108;
t86 = t238 * t109;
t41 = qJD(5) * t162 + t136 * t86 - t141 * t85;
t81 = t108 * t136 + t109 * t141;
t161 = t135 * t41 + t185 * t81;
t243 = -t137 * t182 + t142 * t240;
t127 = -pkin(3) * t143 - pkin(2);
t92 = -pkin(4) * t108 + t127;
t49 = -pkin(5) * t162 - pkin(11) * t81 + t92;
t155 = -t109 * pkin(10) + t90;
t72 = pkin(10) * t108 - t243;
t51 = t136 * t155 + t141 * t72;
t25 = t135 * t49 + t140 * t51;
t191 = t25 * qJD(6);
t61 = t238 * t243;
t147 = t85 * pkin(10) - t61;
t188 = qJD(5) * t136;
t60 = t238 * t90;
t241 = pkin(10) * t86 + qJD(5) * t155 - t60;
t18 = t136 * t147 + t141 * t241 - t188 * t72;
t42 = qJD(5) * t81 - t136 * t85 - t141 * t86;
t129 = qJD(3) * t138 * pkin(3);
t75 = -pkin(4) * t86 + t129;
t22 = pkin(5) * t42 - pkin(11) * t41 + t75;
t3 = -t135 * t18 + t140 * t22 - t191;
t242 = -t3 - t191;
t239 = t245 * t141;
t237 = -mrSges(5,1) * t61 + t60 * mrSges(5,2) - Ifges(5,5) * t85 + Ifges(5,6) * t86;
t105 = pkin(4) * t116 * t188;
t187 = qJD(5) * t141;
t183 = pkin(4) * t187;
t172 = mrSges(7,3) * t183;
t119 = t131 * t172;
t120 = t132 * t172;
t166 = mrSges(7,1) * t135 + mrSges(7,2) * t140;
t112 = t166 * qJD(6);
t125 = -pkin(4) * t141 - pkin(5);
t93 = t125 * t112;
t236 = t105 + t119 + t120 + t93;
t235 = 2 * m(6);
t234 = 2 * m(7);
t233 = -2 * mrSges(6,3);
t19 = t136 * t241 - t141 * t147 + t187 * t72;
t232 = 0.2e1 * t19;
t50 = t136 * t72 - t141 * t155;
t231 = 0.2e1 * t50;
t230 = m(6) / 0.2e1;
t227 = pkin(5) * t112;
t134 = cos(pkin(6));
t133 = sin(pkin(6));
t139 = sin(qJ(2));
t198 = t133 * t139;
t98 = t134 * t143 - t138 * t198;
t99 = t134 * t138 + t143 * t198;
t73 = -t137 * t99 + t142 * t98;
t144 = cos(qJ(2));
t189 = qJD(2) * t144;
t176 = t133 * t189;
t88 = -qJD(3) * t99 - t138 * t176;
t89 = qJD(3) * t98 + t143 * t176;
t35 = qJD(4) * t73 + t137 * t88 + t142 * t89;
t74 = t137 * t98 + t142 * t99;
t36 = -qJD(4) * t74 - t137 * t89 + t142 * t88;
t52 = t136 * t74 - t141 * t73;
t13 = -qJD(5) * t52 + t136 * t36 + t141 * t35;
t190 = qJD(2) * t139;
t177 = t133 * t190;
t197 = t133 * t144;
t53 = t136 * t73 + t141 * t74;
t44 = -t135 * t53 - t140 * t197;
t5 = qJD(6) * t44 + t13 * t140 + t135 * t177;
t226 = t140 * t5;
t225 = t19 * t50;
t224 = t3 * t135;
t126 = pkin(3) * t142 + pkin(4);
t193 = t137 * t141;
t78 = t126 * t188 + (t137 * t187 + (t136 * t142 + t193) * qJD(4)) * pkin(3);
t223 = t50 * t78;
t14 = qJD(5) * t53 + t136 * t35 - t141 * t36;
t222 = t52 * t14;
t221 = t52 * t78;
t159 = t135 * t197 - t140 * t53;
t6 = qJD(6) * t159 - t135 * t13 + t140 * t177;
t220 = t6 * t135;
t195 = t136 * t137;
t77 = t126 * t187 + (-t137 * t188 + (t141 * t142 - t195) * qJD(4)) * pkin(3);
t218 = t77 * mrSges(6,2);
t204 = t140 * t41;
t217 = Ifges(7,5) * t204 + Ifges(7,3) * t42;
t215 = mrSges(7,3) * t140;
t214 = Ifges(7,4) * t135;
t213 = Ifges(7,4) * t140;
t212 = Ifges(7,6) * t135;
t211 = pkin(4) * qJD(5);
t210 = t131 * t77;
t209 = t132 * t77;
t207 = t135 * t81;
t206 = t136 * t50;
t205 = t136 * t52;
t203 = t140 * t81;
t24 = -t135 * t51 + t140 * t49;
t202 = qJD(6) * t24;
t97 = pkin(3) * t193 + t126 * t136;
t95 = pkin(11) + t97;
t201 = qJD(6) * t95;
t196 = t135 * t141;
t192 = t140 * t141;
t186 = qJD(6) * t135;
t184 = 0.2e1 * t138;
t180 = t81 * t186;
t160 = t180 - t204;
t15 = mrSges(7,1) * t161 - mrSges(7,2) * t160;
t178 = m(7) * t19 + t15;
t175 = -t186 / 0.2e1;
t174 = -(2 * Ifges(6,4)) - t212;
t173 = t245 * t77;
t170 = t133 ^ 2 * t139 * t189;
t169 = t14 * t50 + t19 * t52;
t168 = -mrSges(4,1) * t143 + mrSges(4,2) * t138;
t167 = -t136 * mrSges(6,1) - t141 * mrSges(6,2);
t165 = Ifges(7,1) * t140 - t214;
t164 = -Ifges(7,2) * t135 + t213;
t163 = Ifges(7,5) * t135 + Ifges(7,6) * t140;
t96 = -pkin(3) * t195 + t126 * t141;
t114 = t164 * qJD(6);
t115 = t165 * qJD(6);
t117 = Ifges(7,2) * t140 + t214;
t118 = Ifges(7,1) * t135 + t213;
t158 = t114 * t140 + t115 * t135 - t117 * t186 + t118 * t185;
t157 = (-mrSges(5,1) * t137 - mrSges(5,2) * t142) * qJD(4) * pkin(3);
t156 = -t220 + (t135 * t159 - t140 * t44) * qJD(6);
t153 = t156 + t226;
t152 = -t88 * t138 + t89 * t143 + (-t138 * t99 - t143 * t98) * qJD(3);
t67 = t78 * t116;
t70 = mrSges(7,3) * t210;
t71 = mrSges(7,3) * t209;
t76 = t78 * mrSges(6,1);
t94 = -pkin(5) - t96;
t84 = t94 * t112;
t151 = t158 + t67 + t70 + t71 - t76 + t84 - t218;
t150 = -t13 * mrSges(6,2) + mrSges(7,3) * t156 + t52 * t112 + t14 * t244 + t5 * t215;
t10 = -Ifges(7,1) * t160 - Ifges(7,4) * t161 + Ifges(7,5) * t42;
t128 = Ifges(7,5) * t185;
t2 = t135 * t22 + t140 * t18 + t202;
t32 = -Ifges(7,6) * t162 + t164 * t81;
t33 = -Ifges(7,5) * t162 + t165 * t81;
t9 = -Ifges(7,4) * t160 - Ifges(7,2) * t161 + Ifges(7,6) * t42;
t149 = -t18 * mrSges(6,2) + t2 * t215 + t32 * t175 + t33 * t185 / 0.2e1 + Ifges(6,5) * t41 + t50 * t112 - t114 * t207 / 0.2e1 + t115 * t203 / 0.2e1 - t162 * (-Ifges(7,6) * t186 + t128) / 0.2e1 + t135 * t10 / 0.2e1 + t140 * t9 / 0.2e1 + (t163 / 0.2e1 - Ifges(6,6)) * t42 - t161 * t117 / 0.2e1 + t244 * t19 + (t204 / 0.2e1 + t81 * t175) * t118;
t148 = mrSges(5,1) * t36 - t35 * mrSges(5,2) + t150;
t20 = mrSges(7,1) * t42 + mrSges(7,3) * t160;
t21 = -mrSges(7,2) * t42 - mrSges(7,3) * t161;
t55 = mrSges(7,2) * t162 - mrSges(7,3) * t207;
t56 = -mrSges(7,1) * t162 - mrSges(7,3) * t203;
t146 = m(7) * (t140 * t2 - t185 * t24 - t186 * t25 - t224) + t140 * t21 - t135 * t20 - t55 * t186 - t56 * t185;
t145 = (-t224 + (-t135 * t25 - t140 * t24) * qJD(6)) * mrSges(7,3) + t149;
t124 = pkin(4) * t136 + pkin(11);
t113 = (mrSges(4,1) * t138 + mrSges(4,2) * t143) * qJD(3);
t87 = -mrSges(5,1) * t108 + mrSges(5,2) * t109;
t58 = -mrSges(5,1) * t86 - mrSges(5,2) * t85;
t57 = -mrSges(6,1) * t162 + mrSges(6,2) * t81;
t54 = t166 * t81;
t23 = mrSges(6,1) * t42 + mrSges(6,2) * t41;
t1 = [0.2e1 * m(7) * (-t159 * t5 + t44 * t6 + t222) + 0.2e1 * m(6) * (t13 * t53 - t170 + t222) + 0.2e1 * m(5) * (t35 * t74 + t36 * t73 - t170) + 0.2e1 * m(4) * (t88 * t98 + t89 * t99 - t170); t14 * t54 + t52 * t15 + t44 * t20 - t159 * t21 + t5 * t55 + t6 * t56 + (t13 * t162 + t14 * t81 + t41 * t52 - t42 * t53) * mrSges(6,3) + (t108 * t35 - t109 * t36 + t73 * t85 + t74 * t86) * mrSges(5,3) + t152 * mrSges(4,3) + ((-t113 - t23 - t58) * t144 + (-t144 * mrSges(3,2) + (-mrSges(3,1) + t168 + t57 + t87) * t139) * qJD(2)) * t133 + m(5) * (-t243 * t35 + t90 * t36 - t60 * t74 - t61 * t73 + (t127 * t190 - t129 * t144) * t133) + m(6) * (t51 * t13 + t18 * t53 + (-t144 * t75 + t190 * t92) * t133 + t169) + m(7) * (-t159 * t2 + t24 * t6 + t25 * t5 + t3 * t44 + t169) + (-pkin(2) * t177 + pkin(8) * t152) * m(4); t51 * t42 * t233 - 0.2e1 * t85 * t109 * Ifges(5,1) + 0.2e1 * t108 * Ifges(5,2) * t86 - 0.2e1 * pkin(2) * t113 + 0.2e1 * t127 * t58 + t15 * t231 + t54 * t232 + 0.2e1 * t2 * t55 + 0.2e1 * t24 * t20 + 0.2e1 * t25 * t21 + 0.2e1 * t92 * t23 + 0.2e1 * t3 * t56 + 0.2e1 * t75 * t57 + (mrSges(6,3) * t231 - t135 * t32 + t140 * t33) * t41 + 0.2e1 * m(5) * (t127 * t129 + t243 * t60 - t61 * t90) + (t18 * t51 + t75 * t92 + t225) * t235 + (t2 * t25 + t24 * t3 + t225) * t234 - (t18 * t233 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t42 + t174 * t41 + t217) * t162 + (mrSges(6,3) * t232 + 0.2e1 * Ifges(6,1) * t41 + t140 * t10 - t135 * t9 + (Ifges(7,5) * t140 + t174) * t42 + (-t135 * t33 - t140 * t32 + t162 * t163) * qJD(6)) * t81 + ((-Ifges(4,4) * t138 + pkin(3) * t87) * t184 + (0.2e1 * Ifges(4,4) * t143 + (Ifges(4,1) - Ifges(4,2)) * t184) * t143) * qJD(3) + 0.2e1 * (-t108 * t85 + t109 * t86) * Ifges(5,4) + 0.2e1 * (-t108 * t60 + t109 * t61 - t243 * t86 + t85 * t90) * mrSges(5,3); t88 * mrSges(4,1) - t89 * mrSges(4,2) + m(6) * (t13 * t97 - t14 * t96 + t53 * t77 + t221) + m(5) * (t137 * t35 + t142 * t36 + (-t137 * t73 + t142 * t74) * qJD(4)) * pkin(3) + t148 + (t14 * t94 + t153 * t95 + t221 + (-t135 * t44 - t140 * t159) * t77) * m(7); m(6) * (t18 * t97 - t19 * t96 + t51 * t77 + t223) + (t162 * t77 - t41 * t96 - t42 * t97 + t78 * t81) * mrSges(6,3) + t149 + (m(7) * (t2 * t95 - t201 * t24 + t25 * t77) + t95 * t21 + t77 * t55 - mrSges(7,3) * t202 - t56 * t201) * t140 + (-t55 * t201 + t242 * mrSges(7,3) + (-m(7) * t24 - t56) * t77 + (m(7) * t242 - t20) * t95) * t135 + (m(5) * (-t137 * t60 - t142 * t61 + (-t137 * t90 - t142 * t243) * qJD(4)) + (t137 * t86 + t142 * t85 + (t108 * t142 + t109 * t137) * qJD(4)) * mrSges(5,3)) * pkin(3) + t94 * t15 + t78 * t54 + m(7) * (t19 * t94 + t223) + (Ifges(4,5) * t143 - Ifges(4,6) * t138 + pkin(8) * t168) * qJD(3) + t237; -0.2e1 * t218 + 0.2e1 * t67 + 0.2e1 * t70 + 0.2e1 * t71 - 0.2e1 * t76 + 0.2e1 * t84 + 0.2e1 * t157 + (t173 * t95 + t78 * t94) * t234 + (t77 * t97 - t78 * t96) * t235 + t158; m(7) * t125 * t14 + 0.2e1 * ((t13 * t136 - t14 * t141) * t230 + (m(7) * (-t159 * t192 - t44 * t196 + t205) / 0.2e1 + (t141 * t53 + t205) * t230) * qJD(5)) * pkin(4) + t148 + m(7) * (t159 * t186 - t185 * t44 - t220 + t226) * t124; t145 + (m(6) * (t136 * t18 - t141 * t19) + (-t136 * t42 - t141 * t41) * mrSges(6,3) + ((mrSges(6,3) * t81 + t54) * t136 + (mrSges(6,3) * t162 - t135 * t56 + t140 * t55) * t141 + m(7) * (t192 * t25 - t196 * t24 + t206) + m(6) * (t141 * t51 + t206)) * qJD(5)) * pkin(4) + t146 * t124 + t178 * t125 + t237; t157 + t151 + m(7) * (t125 * t78 + (t209 + t210) * t124) + (m(6) * (t136 * t77 - t141 * t78) + (m(7) * (t136 * t94 + t239 * t95) + m(6) * (-t136 * t96 + t141 * t97) + t167) * qJD(5)) * pkin(4) + t236; 0.2e1 * t105 + 0.2e1 * t119 + 0.2e1 * t120 + 0.2e1 * t93 + 0.2e1 * (m(7) * (t124 * t239 + t125 * t136) + t167) * t211 + t158; m(7) * (-pkin(5) * t14 + pkin(11) * t153) + t150; -pkin(5) * t178 + pkin(11) * t146 + t145; m(7) * (-pkin(5) * t78 + pkin(11) * t173) - t227 + t151; -t227 + (m(7) * (-pkin(5) * t136 + pkin(11) * t239) + t167) * t211 + t158 + t236; t158 - 0.2e1 * t227; mrSges(7,1) * t6 - mrSges(7,2) * t5; mrSges(7,1) * t3 - mrSges(7,2) * t2 - Ifges(7,5) * t180 - Ifges(7,6) * t161 + t217; t128 - t166 * t77 + (-t95 * t216 + (mrSges(7,2) * t95 - Ifges(7,6)) * t135) * qJD(6); t128 - t166 * t183 + (t116 * t124 - t212) * qJD(6); t128 + (pkin(11) * t116 - t212) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
