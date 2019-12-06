% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRR9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:21:13
% EndTime: 2019-12-05 17:21:22
% DurationCPUTime: 2.61s
% Computational Cost: add. (9742->316), mult. (19055->453), div. (0->0), fcn. (13840->12), ass. (0->197)
t167 = sin(qJ(5));
t169 = sin(qJ(3));
t199 = qJD(2) * qJD(3);
t154 = t169 * t199;
t173 = cos(qJ(3));
t197 = t173 * qJDD(2);
t136 = -t154 + t197;
t129 = -qJDD(4) + t136;
t126 = -qJDD(5) + t129;
t168 = sin(qJ(4));
t172 = cos(qJ(4));
t203 = qJD(2) * t169;
t130 = -t172 * qJD(3) + t168 * t203;
t132 = t168 * qJD(3) + t172 * t203;
t171 = cos(qJ(5));
t112 = t171 * t130 + t167 * t132;
t114 = -t167 * t130 + t171 * t132;
t86 = t114 * t112;
t230 = -t126 - t86;
t234 = t167 * t230;
t233 = t171 * t230;
t163 = sin(pkin(5));
t164 = cos(pkin(5));
t214 = sin(pkin(10));
t215 = cos(pkin(10));
t180 = t214 * g(1) - t215 * g(2);
t178 = t164 * t180;
t204 = -g(3) + qJDD(1);
t232 = t163 * t204 + t178;
t175 = qJD(2) ^ 2;
t194 = t173 * t199;
t198 = t169 * qJDD(2);
t135 = t194 + t198;
t184 = -t168 * qJDD(3) - t172 * t135;
t108 = -t130 * qJD(4) - t184;
t190 = -t172 * qJDD(3) + t168 * t135;
t182 = t132 * qJD(4) + t190;
t64 = -t112 * qJD(5) + t171 * t108 - t167 * t182;
t151 = t173 * qJD(2) - qJD(4);
t145 = -qJD(5) + t151;
t98 = t112 * t145;
t231 = t98 + t64;
t213 = t132 * t130;
t179 = -t129 - t213;
t229 = t168 * t179;
t228 = t172 * t179;
t123 = t130 * t151;
t91 = t108 - t123;
t191 = t167 * t108 + t171 * t182;
t47 = (qJD(5) + t145) * t114 + t191;
t87 = (qJD(4) + t151) * t132 + t190;
t110 = t112 ^ 2;
t111 = t114 ^ 2;
t227 = t130 ^ 2;
t128 = t132 ^ 2;
t144 = t145 ^ 2;
t149 = t151 ^ 2;
t226 = qJD(3) ^ 2;
t177 = -t163 * t180 + t164 * t204;
t176 = t169 * t177;
t189 = -t173 * pkin(3) - t169 * pkin(8);
t140 = -t215 * g(1) - t214 * g(2);
t170 = sin(qJ(2));
t174 = cos(qJ(2));
t105 = t174 * t140 + t232 * t170;
t95 = -t175 * pkin(2) + qJDD(2) * pkin(7) + t105;
t192 = t175 * t189 + t95;
t68 = -t226 * pkin(3) + qJDD(3) * pkin(8) + t192 * t173 + t176;
t186 = -t136 + t154;
t187 = t135 + t194;
t188 = t170 * t140 - t232 * t174;
t94 = -qJDD(2) * pkin(2) - t175 * pkin(7) + t188;
t73 = t186 * pkin(3) - t187 * pkin(8) + t94;
t36 = t168 * t68 - t172 * t73;
t31 = t179 * pkin(4) - t91 * pkin(9) - t36;
t120 = -t151 * pkin(4) - t132 * pkin(9);
t37 = t168 * t73 + t172 * t68;
t34 = -t227 * pkin(4) - t182 * pkin(9) + t151 * t120 + t37;
t13 = t167 * t34 - t171 * t31;
t14 = t167 * t31 + t171 * t34;
t7 = -t171 * t13 + t167 * t14;
t225 = pkin(4) * t7;
t50 = -t98 + t64;
t26 = -t167 * t47 - t171 * t50;
t224 = pkin(4) * t26;
t223 = t168 * t7;
t222 = t172 * t7;
t119 = t173 * t177;
t67 = -qJDD(3) * pkin(3) - t226 * pkin(8) + t192 * t169 - t119;
t38 = t182 * pkin(4) - t227 * pkin(9) + t132 * t120 + t67;
t221 = t167 * t38;
t74 = t126 - t86;
t220 = t167 * t74;
t219 = t168 * t67;
t218 = t171 * t38;
t217 = t171 * t74;
t216 = t172 * t67;
t212 = t145 * t167;
t211 = t145 * t171;
t210 = t151 * t168;
t101 = t129 - t213;
t209 = t168 * t101;
t150 = t169 * t175 * t173;
t141 = qJDD(3) + t150;
t208 = t169 * t141;
t207 = t172 * t101;
t206 = t172 * t151;
t142 = qJDD(3) - t150;
t205 = t173 * t142;
t202 = qJD(4) - t151;
t196 = t173 * t86;
t195 = t173 * t213;
t8 = t167 * t13 + t171 * t14;
t21 = t168 * t36 + t172 * t37;
t83 = t169 * t95 - t119;
t84 = t173 * t95 + t176;
t52 = t169 * t83 + t173 * t84;
t20 = t168 * t37 - t172 * t36;
t185 = -pkin(2) + t189;
t80 = -t144 - t110;
t40 = t167 * t80 + t233;
t183 = pkin(4) * t40 - t13;
t93 = -t111 - t144;
t54 = t171 * t93 + t220;
t181 = pkin(4) * t54 - t14;
t160 = t173 ^ 2;
t159 = t169 ^ 2;
t158 = t160 * t175;
t156 = t159 * t175;
t148 = -t158 - t226;
t147 = -t156 - t226;
t139 = t156 + t158;
t138 = (t159 + t160) * qJDD(2);
t137 = -0.2e1 * t154 + t197;
t134 = 0.2e1 * t194 + t198;
t122 = -t128 + t149;
t121 = -t149 + t227;
t118 = -t169 * t147 - t205;
t117 = t173 * t148 - t208;
t116 = t128 - t227;
t115 = -t128 - t149;
t109 = -t149 - t227;
t100 = t128 + t227;
t97 = -t111 + t144;
t96 = t110 - t144;
t92 = t202 * t130 + t184;
t90 = t108 + t123;
t88 = -t202 * t132 - t190;
t85 = t111 - t110;
t82 = -t168 * t115 + t207;
t81 = t172 * t115 + t209;
t78 = t172 * t109 - t229;
t77 = t168 * t109 + t228;
t70 = (t112 * t171 - t114 * t167) * t145;
t69 = (t112 * t167 + t114 * t171) * t145;
t65 = -t110 - t111;
t63 = -t114 * qJD(5) - t191;
t62 = t168 * t91 - t172 * t87;
t61 = -t168 * t87 - t172 * t91;
t60 = t171 * t96 + t220;
t59 = -t167 * t97 + t233;
t58 = t167 * t96 - t217;
t57 = t171 * t97 + t234;
t56 = -t169 * t92 + t173 * t82;
t55 = -t167 * t93 + t217;
t53 = -t169 * t88 + t173 * t78;
t46 = (qJD(5) - t145) * t114 + t191;
t45 = t114 * t212 + t171 * t64;
t44 = -t114 * t211 + t167 * t64;
t43 = -t112 * t211 - t167 * t63;
t42 = -t112 * t212 + t171 * t63;
t41 = t171 * t80 - t234;
t39 = -t169 * t100 + t173 * t62;
t33 = -t168 * t54 + t172 * t55;
t32 = t168 * t55 + t172 * t54;
t29 = -t167 * t231 - t171 * t46;
t28 = t167 * t50 - t171 * t47;
t27 = -t167 * t46 + t171 * t231;
t25 = -pkin(9) * t54 + t218;
t24 = -t168 * t40 + t172 * t41;
t23 = t168 * t41 + t172 * t40;
t22 = -pkin(9) * t40 + t221;
t19 = t169 * t231 + t173 * t33;
t18 = -pkin(4) * t231 + pkin(9) * t55 + t221;
t17 = t169 * t67 + t173 * t21;
t16 = t169 * t46 + t173 * t24;
t15 = -pkin(4) * t46 + pkin(9) * t41 - t218;
t11 = -t168 * t26 + t172 * t28;
t10 = t168 * t28 + t172 * t26;
t9 = t173 * t11 + t169 * t65;
t6 = -pkin(4) * t38 + pkin(9) * t8;
t5 = -pkin(9) * t26 - t7;
t4 = -pkin(4) * t65 + pkin(9) * t28 + t8;
t3 = t172 * t8 - t223;
t2 = t168 * t8 + t222;
t1 = t169 * t38 + t173 * t3;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t204, 0, 0, 0, 0, 0, 0, (qJDD(2) * t174 - t170 * t175) * t163, (-qJDD(2) * t170 - t174 * t175) * t163, 0, t164 ^ 2 * t204 + (t170 * t105 - t174 * t188 - t178) * t163, 0, 0, 0, 0, 0, 0, t164 * (t173 * t141 + t169 * t148) + (t170 * t117 + t174 * t137) * t163, t164 * (-t169 * t142 + t173 * t147) + (t170 * t118 - t174 * t134) * t163, (t138 * t170 + t139 * t174) * t163, t164 * (t169 * t84 - t173 * t83) + (t170 * t52 - t174 * t94) * t163, 0, 0, 0, 0, 0, 0, t164 * (t169 * t78 + t173 * t88) + (t170 * t53 - t174 * t77) * t163, t164 * (t169 * t82 + t173 * t92) + (t170 * t56 - t174 * t81) * t163, t164 * (t173 * t100 + t169 * t62) + (t170 * t39 - t174 * t61) * t163, t164 * (t169 * t21 - t173 * t67) + (t170 * t17 - t174 * t20) * t163, 0, 0, 0, 0, 0, 0, t164 * (t169 * t24 - t173 * t46) + (t170 * t16 - t174 * t23) * t163, t164 * (t169 * t33 - t173 * t231) + (t170 * t19 - t174 * t32) * t163, t164 * (t169 * t11 - t173 * t65) + (-t174 * t10 + t170 * t9) * t163, t164 * (t169 * t3 - t173 * t38) + (t170 * t1 - t174 * t2) * t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t188, -t105, 0, 0, t187 * t169, t173 * t134 + t169 * t137, t208 + t173 * (-t156 + t226), -t186 * t173, t169 * (t158 - t226) + t205, 0, pkin(2) * t137 + pkin(7) * t117 - t173 * t94, -pkin(2) * t134 + pkin(7) * t118 + t169 * t94, pkin(2) * t139 + pkin(7) * t138 + t52, -pkin(2) * t94 + pkin(7) * t52, t169 * (t172 * t108 + t132 * t210) - t195, t169 * (-t168 * t90 + t172 * t88) - t173 * t116, t169 * (-t168 * t122 + t228) - t173 * t91, t169 * (-t130 * t206 + t168 * t182) + t195, t169 * (t172 * t121 + t209) + t173 * t87, t173 * t129 + t169 * (t130 * t172 - t132 * t168) * t151, t169 * (-pkin(8) * t77 + t219) + t173 * (-pkin(3) * t77 + t36) - pkin(2) * t77 + pkin(7) * t53, t169 * (-pkin(8) * t81 + t216) + t173 * (-pkin(3) * t81 + t37) - pkin(2) * t81 + pkin(7) * t56, pkin(7) * t39 - t169 * t20 + t185 * t61, pkin(7) * t17 + t185 * t20, t169 * (-t168 * t44 + t172 * t45) - t196, t169 * (-t168 * t27 + t172 * t29) - t173 * t85, t169 * (-t168 * t57 + t172 * t59) - t173 * t50, t169 * (-t168 * t42 + t172 * t43) + t196, t169 * (-t168 * t58 + t172 * t60) + t173 * t47, t169 * (-t168 * t69 + t172 * t70) + t173 * t126, t169 * (-pkin(8) * t23 - t168 * t15 + t172 * t22) + t173 * (-pkin(3) * t23 - t183) - pkin(2) * t23 + pkin(7) * t16, t169 * (-pkin(8) * t32 - t168 * t18 + t172 * t25) + t173 * (-pkin(3) * t32 - t181) - pkin(2) * t32 + pkin(7) * t19, t169 * (-pkin(8) * t10 - t168 * t4 + t172 * t5) + t173 * (-pkin(3) * t10 - t224) - pkin(2) * t10 + pkin(7) * t9, t169 * (-pkin(8) * t2 - pkin(9) * t222 - t168 * t6) + t173 * (-pkin(3) * t2 - t225) - pkin(2) * t2 + pkin(7) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150, t156 - t158, t198, t150, t197, qJDD(3), -t83, -t84, 0, 0, t168 * t108 - t132 * t206, t168 * t88 + t172 * t90, t172 * t122 + t229, -t130 * t210 - t172 * t182, t168 * t121 - t207, (t130 * t168 + t132 * t172) * t151, pkin(3) * t88 + pkin(8) * t78 - t216, pkin(3) * t92 + pkin(8) * t82 + t219, pkin(3) * t100 + pkin(8) * t62 + t21, -pkin(3) * t67 + pkin(8) * t21, t168 * t45 + t172 * t44, t168 * t29 + t172 * t27, t168 * t59 + t172 * t57, t168 * t43 + t172 * t42, t168 * t60 + t172 * t58, t168 * t70 + t172 * t69, -pkin(3) * t46 + pkin(8) * t24 + t172 * t15 + t168 * t22, -pkin(3) * t231 + pkin(8) * t33 + t168 * t25 + t172 * t18, -pkin(3) * t65 + pkin(8) * t11 + t168 * t5 + t172 * t4, -pkin(3) * t38 + pkin(8) * t3 - pkin(9) * t223 + t172 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, t116, t91, -t213, -t87, -t129, -t36, -t37, 0, 0, t86, t85, t50, -t86, -t47, -t126, t183, t181, t224, t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t85, t50, -t86, -t47, -t126, -t13, -t14, 0, 0;];
tauJ_reg = t12;
