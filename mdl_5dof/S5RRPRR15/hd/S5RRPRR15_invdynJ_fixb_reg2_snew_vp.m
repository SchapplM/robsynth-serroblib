% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRR15
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRR15_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:08
% EndTime: 2019-12-31 20:43:17
% DurationCPUTime: 2.66s
% Computational Cost: add. (10397->323), mult. (21650->414), div. (0->0), fcn. (13176->8), ass. (0->204)
t170 = sin(qJ(2));
t164 = t170 ^ 2;
t177 = qJD(1) ^ 2;
t159 = t164 * t177;
t176 = qJD(2) ^ 2;
t147 = -t159 - t176;
t174 = cos(qJ(2));
t216 = t170 * t177;
t204 = t174 * t216;
t142 = qJDD(2) - t204;
t215 = t174 * t142;
t256 = pkin(6) * (t170 * t147 + t215);
t168 = sin(qJ(5));
t209 = qJD(1) * qJD(2);
t156 = t174 * t209;
t158 = t170 * qJDD(1);
t135 = t158 + t156;
t125 = qJDD(4) + t135;
t121 = qJDD(5) + t125;
t169 = sin(qJ(4));
t173 = cos(qJ(4));
t213 = qJD(1) * t174;
t129 = t169 * qJD(2) + t173 * t213;
t131 = t173 * qJD(2) - t169 * t213;
t172 = cos(qJ(5));
t100 = t172 * t129 + t168 * t131;
t102 = -t168 * t129 + t172 * t131;
t75 = t102 * t100;
t251 = -t75 + t121;
t255 = t168 * t251;
t107 = t131 * t129;
t248 = -t107 + t125;
t254 = t169 * t248;
t253 = t172 * t251;
t252 = t173 * t248;
t212 = t170 * qJD(1);
t153 = qJD(4) + t212;
t115 = t153 * t129;
t202 = t170 * t209;
t207 = t174 * qJDD(1);
t136 = -t202 + t207;
t188 = t173 * qJDD(2) - t169 * t136;
t96 = -t129 * qJD(4) + t188;
t250 = -t96 - t115;
t194 = t135 + t156;
t249 = t194 * qJ(3);
t165 = t174 ^ 2;
t223 = t165 * t177;
t247 = t215 + t170 * (-t176 + t223);
t137 = -0.2e1 * t202 + t207;
t149 = -t176 - t223;
t141 = qJDD(2) + t204;
t219 = t170 * t141;
t245 = pkin(6) * (-t174 * t149 + t219) - pkin(1) * t137;
t98 = t100 ^ 2;
t99 = t102 ^ 2;
t123 = t129 ^ 2;
t124 = t131 ^ 2;
t145 = qJD(5) + t153;
t144 = t145 ^ 2;
t150 = t153 ^ 2;
t244 = 2 * qJD(3);
t243 = -pkin(2) - pkin(7);
t143 = pkin(3) * t212 - qJD(2) * pkin(7);
t152 = pkin(2) * t202;
t203 = qJD(3) * t212;
t155 = -0.2e1 * t203;
t171 = sin(qJ(1));
t175 = cos(qJ(1));
t201 = t171 * g(1) - t175 * g(2);
t186 = -qJDD(1) * pkin(1) - t201;
t60 = -t143 * t212 + t152 + t155 + (-pkin(3) * t165 - pkin(6)) * t177 + t243 * t136 - t249 + t186;
t222 = t170 * qJ(3);
t239 = t174 * pkin(2);
t195 = -t222 - t239;
t132 = t195 * qJD(1);
t196 = t175 * g(1) + t171 * g(2);
t228 = qJDD(1) * pkin(6);
t119 = -t177 * pkin(1) - t196 + t228;
t221 = t170 * t119;
t184 = -qJDD(2) * pkin(2) - t176 * qJ(3) + t132 * t212 + qJDD(3) + t221;
t70 = t135 * pkin(3) - qJDD(2) * pkin(7) + (-pkin(3) * t209 - pkin(7) * t216 + g(3)) * t174 + t184;
t32 = t169 * t60 - t173 * t70;
t30 = t248 * pkin(4) + t250 * pkin(8) - t32;
t112 = t153 * pkin(4) - t131 * pkin(8);
t33 = t169 * t70 + t173 * t60;
t198 = t169 * qJDD(2) + t173 * t136;
t95 = -t131 * qJD(4) - t198;
t31 = -t123 * pkin(4) + t95 * pkin(8) - t153 * t112 + t33;
t12 = t168 * t31 - t172 * t30;
t13 = t168 * t30 + t172 * t31;
t7 = -t172 * t12 + t168 * t13;
t241 = t169 * t7;
t161 = t170 * g(3);
t240 = t173 * t7;
t238 = t174 * g(3);
t182 = (qJD(1) * t132 + t119) * t174 - t161 - t176 * pkin(2);
t208 = qJDD(2) * qJ(3);
t69 = t208 + t136 * pkin(3) - pkin(7) * t223 + (t244 + t143) * qJD(2) + t182;
t34 = -t95 * pkin(4) - t123 * pkin(8) + t131 * t112 + t69;
t237 = t168 * t34;
t65 = t75 + t121;
t236 = t168 * t65;
t235 = t169 * t69;
t91 = t107 + t125;
t234 = t169 * t91;
t181 = (-qJD(4) + t153) * t131 - t198;
t55 = t169 * t181 + t173 * t250;
t233 = t170 * t55;
t232 = t172 * t34;
t231 = t172 * t65;
t230 = t173 * t69;
t229 = t173 * t91;
t227 = t145 * t168;
t226 = t145 * t172;
t225 = t153 * t169;
t224 = t153 * t173;
t220 = t170 * t137;
t139 = t159 + t223;
t214 = pkin(1) * t139 + (t164 + t165) * t228;
t211 = qJD(4) + t153;
t210 = qJD(5) + t145;
t206 = t170 * t75;
t205 = t170 * t107;
t8 = t168 * t12 + t172 * t13;
t200 = t168 * t96 - t172 * t95;
t110 = t221 + t238;
t111 = t174 * t119 - t161;
t199 = t170 * t110 + t174 * t111;
t71 = -t144 - t98;
t36 = t168 * t71 + t253;
t197 = pkin(4) * t36 - t12;
t192 = t168 * t95 + t172 * t96;
t16 = t169 * t33 - t173 * t32;
t191 = t169 * t32 + t173 * t33;
t190 = t174 * (-t159 + t176) + t219;
t83 = -t99 - t144;
t49 = t172 * t83 - t236;
t187 = pkin(4) * t49 - t13;
t185 = t174 * t243 - pkin(1) - t222;
t183 = (-qJD(5) + t145) * t102 - t200;
t58 = -t100 * qJD(5) + t192;
t118 = t177 * pkin(6) - t186;
t85 = t184 + t238;
t180 = qJD(2) * t244 + t182;
t179 = t136 * pkin(2) + t118 - t152;
t84 = t180 + t208;
t178 = t179 + 0.2e1 * t203;
t140 = t159 - t223;
t134 = t158 + 0.2e1 * t156;
t114 = -t124 + t150;
t113 = t123 - t150;
t109 = t194 * t170;
t108 = (t136 - t202) * t174;
t105 = t124 - t123;
t104 = -t124 - t150;
t103 = t174 * t134 + t220;
t97 = -t150 - t123;
t89 = -t123 - t124;
t88 = t145 * t100;
t87 = -t99 + t144;
t86 = t98 - t144;
t81 = -t115 + t96;
t80 = -t211 * t129 + t188;
t77 = t211 * t131 + t198;
t74 = t99 - t98;
t72 = t173 * t104 - t234;
t67 = t169 * t97 + t252;
t62 = (-t100 * t172 + t102 * t168) * t145;
t61 = (-t100 * t168 - t102 * t172) * t145;
t59 = -t98 - t99;
t57 = -t102 * qJD(5) - t200;
t54 = t172 * t86 - t236;
t53 = -t168 * t87 + t253;
t52 = t168 * t86 + t231;
t51 = t172 * t87 + t255;
t50 = -t168 * t83 - t231;
t47 = -t210 * t100 + t192;
t46 = t58 + t88;
t45 = t58 - t88;
t42 = t210 * t102 + t200;
t41 = -t102 * t227 + t172 * t58;
t40 = t102 * t226 + t168 * t58;
t39 = t100 * t226 - t168 * t57;
t38 = t100 * t227 + t172 * t57;
t37 = t172 * t71 - t255;
t27 = t169 * t50 + t173 * t49;
t26 = -pkin(8) * t49 + t232;
t25 = t168 * t46 + t172 * t183;
t24 = -t168 * t45 - t172 * t42;
t23 = t168 * t183 - t172 * t46;
t22 = -t168 * t42 + t172 * t45;
t21 = pkin(4) * t23;
t19 = t169 * t37 + t173 * t36;
t18 = -pkin(8) * t36 + t237;
t15 = -pkin(4) * t47 + pkin(8) * t50 + t237;
t14 = -pkin(4) * t42 + pkin(8) * t37 - t232;
t9 = t169 * t25 + t173 * t23;
t6 = pkin(4) * t7;
t5 = -pkin(4) * t34 + pkin(8) * t8;
t4 = -pkin(8) * t23 - t7;
t3 = -pkin(4) * t59 + pkin(8) * t25 + t8;
t1 = t169 * t8 + t240;
t2 = [0, 0, 0, 0, 0, qJDD(1), t201, t196, 0, 0, t109, t103, t190, t108, t247, 0, t174 * t118 - t245, -pkin(1) * t134 - t170 * t118 - t256, t199 + t214, pkin(1) * t118 + pkin(6) * t199, 0, -t190, -t247, t109, t103, t108, t170 * (qJ(3) * t139 + t184) + (pkin(2) * t139 + t161 + t84) * t174 + t214, t174 * (-pkin(2) * t137 + t155 - t179) + (-t174 * t194 - t220) * qJ(3) + t245, t170 * t178 + t256 + (pkin(1) + t239) * t134 + (t134 + t194) * t222, pkin(6) * (t170 * t85 + t174 * t84) + (pkin(1) - t195) * (t178 + t249), t205 + t174 * (-t131 * t224 - t169 * t96), t170 * t105 + t174 * (t169 * t77 - t173 * t81), -t170 * t250 + t174 * (-t173 * t114 - t254), -t205 + t174 * (-t129 * t225 - t173 * t95), t170 * t181 + t174 * (-t169 * t113 - t229), t170 * t125 + t174 * (t129 * t169 + t131 * t173) * t153, t170 * (pkin(3) * t67 - t32) + t174 * (pkin(3) * t77 + t230) + pkin(6) * (t170 * t67 + t174 * t77) + t185 * (t173 * t97 - t254), t170 * (pkin(3) * t72 - t33) + t174 * (pkin(3) * t80 - t235) + pkin(6) * (t170 * t72 + t174 * t80) + t185 * (-t169 * t104 - t229), pkin(3) * t233 + t174 * (pkin(3) * t89 - t191) + pkin(6) * (t174 * t89 + t233) + t185 * (-t169 * t250 + t173 * t181), t185 * t191 + (pkin(3) + pkin(6)) * (t170 * t16 + t174 * t69), t206 + t174 * (-t169 * t41 - t173 * t40), t170 * t74 + t174 * (-t169 * t24 - t173 * t22), t170 * t46 + t174 * (-t169 * t53 - t173 * t51), -t206 + t174 * (-t169 * t39 - t173 * t38), t170 * t183 + t174 * (-t169 * t54 - t173 * t52), t170 * t121 + t174 * (-t169 * t62 - t173 * t61), t170 * (pkin(3) * t19 + t197) + t174 * (pkin(3) * t42 - t173 * t14 - t169 * t18) + pkin(6) * (t170 * t19 + t174 * t42) + t185 * (-t169 * t36 + t173 * t37), t170 * (pkin(3) * t27 + t187) + t174 * (pkin(3) * t47 - t173 * t15 - t169 * t26) + pkin(6) * (t170 * t27 + t174 * t47) + t185 * (-t169 * t49 + t173 * t50), t170 * (pkin(3) * t9 + t21) + t174 * (pkin(3) * t59 - t169 * t4 - t173 * t3) + pkin(6) * (t170 * t9 + t174 * t59) + t185 * (-t169 * t23 + t173 * t25), t170 * (pkin(3) * t1 + t6) + t174 * (pkin(3) * t34 + pkin(8) * t241 - t173 * t5) + pkin(6) * (t170 * t1 + t174 * t34) + t185 * (t173 * t8 - t241); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204, t140, t158, t204, t207, qJDD(2), -t110, -t111, 0, 0, qJDD(2), -t158, -t207, -t204, t140, t204, (-pkin(2) * t170 + qJ(3) * t174) * qJDD(1), -pkin(2) * t141 - qJ(3) * t149 + t85, -pkin(2) * t147 + (qJDD(2) + t142) * qJ(3) + t180, -pkin(2) * t85 + qJ(3) * t84, -t131 * t225 + t173 * t96, -t169 * t81 - t173 * t77, -t169 * t114 + t252, t129 * t224 - t169 * t95, t173 * t113 - t234, (-t129 * t173 + t131 * t169) * t153, qJ(3) * t77 + t243 * t67 + t235, qJ(3) * t80 + t243 * t72 + t230, qJ(3) * t89 + t243 * t55 - t16, qJ(3) * t69 + t243 * t16, -t169 * t40 + t173 * t41, -t169 * t22 + t173 * t24, -t169 * t51 + t173 * t53, -t169 * t38 + t173 * t39, -t169 * t52 + t173 * t54, -t169 * t61 + t173 * t62, qJ(3) * t42 - t169 * t14 + t173 * t18 + t243 * t19, qJ(3) * t47 - t169 * t15 + t173 * t26 + t243 * t27, qJ(3) * t59 - t169 * t3 + t173 * t4 + t243 * t9, -pkin(8) * t240 + qJ(3) * t34 + t243 * t1 - t169 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t141, t147, t85, 0, 0, 0, 0, 0, 0, t67, t72, t55, t16, 0, 0, 0, 0, 0, 0, t19, t27, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, t105, -t250, -t107, t181, t125, -t32, -t33, 0, 0, t75, t74, t46, -t75, t183, t121, t197, t187, t21, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t74, t46, -t75, t183, t121, -t12, -t13, 0, 0;];
tauJ_reg = t2;
