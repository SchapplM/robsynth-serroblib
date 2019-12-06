% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:57
% EndTime: 2019-12-05 17:45:14
% DurationCPUTime: 4.23s
% Computational Cost: add. (14075->371), mult. (39105->547), div. (0->0), fcn. (27680->10), ass. (0->227)
t175 = sin(qJ(5));
t174 = cos(pkin(8));
t219 = t174 * qJDD(1);
t160 = -qJDD(4) + t219;
t154 = -qJDD(5) + t160;
t172 = sin(pkin(8));
t171 = sin(pkin(9));
t173 = cos(pkin(9));
t176 = sin(qJ(4));
t179 = cos(qJ(4));
t195 = t171 * t179 + t173 * t176;
t189 = t172 * t195;
t131 = qJD(1) * t189;
t232 = t172 * t173;
t213 = t179 * t232;
t229 = qJD(1) * t172;
t133 = -t176 * t171 * t229 + qJD(1) * t213;
t178 = cos(qJ(5));
t104 = t178 * t131 + t133 * t175;
t106 = -t131 * t175 + t133 * t178;
t80 = t106 * t104;
t265 = -t154 - t80;
t269 = t175 * t265;
t268 = t178 * t265;
t220 = t172 * qJDD(1);
t209 = t171 * t220;
t190 = qJDD(1) * t213 - t176 * t209;
t113 = -qJD(4) * t131 + t190;
t228 = qJD(1) * t174;
t161 = -qJD(4) + t228;
t244 = t131 * t161;
t91 = t113 - t244;
t181 = qJD(1) ^ 2;
t177 = sin(qJ(1));
t180 = cos(qJ(1));
t200 = g(2) * t180 + g(3) * t177;
t184 = -qJ(2) * t181 + qJDD(2) - t200;
t198 = -t174 * pkin(2) - t172 * qJ(3);
t192 = -pkin(1) + t198;
t267 = -0.2e1 * qJD(3) * t229 + t192 * qJDD(1) + t184;
t199 = -g(2) * t177 + g(3) * t180;
t143 = -pkin(1) * t181 + qJDD(1) * qJ(2) - t199;
t260 = 2 * qJD(2);
t266 = qJD(1) * t260 + t143;
t243 = t133 * t131;
t188 = -t160 - t243;
t264 = t176 * t188;
t263 = t179 * t188;
t262 = (qJD(4) + t161) * t133;
t245 = qJDD(1) * pkin(1);
t139 = -t184 + t245;
t170 = t174 ^ 2;
t235 = t170 * t181;
t168 = t172 ^ 2;
t236 = t168 * t181;
t261 = -t139 - t245 + (t235 + t236) * qJ(2);
t156 = -qJD(5) + t161;
t186 = qJDD(1) * t189;
t112 = -qJD(4) * t133 - t186;
t206 = -t178 * t112 + t113 * t175;
t48 = (qJD(5) + t156) * t106 + t206;
t102 = t104 ^ 2;
t103 = t106 ^ 2;
t129 = t131 ^ 2;
t130 = t133 ^ 2;
t153 = t156 ^ 2;
t159 = t161 ^ 2;
t147 = t198 * qJD(1);
t201 = -g(1) * t172 + t266 * t174;
t109 = t147 * t228 + t201;
t191 = t267 * t173;
t193 = -t174 * pkin(3) - pkin(6) * t232;
t237 = t168 * t173;
t69 = t193 * qJDD(1) + (-t109 + (pkin(6) * t172 * t174 - pkin(3) * t237) * t181) * t171 + t191;
t144 = t193 * qJD(1);
t167 = t171 ^ 2;
t158 = t167 * t236;
t84 = t173 * t109 + t267 * t171;
t70 = -pkin(3) * t158 - pkin(6) * t209 + t144 * t228 + t84;
t37 = t176 * t70 - t179 * t69;
t31 = t188 * pkin(4) - t91 * pkin(7) - t37;
t121 = -pkin(4) * t161 - pkin(7) * t133;
t38 = t176 * t69 + t179 * t70;
t32 = -pkin(4) * t129 + pkin(7) * t112 + t121 * t161 + t38;
t15 = t175 * t32 - t178 * t31;
t16 = t175 * t31 + t178 * t32;
t8 = -t15 * t178 + t16 * t175;
t259 = pkin(4) * t8;
t196 = t112 * t175 + t113 * t178;
t65 = -qJD(5) * t104 + t196;
t96 = t104 * t156;
t51 = t65 - t96;
t26 = -t175 * t48 - t178 * t51;
t258 = pkin(4) * t26;
t257 = g(1) * t174;
t256 = t176 * t8;
t255 = t179 * t8;
t20 = t176 * t38 - t179 * t37;
t254 = t173 * t20;
t207 = qJDD(3) + t257;
t222 = qJDD(1) * t171;
t227 = t260 + t147;
t85 = -pkin(6) * t158 + (pkin(3) * t222 + t143 + (t144 * t173 + t227) * qJD(1)) * t172 + t207;
t46 = -t112 * pkin(4) - t129 * pkin(7) + t121 * t133 + t85;
t253 = t175 * t46;
t76 = t154 - t80;
t252 = t175 * t76;
t251 = t176 * t85;
t99 = t160 - t243;
t250 = t176 * t99;
t249 = t178 * t46;
t248 = t178 * t76;
t247 = t179 * t85;
t246 = t179 * t99;
t242 = t156 * t175;
t241 = t156 * t178;
t240 = t161 * t176;
t239 = t161 * t179;
t169 = t173 ^ 2;
t238 = t168 * t169;
t233 = t171 * t181;
t203 = t233 * t237;
t140 = -t203 + t219;
t234 = t171 * t140;
t141 = -t203 - t219;
t231 = t173 * t141;
t230 = t174 * t181;
t226 = qJD(4) - t161;
t224 = qJD(5) - t156;
t221 = qJDD(1) * t173;
t215 = t169 * t236;
t214 = t171 * t230;
t212 = t173 * t230;
t211 = t174 * t80;
t210 = t174 * t243;
t9 = t15 * t175 + t178 * t16;
t21 = t176 * t37 + t179 * t38;
t205 = t172 * (t266 * t172 + t257) + t174 * t201;
t202 = t171 * t212;
t83 = t109 * t171 - t191;
t54 = t171 * t84 - t173 * t83;
t197 = t168 * t202;
t73 = -t153 - t102;
t40 = t175 * t73 + t268;
t194 = pkin(4) * t40 - t15;
t92 = -t103 - t153;
t55 = t178 * t92 + t252;
t187 = pkin(4) * t55 - t16;
t185 = t195 * t220;
t166 = t170 * qJDD(1);
t165 = t168 * qJDD(1);
t148 = (t168 + t170) * t181;
t146 = (-t170 - t238) * t181;
t145 = -t158 - t235;
t142 = t158 + t215;
t137 = (t214 - t221) * t172;
t136 = (t214 + t221) * t172;
t135 = (-t212 + t222) * t172;
t134 = (t212 + t222) * t172;
t123 = -t130 + t159;
t122 = t129 - t159;
t117 = -t130 - t159;
t116 = t146 * t173 + t234;
t115 = t145 * t171 + t231;
t114 = t130 - t129;
t111 = -t134 * t171 + t137 * t173;
t107 = (t227 * qJD(1) + t143) * t172 + t207;
t97 = -t159 - t129;
t95 = -t103 + t153;
t94 = t102 - t153;
t93 = -t129 - t130;
t90 = t113 + t244;
t89 = -t226 * t131 + t190;
t88 = -t186 - t262;
t87 = t185 + t262;
t86 = t226 * t133 + t185;
t82 = -t117 * t176 + t246;
t81 = t117 * t179 + t250;
t79 = t103 - t102;
t75 = t179 * t97 - t264;
t74 = t176 * t97 + t263;
t72 = (t104 * t178 - t106 * t175) * t156;
t71 = (t104 * t175 + t106 * t178) * t156;
t64 = -qJD(5) * t106 - t206;
t63 = -t102 - t103;
t62 = t176 * t91 + t179 * t88;
t61 = t176 * t88 - t179 * t91;
t60 = t178 * t94 + t252;
t59 = -t175 * t95 + t268;
t58 = t175 * t94 - t248;
t57 = t178 * t95 + t269;
t56 = -t175 * t92 + t248;
t53 = t171 * t82 + t173 * t81;
t52 = -t224 * t104 + t196;
t50 = t65 + t96;
t47 = t224 * t106 + t206;
t45 = t106 * t242 + t178 * t65;
t44 = -t106 * t241 + t175 * t65;
t43 = -t104 * t241 - t175 * t64;
t42 = -t104 * t242 + t178 * t64;
t41 = t178 * t73 - t269;
t39 = t171 * t75 + t173 * t74;
t35 = t171 * t62 + t173 * t61;
t34 = -t176 * t55 + t179 * t56;
t33 = t176 * t56 + t179 * t55;
t29 = -pkin(7) * t55 + t249;
t28 = t175 * t51 - t178 * t48;
t27 = -t175 * t50 - t178 * t47;
t25 = -t175 * t47 + t178 * t50;
t24 = -pkin(7) * t40 + t253;
t23 = -t176 * t40 + t179 * t41;
t22 = t176 * t41 + t179 * t40;
t19 = -pkin(4) * t52 + pkin(7) * t56 + t253;
t18 = -pkin(4) * t47 + pkin(7) * t41 - t249;
t17 = t171 * t34 + t173 * t33;
t13 = -t176 * t26 + t179 * t28;
t12 = t176 * t28 + t179 * t26;
t11 = t171 * t23 + t173 * t22;
t10 = t171 * t21 + t254;
t7 = t173 * t12 + t171 * t13;
t6 = -pkin(4) * t46 + pkin(7) * t9;
t5 = -pkin(7) * t26 - t8;
t4 = -pkin(4) * t63 + pkin(7) * t28 + t9;
t3 = t179 * t9 - t256;
t2 = t176 * t9 + t255;
t1 = t171 * t3 + t173 * t2;
t14 = [0, 0, 0, 0, 0, qJDD(1), t200, t199, 0, 0, t165, 0.2e1 * t172 * t219, 0, t166, 0, 0, -t261 * t174, t261 * t172, pkin(1) * t148 + qJ(2) * (t166 + t165) + t205, pkin(1) * t139 + qJ(2) * t205, -t197 + (qJDD(1) * t169 + t202) * t168, t172 * (-t135 * t173 - t136 * t171) + t174 * (t158 - t215), t172 * (t231 - (t170 - t238) * t233) + t174 * t137, t197 + (qJDD(1) * t167 - t202) * t168, t172 * (t173 * (t158 - t235) + t234) + t174 * t134, t166, t172 * (-qJ(3) * t115 + t107 * t171) + t174 * (-pkin(2) * t115 + t83) - pkin(1) * t115 + qJ(2) * (t174 * (-t141 * t171 + t145 * t173) + t172 * t135), t172 * (-qJ(3) * t116 + t107 * t173) + t174 * (-pkin(2) * t116 + t84) - pkin(1) * t116 + qJ(2) * (t174 * (t140 * t173 - t146 * t171) + t172 * t136), -t172 * t54 + qJ(2) * (t174 * (-t134 * t173 - t137 * t171) - t172 * t142) + t192 * t111, qJ(2) * (t174 * (t171 * t83 + t173 * t84) + t172 * t107) + t192 * t54, t172 * (t173 * (t113 * t179 + t133 * t240) - t171 * (t113 * t176 - t133 * t239)) - t210, t172 * (t173 * (-t176 * t90 - t179 * t86) - t171 * (-t176 * t86 + t179 * t90)) - t174 * t114, t172 * (t173 * (-t123 * t176 + t263) - t171 * (t123 * t179 + t264)) - t174 * t91, t172 * (t173 * (-t112 * t176 - t131 * t239) - t171 * (t112 * t179 - t131 * t240)) + t210, t172 * (t173 * (t122 * t179 + t250) - t171 * (t122 * t176 - t246)) + t174 * t87, t174 * t160 + t172 * (t173 * (t131 * t179 - t133 * t176) - t171 * (t131 * t176 + t133 * t179)) * t161, t172 * (t173 * (-pkin(6) * t74 + t251) - t171 * (-pkin(3) * t86 + pkin(6) * t75 - t247) - qJ(3) * t39) + t174 * (-pkin(2) * t39 - pkin(3) * t74 + t37) - pkin(1) * t39 + qJ(2) * (t174 * (-t171 * t74 + t173 * t75) + t172 * t86), t172 * (t173 * (-pkin(6) * t81 + t247) - t171 * (-pkin(3) * t89 + pkin(6) * t82 + t251) - qJ(3) * t53) + t174 * (-pkin(2) * t53 - pkin(3) * t81 + t38) - pkin(1) * t53 + qJ(2) * (t174 * (-t171 * t81 + t173 * t82) + t172 * t89), t172 * (t173 * (-pkin(6) * t61 - t20) - t171 * (-pkin(3) * t93 + pkin(6) * t62 + t21) - qJ(3) * t35) + t174 * (-pkin(2) * t35 - pkin(3) * t61) - pkin(1) * t35 + qJ(2) * (t174 * (-t171 * t61 + t173 * t62) + t172 * t93), t172 * (-pkin(6) * t254 - t171 * (-pkin(3) * t85 + pkin(6) * t21) - qJ(3) * t10) + t174 * (-pkin(2) * t10 - pkin(3) * t20) - pkin(1) * t10 + qJ(2) * (t174 * (-t171 * t20 + t173 * t21) + t172 * t85), t172 * (t173 * (-t176 * t44 + t179 * t45) - t171 * (t176 * t45 + t179 * t44)) - t211, t172 * (t173 * (-t176 * t25 + t179 * t27) - t171 * (t176 * t27 + t179 * t25)) - t174 * t79, t172 * (t173 * (-t176 * t57 + t179 * t59) - t171 * (t176 * t59 + t179 * t57)) - t174 * t51, t172 * (t173 * (-t176 * t42 + t179 * t43) - t171 * (t176 * t43 + t179 * t42)) + t211, t172 * (t173 * (-t176 * t58 + t179 * t60) - t171 * (t176 * t60 + t179 * t58)) + t174 * t48, t172 * (t173 * (-t176 * t71 + t179 * t72) - t171 * (t176 * t72 + t179 * t71)) + t174 * t154, t172 * (t173 * (-pkin(6) * t22 - t176 * t18 + t179 * t24) - t171 * (-pkin(3) * t47 + pkin(6) * t23 + t176 * t24 + t179 * t18) - qJ(3) * t11) + t174 * (-pkin(2) * t11 - pkin(3) * t22 - t194) - pkin(1) * t11 + qJ(2) * (t174 * (-t171 * t22 + t173 * t23) + t172 * t47), t172 * (t173 * (-pkin(6) * t33 - t176 * t19 + t179 * t29) - t171 * (-pkin(3) * t52 + pkin(6) * t34 + t176 * t29 + t179 * t19) - qJ(3) * t17) + t174 * (-pkin(2) * t17 - pkin(3) * t33 - t187) - pkin(1) * t17 + qJ(2) * (t174 * (-t171 * t33 + t173 * t34) + t172 * t52), t172 * (t173 * (-pkin(6) * t12 - t176 * t4 + t179 * t5) - t171 * (-pkin(3) * t63 + pkin(6) * t13 + t176 * t5 + t179 * t4) - qJ(3) * t7) + t174 * (-pkin(2) * t7 - pkin(3) * t12 - t258) - pkin(1) * t7 + qJ(2) * (t174 * (-t12 * t171 + t13 * t173) + t172 * t63), t172 * (t173 * (-pkin(6) * t2 - pkin(7) * t255 - t176 * t6) - t171 * (-pkin(3) * t46 + pkin(6) * t3 - pkin(7) * t256 + t179 * t6) - qJ(3) * t1) + t174 * (-pkin(2) * t1 - pkin(3) * t2 - t259) - pkin(1) * t1 + qJ(2) * (t174 * (-t171 * t2 + t173 * t3) + t172 * t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t219, t220, -t148, -t139, 0, 0, 0, 0, 0, 0, t115, t116, t111, t54, 0, 0, 0, 0, 0, 0, t39, t53, t35, t10, 0, 0, 0, 0, 0, 0, t11, t17, t7, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t136, -t142, t107, 0, 0, 0, 0, 0, 0, t86, t89, t93, t85, 0, 0, 0, 0, 0, 0, t47, t52, t63, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243, t114, t91, -t243, -t87, -t160, -t37, -t38, 0, 0, t80, t79, t51, -t80, -t48, -t154, t194, t187, t258, t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t79, t51, -t80, -t48, -t154, -t15, -t16, 0, 0;];
tauJ_reg = t14;
