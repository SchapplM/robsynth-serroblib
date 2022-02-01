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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:17:01
% EndTime: 2022-01-23 09:17:07
% DurationCPUTime: 3.93s
% Computational Cost: add. (14075->371), mult. (39105->547), div. (0->0), fcn. (27680->10), ass. (0->227)
t176 = sin(qJ(5));
t175 = cos(pkin(8));
t220 = t175 * qJDD(1);
t160 = -qJDD(4) + t220;
t154 = -qJDD(5) + t160;
t173 = sin(pkin(8));
t172 = sin(pkin(9));
t174 = cos(pkin(9));
t177 = sin(qJ(4));
t180 = cos(qJ(4));
t196 = t172 * t180 + t174 * t177;
t190 = t173 * t196;
t131 = qJD(1) * t190;
t233 = t173 * t174;
t214 = t180 * t233;
t230 = qJD(1) * t173;
t133 = -t177 * t172 * t230 + qJD(1) * t214;
t179 = cos(qJ(5));
t104 = t179 * t131 + t133 * t176;
t106 = -t131 * t176 + t133 * t179;
t80 = t106 * t104;
t266 = -t154 - t80;
t270 = t176 * t266;
t269 = t179 * t266;
t221 = t173 * qJDD(1);
t210 = t172 * t221;
t191 = qJDD(1) * t214 - t177 * t210;
t113 = -qJD(4) * t131 + t191;
t229 = qJD(1) * t175;
t161 = -qJD(4) + t229;
t245 = t131 * t161;
t91 = t113 - t245;
t182 = qJD(1) ^ 2;
t178 = sin(qJ(1));
t181 = cos(qJ(1));
t209 = g(1) * t178 - t181 * g(2);
t187 = -qJ(2) * t182 + qJDD(2) - t209;
t199 = -t175 * pkin(2) - t173 * qJ(3);
t193 = -pkin(1) + t199;
t268 = -0.2e1 * qJD(3) * t230 + qJDD(1) * t193 + t187;
t200 = g(1) * t181 + g(2) * t178;
t143 = -pkin(1) * t182 + qJDD(1) * qJ(2) - t200;
t261 = 2 * qJD(2);
t267 = qJD(1) * t261 + t143;
t244 = t133 * t131;
t189 = -t160 - t244;
t265 = t177 * t189;
t264 = t180 * t189;
t263 = (qJD(4) + t161) * t133;
t246 = qJDD(1) * pkin(1);
t139 = -t187 + t246;
t171 = t175 ^ 2;
t236 = t171 * t182;
t169 = t173 ^ 2;
t237 = t169 * t182;
t262 = -t139 - t246 + (t236 + t237) * qJ(2);
t156 = -qJD(5) + t161;
t186 = qJDD(1) * t190;
t112 = -qJD(4) * t133 - t186;
t206 = -t179 * t112 + t113 * t176;
t48 = (qJD(5) + t156) * t106 + t206;
t102 = t104 ^ 2;
t103 = t106 ^ 2;
t129 = t131 ^ 2;
t130 = t133 ^ 2;
t153 = t156 ^ 2;
t159 = t161 ^ 2;
t147 = t199 * qJD(1);
t201 = -g(3) * t173 + t267 * t175;
t109 = t147 * t229 + t201;
t192 = t268 * t174;
t194 = -t175 * pkin(3) - pkin(6) * t233;
t238 = t169 * t174;
t69 = t194 * qJDD(1) + (-t109 + (pkin(6) * t173 * t175 - pkin(3) * t238) * t182) * t172 + t192;
t144 = t194 * qJD(1);
t168 = t172 ^ 2;
t158 = t168 * t237;
t84 = t174 * t109 + t268 * t172;
t70 = -pkin(3) * t158 - pkin(6) * t210 + t144 * t229 + t84;
t37 = t177 * t70 - t180 * t69;
t31 = t189 * pkin(4) - pkin(7) * t91 - t37;
t121 = -pkin(4) * t161 - pkin(7) * t133;
t38 = t177 * t69 + t180 * t70;
t32 = -pkin(4) * t129 + pkin(7) * t112 + t121 * t161 + t38;
t15 = t176 * t32 - t179 * t31;
t16 = t176 * t31 + t179 * t32;
t8 = -t15 * t179 + t16 * t176;
t260 = pkin(4) * t8;
t197 = t112 * t176 + t113 * t179;
t65 = -qJD(5) * t104 + t197;
t96 = t104 * t156;
t51 = t65 - t96;
t26 = -t176 * t48 - t179 * t51;
t259 = pkin(4) * t26;
t258 = g(3) * t175;
t257 = t177 * t8;
t256 = t180 * t8;
t20 = t177 * t38 - t180 * t37;
t255 = t174 * t20;
t207 = qJDD(3) + t258;
t223 = qJDD(1) * t172;
t228 = t261 + t147;
t85 = -pkin(6) * t158 + (pkin(3) * t223 + t143 + (t144 * t174 + t228) * qJD(1)) * t173 + t207;
t46 = -t112 * pkin(4) - t129 * pkin(7) + t121 * t133 + t85;
t254 = t176 * t46;
t76 = t154 - t80;
t253 = t176 * t76;
t252 = t177 * t85;
t99 = t160 - t244;
t251 = t177 * t99;
t250 = t179 * t46;
t249 = t179 * t76;
t248 = t180 * t85;
t247 = t180 * t99;
t243 = t156 * t176;
t242 = t156 * t179;
t241 = t161 * t177;
t240 = t161 * t180;
t170 = t174 ^ 2;
t239 = t169 * t170;
t234 = t172 * t182;
t203 = t234 * t238;
t140 = -t203 + t220;
t235 = t172 * t140;
t141 = -t203 - t220;
t232 = t174 * t141;
t231 = t175 * t182;
t227 = qJD(4) - t161;
t225 = qJD(5) - t156;
t222 = qJDD(1) * t174;
t216 = t170 * t237;
t215 = t172 * t231;
t213 = t174 * t231;
t212 = t175 * t80;
t211 = t175 * t244;
t9 = t15 * t176 + t179 * t16;
t21 = t177 * t37 + t180 * t38;
t205 = t173 * (t267 * t173 + t258) + t175 * t201;
t202 = t172 * t213;
t83 = t109 * t172 - t192;
t54 = t172 * t84 - t174 * t83;
t198 = t169 * t202;
t73 = -t153 - t102;
t40 = t176 * t73 + t269;
t195 = pkin(4) * t40 - t15;
t92 = -t103 - t153;
t55 = t179 * t92 + t253;
t188 = pkin(4) * t55 - t16;
t185 = t196 * t221;
t166 = t171 * qJDD(1);
t165 = t169 * qJDD(1);
t148 = (t169 + t171) * t182;
t146 = (-t171 - t239) * t182;
t145 = -t158 - t236;
t142 = t158 + t216;
t137 = (t215 - t222) * t173;
t136 = (t215 + t222) * t173;
t135 = (-t213 + t223) * t173;
t134 = (t213 + t223) * t173;
t123 = -t130 + t159;
t122 = t129 - t159;
t117 = -t130 - t159;
t116 = t146 * t174 + t235;
t115 = t145 * t172 + t232;
t114 = t130 - t129;
t111 = -t134 * t172 + t137 * t174;
t107 = (t228 * qJD(1) + t143) * t173 + t207;
t97 = -t159 - t129;
t95 = -t103 + t153;
t94 = t102 - t153;
t93 = -t129 - t130;
t90 = t113 + t245;
t89 = -t227 * t131 + t191;
t88 = -t186 - t263;
t87 = t185 + t263;
t86 = t227 * t133 + t185;
t82 = -t117 * t177 + t247;
t81 = t117 * t180 + t251;
t79 = t103 - t102;
t75 = t180 * t97 - t265;
t74 = t177 * t97 + t264;
t72 = (t104 * t179 - t106 * t176) * t156;
t71 = (t104 * t176 + t106 * t179) * t156;
t64 = -qJD(5) * t106 - t206;
t63 = -t102 - t103;
t62 = t177 * t91 + t180 * t88;
t61 = t177 * t88 - t180 * t91;
t60 = t179 * t94 + t253;
t59 = -t176 * t95 + t269;
t58 = t176 * t94 - t249;
t57 = t179 * t95 + t270;
t56 = -t176 * t92 + t249;
t53 = t172 * t82 + t174 * t81;
t52 = -t225 * t104 + t197;
t50 = t65 + t96;
t47 = t225 * t106 + t206;
t45 = t106 * t243 + t179 * t65;
t44 = -t106 * t242 + t176 * t65;
t43 = -t104 * t242 - t176 * t64;
t42 = -t104 * t243 + t179 * t64;
t41 = t179 * t73 - t270;
t39 = t172 * t75 + t174 * t74;
t35 = t172 * t62 + t174 * t61;
t34 = -t177 * t55 + t180 * t56;
t33 = t177 * t56 + t180 * t55;
t29 = -pkin(7) * t55 + t250;
t28 = t176 * t51 - t179 * t48;
t27 = -t176 * t50 - t179 * t47;
t25 = -t176 * t47 + t179 * t50;
t24 = -pkin(7) * t40 + t254;
t23 = -t177 * t40 + t180 * t41;
t22 = t177 * t41 + t180 * t40;
t19 = -pkin(4) * t52 + pkin(7) * t56 + t254;
t18 = -pkin(4) * t47 + pkin(7) * t41 - t250;
t17 = t172 * t34 + t174 * t33;
t13 = -t177 * t26 + t180 * t28;
t12 = t177 * t28 + t180 * t26;
t11 = t172 * t23 + t174 * t22;
t10 = t172 * t21 + t255;
t7 = t174 * t12 + t172 * t13;
t6 = -pkin(4) * t46 + pkin(7) * t9;
t5 = -pkin(7) * t26 - t8;
t4 = -pkin(4) * t63 + pkin(7) * t28 + t9;
t3 = t180 * t9 - t257;
t2 = t177 * t9 + t256;
t1 = t172 * t3 + t174 * t2;
t14 = [0, 0, 0, 0, 0, qJDD(1), t209, t200, 0, 0, t165, 0.2e1 * t173 * t220, 0, t166, 0, 0, -t262 * t175, t262 * t173, pkin(1) * t148 + qJ(2) * (t166 + t165) + t205, pkin(1) * t139 + qJ(2) * t205, -t198 + (qJDD(1) * t170 + t202) * t169, t173 * (-t135 * t174 - t136 * t172) + t175 * (t158 - t216), t173 * (t232 - (t171 - t239) * t234) + t175 * t137, t198 + (qJDD(1) * t168 - t202) * t169, t173 * (t174 * (t158 - t236) + t235) + t175 * t134, t166, t173 * (-qJ(3) * t115 + t107 * t172) + t175 * (-pkin(2) * t115 + t83) - pkin(1) * t115 + qJ(2) * (t175 * (-t141 * t172 + t145 * t174) + t173 * t135), t173 * (-qJ(3) * t116 + t107 * t174) + t175 * (-pkin(2) * t116 + t84) - pkin(1) * t116 + qJ(2) * (t175 * (t140 * t174 - t146 * t172) + t173 * t136), -t173 * t54 + qJ(2) * (t175 * (-t134 * t174 - t137 * t172) - t173 * t142) + t193 * t111, qJ(2) * (t175 * (t172 * t83 + t174 * t84) + t173 * t107) + t193 * t54, t173 * (t174 * (t113 * t180 + t133 * t241) - t172 * (t113 * t177 - t133 * t240)) - t211, t173 * (t174 * (-t177 * t90 - t180 * t86) - t172 * (-t177 * t86 + t180 * t90)) - t175 * t114, t173 * (t174 * (-t123 * t177 + t264) - t172 * (t123 * t180 + t265)) - t175 * t91, t173 * (t174 * (-t112 * t177 - t131 * t240) - t172 * (t112 * t180 - t131 * t241)) + t211, t173 * (t174 * (t122 * t180 + t251) - t172 * (t122 * t177 - t247)) + t175 * t87, t175 * t160 + t173 * (t174 * (t131 * t180 - t133 * t177) - t172 * (t131 * t177 + t133 * t180)) * t161, t173 * (t174 * (-pkin(6) * t74 + t252) - t172 * (-pkin(3) * t86 + pkin(6) * t75 - t248) - qJ(3) * t39) + t175 * (-pkin(2) * t39 - pkin(3) * t74 + t37) - pkin(1) * t39 + qJ(2) * (t175 * (-t172 * t74 + t174 * t75) + t173 * t86), t173 * (t174 * (-pkin(6) * t81 + t248) - t172 * (-pkin(3) * t89 + pkin(6) * t82 + t252) - qJ(3) * t53) + t175 * (-pkin(2) * t53 - pkin(3) * t81 + t38) - pkin(1) * t53 + qJ(2) * (t175 * (-t172 * t81 + t174 * t82) + t173 * t89), t173 * (t174 * (-pkin(6) * t61 - t20) - t172 * (-pkin(3) * t93 + pkin(6) * t62 + t21) - qJ(3) * t35) + t175 * (-pkin(2) * t35 - pkin(3) * t61) - pkin(1) * t35 + qJ(2) * (t175 * (-t172 * t61 + t174 * t62) + t173 * t93), t173 * (-pkin(6) * t255 - t172 * (-pkin(3) * t85 + pkin(6) * t21) - qJ(3) * t10) + t175 * (-pkin(2) * t10 - pkin(3) * t20) - pkin(1) * t10 + qJ(2) * (t175 * (-t172 * t20 + t174 * t21) + t173 * t85), t173 * (t174 * (-t177 * t44 + t180 * t45) - t172 * (t177 * t45 + t180 * t44)) - t212, t173 * (t174 * (-t177 * t25 + t180 * t27) - t172 * (t177 * t27 + t180 * t25)) - t175 * t79, t173 * (t174 * (-t177 * t57 + t180 * t59) - t172 * (t177 * t59 + t180 * t57)) - t175 * t51, t173 * (t174 * (-t177 * t42 + t180 * t43) - t172 * (t177 * t43 + t180 * t42)) + t212, t173 * (t174 * (-t177 * t58 + t180 * t60) - t172 * (t177 * t60 + t180 * t58)) + t175 * t48, t173 * (t174 * (-t177 * t71 + t180 * t72) - t172 * (t177 * t72 + t180 * t71)) + t175 * t154, t173 * (t174 * (-pkin(6) * t22 - t177 * t18 + t180 * t24) - t172 * (-pkin(3) * t47 + pkin(6) * t23 + t177 * t24 + t18 * t180) - qJ(3) * t11) + t175 * (-pkin(2) * t11 - pkin(3) * t22 - t195) - pkin(1) * t11 + qJ(2) * (t175 * (-t172 * t22 + t174 * t23) + t173 * t47), t173 * (t174 * (-pkin(6) * t33 - t177 * t19 + t180 * t29) - t172 * (-pkin(3) * t52 + pkin(6) * t34 + t177 * t29 + t180 * t19) - qJ(3) * t17) + t175 * (-pkin(2) * t17 - pkin(3) * t33 - t188) - pkin(1) * t17 + qJ(2) * (t175 * (-t172 * t33 + t174 * t34) + t173 * t52), t173 * (t174 * (-pkin(6) * t12 - t177 * t4 + t180 * t5) - t172 * (-pkin(3) * t63 + pkin(6) * t13 + t177 * t5 + t180 * t4) - qJ(3) * t7) + t175 * (-pkin(2) * t7 - pkin(3) * t12 - t259) - pkin(1) * t7 + qJ(2) * (t175 * (-t12 * t172 + t13 * t174) + t173 * t63), t173 * (t174 * (-pkin(6) * t2 - pkin(7) * t256 - t177 * t6) - t172 * (-pkin(3) * t46 + pkin(6) * t3 - pkin(7) * t257 + t180 * t6) - qJ(3) * t1) + t175 * (-pkin(2) * t1 - pkin(3) * t2 - t260) - pkin(1) * t1 + qJ(2) * (t175 * (-t172 * t2 + t174 * t3) + t173 * t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220, t221, -t148, -t139, 0, 0, 0, 0, 0, 0, t115, t116, t111, t54, 0, 0, 0, 0, 0, 0, t39, t53, t35, t10, 0, 0, 0, 0, 0, 0, t11, t17, t7, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t136, -t142, t107, 0, 0, 0, 0, 0, 0, t86, t89, t93, t85, 0, 0, 0, 0, 0, 0, t47, t52, t63, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, t114, t91, -t244, -t87, -t160, -t37, -t38, 0, 0, t80, t79, t51, -t80, -t48, -t154, t195, t188, t259, t260; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t79, t51, -t80, -t48, -t154, -t15, -t16, 0, 0;];
tauJ_reg = t14;
