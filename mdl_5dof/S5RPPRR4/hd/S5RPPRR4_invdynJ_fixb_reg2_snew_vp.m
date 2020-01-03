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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:31:26
% EndTime: 2020-01-03 11:31:44
% DurationCPUTime: 4.89s
% Computational Cost: add. (14075->367), mult. (39105->547), div. (0->0), fcn. (27680->10), ass. (0->228)
t174 = sin(pkin(8));
t168 = t174 ^ 2;
t176 = cos(pkin(8));
t170 = t176 ^ 2;
t270 = t168 + t170;
t177 = sin(qJ(5));
t220 = t176 * qJDD(1);
t160 = -qJDD(4) + t220;
t154 = -qJDD(5) + t160;
t173 = sin(pkin(9));
t175 = cos(pkin(9));
t178 = sin(qJ(4));
t181 = cos(qJ(4));
t196 = t173 * t181 + t175 * t178;
t190 = t174 * t196;
t131 = qJD(1) * t190;
t233 = t174 * t175;
t214 = t181 * t233;
t230 = qJD(1) * t174;
t133 = -t178 * t173 * t230 + qJD(1) * t214;
t180 = cos(qJ(5));
t104 = t180 * t131 + t133 * t177;
t106 = -t131 * t177 + t133 * t180;
t80 = t106 * t104;
t265 = -t154 - t80;
t269 = t177 * t265;
t268 = t180 * t265;
t221 = t174 * qJDD(1);
t210 = t173 * t221;
t191 = qJDD(1) * t214 - t178 * t210;
t113 = -qJD(4) * t131 + t191;
t229 = qJD(1) * t176;
t161 = -qJD(4) + t229;
t245 = t131 * t161;
t91 = t113 - t245;
t183 = qJD(1) ^ 2;
t171 = t183 * qJ(2);
t172 = qJDD(1) * pkin(1);
t179 = sin(qJ(1));
t182 = cos(qJ(1));
t201 = -g(2) * t182 - g(3) * t179;
t139 = qJDD(2) - t171 - t172 - t201;
t199 = -t176 * pkin(2) - t174 * qJ(3);
t267 = -0.2e1 * qJD(3) * t230 + t199 * qJDD(1) + t139;
t200 = g(2) * t179 - g(3) * t182;
t143 = -pkin(1) * t183 + qJDD(1) * qJ(2) - t200;
t260 = 2 * qJD(2);
t266 = qJD(1) * t260 + t143;
t244 = t133 * t131;
t189 = -t160 - t244;
t264 = t178 * t189;
t263 = t181 * t189;
t262 = (qJD(4) + t161) * t133;
t261 = t270 * t171 + t139 - t172;
t156 = -qJD(5) + t161;
t187 = qJDD(1) * t190;
t112 = -qJD(4) * t133 - t187;
t207 = -t180 * t112 + t113 * t177;
t48 = (qJD(5) + t156) * t106 + t207;
t102 = t104 ^ 2;
t103 = t106 ^ 2;
t129 = t131 ^ 2;
t130 = t133 ^ 2;
t153 = t156 ^ 2;
t159 = t161 ^ 2;
t147 = t199 * qJD(1);
t202 = -g(1) * t174 + t266 * t176;
t109 = t147 * t229 + t202;
t192 = t267 * t175;
t194 = -t176 * pkin(3) - pkin(6) * t233;
t238 = t168 * t175;
t69 = t194 * qJDD(1) + (-t109 + (pkin(6) * t174 * t176 - pkin(3) * t238) * t183) * t173 + t192;
t144 = t194 * qJD(1);
t167 = t173 ^ 2;
t237 = t168 * t183;
t158 = t167 * t237;
t84 = t175 * t109 + t267 * t173;
t70 = -pkin(3) * t158 - pkin(6) * t210 + t144 * t229 + t84;
t37 = t178 * t70 - t181 * t69;
t31 = t189 * pkin(4) - t91 * pkin(7) - t37;
t121 = -pkin(4) * t161 - pkin(7) * t133;
t38 = t178 * t69 + t181 * t70;
t32 = -pkin(4) * t129 + pkin(7) * t112 + t121 * t161 + t38;
t15 = t177 * t32 - t180 * t31;
t16 = t177 * t31 + t180 * t32;
t8 = -t15 * t180 + t16 * t177;
t259 = pkin(4) * t8;
t197 = t112 * t177 + t113 * t180;
t65 = -qJD(5) * t104 + t197;
t96 = t104 * t156;
t51 = t65 - t96;
t26 = -t177 * t48 - t180 * t51;
t258 = pkin(4) * t26;
t257 = g(1) * t176;
t256 = t178 * t8;
t255 = t181 * t8;
t20 = t178 * t38 - t181 * t37;
t254 = t175 * t20;
t208 = qJDD(3) + t257;
t223 = qJDD(1) * t173;
t228 = t260 + t147;
t85 = -pkin(6) * t158 + (pkin(3) * t223 + t143 + (t144 * t175 + t228) * qJD(1)) * t174 + t208;
t46 = -t112 * pkin(4) - t129 * pkin(7) + t121 * t133 + t85;
t253 = t177 * t46;
t76 = t154 - t80;
t252 = t177 * t76;
t251 = t178 * t85;
t99 = t160 - t244;
t250 = t178 * t99;
t249 = t180 * t46;
t248 = t180 * t76;
t247 = t181 * t85;
t246 = t181 * t99;
t243 = t156 * t177;
t242 = t156 * t180;
t241 = t161 * t178;
t240 = t161 * t181;
t169 = t175 ^ 2;
t239 = t168 * t169;
t236 = t170 * t183;
t234 = t173 * t183;
t204 = t234 * t238;
t140 = -t204 + t220;
t235 = t173 * t140;
t141 = -t204 - t220;
t232 = t175 * t141;
t231 = t176 * t183;
t227 = qJD(4) - t161;
t225 = qJD(5) - t156;
t222 = qJDD(1) * t175;
t216 = t169 * t237;
t215 = t173 * t231;
t213 = t175 * t231;
t212 = t176 * t80;
t211 = t176 * t244;
t9 = t15 * t177 + t180 * t16;
t21 = t178 * t37 + t181 * t38;
t206 = t174 * (t266 * t174 + t257) + t176 * t202;
t203 = t173 * t213;
t83 = t109 * t173 - t192;
t54 = t173 * t84 - t175 * t83;
t198 = t168 * t203;
t73 = -t153 - t102;
t40 = t177 * t73 + t268;
t195 = pkin(4) * t40 - t15;
t193 = -pkin(1) + t199;
t92 = -t103 - t153;
t55 = t180 * t92 + t252;
t188 = pkin(4) * t55 - t16;
t186 = t196 * t221;
t166 = t170 * qJDD(1);
t165 = t168 * qJDD(1);
t148 = t270 * t183;
t146 = (-t170 - t239) * t183;
t145 = -t158 - t236;
t142 = t158 + t216;
t137 = (t215 - t222) * t174;
t136 = (t215 + t222) * t174;
t135 = (-t213 + t223) * t174;
t134 = (t213 + t223) * t174;
t123 = -t130 + t159;
t122 = t129 - t159;
t117 = -t130 - t159;
t116 = t146 * t175 + t235;
t115 = t145 * t173 + t232;
t114 = t130 - t129;
t111 = -t134 * t173 + t137 * t175;
t107 = (t228 * qJD(1) + t143) * t174 + t208;
t97 = -t159 - t129;
t95 = -t103 + t153;
t94 = t102 - t153;
t93 = -t129 - t130;
t90 = t113 + t245;
t89 = -t227 * t131 + t191;
t88 = -t187 - t262;
t87 = t186 + t262;
t86 = t227 * t133 + t186;
t82 = -t117 * t178 + t246;
t81 = t117 * t181 + t250;
t79 = t103 - t102;
t75 = t181 * t97 - t264;
t74 = t178 * t97 + t263;
t72 = (t104 * t180 - t106 * t177) * t156;
t71 = (t104 * t177 + t106 * t180) * t156;
t64 = -qJD(5) * t106 - t207;
t63 = -t102 - t103;
t62 = t178 * t91 + t181 * t88;
t61 = t178 * t88 - t181 * t91;
t60 = t180 * t94 + t252;
t59 = -t177 * t95 + t268;
t58 = t177 * t94 - t248;
t57 = t180 * t95 + t269;
t56 = -t177 * t92 + t248;
t53 = t173 * t82 + t175 * t81;
t52 = -t225 * t104 + t197;
t50 = t65 + t96;
t47 = t225 * t106 + t207;
t45 = t106 * t243 + t180 * t65;
t44 = -t106 * t242 + t177 * t65;
t43 = -t104 * t242 - t177 * t64;
t42 = -t104 * t243 + t180 * t64;
t41 = t180 * t73 - t269;
t39 = t173 * t75 + t175 * t74;
t35 = t173 * t62 + t175 * t61;
t34 = -t178 * t55 + t181 * t56;
t33 = t178 * t56 + t181 * t55;
t29 = -pkin(7) * t55 + t249;
t28 = t177 * t51 - t180 * t48;
t27 = -t177 * t50 - t180 * t47;
t25 = -t177 * t47 + t180 * t50;
t24 = -pkin(7) * t40 + t253;
t23 = -t178 * t40 + t181 * t41;
t22 = t178 * t41 + t181 * t40;
t19 = -pkin(4) * t52 + pkin(7) * t56 + t253;
t18 = -pkin(4) * t47 + pkin(7) * t41 - t249;
t17 = t173 * t34 + t175 * t33;
t13 = -t178 * t26 + t181 * t28;
t12 = t178 * t28 + t181 * t26;
t11 = t173 * t23 + t175 * t22;
t10 = t173 * t21 + t254;
t7 = t175 * t12 + t173 * t13;
t6 = -pkin(4) * t46 + pkin(7) * t9;
t5 = -pkin(7) * t26 - t8;
t4 = -pkin(4) * t63 + pkin(7) * t28 + t9;
t3 = t181 * t9 - t256;
t2 = t178 * t9 + t255;
t1 = t173 * t3 + t175 * t2;
t14 = [0, 0, 0, 0, 0, qJDD(1), t201, t200, 0, 0, t165, 0.2e1 * t174 * t220, 0, t166, 0, 0, -t261 * t176, t261 * t174, pkin(1) * t148 + qJ(2) * (t166 + t165) + t206, -pkin(1) * t139 + qJ(2) * t206, -t198 + (qJDD(1) * t169 + t203) * t168, t174 * (-t135 * t175 - t136 * t173) + t176 * (t158 - t216), t174 * (t232 - (t170 - t239) * t234) + t176 * t137, t198 + (qJDD(1) * t167 - t203) * t168, t174 * (t175 * (t158 - t236) + t235) + t176 * t134, t166, t174 * (-qJ(3) * t115 + t107 * t173) + t176 * (-pkin(2) * t115 + t83) - pkin(1) * t115 + qJ(2) * (t176 * (-t141 * t173 + t145 * t175) + t174 * t135), t174 * (-qJ(3) * t116 + t107 * t175) + t176 * (-pkin(2) * t116 + t84) - pkin(1) * t116 + qJ(2) * (t176 * (t140 * t175 - t146 * t173) + t174 * t136), -t174 * t54 + qJ(2) * (t176 * (-t134 * t175 - t137 * t173) - t174 * t142) + t193 * t111, qJ(2) * (t176 * (t173 * t83 + t175 * t84) + t174 * t107) + t193 * t54, t174 * (t175 * (t113 * t181 + t133 * t241) - t173 * (t113 * t178 - t133 * t240)) - t211, t174 * (t175 * (-t178 * t90 - t181 * t86) - t173 * (-t178 * t86 + t181 * t90)) - t176 * t114, t174 * (t175 * (-t123 * t178 + t263) - t173 * (t123 * t181 + t264)) - t176 * t91, t174 * (t175 * (-t112 * t178 - t131 * t240) - t173 * (t112 * t181 - t131 * t241)) + t211, t174 * (t175 * (t122 * t181 + t250) - t173 * (t122 * t178 - t246)) + t176 * t87, t176 * t160 + t174 * (t175 * (t131 * t181 - t133 * t178) - t173 * (t131 * t178 + t133 * t181)) * t161, t174 * (t175 * (-pkin(6) * t74 + t251) - t173 * (-pkin(3) * t86 + pkin(6) * t75 - t247) - qJ(3) * t39) + t176 * (-pkin(2) * t39 - pkin(3) * t74 + t37) - pkin(1) * t39 + qJ(2) * (t176 * (-t173 * t74 + t175 * t75) + t174 * t86), t174 * (t175 * (-pkin(6) * t81 + t247) - t173 * (-pkin(3) * t89 + pkin(6) * t82 + t251) - qJ(3) * t53) + t176 * (-pkin(2) * t53 - pkin(3) * t81 + t38) - pkin(1) * t53 + qJ(2) * (t176 * (-t173 * t81 + t175 * t82) + t174 * t89), t174 * (t175 * (-pkin(6) * t61 - t20) - t173 * (-pkin(3) * t93 + pkin(6) * t62 + t21) - qJ(3) * t35) + t176 * (-pkin(2) * t35 - pkin(3) * t61) - pkin(1) * t35 + qJ(2) * (t176 * (-t173 * t61 + t175 * t62) + t174 * t93), t174 * (-pkin(6) * t254 - t173 * (-pkin(3) * t85 + pkin(6) * t21) - qJ(3) * t10) + t176 * (-pkin(2) * t10 - pkin(3) * t20) - pkin(1) * t10 + qJ(2) * (t176 * (-t173 * t20 + t175 * t21) + t174 * t85), t174 * (t175 * (-t178 * t44 + t181 * t45) - t173 * (t178 * t45 + t181 * t44)) - t212, t174 * (t175 * (-t178 * t25 + t181 * t27) - t173 * (t178 * t27 + t181 * t25)) - t176 * t79, t174 * (t175 * (-t178 * t57 + t181 * t59) - t173 * (t178 * t59 + t181 * t57)) - t176 * t51, t174 * (t175 * (-t178 * t42 + t181 * t43) - t173 * (t178 * t43 + t181 * t42)) + t212, t174 * (t175 * (-t178 * t58 + t181 * t60) - t173 * (t178 * t60 + t181 * t58)) + t176 * t48, t174 * (t175 * (-t178 * t71 + t181 * t72) - t173 * (t178 * t72 + t181 * t71)) + t176 * t154, t174 * (t175 * (-pkin(6) * t22 - t178 * t18 + t181 * t24) - t173 * (-pkin(3) * t47 + pkin(6) * t23 + t178 * t24 + t18 * t181) - qJ(3) * t11) + t176 * (-pkin(2) * t11 - pkin(3) * t22 - t195) - pkin(1) * t11 + qJ(2) * (t176 * (-t173 * t22 + t175 * t23) + t174 * t47), t174 * (t175 * (-pkin(6) * t33 - t178 * t19 + t181 * t29) - t173 * (-pkin(3) * t52 + pkin(6) * t34 + t178 * t29 + t181 * t19) - qJ(3) * t17) + t176 * (-pkin(2) * t17 - pkin(3) * t33 - t188) - pkin(1) * t17 + qJ(2) * (t176 * (-t173 * t33 + t175 * t34) + t174 * t52), t174 * (t175 * (-pkin(6) * t12 - t178 * t4 + t181 * t5) - t173 * (-pkin(3) * t63 + pkin(6) * t13 + t178 * t5 + t181 * t4) - qJ(3) * t7) + t176 * (-pkin(2) * t7 - pkin(3) * t12 - t258) - pkin(1) * t7 + qJ(2) * (t176 * (-t12 * t173 + t13 * t175) + t174 * t63), t174 * (t175 * (-pkin(6) * t2 - pkin(7) * t255 - t178 * t6) - t173 * (-pkin(3) * t46 + pkin(6) * t3 - pkin(7) * t256 + t181 * t6) - qJ(3) * t1) + t176 * (-pkin(2) * t1 - pkin(3) * t2 - t259) - pkin(1) * t1 + qJ(2) * (t176 * (-t173 * t2 + t175 * t3) + t174 * t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220, t221, -t148, t139, 0, 0, 0, 0, 0, 0, t115, t116, t111, t54, 0, 0, 0, 0, 0, 0, t39, t53, t35, t10, 0, 0, 0, 0, 0, 0, t11, t17, t7, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t136, -t142, t107, 0, 0, 0, 0, 0, 0, t86, t89, t93, t85, 0, 0, 0, 0, 0, 0, t47, t52, t63, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, t114, t91, -t244, -t87, -t160, -t37, -t38, 0, 0, t80, t79, t51, -t80, -t48, -t154, t195, t188, t258, t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t79, t51, -t80, -t48, -t154, -t15, -t16, 0, 0;];
tauJ_reg = t14;
