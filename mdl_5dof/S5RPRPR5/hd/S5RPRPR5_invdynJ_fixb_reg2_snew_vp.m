% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:43:06
% EndTime: 2020-01-03 11:43:21
% DurationCPUTime: 4.84s
% Computational Cost: add. (16064->372), mult. (41833->537), div. (0->0), fcn. (29414->10), ass. (0->219)
t180 = sin(pkin(8));
t174 = t180 ^ 2;
t182 = cos(pkin(8));
t175 = t182 ^ 2;
t270 = t174 + t175;
t179 = sin(pkin(9));
t181 = cos(pkin(9));
t187 = cos(qJ(3));
t229 = qJD(1) * t187;
t211 = t180 * t229;
t184 = sin(qJ(3));
t230 = qJD(1) * t180;
t212 = t184 * t230;
t139 = -t179 * t212 + t181 * t211;
t200 = -t182 * pkin(2) - t180 * pkin(6);
t154 = t200 * qJD(1);
t189 = qJD(1) ^ 2;
t185 = sin(qJ(1));
t188 = cos(qJ(1));
t198 = t185 * g(2) - t188 * g(3);
t151 = -t189 * pkin(1) + qJDD(1) * qJ(2) - t198;
t257 = 2 * qJD(2);
t265 = qJD(1) * t257 + t151;
t202 = -t180 * g(1) + t265 * t182;
t226 = t182 * qJD(1);
t109 = t154 * t226 + t202;
t176 = t189 * qJ(2);
t178 = qJDD(1) * pkin(1);
t199 = -t188 * g(2) - t185 * g(3);
t148 = qJDD(2) - t176 - t178 - t199;
t129 = t200 * qJDD(1) + t148;
t121 = t187 * t129;
t221 = t180 * qJDD(1);
t164 = t187 * t221;
t147 = -qJD(3) * t212 + t164;
t220 = t182 * qJDD(1);
t167 = -qJDD(3) + t220;
t168 = -qJD(3) + t226;
t233 = t174 * t189;
t216 = t187 * t233;
t69 = -t167 * pkin(3) - t147 * qJ(4) + t121 + (qJ(4) * t168 * t230 - pkin(3) * t216 - t109) * t184;
t143 = -t168 * pkin(3) - qJ(4) * t211;
t146 = (qJD(3) * t229 + qJDD(1) * t184) * t180;
t166 = t184 ^ 2 * t233;
t85 = t187 * t109 + t184 * t129;
t70 = -pkin(3) * t166 - t146 * qJ(4) + t168 * t143 + t85;
t201 = -0.2e1 * qJD(4) * t139 - t179 * t70 + t181 * t69;
t137 = (-t187 * t179 - t184 * t181) * t230;
t239 = t139 * t137;
t194 = -t167 + t239;
t261 = t181 * t194;
t135 = t137 ^ 2;
t258 = t168 ^ 2;
t98 = -t258 - t135;
t73 = t179 * t98 + t261;
t269 = pkin(3) * t73 + t201;
t183 = sin(qJ(5));
t160 = -qJDD(5) + t167;
t186 = cos(qJ(5));
t105 = -t186 * t137 + t183 * t139;
t107 = t183 * t137 + t186 * t139;
t81 = t107 * t105;
t264 = -t160 - t81;
t268 = t183 * t264;
t267 = t186 * t264;
t115 = -t179 * t146 + t181 * t147;
t240 = t137 * t168;
t266 = t115 - t240;
t91 = t115 + t240;
t263 = t179 * t194;
t153 = t168 * t211;
t124 = t153 - t146;
t262 = t180 * t124;
t260 = t270 * t176 + t148 - t178;
t162 = -qJD(5) + t168;
t206 = t181 * t146 + t179 * t147;
t208 = t183 * t115 + t186 * t206;
t48 = (qJD(5) + t162) * t107 + t208;
t259 = (qJD(3) - t168) * t212 - t164;
t102 = t105 ^ 2;
t103 = t107 ^ 2;
t136 = t139 ^ 2;
t159 = t162 ^ 2;
t228 = qJD(4) * t137;
t131 = 0.2e1 * t228;
t250 = t179 * t69 + t181 * t70;
t38 = t131 + t250;
t20 = t179 * t38 + t181 * t201;
t256 = pkin(3) * t20;
t236 = t168 * t139;
t87 = t206 + t236;
t62 = -t179 * t87 - t181 * t91;
t255 = pkin(3) * t62;
t26 = t194 * pkin(4) - t91 * pkin(7) + t201;
t119 = -t168 * pkin(4) - t139 * pkin(7);
t32 = -t135 * pkin(4) - pkin(7) * t206 + t168 * t119 + t38;
t13 = t183 * t32 - t186 * t26;
t14 = t183 * t26 + t186 * t32;
t7 = -t186 * t13 + t183 * t14;
t253 = t179 * t7;
t252 = t181 * t7;
t251 = t182 * g(1);
t225 = t257 + t154;
t80 = t251 + qJDD(4) + t146 * pkin(3) - qJ(4) * t166 + (t151 + (t143 * t187 + t225) * qJD(1)) * t180;
t249 = t179 * t80;
t99 = t167 + t239;
t248 = t179 * t99;
t247 = t181 * t80;
t246 = t181 * t99;
t42 = pkin(4) * t206 - t135 * pkin(7) + t139 * t119 + t80;
t245 = t183 * t42;
t76 = t160 - t81;
t244 = t183 * t76;
t243 = t186 * t42;
t242 = t186 * t76;
t241 = t187 * t20;
t238 = t162 * t183;
t237 = t162 * t186;
t235 = t168 * t179;
t234 = t168 * t181;
t158 = t184 * t216;
t144 = -t158 + t167;
t232 = t184 * t144;
t145 = -t158 - t167;
t231 = t187 * t145;
t223 = qJD(5) - t162;
t215 = t187 ^ 2 * t233;
t214 = t182 * t81;
t213 = t182 * t239;
t8 = t183 * t13 + t186 * t14;
t21 = -t179 * t201 + t181 * t38;
t84 = t184 * t109 - t121;
t207 = t180 * (t265 * t180 + t251) + t182 * t202;
t117 = -t136 - t258;
t82 = t181 * t117 + t248;
t205 = pkin(3) * t82 - t250;
t2 = t179 * t8 + t252;
t204 = pkin(3) * t2 + pkin(4) * t7;
t197 = t186 * t115 - t183 * t206;
t65 = -t105 * qJD(5) + t197;
t97 = t105 * t162;
t51 = t65 - t97;
t29 = -t183 * t48 - t186 * t51;
t31 = t183 * t51 - t186 * t48;
t15 = t179 * t31 + t181 * t29;
t203 = pkin(3) * t15 + pkin(4) * t29;
t60 = t184 * t85 - t187 * t84;
t196 = -pkin(1) + t200;
t75 = -t159 - t102;
t40 = t183 * t75 + t267;
t41 = t186 * t75 - t268;
t22 = t179 * t41 + t181 * t40;
t193 = pkin(3) * t22 + pkin(4) * t40 - t13;
t93 = -t103 - t159;
t54 = t186 * t93 + t244;
t55 = -t183 * t93 + t242;
t33 = t179 * t55 + t181 * t54;
t192 = pkin(3) * t33 + pkin(4) * t54 - t14;
t173 = t175 * qJDD(1);
t172 = t174 * qJDD(1);
t156 = t270 * t189;
t155 = t182 * t167;
t150 = -t166 + t215;
t149 = -t166 - t258;
t134 = -t258 - t215;
t126 = -t136 + t258;
t125 = t135 - t258;
t123 = t153 + t146;
t122 = -t164 + (qJD(3) + t168) * t212;
t116 = t184 * t149 + t231;
t112 = t136 - t135;
t110 = t187 * t134 + t232;
t108 = t251 + (t225 * qJD(1) + t151) * t180;
t96 = -t103 + t159;
t95 = t102 - t159;
t94 = -t135 - t136;
t92 = t187 * t122 - t184 * t123;
t86 = t206 - t236;
t83 = -t179 * t117 + t246;
t79 = t103 - t102;
t74 = t181 * t98 - t263;
t72 = (t105 * t186 - t107 * t183) * t162;
t71 = (t105 * t183 + t107 * t186) * t162;
t64 = -t107 * qJD(5) - t208;
t63 = t179 * t91 - t181 * t87;
t61 = -t102 - t103;
t59 = t186 * t95 + t244;
t58 = -t183 * t96 + t267;
t57 = t183 * t95 - t242;
t56 = t186 * t96 + t268;
t53 = t184 * t83 + t187 * t82;
t52 = -t223 * t105 + t197;
t50 = t65 + t97;
t47 = t223 * t107 + t208;
t46 = t107 * t238 + t186 * t65;
t45 = -t107 * t237 + t183 * t65;
t44 = -t105 * t237 - t183 * t64;
t43 = -t105 * t238 + t186 * t64;
t39 = t184 * t74 + t187 * t73;
t36 = t184 * t63 + t187 * t62;
t34 = -t179 * t54 + t181 * t55;
t30 = -t183 * t50 - t186 * t47;
t28 = -t183 * t47 + t186 * t50;
t27 = -pkin(7) * t54 + t243;
t24 = -pkin(7) * t40 + t245;
t23 = -t179 * t40 + t181 * t41;
t19 = -pkin(4) * t52 + pkin(7) * t55 + t245;
t18 = -pkin(4) * t47 + pkin(7) * t41 - t243;
t17 = t184 * t34 + t187 * t33;
t16 = -t179 * t29 + t181 * t31;
t11 = t184 * t23 + t187 * t22;
t10 = t184 * t21 + t241;
t9 = t187 * t15 + t184 * t16;
t6 = -pkin(4) * t42 + pkin(7) * t8;
t5 = -pkin(7) * t29 - t7;
t4 = -pkin(4) * t61 + pkin(7) * t31 + t8;
t3 = t181 * t8 - t253;
t1 = t184 * t3 + t187 * t2;
t12 = [0, 0, 0, 0, 0, qJDD(1), t199, t198, 0, 0, t172, 0.2e1 * t180 * t220, 0, t173, 0, 0, -t260 * t182, t260 * t180, pkin(1) * t156 + qJ(2) * (t173 + t172) + t207, -pkin(1) * t148 + qJ(2) * t207, (t180 * (t168 * t212 + t147) - t182 * t184 * t233) * t187, t180 * (t187 * t124 + t259 * t184) - t182 * t150, t180 * (t231 - t184 * (t258 - t215)) + t182 * t122, (t182 * t216 - t262) * t184, t180 * (t187 * (t166 - t258) + t232) + t182 * t123, t155, t180 * (-pkin(6) * t116 + t184 * t108) + t182 * (-pkin(2) * t116 + t84) - pkin(1) * t116 + qJ(2) * (t182 * (-t184 * t145 + t187 * t149) - t262), t180 * (-pkin(6) * t110 + t187 * t108) + t182 * (-pkin(2) * t110 + t85) - pkin(1) * t110 + qJ(2) * (t182 * (-t184 * t134 + t187 * t144) - t259 * t180), -t180 * t60 + qJ(2) * (t182 * (-t184 * t122 - t187 * t123) - t180 * (t166 + t215)) + t196 * t92, qJ(2) * (t182 * (t184 * t84 + t187 * t85) + t180 * t108) + t196 * t60, t180 * (t187 * (t181 * t115 + t139 * t235) - t184 * (t179 * t115 - t139 * t234)) + t213, t180 * (t187 * (-t179 * t266 - t181 * t86) - t184 * (-t179 * t86 + t181 * t266)) - t182 * t112, t180 * (t187 * (-t179 * t126 + t261) - t184 * (t181 * t126 + t263)) - t182 * t91, t180 * (t187 * (t137 * t234 + t179 * t206) - t184 * (t137 * t235 - t181 * t206)) - t213, t180 * (t187 * (t181 * t125 + t248) - t184 * (t179 * t125 - t246)) + t182 * t87, t155 + t180 * (t187 * (-t137 * t181 - t139 * t179) - t184 * (-t137 * t179 + t139 * t181)) * t168, t180 * (t187 * (-qJ(4) * t73 + t249) - t184 * (-pkin(3) * t86 + qJ(4) * t74 - t247) - pkin(6) * t39) + t182 * (-pkin(2) * t39 - t269) - pkin(1) * t39 + qJ(2) * (t182 * (-t184 * t73 + t187 * t74) + t180 * t86), t180 * (t187 * (-qJ(4) * t82 + t247) - t184 * (-pkin(3) * t266 + qJ(4) * t83 + t249) - pkin(6) * t53) + t182 * (-pkin(2) * t53 + t131 - t205) - pkin(1) * t53 + qJ(2) * (t182 * (-t184 * t82 + t187 * t83) + t180 * t266), t180 * (t187 * (-qJ(4) * t62 - t20) - t184 * (-pkin(3) * t94 + qJ(4) * t63 + t21) - pkin(6) * t36) + t182 * (-pkin(2) * t36 - t255) - pkin(1) * t36 + qJ(2) * (t182 * (-t184 * t62 + t187 * t63) + t180 * t94), t180 * (-qJ(4) * t241 - t184 * (-pkin(3) * t80 + qJ(4) * t21) - pkin(6) * t10) + t182 * (-pkin(2) * t10 - t256) - pkin(1) * t10 + qJ(2) * (t182 * (-t184 * t20 + t187 * t21) + t180 * t80), t180 * (t187 * (-t179 * t45 + t181 * t46) - t184 * (t179 * t46 + t181 * t45)) - t214, t180 * (t187 * (-t179 * t28 + t181 * t30) - t184 * (t179 * t30 + t181 * t28)) - t182 * t79, t180 * (t187 * (-t179 * t56 + t181 * t58) - t184 * (t179 * t58 + t181 * t56)) - t182 * t51, t180 * (t187 * (-t179 * t43 + t181 * t44) - t184 * (t179 * t44 + t181 * t43)) + t214, t180 * (t187 * (-t179 * t57 + t181 * t59) - t184 * (t179 * t59 + t181 * t57)) + t182 * t48, t180 * (t187 * (-t179 * t71 + t181 * t72) - t184 * (t179 * t72 + t181 * t71)) + t182 * t160, t180 * (t187 * (-qJ(4) * t22 - t179 * t18 + t181 * t24) - t184 * (-pkin(3) * t47 + qJ(4) * t23 + t179 * t24 + t181 * t18) - pkin(6) * t11) + t182 * (-pkin(2) * t11 - t193) - pkin(1) * t11 + qJ(2) * (t182 * (-t184 * t22 + t187 * t23) + t180 * t47), t180 * (t187 * (-qJ(4) * t33 - t179 * t19 + t181 * t27) - t184 * (-pkin(3) * t52 + qJ(4) * t34 + t179 * t27 + t181 * t19) - pkin(6) * t17) + t182 * (-pkin(2) * t17 - t192) - pkin(1) * t17 + qJ(2) * (t182 * (-t184 * t33 + t187 * t34) + t180 * t52), t180 * (t187 * (-qJ(4) * t15 - t179 * t4 + t181 * t5) - t184 * (-pkin(3) * t61 + qJ(4) * t16 + t179 * t5 + t181 * t4) - pkin(6) * t9) + t182 * (-pkin(2) * t9 - t203) - pkin(1) * t9 + qJ(2) * (t182 * (-t184 * t15 + t187 * t16) + t180 * t61), t180 * (t187 * (-pkin(7) * t252 - qJ(4) * t2 - t179 * t6) - t184 * (-pkin(3) * t42 - pkin(7) * t253 + qJ(4) * t3 + t181 * t6) - pkin(6) * t1) + t182 * (-pkin(2) * t1 - t204) - pkin(1) * t1 + qJ(2) * (t182 * (-t184 * t2 + t187 * t3) + t180 * t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220, t221, -t156, t148, 0, 0, 0, 0, 0, 0, t116, t110, t92, t60, 0, 0, 0, 0, 0, 0, t39, t53, t36, t10, 0, 0, 0, 0, 0, 0, t11, t17, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t150, -t122, -t158, -t123, -t167, -t84, -t85, 0, 0, -t239, t112, t91, t239, -t87, -t167, t269, t205 - 0.2e1 * t228, t255, t256, t81, t79, t51, -t81, -t48, -t160, t193, t192, t203, t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t266, t94, t80, 0, 0, 0, 0, 0, 0, t47, t52, t61, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t79, t51, -t81, -t48, -t160, -t13, -t14, 0, 0;];
tauJ_reg = t12;
