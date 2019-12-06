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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:57:16
% EndTime: 2019-12-05 17:57:29
% DurationCPUTime: 4.61s
% Computational Cost: add. (16064->374), mult. (41833->537), div. (0->0), fcn. (29414->10), ass. (0->220)
t178 = sin(pkin(8));
t174 = t178 ^ 2;
t180 = cos(pkin(8));
t175 = t180 ^ 2;
t271 = t174 + t175;
t177 = sin(pkin(9));
t179 = cos(pkin(9));
t185 = cos(qJ(3));
t228 = qJD(1) * t185;
t210 = t178 * t228;
t182 = sin(qJ(3));
t229 = qJD(1) * t178;
t211 = t182 * t229;
t139 = -t177 * t211 + t179 * t210;
t199 = -pkin(2) * t180 - pkin(6) * t178;
t154 = t199 * qJD(1);
t187 = qJD(1) ^ 2;
t183 = sin(qJ(1));
t186 = cos(qJ(1));
t197 = -g(2) * t183 + g(3) * t186;
t151 = -pkin(1) * t187 + qJDD(1) * qJ(2) - t197;
t258 = 2 * qJD(2);
t266 = qJD(1) * t258 + t151;
t201 = -t178 * g(1) + t180 * t266;
t225 = t180 * qJD(1);
t109 = t154 * t225 + t201;
t198 = t186 * g(2) + t183 * g(3);
t230 = t187 * qJ(2);
t191 = qJDD(2) - t198 - t230;
t195 = -pkin(1) + t199;
t129 = qJDD(1) * t195 + t191;
t121 = t185 * t129;
t220 = t178 * qJDD(1);
t164 = t185 * t220;
t147 = -qJD(3) * t211 + t164;
t219 = t180 * qJDD(1);
t167 = -qJDD(3) + t219;
t168 = -qJD(3) + t225;
t233 = t174 * t187;
t215 = t185 * t233;
t69 = -t167 * pkin(3) - t147 * qJ(4) + t121 + (qJ(4) * t168 * t229 - pkin(3) * t215 - t109) * t182;
t143 = -pkin(3) * t168 - qJ(4) * t210;
t146 = (qJD(3) * t228 + qJDD(1) * t182) * t178;
t166 = t182 ^ 2 * t233;
t85 = t109 * t185 + t129 * t182;
t70 = -pkin(3) * t166 - qJ(4) * t146 + t143 * t168 + t85;
t200 = -0.2e1 * qJD(4) * t139 - t177 * t70 + t179 * t69;
t137 = (-t177 * t185 - t179 * t182) * t229;
t239 = t139 * t137;
t193 = -t167 + t239;
t262 = t179 * t193;
t135 = t137 ^ 2;
t259 = t168 ^ 2;
t98 = -t259 - t135;
t73 = t177 * t98 + t262;
t270 = pkin(3) * t73 + t200;
t181 = sin(qJ(5));
t160 = -qJDD(5) + t167;
t184 = cos(qJ(5));
t105 = -t137 * t184 + t139 * t181;
t107 = t137 * t181 + t139 * t184;
t81 = t107 * t105;
t265 = -t160 - t81;
t269 = t181 * t265;
t268 = t184 * t265;
t115 = -t146 * t177 + t147 * t179;
t240 = t137 * t168;
t267 = t115 - t240;
t91 = t115 + t240;
t264 = t177 * t193;
t153 = t168 * t210;
t124 = t153 - t146;
t263 = t178 * t124;
t241 = qJDD(1) * pkin(1);
t148 = -t191 + t241;
t261 = t230 * t271 - t148 - t241;
t162 = -qJD(5) + t168;
t205 = t146 * t179 + t147 * t177;
t207 = t181 * t115 + t184 * t205;
t48 = (qJD(5) + t162) * t107 + t207;
t260 = t211 * (qJD(3) - t168) - t164;
t102 = t105 ^ 2;
t103 = t107 ^ 2;
t136 = t139 ^ 2;
t159 = t162 ^ 2;
t227 = qJD(4) * t137;
t131 = 0.2e1 * t227;
t251 = t177 * t69 + t179 * t70;
t38 = t131 + t251;
t20 = t177 * t38 + t179 * t200;
t257 = pkin(3) * t20;
t236 = t168 * t139;
t87 = t205 + t236;
t62 = -t177 * t87 - t179 * t91;
t256 = pkin(3) * t62;
t26 = t193 * pkin(4) - pkin(7) * t91 + t200;
t119 = -pkin(4) * t168 - pkin(7) * t139;
t32 = -pkin(4) * t135 - pkin(7) * t205 + t119 * t168 + t38;
t13 = t181 * t32 - t184 * t26;
t14 = t181 * t26 + t184 * t32;
t7 = -t13 * t184 + t14 * t181;
t254 = t177 * t7;
t253 = t179 * t7;
t252 = t180 * g(1);
t224 = t258 + t154;
t80 = t252 + qJDD(4) + t146 * pkin(3) - qJ(4) * t166 + (t151 + (t143 * t185 + t224) * qJD(1)) * t178;
t250 = t177 * t80;
t99 = t167 + t239;
t249 = t177 * t99;
t248 = t179 * t80;
t247 = t179 * t99;
t42 = pkin(4) * t205 - pkin(7) * t135 + t139 * t119 + t80;
t246 = t181 * t42;
t76 = t160 - t81;
t245 = t181 * t76;
t244 = t184 * t42;
t243 = t184 * t76;
t242 = t185 * t20;
t238 = t162 * t181;
t237 = t162 * t184;
t235 = t168 * t177;
t234 = t168 * t179;
t158 = t182 * t215;
t144 = -t158 + t167;
t232 = t182 * t144;
t145 = -t158 - t167;
t231 = t185 * t145;
t222 = qJD(5) - t162;
t214 = t185 ^ 2 * t233;
t213 = t180 * t81;
t212 = t180 * t239;
t8 = t13 * t181 + t14 * t184;
t21 = -t177 * t200 + t179 * t38;
t84 = t182 * t109 - t121;
t206 = t178 * (t178 * t266 + t252) + t180 * t201;
t117 = -t136 - t259;
t82 = t117 * t179 + t249;
t204 = pkin(3) * t82 - t251;
t2 = t177 * t8 + t253;
t203 = pkin(3) * t2 + pkin(4) * t7;
t196 = t184 * t115 - t181 * t205;
t65 = -qJD(5) * t105 + t196;
t97 = t105 * t162;
t51 = t65 - t97;
t29 = -t181 * t48 - t184 * t51;
t31 = t181 * t51 - t184 * t48;
t15 = t177 * t31 + t179 * t29;
t202 = pkin(3) * t15 + pkin(4) * t29;
t60 = t182 * t85 - t185 * t84;
t75 = -t159 - t102;
t40 = t181 * t75 + t268;
t41 = t184 * t75 - t269;
t22 = t177 * t41 + t179 * t40;
t192 = pkin(3) * t22 + pkin(4) * t40 - t13;
t93 = -t103 - t159;
t54 = t184 * t93 + t245;
t55 = -t181 * t93 + t243;
t33 = t177 * t55 + t179 * t54;
t190 = pkin(3) * t33 + pkin(4) * t54 - t14;
t173 = t175 * qJDD(1);
t172 = t174 * qJDD(1);
t156 = t271 * t187;
t155 = t180 * t167;
t150 = -t166 + t214;
t149 = -t166 - t259;
t134 = -t259 - t214;
t126 = -t136 + t259;
t125 = t135 - t259;
t123 = t153 + t146;
t122 = -t164 + (qJD(3) + t168) * t211;
t116 = t149 * t182 + t231;
t112 = t136 - t135;
t110 = t134 * t185 + t232;
t108 = t252 + (qJD(1) * t224 + t151) * t178;
t96 = -t103 + t159;
t95 = t102 - t159;
t94 = -t135 - t136;
t92 = t122 * t185 - t123 * t182;
t86 = t205 - t236;
t83 = -t117 * t177 + t247;
t79 = t103 - t102;
t74 = t179 * t98 - t264;
t72 = (t105 * t184 - t107 * t181) * t162;
t71 = (t105 * t181 + t107 * t184) * t162;
t64 = -qJD(5) * t107 - t207;
t63 = t177 * t91 - t179 * t87;
t61 = -t102 - t103;
t59 = t184 * t95 + t245;
t58 = -t181 * t96 + t268;
t57 = t181 * t95 - t243;
t56 = t184 * t96 + t269;
t53 = t182 * t83 + t185 * t82;
t52 = -t105 * t222 + t196;
t50 = t65 + t97;
t47 = t107 * t222 + t207;
t46 = t107 * t238 + t184 * t65;
t45 = -t107 * t237 + t181 * t65;
t44 = -t105 * t237 - t181 * t64;
t43 = -t105 * t238 + t184 * t64;
t39 = t182 * t74 + t185 * t73;
t36 = t182 * t63 + t185 * t62;
t34 = -t177 * t54 + t179 * t55;
t30 = -t181 * t50 - t184 * t47;
t28 = -t181 * t47 + t184 * t50;
t27 = -pkin(7) * t54 + t244;
t24 = -pkin(7) * t40 + t246;
t23 = -t177 * t40 + t179 * t41;
t19 = -pkin(4) * t52 + pkin(7) * t55 + t246;
t18 = -pkin(4) * t47 + pkin(7) * t41 - t244;
t17 = t182 * t34 + t185 * t33;
t16 = -t177 * t29 + t179 * t31;
t11 = t182 * t23 + t185 * t22;
t10 = t182 * t21 + t242;
t9 = t15 * t185 + t16 * t182;
t6 = -pkin(4) * t42 + pkin(7) * t8;
t5 = -pkin(7) * t29 - t7;
t4 = -pkin(4) * t61 + pkin(7) * t31 + t8;
t3 = t179 * t8 - t254;
t1 = t182 * t3 + t185 * t2;
t12 = [0, 0, 0, 0, 0, qJDD(1), t198, t197, 0, 0, t172, 0.2e1 * t178 * t219, 0, t173, 0, 0, -t261 * t180, t261 * t178, pkin(1) * t156 + qJ(2) * (t173 + t172) + t206, pkin(1) * t148 + qJ(2) * t206, (t178 * (t168 * t211 + t147) - t180 * t182 * t233) * t185, t178 * (t185 * t124 + t182 * t260) - t180 * t150, t178 * (t231 - t182 * (t259 - t214)) + t180 * t122, (t180 * t215 - t263) * t182, t178 * (t185 * (t166 - t259) + t232) + t180 * t123, t155, t178 * (-pkin(6) * t116 + t108 * t182) + t180 * (-pkin(2) * t116 + t84) - pkin(1) * t116 + qJ(2) * (t180 * (-t145 * t182 + t149 * t185) - t263), t178 * (-pkin(6) * t110 + t185 * t108) + t180 * (-pkin(2) * t110 + t85) - pkin(1) * t110 + qJ(2) * (t180 * (-t134 * t182 + t144 * t185) - t260 * t178), -t178 * t60 + qJ(2) * (t180 * (-t122 * t182 - t123 * t185) - t178 * (t166 + t214)) + t195 * t92, qJ(2) * (t180 * (t182 * t84 + t185 * t85) + t178 * t108) + t195 * t60, t178 * (t185 * (t115 * t179 + t139 * t235) - t182 * (t115 * t177 - t139 * t234)) + t212, t178 * (t185 * (-t177 * t267 - t179 * t86) - t182 * (-t177 * t86 + t179 * t267)) - t180 * t112, t178 * (t185 * (-t126 * t177 + t262) - t182 * (t126 * t179 + t264)) - t180 * t91, t178 * (t185 * (t137 * t234 + t177 * t205) - t182 * (t137 * t235 - t179 * t205)) - t212, t178 * (t185 * (t125 * t179 + t249) - t182 * (t125 * t177 - t247)) + t180 * t87, t155 + t178 * (t185 * (-t137 * t179 - t139 * t177) - t182 * (-t137 * t177 + t139 * t179)) * t168, t178 * (t185 * (-qJ(4) * t73 + t250) - t182 * (-pkin(3) * t86 + qJ(4) * t74 - t248) - pkin(6) * t39) + t180 * (-pkin(2) * t39 - t270) - pkin(1) * t39 + qJ(2) * (t180 * (-t182 * t73 + t185 * t74) + t178 * t86), t178 * (t185 * (-qJ(4) * t82 + t248) - t182 * (-pkin(3) * t267 + qJ(4) * t83 + t250) - pkin(6) * t53) + t180 * (-pkin(2) * t53 + t131 - t204) - pkin(1) * t53 + qJ(2) * (t180 * (-t182 * t82 + t185 * t83) + t178 * t267), t178 * (t185 * (-qJ(4) * t62 - t20) - t182 * (-pkin(3) * t94 + qJ(4) * t63 + t21) - pkin(6) * t36) + t180 * (-pkin(2) * t36 - t256) - pkin(1) * t36 + qJ(2) * (t180 * (-t182 * t62 + t185 * t63) + t178 * t94), t178 * (-qJ(4) * t242 - t182 * (-pkin(3) * t80 + qJ(4) * t21) - pkin(6) * t10) + t180 * (-pkin(2) * t10 - t257) - pkin(1) * t10 + qJ(2) * (t180 * (-t182 * t20 + t185 * t21) + t178 * t80), t178 * (t185 * (-t177 * t45 + t179 * t46) - t182 * (t177 * t46 + t179 * t45)) - t213, t178 * (t185 * (-t177 * t28 + t179 * t30) - t182 * (t177 * t30 + t179 * t28)) - t180 * t79, t178 * (t185 * (-t177 * t56 + t179 * t58) - t182 * (t177 * t58 + t179 * t56)) - t180 * t51, t178 * (t185 * (-t177 * t43 + t179 * t44) - t182 * (t177 * t44 + t179 * t43)) + t213, t178 * (t185 * (-t177 * t57 + t179 * t59) - t182 * (t177 * t59 + t179 * t57)) + t180 * t48, t178 * (t185 * (-t177 * t71 + t179 * t72) - t182 * (t177 * t72 + t179 * t71)) + t180 * t160, t178 * (t185 * (-qJ(4) * t22 - t177 * t18 + t179 * t24) - t182 * (-pkin(3) * t47 + qJ(4) * t23 + t177 * t24 + t179 * t18) - pkin(6) * t11) + t180 * (-pkin(2) * t11 - t192) - pkin(1) * t11 + qJ(2) * (t180 * (-t182 * t22 + t185 * t23) + t178 * t47), t178 * (t185 * (-qJ(4) * t33 - t177 * t19 + t179 * t27) - t182 * (-pkin(3) * t52 + qJ(4) * t34 + t177 * t27 + t179 * t19) - pkin(6) * t17) + t180 * (-pkin(2) * t17 - t190) - pkin(1) * t17 + qJ(2) * (t180 * (-t182 * t33 + t185 * t34) + t178 * t52), t178 * (t185 * (-qJ(4) * t15 - t177 * t4 + t179 * t5) - t182 * (-pkin(3) * t61 + qJ(4) * t16 + t177 * t5 + t179 * t4) - pkin(6) * t9) + t180 * (-pkin(2) * t9 - t202) - pkin(1) * t9 + qJ(2) * (t180 * (-t15 * t182 + t16 * t185) + t178 * t61), t178 * (t185 * (-pkin(7) * t253 - qJ(4) * t2 - t177 * t6) - t182 * (-pkin(3) * t42 - pkin(7) * t254 + qJ(4) * t3 + t179 * t6) - pkin(6) * t1) + t180 * (-pkin(2) * t1 - t203) - pkin(1) * t1 + qJ(2) * (t180 * (-t182 * t2 + t185 * t3) + t178 * t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t219, t220, -t156, -t148, 0, 0, 0, 0, 0, 0, t116, t110, t92, t60, 0, 0, 0, 0, 0, 0, t39, t53, t36, t10, 0, 0, 0, 0, 0, 0, t11, t17, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t150, -t122, -t158, -t123, -t167, -t84, -t85, 0, 0, -t239, t112, t91, t239, -t87, -t167, t270, t204 - 0.2e1 * t227, t256, t257, t81, t79, t51, -t81, -t48, -t160, t192, t190, t202, t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t267, t94, t80, 0, 0, 0, 0, 0, 0, t47, t52, t61, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t79, t51, -t81, -t48, -t160, -t13, -t14, 0, 0;];
tauJ_reg = t12;
