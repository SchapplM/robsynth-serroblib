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
% Datum: 2022-01-23 09:26
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:26:02
% EndTime: 2022-01-23 09:26:09
% DurationCPUTime: 4.52s
% Computational Cost: add. (16064->374), mult. (41833->537), div. (0->0), fcn. (29414->10), ass. (0->220)
t179 = sin(pkin(8));
t175 = t179 ^ 2;
t181 = cos(pkin(8));
t176 = t181 ^ 2;
t272 = t175 + t176;
t178 = sin(pkin(9));
t180 = cos(pkin(9));
t186 = cos(qJ(3));
t229 = qJD(1) * t186;
t211 = t179 * t229;
t183 = sin(qJ(3));
t230 = qJD(1) * t179;
t212 = t183 * t230;
t139 = -t178 * t212 + t180 * t211;
t199 = -pkin(2) * t181 - pkin(6) * t179;
t154 = t199 * qJD(1);
t188 = qJD(1) ^ 2;
t184 = sin(qJ(1));
t187 = cos(qJ(1));
t198 = g(1) * t187 + g(2) * t184;
t151 = -pkin(1) * t188 + qJDD(1) * qJ(2) - t198;
t259 = 2 * qJD(2);
t267 = qJD(1) * t259 + t151;
t201 = -t179 * g(3) + t181 * t267;
t226 = t181 * qJD(1);
t109 = t154 * t226 + t201;
t210 = t184 * g(1) - g(2) * t187;
t231 = t188 * qJ(2);
t193 = qJDD(2) - t210 - t231;
t196 = -pkin(1) + t199;
t129 = qJDD(1) * t196 + t193;
t121 = t186 * t129;
t221 = t179 * qJDD(1);
t164 = t186 * t221;
t147 = -qJD(3) * t212 + t164;
t220 = t181 * qJDD(1);
t167 = -qJDD(3) + t220;
t168 = -qJD(3) + t226;
t234 = t175 * t188;
t216 = t186 * t234;
t69 = -t167 * pkin(3) - t147 * qJ(4) + t121 + (qJ(4) * t168 * t230 - pkin(3) * t216 - t109) * t183;
t143 = -pkin(3) * t168 - qJ(4) * t211;
t146 = (qJD(3) * t229 + qJDD(1) * t183) * t179;
t166 = t183 ^ 2 * t234;
t85 = t109 * t186 + t129 * t183;
t70 = -pkin(3) * t166 - qJ(4) * t146 + t143 * t168 + t85;
t200 = -0.2e1 * qJD(4) * t139 - t178 * t70 + t180 * t69;
t137 = (-t178 * t186 - t180 * t183) * t230;
t240 = t139 * t137;
t194 = -t167 + t240;
t263 = t180 * t194;
t135 = t137 ^ 2;
t260 = t168 ^ 2;
t98 = -t260 - t135;
t73 = t178 * t98 + t263;
t271 = pkin(3) * t73 + t200;
t182 = sin(qJ(5));
t160 = -qJDD(5) + t167;
t185 = cos(qJ(5));
t105 = -t137 * t185 + t139 * t182;
t107 = t137 * t182 + t139 * t185;
t81 = t107 * t105;
t266 = -t160 - t81;
t270 = t182 * t266;
t269 = t185 * t266;
t115 = -t146 * t178 + t147 * t180;
t241 = t137 * t168;
t268 = t115 - t241;
t91 = t115 + t241;
t265 = t178 * t194;
t153 = t168 * t211;
t124 = t153 - t146;
t264 = t179 * t124;
t242 = qJDD(1) * pkin(1);
t148 = -t193 + t242;
t262 = t231 * t272 - t148 - t242;
t162 = -qJD(5) + t168;
t205 = t146 * t180 + t147 * t178;
t207 = t182 * t115 + t185 * t205;
t48 = (qJD(5) + t162) * t107 + t207;
t261 = (qJD(3) - t168) * t212 - t164;
t102 = t105 ^ 2;
t103 = t107 ^ 2;
t136 = t139 ^ 2;
t159 = t162 ^ 2;
t228 = qJD(4) * t137;
t131 = 0.2e1 * t228;
t252 = t178 * t69 + t180 * t70;
t38 = t131 + t252;
t20 = t178 * t38 + t180 * t200;
t258 = pkin(3) * t20;
t237 = t168 * t139;
t87 = t205 + t237;
t62 = -t178 * t87 - t180 * t91;
t257 = pkin(3) * t62;
t26 = t194 * pkin(4) - pkin(7) * t91 + t200;
t119 = -pkin(4) * t168 - pkin(7) * t139;
t32 = -pkin(4) * t135 - pkin(7) * t205 + t119 * t168 + t38;
t13 = t182 * t32 - t185 * t26;
t14 = t182 * t26 + t185 * t32;
t7 = -t13 * t185 + t14 * t182;
t255 = t178 * t7;
t254 = t180 * t7;
t253 = t181 * g(3);
t225 = t259 + t154;
t80 = t253 + qJDD(4) + t146 * pkin(3) - qJ(4) * t166 + (t151 + (t143 * t186 + t225) * qJD(1)) * t179;
t251 = t178 * t80;
t99 = t167 + t240;
t250 = t178 * t99;
t249 = t180 * t80;
t248 = t180 * t99;
t42 = pkin(4) * t205 - pkin(7) * t135 + t139 * t119 + t80;
t247 = t182 * t42;
t76 = t160 - t81;
t246 = t182 * t76;
t245 = t185 * t42;
t244 = t185 * t76;
t243 = t186 * t20;
t239 = t162 * t182;
t238 = t162 * t185;
t236 = t168 * t178;
t235 = t168 * t180;
t158 = t183 * t216;
t144 = -t158 + t167;
t233 = t183 * t144;
t145 = -t158 - t167;
t232 = t186 * t145;
t223 = qJD(5) - t162;
t215 = t186 ^ 2 * t234;
t214 = t181 * t81;
t213 = t181 * t240;
t8 = t13 * t182 + t14 * t185;
t21 = -t178 * t200 + t180 * t38;
t84 = t183 * t109 - t121;
t206 = t179 * (t179 * t267 + t253) + t181 * t201;
t117 = -t136 - t260;
t82 = t117 * t180 + t250;
t204 = pkin(3) * t82 - t252;
t2 = t178 * t8 + t254;
t203 = pkin(3) * t2 + pkin(4) * t7;
t197 = t185 * t115 - t182 * t205;
t65 = -qJD(5) * t105 + t197;
t97 = t105 * t162;
t51 = t65 - t97;
t29 = -t182 * t48 - t185 * t51;
t31 = t182 * t51 - t185 * t48;
t15 = t178 * t31 + t180 * t29;
t202 = pkin(3) * t15 + pkin(4) * t29;
t60 = t183 * t85 - t186 * t84;
t75 = -t159 - t102;
t40 = t182 * t75 + t269;
t41 = t185 * t75 - t270;
t22 = t178 * t41 + t180 * t40;
t192 = pkin(3) * t22 + pkin(4) * t40 - t13;
t93 = -t103 - t159;
t54 = t185 * t93 + t246;
t55 = -t182 * t93 + t244;
t33 = t178 * t55 + t180 * t54;
t191 = pkin(3) * t33 + pkin(4) * t54 - t14;
t173 = t176 * qJDD(1);
t172 = t175 * qJDD(1);
t156 = t272 * t188;
t155 = t181 * t167;
t150 = -t166 + t215;
t149 = -t166 - t260;
t134 = -t260 - t215;
t126 = -t136 + t260;
t125 = t135 - t260;
t123 = t153 + t146;
t122 = -t164 + (qJD(3) + t168) * t212;
t116 = t149 * t183 + t232;
t112 = t136 - t135;
t110 = t134 * t186 + t233;
t108 = t253 + (qJD(1) * t225 + t151) * t179;
t96 = -t103 + t159;
t95 = t102 - t159;
t94 = -t135 - t136;
t92 = t122 * t186 - t123 * t183;
t86 = t205 - t237;
t83 = -t117 * t178 + t248;
t79 = t103 - t102;
t74 = t180 * t98 - t265;
t72 = (t105 * t185 - t107 * t182) * t162;
t71 = (t105 * t182 + t107 * t185) * t162;
t64 = -qJD(5) * t107 - t207;
t63 = t178 * t91 - t180 * t87;
t61 = -t102 - t103;
t59 = t185 * t95 + t246;
t58 = -t182 * t96 + t269;
t57 = t182 * t95 - t244;
t56 = t185 * t96 + t270;
t53 = t183 * t83 + t186 * t82;
t52 = -t105 * t223 + t197;
t50 = t65 + t97;
t47 = t107 * t223 + t207;
t46 = t107 * t239 + t185 * t65;
t45 = -t107 * t238 + t182 * t65;
t44 = -t105 * t238 - t182 * t64;
t43 = -t105 * t239 + t185 * t64;
t39 = t183 * t74 + t186 * t73;
t36 = t183 * t63 + t186 * t62;
t34 = -t178 * t54 + t180 * t55;
t30 = -t182 * t50 - t185 * t47;
t28 = -t182 * t47 + t185 * t50;
t27 = -pkin(7) * t54 + t245;
t24 = -pkin(7) * t40 + t247;
t23 = -t178 * t40 + t180 * t41;
t19 = -pkin(4) * t52 + pkin(7) * t55 + t247;
t18 = -pkin(4) * t47 + pkin(7) * t41 - t245;
t17 = t183 * t34 + t186 * t33;
t16 = -t178 * t29 + t180 * t31;
t11 = t183 * t23 + t186 * t22;
t10 = t183 * t21 + t243;
t9 = t15 * t186 + t16 * t183;
t6 = -pkin(4) * t42 + pkin(7) * t8;
t5 = -pkin(7) * t29 - t7;
t4 = -pkin(4) * t61 + pkin(7) * t31 + t8;
t3 = t180 * t8 - t255;
t1 = t183 * t3 + t186 * t2;
t12 = [0, 0, 0, 0, 0, qJDD(1), t210, t198, 0, 0, t172, 0.2e1 * t179 * t220, 0, t173, 0, 0, -t262 * t181, t262 * t179, pkin(1) * t156 + qJ(2) * (t173 + t172) + t206, pkin(1) * t148 + qJ(2) * t206, (t179 * (t168 * t212 + t147) - t181 * t183 * t234) * t186, t179 * (t186 * t124 + t183 * t261) - t181 * t150, t179 * (t232 - t183 * (t260 - t215)) + t181 * t122, (t181 * t216 - t264) * t183, t179 * (t186 * (t166 - t260) + t233) + t181 * t123, t155, t179 * (-pkin(6) * t116 + t108 * t183) + t181 * (-pkin(2) * t116 + t84) - pkin(1) * t116 + qJ(2) * (t181 * (-t145 * t183 + t149 * t186) - t264), t179 * (-pkin(6) * t110 + t186 * t108) + t181 * (-pkin(2) * t110 + t85) - pkin(1) * t110 + qJ(2) * (t181 * (-t134 * t183 + t144 * t186) - t261 * t179), -t179 * t60 + qJ(2) * (t181 * (-t122 * t183 - t123 * t186) - t179 * (t166 + t215)) + t196 * t92, qJ(2) * (t181 * (t183 * t84 + t186 * t85) + t179 * t108) + t196 * t60, t179 * (t186 * (t115 * t180 + t139 * t236) - t183 * (t115 * t178 - t139 * t235)) + t213, t179 * (t186 * (-t178 * t268 - t180 * t86) - t183 * (-t178 * t86 + t180 * t268)) - t181 * t112, t179 * (t186 * (-t126 * t178 + t263) - t183 * (t126 * t180 + t265)) - t181 * t91, t179 * (t186 * (t137 * t235 + t178 * t205) - t183 * (t137 * t236 - t180 * t205)) - t213, t179 * (t186 * (t125 * t180 + t250) - t183 * (t125 * t178 - t248)) + t181 * t87, t155 + t179 * (t186 * (-t137 * t180 - t139 * t178) - t183 * (-t137 * t178 + t139 * t180)) * t168, t179 * (t186 * (-qJ(4) * t73 + t251) - t183 * (-pkin(3) * t86 + qJ(4) * t74 - t249) - pkin(6) * t39) + t181 * (-pkin(2) * t39 - t271) - pkin(1) * t39 + qJ(2) * (t181 * (-t183 * t73 + t186 * t74) + t179 * t86), t179 * (t186 * (-qJ(4) * t82 + t249) - t183 * (-pkin(3) * t268 + qJ(4) * t83 + t251) - pkin(6) * t53) + t181 * (-pkin(2) * t53 + t131 - t204) - pkin(1) * t53 + qJ(2) * (t181 * (-t183 * t82 + t186 * t83) + t179 * t268), t179 * (t186 * (-qJ(4) * t62 - t20) - t183 * (-pkin(3) * t94 + qJ(4) * t63 + t21) - pkin(6) * t36) + t181 * (-pkin(2) * t36 - t257) - pkin(1) * t36 + qJ(2) * (t181 * (-t183 * t62 + t186 * t63) + t179 * t94), t179 * (-qJ(4) * t243 - t183 * (-pkin(3) * t80 + qJ(4) * t21) - pkin(6) * t10) + t181 * (-pkin(2) * t10 - t258) - pkin(1) * t10 + qJ(2) * (t181 * (-t183 * t20 + t186 * t21) + t179 * t80), t179 * (t186 * (-t178 * t45 + t180 * t46) - t183 * (t178 * t46 + t180 * t45)) - t214, t179 * (t186 * (-t178 * t28 + t180 * t30) - t183 * (t178 * t30 + t180 * t28)) - t181 * t79, t179 * (t186 * (-t178 * t56 + t180 * t58) - t183 * (t178 * t58 + t180 * t56)) - t181 * t51, t179 * (t186 * (-t178 * t43 + t180 * t44) - t183 * (t178 * t44 + t180 * t43)) + t214, t179 * (t186 * (-t178 * t57 + t180 * t59) - t183 * (t178 * t59 + t180 * t57)) + t181 * t48, t179 * (t186 * (-t178 * t71 + t180 * t72) - t183 * (t178 * t72 + t180 * t71)) + t181 * t160, t179 * (t186 * (-qJ(4) * t22 - t178 * t18 + t180 * t24) - t183 * (-pkin(3) * t47 + qJ(4) * t23 + t178 * t24 + t18 * t180) - pkin(6) * t11) + t181 * (-pkin(2) * t11 - t192) - pkin(1) * t11 + qJ(2) * (t181 * (-t183 * t22 + t186 * t23) + t179 * t47), t179 * (t186 * (-qJ(4) * t33 - t178 * t19 + t180 * t27) - t183 * (-pkin(3) * t52 + qJ(4) * t34 + t178 * t27 + t180 * t19) - pkin(6) * t17) + t181 * (-pkin(2) * t17 - t191) - pkin(1) * t17 + qJ(2) * (t181 * (-t183 * t33 + t186 * t34) + t179 * t52), t179 * (t186 * (-qJ(4) * t15 - t178 * t4 + t180 * t5) - t183 * (-pkin(3) * t61 + qJ(4) * t16 + t178 * t5 + t180 * t4) - pkin(6) * t9) + t181 * (-pkin(2) * t9 - t202) - pkin(1) * t9 + qJ(2) * (t181 * (-t15 * t183 + t16 * t186) + t179 * t61), t179 * (t186 * (-pkin(7) * t254 - qJ(4) * t2 - t178 * t6) - t183 * (-pkin(3) * t42 - pkin(7) * t255 + qJ(4) * t3 + t180 * t6) - pkin(6) * t1) + t181 * (-pkin(2) * t1 - t203) - pkin(1) * t1 + qJ(2) * (t181 * (-t183 * t2 + t186 * t3) + t179 * t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220, t221, -t156, -t148, 0, 0, 0, 0, 0, 0, t116, t110, t92, t60, 0, 0, 0, 0, 0, 0, t39, t53, t36, t10, 0, 0, 0, 0, 0, 0, t11, t17, t9, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t150, -t122, -t158, -t123, -t167, -t84, -t85, 0, 0, -t240, t112, t91, t240, -t87, -t167, t271, t204 - 0.2e1 * t228, t257, t258, t81, t79, t51, -t81, -t48, -t160, t192, t191, t202, t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t268, t94, t80, 0, 0, 0, 0, 0, 0, t47, t52, t61, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t79, t51, -t81, -t48, -t160, -t13, -t14, 0, 0;];
tauJ_reg = t12;
