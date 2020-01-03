% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:19
% EndTime: 2019-12-31 20:18:30
% DurationCPUTime: 6.00s
% Computational Cost: add. (9145->483), mult. (22035->625), div. (0->0), fcn. (16561->14), ass. (0->250)
t197 = cos(qJ(5));
t264 = qJD(5) * t197;
t190 = sin(pkin(9));
t191 = cos(pkin(9));
t198 = cos(qJ(2));
t277 = t191 * t198;
t248 = qJD(1) * t277;
t195 = sin(qJ(2));
t267 = qJD(1) * t195;
t134 = -t190 * t267 + t248;
t194 = sin(qJ(4));
t311 = cos(qJ(4));
t142 = t190 * t198 + t191 * t195;
t314 = t142 * qJD(1);
t86 = t311 * t134 - t194 * t314;
t330 = t197 * t86;
t342 = t264 - t330;
t193 = sin(qJ(5));
t247 = qJD(4) * t311;
t266 = qJD(4) * t194;
t261 = t198 * qJDD(1);
t262 = t195 * qJDD(1);
t226 = t190 * t262 - t191 * t261;
t324 = t314 * qJD(2);
t95 = t226 + t324;
t263 = qJD(1) * qJD(2);
t246 = t195 * t263;
t208 = qJDD(1) * t142 - t190 * t246;
t245 = t198 * t263;
t96 = t191 * t245 + t208;
t216 = t134 * t247 - t194 * t95 - t266 * t314 + t311 * t96;
t218 = -t194 * t134 - t311 * t314;
t260 = qJD(2) + qJD(4);
t234 = t197 * t260;
t259 = qJDD(2) + qJDD(4);
t265 = qJD(5) * t193;
t27 = -qJD(5) * t234 - t193 * t259 - t197 * t216 - t218 * t265;
t72 = t193 * t260 - t197 * t218;
t284 = qJD(5) * t72;
t28 = t193 * t216 - t197 * t259 + t284;
t70 = -t193 * t218 - t234;
t294 = -t193 * t28 - t70 * t264;
t328 = qJD(5) - t86;
t338 = t193 * t328;
t318 = t72 * t338;
t341 = -t197 * t27 + t330 * t70 + t294 - t318;
t240 = t194 * t96 + t311 * t95;
t325 = qJD(4) * t218;
t41 = t240 - t325;
t301 = t198 * pkin(2);
t179 = pkin(1) + t301;
t124 = pkin(2) * t246 - qJDD(1) * t179 + qJDD(3);
t312 = pkin(3) * t95;
t67 = t124 + t312;
t12 = pkin(4) * t41 - pkin(8) * t216 + t67;
t11 = t197 * t12;
t307 = pkin(7) * t314;
t192 = -qJ(3) - pkin(6);
t161 = t192 * t198;
t150 = qJD(1) * t161;
t137 = t190 * t150;
t160 = t192 * t195;
t149 = qJD(1) * t160;
t291 = qJD(2) * pkin(2);
t141 = t149 + t291;
t93 = t191 * t141 + t137;
t64 = qJD(2) * pkin(3) - t307 + t93;
t308 = pkin(7) * t134;
t139 = t191 * t150;
t94 = t190 * t141 - t139;
t69 = t94 + t308;
t45 = t194 * t64 + t311 * t69;
t35 = pkin(8) * t260 + t45;
t152 = -qJD(1) * t179 + qJD(3);
t102 = -pkin(3) * t134 + t152;
t47 = -pkin(4) * t86 + pkin(8) * t218 + t102;
t14 = t193 * t47 + t197 * t35;
t239 = qJD(2) * t192;
t131 = -qJD(3) * t195 + t198 * t239;
t92 = qJDD(2) * pkin(2) + qJD(1) * t131 + qJDD(1) * t160;
t130 = qJD(3) * t198 + t195 * t239;
t99 = qJD(1) * t130 - qJDD(1) * t161;
t56 = -t190 * t99 + t191 * t92;
t43 = qJDD(2) * pkin(3) - pkin(7) * t96 + t56;
t57 = t190 * t92 + t191 * t99;
t46 = -pkin(7) * t95 + t57;
t9 = t194 * t43 + t64 * t247 - t69 * t266 + t311 * t46;
t7 = pkin(8) * t259 + t9;
t3 = -qJD(5) * t14 - t193 * t7 + t11;
t302 = t14 * t328;
t340 = t3 + t302;
t13 = -t193 * t35 + t197 * t47;
t303 = t13 * t328;
t178 = pkin(2) * t191 + pkin(3);
t310 = pkin(2) * t190;
t129 = t194 * t178 + t311 * t310;
t100 = -t149 * t190 + t139;
t211 = t100 - t308;
t101 = t191 * t149 + t137;
t73 = t101 - t307;
t286 = qJD(4) * t129 - t194 * t73 + t311 * t211;
t339 = -t45 - t286;
t319 = t260 * t86;
t337 = t216 - t319;
t296 = t86 ^ 2;
t297 = t218 ^ 2;
t336 = -t296 + t297;
t24 = t27 * t193;
t334 = t342 * t72 - t24;
t298 = t72 * t218;
t39 = qJDD(5) + t41;
t32 = t193 * t39;
t333 = t342 * t328 + t298 + t32;
t44 = -t194 * t69 + t311 * t64;
t34 = -pkin(4) * t260 - t44;
t332 = t34 * t86;
t187 = qJ(2) + pkin(9);
t184 = qJ(4) + t187;
t177 = cos(t184);
t305 = g(3) * t177;
t241 = t194 * t46 - t311 * t43;
t10 = -qJD(4) * t45 - t241;
t8 = -pkin(4) * t259 - t10;
t331 = t8 + t305;
t300 = t70 * t218;
t329 = t328 * t218;
t295 = t86 * t218;
t176 = sin(t184);
t199 = cos(qJ(1));
t280 = t176 * t199;
t196 = sin(qJ(1));
t281 = t176 * t196;
t327 = g(1) * t280 + g(2) * t281;
t272 = t218 * qJD(2);
t326 = -t272 - t240;
t169 = g(3) * t176;
t278 = t177 * t199;
t279 = t177 * t196;
t250 = -g(1) * t278 - g(2) * t279 - t169;
t323 = -t102 * t86 - t250 - t9;
t30 = t34 * t265;
t322 = t13 * t218 + t327 * t197 + t30;
t31 = t34 * t264;
t321 = -t14 * t218 + t331 * t193 + t31;
t320 = t102 * t218 - t241 - t305 + t327;
t55 = -pkin(4) * t218 - pkin(8) * t86;
t128 = t311 * t178 - t194 * t310;
t111 = t128 * qJD(4);
t49 = t194 * t211 + t311 * t73;
t287 = t111 - t49;
t315 = g(1) * t196 - g(2) * t199;
t317 = t315 * t176;
t316 = t177 * pkin(4) + t176 * pkin(8);
t313 = t314 ^ 2;
t200 = qJD(2) ^ 2;
t309 = pkin(2) * t195;
t304 = g(3) * t198;
t2 = qJD(5) * t13 + t12 * t193 + t197 * t7;
t1 = t2 * t197;
t299 = t72 * t70;
t290 = t193 * t70;
t26 = t28 * t197;
t285 = pkin(6) * qJDD(1);
t283 = qJD(5) * t328;
t282 = t314 * t134;
t276 = t193 * t196;
t275 = t193 * t199;
t274 = t196 * t197;
t273 = t197 * t199;
t77 = t191 * t130 + t190 * t131;
t104 = t190 * t160 - t191 * t161;
t188 = t195 ^ 2;
t189 = t198 ^ 2;
t269 = t188 - t189;
t268 = t188 + t189;
t256 = pkin(8) * t283;
t181 = t195 * t291;
t126 = pkin(8) + t129;
t254 = t126 * t283;
t223 = t190 * t195 - t277;
t98 = t311 * t142 - t194 * t223;
t253 = t98 * t265;
t252 = t98 * t264;
t201 = qJD(1) ^ 2;
t251 = t195 * t201 * t198;
t183 = cos(t187);
t244 = pkin(3) * t183 + t301;
t105 = pkin(2) * t267 + pkin(3) * t314;
t182 = sin(t187);
t151 = -pkin(3) * t182 - t309;
t243 = -pkin(4) * t176 + t151;
t103 = t191 * t160 + t161 * t190;
t235 = t1 + t250;
t233 = t195 * t245;
t232 = -pkin(8) * t39 - t332;
t209 = t311 * t223;
t215 = t142 * qJD(2);
t58 = t142 * t266 + t194 * t215 + t209 * t260;
t231 = -t34 * t58 + t8 * t98;
t229 = g(1) * t199 + g(2) * t196;
t227 = -t328 * t58 + t39 * t98;
t225 = t13 * t197 + t14 * t193;
t224 = t13 * t193 - t14 * t197;
t108 = pkin(3) * t223 - t179;
t97 = t142 * t194 + t209;
t51 = t97 * pkin(4) - t98 * pkin(8) + t108;
t78 = -pkin(7) * t142 + t103;
t79 = -pkin(7) * t223 + t104;
t53 = t194 * t78 + t311 * t79;
t21 = -t193 * t53 + t197 * t51;
t22 = t193 * t51 + t197 * t53;
t76 = -t190 * t130 + t191 * t131;
t33 = t197 * t39;
t222 = -t265 * t328 + t338 * t86 + t33;
t221 = -qJD(5) * t47 + t169 - t7;
t220 = t229 * t176;
t219 = -0.2e1 * pkin(1) * t263 - pkin(6) * qJDD(2);
t213 = t223 * qJD(2);
t212 = -t111 * t328 - t126 * t39 - t332;
t106 = pkin(3) * t215 + t181;
t207 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t200 + t315;
t206 = pkin(1) * t201 + t229 - t285;
t205 = t124 - t315;
t204 = -qJD(5) * t225 - t3 * t193 + t1;
t203 = pkin(7) * t213 + t76;
t186 = -pkin(7) + t192;
t154 = pkin(8) * t278;
t153 = pkin(8) * t279;
t148 = pkin(1) + t244;
t140 = t199 * t148;
t132 = t134 ^ 2;
t125 = -pkin(4) - t128;
t123 = t177 * t273 + t276;
t122 = -t177 * t275 + t274;
t121 = -t177 * t274 + t275;
t120 = t177 * t276 + t273;
t62 = -pkin(7) * t215 + t77;
t59 = t98 * qJD(4) - t194 * t213 + t311 * t215;
t52 = t194 * t79 - t311 * t78;
t50 = t105 + t55;
t29 = t59 * pkin(4) + t58 * pkin(8) + t106;
t20 = t193 * t55 + t197 * t44;
t19 = -t193 * t44 + t197 * t55;
t18 = t53 * qJD(4) + t194 * t62 - t311 * t203;
t17 = t194 * t203 + t78 * t247 - t79 * t266 + t311 * t62;
t16 = t193 * t50 + t197 * t49;
t15 = -t193 * t49 + t197 * t50;
t5 = -qJD(5) * t22 - t17 * t193 + t197 * t29;
t4 = qJD(5) * t21 + t17 * t197 + t193 * t29;
t6 = [0, 0, 0, 0, 0, qJDD(1), t315, t229, 0, 0, qJDD(1) * t188 + 0.2e1 * t233, 0.2e1 * t195 * t261 - 0.2e1 * t263 * t269, qJDD(2) * t195 + t198 * t200, qJDD(1) * t189 - 0.2e1 * t233, qJDD(2) * t198 - t195 * t200, 0, t195 * t219 + t198 * t207, -t195 * t207 + t198 * t219, 0.2e1 * t268 * t285 - t229, -g(1) * (-pkin(1) * t196 + pkin(6) * t199) - g(2) * (pkin(1) * t199 + pkin(6) * t196) + (pkin(6) ^ 2 * t268 + pkin(1) ^ 2) * qJDD(1), t96 * t142 - t213 * t314, -t142 * t95 - t96 * t223 + (-t134 * t223 - t142 * t314) * qJD(2), t142 * qJDD(2) - t200 * t223, -t134 * t215 + t223 * t95, -qJDD(2) * t223 - t200 * t142, 0, -t179 * t95 + t124 * t223 + t103 * qJDD(2) + t315 * t183 + (-t134 * t309 + t142 * t152 + t76) * qJD(2), -t104 * qJDD(2) + t124 * t142 - t179 * t96 - t315 * t182 + (-t152 * t223 + t309 * t314 - t77) * qJD(2), t77 * t134 - t104 * t95 - t57 * t223 - t76 * t314 - t103 * t96 - t56 * t142 + (-t142 * t94 + t223 * t93) * qJD(2) - t229, t57 * t104 + t94 * t77 + t56 * t103 + t93 * t76 - t124 * t179 + t152 * t181 - g(1) * (-t179 * t196 - t192 * t199) - g(2) * (t179 * t199 - t192 * t196), t216 * t98 + t218 * t58, -t216 * t97 + t218 * t59 - t41 * t98 - t58 * t86, t259 * t98 - t260 * t58, t41 * t97 - t59 * t86, -t259 * t97 - t260 * t59, 0, t102 * t59 - t106 * t86 + t108 * t41 + t177 * t315 - t18 * t260 - t259 * t52 + t67 * t97, -t102 * t58 - t106 * t218 + t108 * t216 - t17 * t260 - t259 * t53 + t67 * t98 - t317, -t10 * t98 + t17 * t86 - t18 * t218 + t216 * t52 - t41 * t53 + t44 * t58 - t45 * t59 - t9 * t97 - t229, t9 * t53 + t45 * t17 - t10 * t52 - t44 * t18 + t67 * t108 + t102 * t106 - g(1) * (-t148 * t196 - t186 * t199) - g(2) * (-t186 * t196 + t140), -t72 * t253 + (-t27 * t98 - t58 * t72) * t197, (t193 * t72 + t197 * t70) * t58 + (t24 - t26 + (-t197 * t72 + t290) * qJD(5)) * t98, t197 * t227 - t253 * t328 - t27 * t97 + t59 * t72, t70 * t252 + (t28 * t98 - t58 * t70) * t193, -t193 * t227 - t252 * t328 - t28 * t97 - t59 * t70, t328 * t59 + t39 * t97, -g(1) * t121 - g(2) * t123 + t13 * t59 + t18 * t70 + t193 * t231 + t21 * t39 + t28 * t52 + t3 * t97 + t31 * t98 + t328 * t5, -g(1) * t120 - g(2) * t122 - t14 * t59 + t18 * t72 + t197 * t231 - t2 * t97 - t22 * t39 - t27 * t52 - t30 * t98 - t328 * t4, t21 * t27 - t22 * t28 - t4 * t70 - t5 * t72 + t225 * t58 + t317 + (qJD(5) * t224 - t193 * t2 - t197 * t3) * t98, -g(2) * t140 + t13 * t5 + t14 * t4 + t34 * t18 + t2 * t22 + t3 * t21 + t8 * t52 + (g(1) * t186 - g(2) * t316) * t199 + (-g(1) * (-t148 - t316) + g(2) * t186) * t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t251, t269 * t201, t262, t251, t261, qJDD(2), t195 * t206 - t304, g(3) * t195 + t198 * t206, 0, 0, -t282, -t132 + t313, (-t134 + t248) * qJD(2) + t208, t282, -t226, qJDD(2), -g(3) * t183 - qJD(2) * t100 - t314 * t152 + t229 * t182 + (qJDD(2) * t191 + t134 * t267) * pkin(2) + t56, g(3) * t182 + qJD(2) * t101 - t134 * t152 + t229 * t183 + (-qJDD(2) * t190 - t267 * t314) * pkin(2) - t57, (t100 + t94) * t314 + (-t101 + t93) * t134 + (-t190 * t95 - t191 * t96) * pkin(2), -t100 * t93 - t101 * t94 + (-t304 + t190 * t57 + t191 * t56 + (-qJD(1) * t152 + t229) * t195) * pkin(2), t295, t336, t337, -t295, t326, t259, -t286 * qJD(2) + qJD(4) * t339 + t105 * t86 + t128 * t259 + t320, t105 * t218 - t129 * t259 - t260 * t287 + t323, -t128 * t216 - t129 * t41 + (t287 + t44) * t86 + t339 * t218, -g(3) * t244 + t10 * t128 - t102 * t105 + t9 * t129 - t151 * t229 - t286 * t44 + t287 * t45, t334, t341, t333, t338 * t70 - t26, t222 - t300, t329, t125 * t28 - t15 * t328 + t286 * t70 + (-t331 - t254) * t197 + t212 * t193 + t322, -t125 * t27 + t16 * t328 + t286 * t72 + t212 * t197 + (-t220 + t254) * t193 + t321, t15 * t72 + t16 * t70 + (-t111 * t70 - t126 * t28 + t13 * t86 + (t126 * t72 - t13) * qJD(5)) * t197 + (t111 * t72 - t126 * t27 + t14 * t86 - t3 + (t126 * t70 - t14) * qJD(5)) * t193 + t235, t8 * t125 - t14 * t16 - t13 * t15 - g(1) * (t199 * t243 + t154) - g(2) * (t196 * t243 + t153) - g(3) * (t244 + t316) + t286 * t34 - t224 * t111 + t204 * t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226 + 0.2e1 * t324, (t134 + t248) * qJD(2) + t208, -t132 - t313, -t134 * t94 + t314 * t93 + t205, 0, 0, 0, 0, 0, 0, t240 - t272 - 0.2e1 * t325, t216 + t319, -t296 - t297, -t218 * t44 - t45 * t86 + t205 + t312, 0, 0, 0, 0, 0, 0, t222 + t300, -t197 * t328 ^ 2 + t298 - t32, (t70 * t86 + t27) * t197 + t318 + t294, t34 * t218 + t340 * t197 + (t2 - t303) * t193 - t315; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t295, t336, t337, -t295, t326, t259, t45 * qJD(2) + t320, t260 * t44 + t323, 0, 0, t334, t341, t333, t290 * t328 - t26, -t328 * t338 - t300 + t33, t329, -pkin(4) * t28 - t19 * t328 - t45 * t70 + t232 * t193 + (-t331 - t256) * t197 + t322, pkin(4) * t27 + t20 * t328 - t45 * t72 + t232 * t197 + (-t220 + t256) * t193 + t321, t19 * t72 + t20 * t70 + (-t303 + (-t28 + t284) * pkin(8)) * t197 + ((qJD(5) * t70 - t27) * pkin(8) - t340) * t193 + t235, -t8 * pkin(4) - t14 * t20 - t13 * t19 - t34 * t45 - g(1) * (-pkin(4) * t280 + t154) - g(2) * (-pkin(4) * t281 + t153) - g(3) * t316 + t204 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, -t70 ^ 2 + t72 ^ 2, t328 * t70 - t27, -t299, t328 * t72 - t28, t39, -g(1) * t122 + g(2) * t120 + t193 * t221 - t264 * t35 - t34 * t72 + t11 + t302, g(1) * t123 - g(2) * t121 + t303 + t34 * t70 + (qJD(5) * t35 - t12) * t193 + t221 * t197, 0, 0;];
tau_reg = t6;
