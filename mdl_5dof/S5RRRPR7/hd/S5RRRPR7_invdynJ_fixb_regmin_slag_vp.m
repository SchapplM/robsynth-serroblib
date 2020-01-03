% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:17:20
% EndTime: 2019-12-31 21:17:30
% DurationCPUTime: 4.12s
% Computational Cost: add. (4824->390), mult. (10920->531), div. (0->0), fcn. (7999->14), ass. (0->222)
t193 = qJD(2) + qJD(3);
t200 = sin(qJ(3));
t204 = cos(qJ(2));
t301 = cos(qJ(3));
t253 = qJD(1) * t301;
t201 = sin(qJ(2));
t268 = qJD(1) * t201;
t312 = -t200 * t268 + t204 * t253;
t314 = t312 * t193;
t125 = qJD(5) - t312;
t276 = t200 * t204;
t133 = -qJD(1) * t276 - t201 * t253;
t197 = sin(pkin(9));
t198 = cos(pkin(9));
t110 = t198 * t133 - t197 * t193;
t112 = t197 * t133 + t198 * t193;
t199 = sin(qJ(5));
t203 = cos(qJ(5));
t307 = t199 * t110 + t203 * t112;
t313 = t125 * t307;
t311 = qJD(5) - t125;
t141 = t199 * t197 - t203 * t198;
t294 = t125 * t141;
t142 = t203 * t197 + t199 * t198;
t130 = t142 * qJD(5);
t310 = -t142 * t312 + t130;
t233 = t203 * t110 - t199 * t112;
t309 = t125 * t233;
t196 = qJ(2) + qJ(3);
t188 = sin(t196);
t205 = cos(qJ(1));
t283 = t188 * t205;
t202 = sin(qJ(1));
t284 = t188 * t202;
t189 = cos(t196);
t300 = g(3) * t189;
t308 = g(1) * t283 + g(2) * t284 - t300;
t302 = pkin(7) + pkin(6);
t252 = qJD(3) * t301;
t167 = pkin(2) * t252 + qJD(4);
t157 = t302 * t204;
t147 = qJD(1) * t157;
t134 = t200 * t147;
t156 = t302 * t201;
t145 = qJD(1) * t156;
t102 = -t301 * t145 - t134;
t94 = -t133 * pkin(3) - qJ(4) * t312;
t85 = pkin(2) * t268 + t94;
t50 = t198 * t102 + t197 * t85;
t306 = -t198 * t167 + t50;
t296 = qJD(2) * pkin(2);
t138 = -t145 + t296;
t97 = t301 * t138 - t134;
t52 = t197 * t94 + t198 * t97;
t305 = -t198 * qJD(4) + t52;
t241 = g(1) * t202 - g(2) * t205;
t304 = t241 * t188;
t135 = t301 * t147;
t101 = -t200 * t145 + t135;
t267 = qJD(3) * t200;
t244 = pkin(2) * t267 - t101;
t270 = t189 * pkin(3) + t188 * qJ(4);
t191 = qJDD(2) + qJDD(3);
t303 = -t191 * pkin(3) + qJDD(4);
t248 = qJDD(1) * t301;
t261 = t204 * qJDD(1);
t72 = t200 * t261 + t201 * t248 + t314;
t61 = -t198 * t191 + t197 * t72;
t62 = t197 * t191 + t198 * t72;
t15 = -qJD(5) * t233 + t199 * t62 + t203 * t61;
t182 = g(3) * t188;
t298 = t198 * pkin(4);
t190 = t198 * pkin(8);
t297 = t204 * pkin(2);
t185 = pkin(1) + t297;
t263 = qJD(1) * qJD(2);
t251 = t201 * t263;
t126 = pkin(2) * t251 - qJDD(1) * t185;
t144 = t301 * t201 + t276;
t105 = t193 * t144;
t262 = t201 * qJDD(1);
t238 = t200 * t262 - t204 * t248;
t73 = qJD(1) * t105 + t238;
t24 = t73 * pkin(3) - t72 * qJ(4) + t133 * qJD(4) + t126;
t250 = t204 * t263;
t107 = qJDD(2) * pkin(2) + t302 * (-t250 - t262);
t113 = t302 * (-t251 + t261);
t211 = t200 * t107 + t301 * t113 + t138 * t252 - t147 * t267;
t32 = t191 * qJ(4) + t193 * qJD(4) + t211;
t10 = t197 * t24 + t198 * t32;
t224 = -t200 * t201 + t301 * t204;
t104 = t193 * t224;
t259 = t201 * t296;
t48 = t105 * pkin(3) - t104 * qJ(4) - t144 * qJD(4) + t259;
t255 = qJD(2) * t302;
t146 = t201 * t255;
t148 = t204 * t255;
t225 = -t301 * t156 - t200 * t157;
t67 = t225 * qJD(3) - t301 * t146 - t200 * t148;
t21 = t197 * t48 + t198 * t67;
t155 = t185 * qJD(1);
t81 = -pkin(3) * t312 + t133 * qJ(4) - t155;
t98 = t200 * t138 + t135;
t86 = t193 * qJ(4) + t98;
t42 = t197 * t81 + t198 * t86;
t8 = t10 * t198;
t295 = t197 * t73;
t115 = -t200 * t156 + t301 * t157;
t96 = -pkin(3) * t224 - t144 * qJ(4) - t185;
t57 = t198 * t115 + t197 * t96;
t84 = -t193 * pkin(3) + qJD(4) - t97;
t292 = -t167 + t84;
t291 = t104 * t197;
t290 = t125 * t133;
t289 = t312 * t197;
t288 = t312 * t198;
t287 = t133 * t312;
t286 = t144 * t197;
t285 = t144 * t198;
t282 = t189 * t197;
t281 = t189 * t202;
t280 = t189 * t205;
t279 = t198 * qJ(4);
t179 = t200 * pkin(2) + qJ(4);
t277 = t198 * t179;
t192 = pkin(9) + qJ(5);
t186 = sin(t192);
t275 = t202 * t186;
t187 = cos(t192);
t274 = t202 * t187;
t273 = t205 * t186;
t272 = t205 * t187;
t271 = -qJD(4) + t84;
t194 = t201 ^ 2;
t269 = -t204 ^ 2 + t194;
t266 = qJD(5) * t199;
t265 = qJD(5) * t203;
t260 = pkin(8) * t289;
t257 = g(1) * t280 + g(2) * t281 + t182;
t9 = -t197 * t32 + t198 * t24;
t4 = t73 * pkin(4) - t62 * pkin(8) + t9;
t5 = -t61 * pkin(8) + t10;
t256 = -t199 * t5 + t203 * t4;
t20 = -t197 * t67 + t198 * t48;
t41 = -t197 * t86 + t198 * t81;
t51 = -t197 * t97 + t198 * t94;
t49 = -t197 * t102 + t198 * t85;
t56 = -t197 * t115 + t198 * t96;
t247 = t301 * t107 - t200 * t113 - t138 * t267 - t147 * t252;
t184 = -t301 * pkin(2) - pkin(3);
t121 = pkin(4) * t289;
t246 = -t121 + t244;
t245 = -t133 * pkin(4) - pkin(8) * t288;
t243 = -pkin(2) * t201 - pkin(3) * t188;
t242 = g(1) * t205 + g(2) * t202;
t240 = t199 * t4 + t203 * t5;
t239 = -t9 * t197 + t8;
t237 = -t41 * t197 + t42 * t198;
t25 = -pkin(4) * t312 + pkin(8) * t110 + t41;
t26 = pkin(8) * t112 + t42;
t236 = t199 * t26 - t203 * t25;
t7 = t199 * t25 + t203 * t26;
t38 = -pkin(4) * t224 - pkin(8) * t285 + t56;
t47 = -pkin(8) * t286 + t57;
t235 = -t199 * t47 + t203 * t38;
t234 = t199 * t38 + t203 * t47;
t232 = t41 * t288 + t42 * t289 - t257 + t8;
t231 = t185 + t270;
t230 = t242 * t188;
t229 = t241 * t189;
t228 = -0.2e1 * pkin(1) * t263 - pkin(6) * qJDD(2);
t136 = (-pkin(8) - t179) * t197;
t227 = -qJD(5) * t136 - t260 + t306;
t137 = t190 + t277;
t226 = qJD(5) * t137 + t197 * t167 + t245 + t49;
t153 = (-pkin(8) - qJ(4)) * t197;
t223 = -qJD(5) * t153 - t260 + t305;
t154 = t190 + t279;
t222 = t197 * qJD(4) + qJD(5) * t154 + t245 + t51;
t14 = t110 * t266 + t112 * t265 - t199 * t61 + t203 * t62;
t36 = -t247 + t303;
t219 = t247 + t308;
t207 = qJD(2) ^ 2;
t218 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t207 + t241;
t208 = qJD(1) ^ 2;
t217 = pkin(1) * t208 - pkin(6) * qJDD(1) + t242;
t216 = t84 * t104 + t36 * t144 - t242;
t215 = t41 * t133 + (-t36 + t308) * t198;
t214 = -t155 * t133 + t219;
t19 = t61 * pkin(4) + t36;
t55 = -pkin(4) * t112 + t84;
t213 = -t236 * t133 + t19 * t141 + t308 * t187 + t310 * t55;
t212 = g(3) * t282 - t42 * t133 + (-t230 + t36) * t197;
t68 = t115 * qJD(3) - t200 * t146 + t301 * t148;
t210 = -t7 * t133 + t19 * t142 - t294 * t55 + (-t230 + t300) * t186;
t209 = t155 * t312 - t211 + t257;
t180 = -pkin(3) - t298;
t159 = qJ(4) * t280;
t158 = qJ(4) * t281;
t152 = t184 - t298;
t119 = t189 * t272 + t275;
t118 = -t189 * t273 + t274;
t117 = -t189 * t274 + t273;
t116 = t189 * t275 + t272;
t90 = t141 * t144;
t89 = t142 * t144;
t83 = pkin(4) * t286 - t225;
t74 = t133 ^ 2 - t312 ^ 2;
t70 = qJDD(5) + t73;
t69 = t121 + t98;
t54 = -t238 + (-qJD(1) * t144 - t133) * t193;
t53 = t72 - t314;
t40 = pkin(4) * t291 + t68;
t30 = t104 * t142 + t265 * t285 - t266 * t286;
t29 = -t104 * t141 - t144 * t130;
t16 = -pkin(8) * t291 + t21;
t13 = t105 * pkin(4) - t104 * t190 + t20;
t12 = -t125 * t310 + t133 * t307 - t141 * t70;
t11 = -t294 * t125 - t133 * t233 + t142 * t70;
t2 = t14 * t142 + t233 * t294;
t1 = -t14 * t141 - t142 * t15 + t233 * t310 - t294 * t307;
t3 = [qJDD(1), t241, t242, t194 * qJDD(1) + 0.2e1 * t201 * t250, 0.2e1 * t201 * t261 - 0.2e1 * t269 * t263, qJDD(2) * t201 + t207 * t204, qJDD(2) * t204 - t207 * t201, 0, t201 * t228 + t204 * t218, -t201 * t218 + t204 * t228, -t133 * t104 + t72 * t144, t104 * t312 + t133 * t105 - t144 * t73 + t224 * t72, t104 * t193 + t144 * t191, -t105 * t193 + t191 * t224, 0, -t155 * t105 - t126 * t224 - t185 * t73 + t191 * t225 - t68 * t193 - t259 * t312 + t229, -t155 * t104 - t115 * t191 + t126 * t144 - t133 * t259 - t185 * t72 - t67 * t193 - t304, t41 * t105 - t112 * t68 + t197 * t216 + t198 * t229 - t20 * t312 - t224 * t9 - t225 * t61 + t56 * t73, t10 * t224 - t42 * t105 - t110 * t68 + t198 * t216 + t21 * t312 - t225 * t62 - t241 * t282 - t57 * t73, t21 * t112 + t20 * t110 - t56 * t62 - t57 * t61 + t304 + (-t10 * t197 - t198 * t9) * t144 + (-t197 * t42 - t198 * t41) * t104, t10 * t57 - t36 * t225 + t41 * t20 + t42 * t21 + t9 * t56 + t84 * t68 + (-g(1) * t302 - g(2) * t231) * t205 + (g(1) * t231 - g(2) * t302) * t202, -t14 * t90 - t233 * t29, -t14 * t89 + t90 * t15 + t233 * t30 + t29 * t307, -t105 * t233 + t29 * t125 - t14 * t224 - t90 * t70, t105 * t307 - t30 * t125 + t15 * t224 - t89 * t70, t125 * t105 - t224 * t70, (t203 * t13 - t199 * t16) * t125 + t235 * t70 - t256 * t224 - t236 * t105 - t40 * t307 + t83 * t15 + t19 * t89 + t55 * t30 - g(1) * t117 - g(2) * t119 + (-t125 * t234 + t224 * t7) * qJD(5), -(t199 * t13 + t203 * t16) * t125 - t234 * t70 + t240 * t224 - t7 * t105 - t40 * t233 + t83 * t14 - t19 * t90 + t55 * t29 - g(1) * t116 - g(2) * t118 + (-t125 * t235 - t224 * t236) * qJD(5); 0, 0, 0, -t201 * t208 * t204, t269 * t208, t262, t261, qJDD(2), -g(3) * t204 + t201 * t217, g(3) * t201 + t204 * t217, t287, t74, t53, t54, t191, t101 * t193 + (t301 * t191 - t193 * t267 + t268 * t312) * pkin(2) + t214, t102 * t193 + (t133 * t268 - t191 * t200 - t193 * t252) * pkin(2) + t209, -t179 * t295 + t184 * t61 - t244 * t112 - (t292 * t197 - t49) * t312 + t215, -t73 * t277 + t184 * t62 - t244 * t110 - (t292 * t198 + t50) * t312 + t212, -t61 * t277 - t49 * t110 - t306 * t112 + (-t110 * t167 + t179 * t62 - t9) * t197 + t232, t36 * t184 - t42 * t50 - t41 * t49 - g(1) * (t205 * t243 + t159) - g(2) * (t202 * t243 + t158) - g(3) * (t270 + t297) + t244 * t84 + t239 * t179 + t237 * t167, t2, t1, t11, t12, t290, (t203 * t136 - t199 * t137) * t70 + t152 * t15 - t246 * t307 + (t199 * t227 - t203 * t226) * t125 + t213, -(t199 * t136 + t203 * t137) * t70 + t152 * t14 - t246 * t233 + (t199 * t226 + t203 * t227) * t125 + t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t287, t74, t53, t54, t191, t98 * t193 + t214, t97 * t193 + t209, -qJ(4) * t295 - pkin(3) * t61 + t98 * t112 - (t271 * t197 - t51) * t312 + t215, -t73 * t279 - pkin(3) * t62 + t98 * t110 - (t271 * t198 + t52) * t312 + t212, -t61 * t279 - t51 * t110 - t305 * t112 + (qJ(4) * t62 - qJD(4) * t110 - t9) * t197 + t232, -t36 * pkin(3) - t42 * t52 - t41 * t51 - t84 * t98 - g(1) * (-pkin(3) * t283 + t159) - g(2) * (-pkin(3) * t284 + t158) - g(3) * t270 + t237 * qJD(4) + t239 * qJ(4), t2, t1, t11, t12, t290, (t203 * t153 - t199 * t154) * t70 + t180 * t15 + t69 * t307 + (t199 * t223 - t203 * t222) * t125 + t213, -(t199 * t153 + t203 * t154) * t70 + t180 * t14 + t69 * t233 + (t199 * t222 + t203 * t223) * t125 + t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110 * t312 + t61, -t112 * t312 + t62, -t110 ^ 2 - t112 ^ 2, -t41 * t110 - t42 * t112 - t219 + t303, 0, 0, 0, 0, 0, t15 - t309, t14 + t313; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t233 * t307, t233 ^ 2 - t307 ^ 2, t14 - t313, -t15 - t309, t70, -g(1) * t118 + g(2) * t116 + t186 * t182 + t55 * t233 - t311 * t7 + t256, g(1) * t119 - g(2) * t117 + t187 * t182 + t311 * t236 - t307 * t55 - t240;];
tau_reg = t3;
