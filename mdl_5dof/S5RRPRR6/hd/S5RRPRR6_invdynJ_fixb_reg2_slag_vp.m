% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRR6
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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:36:16
% EndTime: 2019-12-05 18:36:23
% DurationCPUTime: 4.04s
% Computational Cost: add. (5893->406), mult. (8918->546), div. (0->0), fcn. (5766->14), ass. (0->247)
t179 = qJDD(1) + qJDD(2);
t182 = qJD(1) + qJD(2);
t187 = sin(pkin(9));
t193 = cos(qJ(5));
t194 = cos(qJ(4));
t289 = t187 * t194;
t247 = t179 * t289;
t189 = sin(qJ(5));
t190 = sin(qJ(4));
t283 = t189 * t190;
t282 = t189 * t194;
t129 = t190 * t193 + t282;
t262 = qJD(4) + qJD(5);
t332 = t262 * t129;
t31 = (t179 * t283 + t182 * t332) * t187 - t193 * t247;
t191 = sin(qJ(2));
t188 = cos(pkin(9));
t267 = qJD(3) * t188;
t195 = cos(qJ(2));
t284 = t188 * t195;
t312 = pkin(1) * qJD(1);
t215 = pkin(3) * t188 + pkin(7) * t187 + pkin(2);
t285 = t188 * t194;
t100 = qJ(3) * t285 - t190 * t215;
t330 = qJD(4) * t100;
t307 = -t190 * t267 - (-t190 * t284 + t191 * t194) * t312 - t330;
t265 = qJD(4) * t194;
t336 = (t190 * t191 + t194 * t284) * t312 - t194 * t267 + t215 * t265;
t266 = qJD(4) * t190;
t242 = t187 * t266;
t153 = pkin(8) * t242;
t334 = t153 + t307;
t286 = t188 * t190;
t253 = qJ(3) * t286;
t259 = pkin(8) * t289;
t333 = -(-t253 - t259) * qJD(4) + t336;
t186 = qJ(1) + qJ(2);
t177 = cos(t186);
t170 = g(2) * t177;
t175 = sin(t186);
t276 = g(3) * t175 + t170;
t248 = t187 * t283;
t331 = -qJD(5) * t248 - t189 * t242;
t271 = qJD(1) * t195;
t256 = pkin(1) * t271;
t222 = qJD(3) - t256;
t258 = t191 * t312;
t132 = qJ(3) * t182 + t258;
t180 = t187 ^ 2;
t329 = t132 * (qJD(4) * t188 + t180 * t182);
t326 = pkin(1) * t195;
t118 = -t215 - t326;
t166 = pkin(1) * t191 + qJ(3);
t73 = t190 * t118 + t166 * t285;
t288 = t187 * (-pkin(8) - pkin(7));
t299 = t177 * t190;
t328 = pkin(4) * t299 + t175 * t288;
t106 = t175 * t286 + t177 * t194;
t108 = -t175 * t194 + t177 * t286;
t327 = -g(2) * t106 + g(3) * t108;
t287 = t188 * t179;
t142 = -qJDD(4) + t287;
t270 = qJD(2) * t191;
t257 = pkin(1) * t270;
t274 = -qJD(1) * t257 + qJDD(1) * t326;
t240 = qJDD(3) - t274;
t65 = -t179 * t215 + t240;
t62 = t194 * t65;
t263 = qJDD(1) * t191;
t269 = qJD(2) * t195;
t90 = t179 * qJ(3) + t182 * qJD(3) + (qJD(1) * t269 + t263) * pkin(1);
t223 = -t286 * t90 + t62;
t246 = t132 * t285;
t295 = t182 * t187;
t81 = -t182 * t215 + t222;
t18 = -pkin(8) * t247 - t142 * pkin(4) + (-t246 + (pkin(8) * t295 - t81) * t190) * qJD(4) + t223;
t243 = t182 * t265;
t298 = t179 * t190;
t213 = t243 + t298;
t241 = t188 * t266;
t261 = t190 * t65 + t81 * t265 + t90 * t285;
t25 = -t132 * t241 + t261;
t19 = -pkin(8) * t187 * t213 + t25;
t264 = qJD(5) * t189;
t294 = t182 * t188;
t145 = -qJD(4) + t294;
t249 = t182 * t289;
t52 = -t132 * t286 + t194 * t81;
t44 = -pkin(8) * t249 + t52;
t34 = -pkin(4) * t145 + t44;
t290 = t187 * t190;
t260 = pkin(8) * t290;
t53 = t190 * t81 + t246;
t45 = -t182 * t260 + t53;
t4 = (qJD(5) * t34 + t19) * t193 + t189 * t18 - t45 * t264;
t325 = g(1) * t187;
t322 = t177 * pkin(2);
t321 = t179 * pkin(2);
t192 = sin(qJ(1));
t320 = t192 * pkin(1);
t196 = cos(qJ(1));
t319 = t196 * pkin(1);
t211 = t182 * t129;
t85 = t187 * t211;
t88 = -t182 * t248 + t193 * t249;
t318 = t88 * t85;
t122 = t194 * t215;
t67 = -t259 - t122 + (-qJ(3) * t190 - pkin(4)) * t188;
t80 = -t260 + t100;
t37 = -t189 * t80 + t193 * t67;
t317 = qJD(5) * t37 + t334 * t189 - t333 * t193;
t38 = t189 * t67 + t193 * t80;
t316 = -qJD(5) * t38 + t333 * t189 + t334 * t193;
t158 = pkin(1) * t269 + qJD(3);
t305 = t132 * t180;
t82 = t180 * t90;
t315 = t158 * t305 + t166 * t82;
t128 = t193 * t194 - t283;
t314 = (t262 - t294) * t128;
t313 = t188 * t211 - t332;
t311 = t189 * t45;
t310 = t193 * t45;
t268 = qJD(3) * t132;
t309 = qJ(3) * t82 + t180 * t268;
t308 = -qJ(3) * t241 - t336;
t306 = qJD(4) * t53;
t133 = -qJDD(5) + t142;
t304 = t133 * t188;
t303 = t142 * t188;
t302 = t158 * t182;
t301 = t175 * t188;
t300 = t177 * t188;
t297 = t179 * t194;
t178 = t182 ^ 2;
t296 = t180 * t178;
t293 = t182 * t190;
t292 = t182 * t194;
t291 = t187 * t179;
t281 = t190 * t194;
t279 = t331 * t182;
t278 = t276 * t187;
t277 = g(2) * t300 + g(3) * t301;
t275 = -g(2) * t175 + g(3) * t177;
t181 = t188 ^ 2;
t273 = t180 + t181;
t183 = t190 ^ 2;
t184 = t194 ^ 2;
t272 = t183 - t184;
t255 = t52 * t242 + t278;
t102 = t240 - t321;
t254 = t102 * t187 - t278;
t252 = t166 * t286;
t250 = t179 * t282;
t245 = t118 * t265 + t158 * t285 + t190 * t257;
t244 = t182 * t270;
t165 = t177 * qJ(3);
t239 = t165 - t320;
t238 = -pkin(2) * t175 + t165;
t235 = t273 * t195;
t234 = t273 * t179;
t233 = t142 + t287;
t232 = t182 * (-qJD(4) - t145);
t231 = t181 * t90 - t275 + t82;
t230 = t262 * t194;
t229 = t182 * t258;
t228 = t281 * t296;
t227 = t274 + t276;
t226 = t190 * t243;
t221 = g(2) * t196 + g(3) * t192;
t220 = -t175 * qJ(3) - t322;
t21 = t189 * t34 + t310;
t114 = t194 * t118;
t56 = -t259 + t114 + (-t166 * t190 - pkin(4)) * t188;
t63 = -t260 + t73;
t29 = -t189 * t63 + t193 * t56;
t30 = t189 * t56 + t193 * t63;
t219 = t190 * t52 - t194 * t53;
t111 = t129 * t187;
t112 = t128 * t187;
t20 = t193 * t34 - t311;
t5 = -qJD(5) * t21 + t193 * t18 - t189 * t19;
t57 = t332 * t187;
t58 = t187 * t193 * t230 + t331;
t218 = -t4 * t111 - t5 * t112 + t20 * t57 - t21 * t58 + t278;
t217 = qJD(4) * (t145 + t294);
t154 = t187 * pkin(4) * t265;
t214 = t222 * t187 + t154;
t212 = -g(2) * t108 - g(3) * t106 + t25 * t188 + t194 * t82;
t210 = -t229 - t321;
t51 = (pkin(4) * t213 + t90) * t187;
t84 = (pkin(4) * t293 + t132) * t187;
t185 = qJ(4) + qJ(5);
t174 = sin(t185);
t176 = cos(t185);
t94 = t174 * t301 + t176 * t177;
t96 = t174 * t300 - t175 * t176;
t209 = -g(2) * t96 - g(3) * t94 + t51 * t112 + t4 * t188 - t84 * t57;
t171 = pkin(4) * t194 + pkin(3);
t208 = -t171 * t300 + t177 * t288 - t322;
t172 = -pkin(2) - t326;
t207 = pkin(1) * t244 + t172 * t179;
t95 = -t177 * t174 + t176 * t301;
t97 = -t174 * t175 - t176 * t300;
t206 = -g(2) * t97 + g(3) * t95 + t51 * t111 - t188 * t5 + t84 * t58;
t205 = -t145 ^ 2 - t296;
t107 = t175 * t285 - t299;
t109 = -t175 * t190 - t177 * t285;
t26 = t223 - t306;
t204 = -g(2) * t109 + g(3) * t107 - t188 * t26 + t190 * t82 + t265 * t305;
t203 = (-g(2) * (-pkin(4) * t190 - qJ(3)) - g(3) * (-t171 * t188 - pkin(2))) * t175;
t50 = -t73 * qJD(4) - t158 * t286 + t194 * t257;
t201 = -g(2) * t95 - g(3) * t97 + t176 * t325 + t84 * t85 - t4;
t200 = -g(2) * t94 + g(3) * t96 + t174 * t325 - t84 * t88 + t5;
t199 = (g(2) * qJ(3) + g(3) * t215) * t175 + t215 * t170;
t163 = pkin(4) * t290;
t162 = t181 * t179;
t161 = t180 * t179;
t135 = -qJD(5) + t145;
t134 = 0.2e1 * t187 * t287;
t126 = qJ(3) * t187 + t163;
t125 = -pkin(2) * t182 + t222;
t115 = t166 * t187 + t163;
t101 = t158 * t187 + t154;
t99 = -t122 - t253;
t92 = (t179 * t184 - 0.2e1 * t226) * t180;
t91 = (t179 * t183 + 0.2e1 * t226) * t180;
t72 = t114 - t252;
t64 = 0.2e1 * (qJD(4) * t182 * t272 - t179 * t281) * t180;
t49 = -t166 * t241 + t245;
t47 = (t190 * t233 + t194 * t217) * t187;
t46 = (t190 * t217 - t194 * t233) * t187;
t43 = t153 + t50;
t42 = (-t252 - t259) * qJD(4) + t245;
t33 = -t85 ^ 2 + t88 ^ 2;
t32 = (t250 + (t182 * t230 + t298) * t193) * t187 + t279;
t28 = -t88 * t135 + (-t250 + (-t262 * t292 - t298) * t193) * t187 - t279;
t27 = -t85 * t135 - t31;
t23 = t193 * t44 - t311;
t22 = -t189 * t44 - t310;
t14 = t111 * t32 + t58 * t85;
t13 = -t112 * t31 - t57 * t88;
t12 = t111 * t133 + t135 * t58 + t188 * t32;
t11 = -t112 * t133 + t135 * t57 + t188 * t31;
t8 = -qJD(5) * t30 - t189 * t42 + t193 * t43;
t7 = qJD(5) * t29 + t189 * t43 + t193 * t42;
t6 = t111 * t31 - t112 * t32 + t57 * t85 - t58 * t88;
t1 = [0, 0, 0, 0, 0, qJDD(1), t221, -g(2) * t192 + g(3) * t196, 0, 0, 0, 0, 0, 0, 0, t179, (t179 * t195 - t244) * pkin(1) + t227, ((-qJDD(1) - t179) * t191 + (-qJD(1) - t182) * t269) * pkin(1) + t275, 0, (t221 + (t191 ^ 2 + t195 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t161, t134, 0, t162, 0, 0, (-t102 - t207) * t188 + t277, t187 * t207 + t254, t166 * t234 + t273 * t302 + t231, t102 * t172 + t125 * t257 - g(2) * (t220 - t319) - g(3) * (t238 - t320) + (t132 * t158 + t166 * t90) * t181 + t315, t92, t64, t46, t91, t47, t303, -t142 * t72 - t145 * t50 + (t158 * t293 + t166 * t213) * t180 + t204, t142 * t73 + t145 * t49 + ((t166 * t179 + t302) * t194 + (-t166 * t182 - t132) * t266) * t180 + t212, ((-t179 * t73 - t25 + (qJD(4) * t72 - t49) * t182) * t190 + (-t179 * t72 - t182 * t50 - t26 + (-t182 * t73 - t53) * qJD(4)) * t194) * t187 + t255, g(2) * t319 - g(3) * t239 + t25 * t73 + t26 * t72 + t53 * t49 + t52 * t50 + t199 + t315, t13, t6, t11, t14, t12, t304, t101 * t85 + t115 * t32 - t133 * t29 - t135 * t8 + t206, t101 * t88 - t115 * t31 + t133 * t30 + t135 * t7 + t209, t29 * t31 - t30 * t32 - t7 * t85 - t8 * t88 + t218, t4 * t30 + t21 * t7 + t5 * t29 + t20 * t8 + t51 * t115 + t84 * t101 - g(2) * (t208 - t319) - g(3) * (t239 + t328) + t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, t227 + t229, (-t263 + (-qJD(2) + t182) * t271) * pkin(1) + t275, 0, 0, t161, t134, 0, t162, 0, 0, (-t102 - t210) * t188 + t277, t187 * t210 + t254, qJ(3) * t234 + (qJD(3) * t273 - t235 * t312) * t182 + t231, -t102 * pkin(2) - g(2) * t220 - g(3) * t238 + (qJ(3) * t90 + t268) * t181 + (-t125 * t191 - t132 * t235) * t312 + t309, t92, t64, t46, t91, t47, t303, -t142 * t99 - t307 * t145 + (qJ(3) * t298 + (qJ(3) * t265 + t190 * t222) * t182) * t180 + t204, t100 * t142 + t308 * t145 + (qJ(3) * t297 - t132 * t266 + (-qJ(3) * t266 + t194 * t222) * t182) * t180 + t212, ((-t100 * t179 - t25) * t190 + (-t179 * t99 - t26 - t306) * t194 + ((-t307 - t330) * t194 + (qJD(4) * t99 - t308) * t190) * t182) * t187 + t255, -g(3) * t165 + t25 * t100 - t256 * t305 + t26 * t99 + t307 * t52 + t308 * t53 + t199 + t309, t13, t6, t11, t14, t12, t304, t126 * t32 - t133 * t37 - t316 * t135 + t214 * t85 + t206, -t126 * t31 + t133 * t38 + t317 * t135 + t214 * t88 + t209, t31 * t37 - t316 * t88 - t317 * t85 - t32 * t38 + t218, t4 * t38 + t5 * t37 + t51 * t126 - g(2) * t208 - g(3) * (t165 + t328) + t214 * t84 + t317 * t21 + t316 * t20 + t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t287, t291, -t273 * t178, -t132 * t182 * t273 + t102 - t276, 0, 0, 0, 0, 0, 0, -t194 * t142 + t190 * t205, t190 * t142 + t194 * t205, (-t183 - t184) * t291, t25 * t190 + t26 * t194 - t219 * qJD(4) + (t188 * t219 - t305) * t182 - t276, 0, 0, 0, 0, 0, 0, -t128 * t133 - t313 * t135 - t85 * t295, t129 * t133 + t314 * t135 - t88 * t295, t128 * t31 - t129 * t32 - t313 * t88 - t314 * t85, t128 * t5 + t129 * t4 + t313 * t20 + t314 * t21 - t84 * t295 - t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, -t272 * t296, (t190 * t232 + t297) * t187, -t228, (t194 * t232 - t298) * t187, -t142, -t145 * t53 + t62 - t194 * t329 + (-qJD(4) * t81 - t188 * t90 + t325) * t190 + t327, g(1) * t289 - g(2) * t107 - g(3) * t109 - t145 * t52 + t190 * t329 - t261, 0, 0, t318, t33, t27, -t318, t28, -t133, t135 * t22 + (-t133 * t193 + t135 * t264 - t249 * t85) * pkin(4) + t200, -t135 * t23 + (qJD(5) * t135 * t193 + t133 * t189 - t249 * t88) * pkin(4) + t201, (t21 + t22) * t88 + (-t20 + t23) * t85 + (-t189 * t32 + t193 * t31 + (t189 * t88 - t193 * t85) * qJD(5)) * pkin(4), -t20 * t22 - t21 * t23 + (t4 * t189 + t5 * t193 + (g(1) * t190 - t292 * t84) * t187 + (-t189 * t20 + t193 * t21) * qJD(5) + t327) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t318, t33, t27, -t318, t28, -t133, -t135 * t21 + t200, -t135 * t20 + t201, 0, 0;];
tau_reg = t1;
