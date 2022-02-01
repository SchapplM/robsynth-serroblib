% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPR3
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:43:16
% EndTime: 2022-01-20 11:43:23
% DurationCPUTime: 3.13s
% Computational Cost: add. (6325->393), mult. (9824->491), div. (0->0), fcn. (6824->16), ass. (0->246)
t238 = cos(qJ(3));
t215 = t238 * qJD(4);
t234 = sin(qJ(3));
t232 = -qJ(4) - pkin(7);
t282 = qJD(3) * t232;
t151 = t234 * t282 + t215;
t152 = -t234 * qJD(4) + t238 * t282;
t230 = sin(pkin(9));
t231 = cos(pkin(9));
t162 = t230 * t238 + t231 * t234;
t239 = cos(qJ(2));
t303 = qJD(1) * t239;
t293 = pkin(1) * t303;
t324 = -t230 * t151 + t231 * t152 + t162 * t293;
t313 = t231 * t238;
t315 = t230 * t234;
t161 = -t313 + t315;
t323 = t231 * t151 + t230 * t152 + t161 * t293;
t233 = sin(qJ(5));
t237 = cos(qJ(5));
t225 = qJD(1) + qJD(2);
t292 = t225 * t313;
t130 = t225 * t315 - t292;
t132 = t162 * t225;
t261 = t233 * t130 - t237 * t132;
t154 = t162 * qJD(3);
t222 = qJDD(1) + qJDD(2);
t311 = t238 * t222;
t312 = t234 * t222;
t263 = t230 * t312 - t231 * t311;
t79 = t154 * t225 + t263;
t300 = qJD(3) * t234;
t288 = t225 * t300;
t250 = t162 * t222 - t230 * t288;
t299 = qJD(3) * t238;
t287 = t225 * t299;
t80 = t231 * t287 + t250;
t246 = qJD(5) * t261 - t233 * t80 - t237 * t79;
t224 = qJD(3) + qJD(5);
t328 = t261 * t224;
t361 = t246 - t328;
t297 = qJD(5) * t237;
t298 = qJD(5) * t233;
t252 = -t130 * t297 - t132 * t298 - t233 * t79 + t237 * t80;
t74 = -t237 * t130 - t233 * t132;
t327 = t74 * t224;
t360 = t252 - t327;
t333 = t261 ^ 2;
t334 = t74 ^ 2;
t359 = t333 - t334;
t332 = t74 * t261;
t155 = t161 * qJD(3);
t342 = t155 * pkin(8);
t358 = t342 + t324;
t147 = t154 * pkin(8);
t357 = -t147 + t323;
t227 = t234 ^ 2;
t228 = t238 ^ 2;
t304 = t227 + t228;
t351 = t239 * t304;
t356 = t225 * t351;
t235 = sin(qJ(2));
t296 = qJDD(1) * t235;
t301 = qJD(2) * t239;
t140 = t222 * pkin(7) + (qJD(1) * t301 + t296) * pkin(1);
t278 = t304 * t140;
t329 = pkin(1) * qJD(1);
t294 = t235 * t329;
t335 = t239 * pkin(1);
t306 = -qJD(2) * t294 + qJDD(1) * t335;
t340 = t222 * pkin(2);
t139 = -t306 - t340;
t229 = qJ(1) + qJ(2);
t217 = cos(t229);
t203 = g(2) * t217;
t355 = t139 + t203;
t216 = sin(t229);
t204 = g(1) * t216;
t349 = t204 - t203;
t226 = qJ(3) + pkin(9);
t214 = qJ(5) + t226;
t196 = sin(t214);
t197 = cos(t214);
t317 = t197 * t217;
t318 = t197 * t216;
t253 = qJ(4) * t222 + qJD(4) * t225 + t140;
t174 = t225 * pkin(7) + t294;
t276 = qJ(4) * t225 + t174;
t259 = qJD(3) * t276;
t55 = qJDD(3) * pkin(3) - t234 * t253 - t238 * t259;
t58 = -t234 * t259 + t238 * t253;
t26 = -t230 * t58 + t231 * t55;
t16 = qJDD(3) * pkin(4) - t80 * pkin(8) + t26;
t27 = t230 * t55 + t231 * t58;
t17 = -t79 * pkin(8) + t27;
t343 = t132 * pkin(8);
t118 = t276 * t238;
t105 = t230 * t118;
t117 = t276 * t234;
t111 = qJD(3) * pkin(3) - t117;
t61 = t231 * t111 - t105;
t43 = qJD(3) * pkin(4) - t343 + t61;
t344 = t130 * pkin(8);
t314 = t231 * t118;
t62 = t230 * t111 + t314;
t45 = t62 - t344;
t4 = (qJD(5) * t43 + t17) * t237 + t233 * t16 - t45 * t298;
t336 = t238 * pkin(3);
t207 = pkin(2) + t336;
t129 = -t207 * t225 + qJD(4) - t293;
t81 = t130 * pkin(4) + t129;
t354 = g(1) * t317 + g(2) * t318 + g(3) * t196 - t81 * t74 - t4;
t319 = t196 * t217;
t320 = t196 * t216;
t19 = t233 * t43 + t237 * t45;
t5 = -qJD(5) * t19 + t237 * t16 - t233 * t17;
t353 = g(1) * t319 + g(2) * t320 - g(3) * t197 + t261 * t81 + t5;
t339 = t225 * pkin(2);
t175 = -t293 - t339;
t352 = t174 * t351 + t175 * t235;
t350 = g(1) * t217 + g(2) * t216;
t348 = t132 ^ 2;
t347 = pkin(3) * t230;
t236 = sin(qJ(1));
t346 = g(1) * t236;
t345 = g(3) * t238;
t341 = t162 * pkin(8);
t338 = t234 * pkin(3);
t337 = t236 * pkin(1);
t186 = t232 * t234;
t218 = t238 * qJ(4);
t187 = t238 * pkin(7) + t218;
t109 = t231 * t186 - t230 * t187;
t84 = t109 - t341;
t110 = t230 * t186 + t231 * t187;
t158 = t161 * pkin(8);
t85 = -t158 + t110;
t38 = -t233 * t85 + t237 * t84;
t331 = qJD(5) * t38 + t358 * t233 + t357 * t237;
t39 = t233 * t84 + t237 * t85;
t330 = -qJD(5) * t39 - t357 * t233 + t358 * t237;
t206 = t235 * pkin(1) + pkin(7);
t310 = -qJ(4) - t206;
t272 = qJD(3) * t310;
t295 = pkin(1) * t301;
t100 = (-qJD(4) - t295) * t234 + t238 * t272;
t99 = t234 * t272 + t238 * t295 + t215;
t57 = t230 * t100 + t231 * t99;
t198 = t231 * pkin(3) + pkin(4);
t149 = t237 * t198 - t233 * t347;
t63 = t230 * t117 - t314;
t46 = t63 + t344;
t64 = -t231 * t117 - t105;
t47 = t64 - t343;
t326 = t149 * qJD(5) - t233 * t46 - t237 * t47;
t150 = t233 * t198 + t237 * t347;
t325 = -t150 * qJD(5) + t233 * t47 - t237 * t46;
t322 = t132 * t130;
t316 = t225 * t234;
t159 = t310 * t234;
t160 = t238 * t206 + t218;
t90 = t230 * t159 + t231 * t160;
t309 = t175 * t300 + t238 * t204;
t308 = t217 * pkin(2) + t216 * pkin(7);
t305 = t227 - t228;
t302 = qJD(2) * t235;
t210 = pkin(3) * t300;
t220 = t225 ^ 2;
t291 = t234 * t220 * t238;
t290 = t175 * t299 + t355 * t234;
t289 = t225 * t302;
t213 = cos(t226);
t286 = pkin(4) * t213 + t336;
t114 = t154 * pkin(4) + t210;
t56 = t231 * t100 - t230 * t99;
t277 = t304 * t222;
t89 = t231 * t159 - t230 * t160;
t166 = pkin(2) + t286;
t223 = pkin(8) - t232;
t275 = t217 * t166 + t216 * t223;
t274 = -t216 * t166 + t223 * t217;
t273 = t217 * t207 - t216 * t232;
t271 = t225 * t294;
t270 = -t350 + t278;
t269 = t306 + t349;
t268 = t234 * t287;
t267 = t114 - t294;
t266 = g(1) * (-t216 * pkin(2) + t217 * pkin(7));
t240 = cos(qJ(1));
t264 = -g(2) * t240 + t346;
t67 = t89 - t341;
t68 = -t158 + t90;
t34 = -t233 * t68 + t237 * t67;
t35 = t233 * t67 + t237 * t68;
t18 = -t233 * t45 + t237 * t43;
t48 = t233 * t154 + t237 * t155 + t161 * t297 + t162 * t298;
t94 = -t233 * t161 + t237 * t162;
t49 = qJD(5) * t94 + t237 * t154 - t233 * t155;
t93 = t237 * t161 + t233 * t162;
t262 = t18 * t48 - t19 * t49 - t4 * t93 - t5 * t94 - t350;
t260 = -t216 * t207 - t217 * t232;
t125 = t161 * pkin(4) - t207;
t258 = -t62 * t154 + t61 * t155 - t27 * t161 - t26 * t162 - t350;
t86 = pkin(3) * t288 - t207 * t222 + qJDD(4) - t306;
t44 = t79 * pkin(4) + t86;
t257 = -g(1) * t320 + g(2) * t319 + t44 * t94 - t81 * t48;
t212 = sin(t226);
t256 = -t129 * t155 + t86 * t162 - t349 * t212;
t255 = g(1) * t318 - g(2) * t317 + t44 * t93 + t81 * t49;
t254 = t129 * t154 + t86 * t161 + t349 * t213;
t251 = -t175 * t225 - t140 + t350;
t241 = qJD(3) ^ 2;
t249 = pkin(7) * t241 - t271 - t340;
t208 = -pkin(2) - t335;
t248 = pkin(1) * t289 + t206 * t241 + t208 * t222;
t247 = -pkin(7) * qJDD(3) + (t293 - t339) * qJD(3);
t245 = -qJDD(3) * t206 + (t208 * t225 - t295) * qJD(3);
t221 = qJDD(3) + qJDD(5);
t219 = t240 * pkin(1);
t211 = pkin(1) * t302;
t184 = -t207 - t335;
t182 = qJDD(3) * t238 - t241 * t234;
t181 = qJDD(3) * t234 + t241 * t238;
t173 = t211 + t210;
t142 = t228 * t222 - 0.2e1 * t268;
t141 = t227 * t222 + 0.2e1 * t268;
t124 = t130 ^ 2;
t113 = t125 - t335;
t104 = -0.2e1 * t305 * t225 * qJD(3) + 0.2e1 * t234 * t311;
t103 = t114 + t211;
t98 = pkin(3) * t316 + t132 * pkin(4);
t92 = -t155 * qJD(3) + t162 * qJDD(3);
t91 = -t154 * qJD(3) - t161 * qJDD(3);
t41 = -t147 + t57;
t40 = t56 + t342;
t37 = -t132 * t155 + t80 * t162;
t36 = t130 * t154 + t79 * t161;
t29 = -t93 * t221 - t49 * t224;
t28 = t94 * t221 - t48 * t224;
t15 = t155 * t130 - t132 * t154 - t80 * t161 - t162 * t79;
t9 = -qJD(5) * t35 - t233 * t41 + t237 * t40;
t8 = qJD(5) * t34 + t233 * t40 + t237 * t41;
t7 = -t246 * t93 - t49 * t74;
t6 = t252 * t94 + t261 * t48;
t1 = t246 * t94 - t252 * t93 + t261 * t49 - t48 * t74;
t2 = [0, 0, 0, 0, 0, qJDD(1), t264, g(1) * t240 + g(2) * t236, 0, 0, 0, 0, 0, 0, 0, t222, (t222 * t239 - t289) * pkin(1) + t269, ((-qJDD(1) - t222) * t235 + (-qJD(1) - t225) * t301) * pkin(1) + t350, 0, (t264 + (t235 ^ 2 + t239 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t141, t104, t181, t142, t182, 0, t245 * t234 + (-t248 - t355) * t238 + t309, t245 * t238 + (t248 - t204) * t234 + t290, pkin(1) * qJD(2) * t356 + t206 * t277 + t270, t139 * t208 - t266 - g(2) * (t219 + t308) + t206 * t278 + (t352 * qJD(2) + t346) * pkin(1), t37, t15, t92, t36, t91, 0, t56 * qJD(3) + t89 * qJDD(3) + t173 * t130 + t184 * t79 + t254, -t57 * qJD(3) - t90 * qJDD(3) + t173 * t132 + t184 * t80 + t256, -t57 * t130 - t56 * t132 - t90 * t79 - t89 * t80 + t258, t27 * t90 + t62 * t57 + t26 * t89 + t61 * t56 + t86 * t184 + t129 * t173 - g(1) * (t260 - t337) - g(2) * (t219 + t273), t6, t1, t28, t7, t29, 0, -t103 * t74 - t113 * t246 + t34 * t221 + t9 * t224 + t255, -t103 * t261 + t113 * t252 - t35 * t221 - t8 * t224 + t257, t246 * t35 - t252 * t34 + t261 * t9 + t74 * t8 + t262, t4 * t35 + t19 * t8 + t5 * t34 + t18 * t9 + t44 * t113 + t81 * t103 - g(1) * (t274 - t337) - g(2) * (t219 + t275); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, t269 + t271, (-t296 + (-qJD(2) + t225) * t303) * pkin(1) + t350, 0, 0, t141, t104, t181, t142, t182, 0, t247 * t234 + (-t249 - t355) * t238 + t309, t247 * t238 + (t249 - t204) * t234 + t290, pkin(7) * t277 - t329 * t356 + t270, -t139 * pkin(2) + pkin(7) * t278 - g(2) * t308 - t352 * t329 - t266, t37, t15, t92, t36, t91, 0, -t130 * t294 + t109 * qJDD(3) - t207 * t79 + (t130 * t338 + t324) * qJD(3) + t254, -t132 * t294 - t110 * qJDD(3) - t207 * t80 + (t132 * t338 - t323) * qJD(3) + t256, -t109 * t80 - t110 * t79 - t323 * t130 - t324 * t132 + t258, t27 * t110 + t26 * t109 - t86 * t207 - g(1) * t260 - g(2) * t273 + t323 * t62 + t324 * t61 + (-t294 + t210) * t129, t6, t1, t28, t7, t29, 0, -t125 * t246 + t38 * t221 + t330 * t224 - t267 * t74 + t255, t125 * t252 - t39 * t221 - t331 * t224 - t261 * t267 + t257, t246 * t39 - t252 * t38 + t261 * t330 + t331 * t74 + t262, -g(1) * t274 - g(2) * t275 + t44 * t125 + t330 * t18 + t331 * t19 + t267 * t81 + t5 * t38 + t4 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t291, t305 * t220, t312, t291, t311, qJDD(3), t234 * t251 - t345, g(3) * t234 + t238 * t251, 0, 0, t322, -t124 + t348, (t130 + t292) * qJD(3) + t250, -t322, -t263, qJDD(3), -g(3) * t213 - t63 * qJD(3) - t129 * t132 + t350 * t212 + (qJDD(3) * t231 - t130 * t316) * pkin(3) + t26, g(3) * t212 + t64 * qJD(3) + t129 * t130 + t350 * t213 + (-qJDD(3) * t230 - t132 * t316) * pkin(3) - t27, (t62 + t63) * t132 + (-t61 + t64) * t130 + (-t230 * t79 - t231 * t80) * pkin(3), -t61 * t63 - t62 * t64 + (-t345 + t230 * t27 + t231 * t26 + (-t129 * t225 + t350) * t234) * pkin(3), t332, t359, t360, -t332, t361, t221, t149 * t221 + t325 * t224 + t74 * t98 + t353, -t150 * t221 - t326 * t224 + t261 * t98 + t354, -t149 * t252 + t150 * t246 + (t18 + t326) * t74 + (-t19 + t325) * t261, t4 * t150 + t5 * t149 - t81 * t98 - g(3) * t286 + t326 * t19 + t325 * t18 - t350 * (-pkin(4) * t212 - t338); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t132 * qJD(3) + t263, (-t130 + t292) * qJD(3) + t250, -t124 - t348, t62 * t130 + t61 * t132 - t349 + t86, 0, 0, 0, 0, 0, 0, -t246 - t328, t252 + t327, -t333 - t334, -t18 * t261 - t19 * t74 - t349 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t332, t359, t360, -t332, t361, t221, t19 * t224 + t353, t18 * t224 + t354, 0, 0;];
tau_reg = t2;
