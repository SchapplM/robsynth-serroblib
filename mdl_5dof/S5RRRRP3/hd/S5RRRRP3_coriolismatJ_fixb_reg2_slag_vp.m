% Calculate inertial parameters regressor of coriolis matrix for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRRP3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:30
% EndTime: 2019-12-31 21:49:37
% DurationCPUTime: 4.20s
% Computational Cost: add. (2898->389), mult. (6272->420), div. (0->0), fcn. (4626->6), ass. (0->300)
t210 = cos(qJ(3));
t356 = cos(qJ(2));
t274 = t356 * t210;
t207 = sin(qJ(3));
t208 = sin(qJ(2));
t320 = t207 * t208;
t164 = (t274 - t320) * pkin(1);
t209 = cos(qJ(4));
t205 = t209 ^ 2;
t323 = t205 * t164;
t206 = sin(qJ(4));
t204 = t206 ^ 2;
t324 = t204 * t164;
t369 = t323 / 0.2e1 + t324 / 0.2e1;
t288 = t356 * pkin(1);
t251 = t288 + pkin(2);
t176 = t210 * t251;
t152 = pkin(1) * t320 - t176;
t303 = t204 + t205;
t70 = t303 * t152;
t192 = t205 - t204;
t290 = -qJD(2) - qJD(3);
t260 = qJD(1) - t290;
t368 = t260 * t192;
t367 = t369 * pkin(8);
t348 = pkin(2) * qJD(2);
t285 = t207 * t348;
t101 = t207 * pkin(2);
t233 = t207 * t251;
t319 = t208 * t210;
t153 = pkin(1) * t319 + t233;
t363 = t153 / 0.2e1;
t250 = t101 / 0.2e1 + t363;
t275 = t356 * t207;
t162 = (t275 + t319) * pkin(1);
t361 = -t162 / 0.2e1;
t229 = t361 + t250;
t59 = t229 * t206;
t339 = t59 * qJD(1) + t206 * t285;
t60 = t229 * t209;
t338 = t60 * qJD(1) + t209 * t285;
t366 = pkin(3) / 0.2e1;
t242 = -t209 * pkin(4) - t206 * qJ(5);
t170 = -pkin(3) + t242;
t92 = t152 + t170;
t365 = t92 / 0.2e1;
t148 = -pkin(3) + t152;
t364 = -t148 / 0.2e1;
t350 = t210 * pkin(2);
t161 = t170 - t350;
t362 = t161 / 0.2e1;
t360 = t170 / 0.2e1;
t196 = -pkin(3) - t350;
t359 = -t196 / 0.2e1;
t358 = t206 / 0.2e1;
t357 = t207 / 0.2e1;
t355 = pkin(1) * t208;
t354 = pkin(3) * t209;
t352 = t162 * pkin(3);
t351 = t206 * pkin(4);
t195 = pkin(8) + t101;
t349 = (t323 + t324) * t195;
t347 = pkin(2) * qJD(3);
t318 = t209 * qJ(5);
t171 = -t318 + t351;
t346 = t92 * t171;
t345 = t92 * t209;
t295 = t153 * qJD(1);
t124 = t209 * t295;
t344 = t60 * qJD(2) + t124;
t147 = t161 * t206;
t131 = t147 / 0.2e1;
t84 = t92 * t206;
t79 = t84 / 0.2e1;
t343 = t131 + t79;
t132 = -t147 / 0.2e1;
t80 = -t84 / 0.2e1;
t342 = t132 + t80;
t294 = t162 * qJD(1);
t144 = t209 * t294;
t341 = -t60 * qJD(3) + t144;
t149 = pkin(8) + t153;
t249 = t149 * t70;
t10 = t92 * t153 - t249;
t337 = t10 * qJD(1);
t73 = t303 * t164;
t248 = t149 * t73;
t12 = t92 * t162 + t248;
t336 = t12 * qJD(1);
t13 = t148 * t153 - t249;
t335 = t13 * qJD(1);
t14 = t148 * t162 + t248;
t334 = t14 * qJD(1);
t333 = t148 * t209;
t332 = t161 * t209;
t331 = t162 * t170;
t330 = t170 * t209;
t329 = t171 * t161;
t328 = t171 * t170;
t327 = t171 * t206;
t326 = t171 * t209;
t325 = t196 * t209;
t322 = t206 * t152;
t321 = t206 * t164;
t317 = t209 * t152;
t316 = t209 * t164;
t37 = t152 * t162 + t153 * t164;
t315 = t37 * qJD(1);
t65 = t327 + t345;
t312 = t65 * qJD(1);
t66 = -t84 + t326;
t311 = t66 * qJD(1);
t310 = t70 * qJD(1);
t309 = t73 * qJD(1);
t146 = t153 * qJD(3);
t122 = t206 * t146;
t154 = t162 * qJD(2);
t142 = t206 * t154;
t308 = t122 + t142;
t307 = pkin(8) * t70;
t306 = -t154 - t146;
t305 = t303 * pkin(8) * t350;
t287 = t207 * t347;
t189 = t206 * t287;
t201 = t204 * qJD(5);
t304 = t201 - t189;
t302 = qJD(1) * t206;
t301 = qJD(2) * t206;
t300 = qJD(3) * t206;
t299 = qJD(4) * qJ(5);
t298 = t101 * qJD(1);
t103 = t176 / 0.2e1 + (-t288 / 0.2e1 + pkin(2) / 0.2e1) * t210;
t297 = t103 * qJD(1);
t296 = t152 * qJD(1);
t293 = t164 * qJD(1);
t292 = t192 * qJD(4);
t202 = t206 * qJD(4);
t203 = t209 * qJD(4);
t291 = t209 * qJD(5);
t286 = pkin(8) * t202;
t284 = pkin(8) * t203;
t283 = -t350 / 0.2e1;
t282 = t350 / 0.2e1;
t281 = t366 + t364;
t280 = t366 + t359;
t279 = qJD(1) * t346;
t278 = t92 * t302;
t277 = t362 + t365;
t276 = t360 + t365;
t273 = t148 * t302;
t272 = qJD(1) * t333;
t271 = t149 * t202;
t270 = t195 * t202;
t269 = t149 * t203;
t268 = t195 * t203;
t265 = t360 + t362;
t264 = t359 + t364;
t263 = t205 / 0.2e1 + t204 / 0.2e1;
t262 = t356 * qJD(1);
t261 = t356 * qJD(2);
t259 = pkin(2) * t290;
t258 = t303 * t210;
t123 = t206 * t295;
t257 = -t59 * qJD(2) - t123;
t143 = t206 * t294;
t256 = t59 * qJD(3) - t143;
t165 = t170 * t206;
t157 = t165 / 0.2e1;
t255 = t157 - t326;
t158 = -t165 / 0.2e1;
t254 = t158 + t326;
t253 = -t330 / 0.2e1 - t327;
t252 = t209 * t287;
t247 = t195 * t258;
t58 = t162 * t358 + t250 * t206;
t246 = t58 * qJD(2) + t122 + t123;
t245 = t58 * qJD(3) + t142 + t143;
t244 = t263 * t152;
t243 = t263 * t210;
t231 = t149 * t243;
t232 = t195 * t244;
t212 = (t92 * t357 + t231) * pkin(2) - t232 + t153 * t362;
t237 = t263 * t164 * pkin(8);
t2 = -t331 / 0.2e1 - t237 + t212;
t63 = (t161 * t207 + t247) * pkin(2);
t241 = t2 * qJD(1) + t63 * qJD(2);
t211 = (t148 * t357 + t231) * pkin(2) - t232 + t196 * t363;
t4 = t352 / 0.2e1 - t237 + t211;
t74 = (t196 * t207 + t247) * pkin(2);
t240 = -t4 * qJD(1) - t74 * qJD(2);
t136 = -t321 / 0.2e1;
t22 = t136 + t326 + t342;
t78 = -t147 + t326;
t239 = t22 * qJD(1) + t78 * qJD(2);
t138 = t316 / 0.2e1;
t23 = t209 * t277 + t138 + t327;
t77 = t327 + t332;
t238 = t23 * qJD(1) + t77 * qJD(2);
t163 = pkin(2) * t258;
t20 = t303 * (t282 - t152 / 0.2e1 - t164 / 0.2e1);
t11 = -t20 * qJD(1) - t163 * qJD(2);
t236 = qJD(4) * t171 - qJD(5) * t206;
t61 = (-t250 + t361) * t209;
t235 = t61 * qJD(3) - t154 * t209 - t144;
t234 = t61 * qJD(2) - t146 * t209 - t124;
t230 = t318 / 0.2e1 - t351 / 0.2e1;
t220 = t230 * t164;
t17 = -t171 * t277 + t220;
t227 = -t17 * qJD(1) + qJD(2) * t329;
t135 = t321 / 0.2e1;
t31 = t135 + t343;
t226 = t31 * qJD(1) + t161 * t301;
t42 = t206 * t264 + t136;
t225 = t42 * qJD(1) - t196 * t301;
t139 = -t316 / 0.2e1;
t43 = t209 * t264 + t139;
t224 = t43 * qJD(1) - qJD(2) * t325;
t223 = -t252 - t338;
t222 = t306 * t209;
t221 = t230 * t152;
t116 = t322 / 0.2e1;
t26 = t116 + t254 + t80;
t179 = t206 * t283;
t38 = t132 + t179 + t254;
t88 = -t165 + t326;
t219 = t26 * qJD(1) + t38 * qJD(2) + t88 * qJD(3);
t118 = -t317 / 0.2e1;
t27 = t209 * t276 + t118 + t327;
t184 = t209 * t282;
t39 = t209 * t265 + t184 + t327;
t87 = t327 + t330;
t218 = t27 * qJD(1) + t39 * qJD(2) + t87 * qJD(3);
t217 = t230 * t350;
t100 = t206 * t280 + t179;
t54 = t206 * t281 + t116;
t216 = pkin(3) * t300 + t54 * qJD(1) + t100 * qJD(2);
t185 = t209 * t283;
t102 = t209 * t280 + t185;
t119 = t317 / 0.2e1;
t55 = t209 * t281 + t119;
t215 = t55 * qJD(1) + t102 * qJD(2) + qJD(3) * t354;
t15 = -t171 * t276 - t221;
t35 = -t171 * t265 + t217;
t214 = -t15 * qJD(1) - t35 * qJD(2) + qJD(3) * t328;
t115 = -t322 / 0.2e1;
t33 = t115 + t157 + t79;
t178 = t206 * t282;
t68 = t178 + t157 + t131;
t213 = t33 * qJD(1) + t68 * qJD(2) + t170 * t300;
t151 = qJD(4) * t242 + t291;
t200 = -t354 / 0.2e1;
t199 = -pkin(3) * t206 / 0.2e1;
t194 = t206 * t203;
t193 = t206 * t291;
t169 = t325 / 0.2e1;
t168 = t196 * t358;
t160 = t260 * t204;
t156 = t164 * qJD(2);
t155 = t163 * qJD(3);
t145 = t152 * qJD(3);
t133 = -t332 / 0.2e1;
t113 = t328 / 0.2e1;
t112 = t333 / 0.2e1;
t111 = t148 * t358;
t106 = t260 * t209 * t206;
t105 = t200 + t169 + t185;
t104 = t199 + t168 + t179;
t93 = t329 / 0.2e1;
t83 = t283 - t176 / 0.2e1 + (t320 - t274 / 0.2e1) * pkin(1);
t82 = -t101 / 0.2e1 - t233 / 0.2e1 + (-t319 - t275 / 0.2e1) * pkin(1);
t81 = -t345 / 0.2e1;
t72 = t346 / 0.2e1;
t71 = t73 * qJD(2);
t69 = t158 + t132 + t178;
t67 = t70 * qJD(3);
t57 = t200 + t112 + t119;
t56 = t199 + t111 + t116;
t45 = t169 + t112 + t139;
t44 = t168 + t111 + t136;
t41 = t131 + t179 + t255;
t40 = t133 + t184 + t253;
t36 = t113 + t93 + t217;
t34 = t158 + t80 + t115;
t32 = t135 + t342;
t29 = t116 + t255 + t79;
t28 = t118 + t253 + t81;
t25 = t136 - t326 + t343;
t24 = t133 + t138 + t81 - t327;
t21 = t71 - t67;
t19 = pkin(2) * t243 - t244 + t369;
t18 = t93 + t72 + t220;
t16 = t113 + t72 - t221;
t9 = -t11 + t155;
t8 = t20 * qJD(3) - t309;
t7 = -t20 * qJD(2) + t310;
t6 = t19 * qJD(3) + t309 + t71;
t5 = t19 * qJD(2) - t310 - t67;
t3 = -t352 / 0.2e1 + t211 + t367;
t1 = t331 / 0.2e1 + t212 + t367;
t30 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t355, -pkin(1) * t261, 0, 0, 0, 0, 0, 0, 0, 0, t306, -t156 + t145, 0, t37 * qJD(2), t194, t292, 0, -t194, 0, 0, t148 * t202 + t222, t148 * t203 + t308, t21, t14 * qJD(2) + t13 * qJD(3), t194, 0, -t292, 0, 0, -t194, -t66 * qJD(4) + t193 + t222, t21, -t65 * qJD(4) + t201 - t308, t12 * qJD(2) + t10 * qJD(3) + t236 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-qJD(1) - qJD(2)) * t355, (-t262 - t261) * pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, t82 * qJD(3) - t154 - t294, t83 * qJD(3) - t156 - t293, 0, t315 + (-t162 * t210 + t164 * t207) * t348, t194, t292, 0, -t194, 0, 0, t44 * qJD(4) + t235, t45 * qJD(4) + t245, t6, t334 + (t162 * t196 + t349) * qJD(2) + t3 * qJD(3), t194, 0, -t292, 0, 0, -t194, t25 * qJD(4) + t193 + t235, t6, t24 * qJD(4) + t201 - t245, t336 + (t162 * t161 + t349) * qJD(2) + t1 * qJD(3) + t18 * qJD(4) + t32 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82 * qJD(2) - t146 - t295, t83 * qJD(2) + t145 + t296, 0, 0, t194, t292, 0, -t194, 0, 0, t56 * qJD(4) + t234, t57 * qJD(4) + t246, t5, t335 + t3 * qJD(2) + (-t153 * pkin(3) - t307) * qJD(3), t194, 0, -t292, 0, 0, -t194, t29 * qJD(4) + t193 + t234, t5, t28 * qJD(4) + t201 - t246, t337 + t1 * qJD(2) + (t153 * t170 - t307) * qJD(3) + t16 * qJD(4) + t34 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t368, t203, -t106, -t202, 0, t44 * qJD(2) + t56 * qJD(3) - t269 + t273, t45 * qJD(2) + t57 * qJD(3) + t271 + t272, 0, 0, t106, t203, -t368, 0, t202, -t106, t25 * qJD(2) + t29 * qJD(3) - t269 - t311, t151, t24 * qJD(2) + t28 * qJD(3) - t271 - t312, t18 * qJD(2) + t16 * qJD(3) + t149 * t151 + t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t203, t160, t32 * qJD(2) + t34 * qJD(3) + t269 - t278; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1) * t355, pkin(1) * t262, 0, 0, 0, 0, 0, 0, 0, 0, -t101 * qJD(3) + t294, -t103 * qJD(3) + t293, 0, -t315, t194, t292, 0, -t194, 0, 0, -t42 * qJD(4) + t341, -t43 * qJD(4) + t256, t8, t4 * qJD(3) - t334, t194, 0, -t292, 0, 0, -t194, -t22 * qJD(4) + t193 + t341, t8, -t23 * qJD(4) + t201 - t256, t2 * qJD(3) - t17 * qJD(4) - t31 * qJD(5) - t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t287, -t210 * t347, 0, 0, t194, t292, 0, -t194, 0, 0, t196 * t202 - t252, t196 * t203 + t189, t155, t74 * qJD(3), t194, 0, -t292, 0, 0, -t194, -t78 * qJD(4) + t193 - t252, t155, -t77 * qJD(4) + t304, t63 * qJD(3) + t161 * t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t207 * t259 - t298, t210 * t259 - t297, 0, 0, t194, t292, 0, -t194, 0, 0, t104 * qJD(4) + t223, t105 * qJD(4) + t189 + t339, t9, (-pkin(3) * t101 + t305) * qJD(3) - t240, t194, 0, -t292, 0, 0, -t194, t41 * qJD(4) + t193 + t223, t9, t40 * qJD(4) + t304 - t339, (t170 * t101 + t305) * qJD(3) + t36 * qJD(4) + t69 * qJD(5) + t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t368, t203, -t106, -t202, 0, t104 * qJD(3) - t225 - t268, t105 * qJD(3) - t224 + t270, 0, 0, t106, t203, -t368, 0, t202, -t106, t41 * qJD(3) - t239 - t268, t151, t40 * qJD(3) - t238 - t270, t36 * qJD(3) + t151 * t195 + t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t203, t160, t69 * qJD(3) - t226 + t268; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101 * qJD(2) + t295, t103 * qJD(2) - t296, 0, 0, t194, t292, 0, -t194, 0, 0, -t54 * qJD(4) + t344, -t55 * qJD(4) + t257, t7, -t4 * qJD(2) - t335, t194, 0, -t292, 0, 0, -t194, -t26 * qJD(4) + t193 + t344, t7, -t27 * qJD(4) + t201 - t257, -t2 * qJD(2) - t15 * qJD(4) - t33 * qJD(5) - t337; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285 + t298, t210 * t348 + t297, 0, 0, t194, t292, 0, -t194, 0, 0, -t100 * qJD(4) + t338, -t102 * qJD(4) - t339, t11, t240, t194, 0, -t292, 0, 0, -t194, -t38 * qJD(4) + t193 + t338, t11, -t39 * qJD(4) + t201 + t339, -t35 * qJD(4) - t68 * qJD(5) - t241; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t194, t292, 0, -t194, 0, 0, -pkin(3) * t202, -pkin(3) * t203, 0, 0, t194, 0, -t292, 0, 0, -t194, -t88 * qJD(4) + t193, 0, -t87 * qJD(4) + t201, t236 * t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t368, t203, -t106, -t202, 0, -t216 - t284, -t215 + t286, 0, 0, t106, t203, -t368, 0, t202, -t106, -t219 - t284, t151, -t218 - t286, pkin(8) * t151 + t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t203, t160, -t213 + t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t368, 0, t106, 0, 0, t42 * qJD(2) + t54 * qJD(3) - t273, t43 * qJD(2) + t55 * qJD(3) - t272, 0, 0, -t106, 0, t368, 0, 0, t106, t22 * qJD(2) + t26 * qJD(3) + t311, 0, t23 * qJD(2) + t27 * qJD(3) + t312, t17 * qJD(2) + t15 * qJD(3) - t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t368, 0, t106, 0, 0, t100 * qJD(3) + t225, t102 * qJD(3) + t224, 0, 0, -t106, 0, t368, 0, 0, t106, t38 * qJD(3) + t239, 0, t39 * qJD(3) + t238, t35 * qJD(3) - t227; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, -t368, 0, t106, 0, 0, t216, t215, 0, 0, -t106, 0, t368, 0, 0, t106, t219, 0, t218, -t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, 0, -t160, t31 * qJD(2) + t33 * qJD(3) + t278; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, 0, -t160, t68 * qJD(3) + t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106, 0, -t160, t213; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t30;
