% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRP11_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP11_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP11_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:48:40
% EndTime: 2019-03-09 12:48:53
% DurationCPUTime: 5.10s
% Computational Cost: add. (5616->511), mult. (11600->645), div. (0->0), fcn. (7493->10), ass. (0->272)
t211 = sin(qJ(2));
t309 = qJD(1) * t211;
t187 = pkin(2) * t309;
t215 = cos(qJ(2));
t332 = qJ(3) * t215;
t254 = pkin(8) * t211 - t332;
t100 = qJD(1) * t254 + t187;
t308 = qJD(1) * t215;
t183 = pkin(7) * t308;
t143 = pkin(3) * t308 + t183;
t214 = cos(qJ(4));
t122 = t214 * t143;
t210 = sin(qJ(4));
t325 = t210 * t211;
t247 = pkin(4) * t215 - pkin(9) * t325;
t301 = qJD(4) * t210;
t217 = -pkin(2) - pkin(8);
t346 = pkin(9) - t217;
t377 = -qJD(1) * t247 + t100 * t210 + t346 * t301 - t122;
t149 = t346 * t214;
t283 = t214 * t309;
t334 = t214 * t100 + t210 * t143;
t373 = pkin(9) * t283 + qJD(4) * t149 + t334;
t294 = t210 * qJD(2);
t133 = -t214 * t308 - t294;
t280 = t210 * t308;
t304 = qJD(2) * t214;
t134 = -t280 + t304;
t209 = sin(qJ(5));
t213 = cos(qJ(5));
t252 = t133 * t209 + t213 * t134;
t358 = t252 ^ 2;
t73 = -t213 * t133 + t134 * t209;
t70 = t73 ^ 2;
t376 = -t70 + t358;
t357 = pkin(3) + pkin(7);
t135 = t209 * t214 + t210 * t213;
t237 = t135 * t211;
t364 = qJD(4) + qJD(5);
t78 = t364 * t135;
t340 = qJD(1) * t237 + t78;
t297 = qJD(5) * t209;
t318 = t213 * t214;
t328 = t209 * t210;
t339 = -t209 * t301 - t210 * t297 + t213 * t283 - t309 * t328 + t318 * t364;
t375 = qJ(6) * t73;
t374 = t252 * t73;
t182 = pkin(7) * t309;
t372 = qJD(3) + t182;
t169 = qJD(4) + t309;
t157 = qJD(5) + t169;
t299 = qJD(4) * t215;
t282 = t210 * t299;
t232 = t211 * t304 + t282;
t290 = t215 * qJDD(1);
t300 = qJD(4) * t214;
t287 = qJD(2) * t300 + t210 * qJDD(2) + t214 * t290;
t221 = qJD(1) * t232 - t287;
t296 = qJD(5) * t213;
t292 = qJD(1) * qJD(2);
t278 = t211 * t292;
t67 = qJD(4) * t133 + t214 * qJDD(2) + (t278 - t290) * t210;
t20 = -t133 * t296 + t134 * t297 - t209 * t221 - t213 * t67;
t371 = t157 * t73 - t20;
t208 = qJ(4) + qJ(5);
t191 = sin(t208);
t201 = g(3) * t215;
t168 = pkin(2) * t278;
t302 = qJD(3) * t211;
t228 = qJD(2) * t254 - t302;
t193 = t211 * qJ(3);
t274 = -pkin(1) - t193;
t234 = t215 * t217 + t274;
t48 = qJD(1) * t228 + qJDD(1) * t234 + t168;
t277 = t215 * t292;
t291 = t211 * qJDD(1);
t236 = t277 + t291;
t167 = pkin(7) * t277;
t179 = pkin(7) * t291;
t276 = qJDD(3) + t167 + t179;
t69 = pkin(3) * t236 + qJDD(2) * t217 + t276;
t293 = pkin(3) * t309 + t372;
t94 = qJD(2) * t217 + t293;
t289 = t210 * t69 + t214 * t48 + t94 * t300;
t88 = t234 * qJD(1);
t11 = pkin(9) * t221 - t301 * t88 + t289;
t50 = -t210 * t88 + t214 * t94;
t40 = -pkin(9) * t134 + t50;
t34 = pkin(4) * t169 + t40;
t51 = t210 * t94 + t214 * t88;
t41 = pkin(9) * t133 + t51;
t132 = qJDD(4) + t236;
t268 = -t210 * t48 + t214 * t69;
t226 = -qJD(4) * t51 + t268;
t7 = pkin(4) * t132 - pkin(9) * t67 + t226;
t275 = -t213 * t11 - t209 * t7 - t34 * t296 + t41 * t297;
t205 = qJD(2) * qJ(3);
t114 = t205 + t143;
t80 = -pkin(4) * t133 + t114;
t192 = cos(t208);
t212 = sin(qJ(1));
t216 = cos(qJ(1));
t322 = t211 * t216;
t96 = t191 * t322 + t192 * t212;
t324 = t211 * t212;
t98 = -t191 * t324 + t192 * t216;
t370 = g(1) * t96 - g(2) * t98 - t191 * t201 + t73 * t80 + t275;
t36 = pkin(5) * t73 + qJD(6) + t80;
t369 = t252 * t36;
t368 = t377 * t213;
t367 = qJ(6) * t252;
t249 = -t318 + t328;
t101 = t249 * t215;
t148 = t346 * t210;
t312 = -t213 * t148 - t209 * t149;
t366 = pkin(4) * t300 + t372;
t365 = -t148 * t297 + t149 * t296 - t209 * t377 + t373 * t213;
t195 = t214 * pkin(4);
t284 = -pkin(3) - t195;
t305 = qJD(2) * t211;
t81 = (-pkin(7) + t284) * t305 - pkin(4) * t282;
t256 = g(1) * t216 + g(2) * t212;
t95 = -t191 * t212 + t192 * t322;
t97 = t191 * t216 + t192 * t324;
t363 = -g(1) * t95 - g(2) * t97 + t192 * t201;
t39 = t213 * t41;
t17 = t209 * t34 + t39;
t279 = -t209 * t11 + t213 * t7;
t227 = -qJD(5) * t17 + t279;
t362 = -t80 * t252 + t227 + t363;
t21 = qJD(5) * t252 + t209 * t67 - t213 * t221;
t361 = t157 * t252 - t21;
t360 = t20 * t249 - t252 * t340;
t124 = qJDD(5) + t132;
t1 = pkin(5) * t124 + qJ(6) * t20 - qJD(6) * t252 + t227;
t13 = t17 - t375;
t2 = -qJ(6) * t21 - qJD(6) * t73 - t275;
t285 = -g(1) * t322 - g(2) * t324 + t201;
t37 = t209 * t41;
t16 = t213 * t34 - t37;
t12 = t16 - t367;
t9 = pkin(5) * t157 + t12;
t359 = -t1 * t249 + t13 * t339 + t135 * t2 - t340 * t9 + t285;
t354 = t12 - t9;
t353 = pkin(4) * t210;
t352 = g(1) * t212;
t349 = g(2) * t216;
t348 = g(3) * t211;
t198 = t215 * pkin(2);
t345 = -qJ(6) * t339 - qJD(6) * t135 - t365;
t344 = -pkin(5) * t308 + qJ(6) * t340 - qJD(5) * t312 + qJD(6) * t249 + t209 * t373 + t368;
t343 = t213 * t40 - t37;
t153 = t357 * t211;
t138 = t214 * t153;
t311 = t198 + t193;
t150 = -pkin(1) - t311;
t125 = -pkin(8) * t215 + t150;
t273 = pkin(9) * t215 - t125;
t58 = pkin(4) * t211 + t210 * t273 + t138;
t137 = t210 * t153;
t313 = t214 * t125 + t137;
t317 = t214 * t215;
t64 = -pkin(9) * t317 + t313;
t341 = t209 * t58 + t213 * t64;
t337 = t67 * t210;
t336 = t67 * t214;
t255 = qJD(1) * t284;
t335 = -t211 * t255 + t366;
t333 = pkin(7) * qJDD(2);
t331 = qJD(4) * t88;
t330 = qJDD(2) * pkin(2);
t202 = -qJ(6) - pkin(9) - pkin(8);
t329 = t202 * t215;
t327 = t210 * t132;
t326 = t210 * t133;
t323 = t211 * t214;
t219 = qJD(1) ^ 2;
t321 = t211 * t219;
t320 = t212 * t214;
t319 = t212 * t215;
t109 = t214 * t132;
t316 = t214 * t216;
t315 = t215 * t216;
t314 = t217 * t132;
t146 = pkin(5) * t192 + t195;
t154 = t357 * t215;
t206 = t211 ^ 2;
t207 = t215 ^ 2;
t310 = t206 - t207;
t307 = qJD(2) * t133;
t306 = qJD(2) * t134;
t303 = qJD(2) * t215;
t298 = qJD(4) * t217;
t295 = t114 * qJD(4);
t288 = t215 * t321;
t180 = pkin(7) * t290;
t203 = qJDD(2) * qJ(3);
t204 = qJD(2) * qJD(3);
t286 = t180 + t203 + t204;
t112 = pkin(4) * t317 + t154;
t281 = t214 * t299;
t175 = qJ(3) + t353;
t144 = t357 * t303;
t186 = pkin(2) * t305;
t85 = t186 + t228;
t267 = t214 * t144 - t210 * t85;
t28 = t247 * qJD(2) + (t214 * t273 - t137) * qJD(4) + t267;
t238 = -t125 * t301 + t210 * t144 + t153 * t300 + t214 * t85;
t30 = pkin(9) * t232 + t238;
t272 = -t209 * t30 + t213 * t28;
t271 = -t209 * t40 - t39;
t270 = -t209 * t64 + t213 * t58;
t266 = -qJD(2) * pkin(2) + qJD(3);
t265 = t148 * t209 - t213 * t149;
t264 = -qJD(1) * t154 - t114;
t263 = t169 + t309;
t261 = t216 * pkin(1) + pkin(2) * t315 + t212 * pkin(7) + qJ(3) * t322;
t260 = -t179 - t285;
t259 = pkin(3) * t290 + t286;
t142 = t357 * t305;
t258 = -t124 * t249 - t157 * t340;
t218 = qJD(2) ^ 2;
t257 = pkin(7) * t218 + t349;
t251 = t134 * t214 + t326;
t147 = t182 + t266;
t152 = -t183 - t205;
t250 = t147 * t215 + t152 * t211;
t248 = t169 * t210;
t246 = t274 - t198;
t115 = t246 * qJD(1);
t245 = t115 * t309 + qJDD(3) - t260;
t244 = -0.2e1 * pkin(1) * t292 - t333;
t243 = t209 * t28 + t213 * t30 + t58 * t296 - t297 * t64;
t242 = -t169 * t300 - t327;
t240 = -t124 * t135 - t157 * t339;
t239 = -qJ(3) * t303 - t302;
t233 = 0.2e1 * qJDD(1) * pkin(1) - t257;
t230 = t333 + (-qJD(1) * t150 - t115) * qJD(2);
t229 = -t215 * t256 - t348;
t105 = t186 + t239;
t68 = qJD(1) * t239 + qJDD(1) * t246 + t168;
t223 = qJD(1) * t105 + qJDD(1) * t150 + t257 + t68;
t91 = pkin(7) * t278 - t286;
t99 = t276 - t330;
t222 = qJD(2) * t250 + t99 * t211 - t91 * t215;
t33 = pkin(4) * t287 + qJD(1) * t81 + t259;
t8 = t21 * pkin(5) + qJDD(6) + t33;
t199 = t216 * pkin(7);
t178 = pkin(4) * t213 + pkin(5);
t173 = g(1) * t319;
t166 = qJ(3) * t315;
t164 = qJ(3) * t319;
t145 = pkin(5) * t191 + t353;
t140 = pkin(3) + t146;
t139 = -qJ(3) * t308 + t187;
t119 = -t210 * t324 + t316;
t118 = t210 * t216 + t211 * t320;
t117 = t210 * t322 + t320;
t116 = -t210 * t212 + t211 * t316;
t102 = t135 * t215;
t71 = -qJD(1) * t142 + t259;
t57 = -qJ(6) * t135 + t312;
t56 = qJ(6) * t249 + t265;
t43 = t215 * t78 - t249 * t305;
t42 = qJD(2) * t237 + t101 * t364;
t25 = qJ(6) * t101 + t341;
t24 = pkin(5) * t211 + qJ(6) * t102 + t270;
t15 = t343 - t367;
t14 = t271 + t375;
t4 = qJ(6) * t43 + qJD(6) * t101 + t243;
t3 = pkin(5) * t303 - qJ(6) * t42 - qJD(5) * t341 + qJD(6) * t102 + t272;
t5 = [qJDD(1), -t349 + t352, t256, qJDD(1) * t206 + 0.2e1 * t211 * t277, 0.2e1 * t211 * t290 - 0.2e1 * t292 * t310, qJDD(2) * t211 + t215 * t218, qJDD(2) * t215 - t211 * t218, 0, t211 * t244 + t215 * t233 + t173, t244 * t215 + (-t233 - t352) * t211 (t206 + t207) * qJDD(1) * pkin(7) + t222 - t256, t211 * t230 + t215 * t223 - t173, t230 * t215 + (-t223 + t352) * t211, pkin(7) * t222 - g(1) * t199 - g(2) * t261 + t115 * t105 + t68 * t150 - t246 * t352, -t215 * t337 + (t211 * t294 - t281) * t134, t251 * t305 + (-t210 * (t214 * t278 - t287) - t336 + (-t214 * t133 + (t134 - t280) * t210) * qJD(4)) * t215 (t169 * t294 + t67) * t211 + (t242 + t306) * t215 (t263 * t304 - t287) * t211 + (t263 * t301 - t109 + t307) * t215, t132 * t211 + t169 * t303, t267 * t169 + (-t125 * t210 + t138) * t132 + t268 * t211 + t142 * t133 + t154 * t287 + t71 * t317 - g(1) * t119 - g(2) * t117 + (t50 * t215 + t264 * t323) * qJD(2) + ((-t125 * t169 - t88 * t211) * t214 + (-t153 * t169 - t94 * t211 + t215 * t264) * t210) * qJD(4), -t238 * t169 - t313 * t132 - t142 * t134 + t154 * t67 + g(1) * t118 - g(2) * t116 + ((qJD(2) * t114 + t331) * t210 - t289) * t211 + (-t51 * qJD(2) - t71 * t210 - t214 * t295) * t215, t102 * t20 + t252 * t42, -t101 * t20 + t102 * t21 + t252 * t43 - t42 * t73, -t102 * t124 + t157 * t42 - t20 * t211 + t252 * t303, t101 * t124 + t157 * t43 - t21 * t211 - t303 * t73, t124 * t211 + t157 * t303, t272 * t157 + t270 * t124 + t279 * t211 + t16 * t303 + t81 * t73 + t112 * t21 - t33 * t101 - t80 * t43 - g(1) * t98 - g(2) * t96 + (-t157 * t341 - t17 * t211) * qJD(5), g(1) * t97 - g(2) * t95 - t33 * t102 - t112 * t20 - t124 * t341 - t157 * t243 - t17 * t303 + t211 * t275 + t252 * t81 + t80 * t42, -g(2) * t315 + t1 * t102 + t101 * t2 + t13 * t43 + t20 * t24 - t21 * t25 - t252 * t3 - t4 * t73 - t42 * t9 + t173, t2 * t25 + t13 * t4 + t1 * t24 + t9 * t3 + t8 * (-pkin(5) * t101 + t112) - g(1) * (t140 * t216 + t199) - g(2) * (t145 * t322 - t202 * t315 + t261) + (-g(1) * (-t145 * t211 + t246 + t329) - g(2) * t140) * t212 + (-pkin(5) * t43 + t81) * t36; 0, 0, 0, -t288, t310 * t219, t291, t290, qJDD(2), pkin(1) * t321 + t260, t348 - t180 + (pkin(1) * t219 + t256) * t215 (-pkin(2) * t211 + t332) * qJDD(1) + ((-t152 - t205) * t211 + (-t147 + t266) * t215) * qJD(1), -t139 * t308 + t245 - 0.2e1 * t330, t180 + 0.2e1 * t203 + 0.2e1 * t204 + (qJD(1) * t139 - g(3)) * t211 + (qJD(1) * t115 - t256) * t215, -t91 * qJ(3) - t152 * qJD(3) - t99 * pkin(2) - t115 * t139 - g(1) * (-pkin(2) * t322 + t166) - g(2) * (-pkin(2) * t324 + t164) - g(3) * t311 - t250 * qJD(1) * pkin(7), -t134 * t248 + t336, -t214 * t287 - t337 - t251 * qJD(4) + (t210 * t281 + (-t326 + (-t134 + t304) * t214) * t211) * qJD(1), -t169 * t301 + t109 + (-t134 * t215 - t169 * t325) * qJD(1) (-t133 * t215 - t169 * t323) * qJD(1) + t242, -t169 * t308, qJ(3) * t287 - t122 * t169 - t50 * t308 - t293 * t133 + (t295 + t314 + (t114 - t205) * t309) * t214 + (-t348 + t71 + (t100 - t298) * t169 + (-qJ(3) * qJD(1) * qJD(4) - t256) * t215) * t210, qJ(3) * t67 + t334 * t169 + t51 * t308 + t293 * t134 + (-t114 * t169 - t314) * t210 + (-t169 * t298 + t229 + t71) * t214, t360, t135 * t20 + t21 * t249 - t252 * t339 + t340 * t73, -t252 * t308 + t258, t308 * t73 + t240, -t157 * t308, t265 * t124 + t175 * t21 + t33 * t135 - t16 * t308 + t339 * t80 + t335 * t73 + (t148 * t296 + (qJD(5) * t149 + t373) * t209 + t368) * t157 + t229 * t191, -t312 * t124 + t157 * t365 + t17 * t308 - t175 * t20 + t229 * t192 - t249 * t33 + t335 * t252 - t340 * t80, t20 * t56 - t21 * t57 - t252 * t344 - t345 * t73 - t359, t2 * t57 + t1 * t56 + t8 * (pkin(5) * t135 + t175) - g(1) * (t145 * t315 + t166) - g(2) * (t145 * t319 + t164) - g(3) * (t311 - t329) + t344 * t9 + (pkin(5) * t339 + t366) * t36 + t345 * t13 + (-g(3) * t145 - t255 * t36 + t256 * (pkin(2) - t202)) * t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t291, qJDD(2) + t288, -t206 * t219 - t218, qJD(2) * t152 + t167 + t245 - t330, 0, 0, 0, 0, 0, -t169 * t248 + t109 + t307, -t169 ^ 2 * t214 - t306 - t327, 0, 0, 0, 0, 0, -qJD(2) * t73 + t258, -qJD(2) * t252 + t240, -t135 * t21 - t339 * t73 - t360, -qJD(2) * t36 + t359; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134 * t133, -t133 ^ 2 + t134 ^ 2, -t133 * t169 + t67, t134 * t169 + t221, t132, -g(1) * t116 - g(2) * t118 + g(3) * t317 - t114 * t134 + t169 * t51 + t226, g(1) * t117 - g(2) * t119 - t114 * t133 + t169 * t50 + (t331 - t201) * t210 - t289, t374, t376, t371, t361, t124, -t271 * t157 + (t124 * t213 - t134 * t73 - t157 * t297) * pkin(4) + t362, t343 * t157 + (-t124 * t209 - t134 * t252 - t157 * t296) * pkin(4) + t370, t13 * t252 + t14 * t252 + t15 * t73 + t178 * t20 - t73 * t9 + (-t209 * t21 + (t209 * t252 - t213 * t73) * qJD(5)) * pkin(4), t1 * t178 - t13 * t15 - t9 * t14 - pkin(5) * t369 - g(1) * (-t145 * t212 + t146 * t322) - g(2) * (t145 * t216 + t146 * t324) + t146 * t201 + (-t36 * t134 + t2 * t209 + (t13 * t213 - t209 * t9) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t374, t376, t371, t361, t124, t157 * t17 + t362, t157 * t16 + t370, pkin(5) * t20 + t354 * t73, -t354 * t13 + (t1 + t363 - t369) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70 - t358, t13 * t73 + t252 * t9 + t229 + t8;];
tau_reg  = t5;
