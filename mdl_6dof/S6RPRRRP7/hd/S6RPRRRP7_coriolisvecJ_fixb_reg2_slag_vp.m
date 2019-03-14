% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRP7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:21:00
% EndTime: 2019-03-09 06:21:18
% DurationCPUTime: 6.78s
% Computational Cost: add. (13402->494), mult. (34197->613), div. (0->0), fcn. (26150->8), ass. (0->225)
t212 = sin(pkin(10));
t213 = cos(pkin(10));
t216 = sin(qJ(3));
t325 = cos(qJ(3));
t185 = t212 * t325 + t216 * t213;
t178 = t185 * qJD(1);
t215 = sin(qJ(4));
t217 = cos(qJ(4));
t152 = qJD(3) * t217 - t215 * t178;
t153 = qJD(3) * t215 + t178 * t217;
t214 = sin(qJ(5));
t324 = cos(qJ(5));
t92 = -t324 * t152 + t153 * t214;
t308 = t178 * t92;
t234 = -t216 * t212 + t325 * t213;
t352 = t234 * qJD(1);
t172 = qJD(4) - t352;
t169 = qJD(5) + t172;
t181 = t185 * qJD(3);
t170 = qJD(1) * t181;
t264 = t324 * t217;
t289 = t214 * t215;
t188 = -t264 + t289;
t265 = t324 * t215;
t288 = t214 * t217;
t189 = t265 + t288;
t332 = qJD(4) + qJD(5);
t143 = t332 * t189;
t344 = -t189 * t352 + t143;
t329 = t344 * t169 + t170 * t188;
t360 = -t308 - t329;
t359 = -t329 + t308;
t261 = t324 * qJD(4);
t250 = t217 * t261;
t260 = t324 * qJD(5);
t280 = -t188 * t352 - t217 * t260 + t289 * t332 - t250;
t278 = qJD(4) * t215;
t222 = qJD(3) * t352;
t342 = qJD(3) * qJD(4) + t222;
t105 = t178 * t278 - t342 * t217;
t233 = t214 * t152 + t153 * t324;
t277 = qJD(4) * t217;
t269 = t178 * t277 + t342 * t215;
t48 = qJD(5) * t233 - t214 * t105 + t324 * t269;
t252 = -t189 * t48 + t280 * t92;
t276 = qJD(5) * t214;
t47 = t324 * t105 - t152 * t260 + t153 * t276 + t214 * t269;
t354 = -t47 * t188 + t344 * t233;
t358 = t252 + t354;
t357 = t354 - t252;
t204 = -pkin(2) * t213 - pkin(1);
t192 = qJD(1) * t204 + qJD(2);
t356 = qJD(2) + t192;
t245 = -t169 * t280 + t189 * t170;
t307 = t178 * t233;
t355 = t245 + t307;
t347 = -t169 * t233 + t48;
t294 = t352 * t215;
t353 = t278 - t294;
t320 = pkin(7) + qJ(2);
t193 = t320 * t212;
t186 = qJD(1) * t193;
t194 = t320 * t213;
t187 = qJD(1) * t194;
t336 = -t325 * t186 - t216 * t187;
t130 = -qJD(3) * pkin(3) - t336;
t87 = -t152 * pkin(4) + t130;
t44 = t92 * pkin(5) - qJ(6) * t233 + t87;
t351 = t44 * t92;
t350 = t87 * t92;
t348 = t92 * t233;
t295 = t172 * t215;
t346 = t153 * t295;
t180 = t234 * qJD(3);
t286 = t215 * t180;
t231 = t185 * t277 + t286;
t328 = t233 ^ 2;
t272 = t92 ^ 2 - t328;
t341 = t169 * t92 - t47;
t62 = pkin(5) * t233 + qJ(6) * t92;
t326 = -pkin(9) - pkin(8);
t195 = t326 * t215;
t196 = t326 * t217;
t151 = t214 * t195 - t196 * t324;
t268 = qJD(4) * t326;
t191 = t215 * t268;
t132 = pkin(3) * t178 - pkin(8) * t352;
t82 = t217 * t132 - t215 * t336;
t67 = -pkin(9) * t217 * t352 + pkin(4) * t178 + t82;
t83 = t215 * t132 + t217 * t336;
t73 = -pkin(9) * t294 + t83;
t317 = qJD(5) * t151 - t250 * t326 + t324 * t67 + (t191 - t73) * t214;
t110 = -pkin(3) * t352 - pkin(8) * t178 + t192;
t118 = t170 * pkin(3) - pkin(8) * t222;
t137 = -t216 * t186 + t187 * t325;
t131 = qJD(3) * pkin(8) + t137;
t223 = t234 * qJD(2);
t97 = qJD(1) * t223 + qJD(3) * t336;
t41 = t110 * t277 + t215 * t118 - t131 * t278 + t217 * t97;
t75 = t217 * t110 - t131 * t215;
t241 = -t172 * t75 + t41;
t76 = t110 * t215 + t131 * t217;
t42 = -qJD(4) * t76 + t217 * t118 - t215 * t97;
t338 = t172 * t76 + t42;
t249 = t353 * pkin(4) - t137;
t134 = -pkin(3) * t234 - pkin(8) * t185 + t204;
t149 = -t216 * t193 + t194 * t325;
t141 = t217 * t149;
t85 = t215 * t134 + t141;
t335 = -t325 * t193 - t194 * t216;
t334 = t178 * qJD(3);
t164 = t170 * pkin(5);
t18 = pkin(4) * t170 + pkin(9) * t105 + t42;
t23 = -pkin(9) * t269 + t41;
t65 = -pkin(9) * t153 + t75;
t55 = pkin(4) * t172 + t65;
t66 = pkin(9) * t152 + t76;
t4 = t324 * t18 - t214 * t23 - t66 * t260 - t55 * t276;
t2 = -t164 - t4;
t227 = t233 * t44 + t2;
t331 = -t87 * t233 + t4;
t291 = t185 * t217;
t84 = t217 * t134 - t149 * t215;
t71 = -pkin(4) * t234 - pkin(9) * t291 + t84;
t292 = t185 * t215;
t77 = -pkin(9) * t292 + t85;
t314 = t214 * t71 + t324 * t77;
t111 = qJD(3) * t335 + t223;
t133 = pkin(3) * t181 - pkin(8) * t180;
t256 = -t111 * t215 + t217 * t133;
t284 = t217 * t180;
t32 = -pkin(9) * t284 + pkin(4) * t181 + (-t141 + (pkin(9) * t185 - t134) * t215) * qJD(4) + t256;
t49 = t217 * t111 + t215 * t133 + t134 * t277 - t149 * t278;
t38 = -pkin(9) * t231 + t49;
t8 = -qJD(5) * t314 - t214 * t38 + t32 * t324;
t327 = t178 ^ 2;
t34 = t214 * t67 + t324 * t73;
t28 = qJ(6) * t178 + t34;
t232 = t195 * t324 + t214 * t196;
t88 = qJD(5) * t232 + t191 * t324 + t268 * t288;
t319 = -t28 + t88;
t318 = -t178 * pkin(5) - t317;
t316 = t34 - t88;
t315 = -t344 * pkin(5) - t280 * qJ(6) + qJD(6) * t189 - t249;
t225 = t185 * qJD(2);
t98 = qJD(1) * t225 + t137 * qJD(3);
t313 = t335 * t98;
t271 = t324 * t66;
t20 = t214 * t55 + t271;
t312 = t169 * t20;
t306 = t214 * t66;
t27 = t324 * t65 - t306;
t304 = t260 * pkin(4) + qJD(6) - t27;
t303 = t105 * t215;
t302 = t232 * t170;
t301 = t151 * t170;
t300 = t152 * t352;
t299 = t152 * t178;
t298 = t153 * t152;
t297 = t153 * t178;
t296 = t169 * t178;
t138 = t170 * t234;
t293 = t178 * t352;
t287 = t215 * t170;
t156 = t217 * t170;
t19 = t324 * t55 - t306;
t283 = qJD(6) - t19;
t101 = t215 * t269;
t275 = t152 * qJD(4);
t282 = t217 * t275 - t101;
t279 = t212 ^ 2 + t213 ^ 2;
t270 = t185 * t289;
t209 = -pkin(4) * t217 - pkin(3);
t263 = t185 * t278;
t259 = -t214 * t18 - t324 * t23 - t55 * t260 + t66 * t276;
t258 = t279 * qJD(1) ^ 2;
t254 = t172 * t217;
t251 = -t151 * t48 + t232 * t47 - t88 * t92;
t26 = t214 * t65 + t271;
t248 = t276 * pkin(4) - t26;
t247 = t269 * t217;
t113 = pkin(4) * t292 - t335;
t124 = t189 * t185;
t57 = t180 * t265 - t214 * t263 - qJD(5) * t270 + (t180 * t214 + (t261 + t260) * t185) * t217;
t243 = t124 * t48 + t57 * t92;
t242 = t76 * t215 + t75 * t217;
t240 = 0.2e1 * t279 * qJD(2) * qJD(1);
t159 = t169 * qJD(6);
t161 = t170 * qJ(6);
t1 = t161 + t159 - t259;
t239 = -t353 * t172 + t156;
t238 = t169 * t19 + t259;
t39 = -t214 * t77 + t324 * t71;
t7 = t214 * t32 + t71 * t260 - t276 * t77 + t324 * t38;
t230 = -t263 + t284;
t229 = t188 * t48 + t344 * t92;
t228 = -pkin(8) * t170 + t130 * t172;
t125 = t185 * t264 - t270;
t56 = t143 * t185 - t180 * t264 + t214 * t286;
t221 = t124 * t47 - t125 * t48 - t233 * t57 + t56 * t92;
t220 = t124 * t170 + t169 * t57 + t181 * t92 - t234 * t48;
t112 = qJD(3) * t149 + t225;
t80 = t231 * pkin(4) + t112;
t70 = pkin(4) * t269 + t98;
t208 = -pkin(4) * t324 - pkin(5);
t203 = pkin(4) * t214 + qJ(6);
t174 = t352 ^ 2;
t135 = pkin(5) * t188 - qJ(6) * t189 + t209;
t86 = t169 * t181 - t138;
t58 = pkin(5) * t124 - qJ(6) * t125 + t113;
t52 = pkin(4) * t153 + t62;
t50 = -qJD(4) * t85 + t256;
t36 = pkin(5) * t234 - t39;
t35 = -qJ(6) * t234 + t314;
t31 = t245 - t307;
t15 = t169 * qJ(6) + t20;
t14 = -t169 * pkin(5) + t283;
t13 = -t189 * t47 - t233 * t280;
t12 = t57 * pkin(5) + t56 * qJ(6) - t125 * qJD(6) + t80;
t11 = -t125 * t47 - t233 * t56;
t10 = t125 * t170 - t169 * t56 + t181 * t233 + t234 * t47;
t9 = t48 * pkin(5) + t47 * qJ(6) - qJD(6) * t233 + t70;
t6 = -t181 * pkin(5) - t8;
t5 = qJ(6) * t181 - qJD(6) * t234 + t7;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, qJ(2) * t240, t178 * t180 + t185 * t222, -t185 * t170 - t178 * t181 + t180 * t352 + t222 * t234, t180 * qJD(3), -t181 * t352 - t138, -t181 * qJD(3), 0, -qJD(3) * t112 + t170 * t204 + t181 * t192, -t111 * qJD(3) + t192 * t180 + t204 * t222, t111 * t352 + t112 * t178 - t137 * t181 - t149 * t170 - t180 * t336 + t98 * t185 - t222 * t335 + t234 * t97, t111 * t137 - t112 * t336 + t149 * t97 - t313, -t105 * t291 + t153 * t230 -(-t217 * t152 + t153 * t215) * t180 + (-t247 + t303 + (-t152 * t215 - t153 * t217) * qJD(4)) * t185, t105 * t234 + t153 * t181 + t156 * t185 + t172 * t230, t101 * t185 - t152 * t231, t152 * t181 - t172 * t231 - t185 * t287 + t234 * t269, t172 * t181 - t138, -t112 * t152 + t130 * t231 + t84 * t170 + t50 * t172 + t75 * t181 - t234 * t42 - t269 * t335 + t292 * t98, t105 * t335 + t112 * t153 + t130 * t230 - t170 * t85 - t172 * t49 - t181 * t76 + t234 * t41 + t291 * t98, t49 * t152 - t85 * t269 - t50 * t153 + t84 * t105 - t242 * t180 + (-t41 * t215 - t42 * t217 + (t215 * t75 - t217 * t76) * qJD(4)) * t185, t112 * t130 + t41 * t85 + t42 * t84 + t49 * t76 + t50 * t75 - t313, t11, t221, t10, t243, -t220, t86, t113 * t48 + t124 * t70 + t169 * t8 + t170 * t39 + t181 * t19 - t234 * t4 + t57 * t87 + t80 * t92, -t113 * t47 + t125 * t70 - t169 * t7 - t170 * t314 - t181 * t20 + t233 * t80 - t234 * t259 - t56 * t87, t124 * t259 - t125 * t4 + t19 * t56 - t20 * t57 - t233 * t8 - t314 * t48 + t39 * t47 - t7 * t92, t113 * t70 + t19 * t8 + t20 * t7 - t259 * t314 + t39 * t4 + t80 * t87, t11, t10, -t221, t86, t220, t243, t12 * t92 + t124 * t9 - t14 * t181 - t169 * t6 - t170 * t36 + t2 * t234 + t44 * t57 + t48 * t58, -t1 * t124 + t125 * t2 - t14 * t56 - t15 * t57 + t233 * t6 - t35 * t48 - t36 * t47 - t5 * t92, -t1 * t234 - t12 * t233 - t125 * t9 + t15 * t181 + t169 * t5 + t170 * t35 + t44 * t56 + t47 * t58, t1 * t35 + t12 * t44 + t14 * t6 + t15 * t5 + t2 * t36 + t58 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t258, -qJ(2) * t258, 0, 0, 0, 0, 0, 0, 0.2e1 * t334, 0.2e1 * t222, -t174 - t327, -t137 * t352 + t178 * t336, 0, 0, 0, 0, 0, 0, t239 + t299, -t172 ^ 2 * t217 - t287 - t297 (t105 - t300) * t217 + t346 + t282, -t130 * t178 + t241 * t215 + t217 * t338, 0, 0, 0, 0, 0, 0, t360, -t355, t358, -t178 * t87 - t188 * t4 - t189 * t259 - t19 * t344 - t20 * t280, 0, 0, 0, 0, 0, 0, t360, t358, t355, t1 * t189 + t14 * t344 - t15 * t280 - t178 * t44 + t188 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t293, -t174 + t327, 0, t293, 0, 0, -t356 * t178, -t356 * t352, 0, 0, t153 * t254 - t303 (-t105 - t300) * t217 - t346 + t282, t172 * t254 + t287 - t297, -t152 * t295 - t247, t239 - t299, -t172 * t178, -pkin(3) * t269 - t98 * t217 - t75 * t178 + t137 * t152 + (-pkin(8) * t277 - t82) * t172 + t228 * t215, pkin(3) * t105 - t137 * t153 + t178 * t76 + t98 * t215 + (pkin(8) * t278 + t83) * t172 + t228 * t217, -t83 * t152 + t82 * t153 + ((t153 * qJD(4) - t269) * pkin(8) + t241) * t217 + ((-t105 - t275) * pkin(8) - t338) * t215, -pkin(3) * t98 - t130 * t137 - t75 * t82 - t76 * t83 + (-qJD(4) * t242 - t42 * t215 + t41 * t217) * pkin(8), t13, -t357, t31, t229, t359, -t296, -t169 * t317 - t178 * t19 + t188 * t70 + t209 * t48 + t249 * t92 + t344 * t87 + t302, t169 * t316 + t178 * t20 + t189 * t70 - t209 * t47 + t233 * t249 - t280 * t87 - t301, t188 * t259 - t189 * t4 + t19 * t280 - t20 * t344 + t233 * t317 + t34 * t92 + t251, -t151 * t259 - t19 * t317 - t20 * t316 + t209 * t70 + t232 * t4 + t249 * t87, t13, t31, t357, -t296, -t359, t229, t135 * t48 + t14 * t178 + t169 * t318 + t188 * t9 - t315 * t92 + t344 * t44 + t302, -t1 * t188 - t14 * t280 - t15 * t344 + t189 * t2 - t233 * t318 + t28 * t92 + t251, t135 * t47 - t15 * t178 + t169 * t319 - t189 * t9 + t233 * t315 + t280 * t44 + t301, t1 * t151 + t135 * t9 - t14 * t318 + t15 * t319 - t2 * t232 - t315 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t298, -t152 ^ 2 + t153 ^ 2, -t152 * t172 - t105, t298, t153 * t172 - t269, t170, -t130 * t153 + t338, -t130 * t152 - t241, 0, 0, t348, -t272, t341, -t348, -t347, t170, t26 * t169 + (-t153 * t92 - t169 * t276 + t170 * t324) * pkin(4) + t331, t27 * t169 + t350 + (-t153 * t233 - t169 * t260 - t170 * t214) * pkin(4) + t259, -t19 * t92 + t20 * t233 - t26 * t233 + t27 * t92 + (t324 * t47 - t214 * t48 + (t214 * t233 - t324 * t92) * qJD(5)) * pkin(4), t19 * t26 - t20 * t27 + (t324 * t4 - t153 * t87 - t214 * t259 + (-t19 * t214 + t324 * t20) * qJD(5)) * pkin(4), t348, t341, t272, t170, t347, -t348, -t169 * t248 - t170 * t208 - t52 * t92 - t227, -t203 * t48 - t208 * t47 + (t15 + t248) * t233 + (t14 - t304) * t92, t169 * t304 + t170 * t203 + t233 * t52 + t1 - t351, t1 * t203 + t14 * t248 + t15 * t304 + t2 * t208 - t44 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t348, -t272, t341, -t348, -t347, t170, t312 + t331, t238 + t350, 0, 0, t348, t341, t272, t170, t347, -t348, -t62 * t92 + t164 - t227 + t312, pkin(5) * t47 - qJ(6) * t48 + (t15 - t20) * t233 + (t14 - t283) * t92, t233 * t62 + 0.2e1 * t159 + 0.2e1 * t161 - t238 - t351, -pkin(5) * t2 + qJ(6) * t1 - t14 * t20 + t15 * t283 - t44 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t348 - t334, t341, -t169 ^ 2 - t328, -t15 * t169 + t227;];
tauc_reg  = t3;