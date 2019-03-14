% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRP12_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:53:32
% EndTime: 2019-03-09 12:53:50
% DurationCPUTime: 6.86s
% Computational Cost: add. (8808->530), mult. (19440->664), div. (0->0), fcn. (12100->6), ass. (0->260)
t212 = cos(qJ(4));
t210 = sin(qJ(4));
t286 = t210 * qJD(2);
t213 = cos(qJ(2));
t294 = qJD(1) * t213;
t146 = -t212 * t294 - t286;
t268 = t210 * t294;
t292 = qJD(2) * t212;
t147 = -t268 + t292;
t209 = sin(qJ(5));
t344 = cos(qJ(5));
t227 = t209 * t146 + t147 * t344;
t148 = t209 * t212 + t210 * t344;
t283 = qJD(1) * qJD(2);
t191 = t213 * t283;
t211 = sin(qJ(2));
t295 = qJD(1) * t211;
t187 = qJD(4) + t295;
t176 = qJD(5) + t187;
t267 = t344 * qJD(5);
t275 = t212 * t295;
t287 = qJD(5) * t209;
t290 = qJD(4) * t210;
t309 = t209 * t210;
t324 = -t275 * t344 + t295 * t309 + t209 * t290 + t210 * t287 - (qJD(4) * t344 + t267) * t212;
t262 = t324 * t176;
t247 = t148 * t191 - t262;
t369 = qJD(2) * t227 + t247;
t288 = qJD(4) * t213;
t270 = t212 * t288;
t274 = t211 * t286;
t223 = -t270 + t274;
t282 = qJD(2) * qJD(4);
t106 = -qJD(1) * t223 + t210 * t282;
t192 = t212 * t282;
t272 = t210 * t288;
t273 = t211 * t292;
t224 = t272 + t273;
t219 = qJD(1) * t224 - t192;
t46 = qJD(5) * t227 - t209 * t106 - t344 * t219;
t87 = -t344 * t146 + t147 * t209;
t228 = t46 * t148 - t324 * t87;
t363 = -t176 * t227 + t46;
t351 = qJD(4) + qJD(5);
t95 = t351 * t148;
t323 = t148 * t295 + t95;
t368 = t323 * t176;
t345 = pkin(3) + pkin(7);
t196 = pkin(7) * t294;
t155 = pkin(3) * t294 + t196;
t206 = qJD(2) * qJ(3);
t135 = t206 + t155;
t97 = -pkin(4) * t146 + t135;
t38 = pkin(5) * t87 - qJ(6) * t227 + t97;
t367 = t38 * t87;
t366 = t87 * t97;
t364 = t87 * t227;
t195 = pkin(7) * t295;
t361 = qJD(3) + t195;
t346 = t227 ^ 2;
t281 = t87 ^ 2 - t346;
t45 = t344 * t106 - t146 * t267 + t147 * t287 - t209 * t219;
t360 = t176 * t87 - t45;
t56 = pkin(5) * t227 + qJ(6) * t87;
t359 = -0.2e1 * t283;
t214 = -pkin(2) - pkin(8);
t339 = pkin(9) - t214;
t160 = t339 * t212;
t145 = qJD(4) * t160;
t248 = t339 * t290;
t307 = t210 * t211;
t235 = pkin(4) * t213 - pkin(9) * t307;
t200 = pkin(2) * t295;
t241 = pkin(8) * t211 - qJ(3) * t213;
t124 = qJD(1) * t241 + t200;
t83 = -t210 * t124 + t212 * t155;
t70 = qJD(1) * t235 + t83;
t84 = t212 * t124 + t210 * t155;
t75 = pkin(9) * t275 + t84;
t159 = t339 * t210;
t99 = -t159 * t344 - t209 * t160;
t336 = qJD(5) * t99 + (-t248 + t70) * t344 + (-t145 - t75) * t209;
t343 = t38 * t227;
t356 = t97 * t227;
t266 = t211 * t283;
t186 = pkin(2) * t266;
t285 = t211 * qJD(3);
t221 = qJD(2) * t241 - t285;
t101 = qJD(1) * t221 + t186;
t305 = t211 * qJ(3);
t264 = -pkin(1) - t305;
t352 = t213 * t214;
t142 = t264 + t352;
t113 = t142 * qJD(1);
t284 = pkin(3) * t295 + t361;
t118 = t214 * qJD(2) + t284;
t185 = pkin(7) * t191;
t141 = pkin(3) * t191 + t185;
t289 = qJD(4) * t212;
t36 = t212 * t101 - t113 * t290 + t118 * t289 + t210 * t141;
t72 = -t113 * t210 + t212 * t118;
t355 = -t187 * t72 + t36;
t73 = t113 * t212 + t118 * t210;
t37 = -qJD(4) * t73 - t210 * t101 + t212 * t141;
t354 = t73 * t187 + t37;
t277 = -pkin(4) * t212 - pkin(3);
t297 = pkin(4) * t289 - t277 * t295 + t361;
t171 = t345 * t211;
t150 = t210 * t171;
t94 = t212 * t142 + t150;
t207 = t211 ^ 2;
t208 = t213 ^ 2;
t296 = t207 - t208;
t276 = t344 * t212;
t149 = t276 - t309;
t11 = -t149 * t45 - t227 * t323;
t251 = t213 * t276;
t306 = t210 * t213;
t125 = t209 * t306 - t251;
t291 = qJD(2) * t213;
t293 = qJD(2) * t211;
t250 = t344 * t293;
t67 = -t209 * t274 + t212 * t250 + t213 * t95;
t350 = t176 * t67 - t211 * t46 + (qJD(1) * t125 - t87) * t291;
t151 = t212 * t171;
t263 = pkin(9) * t213 - t142;
t78 = t211 * pkin(4) + t210 * t263 + t151;
t302 = t212 * t213;
t82 = -pkin(9) * t302 + t94;
t333 = t209 * t78 + t344 * t82;
t199 = pkin(2) * t293;
t108 = t199 + t221;
t156 = t345 * t291;
t254 = -t210 * t108 + t212 * t156;
t40 = t235 * qJD(2) + (t212 * t263 - t150) * qJD(4) + t254;
t51 = t212 * t108 - t142 * t290 + t210 * t156 + t171 * t289;
t43 = pkin(9) * t224 + t51;
t8 = -qJD(5) * t333 - t209 * t43 + t344 * t40;
t162 = t176 * qJD(6);
t179 = qJ(6) * t191;
t22 = pkin(4) * t191 + t106 * pkin(9) + t37;
t27 = pkin(9) * t219 + t36;
t64 = -pkin(9) * t147 + t72;
t55 = pkin(4) * t187 + t64;
t65 = pkin(9) * t146 + t73;
t261 = -t209 * t22 - t55 * t267 - t344 * t27 + t65 * t287;
t1 = t179 + t162 - t261;
t327 = t209 * t65;
t16 = t344 * t55 - t327;
t298 = qJD(6) - t16;
t14 = -t176 * pkin(5) + t298;
t280 = t344 * t65;
t17 = t209 * t55 + t280;
t15 = t176 * qJ(6) + t17;
t253 = pkin(5) * t191;
t260 = t209 * t27 - t344 * t22 + t65 * t267 + t55 * t287;
t2 = -t253 + t260;
t349 = t1 * t148 + t14 * t323 - t2 * t149 - t15 * t324;
t348 = -t148 * t261 - t149 * t260 - t16 * t323 - t17 * t324;
t347 = t45 * t148 - t149 * t46 + t227 * t324 + t323 * t87;
t34 = t209 * t70 + t344 * t75;
t30 = qJ(6) * t294 + t34;
t226 = t209 * t159 - t160 * t344;
t59 = qJD(5) * t226 - t145 * t344 + t209 * t248;
t338 = -t30 + t59;
t337 = -pkin(5) * t294 - t336;
t335 = t34 - t59;
t334 = -t324 * pkin(5) + t323 * qJ(6) - qJD(6) * t149 + t297;
t332 = qJD(2) * pkin(2);
t326 = t210 * t72;
t24 = t344 * t64 - t327;
t322 = pkin(4) * t267 + qJD(6) - t24;
t320 = qJD(2) * t226;
t319 = qJD(2) * t99;
t154 = t345 * t293;
t205 = qJD(2) * qJD(3);
t121 = -qJD(1) * t154 + t205;
t318 = t121 * t210;
t317 = t121 * t212;
t316 = t135 * t211;
t315 = t146 * t212;
t314 = t147 * t146;
t313 = t147 * t213;
t312 = t187 * t211;
t311 = t187 * t214;
t310 = t192 * t210;
t308 = t210 * t146;
t304 = t211 * t212;
t303 = t212 * t106;
t216 = qJD(1) ^ 2;
t301 = t213 * t216;
t215 = qJD(2) ^ 2;
t300 = t215 * t211;
t299 = t215 * t213;
t190 = t210 * pkin(4) + qJ(3);
t172 = t345 * t213;
t161 = -pkin(2) * t213 + t264;
t136 = qJD(1) * t161;
t279 = t146 * t304;
t278 = t187 * t304;
t133 = pkin(4) * t302 + t172;
t271 = t187 * t289;
t269 = t176 * t294;
t259 = qJD(1) * t94 + t73;
t258 = pkin(1) * t359;
t257 = qJD(3) - t332;
t252 = t226 * t45 - t99 * t46 - t59 * t87;
t23 = t209 * t64 + t280;
t249 = pkin(4) * t287 - t23;
t246 = t149 * t191 - t368;
t242 = -t125 * t46 - t67 * t87;
t238 = t147 * t212 + t308;
t237 = -0.2e1 * qJD(2) * t136;
t236 = t187 * t210;
t234 = t212 * t266 - t192;
t233 = t16 * t176 + t261;
t232 = t17 * t176 - t260;
t47 = -t209 * t82 + t344 * t78;
t225 = -qJ(3) * t291 - t285;
t111 = qJD(1) * t225 + t186;
t129 = t199 + t225;
t229 = pkin(7) * t215 + qJD(1) * t129 + t111;
t7 = t209 * t40 + t78 * t267 - t287 * t82 + t344 * t43;
t126 = t148 * t213;
t66 = t210 * t250 + (t351 * t306 + t273) * t209 - t351 * t251;
t222 = t125 * t45 - t126 * t46 - t227 * t67 + t66 * t87;
t220 = -t11 - t228;
t100 = -pkin(4) * t272 + (-pkin(7) + t277) * t293;
t157 = pkin(7) * t266 - t205;
t158 = t195 + t257;
t168 = -t196 - t206;
t217 = -t157 * t213 + (t158 * t213 + (t168 + t196) * t211) * qJD(2);
t77 = t192 * pkin(4) + qJD(1) * t100 + t205;
t194 = -pkin(4) * t344 - pkin(5);
t189 = pkin(4) * t209 + qJ(6);
t182 = t211 * t301;
t174 = t212 * t191;
t173 = t211 * t191;
t170 = -0.2e1 * t173;
t169 = 0.2e1 * t173;
t163 = t296 * t216;
t152 = -qJ(3) * t294 + t200;
t134 = t296 * t359;
t117 = t136 * t295;
t110 = t176 * t291 + t173;
t93 = -t142 * t210 + t151;
t85 = pkin(5) * t148 - qJ(6) * t149 + t190;
t68 = -pkin(5) * t125 + qJ(6) * t126 + t133;
t52 = -t94 * qJD(4) + t254;
t50 = pkin(4) * t147 + t56;
t44 = -t211 * pkin(5) - t47;
t42 = qJ(6) * t211 + t333;
t28 = -t227 * t294 + t246;
t13 = -t67 * pkin(5) - t66 * qJ(6) + t126 * qJD(6) + t100;
t12 = t126 * t45 + t227 * t66;
t10 = t176 * t66 - t211 * t45 + (-qJD(1) * t126 + t227) * t291;
t9 = t46 * pkin(5) + t45 * qJ(6) - qJD(6) * t227 + t77;
t6 = -pkin(5) * t291 - t8;
t5 = qJ(6) * t291 + qJD(6) * t211 + t7;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, t134, t299, t170, -t300, 0, -pkin(7) * t299 + t211 * t258, pkin(7) * t300 + t213 * t258, 0, 0, 0, -t299, t300, t169, t134, t170, t217, t211 * t237 + t213 * t229, -t211 * t229 + t213 * t237, pkin(7) * t217 + t111 * t161 + t136 * t129, t106 * t306 + t147 * t223, t238 * t293 + (-t234 * t210 + t303 + (-t315 + (t147 - t268) * t210) * qJD(4)) * t213, -t187 * t270 - t106 * t211 + (t313 + (-qJD(1) * t208 + t312) * t210) * qJD(2), t146 * t224 - t219 * t302, t187 * t272 + (qJD(4) * t268 - t192) * t211 + (t146 * t213 + (t296 * qJD(1) + t312) * t212) * qJD(2), t187 * t291 + t173, t154 * t146 + t172 * t192 + t52 * t187 + (t37 + (-qJD(1) * t172 - t135) * t292) * t211 + (-t135 * t290 + t72 * qJD(2) + t317 + (qJD(2) * t93 - t172 * t290) * qJD(1)) * t213, -t172 * t106 - t154 * t147 - t51 * t187 + (t135 * t286 - t36) * t211 + (-qJD(2) * t259 - t135 * t289 - t318) * t213, t93 * t106 + t51 * t146 - t52 * t147 - t94 * t192 + (t212 * t259 - t326) * t293 + (t37 * t210 - t36 * t212 + (t210 * t259 + t212 * t72) * qJD(4)) * t213, t121 * t172 - t135 * t154 + t36 * t94 + t37 * t93 + t51 * t73 + t52 * t72, t12, -t222, t10, t242, t350, t110, t100 * t87 - t77 * t125 + t133 * t46 + t8 * t176 - t260 * t211 - t97 * t67 + (qJD(1) * t47 + t16) * t291, t100 * t227 - t77 * t126 - t133 * t45 - t7 * t176 + t261 * t211 + t97 * t66 + (-qJD(1) * t333 - t17) * t291, -t125 * t261 - t126 * t260 - t16 * t66 + t17 * t67 - t227 * t8 - t333 * t46 + t45 * t47 - t7 * t87, t100 * t97 + t133 * t77 + t16 * t8 + t17 * t7 - t260 * t47 - t261 * t333, t12, t10, t222, t110, -t350, t242, -t9 * t125 + t13 * t87 - t6 * t176 - t2 * t211 - t38 * t67 + t68 * t46 + (-qJD(1) * t44 - t14) * t291, t1 * t125 - t126 * t2 + t14 * t66 + t15 * t67 + t227 * t6 - t42 * t46 - t44 * t45 - t5 * t87, t1 * t211 + t9 * t126 - t13 * t227 + t5 * t176 - t38 * t66 + t68 * t45 + (qJD(1) * t42 + t15) * t291, t1 * t42 + t13 * t38 + t14 * t6 + t15 * t5 + t2 * t44 + t68 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, t163, 0, t182, 0, 0, t216 * pkin(1) * t211, pkin(1) * t301, 0, 0, 0, 0, 0, -t182, t163, t182 ((-t168 - t206) * t211 + (-t158 + t257) * t213) * qJD(1), -t152 * t294 + t117, 0.2e1 * t205 + (t136 * t213 + t152 * t211) * qJD(1), -t157 * qJ(3) - t168 * qJD(3) - t136 * t152 + (-t168 * t211 + (-t158 - t332) * t213) * qJD(1) * pkin(7), -t147 * t236 - t303, t106 * t210 - t192 * t212 - t238 * qJD(4) + (t210 * t270 + (-t308 + (-t147 + t292) * t212) * t211) * qJD(1), -t187 * t290 + t174 + (-t187 * t307 - t313) * qJD(1), -t146 * t289 + t310 + (-t210 * t224 - t279) * qJD(1), -t271 + (-t278 + (-t146 - t286) * t213) * qJD(1), -t187 * t294, qJ(3) * t192 + t318 - t83 * t187 - t284 * t146 + (t135 * t212 - t210 * t311) * qJD(4) + ((-qJ(3) * t290 - t72) * t213 + (t316 + (-t305 + t352) * qJD(2)) * t212) * qJD(1), -qJ(3) * t106 + t317 + t84 * t187 + t284 * t147 + (-t135 * t210 - t212 * t311) * qJD(4) + (t213 * t73 + (-t214 * t291 - t316) * t210) * qJD(1), -t84 * t146 + t83 * t147 + (-t73 * t295 + t214 * t106 - t37 + (t146 * t214 - t73) * qJD(4)) * t212 + (t214 * t234 - t36 + t72 * t295 + (t72 + (t147 + t268) * t214) * qJD(4)) * t210, t121 * qJ(3) - t72 * t83 - t73 * t84 + t284 * t135 + (t36 * t210 + t37 * t212 + (t212 * t73 - t326) * qJD(4)) * t214, t11, t347, t28, t228, t262 + (-qJD(2) * t148 + t87) * t294, -t269, t77 * t148 + t190 * t46 - t324 * t97 + t297 * t87 - t336 * t176 + (-t16 + t320) * t294, t149 * t77 - t190 * t45 - t323 * t97 + t297 * t227 + t335 * t176 + (t17 - t319) * t294, t227 * t336 + t34 * t87 + t252 - t348, -t16 * t336 - t17 * t335 + t77 * t190 - t226 * t260 - t261 * t99 + t297 * t97, t11, t28, -t347, -t269, -t294 * t87 + t247, t228, t148 * t9 + t46 * t85 + t334 * t87 - t324 * t38 + t337 * t176 + (t14 + t320) * t294, -t227 * t337 + t30 * t87 + t252 - t349, -t149 * t9 + t45 * t85 - t334 * t227 + t323 * t38 + t338 * t176 + (-t15 + t319) * t294, t1 * t99 - t14 * t337 + t15 * t338 - t2 * t226 + t334 * t38 + t9 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, -t207 * t216 - t215, qJD(2) * t168 + t117 + t185, 0, 0, 0, 0, 0, 0, qJD(2) * t146 - t187 * t236 + t174, -t271 - qJD(2) * t147 + (-t213 * t286 - t278) * qJD(1), t303 - t310 + (t147 * t210 + t315) * qJD(4) + (t279 + (t272 + (t147 + t292) * t211) * t210) * qJD(1), -t135 * qJD(2) + t355 * t210 + t354 * t212, 0, 0, 0, 0, 0, 0, -qJD(2) * t87 + t246, -t369, t220, -t97 * qJD(2) + t348, 0, 0, 0, 0, 0, 0, -t368 + (t149 * t294 - t87) * qJD(2), t220, t369, -t38 * qJD(2) + t349; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t314, -t146 ^ 2 + t147 ^ 2, -t146 * t187 - t106, t314, t147 * t187 + t219, t191, -t135 * t147 + t354, -t135 * t146 - t355, 0, 0, t364, -t281, t360, -t364, -t363, t191, t23 * t176 - t356 + (-t147 * t87 - t176 * t287 + t191 * t344) * pkin(4) - t260, t24 * t176 + t366 + (-t147 * t227 - t176 * t267 - t191 * t209) * pkin(4) + t261, -t16 * t87 + t17 * t227 - t23 * t227 + t24 * t87 + (t344 * t45 - t209 * t46 + (t209 * t227 - t344 * t87) * qJD(5)) * pkin(4), t16 * t23 - t17 * t24 + (-t344 * t260 - t147 * t97 - t209 * t261 + (-t16 * t209 + t17 * t344) * qJD(5)) * pkin(4), t364, t360, t281, t191, t363, -t364, -t343 - t50 * t87 - t249 * t176 + (pkin(5) - t194) * t191 - t260, -t189 * t46 - t194 * t45 + (t15 + t249) * t227 + (t14 - t322) * t87, t176 * t322 + t189 * t191 + t227 * t50 + t1 - t367, t1 * t189 + t14 * t249 + t15 * t322 + t2 * t194 - t38 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t364, -t281, t360, -t364, -t363, t191, t232 - t356, t233 + t366, 0, 0, t364, t360, t281, t191, t363, -t364, -t56 * t87 + t232 + 0.2e1 * t253 - t343, pkin(5) * t45 - t46 * qJ(6) + (t15 - t17) * t227 + (t14 - t298) * t87, t227 * t56 + 0.2e1 * t162 + 0.2e1 * t179 - t233 - t367, -t2 * pkin(5) + t1 * qJ(6) - t14 * t17 + t15 * t298 - t38 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191 + t364, t360, -t176 ^ 2 - t346, -t15 * t176 + t2 + t343;];
tauc_reg  = t3;