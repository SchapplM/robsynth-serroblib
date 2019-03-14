% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:34:26
% EndTime: 2019-03-09 00:34:46
% DurationCPUTime: 8.25s
% Computational Cost: add. (8418->511), mult. (22581->721), div. (0->0), fcn. (18316->12), ass. (0->241)
t178 = sin(pkin(7));
t184 = sin(qJ(3));
t297 = t178 * t184;
t171 = pkin(9) * t297;
t180 = cos(pkin(7));
t188 = cos(qJ(3));
t189 = cos(qJ(2));
t287 = t188 * t189;
t185 = sin(qJ(2));
t292 = t184 * t185;
t206 = -t180 * t292 + t287;
t179 = sin(pkin(6));
t279 = qJD(1) * t179;
t293 = t180 * t188;
t285 = t206 * t279 - (pkin(2) * t293 - t171) * qJD(3);
t227 = pkin(3) * t184 - pkin(10) * t188;
t204 = t227 * qJD(3);
t257 = t185 * t279;
t344 = (t204 - t257) * t178;
t273 = qJD(2) * t188;
t170 = t178 * t273;
t222 = t170 - qJD(4);
t296 = t178 * t188;
t262 = pkin(9) * t296;
t137 = t262 + (pkin(2) * t184 + pkin(10)) * t180;
t228 = -pkin(3) * t188 - pkin(10) * t184;
t138 = (-pkin(2) + t228) * t178;
t183 = sin(qJ(4));
t187 = cos(qJ(4));
t267 = qJD(4) * t187;
t269 = qJD(4) * t183;
t343 = -t137 * t269 + t138 * t267 + t344 * t183 - t285 * t187;
t290 = t185 * t188;
t291 = t184 * t189;
t208 = t180 * t290 + t291;
t294 = t180 * t184;
t283 = -t208 * t279 + (pkin(2) * t294 + t262) * qJD(3);
t272 = qJD(3) * t184;
t255 = t178 * t272;
t342 = -pkin(11) * t255 - t343;
t149 = -t187 * t180 + t183 * t297;
t271 = qJD(3) * t188;
t254 = t178 * t271;
t111 = -qJD(4) * t149 + t187 * t254;
t150 = t180 * t183 + t187 * t297;
t112 = qJD(4) * t150 + t183 * t254;
t341 = -pkin(4) * t112 + pkin(11) * t111 - t283;
t181 = cos(pkin(6));
t278 = qJD(1) * t181;
t258 = t178 * t278;
t277 = qJD(2) * t178;
t154 = pkin(9) * t277 + t257;
t311 = qJD(2) * pkin(2);
t160 = t189 * t279 + t311;
t332 = t188 * t154 + t160 * t294;
t92 = t184 * t258 + t332;
t339 = -t92 - t222 * (pkin(4) * t183 - pkin(11) * t187);
t340 = pkin(10) * t269;
t338 = t137 * t267 + t138 * t269 - t285 * t183 - t344 * t187;
t275 = qJD(2) * t180;
t240 = qJD(3) + t275;
t274 = qJD(2) * t184;
t256 = t178 * t274;
t193 = -t183 * t256 + t187 * t240;
t263 = qJD(2) * qJD(3);
t249 = t178 * t263;
t229 = t188 * t249;
t214 = t187 * t229;
t191 = qJD(4) * t193 + t214;
t337 = -qJD(5) * t222 + t191;
t336 = t183 * t170 - t269;
t125 = qJD(5) - t193;
t182 = sin(qJ(5));
t186 = cos(qJ(5));
t79 = t240 * pkin(10) + t92;
t169 = t180 * t278;
t98 = t169 + (t228 * qJD(2) - t160) * t178;
t44 = t183 * t98 + t187 * t79;
t38 = -t222 * pkin(11) + t44;
t131 = t183 * t240 + t187 * t256;
t145 = t184 * t154;
t299 = t160 * t180;
t91 = (t258 + t299) * t188 - t145;
t78 = -t240 * pkin(3) - t91;
t46 = -pkin(4) * t193 - t131 * pkin(11) + t78;
t18 = t182 * t46 + t186 * t38;
t12 = qJ(6) * t125 + t18;
t215 = t183 * t229;
t104 = t131 * qJD(4) + t215;
t323 = pkin(5) * t104;
t107 = (t204 + t257) * t277;
t196 = t206 * qJD(2);
t233 = t181 * t254;
t61 = (t160 * t293 - t145) * qJD(3) + (t179 * t196 + t233) * qJD(1);
t201 = -t183 * t107 - t187 * t61 - t98 * t267 + t79 * t269;
t230 = t184 * t249;
t15 = pkin(11) * t230 - t201;
t265 = qJD(5) * t186;
t266 = qJD(5) * t182;
t276 = qJD(2) * t179;
t250 = qJD(1) * t276;
t216 = t180 * t185 * t250;
t232 = t181 * t255;
t251 = qJD(3) * t299;
t29 = t104 * pkin(4) - pkin(11) * t191 + qJD(1) * t232 + t154 * t271 + t184 * t251 + t188 * t216 + t250 * t291;
t4 = -t182 * t15 + t186 * t29 - t38 * t265 - t46 * t266;
t2 = -t4 - t323;
t335 = -t12 * t125 + t2;
t313 = -pkin(4) * t255 + t338;
t136 = t171 + (-pkin(2) * t188 - pkin(3)) * t180;
t82 = pkin(4) * t149 - pkin(11) * t150 + t136;
t284 = t187 * t137 + t183 * t138;
t84 = -pkin(11) * t296 + t284;
t334 = t341 * t182 + t186 * t342 - t82 * t265 + t84 * t266;
t333 = t184 * t188;
t234 = t187 * t170;
t331 = t234 - t267;
t168 = -pkin(4) * t187 - pkin(11) * t183 - pkin(3);
t140 = t227 * t277;
t305 = t183 * t140 + t187 * t91;
t59 = pkin(11) * t256 + t305;
t330 = -t168 * t265 - t182 * t339 + t186 * t59;
t100 = t131 * t182 + t186 * t222;
t102 = t186 * t131 - t182 * t222;
t43 = -t183 * t79 + t187 * t98;
t37 = t222 * pkin(4) - t43;
t19 = t100 * pkin(5) - t102 * qJ(6) + t37;
t321 = pkin(11) * t104;
t329 = t125 * t19 - t321;
t218 = t182 * t82 + t186 * t84;
t328 = -qJD(5) * t218 + t182 * t342 - t341 * t186;
t327 = t102 ^ 2;
t326 = t125 ^ 2;
t190 = qJD(2) ^ 2;
t325 = qJ(6) * t112 + qJD(6) * t149 - t334;
t324 = pkin(5) * t112 + t328;
t322 = pkin(10) * t187;
t248 = t187 * t107 - t183 * t61 - t79 * t267 - t98 * t269;
t16 = -pkin(4) * t230 - t248;
t51 = t131 * t266 - t182 * t230 - t337 * t186;
t52 = t131 * t265 + t337 * t182 - t186 * t230;
t5 = pkin(5) * t52 + qJ(6) * t51 - qJD(6) * t102 + t16;
t320 = t182 * t5;
t319 = t186 * t5;
t260 = t182 * t296;
t114 = t150 * t186 - t260;
t288 = t186 * t188;
t113 = t150 * t182 + t178 * t288;
t63 = qJD(5) * t113 - t186 * t111 - t182 * t255;
t64 = -qJD(5) * t260 + t111 * t182 + t150 * t265 - t186 * t255;
t318 = pkin(5) * t64 + qJ(6) * t63 - qJD(6) * t114 + t313;
t264 = qJD(5) * t187;
t268 = qJD(4) * t186;
t317 = qJD(6) * t187 - (-t182 * t264 - t183 * t268) * pkin(10) + t330 + t336 * qJ(6);
t316 = -t168 * t266 + (t59 + t340) * t182 + (-t264 * pkin(10) + t339) * t186 - t336 * pkin(5);
t122 = t182 * t234 - t186 * t256;
t123 = (t182 * t184 + t187 * t288) * t277;
t223 = pkin(5) * t182 - qJ(6) * t186;
t213 = pkin(10) + t223;
t224 = pkin(5) * t186 + qJ(6) * t182;
t87 = t183 * t91;
t58 = -pkin(4) * t256 - t140 * t187 + t87;
t315 = pkin(5) * t122 - qJ(6) * t123 + t58 - (qJD(5) * t224 - qJD(6) * t186) * t183 - t213 * t267;
t90 = pkin(4) * t131 - pkin(11) * t193;
t314 = t182 * t90 + t186 * t43;
t310 = t102 * t19;
t308 = t16 * t182;
t307 = t182 * t51;
t306 = -qJD(6) * t182 + t125 * t223 - t44;
t304 = qJ(6) * t104;
t303 = t100 * t125;
t302 = t102 * t100;
t301 = t102 * t125;
t300 = t104 * t182;
t175 = t178 ^ 2;
t298 = t175 * t190;
t295 = t179 * t190;
t289 = t186 * t104;
t17 = -t182 * t38 + t186 * t46;
t286 = qJD(6) - t17;
t281 = t182 * t168 + t186 * t322;
t280 = t184 ^ 2 - t188 ^ 2;
t270 = qJD(4) * t182;
t261 = pkin(11) * t266;
t259 = t185 * t295;
t253 = t125 * t266;
t252 = t178 * t180 * t190;
t247 = -t183 * t137 + t138 * t187;
t246 = t188 * t222;
t245 = t125 * t186;
t244 = t222 * t178;
t243 = qJD(4) * t222;
t241 = 0.2e1 * t175 * t263;
t239 = qJD(3) + 0.2e1 * t275;
t238 = t175 * t259;
t236 = t178 * t185 * t276;
t83 = pkin(4) * t296 - t247;
t225 = t186 * t267 - t123;
t11 = -pkin(5) * t125 + t286;
t221 = t11 * t186 - t12 * t182;
t219 = -t182 * t84 + t186 * t82;
t209 = t180 * t287 - t292;
t109 = -t179 * t209 - t181 * t296;
t207 = t180 * t291 + t290;
t110 = t179 * t207 + t181 * t297;
t148 = -t178 * t179 * t189 + t180 * t181;
t81 = t110 * t187 + t148 * t183;
t217 = t109 * t186 - t182 * t81;
t50 = t109 * t182 + t186 * t81;
t80 = t110 * t183 - t148 * t187;
t211 = t125 * t18 + t4;
t124 = -t160 * t178 + t169;
t210 = t124 * t178 - t175 * t311;
t3 = t186 * t15 + t182 * t29 + t46 * t265 - t38 * t266;
t202 = -t125 * t265 - t300;
t200 = t125 * t37 - t321;
t74 = t233 + (qJD(3) * t209 + t196) * t179;
t33 = qJD(4) * t81 + t183 * t74 - t187 * t236;
t34 = -qJD(4) * t80 + t183 * t236 + t187 * t74;
t197 = t208 * qJD(2);
t73 = t232 + (qJD(3) * t207 + t197) * t179;
t8 = qJD(5) * t50 + t182 * t34 - t186 * t73;
t199 = t33 * t100 + t104 * t217 - t125 * t8 + t80 * t52;
t195 = qJD(3) * t154 + t216;
t9 = qJD(5) * t217 + t182 * t73 + t186 * t34;
t194 = t102 * t33 - t104 * t50 - t125 * t9 - t51 * t80;
t192 = -t124 * t277 - t251 + (-qJD(3) * t178 * t181 - t189 * t276) * qJD(1);
t163 = -pkin(4) - t224;
t139 = t213 * t183;
t116 = -t168 * t186 + (pkin(10) * t182 + pkin(5)) * t187;
t115 = -qJ(6) * t187 + t281;
t62 = t332 * qJD(3) + (t179 * t197 + t232) * qJD(1);
t55 = pkin(5) * t102 + qJ(6) * t100;
t39 = pkin(5) * t113 - qJ(6) * t114 + t83;
t32 = -pkin(5) * t149 - t219;
t31 = qJ(6) * t149 + t218;
t28 = -t51 + t303;
t21 = -pkin(5) * t131 + t182 * t43 - t186 * t90;
t20 = qJ(6) * t131 + t314;
t1 = qJD(6) * t125 + t3 + t304;
t6 = [0, 0, -t259, -t189 * t295, 0, 0, 0, 0, 0, t148 * t230 - t188 * t238 - t73 * t240, t148 * t229 + t184 * t238 - t74 * t240, 0, 0, 0, 0, 0, t109 * t104 - t193 * t73 + t222 * t33 - t230 * t80, t109 * t191 + t73 * t131 + t222 * t34 - t230 * t81, 0, 0, 0, 0, 0, t199, t194, t199, -t100 * t9 + t102 * t8 + t217 * t51 - t50 * t52, -t194, t1 * t50 + t11 * t8 + t12 * t9 + t19 * t33 - t2 * t217 + t5 * t80; 0, 0, 0, 0, t241 * t333, -t280 * t241, t239 * t254, -t239 * t255, 0 (-t283 * qJD(2) - t62) * t180 + (t184 * t210 - t283) * qJD(3) (t285 * qJD(2) - t61) * t180 + (t188 * t210 + t285) * qJD(3), t131 * t111 + t150 * t191, -t150 * t104 + t111 * t193 - t131 * t112 - t149 * t191, -t111 * t222 + t131 * t255 + t150 * t230 - t191 * t296, t112 * t222 + (t104 * t188 + (-qJD(2) * t149 + t193) * t272) * t178 (-t175 * t273 - t244) * t272, t136 * t104 + t62 * t149 + t78 * t112 - t283 * t193 + (-t248 * t188 + (t247 * qJD(2) + t43) * t272) * t178 + t338 * t222, t78 * t111 + t283 * t131 + t136 * t191 + t62 * t150 - t201 * t296 + t343 * t222 - t230 * t284 - t255 * t44, -t102 * t63 - t114 * t51, t100 * t63 - t102 * t64 + t113 * t51 - t114 * t52, t102 * t112 + t104 * t114 - t125 * t63 - t149 * t51, -t100 * t112 - t104 * t113 - t125 * t64 - t149 * t52, t104 * t149 + t112 * t125, t313 * t100 + t219 * t104 + t17 * t112 + t16 * t113 + t328 * t125 + t4 * t149 + t37 * t64 + t83 * t52, t313 * t102 - t218 * t104 - t18 * t112 + t16 * t114 + t334 * t125 - t3 * t149 - t37 * t63 - t83 * t51, t318 * t100 - t104 * t32 - t11 * t112 + t113 * t5 + t324 * t125 - t149 * t2 + t19 * t64 + t39 * t52, -t1 * t113 - t325 * t100 - t324 * t102 - t11 * t63 + t114 * t2 - t12 * t64 - t31 * t52 - t32 * t51, t1 * t149 - t318 * t102 + t104 * t31 + t112 * t12 - t114 * t5 + t325 * t125 + t19 * t63 + t39 * t51, t1 * t31 - t324 * t11 + t325 * t12 + t318 * t19 + t2 * t32 + t39 * t5; 0, 0, 0, 0, -t298 * t333, t280 * t298, -t188 * t252, t184 * t252, 0, t192 * t184 - t195 * t188 + t92 * t240, t195 * t184 + t192 * t188 + t91 * t240, -qJD(4) * t183 ^ 2 * t256 + ((qJD(4) * t240 + t229) * t183 - t222 * t131) * t187, -t183 * t104 + t336 * t131 + t187 * t191 - t193 * t331, -t187 * t243 + (t187 * t246 + (t183 * qJD(3) - t131) * t184) * t277, t183 * t243 + (-t183 * t246 + (t187 * qJD(3) - t193) * t184) * t277, t244 * t274, -pkin(3) * t104 + t78 * t269 - t87 * t222 + t92 * t193 + (pkin(10) * t243 + t140 * t222 - t62) * t187 + (-t184 * t43 + (-pkin(10) * t272 - t188 * t78) * t183) * t277, -pkin(3) * t191 - t92 * t131 + t62 * t183 - t230 * t322 + t44 * t256 - t331 * t78 + (-t305 - t340) * t222, -t183 * t186 * t51 + (-t183 * t266 + t225) * t102, t100 * t123 + t102 * t122 + (-t100 * t186 - t102 * t182) * t267 + (t307 - t186 * t52 + (t100 * t182 - t102 * t186) * qJD(5)) * t183, t187 * t51 + t225 * t125 + (-t102 * t222 - t253 + t289) * t183, t187 * t52 + (-t182 * t267 + t122) * t125 + (t100 * t222 + t202) * t183, -t125 * t183 * t222 - t104 * t187, t168 * t289 - t58 * t100 - t37 * t122 + (t339 * t186 + (-qJD(5) * t168 + t59) * t182) * t125 + (t37 * t270 - t4 + (qJD(4) * t100 + t202) * pkin(10)) * t187 + (t37 * t265 + t308 - t222 * t17 + (t125 * t270 + t52) * pkin(10)) * t183, -t281 * t104 - t58 * t102 - t37 * t123 + t330 * t125 + (t37 * t268 + t3 + (qJD(4) * t102 + t253) * pkin(10)) * t187 + (-t37 * t266 + t16 * t186 + t222 * t18 + (t125 * t268 - t51) * pkin(10)) * t183, -t104 * t116 - t122 * t19 + t139 * t52 + (t19 * t270 + t2) * t187 + t316 * t125 - t315 * t100 + (t11 * t222 + t19 * t265 + t320) * t183, -t11 * t123 - t115 * t52 - t116 * t51 + t12 * t122 - t316 * t102 + t317 * t100 + t221 * t267 + (-t1 * t182 + t186 * t2 + (-t11 * t182 - t12 * t186) * qJD(5)) * t183, t104 * t115 + t123 * t19 + t139 * t51 + (-t19 * t268 - t1) * t187 - t317 * t125 + t315 * t102 + (-t12 * t222 + t19 * t266 - t319) * t183, t1 * t115 - t316 * t11 + t116 * t2 - t317 * t12 + t139 * t5 - t315 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131 * t193, t131 ^ 2 - t193 ^ 2, t193 * t170 + t214, -t131 * t170 - t215, t230, -t78 * t131 - t222 * t44 + t248, -t193 * t78 - t222 * t43 + t201, t102 * t245 - t307 (-t51 - t303) * t186 + (-t52 - t301) * t182, -t102 * t131 + t125 * t245 + t300, t100 * t131 - t326 * t182 + t289, -t125 * t131, -pkin(4) * t52 - t44 * t100 - t17 * t131 + (-t16 + (-pkin(11) * qJD(5) - t90) * t125) * t186 + (t43 * t125 + t200) * t182, pkin(4) * t51 - t44 * t102 + t18 * t131 + t308 + (t261 + t314) * t125 + t200 * t186, t11 * t131 + t163 * t52 - t319 + (-pkin(11) * t265 + t21) * t125 + t306 * t100 + t329 * t182, t100 * t20 - t102 * t21 + (t1 + t125 * t11 + (qJD(5) * t102 - t52) * pkin(11)) * t186 + ((qJD(5) * t100 - t51) * pkin(11) + t335) * t182, -t12 * t131 + t163 * t51 - t320 + (-t20 - t261) * t125 - t306 * t102 - t329 * t186, -t11 * t21 - t12 * t20 + t163 * t5 + t306 * t19 + (qJD(5) * t221 + t1 * t186 + t182 * t2) * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, -t100 ^ 2 + t327, t28, -t52 + t301, t104, -t102 * t37 + t211, t100 * t37 + t125 * t17 - t3, -t100 * t55 + t211 - t310 + 0.2e1 * t323, pkin(5) * t51 - qJ(6) * t52 + (t12 - t18) * t102 + (t11 - t286) * t100, 0.2e1 * t304 - t100 * t19 + t102 * t55 + (0.2e1 * qJD(6) - t17) * t125 + t3, -pkin(5) * t2 + qJ(6) * t1 - t11 * t18 + t12 * t286 - t19 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104 + t302, t28, -t326 - t327, t310 + t335;];
tauc_reg  = t6;