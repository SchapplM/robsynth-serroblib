% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x31]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR15_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:26:25
% EndTime: 2024-09-27 22:26:32
% DurationCPUTime: 5.02s
% Computational Cost: add. (10374->381), mult. (35013->566), div. (0->0), fcn. (28421->12), ass. (0->256)
t204 = sin(qJ(3));
t205 = sin(qJ(2));
t199 = sin(pkin(5));
t301 = qJD(1) * t199;
t343 = cos(qJ(3));
t257 = t343 * t301;
t207 = cos(qJ(2));
t280 = t207 * t301;
t164 = -t204 * t280 - t205 * t257;
t203 = sin(qJ(4));
t342 = cos(qJ(4));
t281 = t205 * t301;
t367 = -t204 * t281 + t207 * t257;
t123 = t164 * t203 + t342 * t367;
t198 = sin(pkin(6));
t206 = cos(qJ(5));
t231 = -t342 * t164 + t203 * t367;
t294 = qJD(5) * t206;
t200 = cos(pkin(6));
t202 = sin(qJ(5));
t308 = t200 * t202;
t368 = -t123 * t206 + t198 * t294 + t231 * t308;
t275 = qJD(4) * t342;
t297 = qJD(4) * t203;
t234 = t204 * t207 + t343 * t205;
t351 = qJD(2) + qJD(3);
t214 = t351 * t234;
t211 = qJD(1) * t214;
t210 = t199 * t211;
t201 = cos(pkin(5));
t340 = pkin(1) * t207;
t290 = t201 * t340;
t187 = qJD(1) * t290;
t345 = pkin(8) + pkin(9);
t286 = t199 * t345;
t260 = t205 * t286;
t152 = -qJD(1) * t260 + t187;
t300 = qJD(1) * t201;
t189 = qJD(2) + t300;
t135 = pkin(2) * t189 + t152;
t180 = qJD(2) * t187;
t246 = qJD(2) * t260;
t143 = -qJD(1) * t246 + t180;
t341 = pkin(1) * t205;
t291 = t201 * t341;
t155 = (-t207 * t286 - t291) * qJD(2);
t144 = qJD(1) * t155;
t310 = t199 * t207;
t160 = t345 * t310 + t291;
t153 = t160 * qJD(1);
t276 = qJD(3) * t343;
t298 = qJD(3) * t204;
t259 = t135 * t276 + t343 * t143 + t204 * t144 - t153 * t298;
t55 = -pkin(10) * t210 + t259;
t125 = t367 * t351;
t148 = t343 * t153;
t236 = -t204 * t135 - t148;
t265 = -t204 * t143 + t343 * t144;
t217 = t236 * qJD(3) + t265;
t56 = -t125 * pkin(10) + t217;
t184 = qJD(3) + t189;
t159 = t164 * pkin(10);
t145 = t204 * t153;
t266 = t343 * t135 - t145;
t93 = t159 + t266;
t88 = pkin(3) * t184 + t93;
t337 = t367 * pkin(10);
t94 = -t236 + t337;
t269 = -t203 * t56 - t88 * t275 + t94 * t297 - t342 * t55;
t213 = t199 * t214;
t209 = t342 * t213;
t267 = qJD(1) * t209 + t203 * t125;
t66 = t231 * qJD(4) + t267;
t329 = t200 * t66;
t13 = -pkin(11) * t329 - t269;
t92 = t342 * t94;
t242 = -t203 * t88 - t92;
t272 = -t203 * t55 + t342 * t56;
t219 = t242 * qJD(4) + t272;
t338 = pkin(11) * t200;
t65 = t342 * t125 + t164 * t297 - t203 * t210 + t275 * t367;
t14 = -t65 * t338 + t219;
t238 = qJD(4) + t184;
t114 = t231 * t338;
t90 = t203 * t94;
t271 = t342 * t88 - t90;
t41 = t271 - t114;
t40 = t238 * pkin(4) + t41;
t173 = (-pkin(2) * t207 - pkin(1)) * t199;
t171 = qJD(1) * t173;
t132 = -pkin(3) * t367 + t171;
t355 = t231 * t198;
t67 = -pkin(4) * t123 - pkin(11) * t355 + t132;
t249 = t198 * t67 + t200 * t40;
t224 = t198 * t238;
t309 = t200 * t123;
t220 = -t224 - t309;
t39 = -t220 * pkin(11) - t242;
t15 = -t202 * t39 + t249 * t206;
t299 = qJD(2) * t199;
t279 = t205 * t299;
t255 = qJD(1) * t279;
t179 = pkin(2) * t255;
t108 = pkin(3) * t210 + t179;
t193 = t198 * pkin(11);
t33 = t66 * pkin(4) - t65 * t193 + t108;
t2 = (t14 * t200 + t198 * t33) * t202 + t13 * t206 + t15 * qJD(5);
t28 = -t198 * t40 + t200 * t67;
t312 = t198 * t202;
t6 = -t14 * t198 + t200 * t33;
t237 = -t2 * t200 + t368 * t28 + t6 * t312;
t295 = qJD(5) * t202;
t307 = t200 * t206;
t364 = t123 * t307 + t206 * t224;
t62 = t66 * t308;
t24 = t364 * qJD(5) + t206 * t65 - t231 * t295 - t62;
t304 = t206 * t231;
t349 = t220 * t202 - t304;
t25 = -t349 * qJD(5) + t202 * t65 + t66 * t307;
t311 = t198 * t206;
t76 = t202 * t231 - t364;
t366 = t24 * t311 - t25 * t312 - t368 * t76;
t320 = t123 * t198;
t106 = -t200 * t238 - qJD(5) + t320;
t194 = t198 ^ 2;
t327 = t202 * t66;
t8 = -t368 * t106 + t194 * t327 + t24 * t200 + t349 * t355;
t10 = t24 * t312 - t368 * t349;
t365 = t123 * t202;
t319 = t231 * t123;
t51 = -t123 ^ 2 + t231 ^ 2;
t49 = -t123 * t238 + t65;
t362 = t194 * t206 * t66 - t25 * t200 + t355 * t76;
t241 = -t132 * t123 + t269;
t359 = pkin(4) * t231;
t16 = t249 * t202 + t206 * t39;
t358 = t16 * t231;
t270 = -t203 * t93 - t92;
t289 = t123 * t338;
t43 = t270 - t289;
t339 = pkin(3) * t164;
t74 = -t123 * t193 - t339 + t359;
t357 = -t200 * t74 + (t297 * pkin(3) + t43) * t198;
t303 = t343 * t152 - t145;
t100 = t159 + t303;
t192 = t343 * pkin(2) + pkin(3);
t306 = t203 * t204;
t264 = -t152 * t204 - t148;
t99 = t264 - t337;
t353 = -t192 * t275 - (-t204 * t297 + (t343 * t342 - t306) * qJD(3)) * pkin(2) + t342 * t100 + t203 * t99;
t352 = -t204 * t205 + t343 * t207;
t350 = -t132 * t231 + t272;
t50 = t231 * t184 - t267;
t151 = (pkin(2) + t340) * t201 - t260;
t167 = t234 * t199;
t101 = t201 * pkin(3) - t167 * pkin(10) + t343 * t151 - t204 * t160;
t166 = t352 * t199;
t235 = -t204 * t151 - t343 * t160;
t105 = t166 * pkin(10) - t235;
t128 = t203 * t166 + t342 * t167;
t48 = t201 * pkin(4) + t342 * t101 - t203 * t105 - t128 * t338;
t139 = -pkin(3) * t166 + t173;
t230 = t342 * t166 - t203 * t167;
t75 = -pkin(4) * t230 - t128 * t193 + t139;
t247 = t198 * t75 + t200 * t48;
t232 = -t203 * t101 - t342 * t105;
t245 = t198 * t201 + t200 * t230;
t47 = t245 * pkin(11) - t232;
t348 = t247 * t202 + t206 * t47;
t346 = t199 * t351;
t38 = (t106 * t231 + t329) * t198;
t278 = t198 * t295;
t3 = -t16 * qJD(5) - t13 * t202 + t14 * t307 + t33 * t311;
t344 = t3 * t200 + t28 * t278;
t335 = t342 * t93 - t90;
t333 = pkin(3) * qJD(4);
t331 = t198 * t66;
t129 = t352 * t346;
t71 = t128 * qJD(4) + t203 * t129 + t209;
t330 = t198 * t71;
t324 = t114 - t353;
t282 = t342 * t204;
t134 = -t192 * t297 + (-t204 * t275 + (-t343 * t203 - t282) * qJD(3)) * pkin(2);
t268 = -t100 * t203 + t342 * t99;
t45 = t268 - t289;
t323 = t134 - t45;
t322 = t106 * t202;
t317 = t128 * t202;
t316 = t164 * t367;
t315 = t171 * t164;
t191 = t342 * pkin(3) + pkin(4);
t314 = t191 * t198;
t195 = t199 ^ 2;
t208 = qJD(1) ^ 2;
t313 = t195 * t208;
t302 = t205 ^ 2 - t207 ^ 2;
t296 = qJD(5) * t106;
t293 = qJD(2) - t189;
t288 = pkin(11) * t312;
t287 = t195 * t341;
t285 = t207 * t313;
t274 = qJD(1) * qJD(2) * t195;
t263 = t189 + t300;
t186 = pkin(2) * t279;
t185 = pkin(2) * t281;
t256 = t207 * t274;
t79 = t200 * t304 + t365;
t253 = -t79 + t278;
t80 = t231 * t307 + t365;
t252 = -t80 + t278;
t251 = -t15 * t355 - t28 * t80 + t344;
t248 = t198 * t74 + t200 * t43;
t181 = pkin(3) * t203 + t193;
t240 = -t181 * t202 + t191 * t307;
t239 = t181 * t206 + t191 * t308;
t188 = qJD(2) * t290;
t154 = t188 - t246;
t225 = t151 * t276 + t343 * t154 + t204 * t155 - t160 * t298;
t60 = -pkin(10) * t213 + t225;
t216 = t235 * qJD(3) - t204 * t154 + t343 * t155;
t61 = -t129 * pkin(10) + t216;
t229 = t101 * t275 - t105 * t297 + t203 * t61 + t342 * t60;
t228 = -pkin(8) * t310 - t291;
t227 = -pkin(8) * t255 + t180;
t226 = -t171 * t367 - t259;
t223 = t228 * t189;
t221 = -t202 * t47 + t247 * t206;
t84 = t128 * t206 + t245 * t202;
t218 = t232 * qJD(4) - t203 * t60 + t342 * t61;
t215 = -t94 * t275 - t88 * t297 + t350;
t115 = pkin(3) * t213 + t186;
t169 = -pkin(2) * t306 + t342 * t192 + pkin(4);
t165 = pkin(2) * t282 + t203 * t192 + t193;
t140 = t185 - t339;
t113 = t198 * t230 - t200 * t201;
t107 = t164 ^ 2 - t367 ^ 2;
t103 = -t164 * t184 - t210;
t102 = -t184 * t367 + t125;
t85 = -pkin(11) * t320 + t359;
t83 = -t201 * t311 - t230 * t307 + t317;
t70 = t230 * qJD(4) + t342 * t129 - t203 * t213;
t69 = t185 + t74;
t44 = -t114 + t335;
t42 = -pkin(11) * t309 + t242;
t37 = t71 * pkin(4) - t70 * t193 + t115;
t36 = -t198 * t48 + t200 * t75;
t35 = -t198 * t42 + t200 * t85;
t34 = -t198 * t45 + t200 * t69;
t30 = t84 * qJD(5) + t202 * t70 + t71 * t307;
t29 = -t71 * t308 + t206 * t70 + (t245 * t206 - t317) * qJD(5);
t18 = -t70 * t338 + t218;
t17 = -t71 * t338 + t229;
t9 = -t18 * t198 + t200 * t37;
t7 = t252 * t106 + t362;
t4 = t252 * t349 + t366;
t1 = [0, 0, 0, 0.2e1 * t205 * t256, -0.2e1 * t302 * t274, t263 * t207 * t299, -t263 * t279, 0, (t223 + (t228 * t201 - 0.2e1 * t287) * qJD(1)) * qJD(2), -0.2e1 * pkin(1) * t256 - (-pkin(8) * t279 + t188) * t189 - t227 * t201, t125 * t167 - t129 * t164, t125 * t166 + t129 * t367 + (t164 * t214 - t167 * t211) * t199, t125 * t201 + t129 * t184, (-t214 * t184 - t201 * t211) * t199, 0, 0.2e1 * t171 * t234 * t346 - t166 * t179 + t216 * t184 - t186 * t367 + t217 * t201, -t225 * t184 - t259 * t201 + t173 * t125 + t171 * t129 + (qJD(1) * t167 - t164) * t186, t128 * t65 + t231 * t70, t123 * t70 - t128 * t66 + t230 * t65 - t231 * t71, t65 * t201 + t70 * t238, -t66 * t201 - t71 * t238, 0, -t108 * t230 - t115 * t123 + t132 * t71 + t139 * t66 + t219 * t201 + t218 * t238, t108 * t128 + t115 * t231 + t132 * t70 + t139 * t65 + t269 * t201 - t229 * t238, t24 * t84 - t29 * t349, -t24 * t83 - t25 * t84 - t29 * t76 + t30 * t349, -t106 * t29 - t113 * t24 + (-t349 * t71 + t66 * t84) * t198, t106 * t30 + t113 * t25 + (-t66 * t83 - t71 * t76) * t198, (-t106 * t71 - t113 * t66) * t198, -(-t17 * t202 + t18 * t307 + t37 * t311) * t106 + t221 * t331 - t3 * t113 + t15 * t330 + t9 * t76 + t36 * t25 + t6 * t83 + t28 * t30 + t348 * t296, (t17 * t206 + t18 * t308 + t37 * t312) * t106 - t348 * t331 + t2 * t113 - t16 * t330 - t9 * t349 + t36 * t24 + t6 * t84 + t28 * t29 + t221 * t296; 0, 0, 0, -t205 * t285, t302 * t313, t293 * t280, -t293 * t281, 0, t208 * t287 + (t228 * qJD(2) - t223) * qJD(1), pkin(1) * t285 + (-pkin(8) * t281 + t187) * t189 - t227, t316, t107, t102, t103, 0, -t264 * t184 + t367 * t185 + t315 + (-t148 + (-pkin(2) * t184 - t135) * t204) * qJD(3) + t265, t303 * t184 + (t164 * t281 - t184 * t276) * pkin(2) + t226, -t319, t51, t49, t50, 0, t140 * t123 + t215 + (t134 - t268) * t238, -t140 * t231 + t353 * t238 + t241, t10, t4, t8, t7, t38, -t34 * t76 + (t165 * t294 + t324 * t202 + (t169 * t295 - t323 * t206) * t200) * t106 + (-t165 * t327 - t134 * t76 - t169 * t25 + (t106 * t69 + t169 * t329 - t6) * t206) * t198 + t251, t34 * t349 + (-t165 * t295 + t324 * t206 + (t169 * t294 + t323 * t202) * t200) * t106 + (-(t165 * t206 + t169 * t308) * t66 + t134 * t349 - t169 * t24 - t69 * t322 + t358) * t198 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t316, t107, t102, t103, 0, -t236 * t184 + t217 + t315, t266 * t184 + t226, -t319, t51, t49, t50, 0, -t270 * t238 + (-t123 * t164 - t238 * t297) * pkin(3) + t215, t335 * t238 + (t164 * t231 - t238 * t275) * pkin(3) + t241, t10, t4, t8, t7, t38, t240 * t331 - t25 * t314 - t6 * t311 + t251 + t357 * t76 + (t239 * qJD(5) - (-t342 * t202 - t203 * t307) * t333 - t202 * t44 + t248 * t206) * t106, t16 * t355 - t239 * t331 - t24 * t314 + t237 - t357 * t349 + (t240 * qJD(5) + (-t203 * t308 + t206 * t342) * t333 - t248 * t202 - t206 * t44) * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t319, t51, t49, t50, 0, -t184 * t242 + t350, t238 * t271 + t241, t10, t253 * t349 + t366, t8, t106 * t253 + t362, t38, -t28 * t79 - t35 * t76 + (-t202 * t41 + (pkin(4) * t295 + t206 * t42) * t200) * t106 + (-t66 * t288 - pkin(4) * t25 - t15 * t231 + (pkin(4) * t329 - t6 + (pkin(11) * qJD(5) + t85) * t106) * t206) * t198 + t344, -(t206 * t41 + t308 * t42) * t106 + t35 * t349 + (-pkin(11) * t66 * t311 - t85 * t322 + t358 + (-t24 - t62) * pkin(4)) * t198 + (pkin(4) * t307 - t288) * t296 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t349 * t76, t349 ^ 2 - t76 ^ 2, -t106 * t76 + t24, t106 * t349 - t25, t331, -t106 * t16 + t28 * t349 + t3, -t106 * t15 + t28 * t76 - t2;];
tauc_reg = t1;
