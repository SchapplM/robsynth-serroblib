% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRP4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:12:47
% EndTime: 2021-01-16 03:13:08
% DurationCPUTime: 5.71s
% Computational Cost: add. (3396->493), mult. (7472->662), div. (0->0), fcn. (5385->10), ass. (0->275)
t168 = sin(pkin(6));
t173 = sin(qJ(2));
t176 = cos(qJ(2));
t275 = qJD(1) * qJD(2);
t253 = t176 * t275;
t170 = cos(pkin(6));
t286 = qJD(3) * t170;
t325 = qJDD(2) * pkin(8);
t389 = qJD(1) * t286 + t325 + (qJDD(1) * t173 + t253) * t168;
t174 = cos(qJ(5));
t171 = sin(qJ(5));
t285 = qJD(3) * t171;
t175 = cos(qJ(3));
t289 = qJD(2) * t175;
t112 = t174 * t289 + t285;
t172 = sin(qJ(3));
t291 = qJD(2) * t172;
t148 = qJD(5) + t291;
t323 = t112 * t148;
t274 = qJD(2) * qJD(3);
t251 = t172 * t274;
t270 = t175 * qJDD(2);
t384 = -t251 + t270;
t42 = qJD(5) * t112 - t174 * qJDD(3) + t384 * t171;
t388 = t42 - t323;
t387 = t42 + t323;
t294 = qJD(1) * t168;
t261 = t173 * t294;
t121 = qJD(2) * pkin(8) + t261;
t293 = qJD(1) * t170;
t297 = -t172 * t121 + t175 * t293;
t382 = qJD(4) - t297;
t177 = pkin(4) + pkin(8);
t282 = qJD(3) * t175;
t120 = t177 * t282;
t303 = t174 * t176;
t81 = (-t171 * t173 + t172 * t303) * t168;
t386 = -qJD(1) * t81 + t174 * t120;
t178 = -pkin(3) - pkin(9);
t328 = qJ(4) * t172;
t108 = t178 * t175 - pkin(2) - t328;
t129 = t177 * t172;
t279 = qJD(5) * t174;
t280 = qJD(5) * t171;
t284 = qJD(3) * t172;
t154 = pkin(3) * t284;
t327 = qJ(4) * t175;
t225 = pkin(9) * t172 - t327;
t281 = qJD(4) * t172;
t191 = t225 * qJD(3) - t281;
t78 = t154 + t191;
t306 = t173 * t174;
t307 = t172 * t176;
t82 = (t171 * t307 + t306) * t168;
t385 = qJD(1) * t82 + t108 * t280 - t171 * t120 - t129 * t279 - t174 * t78;
t68 = -qJD(3) * pkin(3) + t382;
t167 = sin(pkin(10));
t169 = cos(pkin(10));
t318 = t168 * t175;
t309 = t172 * t173;
t267 = t168 * t309;
t101 = -t170 * t175 + t267;
t345 = g(3) * t101;
t314 = t170 * t173;
t366 = t169 * t176;
t93 = t167 * t314 - t366;
t370 = t167 * t176;
t95 = t169 * t314 + t370;
t383 = g(1) * (t167 * t318 + t172 * t93) - g(2) * (t169 * t318 + t172 * t95) - t345;
t305 = t173 * t175;
t102 = t168 * t305 + t170 * t172;
t180 = qJD(2) ^ 2;
t217 = qJDD(2) * t176 - t173 * t180;
t249 = t176 * t274;
t317 = t168 * t176;
t259 = qJD(2) * t317;
t59 = -qJD(3) * t267 + (t259 + t286) * t175;
t381 = t59 * qJD(3) + t102 * qJDD(3) + (t172 * t217 + t175 * t249) * t168;
t60 = qJD(3) * t102 + t172 * t259;
t380 = -t60 * qJD(3) - t101 * qJDD(3) + (-t172 * t249 + t175 * t217) * t168;
t255 = t171 * t289;
t283 = qJD(3) * t174;
t114 = -t255 + t283;
t299 = pkin(4) * t291 + t382;
t44 = t178 * qJD(3) + t299;
t256 = t176 * t294;
t65 = qJD(2) * t108 - t256;
t18 = t171 * t44 + t174 * t65;
t273 = qJDD(1) * t170;
t248 = t175 * t273;
t359 = t121 * t282 + t389 * t172;
t193 = qJDD(4) - t248 + t359;
t250 = t175 * t274;
t271 = t172 * qJDD(2);
t197 = t250 + t271;
t13 = t197 * pkin(4) + t178 * qJDD(3) + t193;
t254 = t173 * t275;
t224 = -qJDD(1) * t317 + t168 * t254;
t213 = pkin(3) * t251 + t224;
t21 = qJD(2) * t191 + qJDD(2) * t108 + t213;
t246 = t174 * t13 - t171 * t21;
t188 = -t18 * qJD(5) + t246;
t332 = t42 * qJ(6);
t111 = qJDD(5) + t197;
t342 = t111 * pkin(5);
t1 = -t114 * qJD(6) + t188 + t332 + t342;
t245 = -t171 * t13 - t174 * t21 - t44 * t279 + t65 * t280;
t43 = -qJD(5) * t255 + t171 * qJDD(3) + (qJD(3) * qJD(5) + t384) * t174;
t338 = qJ(6) * t43;
t2 = -qJD(6) * t112 - t245 - t338;
t10 = -qJ(6) * t112 + t18;
t336 = t10 * t148;
t17 = -t171 * t65 + t174 * t44;
t9 = -qJ(6) * t114 + t17;
t8 = pkin(5) * t148 + t9;
t379 = -(t148 * t8 - t2) * t171 + (t1 + t336) * t174;
t378 = t8 - t9;
t242 = qJ(6) * t175 - t108;
t349 = pkin(5) * t282 + t242 * t279 + (-qJ(6) * t284 - qJD(5) * t129 + qJD(6) * t175 - t78) * t171 + t386;
t278 = qJD(5) * t175;
t257 = t171 * t278;
t276 = t174 * qJD(6);
t377 = -t175 * t276 + (t172 * t283 + t257) * qJ(6) - t385;
t311 = t171 * t172;
t363 = qJ(6) - t178;
t76 = t175 * t121 + t172 * t293;
t64 = pkin(4) * t289 + t76;
t48 = t174 * t64;
t155 = pkin(3) * t291;
t89 = t225 * qJD(2) + t155;
t339 = t363 * t280 - t276 + t171 * t89 - t48 - (pkin(5) * t175 - qJ(6) * t311) * qJD(2);
t124 = t363 * t174;
t340 = t171 * t64 + t174 * t89;
t375 = -qJ(6) * t174 * t291 - qJD(5) * t124 - t171 * qJD(6) - t340;
t347 = pkin(5) * t174;
t262 = -pkin(4) - t347;
t331 = pkin(5) * t279 - t262 * t291 + t382;
t164 = qJD(3) * qJ(4);
t69 = -t164 - t76;
t322 = t114 * t148;
t374 = -t43 + t322;
t373 = t43 + t322;
t205 = -qJ(4) * t282 - t281;
t372 = t154 + t205 - t261;
t92 = t174 * t111;
t371 = -t148 * t280 + t92;
t151 = pkin(5) * t171 + qJ(4);
t221 = t151 * t172 + t175 * t363;
t369 = t168 * t221;
t228 = pkin(3) * t175 + t328;
t368 = t168 * t228;
t367 = t169 * t173;
t247 = t177 + t347;
t365 = t173 * t247;
t364 = t176 * (pkin(2) + t221);
t162 = qJDD(3) * qJ(4);
t163 = qJD(3) * qJD(4);
t239 = t121 * t284 - t172 * t273 - t389 * t175;
t22 = -t162 - t163 + t239;
t324 = qJDD(3) * pkin(3);
t23 = t193 - t324;
t362 = t23 * t172 - t22 * t175;
t264 = t170 * t305;
t319 = t168 * t172;
t103 = t264 - t319;
t302 = t175 * t176;
t344 = g(3) * t102;
t361 = -t344 - g(2) * (t103 * t169 + t167 * t302);
t360 = -t76 * qJD(3) + t359 + t383;
t46 = t164 + t64;
t357 = t178 * t111 + t148 * t46;
t313 = t170 * t176;
t210 = -t167 * t173 + t169 * t313;
t96 = t167 * t313 + t367;
t231 = g(1) * t96 - g(2) * t210;
t179 = qJD(3) ^ 2;
t346 = pkin(8) * t179;
t356 = 0.2e1 * qJDD(2) * pkin(2) + t168 * (-g(3) * t176 + t254) - t224 + t231 - t346;
t194 = -g(3) * t317 + t231;
t125 = pkin(2) + t228;
t272 = qJDD(2) * t125;
t27 = t205 * qJD(2) + t213 - t272;
t353 = -qJD(2) * t372 + t194 - t27 + t272 - t346;
t352 = t114 ^ 2;
t343 = g(3) * t168;
t337 = qJD(2) * pkin(2);
t335 = t174 * t42;
t330 = t174 * t108 + t171 * t129;
t329 = pkin(8) * qJDD(3);
t321 = pkin(2) * t367;
t320 = t167 * t170;
t316 = t169 * t170;
t315 = t170 * t171;
t312 = t171 * t111;
t310 = t171 * t176;
t308 = t172 * t174;
t304 = t174 * t175;
t298 = qJDD(1) - g(3);
t130 = t177 * t175;
t165 = t172 ^ 2;
t166 = t175 ^ 2;
t296 = t165 - t166;
t295 = t165 + t166;
t292 = qJD(2) * t125;
t290 = qJD(2) * t173;
t288 = qJD(3) * t112;
t287 = qJD(3) * t114;
t277 = qJD(5) * t178;
t244 = -pkin(5) * t112 - qJD(6);
t32 = -t244 + t46;
t269 = t32 * t280;
t268 = t32 * t279;
t266 = t168 * t304;
t265 = t170 * t308;
t263 = t172 * t180 * t175;
t260 = t168 * t290;
t100 = t170 * t309 + t318;
t241 = t100 * t169 + t167 * t307;
t240 = -t100 * t167 + t169 * t307;
t238 = pkin(8) - t262;
t237 = t112 * t256;
t236 = t114 * t256;
t235 = t175 * t256;
t232 = -g(1) * t93 + g(2) * t95;
t227 = pkin(3) * t172 - t327;
t226 = pkin(8) * t173 + t125 * t176;
t222 = t151 * t175 - t172 * t363;
t220 = t148 * t171;
t62 = -t101 * t171 + t168 * t303;
t61 = t101 * t174 + t168 * t310;
t209 = -t148 * t279 - t312;
t208 = t173 * t227;
t207 = t176 * t227;
t206 = t228 * t170;
t79 = t210 * t172;
t80 = t96 * t172;
t204 = -g(1) * (-t171 * t80 - t174 * t93) - g(2) * (t171 * t79 + t174 * t95) - g(3) * t82;
t203 = -g(1) * (t171 * t93 - t174 * t80) - g(2) * (-t171 * t95 + t174 * t79) - g(3) * t81;
t202 = -g(1) * t240 - g(2) * t241 - t345;
t137 = t169 * t302;
t201 = -g(1) * (-t103 * t167 + t137) + t361;
t199 = t173 * t222;
t198 = t176 * t222;
t14 = t384 * pkin(4) - t22;
t5 = t43 * pkin(5) + qJDD(6) + t14;
t196 = t201 + t5;
t195 = t14 + t201;
t190 = -g(1) * (-t171 * t240 - t174 * t96) - g(2) * (-t171 * t241 + t174 * t210) - g(3) * t62 + t245;
t122 = -t256 - t337;
t187 = -t329 + (t122 + t256 - t337) * qJD(3);
t77 = -t256 - t292;
t186 = t329 + (-t256 - t77 + t292) * qJD(3);
t185 = -t221 * t173 + t247 * t176;
t127 = t167 * t319;
t184 = g(1) * (t175 * t93 - t127) - g(2) * (-t169 * t319 + t175 * t95) - qJD(3) * t297 - t239 - t344;
t181 = -g(1) * (-t171 * t96 + t174 * t240) - g(2) * (t171 * t210 + t174 * t241) - g(3) * t61 + t188;
t157 = t167 * pkin(2);
t145 = pkin(2) * t316;
t123 = t363 * t171;
t119 = t177 * t284;
t118 = -qJ(4) * t289 + t155;
t117 = t174 * t129;
t110 = t112 ^ 2;
t99 = pkin(5) * t304 + t130;
t67 = t77 * t291;
t66 = -pkin(5) * t257 - t238 * t284;
t40 = -qJ(6) * t304 + t330;
t37 = t172 * pkin(5) + t242 * t171 + t117;
t26 = -t148 ^ 2 * t174 - t287 - t312;
t25 = -t148 * t220 - t288 + t92;
t16 = qJD(5) * t61 + t60 * t171 + t174 * t260;
t15 = qJD(5) * t62 - t171 * t260 + t60 * t174;
t4 = -t102 * t42 + t111 * t62 + t114 * t59 - t148 * t16;
t3 = t102 * t43 + t111 * t61 + t112 * t59 + t148 * t15;
t6 = [t298, 0, t217 * t168, (-qJDD(2) * t173 - t176 * t180) * t168, 0, 0, 0, 0, 0, t380, -t381, (t101 * t172 + t102 * t175) * qJDD(2) + (t172 * t60 + t175 * t59 + (t101 * t175 - t102 * t172) * qJD(3)) * qJD(2), -t380, t381, t23 * t101 - t22 * t102 - t69 * t59 + t68 * t60 - g(3) + (-t176 * t27 + t77 * t290) * t168, 0, 0, 0, 0, 0, t3, t4, t3, t4, -t112 * t16 - t114 * t15 + t42 * t61 + t43 * t62, t1 * t61 + t10 * t16 + t102 * t5 + t15 * t8 - t2 * t62 + t32 * t59 - g(3); 0, qJDD(2), t298 * t317 + t231, -t298 * t173 * t168 + t232, qJDD(2) * t165 + 0.2e1 * t172 * t250, 0.2e1 * t172 * t270 - 0.2e1 * t296 * t274, qJDD(3) * t172 + t175 * t179, qJDD(3) * t175 - t172 * t179, 0, t187 * t172 + t356 * t175, -t356 * t172 + t187 * t175, (t69 * t172 + t68 * t175) * qJD(3) + t295 * t325 + (-g(3) * t173 - t295 * t253) * t168 - t232 + t362, t186 * t172 - t353 * t175, t353 * t172 + t186 * t175, -t27 * t125 + t69 * t235 - t68 * t172 * t256 - g(1) * (-t226 * t320 - t228 * t367 - t321) - g(2) * (-(t167 * t228 + t157) * t173 + (t169 * t206 + t145) * t176) - t226 * t343 + t372 * t77 + (t69 * t284 + t68 * t282 - g(1) * t366 - g(2) * (t316 * t173 + t370) + t362) * pkin(8), t42 * t171 * t175 + (t171 * t284 - t174 * t278) * t114, (-t112 * t171 + t114 * t174) * t284 + (t171 * t43 + t335 + (t112 * t174 + t114 * t171) * qJD(5)) * t175, (t148 * t285 - t42) * t172 + (t209 + t287) * t175, (t148 * t283 - t43) * t172 + (-t288 - t371) * t175, t111 * t172 + t148 * t282, (-t108 * t171 + t117) * t111 - t119 * t112 + t130 * t43 + (-t283 * t46 + t246) * t172 + (-t171 * t78 + t386) * t148 + (-t330 * t148 - t18 * t172) * qJD(5) + (qJD(3) * t17 + t14 * t174 - t280 * t46 - t237) * t175 + t204, -t330 * t111 - t119 * t114 - t130 * t42 + (t285 * t46 + t245) * t172 + t385 * t148 + (-qJD(3) * t18 - t14 * t171 - t279 * t46 - t236) * t175 + t203, t111 * t37 + t112 * t66 + t43 * t99 + (-t283 * t32 + t1) * t172 + t349 * t148 + (qJD(3) * t8 + t174 * t5 - t237 - t269) * t175 + t204, -t111 * t40 + t114 * t66 - t42 * t99 + (t285 * t32 - t2) * t172 - t377 * t148 + (-qJD(3) * t10 - t171 * t5 - t236 - t268) * t175 + t203, t37 * t42 - t40 * t43 - t349 * t114 - t377 * t112 + (t10 * t174 - t171 * t8) * t284 + (t1 * t171 - t174 * t2 + (t10 * t171 + t174 * t8) * qJD(5) + t194) * t175, t2 * t40 + t1 * t37 + t5 * t99 - g(1) * (-t321 + t185 * t169 + (-t364 - t365) * t320) - g(2) * (t145 * t176 - t157 * t173 + (t176 * t221 + t365) * t316 + t185 * t167) - (t173 * t238 + t364) * t343 + t349 * t8 + (t66 - t235) * t32 + t377 * t10; 0, 0, 0, 0, -t263, t296 * t180, t271, t270, qJDD(3), -t122 * t291 + t248 - t360, -t122 * t289 - t184, -t227 * qJDD(2), -0.2e1 * t324 + qJDD(4) + t67 + (-qJD(2) * t118 - t273) * t175 + t360, 0.2e1 * t162 + 0.2e1 * t163 + (t118 * t172 + t175 * t77) * qJD(2) + t184, -t22 * qJ(4) - t23 * pkin(3) - t77 * t118 - t68 * t76 - g(1) * (-t169 * t207 + (t227 * t314 + t368) * t167) - g(2) * (-t167 * t207 + (-t170 * t208 - t368) * t169) - g(3) * (-t168 * t208 + t206) - t382 * t69, -t114 * t220 - t335, t387 * t171 - t373 * t174, (-t114 * t175 - t148 * t311) * qJD(2) + t371, (t112 * t175 - t148 * t308) * qJD(2) + t209, -t148 * t289, -t17 * t289 + qJ(4) * t43 - t48 * t148 + t299 * t112 + t357 * t174 + ((t89 - t277) * t148 + t195) * t171, -qJ(4) * t42 + t340 * t148 + t18 * t289 + t299 * t114 - t357 * t171 + (-t148 * t277 + t195) * t174, t268 - t111 * t124 + t151 * t43 + t339 * t148 + t331 * t112 + (-t175 * t8 + t308 * t32) * qJD(2) + t196 * t171, -t269 + t111 * t123 - t151 * t42 - t375 * t148 + t331 * t114 + (t10 * t175 - t311 * t32) * qJD(2) + t196 * t174, -t112 * t375 - t339 * t114 + t123 * t43 - t124 * t42 - t379 - t383, -t2 * t123 - t1 * t124 + t5 * t151 - g(1) * (t169 * t198 + (-t222 * t314 + t369) * t167) - g(2) * (t167 * t198 + (t170 * t199 - t369) * t169) - g(3) * (t168 * t199 + t170 * t221) + t339 * t8 + t331 * t32 + t375 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t271, qJDD(3) + t263, -t165 * t180 - t179, qJD(3) * t69 + t202 + t23 + t67, 0, 0, 0, 0, 0, t25, t26, t25, t26, t374 * t171 + t388 * t174, -qJD(3) * t32 + t202 + t379; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114 * t112, -t110 + t352, -t388, t374, t111, -t46 * t114 + t18 * t148 + t181, t112 * t46 + t148 * t17 + t190, 0.2e1 * t342 + t332 + t336 + (t244 - t32) * t114 + t181, -pkin(5) * t352 + t338 + t148 * t9 + (qJD(6) + t32) * t112 + t190, pkin(5) * t42 - t378 * t112, t378 * t10 + (t1 - t32 * t114 - g(1) * ((-t167 * t315 + t169 * t308) * t176 + (-t167 * t265 - t169 * t171) * t173 - t167 * t266) - g(2) * ((t167 * t308 + t169 * t315) * t176 + (-t167 * t171 + t169 * t265) * t173 + t169 * t266) - g(3) * (-t170 * t304 + (t172 * t306 + t310) * t168)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t373, -t387, -t110 - t352, -g(1) * (-t167 * t264 + t127 + t137) + t10 * t112 + t8 * t114 + t5 + t361;];
tau_reg = t6;
