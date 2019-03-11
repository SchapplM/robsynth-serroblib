% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:26:12
% EndTime: 2019-03-09 03:26:20
% DurationCPUTime: 4.97s
% Computational Cost: add. (7190->535), mult. (14250->625), div. (0->0), fcn. (9676->10), ass. (0->259)
t184 = sin(qJ(3));
t187 = cos(qJ(3));
t316 = sin(pkin(9));
t317 = cos(pkin(9));
t120 = t184 * t317 + t187 * t316;
t111 = t120 * qJD(1);
t360 = qJD(5) + t111;
t252 = t360 ^ 2;
t183 = sin(qJ(5));
t186 = cos(qJ(5));
t189 = -pkin(1) - pkin(7);
t140 = qJD(1) * t189 + qJD(2);
t251 = -qJ(4) * qJD(1) + t140;
t102 = t251 * t184;
t262 = t317 * t102;
t284 = qJD(1) * t187;
t103 = -qJ(4) * t284 + t140 * t187;
t99 = qJD(3) * pkin(3) + t103;
t61 = t316 * t99 + t262;
t54 = qJD(3) * pkin(8) + t61;
t260 = t316 * t184;
t246 = qJD(1) * t260;
t261 = t317 * t187;
t114 = qJD(1) * t261 - t246;
t171 = t184 * pkin(3);
t155 = qJ(2) + t171;
t125 = qJD(1) * t155 + qJD(4);
t62 = pkin(4) * t111 - pkin(8) * t114 + t125;
t25 = t183 * t62 + t186 * t54;
t16 = qJ(6) * t360 + t25;
t365 = t16 * t360;
t24 = -t183 * t54 + t186 * t62;
t364 = t24 * t360;
t333 = t360 * t25;
t278 = t186 * qJD(3);
t88 = t114 * t183 - t278;
t256 = t360 * t88;
t258 = qJD(3) * t316;
t238 = qJD(1) * t258;
t254 = qJDD(1) * t317;
t259 = qJD(3) * t317;
t363 = qJD(1) * t259 + qJDD(1) * t316;
t209 = (t238 - t254) * t187 + t363 * t184;
t281 = qJD(5) * t183;
t44 = -qJD(5) * t278 - qJDD(3) * t183 + t114 * t281 + t186 * t209;
t90 = qJD(3) * t183 + t114 * t186;
t86 = t90 * t281;
t356 = -t186 * t44 - t86;
t175 = qJ(3) + pkin(9);
t164 = cos(t175);
t163 = sin(t175);
t341 = g(3) * t163;
t188 = cos(qJ(1));
t173 = g(2) * t188;
t185 = sin(qJ(1));
t353 = g(1) * t185 - t173;
t206 = -t164 * t353 + t341;
t153 = pkin(3) * t316 + pkin(8);
t265 = qJD(5) * t360 * t153;
t362 = t265 - t206;
t191 = qJD(1) ^ 2;
t361 = -qJ(2) * t191 - t353;
t177 = qJDD(1) * qJ(2);
t116 = t120 * qJD(3);
t121 = t261 - t260;
t248 = qJD(5) * t120 + qJD(1);
t306 = t120 * t186;
t113 = t184 * t258 - t187 * t259;
t311 = t113 * t186;
t268 = t184 * t254 + t187 * t363;
t82 = t184 * t238 - t268;
t80 = -qJDD(5) + t82;
t359 = t360 * (t183 * t248 + t311) + t116 * t90 + t121 * t44 + t306 * t80;
t321 = t186 * t88;
t324 = t183 * t90;
t228 = t321 + t324;
t320 = t186 * t90;
t45 = qJD(5) * t90 - qJDD(3) * t186 - t183 * t209;
t322 = t186 * t45;
t325 = t183 * t44;
t358 = t121 * (qJD(5) * (t183 * t88 - t320) - t322 + t325) + t228 * t116;
t357 = -t360 * t90 + t45;
t280 = qJD(5) * t186;
t336 = -t183 * t45 - t280 * t88;
t355 = (-t321 + t324) * t111 + t336 - t356;
t139 = qJDD(1) * t189 + qJDD(2);
t180 = t184 ^ 2;
t181 = t187 ^ 2;
t287 = t180 + t181;
t255 = t287 * t139;
t242 = g(1) * t188 + g(2) * t185;
t178 = qJD(1) * qJD(2);
t271 = 0.2e1 * t178;
t352 = 0.2e1 * t177 + t271 - t242;
t96 = t316 * t102;
t60 = t317 * t99 - t96;
t53 = -qJD(3) * pkin(4) - t60;
t26 = pkin(5) * t88 - qJ(6) * t90 + t53;
t328 = t153 * t80;
t351 = t26 * t360 + t328;
t308 = t116 * t183;
t75 = t183 * t80;
t350 = t360 * (t121 * t280 - t308) - t113 * t88 + t120 * t45 - t121 * t75;
t349 = -t114 * t116 - t121 * t209;
t81 = pkin(4) * t120 - pkin(8) * t121 + t155;
t297 = qJ(4) - t189;
t126 = t297 * t184;
t257 = t297 * t187;
t85 = -t126 * t317 - t257 * t316;
t335 = t183 * t81 + t186 * t85;
t100 = -qJD(3) * t257 - qJD(4) * t184;
t283 = qJD(3) * t184;
t210 = -qJD(4) * t187 + t283 * t297;
t65 = t100 * t317 + t210 * t316;
t282 = qJD(3) * t187;
t144 = pkin(3) * t282 + qJD(2);
t71 = -pkin(4) * t113 + pkin(8) * t116 + t144;
t11 = -qJD(5) * t335 - t183 * t65 + t186 * t71;
t346 = t90 ^ 2;
t345 = t114 ^ 2;
t343 = pkin(5) * t80;
t342 = pkin(3) * t187;
t340 = g(3) * t164;
t339 = g(3) * t184;
t338 = g(3) * t186;
t337 = t90 * t88;
t122 = t187 * t139;
t272 = t187 * qJDD(1);
t275 = qJD(1) * qJD(4);
t276 = qJD(1) * qJD(3);
t58 = -t187 * t275 - t140 * t283 + qJDD(3) * pkin(3) + t122 + (t184 * t276 - t272) * qJ(4);
t66 = t251 * t282 + (-qJ(4) * qJDD(1) + t139 - t275) * t184;
t27 = -t316 * t66 + t317 * t58;
t28 = t316 * t58 + t317 * t66;
t70 = t103 * t317 - t96;
t72 = pkin(3) * t284 + pkin(4) * t114 + pkin(8) * t111;
t33 = t183 * t72 + t186 * t70;
t334 = qJ(6) * t80;
t332 = t113 * t90;
t331 = t114 * t88;
t330 = t114 * t90;
t329 = t153 * t44;
t327 = t153 * t88;
t326 = t153 * t90;
t76 = t186 * t80;
t236 = pkin(5) * t183 - qJ(6) * t186;
t69 = t103 * t316 + t262;
t318 = -qJD(6) * t183 + t236 * t360 - t69;
t315 = pkin(1) * qJDD(1);
t313 = t360 * t111;
t312 = t360 * t114;
t310 = t114 * t111;
t307 = t116 * t186;
t305 = t163 * t185;
t304 = t163 * t188;
t303 = t164 * t185;
t302 = t164 * t188;
t301 = t183 * t185;
t300 = t183 * t188;
t299 = t185 * t186;
t298 = t188 * t186;
t296 = qJD(6) - t24;
t295 = -qJD(3) * t116 + qJDD(3) * t121;
t294 = g(2) * t164 * t298 + t163 * t338;
t293 = g(1) * t302 + g(2) * t303;
t292 = pkin(8) * t164 - t171;
t291 = (t271 + t177) * qJ(2);
t290 = pkin(1) * t188 + qJ(2) * t185;
t288 = t180 - t181;
t190 = qJD(3) ^ 2;
t286 = -t190 - t191;
t285 = qJD(1) * t125;
t279 = t111 * qJD(3);
t274 = qJDD(3) * t184;
t273 = t184 * qJDD(1);
t270 = t88 ^ 2 - t346;
t269 = t187 * t191 * t184;
t267 = pkin(4) * t303 + pkin(8) * t305 + t185 * t342;
t266 = t317 * pkin(3);
t264 = t187 * t276;
t170 = t188 * qJ(2);
t263 = -pkin(1) * t185 + t170;
t23 = qJDD(3) * pkin(8) + t28;
t95 = qJDD(4) + t177 + t178 + (t264 + t273) * pkin(3);
t31 = -pkin(4) * t82 + pkin(8) * t209 + t95;
t4 = -t183 * t23 + t186 * t31 - t280 * t54 - t281 * t62;
t250 = t287 * qJDD(1);
t249 = qJDD(2) - t315;
t247 = t184 * t264;
t22 = -qJDD(3) * pkin(4) - t27;
t245 = -pkin(8) * t163 - t342;
t105 = t163 * t301 - t298;
t107 = t163 * t300 + t299;
t244 = -g(1) * t107 - g(2) * t105;
t106 = t163 * t299 + t300;
t108 = t163 * t298 - t301;
t243 = -g(1) * t108 - g(2) * t106;
t64 = t100 * t316 - t210 * t317;
t237 = pkin(5) * t186 + qJ(6) * t183;
t15 = -pkin(5) * t360 + t296;
t235 = t15 * t186 - t16 * t183;
t234 = -t15 * t183 - t16 * t186;
t233 = t183 * t25 + t186 * t24;
t232 = t183 * t24 - t186 * t25;
t32 = -t183 * t70 + t186 * t72;
t42 = -t183 * t85 + t186 * t81;
t84 = -t126 * t316 + t257 * t317;
t226 = -t111 * t113 - t120 * t82;
t225 = t186 * t313 + t280 * t360 - t75;
t182 = -qJ(4) - pkin(7);
t224 = t171 * t188 + t182 * t185 + t263;
t223 = -t183 * t313 - t281 * t360 - t76;
t222 = qJD(3) * t113 - qJDD(3) * t120;
t221 = -pkin(4) - t237;
t220 = t171 * t185 - t182 * t188 + t290;
t3 = t183 * t31 + t186 * t23 + t280 * t62 - t281 * t54;
t10 = t183 * t71 + t186 * t65 + t280 * t81 - t281 * t85;
t218 = t360 * t53 + t328;
t217 = -g(1) * t305 + g(2) * t304 - t340;
t216 = 0.2e1 * qJ(2) * t276 + qJDD(3) * t189;
t215 = -g(1) * t303 - t265;
t213 = -t153 * t322 + t217;
t212 = t183 * t256 - t322;
t211 = -t223 - t331;
t205 = g(1) * t105 - g(2) * t107 + t183 * t340 + t4;
t204 = pkin(4) * t304 - pkin(8) * t302 + t224;
t203 = pkin(4) * t305 - pkin(8) * t303 + t220;
t200 = -t121 * t336 - t308 * t88;
t1 = qJD(6) * t360 + t3 - t334;
t2 = qJDD(6) - t4 + t343;
t199 = qJD(5) * t235 + t1 * t186 + t2 * t183;
t198 = -qJD(5) * t233 - t4 * t183 + t3 * t186;
t197 = t113 * t61 + t116 * t60 - t120 * t28 - t121 * t27 + t353;
t196 = t26 * t90 + qJDD(6) - t205;
t195 = -t189 * t190 + t352;
t194 = -g(1) * t106 + g(2) * t108 - t164 * t338 + t3;
t193 = t120 * t75 - t121 * t45 + t116 * t88 + (t113 * t183 - t186 * t248) * t360;
t192 = -t45 * t306 + t88 * t311 + t248 * t320 + (qJD(1) * t88 - t332 + (qJD(5) * t88 - t44) * t120) * t183;
t167 = qJDD(3) * t187;
t154 = -t266 - pkin(4);
t118 = -t266 + t221;
t110 = t111 ^ 2;
t48 = pkin(5) * t90 + qJ(6) * t88;
t46 = t121 * t236 + t84;
t37 = -t113 * t360 - t120 * t80;
t35 = -pkin(5) * t120 - t42;
t34 = qJ(6) * t120 + t335;
t19 = -pkin(5) * t114 - t32;
t18 = qJ(6) * t114 + t33;
t17 = -t44 + t256;
t14 = t225 - t330;
t13 = -t236 * t116 + (qJD(5) * t237 - qJD(6) * t186) * t121 + t64;
t12 = t320 * t360 - t325;
t9 = t121 * t356 - t307 * t90;
t8 = pkin(5) * t113 - t11;
t7 = -qJ(6) * t113 + qJD(6) * t120 + t10;
t6 = -t121 * t76 - t332 - t120 * t44 + (-t121 * t281 - t307) * t360;
t5 = pkin(5) * t45 + qJ(6) * t44 - qJD(6) * t90 + t22;
t20 = [0, 0, 0, 0, 0, qJDD(1), t353, t242, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t353 - 0.2e1 * t315, t352, -pkin(1) * t249 - g(1) * t263 - g(2) * t290 + t291, qJDD(1) * t181 - 0.2e1 * t247, -0.2e1 * t184 * t272 + 0.2e1 * t276 * t288, -t184 * t190 + t167, qJDD(1) * t180 + 0.2e1 * t247, -t187 * t190 - t274, 0, t184 * t195 + t187 * t216, -t184 * t216 + t187 * t195, -t189 * t250 - t255 + t353, -g(1) * (t185 * t189 + t170) - g(2) * (pkin(7) * t188 + t290) + t189 * t255 + t291, t349, t111 * t116 + t113 * t114 + t120 * t209 + t121 * t82, t295, t226, t222, 0, -qJD(3) * t64 - qJDD(3) * t84 + t111 * t144 - t113 * t125 + t120 * t95 - t155 * t82 - t163 * t242, -qJD(3) * t65 - qJDD(3) * t85 + t114 * t144 - t116 * t125 + t121 * t95 - t155 * t209 - t293, -t111 * t65 + t114 * t64 - t209 * t84 + t82 * t85 + t197, -g(1) * t224 - g(2) * t220 + t125 * t144 + t155 * t95 - t27 * t84 + t28 * t85 - t60 * t64 + t61 * t65, t9, t358, t6, t200, -t350, t37, -t53 * t308 + t360 * t11 - t113 * t24 + t120 * t4 - t42 * t80 + t45 * t84 + t64 * t88 + (t183 * t22 + t280 * t53) * t121 + t243, -t53 * t307 - t10 * t360 + t113 * t25 - t120 * t3 + t335 * t80 - t44 * t84 + t64 * t90 + (t186 * t22 - t281 * t53) * t121 - t244, -t10 * t88 - t11 * t90 + t42 * t44 - t335 * t45 + t233 * t116 + (qJD(5) * t232 - t183 * t3 - t186 * t4) * t121 + t293, -g(1) * t204 - g(2) * t203 + t25 * t10 + t24 * t11 + t22 * t84 + t3 * t335 + t4 * t42 + t53 * t64, t9, t6, -t358, t37, t350, t200, -t26 * t308 - t360 * t8 + t113 * t15 - t120 * t2 + t13 * t88 + t35 * t80 + t45 * t46 + (t183 * t5 + t26 * t280) * t121 + t243, -t34 * t45 - t35 * t44 - t7 * t88 + t8 * t90 - t235 * t116 + (qJD(5) * t234 - t1 * t183 + t186 * t2) * t121 + t293, t26 * t307 + t1 * t120 + t360 * t7 - t113 * t16 - t13 * t90 - t34 * t80 + t44 * t46 + (-t186 * t5 + t26 * t281) * t121 + t244, t1 * t34 + t16 * t7 + t5 * t46 + t26 * t13 + t2 * t35 + t15 * t8 - g(1) * (pkin(5) * t108 + qJ(6) * t107 + t204) - g(2) * (pkin(5) * t106 + qJ(6) * t105 + t203); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t191, t361 + t249, 0, 0, 0, 0, 0, 0, t184 * t286 + t167, t187 * t286 - t274, -t250, t255 + t361, 0, 0, 0, 0, 0, 0, -qJD(1) * t111 + t295, -qJD(1) * t114 + t222, -t226 - t349, -t197 - t285, 0, 0, 0, 0, 0, 0, t193, t359, t192, -qJD(1) * t233 + t113 * t232 + t116 * t53 + t120 * t198 - t121 * t22 - t353, 0, 0, 0, 0, 0, 0, t193, t192, -t359, qJD(1) * t235 + t113 * t234 + t116 * t26 + t120 * t199 - t121 * t5 - t353; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t269, -t288 * t191, t272, -t269, -t273, qJDD(3), t187 * t361 + t122 + t339, g(3) * t187 + (-t139 - t361) * t184, 0, 0, t310, -t110 + t345, -t209 + t279, -t310 (t246 + t114) * qJD(3) - t268, qJDD(3), t69 * qJD(3) - t125 * t114 + (qJDD(3) * t317 - t111 * t284) * pkin(3) + t206 + t27, t70 * qJD(3) + t125 * t111 + (-qJDD(3) * t316 - t114 * t284) * pkin(3) - t217 - t28 (t61 - t69) * t114 + (t70 - t60) * t111 + (t209 * t317 + t316 * t82) * pkin(3), t60 * t69 - t61 * t70 + (t317 * t27 + t316 * t28 + t339 + (-t353 - t285) * t187) * pkin(3), t12, -t111 * t228 + t336 + t356, t14, t212, -t211, -t312, -t360 * t32 - t114 * t24 + t154 * t45 - t69 * t88 + (t215 - t22) * t186 + t218 * t183 + t294, t360 * t33 + t114 * t25 - t154 * t44 - t69 * t90 + t218 * t186 + (t22 + t362) * t183, t32 * t90 + t33 * t88 + (-t111 * t24 + t3 + (-t24 + t326) * qJD(5)) * t186 + (-t111 * t25 - t329 - t4 + (-t25 + t327) * qJD(5)) * t183 + t213, t22 * t154 - t25 * t33 - t24 * t32 - t53 * t69 - g(1) * t267 - g(3) * (-pkin(4) * t163 + t292) - (-pkin(4) * t164 + t245) * t173 + t198 * t153, t12, t14, t86 + (t111 * t90 + t45) * t183 + (t44 + t256) * t186, -t312, t211, t212, t360 * t19 + t114 * t15 + t118 * t45 + t318 * t88 + (t215 - t5) * t186 + t351 * t183 + t294, t18 * t88 - t19 * t90 + (t111 * t15 + t1 + (t15 + t326) * qJD(5)) * t186 + (-t111 * t16 - t329 + t2 + (-t16 + t327) * qJD(5)) * t183 + t213, -t360 * t18 - t114 * t16 + t118 * t44 - t318 * t90 - t351 * t186 + (-t5 - t362) * t183, t5 * t118 - t16 * t18 - t15 * t19 - g(1) * (t237 * t303 + t267) - g(3) * t292 + t318 * t26 - t221 * t341 + t199 * t153 - (t164 * t221 + t245) * t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t246 + t114) * qJD(3) + t268, -t209 - t279, -t110 - t345, t111 * t61 + t114 * t60 - t242 + t95, 0, 0, 0, 0, 0, 0, t223 - t331, -t186 * t252 - t330 + t75, t355, -t114 * t53 + (t4 + t333) * t186 + (t3 - t364) * t183 - t242, 0, 0, 0, 0, 0, 0, -t183 * t252 - t331 - t76, t355, t225 + t330, -t114 * t26 + (-t2 + t365) * t186 + (t15 * t360 + t1) * t183 - t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t337, -t270, t17, -t337, -t357, -t80, -t53 * t90 + t205 + t333, t53 * t88 - t194 + t364, 0, 0, t337, t17, t270, -t80, t357, -t337, -t48 * t88 - t196 + t333 - 0.2e1 * t343, pkin(5) * t44 - qJ(6) * t45 + (t16 - t25) * t90 + (t15 - t296) * t88, -0.2e1 * t334 - t26 * t88 + t48 * t90 + (0.2e1 * qJD(6) - t24) * t360 + t194, t1 * qJ(6) - t2 * pkin(5) - t26 * t48 - t15 * t25 - g(1) * (-pkin(5) * t105 + qJ(6) * t106) - g(2) * (pkin(5) * t107 - qJ(6) * t108) + t236 * t340 + t296 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80 + t337, t17, -t346 - t252, t196 + t343 - t365;];
tau_reg  = t20;
