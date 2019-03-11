% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% tau_reg [6x31]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRP8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:24:37
% EndTime: 2019-03-09 06:24:47
% DurationCPUTime: 4.78s
% Computational Cost: add. (6991->476), mult. (13342->585), div. (0->0), fcn. (8989->10), ass. (0->251)
t186 = cos(qJ(3));
t329 = cos(qJ(4));
t252 = t329 * qJD(4);
t212 = t329 * qJD(3) + t252;
t183 = sin(qJ(3));
t328 = sin(qJ(4));
t256 = t328 * t183;
t238 = qJD(1) * t256;
t258 = t329 * t183;
t269 = qJD(3) + qJD(4);
t57 = qJDD(1) * t258 + t186 * (t212 * qJD(1) + t328 * qJDD(1)) - t269 * t238;
t56 = qJDD(5) + t57;
t118 = t328 * t186 + t258;
t112 = t118 * qJD(1);
t345 = qJD(5) + t112;
t351 = t345 ^ 2;
t257 = t329 * t186;
t111 = -qJD(1) * t257 + t238;
t182 = sin(qJ(5));
t185 = cos(qJ(5));
t245 = t185 * t269;
t87 = -t111 * t182 - t245;
t350 = t345 * t87;
t217 = t185 * t111 - t182 * t269;
t349 = t217 * t345;
t275 = qJD(5) * t185;
t293 = t112 * t185;
t348 = t275 + t293;
t276 = qJD(5) * t182;
t347 = t112 * t182 + t276;
t187 = cos(qJ(1));
t175 = g(2) * t187;
t184 = sin(qJ(1));
t326 = g(1) * t184;
t250 = -t175 + t326;
t273 = qJD(1) * qJD(3);
t251 = t186 * t273;
t271 = t183 * qJDD(1);
t346 = t251 + t271;
t117 = t256 - t257;
t243 = qJD(5) * t118 + qJD(1);
t305 = t185 * t56;
t341 = t112 * t269;
t193 = -t117 * qJDD(1) - t341;
t268 = qJDD(3) + qJDD(4);
t36 = -qJD(5) * t245 - t111 * t276 - t182 * t268 - t185 * t193;
t254 = qJD(4) * t328;
t338 = -qJD(3) * t328 - t254;
t83 = -t183 * t212 + t186 * t338;
t84 = t183 * t338 + t186 * t212;
t344 = t345 * (t182 * t243 - t185 * t84) - t117 * t36 - t118 * t305 + t217 * t83;
t189 = -pkin(1) - pkin(7);
t142 = t189 * qJD(1) + qJD(2);
t280 = qJD(1) * t183;
t101 = -pkin(8) * t280 + t142 * t183;
t140 = t189 * qJDD(1) + qJDD(2);
t119 = t186 * t140;
t270 = t186 * qJDD(1);
t278 = qJD(3) * t183;
t73 = -t142 * t278 + qJDD(3) * pkin(3) + t119 + (t183 * t273 - t270) * pkin(8);
t277 = qJD(3) * t186;
t76 = -t346 * pkin(8) + t183 * t140 + t142 * t277;
t279 = qJD(1) * t186;
t102 = -pkin(8) * t279 + t186 * t142;
t98 = qJD(3) * pkin(3) + t102;
t202 = -t101 * t254 + t98 * t252 + t328 * t73 + t329 * t76;
t16 = pkin(9) * t268 + t202;
t176 = qJDD(1) * qJ(2);
t177 = qJD(1) * qJD(2);
t99 = t346 * pkin(3) + t176 + t177;
t24 = t57 * pkin(4) - pkin(9) * t193 + t99;
t97 = t329 * t101;
t66 = t328 * t98 + t97;
t59 = pkin(9) * t269 + t66;
t129 = pkin(3) * t280 + qJD(1) * qJ(2);
t68 = pkin(4) * t112 + pkin(9) * t111 + t129;
t221 = t185 * t16 + t182 * t24 + t68 * t275 - t276 * t59;
t315 = qJ(6) * t56;
t2 = qJD(6) * t345 + t221 + t315;
t248 = t182 * t16 - t185 * t24 + t59 * t275 + t68 * t276;
t331 = pkin(5) * t56;
t4 = qJDD(6) + t248 - t331;
t343 = t4 * t182 + t2 * t185;
t342 = (t276 * t345 - t305) * pkin(9);
t70 = t328 * t102 + t97;
t224 = pkin(3) * t254 - t70;
t340 = t347 * pkin(5) - t348 * qJ(6) - t182 * qJD(6);
t181 = qJ(3) + qJ(4);
t170 = cos(t181);
t156 = t170 * pkin(9);
t173 = t183 * pkin(3);
t339 = t156 - t173;
t191 = qJD(1) ^ 2;
t214 = -qJ(2) * t191 - t250;
t235 = g(1) * t187 + g(2) * t184;
t265 = 0.2e1 * t177;
t337 = 0.2e1 * t176 + t265 - t235;
t37 = -qJD(5) * t217 + t182 * t193 - t185 * t268;
t334 = -t118 * t268 - t269 * t84;
t333 = t217 ^ 2;
t330 = pkin(9) * t56;
t327 = pkin(5) * t111;
t169 = sin(t181);
t157 = g(3) * t169;
t158 = g(3) * t170;
t325 = g(3) * t182;
t324 = g(3) * t185;
t247 = t101 * t252 + t98 * t254 + t328 * t76 - t329 * t73;
t17 = -pkin(4) * t268 + t247;
t5 = t37 * pkin(5) + t36 * qJ(6) + qJD(6) * t217 + t17;
t323 = t5 * t182;
t322 = t217 * t87;
t321 = pkin(8) - t189;
t320 = -t340 - t224;
t96 = t328 * t101;
t65 = t329 * t98 - t96;
t79 = -pkin(4) * t111 + pkin(9) * t112;
t319 = t182 * t79 + t185 * t65;
t71 = t329 * t102 - t96;
t74 = pkin(3) * t279 + t79;
t318 = t182 * t74 + t185 * t71;
t155 = qJ(2) + t173;
t80 = pkin(4) * t118 + pkin(9) * t117 + t155;
t126 = t321 * t183;
t127 = t321 * t186;
t86 = -t329 * t126 - t328 * t127;
t317 = t182 * t80 + t185 * t86;
t316 = pkin(9) * qJD(5);
t35 = t182 * t68 + t185 * t59;
t314 = t345 * t35;
t58 = -pkin(4) * t269 - t65;
t30 = t87 * pkin(5) + qJ(6) * t217 + t58;
t313 = t112 * t30;
t312 = t112 * t58;
t161 = t328 * pkin(3) + pkin(9);
t311 = t161 * t56;
t310 = t182 * t56;
t309 = t182 * t83;
t308 = t182 * t87;
t307 = t182 * t217;
t306 = t185 * t37;
t304 = t185 * t83;
t302 = t185 * t87;
t301 = t185 * t217;
t300 = t36 * t182;
t299 = -t66 + t340;
t298 = -t117 * t268 + t83 * t269;
t297 = pkin(1) * qJDD(1);
t296 = t345 * t111;
t295 = t111 * t112;
t292 = t169 * t184;
t291 = t170 * t184;
t290 = t170 * t187;
t289 = t182 * t184;
t288 = t182 * t187;
t287 = t184 * t185;
t286 = t187 * t185;
t34 = -t182 * t59 + t185 * t68;
t285 = qJD(6) - t34;
t138 = g(2) * t290;
t284 = t185 * t138 + t169 * t324;
t283 = t187 * pkin(1) + t184 * qJ(2);
t180 = t186 ^ 2;
t282 = t183 ^ 2 - t180;
t190 = qJD(3) ^ 2;
t281 = -t190 - t191;
t143 = pkin(3) * t277 + qJD(2);
t272 = qJDD(3) * t183;
t267 = g(1) * t291;
t137 = t169 * t175;
t266 = t329 * pkin(3);
t264 = t170 * t289;
t27 = t30 * t276;
t263 = t30 * t275;
t52 = t58 * t276;
t262 = g(1) * t264 - t182 * t138 - t169 * t325;
t261 = g(1) * t292 - t137 + t158;
t260 = t182 * t329;
t259 = t185 * t329;
t246 = t185 * t345;
t244 = qJDD(2) - t297;
t242 = pkin(3) * t252;
t240 = -t5 - t267;
t239 = -t17 - t267;
t104 = t169 * t289 - t286;
t106 = t169 * t288 + t287;
t237 = g(1) * t106 + g(2) * t104;
t105 = t169 * t287 + t288;
t107 = t169 * t286 - t289;
t236 = -g(1) * t107 - g(2) * t105;
t233 = -t185 * pkin(5) - t182 * qJ(6);
t232 = pkin(5) * t182 - qJ(6) * t185;
t230 = -t311 + t312;
t25 = -pkin(5) * t345 + t285;
t26 = qJ(6) * t345 + t35;
t229 = -t182 * t26 + t185 * t25;
t228 = t182 * t25 + t185 * t26;
t227 = -t301 + t308;
t226 = -t25 * t111 + t27 + t284;
t225 = t34 * t111 + t284 + t52;
t223 = pkin(4) - t233;
t43 = pkin(4) * t84 - pkin(9) * t83 + t143;
t115 = t321 * t278;
t116 = qJD(3) * t127;
t85 = -t328 * t126 + t329 * t127;
t47 = -t85 * qJD(4) + t328 * t115 - t329 * t116;
t220 = t182 * t43 + t185 * t47 + t80 * t275 - t276 * t86;
t219 = 0.2e1 * qJ(2) * t273 + qJDD(3) * t189;
t218 = pkin(4) * t169 - t339;
t215 = t223 * t170;
t213 = -t35 * t111 + t17 * t182 + t58 * t275 + t262;
t210 = t26 * t111 - t30 * t293 - t262 - t323;
t209 = g(1) * t104 - g(2) * t106 + t170 * t325 - t248;
t207 = -t161 * t276 + t185 * t242;
t206 = t129 * t111 + t138 + t157 - t247 - t267;
t205 = qJD(5) * t229 + t343;
t204 = qJD(5) * t227 - t300 - t306;
t203 = -t217 * t30 + qJDD(6) - t209;
t201 = -t189 * t190 + t337;
t200 = t348 * t25 - t347 * t26 - t261 + t343;
t199 = -g(1) * t105 + g(2) * t107 - t170 * t324 + t221;
t198 = -g(1) * (t170 * pkin(5) * t287 + pkin(4) * t291 + pkin(9) * t292 + qJ(6) * t264) + t223 * t157;
t197 = -t118 * t310 + t117 * t37 - t83 * t87 + (-t182 * t84 - t185 * t243) * t345;
t48 = t86 * qJD(4) - t329 * t115 - t328 * t116;
t196 = t129 * t112 - t202 + t261;
t188 = -pkin(8) - pkin(7);
t172 = t187 * qJ(2);
t168 = qJDD(3) * t186;
t162 = -t266 - pkin(4);
t114 = -t266 - t223;
t103 = t111 * qJ(6);
t61 = t111 ^ 2 - t112 ^ 2;
t51 = -pkin(5) * t217 + qJ(6) * t87;
t49 = -t117 * t232 + t85;
t45 = -t111 * t269 - t57;
t44 = t193 + t341;
t40 = -pkin(5) * t118 + t182 * t86 - t185 * t80;
t39 = qJ(6) * t118 + t317;
t32 = t182 * t65 - t185 * t79 + t327;
t31 = -t103 + t319;
t29 = t182 * t71 - t185 * t74 + t327;
t28 = -t103 + t318;
t19 = -t36 + t350;
t12 = -t111 * t217 + t246 * t345 + t310;
t11 = -t87 * t111 - t182 * t351 + t305;
t10 = -t217 * t246 - t300;
t9 = t232 * t83 + (qJD(5) * t233 + qJD(6) * t185) * t117 + t48;
t8 = -t84 * pkin(5) + t317 * qJD(5) + t182 * t47 - t185 * t43;
t7 = qJ(6) * t84 + qJD(6) * t118 + t220;
t6 = (-t36 - t350) * t185 + (-t37 + t349) * t182;
t1 = [qJDD(1), t250, t235, qJDD(2) - t250 - 0.2e1 * t297, t337, -t244 * pkin(1) - g(1) * (-pkin(1) * t184 + t172) - g(2) * t283 + (t265 + t176) * qJ(2), qJDD(1) * t180 - 0.2e1 * t183 * t251, -0.2e1 * t183 * t270 + 0.2e1 * t273 * t282, -t183 * t190 + t168, -t186 * t190 - t272, 0, t183 * t201 + t186 * t219, -t183 * t219 + t186 * t201, -t111 * t83 - t117 * t193, t111 * t84 - t83 * t112 + t117 * t57 - t118 * t193, t298, t334, 0, t143 * t112 + t99 * t118 + t129 * t84 + t155 * t57 - t169 * t235 - t268 * t85 - t269 * t48, -g(1) * t290 - g(2) * t291 - t143 * t111 - t99 * t117 + t129 * t83 + t155 * t193 - t268 * t86 - t269 * t47, -t83 * t301 + (t36 * t185 - t217 * t276) * t117 (-t302 + t307) * t83 + (-t300 + t306 + (-t301 - t308) * qJD(5)) * t117, -t117 * t305 - t36 * t118 - t217 * t84 + (t117 * t276 + t304) * t345, t117 * t310 - t37 * t118 - t87 * t84 + (t117 * t275 - t309) * t345, t118 * t56 + t345 * t84, -t248 * t118 + t34 * t84 + t48 * t87 + t85 * t37 + ((-qJD(5) * t86 + t43) * t345 + t80 * t56 - t58 * qJD(5) * t117) * t185 + ((-qJD(5) * t80 - t47) * t345 - t86 * t56 - t17 * t117 + t58 * t83) * t182 + t236, -t220 * t345 - t317 * t56 - t221 * t118 - t35 * t84 - t48 * t217 - t85 * t36 + t58 * t304 + (-t17 * t185 + t52) * t117 + t237, t30 * t309 - t345 * t8 - t118 * t4 - t25 * t84 + t37 * t49 - t40 * t56 + t87 * t9 + (-t263 - t323) * t117 + t236, -t36 * t40 - t37 * t39 - t7 * t87 - t8 * t217 + t229 * t83 + t235 * t170 + (qJD(5) * t228 + t182 * t2 - t185 * t4) * t117, -t30 * t304 + t345 * t7 + t118 * t2 + t26 * t84 + t36 * t49 + t39 * t56 + t217 * t9 + (t5 * t185 - t27) * t117 - t237, t2 * t39 + t26 * t7 + t5 * t49 + t30 * t9 + t4 * t40 + t25 * t8 - g(1) * (pkin(5) * t107 + qJ(6) * t106 + t172) - g(2) * (pkin(5) * t105 + qJ(6) * t104 + t283) + (-g(1) * t218 + g(2) * t188) * t187 + (-g(1) * (-pkin(1) + t188) - g(2) * t218) * t184; 0, 0, 0, qJDD(1), -t191, t244 + t214, 0, 0, 0, 0, 0, t183 * t281 + t168, t186 * t281 - t272, 0, 0, 0, 0, 0, -qJD(1) * t112 + t298, qJD(1) * t111 + t334, 0, 0, 0, 0, 0, t197, t344, t197 (-t302 - t307) * t84 + t227 * qJD(1) + t204 * t118, -t344, t229 * qJD(1) + t117 * t5 + t118 * t205 + t228 * t84 - t30 * t83 - t250; 0, 0, 0, 0, 0, 0, t186 * t191 * t183, -t282 * t191, t270, -t271, qJDD(3), g(3) * t183 + t186 * t214 + t119, g(3) * t186 + (-t140 - t214) * t183, -t295, t61, t44, t45, t268, t70 * t269 + (-t112 * t279 - t269 * t254 + t329 * t268) * pkin(3) + t206, t71 * t269 + (t111 * t279 - t269 * t252 - t328 * t268) * pkin(3) + t196, t10, t6, t12, t11, t296, t162 * t37 + t224 * t87 + t239 * t185 + t230 * t182 + ((-qJD(5) * t161 - t74) * t185 + (-t242 + t71) * t182) * t345 + t225, -t162 * t36 - t224 * t217 + t230 * t185 + (-t207 + t318) * t345 + t213, t114 * t37 - t320 * t87 + t240 * t185 + (-t311 + t313) * t182 + (-t161 * t275 - t182 * t242 + t29) * t345 + t226, t28 * t87 + t29 * t217 + (-t217 * t260 - t259 * t87) * qJD(4) * pkin(3) + t204 * t161 + t200, t114 * t36 - t320 * t217 + (-qJD(5) * t30 + t311) * t185 + (-t28 + t207) * t345 + t210, t5 * t114 - t26 * t28 - t25 * t29 - g(3) * t339 - t320 * t30 + (-t186 * t326 + (t25 * t260 + t26 * t259) * qJD(4)) * pkin(3) + t205 * t161 - (-pkin(3) * t186 - pkin(9) * t169 - t215) * t175 + t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t295, t61, t44, t45, t268, t269 * t66 + t206, t269 * t65 + t196, t10, t6, t12, t11, t296, -pkin(4) * t37 - t66 * t87 + (t345 * t65 + t312 - t330) * t182 + ((-t79 - t316) * t345 + t239) * t185 + t225, pkin(4) * t36 + t217 * t66 + t58 * t293 + t319 * t345 + t213 + t342, t345 * t32 - t223 * t37 + t299 * t87 + (t313 - t330) * t182 + (-t316 * t345 + t240) * t185 + t226, pkin(9) * t204 + t217 * t32 + t31 * t87 + t200, t217 * t299 - t223 * t36 - t31 * t345 + t210 - t263 - t342, -t5 * t223 - t26 * t31 - t25 * t32 - g(3) * t156 + t299 * t30 + t215 * t175 + (t205 + t137) * pkin(9) + t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t322, -t87 ^ 2 + t333, t19, -t37 - t349, t56, t217 * t58 + t209 + t314, t34 * t345 + t58 * t87 - t199, -t51 * t87 - t203 + t314 + 0.2e1 * t331, pkin(5) * t36 - t37 * qJ(6) - (t26 - t35) * t217 + (t25 - t285) * t87, 0.2e1 * t315 - t30 * t87 - t51 * t217 + (0.2e1 * qJD(6) - t34) * t345 + t199, t2 * qJ(6) - t4 * pkin(5) - t30 * t51 - t25 * t35 - g(1) * (-pkin(5) * t104 + qJ(6) * t105) - g(2) * (pkin(5) * t106 - qJ(6) * t107) + t285 * t26 + t232 * t158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t322 - t56, t19, -t333 - t351, -t26 * t345 + t203 - t331;];
tau_reg  = t1;
