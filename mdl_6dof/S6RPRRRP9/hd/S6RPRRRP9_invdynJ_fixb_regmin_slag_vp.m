% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRP9
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
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRP9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:28:34
% EndTime: 2019-03-09 06:28:45
% DurationCPUTime: 4.37s
% Computational Cost: add. (5075->464), mult. (10100->597), div. (0->0), fcn. (6728->10), ass. (0->235)
t184 = sin(qJ(4));
t185 = sin(qJ(3));
t274 = qJD(1) * t185;
t245 = t184 * t274;
t320 = pkin(8) + pkin(9);
t246 = qJD(4) * t320;
t189 = cos(qJ(3));
t222 = pkin(3) * t189 + pkin(8) * t185;
t132 = t222 * qJD(1);
t191 = -pkin(1) - pkin(7);
t155 = t191 * qJD(1) + qJD(2);
t188 = cos(qJ(4));
t282 = t188 * t189;
t279 = t184 * t132 + t155 * t282;
t335 = pkin(9) * t245 + t184 * t246 + t279;
t111 = t188 * t132;
t286 = t185 * t188;
t251 = pkin(9) * t286;
t289 = t184 * t189;
t340 = t188 * t246 - t155 * t289 + t111 + (pkin(4) * t189 + t251) * qJD(1);
t273 = qJD(1) * t189;
t239 = t184 * t273;
t269 = qJD(3) * t188;
t123 = -t239 + t269;
t271 = qJD(3) * t184;
t124 = t188 * t273 + t271;
t183 = sin(qJ(5));
t187 = cos(qJ(5));
t217 = t123 * t183 + t187 * t124;
t321 = t217 ^ 2;
t68 = -t187 * t123 + t124 * t183;
t65 = t68 ^ 2;
t339 = -t65 + t321;
t162 = qJD(4) + t274;
t314 = g(3) * t185;
t186 = sin(qJ(1));
t190 = cos(qJ(1));
t326 = -g(1) * t186 + g(2) * t190;
t203 = -t189 * t326 - t314;
t270 = qJD(3) * t185;
t297 = qJDD(3) * pkin(3);
t228 = t155 * t270 - t297;
t151 = t191 * qJDD(1) + qJDD(2);
t280 = t189 * t151;
t79 = t228 - t280;
t338 = qJD(4) * pkin(8) * t162 + t203 + t79;
t253 = qJD(4) + qJD(5);
t262 = qJD(5) * t187;
t265 = qJD(4) * t188;
t283 = t187 * t188;
t291 = t183 * t184;
t305 = t183 * t245 - t187 * t265 - t188 * t262 + t253 * t291 - t274 * t283;
t127 = t183 * t188 + t184 * t187;
t101 = t127 * qJD(1);
t76 = t253 * t127;
t304 = t185 * t101 + t76;
t337 = qJ(6) * t68;
t336 = t217 * t68;
t264 = qJD(4) * t189;
t240 = t188 * t264;
t244 = t184 * t270;
t334 = t240 - t244;
t254 = t189 * qJDD(1);
t333 = qJD(3) * qJD(4) + t254;
t153 = qJD(5) + t162;
t247 = qJD(1) * t240 + t333 * t184;
t213 = t188 * qJDD(3) - t247;
t198 = qJD(1) * t244 + t213;
t263 = qJD(5) * t183;
t243 = t185 * t269;
t206 = -t184 * t264 - t243;
t61 = qJD(1) * t206 + t184 * qJDD(3) + t333 * t188;
t19 = -t123 * t262 + t124 * t263 - t183 * t198 - t187 * t61;
t332 = t153 * t68 - t19;
t182 = qJ(4) + qJ(5);
t173 = cos(t182);
t135 = t185 * t155;
t112 = qJD(3) * pkin(8) + t135;
t138 = pkin(3) * t185 - pkin(8) * t189 + qJ(2);
t102 = t138 * qJD(1);
t121 = qJD(3) * t222 + qJD(2);
t66 = qJD(1) * t121 + qJDD(1) * t138;
t268 = qJD(3) * t189;
t80 = qJDD(3) * pkin(8) + t151 * t185 + t155 * t268;
t252 = -t102 * t265 - t184 * t66 - t188 * t80;
t267 = qJD(4) * t184;
t209 = -t112 * t267 - t252;
t13 = pkin(9) * t198 + t209;
t57 = t188 * t102 - t112 * t184;
t44 = -pkin(9) * t124 + t57;
t39 = pkin(4) * t162 + t44;
t58 = t102 * t184 + t112 * t188;
t45 = pkin(9) * t123 + t58;
t259 = qJD(1) * qJD(3);
t236 = t189 * t259;
t255 = t185 * qJDD(1);
t120 = qJDD(4) + t236 + t255;
t63 = t188 * t66;
t9 = pkin(4) * t120 - pkin(9) * t61 - qJD(4) * t58 - t184 * t80 + t63;
t235 = -t187 * t13 - t183 * t9 - t39 * t262 + t45 * t263;
t313 = g(3) * t189;
t293 = t155 * t189;
t303 = qJD(3) * pkin(3);
t113 = -t293 - t303;
t78 = -pkin(4) * t123 + t113;
t172 = sin(t182);
t287 = t185 * t186;
t87 = t172 * t190 + t173 * t287;
t285 = t185 * t190;
t89 = -t172 * t186 + t173 * t285;
t331 = g(1) * t87 - g(2) * t89 + t173 * t313 + t68 * t78 + t235;
t32 = pkin(5) * t68 + qJD(6) + t78;
t330 = t217 * t32;
t329 = t340 * t187;
t327 = qJ(6) * t217;
t126 = -t283 + t291;
t93 = t126 * t185;
t91 = t127 * t185;
t146 = t320 * t184;
t147 = t320 * t188;
t278 = -t183 * t146 + t187 * t147;
t325 = t146 * t262 + t147 * t263 + t340 * t183 + t335 * t187;
t86 = -t172 * t287 + t173 * t190;
t88 = t172 * t285 + t173 * t186;
t324 = -g(1) * t86 - g(2) * t88 + t172 * t313;
t43 = t187 * t45;
t17 = t183 * t39 + t43;
t238 = -t183 * t13 + t187 * t9;
t201 = -qJD(5) * t17 + t238;
t323 = -t78 * t217 + t201 + t324;
t20 = qJD(5) * t217 + t183 * t61 - t187 * t198;
t322 = t153 * t217 - t20;
t41 = t183 * t45;
t16 = t187 * t39 - t41;
t10 = t16 - t327;
t8 = pkin(5) * t153 + t10;
t317 = t10 - t8;
t316 = pkin(8) * t120;
t176 = t188 * pkin(4);
t312 = pkin(3) + t176;
t136 = pkin(4) * t184 + pkin(5) * t172;
t311 = pkin(7) + t136;
t310 = -t304 * qJ(6) - qJD(6) * t126 - t325;
t309 = -pkin(5) * t273 + t305 * qJ(6) - t278 * qJD(5) - qJD(6) * t127 + t335 * t183 - t329;
t308 = t187 * t44 - t41;
t119 = t188 * t138;
t288 = t184 * t191;
t234 = pkin(4) - t288;
t64 = -pkin(9) * t282 + t185 * t234 + t119;
t152 = t191 * t286;
t277 = t184 * t138 + t152;
t77 = -pkin(9) * t289 + t277;
t306 = t183 * t64 + t187 * t77;
t302 = t61 * t184;
t301 = -t126 * qJD(1) + t127 * t268 - t253 * t93;
t94 = t126 * t189;
t300 = qJD(3) * t94 + t253 * t91 + t101;
t299 = pkin(1) * qJDD(1);
t193 = qJD(1) ^ 2;
t298 = qJ(2) * t193;
t296 = t120 * t188;
t295 = t123 * t162;
t294 = t124 * t188;
t292 = t162 * t184;
t290 = t184 * t120;
t284 = t186 * t188;
t281 = t188 * t190;
t137 = pkin(5) * t173 + t176;
t181 = t189 ^ 2;
t276 = t185 ^ 2 - t181;
t192 = qJD(3) ^ 2;
t275 = -t192 - t193;
t272 = qJD(3) * t124;
t266 = qJD(4) * t185;
t261 = t113 * qJD(4);
t260 = t123 * qJD(3);
t257 = qJDD(1) * qJ(2);
t256 = qJDD(3) * t185;
t249 = 0.2e1 * qJD(1) * qJD(2);
t242 = t188 * t268;
t248 = t184 * t121 + t138 * t265 + t191 * t242;
t237 = g(2) * (t190 * pkin(1) + t186 * qJ(2));
t97 = t188 * t121;
t28 = t97 + (-t152 + (pkin(9) * t189 - t138) * t184) * qJD(4) + (t189 * t234 + t251) * qJD(3);
t30 = -t334 * pkin(9) - t266 * t288 + t248;
t232 = -t183 * t30 + t187 * t28;
t231 = -t183 * t44 - t43;
t229 = -t183 * t77 + t187 * t64;
t227 = -qJD(4) * t102 - t80;
t226 = t162 * t191 + t112;
t225 = -t187 * t146 - t147 * t183;
t122 = pkin(4) * t289 - t189 * t191;
t224 = qJD(1) + t266;
t223 = -t135 + (t245 + t267) * pkin(4);
t221 = g(1) * t190 + g(2) * t186;
t219 = qJDD(2) + t326;
t131 = pkin(3) + t137;
t179 = -qJ(6) - t320;
t215 = t131 * t185 + t179 * t189;
t214 = -t151 - t326;
t212 = t183 * t28 + t187 * t30 + t64 * t262 - t263 * t77;
t211 = t162 * t265 + t290;
t210 = -t162 * t267 + t296;
t81 = t334 * pkin(4) + t191 * t270;
t208 = 0.2e1 * qJ(2) * t259 + qJDD(3) * t191;
t205 = t214 + t298;
t200 = -t221 + t249 + 0.2e1 * t257;
t197 = -t191 * t192 + t200;
t195 = -pkin(4) * t198 + t228;
t194 = t20 * pkin(5) + qJDD(6) + t195;
t175 = t190 * qJ(2);
t171 = qJDD(3) * t189;
t168 = pkin(4) * t187 + pkin(5);
t115 = qJDD(5) + t120;
t106 = -t184 * t186 + t185 * t281;
t105 = t184 * t285 + t284;
t104 = t184 * t190 + t185 * t284;
t103 = -t184 * t287 + t281;
t92 = t127 * t189;
t51 = -qJ(6) * t126 + t278;
t50 = -qJ(6) * t127 + t225;
t38 = -t263 * t289 + (t253 * t282 - t244) * t187 + t206 * t183;
t36 = -t183 * t244 + t187 * t243 + t189 * t76;
t31 = t195 - t280;
t25 = -qJ(6) * t92 + t306;
t24 = pkin(5) * t185 + qJ(6) * t94 + t229;
t15 = t308 - t327;
t14 = t231 + t337;
t11 = t17 - t337;
t5 = t194 - t280;
t4 = -qJ(6) * t38 - qJD(6) * t92 + t212;
t3 = pkin(5) * t268 + qJ(6) * t36 - qJD(5) * t306 + qJD(6) * t94 + t232;
t2 = -qJ(6) * t20 - qJD(6) * t68 - t235;
t1 = pkin(5) * t115 + qJ(6) * t19 - qJD(6) * t217 + t201;
t6 = [qJDD(1), -t326, t221, t219 - 0.2e1 * t299, t200 -(qJDD(2) - t299) * pkin(1) - g(1) * (-pkin(1) * t186 + t175) - t237 + (t249 + t257) * qJ(2), qJDD(1) * t181 - 0.2e1 * t185 * t236, -0.2e1 * t185 * t254 + 0.2e1 * t259 * t276, -t185 * t192 + t171, -t189 * t192 - t256, 0, t185 * t197 + t189 * t208, -t185 * t208 + t189 * t197, t124 * t206 + t282 * t61 (-t123 * t188 + t124 * t184) * t270 + (t188 * t198 - t302 + (-t123 * t184 - t294) * qJD(4)) * t189 (-t162 * t269 + t61) * t185 + (t210 + t272) * t189 ((t162 + t274) * t271 + t213) * t185 + (-t211 + t260) * t189, t120 * t185 + t162 * t268 (-t138 * t267 + t97) * t162 + t119 * t120 - g(1) * t106 - g(2) * t104 + (-t191 * t260 + t63 - t226 * t265 + (-qJD(3) * t113 - t120 * t191 + t227) * t184) * t185 + (t191 * t213 + t79 * t184 + t188 * t261 + (t57 + (-t162 + t274) * t288) * qJD(3)) * t189, -t248 * t162 - t277 * t120 + g(1) * t105 - g(2) * t103 + (t226 * t267 + (-t113 * t188 + t124 * t191) * qJD(3) + t252) * t185 + (-qJD(3) * t58 - t184 * t261 + t79 * t188 - t191 * t61) * t189, t19 * t94 - t217 * t36, t19 * t92 + t20 * t94 - t217 * t38 + t36 * t68, -t115 * t94 - t153 * t36 - t185 * t19 + t217 * t268, -t115 * t92 - t153 * t38 - t185 * t20 - t268 * t68, t115 * t185 + t153 * t268, t232 * t153 + t229 * t115 + t238 * t185 + t16 * t268 + t81 * t68 + t122 * t20 + t31 * t92 + t78 * t38 - g(1) * t89 - g(2) * t87 + (-t153 * t306 - t17 * t185) * qJD(5), g(1) * t88 - g(2) * t86 - t115 * t306 - t122 * t19 - t153 * t212 - t17 * t268 + t185 * t235 + t217 * t81 - t31 * t94 - t78 * t36, t1 * t94 - t11 * t38 + t189 * t221 + t19 * t24 - t2 * t92 - t20 * t25 - t217 * t3 + t36 * t8 - t4 * t68, t2 * t25 + t11 * t4 + t1 * t24 + t8 * t3 + t5 * (pkin(5) * t92 + t122) + t32 * (pkin(5) * t38 + t81) - g(1) * t175 - t237 + (-g(1) * t215 - g(2) * t311) * t190 + (-g(1) * (-pkin(1) - t311) - g(2) * t215) * t186; 0, 0, 0, qJDD(1), -t193, t219 - t298 - t299, 0, 0, 0, 0, 0, t185 * t275 + t171, t189 * t275 - t256, 0, 0, 0, 0, 0, t189 * t213 + (-t290 + (-t123 + t239) * qJD(3)) * t185 + (-t184 * t268 - t188 * t224) * t162, -t189 * t61 + (t272 - t296) * t185 + (t184 * t224 - t242) * t162, 0, 0, 0, 0, 0, -t115 * t91 - t153 * t301 - t189 * t20 + t270 * t68, t115 * t93 + t153 * t300 + t189 * t19 + t217 * t270, -t19 * t91 + t20 * t93 + t217 * t301 + t300 * t68, -t1 * t91 - t11 * t300 - t189 * t5 - t2 * t93 + t270 * t32 - t301 * t8 + t326; 0, 0, 0, 0, 0, 0, t189 * t193 * t185, -t276 * t193, t254, -t255, qJDD(3), -t189 * t205 + t314, t185 * t205 + t313, t162 * t294 + t302 (t61 + t295) * t188 + (-t124 * qJD(4) + (-t124 + t271) * t274 + t213) * t184 (-t124 * t189 + t162 * t286) * qJD(1) + t211 (-t123 * t189 - t185 * t292) * qJD(1) + t210, -t162 * t273, -pkin(3) * t247 - t111 * t162 - t57 * t273 + t123 * t135 + (t162 * t293 - t316 + t261 + (t113 + t303) * t274) * t184 + (t297 - t338) * t188, -pkin(3) * t61 + t279 * t162 + t58 * t273 - t124 * t135 + (t113 * t162 - t316) * t188 + t338 * t184, -t127 * t19 - t217 * t305, t126 * t19 - t127 * t20 - t217 * t304 + t305 * t68, t115 * t127 - t153 * t305 - t217 * t273, -t115 * t126 - t153 * t304 + t273 * t68, -t153 * t273, t225 * t115 - t312 * t20 + t31 * t126 - t16 * t273 + t304 * t78 + t223 * t68 + (-t147 * t262 + (qJD(5) * t146 + t335) * t183 - t329) * t153 - t203 * t173, -t278 * t115 + t31 * t127 + t153 * t325 + t17 * t273 + t203 * t172 + t19 * t312 + t223 * t217 - t305 * t78, -t1 * t127 - t11 * t304 - t126 * t2 + t185 * t326 + t19 * t50 - t20 * t51 - t217 * t309 + t305 * t8 - t310 * t68 - t313, t2 * t51 + t1 * t50 + t5 * (pkin(5) * t126 - t312) + g(3) * t215 + t309 * t8 + (pkin(4) * t292 + pkin(5) * t304 - t135) * t32 + t310 * t11 + t326 * (t131 * t189 - t179 * t185); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124 * t123, -t123 ^ 2 + t124 ^ 2, t61 - t295, t124 * t162 + t198, t120, -t112 * t265 - g(1) * t103 - g(2) * t105 - t113 * t124 + t162 * t58 + t63 + (t227 + t313) * t184, g(1) * t104 - g(2) * t106 + g(3) * t282 - t113 * t123 + t162 * t57 - t209, t336, t339, t332, t322, t115, -t231 * t153 + (t187 * t115 - t124 * t68 - t153 * t263) * pkin(4) + t323, t308 * t153 + (-t183 * t115 - t124 * t217 - t153 * t262) * pkin(4) + t331, t11 * t217 + t14 * t217 + t15 * t68 + t168 * t19 - t68 * t8 + (-t183 * t20 + (t183 * t217 - t187 * t68) * qJD(5)) * pkin(4), t1 * t168 - t11 * t15 - t8 * t14 - pkin(5) * t330 - g(1) * (-t136 * t287 + t137 * t190) - g(2) * (t136 * t285 + t137 * t186) + t136 * t313 + (-t32 * t124 + t2 * t183 + (t11 * t187 - t183 * t8) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t336, t339, t332, t322, t115, t153 * t17 + t323, t153 * t16 + t331, pkin(5) * t19 + t317 * t68, -t317 * t11 + (t1 + t324 - t330) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65 - t321, t11 * t68 + t189 * t214 + t217 * t8 + t194 - t314;];
tau_reg  = t6;
