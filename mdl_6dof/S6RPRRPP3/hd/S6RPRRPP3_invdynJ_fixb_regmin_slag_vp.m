% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:37:02
% EndTime: 2019-03-09 04:37:13
% DurationCPUTime: 4.51s
% Computational Cost: add. (4139->486), mult. (8160->570), div. (0->0), fcn. (5100->10), ass. (0->232)
t173 = cos(qJ(3));
t273 = qJD(1) * t173;
t139 = -qJD(4) + t273;
t172 = cos(qJ(4));
t265 = t172 * qJD(3);
t242 = t173 * t265;
t169 = sin(qJ(4));
t170 = sin(qJ(3));
t267 = qJD(4) * t170;
t245 = t169 * t267;
t285 = t170 * t172;
t156 = t173 * qJDD(1);
t264 = qJD(1) * qJD(3);
t93 = t170 * t264 + qJDD(4) - t156;
t204 = -t139 * (t242 - t245) + t93 * t285;
t270 = qJD(3) * t170;
t241 = t173 * t264;
t262 = t170 * qJDD(1);
t336 = qJD(1) * t267 - qJDD(3);
t40 = -qJD(4) * t265 + (-t241 - t262) * t172 + t336 * t169;
t271 = qJD(3) * t169;
t274 = qJD(1) * t170;
t99 = t172 * t274 + t271;
t332 = -t173 * t40 - t270 * t99;
t342 = t204 + t332;
t311 = pkin(4) + qJ(6);
t247 = t311 * t93;
t166 = sin(pkin(9));
t147 = pkin(1) * t166 + pkin(7);
t117 = t147 * qJD(1);
t76 = qJD(2) * t173 - t117 * t170;
t341 = qJD(3) * t76;
t115 = t147 * qJDD(1);
t340 = -qJD(2) * qJD(3) - t115;
t287 = t169 * t170;
t266 = qJD(4) * t172;
t269 = qJD(3) * t173;
t338 = t169 * t269 + t170 * t266;
t183 = t139 * t338 - t287 * t93;
t41 = ((qJD(4) + t273) * qJD(3) + t262) * t169 + t336 * t172;
t97 = t169 * t274 - t265;
t235 = t173 * t41 - t270 * t97;
t180 = t183 - t235;
t163 = qJ(1) + pkin(9);
t153 = cos(t163);
t289 = t153 * t170;
t152 = sin(t163);
t291 = t152 * t170;
t337 = g(1) * t289 + g(2) * t291;
t77 = qJD(2) * t170 + t117 * t173;
t60 = qJD(3) * pkin(8) + t77;
t159 = t170 * pkin(8);
t161 = t173 * pkin(3);
t249 = -pkin(2) - t161;
t213 = t249 - t159;
t167 = cos(pkin(9));
t319 = pkin(1) * t167;
t197 = t213 - t319;
t65 = t197 * qJD(1);
t24 = t169 * t60 - t172 * t65;
t212 = pkin(5) * t99 + t24;
t282 = qJD(5) + t212;
t127 = qJD(5) * t139;
t86 = t93 * qJ(5);
t335 = t127 - t86;
t136 = t139 ^ 2;
t92 = t99 ^ 2;
t334 = -t92 - t136;
t268 = qJD(4) * t169;
t333 = -qJD(5) * t169 - t77 + (-t169 * t273 + t268) * pkin(4);
t331 = g(1) * t153 + g(2) * t152;
t330 = -pkin(5) * t41 + qJDD(6);
t327 = t97 ^ 2;
t326 = 0.2e1 * t86;
t325 = pkin(5) + pkin(8);
t324 = pkin(4) * t93;
t323 = pkin(5) * t40;
t321 = pkin(5) * t97;
t320 = pkin(8) * t93;
t318 = g(1) * t152;
t315 = g(2) * t153;
t314 = g(3) * t170;
t313 = g(3) * t173;
t312 = t99 * t97;
t294 = qJ(5) * t172;
t214 = qJ(6) * t169 - t294;
t202 = t214 * t173;
t310 = -qJD(1) * t202 + qJD(4) * t214 - qJD(6) * t172 + t333;
t25 = t169 * t65 + t172 * t60;
t222 = pkin(3) * t170 - pkin(8) * t173;
t103 = t222 * qJD(1);
t309 = t103 * t169 + t172 * t76;
t106 = t222 * qJD(3);
t148 = -pkin(2) - t319;
t276 = t159 + t161;
t88 = t148 - t276;
t308 = t106 * t169 + t266 * t88;
t307 = -qJ(5) * t266 + t273 * t294 + t333;
t306 = pkin(8) * qJD(4);
t305 = qJ(5) * t41;
t304 = qJ(5) * t97;
t18 = qJ(5) * t139 - t25;
t303 = t139 * t18;
t302 = t139 * t25;
t301 = t139 * t97;
t300 = t139 * t99;
t284 = t172 * t173;
t102 = t147 * t284;
t297 = t169 * t88 + t102;
t286 = t169 * t173;
t260 = pkin(5) * t286;
t296 = -t325 * t268 - (qJ(5) * t170 - t260) * qJD(1) - t309;
t121 = t325 * t172;
t57 = t169 * t76;
t231 = -t103 * t172 + t57;
t259 = pkin(5) * t284;
t295 = qJD(4) * t121 - (-t170 * t311 + t259) * qJD(1) - t231;
t293 = t147 * t169;
t292 = t147 * t172;
t290 = t152 * t173;
t288 = t153 * t173;
t283 = -qJD(5) - t24;
t16 = t25 - t321;
t281 = -qJD(6) - t16;
t280 = qJDD(2) - g(3);
t279 = t337 * t169;
t278 = t337 * t172;
t277 = pkin(4) * t287 - qJ(5) * t285;
t164 = t170 ^ 2;
t275 = -t173 ^ 2 + t164;
t118 = qJD(1) * t148;
t261 = -t40 * t287 + t338 * t99;
t256 = -t117 * t269 + t170 * t340;
t59 = -qJD(3) * pkin(3) - t76;
t193 = -qJ(5) * t99 + t59;
t14 = t311 * t97 + t193;
t252 = t14 * t268;
t251 = t14 * t266;
t101 = t147 * t286;
t36 = -qJDD(3) * pkin(3) - qJDD(2) * t173 - t256;
t182 = qJ(5) * t40 - qJD(5) * t99 + t36;
t3 = qJD(6) * t97 + t311 * t41 + t182;
t248 = -t3 - t313;
t246 = t147 * t268;
t239 = -qJ(5) * t169 - pkin(3);
t238 = -pkin(4) - t293;
t69 = t152 * t286 + t153 * t172;
t70 = t152 * t284 - t153 * t169;
t237 = -pkin(4) * t69 + qJ(5) * t70;
t71 = -t152 * t172 + t153 * t286;
t72 = t152 * t169 + t153 * t284;
t236 = -pkin(4) * t71 + qJ(5) * t72;
t35 = qJDD(3) * pkin(8) + qJDD(2) * t170 + t115 * t173 + t341;
t42 = qJD(1) * t106 + qJDD(1) * t197;
t234 = t169 * t35 - t172 * t42 + t266 * t60 + t268 * t65;
t233 = -t169 * t42 - t172 * t35 - t266 * t65 + t268 * t60;
t232 = t172 * t88 - t101;
t230 = -qJ(5) + t292;
t66 = t97 * t242;
t227 = pkin(4) * t338 + qJ(5) * t245 + t147 * t269;
t226 = pkin(4) * t284 + qJ(5) * t286 + t276;
t225 = g(1) * t69 - g(2) * t71;
t224 = g(1) * t70 - g(2) * t72;
t223 = -qJD(4) * t102 + t106 * t172 - t268 * t88;
t44 = qJ(5) * t173 - t297;
t171 = sin(qJ(1));
t174 = cos(qJ(1));
t220 = g(1) * t171 - g(2) * t174;
t219 = -qJDD(5) - t234;
t11 = qJD(6) - t18 - t321;
t9 = t139 * t311 + t282;
t218 = t11 * t172 + t169 * t9;
t217 = t11 * t169 - t172 * t9;
t17 = pkin(4) * t139 - t283;
t216 = t169 * t18 + t17 * t172;
t215 = t169 * t17 - t172 * t18;
t4 = t233 + t335;
t211 = pkin(4) * t172 - t239;
t210 = -qJ(6) * t287 - t277;
t209 = -t139 * t306 + t313;
t206 = t245 * t97 - t285 * t41 - t66;
t6 = pkin(4) * t41 + t182;
t205 = -t209 - t6;
t203 = -pkin(1) * t171 - pkin(4) * t70 + pkin(7) * t153 - qJ(5) * t69;
t201 = t172 * t311 - t239;
t200 = -qJD(1) * t118 + t331;
t196 = -t139 * t59 - t320;
t23 = pkin(4) * t97 + t193;
t195 = t139 * t23 + t320;
t194 = pkin(1) * t174 + pkin(2) * t153 + pkin(3) * t288 + pkin(4) * t72 + pkin(7) * t152 + pkin(8) * t289 + qJ(5) * t71;
t192 = t93 - t312;
t175 = qJD(3) ^ 2;
t191 = 0.2e1 * qJDD(1) * t148 + t147 * t175 + t315;
t190 = g(1) * t71 + g(2) * t69 + g(3) * t287 - t234;
t189 = 0.2e1 * qJD(3) * t118 - qJDD(3) * t147;
t188 = -qJDD(5) + t190;
t187 = g(1) * t72 + g(2) * t70 + g(3) * t285 + t233;
t20 = -t40 - t301;
t1 = qJD(6) * t139 - t219 - t247 - t323;
t2 = -t4 + t330;
t186 = -qJD(4) * t217 + t1 * t169 + t172 * t2;
t5 = -t219 - t324;
t185 = qJD(4) * t216 + t5 * t169 - t4 * t172;
t184 = t23 * t99 - t188;
t181 = t14 * t99 - t188 - t323;
t179 = -t14 * t97 - t187 + t330;
t178 = -t300 - t41;
t176 = qJD(1) ^ 2;
t160 = t173 * pkin(4);
t126 = g(1) * t291;
t124 = pkin(8) * t288;
t122 = pkin(8) * t290;
t120 = t325 * t169;
t119 = t170 * t147;
t114 = qJDD(3) * t173 - t170 * t175;
t113 = qJDD(3) * t170 + t173 * t175;
t54 = t119 + t277;
t47 = pkin(4) * t99 + t304;
t46 = t119 - t210;
t45 = t160 - t232;
t37 = -pkin(5) * t287 - t44;
t32 = -pkin(4) * t274 + t231;
t31 = -qJ(5) * t274 - t309;
t29 = t311 * t99 + t304;
t26 = qJ(6) * t173 + t101 + t160 + (pkin(5) * t170 - t88) * t172;
t21 = (-qJ(5) * t269 - qJD(5) * t170) * t172 + t227;
t13 = t238 * t270 - t223;
t12 = (qJD(5) + t246) * t173 + t230 * t270 - t308;
t10 = qJD(3) * t202 + (qJD(6) * t169 + (qJ(6) * qJD(4) - qJD(5)) * t172) * t170 + t227;
t8 = -qJD(5) * t173 + (-pkin(5) * t285 - t101) * qJD(4) + (-t170 * t230 - t260) * qJD(3) + t308;
t7 = -pkin(5) * t245 + qJD(6) * t173 + (t259 + (-qJ(6) + t238) * t170) * qJD(3) - t223;
t15 = [qJDD(1), t220, g(1) * t174 + g(2) * t171 (t220 + (t166 ^ 2 + t167 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t164 + 0.2e1 * t170 * t241, 0.2e1 * t156 * t170 - 0.2e1 * t264 * t275, t113, t114, 0, t189 * t170 + (-t191 + t318) * t173, t170 * t191 + t173 * t189 - t126, t99 * t242 + (-t172 * t40 - t268 * t99) * t170, t206 - t261, t204 - t332, t183 + t235, -t139 * t270 - t173 * t93, -t223 * t139 + t232 * t93 + ((t147 * t97 + t169 * t59) * qJD(3) + t234) * t173 + (t59 * t266 + t147 * t41 + t36 * t169 + (-t139 * t293 - t24) * qJD(3)) * t170 + t224, t308 * t139 - t297 * t93 + (-t139 * t246 + (t147 * t99 + t172 * t59) * qJD(3) - t233) * t173 + (-t59 * t268 - t147 * t40 + t36 * t172 + (-t139 * t292 - t25) * qJD(3)) * t170 - t225, t12 * t97 + t13 * t99 - t40 * t45 + t41 * t44 + t126 + t216 * t269 + (-qJD(4) * t215 + t169 * t4 + t172 * t5 - t315) * t170, -t13 * t139 - t21 * t97 - t41 * t54 + t45 * t93 + (-t23 * t271 - t5) * t173 + (qJD(3) * t17 - t169 * t6 - t23 * t266) * t170 - t224, t12 * t139 - t21 * t99 + t40 * t54 - t44 * t93 + (-t23 * t265 + t4) * t173 + (-qJD(3) * t18 - t172 * t6 + t23 * t268) * t170 + t225, -g(1) * t203 - g(2) * t194 + t18 * t12 + t17 * t13 + t23 * t21 - t213 * t318 + t4 * t44 + t5 * t45 + t6 * t54, -t26 * t40 - t37 * t41 + t7 * t99 - t8 * t97 + t126 - t217 * t269 + (-qJD(4) * t218 + t1 * t172 - t169 * t2 - t315) * t170, -t10 * t99 - t139 * t8 + t37 * t93 + t40 * t46 + (-t14 * t265 - t2) * t173 + (qJD(3) * t11 - t172 * t3 + t252) * t170 + t225, t10 * t97 + t139 * t7 - t26 * t93 + t41 * t46 + (t14 * t271 + t1) * t173 + (-qJD(3) * t9 + t169 * t3 + t251) * t170 + t224, t3 * t46 + t14 * t10 + t1 * t26 + t9 * t7 + t2 * t37 + t11 * t8 - g(1) * (-qJ(6) * t70 + t203) - g(2) * (pkin(5) * t289 + qJ(6) * t72 + t194) - (-t170 * t325 + t249) * t318; 0, 0, 0, t280, 0, 0, 0, 0, 0, t114, -t113, 0, 0, 0, 0, 0, t180, -t342, -t66 + (-t172 * t41 + t268 * t97) * t170 + t261, -t180, t342, -g(3) + (qJD(3) * t215 - t6) * t173 + (qJD(3) * t23 + t185) * t170, t206 + t261, t342, t180, -g(3) + (qJD(3) * t218 - t3) * t173 + (qJD(3) * t14 + t186) * t170; 0, 0, 0, 0, -t170 * t176 * t173, t275 * t176, t262, t156, qJDD(3), qJD(3) * t77 + t170 * t200 + t173 * t280 + t256, t341 + (qJD(3) * t117 - t280) * t170 + (t200 + t340) * t173, -t169 * t40 - t172 * t300 (-t40 + t301) * t172 + (-t41 + t300) * t169, -t139 * t266 + t169 * t93 + (t139 * t284 - t170 * t99) * qJD(1), t139 * t268 + t172 * t93 + (-t139 * t286 + t170 * t97) * qJD(1), t139 * t274, t24 * t274 - pkin(3) * t41 - t57 * t139 - t77 * t97 + (-t313 - t36 + (t103 + t306) * t139) * t172 + t196 * t169 + t278, pkin(3) * t40 - t309 * t139 + t25 * t274 - t77 * t99 + t196 * t172 + (t209 + t36) * t169 - t279, -t314 - t31 * t97 - t32 * t99 - t331 * t173 + (-t4 - t139 * t17 + (qJD(4) * t99 - t41) * pkin(8)) * t172 + (t5 - t303 + (qJD(4) * t97 - t40) * pkin(8)) * t169, t139 * t32 + t169 * t195 - t17 * t274 - t172 * t205 + t211 * t41 - t307 * t97 - t278, -t139 * t31 + t169 * t205 + t172 * t195 + t18 * t274 - t211 * t40 - t307 * t99 + t279, t185 * pkin(8) - g(1) * t124 - g(2) * t122 - g(3) * t226 - t17 * t32 - t18 * t31 + t307 * t23 + (t170 * t331 - t6) * t211, -t314 - t120 * t40 - t121 * t41 + t295 * t99 - t296 * t97 + (qJD(1) * t217 - t331) * t173 + t186, -t251 + t121 * t93 - t40 * t201 - t310 * t99 + t248 * t169 - t296 * t139 + (-t11 * t170 + t14 * t284) * qJD(1) + t279, t252 - t120 * t93 - t41 * t201 + t310 * t97 + t248 * t172 + t295 * t139 + (-t14 * t286 + t170 * t9) * qJD(1) + t278, -t3 * t201 + t1 * t120 + t2 * t121 - g(1) * (pkin(5) * t288 + t124) - g(2) * (pkin(5) * t290 + t122) - g(3) * (qJ(6) * t284 + t226) + t295 * t9 + t310 * t14 + t296 * t11 + (-g(3) * pkin(5) + t201 * t331) * t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t312, t92 - t327, t20, t178, t93, -t59 * t99 + t190 - t302, t139 * t24 + t59 * t97 + t187, pkin(4) * t40 - t305 + (-t18 - t25) * t99 + (t17 + t283) * t97, t47 * t97 + t184 + t302 - 0.2e1 * t324, t139 * t283 - t23 * t97 + t47 * t99 - t127 - t187 + t326, -t5 * pkin(4) - g(1) * t236 - g(2) * t237 + g(3) * t277 - t4 * qJ(5) - t17 * t25 + t18 * t283 - t23 * t47, -t305 + t311 * t40 + (t11 + t281) * t99 + (t9 - t282) * t97, -t139 * t212 + t29 * t99 - 0.2e1 * t127 + t179 + t326, -t29 * t97 + (-0.2e1 * qJD(6) - t16) * t139 + 0.2e1 * t247 - t181, -t1 * t311 + t2 * qJ(5) - t14 * t29 - g(1) * (-qJ(6) * t71 + t236) - g(2) * (-qJ(6) * t69 + t237) - g(3) * t210 + t281 * t9 + t282 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t192, t334, t184 - t303 - t324, t20, t334, -t192, -t247 + (qJD(6) + t11) * t139 + t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, t93 + t312, -t136 - t327, -t139 * t9 + t179 - t335;];
tau_reg  = t15;
