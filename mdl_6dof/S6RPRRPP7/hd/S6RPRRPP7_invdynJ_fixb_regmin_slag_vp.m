% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPP7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:52:15
% EndTime: 2019-03-09 04:52:25
% DurationCPUTime: 4.56s
% Computational Cost: add. (3882->496), mult. (7245->582), div. (0->0), fcn. (4332->6), ass. (0->236)
t158 = sin(qJ(3));
t264 = qJD(1) * t158;
t125 = qJD(4) + t264;
t157 = sin(qJ(4));
t160 = cos(qJ(4));
t253 = t160 * qJD(3);
t161 = cos(qJ(3));
t263 = qJD(1) * t161;
t94 = t157 * t263 - t253;
t297 = t94 * t125;
t257 = qJD(4) * t161;
t230 = t157 * t257;
t247 = t161 * qJDD(1);
t32 = (t158 * t253 + t230) * qJD(1) - qJD(4) * t253 - t157 * qJDD(3) - t160 * t247;
t344 = -t32 - t297;
t343 = t32 - t297;
t262 = qJD(3) * t157;
t96 = t160 * t263 + t262;
t296 = t96 * t125;
t234 = t157 * t264;
t290 = qJD(4) * t96;
t33 = -qJD(3) * t234 - t160 * qJDD(3) + t157 * t247 + t290;
t342 = -t33 - t296;
t216 = qJD(4) * t158 + qJD(1);
t341 = (-t157 * t216 + t161 * t253) * t125;
t251 = qJD(1) * qJD(3);
t227 = t161 * t251;
t248 = t158 * qJDD(1);
t92 = qJDD(4) + t227 + t248;
t340 = t92 * qJ(5) + t125 * qJD(5);
t261 = qJD(3) * t158;
t287 = t157 * t158;
t339 = -t161 * t33 + t94 * t261 - t92 * t287;
t320 = pkin(4) + pkin(5);
t240 = t320 * t157;
t292 = qJ(5) * t160;
t194 = t240 - t292;
t219 = -t33 + t290;
t220 = qJD(4) * t94 - t32;
t260 = qJD(3) * t161;
t337 = (qJD(1) * t94 + t158 * t220 + t96 * t260) * t157 + (qJD(1) * t96 + t158 * t219 - t94 * t260) * t160;
t288 = t157 * qJ(5);
t193 = -t320 * t160 - t288;
t336 = -pkin(3) + t193;
t335 = t33 * qJ(6) + t94 * qJD(6);
t334 = 0.2e1 * t340;
t91 = t96 ^ 2;
t333 = -t125 ^ 2 - t91;
t258 = qJD(4) * t160;
t228 = t125 * t258;
t233 = t157 * t260;
t289 = t125 * t160;
t332 = -qJD(1) * t289 - t125 * t233 - t158 * t228 + t339;
t148 = t161 * pkin(8);
t328 = t158 * pkin(3) + qJ(2) - t148;
t70 = t328 * qJD(1);
t164 = -pkin(1) - pkin(7);
t121 = t164 * qJD(1) + qJD(2);
t102 = t158 * t121;
t77 = qJD(3) * pkin(8) + t102;
t29 = -t157 * t77 + t160 * t70;
t272 = qJD(5) - t29;
t329 = t157 * qJD(5) + t102;
t162 = cos(qJ(1));
t275 = t161 * t162;
t133 = g(2) * t275;
t327 = g(3) * t158 + t133;
t81 = t92 * pkin(4);
t326 = t81 - qJDD(5);
t276 = t161 * t121;
t78 = -qJD(3) * pkin(3) - t276;
t183 = t96 * qJ(5) - t78;
t23 = pkin(4) * t94 - t183;
t319 = pkin(8) * t92;
t325 = t125 * t23 - t319;
t13 = -t320 * t94 + qJD(6) + t183;
t259 = qJD(4) * t157;
t211 = pkin(3) * t161 + pkin(8) * t158;
t93 = t211 * qJD(3) + qJD(2);
t39 = t93 * qJD(1) + qJDD(1) * t328;
t118 = t164 * qJDD(1) + qJDD(2);
t45 = qJDD(3) * pkin(8) + t118 * t158 + t121 * t260;
t221 = -t157 * t45 + t160 * t39 - t77 * t258 - t70 * t259;
t286 = t157 * t161;
t274 = t162 * t160;
t159 = sin(qJ(1));
t281 = t159 * t157;
t71 = t158 * t281 - t274;
t280 = t159 * t160;
t283 = t158 * t162;
t73 = t157 * t283 + t280;
t179 = g(1) * t71 - g(2) * t73 + g(3) * t286 + t221;
t176 = t179 + t326;
t300 = t32 * qJ(6);
t324 = (qJD(6) + t13) * t96 + t176 - t300;
t321 = t94 ^ 2;
t318 = t92 * pkin(5);
t317 = pkin(4) * t157;
t316 = g(1) * t159;
t153 = g(2) * t162;
t315 = t96 * t94;
t314 = pkin(8) - qJ(6);
t313 = -t125 * t194 + t329;
t30 = t157 * t70 + t160 * t77;
t252 = t160 * qJD(6);
t99 = t211 * qJD(1);
t309 = t157 * t99 + t160 * t276;
t36 = qJ(5) * t263 + t309;
t312 = qJ(6) * t234 - t314 * t259 - t252 - t36;
t205 = -t292 + t317;
t311 = t125 * t205 - t329;
t111 = t314 * t160;
t88 = t157 * t276;
t222 = -t160 * t99 + t88;
t239 = t320 * t161;
t284 = t158 * t160;
t310 = qJD(4) * t111 - t157 * qJD(6) - (qJ(6) * t284 - t239) * qJD(1) - t222;
t308 = pkin(8) * qJD(4);
t307 = qJ(5) * t33;
t117 = t125 * qJ(5);
t17 = qJ(6) * t94 + t30;
t12 = t117 + t17;
t306 = t12 * t125;
t21 = t117 + t30;
t305 = t125 * t21;
t304 = t125 * t30;
t303 = t160 * t92;
t302 = t161 * t32;
t299 = t32 * t157;
t298 = t94 * qJ(5);
t282 = t158 * t164;
t295 = t157 * t328 + t160 * t282;
t294 = pkin(1) * qJDD(1);
t166 = qJD(1) ^ 2;
t293 = qJ(2) * t166;
t291 = qJD(3) * t96;
t285 = t158 * t159;
t279 = t159 * t161;
t278 = t160 * t161;
t277 = t161 * qJ(6);
t16 = t96 * qJ(6) + t29;
t273 = qJD(5) - t16;
t270 = g(3) * t287 + t157 * t133;
t269 = g(3) * t284 + t160 * t133;
t268 = g(1) * t285 + g(3) * t161;
t267 = t162 * pkin(1) + t159 * qJ(2);
t155 = t161 ^ 2;
t266 = t158 ^ 2 - t155;
t165 = qJD(3) ^ 2;
t265 = -t165 - t166;
t256 = qJD(4) * t164;
t255 = qJD(5) * t160;
t250 = qJDD(1) * qJ(2);
t249 = qJDD(3) * t158;
t246 = -t157 * t39 - t160 * t45 - t70 * t258;
t245 = g(1) * t279;
t244 = g(2) * t283;
t231 = t164 * t260;
t243 = t157 * t93 + t160 * t231 + t258 * t328;
t242 = 0.2e1 * qJD(1) * qJD(2);
t100 = t161 * t118;
t44 = -qJDD(3) * pkin(3) + t121 * t261 - t100;
t229 = t158 * t256;
t241 = t157 * t231 + t160 * t229 + t259 * t328;
t47 = t158 * qJ(5) + t295;
t238 = t13 * t259;
t237 = t13 * t258;
t236 = t159 * t278;
t225 = -t153 + t316;
t72 = t157 * t162 + t158 * t280;
t224 = -t71 * pkin(4) + qJ(5) * t72;
t74 = t158 * t274 - t281;
t223 = t73 * pkin(4) - qJ(5) * t74;
t119 = t157 * t282;
t218 = t160 * t328 - t119;
t217 = qJDD(2) - t294;
t215 = pkin(4) * t236 + pkin(8) * t285 + (pkin(3) + t288) * t279;
t6 = t33 * pkin(4) + t32 * qJ(5) - t96 * qJD(5) + t44;
t3 = -pkin(5) * t33 + qJDD(6) - t6;
t214 = t3 - t245;
t213 = g(1) * t73 + g(2) * t71;
t212 = -g(1) * t74 - g(2) * t72;
t210 = g(1) * t162 + g(2) * t159;
t209 = qJ(5) * t260 + t158 * qJD(5) + t243;
t9 = -t320 * t125 + t273;
t208 = t12 * t160 + t157 * t9;
t207 = -t12 * t157 + t160 * t9;
t206 = pkin(4) * t160 + t288;
t204 = -t293 - t316;
t20 = -pkin(4) * t125 + t272;
t203 = t157 * t21 - t160 * t20;
t202 = t157 * t20 + t160 * t21;
t201 = t244 - t268;
t200 = t160 * t93 - t241;
t5 = -t221 - t326;
t199 = pkin(3) + t206;
t196 = -t96 * t261 - t302;
t190 = t157 * t92 + t228;
t189 = -t125 * t259 + t303;
t188 = -t77 * t259 - t246;
t187 = 0.2e1 * qJ(2) * t251 + qJDD(3) * t164;
t186 = t125 * t308 + t245;
t185 = pkin(3) * t285 + t72 * pkin(4) + t162 * pkin(7) + t71 * qJ(5) + t267;
t182 = -t186 - t6;
t181 = t125 * t78 - t319;
t4 = t188 + t340;
t147 = t162 * qJ(2);
t180 = pkin(3) * t283 + t74 * pkin(4) - pkin(8) * t275 + t73 * qJ(5) + t147;
t178 = -t210 + t242 + 0.2e1 * t250;
t177 = -t92 + t315;
t175 = -t203 * qJD(4) + t5 * t157 + t4 * t160;
t173 = -t164 * t165 + t178;
t172 = t23 * t96 - t176;
t170 = g(1) * t72 - g(2) * t74 + g(3) * t278 - t188;
t168 = t92 * t284 + t196 + t341;
t167 = t125 * t29 + t170;
t142 = qJDD(3) * t161;
t124 = qJ(5) * t278;
t110 = t314 * t157;
t50 = -t124 + (-t164 + t317) * t161;
t48 = -t158 * pkin(4) - t218;
t46 = t124 + (t164 - t240) * t161;
t42 = pkin(4) * t96 + t298;
t38 = t157 * t277 + t47;
t37 = -pkin(4) * t263 + t222;
t27 = t119 + (-t328 - t277) * t160 - t320 * t158;
t24 = -t320 * t96 - t298;
t19 = (t206 * qJD(4) - t255) * t161 + (t164 - t205) * t261;
t14 = -pkin(4) * t260 - t200;
t11 = -t157 * t229 + t209;
t10 = (t193 * qJD(4) + t255) * t161 + (-t164 + t194) * t261;
t8 = t160 * qJ(6) * t257 + (qJD(6) * t161 + (-qJ(6) * qJD(3) - t256) * t158) * t157 + t209;
t7 = (qJ(6) * t261 - t93) * t160 + (qJ(6) * t259 - t320 * qJD(3) - t252) * t161 + t241;
t2 = t4 + t335;
t1 = -t96 * qJD(6) + t300 - t318 + t5;
t15 = [qJDD(1), t225, t210, qJDD(2) - t225 - 0.2e1 * t294, t178, -t217 * pkin(1) - g(1) * (-t159 * pkin(1) + t147) - g(2) * t267 + (t242 + t250) * qJ(2), qJDD(1) * t155 - 0.2e1 * t158 * t227, -0.2e1 * t158 * t247 + 0.2e1 * t266 * t251, -t158 * t165 + t142, -t161 * t165 - t249, 0, t173 * t158 + t187 * t161, -t187 * t158 + t161 * t173, t196 * t160 - t230 * t96 (t157 * t96 + t160 * t94) * t261 + (t299 - t160 * t33 + (t157 * t94 - t160 * t96) * qJD(4)) * t161 (-t125 * t253 - t32) * t158 + (t189 + t291) * t161 (t125 * t262 - t33) * t158 + (-qJD(3) * t94 - t190) * t161, t125 * t260 + t158 * t92, t200 * t125 + t218 * t92 + ((-t157 * t78 + t164 * t94) * qJD(3) + t221) * t158 + (qJD(3) * t29 + t44 * t157 - t164 * t33 + t258 * t78) * t161 + t212, -t243 * t125 - t295 * t92 + ((t125 * t164 + t77) * t259 + (-t160 * t78 + t164 * t96) * qJD(3) + t246) * t158 + (-qJD(3) * t30 + t44 * t160 + t164 * t32 - t259 * t78) * t161 + t213, -t14 * t125 + t19 * t94 + t50 * t33 - t48 * t92 + (-t23 * t262 - t5) * t158 + (-qJD(3) * t20 + t6 * t157 + t23 * t258) * t161 + t212, -t11 * t94 + t14 * t96 - t48 * t32 - t47 * t33 + t203 * t261 + (-qJD(4) * t202 - t157 * t4 + t160 * t5 + t210) * t161, t11 * t125 - t19 * t96 + t50 * t32 + t47 * t92 + (t23 * t253 + t4) * t158 + (qJD(3) * t21 - t6 * t160 + t23 * t259) * t161 - t213, t4 * t47 + t21 * t11 + t6 * t50 + t23 * t19 + t5 * t48 + t20 * t14 - g(1) * (t164 * t159 + t180) - g(2) * (-pkin(8) * t279 + t185) -t10 * t94 - t7 * t125 - t27 * t92 - t46 * t33 + (t13 * t262 - t1) * t158 + (-qJD(3) * t9 - t3 * t157 - t237) * t161 + t212, t10 * t96 + t8 * t125 - t46 * t32 + t38 * t92 + (-t13 * t253 + t2) * t158 + (qJD(3) * t12 + t3 * t160 - t238) * t161 - t213, t27 * t32 + t38 * t33 - t7 * t96 + t8 * t94 + t207 * t261 + (qJD(4) * t208 - t1 * t160 + t157 * t2 - t210) * t161, t2 * t38 + t12 * t8 + t1 * t27 + t9 * t7 + t3 * t46 + t13 * t10 - g(1) * (t74 * pkin(5) + qJ(6) * t275 + t180) - g(2) * (t72 * pkin(5) + t185) + (g(2) * t314 * t161 - g(1) * t164) * t159; 0, 0, 0, qJDD(1), -t166, t153 + t204 + t217, 0, 0, 0, 0, 0, t265 * t158 + t142, t161 * t265 - t249, 0, 0, 0, 0, 0 (-t160 * t216 - t233) * t125 + t339, t302 + (t291 - t303) * t158 - t341, t332, t337, t168, -t203 * qJD(1) + (qJD(3) * t202 - t6) * t161 + (qJD(3) * t23 + t175) * t158 - t225, t332, t168, -t337, t207 * qJD(1) + (qJD(3) * t208 + t3) * t161 + (-qJD(3) * t13 + qJD(4) * t207 + t1 * t157 + t2 * t160) * t158 - t225; 0, 0, 0, 0, 0, 0, t161 * t166 * t158, -t266 * t166, t247, -t248, qJDD(3), t204 * t161 + t100 + t327 (-t118 + t293 - t153) * t158 + t268, t96 * t289 - t299, t342 * t157 + t160 * t344 (t125 * t284 - t161 * t96) * qJD(1) + t190 (-t125 * t287 + t161 * t94) * qJD(1) + t189, -t125 * t263, -t29 * t263 - t94 * t102 - pkin(3) * t33 + t88 * t125 + (-t245 - t44 + (-t99 - t308) * t125) * t160 + t181 * t157 + t269, pkin(3) * t32 + t309 * t125 + t30 * t263 - t96 * t102 + t181 * t160 + (t186 + t44) * t157 - t270, t37 * t125 + t325 * t157 + t182 * t160 - t199 * t33 + t20 * t263 + t311 * t94 + t269, t36 * t94 - t37 * t96 + (pkin(8) * t219 + t125 * t20 + t4) * t160 + (pkin(8) * t220 - t305 + t5) * t157 + t201, -t36 * t125 + t182 * t157 - t325 * t160 - t199 * t32 - t21 * t263 - t311 * t96 + t270, -t21 * t36 - t20 * t37 - g(1) * t215 - g(3) * t148 + t311 * t23 + (t175 + t244) * pkin(8) + (-t6 + t327) * t199, -t238 - t110 * t92 + t336 * t33 - t313 * t94 + t214 * t160 - t310 * t125 + (-t13 * t287 + t161 * t9) * qJD(1) + t269, t237 + t111 * t92 + t336 * t32 + t313 * t96 + t214 * t157 + t312 * t125 + (-t12 * t161 + t13 * t284) * qJD(1) + t270, t110 * t32 + t111 * t33 - t310 * t96 + t312 * t94 + (-t125 * t9 - t2) * t160 + (-t1 + t306) * t157 - t201, t2 * t111 + t1 * t110 - t3 * t336 - g(1) * (pkin(5) * t236 + t215) - g(3) * (t148 - t277) + t310 * t9 + t313 * t13 + t312 * t12 + (qJ(6) * t316 - g(3) * (-pkin(5) * t160 - t199)) * t158 - (-t314 * t158 + t336 * t161) * t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t315, t91 - t321, -t343, -t33 + t296, t92, -t78 * t96 + t179 + t304, t78 * t94 + t167, -t42 * t94 - t172 + t304 + t81, pkin(4) * t32 - t307 + (t21 - t30) * t96 + (t20 - t272) * t94, -t23 * t94 + t42 * t96 - t167 + t334, t4 * qJ(5) - t5 * pkin(4) - t23 * t42 - t20 * t30 - g(1) * t224 - g(2) * t223 - g(3) * (-pkin(4) * t286 + t124) + t272 * t21 (pkin(5) + t320) * t92 + t24 * t94 + t17 * t125 + t324, -t125 * t16 + t13 * t94 - t24 * t96 - t170 + t334 + t335, t307 - t320 * t32 + (-t12 + t17) * t96 + (-t9 + t273) * t94, t2 * qJ(5) - t1 * t320 - t9 * t17 - t13 * t24 - g(1) * (-pkin(5) * t71 + t224) - g(2) * (pkin(5) * t73 + t223) - g(3) * (-t157 * t239 + t124) + t273 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, -t343, t333, t172 - t305, t177, t333, t343, -t306 - t318 - t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t342, t344, -t91 - t321, -t12 * t94 + t9 * t96 + t214 + t327;];
tau_reg  = t15;
