% Calculate minimal parameter regressor of coriolis joint torque vector for
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
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRP12_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:53:32
% EndTime: 2019-03-09 12:53:46
% DurationCPUTime: 5.18s
% Computational Cost: add. (5743->432), mult. (12622->561), div. (0->0), fcn. (7884->6), ass. (0->227)
t179 = sin(qJ(2));
t259 = qJD(1) * t179;
t155 = qJD(4) + t259;
t296 = cos(qJ(5));
t231 = t296 * qJD(5);
t322 = t155 * t296 + t231;
t181 = cos(qJ(2));
t178 = sin(qJ(4));
t266 = t178 * t179;
t207 = pkin(4) * t181 - pkin(9) * t266;
t258 = qJD(1) * t181;
t164 = pkin(7) * t258;
t123 = pkin(3) * t258 + t164;
t180 = cos(qJ(4));
t168 = pkin(2) * t259;
t211 = pkin(8) * t179 - qJ(3) * t181;
t94 = qJD(1) * t211 + t168;
t223 = t123 * t180 - t178 * t94;
t253 = qJD(4) * t178;
t182 = -pkin(2) - pkin(8);
t292 = pkin(9) - t182;
t321 = -qJD(1) * t207 + t253 * t292 - t223;
t128 = t292 * t180;
t240 = t180 * t259;
t277 = t123 * t178 + t180 * t94;
t320 = pkin(9) * t240 + qJD(4) * t128 + t277;
t141 = qJD(5) + t155;
t177 = sin(qJ(5));
t250 = qJD(5) * t177;
t252 = qJD(4) * t180;
t287 = -t180 * t250 - t322 * t178 + (-t240 - t252) * t177;
t319 = t287 * t141;
t297 = pkin(3) + pkin(7);
t234 = t180 * t258;
t249 = t178 * qJD(2);
t114 = -t234 - t249;
t256 = qJD(2) * t180;
t115 = -t178 * t258 + t256;
t197 = t114 * t177 + t115 * t296;
t62 = -t114 * t296 + t115 * t177;
t174 = qJD(2) * qJ(3);
t103 = t174 + t123;
t70 = -pkin(4) * t114 + t103;
t23 = pkin(5) * t62 - qJ(6) * t197 + t70;
t318 = t23 * t62;
t317 = t62 * t70;
t267 = t177 * t178;
t286 = -t177 * t253 - t178 * t250 + t180 * t322 - t259 * t267;
t293 = t197 * t62;
t117 = t177 * t180 + t178 * t296;
t233 = t286 * t141;
t316 = (t117 * t258 + t197) * qJD(2) + t233;
t163 = pkin(7) * t259;
t315 = qJD(3) + t163;
t299 = t197 ^ 2;
t314 = -t62 ^ 2 + t299;
t251 = qJD(4) * t181;
t235 = t180 * t251;
t239 = t179 * t249;
t191 = -t235 + t239;
t246 = qJD(2) * qJD(4);
t188 = -qJD(1) * t191 + t178 * t246;
t160 = t180 * t246;
t236 = t178 * t251;
t238 = t179 * t256;
t192 = t236 + t238;
t189 = qJD(1) * t192 - t160;
t29 = -t114 * t231 + t115 * t250 - t177 * t189 + t188 * t296;
t13 = t141 * t62 - t29;
t36 = pkin(5) * t197 + qJ(6) * t62;
t247 = qJD(1) * qJD(2);
t313 = -0.2e1 * t247;
t294 = t23 * t197;
t311 = t70 * t197;
t242 = -pkin(4) * t180 - pkin(3);
t278 = pkin(4) * t252 - t242 * t259 + t315;
t127 = t292 * t178;
t196 = t127 * t177 - t128 * t296;
t310 = -t196 * qJD(5) - t177 * t321 + t296 * t320;
t72 = -t127 * t296 - t128 * t177;
t309 = -t72 * qJD(5) + t177 * t320 + t296 * t321;
t275 = qJ(3) * t179;
t228 = -pkin(1) - t275;
t308 = t181 * t182;
t110 = t228 + t308;
t84 = t110 * qJD(1);
t248 = pkin(3) * t259 + t315;
t89 = qJD(2) * t182 + t248;
t51 = t178 * t89 + t180 * t84;
t43 = pkin(9) * t114 + t51;
t175 = t179 ^ 2;
t176 = t181 ^ 2;
t260 = t175 - t176;
t30 = qJD(5) * t197 - t177 * t188 - t189 * t296;
t305 = t141 * t197 - t30;
t241 = t296 * t180;
t116 = -t241 + t267;
t304 = t116 * t29 + t197 * t287;
t272 = t103 * t179;
t303 = qJD(2) * (t275 - t308) - t272;
t136 = t297 * t179;
t119 = t180 * t136;
t227 = pkin(9) * t181 - t110;
t56 = pkin(4) * t179 + t178 * t227 + t119;
t264 = t180 * t181;
t118 = t178 * t136;
t276 = t110 * t180 + t118;
t59 = -pkin(9) * t264 + t276;
t200 = t177 * t56 + t296 * t59;
t193 = t207 * qJD(2);
t255 = qJD(2) * t181;
t124 = t297 * t255;
t257 = qJD(2) * t179;
t167 = pkin(2) * t257;
t254 = qJD(3) * t179;
t190 = qJD(2) * t211 - t254;
t80 = t167 + t190;
t222 = t124 * t180 - t178 * t80;
t24 = t193 + (t180 * t227 - t118) * qJD(4) + t222;
t194 = -t110 * t253 + t124 * t178 + t136 * t252 + t180 * t80;
t27 = pkin(9) * t192 + t194;
t302 = -qJD(5) * t200 - t177 * t27 + t24 * t296;
t301 = (qJD(4) + qJD(5)) * t181;
t130 = t141 * qJD(6);
t159 = t181 * t247;
t144 = qJ(6) * t159;
t153 = pkin(7) * t159;
t109 = pkin(3) * t159 + t153;
t230 = t179 * t247;
t154 = pkin(2) * t230;
t74 = qJD(1) * t190 + t154;
t226 = t109 * t180 - t178 * t74;
t14 = qJD(1) * t193 - qJD(4) * t43 + t226;
t199 = t109 * t178 + t180 * t74 + t252 * t89 - t253 * t84;
t19 = pkin(9) * t189 + t199;
t50 = -t178 * t84 + t180 * t89;
t42 = -pkin(9) * t115 + t50;
t35 = pkin(4) * t155 + t42;
t225 = -t14 * t177 - t19 * t296 - t231 * t35 + t250 * t43;
t1 = t144 + t130 - t225;
t217 = pkin(5) * t159;
t224 = -t14 * t296 + t177 * t19 + t231 * t43 + t250 * t35;
t2 = -t217 + t224;
t283 = t177 * t43;
t9 = t296 * t35 - t283;
t280 = qJD(6) - t9;
t7 = -pkin(5) * t141 + t280;
t245 = t296 * t43;
t10 = t177 * t35 + t245;
t8 = qJ(6) * t141 + t10;
t300 = t1 * t117 + t116 * t2 + t286 * t8 - t287 * t7;
t298 = t180 ^ 2;
t291 = qJ(6) * t258 + t310;
t290 = -pkin(5) * t258 + t309;
t289 = pkin(5) * t286 - qJ(6) * t287 + qJD(6) * t116 + t278;
t285 = qJD(2) * pkin(2);
t122 = t297 * t257;
t173 = qJD(2) * qJD(3);
t91 = -qJD(1) * t122 + t173;
t282 = t91 * t178;
t281 = t91 * t180;
t16 = t296 * t42 - t283;
t279 = pkin(4) * t231 + qJD(6) - t16;
t274 = qJD(2) * t196;
t273 = qJD(2) * t72;
t271 = t115 * t155;
t270 = t115 * t181;
t269 = t155 * t179;
t268 = t155 * t182;
t265 = t178 * t181;
t184 = qJD(1) ^ 2;
t263 = t181 * t184;
t183 = qJD(2) ^ 2;
t262 = t183 * t179;
t261 = t183 * t181;
t158 = pkin(4) * t178 + qJ(3);
t137 = t297 * t181;
t129 = -pkin(2) * t181 + t228;
t104 = qJD(1) * t129;
t244 = t179 * t263;
t243 = t180 * t269;
t102 = pkin(4) * t264 + t137;
t237 = t155 * t252;
t229 = qJD(1) * t251;
t221 = pkin(1) * t313;
t220 = qJD(3) - t285;
t219 = -t114 - t249;
t218 = -qJD(4) + t259;
t216 = -t116 * t159 + t319;
t15 = t177 * t42 + t245;
t214 = pkin(4) * t250 - t15;
t209 = (-qJD(1) * t137 - t103) * t180;
t208 = -0.2e1 * qJD(2) * t104;
t195 = -qJ(3) * t255 - t254;
t82 = qJD(1) * t195 + t154;
t99 = t167 + t195;
t206 = pkin(7) * t183 + qJD(1) * t99 + t82;
t205 = t141 * t9 + t225;
t204 = t10 * t141 - t224;
t202 = -t177 * t59 + t296 * t56;
t198 = t177 * t24 + t231 * t56 - t250 * t59 + t27 * t296;
t73 = -pkin(4) * t236 + (-pkin(7) + t242) * t257;
t125 = pkin(7) * t230 - t173;
t126 = t163 + t220;
t135 = -t164 - t174;
t187 = -t125 * t181 + (t126 * t181 + (t135 + t164) * t179) * qJD(2);
t55 = t160 * pkin(4) + qJD(1) * t73 + t173;
t162 = -pkin(4) * t296 - pkin(5);
t157 = pkin(4) * t177 + qJ(6);
t139 = t180 * t159;
t138 = t179 * t159;
t120 = -qJ(3) * t258 + t168;
t96 = t117 * t181;
t95 = t177 * t265 - t181 * t241;
t88 = t104 * t259;
t60 = pkin(5) * t117 + qJ(6) * t116 + t158;
t46 = -pkin(5) * t95 + qJ(6) * t96 + t102;
t45 = t117 * t301 - t177 * t239 + t296 * t238;
t44 = t116 * t301 + t117 * t257;
t32 = pkin(4) * t115 + t36;
t28 = -pkin(5) * t179 - t202;
t26 = qJ(6) * t179 + t200;
t6 = -pkin(5) * t45 - qJ(6) * t44 + qJD(6) * t96 + t73;
t5 = t30 * pkin(5) + t29 * qJ(6) - qJD(6) * t197 + t55;
t4 = -pkin(5) * t255 - t302;
t3 = qJ(6) * t255 + qJD(6) * t179 + t198;
t11 = [0, 0, 0, 0.2e1 * t138, t260 * t313, t261, -t262, 0, -pkin(7) * t261 + t179 * t221, pkin(7) * t262 + t181 * t221, t187, t179 * t208 + t181 * t206, -t179 * t206 + t181 * t208, pkin(7) * t187 + t104 * t99 + t129 * t82, t115 * t191 + t188 * t265 (t178 * t160 + (-t180 * t114 + (t115 + t256) * t178 + (-t178 ^ 2 + t298) * t258) * qJD(4)) * t181 + (t115 * t180 + (t114 - 0.2e1 * t234) * t178) * t257 (-t155 - t259) * t235 + (t270 + (-qJD(1) * t176 + (t155 + t218) * t179) * t178) * qJD(2), t155 * t236 + (t178 * t229 - t160) * t179 + (t114 * t181 + (qJD(1) * t260 + t269) * t180) * qJD(2), t155 * t255 + t138, t222 * t155 + t122 * t114 + t137 * t160 + (qJD(2) * t209 + t226) * t179 + (-t155 * t276 - t179 * t51) * qJD(4) + (-t103 * t253 + t50 * qJD(2) + t281 + ((-t110 * t178 + t119) * qJD(2) - t137 * t253) * qJD(1)) * t181, -t194 * t155 - t199 * t179 - t122 * t115 + (qJD(4) * t209 - t282) * t181 + ((-qJD(1) * t276 - t51) * t181 + (t137 * t218 + t272) * t178) * qJD(2), t197 * t44 + t29 * t96, t197 * t45 - t29 * t95 + t30 * t96 - t44 * t62, t141 * t44 - t179 * t29 + (-qJD(1) * t96 + t197) * t255, t141 * t45 - t179 * t30 + (qJD(1) * t95 - t62) * t255, t141 * t255 + t138, t102 * t30 + t141 * t302 + t159 * t202 - t179 * t224 + t255 * t9 - t70 * t45 - t55 * t95 + t73 * t62, -t198 * t141 + t225 * t179 + t73 * t197 - t102 * t29 - t55 * t96 + t70 * t44 + (-qJD(1) * t200 - t10) * t255, -t141 * t4 - t179 * t2 - t23 * t45 + t30 * t46 - t5 * t95 + t6 * t62 + (-qJD(1) * t28 - t7) * t255, t1 * t95 + t197 * t4 - t2 * t96 - t26 * t30 - t28 * t29 - t3 * t62 + t44 * t7 + t45 * t8, t1 * t179 + t141 * t3 - t23 * t44 + t29 * t46 + t5 * t96 - t6 * t197 + (qJD(1) * t26 + t8) * t255, t1 * t26 + t2 * t28 + t23 * t6 + t3 * t8 + t4 * t7 + t46 * t5; 0, 0, 0, -t244, t260 * t184, 0, 0, 0, t184 * pkin(1) * t179, pkin(1) * t263 ((-t135 - t174) * t179 + (-t126 + t220) * t181) * qJD(1), -t120 * t258 + t88, 0.2e1 * t173 + (t104 * t181 + t120 * t179) * qJD(1), -qJ(3) * t125 - qJD(3) * t135 - t104 * t120 + (-t135 * t179 + (-t126 - t285) * t181) * qJD(1) * pkin(7), -t298 * t229 + (t218 * t256 - t271) * t178 (-t115 * qJD(4) - t160 + (-t115 + t256) * t259) * t180 + ((-t114 + t249) * qJD(4) + (t179 * t219 + 0.2e1 * t235) * qJD(1)) * t178, -t155 * t253 + t139 + (-t155 * t266 - t270) * qJD(1), -t237 + (t181 * t219 - t243) * qJD(1), -t155 * t258, qJ(3) * t160 + t282 - t223 * t155 - t248 * t114 + (t103 * t180 - t178 * t268) * qJD(4) + ((-qJ(3) * t253 - t50) * t181 - t303 * t180) * qJD(1), t281 + t277 * t155 + t248 * t115 + (-t180 * t268 + (-t103 - t174) * t178) * qJD(4) + ((-qJ(3) * t252 + t51) * t181 + t303 * t178) * qJD(1), t304, t116 * t30 + t117 * t29 - t197 * t286 - t287 * t62, -t197 * t258 + t216, -t233 + (-qJD(2) * t117 + t62) * t258, -t141 * t258, t55 * t117 + t158 * t30 + t286 * t70 + t278 * t62 + t309 * t141 + (-t9 + t274) * t258, -t55 * t116 - t158 * t29 + t287 * t70 + t278 * t197 + t310 * t141 + (t10 - t273) * t258, t117 * t5 + t30 * t60 + t289 * t62 + t286 * t23 + t290 * t141 + (t7 + t274) * t258, t196 * t29 - t197 * t290 + t291 * t62 - t30 * t72 - t300, t116 * t5 + t29 * t60 - t289 * t197 - t287 * t23 - t291 * t141 + (-t8 + t273) * t258, t1 * t72 - t196 * t2 + t23 * t289 - t290 * t7 - t291 * t8 + t5 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, -t175 * t184 - t183, qJD(2) * t135 + t153 + t88, 0, 0, 0, 0, 0, -t155 ^ 2 * t178 + qJD(2) * t114 + t139, -t237 - qJD(2) * t115 + (-t181 * t249 - t243) * qJD(1), 0, 0, 0, 0, 0, -qJD(2) * t62 + t216, -t316, t319 + (-t116 * t258 - t62) * qJD(2), -t117 * t30 - t286 * t62 - t304, t316, -qJD(2) * t23 + t300; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115 * t114, -t114 ^ 2 + t115 ^ 2, -t114 * t155 - t188, t189 + t271, t159, -t103 * t115 + t226 + (-qJD(4) + t155) * t51, -t103 * t114 + t155 * t50 - t199, t293, t314, t13, t305, t159, t15 * t141 - t311 + (-t115 * t62 - t141 * t250 + t159 * t296) * pkin(4) - t224, t16 * t141 + t317 + (-t115 * t197 - t141 * t231 - t159 * t177) * pkin(4) + t225, -t294 - t32 * t62 - t214 * t141 + (pkin(5) - t162) * t159 - t224, -t157 * t30 - t162 * t29 + (t214 + t8) * t197 + (-t279 + t7) * t62, t141 * t279 + t157 * t159 + t197 * t32 + t1 - t318, t1 * t157 + t162 * t2 + t214 * t7 - t23 * t32 + t279 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t293, t314, t13, t305, t159, t204 - t311, t205 + t317, -t36 * t62 + t204 + 0.2e1 * t217 - t294, pkin(5) * t29 - qJ(6) * t30 + (-t10 + t8) * t197 + (t7 - t280) * t62, t197 * t36 + 0.2e1 * t130 + 0.2e1 * t144 - t205 - t318, -pkin(5) * t2 + qJ(6) * t1 - t10 * t7 - t23 * t36 + t280 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t159 + t293, t13, -t141 ^ 2 - t299, -t141 * t8 + t2 + t294;];
tauc_reg  = t11;
