% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:45:13
% EndTime: 2019-03-09 11:45:25
% DurationCPUTime: 4.27s
% Computational Cost: add. (10053->408), mult. (25838->525), div. (0->0), fcn. (19779->8), ass. (0->224)
t182 = sin(pkin(10));
t183 = cos(pkin(10));
t186 = sin(qJ(2));
t188 = cos(qJ(2));
t160 = -t182 * t186 + t183 * t188;
t150 = t160 * qJD(1);
t161 = t182 * t188 + t183 * t186;
t152 = t161 * qJD(1);
t185 = sin(qJ(4));
t291 = cos(qJ(4));
t297 = t291 * t150 - t185 * t152;
t106 = qJD(5) - t297;
t187 = cos(qJ(5));
t241 = qJD(5) * t187;
t184 = sin(qJ(5));
t236 = qJD(1) * qJD(2);
t224 = t188 * t236;
t225 = t186 * t236;
t141 = -t182 * t225 + t183 * t224;
t151 = t161 * qJD(2);
t196 = qJD(1) * t151;
t218 = t185 * t141 + t291 * t196;
t200 = -t185 * t150 - t291 * t152;
t312 = t200 * qJD(4);
t70 = t218 - t312;
t67 = t184 * t70;
t262 = t106 * t241 + t67;
t316 = t297 * t187;
t206 = -t106 * t316 + t262;
t235 = qJD(2) + qJD(4);
t96 = t184 * t235 - t187 * t200;
t272 = t200 * t96;
t322 = t206 + t272;
t231 = t291 * t141;
t193 = -t185 * t196 + t231;
t313 = qJD(4) * t297;
t192 = t193 + t313;
t216 = qJD(5) * t235;
t242 = qJD(5) * t184;
t42 = -t200 * t242 + (-t192 - t216) * t187;
t40 = t42 * t184;
t321 = -t40 + (t241 - t316) * t96;
t320 = pkin(5) * t200;
t93 = -t184 * t200 - t187 * t235;
t273 = t200 * t93;
t319 = qJ(6) * t200;
t318 = t106 * t200;
t317 = t200 * t297;
t239 = t297 * qJD(2);
t315 = -t239 + t193;
t240 = t200 * qJD(2);
t314 = -t240 - t218;
t311 = t200 ^ 2 - t297 ^ 2;
t283 = -qJ(3) - pkin(7);
t170 = t283 * t188;
t166 = qJD(1) * t170;
t155 = t182 * t166;
t169 = t283 * t186;
t165 = qJD(1) * t169;
t275 = qJD(2) * pkin(2);
t159 = t165 + t275;
t113 = t183 * t159 + t155;
t288 = pkin(8) * t152;
t88 = qJD(2) * pkin(3) + t113 - t288;
t253 = t183 * t166;
t114 = t182 * t159 - t253;
t289 = pkin(8) * t150;
t91 = t114 + t289;
t51 = t185 * t88 + t291 * t91;
t48 = t235 * pkin(9) + t51;
t233 = -pkin(2) * t188 - pkin(1);
t214 = t233 * qJD(1);
t167 = qJD(3) + t214;
t119 = -pkin(3) * t150 + t167;
t55 = -pkin(4) * t297 + pkin(9) * t200 + t119;
t21 = -t184 * t48 + t187 * t55;
t246 = qJD(6) - t21;
t10 = -pkin(5) * t106 + t246;
t50 = -t185 * t91 + t291 * t88;
t47 = -t235 * pkin(4) - t50;
t25 = t93 * pkin(5) - t96 * qJ(6) + t47;
t226 = qJD(4) * t291;
t243 = qJD(4) * t185;
t222 = qJD(2) * t283;
t146 = qJD(3) * t188 + t186 * t222;
t130 = t146 * qJD(1);
t147 = -qJD(3) * t186 + t188 * t222;
t131 = t147 * qJD(1);
t92 = -t130 * t182 + t183 * t131;
t78 = -pkin(8) * t141 + t92;
t95 = t183 * t130 + t182 * t131;
t79 = -pkin(8) * t196 + t95;
t16 = t185 * t79 + t91 * t226 + t88 * t243 - t291 * t78;
t191 = t184 * t192;
t258 = qJD(5) * t96;
t43 = t191 + t258;
t7 = pkin(5) * t43 + qJ(6) * t42 - qJD(6) * t96 + t16;
t227 = -t187 * t7 + t25 * t242;
t310 = -t10 * t200 + t227;
t22 = t184 * t55 + t187 * t48;
t11 = qJ(6) * t106 + t22;
t287 = t184 * t7;
t309 = t200 * t11 - t287;
t308 = t119 * t200 - t16;
t223 = -t16 * t187 + t47 * t242;
t307 = t21 * t200 + t223;
t306 = t16 * t184 - t22 * t200 + t47 * t241;
t195 = -t185 * t78 - t88 * t226 + t91 * t243 - t291 * t79;
t305 = -t119 * t297 + t195;
t212 = pkin(5) * t184 - qJ(6) * t187;
t304 = pkin(5) * t242 - qJ(6) * t241 - t184 * qJD(6) - t212 * t297;
t264 = t187 * t93;
t267 = t184 * t96;
t208 = t264 + t267;
t215 = -t42 * t187 - t96 * t242;
t281 = -t184 * t43 - t93 * t241;
t303 = t208 * t297 + t215 + t281;
t73 = -pkin(4) * t200 - pkin(9) * t297;
t302 = -0.2e1 * t236;
t174 = pkin(2) * t225;
t251 = t185 * t161;
t33 = t174 - (t150 * t226 - t152 * t243 + t231) * pkin(9) + t70 * pkin(4) + (t161 * pkin(3) + pkin(9) * t251) * t236;
t204 = t184 * t33 - t187 * t195 + t55 * t241 - t48 * t242;
t276 = qJ(6) * t70;
t2 = qJD(6) * t106 + t204 + t276;
t221 = -t184 * t195 - t187 * t33 + t48 * t241 + t55 * t242;
t292 = pkin(5) * t70;
t4 = t221 - t292;
t301 = t4 * t184 + t2 * t187;
t300 = t106 * t11 - t4;
t176 = pkin(2) * t183 + pkin(3);
t290 = pkin(2) * t182;
t296 = t291 * t176 - t185 * t290;
t136 = t296 * qJD(4);
t117 = -t165 * t182 + t253;
t98 = t117 - t289;
t118 = t183 * t165 + t155;
t99 = t118 - t288;
t58 = t185 * t98 + t291 * t99;
t299 = -t136 + t58;
t198 = t185 * t176 + t291 * t290;
t260 = t198 * qJD(4) - t185 * t99 + t291 * t98;
t69 = t187 * t70;
t298 = t106 * t242 - t69;
t122 = t183 * t169 + t170 * t182;
t104 = -pkin(8) * t161 + t122;
t123 = t182 * t169 - t183 * t170;
t105 = pkin(8) * t160 + t123;
t295 = t291 * t104 - t185 * t105;
t294 = t96 ^ 2;
t293 = t106 ^ 2;
t285 = t25 * t96;
t284 = t96 * t93;
t280 = t184 * t73 + t187 * t50;
t244 = qJD(1) * t186;
t127 = pkin(2) * t244 + pkin(3) * t152;
t59 = t127 + t73;
t279 = t184 * t59 + t187 * t58;
t116 = t185 * t160 + t291 * t161;
t132 = -pkin(3) * t160 + t233;
t199 = t291 * t160 - t251;
t64 = -pkin(4) * t199 - pkin(9) * t116 + t132;
t66 = t185 * t104 + t291 * t105;
t278 = t184 * t64 + t187 * t66;
t145 = pkin(9) + t198;
t271 = t145 * t70;
t154 = t160 * qJD(2);
t74 = t199 * qJD(4) - t185 * t151 + t291 * t154;
t269 = t184 * t74;
t268 = t184 * t93;
t266 = t187 * t43;
t265 = t187 * t74;
t263 = t187 * t96;
t261 = t260 + t304;
t259 = t304 - t51;
t257 = t106 * t184;
t256 = t297 * t184;
t250 = t187 * t136;
t190 = qJD(1) ^ 2;
t249 = t188 * t190;
t189 = qJD(2) ^ 2;
t248 = t189 * t186;
t247 = t189 * t188;
t103 = t183 * t146 + t182 * t147;
t245 = t186 ^ 2 - t188 ^ 2;
t237 = t10 * t241 + t301;
t179 = t186 * t275;
t228 = t145 * t242;
t128 = pkin(3) * t151 + t179;
t18 = t184 * t58 - t187 * t59 + t320;
t220 = t184 * t136 - t18;
t219 = pkin(1) * t302;
t102 = -t146 * t182 + t183 * t147;
t213 = t187 * pkin(5) + t184 * qJ(6);
t211 = t10 * t187 - t11 * t184;
t210 = -t297 * t47 - t271;
t209 = -t184 * t50 + t187 * t73;
t207 = t257 * t297 - t298;
t168 = -pkin(4) - t213;
t205 = t106 * t22 - t221;
t81 = -pkin(8) * t154 + t102;
t82 = -pkin(8) * t151 + t103;
t29 = qJD(4) * t295 + t185 * t81 + t291 * t82;
t75 = t116 * qJD(4) + t291 * t151 + t185 * t154;
t37 = pkin(4) * t75 - pkin(9) * t74 + t128;
t203 = t184 * t37 + t187 * t29 + t64 * t241 - t66 * t242;
t197 = t262 * pkin(9);
t194 = t211 * qJD(5) + t301;
t30 = t66 * qJD(4) + t185 * t82 - t291 * t81;
t144 = -pkin(4) - t296;
t121 = t168 - t296;
t120 = pkin(3) * t196 + t174;
t60 = pkin(5) * t96 + qJ(6) * t93;
t38 = t212 * t116 - t295;
t27 = pkin(5) * t199 + t184 * t66 - t187 * t64;
t26 = -qJ(6) * t199 + t278;
t23 = t106 * t93 - t42;
t20 = -t209 + t320;
t19 = t280 - t319;
t17 = t279 - t319;
t8 = t212 * t74 + (t213 * qJD(5) - qJD(6) * t187) * t116 + t30;
t6 = -pkin(5) * t75 + t278 * qJD(5) + t184 * t29 - t187 * t37;
t5 = qJ(6) * t75 - qJD(6) * t199 + t203;
t1 = [0, 0, 0, 0.2e1 * t186 * t224, t245 * t302, t247, -t248, 0, -pkin(7) * t247 + t186 * t219, pkin(7) * t248 + t188 * t219, -t102 * t152 + t103 * t150 - t113 * t154 - t114 * t151 - t122 * t141 - t123 * t196 + t95 * t160 - t92 * t161, t113 * t102 + t114 * t103 + t92 * t122 + t95 * t123 + (t167 + t214) * t179, t192 * t116 - t200 * t74, -t116 * t70 + t192 * t199 + t200 * t75 + t297 * t74, t74 * t235, -t75 * t235, 0, t119 * t75 - t120 * t199 - t128 * t297 + t132 * t70 - t235 * t30, t120 * t116 + t119 * t74 - t128 * t200 + t132 * t192 - t235 * t29, t116 * t215 + t74 * t263, -t208 * t74 + (t40 - t266 + (-t263 + t268) * qJD(5)) * t116, t116 * t69 + t199 * t42 + t75 * t96 + (-t116 * t242 + t265) * t106, -t116 * t67 + t199 * t43 - t75 * t93 + (-t116 * t241 - t269) * t106, t106 * t75 - t199 * t70, t221 * t199 + t21 * t75 + t30 * t93 - t295 * t43 + ((-qJD(5) * t66 + t37) * t106 + t64 * t70 + t47 * qJD(5) * t116) * t187 + ((-qJD(5) * t64 - t29) * t106 - t66 * t70 + t16 * t116 + t47 * t74) * t184, -t203 * t106 - t116 * t223 + t199 * t204 - t22 * t75 + t47 * t265 - t278 * t70 + t295 * t42 + t30 * t96, t25 * t269 - t10 * t75 - t106 * t6 + t199 * t4 - t27 * t70 + t38 * t43 + t8 * t93 + (t25 * t241 + t287) * t116, -t26 * t43 - t27 * t42 - t5 * t93 + t6 * t96 + t211 * t74 + (-t184 * t2 + t187 * t4 + (-t10 * t184 - t11 * t187) * qJD(5)) * t116, t106 * t5 + t11 * t75 + t116 * t227 - t199 * t2 - t25 * t265 + t26 * t70 + t38 * t42 - t8 * t96, t10 * t6 + t11 * t5 + t2 * t26 + t25 * t8 + t27 * t4 + t38 * t7; 0, 0, 0, -t186 * t249, t245 * t190, 0, 0, 0, t190 * pkin(1) * t186, pkin(1) * t249 (t114 + t117) * t152 + (-t118 + t113) * t150 + (-t183 * t141 - t182 * t196) * pkin(2), -t113 * t117 - t114 * t118 + (-t167 * t244 + t182 * t95 + t183 * t92) * pkin(2), t317, t311, t315, t314, 0, t127 * t297 - t235 * t260 + t308, t127 * t200 + t235 * t299 + t305, t321, t303, t322, t207 - t273, t318, t144 * t43 + t260 * t93 + t210 * t184 + ((-qJD(5) * t145 - t59) * t187 + t299 * t184) * t106 + t307, -t144 * t42 + t260 * t96 + t210 * t187 + (t228 - t250 + t279) * t106 + t306, t121 * t43 + t261 * t93 + (-t25 * t297 - t271) * t184 + (-t145 * t241 - t220) * t106 + t310, t17 * t93 - t18 * t96 + (-t10 * t297 - t136 * t93 + (-t43 + t258) * t145) * t187 + (t11 * t297 + t136 * t96 - t145 * t42 + (t145 * t93 - t11) * qJD(5)) * t184 + t237, t121 * t42 + t271 * t187 - t261 * t96 + (-t17 - t228 + (t136 - t25) * t187) * t106 + t309, t121 * t7 + t261 * t25 + (-t17 + t250) * t11 + t220 * t10 + t194 * t145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150 ^ 2 - t152 ^ 2, t113 * t152 - t114 * t150 + t174, 0, 0, 0, 0, 0, t218 - t240 - 0.2e1 * t312, t193 + t239 + 0.2e1 * t313, 0, 0, 0, 0, 0, t207 + t273, -t187 * t293 + t272 - t67, -t184 * t293 + t273 + t69 (t264 - t267) * t297 - t215 + t281, t206 - t272, t200 * t25 + t300 * t187 + (t10 * t106 + t2) * t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t317, t311, t315, t314, 0, t235 * t51 + t308, t235 * t50 + t305, t321, t303, t322, -t106 * t257 - t273 + t69, t318, -pkin(4) * t43 - t106 * t209 - t47 * t256 - t51 * t93 - t197 + t307, pkin(4) * t42 + pkin(9) * t298 + t280 * t106 - t316 * t47 - t51 * t96 + t306, t106 * t20 + t168 * t43 - t25 * t256 + t259 * t93 - t197 + t310, -t11 * t242 + t19 * t93 - t20 * t96 - t211 * t297 + (-t40 - t266 + (t263 + t268) * qJD(5)) * pkin(9) + t237, t168 * t42 - t259 * t96 + (-pkin(9) * t242 - t19) * t106 + (pkin(9) * t70 - t106 * t25) * t187 + t309, t194 * pkin(9) - t10 * t20 - t11 * t19 + t168 * t7 + t259 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t284, -t93 ^ 2 + t294, t23, t96 * t106 - t184 * t216 + t200 * t241 - t191, t70, -t47 * t96 + t205, t106 * t21 + t47 * t93 - t204, -t60 * t93 + t205 - t285 + 0.2e1 * t292, pkin(5) * t42 - qJ(6) * t43 + (t11 - t22) * t96 + (t10 - t246) * t93, 0.2e1 * t276 - t25 * t93 + t60 * t96 + (0.2e1 * qJD(6) - t21) * t106 + t204, -pkin(5) * t4 + qJ(6) * t2 - t10 * t22 + t11 * t246 - t25 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t284 - t70, t23, -t293 - t294, t285 - t300;];
tauc_reg  = t1;
