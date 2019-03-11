% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x38]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:29:48
% EndTime: 2019-03-10 03:30:05
% DurationCPUTime: 6.59s
% Computational Cost: add. (14657->415), mult. (38123->547), div. (0->0), fcn. (30100->10), ass. (0->251)
t217 = cos(qJ(6));
t291 = qJD(6) * t217;
t213 = sin(qJ(5));
t218 = cos(qJ(5));
t220 = cos(qJ(3));
t221 = cos(qJ(2));
t298 = qJD(1) * t221;
t286 = t220 * t298;
t215 = sin(qJ(3));
t216 = sin(qJ(2));
t299 = qJD(1) * t216;
t287 = t215 * t299;
t169 = -t286 + t287;
t171 = -t215 * t298 - t220 * t299;
t214 = sin(qJ(4));
t219 = cos(qJ(4));
t245 = t214 * t169 + t219 * t171;
t246 = -t219 * t169 + t214 * t171;
t99 = t213 * t245 + t218 * t246;
t391 = t217 * t99;
t397 = t291 - t391;
t165 = t171 * pkin(9);
t344 = pkin(7) + pkin(8);
t192 = t344 * t221;
t186 = qJD(1) * t192;
t172 = t215 * t186;
t191 = t344 * t216;
t184 = qJD(1) * t191;
t331 = qJD(2) * pkin(2);
t178 = -t184 + t331;
t260 = t220 * t178 - t172;
t118 = t165 + t260;
t209 = qJD(2) + qJD(3);
t107 = t209 * pkin(3) + t118;
t105 = t219 * t107;
t133 = pkin(10) * t245;
t290 = qJD(1) * qJD(2);
t283 = t221 * t290;
t147 = qJD(3) * t286 - t209 * t287 + t220 * t283;
t183 = t215 * t221 + t220 * t216;
t153 = t209 * t183;
t148 = t153 * qJD(1);
t261 = t214 * t147 + t219 * t148;
t176 = t220 * t186;
t244 = -t215 * t178 - t176;
t340 = t169 * pkin(9);
t119 = -t244 - t340;
t296 = qJD(4) * t214;
t288 = qJD(2) * t344;
t254 = qJD(1) * t288;
t179 = t216 * t254;
t180 = t221 * t254;
t259 = t215 * t179 - t220 * t180;
t229 = qJD(3) * t244 + t259;
t78 = -t147 * pkin(9) + t229;
t272 = -t119 * t296 + t214 * t78;
t297 = qJD(3) * t215;
t258 = -t215 * t180 - t186 * t297;
t352 = (qJD(3) * t178 - t179) * t220;
t77 = -t148 * pkin(9) + t258 + t352;
t16 = t219 * t77 - t261 * pkin(10) + (t105 + t133) * qJD(4) + t272;
t113 = t219 * t119;
t251 = -t107 * t214 - t113;
t275 = -t214 * t77 + t219 * t78;
t230 = qJD(4) * t251 + t275;
t295 = qJD(4) * t219;
t75 = t219 * t147 - t214 * t148 - t169 * t295 + t171 * t296;
t17 = -t75 * pkin(10) + t230;
t294 = qJD(5) * t213;
t341 = pkin(10) * t246;
t62 = -t251 + t341;
t277 = t213 * t17 - t62 * t294;
t208 = qJD(4) + t209;
t111 = t214 * t119;
t264 = t105 - t111;
t61 = t264 + t133;
t58 = pkin(4) * t208 + t61;
t3 = (qJD(5) * t58 + t16) * t218 + t277;
t205 = -t221 * pkin(2) - pkin(1);
t190 = t205 * qJD(1);
t154 = t169 * pkin(3) + t190;
t103 = -pkin(4) * t246 + t154;
t370 = t103 * t99;
t225 = -t3 - t370;
t199 = qJD(5) + t208;
t212 = sin(qJ(6));
t292 = qJD(6) * t212;
t227 = t245 * qJD(4) - t261;
t293 = qJD(5) * t218;
t33 = t213 * t227 + t218 * t75 + t245 * t294 + t246 * t293;
t361 = t213 * t246 - t218 * t245;
t21 = t199 * t291 + t217 * t33 - t292 * t361;
t19 = t21 * t212;
t86 = t199 * t212 + t217 * t361;
t396 = t397 * t86 + t19;
t304 = -qJD(6) + t99;
t34 = qJD(5) * t361 + t213 * t75 - t218 * t227;
t334 = t212 * t34 - t291 * t304;
t395 = t304 * t391 - t361 * t86 + t334;
t32 = t217 * t34;
t389 = t212 * t304;
t84 = -t217 * t199 + t212 * t361;
t394 = -t304 * t389 + t361 * t84 + t32;
t20 = t21 * t217;
t22 = qJD(6) * t86 + t212 * t33;
t393 = -t212 * t22 + t389 * t86 - t397 * t84 + t20;
t326 = t213 * t62;
t37 = t218 * t58 - t326;
t35 = -pkin(5) * t199 - t37;
t392 = t35 * t99;
t355 = t103 * t361;
t278 = t213 * t16 - t218 * t17;
t321 = t218 * t62;
t38 = t213 * t58 + t321;
t4 = qJD(5) * t38 + t278;
t224 = -t4 - t355;
t390 = t361 * t99;
t388 = t361 ^ 2 - t99 ^ 2;
t54 = pkin(5) * t361 - pkin(11) * t99;
t387 = -t199 * t99 + t33;
t384 = t304 * t361;
t257 = t215 * t184 - t176;
t121 = t257 + t340;
t301 = -t220 * t184 - t172;
t122 = t165 + t301;
t204 = t220 * pkin(2) + pkin(3);
t311 = t214 * t215;
t383 = -t204 * t295 - (-t215 * t296 + (t219 * t220 - t311) * qJD(3)) * pkin(2) + t214 * t121 + t219 * t122;
t309 = t215 * t219;
t382 = -t204 * t296 + (-t215 * t295 + (-t214 * t220 - t309) * qJD(3)) * pkin(2) - t219 * t121 + t122 * t214;
t316 = t154 * t245;
t381 = t230 + t316;
t380 = -t154 * t246 - t272;
t379 = qJD(6) + t304;
t36 = pkin(11) * t199 + t38;
t49 = -pkin(5) * t99 - pkin(11) * t361 + t103;
t14 = t212 * t49 + t217 * t36;
t378 = t14 * t361 + t4 * t212 + t35 * t291;
t377 = t199 * t361 - t34;
t252 = t212 * t36 - t217 * t49;
t376 = t252 * t361 + t35 * t292;
t342 = pkin(4) * t245;
t367 = -t133 - t383;
t366 = -t341 - t382;
t364 = (-qJD(4) * t107 - t77) * t219 + t380;
t60 = t245 ^ 2 - t246 ^ 2;
t57 = -t208 * t245 + t227;
t358 = -0.2e1 * t290;
t317 = t245 * t246;
t351 = qJD(1) * t183;
t56 = -t208 * t246 + t75;
t182 = t215 * t216 - t220 * t221;
t150 = t219 * t182 + t214 * t183;
t151 = -t214 * t182 + t219 * t183;
t101 = t218 * t150 + t213 * t151;
t102 = -t213 * t150 + t218 * t151;
t152 = t209 * t182;
t80 = -qJD(4) * t150 - t219 * t152 - t214 * t153;
t81 = qJD(4) * t151 - t214 * t152 + t219 * t153;
t45 = -qJD(5) * t101 - t213 * t81 + t218 * t80;
t308 = t220 * t191;
t131 = -t183 * pkin(9) - t215 * t192 - t308;
t243 = t215 * t191 - t220 * t192;
t132 = -t182 * pkin(9) - t243;
t73 = -pkin(10) * t151 + t131 * t219 - t132 * t214;
t250 = -t131 * t214 - t132 * t219;
t74 = -pkin(10) * t150 - t250;
t48 = t213 * t73 + t218 * t74;
t185 = t216 * t288;
t187 = t221 * t288;
t238 = -qJD(3) * t308 - t220 * t185 - t215 * t187 - t192 * t297;
t90 = -pkin(9) * t153 + t238;
t226 = qJD(3) * t243 + t215 * t185 - t220 * t187;
t91 = t152 * pkin(9) + t226;
t239 = t131 * t295 - t132 * t296 + t214 * t91 + t219 * t90;
t25 = -pkin(10) * t81 + t239;
t274 = -t214 * t90 + t219 * t91;
t26 = -t80 * pkin(10) + t250 * qJD(4) + t274;
t47 = t213 * t74 - t218 * t73;
t5 = -qJD(5) * t47 + t213 * t26 + t218 * t25;
t157 = t182 * pkin(3) + t205;
t114 = t150 * pkin(4) + t157;
t53 = pkin(5) * t101 - pkin(11) * t102 + t114;
t345 = -(qJD(6) * t49 + t3) * t101 + t4 * t102 + (qJD(6) * t53 + t5) * t304 - t48 * t34 + t35 * t45;
t339 = t171 * pkin(3);
t337 = t53 * t34;
t163 = -pkin(2) * t311 + t219 * t204 + pkin(4);
t166 = pkin(2) * t309 + t214 * t204;
t248 = t218 * t163 - t213 * t166;
t333 = -qJD(5) * t248 + t213 * t366 - t218 * t367;
t247 = t213 * t163 + t218 * t166;
t332 = qJD(5) * t247 + t213 * t367 + t218 * t366;
t329 = t102 * t35;
t325 = t217 * t86;
t203 = t219 * pkin(3) + pkin(4);
t312 = t213 * t214;
t263 = -t118 * t214 - t113;
t63 = t263 - t341;
t303 = t219 * t118 - t111;
t64 = t133 + t303;
t320 = t213 * t63 + t218 * t64 - t203 * t293 - (-t214 * t294 + (t218 * t219 - t312) * qJD(4)) * pkin(3);
t310 = t214 * t218;
t319 = -t213 * t64 + t218 * t63 + t203 * t294 + (t214 * t293 + (t213 * t219 + t310) * qJD(4)) * pkin(3);
t314 = t171 * t169;
t313 = t190 * t171;
t223 = qJD(1) ^ 2;
t307 = t221 * t223;
t222 = qJD(2) ^ 2;
t306 = t222 * t216;
t305 = t222 * t221;
t300 = t216 ^ 2 - t221 ^ 2;
t207 = t216 * t331;
t206 = pkin(2) * t299;
t285 = -pkin(4) * t199 - t58;
t282 = -pkin(2) * t209 - t178;
t281 = -pkin(3) * t208 - t107;
t123 = t148 * pkin(3) + qJD(2) * t206;
t141 = t153 * pkin(3) + t207;
t268 = pkin(1) * t358;
t127 = pkin(11) + t247;
t106 = -t339 - t342;
t51 = t106 + t54;
t267 = qJD(6) * t127 + t206 + t51;
t164 = pkin(3) * t310 + t213 * t203 + pkin(11);
t266 = qJD(6) * t164 + t51;
t201 = t213 * pkin(4) + pkin(11);
t265 = qJD(6) * t201 - t342 + t54;
t39 = t213 * t61 + t321;
t253 = pkin(4) * t294 - t39;
t67 = pkin(4) * t81 + t141;
t241 = -t292 * t304 - t32;
t240 = t190 * t169 - t258;
t236 = -t127 * t34 - t304 * t333 - t392;
t235 = -t164 * t34 - t304 * t320 - t392;
t40 = t218 * t61 - t326;
t228 = -t201 * t34 - t392 - (-pkin(4) * t293 + t40) * t304;
t55 = -pkin(4) * t227 + t123;
t202 = -t218 * pkin(4) - pkin(5);
t162 = pkin(3) * t312 - t218 * t203 - pkin(5);
t155 = t206 - t339;
t126 = -pkin(5) - t248;
t120 = -t169 ^ 2 + t171 ^ 2;
t110 = (-t171 - t351) * t209;
t109 = t169 * t209 + t147;
t104 = t106 + t206;
t46 = qJD(5) * t102 + t213 * t80 + t218 * t81;
t12 = pkin(5) * t46 - pkin(11) * t45 + t67;
t8 = t34 * pkin(5) - t33 * pkin(11) + t55;
t7 = t217 * t8;
t6 = qJD(5) * t48 + t213 * t25 - t218 * t26;
t1 = [0, 0, 0, 0.2e1 * t216 * t283, t300 * t358, t305, -t306, 0, -pkin(7) * t305 + t216 * t268, pkin(7) * t306 + t221 * t268, t147 * t183 + t152 * t171, -t147 * t182 - t148 * t183 + t152 * t169 + t153 * t171, -t152 * t209, -t153 * t209, 0, t205 * t148 + t190 * t153 + t226 * t209 + (qJD(1) * t182 + t169) * t207, t205 * t147 - t190 * t152 - t238 * t209 + (-t171 + t351) * t207, t151 * t75 - t245 * t80, -t75 * t150 + t151 * t227 + t245 * t81 + t246 * t80, t80 * t208, -t81 * t208, 0, -t141 * t246 + t157 * t261 + t123 * t150 + t154 * t81 + t274 * t208 + (-t157 * t245 + t208 * t250) * qJD(4), t123 * t151 - t141 * t245 + t154 * t80 + t157 * t75 - t208 * t239, t102 * t33 + t361 * t45, -t101 * t33 - t102 * t34 - t361 * t46 + t45 * t99, t45 * t199, -t46 * t199, 0, t101 * t55 + t103 * t46 + t114 * t34 - t199 * t6 - t67 * t99, t102 * t55 + t103 * t45 + t114 * t33 - t199 * t5 + t361 * t67, t45 * t325 + (-t292 * t86 + t20) * t102 (-t212 * t86 - t217 * t84) * t45 + (-t19 - t217 * t22 + (t212 * t84 - t325) * qJD(6)) * t102, -t217 * t304 * t45 + t21 * t101 - t102 * t241 + t86 * t46, -t22 * t101 - t102 * t334 + t389 * t45 - t84 * t46, t101 * t34 - t304 * t46, t7 * t101 - t252 * t46 + t47 * t22 + t6 * t84 + (-t12 * t304 + t337 + (-t101 * t36 + t304 * t48 + t329) * qJD(6)) * t217 + t345 * t212, -t14 * t46 + t47 * t21 + t6 * t86 + ((-qJD(6) * t48 + t12) * t304 - t337 - (-qJD(6) * t36 + t8) * t101 - qJD(6) * t329) * t212 + t345 * t217; 0, 0, 0, -t216 * t307, t300 * t223, 0, 0, 0, t223 * pkin(1) * t216, pkin(1) * t307, -t314, t120, t109, t110, 0, -t169 * t206 + t313 - t257 * t209 + (t215 * t282 - t176) * qJD(3) + t259, t171 * t206 + t301 * t209 + (qJD(3) * t282 + t179) * t220 + t240, t317, t60, t56, t57, 0, t155 * t246 + t208 * t382 + t381, t155 * t245 + t208 * t383 + t364, -t390, t388, t387, t377, 0, t104 * t99 - t199 * t332 + t224, -t104 * t361 + t199 * t333 + t225, t396, t393, t395, t394, t384, t126 * t22 + t332 * t84 + (t267 * t304 - t4) * t217 + t236 * t212 + t376, t126 * t21 + t217 * t236 - t267 * t389 + t332 * t86 + t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t314, t120, t109, t110, 0, -t209 * t244 + t229 + t313, t209 * t260 + t240 - t352, t317, t60, t56, t57, 0, -t246 * t339 + t316 - t263 * t208 + (t214 * t281 - t113) * qJD(4) + t275, -t245 * t339 + t303 * t208 + (qJD(4) * t281 - t77) * t219 + t380, -t390, t388, t387, t377, 0, t106 * t99 - t199 * t319 + t224, -t106 * t361 + t199 * t320 + t225, t396, t393, t395, t394, t384, t162 * t22 + t319 * t84 + (t266 * t304 - t4) * t217 + t235 * t212 + t376, t162 * t21 + t217 * t235 - t266 * t389 + t319 * t86 + t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t317, t60, t56, t57, 0, -t208 * t251 + t381, t208 * t264 + t364, -t390, t388, t387, t377, 0, -t99 * t342 - t355 + t39 * t199 + (t213 * t285 - t321) * qJD(5) - t278, t361 * t342 - t370 + t40 * t199 + (qJD(5) * t285 - t16) * t218 - t277, t396, t393, t395, t394, t384, t202 * t22 + t253 * t84 + (t265 * t304 - t4) * t217 + t228 * t212 + t376, t202 * t21 + t217 * t228 + t253 * t86 - t265 * t389 + t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t390, t388, t387, t377, 0, t38 * t199 + t224, t37 * t199 + t225, t396, t393, t395, t394, t384, -pkin(5) * t22 - t4 * t217 + (-t212 * t37 + t217 * t54) * t304 - t38 * t84 - t212 * t392 - t334 * pkin(11) + t376, -pkin(5) * t21 - (t212 * t54 + t217 * t37) * t304 - t38 * t86 - t35 * t391 + t241 * pkin(11) + t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86 * t84, -t84 ^ 2 + t86 ^ 2, -t304 * t84 + t21, -t304 * t86 - t22, t34, -t14 * t379 - t212 * t3 - t35 * t86 + t7, -t212 * t8 - t217 * t3 + t252 * t379 + t35 * t84;];
tauc_reg  = t1;
