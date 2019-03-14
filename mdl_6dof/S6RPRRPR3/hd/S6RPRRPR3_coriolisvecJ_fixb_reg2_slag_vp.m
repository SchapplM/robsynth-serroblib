% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPR3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:07:25
% EndTime: 2019-03-09 05:07:34
% DurationCPUTime: 4.94s
% Computational Cost: add. (6187->492), mult. (13866->624), div. (0->0), fcn. (8706->8), ass. (0->253)
t174 = sin(pkin(10)) * pkin(1) + pkin(7);
t157 = t174 * qJD(1);
t192 = cos(qJ(3));
t189 = sin(qJ(3));
t269 = t189 * qJD(2);
t119 = t192 * t157 + t269;
t104 = qJD(3) * pkin(8) + t119;
t175 = -cos(pkin(10)) * pkin(1) - pkin(2);
t133 = -pkin(3) * t192 - pkin(8) * t189 + t175;
t108 = t133 * qJD(1);
t188 = sin(qJ(4));
t191 = cos(qJ(4));
t57 = -t188 * t104 + t191 * t108;
t285 = qJD(5) - t57;
t278 = qJD(3) * t188;
t280 = qJD(1) * t189;
t140 = t191 * t280 + t278;
t279 = qJD(1) * t192;
t172 = -qJD(4) + t279;
t296 = t140 * t172;
t268 = t191 * qJD(3);
t138 = t188 * t280 - t268;
t298 = t138 * t172;
t267 = qJD(1) * qJD(3);
t246 = t191 * t267;
t275 = qJD(4) * t188;
t253 = t189 * t275;
t266 = qJD(3) * qJD(4);
t96 = qJD(1) * t253 - t191 * t266 - t192 * t246;
t274 = qJD(4) * t191;
t252 = t189 * t274;
t276 = qJD(3) * t192;
t256 = t188 * t276;
t203 = t252 + t256;
t97 = qJD(1) * t203 + t188 * t266;
t357 = t188 * (t97 - t296) + t191 * (t96 - t298);
t341 = t192 * qJD(2) - t189 * t157;
t234 = qJD(3) * pkin(3) + t341;
t207 = qJ(5) * t140 + t234;
t334 = pkin(4) + pkin(5);
t41 = -t334 * t138 + t207;
t187 = sin(qJ(6));
t190 = cos(qJ(6));
t75 = t138 * t187 + t140 * t190;
t356 = t41 * t75;
t222 = -t190 * t138 + t140 * t187;
t328 = t75 * t222;
t257 = t188 * t279;
t333 = pkin(8) - pkin(9);
t233 = pkin(3) * t189 - pkin(8) * t192;
t146 = t233 * qJD(1);
t67 = t188 * t146 + t191 * t341;
t62 = qJ(5) * t280 + t67;
t355 = pkin(9) * t257 + t333 * t275 + t62;
t160 = t333 * t191;
t288 = t191 * t192;
t263 = pkin(9) * t288;
t66 = t146 * t191 - t188 * t341;
t354 = -qJD(4) * t160 + (-t334 * t189 - t263) * qJD(1) - t66;
t351 = -pkin(9) * t140 + t285;
t106 = t119 * qJD(3);
t350 = qJ(5) * t96 - qJD(5) * t140 + t106;
t105 = t341 * qJD(3);
t149 = t233 * qJD(3);
t132 = qJD(1) * t149;
t241 = t104 * t274 + t188 * t105 + t108 * t275 - t191 * t132;
t58 = t104 * t191 + t108 * t188;
t214 = -t172 * t58 - t241;
t349 = -t222 ^ 2 + t75 ^ 2;
t166 = qJD(6) + t172;
t27 = qJD(6) * t75 - t187 * t96 - t190 * t97;
t343 = t166 * t75 - t27;
t254 = t172 * t275;
t290 = t188 * t192;
t348 = qJD(1) * (t189 * (t138 + t268) - t172 * t290) + t254;
t161 = t172 * qJD(5);
t176 = t189 * t267;
t167 = qJ(5) * t176;
t212 = t104 * t275 - t191 * t105 - t108 * t274 - t188 * t132;
t17 = -t161 + t167 - t212;
t10 = pkin(9) * t97 + t17;
t13 = pkin(9) * t96 - t334 * t176 + t241;
t33 = t334 * t172 + t351;
t164 = t172 * qJ(5);
t43 = pkin(9) * t138 + t58;
t35 = -t164 + t43;
t6 = t187 * t33 + t190 * t35;
t2 = -t6 * qJD(6) - t187 * t10 + t190 * t13;
t346 = t166 * t6 + t2;
t237 = pkin(4) * t176;
t19 = -t237 + t241;
t51 = -t164 + t58;
t345 = t172 * t51 + t19;
t271 = qJD(6) * t190;
t272 = qJD(6) * t187;
t26 = -t138 * t271 + t140 * t272 - t187 * t97 + t190 * t96;
t344 = -t166 * t222 + t26;
t342 = t97 + t296;
t340 = t188 * qJD(5) + t269;
t142 = t187 * t188 + t190 * t191;
t339 = t142 * qJD(6) - t187 * t275 - t190 * t274;
t1 = t190 * t10 + t187 * t13 + t33 * t271 - t35 * t272;
t338 = t222 * t41 - t1;
t293 = t188 * qJ(5);
t336 = -t334 * t191 - t293;
t335 = t140 ^ 2;
t332 = pkin(8) * t140;
t331 = pkin(8) * t172;
t5 = -t187 * t35 + t190 * t33;
t330 = t166 * t5;
t159 = t333 * t188;
t90 = t159 * t187 + t160 * t190;
t327 = qJD(6) * t90 - t187 * t355 + t190 * t354;
t89 = t159 * t190 - t160 * t187;
t326 = -qJD(6) * t89 + t187 * t354 + t190 * t355;
t123 = t142 * t189;
t210 = t192 * t142;
t291 = t188 * t190;
t294 = t187 * t191;
t221 = -t291 + t294;
t47 = qJD(3) * t210 + (qJD(4) - qJD(6)) * t189 * t221;
t325 = -t123 * t27 - t222 * t47;
t56 = pkin(4) * t138 - t207;
t324 = t140 * t56;
t30 = pkin(4) * t97 + t350;
t319 = t188 * t30;
t318 = t191 * t30;
t317 = t192 * t27;
t316 = t192 * t96;
t315 = t192 * t97;
t314 = t47 * t166;
t313 = t97 * t191;
t305 = qJ(5) * t191;
t213 = -t334 * t188 + t305;
t312 = qJD(4) * t213 - (qJD(1) * t213 - t157) * t192 + t340;
t311 = -t187 * t274 - t188 * t271 + t191 * t272 + t279 * t294 + (-t257 + t275) * t190;
t151 = qJ(5) * t190 - t187 * t334;
t310 = qJD(6) * t151 + t187 * t351 + t190 * t43;
t309 = qJD(1) * t210 + t339;
t255 = t192 * t268;
t289 = t189 * t191;
t79 = t97 * t289;
t308 = -t138 * t255 - t79;
t150 = -qJ(5) * t187 - t190 * t334;
t307 = qJD(6) * t150 - t187 * t43 + t190 * t351;
t228 = pkin(4) * t188 - t305;
t306 = qJD(4) * t228 - (qJD(1) * t228 + t157) * t192 - t340;
t304 = t234 * t188;
t303 = t234 * t191;
t302 = t106 * t188;
t301 = t106 * t189;
t300 = t106 * t191;
t299 = t138 * qJ(5);
t297 = t140 * t138;
t295 = t174 * t188;
t292 = t188 * t189;
t194 = qJD(3) ^ 2;
t287 = t194 * t189;
t286 = t194 * t192;
t284 = t133 * t274 + t188 * t149;
t145 = t174 * t288;
t78 = t188 * t133 + t145;
t183 = t189 ^ 2;
t283 = -t192 ^ 2 + t183;
t122 = t187 * t289 - t189 * t291;
t282 = qJD(1) * t122;
t281 = qJD(1) * t123;
t158 = qJD(1) * t175;
t277 = qJD(3) * t189;
t273 = qJD(5) * t191;
t265 = t188 * t331;
t264 = t191 * t331;
t262 = pkin(8) * t277;
t261 = pkin(8) * t268;
t113 = t140 * t252;
t260 = -t140 * t256 + t96 * t292 - t113;
t195 = qJD(1) ^ 2;
t258 = t189 * t195 * t192;
t251 = t172 * t280;
t248 = t138 ^ 2 - t335;
t244 = -pkin(4) - t295;
t162 = t183 * t246;
t242 = t162 - t316;
t240 = qJD(4) * t138 - t96;
t144 = t174 * t290;
t77 = t133 * t191 - t144;
t239 = t166 ^ 2;
t236 = t172 * t252;
t235 = t192 * t176;
t69 = -qJ(5) * t192 + t78;
t232 = qJD(4) * t145 + t133 * t275 - t149 * t191;
t46 = t187 * t255 + t189 * t339 - t190 * t256;
t231 = -t122 * t26 + t46 * t75;
t230 = t240 * pkin(8);
t229 = pkin(4) * t191 + t293;
t182 = t192 * pkin(4);
t59 = pkin(5) * t192 + t144 + t182 + (-pkin(9) * t189 - t133) * t191;
t65 = pkin(9) * t292 + t69;
t20 = -t187 * t65 + t190 * t59;
t21 = t187 * t59 + t190 * t65;
t48 = pkin(4) * t172 + t285;
t227 = -t188 * t51 + t191 * t48;
t226 = t188 * t48 + t191 * t51;
t225 = -t188 * t58 - t191 * t57;
t224 = t188 * t57 - t191 * t58;
t220 = qJD(1) * t183 - t172 * t192;
t218 = 0.2e1 * qJD(3) * t158;
t215 = t174 + t228;
t211 = t138 * t253 + t308;
t209 = t220 * t188;
t206 = -t174 + t213;
t205 = t236 + t315;
t204 = -t253 + t255;
t202 = -t172 * t57 + t212;
t201 = -t188 * t298 - t313;
t44 = (-t189 * t268 - t192 * t275) * t174 + t284;
t200 = qJD(4) * t227 + t17 * t191 + t188 * t19;
t199 = qJD(4) * t225 + t188 * t241 - t191 * t212;
t198 = t105 * t192 + t301 + (-t119 * t189 - t192 * t341) * qJD(3);
t197 = t138 * t203 + t97 * t292;
t127 = t138 * t277;
t196 = -qJD(3) * t209 + t127 + t236 - t315;
t179 = qJ(5) * t277;
t152 = -pkin(3) - t229;
t135 = t172 * t255;
t134 = pkin(3) - t336;
t125 = t140 * t277;
t107 = (-t172 - t279) * t277;
t98 = t215 * t189;
t84 = pkin(8) * t313;
t83 = pkin(4) * t140 + t299;
t76 = t206 * t189;
t70 = t182 - t77;
t64 = -pkin(4) * t280 - t66;
t61 = -t334 * t140 - t299;
t60 = -t96 - t298;
t54 = (qJD(4) * t229 - t273) * t189 + t215 * t276;
t53 = -t172 * t274 + (t172 * t288 + (-t140 + t278) * t189) * qJD(1);
t45 = t277 * t295 - t232;
t40 = -t96 * t188 - t191 * t296;
t39 = t244 * t277 + t232;
t38 = (qJD(4) * t336 + t273) * t189 + t206 * t276;
t37 = t46 * t166;
t36 = t140 * t204 - t96 * t289;
t34 = -qJD(5) * t192 + t179 + t44;
t32 = t172 * t253 + t125 - t135 + t162 + t316;
t29 = t179 + (pkin(9) * qJD(4) - qJD(3) * t174) * t289 + (-qJD(5) + (pkin(9) * qJD(3) - qJD(4) * t174) * t188) * t192 + t284;
t28 = pkin(9) * t253 + (-t263 + (-pkin(5) + t244) * t189) * qJD(3) + t232;
t24 = t26 * t192;
t16 = -t334 * t97 - t350;
t4 = -qJD(6) * t21 - t187 * t29 + t190 * t28;
t3 = qJD(6) * t20 + t187 * t28 + t190 * t29;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t235, -0.2e1 * t283 * t267, t286, -0.2e1 * t235, -t287, 0, -t174 * t286 + t189 * t218, t174 * t287 + t192 * t218, t198, t198 * t174, t36, t211 + t260, t32, t197 (-t138 * t189 - t209) * qJD(3) + t205, t107, -t172 * t45 + (t241 + (t138 * t174 - t304) * qJD(3)) * t192 + (-t234 * t274 + t302 + t174 * t97 + (qJD(1) * t77 + t57) * qJD(3)) * t189, t172 * t44 + (-t212 + (t140 * t174 - t303) * qJD(3)) * t192 + (t234 * t275 + t300 - t174 * t96 + (-qJD(1) * t78 - t58) * qJD(3)) * t189, -t138 * t44 - t140 * t45 + t77 * t96 - t78 * t97 + t225 * t276 + (qJD(4) * t224 + t188 * t212 + t191 * t241) * t189, -t212 * t78 - t241 * t77 + t44 * t58 + t45 * t57 + (-t234 * t276 + t301) * t174, t36, t32, t138 * t204 - t260 + t79, t107, t220 * t278 + t127 - t205, t197, t138 * t54 + t172 * t39 + t97 * t98 + (t278 * t56 + t19) * t192 + (t56 * t274 + t319 + (-qJD(1) * t70 - t48) * qJD(3)) * t189, -t138 * t34 + t140 * t39 - t69 * t97 - t70 * t96 + t227 * t276 + (-qJD(4) * t226 - t17 * t188 + t19 * t191) * t189, -t140 * t54 - t172 * t34 + t96 * t98 + (-t268 * t56 - t17) * t192 + (t56 * t275 - t318 + (qJD(1) * t69 + t51) * qJD(3)) * t189, t17 * t69 + t19 * t70 + t30 * t98 + t34 * t51 + t39 * t48 + t54 * t56, -t123 * t26 + t47 * t75, -t231 + t325, t314 - t24 + (-t75 - t281) * t277, t122 * t27 + t222 * t46, -t317 - t37 + (t222 + t282) * t277 (-t166 - t279) * t277, t122 * t16 + t166 * t4 + t192 * t2 + t27 * t76 + t38 * t222 + t41 * t46 + (-qJD(1) * t20 - t5) * t277, -t1 * t192 + t123 * t16 - t166 * t3 - t26 * t76 + t38 * t75 + t41 * t47 + (qJD(1) * t21 + t6) * t277, -t1 * t122 - t123 * t2 + t20 * t26 - t21 * t27 - t222 * t3 - t4 * t75 - t46 * t6 - t47 * t5, t1 * t21 + t16 * t76 + t2 * t20 + t3 * t6 + t38 * t41 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t287, -t286, 0, t105 * t189 - t106 * t192 + (t119 * t192 - t189 * t341) * qJD(3), 0, 0, 0, 0, 0, 0, t196, t172 * t204 + t125 - t242, t113 + (t140 * t276 + t189 * t240) * t188 + t308 (-qJD(3) * t224 - t106) * t192 + (-qJD(3) * t234 + t199) * t189, 0, 0, 0, 0, 0, 0, t196, t211 - t260, -t135 + (-qJD(3) * t140 + t254) * t189 + t242 (qJD(3) * t226 - t30) * t192 + (qJD(3) * t56 + t200) * t189, 0, 0, 0, 0, 0, 0, t317 - t37 + (-t222 + t282) * t277, -t314 - t24 + (-t75 + t281) * t277, t231 + t325, t1 * t123 - t122 * t2 + t16 * t192 - t277 * t41 - t46 * t5 + t47 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t258, t283 * t195, 0, t258, 0, 0, -t158 * t280, -t158 * t279, 0, 0, t40, -t357, t53, t201, t348, t251, -pkin(3) * t97 - t300 - t119 * t138 + t172 * t66 + (t264 - t304) * qJD(4) + (-t189 * t57 + (t192 * t234 - t262) * t188) * qJD(1), pkin(3) * t96 + t302 - t119 * t140 - t172 * t67 + (-t265 - t303) * qJD(4) + (t234 * t288 + (t58 - t261) * t189) * qJD(1), t138 * t67 + t140 * t66 - t84 + (t57 * t279 - t212 + (-t57 + t332) * qJD(4)) * t191 + (t230 - t214) * t188, -pkin(3) * t106 + pkin(8) * t199 + t119 * t234 - t57 * t66 - t58 * t67, t40, t53, t357, t251, -t348, t201, t152 * t97 - t172 * t64 - t318 + t306 * t138 + (t188 * t56 + t264) * qJD(4) + (t189 * t48 + (-t192 * t56 - t262) * t188) * qJD(1), t138 * t62 - t140 * t64 - t84 + (-t48 * t279 + t17 + (t48 + t332) * qJD(4)) * t191 + (t230 + t345) * t188, t152 * t96 + t172 * t62 - t319 - t306 * t140 + (-t191 * t56 + t265) * qJD(4) + (t56 * t288 + (-t51 + t261) * t189) * qJD(1), pkin(8) * t200 + t152 * t30 + t306 * t56 - t48 * t64 - t51 * t62, t221 * t26 - t309 * t75, t142 * t26 + t221 * t27 + t222 * t309 + t311 * t75, -t309 * t166 + (qJD(3) * t221 + t75) * t280, t27 * t142 - t222 * t311, t311 * t166 + (qJD(3) * t142 - t222) * t280, t166 * t280, t134 * t27 + t142 * t16 + t312 * t222 - t311 * t41 - t327 * t166 + (-qJD(3) * t89 + t5) * t280, -t134 * t26 - t221 * t16 + t312 * t75 - t309 * t41 + t326 * t166 + (qJD(3) * t90 - t6) * t280, -t1 * t142 + t2 * t221 + t222 * t326 + t26 * t89 - t27 * t90 + t309 * t5 + t311 * t6 + t327 * t75, t1 * t90 + t134 * t16 + t2 * t89 + t312 * t41 - t326 * t6 - t327 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t297, -t248, t60, -t297, -t342, t176, t140 * t234 + t214, -t138 * t234 + t202, 0, 0, t297, t60, t248, t176, t342, -t297, -t138 * t83 + t214 + 0.2e1 * t237 - t324, pkin(4) * t96 - qJ(5) * t97 + (t51 - t58) * t140 + (t48 - t285) * t138, -t138 * t56 + t140 * t83 - 0.2e1 * t161 + 0.2e1 * t167 - t202, -pkin(4) * t19 + qJ(5) * t17 + t285 * t51 - t48 * t58 - t56 * t83, -t328, -t349, t344, t328, -t343, t176, -t150 * t176 - t166 * t310 - t222 * t61 - t2 + t356, t151 * t176 - t166 * t307 - t61 * t75 - t338, t150 * t26 - t151 * t27 + (-t307 + t5) * t222 + (t310 - t6) * t75, t1 * t151 + t150 * t2 + t307 * t6 - t310 * t5 - t41 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176 + t297, t60, -t172 ^ 2 - t335, t324 + t345, 0, 0, 0, 0, 0, 0, -t140 * t222 - t176 * t190 - t187 * t239, -t140 * t75 + t176 * t187 - t190 * t239, t187 * t343 + t190 * t344, -t41 * t140 + t346 * t190 + (t1 - t330) * t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t328, t349, -t344, -t328, t343, -t176, t346 - t356, t330 + t338, 0, 0;];
tauc_reg  = t7;