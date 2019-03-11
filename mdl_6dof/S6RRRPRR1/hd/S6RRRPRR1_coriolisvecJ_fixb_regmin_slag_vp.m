% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% tauc_reg [6x33]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRPRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:04:06
% EndTime: 2019-03-09 18:04:21
% DurationCPUTime: 5.48s
% Computational Cost: add. (11185->365), mult. (29591->494), div. (0->0), fcn. (23376->10), ass. (0->233)
t213 = cos(qJ(6));
t273 = qJD(6) * t213;
t214 = cos(qJ(5));
t215 = cos(qJ(3));
t216 = cos(qJ(2));
t277 = qJD(1) * t216;
t267 = t215 * t277;
t211 = sin(qJ(3));
t212 = sin(qJ(2));
t278 = qJD(1) * t212;
t268 = t211 * t278;
t166 = -t267 + t268;
t168 = -t211 * t277 - t215 * t278;
t207 = sin(pkin(11));
t208 = cos(pkin(11));
t239 = t208 * t166 - t207 * t168;
t138 = t214 * t239;
t142 = t207 * t166 + t208 * t168;
t210 = sin(qJ(5));
t96 = -t142 * t210 + t138;
t339 = t213 * t96;
t341 = t273 + t339;
t204 = qJD(2) + qJD(3);
t203 = qJD(5) + t204;
t305 = t203 * t96;
t272 = qJD(1) * qJD(2);
t265 = t216 * t272;
t146 = qJD(3) * t267 - t204 * t268 + t215 * t265;
t180 = t211 * t216 + t215 * t212;
t151 = t204 * t180;
t147 = t151 * qJD(1);
t103 = -t146 * t207 - t147 * t208;
t104 = t146 * t208 - t147 * t207;
t275 = qJD(5) * t210;
t43 = -qJD(5) * t138 + t210 * t103 + t214 * t104 + t142 * t275;
t23 = t43 + t305;
t209 = sin(qJ(6));
t274 = qJD(6) * t209;
t329 = -t214 * t142 - t210 * t239;
t18 = t203 * t273 + t213 * t43 - t274 * t329;
t17 = t18 * t213;
t88 = t203 * t209 + t213 * t329;
t19 = qJD(6) * t88 + t209 * t43;
t86 = -t213 * t203 + t209 * t329;
t340 = -t209 * t19 - t341 * t86 + t17;
t16 = t18 * t209;
t9 = t341 * t88 + t16;
t338 = -qJD(6) - t96;
t44 = qJD(5) * t329 - t214 * t103 + t210 * t104;
t40 = t209 * t44;
t310 = -t273 * t338 + t40;
t315 = t88 * t329;
t8 = -t338 * t339 + t310 - t315;
t139 = t142 * pkin(9);
t164 = t168 * qJ(4);
t323 = pkin(7) + pkin(8);
t189 = t323 * t216;
t183 = qJD(1) * t189;
t169 = t211 * t183;
t188 = t323 * t212;
t181 = qJD(1) * t188;
t307 = qJD(2) * pkin(2);
t175 = -t181 + t307;
t251 = t215 * t175 - t169;
t126 = t164 + t251;
t115 = t204 * pkin(3) + t126;
t173 = t215 * t183;
t238 = -t211 * t175 - t173;
t292 = t166 * qJ(4);
t127 = -t238 - t292;
t118 = t207 * t127;
t72 = t208 * t115 - t118;
t57 = pkin(4) * t204 + t139 + t72;
t320 = pkin(9) * t239;
t288 = t208 * t127;
t73 = t207 * t115 + t288;
t59 = t73 - t320;
t30 = -t210 * t59 + t214 * t57;
t28 = -pkin(5) * t203 - t30;
t318 = t28 * t96;
t313 = t329 * t96;
t25 = t329 ^ 2 - t96 ^ 2;
t56 = pkin(5) * t329 + pkin(10) * t96;
t200 = -t216 * pkin(2) - pkin(1);
t187 = t200 * qJD(1);
t152 = t166 * pkin(3) + qJD(4) + t187;
t110 = pkin(4) * t239 + t152;
t269 = qJD(2) * t323;
t245 = qJD(1) * t269;
t177 = t216 * t245;
t276 = qJD(3) * t211;
t249 = -t211 * t177 - t183 * t276;
t176 = t212 * t245;
t326 = (qJD(3) * t175 - t176) * t215;
t67 = -t147 * qJ(4) - t166 * qJD(4) + t249 + t326;
t250 = t211 * t176 - t215 * t177;
t223 = qJD(3) * t238 + t250;
t68 = -t146 * qJ(4) + t168 * qJD(4) + t223;
t38 = -t207 * t67 + t208 * t68;
t21 = -pkin(9) * t104 + t38;
t39 = t207 * t68 + t208 * t67;
t22 = pkin(9) * t103 + t39;
t3 = (qJD(5) * t57 + t22) * t214 + t210 * t21 - t59 * t275;
t221 = t110 * t96 - t3;
t316 = t86 * t329;
t248 = t211 * t181 - t173;
t128 = t248 + t292;
t280 = -t215 * t181 - t169;
t129 = t164 + t280;
t287 = t208 * t211;
t308 = pkin(2) * qJD(3);
t294 = -t208 * t128 + t129 * t207 + (-t207 * t215 - t287) * t308;
t289 = t207 * t211;
t327 = -t207 * t128 - t208 * t129 + (t208 * t215 - t289) * t308;
t297 = t329 * t203;
t24 = -t44 + t297;
t314 = t338 * t329;
t334 = qJD(6) + t338;
t31 = t210 * t57 + t214 * t59;
t4 = qJD(5) * t31 - t214 * t21 + t210 * t22;
t220 = -t110 * t329 - t4;
t29 = pkin(10) * t203 + t31;
t47 = t96 * pkin(5) - pkin(10) * t329 + t110;
t13 = t209 * t47 + t213 * t29;
t246 = t13 * t329 + t4 * t209 + t28 * t273;
t243 = t209 * t29 - t213 * t47;
t266 = t243 * t329 + t28 * t274;
t333 = t320 - t294;
t332 = -t139 + t327;
t330 = t209 * t338;
t328 = -0.2e1 * t272;
t42 = t213 * t44;
t231 = -t274 * t338 - t42;
t325 = qJD(1) * t180;
t179 = t211 * t212 - t215 * t216;
t148 = -t208 * t179 - t207 * t180;
t149 = -t207 * t179 + t208 * t180;
t106 = t148 * t210 + t149 * t214;
t242 = t214 * t148 - t149 * t210;
t286 = t215 * t188;
t140 = -t180 * qJ(4) - t211 * t189 - t286;
t237 = t211 * t188 - t215 * t189;
t141 = -t179 * qJ(4) - t237;
t90 = t208 * t140 - t141 * t207;
t70 = -pkin(9) * t149 + t90;
t91 = t207 * t140 + t208 * t141;
t71 = pkin(9) * t148 + t91;
t46 = t210 * t70 + t214 * t71;
t150 = t204 * t179;
t109 = -t150 * t208 - t151 * t207;
t182 = t212 * t269;
t184 = t216 * t269;
t229 = -qJD(3) * t286 - t215 * t182 - t211 * t184 - t189 * t276;
t80 = -qJ(4) * t151 - qJD(4) * t179 + t229;
t222 = qJD(3) * t237 + t211 * t182 - t215 * t184;
t81 = t150 * qJ(4) - t180 * qJD(4) + t222;
t52 = -t207 * t80 + t208 * t81;
t36 = -pkin(9) * t109 + t52;
t108 = t150 * t207 - t151 * t208;
t53 = t207 * t81 + t208 * t80;
t37 = pkin(9) * t108 + t53;
t45 = t210 * t71 - t214 * t70;
t5 = -qJD(5) * t45 + t210 * t36 + t214 * t37;
t50 = qJD(5) * t242 + t210 * t108 + t214 * t109;
t236 = t179 * pkin(3) + t200;
t122 = -t148 * pkin(4) + t236;
t54 = -pkin(5) * t242 - pkin(10) * t106 + t122;
t324 = (qJD(6) * t47 + t3) * t242 + t4 * t106 + (qJD(6) * t54 + t5) * t338 + t28 * t50 - t46 * t44;
t321 = pkin(3) * t207;
t319 = t168 * pkin(3);
t317 = t54 * t44;
t199 = t215 * pkin(2) + pkin(3);
t162 = -pkin(2) * t289 + t208 * t199;
t155 = pkin(4) + t162;
t163 = pkin(2) * t287 + t207 * t199;
t241 = t214 * t155 - t210 * t163;
t311 = -qJD(5) * t241 + t333 * t210 - t332 * t214;
t240 = t210 * t155 + t214 * t163;
t309 = qJD(5) * t240 + t332 * t210 + t333 * t214;
t306 = t106 * t28;
t303 = t209 * t88;
t301 = t213 * t88;
t300 = t213 * t338;
t197 = t208 * pkin(3) + pkin(4);
t233 = t214 * t197 - t210 * t321;
t74 = -t126 * t207 - t288;
t60 = t74 + t320;
t75 = t208 * t126 - t118;
t61 = t139 + t75;
t296 = -t233 * qJD(5) + t210 * t60 + t214 * t61;
t234 = t210 * t197 + t214 * t321;
t295 = t234 * qJD(5) - t210 * t61 + t214 * t60;
t291 = t168 * t166;
t290 = t187 * t168;
t218 = qJD(1) ^ 2;
t285 = t216 * t218;
t217 = qJD(2) ^ 2;
t284 = t217 * t212;
t283 = t217 * t216;
t279 = t212 ^ 2 - t216 ^ 2;
t202 = t212 * t307;
t201 = pkin(2) * t278;
t264 = -pkin(2) * t204 - t175;
t263 = pkin(3) * t147 + qJD(2) * t201;
t262 = pkin(3) * t151 + t202;
t255 = pkin(1) * t328;
t132 = pkin(10) + t240;
t114 = -t142 * pkin(4) - t319;
t49 = t114 + t56;
t254 = qJD(6) * t132 + t201 + t49;
t157 = pkin(10) + t234;
t253 = qJD(6) * t157 + t49;
t244 = -t73 * t142 - t239 * t72;
t235 = t330 * t96 - t231;
t69 = -pkin(4) * t103 + t263;
t79 = -pkin(4) * t108 + t262;
t230 = t187 * t166 - t249;
t227 = -t132 * t44 - t311 * t338 + t318;
t226 = -t157 * t44 - t296 * t338 + t318;
t156 = -pkin(5) - t233;
t131 = -pkin(5) - t241;
t130 = -t166 ^ 2 + t168 ^ 2;
t117 = (-t168 - t325) * t204;
t116 = t166 * t204 + t146;
t111 = t114 + t201;
t51 = qJD(5) * t106 - t214 * t108 + t210 * t109;
t14 = pkin(5) * t51 - pkin(10) * t50 + t79;
t11 = pkin(5) * t44 - pkin(10) * t43 + t69;
t10 = t213 * t11;
t7 = t235 + t316;
t6 = qJD(5) * t46 + t210 * t37 - t214 * t36;
t1 = t303 * t338 + t340;
t2 = [0, 0, 0, 0.2e1 * t212 * t265, t279 * t328, t283, -t284, 0, -pkin(7) * t283 + t212 * t255, pkin(7) * t284 + t216 * t255, t146 * t180 + t150 * t168, -t146 * t179 - t147 * t180 + t150 * t166 + t151 * t168, -t150 * t204, -t151 * t204, 0, t200 * t147 + t187 * t151 + t222 * t204 + (qJD(1) * t179 + t166) * t202, t200 * t146 - t187 * t150 - t229 * t204 + (-t168 + t325) * t202, t91 * t103 - t90 * t104 + t73 * t108 - t72 * t109 + t142 * t52 + t39 * t148 - t38 * t149 - t239 * t53, t152 * t262 + t236 * t263 + t38 * t90 + t39 * t91 + t72 * t52 + t73 * t53, t106 * t43 + t329 * t50, -t106 * t44 + t242 * t43 - t329 * t51 - t50 * t96, t50 * t203, -t51 * t203, 0, t110 * t51 + t122 * t44 - t203 * t6 - t242 * t69 + t79 * t96, t106 * t69 + t110 * t50 + t122 * t43 - t203 * t5 + t329 * t79, t50 * t301 + (-t274 * t88 + t17) * t106 (-t213 * t86 - t303) * t50 + (-t16 - t19 * t213 + (t209 * t86 - t301) * qJD(6)) * t106, -t106 * t231 - t18 * t242 - t50 * t300 + t88 * t51, -t106 * t310 + t19 * t242 + t330 * t50 - t86 * t51, -t242 * t44 - t338 * t51, -t10 * t242 - t243 * t51 + t45 * t19 + t6 * t86 + (-t14 * t338 + t317 + (t242 * t29 + t338 * t46 + t306) * qJD(6)) * t213 + t324 * t209, -t13 * t51 + t45 * t18 + t6 * t88 + ((-qJD(6) * t46 + t14) * t338 - t317 + (-qJD(6) * t29 + t11) * t242 - qJD(6) * t306) * t209 + t324 * t213; 0, 0, 0, -t212 * t285, t279 * t218, 0, 0, 0, t218 * pkin(1) * t212, pkin(1) * t285, -t291, t130, t116, t117, 0, -t166 * t201 + t290 - t248 * t204 + (t211 * t264 - t173) * qJD(3) + t250, t168 * t201 + t280 * t204 + (qJD(3) * t264 + t176) * t215 + t230, t163 * t103 - t162 * t104 + t142 * t294 - t239 * t327 + t244, t39 * t163 + t38 * t162 - t152 * (t201 - t319) + t327 * t73 + t294 * t72, t313, t25, t23, t24, 0, -t111 * t96 - t309 * t203 + t220, -t111 * t329 + t311 * t203 + t221, t9, t1, t8, t7, t314, t131 * t19 + t309 * t86 + (t254 * t338 - t4) * t213 + t227 * t209 + t266, t131 * t18 + t227 * t213 - t254 * t330 + t309 * t88 + t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t291, t130, t116, t117, 0, -t204 * t238 + t223 + t290, t204 * t251 + t230 - t326, t75 * t239 - t74 * t142 + (t103 * t207 - t104 * t208) * pkin(3) + t244, -t72 * t74 - t73 * t75 + (t152 * t168 + t207 * t39 + t208 * t38) * pkin(3), t313, t25, t23, t24, 0, -t114 * t96 - t295 * t203 + t220, -t114 * t329 + t296 * t203 + t221, t9, t1, t8, t7, t314, t156 * t19 + t295 * t86 + (t253 * t338 - t4) * t213 + t226 * t209 + t266, t156 * t18 + t213 * t226 - t253 * t330 + t295 * t88 + t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142 ^ 2 - t239 ^ 2, -t142 * t72 + t239 * t73 + t263, 0, 0, 0, 0, 0, t44 + t297, t43 - t305, 0, 0, 0, 0, 0, t235 - t316, -t300 * t338 - t315 - t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t313, t25, t23, t24, 0, t31 * t203 + t220, t30 * t203 + t221, t9, t330 * t88 + t340, t8, -t330 * t338 + t316 + t42, t314, -pkin(5) * t19 - t4 * t213 + (-t209 * t30 + t213 * t56) * t338 - t31 * t86 + t209 * t318 - t310 * pkin(10) + t266, -pkin(5) * t18 - (t209 * t56 + t213 * t30) * t338 - t31 * t88 + t28 * t339 + t231 * pkin(10) + t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88 * t86, -t86 ^ 2 + t88 ^ 2, -t338 * t86 + t18, -t338 * t88 - t19, t44, -t334 * t13 - t209 * t3 - t28 * t88 + t10, -t209 * t11 - t213 * t3 + t334 * t243 + t28 * t86;];
tauc_reg  = t2;
