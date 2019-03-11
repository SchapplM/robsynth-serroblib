% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRR11_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:42:44
% EndTime: 2019-03-09 09:43:04
% DurationCPUTime: 7.81s
% Computational Cost: add. (6296->452), mult. (16531->646), div. (0->0), fcn. (12968->10), ass. (0->235)
t211 = sin(qJ(6));
t209 = cos(pkin(6));
t300 = qJD(1) * t209;
t192 = qJD(2) + t300;
t206 = sin(pkin(11));
t208 = cos(pkin(11));
t215 = cos(qJ(2));
t207 = sin(pkin(6));
t301 = qJD(1) * t207;
t282 = t215 * t301;
t306 = t206 * t192 + t208 * t282;
t327 = cos(qJ(5));
t131 = t327 * t306;
t137 = t192 * t208 - t206 * t282;
t212 = sin(qJ(5));
t75 = t137 * t212 + t131;
t71 = qJD(6) + t75;
t342 = t71 ^ 2;
t214 = cos(qJ(6));
t213 = sin(qJ(2));
t312 = t207 * t213;
t262 = t327 * t312;
t244 = qJD(2) * t262;
t228 = qJD(1) * t244;
t290 = qJD(1) * qJD(2);
t276 = t207 * t290;
t260 = t213 * t276;
t243 = t212 * t260;
t223 = t206 * t243 - t208 * t228;
t338 = t327 * t137 - t212 * t306;
t53 = qJD(5) * t338 + t223;
t47 = t214 * t53;
t343 = -t211 * t342 + t47;
t299 = qJD(1) * t213;
t283 = t207 * t299;
t176 = qJD(5) + t283;
t341 = t176 * t75;
t61 = -t214 * t176 + t211 * t338;
t340 = t338 * t61;
t63 = t176 * t211 + t214 * t338;
t339 = t338 * t63;
t289 = pkin(1) * t300;
t144 = pkin(8) * t283 - t215 * t289;
t291 = qJD(3) + t144;
t245 = qJD(1) * t262;
t261 = t212 * t283;
t122 = -t206 * t245 - t208 * t261;
t278 = qJD(5) * t327;
t295 = qJD(5) * t212;
t149 = -t206 * t278 - t208 * t295;
t308 = t149 + t122;
t307 = (t245 + t278) * t208 + (-t261 - t295) * t206;
t292 = pkin(3) * t283 + t291;
t230 = -t327 * t206 - t212 * t208;
t210 = -pkin(2) - qJ(4);
t323 = -pkin(9) + t210;
t162 = t323 * t206;
t163 = t323 * t208;
t231 = -t212 * t162 + t327 * t163;
t313 = t206 * t213;
t225 = (pkin(4) * t215 - pkin(9) * t313) * t207;
t188 = pkin(2) * t283;
t247 = -qJ(3) * t215 + qJ(4) * t213;
t124 = t247 * t301 + t188;
t145 = pkin(8) * t282 + t213 * t289;
t129 = pkin(3) * t282 + t145;
t66 = -t124 * t206 + t208 * t129;
t56 = qJD(1) * t225 + t66;
t310 = t208 * t213;
t268 = pkin(9) * t207 * t310;
t67 = t208 * t124 + t206 * t129;
t58 = qJD(1) * t268 + t67;
t335 = -qJD(4) * t230 - qJD(5) * t231 + t212 * t56 + t327 * t58;
t113 = t327 * t162 + t212 * t163;
t314 = t206 * t212;
t229 = t327 * t208 - t314;
t334 = -qJD(4) * t229 - qJD(5) * t113 + t212 * t58 - t327 * t56;
t284 = -pkin(4) * t208 - pkin(3);
t226 = t284 * t283;
t309 = -t226 + t291;
t326 = pkin(1) * t213;
t197 = t209 * t326;
t311 = t207 * t215;
t333 = pkin(8) * t311 + t197;
t178 = t215 * t276;
t275 = -qJ(3) * t213 - pkin(1);
t120 = (t210 * t215 + t275) * t207;
t106 = qJD(1) * t120;
t85 = t210 * t192 + t292;
t48 = -t106 * t206 + t208 * t85;
t30 = pkin(4) * t283 - pkin(9) * t137 + t48;
t49 = t208 * t106 + t206 * t85;
t32 = -t306 * pkin(9) + t49;
t12 = t212 * t30 + t327 * t32;
t221 = qJD(2) * t225;
t171 = pkin(2) * t260;
t296 = qJD(3) * t213;
t217 = (qJD(2) * t247 - qJD(4) * t215 - t296) * t207;
t79 = qJD(1) * t217 + t171;
t288 = pkin(1) * qJD(2) * t209;
t264 = qJD(1) * t288;
t138 = pkin(8) * t178 + t213 * t264;
t96 = pkin(3) * t178 - qJD(4) * t192 + t138;
t37 = -t206 * t79 + t208 * t96;
t29 = qJD(1) * t221 + t37;
t253 = qJD(2) * t268;
t38 = t206 * t96 + t208 * t79;
t33 = qJD(1) * t253 + t38;
t219 = -qJD(5) * t12 - t212 * t33 + t327 * t29;
t4 = -pkin(5) * t178 - t219;
t332 = (pkin(5) * t338 + t71 * pkin(10)) * t71 + t4;
t198 = t206 * pkin(4) + qJ(3);
t101 = -pkin(5) * t230 - pkin(10) * t229 + t198;
t11 = -t212 * t32 + t327 * t30;
t9 = -t176 * pkin(5) - t11;
t331 = -t9 * qJD(6) * t229 - t101 * t53 + (-t307 * pkin(5) + t308 * pkin(10) + qJD(6) * t113 - t309) * t71;
t52 = -qJD(5) * t131 - t137 * t295 + t206 * t228 + t208 * t243;
t19 = qJD(6) * t63 - t214 * t178 + t211 * t52;
t297 = qJD(2) * t215;
t302 = qJ(3) * qJD(2);
t179 = t192 * qJ(3);
t99 = t179 + qJD(4) + t129;
t330 = t213 * (qJD(4) - t99 + t302) - t210 * t297;
t148 = -t206 * t311 + t208 * t209;
t194 = pkin(8) * t312;
t285 = -pkin(1) * t215 - pkin(2);
t108 = pkin(3) * t312 + t194 + (-qJ(4) + t285) * t209;
t64 = t208 * t108 - t120 * t206;
t44 = pkin(4) * t312 - pkin(9) * t148 + t64;
t147 = t206 * t209 + t208 * t311;
t65 = t206 * t108 + t208 * t120;
t50 = -pkin(9) * t147 + t65;
t236 = t212 * t44 + t327 * t50;
t328 = pkin(3) + pkin(8);
t110 = -qJD(4) * t209 + (t328 * t311 + t197) * qJD(2);
t298 = qJD(2) * t213;
t281 = t207 * t298;
t183 = pkin(2) * t281;
t90 = t183 + t217;
t54 = t208 * t110 - t206 * t90;
t36 = t221 + t54;
t55 = t206 * t110 + t208 * t90;
t43 = t253 + t55;
t329 = -qJD(5) * t236 - t212 * t43 + t327 * t36;
t305 = pkin(8) * t260 - t215 * t264;
t118 = -t192 * qJD(3) + t305;
t73 = qJD(2) * t226 - t118;
t14 = pkin(5) * t53 - pkin(10) * t52 + t73;
t10 = t176 * pkin(10) + t12;
t68 = t306 * pkin(4) + t99;
t22 = t75 * pkin(5) - pkin(10) * t338 + t68;
t252 = t10 * t211 - t214 * t22;
t234 = t212 * t29 + t30 * t278 - t32 * t295 + t327 * t33;
t3 = pkin(10) * t178 + t234;
t1 = -t252 * qJD(6) + t14 * t211 + t214 * t3;
t325 = t61 * t71;
t324 = t63 * t71;
t321 = pkin(5) * t282 - t334;
t293 = qJD(6) * t214;
t294 = qJD(6) * t211;
t18 = t176 * t293 + t211 * t178 + t214 * t52 - t294 * t338;
t320 = t18 * t211;
t319 = t211 * t53;
t318 = t215 * t48;
t317 = t215 * t49;
t316 = t229 * t214;
t203 = t207 ^ 2;
t315 = t203 * qJD(1) ^ 2;
t187 = t215 * t288;
t200 = t209 * qJD(3);
t304 = t187 + t200;
t204 = t213 ^ 2;
t303 = -t215 ^ 2 + t204;
t287 = t204 * t315;
t286 = t215 * t315;
t139 = -t209 * qJ(3) - t333;
t280 = t207 * t297;
t277 = t203 * t290;
t270 = t214 * t71;
t269 = t307 * t176;
t267 = t192 + t300;
t266 = -qJD(6) * t230 + t192;
t265 = 0.2e1 * t277;
t263 = t213 * t286;
t123 = pkin(3) * t311 - t139;
t258 = -0.2e1 * pkin(1) * t277;
t257 = t308 * t176 + t229 * t178;
t256 = t145 * t192 - t138;
t6 = t10 * t214 + t211 * t22;
t16 = pkin(10) * t312 + t236;
t232 = -t327 * t147 - t212 * t148;
t86 = pkin(4) * t147 + t123;
t94 = -t212 * t147 + t327 * t148;
t25 = -pkin(5) * t232 - pkin(10) * t94 + t86;
t251 = t16 * t214 + t211 * t25;
t250 = -t16 * t211 + t214 * t25;
t249 = t38 * t206 + t37 * t208;
t248 = -t206 * t48 + t208 * t49;
t146 = t333 * qJD(2);
t246 = t138 * t209 + t146 * t192;
t242 = -pkin(8) * t281 + t187;
t241 = -t192 * t282 + t178;
t239 = -t212 * t50 + t327 * t44;
t235 = -t211 * t94 + t214 * t312;
t70 = t211 * t312 + t214 * t94;
t233 = t212 * t36 + t44 * t278 - t50 * t295 + t327 * t43;
t140 = (-pkin(2) * t215 + t275) * t207;
t92 = -t122 * t214 + t211 * t282;
t224 = t214 * t149 - t229 * t294 - t92;
t222 = (-qJ(3) * t297 - t296) * t207;
t220 = -pkin(10) * t53 + (t11 + t9) * t71;
t95 = -pkin(3) * t260 - t118;
t89 = (-pkin(8) + t284) * t281 + t304;
t2 = -t6 * qJD(6) + t214 * t14 - t211 * t3;
t218 = -t113 * t53 + t9 * t149 + t4 * t229 + (pkin(10) * t282 - qJD(6) * t101 + t335) * t71;
t143 = -qJ(3) * t282 + t188;
t142 = t285 * t209 + t194;
t133 = -t200 - t242;
t132 = qJD(1) * t140;
t130 = t183 + t222;
t125 = -t179 - t145;
t119 = -pkin(2) * t192 + t291;
t115 = qJD(1) * t222 + t171;
t109 = -t328 * t281 + t304;
t107 = t132 * t283;
t91 = -t122 * t211 - t214 * t282;
t60 = qJD(5) * t94 - t208 * t244 + t281 * t314;
t59 = qJD(5) * t232 - t230 * t281;
t24 = qJD(6) * t70 + t211 * t59 - t214 * t280;
t23 = qJD(6) * t235 + t211 * t280 + t214 * t59;
t17 = pkin(5) * t60 - pkin(10) * t59 + t89;
t15 = -pkin(5) * t312 - t239;
t8 = -pkin(5) * t280 - t329;
t7 = pkin(10) * t280 + t233;
t5 = [0, 0, 0, t213 * t215 * t265, -t303 * t265, t267 * t280, -t267 * t281, 0, t213 * t258 - t246, -t242 * t192 + t305 * t209 + t215 * t258 (-t118 * t215 + t138 * t213 + (t119 * t215 + t125 * t213) * qJD(2) + (-t133 * t215 + t146 * t213 + (t139 * t213 + t142 * t215) * qJD(2)) * qJD(1)) * t207 (-t132 * t298 + t115 * t215 + (t130 * t215 - t140 * t298) * qJD(1)) * t207 + t246, -t118 * t209 - t133 * t192 + (-t132 * t297 - t115 * t213 + (-t130 * t213 - t140 * t297) * qJD(1)) * t207, t115 * t140 + t118 * t139 + t119 * t146 + t125 * t133 + t130 * t132 + t138 * t142, t109 * t306 + t95 * t147 + ((qJD(1) * t54 + t37) * t213 + (-t99 * t310 + t318 + (-t123 * t310 + t215 * t64) * qJD(1)) * qJD(2)) * t207, t109 * t137 + t148 * t95 + ((-qJD(1) * t55 - t38) * t213 + (t99 * t313 - t317 + (t123 * t313 - t215 * t65) * qJD(1)) * qJD(2)) * t207, -t55 * t306 - t38 * t147 - t54 * t137 - t37 * t148 + ((-t206 * t64 + t208 * t65) * qJD(1) + t248) * t281, t109 * t99 + t123 * t95 + t37 * t64 + t38 * t65 + t48 * t54 + t49 * t55, t338 * t59 + t52 * t94, t232 * t52 - t338 * t60 - t53 * t94 - t59 * t75, t176 * t59 + (t213 * t52 + (qJD(1) * t94 + t338) * t297) * t207, -t176 * t60 + (-t213 * t53 + (qJD(1) * t232 - t75) * t297) * t207 (t176 * t207 + t203 * t299) * t297, t11 * t280 + t176 * t329 + t239 * t178 + t219 * t312 - t232 * t73 + t86 * t53 + t68 * t60 + t89 * t75, -t233 * t176 + t89 * t338 + t86 * t52 + t73 * t94 + t68 * t59 + (-t234 * t213 + (-t236 * qJD(1) - t12) * t297) * t207, t18 * t70 + t23 * t63, t18 * t235 - t19 * t70 - t23 * t61 - t24 * t63, -t18 * t232 + t23 * t71 + t53 * t70 + t60 * t63, t19 * t232 + t235 * t53 - t24 * t71 - t60 * t61, -t232 * t53 + t60 * t71 (-qJD(6) * t251 + t17 * t214 - t211 * t7) * t71 + t250 * t53 - t2 * t232 - t252 * t60 + t8 * t61 + t15 * t19 - t4 * t235 + t9 * t24 -(qJD(6) * t250 + t17 * t211 + t214 * t7) * t71 - t251 * t53 + t1 * t232 - t6 * t60 + t8 * t63 + t15 * t18 + t4 * t70 + t9 * t23; 0, 0, 0, -t263, t303 * t315, t241 (-qJD(2) + t192) * t283, 0, t315 * t326 + t256, pkin(1) * t286 - t144 * t192 + t305 ((-t125 - t145 - t302) * t213 + (-pkin(2) * qJD(2) - t119 + t291) * t215) * t301, -t143 * t282 + t107 - t256, t291 * t192 + (t132 * t215 + t143 * t213) * t301 - t118, -pkin(2) * t138 - qJ(3) * t118 - t119 * t145 - t125 * t291 - t132 * t143, t95 * t206 + (-t208 * t330 - t213 * t66 - t318) * t301 + t292 * t306, t208 * t95 + t292 * t137 + (t206 * t330 + t213 * t67 + t317) * t301, t67 * t306 + t66 * t137 + (qJD(4) * t137 - t283 * t49 - t37) * t208 + (qJD(4) * t306 + t283 * t48 - t38) * t206, qJ(3) * t95 - t48 * t66 - t49 * t67 + t292 * t99 + t249 * t210 + (-t206 * t49 - t208 * t48) * qJD(4), t229 * t52 + t308 * t338, -t229 * t53 + t230 * t52 - t307 * t338 - t308 * t75, -t282 * t338 + t257, -t269 + (qJD(2) * t230 + t75) * t282, -t176 * t282, -t73 * t230 + t198 * t53 + t309 * t75 + t307 * t68 + t334 * t176 + (qJD(2) * t231 - t11) * t282, t73 * t229 + t198 * t52 + t309 * t338 + t308 * t68 + t335 * t176 + (-qJD(2) * t113 + t12) * t282, t18 * t316 + t224 * t63, t61 * t92 + t63 * t91 + (-t211 * t63 - t214 * t61) * t149 - (t320 + t19 * t214 + (-t211 * t61 + t214 * t63) * qJD(6)) * t229, -t18 * t230 + t224 * t71 + t307 * t63 + t316 * t53, -t229 * t319 + t230 * t19 - t307 * t61 + (-t211 * t149 - t229 * t293 + t91) * t71, -t230 * t53 + t307 * t71, -t19 * t231 - t2 * t230 + t218 * t211 - t214 * t331 - t252 * t307 + t321 * t61 - t9 * t91, t1 * t230 - t18 * t231 + t211 * t331 + t214 * t218 - t307 * t6 + t321 * t63 - t9 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, t263, -t192 ^ 2 - t287, t125 * t192 + t107 + t138, t178 * t208 - t192 * t306 - t206 * t287, -t137 * t192 - t178 * t206 - t208 * t287 (t137 * t206 - t208 * t306) * t283, -t192 * t99 + t248 * t283 + t249, 0, 0, 0, 0, 0, -t192 * t75 + t257, t178 * t230 - t192 * t338 - t269, 0, 0, 0, 0, 0, t230 * t319 - t229 * t19 - t308 * t61 + (-t211 * t307 - t214 * t266) * t71, t230 * t47 - t229 * t18 - t308 * t63 + (t211 * t266 - t214 * t307) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-qJD(2) * t208 + t137) * t283 (qJD(2) * t206 - t306) * t283, -t137 ^ 2 - t306 ^ 2, t137 * t48 + t306 * t49 + t95, 0, 0, 0, 0, 0, t176 * t338 + t53, t52 - t341, 0, 0, 0, 0, 0, -t340 + t343, -t214 * t342 - t319 - t339; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t338 * t75, t338 ^ 2 - t75 ^ 2, t52 + t341, -t223 + (-qJD(5) + t176) * t338, t178, t12 * t176 - t338 * t68 + t219, t11 * t176 + t68 * t75 - t234, t270 * t63 + t320 (t18 - t325) * t214 + (-t19 - t324) * t211, t270 * t71 + t319 - t339, t340 + t343, -t71 * t338, -pkin(5) * t19 - t12 * t61 + t211 * t220 - t214 * t332 + t252 * t338, -pkin(5) * t18 - t12 * t63 + t211 * t332 + t214 * t220 + t338 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63 * t61, -t61 ^ 2 + t63 ^ 2, t18 + t325, -t19 + t324, t53, t6 * t71 - t63 * t9 + t2, -t252 * t71 + t61 * t9 - t1;];
tauc_reg  = t5;
