% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:52:36
% EndTime: 2019-03-09 08:52:49
% DurationCPUTime: 4.97s
% Computational Cost: add. (7125->402), mult. (18546->562), div. (0->0), fcn. (14776->10), ass. (0->206)
t221 = cos(qJ(2));
t283 = cos(pkin(10));
t245 = t283 * t221;
t203 = qJD(1) * t245;
t214 = sin(pkin(10));
t218 = sin(qJ(2));
t264 = qJD(1) * t218;
t177 = t214 * t264 - t203;
t172 = qJD(5) + t177;
t164 = qJD(6) + t172;
t216 = sin(qJ(6));
t219 = cos(qJ(6));
t193 = t214 * t221 + t283 * t218;
t180 = t193 * qJD(1);
t213 = sin(pkin(11));
t215 = cos(pkin(11));
t149 = t213 * qJD(2) + t215 * t180;
t150 = t215 * qJD(2) - t213 * t180;
t217 = sin(qJ(5));
t220 = cos(qJ(5));
t301 = -t220 * t149 - t217 * t150;
t95 = t217 * t149 - t220 * t150;
t41 = -t216 * t301 + t219 * t95;
t314 = t41 * t164;
t313 = t95 * t172;
t236 = t216 * t95 + t219 * t301;
t312 = t164 * t236;
t273 = t220 * t215;
t192 = t217 * t213 - t273;
t267 = t172 * t192;
t194 = t220 * t213 + t217 * t215;
t184 = t194 * qJD(5);
t266 = t194 * t177 + t184;
t280 = t177 * t213;
t121 = pkin(2) * t264 + t180 * pkin(3) + t177 * qJ(4);
t294 = -qJ(3) - pkin(7);
t201 = t294 * t221;
t197 = qJD(1) * t201;
t185 = t214 * t197;
t200 = t294 * t218;
t196 = qJD(1) * t200;
t140 = t283 * t196 + t185;
t77 = t213 * t121 + t215 * t140;
t61 = pkin(8) * t280 + t77;
t311 = -qJD(4) * t215 + t61;
t310 = t172 * t301;
t136 = -t216 * t192 + t219 * t194;
t290 = qJD(6) * t136 - t216 * t267 + t219 * t266;
t307 = t236 * t41;
t306 = t236 ^ 2 - t41 ^ 2;
t256 = -t221 * pkin(2) - pkin(1);
t238 = t256 * qJD(1);
t199 = qJD(3) + t238;
t112 = t177 * pkin(3) - t180 * qJ(4) + t199;
t289 = qJD(2) * pkin(2);
t190 = t196 + t289;
t246 = t283 * t197;
t133 = t214 * t190 - t246;
t129 = qJD(2) * qJ(4) + t133;
t66 = t215 * t112 - t213 * t129;
t34 = t177 * pkin(4) - t149 * pkin(8) + t66;
t67 = t213 * t112 + t215 * t129;
t48 = pkin(8) * t150 + t67;
t19 = t217 * t34 + t220 * t48;
t13 = -t95 * pkin(9) + t19;
t260 = qJD(6) * t216;
t11 = t13 * t260;
t132 = t283 * t190 + t185;
t124 = -qJD(2) * pkin(3) + qJD(4) - t132;
t88 = -pkin(4) * t150 + t124;
t39 = t95 * pkin(5) + t88;
t305 = t39 * t41 + t11;
t259 = qJD(6) * t219;
t258 = qJD(1) * qJD(2);
t253 = t218 * t258;
t166 = qJD(2) * t203 - t214 * t253;
t261 = qJD(5) * t220;
t275 = t213 * t166;
t52 = t150 * t261 + t166 * t273 + (-qJD(5) * t149 - t275) * t217;
t53 = -qJD(5) * t301 + t166 * t194;
t8 = -t216 * t53 + t219 * t52 - t95 * t259 + t260 * t301;
t304 = t8 + t314;
t179 = t193 * qJD(2);
t165 = qJD(1) * t179;
t274 = t215 * t166;
t247 = qJD(2) * t294;
t173 = t221 * qJD(3) + t218 * t247;
t156 = t173 * qJD(1);
t174 = -t218 * qJD(3) + t221 * t247;
t157 = t174 * qJD(1);
t108 = t283 * t156 + t214 * t157;
t102 = qJD(2) * qJD(4) + t108;
t205 = pkin(2) * t253;
t87 = t165 * pkin(3) - t166 * qJ(4) - t180 * qJD(4) + t205;
t36 = -t213 * t102 + t215 * t87;
t25 = t165 * pkin(4) - pkin(8) * t274 + t36;
t37 = t215 * t102 + t213 * t87;
t27 = -pkin(8) * t275 + t37;
t250 = -t217 * t27 + t220 * t25;
t225 = -qJD(5) * t19 + t250;
t2 = t165 * pkin(5) - t52 * pkin(9) + t225;
t262 = qJD(5) * t217;
t231 = t217 * t25 + t220 * t27 + t34 * t261 - t48 * t262;
t3 = -t53 * pkin(9) + t231;
t255 = t219 * t2 - t216 * t3;
t18 = -t217 * t48 + t220 * t34;
t12 = pkin(9) * t301 + t18;
t10 = t172 * pkin(5) + t12;
t286 = t219 * t13;
t5 = t216 * t10 + t286;
t303 = -qJD(6) * t5 + t39 * t236 + t255;
t224 = qJD(6) * t236 - t216 * t52 - t219 * t53;
t302 = t224 - t312;
t300 = -0.2e1 * t258;
t206 = t214 * pkin(2) + qJ(4);
t295 = pkin(8) + t206;
t188 = t295 * t213;
t189 = t295 * t215;
t268 = -t217 * t188 + t220 * t189;
t135 = t219 * t192 + t216 * t194;
t291 = -qJD(6) * t135 - t216 * t266 - t219 * t267;
t299 = -t136 * t165 - t291 * t164;
t298 = -t194 * t165 + t172 * t267;
t232 = qJD(4) * t213 + qJD(5) * t189;
t296 = pkin(8) * t215;
t76 = t215 * t121 - t213 * t140;
t47 = t180 * pkin(4) + t177 * t296 + t76;
t297 = t188 * t261 + t311 * t220 + (t232 + t47) * t217;
t175 = t177 ^ 2;
t229 = -t214 * t218 + t245;
t277 = t193 * t215;
t131 = -pkin(3) * t229 - t193 * qJ(4) + t256;
t145 = t214 * t200 - t283 * t201;
t83 = t215 * t131 - t213 * t145;
t60 = -pkin(4) * t229 - pkin(8) * t277 + t83;
t278 = t193 * t213;
t84 = t213 * t131 + t215 * t145;
t70 = -pkin(8) * t278 + t84;
t292 = t217 * t60 + t220 * t70;
t288 = t180 * t41;
t287 = t180 * t95;
t285 = t236 * t180;
t284 = t301 * t180;
t182 = t229 * qJD(2);
t257 = t218 * t289;
t101 = t179 * pkin(3) - t182 * qJ(4) - t193 * qJD(4) + t257;
t123 = t283 * t173 + t214 * t174;
t56 = t213 * t101 + t215 * t123;
t107 = t214 * t156 - t283 * t157;
t144 = -t283 * t200 - t214 * t201;
t282 = t107 * t144;
t279 = t182 * t213;
t223 = qJD(1) ^ 2;
t272 = t221 * t223;
t222 = qJD(2) ^ 2;
t271 = t222 * t218;
t270 = t222 * t221;
t265 = t218 ^ 2 - t221 ^ 2;
t139 = t214 * t196 - t246;
t104 = -pkin(4) * t280 + t139;
t254 = pkin(5) * t266 - t104;
t252 = qJD(6) * t10 + t3;
t55 = t215 * t101 - t213 * t123;
t31 = t179 * pkin(4) - t182 * t296 + t55;
t38 = -pkin(8) * t279 + t56;
t249 = -t217 * t38 + t220 * t31;
t248 = -t217 * t70 + t220 * t60;
t244 = pkin(1) * t300;
t122 = t214 * t173 - t283 * t174;
t243 = -t220 * t188 - t217 * t189;
t242 = -t135 * t165 - t164 * t290;
t209 = -t283 * pkin(2) - pkin(3);
t105 = -t194 * pkin(9) + t243;
t241 = pkin(9) * t266 - qJD(6) * t105 + t297;
t106 = -t192 * pkin(9) + t268;
t46 = t220 * t47;
t240 = t180 * pkin(5) - pkin(9) * t267 + t194 * qJD(4) + t268 * qJD(5) + qJD(6) * t106 - t217 * t61 + t46;
t239 = -t192 * t165 - t172 * t266;
t82 = pkin(4) * t275 + t107;
t90 = pkin(4) * t279 + t122;
t118 = pkin(4) * t278 + t144;
t237 = -t213 * t66 + t215 * t67;
t234 = t107 * t193 + t144 * t166;
t125 = t194 * t193;
t126 = t192 * t193;
t73 = t219 * t125 - t216 * t126;
t74 = -t216 * t125 - t219 * t126;
t230 = t217 * t31 + t220 * t38 + t60 * t261 - t70 * t262;
t198 = -t215 * pkin(4) + t209;
t228 = t124 * t182 + t234;
t227 = -t165 * t206 + t166 * t209 + (-qJD(4) + t124) * t177;
t146 = t192 * pkin(5) + t198;
t137 = t165 * t229;
t75 = t125 * pkin(5) + t118;
t72 = t182 * t194 + t261 * t277 - t262 * t278;
t71 = -t182 * t192 - t193 * t184;
t28 = t72 * pkin(5) + t90;
t22 = t53 * pkin(5) + t82;
t21 = -t125 * pkin(9) + t292;
t20 = -pkin(5) * t229 + t126 * pkin(9) + t248;
t17 = qJD(6) * t74 + t216 * t71 + t219 * t72;
t16 = -qJD(6) * t73 - t216 * t72 + t219 * t71;
t7 = -t72 * pkin(9) + t230;
t6 = t179 * pkin(5) - t71 * pkin(9) - qJD(5) * t292 + t249;
t4 = t219 * t10 - t216 * t13;
t1 = [0, 0, 0, 0.2e1 * t221 * t253, t265 * t300, t270, -t271, 0, -pkin(7) * t270 + t218 * t244, pkin(7) * t271 + t221 * t244, t108 * t229 + t122 * t180 - t123 * t177 - t132 * t182 - t133 * t179 - t145 * t165 + t234, t282 + t108 * t145 - t132 * t122 + t133 * t123 + (t199 + t238) * t257, -t122 * t150 + t83 * t165 + t55 * t177 + t66 * t179 + t213 * t228 - t229 * t36, t122 * t149 - t84 * t165 - t56 * t177 - t67 * t179 + t215 * t228 + t229 * t37, t56 * t150 - t55 * t149 + (-t166 * t83 - t182 * t66 - t193 * t36) * t215 + (-t166 * t84 - t182 * t67 - t193 * t37) * t213, t124 * t122 + t36 * t83 + t37 * t84 + t66 * t55 + t67 * t56 + t282, -t52 * t126 - t301 * t71, -t52 * t125 + t126 * t53 + t301 * t72 - t71 * t95, -t126 * t165 + t71 * t172 - t179 * t301 - t229 * t52, -t125 * t165 - t72 * t172 - t95 * t179 + t229 * t53, t172 * t179 - t137, t249 * t172 + t248 * t165 - t250 * t229 + t18 * t179 + t90 * t95 + t118 * t53 + t82 * t125 + t88 * t72 + (-t172 * t292 + t19 * t229) * qJD(5), t118 * t52 - t82 * t126 - t292 * t165 - t230 * t172 - t19 * t179 + t229 * t231 - t301 * t90 + t88 * t71, -t16 * t236 + t8 * t74, -t16 * t41 + t17 * t236 + t224 * t74 - t8 * t73, t16 * t164 + t74 * t165 - t179 * t236 - t229 * t8, -t17 * t164 - t73 * t165 - t41 * t179 - t224 * t229, t164 * t179 - t137 (-t216 * t7 + t219 * t6) * t164 + (t219 * t20 - t216 * t21) * t165 - t255 * t229 + t4 * t179 + t28 * t41 - t75 * t224 + t22 * t73 + t39 * t17 + ((-t216 * t20 - t219 * t21) * t164 + t5 * t229) * qJD(6), -t11 * t229 + t39 * t16 - t5 * t179 + t22 * t74 - t28 * t236 + t75 * t8 + (-(-qJD(6) * t21 + t6) * t164 - t20 * t165 + t2 * t229) * t216 + (-(qJD(6) * t20 + t7) * t164 - t21 * t165 + t252 * t229) * t219; 0, 0, 0, -t218 * t272, t265 * t223, 0, 0, 0, t223 * pkin(1) * t218, pkin(1) * t272 (t133 - t139) * t180 + (-t132 + t140) * t177 + (-t165 * t214 - t283 * t166) * pkin(2), t132 * t139 - t133 * t140 + (-t283 * t107 + t108 * t214 - t199 * t264) * pkin(2), -t107 * t215 + t139 * t150 - t76 * t177 - t66 * t180 + t213 * t227, t107 * t213 - t139 * t149 + t77 * t177 + t67 * t180 + t215 * t227, -t77 * t150 + t76 * t149 + (qJD(4) * t150 - t177 * t66 + t37) * t215 + (qJD(4) * t149 - t177 * t67 - t36) * t213, t107 * t209 - t124 * t139 - t66 * t76 - t67 * t77 + (-t36 * t213 + t37 * t215) * t206 + t237 * qJD(4), t52 * t194 + t267 * t301, -t52 * t192 - t194 * t53 + t266 * t301 + t267 * t95, t284 - t298, t239 + t287, -t172 * t180, t243 * t165 + t198 * t53 + t82 * t192 - t18 * t180 - t104 * t95 + t266 * t88 + (-t46 - t232 * t220 + (qJD(5) * t188 + t311) * t217) * t172, t104 * t301 - t268 * t165 + t297 * t172 + t19 * t180 + t82 * t194 + t198 * t52 - t267 * t88, t8 * t136 - t236 * t291, -t8 * t135 + t136 * t224 + t236 * t290 - t291 * t41, t285 - t299, t242 + t288, -t164 * t180 (t219 * t105 - t216 * t106) * t165 - t146 * t224 + t22 * t135 - t4 * t180 + t254 * t41 + t290 * t39 + (t216 * t241 - t219 * t240) * t164 -(t216 * t105 + t219 * t106) * t165 + t146 * t8 + t22 * t136 + t5 * t180 - t254 * t236 + t291 * t39 + (t216 * t240 + t219 * t241) * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180 ^ 2 - t175, t132 * t180 + t133 * t177 + t205, t150 * t180 + t215 * t165 - t213 * t175, -t180 * t149 - t213 * t165 - t215 * t175 (t149 * t213 + t150 * t215) * t177 + (-t213 ^ 2 - t215 ^ 2) * t166, -t124 * t180 + t177 * t237 + t37 * t213 + t36 * t215, 0, 0, 0, 0, 0, t239 - t287, t284 + t298, 0, 0, 0, 0, 0, t242 - t288, t285 + t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149 * t177 + t275, t150 * t177 + t274, -t149 ^ 2 - t150 ^ 2, t149 * t66 - t67 * t150 + t107, 0, 0, 0, 0, 0, t53 - t310, t52 - t313, 0, 0, 0, 0, 0, -t224 - t312, t8 - t314; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301 * t95, t301 ^ 2 - t95 ^ 2, t52 + t313, -t53 - t310, t165, t19 * t172 + t301 * t88 + t225, t18 * t172 + t88 * t95 - t231, -t307, t306, t304, t302, t165 -(-t216 * t12 - t286) * t164 + (-t164 * t260 + t219 * t165 + t301 * t41) * pkin(5) + t303 (-t13 * t164 - t2) * t216 + (t12 * t164 - t252) * t219 + (-t164 * t259 - t216 * t165 - t236 * t301) * pkin(5) + t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t307, t306, t304, t302, t165, t5 * t164 + t303, t4 * t164 - t216 * t2 - t219 * t252 + t305;];
tauc_reg  = t1;
