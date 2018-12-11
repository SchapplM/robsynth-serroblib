% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPRRR14_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_inertiaDJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:32:34
% EndTime: 2018-12-10 18:33:07
% DurationCPUTime: 8.99s
% Computational Cost: add. (21800->540), mult. (68696->1006), div. (0->0), fcn. (74396->16), ass. (0->243)
t197 = sin(pkin(8));
t201 = cos(pkin(8));
t198 = sin(pkin(7));
t199 = sin(pkin(6));
t207 = sin(qJ(2));
t307 = qJD(2) * t207;
t275 = t199 * t307;
t260 = t198 * t275;
t196 = sin(pkin(14));
t200 = cos(pkin(14));
t202 = cos(pkin(7));
t259 = t202 * t275;
t210 = cos(qJ(2));
t306 = qJD(2) * t210;
t274 = t199 * t306;
t310 = t196 * t274 + t200 * t259;
t114 = -t197 * t310 - t201 * t260;
t206 = sin(qJ(4));
t334 = cos(qJ(4));
t283 = t201 * t334;
t284 = t197 * t334;
t217 = t197 * t260 - t310 * t201;
t203 = cos(pkin(6));
t333 = pkin(1) * t203;
t183 = t306 * t333;
t269 = -qJ(3) * t202 - pkin(10);
t250 = t269 * t207;
t303 = qJD(3) * t202;
t305 = qJD(3) * t198;
t108 = t203 * t305 + t183 + (qJD(2) * t250 + t210 * t303) * t199;
t290 = t207 * t333;
t302 = qJD(3) * t207;
t314 = t199 * t210;
t115 = -t199 * t202 * t302 + (t269 * t314 - t290) * qJD(2);
t322 = qJ(3) * t198;
t134 = (-t198 * t302 + (pkin(2) * t207 - t210 * t322) * qJD(2)) * t199;
t320 = t196 * t202;
t321 = t196 * t198;
t73 = t200 * t108 + t115 * t320 + t134 * t321;
t56 = pkin(11) * t217 + t73;
t152 = -t196 * t259 + t200 * t274;
t331 = pkin(11) * t201;
t313 = t200 * t202;
t318 = t198 * t200;
t72 = -t108 * t196 + t115 * t313 + t134 * t318;
t58 = pkin(3) * t260 - t152 * t331 + t72;
t229 = t198 * t203 + t202 * t314;
t315 = t199 * t207;
t135 = t196 * t315 - t200 * t229;
t267 = -t198 * t314 + t202 * t203;
t242 = t197 * t267;
t220 = -t201 * t135 + t242;
t224 = -pkin(10) * t314 - t290;
t132 = qJ(3) * t229 - t224;
t142 = (pkin(1) * t210 + pkin(2)) * t203 + t199 * t250;
t154 = (-pkin(2) * t210 - t207 * t322 - pkin(1)) * t199;
t88 = t200 * t132 + t142 * t320 + t154 * t321;
t65 = pkin(11) * t220 + t88;
t64 = t334 * t65;
t136 = t196 * t229 + t200 * t315;
t87 = -t196 * t132 + t142 * t313 + t154 * t318;
t68 = pkin(3) * t267 - t136 * t331 + t87;
t332 = pkin(11) * t197;
t98 = -t198 * t115 + t202 * t134;
t76 = t310 * pkin(3) - t152 * t332 + t98;
t102 = -t142 * t198 + t202 * t154;
t82 = pkin(3) * t135 - t136 * t332 + t102;
t20 = -((t197 * t82 + t201 * t68) * t206 + t64) * qJD(4) - t206 * t56 + t58 * t283 + t76 * t284;
t16 = t114 * pkin(4) - t20;
t205 = sin(qJ(5));
t209 = cos(qJ(5));
t92 = t135 * t283 + t206 * t136 - t334 * t242;
t59 = -t92 * qJD(4) + t334 * t152 + t217 * t206;
t104 = -t135 * t197 - t201 * t267;
t93 = t334 * t136 + t220 * t206;
t66 = t209 * t104 + t205 * t93;
t39 = -qJD(5) * t66 - t205 * t114 + t209 * t59;
t67 = -t104 * t205 + t209 * t93;
t40 = qJD(5) * t67 + t209 * t114 + t205 * t59;
t341 = t40 * pkin(5) - t39 * pkin(13) + t16;
t340 = -0.4e1 * t205;
t208 = cos(qJ(6));
t194 = t208 ^ 2;
t204 = sin(qJ(6));
t309 = t204 ^ 2 - t194;
t266 = t309 * qJD(6);
t187 = pkin(2) * t313;
t140 = t202 * pkin(3) + t187 + (-qJ(3) - t331) * t321;
t153 = (-pkin(3) * t200 - t196 * t332 - pkin(2)) * t198;
t101 = -t140 * t197 + t201 * t153;
t261 = t200 * t283;
t262 = t202 * t284;
t316 = t198 * t206;
t137 = t196 * t316 - t198 * t261 - t262;
t282 = t334 * t196;
t312 = t201 * t206;
t319 = t197 * t206;
t138 = t202 * t319 + (t200 * t312 + t282) * t198;
t83 = pkin(4) * t137 - pkin(12) * t138 + t101;
t162 = t197 * t318 - t201 * t202;
t165 = pkin(2) * t320 + qJ(3) * t318;
t131 = (t197 * t202 + t201 * t318) * pkin(11) + t165;
t123 = t334 * t131;
t285 = t140 * t312 + t153 * t319 + t123;
t86 = -pkin(12) * t162 + t285;
t327 = t205 * t83 + t209 * t86;
t271 = qJD(4) * t334;
t254 = t197 * t271;
t255 = t201 * t271;
t256 = t334 * t305;
t280 = t196 * t305;
t78 = -t140 * t255 - t153 * t254 - t200 * t256 + (qJD(4) * t131 + t201 * t280) * t206;
t127 = (t262 + (-t196 * t206 + t261) * t198) * qJD(4);
t128 = t138 * qJD(4);
t258 = t197 * t280;
t97 = pkin(4) * t128 - pkin(12) * t127 + t258;
t38 = -qJD(5) * t327 + t205 * t78 + t209 * t97;
t301 = qJD(4) * t206;
t19 = -t82 * t254 - t68 * t255 + t65 * t301 - t58 * t312 - t76 * t319 - t334 * t56;
t15 = -pkin(12) * t114 - t19;
t42 = -t197 * t58 + t201 * t76;
t60 = t93 * qJD(4) + t206 * t152 - t217 * t334;
t26 = pkin(4) * t60 - pkin(12) * t59 + t42;
t291 = t68 * t312 + t82 * t319 + t64;
t33 = -pkin(12) * t104 + t291;
t43 = -t197 * t68 + t201 * t82;
t36 = pkin(4) * t92 - pkin(12) * t93 + t43;
t328 = t205 * t36 + t209 * t33;
t6 = -qJD(5) * t328 - t205 * t15 + t209 * t26;
t337 = 0.2e1 * t128;
t330 = pkin(12) * t204;
t329 = pkin(12) * t205;
t46 = t204 * t67 - t208 * t92;
t23 = -qJD(6) * t46 + t204 * t60 + t208 * t39;
t326 = t23 * t204;
t325 = t23 * t208;
t105 = t138 * t205 + t209 * t162;
t90 = -qJD(5) * t105 + t209 * t127;
t106 = t138 * t209 - t162 * t205;
t94 = t106 * t204 - t208 * t137;
t52 = -qJD(6) * t94 + t204 * t128 + t208 * t90;
t324 = t52 * t204;
t323 = t52 * t208;
t311 = t208 * t209;
t193 = t205 ^ 2;
t308 = -t209 ^ 2 + t193;
t304 = qJD(3) * t200;
t300 = qJD(5) * t204;
t299 = qJD(5) * t205;
t298 = qJD(5) * t208;
t297 = qJD(5) * t209;
t296 = qJD(6) * t204;
t295 = qJD(6) * t208;
t294 = qJD(6) * t209;
t293 = -0.2e1 * pkin(4) * qJD(5);
t292 = -0.2e1 * pkin(5) * qJD(6);
t289 = t209 * t330;
t288 = pkin(12) * t311;
t287 = t205 * t319;
t191 = t199 ^ 2;
t281 = t191 * t306;
t279 = t208 * t297;
t278 = t204 * t294;
t277 = t208 * t294;
t276 = t197 * t301;
t273 = t204 * t295;
t272 = t205 * t297;
t270 = qJD(5) * t334;
t265 = t308 * qJD(5);
t264 = -0.2e1 * t198 * t303;
t263 = 0.2e1 * t272;
t257 = t204 * t279;
t91 = qJD(5) * t106 + t205 * t127;
t253 = t91 * pkin(5) - t90 * pkin(13);
t252 = -pkin(5) * t209 - pkin(13) * t205;
t251 = pkin(5) * t205 - pkin(13) * t209;
t13 = pkin(13) * t92 + t328;
t216 = -t206 * t65 + t68 * t283 + t82 * t284;
t32 = t104 * pkin(4) - t216;
t22 = t66 * pkin(5) - t67 * pkin(13) + t32;
t8 = t13 * t208 + t204 * t22;
t45 = pkin(13) * t137 + t327;
t215 = -t206 * t131 + t140 * t283 + t153 * t284;
t85 = t162 * pkin(4) - t215;
t51 = t105 * pkin(5) - t106 * pkin(13) + t85;
t28 = t204 * t51 + t208 * t45;
t47 = t204 * t92 + t208 * t67;
t247 = -t204 * t47 - t208 * t46;
t95 = t106 * t208 + t137 * t204;
t246 = -t204 * t95 - t208 * t94;
t17 = -t205 * t33 + t209 * t36;
t48 = -t205 * t86 + t209 * t83;
t179 = -pkin(4) + t252;
t157 = t179 * t204 + t288;
t12 = -pkin(5) * t92 - t17;
t4 = -t60 * pkin(5) - t6;
t240 = t12 * t295 + t4 * t204;
t239 = t12 * t296 - t4 * t208;
t238 = t205 * t60 + t92 * t297;
t237 = -t209 * t60 + t92 * t299;
t31 = -t128 * pkin(5) - t38;
t44 = -pkin(5) * t137 - t48;
t236 = t31 * t204 + t44 * t295;
t235 = -t31 * t208 + t44 * t296;
t234 = t204 * t40 + t66 * t295;
t233 = -t208 * t40 + t66 * t296;
t232 = t251 * t204;
t231 = t105 * t295 + t204 * t91;
t230 = t105 * t296 - t208 * t91;
t167 = t201 * t205 + t209 * t319;
t5 = -t209 * t15 - t205 * t26 - t36 * t297 + t33 * t299;
t37 = -t205 * t97 + t209 * t78 - t83 * t297 + t86 * t299;
t228 = t128 * t205 + t137 * t297;
t227 = -t128 * t209 + t137 * t299;
t147 = qJD(5) * t167 + t205 * t254;
t166 = -t209 * t201 + t287;
t226 = t147 * t204 + t166 * t295;
t225 = -t147 * t208 + t166 * t296;
t223 = -t208 * t167 + t204 * t284;
t222 = t204 * t167 + t208 * t284;
t221 = t205 * t298 + t278;
t219 = -pkin(13) * t60 + t5;
t218 = -pkin(13) * t128 + t37;
t212 = (t123 + (t140 * t201 + t153 * t197) * t206) * qJD(4);
t79 = (t200 * t206 + t201 * t282) * t305 + t212;
t164 = -qJ(3) * t321 + t187;
t161 = t224 * qJD(2);
t160 = pkin(10) * t275 - t183;
t156 = t179 * t208 - t289;
t146 = qJD(5) * t287 - t201 * t297 - t209 * t254;
t110 = -t157 * qJD(6) + (t204 * t329 + t208 * t251) * qJD(5);
t109 = pkin(12) * t221 - qJD(5) * t232 - t179 * t295;
t100 = qJD(6) * t223 + t204 * t146 + t208 * t276;
t99 = qJD(6) * t222 + t208 * t146 - t204 * t276;
t53 = qJD(6) * t95 - t208 * t128 + t204 * t90;
t27 = -t204 * t45 + t208 * t51;
t24 = qJD(6) * t47 + t204 * t39 - t208 * t60;
t11 = -t28 * qJD(6) + t204 * t218 + (t196 * t201 * t256 + t304 * t316 + t212 + t253) * t208;
t10 = t45 * t296 + t208 * t218 - t204 * (t253 + t79) - t51 * t295;
t7 = -t13 * t204 + t208 * t22;
t2 = -t8 * qJD(6) + t204 * t219 + t208 * t341;
t1 = t13 * t296 - t204 * t341 + t208 * t219 - t22 * t295;
t3 = [0, 0, 0, 0.2e1 * t207 * t281, 0.2e1 * (-t207 ^ 2 + t210 ^ 2) * t191 * qJD(2), 0.2e1 * t203 * t274, -0.2e1 * t203 * t275, 0, -0.2e1 * pkin(1) * t191 * t307 + 0.2e1 * t161 * t203, -0.2e1 * pkin(1) * t281 + 0.2e1 * t160 * t203, 0.2e1 * t102 * t310 + 0.2e1 * t98 * t135 + 0.2e1 * t87 * t260 + 0.2e1 * t72 * t267, 0.2e1 * t102 * t152 + 0.2e1 * t98 * t136 - 0.2e1 * t260 * t88 - 0.2e1 * t267 * t73, -0.2e1 * t73 * t135 - 0.2e1 * t72 * t136 - 0.2e1 * t87 * t152 - 0.2e1 * t88 * t310, 0.2e1 * t102 * t98 + 0.2e1 * t72 * t87 + 0.2e1 * t73 * t88, 0.2e1 * t93 * t59, -0.2e1 * t59 * t92 - 0.2e1 * t60 * t93, -0.2e1 * t104 * t59 - 0.2e1 * t114 * t93, 0.2e1 * t104 * t60 + 0.2e1 * t114 * t92, 0.2e1 * t104 * t114, -0.2e1 * t20 * t104 - 0.2e1 * t114 * t216 + 0.2e1 * t42 * t92 + 0.2e1 * t43 * t60, -0.2e1 * t19 * t104 + 0.2e1 * t114 * t291 + 0.2e1 * t42 * t93 + 0.2e1 * t43 * t59, 0.2e1 * t67 * t39, -0.2e1 * t39 * t66 - 0.2e1 * t40 * t67, 0.2e1 * t39 * t92 + 0.2e1 * t60 * t67, -0.2e1 * t40 * t92 - 0.2e1 * t60 * t66, 0.2e1 * t92 * t60, 0.2e1 * t16 * t66 + 0.2e1 * t17 * t60 + 0.2e1 * t32 * t40 + 0.2e1 * t6 * t92, 0.2e1 * t16 * t67 + 0.2e1 * t32 * t39 - 0.2e1 * t328 * t60 + 0.2e1 * t5 * t92, 0.2e1 * t47 * t23, -0.2e1 * t23 * t46 - 0.2e1 * t24 * t47, 0.2e1 * t23 * t66 + 0.2e1 * t40 * t47, -0.2e1 * t24 * t66 - 0.2e1 * t40 * t46, 0.2e1 * t66 * t40, 0.2e1 * t12 * t24 + 0.2e1 * t2 * t66 + 0.2e1 * t4 * t46 + 0.2e1 * t40 * t7, 0.2e1 * t1 * t66 + 0.2e1 * t12 * t23 + 0.2e1 * t4 * t47 - 0.2e1 * t40 * t8; 0, 0, 0, 0, 0, t274, -t275, 0, t161, t160, t72 * t202 + (-t196 * qJD(3) * t267 - pkin(2) * t310 + t164 * t275 - t98 * t200) * t198, -t73 * t202 + (-pkin(2) * t152 - t165 * t275 + t98 * t196 - t267 * t304) * t198, -t164 * t152 - t165 * t310 + (-t72 * t196 + t73 * t200 + (-t135 * t200 + t136 * t196) * qJD(3)) * t198, t72 * t164 + t73 * t165 + (-pkin(2) * t98 + (-t196 * t87 + t200 * t88) * qJD(3)) * t198, t127 * t93 + t138 * t59, -t127 * t92 - t128 * t93 - t137 * t59 - t138 * t60, -t104 * t127 - t114 * t138 - t162 * t59, t104 * t128 + t114 * t137 + t162 * t60, t114 * t162, t101 * t60 + t79 * t104 - t114 * t215 + t43 * t128 + t42 * t137 - t20 * t162 + t258 * t92, t101 * t59 - t78 * t104 + t114 * t285 + t43 * t127 + t42 * t138 - t19 * t162 + t258 * t93, t106 * t39 + t67 * t90, -t105 * t39 - t106 * t40 - t66 * t90 - t67 * t91, t106 * t60 + t128 * t67 + t137 * t39 + t90 * t92, -t105 * t60 - t128 * t66 - t137 * t40 - t91 * t92, t128 * t92 + t137 * t60, t105 * t16 + t128 * t17 + t137 * t6 + t32 * t91 + t38 * t92 + t40 * t85 + t48 * t60 + t66 * t79, t106 * t16 - t128 * t328 + t137 * t5 + t32 * t90 - t327 * t60 + t37 * t92 + t39 * t85 + t67 * t79, t23 * t95 + t47 * t52, -t23 * t94 - t24 * t95 - t46 * t52 - t47 * t53, t105 * t23 + t40 * t95 + t47 * t91 + t52 * t66, -t105 * t24 - t40 * t94 - t46 * t91 - t53 * t66, t105 * t40 + t66 * t91, t105 * t2 + t11 * t66 + t12 * t53 + t24 * t44 + t27 * t40 + t31 * t46 + t4 * t94 + t7 * t91, t1 * t105 + t10 * t66 + t12 * t52 + t23 * t44 - t28 * t40 + t31 * t47 + t4 * t95 - t8 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196 * t264, t200 * t264, 0.2e1 * (t196 ^ 2 + t200 ^ 2) * t198 ^ 2 * qJD(3), 0.2e1 * (-t164 * t196 + t165 * t200) * t305, 0.2e1 * t138 * t127, -0.2e1 * t127 * t137 - 0.2e1 * t128 * t138, -0.2e1 * t162 * t127, t162 * t337, 0, 0.2e1 * t101 * t128 + 0.2e1 * t137 * t258 + 0.2e1 * t162 * t79, 0.2e1 * t101 * t127 + 0.2e1 * t138 * t258 - 0.2e1 * t162 * t78, 0.2e1 * t106 * t90, -0.2e1 * t105 * t90 - 0.2e1 * t106 * t91, 0.2e1 * t106 * t128 + 0.2e1 * t137 * t90, -0.2e1 * t105 * t128 - 0.2e1 * t137 * t91, t137 * t337, 0.2e1 * t105 * t79 + 0.2e1 * t128 * t48 + 0.2e1 * t137 * t38 + 0.2e1 * t85 * t91, 0.2e1 * t106 * t79 - 0.2e1 * t128 * t327 + 0.2e1 * t137 * t37 + 0.2e1 * t85 * t90, 0.2e1 * t95 * t52, -0.2e1 * t52 * t94 - 0.2e1 * t53 * t95, 0.2e1 * t105 * t52 + 0.2e1 * t91 * t95, -0.2e1 * t105 * t53 - 0.2e1 * t91 * t94, 0.2e1 * t105 * t91, 0.2e1 * t105 * t11 + 0.2e1 * t27 * t91 + 0.2e1 * t31 * t94 + 0.2e1 * t44 * t53, 0.2e1 * t10 * t105 - 0.2e1 * t28 * t91 + 0.2e1 * t31 * t95 + 0.2e1 * t44 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t310, t152, 0, t98, 0, 0, 0, 0, 0, t201 * t60 + (t104 * t301 - t334 * t114) * t197, t201 * t59 + (t104 * t271 + t114 * t206) * t197, 0, 0, 0, 0, 0, -t147 * t92 - t166 * t60 + (t66 * t301 - t334 * t40) * t197, t146 * t92 - t167 * t60 + (t67 * t301 - t334 * t39) * t197, 0, 0, 0, 0, 0, t100 * t66 + t147 * t46 + t166 * t24 - t222 * t40, t147 * t47 + t166 * t23 + t223 * t40 + t66 * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128 * t201 + t162 * t276, t201 * t127 + t162 * t254, 0, 0, 0, 0, 0, -t166 * t128 - t147 * t137 + (t105 * t301 - t334 * t91) * t197, -t167 * t128 + t146 * t137 + (t106 * t301 - t334 * t90) * t197, 0, 0, 0, 0, 0, t100 * t105 + t147 * t94 + t166 * t53 - t222 * t91, t105 * t99 + t147 * t95 + t166 * t52 + t223 * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t60, -t114, t20, t19, t205 * t39 + t297 * t67, -t40 * t205 + t209 * t39 + (-t205 * t67 - t209 * t66) * qJD(5), t238, -t237, 0, -pkin(4) * t40 - pkin(12) * t238 - t16 * t209 + t299 * t32, -pkin(4) * t39 + pkin(12) * t237 + t16 * t205 + t297 * t32, t47 * t279 + (-t296 * t47 + t325) * t205, t247 * t297 + (-t326 - t208 * t24 + (t204 * t46 - t208 * t47) * qJD(6)) * t205 (t298 * t66 - t23) * t209 + (qJD(5) * t47 - t233) * t205 (-t300 * t66 + t24) * t209 + (-qJD(5) * t46 - t234) * t205, -t209 * t40 + t299 * t66, t110 * t66 + t156 * t40 + (-t2 + (pkin(12) * t46 + t12 * t204) * qJD(5)) * t209 + (pkin(12) * t24 + qJD(5) * t7 + t240) * t205, t109 * t66 - t157 * t40 + (-t1 + (pkin(12) * t47 + t12 * t208) * qJD(5)) * t209 + (pkin(12) * t23 - qJD(5) * t8 - t239) * t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, -t128, 0, -t79, t78, t106 * t297 + t205 * t90, -t91 * t205 + t209 * t90 + (-t105 * t209 - t106 * t205) * qJD(5), t228, -t227, 0, -pkin(4) * t91 - pkin(12) * t228 - t79 * t209 + t299 * t85, -pkin(4) * t90 + pkin(12) * t227 + t79 * t205 + t297 * t85, t95 * t279 + (-t296 * t95 + t323) * t205, t246 * t297 + (-t324 - t208 * t53 + (t204 * t94 - t208 * t95) * qJD(6)) * t205 (t105 * t298 - t52) * t209 + (qJD(5) * t95 - t230) * t205 (-t105 * t300 + t53) * t209 + (-qJD(5) * t94 - t231) * t205, t105 * t299 - t209 * t91, t110 * t105 + t156 * t91 + (-t11 + (pkin(12) * t94 + t204 * t44) * qJD(5)) * t209 + (pkin(12) * t53 + qJD(5) * t27 + t236) * t205, t109 * t105 - t157 * t91 + (-t10 + (pkin(12) * t95 + t208 * t44) * qJD(5)) * t209 + (pkin(12) * t52 - qJD(5) * t28 - t235) * t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t276, -t254, 0, 0, 0, 0, 0 (-t205 * t270 - t209 * t301) * t197 (t205 * t301 - t209 * t270) * t197, 0, 0, 0, 0, 0 (t166 * t300 - t100) * t209 + (-qJD(5) * t222 + t226) * t205 (t166 * t298 - t99) * t209 + (qJD(5) * t223 - t225) * t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t263, -0.2e1 * t265, 0, 0, 0, t205 * t293, t209 * t293, -0.2e1 * t193 * t273 + 0.2e1 * t194 * t272, 0.2e1 * t193 * t266 + t257 * t340, 0.2e1 * t205 * t278 + 0.2e1 * t298 * t308, -0.2e1 * t204 * t265 + 0.2e1 * t205 * t277, -0.2e1 * t272, 0.2e1 * t156 * t299 - 0.2e1 * t110 * t209 + 0.2e1 * (t193 * t295 + t204 * t263) * pkin(12), -0.2e1 * t157 * t299 - 0.2e1 * t109 * t209 + 0.2e1 * (-t193 * t296 + t208 * t263) * pkin(12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t40, t60, t6, t5, t295 * t47 + t326, qJD(6) * t247 - t204 * t24 + t325, t234, -t233, 0, -pkin(5) * t24 - pkin(13) * t234 + t239, -pkin(5) * t23 + pkin(13) * t233 + t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, -t91, t128, t38, t37, t295 * t95 + t324, qJD(6) * t246 - t204 * t53 + t323, t231, -t230, 0, -pkin(5) * t53 - pkin(13) * t231 + t235, -pkin(5) * t52 + pkin(13) * t230 + t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t146, 0, 0, 0, 0, 0, t225, t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t297, -t299, 0, -pkin(12) * t297, pkin(12) * t299, -t205 * t266 + t257, t273 * t340 - t297 * t309, t204 * t299 - t277, t221, 0 (pkin(13) * t311 + (-pkin(5) * t208 + t330) * t205) * qJD(6) + (t204 * t252 - t288) * qJD(5) (t208 * t329 + t232) * qJD(6) + (t208 * t252 + t289) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t273, -0.2e1 * t266, 0, 0, 0, t204 * t292, t208 * t292; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t24, t40, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, -t53, t91, t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205 * t296 + t279, -t204 * t297 - t205 * t295, t299, t110, t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t295, -t296, 0, -pkin(13) * t295, pkin(13) * t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
