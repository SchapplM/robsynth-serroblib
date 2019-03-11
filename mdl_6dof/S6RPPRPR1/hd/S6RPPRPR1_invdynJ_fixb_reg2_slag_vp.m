% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:40:09
% EndTime: 2019-03-09 01:40:20
% DurationCPUTime: 6.25s
% Computational Cost: add. (8433->472), mult. (18568->584), div. (0->0), fcn. (14074->18), ass. (0->238)
t204 = sin(pkin(10));
t207 = cos(pkin(10));
t212 = sin(qJ(4));
t328 = cos(qJ(4));
t164 = t328 * t204 + t212 * t207;
t151 = t164 * qJD(1);
t203 = sin(pkin(11));
t206 = cos(pkin(11));
t132 = t203 * qJD(4) + t206 * t151;
t211 = sin(qJ(6));
t327 = cos(qJ(6));
t285 = t203 * t151;
t343 = t206 * qJD(4) - t285;
t225 = t327 * t343;
t66 = -t211 * t132 + t225;
t347 = t66 ^ 2;
t267 = t328 * t207;
t176 = qJD(1) * t267;
t282 = t212 * t204;
t264 = qJD(1) * t282;
t149 = -t176 + t264;
t145 = qJD(6) + t149;
t346 = t66 * t145;
t65 = t327 * t132 + t211 * t343;
t345 = t65 ^ 2;
t201 = pkin(10) + qJ(4);
t192 = sin(t201);
t195 = cos(t201);
t202 = qJ(1) + pkin(9);
t193 = sin(t202);
t196 = cos(t202);
t253 = g(1) * t196 + g(2) * t193;
t222 = -g(3) * t195 + t253 * t192;
t205 = sin(pkin(9));
t179 = t205 * pkin(1) + qJ(3);
t159 = qJD(1) * qJD(3) + qJDD(1) * t179;
t188 = t207 * qJDD(2);
t121 = t188 + (-pkin(7) * qJDD(1) - t159) * t204;
t137 = t204 * qJDD(2) + t207 * t159;
t270 = t207 * qJDD(1);
t122 = pkin(7) * t270 + t137;
t171 = t179 * qJD(1);
t190 = t207 * qJD(2);
t309 = pkin(7) * qJD(1);
t134 = t190 + (-t171 - t309) * t204;
t142 = t204 * qJD(2) + t207 * t171;
t135 = t207 * t309 + t142;
t263 = qJD(4) * t328;
t276 = qJD(4) * t212;
t32 = t328 * t121 - t212 * t122 - t134 * t276 - t135 * t263;
t30 = -qJDD(4) * pkin(4) + qJDD(5) - t32;
t217 = t30 - t222;
t208 = cos(pkin(9));
t321 = t208 * pkin(1);
t186 = -pkin(2) - t321;
t272 = qJDD(1) * t186;
t166 = qJDD(3) + t272;
t326 = g(1) * t193;
t261 = g(2) * t196 - t326;
t344 = -t166 - t261;
t266 = t327 * t206;
t226 = -t211 * t203 + t266;
t275 = qJD(6) * t211;
t335 = qJD(6) * t266 - t203 * t275;
t302 = -t226 * t149 - t335;
t163 = t327 * t203 + t211 * t206;
t154 = t163 * qJD(6);
t301 = t163 * t149 + t154;
t300 = pkin(1) * qJDD(1);
t340 = t151 * t343;
t155 = t204 * t276 - t207 * t263;
t259 = qJDD(1) * t328;
t268 = qJD(4) * t176 + t204 * t259 + t212 * t270;
t107 = qJD(4) * t264 - t268;
t93 = t203 * qJDD(4) - t206 * t107;
t242 = t132 * t155 - t164 * t93;
t339 = t203 * t242;
t156 = t164 * qJD(4);
t271 = t204 * qJDD(1);
t248 = -t207 * t259 + t212 * t271;
t108 = qJD(1) * t156 + t248;
t303 = -t164 * t108 + t155 * t149;
t338 = t206 * t303;
t228 = t261 * t195;
t318 = pkin(7) + t179;
t157 = t318 * t204;
t158 = t318 * t207;
t336 = -t328 * t157 - t212 * t158;
t92 = -t206 * qJDD(4) - t203 * t107;
t21 = -qJD(6) * t225 + t132 * t275 + t211 * t92 - t327 * t93;
t334 = -t21 * t226 - t301 * t65;
t105 = qJDD(6) + t108;
t333 = t163 * t105 - t145 * t302;
t148 = t149 ^ 2;
t286 = t203 * t108;
t332 = -t206 * t148 - t286;
t227 = t267 - t282;
t331 = -t156 * t343 - t227 * t92;
t323 = g(3) * t192;
t220 = -t253 * t195 - t323;
t69 = t212 * t134 + t328 * t135;
t59 = qJD(4) * qJ(5) + t69;
t185 = t207 * pkin(3) + pkin(2);
t168 = -t185 - t321;
t147 = qJD(1) * t168 + qJD(3);
t77 = t149 * pkin(4) - t151 * qJ(5) + t147;
t34 = -t203 * t59 + t206 * t77;
t20 = t149 * pkin(5) - t132 * pkin(8) + t34;
t35 = t203 * t77 + t206 * t59;
t24 = pkin(8) * t343 + t35;
t231 = -t327 * t20 + t211 * t24;
t124 = t212 * t135;
t269 = t212 * t121 + t328 * t122 + t134 * t263;
t29 = qJDD(4) * qJ(5) + (qJD(5) - t124) * qJD(4) + t269;
t144 = qJDD(1) * t168 + qJDD(3);
t41 = t108 * pkin(4) + t107 * qJ(5) - t151 * qJD(5) + t144;
t12 = -t203 * t29 + t206 * t41;
t6 = t108 * pkin(5) - t93 * pkin(8) + t12;
t13 = t203 * t41 + t206 * t29;
t9 = -t92 * pkin(8) + t13;
t1 = -qJD(6) * t231 + t211 * t6 + t327 * t9;
t330 = t151 ^ 2;
t329 = t92 * pkin(5);
t322 = t195 * pkin(4);
t213 = sin(qJ(1));
t320 = t213 * pkin(1);
t319 = t65 * t66;
t317 = pkin(8) + qJ(5);
t100 = t226 * t164;
t260 = t211 * t93 + t327 * t92;
t22 = qJD(6) * t65 + t260;
t47 = t154 * t164 + t226 * t155;
t316 = -t100 * t22 - t47 * t66;
t48 = -t155 * t163 + t335 * t164;
t99 = t163 * t164;
t315 = -t99 * t105 - t48 * t145;
t314 = t65 * t156 + t21 * t227;
t104 = t151 * pkin(4) + t149 * qJ(5);
t68 = t328 * t134 - t124;
t45 = t203 * t104 + t206 * t68;
t233 = t343 * t155;
t291 = t164 * t206;
t313 = -t206 * t233 - t92 * t291;
t70 = t227 * qJD(3) + t336 * qJD(4);
t78 = t156 * pkin(4) + t155 * qJ(5) - t164 * qJD(5);
t38 = t203 * t78 + t206 * t70;
t312 = t132 * t156 - t227 * t93;
t169 = t317 * t203;
t170 = t317 * t206;
t125 = -t327 * t169 - t211 * t170;
t295 = t149 * t206;
t44 = t206 * t104 - t203 * t68;
t26 = t151 * pkin(5) + pkin(8) * t295 + t44;
t296 = t149 * t203;
t33 = pkin(8) * t296 + t45;
t311 = qJD(5) * t226 + qJD(6) * t125 - t211 * t26 - t327 * t33;
t126 = -t211 * t169 + t327 * t170;
t310 = -qJD(5) * t163 - qJD(6) * t126 + t211 * t33 - t327 * t26;
t103 = -t212 * t157 + t328 * t158;
t96 = -pkin(4) * t227 - t164 * qJ(5) + t168;
t50 = t206 * t103 + t203 * t96;
t308 = t151 * t66;
t305 = t65 * t151;
t304 = t93 * t206;
t299 = t132 * t151;
t298 = t132 * t203;
t297 = t149 * t151;
t294 = t155 * t203;
t292 = t164 * t203;
t290 = t192 * t196;
t289 = t193 * t195;
t200 = pkin(11) + qJ(6);
t191 = sin(t200);
t288 = t196 * t191;
t194 = cos(t200);
t287 = t196 * t194;
t101 = t206 * t108;
t281 = t69 * qJD(4);
t57 = -qJD(4) * pkin(4) + qJD(5) - t68;
t280 = -qJD(5) + t57;
t279 = -t203 * t148 + t101;
t214 = cos(qJ(1));
t197 = t214 * pkin(1);
t278 = t196 * t185 + t197;
t198 = t204 ^ 2;
t199 = t207 ^ 2;
t277 = t198 + t199;
t210 = -pkin(7) - qJ(3);
t262 = pkin(5) * t203 - t210;
t37 = -t203 * t70 + t206 * t78;
t49 = -t203 * t103 + t206 * t96;
t257 = -t163 * t22 - t302 * t66;
t255 = g(2) * t290 - t192 * t326;
t254 = t226 * t105 - t301 * t145;
t251 = g(1) * t213 - g(2) * t214;
t250 = -t99 * t21 + t65 * t48;
t249 = -t196 * t210 - t320;
t247 = -t12 * t206 - t13 * t203;
t246 = -t12 * t203 + t13 * t206;
t245 = t156 * t66 + t22 * t227;
t244 = t34 * t203 - t35 * t206;
t243 = -t100 * t105 + t47 * t145;
t241 = t107 * t227 + t151 * t156;
t239 = -t108 * t227 + t149 * t156;
t136 = -t204 * t159 + t188;
t238 = -t136 * t204 + t137 * t207;
t237 = (-t204 * t171 + t190) * t204 - t142 * t207;
t184 = t206 * pkin(5) + pkin(4);
t236 = t195 * t184 + t192 * t317;
t234 = t206 * t343;
t36 = -pkin(5) * t227 - pkin(8) * t291 + t49;
t43 = -pkin(8) * t292 + t50;
t14 = -t211 * t43 + t327 * t36;
t8 = t211 * t20 + t327 * t24;
t15 = t211 * t36 + t327 * t43;
t224 = t303 * t203;
t223 = -t272 + t344;
t219 = -t57 * t155 + t30 * t164 - t253;
t2 = -qJD(6) * t8 - t211 * t9 + t327 * t6;
t71 = qJD(3) * t164 + qJD(4) * t103;
t130 = t193 * t191 + t195 * t287;
t129 = t193 * t194 - t195 * t288;
t128 = -t194 * t289 + t288;
t127 = t191 * t289 + t287;
t110 = -t156 * qJD(4) + qJDD(4) * t227;
t109 = -t155 * qJD(4) + t164 * qJDD(4);
t81 = t203 * t92;
t74 = pkin(5) * t292 - t336;
t52 = -pkin(5) * t294 + t71;
t51 = -pkin(5) * t296 + t69;
t46 = -pkin(5) * t343 + t57;
t31 = -t135 * t276 + t269;
t28 = pkin(8) * t294 + t38;
t23 = t206 * t155 * pkin(8) + t156 * pkin(5) + t37;
t17 = t30 + t329;
t4 = -t15 * qJD(6) - t211 * t28 + t327 * t23;
t3 = t14 * qJD(6) + t211 * t23 + t327 * t28;
t5 = [0, 0, 0, 0, 0, qJDD(1), t251, g(1) * t214 + g(2) * t213, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t208 * t300 - t261, -0.2e1 * t205 * t300 + t253, 0 (t251 + (t205 ^ 2 + t208 ^ 2) * t300) * pkin(1), t198 * qJDD(1), 0.2e1 * t204 * t270, 0, t199 * qJDD(1), 0, 0, t223 * t207, -t223 * t204, t159 * t277 + t238 - t253, t166 * t186 - g(1) * (-t193 * pkin(2) + t196 * qJ(3) - t320) - g(2) * (t196 * pkin(2) + t193 * qJ(3) + t197) + t238 * t179 - t237 * qJD(3), -t107 * t164 - t151 * t155, -t241 + t303, t109, t239, t110, 0, -t71 * qJD(4) + qJDD(4) * t336 + t168 * t108 - t144 * t227 + t147 * t156 - t228, -t70 * qJD(4) - t103 * qJDD(4) - t168 * t107 + t144 * t164 - t147 * t155 + t255, -t103 * t108 + t107 * t336 - t70 * t149 + t71 * t151 + t68 * t155 - t69 * t156 - t32 * t164 + t227 * t31 - t253, t31 * t103 + t69 * t70 + t32 * t336 - t68 * t71 + t144 * t168 - g(1) * (-t193 * t185 + t249) - g(2) * (-t193 * t210 + t278) -t242 * t206, t313 + t339, t312 - t338 (t92 * t164 + t233) * t203, t224 - t331, t239, t49 * t108 - t12 * t227 + t37 * t149 + t34 * t156 + t203 * t219 - t206 * t228 - t336 * t92 - t343 * t71, -t50 * t108 + t13 * t227 + t71 * t132 - t38 * t149 - t35 * t156 + t203 * t228 + t206 * t219 - t336 * t93, t38 * t343 - t50 * t92 - t37 * t132 - t49 * t93 + t247 * t164 + (t35 * t203 + t34 * t206) * t155 - t255, t13 * t50 + t35 * t38 + t12 * t49 + t34 * t37 - t30 * t336 + t57 * t71 - g(1) * t249 - g(2) * (qJ(5) * t290 + t196 * t322 + t278) + (-g(1) * (-t192 * qJ(5) - t185 - t322) + g(2) * t210) * t193, -t21 * t100 - t65 * t47, -t250 + t316, -t243 + t314, t22 * t99 - t48 * t66, t245 + t315, -t105 * t227 + t145 * t156, -g(1) * t128 - g(2) * t130 + t14 * t105 + t4 * t145 - t156 * t231 + t17 * t99 - t2 * t227 + t74 * t22 + t46 * t48 - t52 * t66, -g(1) * t127 - g(2) * t129 + t1 * t227 + t17 * t100 - t15 * t105 - t3 * t145 - t8 * t156 - t74 * t21 - t46 * t47 + t52 * t65, -t1 * t99 - t2 * t100 + t14 * t21 - t15 * t22 - t231 * t47 + t3 * t66 - t4 * t65 - t8 * t48 - t255, t1 * t15 + t8 * t3 + t2 * t14 - t231 * t4 + t17 * t74 + t46 * t52 + g(1) * t320 - g(2) * t278 + (-g(1) * t262 - g(2) * t236) * t196 + (-g(1) * (-t185 - t236) - g(2) * t262) * t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t136 * t207 + t137 * t204 - g(3), 0, 0, 0, 0, 0, 0, t110, -t109, t241 + t303, -t69 * t155 - t68 * t156 + t31 * t164 + t227 * t32 - g(3), 0, 0, 0, 0, 0, 0, t224 + t331, t312 + t338, t313 - t339, t155 * t244 + t57 * t156 + t164 * t246 - t227 * t30 - g(3), 0, 0, 0, 0, 0, 0, -t245 + t315, t243 + t314, t250 + t316, t1 * t100 + t46 * t156 - t17 * t227 - t2 * t99 + t231 * t48 - t8 * t47 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t270, t271, -t277 * qJD(1) ^ 2, qJD(1) * t237 - t344, 0, 0, 0, 0, 0, 0, 0.2e1 * t151 * qJD(4) + t248 (-t149 - t264) * qJD(4) + t268, -t148 - t330, t69 * t149 + t68 * t151 + t144 + t261, 0, 0, 0, 0, 0, 0, t279 + t340, -t299 + t332, -t304 - t81 + (t234 + t298) * t149, -t149 * t244 - t57 * t151 - t247 + t261, 0, 0, 0, 0, 0, 0, t254 + t308, -t305 - t333, t257 - t334, t1 * t163 - t46 * t151 + t2 * t226 + t231 * t301 - t302 * t8 + t261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t297, -t148 + t330 (t149 - t264) * qJD(4) + t268, -t297, -t248, qJDD(4), -t147 * t151 + t222 + t281 + t32, t147 * t149 + (t68 + t124) * qJD(4) - t269 - t220, 0, 0, t132 * t295 + t93 * t203, t304 - t81 + (t234 - t298) * t149, -t299 - t332, -t92 * t206 - t296 * t343, t279 - t340, -t297, -qJ(5) * t286 - pkin(4) * t92 - t69 * t285 - t34 * t151 + (t203 * t280 - t44) * t149 + (-t217 + t281) * t206, -qJ(5) * t101 - pkin(4) * t93 - t69 * t132 + t35 * t151 + (t206 * t280 + t45) * t149 + t217 * t203, t44 * t132 + t45 * t285 + (-qJ(5) * t92 - qJD(5) * t285 - t34 * t149 + t13 + (qJD(5) * t206 - t45) * qJD(4)) * t206 + (qJ(5) * t93 + qJD(5) * t132 - t35 * t149 - t12) * t203 + t220, -t34 * t44 - t35 * t45 - t57 * t69 - t244 * qJD(5) - t217 * pkin(4) + (t220 + t246) * qJ(5), -t21 * t163 - t302 * t65, t257 + t334, -t305 + t333, -t22 * t226 - t301 * t66, t254 - t308, -t145 * t151, t125 * t105 + t145 * t310 + t151 * t231 - t17 * t226 - t184 * t22 + t194 * t222 + t301 * t46 + t51 * t66, -t126 * t105 - t145 * t311 + t8 * t151 + t17 * t163 + t184 * t21 - t191 * t222 - t302 * t46 - t51 * t65, t1 * t226 + t125 * t21 - t126 * t22 - t2 * t163 - t231 * t302 - t301 * t8 - t310 * t65 + t311 * t66 + t220, -g(3) * t236 + t1 * t126 + t2 * t125 - t17 * t184 - t310 * t231 + t311 * t8 - t46 * t51 + t253 * (t184 * t192 - t195 * t317); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132 * t149 + t92, t149 * t343 + t93, -t132 ^ 2 - t343 ^ 2, t132 * t34 - t343 * t35 + t217, 0, 0, 0, 0, 0, 0, t65 * t145 + t22, -t21 + t346, -t345 - t347, -t231 * t65 - t8 * t66 + t217 + t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t319, t345 - t347, -t21 - t346, t319, -t260 + (-qJD(6) + t145) * t65, t105, -g(1) * t129 + g(2) * t127 + t8 * t145 + t191 * t323 - t46 * t65 + t2, g(1) * t130 - g(2) * t128 - t145 * t231 + t194 * t323 - t46 * t66 - t1, 0, 0;];
tau_reg  = t5;
