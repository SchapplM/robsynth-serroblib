% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:45:23
% EndTime: 2019-03-09 00:45:36
% DurationCPUTime: 4.69s
% Computational Cost: add. (5613->383), mult. (13688->560), div. (0->0), fcn. (10832->12), ass. (0->233)
t190 = qJD(3) + qJD(4);
t197 = sin(qJ(4));
t202 = cos(qJ(3));
t321 = cos(qJ(4));
t254 = qJD(2) * t321;
t198 = sin(qJ(3));
t280 = qJD(2) * t198;
t342 = -t197 * t280 + t202 * t254;
t110 = t342 * t190;
t326 = qJD(5) + qJD(6);
t345 = t342 - t326;
t199 = sin(qJ(2));
t193 = sin(pkin(6));
t283 = qJD(1) * t193;
t263 = t199 * t283;
t315 = qJD(3) * pkin(3);
t211 = t198 * t315 - t263;
t219 = -t197 * t198 + t321 * t202;
t203 = cos(qJ(2));
t262 = t203 * t283;
t133 = t219 * t262;
t324 = pkin(8) + pkin(9);
t266 = qJD(3) * t324;
t161 = t198 * t266;
t162 = t202 * t266;
t170 = t324 * t198;
t172 = t324 * t202;
t328 = -t321 * t170 - t197 * t172;
t76 = qJD(4) * t328 - t321 * t161 - t197 * t162;
t344 = t133 - t76;
t253 = qJD(4) * t321;
t247 = qJD(2) * t324 + t263;
t194 = cos(pkin(6));
t282 = qJD(1) * t194;
t131 = t198 * t282 + t247 * t202;
t117 = t197 * t131;
t130 = -t247 * t198 + t202 * t282;
t69 = t321 * t130 - t117;
t337 = pkin(3) * t253 - t69;
t123 = t190 * t219;
t290 = t197 * t202;
t159 = t321 * t198 + t290;
t124 = t190 * t159;
t343 = -pkin(4) * t124 + pkin(10) * t123 - t211;
t196 = sin(qJ(5));
t201 = cos(qJ(5));
t153 = -qJD(2) * t290 - t198 * t254;
t113 = -pkin(4) * t153 - pkin(10) * t342;
t98 = pkin(3) * t280 + t113;
t341 = -t196 * t337 - t201 * t98;
t134 = -t153 * t196 - t201 * t190;
t119 = t130 + t315;
t66 = t321 * t119 - t117;
t56 = -t190 * pkin(4) - t66;
t44 = t134 * pkin(5) + t56;
t195 = sin(qJ(6));
t200 = cos(qJ(6));
t224 = t153 * t201 - t190 * t196;
t72 = t200 * t134 - t195 * t224;
t340 = t44 * t72;
t225 = t134 * t195 + t200 * t224;
t339 = t225 * t72;
t156 = t195 * t196 - t200 * t201;
t311 = t345 * t156;
t158 = t195 * t201 + t196 * t200;
t338 = t345 * t158;
t336 = t225 ^ 2 - t72 ^ 2;
t147 = qJD(5) - t342;
t144 = qJD(6) + t147;
t274 = qJD(6) * t200;
t275 = qJD(6) * t195;
t276 = qJD(5) * t201;
t277 = qJD(5) * t196;
t63 = t201 * t110 + t153 * t277 + t190 * t276;
t64 = -qJD(5) * t224 + t110 * t196;
t16 = -t134 * t274 - t195 * t64 + t200 * t63 + t224 * t275;
t335 = t144 * t72 + t16;
t111 = t124 * qJD(2);
t278 = qJD(4) * t197;
t281 = qJD(2) * t193;
t252 = qJD(1) * t281;
t240 = t203 * t252;
t89 = qJD(3) * t130 + t202 * t240;
t90 = -qJD(3) * t131 - t198 * t240;
t23 = t119 * t253 - t131 * t278 + t197 * t90 + t321 * t89;
t118 = t321 * t131;
t67 = t197 * t119 + t118;
t57 = pkin(10) * t190 + t67;
t188 = -pkin(3) * t202 - pkin(2);
t145 = t188 * qJD(2) - t262;
t84 = -pkin(4) * t342 + pkin(10) * t153 + t145;
t39 = t196 * t84 + t201 * t57;
t273 = qJD(2) * qJD(3);
t251 = t198 * t273;
t148 = pkin(3) * t251 + t199 * t252;
t48 = pkin(4) * t111 - pkin(10) * t110 + t148;
t47 = t201 * t48;
t208 = -qJD(5) * t39 - t196 * t23 + t47;
t3 = pkin(5) * t111 - pkin(11) * t63 + t208;
t221 = t196 * t48 + t201 * t23 + t84 * t276 - t57 * t277;
t4 = -pkin(11) * t64 + t221;
t267 = -t195 * t4 + t200 * t3;
t38 = -t196 * t57 + t201 * t84;
t26 = pkin(11) * t224 + t38;
t18 = pkin(5) * t147 + t26;
t27 = -pkin(11) * t134 + t39;
t313 = t200 * t27;
t7 = t18 * t195 + t313;
t334 = -qJD(6) * t7 + t44 * t225 + t267;
t207 = qJD(6) * t225 - t195 * t63 - t200 * t64;
t333 = -t144 * t225 + t207;
t332 = t133 * t196 - t343 * t201;
t114 = -pkin(4) * t219 - pkin(10) * t159 + t188;
t138 = -t197 * t170 + t321 * t172;
t331 = -t114 * t276 + t138 * t277 + t343 * t196 + t344 * t201;
t24 = t119 * t278 + t131 * t253 + t197 * t89 - t321 * t90;
t330 = t24 * t201 - t56 * t277;
t68 = t197 * t130 + t118;
t239 = pkin(3) * t278 - t68;
t306 = t138 * qJD(4) - t159 * t262 - t197 * t161 + t321 * t162;
t105 = t158 * t159;
t300 = t342 * t196;
t329 = (t277 - t300) * pkin(5);
t288 = t201 * t123;
t217 = -t159 * t277 + t288;
t327 = t196 * t98 - t337 * t201;
t25 = t27 * t275;
t250 = qJD(6) * t18 + t4;
t325 = t195 * t3 + t200 * t250 - t25;
t323 = -pkin(10) - pkin(11);
t129 = t201 * t138;
t322 = -pkin(11) * t288 + pkin(5) * t124 - t196 * t76 + (-t129 + (pkin(11) * t159 - t114) * t196) * qJD(5) + t332;
t320 = t201 * pkin(5);
t185 = pkin(3) * t197 + pkin(10);
t319 = -pkin(11) - t185;
t256 = t159 * t276;
t291 = t196 * t123;
t218 = t256 + t291;
t318 = pkin(11) * t218 + t331;
t316 = qJD(2) * pkin(2);
t314 = t196 * t63;
t309 = t329 + t239;
t308 = pkin(5) * t218 + t306;
t307 = t196 * t113 + t201 * t66;
t304 = t134 * t147;
t303 = t224 * t147;
t302 = t144 * t153;
t301 = t147 * t153;
t299 = t342 * t201;
t298 = t153 * t342;
t297 = t159 * t196;
t296 = t159 * t201;
t295 = t193 * t199;
t294 = t193 * t203;
t205 = qJD(2) ^ 2;
t293 = t193 * t205;
t292 = t196 * t111;
t289 = t201 * t111;
t204 = qJD(3) ^ 2;
t287 = t204 * t198;
t286 = t204 * t202;
t285 = t196 * t114 + t129;
t284 = t198 ^ 2 - t202 ^ 2;
t279 = qJD(2) * t199;
t272 = pkin(11) * t300;
t268 = t199 * t293;
t265 = qJD(5) * t323;
t260 = t193 * t279;
t259 = t203 * t281;
t248 = qJD(5) * t319;
t246 = t201 * t113 - t196 * t66;
t245 = t147 * t201;
t243 = t198 * t259;
t242 = t202 * t259;
t186 = -t321 * pkin(3) - pkin(4);
t241 = -t39 * t153 + t24 * t196 + t56 * t276;
t238 = -t67 + t329;
t237 = -t153 * pkin(5) - pkin(11) * t299;
t154 = t319 * t196;
t236 = -qJD(6) * t154 - t196 * t248 - t272 + t327;
t189 = t201 * pkin(11);
t155 = t185 * t201 + t189;
t235 = qJD(6) * t155 - t201 * t248 + t237 - t341;
t169 = t323 * t196;
t234 = -qJD(6) * t169 - t196 * t265 - t272 + t307;
t171 = pkin(10) * t201 + t189;
t233 = qJD(6) * t171 - t201 * t265 + t237 + t246;
t108 = t201 * t114;
t45 = -pkin(5) * t219 - pkin(11) * t296 - t138 * t196 + t108;
t51 = -pkin(11) * t297 + t285;
t231 = t195 * t45 + t200 * t51;
t149 = t194 * t202 - t198 * t295;
t150 = t194 * t198 + t202 * t295;
t100 = t197 * t149 + t321 * t150;
t222 = -t100 * t201 + t196 * t294;
t82 = -t100 * t196 - t201 * t294;
t230 = t195 * t222 + t200 * t82;
t229 = t195 * t82 - t200 * t222;
t226 = -t185 * t111 - t342 * t56;
t223 = t38 * t153 - t330;
t220 = t321 * t149 - t197 * t150;
t15 = pkin(5) * t64 + t24;
t6 = t18 * t200 - t195 * t27;
t216 = t15 * t156 + t6 * t153 - t338 * t44;
t215 = t15 * t158 - t7 * t153 + t311 * t44;
t214 = t145 * t153 - t24;
t213 = t316 * qJD(2);
t210 = -0.2e1 * qJD(3) * t316;
t206 = -t145 * t342 - t23;
t187 = -pkin(4) - t320;
t168 = t186 - t320;
t128 = -qJD(3) * t150 - t243;
t127 = qJD(3) * t149 + t242;
t106 = t156 * t159;
t95 = pkin(5) * t297 - t328;
t86 = t153 ^ 2 - t342 ^ 2;
t85 = t111 * t219;
t80 = (-qJD(2) * t159 - t153) * t190;
t43 = t100 * qJD(4) + t197 * t127 - t321 * t128;
t42 = t220 * qJD(4) + t321 * t127 + t197 * t128;
t34 = t147 * t245 - t153 * t224 + t292;
t33 = -t147 ^ 2 * t196 - t134 * t153 + t289;
t31 = -t224 * t245 + t314;
t29 = -t275 * t297 + (t296 * t326 + t291) * t200 + t217 * t195;
t28 = -t105 * t326 - t156 * t123;
t20 = qJD(5) * t222 - t196 * t42 + t201 * t260;
t19 = qJD(5) * t82 + t196 * t260 + t201 * t42;
t14 = -t111 * t156 + t144 * t338 - t153 * t72;
t13 = t111 * t158 + t311 * t144 - t153 * t225;
t8 = (t63 - t304) * t201 + (-t64 + t303) * t196;
t5 = t158 * t16 - t225 * t311;
t1 = -t156 * t16 + t158 * t207 - t225 * t338 - t311 * t72;
t2 = [0, 0, -t268, -t203 * t293, 0, 0, 0, 0, 0, -t202 * t268 + (t128 - t243) * qJD(3), t198 * t268 + (-t127 - t242) * qJD(3), 0, 0, 0, 0, 0, -t190 * t43 + (-t111 * t203 - t279 * t342) * t193, -t190 * t42 + (-t110 * t203 - t153 * t279) * t193, 0, 0, 0, 0, 0, t82 * t111 + t134 * t43 + t147 * t20 - t220 * t64, t111 * t222 - t147 * t19 - t220 * t63 - t224 * t43, 0, 0, 0, 0, 0 (-qJD(6) * t229 - t19 * t195 + t20 * t200) * t144 + t230 * t111 + t43 * t72 + t220 * t207 -(qJD(6) * t230 + t19 * t200 + t195 * t20) * t144 - t229 * t111 - t43 * t225 - t220 * t16; 0, 0, 0, 0, 0.2e1 * t202 * t251, -0.2e1 * t284 * t273, t286, -t287, 0, -pkin(8) * t286 + t198 * t210, pkin(8) * t287 + t202 * t210, t110 * t159 - t123 * t153, t110 * t219 - t111 * t159 + t123 * t342 + t124 * t153, t123 * t190, -t124 * t190, 0, t111 * t188 + t124 * t145 - t148 * t219 - t306 * t190 - t211 * t342, t110 * t188 + t123 * t145 + t148 * t159 - t211 * t153 + t344 * t190, -t217 * t224 + t296 * t63 (-t134 * t201 + t196 * t224) * t123 + (-t314 - t201 * t64 + (t134 * t196 + t201 * t224) * qJD(5)) * t159, -t124 * t224 + t147 * t217 + t159 * t289 - t219 * t63, -t124 * t134 - t147 * t218 - t159 * t292 + t219 * t64, t124 * t147 - t85, t108 * t111 - (-t276 * t57 + t47) * t219 + t38 * t124 - t328 * t64 + t56 * t256 + (-t138 * t276 + t332) * t147 + t306 * t134 + ((-qJD(5) * t114 - t76) * t147 - t138 * t111 - (-qJD(5) * t84 - t23) * t219 + t24 * t159 + t56 * t123) * t196, -t285 * t111 - t39 * t124 + t331 * t147 + t330 * t159 + t219 * t221 - t224 * t306 + t56 * t288 - t328 * t63, -t106 * t16 - t225 * t28, -t105 * t16 - t106 * t207 + t225 * t29 - t28 * t72, -t106 * t111 - t124 * t225 + t144 * t28 - t16 * t219, -t105 * t111 - t124 * t72 - t144 * t29 - t207 * t219, t124 * t144 - t85 (-t195 * t51 + t200 * t45) * t111 - t267 * t219 + t6 * t124 - t95 * t207 + t15 * t105 + t44 * t29 + t308 * t72 + (t318 * t195 + t322 * t200) * t144 + (-t144 * t231 + t219 * t7) * qJD(6), -t231 * t111 + t325 * t219 - t7 * t124 + t95 * t16 - t15 * t106 + t44 * t28 - t308 * t225 + ((-qJD(6) * t45 + t318) * t200 + (qJD(6) * t51 - t322) * t195) * t144; 0, 0, 0, 0, -t198 * t205 * t202, t284 * t205, 0, 0, 0, t198 * t213, t202 * t213, t298, t86, 0, t80, 0, t190 * t68 + (-t190 * t278 + t280 * t342) * pkin(3) + t214, t69 * t190 + (t153 * t280 - t190 * t253) * pkin(3) + t206, t31, t8, t34, t33, t301, t186 * t64 + t226 * t196 + t239 * t134 + (-t185 * t276 + t341) * t147 + t223, t186 * t63 + t226 * t201 - t239 * t224 + (t185 * t277 + t327) * t147 + t241, t5, t1, t13, t14, t302 (t154 * t200 - t155 * t195) * t111 - t168 * t207 + t309 * t72 + (t195 * t236 - t200 * t235) * t144 + t216 -(t154 * t195 + t155 * t200) * t111 + t168 * t16 - t309 * t225 + (t195 * t235 + t200 * t236) * t144 + t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t298, t86, 0, t80, 0, t190 * t67 + t214, t66 * t190 + t206, t31, t8, t34, t33, t301, -pkin(4) * t64 - t246 * t147 - t67 * t134 - t56 * t300 + (-t147 * t276 - t292) * pkin(10) + t223, -pkin(4) * t63 + t307 * t147 + t67 * t224 - t56 * t299 + (t147 * t277 - t289) * pkin(10) + t241, t5, t1, t13, t14, t302 (t169 * t200 - t171 * t195) * t111 - t187 * t207 + t238 * t72 + (t195 * t234 - t200 * t233) * t144 + t216 -(t169 * t195 + t171 * t200) * t111 + t187 * t16 - t238 * t225 + (t195 * t233 + t200 * t234) * t144 + t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t224 * t134, -t134 ^ 2 + t224 ^ 2, t63 + t304, -t64 - t303, t111, t147 * t39 + t224 * t56 + t208, t134 * t56 + t147 * t38 - t221, -t339, t336, t335, t333, t111 -(-t195 * t26 - t313) * t144 + (t111 * t200 - t144 * t275 + t224 * t72) * pkin(5) + t334, t340 + t25 + (-t144 * t27 - t3) * t195 + (t144 * t26 - t250) * t200 + (-t111 * t195 - t144 * t274 - t224 * t225) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t339, t336, t335, t333, t111, t144 * t7 + t334, t144 * t6 - t325 + t340;];
tauc_reg  = t2;
