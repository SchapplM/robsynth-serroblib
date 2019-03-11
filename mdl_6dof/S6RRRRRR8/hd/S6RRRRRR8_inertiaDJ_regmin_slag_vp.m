% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x38]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRRRR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_inertiaDJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:06:45
% EndTime: 2019-03-10 05:07:10
% DurationCPUTime: 8.73s
% Computational Cost: add. (19125->592), mult. (56894->1075), div. (0->0), fcn. (59595->14), ass. (0->261)
t229 = sin(pkin(6));
t238 = cos(qJ(3));
t228 = sin(pkin(7));
t341 = cos(pkin(6));
t286 = t341 * t228;
t230 = cos(pkin(7));
t239 = cos(qJ(2));
t326 = t238 * t239;
t234 = sin(qJ(3));
t235 = sin(qJ(2));
t329 = t234 * t235;
t359 = t230 * t326 - t329;
t149 = -t229 * t359 - t238 * t286;
t232 = sin(qJ(5));
t353 = cos(qJ(5));
t327 = t235 * t238;
t328 = t234 * t239;
t260 = t230 * t328 + t327;
t150 = t229 * t260 + t234 * t286;
t333 = t229 * t239;
t186 = t228 * t333 - t230 * t341;
t233 = sin(qJ(4));
t237 = cos(qJ(4));
t123 = t150 * t237 - t186 * t233;
t352 = pkin(4) * t149;
t274 = t229 * (-pkin(10) * t230 - pkin(9));
t264 = t235 * t274;
t301 = pkin(1) * t341;
t157 = pkin(2) * t341 + t239 * t301 + t264;
t351 = pkin(10) * t228;
t172 = (-pkin(2) * t239 - t235 * t351 - pkin(1)) * t229;
t115 = -t157 * t228 + t172 * t230;
t86 = pkin(3) * t149 - pkin(11) * t150 + t115;
t332 = t230 * t234;
t138 = t157 * t332;
t278 = t235 * t301;
t247 = -pkin(9) * t333 - t278;
t145 = (t230 * t333 + t286) * pkin(10) - t247;
t335 = t228 * t234;
t302 = t145 * t238 + t172 * t335 + t138;
t91 = -pkin(11) * t186 + t302;
t51 = -t233 * t91 + t237 * t86;
t45 = -pkin(12) * t123 + t352 + t51;
t122 = t150 * t233 + t186 * t237;
t52 = t233 * t86 + t237 * t91;
t47 = -pkin(12) * t122 + t52;
t24 = -t232 * t47 + t353 * t45;
t22 = -pkin(5) * t149 - t24;
t236 = cos(qJ(6));
t231 = sin(qJ(6));
t313 = qJD(6) * t231;
t272 = qJD(3) * t286;
t118 = t234 * t272 + (t260 * qJD(3) + (t230 * t327 + t328) * qJD(2)) * t229;
t119 = t238 * t272 + (t359 * qJD(3) + (-t230 * t329 + t326) * qJD(2)) * t229;
t159 = (t239 * t274 - t278) * qJD(2);
t177 = (pkin(2) * t235 - t239 * t351) * t229 * qJD(2);
t120 = -t159 * t228 + t177 * t230;
t245 = pkin(3) * t118 - pkin(11) * t119 + t120;
t321 = qJD(2) * t235;
t293 = t229 * t321;
t275 = t228 * t293;
t285 = qJD(2) * t341;
t273 = t239 * t285;
t210 = pkin(1) * t273;
t158 = qJD(2) * t264 + t210;
t318 = qJD(3) * t238;
t294 = t228 * t318;
t295 = t230 * t318;
t319 = qJD(3) * t234;
t59 = t145 * t319 - t157 * t295 - t158 * t238 - t159 * t332 - t172 * t294 - t177 * t335;
t55 = pkin(11) * t275 - t59;
t27 = -qJD(4) * t52 - t233 * t55 + t237 * t245;
t80 = -qJD(4) * t122 + t119 * t237 + t233 * t275;
t16 = pkin(4) * t118 - pkin(12) * t80 + t27;
t316 = qJD(4) * t237;
t317 = qJD(4) * t233;
t26 = -t233 * t245 - t237 * t55 - t316 * t86 + t317 * t91;
t79 = t119 * t233 + t150 * t316 - t186 * t317 - t237 * t275;
t18 = -pkin(12) * t79 - t26;
t289 = -t16 * t353 + t232 * t18;
t46 = t353 * t47;
t347 = t232 * t45 + t46;
t7 = -qJD(5) * t347 - t289;
t5 = -t118 * pkin(5) - t7;
t291 = t22 * t313 - t5 * t236;
t213 = t228 * t319;
t187 = -t230 * t237 + t233 * t335;
t160 = -qJD(4) * t187 + t237 * t294;
t282 = pkin(4) * t213;
t334 = t228 * t238;
t323 = pkin(2) * t332 + pkin(10) * t334;
t179 = pkin(11) * t230 + t323;
t180 = (-pkin(3) * t238 - pkin(11) * t234 - pkin(2)) * t228;
t130 = t179 * t237 + t180 * t233;
t182 = -pkin(2) * t295 + pkin(10) * t213;
t250 = qJD(3) * t228 * (pkin(3) * t234 - pkin(11) * t238);
t99 = -qJD(4) * t130 + t233 * t182 + t237 * t250;
t81 = -pkin(12) * t160 + t282 + t99;
t304 = t237 * t335;
t161 = qJD(4) * t304 + t230 * t317 + t233 * t294;
t98 = t179 * t317 - t180 * t316 + t182 * t237 - t233 * t250;
t89 = -pkin(12) * t161 - t98;
t288 = t232 * t89 - t353 * t81;
t129 = -t179 * t233 + t180 * t237;
t188 = t230 * t233 + t304;
t311 = pkin(4) * t334;
t106 = -pkin(12) * t188 + t129 - t311;
t110 = -pkin(12) * t187 + t130;
t107 = t353 * t110;
t360 = t106 * t232 + t107;
t39 = -qJD(5) * t360 - t288;
t36 = -pkin(5) * t213 - t39;
t74 = t106 * t353 - t110 * t232;
t72 = pkin(5) * t334 - t74;
t287 = -t36 * t236 + t313 * t72;
t227 = t236 ^ 2;
t322 = t231 ^ 2 - t227;
t283 = t322 * qJD(6);
t358 = qJD(4) + qJD(5);
t357 = -t234 * t145 + t238 * (t157 * t230 + t172 * t228);
t356 = 0.2e1 * t228;
t355 = pkin(11) + pkin(12);
t223 = qJD(6) * t236;
t354 = t22 * t223 + t231 * t5;
t350 = pkin(11) * t228;
t348 = t223 * t72 + t231 * t36;
t346 = pkin(4) * qJD(5);
t94 = -t122 * t232 + t123 * t353;
t270 = t149 * t236 - t231 * t94;
t258 = -t122 * t353 - t123 * t232;
t40 = qJD(5) * t258 - t232 * t79 + t353 * t80;
t28 = qJD(6) * t270 + t118 * t231 + t236 * t40;
t345 = t28 * t231;
t344 = t28 * t236;
t132 = -t187 * t232 + t188 * t353;
t124 = t132 * t231 + t236 * t334;
t257 = -t187 * t353 - t188 * t232;
t95 = qJD(5) * t257 + t160 * t353 - t161 * t232;
t64 = -qJD(6) * t124 + t213 * t231 + t236 * t95;
t342 = t64 * t231;
t256 = -t232 * t233 + t237 * t353;
t155 = t358 * t256;
t340 = t155 * t231;
t339 = t155 * t236;
t196 = t232 * t237 + t233 * t353;
t338 = t196 * t155;
t337 = t196 * t231;
t336 = t196 * t236;
t331 = t230 * t238;
t330 = t231 * t236;
t207 = t355 * t233;
t208 = t355 * t237;
t166 = -t207 * t232 + t208 * t353;
t300 = qJD(4) * t355;
t266 = t353 * t300;
t279 = t232 * t300;
t117 = qJD(5) * t166 - t233 * t279 + t237 * t266;
t165 = t207 * t353 + t208 * t232;
t325 = t117 * t231 + t165 * t223;
t310 = t353 * pkin(4);
t221 = -t310 - pkin(5);
t314 = qJD(5) * t232;
t306 = pkin(4) * t314;
t324 = t221 * t223 + t231 * t306;
t320 = qJD(2) * t239;
t315 = qJD(4) * t238;
t312 = -0.2e1 * pkin(3) * qJD(4);
t309 = pkin(4) * t317;
t308 = pkin(5) * t313;
t307 = pkin(5) * t223;
t305 = t231 * t334;
t222 = -pkin(4) * t237 - pkin(3);
t299 = t231 * t353;
t298 = t236 * t353;
t225 = t229 ^ 2;
t297 = t225 * t320;
t224 = t228 ^ 2;
t296 = t224 * t318;
t292 = t231 * t223;
t290 = qJD(5) * t353;
t284 = -0.4e1 * t196 * t330;
t281 = pkin(4) * t290;
t280 = -qJD(3) * t138 - t145 * t318 - t158 * t234 - t172 * t213;
t277 = t224 * t293;
t276 = t234 * t296;
t23 = pkin(13) * t149 + t347;
t90 = t186 * pkin(3) - t357;
t63 = t122 * pkin(4) + t90;
t37 = -pkin(5) * t258 - pkin(13) * t94 + t63;
t11 = t23 * t236 + t231 * t37;
t67 = t149 * t231 + t236 * t94;
t271 = -t231 * t67 + t236 * t270;
t73 = -pkin(13) * t334 + t360;
t215 = pkin(10) * t335;
t178 = t215 + (-pkin(2) * t238 - pkin(3)) * t230;
t137 = t187 * pkin(4) + t178;
t92 = -pkin(5) * t257 - t132 * pkin(13) + t137;
t49 = t231 * t92 + t236 * t73;
t125 = t132 * t236 - t305;
t269 = -t124 * t236 - t125 * t231;
t143 = -pkin(5) * t256 - pkin(13) * t196 + t222;
t109 = t143 * t231 + t166 * t236;
t220 = pkin(4) * t232 + pkin(13);
t267 = -t196 * t221 - t220 * t256;
t265 = t221 * t313 - t236 * t306;
t41 = qJD(5) * t94 + t232 * t80 + t353 * t79;
t31 = -t223 * t258 + t231 * t41;
t262 = -t236 * t41 - t258 * t313;
t96 = qJD(5) * t132 + t160 * t232 + t161 * t353;
t69 = -t223 * t257 + t231 * t96;
t261 = -t236 * t96 - t257 * t313;
t6 = -t16 * t232 - t18 * t353 - t290 * t45 + t314 * t47;
t255 = t118 * t233 + t149 * t316;
t254 = -t118 * t237 + t149 * t317;
t253 = t196 * t223 + t340;
t252 = t196 * t313 - t339;
t156 = t358 * t196;
t251 = -t156 * t236 - t256 * t313;
t38 = -t106 * t290 + t110 * t314 - t232 * t81 - t353 * t89;
t249 = t233 * t315 + t237 * t319;
t248 = t233 * t319 - t237 * t315;
t246 = pkin(5) * t156 - pkin(13) * t155 + t309;
t183 = t323 * qJD(3);
t244 = -pkin(13) * t118 + t6;
t243 = pkin(13) * t213 - t38;
t128 = t161 * pkin(4) + t183;
t242 = t155 * t221 - t156 * t220 + (t196 * t232 + t256 * t353) * t346;
t241 = t96 * pkin(5) - t95 * pkin(13) + t128;
t56 = -t159 * t331 + (-pkin(3) * t293 - t177 * t238) * t228 - t280;
t42 = pkin(4) * t79 + t56;
t240 = pkin(5) * t41 - pkin(13) * t40 + t42;
t212 = 0.2e1 * t292;
t199 = -0.2e1 * t276;
t194 = -0.2e1 * t283;
t193 = t196 ^ 2;
t185 = t247 * qJD(2);
t184 = pkin(9) * t293 - t210;
t153 = t165 * t313;
t127 = t156 * t231 - t223 * t256;
t116 = t207 * t290 + t208 * t314 + t233 * t266 + t237 * t279;
t108 = t143 * t236 - t166 * t231;
t105 = t155 * t330 - t196 * t283;
t101 = qJD(6) * t284 - t155 * t322;
t100 = 0.2e1 * t149 * t118;
t97 = (-t118 * t238 + t149 * t319) * t228;
t65 = -qJD(6) * t305 + t132 * t223 - t213 * t236 + t231 * t95;
t60 = (t159 * t230 + t177 * t228) * t238 + t280;
t58 = -qJD(6) * t109 + t231 * t116 + t236 * t246;
t57 = t116 * t236 - t143 * t223 + t166 * t313 - t231 * t246;
t53 = t125 * t223 + t342;
t48 = -t231 * t73 + t236 * t92;
t32 = qJD(6) * t269 - t231 * t65 + t64 * t236;
t29 = qJD(6) * t67 - t118 * t236 + t231 * t40;
t19 = t223 * t67 + t345;
t13 = -qJD(6) * t49 - t231 * t243 + t236 * t241;
t12 = -t223 * t92 - t231 * t241 - t236 * t243 + t313 * t73;
t10 = -t23 * t231 + t236 * t37;
t8 = qJD(6) * t271 - t231 * t29 + t344;
t2 = -qJD(6) * t11 + t231 * t244 + t236 * t240;
t1 = -t223 * t37 + t23 * t313 - t231 * t240 + t236 * t244;
t3 = [0, 0, 0, 0.2e1 * t235 * t297, 0.2e1 * (-t235 ^ 2 + t239 ^ 2) * t225 * qJD(2), 0.2e1 * t229 * t273, -0.2e1 * t229 * t235 * t285, 0, -0.2e1 * pkin(1) * t225 * t321 + 0.2e1 * t185 * t341, -0.2e1 * pkin(1) * t297 + 0.2e1 * t184 * t341, 0.2e1 * t150 * t119, -0.2e1 * t118 * t150 - 0.2e1 * t119 * t149, -0.2e1 * t119 * t186 + 0.2e1 * t150 * t275, 0.2e1 * t118 * t186 - 0.2e1 * t149 * t275, -0.2e1 * t186 * t275, 0.2e1 * t115 * t118 + 0.2e1 * t120 * t149 - 0.2e1 * t60 * t186 + 0.2e1 * t275 * t357, 0.2e1 * t115 * t119 + 0.2e1 * t120 * t150 - 0.2e1 * t186 * t59 - 0.2e1 * t275 * t302, 0.2e1 * t123 * t80, -0.2e1 * t122 * t80 - 0.2e1 * t123 * t79, 0.2e1 * t118 * t123 + 0.2e1 * t149 * t80, -0.2e1 * t118 * t122 - 0.2e1 * t149 * t79, t100, 0.2e1 * t118 * t51 + 0.2e1 * t122 * t56 + 0.2e1 * t149 * t27 + 0.2e1 * t79 * t90, -0.2e1 * t118 * t52 + 0.2e1 * t123 * t56 + 0.2e1 * t149 * t26 + 0.2e1 * t80 * t90, 0.2e1 * t94 * t40, 0.2e1 * t258 * t40 - 0.2e1 * t41 * t94, 0.2e1 * t118 * t94 + 0.2e1 * t149 * t40, 0.2e1 * t118 * t258 - 0.2e1 * t149 * t41, t100, 0.2e1 * t118 * t24 + 0.2e1 * t149 * t7 - 0.2e1 * t258 * t42 + 0.2e1 * t41 * t63, -0.2e1 * t118 * t347 + 0.2e1 * t149 * t6 + 0.2e1 * t40 * t63 + 0.2e1 * t42 * t94, 0.2e1 * t67 * t28, 0.2e1 * t270 * t28 - 0.2e1 * t29 * t67, -0.2e1 * t258 * t28 + 0.2e1 * t41 * t67, 0.2e1 * t258 * t29 + 0.2e1 * t270 * t41, -0.2e1 * t258 * t41, 0.2e1 * t10 * t41 - 0.2e1 * t2 * t258 + 0.2e1 * t22 * t29 - 0.2e1 * t270 * t5, -0.2e1 * t1 * t258 - 0.2e1 * t11 * t41 + 0.2e1 * t22 * t28 + 0.2e1 * t5 * t67; 0, 0, 0, 0, 0, t229 * t320, -t293, 0, t185, t184 (t119 * t234 + t150 * t318) * t228 (-t118 * t234 + t119 * t238 + (-t149 * t238 - t150 * t234) * qJD(3)) * t228, t119 * t230 - t186 * t294 + t234 * t277, -t118 * t230 + t186 * t213 + t238 * t277, t230 * t275, t183 * t186 + t60 * t230 + ((pkin(2) * t331 - t215) * t293 - pkin(2) * t118 - t120 * t238 + t115 * t319) * t228, -t182 * t186 + t59 * t230 + (-pkin(2) * t119 + t115 * t318 + t120 * t234 - t293 * t323) * t228, t123 * t160 + t188 * t80, -t122 * t160 - t123 * t161 - t187 * t80 - t188 * t79, t188 * t118 + t160 * t149 + (t123 * t319 - t238 * t80) * t228, -t187 * t118 - t161 * t149 + (-t122 * t319 + t238 * t79) * t228, t97, t118 * t129 + t122 * t183 + t149 * t99 + t161 * t90 + t178 * t79 + t187 * t56 + (-t238 * t27 + t319 * t51) * t228, -t118 * t130 + t123 * t183 + t149 * t98 + t160 * t90 + t178 * t80 + t188 * t56 + (-t238 * t26 - t319 * t52) * t228, t132 * t40 + t94 * t95, -t132 * t41 + t257 * t40 + t258 * t95 - t94 * t96, t132 * t118 + t95 * t149 + (-t238 * t40 + t319 * t94) * t228, t257 * t118 - t96 * t149 + (t238 * t41 + t258 * t319) * t228, t97, t118 * t74 - t128 * t258 - t257 * t42 + t137 * t41 + t149 * t39 + t63 * t96 + (-t238 * t7 + t24 * t319) * t228, -t118 * t360 + t128 * t94 + t132 * t42 + t137 * t40 + t149 * t38 + t63 * t95 + (-t238 * t6 - t319 * t347) * t228, t125 * t28 + t64 * t67, -t124 * t28 - t125 * t29 + t270 * t64 - t65 * t67, t125 * t41 - t257 * t28 - t258 * t64 + t67 * t96, -t124 * t41 + t257 * t29 + t258 * t65 + t270 * t96, -t257 * t41 - t258 * t96, t10 * t96 + t124 * t5 - t13 * t258 - t2 * t257 + t22 * t65 - t270 * t36 + t29 * t72 + t41 * t48, -t1 * t257 - t11 * t96 - t12 * t258 + t125 * t5 + t22 * t64 + t28 * t72 + t36 * t67 - t41 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t276, 0.2e1 * (-t234 ^ 2 + t238 ^ 2) * t224 * qJD(3), 0.2e1 * t230 * t294, -0.2e1 * t230 * t213, 0, -0.2e1 * pkin(2) * t224 * t319 - 0.2e1 * t183 * t230, -0.2e1 * pkin(2) * t296 + 0.2e1 * t182 * t230, 0.2e1 * t188 * t160, -0.2e1 * t160 * t187 - 0.2e1 * t161 * t188 (-t160 * t238 + t188 * t319) * t356 (t161 * t238 - t187 * t319) * t356, t199, 0.2e1 * t178 * t161 + 0.2e1 * t183 * t187 + 0.2e1 * (t129 * t319 - t238 * t99) * t228, 0.2e1 * t178 * t160 + 0.2e1 * t183 * t188 + 0.2e1 * (-t130 * t319 - t238 * t98) * t228, 0.2e1 * t132 * t95, -0.2e1 * t132 * t96 + 0.2e1 * t257 * t95 (t132 * t319 - t238 * t95) * t356 (t238 * t96 + t257 * t319) * t356, t199, -0.2e1 * t128 * t257 + 0.2e1 * t137 * t96 + 0.2e1 * (-t238 * t39 + t319 * t74) * t228, 0.2e1 * t128 * t132 + 0.2e1 * t137 * t95 + 0.2e1 * (-t238 * t38 - t319 * t360) * t228, 0.2e1 * t125 * t64, -0.2e1 * t124 * t64 - 0.2e1 * t125 * t65, 0.2e1 * t125 * t96 - 0.2e1 * t257 * t64, -0.2e1 * t124 * t96 + 0.2e1 * t257 * t65, -0.2e1 * t257 * t96, 0.2e1 * t124 * t36 - 0.2e1 * t13 * t257 + 0.2e1 * t48 * t96 + 0.2e1 * t65 * t72, -0.2e1 * t12 * t257 + 0.2e1 * t125 * t36 - 0.2e1 * t49 * t96 + 0.2e1 * t64 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, -t118, t275, t60, t59, t123 * t316 + t233 * t80, -t233 * t79 + t80 * t237 + (-t122 * t237 - t123 * t233) * qJD(4), t255, -t254, 0, -pkin(3) * t79 - pkin(11) * t255 - t56 * t237 + t317 * t90, -pkin(3) * t80 + pkin(11) * t254 + t56 * t233 + t316 * t90, t155 * t94 + t196 * t40, t155 * t258 - t156 * t94 - t196 * t41 + t256 * t40, t118 * t196 + t149 * t155, t118 * t256 - t149 * t156, 0, -t117 * t149 - t118 * t165 + t156 * t63 + t222 * t41 - t256 * t42 - t258 * t309, t116 * t149 - t118 * t166 + t155 * t63 + t196 * t42 + t222 * t40 + t309 * t94, t67 * t339 + (-t313 * t67 + t344) * t196, t271 * t155 + (-t345 - t236 * t29 + (-t231 * t270 - t236 * t67) * qJD(6)) * t196, t67 * t156 - t196 * t262 - t256 * t28 - t258 * t339, t156 * t270 - t196 * t31 + t256 * t29 + t258 * t340, -t156 * t258 - t256 * t41, t10 * t156 + t108 * t41 - t117 * t270 + t165 * t29 + t196 * t354 - t2 * t256 + t22 * t340 - t258 * t58, -t1 * t256 - t109 * t41 - t11 * t156 + t117 * t67 + t165 * t28 - t196 * t291 + t22 * t339 - t258 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t294, -t213, 0, -t183, t182, t160 * t233 + t188 * t316, t160 * t237 - t233 * t161 + (-t187 * t237 - t188 * t233) * qJD(4), t248 * t228, t249 * t228, 0, -pkin(3) * t161 + t178 * t317 - t183 * t237 - t248 * t350, -pkin(3) * t160 + t178 * t316 + t183 * t233 - t249 * t350, t132 * t155 + t196 * t95, -t132 * t156 + t155 * t257 - t196 * t96 + t256 * t95 (-t155 * t238 + t196 * t319) * t228 (t156 * t238 + t256 * t319) * t228, 0, -t257 * t309 - t128 * t256 + t137 * t156 + t222 * t96 + (t117 * t238 - t165 * t319) * t228, t132 * t309 + t128 * t196 + t137 * t155 + t222 * t95 + (-t116 * t238 - t166 * t319) * t228, -t125 * t252 + t336 * t64, t269 * t155 + (-t342 - t236 * t65 + (t124 * t231 - t125 * t236) * qJD(6)) * t196, t125 * t156 + t252 * t257 - t256 * t64 + t336 * t96, -t124 * t156 + t253 * t257 + t256 * t65 - t337 * t96, -t156 * t257 - t256 * t96, t108 * t96 + t117 * t124 - t13 * t256 + t156 * t48 + t165 * t65 + t196 * t348 - t257 * t58 + t340 * t72, -t109 * t96 + t117 * t125 - t12 * t256 - t156 * t49 + t165 * t64 - t196 * t287 - t257 * t57 + t339 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t233 * t316, 0.2e1 * (-t233 ^ 2 + t237 ^ 2) * qJD(4), 0, 0, 0, t233 * t312, t237 * t312, 0.2e1 * t338, 0.2e1 * t155 * t256 - 0.2e1 * t156 * t196, 0, 0, 0, 0.2e1 * t156 * t222 - 0.2e1 * t256 * t309, 0.2e1 * t155 * t222 + 0.2e1 * t196 * t309, -0.2e1 * t193 * t292 + 0.2e1 * t227 * t338, t155 * t284 + 0.2e1 * t193 * t283, 0.2e1 * t156 * t336 + 0.2e1 * t252 * t256, -0.2e1 * t156 * t337 + 0.2e1 * t253 * t256, -0.2e1 * t256 * t156, 0.2e1 * t108 * t156 + 0.2e1 * t117 * t337 + 0.2e1 * t165 * t253 - 0.2e1 * t256 * t58, -0.2e1 * t109 * t156 + 0.2e1 * t117 * t336 - 0.2e1 * t165 * t252 - 0.2e1 * t256 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t79, t118, t27, t26, 0, 0, t40, -t41, t118, t118 * t310 + (-t46 + (-t45 - t352) * t232) * qJD(5) - t289 (-t118 * t232 - t149 * t290) * pkin(4) + t6, t19, t8, t31, -t262, 0, t221 * t29 - t31 * t220 + (-t232 * t270 + t258 * t299) * t346 + t291, t221 * t28 + t262 * t220 + (t232 * t67 + t258 * t298) * t346 + t354; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t160, -t161, t213, t99, t98, 0, 0, t95, -t96, t213, t353 * t282 + (-t107 + (-t106 + t311) * t232) * qJD(5) - t288 (-t232 * t319 + t238 * t290) * t228 * pkin(4) + t38, t53, t32, t69, -t261, 0, t221 * t65 - t69 * t220 + (t124 * t232 + t257 * t299) * t346 + t287, t221 * t64 + t261 * t220 + (t125 * t232 + t257 * t298) * t346 + t348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t316, -t317, 0, -pkin(11) * t316, pkin(11) * t317, 0, 0, t155, -t156, 0, -t117, t116, t105, t101, t127, -t251, 0, t153 + (-qJD(6) * t267 - t117) * t236 + t242 * t231, t236 * t242 + t267 * t313 + t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t306, -0.2e1 * t281, t212, t194, 0, 0, 0, 0.2e1 * t265, 0.2e1 * t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t41, t118, t7, t6, t19, t8, t31, -t262, 0, -pkin(5) * t29 - pkin(13) * t31 + t291, -pkin(5) * t28 + pkin(13) * t262 + t354; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, -t96, t213, t39, t38, t53, t32, t69, -t261, 0, -pkin(5) * t65 - pkin(13) * t69 + t287, -pkin(5) * t64 + pkin(13) * t261 + t348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, -t156, 0, -t117, t116, t105, t101, t127, -t251, 0, t153 + (-pkin(5) * t155 - pkin(13) * t156) * t231 + (-t117 + (-pkin(5) * t196 + pkin(13) * t256) * qJD(6)) * t236, pkin(5) * t252 + pkin(13) * t251 + t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t306, -t281, t212, t194, 0, 0, 0, t265 - t308, -t307 + t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t212, t194, 0, 0, 0, -0.2e1 * t308, -0.2e1 * t307; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t29, t41, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, -t65, t96, t13, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t252, -t253, t156, t58, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, -t313, 0, -t220 * t223 - t231 * t281, t220 * t313 - t236 * t281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223, -t313, 0, -pkin(13) * t223, pkin(13) * t313; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
