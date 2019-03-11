% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:52:57
% EndTime: 2019-03-09 03:53:11
% DurationCPUTime: 6.43s
% Computational Cost: add. (8478->479), mult. (20497->626), div. (0->0), fcn. (16880->18), ass. (0->239)
t239 = cos(pkin(10));
t353 = cos(qJ(3));
t298 = t353 * t239;
t216 = qJD(1) * t298;
t237 = sin(pkin(10));
t243 = sin(qJ(3));
t328 = t243 * t237;
t297 = qJD(1) * t328;
t181 = -t216 + t297;
t175 = qJD(5) + t181;
t171 = qJD(6) + t175;
t241 = sin(qJ(6));
t245 = cos(qJ(6));
t198 = t353 * t237 + t243 * t239;
t183 = t198 * qJD(1);
t236 = sin(pkin(11));
t238 = cos(pkin(11));
t151 = -t238 * qJD(3) + t236 * t183;
t153 = t236 * qJD(3) + t238 * t183;
t242 = sin(qJ(5));
t246 = cos(qJ(5));
t365 = t242 * t151 - t246 * t153;
t97 = t246 * t151 + t242 * t153;
t374 = t241 * t97 + t245 * t365;
t377 = t374 * t171;
t47 = -t241 * t365 + t245 * t97;
t376 = t47 * t171;
t375 = t97 * t175;
t195 = t242 * t236 - t246 * t238;
t318 = t175 * t195;
t197 = t246 * t236 + t242 * t238;
t185 = t197 * qJD(5);
t317 = t197 * t181 + t185;
t335 = t181 * t236;
t136 = t183 * pkin(3) + t181 * qJ(4);
t348 = pkin(7) + qJ(2);
t206 = t348 * t237;
t199 = qJD(1) * t206;
t208 = t348 * t239;
t200 = qJD(1) * t208;
t264 = t353 * t199 + t243 * t200;
t84 = t236 * t136 - t238 * t264;
t66 = pkin(8) * t335 + t84;
t373 = -t238 * qJD(4) + t66;
t235 = pkin(10) + qJ(3);
t227 = sin(t235);
t229 = cos(t235);
t244 = sin(qJ(1));
t247 = cos(qJ(1));
t280 = g(1) * t247 + g(2) * t244;
t259 = -g(3) * t229 + t227 * t280;
t306 = qJD(1) * qJD(2);
t355 = t348 * qJDD(1) + t306;
t168 = t355 * t237;
t169 = t355 * t239;
t295 = qJD(3) * t353;
t313 = qJD(3) * t243;
t261 = t353 * t168 + t243 * t169 - t199 * t313 + t200 * t295;
t81 = -qJDD(3) * pkin(3) + qJDD(4) + t261;
t254 = t81 - t259;
t372 = t175 * t365;
t143 = -t241 * t195 + t245 * t197;
t343 = qJD(6) * t143 - t318 * t241 + t317 * t245;
t369 = t374 * t47;
t338 = qJDD(1) * pkin(1);
t362 = g(1) * t244 - g(2) * t247;
t270 = -qJDD(2) + t338 + t362;
t368 = t374 ^ 2 - t47 ^ 2;
t309 = qJD(6) * t245;
t310 = qJD(6) * t241;
t288 = qJDD(1) * t353;
t304 = t239 * qJDD(1);
t303 = qJD(3) * t216 + t237 * t288 + t243 * t304;
t139 = -qJD(3) * t297 + t303;
t120 = -t238 * qJDD(3) + t236 * t139;
t121 = t236 * qJDD(3) + t238 * t139;
t311 = qJD(5) * t246;
t312 = qJD(5) * t242;
t33 = -t242 * t120 + t246 * t121 - t151 * t311 - t153 * t312;
t34 = -qJD(5) * t365 + t246 * t120 + t242 * t121;
t8 = -t241 * t34 + t245 * t33 - t97 * t309 + t310 * t365;
t367 = t8 + t376;
t221 = t239 * pkin(2) + pkin(1);
t204 = -qJD(1) * t221 + qJD(2);
t114 = t181 * pkin(3) - t183 * qJ(4) + t204;
t145 = -t243 * t199 + t353 * t200;
t135 = qJD(3) * qJ(4) + t145;
t72 = t238 * t114 - t236 * t135;
t41 = t181 * pkin(4) - t153 * pkin(8) + t72;
t73 = t236 * t114 + t238 * t135;
t53 = -t151 * pkin(8) + t73;
t23 = t242 * t41 + t246 * t53;
t13 = -t97 * pkin(9) + t23;
t11 = t13 * t310;
t234 = pkin(11) + qJ(5);
t230 = qJ(6) + t234;
t219 = sin(t230);
t322 = t247 * t219;
t220 = cos(t230);
t325 = t244 * t220;
t156 = -t229 * t325 + t322;
t321 = t247 * t220;
t326 = t244 * t219;
t158 = t229 * t321 + t326;
t350 = g(3) * t227;
t133 = -qJD(3) * pkin(3) + qJD(4) + t264;
t91 = t151 * pkin(4) + t133;
t45 = t97 * pkin(5) + t91;
t366 = g(1) * t158 - g(2) * t156 + t220 * t350 + t45 * t47 + t11;
t363 = t362 * t227;
t268 = t362 * t229;
t347 = pkin(8) + qJ(4);
t205 = t347 * t236;
t207 = t347 * t238;
t315 = -t242 * t205 + t246 * t207;
t361 = -t353 * t206 - t243 * t208;
t360 = qJ(2) * qJDD(1);
t155 = t229 * t326 + t321;
t157 = -t229 * t322 + t325;
t187 = t198 * qJD(3);
t305 = t237 * qJDD(1);
t277 = -t239 * t288 + t243 * t305;
t140 = qJD(1) * t187 + t277;
t137 = qJDD(5) + t140;
t203 = -qJDD(1) * t221 + qJDD(2);
t63 = t140 * pkin(3) - t139 * qJ(4) - t183 * qJD(4) + t203;
t265 = -t243 * t168 + t353 * t169;
t76 = qJDD(3) * qJ(4) + (qJD(4) - t264) * qJD(3) + t265;
t28 = -t236 * t76 + t238 * t63;
t18 = t140 * pkin(4) - t121 * pkin(8) + t28;
t29 = t236 * t63 + t238 * t76;
t26 = -t120 * pkin(8) + t29;
t292 = t246 * t18 - t242 * t26;
t253 = -qJD(5) * t23 + t292;
t2 = t137 * pkin(5) - t33 * pkin(9) + t253;
t267 = t242 * t18 + t246 * t26 + t41 * t311 - t312 * t53;
t3 = -t34 * pkin(9) + t267;
t302 = t245 * t2 - t241 * t3;
t22 = -t242 * t53 + t246 * t41;
t12 = pkin(9) * t365 + t22;
t10 = t175 * pkin(5) + t12;
t340 = t245 * t13;
t5 = t241 * t10 + t340;
t359 = -g(1) * t157 + g(2) * t155 - qJD(6) * t5 + t219 * t350 + t374 * t45 + t302;
t9 = -qJD(6) * t374 + t241 * t33 + t245 * t34;
t358 = -t9 - t377;
t132 = qJDD(6) + t137;
t142 = t245 * t195 + t241 * t197;
t344 = -qJD(6) * t142 - t317 * t241 - t318 * t245;
t357 = -t143 * t132 - t171 * t344;
t356 = -t197 * t137 + t175 * t318;
t257 = -t280 * t229 - t350;
t271 = t236 * qJD(4) + qJD(5) * t207;
t352 = pkin(8) * t238;
t83 = t238 * t136 + t236 * t264;
t54 = t183 * pkin(4) + t181 * t352 + t83;
t354 = t205 * t311 + t373 * t246 + (t271 + t54) * t242;
t176 = t181 ^ 2;
t263 = t298 - t328;
t331 = t198 * t238;
t138 = -pkin(3) * t263 - t198 * qJ(4) - t221;
t150 = -t243 * t206 + t353 * t208;
t88 = t238 * t138 - t236 * t150;
t62 = -pkin(4) * t263 - pkin(8) * t331 + t88;
t332 = t198 * t236;
t89 = t236 * t138 + t238 * t150;
t75 = -pkin(8) * t332 + t89;
t345 = t242 * t62 + t246 * t75;
t342 = t183 * t47;
t341 = t183 * t97;
t339 = t374 * t183;
t337 = t365 * t183;
t186 = t263 * qJD(3);
t334 = t186 * t236;
t330 = t236 * t140;
t329 = t238 * t140;
t226 = sin(t234);
t324 = t244 * t226;
t228 = cos(t234);
t323 = t244 * t228;
t320 = t247 * t226;
t319 = t247 * t228;
t105 = t187 * pkin(3) - t186 * qJ(4) - t198 * qJD(4);
t115 = t263 * qJD(2) + qJD(3) * t361;
t61 = t236 * t105 + t238 * t115;
t314 = t237 ^ 2 + t239 ^ 2;
t307 = -qJD(4) + t133;
t104 = -pkin(4) * t335 + t145;
t222 = -t238 * pkin(4) - pkin(3);
t296 = t317 * pkin(5) - t104;
t293 = qJD(6) * t10 + t3;
t60 = t238 * t105 - t236 * t115;
t38 = t187 * pkin(4) - t186 * t352 + t60;
t44 = -pkin(8) * t334 + t61;
t290 = -t242 * t44 + t246 * t38;
t289 = -t242 * t75 + t246 * t62;
t287 = t314 * qJD(1) ^ 2;
t285 = -t246 * t205 - t242 * t207;
t116 = t198 * qJD(2) - t206 * t313 + t208 * t295;
t284 = 0.2e1 * t314;
t283 = -t142 * t132 - t343 * t171;
t123 = -t197 * pkin(9) + t285;
t282 = t317 * pkin(9) - qJD(6) * t123 + t354;
t124 = -t195 * pkin(9) + t315;
t52 = t246 * t54;
t281 = t183 * pkin(5) - t318 * pkin(9) + t197 * qJD(4) + t315 * qJD(5) + qJD(6) * t124 - t242 * t66 + t52;
t278 = -t195 * t137 - t317 * t175;
t276 = -t29 * t236 - t28 * t238;
t275 = -t72 * t236 + t73 * t238;
t90 = pkin(4) * t334 + t116;
t128 = t197 * t198;
t129 = t195 * t198;
t79 = t245 * t128 - t241 * t129;
t80 = -t241 * t128 - t245 * t129;
t117 = pkin(4) * t332 - t361;
t269 = t229 * pkin(3) + t227 * qJ(4) + t221;
t266 = t242 * t38 + t246 * t44 + t62 * t311 - t312 * t75;
t260 = t270 + t338;
t256 = t133 * t186 + t81 * t198 - t280;
t251 = t284 * t306 - t280;
t42 = t120 * pkin(4) + t81;
t164 = t195 * pkin(5) + t222;
t163 = t229 * t319 + t324;
t162 = -t229 * t320 + t323;
t161 = -t229 * t323 + t320;
t160 = t229 * t324 + t319;
t82 = t128 * pkin(5) + t117;
t78 = t186 * t197 + t311 * t331 - t312 * t332;
t77 = -t185 * t198 - t186 * t195;
t35 = t78 * pkin(5) + t90;
t27 = -t128 * pkin(9) + t345;
t24 = -pkin(5) * t263 + t129 * pkin(9) + t289;
t20 = qJD(6) * t80 + t241 * t77 + t245 * t78;
t19 = -qJD(6) * t79 - t241 * t78 + t245 * t77;
t14 = t34 * pkin(5) + t42;
t7 = -t78 * pkin(9) + t266;
t6 = t187 * pkin(5) - t77 * pkin(9) - qJD(5) * t345 + t290;
t4 = t245 * t10 - t241 * t13;
t1 = [qJDD(1), t362, t280, t260 * t239, -t260 * t237, t284 * t360 + t251, pkin(1) * t270 + (t314 * t360 + t251) * qJ(2), t139 * t198 + t183 * t186, t139 * t263 - t198 * t140 - t186 * t181 - t183 * t187, t186 * qJD(3) + t198 * qJDD(3), -t187 * qJD(3) + qJDD(3) * t263, 0, -t116 * qJD(3) + qJDD(3) * t361 - t221 * t140 + t204 * t187 - t203 * t263 + t268, -t115 * qJD(3) - t150 * qJDD(3) - t221 * t139 + t204 * t186 + t203 * t198 - t363, t116 * t151 - t120 * t361 + t88 * t140 + t60 * t181 + t72 * t187 + t236 * t256 + t238 * t268 - t263 * t28, t116 * t153 - t121 * t361 - t89 * t140 - t61 * t181 - t73 * t187 - t236 * t268 + t238 * t256 + t263 * t29, -t89 * t120 - t88 * t121 - t61 * t151 - t60 * t153 + t363 + t276 * t198 + (-t236 * t73 - t238 * t72) * t186, t133 * t116 - t81 * t361 + t28 * t88 + t29 * t89 + t72 * t60 + t73 * t61 + (-g(1) * t348 - g(2) * t269) * t247 + (g(1) * t269 - g(2) * t348) * t244, -t33 * t129 - t365 * t77, -t33 * t128 + t129 * t34 + t365 * t78 - t77 * t97, -t129 * t137 + t77 * t175 - t187 * t365 - t263 * t33, -t128 * t137 - t78 * t175 - t97 * t187 + t263 * t34, -t137 * t263 + t175 * t187, t290 * t175 + t289 * t137 - t292 * t263 + t22 * t187 + t90 * t97 + t117 * t34 + t42 * t128 + t91 * t78 - g(1) * t161 - g(2) * t163 + (-t175 * t345 + t23 * t263) * qJD(5), -g(1) * t160 - g(2) * t162 + t117 * t33 - t42 * t129 - t345 * t137 - t266 * t175 - t23 * t187 + t263 * t267 - t365 * t90 + t91 * t77, -t19 * t374 + t8 * t80, -t19 * t47 + t20 * t374 - t8 * t79 - t80 * t9, t80 * t132 + t19 * t171 - t187 * t374 - t263 * t8, -t79 * t132 - t20 * t171 - t47 * t187 + t263 * t9, -t132 * t263 + t171 * t187 (-t241 * t7 + t245 * t6) * t171 + (t245 * t24 - t241 * t27) * t132 - t302 * t263 + t4 * t187 + t35 * t47 + t82 * t9 + t14 * t79 + t45 * t20 - g(1) * t156 - g(2) * t158 + ((-t241 * t24 - t245 * t27) * t171 + t5 * t263) * qJD(6), -g(1) * t155 - g(2) * t157 - t11 * t263 + t14 * t80 - t5 * t187 + t45 * t19 - t35 * t374 + t82 * t8 + (-(-qJD(6) * t27 + t6) * t171 - t24 * t132 + t2 * t263) * t241 + (-(qJD(6) * t24 + t7) * t171 - t27 * t132 + t293 * t263) * t245; 0, 0, 0, -t304, t305, -t287, -qJ(2) * t287 - t270, 0, 0, 0, 0, 0, 0.2e1 * t183 * qJD(3) + t277 (-t181 - t297) * qJD(3) + t303, -t183 * t151 - t236 * t176 + t329, -t183 * t153 - t238 * t176 - t330, -t236 * t120 - t238 * t121 + (-t151 * t238 + t153 * t236) * t181, -t133 * t183 + t181 * t275 - t276 - t362, 0, 0, 0, 0, 0, t278 - t341, t337 + t356, 0, 0, 0, 0, 0, t283 - t342, t339 + t357; 0, 0, 0, 0, 0, 0, 0, t183 * t181, t183 ^ 2 - t176 (t181 - t297) * qJD(3) + t303, -t277, qJDD(3), t145 * qJD(3) - t204 * t183 + t259 - t261, t204 * t181 - t257 - t265, -qJ(4) * t330 - pkin(3) * t120 - t145 * t151 - t72 * t183 + (t236 * t307 - t83) * t181 - t254 * t238, -qJ(4) * t329 - pkin(3) * t121 - t145 * t153 + t73 * t183 + (t238 * t307 + t84) * t181 + t254 * t236, t84 * t151 + t83 * t153 + (-qJ(4) * t120 - qJD(4) * t151 - t181 * t72 + t29) * t238 + (qJ(4) * t121 + qJD(4) * t153 - t181 * t73 - t28) * t236 + t257, -t133 * t145 - t72 * t83 - t73 * t84 + t275 * qJD(4) - t254 * pkin(3) + (-t28 * t236 + t29 * t238 + t257) * qJ(4), t33 * t197 + t318 * t365, -t33 * t195 - t197 * t34 + t317 * t365 + t318 * t97, t337 - t356, t278 + t341, -t175 * t183, t285 * t137 + t222 * t34 + t42 * t195 - t22 * t183 - t104 * t97 + t317 * t91 + (-t52 - t271 * t246 + (qJD(5) * t205 + t373) * t242) * t175 + t259 * t228, t104 * t365 - t315 * t137 + t175 * t354 + t23 * t183 + t42 * t197 + t222 * t33 - t226 * t259 - t318 * t91, t8 * t143 - t344 * t374, -t8 * t142 - t143 * t9 + t343 * t374 - t344 * t47, t339 - t357, t283 + t342, -t171 * t183 (t245 * t123 - t241 * t124) * t132 + t164 * t9 + t14 * t142 - t4 * t183 + t296 * t47 + t343 * t45 + (t241 * t282 - t245 * t281) * t171 + t259 * t220 -(t241 * t123 + t245 * t124) * t132 + t164 * t8 + t14 * t143 + t5 * t183 - t296 * t374 + t344 * t45 + (t241 * t281 + t245 * t282) * t171 - t259 * t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153 * t181 + t120, -t151 * t181 + t121, -t151 ^ 2 - t153 ^ 2, t151 * t73 + t153 * t72 + t254, 0, 0, 0, 0, 0, t34 - t372, t33 - t375, 0, 0, 0, 0, 0, t9 - t377, t8 - t376; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t365 * t97, t365 ^ 2 - t97 ^ 2, t33 + t375, -t34 - t372, t137, -g(1) * t162 + g(2) * t160 + t23 * t175 + t226 * t350 + t365 * t91 + t253, g(1) * t163 - g(2) * t161 + t22 * t175 + t228 * t350 + t91 * t97 - t267, -t369, t368, t367, t358, t132 -(-t241 * t12 - t340) * t171 + (t245 * t132 - t171 * t310 + t365 * t47) * pkin(5) + t359 (-t13 * t171 - t2) * t241 + (t12 * t171 - t293) * t245 + (-t241 * t132 - t171 * t309 - t365 * t374) * pkin(5) + t366; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t369, t368, t367, t358, t132, t5 * t171 + t359, t4 * t171 - t241 * t2 - t245 * t293 + t366;];
tau_reg  = t1;
