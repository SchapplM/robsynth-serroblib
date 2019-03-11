% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:16:31
% EndTime: 2019-03-09 02:16:40
% DurationCPUTime: 4.79s
% Computational Cost: add. (6819->484), mult. (13554->558), div. (0->0), fcn. (9527->10), ass. (0->240)
t169 = cos(pkin(9));
t168 = sin(pkin(9));
t175 = cos(qJ(4));
t308 = sin(qJ(4));
t114 = t168 * t175 + t169 * t308;
t244 = qJD(4) * t308;
t231 = qJD(1) * t244;
t190 = qJDD(1) * t114 - t168 * t231;
t261 = qJD(4) * t175;
t243 = qJD(1) * t261;
t74 = t169 * t243 + t190;
t71 = qJDD(5) + t74;
t105 = t114 * qJD(1);
t334 = qJD(5) + t105;
t238 = t334 ^ 2;
t172 = sin(qJ(5));
t174 = cos(qJ(5));
t171 = -pkin(1) - qJ(3);
t135 = qJD(1) * t171 + qJD(2);
t240 = -pkin(7) * qJD(1) + t135;
t95 = t240 * t168;
t96 = t240 * t169;
t59 = t175 * t95 + t308 * t96;
t56 = qJD(4) * pkin(8) + t59;
t273 = t169 * t175;
t246 = qJD(1) * t273;
t247 = t308 * t168;
t107 = -qJD(1) * t247 + t246;
t150 = qJD(1) * qJ(2) + qJD(3);
t156 = t168 * pkin(3);
t120 = qJD(1) * t156 + t150;
t57 = pkin(4) * t105 - pkin(8) * t107 + t120;
t25 = t172 * t57 + t174 * t56;
t16 = qJ(6) * t334 + t25;
t336 = t16 * t334;
t24 = -t172 * t56 + t174 * t57;
t335 = t24 * t334;
t304 = t25 * t334;
t258 = t174 * qJD(4);
t84 = t107 * t172 - t258;
t237 = t334 * t84;
t254 = t169 * qJDD(1);
t201 = -qJDD(1) * t247 - t168 * t243 - t169 * t231 + t175 * t254;
t260 = qJD(5) * t172;
t44 = -qJD(5) * t258 - qJDD(4) * t172 + t107 * t260 - t174 * t201;
t86 = qJD(4) * t172 + t107 * t174;
t78 = t86 * t260;
t329 = -t174 * t44 - t78;
t115 = -t247 + t273;
t109 = t114 * qJD(4);
t233 = qJD(5) * t114 + qJD(1);
t278 = t114 * t174;
t110 = -t168 * t244 + t169 * t261;
t279 = t110 * t174;
t333 = t109 * t86 + t115 * t44 - t71 * t278 + (t172 * t233 - t279) * t334;
t289 = t174 * t84;
t292 = t172 * t86;
t216 = t289 + t292;
t288 = t174 * t86;
t45 = qJD(5) * t86 - qJDD(4) * t174 + t172 * t201;
t290 = t174 * t45;
t293 = t172 * t44;
t332 = (qJD(5) * (t172 * t84 - t288) - t290 + t293) * t115 + t216 * t109;
t164 = pkin(9) + qJ(4);
t152 = cos(t164);
t306 = g(3) * t152;
t259 = qJD(5) * t174;
t64 = t172 * t71;
t207 = t259 * t334 + t64;
t66 = t174 * t71;
t331 = -t260 * t334 + t66;
t330 = -t334 * t86 + t45;
t300 = -t172 * t45 - t259 * t84;
t328 = (-t289 + t292) * t105 + t300 - t329;
t162 = t168 ^ 2;
t163 = t169 ^ 2;
t263 = t162 + t163;
t327 = t135 * t263;
t176 = cos(qJ(1));
t160 = g(2) * t176;
t173 = sin(qJ(1));
t161 = g(1) * t173;
t324 = t161 - t160;
t326 = t152 * t324;
t165 = qJDD(1) * qJ(2);
t166 = qJD(1) * qJD(2);
t323 = t165 + t166;
t122 = qJDD(3) + t323;
t228 = g(1) * t176 + g(2) * t173;
t197 = -t228 + t122;
t301 = -pkin(7) + t171;
t117 = t301 * t168;
t118 = t301 * t169;
t76 = t117 * t308 - t118 * t175;
t315 = -qJD(1) * qJD(3) + qJDD(1) * t171;
t116 = qJDD(2) + t315;
t235 = -pkin(7) * qJDD(1) + t116;
t91 = t235 * t168;
t92 = t235 * t169;
t253 = t175 * t91 + t261 * t96 + t308 * t92;
t27 = -t244 * t95 + t253;
t22 = qJDD(4) * pkin(8) + t27;
t255 = t168 * qJDD(1);
t145 = pkin(3) * t255;
t112 = t145 + t122;
t35 = pkin(4) * t74 - pkin(8) * t201 + t112;
t241 = t172 * t22 - t174 * t35 + t259 * t56 + t260 * t57;
t311 = pkin(5) * t71;
t2 = qJDD(6) + t241 - t311;
t322 = -t2 + t336;
t321 = -t241 + t304;
t94 = t308 * t95;
t58 = t175 * t96 - t94;
t55 = -qJD(4) * pkin(4) - t58;
t26 = pkin(5) * t84 - qJ(6) * t86 + t55;
t310 = pkin(8) * t71;
t320 = t26 * t334 - t310;
t281 = t109 * t172;
t319 = t110 * t84 + t114 * t45 + t115 * t207 - t281 * t334;
t318 = -t107 * t109 + t115 * t201;
t142 = qJ(2) + t156;
t72 = pkin(4) * t114 - pkin(8) * t115 + t142;
t77 = t117 * t175 + t118 * t308;
t299 = t172 * t72 + t174 * t77;
t52 = -qJD(3) * t114 - qJD(4) * t76;
t69 = pkin(4) * t110 + pkin(8) * t109 + qJD(2);
t11 = -qJD(5) * t299 - t172 * t52 + t174 * t69;
t314 = t86 ^ 2;
t313 = t107 ^ 2;
t312 = 0.2e1 * t166;
t309 = pkin(8) * t86;
t151 = sin(t164);
t307 = g(3) * t151;
t305 = g(3) * t174;
t303 = t86 * t84;
t73 = pkin(4) * t107 + pkin(8) * t105;
t34 = t172 * t73 + t174 * t58;
t298 = qJ(6) * t71;
t297 = t105 * t334;
t296 = t107 * t84;
t295 = t107 * t86;
t294 = t110 * t86;
t287 = t334 * t107;
t224 = pkin(5) * t172 - qJ(6) * t174;
t286 = -qJD(6) * t172 + t224 * t334 - t59;
t285 = -qJD(4) * t109 + qJDD(4) * t115;
t284 = pkin(1) * qJDD(1);
t283 = t107 * t105;
t280 = t109 * t174;
t277 = t151 * t173;
t276 = t151 * t176;
t275 = t152 * t173;
t274 = t152 * t176;
t272 = t172 * t173;
t271 = t172 * t176;
t270 = t173 * t174;
t269 = t176 * t174;
t268 = qJD(6) - t24;
t251 = g(2) * t274;
t267 = t151 * t305 + t174 * t251;
t266 = g(1) * t274 + g(2) * t275;
t265 = pkin(1) * t176 + qJ(2) * t173;
t262 = qJD(4) * t105;
t252 = pkin(8) * qJD(5) * t334;
t133 = g(2) * t276;
t250 = t84 ^ 2 - t314;
t28 = t175 * t92 - t244 * t96 - t261 * t95 - t308 * t91;
t23 = -qJDD(4) * pkin(4) - t28;
t245 = t23 - t307;
t158 = t176 * qJ(2);
t242 = -pkin(1) * t173 + t158;
t239 = qJD(5) * t84 - t44;
t236 = t263 * t116;
t234 = qJDD(2) - t284;
t101 = t151 * t271 + t270;
t99 = t151 * t272 - t269;
t230 = -g(1) * t101 - g(2) * t99;
t100 = t151 * t270 + t271;
t102 = t151 * t269 - t272;
t229 = -g(1) * t102 - g(2) * t100;
t226 = t239 * pkin(8);
t225 = pkin(5) * t174 + qJ(6) * t172;
t15 = -pkin(5) * t334 + t268;
t223 = t15 * t174 - t16 * t172;
t222 = t15 * t172 + t16 * t174;
t221 = t172 * t25 + t174 * t24;
t220 = t172 * t24 - t174 * t25;
t33 = -t172 * t58 + t174 * t73;
t41 = -t172 * t77 + t174 * t72;
t214 = t105 * t110 + t114 * t74;
t213 = t174 * t297 + t207;
t212 = -t172 * t297 + t331;
t170 = -pkin(7) - qJ(3);
t211 = t156 * t176 + t170 * t173 + t242;
t210 = -qJD(4) * t110 - qJDD(4) * t114;
t209 = pkin(4) + t225;
t208 = t156 * t173 - t170 * t176 + t265;
t206 = -g(1) * (pkin(4) * t275 + pkin(8) * t277) - pkin(8) * t306;
t3 = t172 * t35 + t174 * t22 + t259 * t57 - t260 * t56;
t10 = t172 * t69 + t174 * t52 + t259 * t72 - t260 * t77;
t205 = t334 * t55 - t310;
t204 = -g(1) * t275 - t252;
t203 = -g(1) * t277 + t133 - t306;
t200 = -pkin(8) * t290 + t203;
t199 = t172 * t237 - t290;
t198 = -t212 - t296;
t194 = g(1) * t99 - g(2) * t101 + t172 * t306 - t241;
t193 = t307 - t326;
t192 = pkin(4) * t276 - pkin(8) * t274 + t211;
t189 = pkin(4) * t277 - pkin(8) * t275 + t208;
t186 = t197 + t323;
t185 = -t115 * t300 - t281 * t84;
t1 = qJD(6) * t334 + t298 + t3;
t184 = qJD(5) * t223 + t1 * t174 + t2 * t172;
t183 = -qJD(5) * t221 + t172 * t241 + t3 * t174;
t182 = t109 * t58 - t110 * t59 - t114 * t27 - t115 * t28 + t324;
t181 = t26 * t86 + qJDD(6) - t194;
t53 = qJD(3) * t115 + qJD(4) * t77;
t180 = -g(1) * t100 + g(2) * t102 - t152 * t305 + t3;
t179 = -t114 * t64 - t115 * t45 + t109 * t84 + (-t110 * t172 - t174 * t233) * t334;
t178 = -t45 * t278 - t84 * t279 + t233 * t288 + (qJD(1) * t84 + t114 * t239 + t294) * t172;
t177 = qJD(1) ^ 2;
t104 = t105 ^ 2;
t47 = pkin(5) * t86 + qJ(6) * t84;
t46 = t115 * t224 + t76;
t38 = t110 * t334 + t114 * t71;
t32 = -pkin(5) * t114 - t41;
t31 = qJ(6) * t114 + t299;
t21 = -pkin(5) * t107 - t33;
t20 = qJ(6) * t107 + t34;
t17 = -t44 + t237;
t14 = t213 - t295;
t13 = -t224 * t109 + (qJD(5) * t225 - qJD(6) * t174) * t115 + t53;
t12 = t288 * t334 - t293;
t9 = t115 * t329 - t280 * t86;
t8 = -pkin(5) * t110 - t11;
t7 = qJ(6) * t110 + qJD(6) * t114 + t10;
t6 = -t114 * t44 + t115 * t331 - t280 * t334 + t294;
t5 = pkin(5) * t45 + qJ(6) * t44 - qJD(6) * t86 + t23;
t4 = [0, 0, 0, 0, 0, qJDD(1), t324, t228, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t324 - 0.2e1 * t284, 0.2e1 * t165 + t312 - t228, -t234 * pkin(1) - g(1) * t242 - g(2) * t265 + (t165 + t312) * qJ(2), t163 * qJDD(1), -0.2e1 * t168 * t254, 0, t162 * qJDD(1), 0, 0, t186 * t168, t186 * t169, t324 + t263 * (-t116 - t315) t122 * qJ(2) + t150 * qJD(2) - g(1) * (t171 * t173 + t158) - g(2) * (qJ(3) * t176 + t265) + t171 * t236 - qJD(3) * t327, t318, t105 * t109 - t107 * t110 - t114 * t201 - t115 * t74, t285, t214, t210, 0, qJD(2) * t105 - qJD(4) * t53 - qJDD(4) * t76 + t110 * t120 + t112 * t114 + t142 * t74 - t151 * t228, qJD(2) * t107 - qJD(4) * t52 - qJDD(4) * t77 - t109 * t120 + t112 * t115 + t142 * t201 - t266, -t105 * t52 + t107 * t53 + t201 * t76 - t74 * t77 + t182, -g(1) * t211 - g(2) * t208 + qJD(2) * t120 + t112 * t142 + t27 * t77 - t28 * t76 + t52 * t59 - t53 * t58, t9, t332, t6, t185, -t319, t38, -t55 * t281 + t11 * t334 + t110 * t24 - t114 * t241 + t41 * t71 + t45 * t76 + t53 * t84 + (t172 * t23 + t259 * t55) * t115 + t229, -t55 * t280 - t10 * t334 - t110 * t25 - t114 * t3 - t299 * t71 - t44 * t76 + t53 * t86 + (t174 * t23 - t260 * t55) * t115 - t230, -t10 * t84 - t11 * t86 + t41 * t44 - t299 * t45 + t221 * t109 + (qJD(5) * t220 - t172 * t3 + t174 * t241) * t115 + t266, -g(1) * t192 - g(2) * t189 + t25 * t10 + t24 * t11 + t23 * t76 - t241 * t41 + t299 * t3 + t55 * t53, t9, t6, -t332, t38, t319, t185, -t26 * t281 - t110 * t15 - t114 * t2 + t13 * t84 - t32 * t71 + t45 * t46 - t8 * t334 + (t172 * t5 + t259 * t26) * t115 + t229, -t31 * t45 - t32 * t44 - t7 * t84 + t8 * t86 - t223 * t109 + (-qJD(5) * t222 - t1 * t172 + t174 * t2) * t115 + t266, t26 * t280 + t1 * t114 + t110 * t16 - t13 * t86 + t31 * t71 + t44 * t46 + t7 * t334 + (-t174 * t5 + t26 * t260) * t115 + t230, t1 * t31 + t16 * t7 + t5 * t46 + t26 * t13 + t2 * t32 + t15 * t8 - g(1) * (pkin(5) * t102 + qJ(6) * t101 + t192) - g(2) * (pkin(5) * t100 + qJ(6) * t99 + t189); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t177, -qJ(2) * t177 + t234 - t324, 0, 0, 0, 0, 0, 0, -t177 * t168, -t177 * t169, -t263 * qJDD(1), -qJD(1) * t150 + t236 - t324, 0, 0, 0, 0, 0, 0, -qJD(1) * t105 + t285, -qJD(1) * t107 + t210, -t214 - t318, -qJD(1) * t120 - t182, 0, 0, 0, 0, 0, 0, t179, t333, t178, -qJD(1) * t221 + t109 * t55 - t110 * t220 + t114 * t183 - t115 * t23 - t324, 0, 0, 0, 0, 0, 0, t179, t178, -t333, qJD(1) * t223 + t109 * t26 + t110 * t222 + t114 * t184 - t115 * t5 - t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t255, t254, -t263 * t177, qJD(1) * t327 + t197, 0, 0, 0, 0, 0, 0 (t107 + t246) * qJD(4) + t190, t201 - t262, -t104 - t313, t105 * t59 + t107 * t58 + t145 + t197, 0, 0, 0, 0, 0, 0, t212 - t296, -t174 * t238 - t295 - t64, t328, -t107 * t55 + t321 * t174 + (t3 - t335) * t172 - t228, 0, 0, 0, 0, 0, 0, -t172 * t238 - t296 + t66, t328, t213 + t295, -t107 * t26 + t322 * t174 + (t15 * t334 + t1) * t172 - t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t283, -t104 + t313, t201 + t262, -t283 (t107 - t246) * qJD(4) - t190, qJDD(4), qJD(4) * t59 - t107 * t120 + t193 + t28, t120 * t105 + (t94 + t58) * qJD(4) - t203 - t253, 0, 0, t12, -t105 * t216 + t300 + t329, t14, t199, -t198, -t287, -pkin(4) * t45 - t107 * t24 - t33 * t334 - t59 * t84 + (t204 - t23) * t174 + t205 * t172 + t267, pkin(4) * t44 + t107 * t25 + t34 * t334 - t59 * t86 + t205 * t174 + (t245 + t252 + t326) * t172, t33 * t86 + t34 * t84 + (-t105 * t24 + t3 + (-t24 + t309) * qJD(5)) * t174 + (t226 - t321) * t172 + t200, -t25 * t34 - t24 * t33 - t55 * t59 + (-t245 + t251) * pkin(4) + (t183 + t133) * pkin(8) + t206, t12, t14, t78 + (t105 * t86 + t45) * t172 + (t44 + t237) * t174, -t287, t198, t199, t107 * t15 - t209 * t45 + t21 * t334 + t286 * t84 + (t204 - t5) * t174 + t320 * t172 + t267, t20 * t84 - t21 * t86 + (t105 * t15 + t1 + (t15 + t309) * qJD(5)) * t174 + (t226 - t322) * t172 + t200, -t107 * t16 - t209 * t44 - t20 * t334 - t286 * t86 - t320 * t174 + (t193 - t5 - t252) * t172, -t16 * t20 - t15 * t21 + t286 * t26 - t161 * t225 * t152 + (t184 + t133) * pkin(8) + t206 + (t152 * t160 + t307 - t5) * t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t303, -t250, t17, -t303, -t330, t71, -t55 * t86 + t194 + t304, t55 * t84 - t180 + t335, 0, 0, t303, t17, t250, t71, t330, -t303, -t47 * t84 - t181 + t304 + 0.2e1 * t311, pkin(5) * t44 - qJ(6) * t45 + (t16 - t25) * t86 + (t15 - t268) * t84, 0.2e1 * t298 - t26 * t84 + t47 * t86 + (0.2e1 * qJD(6) - t24) * t334 + t180, t1 * qJ(6) - t2 * pkin(5) - t26 * t47 - t15 * t25 - g(1) * (-pkin(5) * t99 + qJ(6) * t100) - g(2) * (pkin(5) * t101 - qJ(6) * t102) + t268 * t16 + t224 * t306; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t303 - t71, t17, -t314 - t238, t181 - t311 - t336;];
tau_reg  = t4;
