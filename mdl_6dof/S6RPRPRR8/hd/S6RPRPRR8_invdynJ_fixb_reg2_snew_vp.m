% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 19:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPRR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:21:42
% EndTime: 2019-05-05 19:21:59
% DurationCPUTime: 7.26s
% Computational Cost: add. (31901->458), mult. (70682->645), div. (0->0), fcn. (49162->10), ass. (0->286)
t237 = sin(pkin(10));
t238 = cos(pkin(10));
t242 = sin(qJ(3));
t246 = cos(qJ(3));
t295 = t246 * t237;
t216 = (t242 * t238 + t295) * qJD(1);
t290 = qJD(1) * t246;
t218 = -t237 * t242 * qJD(1) + t238 * t290;
t310 = t218 * t216;
t334 = qJDD(3) - t310;
t336 = t237 * t334;
t335 = t238 * t334;
t281 = qJD(1) * qJD(3);
t270 = t246 * t281;
t279 = t242 * qJDD(1);
t222 = -t270 - t279;
t230 = t246 * qJDD(1);
t271 = t242 * t281;
t223 = t230 - t271;
t190 = t237 * t222 + t238 * t223;
t289 = qJD(3) * t216;
t175 = t190 - t289;
t240 = sin(qJ(6));
t241 = sin(qJ(5));
t245 = cos(qJ(5));
t196 = -t245 * qJD(3) + t241 * t218;
t198 = t241 * qJD(3) + t245 * t218;
t244 = cos(qJ(6));
t167 = t244 * t196 + t240 * t198;
t169 = -t240 * t196 + t244 * t198;
t131 = t169 * t167;
t266 = -t238 * t222 + t237 * t223;
t187 = qJDD(5) + t266;
t186 = qJDD(6) + t187;
t329 = -t131 + t186;
t333 = t240 * t329;
t172 = t198 * t196;
t327 = -t172 + t187;
t332 = t241 * t327;
t331 = t244 * t329;
t330 = t245 * t327;
t288 = qJD(3) * t218;
t173 = t266 + t288;
t261 = -t241 * qJDD(3) - t245 * t190;
t155 = -t196 * qJD(5) - t261;
t267 = -t245 * qJDD(3) + t241 * t190;
t256 = t198 * qJD(5) + t267;
t106 = -t167 * qJD(6) + t244 * t155 - t240 * t256;
t212 = qJD(5) + t216;
t207 = qJD(6) + t212;
t150 = t207 * t167;
t328 = -t150 + t106;
t180 = t212 * t196;
t136 = t155 + t180;
t243 = sin(qJ(1));
t247 = cos(qJ(1));
t263 = t243 * g(1) - t247 * g(2);
t258 = qJDD(2) - t263;
t249 = qJD(1) ^ 2;
t293 = t249 * qJ(2);
t254 = t258 - t293;
t323 = pkin(7) + pkin(1);
t252 = -qJDD(1) * t323 + t254;
t192 = t246 * g(3) - t242 * t252;
t257 = qJD(3) * pkin(3) - qJ(4) * t290;
t234 = t242 ^ 2;
t309 = t234 * t249;
t163 = -pkin(3) * t309 + t222 * qJ(4) - qJD(3) * t257 - t192;
t251 = t246 * t252;
t294 = t246 * t249;
t250 = t251 - t223 * qJ(4) + qJDD(3) * pkin(3) + (-pkin(3) * t294 - qJ(4) * t281 + g(3)) * t242;
t123 = -0.2e1 * qJD(4) * t216 + t238 * t163 + t237 * t250;
t269 = t240 * t155 + t244 * t256;
t85 = (qJD(6) - t207) * t169 + t269;
t132 = (qJD(5) - t212) * t198 + t267;
t233 = qJDD(1) * qJ(2);
t264 = t247 * g(1) + t243 * g(2);
t259 = -t233 + t264;
t326 = -t222 * pkin(3) - (qJ(4) * t234 + t323) * t249 + t257 * t290 + qJDD(4) - t259;
t164 = t167 ^ 2;
t165 = t169 ^ 2;
t325 = t196 ^ 2;
t195 = t198 ^ 2;
t206 = t207 ^ 2;
t211 = t212 ^ 2;
t214 = t216 ^ 2;
t215 = t218 ^ 2;
t324 = 0.2e1 * qJD(4);
t181 = t216 * pkin(4) - t218 * pkin(8);
t248 = qJD(3) ^ 2;
t104 = -t248 * pkin(4) + qJDD(3) * pkin(8) - t216 * t181 + t123;
t280 = qJD(2) * qJD(1);
t232 = 0.2e1 * t280;
t111 = t173 * pkin(4) - t175 * pkin(8) + t232 + t326;
t66 = t241 * t104 - t245 * t111;
t51 = pkin(5) * t327 - pkin(9) * t136 - t66;
t177 = t212 * pkin(5) - t198 * pkin(9);
t67 = t245 * t104 + t241 * t111;
t59 = -pkin(5) * t325 - pkin(9) * t256 - t212 * t177 + t67;
t29 = t240 * t59 - t244 * t51;
t30 = t240 * t51 + t244 * t59;
t13 = t240 * t30 - t244 * t29;
t322 = pkin(5) * t13;
t88 = t150 + t106;
t54 = -t240 * t85 - t244 * t88;
t321 = pkin(5) * t54;
t268 = t237 * t163 - t238 * t250;
t103 = -qJDD(3) * pkin(4) - t248 * pkin(8) + (t324 + t181) * t218 + t268;
t68 = pkin(5) * t256 - pkin(9) * t325 + t198 * t177 + t103;
t320 = t240 * t68;
t319 = t241 * t13;
t318 = t244 * t68;
t317 = t245 * t13;
t122 = t218 * t324 + t268;
t74 = -t238 * t122 + t237 * t123;
t316 = t246 * t74;
t315 = qJDD(1) * pkin(1);
t314 = t207 * t240;
t313 = t207 * t244;
t312 = t212 * t241;
t311 = t212 * t245;
t235 = t246 ^ 2;
t308 = t235 * t249;
t278 = -0.2e1 * t280;
t166 = t278 - t326;
t307 = t237 * t166;
t184 = qJDD(3) + t310;
t306 = t237 * t184;
t305 = t238 * t166;
t304 = t238 * t184;
t120 = t131 + t186;
t303 = t240 * t120;
t302 = t241 * t103;
t143 = t172 + t187;
t301 = t241 * t143;
t273 = t242 * t294;
t300 = t242 * (qJDD(3) + t273);
t299 = t244 * t120;
t298 = t245 * t103;
t297 = t245 * t143;
t296 = t246 * (qJDD(3) - t273);
t291 = t234 + t235;
t287 = qJD(3) * t237;
t286 = qJD(3) * t238;
t283 = qJD(5) + t212;
t277 = t237 * t131;
t276 = t237 * t172;
t275 = t238 * t131;
t274 = t238 * t172;
t272 = -pkin(4) * t238 - pkin(3);
t14 = t240 * t29 + t244 * t30;
t41 = t241 * t66 + t245 * t67;
t75 = t237 * t122 + t238 * t123;
t40 = t241 * t67 - t245 * t66;
t32 = -t238 * t103 + t237 * t41;
t18 = t242 * (t237 * t103 + t238 * t41) + t246 * t32;
t191 = t242 * g(3) + t251;
t158 = t246 * t191 - t242 * t192;
t126 = -t206 - t164;
t76 = t240 * t126 + t331;
t260 = pkin(5) * t76 - t29;
t174 = -t266 + t288;
t140 = -t165 - t206;
t92 = t244 * t140 - t303;
t255 = pkin(5) * t92 - t30;
t225 = t291 * qJDD(1);
t224 = t230 - 0.2e1 * t271;
t221 = 0.2e1 * t270 + t279;
t213 = -t254 + t315;
t205 = -t215 - t248;
t204 = -t215 + t248;
t203 = t214 - t248;
t202 = t249 * t323 + t259 + t278;
t200 = -t300 + t246 * (-t248 - t308);
t199 = t242 * (-t248 - t309) + t296;
t182 = -t248 - t214;
t179 = -t195 + t211;
t178 = -t211 + t325;
t176 = t190 + t289;
t171 = -t214 - t215;
t170 = t195 - t325;
t162 = -t195 - t211;
t161 = -t237 * t205 - t304;
t160 = t238 * t205 - t306;
t153 = -t211 - t325;
t149 = t195 + t325;
t148 = -t165 + t206;
t147 = t164 - t206;
t146 = t238 * t182 - t336;
t145 = t237 * t182 + t335;
t141 = (-t196 * t245 + t198 * t241) * t212;
t139 = t238 * t174 + t237 * t176;
t138 = t237 * t174 - t238 * t176;
t137 = t196 * t283 + t261;
t135 = t155 - t180;
t133 = -t198 * t283 - t267;
t130 = t165 - t164;
t129 = t245 * t155 - t198 * t312;
t128 = t196 * t311 + t241 * t256;
t127 = t246 * t160 + t242 * t161;
t125 = t245 * t178 - t301;
t124 = -t241 * t179 + t330;
t117 = -t241 * t162 - t297;
t116 = t245 * t162 - t301;
t115 = (-t167 * t244 + t169 * t240) * t207;
t114 = (-t167 * t240 - t169 * t244) * t207;
t113 = t245 * t153 - t332;
t112 = t241 * t153 + t330;
t108 = t246 * t145 + t242 * t146;
t107 = -t164 - t165;
t105 = -t169 * qJD(6) - t269;
t101 = t246 * t138 + t242 * t139;
t100 = -t132 * t245 + t241 * t136;
t99 = t245 * t133 - t241 * t135;
t98 = -t132 * t241 - t245 * t136;
t97 = t244 * t147 - t303;
t96 = -t240 * t148 + t331;
t95 = t240 * t147 + t299;
t94 = t244 * t148 + t333;
t93 = -t240 * t140 - t299;
t91 = t238 * t117 - t237 * t137;
t90 = t237 * t117 + t238 * t137;
t84 = (qJD(6) + t207) * t169 + t269;
t83 = t238 * t113 - t237 * t133;
t82 = t237 * t113 + t238 * t133;
t81 = t244 * t106 - t169 * t314;
t80 = t240 * t106 + t169 * t313;
t79 = -t240 * t105 + t167 * t313;
t78 = t244 * t105 + t167 * t314;
t77 = t244 * t126 - t333;
t73 = t238 * t100 - t237 * t149;
t72 = t237 * t100 + t238 * t149;
t71 = -t241 * t114 + t245 * t115;
t70 = -pkin(8) * t116 + t298;
t69 = -pkin(8) * t112 + t302;
t64 = -t241 * t95 + t245 * t97;
t63 = -t241 * t94 + t245 * t96;
t62 = -t241 * t92 + t245 * t93;
t61 = t241 * t93 + t245 * t92;
t60 = -pkin(4) * t116 + t67;
t58 = -pkin(4) * t112 + t66;
t57 = t242 * t91 + t246 * t90;
t56 = t240 * t88 - t244 * t85;
t55 = -t240 * t328 - t244 * t84;
t53 = -t240 * t84 + t244 * t328;
t52 = t242 * t83 + t246 * t82;
t50 = -t241 * t80 + t245 * t81;
t49 = -t241 * t78 + t245 * t79;
t47 = -t241 * t76 + t245 * t77;
t46 = t241 * t77 + t245 * t76;
t45 = t242 * t75 + t316;
t44 = -pkin(9) * t92 + t318;
t43 = t242 * t73 + t246 * t72;
t42 = -pkin(9) * t76 + t320;
t39 = t237 * t328 + t238 * t62;
t38 = t237 * t62 - t238 * t328;
t37 = t237 * t84 + t238 * t47;
t36 = t237 * t47 - t238 * t84;
t35 = -pkin(5) * t328 + pkin(9) * t93 + t320;
t34 = -pkin(5) * t84 + pkin(9) * t77 - t318;
t31 = -pkin(8) * t98 - t40;
t27 = -t241 * t54 + t245 * t56;
t26 = -t241 * t53 + t245 * t55;
t25 = t241 * t56 + t245 * t54;
t24 = t237 * t107 + t238 * t27;
t23 = -t238 * t107 + t237 * t27;
t22 = t242 * t39 + t246 * t38;
t21 = -pkin(4) * t25 - t321;
t20 = t242 * t37 + t246 * t36;
t19 = -pkin(4) * t61 - t255;
t17 = -pkin(4) * t46 - t260;
t16 = -pkin(8) * t61 - t241 * t35 + t245 * t44;
t15 = -pkin(8) * t46 - t241 * t34 + t245 * t42;
t12 = -pkin(5) * t68 + pkin(9) * t14;
t11 = t246 * t23 + t242 * t24;
t10 = -pkin(9) * t54 - t13;
t9 = -pkin(5) * t107 + pkin(9) * t56 + t14;
t8 = t245 * t14 - t319;
t7 = t241 * t14 + t317;
t6 = t237 * t68 + t238 * t8;
t5 = t237 * t8 - t238 * t68;
t4 = -pkin(4) * t7 - t322;
t3 = -pkin(8) * t25 + t245 * t10 - t241 * t9;
t2 = -pkin(8) * t7 - pkin(9) * t317 - t241 * t12;
t1 = t242 * t6 + t246 * t5;
t28 = [0, 0, 0, 0, 0, qJDD(1), t263, t264, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t258 - 0.2e1 * t315, t232 + 0.2e1 * t233 - t264, pkin(1) * t213 + qJ(2) * (-t249 * pkin(1) + t232 - t259), (t223 - t271) * t246, -t246 * t221 - t242 * t224, t296 - t242 * (t248 - t308), (-t222 + t270) * t242, t246 * (-t248 + t309) - t300, 0, qJ(2) * t221 - t199 * t323 - t242 * t202, qJ(2) * t224 - t200 * t323 - t246 * t202, t225 * t323 - t291 * t293 - t158, -qJ(2) * t202 - t158 * t323, t246 * (t238 * t190 - t218 * t287) - t242 * (t237 * t190 + t218 * t286), t246 * (-t238 * t173 - t237 * t175) - t242 * (-t237 * t173 + t238 * t175), t246 * (-t237 * t204 + t335) - t242 * (t238 * t204 + t336), t246 * (t216 * t286 + t237 * t266) - t242 * (t216 * t287 - t238 * t266), t246 * (t238 * t203 - t306) - t242 * (t237 * t203 + t304), (t246 * (-t216 * t238 + t218 * t237) - t242 * (-t216 * t237 - t218 * t238)) * qJD(3), t246 * (-qJ(4) * t145 - t307) - t242 * (-pkin(3) * t173 + qJ(4) * t146 + t305) + qJ(2) * t173 - t323 * t108, t246 * (-qJ(4) * t160 - t305) - t242 * (-pkin(3) * t175 + qJ(4) * t161 - t307) + qJ(2) * t175 - t323 * t127, t246 * (-qJ(4) * t138 - t74) - t242 * (-pkin(3) * t171 + qJ(4) * t139 + t75) + qJ(2) * t171 - t323 * t101, -qJ(4) * t316 - t242 * (pkin(3) * t166 + qJ(4) * t75) - qJ(2) * t166 - t323 * t45, t246 * (t238 * t129 + t276) - t242 * (t237 * t129 - t274), t246 * (t237 * t170 + t238 * t99) - t242 * (-t238 * t170 + t237 * t99), t246 * (t238 * t124 + t237 * t136) - t242 * (t237 * t124 - t238 * t136), t246 * (t238 * t128 - t276) - t242 * (t237 * t128 + t274), t246 * (t238 * t125 - t237 * t132) - t242 * (t237 * t125 + t238 * t132), t246 * (t238 * t141 + t237 * t187) - t242 * (t237 * t141 - t238 * t187), t246 * (-qJ(4) * t82 - t237 * t58 + t238 * t69) - t242 * (-pkin(3) * t112 + qJ(4) * t83 + t237 * t69 + t238 * t58) + qJ(2) * t112 - t323 * t52, t246 * (-qJ(4) * t90 - t237 * t60 + t238 * t70) - t242 * (-pkin(3) * t116 + qJ(4) * t91 + t237 * t70 + t238 * t60) + qJ(2) * t116 - t323 * t57, t246 * (-qJ(4) * t72 + t238 * t31) - t242 * (qJ(4) * t73 + t237 * t31) + (pkin(4) * t295 - t242 * t272 + qJ(2)) * t98 - t323 * t43, (t246 * (pkin(4) * t237 - pkin(8) * t238) - t242 * (-pkin(8) * t237 + t272) + qJ(2)) * t40 + (-t323 - qJ(4)) * t18, t246 * (t238 * t50 + t277) - t242 * (t237 * t50 - t275), t246 * (t237 * t130 + t238 * t26) - t242 * (-t238 * t130 + t237 * t26), t246 * (t237 * t88 + t238 * t63) - t242 * (t237 * t63 - t238 * t88), t246 * (t238 * t49 - t277) - t242 * (t237 * t49 + t275), t246 * (-t237 * t85 + t238 * t64) - t242 * (t237 * t64 + t238 * t85), t246 * (t237 * t186 + t238 * t71) - t242 * (-t238 * t186 + t237 * t71), t246 * (-qJ(4) * t36 + t238 * t15 - t237 * t17) - t242 * (-pkin(3) * t46 + qJ(4) * t37 + t237 * t15 + t238 * t17) + qJ(2) * t46 - t323 * t20, t246 * (-qJ(4) * t38 + t238 * t16 - t237 * t19) - t242 * (-pkin(3) * t61 + qJ(4) * t39 + t237 * t16 + t238 * t19) + qJ(2) * t61 - t323 * t22, t246 * (-qJ(4) * t23 - t237 * t21 + t238 * t3) - t242 * (-pkin(3) * t25 + qJ(4) * t24 + t238 * t21 + t237 * t3) + qJ(2) * t25 - t323 * t11, t246 * (-qJ(4) * t5 + t238 * t2 - t237 * t4) - t242 * (-pkin(3) * t7 + qJ(4) * t6 + t237 * t2 + t238 * t4) + qJ(2) * t7 - t323 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t249, -t213, 0, 0, 0, 0, 0, 0, t199, t200, -t225, t158, 0, 0, 0, 0, 0, 0, t108, t127, t101, t45, 0, 0, 0, 0, 0, 0, t52, t57, t43, t18, 0, 0, 0, 0, 0, 0, t20, t22, t11, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t273, (-t234 + t235) * t249, t230, -t273, -t279, qJDD(3), t191, t192, 0, 0, t310, t215 - t214, t176, -t310, t174, qJDD(3), pkin(3) * t145 - t122, pkin(3) * t160 - t123, pkin(3) * t138, pkin(3) * t74, t241 * t155 + t198 * t311, t241 * t133 + t245 * t135, t245 * t179 + t332, t196 * t312 - t245 * t256, t241 * t178 + t297, (-t196 * t241 - t198 * t245) * t212, pkin(3) * t82 + pkin(4) * t133 + pkin(8) * t113 - t298, pkin(3) * t90 + pkin(4) * t137 + pkin(8) * t117 + t302, pkin(3) * t72 + pkin(4) * t149 + pkin(8) * t100 + t41, pkin(3) * t32 - pkin(4) * t103 + pkin(8) * t41, t241 * t81 + t245 * t80, t241 * t55 + t245 * t53, t241 * t96 + t245 * t94, t241 * t79 + t245 * t78, t241 * t97 + t245 * t95, t245 * t114 + t241 * t115, pkin(3) * t36 - pkin(4) * t84 + pkin(8) * t47 + t241 * t42 + t245 * t34, pkin(3) * t38 - pkin(4) * t328 + pkin(8) * t62 + t241 * t44 + t245 * t35, pkin(3) * t23 - pkin(4) * t107 + pkin(8) * t27 + t241 * t10 + t245 * t9, pkin(3) * t5 - pkin(4) * t68 + pkin(8) * t8 - pkin(9) * t319 + t245 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t175, t171, -t166, 0, 0, 0, 0, 0, 0, t112, t116, t98, t40, 0, 0, 0, 0, 0, 0, t46, t61, t25, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, t170, t136, -t172, -t132, t187, -t66, -t67, 0, 0, t131, t130, t88, -t131, -t85, t186, t260, t255, t321, t322; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, t130, t88, -t131, -t85, t186, -t29, -t30, 0, 0;];
tauJ_reg  = t28;
