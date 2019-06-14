% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 15:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRRR2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:21:32
% EndTime: 2019-05-05 15:21:45
% DurationCPUTime: 6.45s
% Computational Cost: add. (31616->465), mult. (71431->685), div. (0->0), fcn. (52633->12), ass. (0->271)
t236 = sin(pkin(11));
t231 = t236 ^ 2;
t238 = cos(pkin(11));
t232 = t238 ^ 2;
t313 = qJD(1) ^ 2;
t219 = (t231 + t232) * t313;
t241 = sin(qJ(6));
t243 = sin(qJ(4));
t246 = cos(qJ(4));
t260 = t236 * t246 + t238 * t243;
t213 = t260 * qJD(1);
t242 = sin(qJ(5));
t245 = cos(qJ(5));
t196 = -qJD(4) * t245 + t213 * t242;
t198 = qJD(4) * t242 + t213 * t245;
t244 = cos(qJ(6));
t166 = t196 * t244 + t198 * t241;
t168 = -t196 * t241 + t198 * t244;
t132 = t168 * t166;
t275 = t238 * qJDD(1);
t276 = t236 * qJDD(1);
t263 = t243 * t276 - t246 * t275;
t280 = t213 * qJD(4);
t189 = -t263 - t280;
t182 = qJDD(5) - t189;
t181 = qJDD(6) + t182;
t320 = -t132 + t181;
t329 = t241 * t320;
t174 = t198 * t196;
t318 = -t174 + t182;
t328 = t242 * t318;
t283 = qJD(1) * t238;
t288 = t236 * t243;
t211 = qJD(1) * t288 - t246 * t283;
t192 = t213 * t211;
t316 = qJDD(4) - t192;
t327 = t243 * t316;
t326 = t244 * t320;
t325 = t245 * t318;
t324 = t246 * t316;
t237 = sin(pkin(10));
t307 = sin(qJ(1));
t308 = cos(qJ(1));
t255 = g(1) * t307 - g(2) * t308;
t253 = qJDD(1) * pkin(1) + t255;
t256 = g(1) * t308 + g(2) * t307;
t217 = -pkin(1) * t313 - t256;
t239 = cos(pkin(10));
t287 = t239 * t217;
t251 = qJDD(1) * qJ(3) + t237 * t253 + t287;
t284 = -g(3) + qJDD(2);
t311 = 2 * qJD(3);
t170 = t238 * (-pkin(2) * t313 + t251) + t236 * t284 + t283 * t311;
t282 = t313 * t232;
t158 = -pkin(3) * t282 + pkin(7) * t275 + t170;
t252 = -t237 * t255 - t287;
t267 = t238 * t284;
t268 = pkin(1) * t237 + qJ(3);
t322 = t268 + pkin(7);
t249 = t267 + (-t322 * qJDD(1) + (-(2 * qJD(3)) + (t238 * pkin(3) + pkin(2)) * qJD(1)) * qJD(1) + t252) * t236;
t123 = t246 * t158 + t243 * t249;
t214 = t239 * t253;
t264 = -t237 * t217 + t214;
t180 = -qJDD(1) * pkin(2) - qJ(3) * t313 + qJDD(3) - t264;
t270 = pkin(1) * t239 + pkin(2);
t323 = -qJDD(1) * t270 + t219 * t268 + t180;
t210 = t260 * qJDD(1);
t281 = t211 * qJD(4);
t191 = t210 - t281;
t259 = -qJDD(4) * t242 - t245 * t191;
t156 = -qJD(5) * t196 - t259;
t265 = -qJDD(4) * t245 + t191 * t242;
t257 = qJD(5) * t198 + t265;
t103 = -qJD(6) * t166 + t156 * t244 - t241 * t257;
t207 = qJD(5) + t211;
t203 = qJD(6) + t207;
t151 = t203 * t166;
t319 = -t151 + t103;
t178 = t207 * t196;
t137 = t156 + t178;
t266 = t156 * t241 + t244 * t257;
t83 = (qJD(6) - t203) * t168 + t266;
t133 = (qJD(5) - t207) * t198 + t265;
t164 = t166 ^ 2;
t165 = t168 ^ 2;
t314 = t196 ^ 2;
t195 = t198 ^ 2;
t202 = t203 ^ 2;
t206 = t207 ^ 2;
t208 = t211 ^ 2;
t209 = t213 ^ 2;
t312 = qJD(4) ^ 2;
t183 = pkin(4) * t211 - pkin(8) * t213;
t105 = -pkin(4) * t312 + qJDD(4) * pkin(8) - t211 * t183 + t123;
t163 = -pkin(3) * t275 + t180 + (-t231 * t313 - t282) * pkin(7);
t109 = (-t191 + t281) * pkin(8) + (-t189 + t280) * pkin(4) + t163;
t66 = t105 * t242 - t109 * t245;
t51 = pkin(5) * t318 - pkin(9) * t137 - t66;
t175 = pkin(5) * t207 - pkin(9) * t198;
t67 = t105 * t245 + t109 * t242;
t59 = -pkin(5) * t314 - pkin(9) * t257 - t175 * t207 + t67;
t29 = t241 * t59 - t244 * t51;
t30 = t241 * t51 + t244 * t59;
t13 = t241 * t30 - t244 * t29;
t310 = pkin(5) * t13;
t86 = t151 + t103;
t54 = -t241 * t83 - t244 * t86;
t309 = pkin(5) * t54;
t306 = t13 * t242;
t305 = t13 * t245;
t122 = t158 * t243 - t246 * t249;
t74 = -t122 * t246 + t123 * t243;
t304 = t236 * t74;
t104 = -qJDD(4) * pkin(4) - pkin(8) * t312 + t183 * t213 + t122;
t68 = pkin(5) * t257 - pkin(9) * t314 + t175 * t198 + t104;
t303 = t241 * t68;
t302 = t244 * t68;
t301 = t104 * t242;
t300 = t104 * t245;
t119 = t132 + t181;
t299 = t119 * t241;
t298 = t119 * t244;
t142 = t174 + t182;
t297 = t142 * t242;
t296 = t142 * t245;
t295 = t163 * t243;
t294 = t163 * t246;
t186 = qJDD(4) + t192;
t293 = t186 * t246;
t292 = t203 * t241;
t291 = t203 * t244;
t290 = t207 * t242;
t289 = t207 * t245;
t286 = t243 * t186;
t278 = qJD(5) + t207;
t274 = t246 * t132;
t273 = t246 * t174;
t272 = t243 * t132;
t271 = t243 * t174;
t269 = -pkin(4) * t246 - pkin(3);
t14 = t241 * t29 + t244 * t30;
t41 = t242 * t66 + t245 * t67;
t75 = t122 * t243 + t123 * t246;
t169 = -t267 + ((-pkin(2) * qJD(1) + t311) * qJD(1) + t251) * t236;
t128 = t236 * t169 + t170 * t238;
t40 = t242 * t67 - t245 * t66;
t126 = -t202 - t164;
t76 = t126 * t241 + t326;
t258 = pkin(5) * t76 - t29;
t139 = -t165 - t202;
t92 = t139 * t244 - t299;
t254 = pkin(5) * t92 - t30;
t228 = t232 * qJDD(1);
t227 = t231 * qJDD(1);
t218 = t228 + t227;
t201 = -t209 - t312;
t200 = -t209 + t312;
t199 = t208 - t312;
t190 = t210 - 0.2e1 * t281;
t188 = t263 + 0.2e1 * t280;
t184 = -t312 - t208;
t177 = -t195 + t206;
t176 = -t206 + t314;
t173 = -t208 - t209;
t172 = t195 - t314;
t161 = -t195 - t206;
t160 = -t201 * t243 - t293;
t159 = t201 * t246 - t286;
t157 = -t206 - t314;
t150 = t195 + t314;
t149 = t210 * t243 - t246 * t263;
t148 = -t210 * t246 - t243 * t263;
t147 = -t165 + t202;
t146 = t164 - t202;
t145 = t184 * t246 - t327;
t144 = t184 * t243 + t324;
t140 = (-t196 * t245 + t198 * t242) * t207;
t138 = t196 * t278 + t259;
t136 = t156 - t178;
t134 = -t198 * t278 - t265;
t131 = t165 - t164;
t130 = t156 * t245 - t198 * t290;
t129 = t196 * t289 + t242 * t257;
t127 = -t159 * t236 + t160 * t238;
t125 = t176 * t245 - t297;
t124 = -t177 * t242 + t325;
t117 = -t161 * t242 - t296;
t116 = t161 * t245 - t297;
t115 = (-t166 * t244 + t168 * t241) * t203;
t114 = (-t166 * t241 - t168 * t244) * t203;
t113 = -t148 * t236 + t149 * t238;
t112 = t157 * t245 - t328;
t111 = t157 * t242 + t325;
t110 = -t144 * t236 + t145 * t238;
t106 = -t164 - t165;
t102 = -qJD(6) * t168 - t266;
t100 = -t133 * t245 + t137 * t242;
t99 = t134 * t245 - t136 * t242;
t98 = -t133 * t242 - t137 * t245;
t97 = t146 * t244 - t299;
t96 = -t147 * t241 + t326;
t95 = t146 * t241 + t298;
t94 = t147 * t244 + t329;
t93 = -t139 * t241 - t298;
t91 = t117 * t246 - t138 * t243;
t90 = t117 * t243 + t138 * t246;
t89 = t112 * t246 - t134 * t243;
t88 = t112 * t243 + t134 * t246;
t82 = (qJD(6) + t203) * t168 + t266;
t81 = t103 * t244 - t168 * t292;
t80 = t103 * t241 + t168 * t291;
t79 = -t102 * t241 + t166 * t291;
t78 = t102 * t244 + t166 * t292;
t77 = t126 * t244 - t329;
t73 = t100 * t246 - t150 * t243;
t72 = t100 * t243 + t150 * t246;
t71 = -t114 * t242 + t115 * t245;
t70 = -pkin(8) * t116 + t300;
t69 = -pkin(8) * t111 + t301;
t64 = -t242 * t95 + t245 * t97;
t63 = -t242 * t94 + t245 * t96;
t62 = -t242 * t92 + t245 * t93;
t61 = t242 * t93 + t245 * t92;
t60 = -pkin(4) * t116 + t67;
t58 = -pkin(4) * t111 + t66;
t57 = -t236 * t90 + t238 * t91;
t56 = t241 * t86 - t244 * t83;
t55 = -t241 * t319 - t244 * t82;
t53 = -t241 * t82 + t244 * t319;
t52 = -t236 * t88 + t238 * t89;
t50 = -t242 * t80 + t245 * t81;
t49 = -t242 * t78 + t245 * t79;
t47 = -t242 * t76 + t245 * t77;
t46 = t242 * t77 + t245 * t76;
t45 = t238 * t75 - t304;
t44 = -pkin(9) * t92 + t302;
t42 = -pkin(9) * t76 + t303;
t39 = t243 * t319 + t246 * t62;
t38 = t243 * t62 - t246 * t319;
t37 = t243 * t82 + t246 * t47;
t36 = t243 * t47 - t246 * t82;
t35 = -pkin(5) * t319 + pkin(9) * t93 + t303;
t34 = -pkin(5) * t82 + pkin(9) * t77 - t302;
t33 = t104 * t243 + t246 * t41;
t32 = -t104 * t246 + t243 * t41;
t31 = -pkin(8) * t98 - t40;
t27 = -t242 * t54 + t245 * t56;
t26 = -t242 * t53 + t245 * t55;
t25 = t242 * t56 + t245 * t54;
t24 = t106 * t243 + t246 * t27;
t23 = -t106 * t246 + t243 * t27;
t22 = -t236 * t38 + t238 * t39;
t21 = -pkin(4) * t25 - t309;
t20 = -t236 * t36 + t238 * t37;
t19 = -pkin(4) * t61 - t254;
t17 = -pkin(4) * t46 - t258;
t16 = -pkin(8) * t61 - t242 * t35 + t245 * t44;
t15 = -pkin(8) * t46 - t242 * t34 + t245 * t42;
t12 = -pkin(5) * t68 + pkin(9) * t14;
t11 = -t23 * t236 + t238 * t24;
t10 = -pkin(9) * t54 - t13;
t9 = -pkin(5) * t106 + pkin(9) * t56 + t14;
t8 = t14 * t245 - t306;
t7 = t14 * t242 + t305;
t6 = t243 * t68 + t246 * t8;
t5 = t243 * t8 - t246 * t68;
t4 = -pkin(4) * t7 - t310;
t3 = -pkin(8) * t25 + t10 * t245 - t242 * t9;
t2 = -pkin(8) * t7 - pkin(9) * t305 - t12 * t242;
t1 = -t236 * t5 + t238 * t6;
t18 = [0, 0, 0, 0, 0, qJDD(1), t255, t256, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t239 - t237 * t313) + t264, (-0.2e1 * qJDD(1) * t237 - t239 * t313) * pkin(1) + t252, 0, pkin(1) * (t237 ^ 2 * t253 + t239 * t214), t227, 0.2e1 * t236 * t275, 0, t228, 0, 0, -t323 * t238, t323 * t236, pkin(2) * t219 + qJ(3) * t218 + pkin(1) * (t218 * t237 + t219 * t239) + t128, -pkin(2) * t180 + qJ(3) * t128 + pkin(1) * (t128 * t237 - t180 * t239), t236 * (t191 * t246 - t243 * t280) + t238 * (t191 * t243 + t246 * t280), t236 * (-t188 * t246 - t190 * t243) + t238 * (-t188 * t243 + t190 * t246), t236 * (-t200 * t243 + t324) + t238 * (t200 * t246 + t327), t236 * (-t189 * t243 + t246 * t281) + t238 * (t189 * t246 + t243 * t281), t236 * (t199 * t246 - t286) + t238 * (t199 * t243 + t293), (t236 * (-t211 * t246 + t213 * t243) + t238 * (-t211 * t243 - t213 * t246)) * qJD(4), t236 * (-pkin(7) * t144 + t295) + t238 * (-pkin(3) * t188 + pkin(7) * t145 - t294) - pkin(2) * t188 + qJ(3) * t110 + pkin(1) * (t110 * t237 - t188 * t239), t236 * (-pkin(7) * t159 + t294) + t238 * (-pkin(3) * t190 + pkin(7) * t160 + t295) - pkin(2) * t190 + qJ(3) * t127 + pkin(1) * (t127 * t237 - t190 * t239), t236 * (-pkin(7) * t148 - t74) + t238 * (-pkin(3) * t173 + pkin(7) * t149 + t75) - pkin(2) * t173 + qJ(3) * t113 + pkin(1) * (t113 * t237 - t173 * t239), -pkin(7) * t304 + t238 * (-pkin(3) * t163 + pkin(7) * t75) - pkin(2) * t163 + qJ(3) * t45 + pkin(1) * (-t163 * t239 + t237 * t45), t236 * (t130 * t246 + t271) + t238 * (t130 * t243 - t273), t236 * (t172 * t243 + t246 * t99) + t238 * (-t172 * t246 + t243 * t99), t236 * (t124 * t246 + t137 * t243) + t238 * (t124 * t243 - t137 * t246), t236 * (t129 * t246 - t271) + t238 * (t129 * t243 + t273), t236 * (t125 * t246 - t133 * t243) + t238 * (t125 * t243 + t133 * t246), t236 * (t140 * t246 + t182 * t243) + t238 * (t140 * t243 - t182 * t246), t236 * (-pkin(7) * t88 - t243 * t58 + t246 * t69) + t238 * (-pkin(3) * t111 + pkin(7) * t89 + t243 * t69 + t246 * t58) - pkin(2) * t111 + qJ(3) * t52 + pkin(1) * (-t111 * t239 + t237 * t52), t236 * (-pkin(7) * t90 - t243 * t60 + t246 * t70) + t238 * (-pkin(3) * t116 + pkin(7) * t91 + t243 * t70 + t246 * t60) - pkin(2) * t116 + qJ(3) * t57 + pkin(1) * (-t116 * t239 + t237 * t57), t236 * (-pkin(7) * t72 + t246 * t31) + t238 * (pkin(7) * t73 + t243 * t31) + t268 * (-t236 * t72 + t238 * t73) + (pkin(4) * t288 + t238 * t269 - t270) * t98, (t236 * (pkin(4) * t243 - pkin(8) * t246) + t238 * (-pkin(8) * t243 + t269) - t270) * t40 + t322 * (-t236 * t32 + t238 * t33), t236 * (t246 * t50 + t272) + t238 * (t243 * t50 - t274), t236 * (t131 * t243 + t246 * t26) + t238 * (-t131 * t246 + t243 * t26), t236 * (t243 * t86 + t246 * t63) + t238 * (t243 * t63 - t246 * t86), t236 * (t246 * t49 - t272) + t238 * (t243 * t49 + t274), t236 * (-t243 * t83 + t246 * t64) + t238 * (t243 * t64 + t246 * t83), t236 * (t181 * t243 + t246 * t71) + t238 * (-t181 * t246 + t243 * t71), t236 * (-pkin(7) * t36 + t15 * t246 - t17 * t243) + t238 * (-pkin(3) * t46 + pkin(7) * t37 + t15 * t243 + t17 * t246) - pkin(2) * t46 + qJ(3) * t20 + pkin(1) * (t20 * t237 - t239 * t46), t236 * (-pkin(7) * t38 + t16 * t246 - t19 * t243) + t238 * (-pkin(3) * t61 + pkin(7) * t39 + t16 * t243 + t19 * t246) - pkin(2) * t61 + qJ(3) * t22 + pkin(1) * (t22 * t237 - t239 * t61), t236 * (-pkin(7) * t23 - t21 * t243 + t246 * t3) + t238 * (-pkin(3) * t25 + pkin(7) * t24 + t21 * t246 + t243 * t3) - pkin(2) * t25 + qJ(3) * t11 + pkin(1) * (t11 * t237 - t239 * t25), t236 * (-pkin(7) * t5 + t2 * t246 - t243 * t4) + t238 * (-pkin(3) * t7 + pkin(7) * t6 + t2 * t243 + t246 * t4) - pkin(2) * t7 + qJ(3) * t1 + pkin(1) * (t1 * t237 - t239 * t7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t284, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169 * t238 + t170 * t236, 0, 0, 0, 0, 0, 0, t144 * t238 + t145 * t236, t159 * t238 + t160 * t236, t148 * t238 + t149 * t236, t236 * t75 + t238 * t74, 0, 0, 0, 0, 0, 0, t236 * t89 + t238 * t88, t236 * t91 + t238 * t90, t236 * t73 + t238 * t72, t236 * t33 + t238 * t32, 0, 0, 0, 0, 0, 0, t236 * t37 + t238 * t36, t236 * t39 + t238 * t38, t23 * t238 + t236 * t24, t236 * t6 + t238 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t275, t276, -t219, t180, 0, 0, 0, 0, 0, 0, t188, t190, t173, t163, 0, 0, 0, 0, 0, 0, t111, t116, t98, t40, 0, 0, 0, 0, 0, 0, t46, t61, t25, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, t209 - t208, t210, -t192, -t263, qJDD(4), -t122, -t123, 0, 0, t156 * t242 + t198 * t289, t134 * t242 + t136 * t245, t177 * t245 + t328, t196 * t290 - t245 * t257, t176 * t242 + t296, (-t196 * t242 - t198 * t245) * t207, pkin(4) * t134 + pkin(8) * t112 - t300, pkin(4) * t138 + pkin(8) * t117 + t301, pkin(4) * t150 + pkin(8) * t100 + t41, -pkin(4) * t104 + pkin(8) * t41, t242 * t81 + t245 * t80, t242 * t55 + t245 * t53, t242 * t96 + t245 * t94, t242 * t79 + t245 * t78, t242 * t97 + t245 * t95, t114 * t245 + t115 * t242, -pkin(4) * t82 + pkin(8) * t47 + t242 * t42 + t245 * t34, -pkin(4) * t319 + pkin(8) * t62 + t242 * t44 + t245 * t35, -pkin(4) * t106 + pkin(8) * t27 + t10 * t242 + t245 * t9, -pkin(4) * t68 + pkin(8) * t8 - pkin(9) * t306 + t12 * t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, t172, t137, -t174, -t133, t182, -t66, -t67, 0, 0, t132, t131, t86, -t132, -t83, t181, t258, t254, t309, t310; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t131, t86, -t132, -t83, t181, -t29, -t30, 0, 0;];
tauJ_reg  = t18;
