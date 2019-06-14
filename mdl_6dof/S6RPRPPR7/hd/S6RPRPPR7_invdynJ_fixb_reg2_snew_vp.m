% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 17:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPPR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:16:30
% EndTime: 2019-05-05 17:16:42
% DurationCPUTime: 5.91s
% Computational Cost: add. (10548->364), mult. (23751->457), div. (0->0), fcn. (15469->8), ass. (0->224)
t272 = pkin(7) + pkin(1);
t192 = sin(qJ(3));
t195 = cos(qJ(3));
t197 = qJD(3) ^ 2;
t188 = sin(pkin(9));
t189 = cos(pkin(9));
t241 = qJD(1) * t195;
t167 = -t188 * t192 * qJD(1) + t189 * t241;
t276 = t167 ^ 2;
t147 = t276 + t197;
t244 = t195 * t188;
t165 = (-t192 * t189 - t244) * qJD(1);
t249 = t167 * t165;
t285 = qJDD(3) - t249;
t301 = t285 * t188;
t87 = t147 * t189 + t301;
t300 = t285 * t189;
t89 = -t147 * t188 + t300;
t57 = t192 * t89 + t195 * t87;
t329 = t272 * t57;
t328 = pkin(3) * t87;
t327 = qJ(4) * t87;
t326 = qJ(4) * t89;
t148 = t276 - t197;
t284 = qJDD(3) + t249;
t302 = t189 * t284;
t303 = t188 * t284;
t325 = t192 * (-t148 * t189 + t303) - t195 * (t148 * t188 + t302);
t277 = t165 ^ 2;
t111 = -t197 - t277;
t75 = t111 * t188 + t302;
t78 = -t111 * t189 + t303;
t47 = t192 * t78 - t195 * t75;
t324 = t272 * t47;
t280 = -t276 - t277;
t234 = qJD(1) * qJD(3);
t225 = t195 * t234;
t231 = t192 * qJDD(1);
t172 = -t225 - t231;
t179 = t195 * qJDD(1);
t226 = t192 * t234;
t173 = t179 - t226;
t128 = -t189 * t172 + t173 * t188;
t240 = qJD(3) * t167;
t101 = -t128 + t240;
t129 = t172 * t188 + t173 * t189;
t237 = t165 * qJD(3);
t288 = t129 - t237;
t298 = t101 * t188 - t189 * t288;
t299 = t101 * t189 + t188 * t288;
t313 = t192 * t299 + t195 * t298;
t321 = qJ(2) * t280 - t272 * t313;
t143 = t277 - t197;
t320 = t192 * (t143 * t188 + t300) + t195 * (-t143 * t189 + t301);
t319 = pkin(3) * t298;
t318 = qJ(4) * t298;
t287 = t129 + t237;
t317 = qJ(5) * t287;
t314 = -pkin(3) * t280 + qJ(4) * t299;
t312 = pkin(3) * t75;
t311 = qJ(4) * t75;
t310 = qJ(4) * t78;
t99 = t128 + t240;
t297 = t195 * (-t188 * t287 - t189 * t99) - t192 * (-t188 * t99 + t189 * t287);
t275 = 2 * qJD(4);
t274 = 2 * qJD(5);
t191 = sin(qJ(6));
t122 = qJDD(6) + t129;
t194 = cos(qJ(6));
t134 = qJD(3) * t191 + t194 * t165;
t136 = qJD(3) * t194 - t165 * t191;
t98 = t136 * t134;
t282 = t122 - t98;
t290 = t191 * t282;
t289 = t194 * t282;
t198 = qJD(1) ^ 2;
t243 = t195 * t198;
t193 = sin(qJ(1));
t196 = cos(qJ(1));
t223 = g(1) * t193 - t196 * g(2);
t215 = qJDD(2) - t223;
t261 = qJ(2) * t198;
t208 = t215 - t261;
t142 = -t272 * qJDD(1) + t208;
t252 = t142 * t195;
t92 = qJDD(3) * pkin(3) - qJ(4) * t173 + t252 + (-pkin(3) * t243 - qJ(4) * t234 + g(3)) * t192;
t131 = t195 * g(3) - t192 * t142;
t211 = qJD(3) * pkin(3) - qJ(4) * t241;
t185 = t192 ^ 2;
t248 = t185 * t198;
t93 = -pkin(3) * t248 + t172 * qJ(4) - qJD(3) * t211 - t131;
t222 = t188 * t93 - t189 * t92;
t209 = -qJDD(3) * pkin(4) - t197 * qJ(5) + qJDD(5) + t222;
t206 = -qJDD(3) * pkin(8) + t209;
t286 = pkin(5) * t288 + t206;
t281 = -pkin(4) * t240 + t167 * t274;
t54 = t165 * t275 + t188 * t92 + t189 * t93;
t279 = t276 - t277;
t183 = qJDD(1) * qJ(2);
t216 = g(1) * t196 + g(2) * t193;
t213 = -t183 + t216;
t278 = pkin(3) * t172 + (qJ(4) * t185 + t272) * t198 - t211 * t241 - qJDD(4) + t213;
t132 = t134 ^ 2;
t133 = t136 ^ 2;
t159 = qJD(6) + t167;
t157 = t159 ^ 2;
t273 = pkin(4) + pkin(8);
t271 = pkin(4) * t128;
t270 = pkin(4) * t189;
t233 = qJD(2) * qJD(1);
t230 = -0.2e1 * t233;
t94 = t230 + t278;
t268 = t188 * t94;
t267 = t189 * t94;
t141 = pkin(5) * t167 - qJD(3) * pkin(8);
t110 = -pkin(4) * t165 - qJ(5) * t167;
t210 = -t197 * pkin(4) + t110 * t165 + t54;
t232 = qJDD(3) * qJ(5);
t28 = t232 - t128 * pkin(5) - t277 * pkin(8) + (t274 + t141) * qJD(3) + t210;
t266 = t191 * t28;
t73 = t122 + t98;
t265 = t191 * t73;
t264 = t194 * t28;
t263 = t194 * t73;
t53 = t167 * t275 + t222;
t25 = t188 * t54 - t189 * t53;
t262 = t195 * t25;
t260 = qJ(5) * t189;
t259 = qJDD(1) * pkin(1);
t251 = t159 * t191;
t250 = t159 * t194;
t186 = t195 ^ 2;
t247 = t186 * t198;
t227 = t192 * t243;
t246 = t192 * (qJDD(3) + t227);
t245 = t195 * (qJDD(3) - t227);
t242 = t185 + t186;
t236 = t275 + t110;
t235 = qJD(6) + t159;
t229 = t188 * t98;
t228 = t189 * t98;
t224 = -qJ(5) * t188 - pkin(3);
t26 = t188 * t53 + t189 * t54;
t212 = (-pkin(8) * t165 + t236) * t167;
t181 = 0.2e1 * t233;
t199 = t181 - t278 - t281 - t317;
t33 = -pkin(5) * t277 + t273 * t128 - t141 * t167 + t199;
t16 = t191 * t33 - t194 * (t212 + t286);
t219 = qJDD(3) * t191 - t194 * t128;
t204 = t191 * t212 + t194 * t33;
t17 = t191 * t286 + t204;
t7 = -t16 * t194 + t17 * t191;
t8 = t16 * t191 + t17 * t194;
t205 = qJD(3) * t274 + t210;
t43 = t205 + t232;
t44 = t236 * t167 + t209;
t21 = t188 * t43 - t189 * t44;
t11 = t192 * (t188 * t44 + t189 * t43) + t195 * t21;
t130 = g(3) * t192 + t252;
t85 = t195 * t130 - t192 * t131;
t214 = qJDD(3) * t194 + t128 * t191;
t207 = (-qJD(6) + t159) * t136 - t219;
t82 = -qJD(6) * t134 + t214;
t202 = t195 * (t189 * t129 - t188 * t240) - t192 * (t188 * t129 + t189 * t240);
t201 = t195 * (t128 * t188 - t189 * t237) - t192 * (-t189 * t128 - t188 * t237);
t200 = (t195 * (t165 * t189 + t167 * t188) - t192 * (t165 * t188 - t167 * t189)) * qJD(3);
t175 = t242 * qJDD(1);
t174 = t179 - 0.2e1 * t226;
t171 = 0.2e1 * t225 + t231;
t160 = -t208 + t259;
t140 = t272 * t198 + t213 + t230;
t138 = -t246 + t195 * (-t197 - t247);
t137 = t192 * (-t197 - t248) + t245;
t109 = t159 * t134;
t108 = -t133 + t157;
t107 = t132 - t157;
t95 = t133 - t132;
t91 = -t133 - t157;
t81 = -qJD(6) * t136 - t219;
t80 = -t157 - t132;
t79 = -t132 - t133;
t71 = (t134 * t191 + t136 * t194) * t159;
t66 = -t235 * t134 + t214;
t65 = t109 + t82;
t64 = -t109 + t82;
t61 = t235 * t136 + t219;
t60 = -t136 * t250 - t191 * t82;
t59 = -t134 * t251 - t194 * t81;
t56 = -t107 * t191 - t263;
t55 = -t108 * t194 - t290;
t51 = -t191 * t91 - t263;
t50 = t194 * t91 - t265;
t49 = t194 * t80 - t290;
t48 = t191 * t80 + t289;
t45 = t199 + t271;
t40 = t191 * t65 + t194 * t207;
t39 = t191 * t207 - t194 * t65;
t38 = t191 * t61 - t194 * t64;
t37 = (t99 + t128) * pkin(4) + t199;
t36 = -t271 + t281 + t94 + 0.2e1 * t317;
t35 = -qJ(5) * t280 + t44;
t34 = -pkin(4) * t280 + t43;
t32 = t188 * t50 + t189 * t66;
t31 = t188 * t66 - t189 * t50;
t30 = t188 * t48 + t189 * t61;
t29 = t188 * t61 - t189 * t48;
t24 = t188 * t39 + t189 * t79;
t23 = t188 * t79 - t189 * t39;
t20 = pkin(5) * t39 - qJ(5) * t40;
t19 = t192 * t32 + t195 * t31;
t18 = t192 * t30 + t195 * t29;
t15 = t192 * t26 + t262;
t14 = t192 * t24 + t195 * t23;
t13 = pkin(5) * t66 - t273 * t51 - t266;
t12 = pkin(5) * t61 - t273 * t49 + t264;
t10 = -t191 * t206 - qJ(5) * t51 + (-t191 * t288 + t50) * pkin(5) - t204;
t9 = pkin(5) * t48 - qJ(5) * t49 - t16;
t6 = t188 * t7 + t189 * t28;
t5 = t188 * t28 - t189 * t7;
t4 = pkin(5) * t79 - t273 * t40 - t8;
t3 = pkin(5) * t7 - qJ(5) * t8;
t2 = pkin(5) * t28 - t273 * t8;
t1 = t192 * t6 + t195 * t5;
t22 = [0, 0, 0, 0, 0, qJDD(1), t223, t216, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t215 - 0.2e1 * t259, t181 + 0.2e1 * t183 - t216, pkin(1) * t160 + qJ(2) * (-pkin(1) * t198 + t181 - t213), (t173 - t226) * t195, -t171 * t195 - t174 * t192, t245 - t192 * (t197 - t247), (-t172 + t225) * t192, t195 * (-t197 + t248) - t246, 0, qJ(2) * t171 - t272 * t137 - t140 * t192, qJ(2) * t174 - t272 * t138 - t140 * t195, t272 * t175 - t242 * t261 - t85, -qJ(2) * t140 - t272 * t85, t202, t297, -t325, t201, -t320, t200, t195 * (-t268 - t311) - t192 * (-pkin(3) * t99 + t267 - t310) + qJ(2) * t99 + t324, t195 * (-t267 + t327) - t192 * (-pkin(3) * t287 - t268 - t326) + qJ(2) * t287 + t329, t195 * (-t25 - t318) - t192 * (t26 + t314) + t321, -qJ(4) * t262 - t192 * (pkin(3) * t94 + qJ(4) * t26) - qJ(2) * t94 - t272 * t15, t200, t325, t320, t202, t297, t201, t195 * (-t188 * t34 + t189 * t35 - t318) - t192 * (t188 * t35 + t189 * t34 + t314) + t321, t195 * (-t188 * t37 + t311) - t192 * (t189 * t37 + t310) - t324 - (-t192 * t224 - t195 * t260 + qJ(2)) * t99, t195 * (t189 * t36 - t327) - t192 * (t188 * t36 + t326) - t329 + (-pkin(4) * t244 - t192 * (pkin(3) + t270) - qJ(2)) * t287, (t195 * (pkin(4) * t188 - t260) - t192 * (t224 - t270) + qJ(2)) * t45 + (-t272 - qJ(4)) * t11, t195 * (-t188 * t60 + t228) - t192 * (t189 * t60 + t229), t195 * (-t188 * t38 + t189 * t95) - t192 * (t188 * t95 + t189 * t38), t195 * (-t188 * t55 + t189 * t65) - t192 * (t188 * t65 + t189 * t55), t195 * (-t188 * t59 - t228) - t192 * (t189 * t59 - t229), t195 * (-t188 * t56 + t189 * t207) - t192 * (t188 * t207 + t189 * t56), t195 * (t122 * t189 - t188 * t71) - t192 * (t122 * t188 + t189 * t71), t195 * (-qJ(4) * t29 - t12 * t188 + t189 * t9) - t192 * (-pkin(3) * t49 + qJ(4) * t30 + t12 * t189 + t188 * t9) + qJ(2) * t49 - t272 * t18, t195 * (-qJ(4) * t31 + t10 * t189 - t13 * t188) - t192 * (-pkin(3) * t51 + qJ(4) * t32 + t10 * t188 + t13 * t189) + qJ(2) * t51 - t272 * t19, t195 * (-qJ(4) * t23 - t188 * t4 + t189 * t20) - t192 * (-pkin(3) * t40 + qJ(4) * t24 + t188 * t20 + t189 * t4) + qJ(2) * t40 - t272 * t14, t195 * (-qJ(4) * t5 - t188 * t2 + t189 * t3) - t192 * (-pkin(3) * t8 + qJ(4) * t6 + t188 * t3 + t189 * t2) + qJ(2) * t8 - t272 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t198, -t160, 0, 0, 0, 0, 0, 0, t137, t138, -t175, t85, 0, 0, 0, 0, 0, 0, -t47, -t57, t313, t15, 0, 0, 0, 0, 0, 0, t313, t47, t57, t11, 0, 0, 0, 0, 0, 0, t18, t19, t14, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, (-t185 + t186) * t198, t179, -t227, -t231, qJDD(3), t130, t131, 0, 0, -t249, t279, t288, t249, t101, qJDD(3), -t53 + t312, -t54 - t328, t319, pkin(3) * t25, qJDD(3), -t288, -t101, -t249, t279, t249, -pkin(4) * t288 + qJ(5) * t101 + t319, -pkin(4) * t284 - qJ(5) * t111 - t312 + t44, t328 + pkin(4) * t147 + (qJDD(3) + t285) * qJ(5) + t205, pkin(3) * t21 - pkin(4) * t44 + qJ(5) * t43, -t136 * t251 + t194 * t82, -t191 * t64 - t194 * t61, -t108 * t191 + t289, t134 * t250 - t191 * t81, t107 * t194 - t265, (-t134 * t194 + t136 * t191) * t159, pkin(3) * t29 + qJ(5) * t61 - t273 * t48 + t266, pkin(3) * t31 + qJ(5) * t66 - t273 * t50 + t264, pkin(3) * t23 + qJ(5) * t79 - t273 * t39 - t7, pkin(3) * t5 + qJ(5) * t28 - t273 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t287, t280, -t94, 0, 0, 0, 0, 0, 0, t280, -t99, -t287, t45, 0, 0, 0, 0, 0, 0, t49, t51, t40, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288, t284, -t147, t44, 0, 0, 0, 0, 0, 0, t48, t50, t39, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t95, t65, -t98, t207, t122, -t16, -t17, 0, 0;];
tauJ_reg  = t22;
