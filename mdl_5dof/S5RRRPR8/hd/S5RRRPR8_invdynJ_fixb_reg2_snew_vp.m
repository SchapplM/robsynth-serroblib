% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:26
% EndTime: 2019-12-31 21:20:39
% DurationCPUTime: 5.25s
% Computational Cost: add. (10568->360), mult. (22290->465), div. (0->0), fcn. (15150->8), ass. (0->219)
t205 = cos(qJ(3));
t197 = qJD(2) + qJD(3);
t195 = t197 ^ 2;
t201 = sin(qJ(3));
t206 = cos(qJ(2));
t202 = sin(qJ(2));
t253 = t202 * t205;
t172 = (t201 * t206 + t253) * qJD(1);
t284 = t172 ^ 2;
t289 = -t284 - t195;
t196 = qJDD(2) + qJDD(3);
t250 = qJD(1) * t202;
t170 = -t205 * t206 * qJD(1) + t201 * t250;
t259 = t172 * t170;
t302 = t259 + t196;
t312 = t302 * t201;
t101 = -t205 * t289 + t312;
t311 = t302 * t205;
t103 = t201 * t289 + t311;
t326 = pkin(6) * (t101 * t202 - t103 * t206);
t325 = pkin(2) * t101;
t324 = pkin(7) * t101;
t323 = pkin(7) * t103;
t155 = -t284 + t195;
t303 = t196 - t259;
t313 = t205 * t303;
t314 = t201 * t303;
t321 = t202 * (t155 * t201 - t313) - t206 * (t155 * t205 + t314);
t285 = t170 ^ 2;
t154 = t285 - t195;
t318 = t202 * (-t154 * t205 + t312) - t206 * (t154 * t201 + t311);
t129 = -t195 - t285;
t84 = t129 * t201 + t313;
t87 = -t129 * t205 + t314;
t315 = pkin(6) * (t202 * t84 + t206 * t87);
t310 = pkin(2) * t84;
t309 = pkin(7) * t84;
t308 = pkin(7) * t87;
t283 = 2 * qJD(4);
t282 = pkin(3) + pkin(8);
t287 = -t285 - t284;
t300 = pkin(1) * t287;
t299 = pkin(2) * t287;
t200 = sin(qJ(5));
t204 = cos(qJ(5));
t145 = -t204 * t170 + t197 * t200;
t147 = t170 * t200 + t197 * t204;
t110 = t147 * t145;
t190 = t202 * qJDD(1);
t245 = qJD(1) * qJD(2);
t240 = t206 * t245;
t179 = t190 + t240;
t191 = t206 * qJDD(1);
t241 = t202 * t245;
t180 = t191 - t241;
t228 = t205 * t179 + t201 * t180;
t119 = -qJD(3) * t170 + t228;
t116 = qJDD(5) + t119;
t292 = -t110 + t116;
t298 = t200 * t292;
t296 = t204 * t292;
t199 = t206 ^ 2;
t208 = qJD(1) ^ 2;
t203 = sin(qJ(1));
t281 = cos(qJ(1));
t238 = t203 * g(1) - t281 * g(2);
t223 = qJDD(1) * pkin(1) + t238;
t224 = qJD(2) * pkin(2) - pkin(7) * t250;
t122 = t180 * pkin(2) + (pkin(7) * t199 + pkin(6)) * t208 - t224 * t250 + t223;
t258 = t197 * t172;
t294 = pkin(3) * t258 - t172 * t283 - t122;
t136 = pkin(3) * t170 - qJ(4) * t172;
t252 = t202 * t208;
t225 = t281 * g(1) + t203 * g(2);
t270 = qJDD(1) * pkin(6);
t175 = -t208 * pkin(1) - t225 + t270;
t255 = t202 * t175;
t107 = qJDD(2) * pkin(2) - t179 * pkin(7) - t255 + (pkin(2) * t252 + pkin(7) * t245 - g(3)) * t206;
t152 = -t202 * g(3) + t206 * t175;
t193 = t199 * t208;
t109 = -pkin(2) * t193 + t180 * pkin(7) - qJD(2) * t224 + t152;
t73 = -t205 * t107 + t201 * t109;
t48 = -t196 * pkin(3) - t195 * qJ(4) + t172 * t136 + qJDD(4) + t73;
t286 = -pkin(8) * t303 + t48;
t159 = t170 * t197;
t291 = t119 + t159;
t211 = pkin(4) * t291 + t286;
t233 = t201 * t179 - t205 * t180;
t118 = qJD(3) * t172 + t233;
t153 = pkin(4) * t172 - pkin(8) * t197;
t99 = t119 - t159;
t293 = qJ(4) * t99;
t209 = -t293 + t294;
t23 = -pkin(4) * t285 + t282 * t118 - t172 * t153 + t209;
t14 = t200 * t23 - t204 * t211;
t272 = t204 * t23;
t15 = t200 * t211 + t272;
t7 = -t204 * t14 + t15 * t200;
t288 = t284 - t285;
t142 = t145 ^ 2;
t143 = t147 ^ 2;
t165 = qJD(5) + t172;
t163 = t165 ^ 2;
t280 = pkin(3) * t201;
t279 = pkin(3) * t205;
t74 = t201 * t107 + t205 * t109;
t216 = -t195 * pkin(3) + t196 * qJ(4) - t170 * t136 + t74;
t43 = t197 * t283 + t216;
t278 = -pkin(3) * t48 + qJ(4) * t43;
t248 = -qJD(3) + t197;
t94 = t248 * t172 - t233;
t95 = t248 * t170 + t228;
t277 = -pkin(3) * t95 + qJ(4) * t94;
t29 = -t118 * pkin(4) - t285 * pkin(8) + (t283 + t153) * t197 + t216;
t275 = t200 * t29;
t76 = t110 + t116;
t274 = t200 * t76;
t37 = t201 * t74 - t205 * t73;
t273 = t202 * t37;
t27 = t204 * t29;
t271 = t204 * t76;
t269 = t122 * t201;
t268 = t122 * t205;
t261 = t165 * t200;
t260 = t165 * t204;
t257 = t197 * t201;
t256 = t197 * t205;
t185 = t206 * t252;
t254 = t202 * (qJDD(2) + t185);
t251 = t206 * (qJDD(2) - t185);
t247 = qJD(3) + t197;
t246 = qJD(5) + t165;
t244 = t201 * t110;
t243 = t205 * t110;
t239 = qJ(4) * t201 + pkin(2);
t38 = t201 * t73 + t205 * t74;
t237 = -t204 * t118 + t200 * t196;
t151 = t206 * g(3) + t255;
t234 = t202 * t151 + t206 * t152;
t232 = qJ(4) * t29 - t282 * t7;
t105 = -t143 - t163;
t50 = t105 * t204 - t274;
t229 = t200 * t118 + t204 * t196;
t71 = -t246 * t145 + t229;
t231 = qJ(4) * t71 - t282 * t50 + t27;
t8 = t200 * t14 + t204 * t15;
t88 = -t163 - t142;
t45 = t200 * t88 + t296;
t66 = t246 * t147 + t237;
t222 = qJ(4) * t66 - t282 * t45 + t275;
t220 = (-qJD(5) + t165) * t147 - t237;
t125 = t165 * t145;
t81 = -qJD(5) * t145 + t229;
t70 = t125 + t81;
t34 = t200 * t220 - t204 * t70;
t83 = -t142 - t143;
t221 = qJ(4) * t83 - t282 * t34 - t7;
t218 = -t247 * t170 + t228;
t217 = -pkin(3) * t303 - qJ(4) * t129 + t48;
t215 = t202 * (t205 * t119 - t172 * t257) + t206 * (t201 * t119 + t172 * t256);
t214 = t202 * (t118 * t201 + t170 * t256) + t206 * (-t205 * t118 + t170 * t257);
t213 = -pkin(3) * t289 + qJ(4) * t302 + t43;
t212 = (t202 * (-t170 * t205 + t172 * t201) + t206 * (-t170 * t201 - t172 * t205)) * t197;
t210 = -t118 * pkin(3) - t294;
t207 = qJD(2) ^ 2;
t198 = t202 ^ 2;
t192 = t198 * t208;
t181 = t191 - 0.2e1 * t241;
t178 = t190 + 0.2e1 * t240;
t174 = t208 * pkin(6) + t223;
t124 = -t143 + t163;
t123 = t142 - t163;
t108 = t143 - t142;
t93 = t118 + t258;
t92 = t118 - t258;
t91 = t247 * t172 + t233;
t80 = -qJD(5) * t147 - t237;
t79 = (-t145 * t204 + t147 * t200) * t165;
t78 = (t145 * t200 + t147 * t204) * t165;
t69 = -t125 + t81;
t63 = -t147 * t261 + t204 * t81;
t62 = -t147 * t260 - t200 * t81;
t61 = t145 * t260 - t200 * t80;
t60 = -t145 * t261 - t204 * t80;
t59 = t201 * t291 - t205 * t92;
t58 = t201 * t95 + t205 * t94;
t57 = -t201 * t92 - t205 * t291;
t56 = t201 * t94 - t205 * t95;
t55 = t123 * t204 - t274;
t54 = -t124 * t200 + t296;
t53 = -t123 * t200 - t271;
t52 = -t124 * t204 - t298;
t51 = -t105 * t200 - t271;
t46 = t204 * t88 - t298;
t40 = -qJ(4) * t287 + t48;
t39 = -pkin(3) * t287 + t43;
t36 = t200 * t70 + t204 * t220;
t35 = -t200 * t69 - t204 * t66;
t33 = t200 * t66 - t204 * t69;
t31 = (t118 + t93) * pkin(3) + t209;
t30 = (t99 + t218) * qJ(4) + t210;
t25 = t201 * t50 + t205 * t71;
t24 = t201 * t71 - t205 * t50;
t22 = t201 * t45 + t205 * t66;
t21 = t201 * t66 - t205 * t45;
t20 = t201 * t34 + t205 * t83;
t19 = t201 * t83 - t205 * t34;
t17 = t201 * t43 - t205 * t48;
t16 = pkin(4) * t34 - qJ(4) * t36;
t12 = pkin(4) * t71 - t282 * t51 - t275;
t11 = pkin(4) * t66 - t282 * t46 + t27;
t10 = -t272 - t200 * t286 - qJ(4) * t51 + (-t200 * t291 + t50) * pkin(4);
t9 = pkin(4) * t45 - qJ(4) * t46 - t14;
t5 = t201 * t7 + t205 * t29;
t4 = t201 * t29 - t205 * t7;
t3 = pkin(4) * t83 - t282 * t36 - t8;
t2 = pkin(4) * t7 - qJ(4) * t8;
t1 = pkin(4) * t29 - t282 * t8;
t6 = [0, 0, 0, 0, 0, qJDD(1), t238, t225, 0, 0, (t179 + t240) * t202, t178 * t206 + t181 * t202, t254 + t206 * (-t192 + t207), (t180 - t241) * t206, t202 * (t193 - t207) + t251, 0, t206 * t174 + pkin(1) * t181 + pkin(6) * (t206 * (-t193 - t207) - t254), -t202 * t174 - pkin(1) * t178 + pkin(6) * (-t251 - t202 * (-t192 - t207)), pkin(1) * (t192 + t193) + (t198 + t199) * t270 + t234, pkin(1) * t174 + pkin(6) * t234, t215, t202 * (-t201 * t99 - t205 * t91) + t206 * (-t201 * t91 + t205 * t99), -t321, t214, -t318, t212, t202 * (-t269 - t309) + t206 * (-pkin(2) * t91 + t268 - t308) - pkin(1) * t91 - t315, t202 * (-t268 + t324) + t206 * (-pkin(2) * t218 - t269 - t323) - pkin(1) * t218 + t326, t202 * (-pkin(7) * t57 - t37) + t206 * (pkin(7) * t59 - t299 + t38) - t300 + pkin(6) * (-t202 * t57 + t206 * t59), -pkin(7) * t273 + t206 * (pkin(2) * t122 + pkin(7) * t38) + pkin(1) * t122 + pkin(6) * (t206 * t38 - t273), t212, t321, t318, t215, t202 * (-t201 * t218 - t205 * t93) + t206 * (-t201 * t93 + t205 * t218), t214, t202 * (-pkin(7) * t56 - t201 * t39 + t205 * t40) + t206 * (pkin(7) * t58 + t201 * t40 + t205 * t39 - t299) - t300 + pkin(6) * (-t202 * t56 + t206 * t58), t202 * (-t201 * t31 + t309) + t206 * (t205 * t31 + t308) + t315 + (qJ(4) * t253 + t206 * t239 + pkin(1)) * t93, t202 * (t205 * t30 - t324) + t206 * (t201 * t30 + t323) - t326 + (-t202 * t280 + t206 * (pkin(2) + t279) + pkin(1)) * t218, (t202 * (qJ(4) * t205 - t280) + t206 * (t239 + t279) + pkin(1)) * (t210 + t293) + (pkin(6) + pkin(7)) * (-t17 * t202 + (t201 * t48 + t205 * t43) * t206), t202 * (-t201 * t62 + t243) + t206 * (t205 * t62 + t244), t202 * (t108 * t205 - t201 * t33) + t206 * (t108 * t201 + t205 * t33), t202 * (-t201 * t52 + t205 * t70) + t206 * (t201 * t70 + t205 * t52), t202 * (-t201 * t60 - t243) + t206 * (t205 * t60 - t244), t202 * (-t201 * t53 + t205 * t220) + t206 * (t201 * t220 + t205 * t53), t202 * (t116 * t205 - t201 * t78) + t206 * (t116 * t201 + t205 * t78), t202 * (-pkin(7) * t21 - t11 * t201 + t205 * t9) + t206 * (-pkin(2) * t46 + pkin(7) * t22 + t11 * t205 + t201 * t9) - pkin(1) * t46 + pkin(6) * (-t202 * t21 + t206 * t22), t202 * (-pkin(7) * t24 + t10 * t205 - t12 * t201) + t206 * (-pkin(2) * t51 + pkin(7) * t25 + t10 * t201 + t12 * t205) - pkin(1) * t51 + pkin(6) * (-t202 * t24 + t206 * t25), t202 * (-pkin(7) * t19 + t16 * t205 - t201 * t3) + t206 * (-pkin(2) * t36 + pkin(7) * t20 + t16 * t201 + t205 * t3) - pkin(1) * t36 + pkin(6) * (-t19 * t202 + t20 * t206), t202 * (-pkin(7) * t4 - t1 * t201 + t2 * t205) + t206 * (-pkin(2) * t8 + pkin(7) * t5 + t1 * t205 + t2 * t201) - pkin(1) * t8 + pkin(6) * (-t202 * t4 + t206 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185, t192 - t193, t190, t185, t191, qJDD(2), -t151, -t152, 0, 0, t259, t288, t291, -t259, t94, t196, -t73 + t310, -t74 - t325, pkin(2) * t57, pkin(2) * t37, t196, -t95, t92, t259, t288, -t259, pkin(2) * t56 + t277, t217 - t310, t213 + t325, pkin(2) * t17 + t278, t63, t35, t54, t61, t55, t79, pkin(2) * t21 + t222, pkin(2) * t24 + t231, pkin(2) * t19 + t221, pkin(2) * t4 + t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t259, t288, t291, -t259, t94, t196, -t73, -t74, 0, 0, t196, -t95, t92, t259, t288, -t259, t277, t217, t213, t278, t63, t35, t54, t61, t55, t79, t222, t231, t221, t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, t303, t289, t48, 0, 0, 0, 0, 0, 0, t45, t50, t34, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, t108, t70, -t110, t220, t116, -t14, -t15, 0, 0;];
tauJ_reg = t6;
