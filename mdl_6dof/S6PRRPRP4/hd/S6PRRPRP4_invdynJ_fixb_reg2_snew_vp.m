% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 04:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRPRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:04:32
% EndTime: 2019-05-05 04:04:40
% DurationCPUTime: 2.50s
% Computational Cost: add. (6305->283), mult. (12734->364), div. (0->0), fcn. (8080->10), ass. (0->184)
t189 = sin(qJ(5));
t192 = cos(qJ(5));
t193 = cos(qJ(3));
t235 = qJD(2) * t193;
t140 = qJD(3) * t189 + t192 * t235;
t142 = qJD(3) * t192 - t189 * t235;
t106 = t142 * t140;
t232 = qJD(2) * qJD(3);
t225 = t193 * t232;
t190 = sin(qJ(3));
t230 = t190 * qJDD(2);
t146 = t225 + t230;
t135 = qJDD(5) + t146;
t273 = t106 - t135;
t281 = t273 * pkin(5);
t178 = t190 ^ 2;
t196 = qJD(2) ^ 2;
t172 = t178 * t196;
t195 = qJD(3) ^ 2;
t160 = -t172 - t195;
t238 = t193 * t196;
t227 = t190 * t238;
t156 = qJDD(3) - t227;
t239 = t193 * t156;
t111 = t160 * t190 + t239;
t145 = 0.2e1 * t225 + t230;
t184 = sin(pkin(6));
t186 = cos(pkin(6));
t191 = sin(qJ(2));
t194 = cos(qJ(2));
t280 = t186 * (t156 * t190 - t160 * t193) + (t111 * t191 + t145 * t194) * t184;
t279 = pkin(8) * t111;
t179 = t193 ^ 2;
t173 = t179 * t196;
t162 = -t173 - t195;
t155 = qJDD(3) + t227;
t247 = t155 * t190;
t110 = -t162 * t193 + t247;
t226 = t190 * t232;
t229 = t193 * qJDD(2);
t148 = -0.2e1 * t226 + t229;
t277 = (t110 * t191 - t148 * t194) * t184 - t186 * (t155 * t193 + t162 * t190);
t262 = -pkin(3) - pkin(9);
t276 = pkin(8) * t110;
t275 = t189 * t273;
t274 = t192 * t273;
t272 = t146 + t225;
t183 = sin(pkin(10));
t185 = cos(pkin(10));
t153 = g(1) * t183 - g(2) * t185;
t180 = -g(3) + qJDD(1);
t271 = t153 * t186 + t180 * t184;
t133 = t140 ^ 2;
t236 = qJD(2) * t190;
t166 = qJD(5) + t236;
t163 = t166 ^ 2;
t97 = -t163 - t133;
t55 = t189 * t97 - t274;
t56 = t192 * t97 + t275;
t269 = pkin(4) * t55 - qJ(4) * t56;
t134 = t142 ^ 2;
t103 = -t134 - t163;
t89 = t106 + t135;
t255 = t189 * t89;
t61 = t103 * t192 - t255;
t253 = t192 * t89;
t62 = -t103 * t189 - t253;
t268 = pkin(4) * t61 - qJ(4) * t62;
t119 = -t153 * t184 + t180 * t186;
t114 = t193 * t119;
t217 = -qJDD(3) * pkin(3) - qJ(4) * t195 + qJDD(4) - t114;
t240 = t190 * qJ(4);
t213 = -pkin(3) * t193 - t240;
t143 = t213 * qJD(2);
t154 = -g(1) * t185 - g(2) * t183;
t92 = t194 * t154 + t191 * t271;
t83 = -t196 * pkin(2) + qJDD(2) * pkin(8) + t92;
t222 = qJD(2) * t143 + t83;
t37 = -qJDD(3) * pkin(9) + (t146 - t225) * pkin(4) + (-pkin(9) * t238 + t222) * t190 + t217;
t147 = -t226 + t229;
t158 = pkin(4) * t236 - qJD(3) * pkin(9);
t215 = t191 * t154 - t194 * t271;
t82 = -qJDD(2) * pkin(2) - pkin(8) * t196 + t215;
t199 = -t147 * pkin(3) - qJ(4) * t272 + t82;
t223 = pkin(3) * qJD(3) - (2 * qJD(4));
t39 = -pkin(4) * t173 - t147 * pkin(9) + (-t158 + t223) * t236 + t199;
t19 = t189 * t39 - t192 * t37;
t244 = t166 * t140;
t264 = qJ(6) * t244 + 0.2e1 * qJD(6) * t142 + t19 + t281;
t115 = pkin(5) * t166 - qJ(6) * t142;
t20 = t189 * t37 + t192 * t39;
t221 = t189 * qJDD(3) + t147 * t192;
t95 = -qJD(5) * t142 - t221;
t16 = -pkin(5) * t133 + qJ(6) * t95 - 0.2e1 * qJD(6) * t140 - t115 * t166 + t20;
t263 = t239 + (t173 - t195) * t190;
t201 = (-qJD(5) + t166) * t142 - t221;
t96 = -qJD(5) * t140 + qJDD(3) * t192 - t147 * t189;
t78 = -t96 - t244;
t43 = t189 * t201 + t192 * t78;
t87 = -t133 - t134;
t25 = t190 * t43 + t193 * t87;
t44 = -t189 * t78 + t192 * t201;
t261 = -pkin(2) * t44 + pkin(8) * t25;
t74 = (qJD(5) + t166) * t142 + t221;
t29 = t190 * t55 + t193 * t74;
t260 = -pkin(2) * t56 + pkin(8) * t29;
t205 = t96 - t244;
t31 = t190 * t61 + t193 * t205;
t259 = -pkin(2) * t62 + pkin(8) * t31;
t231 = qJD(4) * qJD(3);
t174 = 0.2e1 * t231;
t64 = t119 * t190 + t193 * t83;
t214 = -pkin(3) * t195 + qJDD(3) * qJ(4) + t143 * t235 + t64;
t200 = pkin(4) * t147 - pkin(9) * t173 + qJD(3) * t158 + t214;
t36 = t174 + t200;
t256 = t189 * t36;
t254 = t192 * t36;
t252 = qJ(6) * t189;
t251 = qJ(6) * t192;
t243 = t166 * t189;
t242 = t166 * t192;
t150 = (t178 + t179) * qJDD(2);
t151 = t172 + t173;
t237 = pkin(2) * t151 + pkin(8) * t150;
t228 = t190 * t106;
t224 = pkin(4) * t43 - qJ(4) * t44;
t63 = t190 * t83 - t114;
t26 = t190 * t63 + t193 * t64;
t220 = qJ(4) * t87 + t262 * t43;
t219 = qJ(4) * t74 + t262 * t55;
t218 = qJ(4) * t205 + t262 * t61;
t6 = t189 * t20 - t19 * t192;
t4 = t190 * t6 + t193 * t36;
t7 = t189 * t19 + t192 * t20;
t212 = (-t172 + t195) * t193 + t247;
t209 = pkin(4) * t74 + t262 * t56;
t208 = pkin(4) * t205 + t262 * t62;
t207 = pkin(4) * t87 + t262 * t44;
t206 = pkin(2) - t213;
t46 = t174 + t214;
t203 = pkin(5) * t103 - t16;
t202 = t193 * t262 - pkin(2) - t240;
t14 = -qJ(6) * t96 - t264;
t48 = t222 * t190 + t217;
t198 = t14 - t281;
t197 = -pkin(5) * t95 - qJ(6) * t133 + t115 * t142 + qJDD(6) + t200;
t22 = t174 + t197;
t51 = t223 * t236 + t199;
t152 = t172 - t173;
t117 = -t134 + t163;
t116 = t133 - t163;
t108 = t272 * t190;
t107 = (t147 - t226) * t193;
t104 = t134 - t133;
t101 = t145 * t193 + t148 * t190;
t98 = (t150 * t191 + t151 * t194) * t184;
t81 = (-t140 * t192 + t142 * t189) * t166;
t73 = pkin(5) * t78;
t69 = -t142 * t243 + t192 * t96;
t68 = t140 * t242 - t189 * t95;
t67 = t190 * t135 + t193 * (t140 * t189 + t142 * t192) * t166;
t66 = t116 * t192 - t255;
t65 = -t117 * t189 - t274;
t50 = t228 + t193 * (-t142 * t242 - t189 * t96);
t49 = -t228 + t193 * (-t140 * t243 - t192 * t95);
t47 = -pkin(5) * t205 - qJ(6) * t89;
t45 = -t189 * t205 - t192 * t74;
t35 = t190 * t201 + t193 * (-t116 * t189 - t253);
t34 = -t190 * t78 + t193 * (-t117 * t192 + t275);
t27 = t190 * t104 + t193 * (t189 * t74 - t192 * t205);
t23 = t190 * t48 + t193 * t46;
t21 = -qJ(6) * t103 + t22;
t17 = -pkin(5) * t74 + qJ(6) * t97 - t197 - 0.2e1 * t231;
t15 = t186 * (t190 * t205 - t193 * t61) + (t191 * t31 - t194 * t62) * t184;
t13 = pkin(5) * t14;
t11 = t186 * (t190 * t74 - t193 * t55) + (t191 * t29 - t194 * t56) * t184;
t10 = (-t78 + t96) * qJ(6) + t264;
t9 = -pkin(5) * t87 + qJ(6) * t201 + t16;
t8 = t186 * (t190 * t87 - t193 * t43) + (t191 * t25 - t194 * t44) * t184;
t5 = -pkin(5) * t22 + qJ(6) * t16;
t3 = -t14 * t189 + t16 * t192;
t2 = t14 * t192 + t16 * t189;
t1 = t190 * t2 + t193 * t22;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t180, 0, 0, 0, 0, 0, 0, (qJDD(2) * t194 - t191 * t196) * t184, (-qJDD(2) * t191 - t194 * t196) * t184, 0, t186 * t119 + (t191 * t92 - t194 * t215) * t184, 0, 0, 0, 0, 0, 0, -t277, -t280, t98, t186 * (t190 * t64 - t193 * t63) + (t191 * t26 - t194 * t82) * t184, 0, 0, 0, 0, 0, 0, t98, t277, t280, t186 * (t190 * t46 - t193 * t48) + (t191 * t23 - t194 * t51) * t184, 0, 0, 0, 0, 0, 0, t11, t15, t8, t186 * (t190 * t36 - t193 * t6) + (t191 * t4 - t194 * t7) * t184, 0, 0, 0, 0, 0, 0, t11, t15, t8, t186 * (t190 * t22 - t193 * t2) + (t1 * t191 - t194 * t3) * t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t215, -t92, 0, 0, t108, t101, t212, t107, t263, 0, pkin(2) * t148 - t193 * t82 - t276, -pkin(2) * t145 + t190 * t82 - t279, t26 + t237, -pkin(2) * t82 + pkin(8) * t26, 0, -t212, -t263, t108, t101, t107, t193 * (pkin(3) * t151 + t46) + (qJ(4) * t151 + t48) * t190 + t237, -t148 * t206 + t193 * t51 + t276, t190 * (-pkin(3) * t226 + 0.2e1 * qJD(4) * t236 - t199) + t279 + t206 * t145, pkin(8) * t23 - t206 * t51, t50, t27, t34, t49, t35, t67, t190 * (-t19 + t269) + t193 * (t209 + t254) + t260, t190 * (-t20 + t268) + t193 * (t208 - t256) + t259, t190 * t224 + t193 * (t207 - t7) + t261, t202 * t7 + (pkin(4) + pkin(8)) * t4, t50, t27, t34, t49, t35, t67, t190 * (t198 + t269) + t193 * (-t192 * t17 - t252 * t273 + t209) + t260, t190 * (t203 + t268) + t193 * (-t189 * t21 - t192 * t47 + t208) + t259, t190 * (t224 + t73) + t193 * (-t189 * t10 - t192 * t9 + t207) + t261, t190 * (pkin(4) * t2 + t13) + t193 * (pkin(4) * t22 + t14 * t252 - t192 * t5) + pkin(8) * t1 + t202 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t227, t152, t230, t227, t229, qJDD(3), -t63, -t64, 0, 0, qJDD(3), -t230, -t229, -t227, t152, t227, (-pkin(3) * t190 + qJ(4) * t193) * qJDD(2), -pkin(3) * t155 - qJ(4) * t162 + t48, -pkin(3) * t160 + qJ(4) * t156 + t46, -pkin(3) * t48 + qJ(4) * t46, t69, t45, t65, t68, t66, t81, t219 + t256, t218 + t254, t220 - t6, qJ(4) * t36 + t262 * t6, t69, t45, t65, t68, t66, t81, -t17 * t189 + t251 * t273 + t219, -t189 * t47 + t192 * t21 + t218, t10 * t192 - t189 * t9 + t220, qJ(4) * t22 - t14 * t251 - t189 * t5 + t2 * t262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230, t155, t160, t48, 0, 0, 0, 0, 0, 0, t55, t61, t43, t6, 0, 0, 0, 0, 0, 0, t55, t61, t43, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t104, -t78, -t106, t201, t135, -t19, -t20, 0, 0, t106, t104, -t78, -t106, t201, t135, t198, t203, t73, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t205, t87, t22;];
tauJ_reg  = t12;
