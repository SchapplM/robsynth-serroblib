% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRPR6
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:55
% EndTime: 2019-12-05 16:33:08
% DurationCPUTime: 5.00s
% Computational Cost: add. (4020->473), mult. (9701->684), div. (0->0), fcn. (7593->14), ass. (0->220)
t169 = sin(qJ(3));
t172 = cos(qJ(3));
t196 = pkin(3) * t169 - qJ(4) * t172;
t101 = qJD(3) * t196 - t169 * qJD(4);
t163 = sin(pkin(10));
t165 = cos(pkin(10));
t170 = sin(qJ(2));
t241 = qJD(3) * t169;
t232 = pkin(7) * t241;
t164 = sin(pkin(5));
t250 = qJD(1) * t164;
t173 = cos(qJ(2));
t254 = t172 * t173;
t278 = t165 * t101 + t163 * t232 - (-t163 * t254 + t165 * t170) * t250;
t261 = t163 * t170;
t298 = t163 * t101 - (t165 * t254 + t261) * t250;
t242 = qJD(3) * t165;
t246 = qJD(2) * t169;
t119 = -t163 * t246 + t242;
t243 = qJD(3) * t163;
t120 = t165 * t246 + t243;
t168 = sin(qJ(5));
t171 = cos(qJ(5));
t51 = -t171 * t119 + t120 * t168;
t297 = t51 ^ 2;
t54 = t119 * t168 + t120 * t171;
t296 = t54 ^ 2;
t256 = t165 * t172;
t189 = pkin(4) * t169 - pkin(8) * t256;
t295 = qJD(3) * t189 + t278;
t257 = t165 * t169;
t260 = t163 * t172;
t294 = (-pkin(7) * t257 - pkin(8) * t260) * qJD(3) + t298;
t244 = qJD(2) * t172;
t152 = -qJD(5) + t244;
t293 = t152 * t51;
t237 = qJD(1) * qJD(2);
t218 = t173 * t237;
t166 = cos(pkin(5));
t249 = qJD(1) * t166;
t267 = qJDD(2) * pkin(7);
t292 = qJD(3) * t249 + t267 + (qJDD(1) * t170 + t218) * t164;
t158 = t172 * qJDD(2);
t236 = qJD(2) * qJD(3);
t291 = t169 * t236 - t158;
t238 = qJD(5) * t171;
t239 = qJD(5) * t168;
t290 = -t163 * t239 + t165 * t238;
t123 = t163 * t168 - t171 * t165;
t197 = pkin(3) * t172 + qJ(4) * t169;
t161 = t169 ^ 2;
t162 = t172 ^ 2;
t252 = t161 - t162;
t203 = qJD(2) * t252;
t289 = -qJDD(3) * pkin(3) + qJDD(4);
t174 = qJD(3) ^ 2;
t270 = sin(pkin(9));
t205 = t270 * t170;
t271 = cos(pkin(9));
t206 = t271 * t173;
t103 = -t166 * t206 + t205;
t204 = t270 * t173;
t207 = t271 * t170;
t105 = t166 * t204 + t207;
t200 = g(1) * t105 + g(2) * t103;
t219 = t170 * t237;
t268 = qJDD(2) * pkin(2);
t258 = t164 * t173;
t195 = -qJDD(1) * t258 + t164 * t219;
t92 = t195 - t268;
t288 = -pkin(7) * t174 + t164 * (-g(3) * t173 + t219) + t200 + t268 - t92;
t127 = qJD(2) * pkin(7) + t170 * t250;
t115 = t169 * t127;
t235 = qJDD(1) * t166;
t231 = t169 * t235 + t292 * t172;
t26 = qJDD(3) * qJ(4) + (qJD(4) - t115) * qJD(3) + t231;
t129 = -pkin(2) - t197;
t37 = qJD(2) * t101 + qJDD(2) * t129 + t195;
t13 = t163 * t37 + t165 * t26;
t233 = t169 * qJDD(2);
t143 = t163 * t233;
t215 = t172 * t236;
t86 = -qJDD(3) * t165 + t163 * t215 + t143;
t11 = -pkin(8) * t86 + t13;
t248 = qJD(1) * t169;
t221 = t166 * t248;
t84 = t127 * t172 + t221;
t75 = qJD(3) * qJ(4) + t84;
t247 = qJD(1) * t173;
t229 = t164 * t247;
t85 = qJD(2) * t129 - t229;
t31 = -t163 * t75 + t165 * t85;
t19 = -pkin(4) * t244 - pkin(8) * t120 + t31;
t32 = t163 * t85 + t165 * t75;
t20 = pkin(8) * t119 + t32;
t194 = t168 * t20 - t171 * t19;
t12 = -t163 * t26 + t165 * t37;
t234 = t163 * qJDD(3);
t87 = t234 + (t215 + t233) * t165;
t6 = t291 * pkin(4) - pkin(8) * t87 + t12;
t1 = -t194 * qJD(5) + t171 * t11 + t168 * t6;
t287 = pkin(4) * t86;
t118 = t165 * t129;
t55 = -pkin(8) * t257 + t118 + (-pkin(7) * t163 - pkin(4)) * t172;
t89 = pkin(7) * t256 + t163 * t129;
t68 = -pkin(8) * t163 * t169 + t89;
t21 = -t168 * t68 + t171 * t55;
t286 = qJD(5) * t21 + t295 * t168 + t294 * t171;
t22 = t168 * t55 + t171 * t68;
t285 = -qJD(5) * t22 - t294 * t168 + t295 * t171;
t283 = pkin(4) * t163;
t282 = t54 * t51;
t281 = pkin(8) + qJ(4);
t125 = t196 * qJD(2);
t83 = t172 * t249 - t115;
t44 = t165 * t125 - t163 * t83;
t30 = qJD(2) * t189 + t44;
t226 = t163 * t244;
t45 = t163 * t125 + t165 * t83;
t38 = -pkin(8) * t226 + t45;
t132 = t281 * t163;
t133 = t281 * t165;
t69 = -t132 * t171 - t133 * t168;
t280 = -qJD(4) * t123 + qJD(5) * t69 - t168 * t30 - t171 * t38;
t124 = t163 * t171 + t165 * t168;
t70 = -t132 * t168 + t133 * t171;
t279 = -qJD(4) * t124 - qJD(5) * t70 + t168 * t38 - t171 * t30;
t277 = -t165 * t232 + t298;
t276 = qJD(2) * pkin(2);
t275 = t163 * t87;
t274 = t165 * t86;
t273 = -t123 * t244 - t290;
t108 = t124 * qJD(5);
t185 = t124 * t172;
t272 = -qJD(2) * t185 + t108;
t160 = pkin(10) + qJ(5);
t156 = sin(t160);
t265 = t156 * t172;
t157 = cos(t160);
t264 = t157 * t172;
t175 = qJD(2) ^ 2;
t263 = t162 * t175;
t259 = t164 * t170;
t253 = qJDD(1) - g(3);
t251 = t161 + t162;
t245 = qJD(2) * t170;
t240 = qJD(3) * t172;
t230 = pkin(7) + t283;
t228 = t169 * t247;
t227 = t120 * t244;
t225 = t164 * t245;
t224 = qJD(2) * t258;
t220 = g(3) * (pkin(2) * t258 + pkin(7) * t259);
t214 = t173 * t236;
t213 = t163 * t158;
t212 = t165 * t233;
t211 = t165 * t158;
t210 = t168 * t87 + t171 * t86;
t209 = t164 * t271;
t208 = t164 * t270;
t202 = t127 * t240 + t292 * t169 - t172 * t235;
t201 = t169 * t215;
t104 = t166 * t207 + t204;
t106 = -t166 * t205 + t206;
t199 = g(1) * t106 + g(2) * t104;
t10 = t168 * t19 + t171 * t20;
t111 = t166 * t169 + t172 * t259;
t60 = -t111 * t163 - t165 * t258;
t61 = t111 * t165 - t163 * t258;
t24 = -t168 * t61 + t171 * t60;
t25 = t168 * t60 + t171 * t61;
t193 = t119 * t165 - t120 * t163;
t155 = pkin(4) * t165 + pkin(3);
t192 = t155 * t172 + t169 * t281;
t191 = qJD(2) * (-t119 + t242);
t190 = qJDD(2) * t173 - t170 * t175;
t110 = -t166 * t172 + t169 * t259;
t16 = -t119 * t238 + t120 * t239 + t168 * t86 - t171 * t87;
t62 = t104 * t169 + t172 * t209;
t64 = t106 * t169 - t172 * t208;
t187 = g(1) * t64 + g(2) * t62 + g(3) * t110;
t63 = t104 * t172 - t169 * t209;
t65 = t106 * t172 + t169 * t208;
t186 = g(1) * t65 + g(2) * t63 + g(3) * t111;
t27 = t202 + t289;
t184 = t187 - t27;
t71 = -qJD(3) * pkin(3) + qJD(4) - t83;
t183 = -qJ(4) * t241 + (qJD(4) - t71) * t172;
t182 = g(3) * t258 - t200;
t181 = -g(3) * t259 - t199;
t2 = -qJD(5) * t10 - t168 * t11 + t171 * t6;
t179 = t187 - t202;
t17 = qJD(5) * t54 + t210;
t128 = -t229 - t276;
t178 = -pkin(7) * qJDD(3) + (t128 + t229 - t276) * qJD(3);
t177 = -t179 + t289;
t28 = -t127 * t241 + t231;
t176 = t202 * t169 + t28 * t172 + (-t169 * t84 - t172 * t83) * qJD(3) - t199;
t148 = t169 * t175 * t172;
t126 = t230 * t169;
t122 = qJDD(5) + t291;
t116 = qJDD(2) * t162 - 0.2e1 * t201;
t114 = t230 * t240;
t100 = t105 * pkin(2);
t99 = t103 * pkin(2);
t97 = t123 * t169;
t96 = t124 * t169;
t88 = -pkin(7) * t260 + t118;
t67 = qJD(3) * t111 + t169 * t224;
t66 = -qJD(3) * t110 + t172 * t224;
t56 = t221 + (qJD(2) * t283 + t127) * t172;
t48 = qJD(3) * t185 + t290 * t169;
t47 = t108 * t169 + t123 * t240;
t46 = -pkin(4) * t119 + t71;
t43 = t163 * t225 + t165 * t66;
t42 = -t163 * t66 + t165 * t225;
t18 = t27 + t287;
t5 = -qJD(5) * t25 - t168 * t43 + t171 * t42;
t4 = qJD(5) * t24 + t168 * t42 + t171 * t43;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t253, 0, 0, 0, 0, 0, 0, t190 * t164, (-qJDD(2) * t170 - t173 * t175) * t164, 0, -g(3) + (t166 ^ 2 + (t170 ^ 2 + t173 ^ 2) * t164 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -qJD(3) * t67 - qJDD(3) * t110 + (-t169 * t214 + t172 * t190) * t164, -qJD(3) * t66 - qJDD(3) * t111 + (-t169 * t190 - t172 * t214) * t164, (t110 * t169 + t111 * t172) * qJDD(2) + (t169 * t67 + t172 * t66 + (t110 * t172 - t111 * t169) * qJD(3)) * qJD(2), t110 * t202 + t111 * t28 + t66 * t84 - t67 * t83 - g(3) + (t128 * t245 - t173 * t92) * t164, 0, 0, 0, 0, 0, 0, -t60 * t158 + t110 * t86 - t119 * t67 + (-t172 * t42 + t241 * t60) * qJD(2), t61 * t158 + t110 * t87 + t120 * t67 + (t172 * t43 - t241 * t61) * qJD(2), t119 * t43 - t120 * t42 - t60 * t87 - t61 * t86, t110 * t27 + t12 * t60 + t13 * t61 + t31 * t42 + t32 * t43 + t67 * t71 - g(3), 0, 0, 0, 0, 0, 0, t110 * t17 + t122 * t24 - t152 * t5 + t51 * t67, -t110 * t16 - t122 * t25 + t152 * t4 + t54 * t67, t16 * t24 - t17 * t25 - t4 * t51 - t5 * t54, t1 * t25 + t10 * t4 + t110 * t18 - t194 * t5 + t2 * t24 + t46 * t67 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t253 * t258 + t200, -t253 * t259 + t199, 0, 0, qJDD(2) * t161 + 0.2e1 * t201, -0.2e1 * qJD(3) * t203 + 0.2e1 * t158 * t169, qJDD(3) * t169 + t172 * t174, t116, qJDD(3) * t172 - t169 * t174, 0, t178 * t169 + t288 * t172, -t288 * t169 + t178 * t172, t251 * t267 + (-g(3) * t170 - t218 * t251) * t164 + t176, -t92 * pkin(2) + g(1) * t100 + g(2) * t99 - t220 + (-t128 * t170 + (t169 * t83 - t172 * t84) * t173) * t250 + t176 * pkin(7), (t120 * t240 + t169 * t87) * t165, (-t274 - t275) * t169 + t193 * t240, (-t87 - t212) * t172 + (t120 * t169 + t165 * t203) * qJD(3), (-t119 * t240 + t169 * t86) * t163, (t86 + t143) * t172 + (t119 * t169 - t163 * t203) * qJD(3), t116, t181 * t163 + (t119 * t229 + pkin(7) * t86 + t27 * t163 + (qJD(2) * t88 + t31) * qJD(3)) * t169 + (-t88 * qJDD(2) - t12 + (-pkin(7) * t119 + t163 * t71) * qJD(3) - t278 * qJD(2) - t182 * t165) * t172, t181 * t165 + (-t120 * t229 + pkin(7) * t87 + t27 * t165 + (-qJD(2) * t89 - t32) * qJD(3)) * t169 + (t89 * qJDD(2) + t13 + (pkin(7) * t120 + t165 * t71) * qJD(3) + t277 * qJD(2) + t182 * t163) * t172, -t86 * t89 - t87 * t88 - t278 * t120 + t277 * t119 + (-t163 * t32 - t165 * t31) * t240 + (-t12 * t165 - t13 * t163 - t182) * t169, t13 * t89 + t12 * t88 - g(1) * (-t197 * t105 - t100) - g(2) * (-t197 * t103 - t99) - t220 + t277 * t32 + t278 * t31 + (-g(3) * t197 - t248 * t71) * t258 + (t169 * t27 + t240 * t71 - t199) * pkin(7), t16 * t97 - t47 * t54, t16 * t96 + t17 * t97 + t47 * t51 - t48 * t54, -t122 * t97 + t152 * t47 + t16 * t172 + t241 * t54, t17 * t96 + t48 * t51, -t122 * t96 + t152 * t48 + t17 * t172 - t241 * t51, -t122 * t172 - t152 * t241, t21 * t122 - t2 * t172 - t194 * t241 + t114 * t51 + t126 * t17 + t18 * t96 + t46 * t48 - g(1) * (-t105 * t264 + t106 * t156) - g(2) * (-t103 * t264 + t104 * t156) - t285 * t152 + (-t51 * t228 - g(3) * (t156 * t170 + t157 * t254)) * t164, t114 * t54 - t126 * t16 - t18 * t97 - t46 * t47 - t22 * t122 + t1 * t172 - t10 * t241 - g(1) * (t105 * t265 + t106 * t157) - g(2) * (t103 * t265 + t104 * t157) + t286 * t152 + (-t54 * t228 - g(3) * (-t156 * t254 + t157 * t170)) * t164, -t1 * t96 - t10 * t48 + t16 * t21 - t182 * t169 - t17 * t22 - t194 * t47 + t2 * t97 - t285 * t54 - t286 * t51, t1 * t22 + t2 * t21 + t18 * t126 + t46 * t114 - g(1) * (-t105 * t192 + t106 * t230 - t100) - g(2) * (-t103 * t192 + t104 * t230 - t99) - t220 - t285 * t194 + t286 * t10 + (-g(3) * pkin(4) * t261 + (-g(3) * t192 - t248 * t46) * t173) * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, t252 * t175, t233, t148, t158, qJDD(3), qJD(3) * t84 - t128 * t246 + t179, -t128 * t244 + (t83 + t115) * qJD(3) + t186 - t231, 0, 0, -t165 * t227 + t275, -t163 * t86 + t165 * t87 - t193 * t244, -t213 + t165 * t263 + (-t120 + t243) * t246, t119 * t226 - t274, -t163 * t263 + t169 * t191 - t211, t148, qJ(4) * t213 - pkin(3) * t86 + t119 * t84 + t184 * t165 + (t163 * t183 - t169 * t31 + t172 * t44) * qJD(2), qJ(4) * t211 - pkin(3) * t87 - t120 * t84 - t184 * t163 + (t165 * t183 + t169 * t32 - t172 * t45) * qJD(2), -t119 * t45 + t120 * t44 + (-qJ(4) * t86 + qJD(4) * t119 + t244 * t31 + t13) * t165 + (qJ(4) * t87 + qJD(4) * t120 + t244 * t32 - t12) * t163 - t186, -t31 * t44 - t32 * t45 - t71 * t84 + (-t163 * t31 + t165 * t32) * qJD(4) + t184 * pkin(3) + (-t12 * t163 + t13 * t165 - t186) * qJ(4), -t16 * t124 - t273 * t54, t123 * t16 - t124 * t17 - t272 * t54 + t273 * t51, t122 * t124 + t152 * t273 - t246 * t54, t17 * t123 + t272 * t51, -t122 * t123 + t152 * t272 + t246 * t51, t152 * t246, t122 * t69 + t123 * t18 - t152 * t279 - t155 * t17 + t157 * t187 + t194 * t246 + t272 * t46 - t51 * t56, t10 * t246 - t122 * t70 + t124 * t18 + t152 * t280 + t155 * t16 - t156 * t187 - t273 * t46 - t54 * t56, -t1 * t123 - t10 * t272 - t124 * t2 + t16 * t69 - t17 * t70 - t194 * t273 - t279 * t54 - t280 * t51 - t186, t1 * t70 + t2 * t69 - t18 * t155 - t46 * t56 - g(1) * (-t155 * t64 + t281 * t65) - g(2) * (-t155 * t62 + t281 * t63) - g(3) * (-t110 * t155 + t111 * t281) - t279 * t194 + t280 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86 - t227, t172 * t191 + t212 + t234, -t119 ^ 2 - t120 ^ 2, -t119 * t32 + t120 * t31 + t177, 0, 0, 0, 0, 0, 0, -t152 * t54 + t17, -t16 + t293, -t296 - t297, t10 * t51 - t194 * t54 + t177 + t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t282, t296 - t297, -t16 - t293, -t282, -t210 + (-qJD(5) - t152) * t54, t122, -t10 * t152 - t46 * t54 - g(1) * (t105 * t157 - t156 * t65) - g(2) * (t103 * t157 - t156 * t63) - g(3) * (-t111 * t156 - t157 * t258) + t2, t46 * t51 + t194 * t152 - g(1) * (-t105 * t156 - t157 * t65) - g(2) * (-t103 * t156 - t157 * t63) - g(3) * (-t111 * t157 + t156 * t258) - t1, 0, 0;];
tau_reg = t3;
