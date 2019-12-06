% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:21
% EndTime: 2019-12-05 18:12:32
% DurationCPUTime: 4.43s
% Computational Cost: add. (8188->379), mult. (20755->482), div. (0->0), fcn. (16432->16), ass. (0->195)
t201 = cos(pkin(9));
t209 = cos(qJ(3));
t254 = t209 * t201;
t240 = qJD(1) * t254;
t200 = sin(pkin(9));
t205 = sin(qJ(3));
t255 = t205 * t200;
t241 = qJD(1) * t255;
t139 = -t240 + t241;
t151 = t209 * t200 + t205 * t201;
t141 = t151 * qJD(1);
t204 = sin(qJ(4));
t208 = cos(qJ(4));
t101 = -t208 * t139 - t204 * t141;
t203 = sin(qJ(5));
t207 = cos(qJ(5));
t244 = t201 * qJDD(1);
t245 = t200 * qJDD(1);
t242 = qJD(3) * t240 + t205 * t244 + t209 * t245;
t110 = qJD(3) * t241 - t242;
t144 = t151 * qJD(3);
t228 = t205 * t245 - t209 * t244;
t111 = qJD(1) * t144 + t228;
t226 = t204 * t139 - t208 * t141;
t218 = qJD(4) * t226 + t204 * t110 - t208 * t111;
t249 = qJD(4) * t208;
t250 = qJD(4) * t204;
t223 = -t208 * t110 - t204 * t111 - t139 * t249 - t141 * t250;
t247 = qJD(5) * t207;
t248 = qJD(5) * t203;
t10 = -t101 * t247 - t203 * t218 - t207 * t223 - t226 * t248;
t199 = qJD(3) + qJD(4);
t190 = qJD(5) + t199;
t57 = t207 * t101 + t203 * t226;
t268 = t190 * t57;
t294 = -t10 - t268;
t274 = t57 ^ 2;
t290 = t203 * t101 - t207 * t226;
t275 = t290 ^ 2;
t295 = -t274 + t275;
t273 = t290 * t57;
t297 = t101 * pkin(8);
t272 = pkin(6) + qJ(2);
t165 = t272 * t200;
t152 = qJD(1) * t165;
t166 = t272 * t201;
t153 = qJD(1) * t166;
t115 = -t205 * t152 + t209 * t153;
t84 = -t139 * pkin(7) + t115;
t81 = t208 * t84;
t256 = t205 * t153;
t114 = -t209 * t152 - t256;
t83 = -t141 * pkin(7) + t114;
t82 = qJD(3) * pkin(3) + t83;
t46 = t204 * t82 + t81;
t30 = t46 + t297;
t267 = t203 * t30;
t298 = pkin(8) * t226;
t79 = t204 * t84;
t45 = t208 * t82 - t79;
t29 = t45 + t298;
t27 = t199 * pkin(4) + t29;
t12 = t207 * t27 - t267;
t266 = t207 * t30;
t13 = t203 * t27 + t266;
t304 = t12 * t57 + t13 * t290;
t194 = qJDD(3) + qJDD(4);
t246 = qJD(1) * qJD(2);
t284 = t272 * qJDD(1) + t246;
t125 = t284 * t200;
t126 = t284 * t201;
t232 = -t209 * t125 - t205 * t126;
t60 = -t115 * qJD(3) + t232;
t38 = qJDD(3) * pkin(3) + t110 * pkin(7) + t60;
t251 = qJD(3) * t209;
t243 = t205 * t125 - t209 * t126 + t152 * t251;
t252 = qJD(3) * t205;
t59 = -t153 * t252 - t243;
t42 = -t111 * pkin(7) + t59;
t9 = -qJD(4) * t46 - t204 * t42 + t208 * t38;
t6 = t194 * pkin(4) - pkin(8) * t223 + t9;
t235 = -t204 * t38 - t208 * t42 - t82 * t249 + t84 * t250;
t7 = pkin(8) * t218 - t235;
t1 = (qJD(5) * t27 + t7) * t207 + t203 * t6 - t30 * t248;
t198 = pkin(9) + qJ(3);
t191 = qJ(4) + t198;
t184 = qJ(5) + t191;
t177 = sin(t184);
t178 = cos(t184);
t206 = sin(qJ(1));
t210 = cos(qJ(1));
t230 = g(1) * t210 + g(2) * t206;
t183 = t201 * pkin(2) + pkin(1);
t160 = -t183 * qJD(1) + qJD(2);
t116 = t139 * pkin(3) + t160;
t68 = -pkin(4) * t101 + t116;
t293 = g(3) * t177 + t230 * t178 - t68 * t57 - t1;
t11 = qJD(5) * t290 + t203 * t223 - t207 * t218;
t265 = t290 * t190;
t285 = -t11 + t265;
t264 = t226 * t199;
t302 = t218 - t264;
t260 = t101 * t199;
t301 = t223 - t260;
t270 = t226 ^ 2;
t271 = t101 ^ 2;
t300 = t270 - t271;
t181 = sin(t191);
t182 = cos(t191);
t299 = -g(3) * t182 + t230 * t181;
t2 = -qJD(5) * t13 - t203 * t7 + t207 * t6;
t286 = -g(3) * t178 + t230 * t177 - t290 * t68 + t2;
t269 = t101 * t226;
t261 = qJDD(1) * pkin(1);
t288 = g(1) * t206 - g(2) * t210;
t225 = -qJDD(2) + t261 + t288;
t292 = g(3) * t181 - t116 * t101 + t230 * t182 + t235;
t291 = t116 * t226 + t299 + t9;
t118 = -t205 * t165 + t209 * t166;
t287 = qJ(2) * qJDD(1);
t283 = t141 ^ 2;
t282 = t218 * pkin(4);
t277 = t111 * pkin(3);
t276 = t144 * pkin(3);
t48 = t208 * t83 - t79;
t117 = -t209 * t165 - t205 * t166;
t96 = -t151 * pkin(7) + t117;
t150 = -t254 + t255;
t97 = -t150 * pkin(7) + t118;
t51 = t204 * t96 + t208 * t97;
t185 = t208 * pkin(3) + pkin(4);
t258 = t203 * t204;
t47 = -t204 * t83 - t81;
t31 = t47 - t297;
t32 = t48 + t298;
t263 = -t203 * t31 - t207 * t32 + t185 * t247 + (-t204 * t248 + (t207 * t208 - t258) * qJD(4)) * pkin(3);
t257 = t204 * t207;
t262 = t203 * t32 - t207 * t31 - t185 * t248 + (-t204 * t247 + (-t203 * t208 - t257) * qJD(4)) * pkin(3);
t259 = t141 * t139;
t196 = t200 ^ 2;
t197 = t201 ^ 2;
t253 = t196 + t197;
t195 = -pkin(7) - t272;
t50 = -t204 * t97 + t208 * t96;
t234 = t253 * qJD(1) ^ 2;
t189 = cos(t198);
t179 = pkin(3) * t189;
t154 = t179 + t183;
t231 = 0.2e1 * t253;
t113 = -t204 * t150 + t208 * t151;
t33 = -t113 * pkin(8) + t50;
t112 = t208 * t150 + t204 * t151;
t34 = -t112 * pkin(8) + t51;
t20 = -t203 * t34 + t207 * t33;
t21 = t203 * t33 + t207 * t34;
t67 = -t203 * t112 + t207 * t113;
t122 = t150 * pkin(3) - t183;
t87 = -t165 * t251 + qJD(2) * t254 + (-qJD(2) * t200 - qJD(3) * t166) * t205;
t72 = -t144 * pkin(7) + t87;
t143 = t200 * t252 - t201 * t251;
t88 = -t151 * qJD(2) - t118 * qJD(3);
t73 = t143 * pkin(7) + t88;
t24 = t204 * t73 + t208 * t72 + t96 * t249 - t97 * t250;
t155 = -t183 * qJDD(1) + qJDD(2);
t222 = t225 + t261;
t188 = sin(t198);
t221 = -g(3) * t189 + t188 * t230;
t85 = t155 + t277;
t220 = t155 - t288;
t25 = -t51 * qJD(4) - t204 * t72 + t208 * t73;
t217 = t231 * t246 - t230;
t215 = t220 + t277;
t192 = -pkin(8) + t195;
t187 = qJDD(5) + t194;
t176 = pkin(4) * t182;
t136 = pkin(3) * t257 + t203 * t185;
t135 = -pkin(3) * t258 + t207 * t185;
t134 = t139 ^ 2;
t124 = t154 + t176;
t77 = t112 * pkin(4) + t122;
t74 = t141 * pkin(3) - pkin(4) * t226;
t66 = t207 * t112 + t203 * t113;
t65 = qJD(4) * t113 - t204 * t143 + t208 * t144;
t64 = t208 * t143 + t204 * t144 + t150 * t249 + t151 * t250;
t49 = t65 * pkin(4) + t276;
t26 = t85 - t282;
t23 = qJD(5) * t67 - t203 * t64 + t207 * t65;
t22 = t112 * t247 + t113 * t248 + t203 * t65 + t207 * t64;
t19 = t64 * pkin(8) + t25;
t18 = -t65 * pkin(8) + t24;
t15 = t207 * t29 - t267;
t14 = -t203 * t29 - t266;
t4 = -qJD(5) * t21 - t203 * t18 + t207 * t19;
t3 = qJD(5) * t20 + t207 * t18 + t203 * t19;
t5 = [0, 0, 0, 0, 0, qJDD(1), t288, t230, 0, 0, t196 * qJDD(1), 0.2e1 * t200 * t244, 0, t197 * qJDD(1), 0, 0, t222 * t201, -t222 * t200, t231 * t287 + t217, pkin(1) * t225 + (t253 * t287 + t217) * qJ(2), -t110 * t151 - t141 * t143, t110 * t150 - t151 * t111 + t143 * t139 - t141 * t144, -t143 * qJD(3) + t151 * qJDD(3), t111 * t150 + t139 * t144, -t144 * qJD(3) - t150 * qJDD(3), 0, t88 * qJD(3) + t117 * qJDD(3) - t183 * t111 + t160 * t144 + t155 * t150 + t189 * t288, -t87 * qJD(3) - t118 * qJDD(3) + t183 * t110 - t160 * t143 + t155 * t151 - t188 * t288, t117 * t110 - t118 * t111 + t114 * t143 - t115 * t144 - t87 * t139 - t88 * t141 - t59 * t150 - t60 * t151 - t230, t59 * t118 + t115 * t87 + t60 * t117 + t114 * t88 - t155 * t183 - g(1) * (-t206 * t183 + t210 * t272) - g(2) * (t210 * t183 + t206 * t272), t113 * t223 + t226 * t64, -t101 * t64 - t112 * t223 + t113 * t218 + t226 * t65, t113 * t194 - t64 * t199, -t101 * t65 - t112 * t218, -t112 * t194 - t65 * t199, 0, -t101 * t276 + t85 * t112 + t116 * t65 - t122 * t218 + t182 * t288 + t50 * t194 + t25 * t199, t85 * t113 - t116 * t64 + t122 * t223 - t181 * t288 - t51 * t194 - t24 * t199 - t226 * t276, t101 * t24 + t112 * t235 - t9 * t113 + t218 * t51 - t223 * t50 + t226 * t25 + t45 * t64 - t46 * t65 - t230, -t235 * t51 + t46 * t24 + t9 * t50 + t45 * t25 + t85 * t122 + t116 * t276 - g(1) * (-t206 * t154 - t210 * t195) - g(2) * (t210 * t154 - t206 * t195), -t10 * t67 - t22 * t290, t10 * t66 - t67 * t11 - t22 * t57 - t23 * t290, t67 * t187 - t22 * t190, t11 * t66 - t23 * t57, -t66 * t187 - t23 * t190, 0, t77 * t11 + t178 * t288 + t20 * t187 + t4 * t190 + t68 * t23 + t26 * t66 - t49 * t57, -t77 * t10 - t177 * t288 - t21 * t187 - t3 * t190 - t68 * t22 + t26 * t67 + t290 * t49, -t1 * t66 + t20 * t10 - t21 * t11 + t12 * t22 - t13 * t23 - t2 * t67 - t290 * t4 + t3 * t57 - t230, t1 * t21 + t13 * t3 + t2 * t20 + t12 * t4 + t26 * t77 + t68 * t49 - g(1) * (-t206 * t124 - t210 * t192) - g(2) * (t210 * t124 - t206 * t192); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244, t245, -t234, -qJ(2) * t234 - t225, 0, 0, 0, 0, 0, 0, 0.2e1 * t141 * qJD(3) + t228, (-t139 - t241) * qJD(3) + t242, -t134 - t283, t114 * t141 + t115 * t139 + t220, 0, 0, 0, 0, 0, 0, -t218 - t264, t223 + t260, -t270 - t271, -t46 * t101 - t226 * t45 + t215, 0, 0, 0, 0, 0, 0, t11 + t265, -t10 + t268, -t274 - t275, t12 * t290 - t13 * t57 + t215 - t282; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t259, -t134 + t283, (t139 - t241) * qJD(3) + t242, -t259, -t228, qJDD(3), -t160 * t141 + t221 + t232, g(3) * t188 + t160 * t139 + t230 * t189 + (t114 + t256) * qJD(3) + t243, 0, 0, t269, t300, t301, -t269, t302, t194, -t47 * t199 + (t101 * t141 + t194 * t208 - t199 * t250) * pkin(3) + t291, t48 * t199 + (t141 * t226 - t194 * t204 - t199 * t249) * pkin(3) + t292, t45 * t101 - t47 * t226 - t46 * t226 - t48 * t101 + (t204 * t218 - t208 * t223 + (t101 * t208 - t204 * t226) * qJD(4)) * pkin(3), -t45 * t47 - t46 * t48 + (-t116 * t141 - t204 * t235 + t208 * t9 + (-t204 * t45 + t208 * t46) * qJD(4) + t221) * pkin(3), -t273, t295, t294, t273, t285, t187, t135 * t187 + t190 * t262 + t57 * t74 + t286, -t136 * t187 - t190 * t263 - t290 * t74 + t293, t135 * t10 - t136 * t11 - t262 * t290 + t263 * t57 + t304, t1 * t136 + t2 * t135 - t68 * t74 - g(3) * (t176 + t179) - t230 * (-pkin(3) * t188 - pkin(4) * t181) + t263 * t13 + t262 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t269, t300, t301, -t269, t302, t194, t46 * t199 + t291, t199 * t45 + t292, 0, 0, -t273, t295, t294, t273, t285, t187, -t14 * t190 + (t187 * t207 - t190 * t248 - t226 * t57) * pkin(4) + t286, t15 * t190 + (-t187 * t203 - t190 * t247 + t226 * t290) * pkin(4) + t293, t14 * t290 - t15 * t57 + (t10 * t207 - t11 * t203 + (t203 * t290 + t207 * t57) * qJD(5)) * pkin(4) + t304, -t12 * t14 - t13 * t15 + (t1 * t203 + t226 * t68 + t2 * t207 + (-t12 * t203 + t13 * t207) * qJD(5) + t299) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t273, t295, t294, t273, t285, t187, t13 * t190 + t286, t12 * t190 + t293, 0, 0;];
tau_reg = t5;
