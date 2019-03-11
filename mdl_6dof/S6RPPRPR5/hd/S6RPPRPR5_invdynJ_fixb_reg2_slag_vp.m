% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:49:28
% EndTime: 2019-03-09 01:49:35
% DurationCPUTime: 3.69s
% Computational Cost: add. (4362->472), mult. (8003->578), div. (0->0), fcn. (5001->10), ass. (0->240)
t160 = sin(qJ(6));
t299 = cos(qJ(6));
t155 = sin(pkin(9));
t163 = cos(qJ(4));
t251 = qJD(1) * t163;
t218 = t155 * t251;
t156 = cos(pkin(9));
t244 = t156 * qJD(4);
t94 = -t218 + t244;
t217 = t156 * t251;
t245 = t155 * qJD(4);
t95 = t217 + t245;
t43 = t160 * t95 - t299 * t94;
t314 = t43 ^ 2;
t46 = t160 * t94 + t299 * t95;
t313 = t46 ^ 2;
t161 = sin(qJ(4));
t243 = t161 * qJD(1);
t125 = qJD(6) + t243;
t312 = t43 * t125;
t240 = qJD(1) * qJD(4);
t214 = t163 * t240;
t234 = t161 * qJDD(1);
t173 = t214 + t234;
t158 = -pkin(7) + qJ(2);
t248 = qJD(4) * t163;
t250 = qJD(2) * t161;
t311 = t158 * t248 + t250;
t174 = -t160 * t155 + t299 * t156;
t76 = t174 * t161;
t162 = sin(qJ(1));
t164 = cos(qJ(1));
t307 = g(1) * t164 + g(2) * t162;
t310 = t307 * t156;
t131 = qJ(2) * qJD(1) + qJD(3);
t116 = -pkin(7) * qJD(1) + t131;
t152 = t161 ^ 2;
t153 = t163 ^ 2;
t253 = t152 + t153;
t309 = t253 * t116;
t144 = g(2) * t164;
t145 = g(1) * t162;
t308 = t145 - t144;
t216 = qJD(6) * t299;
t247 = qJD(6) * t160;
t306 = -t155 * t247 + t156 * t216;
t159 = pkin(1) + qJ(3);
t305 = qJD(1) * t159;
t205 = qJD(1) * (t152 - t153);
t233 = t163 * qJDD(1);
t119 = t155 * t233;
t203 = -t156 * qJDD(4) + t119;
t215 = t161 * t240;
t60 = t155 * t215 - t203;
t120 = t156 * t233;
t258 = t155 * qJDD(4) + t120;
t61 = t156 * t215 - t258;
t208 = -t160 * t61 - t299 * t60;
t13 = qJD(6) * t46 + t208;
t222 = t156 * t243;
t223 = t155 * t243;
t282 = -t160 * t223 + t299 * t222 + t306;
t99 = t299 * t155 + t160 * t156;
t304 = -t13 * t99 - t282 * t43;
t83 = t99 * qJD(1);
t85 = t99 * qJD(6);
t283 = t161 * t83 + t85;
t97 = qJDD(6) + t173;
t303 = -t283 * t125 + t174 * t97;
t295 = g(3) * t161;
t302 = t307 * t163 - t295;
t117 = -qJD(2) + t305;
t301 = qJD(4) * (qJD(2) + t117 + t305) + qJDD(4) * t158;
t292 = t161 * pkin(4);
t190 = -t163 * qJ(5) + t292;
t105 = t190 + t159;
t68 = qJD(1) * t105 - qJD(2);
t104 = t161 * t116;
t88 = qJD(4) * qJ(5) + t104;
t29 = -t155 * t88 + t156 * t68;
t18 = pkin(5) * t243 - t95 * pkin(8) + t29;
t30 = t155 * t68 + t156 * t88;
t20 = pkin(8) * t94 + t30;
t177 = t160 * t20 - t299 * t18;
t148 = qJDD(1) * qJ(3);
t154 = qJDD(1) * pkin(1);
t232 = t154 - qJDD(2);
t211 = t148 + t232;
t191 = pkin(4) * t163 + qJ(5) * t161;
t78 = qJD(4) * t191 - t163 * qJD(5) + qJD(3);
t28 = qJD(1) * t78 + qJDD(1) * t190 + t211;
t149 = qJD(1) * qJD(2);
t150 = qJ(2) * qJDD(1);
t212 = qJDD(3) + t149 + t150;
t101 = -pkin(7) * qJDD(1) + t212;
t263 = t163 * t116;
t39 = qJDD(4) * qJ(5) + t161 * t101 + (qJD(5) + t263) * qJD(4);
t14 = -t155 * t39 + t156 * t28;
t8 = pkin(5) * t173 + t61 * pkin(8) + t14;
t15 = t155 * t28 + t156 * t39;
t9 = pkin(8) * t60 + t15;
t1 = -qJD(6) * t177 + t160 * t8 + t299 * t9;
t140 = 0.2e1 * t149;
t300 = t60 * pkin(5);
t298 = pkin(5) * t155;
t297 = pkin(8) * t161;
t296 = pkin(8) * t163;
t294 = g(3) * t163;
t291 = t164 * pkin(7);
t290 = t46 * t43;
t288 = pkin(8) + qJ(5);
t231 = t156 * t297;
t103 = t191 * qJD(1);
t49 = t156 * t103 - t155 * t263;
t27 = (pkin(5) * t163 + t231) * qJD(1) + t49;
t50 = t155 * t103 + t156 * t263;
t38 = pkin(8) * t223 + t50;
t112 = t288 * t155;
t113 = t288 * t156;
t51 = -t299 * t112 - t160 * t113;
t287 = qJD(5) * t174 + qJD(6) * t51 - t160 * t27 - t299 * t38;
t52 = -t160 * t112 + t299 * t113;
t286 = -qJD(5) * t99 - qJD(6) * t52 + t160 * t38 - t299 * t27;
t74 = t99 * t161;
t77 = t174 * t163;
t285 = qJD(4) * t77 - qJD(6) * t74 - t83;
t284 = -t174 * qJD(1) - qJD(6) * t76 - t99 * t248;
t281 = t155 * t60;
t280 = t155 * t95;
t279 = t156 * t61;
t278 = t156 * t94;
t277 = t163 * t60;
t276 = t163 * t61;
t249 = qJD(4) * t161;
t41 = -qJDD(4) * pkin(4) - t163 * t101 + t116 * t249 + qJDD(5);
t275 = t41 * t155;
t274 = t41 * t156;
t273 = t41 * t163;
t272 = t60 * t156;
t271 = t61 * t155;
t267 = t158 * t161;
t57 = t155 * t105 + t156 * t267;
t166 = qJD(1) ^ 2;
t270 = t152 * t166;
t269 = t155 * t166;
t268 = t156 * t166;
t147 = pkin(9) + qJ(6);
t133 = sin(t147);
t266 = t162 * t133;
t134 = cos(t147);
t265 = t162 * t134;
t264 = t162 * t163;
t262 = t163 * t164;
t261 = t164 * t133;
t260 = t164 * t134;
t80 = -qJD(4) * pkin(4) + qJD(5) - t263;
t259 = -qJD(5) + t80;
t257 = t164 * pkin(1) + t162 * qJ(2);
t165 = qJD(4) ^ 2;
t252 = -t165 - t166;
t246 = t117 * qJD(1);
t242 = t163 * qJD(2);
t241 = qJ(5) * qJDD(1);
t239 = qJD(3) * qJD(1);
t237 = qJDD(4) * t161;
t236 = t152 * qJDD(1);
t235 = t159 * qJDD(1);
t36 = t155 * t78 + t311 * t156;
t230 = t43 * t251;
t229 = t95 * t251;
t228 = t163 * t166 * t161;
t227 = t46 * t251;
t226 = t152 * t268 + t173 * t155;
t225 = t164 * qJ(3) + t257;
t224 = pkin(7) + t298;
t220 = t160 * t249;
t213 = qJDD(2) - t308;
t210 = -t155 * t158 + pkin(5);
t209 = -t158 + t298;
t207 = t253 * t101;
t206 = -t94 + t244;
t204 = t101 - t246;
t202 = 0.2e1 * t150 + t140 - t307;
t201 = g(2) * t225;
t200 = t299 * t249;
t198 = t161 * t214;
t197 = -t154 + t213;
t196 = -g(1) * t262 - g(2) * t264 + t295;
t193 = -t246 - t307;
t139 = t164 * qJ(2);
t192 = -t159 * t162 + t139;
t189 = -t14 * t156 - t15 * t155;
t188 = -t14 * t155 + t15 * t156;
t187 = t155 * t30 + t156 * t29;
t186 = t29 * t155 - t30 * t156;
t185 = -t278 + t280;
t130 = t156 * pkin(5) + pkin(4);
t182 = t161 * t130 - t163 * t288;
t181 = -t148 + t197;
t102 = t211 + t239;
t180 = t117 * qJD(3) + t102 * t159;
t12 = -t160 * t60 - t94 * t216 + t95 * t247 + t299 * t61;
t178 = -t12 * t174 - t283 * t46;
t90 = t156 * t105;
t40 = -t156 * t296 + t210 * t161 + t90;
t47 = -t155 * t296 + t57;
t16 = -t160 * t47 + t299 * t40;
t6 = t160 * t18 + t299 * t20;
t17 = t160 * t40 + t299 * t47;
t176 = t307 * t155;
t175 = t282 * t125 + t99 * t97;
t171 = t152 * t269 - t156 * t234;
t169 = -t161 * t307 - t294;
t2 = -qJD(6) * t6 - t160 * t9 + t299 * t8;
t168 = -t158 * t165 + t102 + t145 + t235 + t239;
t167 = -t196 + t41;
t136 = t153 * t166;
t135 = qJDD(4) * t163;
t127 = g(2) * t262;
t93 = t209 * t163;
t92 = 0.2e1 * t198 + t236;
t75 = t99 * t163;
t72 = t161 * t260 - t266;
t71 = -t161 * t261 - t265;
t70 = -t161 * t265 - t261;
t69 = t161 * t266 - t260;
t67 = -pkin(5) * t223 + t104;
t64 = -t209 * t249 - t242;
t63 = t156 * t78;
t56 = -t155 * t267 + t90;
t48 = -t94 * pkin(5) + t80;
t35 = -t311 * t155 + t63;
t34 = -t155 * t200 - t156 * t220 + t306 * t163;
t32 = -t155 * t220 + t156 * t200 + t163 * t85;
t22 = t245 * t297 + t36;
t21 = -t155 * t250 + t63 + (t210 * t163 + t231) * qJD(4);
t19 = t41 - t300;
t4 = -t17 * qJD(6) - t160 * t22 + t299 * t21;
t3 = t16 * qJD(6) + t160 * t21 + t299 * t22;
t5 = [0, 0, 0, 0, 0, qJDD(1), t308, t307, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t154 + t213, t202, t232 * pkin(1) - g(1) * (-t162 * pkin(1) + t139) - g(2) * t257 + (t140 + t150) * qJ(2), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) + t202, -t181 + t235 + 0.2e1 * t239, -g(1) * t192 + t212 * qJ(2) + t131 * qJD(2) + t180 - t201, t153 * qJDD(1) - 0.2e1 * t198, 0.2e1 * qJD(4) * t205 - 0.2e1 * t161 * t233, -t165 * t161 + t135, t92, -t165 * t163 - t237, 0, t301 * t163 + (t168 - t144) * t161, -t301 * t161 + t168 * t163 - t127, t307 + t253 * (-qJDD(1) * t158 - t101 - t149) -g(1) * (t192 - t291) - g(2) * (-t162 * pkin(7) + t225) + t158 * t207 + qJD(2) * t309 + t180 (-t95 * t249 - t276) * t156 (t271 + t272) * t163 + t185 * t249 (-t61 + t120) * t161 + (-t156 * t205 + t163 * t95) * qJD(4) (t94 * t249 - t277) * t155 (t60 - t119) * t161 + (t155 * t205 + t163 * t94) * qJD(4), t92, t176 + (qJD(2) * t94 + t275 + t158 * t60 + (qJD(1) * t56 + t29) * qJD(4)) * t163 + (t35 * qJD(1) + t56 * qJDD(1) + t14 + t308 * t156 + (-t155 * t80 - t158 * t94) * qJD(4)) * t161, t310 + (-qJD(2) * t95 + t274 + t158 * t61 + (-qJD(1) * t57 - t30) * qJD(4)) * t163 + (-t36 * qJD(1) - t57 * qJDD(1) - t15 - t308 * t155 + (-t156 * t80 + t158 * t95) * qJD(4)) * t161, -t35 * t95 + t36 * t94 + t56 * t61 + t57 * t60 + t127 + t187 * t249 + (t189 - t145) * t163, t15 * t57 + t30 * t36 + t14 * t56 + t29 * t35 - t80 * t242 - g(1) * (t139 - t291) - g(2) * (-qJ(5) * t262 + t164 * t292 + t225) + (t80 * t249 - t273) * t158 + (g(2) * pkin(7) + g(1) * t105) * t162, -t12 * t77 - t32 * t46, t12 * t75 - t13 * t77 + t32 * t43 - t34 * t46, -t12 * t161 - t32 * t125 + t46 * t248 + t77 * t97, t13 * t75 + t34 * t43, -t34 * t125 - t13 * t161 - t248 * t43 - t75 * t97, t125 * t248 + t97 * t161, -g(1) * t70 - g(2) * t72 + t4 * t125 + t93 * t13 + t16 * t97 + t2 * t161 - t177 * t248 + t19 * t75 + t48 * t34 + t64 * t43, -g(1) * t69 - g(2) * t71 - t1 * t161 - t93 * t12 - t3 * t125 - t17 * t97 + t19 * t77 - t248 * t6 - t48 * t32 + t64 * t46, -g(1) * t264 - t1 * t75 + t16 * t12 - t17 * t13 - t177 * t32 - t2 * t77 - t3 * t43 - t6 * t34 - t4 * t46 + t127, t1 * t17 + t6 * t3 + t2 * t16 - t177 * t4 + t19 * t93 + t48 * t64 - g(1) * t139 - t201 + (g(1) * t224 - g(2) * t182) * t164 + (-g(1) * (-t182 - t159) + g(2) * t224) * t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t166, -t166 * qJ(2) + t197, 0, 0, 0, 0, 0, 0, 0, -t166, -qJDD(1) (-qJD(3) - t131) * qJD(1) + t181, 0, 0, 0, 0, 0, 0, -0.2e1 * t214 - t234, 0.2e1 * t215 - t233, t136 + t270 (-qJD(3) - t309) * qJD(1) + t181, 0, 0, 0, 0, 0, 0 (-t94 - t244) * t251 + t171, t226 + t229, -t281 - t279 + (-t278 - t280) * t243 (t161 * t186 + t163 * t80) * qJD(1) + t189 - t308, 0, 0, 0, 0, 0, 0, t230 - t303, t175 + t227, t178 - t304, -t1 * t99 - t174 * t2 - t177 * t283 + t48 * t251 - t282 * t6 - t308; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t166, t193 + t212, 0, 0, 0, 0, 0, 0, t252 * t161 + t135, t252 * t163 - t237, -t253 * qJDD(1), t207 + t193, 0, 0, 0, 0, 0, 0, -t155 * t236 + t277 + (-t268 + (-t94 - 0.2e1 * t218) * qJD(4)) * t161, -t156 * t236 + t276 + (t269 + (t95 - 0.2e1 * t217) * qJD(4)) * t161 (qJD(1) * t95 + t161 * t60 + t94 * t248) * t156 + (-qJD(1) * t94 - t161 * t61 + t95 * t248) * t155, -t273 + t188 * t161 - t187 * qJD(1) + (t161 * t80 - t163 * t186) * qJD(4) - t307, 0, 0, 0, 0, 0, 0, t284 * t125 - t163 * t13 + t43 * t249 - t74 * t97, t163 * t12 - t285 * t125 + t46 * t249 - t76 * t97, -t12 * t74 - t13 * t76 - t284 * t46 - t285 * t43, t1 * t76 - t19 * t163 - t177 * t284 - t2 * t74 + t48 * t249 + t285 * t6 - t307; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, t136 - t270, t233, -t228, -t234, qJDD(4), t163 * t204 + t196, t294 + (t307 - t204) * t161, 0, 0, t95 * t222 - t271, -t185 * t243 - t279 + t281, t226 - t229, -t94 * t223 + t272, t206 * t251 - t171, -t228, pkin(4) * t60 - t274 + (-t310 + (-qJ(5) * t245 - t29) * qJD(1)) * t163 + (-t155 * t241 + g(3) * t156 + t116 * t94 + (t259 * t155 - t49) * qJD(1)) * t161, pkin(4) * t61 + t275 + (t176 + (-qJ(5) * t244 + t30) * qJD(1)) * t163 + (-t156 * t241 - g(3) * t155 - t116 * t95 + (t259 * t156 + t50) * qJD(1)) * t161, t49 * t95 - t50 * t94 + (qJ(5) * t60 + qJD(5) * t94 - t29 * t243 + t15) * t156 + (-qJ(5) * t61 + qJD(5) * t95 - t30 * t243 - t14) * t155 + t169, -t80 * t104 - t29 * t49 - t30 * t50 - t186 * qJD(5) + (-t41 - t302) * pkin(4) + (t169 + t188) * qJ(5), -t12 * t99 + t282 * t46, t178 + t304, t175 - t227, -t13 * t174 + t283 * t43, t230 + t303, -t125 * t251, t286 * t125 - t130 * t13 - t134 * t302 - t174 * t19 + t177 * t251 + t283 * t48 - t67 * t43 + t51 * t97, t130 * t12 - t287 * t125 + t302 * t133 + t19 * t99 + t6 * t251 + t282 * t48 - t67 * t46 - t52 * t97, t1 * t174 + t51 * t12 - t52 * t13 + t177 * t282 - t2 * t99 - t283 * t6 - t286 * t46 - t287 * t43 + t169, g(3) * t182 + t1 * t52 - t19 * t130 + t2 * t51 - t286 * t177 + t287 * t6 - t48 * t67 - t307 * (t130 * t163 + t161 * t288); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t95 - t245) * t243 + t203, -t206 * t243 + t258, -t94 ^ 2 - t95 ^ 2, t29 * t95 - t30 * t94 + t167, 0, 0, 0, 0, 0, 0, t46 * t125 + t13, -t12 - t312, -t313 - t314, -t177 * t46 + t43 * t6 + t167 - t300; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290, t313 - t314, -t12 + t312, -t290, -t208 + (-qJD(6) + t125) * t46, t97, -g(1) * t71 + g(2) * t69 + t6 * t125 + t133 * t294 - t48 * t46 + t2, g(1) * t72 - g(2) * t70 - t125 * t177 + t134 * t294 + t48 * t43 - t1, 0, 0;];
tau_reg  = t5;
