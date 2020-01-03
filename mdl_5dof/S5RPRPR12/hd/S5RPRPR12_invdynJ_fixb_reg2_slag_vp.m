% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR12_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:30:21
% EndTime: 2019-12-31 18:30:29
% DurationCPUTime: 4.95s
% Computational Cost: add. (5978->412), mult. (14482->527), div. (0->0), fcn. (11184->14), ass. (0->202)
t171 = sin(pkin(8));
t173 = cos(pkin(8));
t177 = sin(qJ(3));
t276 = cos(qJ(3));
t137 = t171 * t276 + t177 * t173;
t127 = t137 * qJD(1);
t170 = sin(pkin(9));
t172 = cos(pkin(9));
t104 = qJD(3) * t170 + t127 * t172;
t176 = sin(qJ(5));
t275 = cos(qJ(5));
t241 = t170 * t127;
t292 = qJD(3) * t172 - t241;
t191 = t275 * t292;
t53 = -t104 * t176 + t191;
t295 = t53 ^ 2;
t227 = t276 * t173;
t152 = qJD(1) * t227;
t240 = t171 * t177;
t225 = qJD(1) * t240;
t125 = -t152 + t225;
t122 = qJD(5) + t125;
t294 = t122 * t53;
t52 = t275 * t104 + t176 * t292;
t293 = t52 ^ 2;
t169 = pkin(8) + qJ(3);
t162 = sin(t169);
t164 = cos(t169);
t178 = sin(qJ(1));
t179 = cos(qJ(1));
t210 = g(1) * t179 + g(2) * t178;
t188 = -g(3) * t164 + t162 * t210;
t232 = qJD(1) * qJD(2);
t267 = pkin(6) + qJ(2);
t279 = qJDD(1) * t267 + t232;
t116 = t279 * t171;
t117 = t279 * t173;
t143 = t267 * t171;
t138 = qJD(1) * t143;
t145 = t267 * t173;
t139 = qJD(1) * t145;
t222 = qJD(3) * t276;
t234 = qJD(3) * t177;
t41 = -t116 * t276 - t177 * t117 + t138 * t234 - t139 * t222;
t39 = -qJDD(3) * pkin(3) + qJDD(4) - t41;
t182 = t39 - t188;
t226 = t275 * t172;
t192 = -t176 * t170 + t226;
t233 = qJD(5) * t176;
t284 = qJD(5) * t226 - t170 * t233;
t254 = -t125 * t192 - t284;
t136 = t170 * t275 + t176 * t172;
t130 = t136 * qJD(5);
t253 = t125 * t136 + t130;
t252 = qJDD(1) * pkin(1);
t289 = -g(1) * t178 + g(2) * t179;
t199 = -qJDD(2) + t252 - t289;
t194 = t289 * t164;
t287 = t292 * t127;
t286 = -t143 * t276 - t145 * t177;
t283 = qJ(2) * qJDD(1);
t224 = t171 * t234;
t218 = qJDD(1) * t276;
t230 = t173 * qJDD(1);
t228 = qJD(3) * t152 + t171 * t218 + t177 * t230;
t91 = qJD(1) * t224 - t228;
t74 = -qJDD(3) * t172 - t170 * t91;
t75 = qJDD(3) * t170 - t172 * t91;
t17 = -qJD(5) * t191 + t104 * t233 + t176 * t74 - t275 * t75;
t282 = -t17 * t192 - t253 * t52;
t132 = t137 * qJD(3);
t231 = t171 * qJDD(1);
t207 = -t173 * t218 + t177 * t231;
t92 = qJD(1) * t132 + t207;
t88 = qJDD(5) + t92;
t281 = -t122 * t254 + t136 * t88;
t123 = t125 ^ 2;
t258 = t170 * t92;
t280 = -t123 * t172 - t258;
t269 = g(3) * t162;
t186 = -t164 * t210 - t269;
t159 = pkin(2) * t173 + pkin(1);
t141 = -qJD(1) * t159 + qJD(2);
t66 = pkin(3) * t125 - qJ(4) * t127 + t141;
t95 = -t177 * t138 + t139 * t276;
t86 = qJD(3) * qJ(4) + t95;
t33 = -t170 * t86 + t172 * t66;
t20 = pkin(4) * t125 - pkin(7) * t104 + t33;
t34 = t170 * t66 + t172 * t86;
t23 = pkin(7) * t292 + t34;
t197 = t176 * t23 - t20 * t275;
t140 = -qJDD(1) * t159 + qJDD(2);
t30 = pkin(3) * t92 + qJ(4) * t91 - qJD(4) * t127 + t140;
t124 = t177 * t139;
t229 = -t177 * t116 + t117 * t276 - t138 * t222;
t36 = qJDD(3) * qJ(4) + (qJD(4) - t124) * qJD(3) + t229;
t14 = -t170 * t36 + t172 * t30;
t6 = pkin(4) * t92 - pkin(7) * t75 + t14;
t15 = t170 * t30 + t172 * t36;
t9 = -pkin(7) * t74 + t15;
t1 = -qJD(5) * t197 + t176 * t6 + t275 * t9;
t278 = t127 ^ 2;
t277 = pkin(4) * t74;
t274 = pkin(7) * t172;
t146 = t179 * t159;
t271 = g(2) * t146;
t268 = t52 * t53;
t266 = pkin(7) + qJ(4);
t131 = -t173 * t222 + t224;
t57 = pkin(3) * t132 + qJ(4) * t131 - qJD(4) * t137;
t193 = t227 - t240;
t69 = qJD(2) * t193 + qJD(3) * t286;
t28 = t170 * t57 + t172 * t69;
t247 = t125 * t172;
t87 = pkin(3) * t127 + qJ(4) * t125;
t94 = -t138 * t276 - t124;
t42 = -t170 * t94 + t172 * t87;
t24 = pkin(4) * t127 + pkin(7) * t247 + t42;
t248 = t125 * t170;
t43 = t170 * t87 + t172 * t94;
t31 = pkin(7) * t248 + t43;
t142 = t266 * t170;
t144 = t266 * t172;
t99 = -t142 * t275 - t176 * t144;
t265 = qJD(4) * t192 + qJD(5) * t99 - t176 * t24 - t275 * t31;
t101 = -t176 * t142 + t144 * t275;
t264 = -qJD(4) * t136 - qJD(5) * t101 + t176 * t31 - t24 * t275;
t83 = t172 * t92;
t263 = -t123 * t170 + t83;
t102 = -t177 * t143 + t145 * t276;
t89 = -pkin(3) * t193 - qJ(4) * t137 - t159;
t46 = t172 * t102 + t170 * t89;
t262 = t127 * t53;
t261 = t127 * t52;
t257 = t74 * t172;
t256 = t75 * t170;
t255 = t75 * t172;
t251 = t104 * t127;
t250 = t104 * t170;
t249 = t125 * t127;
t246 = t131 * t170;
t245 = t137 * t170;
t243 = t164 * t178;
t242 = t164 * t179;
t238 = t95 * qJD(3);
t84 = -qJD(3) * pkin(3) + qJD(4) - t94;
t237 = -qJD(4) + t84;
t166 = t171 ^ 2;
t167 = t173 ^ 2;
t236 = t166 + t167;
t221 = pkin(4) * t170 + t267;
t27 = -t170 * t69 + t172 * t57;
t219 = t176 * t75 + t275 * t74;
t45 = -t102 * t170 + t172 * t89;
t217 = t236 * qJD(1) ^ 2;
t18 = qJD(5) * t52 + t219;
t215 = -t136 * t18 - t254 * t53;
t214 = 0.2e1 * t236;
t213 = -t122 * t253 + t192 * t88;
t211 = t289 * t162;
t208 = pkin(3) * t164 + qJ(4) * t162;
t206 = -t14 * t172 - t15 * t170;
t205 = -t33 * t170 + t34 * t172;
t204 = t131 * t125 - t137 * t92;
t203 = t125 * t132 - t193 * t92;
t158 = t172 * pkin(4) + pkin(3);
t202 = t158 * t164 + t162 * t266;
t200 = t172 * t292;
t29 = -pkin(4) * t193 - t137 * t274 + t45;
t35 = -pkin(7) * t245 + t46;
t12 = -t176 * t35 + t275 * t29;
t8 = t176 * t20 + t23 * t275;
t13 = t176 * t29 + t275 * t35;
t190 = t199 + t252;
t189 = t200 - t250;
t185 = -t84 * t131 + t39 * t137 - t210;
t2 = -qJD(5) * t8 - t176 * t9 + t275 * t6;
t183 = t214 * t232 - t210;
t70 = qJD(2) * t137 + qJD(3) * t102;
t168 = pkin(9) + qJ(5);
t163 = cos(t168);
t161 = sin(t168);
t114 = t161 * t178 + t163 * t242;
t113 = -t161 * t242 + t163 * t178;
t112 = t161 * t179 - t163 * t243;
t111 = t161 * t243 + t163 * t179;
t80 = t192 * t137;
t79 = t136 * t137;
t71 = pkin(4) * t245 - t286;
t62 = t170 * t74;
t56 = -pkin(4) * t248 + t95;
t48 = -pkin(4) * t292 + t84;
t47 = -pkin(4) * t246 + t70;
t40 = -t139 * t234 + t229;
t38 = -t136 * t131 + t137 * t284;
t37 = t137 * t130 + t131 * t192;
t22 = pkin(7) * t246 + t28;
t21 = t39 + t277;
t19 = pkin(4) * t132 + t131 * t274 + t27;
t4 = -qJD(5) * t13 - t176 * t22 + t19 * t275;
t3 = qJD(5) * t12 + t176 * t19 + t22 * t275;
t5 = [0, 0, 0, 0, 0, qJDD(1), -t289, t210, 0, 0, t166 * qJDD(1), 0.2e1 * t171 * t230, 0, t167 * qJDD(1), 0, 0, t190 * t173, -t190 * t171, t214 * t283 + t183, t199 * pkin(1) + (t236 * t283 + t183) * qJ(2), -t127 * t131 - t137 * t91, -t127 * t132 - t193 * t91 + t204, -t131 * qJD(3) + qJDD(3) * t137, t203, -qJD(3) * t132 + qJDD(3) * t193, 0, -t70 * qJD(3) + qJDD(3) * t286 + t132 * t141 - t140 * t193 - t159 * t92 - t194, -qJD(3) * t69 - qJDD(3) * t102 - t131 * t141 + t137 * t140 + t159 * t91 + t211, -t102 * t92 - t125 * t69 + t127 * t70 + t131 * t94 - t132 * t95 - t137 * t41 + t193 * t40 + t286 * t91 - t210, t40 * t102 + t95 * t69 + t41 * t286 - t94 * t70 - t140 * t159 - g(1) * (-t159 * t178 + t179 * t267) - g(2) * (t178 * t267 + t146), (-t104 * t131 + t75 * t137) * t172, (-t256 - t257) * t137 - t189 * t131, t104 * t132 - t172 * t204 - t193 * t75, (t131 * t292 + t74 * t137) * t170, t132 * t292 + t170 * t204 + t193 * t74, t203, t27 * t125 + t33 * t132 - t14 * t193 + t170 * t185 - t172 * t194 - t286 * t74 - t292 * t70 + t45 * t92, t70 * t104 - t28 * t125 - t34 * t132 + t15 * t193 + t170 * t194 + t172 * t185 - t286 * t75 - t46 * t92, t28 * t292 - t46 * t74 - t27 * t104 - t45 * t75 + t206 * t137 + (t34 * t170 + t33 * t172) * t131 - t211, -t271 - t39 * t286 + t14 * t45 + t15 * t46 + t33 * t27 + t34 * t28 + t84 * t70 + (-g(1) * t267 - g(2) * t208) * t179 + (-g(1) * (-t159 - t208) - g(2) * t267) * t178, -t17 * t80 - t37 * t52, t17 * t79 - t18 * t80 - t37 * t53 - t38 * t52, -t122 * t37 + t132 * t52 + t17 * t193 + t80 * t88, t18 * t79 - t38 * t53, -t122 * t38 + t132 * t53 + t18 * t193 - t79 * t88, t122 * t132 - t193 * t88, -g(1) * t112 - g(2) * t114 + t12 * t88 + t122 * t4 - t132 * t197 + t18 * t71 - t193 * t2 + t21 * t79 + t38 * t48 - t47 * t53, -g(1) * t111 - g(2) * t113 + t1 * t193 - t122 * t3 - t13 * t88 - t132 * t8 - t17 * t71 + t21 * t80 - t37 * t48 + t47 * t52, -t1 * t79 + t12 * t17 - t13 * t18 - t197 * t37 - t2 * t80 + t3 * t53 - t38 * t8 - t4 * t52 - t211, -t271 + t1 * t13 + t2 * t12 + t21 * t71 + t8 * t3 - t197 * t4 + t48 * t47 + (-g(1) * t221 - g(2) * t202) * t179 + (-g(1) * (-t159 - t202) - g(2) * t221) * t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t230, t231, -t217, -qJ(2) * t217 - t199, 0, 0, 0, 0, 0, 0, 0.2e1 * t127 * qJD(3) + t207, (-t125 - t225) * qJD(3) + t228, -t123 - t278, t125 * t95 + t127 * t94 + t140 + t289, 0, 0, 0, 0, 0, 0, t263 + t287, -t251 + t280, -t255 - t62 + (t200 + t250) * t125, t125 * t205 - t127 * t84 - t206 + t289, 0, 0, 0, 0, 0, 0, t213 + t262, -t261 - t281, t215 - t282, t1 * t136 - t127 * t48 + t192 * t2 + t197 * t253 - t254 * t8 + t289; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t249, -t123 + t278, (t125 - t225) * qJD(3) + t228, -t249, -t207, qJDD(3), -t127 * t141 + t188 + t238 + t41, t125 * t141 + (t94 + t124) * qJD(3) - t229 - t186, 0, 0, t104 * t247 + t256, t125 * t189 + t255 - t62, -t251 - t280, -t248 * t292 - t257, t263 - t287, -t249, -qJ(4) * t258 - pkin(3) * t74 - t95 * t241 - t33 * t127 + (t170 * t237 - t42) * t125 + (-t182 + t238) * t172, -qJ(4) * t83 - pkin(3) * t75 - t104 * t95 + t127 * t34 + (t172 * t237 + t43) * t125 + t182 * t170, t42 * t104 + t43 * t241 + (-qJ(4) * t74 - qJD(4) * t241 - t33 * t125 + t15 + (qJD(4) * t172 - t43) * qJD(3)) * t172 + (qJ(4) * t75 + qJD(4) * t104 - t34 * t125 - t14) * t170 + t186, -t33 * t42 - t34 * t43 - t84 * t95 + t205 * qJD(4) - t182 * pkin(3) + (-t14 * t170 + t15 * t172 + t186) * qJ(4), -t136 * t17 - t254 * t52, t215 + t282, -t261 + t281, -t18 * t192 - t253 * t53, t213 - t262, -t122 * t127, t122 * t264 + t127 * t197 - t158 * t18 + t163 * t188 - t192 * t21 + t253 * t48 + t53 * t56 + t88 * t99, -t101 * t88 - t122 * t265 + t127 * t8 + t136 * t21 + t158 * t17 - t161 * t188 - t254 * t48 - t52 * t56, t1 * t192 - t101 * t18 - t136 * t2 + t17 * t99 - t197 * t254 - t253 * t8 - t264 * t52 + t265 * t53 + t186, -g(3) * t202 + t1 * t101 - t21 * t158 + t2 * t99 - t264 * t197 + t265 * t8 - t48 * t56 + t210 * (t158 * t162 - t164 * t266); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104 * t125 + t74, t125 * t292 + t75, -t104 ^ 2 - t292 ^ 2, t104 * t33 - t292 * t34 + t182, 0, 0, 0, 0, 0, 0, t52 * t122 + t18, -t17 + t294, -t293 - t295, -t197 * t52 - t53 * t8 + t182 + t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t268, t293 - t295, -t17 - t294, t268, -t219 + (-qJD(5) + t122) * t52, t88, -g(1) * t113 + g(2) * t111 + t8 * t122 + t161 * t269 - t48 * t52 + t2, g(1) * t114 - g(2) * t112 - t122 * t197 + t163 * t269 - t48 * t53 - t1, 0, 0;];
tau_reg = t5;
