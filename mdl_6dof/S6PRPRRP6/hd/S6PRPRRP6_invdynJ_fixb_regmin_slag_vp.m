% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% tau_reg [6x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:20:59
% EndTime: 2019-03-08 20:21:07
% DurationCPUTime: 2.90s
% Computational Cost: add. (3031->409), mult. (6221->535), div. (0->0), fcn. (4617->10), ass. (0->208)
t131 = sin(qJ(4));
t134 = cos(qJ(4));
t203 = t134 * qJDD(2);
t210 = qJD(2) * qJD(4);
t286 = -t131 * t210 + t203;
t118 = qJD(2) * t131 + qJD(5);
t130 = sin(qJ(5));
t181 = qJD(5) * t131 + qJD(2);
t133 = cos(qJ(5));
t211 = t133 * qJD(4);
t195 = t134 * t211;
t216 = qJD(4) * t130;
t217 = qJD(2) * t134;
t95 = t133 * t217 + t216;
t244 = qJD(4) * t95;
t213 = qJD(5) * t130;
t194 = t134 * t213;
t35 = qJD(2) * (t131 * t211 + t194) - qJD(5) * t211 - t130 * qJDD(4) - t133 * t203;
t251 = t134 * t35;
t190 = t134 * t210;
t204 = t131 * qJDD(2);
t90 = qJDD(5) + t190 + t204;
t255 = t133 * t90;
t285 = t118 * (t130 * t181 - t195) + (t244 - t255) * t131 + t251;
t136 = -pkin(2) - pkin(8);
t231 = t131 * t133;
t177 = pkin(4) * t131 - pkin(9) * t134;
t99 = qJ(3) + t177;
t247 = t130 * t99 + t136 * t231;
t127 = sin(pkin(6));
t222 = qJD(1) * t127;
t135 = cos(qJ(2));
t132 = sin(qJ(2));
t233 = t130 * t132;
t281 = t131 * t233 - t133 * t135;
t178 = pkin(4) * t134 + pkin(9) * t131;
t91 = qJD(4) * t178 + qJD(3);
t284 = -t133 * t91 - t281 * t222;
t212 = qJD(5) * t133;
t230 = t132 * t133;
t64 = (t130 * t135 + t131 * t230) * t127;
t283 = qJD(1) * t64 - t130 * t91 - t136 * t195 - t99 * t212;
t129 = cos(pkin(6));
t237 = t129 * t131;
t192 = t135 * t222;
t172 = qJD(3) - t192;
t86 = t136 * qJD(2) + t172;
t282 = -qJD(1) * t237 + t134 * t86;
t243 = qJD(5) * t95;
t36 = -t133 * qJDD(4) + t286 * t130 + t243;
t126 = sin(pkin(10));
t128 = cos(pkin(10));
t236 = t129 * t132;
t81 = t126 * t135 + t128 * t236;
t83 = -t126 * t236 + t128 * t135;
t280 = -g(1) * t83 - g(2) * t81;
t241 = t127 * t131;
t235 = t129 * t135;
t80 = t126 * t132 - t128 * t235;
t82 = t126 * t235 + t128 * t132;
t238 = t127 * t135;
t84 = t134 * t238 + t237;
t156 = g(3) * t84 - g(2) * (t128 * t241 + t134 * t80) - g(1) * (-t126 * t241 + t134 * t82);
t39 = -qJD(4) * pkin(4) - t282;
t93 = t130 * t217 - t211;
t15 = pkin(5) * t93 - qJ(6) * t95 + t39;
t272 = pkin(9) * t90;
t279 = t118 * t15 - t272;
t226 = qJDD(1) - g(3);
t240 = t127 * t132;
t276 = -t226 * t240 - t280;
t193 = t132 * t222;
t220 = qJD(2) * qJ(3);
t98 = t193 + t220;
t275 = qJD(4) * (t193 - t98 - t220) - qJDD(4) * t136;
t274 = t95 ^ 2;
t273 = pkin(5) * t90;
t266 = t95 * t93;
t214 = qJD(4) * t134;
t265 = qJ(6) * t214 + (-t136 * t213 + qJD(6)) * t131 - t283;
t232 = t130 * t136;
t185 = -pkin(5) + t232;
t264 = qJD(5) * t247 + t185 * t214 + t284;
t97 = t178 * qJD(2);
t263 = t130 * t97 + t133 * t282;
t173 = pkin(5) * t130 - qJ(6) * t133;
t221 = qJD(1) * t134;
t113 = t129 * t221;
t53 = t131 * t86 + t113;
t262 = -t130 * qJD(6) + t118 * t173 - t53;
t261 = pkin(9) * qJD(5);
t260 = qJ(6) * t90;
t40 = qJD(4) * pkin(9) + t53;
t61 = qJD(2) * t99 + t193;
t14 = t130 * t61 + t133 * t40;
t11 = qJ(6) * t118 + t14;
t259 = t11 * t118;
t258 = t118 * t14;
t257 = t118 * t93;
t256 = t118 * t95;
t253 = t133 * t95;
t252 = t133 * t99;
t249 = t134 * t93;
t248 = t35 * t130;
t246 = qJD(2) * t98;
t245 = qJD(4) * t93;
t242 = qJDD(2) * pkin(2);
t239 = t127 * t134;
t234 = t130 * t131;
t138 = qJD(2) ^ 2;
t228 = t135 * t138;
t13 = -t130 * t40 + t133 * t61;
t227 = qJD(6) - t13;
t225 = pkin(2) * t238 + qJ(3) * t240;
t125 = t134 ^ 2;
t224 = t131 ^ 2 - t125;
t137 = qJD(4) ^ 2;
t223 = -t137 - t138;
t219 = qJD(2) * t127;
t215 = qJD(4) * t131;
t209 = qJDD(1) * t127;
t208 = qJDD(1) * t129;
t207 = qJDD(2) * qJ(3);
t206 = qJDD(4) * t131;
t189 = t132 * t209;
t27 = t189 + t99 * qJDD(2) + (t91 + t192) * qJD(2);
t187 = t134 * t208;
t104 = qJD(2) * t193;
t188 = t135 * t209;
t159 = qJDD(3) + t104 - t188;
t56 = t136 * qJDD(2) + t159;
t9 = qJDD(4) * pkin(9) + qJD(4) * t282 + t131 * t56 + t187;
t202 = -t130 * t27 - t133 * t9 - t61 * t212;
t200 = -qJD(4) * t113 - t131 * t208 - t86 * t215;
t197 = t132 * t219;
t196 = t135 * t219;
t186 = t130 * t9 - t133 * t27 + t40 * t212 + t61 * t213;
t184 = -t56 + t246;
t183 = qJD(5) * t93 - t35;
t182 = -t36 + t243;
t180 = t95 * t193;
t8 = -pkin(5) * t118 + t227;
t176 = t11 * t133 + t130 * t8;
t175 = t11 * t130 - t133 * t8;
t174 = pkin(5) * t133 + qJ(6) * t130;
t171 = t132 * (-qJD(2) * pkin(2) + t172) + t135 * t98;
t169 = -g(1) * t82 - g(2) * t80 + g(3) * t238;
t168 = qJDD(2) * t132 + t228;
t167 = pkin(4) + t174;
t166 = -t136 + t173;
t85 = t129 * t134 - t131 * t238;
t164 = t127 * t230 - t130 * t85;
t49 = t127 * t233 + t133 * t85;
t163 = -t213 * t40 - t202;
t161 = t118 * t212 + t130 * t90;
t160 = -t118 * t213 + t255;
t46 = -qJD(4) * t84 + t131 * t197;
t4 = qJD(5) * t49 + t46 * t130 - t133 * t196;
t47 = qJD(4) * t85 - t134 * t197;
t158 = -t118 * t4 + t164 * t90 + t84 * t36 + t47 * t93;
t30 = t80 * t133 + t81 * t234;
t32 = t82 * t133 + t83 * t234;
t63 = t281 * t127;
t157 = g(1) * t32 + g(2) * t30 + g(3) * t63;
t43 = t126 * t239 + t131 * t82;
t45 = t128 * t239 - t80 * t131;
t155 = g(1) * t43 - g(2) * t45 + g(3) * t85;
t31 = -t130 * t80 + t81 * t231;
t33 = -t130 * t82 + t83 * t231;
t153 = -g(1) * t33 - g(2) * t31 - g(3) * t64 + t193 * t249;
t151 = g(3) * t240 - t280;
t10 = -qJDD(4) * pkin(4) - t134 * t56 - t200;
t150 = t118 * t39 - t272;
t149 = -t169 + t188;
t5 = qJD(5) * t164 + t130 * t196 + t46 * t133;
t148 = t118 * t5 + t35 * t84 - t47 * t95 + t49 * t90;
t21 = t130 * t43 - t83 * t133;
t23 = -t130 * t45 - t81 * t133;
t147 = g(1) * t21 + g(2) * t23 - g(3) * t164 - t186;
t146 = t118 * t261 - t156;
t145 = qJDD(3) - t149;
t3 = pkin(5) * t36 + qJ(6) * t35 - qJD(6) * t95 + t10;
t144 = -t146 - t3;
t1 = qJD(6) * t118 + t163 + t260;
t2 = qJDD(6) + t186 - t273;
t143 = -qJD(5) * t175 + t1 * t133 + t2 * t130;
t22 = t130 * t83 + t133 * t43;
t24 = t130 * t81 - t133 * t45;
t142 = -g(1) * t22 - g(2) * t24 - g(3) * t49 + t163;
t141 = t15 * t95 + qJDD(6) - t147;
t140 = -t90 * t234 - t134 * t36 + t93 * t215 + (-t130 * t214 - t133 * t181) * t118;
t57 = t189 + t207 + (qJD(3) + t192) * qJD(2);
t139 = qJD(2) * t172 - t136 * t137 - t151 + t207 + t57;
t122 = qJDD(4) * t134;
t79 = t168 * t127;
t78 = (-qJDD(2) * t135 + t132 * t138) * t127;
t74 = t82 * pkin(2);
t73 = t80 * pkin(2);
t65 = t166 * t134;
t62 = t159 - t242;
t55 = t131 * t185 - t252;
t54 = qJ(6) * t131 + t247;
t41 = pkin(5) * t95 + qJ(6) * t93;
t28 = (qJD(5) * t174 - qJD(6) * t133) * t134 - t166 * t215;
t26 = -pkin(5) * t217 + t130 * t282 - t133 * t97;
t25 = qJ(6) * t217 + t263;
t17 = t257 - t35;
t6 = [t226, 0, -t78, -t79, t78, t79, qJDD(1) * t129 ^ 2 - g(3) + (qJD(2) * t171 + t132 * t57 - t135 * t62) * t127, 0, 0, 0, 0, 0, -t47 * qJD(4) - t84 * qJDD(4) + (t131 * t168 + t132 * t190) * t127, -t46 * qJD(4) - t85 * qJDD(4) + (t286 * t132 + t134 * t228) * t127, 0, 0, 0, 0, 0, t158, -t148, t158, t164 * t35 - t36 * t49 + t4 * t95 - t5 * t93, t148, t1 * t49 + t11 * t5 + t15 * t47 - t164 * t2 + t3 * t84 + t4 * t8 - g(3); 0, qJDD(2), t149, t276, t145 - 0.2e1 * t242, 0.2e1 * qJD(2) * qJD(3) + 0.2e1 * t207 - t276, t57 * qJ(3) + t98 * qJD(3) - t62 * pkin(2) - g(1) * (qJ(3) * t83 - t74) - g(2) * (qJ(3) * t81 - t73) - g(3) * t225 - t171 * t222, qJDD(2) * t125 - 0.2e1 * t131 * t190, -0.2e1 * t131 * t203 + 0.2e1 * t210 * t224, -t131 * t137 + t122, -t134 * t137 - t206, 0, t139 * t131 - t134 * t275, t131 * t275 + t139 * t134, -t95 * t194 + (-t215 * t95 - t251) * t133 (t130 * t95 + t133 * t93) * t215 + (t248 - t133 * t36 + (t130 * t93 - t253) * qJD(5)) * t134 (-t118 * t211 - t35) * t131 + (t160 + t244) * t134 (t118 * t216 - t36) * t131 + (-t161 - t245) * t134, t118 * t214 + t131 * t90, t90 * t252 + (-t213 * t99 - t284) * t118 + (-t39 * t216 + (-t161 + t245) * t136 - t186) * t131 + (t39 * t212 + t10 * t130 - t136 * t36 + (-t118 * t232 + t13) * qJD(4)) * t134 + t153, -t247 * t90 + t283 * t118 + ((t118 * t136 + t40) * t213 + (-t133 * t39 + t136 * t95) * qJD(4) + t202) * t131 + (-qJD(4) * t14 + t10 * t133 + t136 * t35 - t213 * t39 + t180) * t134 + t157, t28 * t93 + t36 * t65 - t55 * t90 + (-t15 * t216 - t2) * t131 - t264 * t118 + (-qJD(4) * t8 + t130 * t3 + t15 * t212) * t134 + t153, -t35 * t55 - t36 * t54 + t264 * t95 - t265 * t93 + t175 * t215 + (-qJD(5) * t176 - t1 * t130 + t133 * t2 + t151) * t134, -t28 * t95 + t35 * t65 + t54 * t90 + (t15 * t211 + t1) * t131 + t265 * t118 + (qJD(4) * t11 - t133 * t3 + t15 * t213 - t180) * t134 - t157, t1 * t54 + t3 * t65 + t15 * t28 + t2 * t55 - g(1) * (pkin(5) * t33 - pkin(8) * t82 + qJ(6) * t32 - t74) - g(2) * (pkin(5) * t31 - pkin(8) * t80 + qJ(6) * t30 - t73) - g(3) * (pkin(5) * t64 + qJ(6) * t63 + t225) + t264 * t8 + t265 * t11 + (-g(3) * pkin(8) * t135 + (-g(3) * t177 + t15 * t221) * t132) * t127 + t280 * t99; 0, 0, 0, 0, qJDD(2), -t138, t104 + t145 - t242 - t246, 0, 0, 0, 0, 0, t131 * t223 + t122, t134 * t223 - t206, 0, 0, 0, 0, 0, t140, t285, t140 (qJD(2) * t95 + t131 * t182 - t214 * t93) * t133 + (qJD(2) * t93 + t131 * t183 + t214 * t95) * t130, -t285, -t175 * qJD(2) + (qJD(4) * t176 - t3) * t134 + (qJD(4) * t15 + t143) * t131 + t169; 0, 0, 0, 0, 0, 0, 0, t134 * t138 * t131, -t224 * t138, t203, -t204, qJDD(4), qJD(4) * t53 - t134 * t184 + t156 + t200, t184 * t131 + t155 - t187, t118 * t253 - t248 (-t35 - t257) * t133 + (-t36 - t256) * t130 (t118 * t231 - t134 * t95) * qJD(2) + t161 (-t118 * t234 + t249) * qJD(2) + t160, -t118 * t217, -t13 * t217 - pkin(4) * t36 - t53 * t93 + (t118 * t282 + t150) * t130 + (-t10 + (-t97 - t261) * t118 + t156) * t133, pkin(4) * t35 + t263 * t118 + t14 * t217 - t53 * t95 + t150 * t133 + (t10 + t146) * t130, t118 * t26 + t130 * t279 + t144 * t133 - t167 * t36 + t8 * t217 + t262 * t93, t25 * t93 - t26 * t95 + (pkin(9) * t182 + t118 * t8 + t1) * t133 + (pkin(9) * t183 + t2 - t259) * t130 - t155, -t11 * t217 - t118 * t25 + t144 * t130 - t133 * t279 - t167 * t35 - t262 * t95, -t11 * t25 - t8 * t26 + t262 * t15 + (t143 - t155) * pkin(9) + (-t3 + t156) * t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t266, -t93 ^ 2 + t274, t17, t256 - t36, t90, -t39 * t95 + t147 + t258, t118 * t13 + t39 * t93 - t142, -t41 * t93 - t141 + t258 + 0.2e1 * t273, pkin(5) * t35 - qJ(6) * t36 + (t11 - t14) * t95 + (t8 - t227) * t93, 0.2e1 * t260 - t15 * t93 + t41 * t95 + (0.2e1 * qJD(6) - t13) * t118 + t142, t1 * qJ(6) - t2 * pkin(5) - t15 * t41 - t8 * t14 - g(1) * (-pkin(5) * t21 + qJ(6) * t22) - g(2) * (-pkin(5) * t23 + qJ(6) * t24) - g(3) * (pkin(5) * t164 + qJ(6) * t49) + t227 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t90 + t266, t17, -t118 ^ 2 - t274, t141 - t259 - t273;];
tau_reg  = t6;
