% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPPRR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:20:07
% EndTime: 2019-03-08 19:20:13
% DurationCPUTime: 3.75s
% Computational Cost: add. (3581->397), mult. (7812->525), div. (0->0), fcn. (6345->12), ass. (0->212)
t123 = sin(pkin(11));
t126 = cos(pkin(11));
t131 = sin(qJ(2));
t134 = cos(qJ(2));
t286 = -t131 * t123 + t134 * t126;
t130 = sin(qJ(5));
t133 = cos(qJ(5));
t161 = t130 * pkin(5) - t133 * pkin(9) + qJ(4);
t125 = sin(pkin(6));
t222 = qJD(1) * t125;
t188 = qJD(2) * t222;
t208 = t125 * qJDD(1);
t281 = t131 * t208 + t134 * t188;
t104 = t134 * t208;
t247 = qJDD(2) * pkin(2);
t73 = -t131 * t188 + t104 + t247;
t37 = t123 * t73 + t281 * t126;
t174 = pkin(5) * t133 + pkin(9) * t130;
t81 = qJD(5) * t174 + qJD(4);
t12 = qJD(2) * t81 + qJDD(2) * t161 + t37;
t129 = sin(qJ(6));
t132 = cos(qJ(6));
t128 = cos(pkin(6));
t108 = t128 * qJD(1) + qJD(3);
t190 = t134 * t222;
t90 = qJD(2) * pkin(2) + t190;
t193 = t131 * t222;
t95 = t123 * t193;
t61 = t126 * t90 - t95;
t175 = qJD(4) - t61;
t274 = -pkin(3) - pkin(8);
t53 = t274 * qJD(2) + t175;
t35 = t133 * t108 + t130 * t53;
t33 = qJD(5) * pkin(9) + t35;
t62 = t123 * t90 + t126 * t193;
t42 = qJD(2) * t161 + t62;
t169 = t129 * t33 - t132 * t42;
t105 = t128 * qJDD(1) + qJDD(3);
t219 = qJD(5) * t133;
t204 = t281 * t123 - t126 * t73;
t181 = qJDD(4) + t204;
t23 = t274 * qJDD(2) + t181;
t205 = -t133 * t105 - t130 * t23 - t53 * t219;
t220 = qJD(5) * t130;
t7 = -t108 * t220 - t205;
t5 = qJDD(5) * pkin(9) + t7;
t1 = -t169 * qJD(6) + t129 * t12 + t132 * t5;
t214 = t130 * qJD(2);
t111 = qJD(6) + t214;
t285 = -t169 * t111 - t1;
t72 = t126 * t190 - t95;
t284 = -t72 + t81;
t124 = sin(pkin(10));
t127 = cos(pkin(10));
t239 = t128 * t134;
t158 = -t124 * t239 - t127 * t131;
t151 = t286 * t128;
t165 = t134 * t123 + t131 * t126;
t47 = -t124 * t151 - t127 * t165;
t283 = pkin(2) * t158 + t47 * pkin(3);
t237 = t130 * t108;
t34 = t133 * t53 - t237;
t183 = t130 * t105 - t133 * t23;
t8 = -t35 * qJD(5) - t183;
t139 = -(t130 * t34 - t133 * t35) * qJD(5) + t7 * t130 + t8 * t133;
t278 = t286 * t125;
t44 = -t124 * t165 + t127 * t151;
t153 = g(1) * t47 + g(2) * t44 + g(3) * t278;
t282 = t139 + t153;
t10 = t129 * t42 + t132 * t33;
t2 = -qJD(6) * t10 + t132 * t12 - t129 * t5;
t279 = t10 * t111 + t2;
t215 = t129 * qJD(5);
t189 = t130 * t215;
t206 = t133 * qJDD(2);
t221 = qJD(2) * t133;
t86 = t132 * t221 + t215;
t248 = qJD(6) * t86;
t52 = -qJD(2) * t189 - t132 * qJDD(5) + t129 * t206 + t248;
t226 = t165 * t128;
t43 = -t124 * t286 - t127 * t226;
t48 = -t124 * t226 + t127 * t286;
t76 = t165 * t125;
t154 = -g(1) * t48 + g(2) * t43 - g(3) * t76;
t172 = t10 * t132 + t129 * t169;
t6 = -qJDD(5) * pkin(5) - t8;
t277 = qJD(5) * t172 - t6;
t216 = qJD(6) * t133;
t191 = t129 * t216;
t213 = t132 * qJD(5);
t149 = t130 * t213 + t191;
t233 = t132 * t133;
t212 = qJD(2) * qJD(5);
t187 = t133 * t212;
t207 = t130 * qJDD(2);
t80 = qJDD(6) + t187 + t207;
t276 = -t111 * t149 + t80 * t233;
t51 = qJD(2) * t149 - qJD(6) * t213 - t129 * qJDD(5) - t132 * t206;
t114 = -t126 * pkin(2) - pkin(3);
t110 = -pkin(8) + t114;
t267 = t123 * pkin(2);
t112 = qJ(4) + t267;
t57 = qJD(2) * qJ(4) + t62;
t69 = qJD(1) * t76;
t275 = qJDD(5) * t110 + (qJD(2) * t112 + t57 - t69) * qJD(5);
t270 = t44 * pkin(8);
t269 = t47 * pkin(8);
t268 = t278 * pkin(8);
t84 = t129 * t221 - t213;
t266 = t86 * t84;
t192 = t110 * t219;
t236 = t130 * t132;
t238 = t129 * t130;
t78 = t161 + t267;
t54 = -t110 * t238 + t132 * t78;
t264 = qJD(6) * t54 + t284 * t129 + t132 * t192 - t69 * t236;
t55 = t110 * t236 + t129 * t78;
t263 = -qJD(6) * t55 - t129 * t192 + t284 * t132 + t69 * t238;
t77 = t86 * t219;
t262 = -t51 * t130 + t77;
t258 = t129 * t80;
t257 = t130 * t52;
t256 = t132 * t80;
t255 = t133 * t51;
t254 = t133 * t52;
t253 = t84 * t111;
t252 = t86 * t111;
t242 = t125 * t134;
t251 = pkin(2) * t242 + pkin(3) * t278;
t250 = qJD(5) * t84;
t249 = qJD(6) * t84;
t246 = qJDD(2) * pkin(3);
t245 = t124 * t131;
t244 = t125 * t130;
t243 = t125 * t133;
t240 = t128 * t131;
t230 = t57 * qJD(2);
t229 = t69 * qJD(2);
t228 = qJD(4) - t72;
t227 = qJDD(1) - g(3);
t121 = t130 ^ 2;
t122 = t133 ^ 2;
t225 = t121 - t122;
t224 = t121 + t122;
t135 = qJD(5) ^ 2;
t136 = qJD(2) ^ 2;
t223 = -t135 - t136;
t218 = qJD(6) * t129;
t217 = qJD(6) * t132;
t210 = qJDD(5) * t130;
t209 = qJDD(5) * t133;
t202 = t84 * t220;
t201 = t86 * t220;
t200 = t86 * t216;
t196 = t127 * t239;
t195 = t133 * t136 * t130;
t185 = -t51 + t249;
t184 = -t52 + t248;
t182 = t224 * qJDD(2);
t180 = qJD(6) * t130 + qJD(2);
t178 = t130 * t187;
t177 = t76 * qJ(4) + t251;
t171 = t10 * t129 - t132 * t169;
t59 = t128 * t133 - t130 * t278;
t22 = t76 * t129 + t59 * t132;
t21 = -t59 * t129 + t76 * t132;
t168 = -t35 * t130 - t34 * t133;
t58 = t128 * t130 + t133 * t278;
t99 = pkin(2) * t196;
t166 = -pkin(2) * t245 + t44 * pkin(3) + t99;
t70 = qJD(2) * t76;
t164 = qJD(2) * t70 - qJDD(2) * t278;
t71 = t278 * qJD(2);
t163 = qJD(2) * t71 + qJDD(2) * t76;
t24 = qJDD(2) * qJ(4) + qJD(4) * qJD(2) + t37;
t162 = t24 * t76 + t57 * t71 - g(3);
t159 = t111 * t217 + t258;
t157 = t24 * t112 + t228 * t57;
t156 = -g(1) * (-t124 * t244 - t47 * t133) - g(2) * (t127 * t244 - t44 * t133) + g(3) * t58;
t28 = t124 * t243 - t47 * t130;
t30 = t127 * t243 + t44 * t130;
t155 = -g(1) * t28 + g(2) * t30 - g(3) * t59;
t152 = -g(3) * t128 + (-g(1) * t124 + g(2) * t127) * t125;
t150 = t156 - t6;
t148 = -t43 * qJ(4) + t166;
t32 = -qJD(5) * pkin(5) - t34;
t147 = -pkin(9) * t80 + t111 * t32;
t146 = t154 + t37;
t145 = -t153 - t204;
t144 = pkin(9) * qJD(6) * t111 - t150;
t143 = -g(1) * t158 - g(3) * t242;
t142 = qJ(4) * t48 + t283;
t141 = -qJD(6) * t171 + t1 * t132 - t2 * t129;
t140 = -t145 - t229;
t138 = qJD(5) * t32 + t141;
t137 = t228 * qJD(2) + t112 * qJDD(2) - t110 * t135 + t154 + t24;
t98 = -t135 * t130 + t209;
t97 = -t135 * t133 - t210;
t91 = t105 * t128;
t89 = t174 * qJD(2);
t79 = t111 * t189;
t60 = t152 + t105;
t56 = -qJD(2) * pkin(3) + t175;
t39 = t52 * t233;
t31 = t181 - t246;
t20 = qJD(5) * t59 - t70 * t133;
t19 = -qJD(5) * t58 + t70 * t130;
t14 = t129 * t89 + t132 * t34;
t13 = -t129 * t34 + t132 * t89;
t4 = qJD(6) * t21 + t71 * t129 + t19 * t132;
t3 = -qJD(6) * t22 - t19 * t129 + t71 * t132;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t227, 0, 0, 0, 0, 0, 0 (qJDD(2) * t134 - t131 * t136) * t125 (-qJDD(2) * t131 - t134 * t136) * t125, 0, -g(3) + (t128 ^ 2 + (t131 ^ 2 + t134 ^ 2) * t125 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t164, -t163, 0, -t204 * t278 + t37 * t76 - t61 * t70 + t62 * t71 - g(3) + t91, 0, 0, 0, 0, 0, 0, 0, t164, t163, -t278 * t31 + t56 * t70 + t162 + t91, 0, 0, 0, 0, 0, 0, t76 * t207 - t20 * qJD(5) - t58 * qJDD(5) + (t130 * t71 + t219 * t76) * qJD(2), t76 * t206 - t19 * qJD(5) - t59 * qJDD(5) + (t133 * t71 - t220 * t76) * qJD(2) (-t130 * t59 + t133 * t58) * qJDD(2) + (-t130 * t19 + t133 * t20 + (-t130 * t58 - t133 * t59) * qJD(5)) * qJD(2), t19 * t35 - t20 * t34 - t58 * t8 + t59 * t7 + t162, 0, 0, 0, 0, 0, 0, t3 * t111 + t20 * t84 + t21 * t80 + t58 * t52, -t4 * t111 + t20 * t86 - t22 * t80 - t58 * t51, t21 * t51 - t22 * t52 - t3 * t86 - t4 * t84, t1 * t22 + t10 * t4 - t169 * t3 + t2 * t21 + t20 * t32 + t58 * t6 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t104 - g(2) * (t196 - t245) + t143, -g(1) * (t124 * t240 - t127 * t134) - g(2) * (-t124 * t134 - t127 * t240) - t227 * t131 * t125, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t126 * t247 - t140, t72 * qJD(2) - t123 * t247 - t146, 0, -g(2) * t99 + t61 * t69 - t62 * t72 + (g(2) * t245 + t37 * t123 - t126 * t204 + t143) * pkin(2), qJDD(2), 0, 0, 0, 0, 0, 0, qJDD(4) + (-pkin(3) + t114) * qJDD(2) + t140 (qJ(4) + t112) * qJDD(2) + (0.2e1 * qJD(4) - t72) * qJD(2) + t146, -g(1) * t142 - g(2) * t148 - g(3) * t177 + t31 * t114 - t56 * t69 + t157, t122 * qJDD(2) - 0.2e1 * t178, -0.2e1 * t130 * t206 + 0.2e1 * t212 * t225, t98, t121 * qJDD(2) + 0.2e1 * t178, t97, 0, t137 * t130 + t275 * t133, -t275 * t130 + t137 * t133, -t110 * t182 + t224 * t229 - t282, -g(1) * (t142 + t269) - g(2) * (t148 + t270) - g(3) * (t177 + t268) + t168 * t69 + t139 * t110 + t157, -t86 * t191 + (-t201 - t255) * t132, -t39 + (-t200 + t202) * t132 + (t201 + (t51 + t249) * t133) * t129, t262 + t276, t84 * t132 * t216 + (-t202 + t254) * t129, -t257 + t79 + (-t159 - t250) * t133, t111 * t219 + t80 * t130, t54 * t80 + t263 * t111 - t153 * t129 + (t2 + (t110 * t84 - t129 * t32) * qJD(5) + t154 * t132) * t130 + (-qJD(5) * t169 - t110 * t52 + t6 * t129 + t217 * t32 + t69 * t84) * t133, -t55 * t80 - t264 * t111 - t153 * t132 + (-t1 + (t110 * t86 - t132 * t32) * qJD(5) - t154 * t129) * t130 + (-t10 * qJD(5) + t110 * t51 + t6 * t132 - t218 * t32 + t69 * t86) * t133, t54 * t51 - t55 * t52 - t263 * t86 - t264 * t84 + t171 * t220 + (-qJD(6) * t172 - t1 * t129 - t132 * t2 - t154) * t133, t1 * t55 + t2 * t54 + t32 * t133 * t69 - g(1) * (t269 + t283) - g(2) * (t166 + t270) - g(3) * (t251 + t268) - t263 * t169 + (-t6 * t133 + t220 * t32) * t110 + t264 * t10 + t154 * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, 0, 0, 0, 0, 0, t97, -t98, 0, qJD(5) * t168 - t8 * t130 + t7 * t133 + t152, 0, 0, 0, 0, 0, 0, t257 + t79 + (-t159 + t250) * t133, t262 - t276, -t39 + (t200 + t202) * t132 + (t133 * t185 - t201) * t129, -t277 * t130 + t138 * t133 + t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t136, qJDD(4) - t145 - t230 - t246, 0, 0, 0, 0, 0, 0, t130 * t223 + t209, t133 * t223 - t210, -t182, -t230 + t282, 0, 0, 0, 0, 0, 0, -t254 + (t250 - t258) * t130 + (-t132 * t180 - t133 * t215) * t111, t255 + (qJD(5) * t86 - t256) * t130 + (t129 * t180 - t133 * t213) * t111 (qJD(2) * t86 + t130 * t184 - t219 * t84) * t132 + (qJD(2) * t84 + t130 * t185 + t77) * t129, -t171 * qJD(2) + t138 * t130 + t277 * t133 + t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, -t225 * t136, t206, -t195, -t207, qJDD(5), -t221 * t57 + t156 - t183, t57 * t214 + (t34 + t237) * qJD(5) - t155 + t205, 0, 0, -t51 * t129 + t132 * t252 (-t51 - t253) * t132 + (-t52 - t252) * t129 (t111 * t236 - t133 * t86) * qJD(2) + t159, t129 * t253 - t52 * t132, -t111 * t218 + t256 + (-t111 * t238 + t133 * t84) * qJD(2), -t111 * t221, -pkin(5) * t52 - t13 * t111 + t129 * t147 - t132 * t144 + t169 * t221 - t35 * t84, pkin(5) * t51 + t10 * t221 + t14 * t111 + t129 * t144 + t132 * t147 - t35 * t86, t13 * t86 + t14 * t84 + (pkin(9) * t184 - t285) * t132 + (pkin(9) * t185 - t279) * t129 + t155, -t10 * t14 + t169 * t13 - t32 * t35 + t150 * pkin(5) + (t141 + t155) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t266, -t84 ^ 2 + t86 ^ 2, t253 - t51, -t266, t252 - t52, t80, -t32 * t86 - g(1) * (-t28 * t129 + t48 * t132) - g(2) * (t30 * t129 - t132 * t43) - g(3) * t21 + t279, t32 * t84 - g(1) * (-t48 * t129 - t28 * t132) - g(2) * (t129 * t43 + t30 * t132) + g(3) * t22 + t285, 0, 0;];
tau_reg  = t9;
