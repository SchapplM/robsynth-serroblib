% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PPRRRP1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tau_reg [6x21]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PPRRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:54:32
% EndTime: 2019-03-08 18:54:39
% DurationCPUTime: 3.31s
% Computational Cost: add. (3457->368), mult. (8803->541), div. (0->0), fcn. (8067->14), ass. (0->194)
t136 = sin(qJ(4));
t139 = cos(qJ(4));
t182 = pkin(4) * t136 - pkin(10) * t139;
t105 = t182 * qJD(4);
t237 = cos(pkin(6));
t117 = t237 * qJD(1) + qJD(2);
t132 = sin(pkin(6));
t130 = sin(pkin(12));
t137 = sin(qJ(3));
t140 = cos(qJ(3));
t133 = cos(pkin(12));
t236 = cos(pkin(7));
t190 = t133 * t236;
t165 = t130 * t140 + t137 * t190;
t157 = t165 * t132;
t131 = sin(pkin(7));
t227 = t131 * t137;
t57 = qJD(1) * t157 + t117 * t227;
t272 = t105 - t57;
t114 = t237 * qJDD(1) + qJDD(2);
t186 = t132 * t190;
t176 = qJD(1) * t186;
t219 = qJD(3) * t137;
t201 = t131 * t219;
t221 = qJD(1) * t132;
t202 = t130 * t221;
t228 = t130 * t137;
t205 = t132 * t228;
t218 = qJD(3) * t140;
t146 = -t140 * (qJDD(1) * t186 + t114 * t131) + qJDD(1) * t205 + t117 * t201 + t176 * t219 + t202 * t218;
t235 = cos(pkin(11));
t181 = t237 * t235;
t234 = sin(pkin(11));
t153 = t234 * t130 - t133 * t181;
t192 = t132 * t235;
t264 = t131 * t192 + t153 * t236;
t81 = t130 * t181 + t234 * t133;
t43 = t81 * t137 + t140 * t264;
t180 = t237 * t234;
t154 = t235 * t130 + t133 * t180;
t191 = t132 * t234;
t263 = -t131 * t191 + t154 * t236;
t82 = -t130 * t180 + t235 * t133;
t45 = t82 * t137 + t140 * t263;
t184 = t140 * t190;
t193 = t131 * t237;
t185 = t140 * t193;
t65 = -t132 * t184 - t185 + t205;
t167 = g(1) * t45 + g(2) * t43 + g(3) * t65;
t271 = t57 * qJD(3) - t146 + t167;
t220 = qJD(3) * t136;
t270 = qJD(5) * t220 - qJDD(4);
t135 = sin(qJ(5));
t138 = cos(qJ(5));
t207 = t136 * qJDD(3);
t209 = t139 * qJD(3);
t59 = t135 * ((qJD(5) + t209) * qJD(4) + t207) + t270 * t138;
t121 = t138 * pkin(5) + pkin(4);
t269 = -t121 * t139 - pkin(3);
t203 = pkin(5) * t135 + pkin(9);
t212 = t135 * qJD(4);
t225 = t135 * t139;
t56 = -t137 * t202 + t140 * (t117 * t131 + t176);
t268 = -t136 * pkin(9) * t212 - t272 * t138 - t56 * t225;
t108 = -t139 * pkin(4) - t136 * pkin(10) - pkin(3);
t213 = qJD(5) * t138;
t223 = t138 * t139;
t267 = t108 * t213 + t272 * t135 - t56 * t223;
t55 = qJD(3) * pkin(9) + t57;
t204 = t132 * t133 * t131;
t77 = -qJD(1) * t204 + t236 * t117;
t266 = -t136 * t55 + t139 * t77;
t143 = t153 * t131 - t236 * t192;
t44 = -t137 * t264 + t81 * t140;
t22 = t136 * t143 + t44 * t139;
t144 = t154 * t131 + t236 * t191;
t46 = -t137 * t263 + t82 * t140;
t24 = t136 * t144 + t46 * t139;
t159 = t237 * t236 - t204;
t66 = t137 * t193 + t157;
t49 = t136 * t159 + t66 * t139;
t25 = -t49 * t135 + t65 * t138;
t265 = -g(1) * (-t24 * t135 + t45 * t138) - g(2) * (-t22 * t135 + t43 * t138) - g(3) * t25;
t262 = qJDD(1) * t165;
t261 = t59 * pkin(5) + qJDD(6);
t118 = -qJD(5) + t209;
t41 = t136 * t77 + t139 * t55;
t38 = qJD(4) * pkin(10) + t41;
t259 = (pkin(9) * t118 + t38) * qJD(5) + t167;
t98 = t138 * t220 + t212;
t258 = t98 ^ 2;
t50 = qJD(3) * t108 - t56;
t12 = -t135 * t38 + t138 * t50;
t10 = -t98 * qJ(6) + t12;
t253 = t118 * pkin(5);
t9 = t10 - t253;
t255 = t10 - t9;
t252 = qJ(6) + pkin(10);
t119 = pkin(9) * t223;
t175 = pkin(5) * t136 - qJ(6) * t223;
t210 = t138 * qJD(6);
t233 = qJ(6) * t136;
t251 = -t136 * t210 + t175 * qJD(4) + (-t119 + (-t108 + t233) * t135) * qJD(5) - t268;
t104 = t182 * qJD(3);
t250 = t135 * t104 + t138 * t266;
t224 = t136 * t138;
t249 = (-pkin(9) * qJD(4) - qJ(6) * qJD(5)) * t224 + (-qJD(6) * t136 + (-pkin(9) * qJD(5) - qJ(6) * qJD(4)) * t139) * t135 + t267;
t194 = qJD(5) * t252;
t248 = t210 - t250 + (qJ(6) * t209 - t194) * t135;
t87 = t138 * t104;
t247 = -qJD(3) * t175 - t138 * t194 - t87 + (-qJD(6) + t266) * t135;
t245 = qJD(3) * pkin(3);
t76 = -qJDD(1) * t204 + t236 * t114;
t244 = t136 * t76;
t208 = qJD(3) * qJD(4);
t197 = t139 * t208;
t211 = t138 * qJD(4);
t58 = -qJD(5) * t211 + (-t197 - t207) * t138 + t270 * t135;
t242 = t58 * t135;
t96 = t135 * t220 - t211;
t241 = t96 * t118;
t240 = t98 * t118;
t238 = t135 * t108 + t119;
t232 = qJD(4) * t96;
t230 = qJDD(4) * pkin(4);
t226 = t131 * t140;
t128 = t136 ^ 2;
t222 = -t139 ^ 2 + t128;
t217 = qJD(4) * t136;
t216 = qJD(4) * t139;
t215 = qJD(5) * t118;
t214 = qJD(5) * t135;
t125 = t139 * qJDD(3);
t16 = qJD(3) * t105 + qJDD(3) * t108 + t146;
t164 = t184 - t228;
t31 = qJDD(3) * pkin(9) + (t114 * t137 + t117 * t218) * t131 + (qJD(1) * qJD(3) * t164 + t262) * t132;
t7 = qJDD(4) * pkin(10) + qJD(4) * t266 + t139 * t31 + t244;
t206 = t135 * t16 + t138 * t7 + t50 * t213;
t200 = t131 * t218;
t199 = t118 * t211;
t195 = t136 * t208;
t13 = t135 * t50 + t138 * t38;
t26 = t65 * t135 + t49 * t138;
t178 = t136 * t31 - t139 * t76 + t55 * t216 + t77 * t217;
t142 = qJD(3) ^ 2;
t177 = qJDD(3) * t140 - t137 * t142;
t37 = -qJD(4) * pkin(4) - t266;
t84 = t136 * t236 + t139 * t227;
t69 = -t135 * t84 - t138 * t226;
t173 = t135 * t226 - t138 * t84;
t172 = t38 * t214 - t206;
t93 = qJDD(5) - t125 + t195;
t171 = -t118 * t213 + t135 * t93;
t170 = t118 * t214 + t138 * t93;
t21 = t44 * t136 - t139 * t143;
t23 = t46 * t136 - t139 * t144;
t48 = t66 * t136 - t139 * t159;
t169 = g(1) * t23 + g(2) * t21 + g(3) * t48;
t168 = g(1) * t24 + g(2) * t22 + g(3) * t49;
t166 = g(1) * t46 + g(2) * t44 + g(3) * t66;
t83 = t136 * t227 - t139 * t236;
t8 = t178 - t230;
t163 = -pkin(10) * t93 - t118 * t37;
t54 = -t56 - t245;
t160 = -pkin(9) * qJDD(4) + (t54 + t56 - t245) * qJD(4);
t156 = -g(1) * t191 + g(2) * t192 - g(3) * t237;
t15 = t138 * t16;
t155 = -qJD(5) * t13 - t135 * t7 + t15;
t152 = pkin(10) * t215 + t169 - t8;
t151 = t169 - t178;
t141 = qJD(4) ^ 2;
t145 = 0.2e1 * qJDD(3) * pkin(3) - pkin(9) * t141 + t271;
t110 = t252 * t138;
t109 = t252 * t135;
t95 = t138 * t108;
t92 = t96 ^ 2;
t71 = -t135 * t233 + t238;
t68 = qJD(4) * t84 + t136 * t200;
t67 = -qJD(4) * t83 + t139 * t200;
t64 = -qJ(6) * t224 + t95 + (-pkin(9) * t135 - pkin(5)) * t139;
t61 = t66 * qJD(3);
t60 = (t164 * t132 + t185) * qJD(3);
t36 = qJD(5) * t173 - t135 * t67 + t138 * t201;
t35 = qJD(5) * t69 + t135 * t201 + t138 * t67;
t27 = t96 * pkin(5) + qJD(6) + t37;
t20 = -qJD(4) * t48 + t60 * t139;
t19 = qJD(4) * t49 + t60 * t136;
t11 = -t96 * qJ(6) + t13;
t5 = t8 + t261;
t4 = qJD(5) * t25 + t61 * t135 + t20 * t138;
t3 = -qJD(5) * t26 - t20 * t135 + t61 * t138;
t2 = -t59 * qJ(6) - t96 * qJD(6) - t172;
t1 = t93 * pkin(5) + t58 * qJ(6) - t98 * qJD(6) + t155;
t6 = [qJDD(1) - g(3), t114 * t237 - g(3) + (t130 ^ 2 + t133 ^ 2) * t132 ^ 2 * qJDD(1), 0, -t61 * qJD(3) - t65 * qJDD(3), -t60 * qJD(3) - t66 * qJDD(3), 0, 0, 0, 0, 0, -t65 * t125 - t19 * qJD(4) - t48 * qJDD(4) + (-t139 * t61 + t217 * t65) * qJD(3), t65 * t207 - t20 * qJD(4) - t49 * qJDD(4) + (t136 * t61 + t216 * t65) * qJD(3), 0, 0, 0, 0, 0, -t3 * t118 + t19 * t96 + t25 * t93 + t48 * t59, t4 * t118 + t19 * t98 - t26 * t93 - t48 * t58, t25 * t58 - t26 * t59 - t3 * t98 - t4 * t96, t1 * t25 + t11 * t4 + t27 * t19 + t2 * t26 + t9 * t3 + t5 * t48 - g(3); 0, t156 + t114, 0, t177 * t131 (-qJDD(3) * t137 - t140 * t142) * t131, 0, 0, 0, 0, 0, -t68 * qJD(4) - t83 * qJDD(4) + (t139 * t177 - t140 * t195) * t131, -t67 * qJD(4) - t84 * qJDD(4) + (-t136 * t177 - t140 * t197) * t131, 0, 0, 0, 0, 0, -t36 * t118 + t83 * t59 + t68 * t96 + t69 * t93, t35 * t118 + t173 * t93 - t83 * t58 + t68 * t98, t173 * t59 - t35 * t96 - t36 * t98 + t69 * t58, t1 * t69 + t11 * t35 - t173 * t2 + t27 * t68 + t9 * t36 + t5 * t83 + t156; 0, 0, qJDD(3), t271, -t114 * t227 - t132 * t262 + (-t117 * t226 - t164 * t221 + t56) * qJD(3) + t166, t128 * qJDD(3) + 0.2e1 * t139 * t195, 0.2e1 * t136 * t125 - 0.2e1 * t222 * t208, qJDD(4) * t136 + t141 * t139, qJDD(4) * t139 - t141 * t136, 0, t136 * t160 + t139 * t145, -t136 * t145 + t139 * t160, t98 * t139 * t211 + (-t58 * t138 - t214 * t98) * t136 (-t135 * t98 - t138 * t96) * t216 + (t242 - t138 * t59 + (t135 * t96 - t138 * t98) * qJD(5)) * t136 (t58 - t199) * t139 + (qJD(4) * t98 + t170) * t136 (t118 * t212 + t59) * t139 + (-t171 - t232) * t136, -t118 * t217 - t93 * t139, t95 * t93 + t268 * t118 + (t108 * t215 - t166) * t135 + (pkin(9) * t232 - t15 + (-pkin(9) * t93 + qJD(4) * t37 + qJD(5) * t50 + t7) * t135 + t259 * t138) * t139 + (pkin(9) * t59 + t12 * qJD(4) + t8 * t135 + t213 * t37 - t56 * t96) * t136, -t238 * t93 + t267 * t118 - t166 * t138 + ((pkin(9) * t98 + t37 * t138) * qJD(4) - t259 * t135 + t206) * t139 + (-t37 * t214 - t13 * qJD(4) + t8 * t138 - t56 * t98 + (-t58 - t199) * pkin(9)) * t136, t64 * t58 - t71 * t59 - t251 * t98 - t249 * t96 + (-t11 * t135 - t138 * t9) * t216 + (-t1 * t138 - t135 * t2 + (-t11 * t138 + t135 * t9) * qJD(5) + t167) * t136, t2 * t71 + t1 * t64 - g(1) * (t203 * t46 + t269 * t45) - g(2) * (t203 * t44 + t269 * t43) - g(3) * (t203 * t66 + t269 * t65) + t251 * t9 + t249 * t11 + t27 * t203 * t216 + (t5 * t203 + (pkin(5) * t213 - t56) * t27 + t167 * t252) * t136; 0, 0, 0, 0, 0, -t136 * t142 * t139, t222 * t142, t207, t125, qJDD(4), t41 * qJD(4) - t220 * t54 + t151, -t244 + (-qJD(3) * t54 - t31) * t139 + t168, -t138 * t240 - t242 (-t58 + t241) * t138 + (-t59 + t240) * t135 (t118 * t223 - t136 * t98) * qJD(3) + t171 (-t118 * t225 + t136 * t96) * qJD(3) + t170, t118 * t220, -t12 * t220 - pkin(4) * t59 + t87 * t118 - t41 * t96 + (-t118 * t266 + t163) * t135 + t152 * t138, pkin(4) * t58 - t250 * t118 + t13 * t220 - t152 * t135 + t163 * t138 - t41 * t98, -t109 * t58 - t110 * t59 - t247 * t98 - t248 * t96 + (t118 * t9 + t2) * t138 + (t11 * t118 - t1) * t135 - t168, t2 * t110 - t1 * t109 - t5 * t121 - g(1) * (-t23 * t121 + t24 * t252) - g(2) * (-t21 * t121 + t22 * t252) - g(3) * (-t48 * t121 + t252 * t49) + t247 * t9 + (-t135 * t253 - t41) * t27 + t248 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98 * t96, -t92 + t258, -t58 - t241, -t240 - t59, t93, -t13 * t118 - t37 * t98 + t155 + t265, -t12 * t118 + t37 * t96 - g(1) * (-t45 * t135 - t24 * t138) - g(2) * (-t43 * t135 - t22 * t138) + g(3) * t26 + t172, pkin(5) * t58 + t255 * t96, -t255 * t11 + (-t27 * t98 + t1 + t265) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92 - t258, t11 * t96 + t9 * t98 - t151 - t230 + t261;];
tau_reg  = t6;
