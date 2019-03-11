% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRP10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% 
% Output:
% tau_reg [6x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:32:35
% EndTime: 2019-03-09 03:32:43
% DurationCPUTime: 2.94s
% Computational Cost: add. (2969->430), mult. (5275->501), div. (0->0), fcn. (2904->6), ass. (0->207)
t133 = cos(qJ(3));
t205 = qJD(1) * qJD(3);
t188 = t133 * t205;
t130 = sin(qJ(3));
t201 = t130 * qJDD(1);
t276 = t188 + t201;
t134 = cos(qJ(1));
t119 = g(2) * t134;
t131 = sin(qJ(1));
t120 = g(1) * t131;
t267 = t120 - t119;
t214 = qJD(1) * t133;
t216 = qJD(1) * t130;
t222 = pkin(3) * t216 + qJD(1) * qJ(2);
t54 = -qJ(4) * t214 + t222;
t275 = -qJD(1) * t54 - t267;
t138 = qJD(1) ^ 2;
t237 = qJ(2) * t138;
t274 = t267 + t237;
t136 = -pkin(1) - pkin(7);
t248 = pkin(4) - t136;
t129 = sin(qJ(5));
t209 = qJD(5) * t129;
t132 = cos(qJ(5));
t111 = t133 * qJDD(1);
t189 = t130 * t205;
t266 = -t189 + t111;
t67 = -qJDD(5) - t266;
t52 = t132 * t67;
t97 = qJD(5) + t214;
t273 = -t97 * t209 - t52;
t215 = qJD(1) * t132;
t191 = t130 * t215;
t213 = qJD(3) * t129;
t68 = -t191 + t213;
t251 = t68 * t97;
t26 = qJD(3) * t209 - qJD(5) * t191 - t132 * qJDD(3) - t276 * t129;
t272 = t26 - t251;
t211 = qJD(3) * t132;
t70 = t129 * t216 + t211;
t249 = t70 * t97;
t231 = qJD(5) * t70;
t27 = qJDD(3) * t129 - t276 * t132 + t231;
t271 = -t27 + t249;
t234 = qJD(3) * t68;
t242 = t129 * t97;
t270 = -t242 * t97 - t234 - t52;
t269 = t130 * t267;
t89 = t136 * qJD(1) + qJD(2);
t77 = t133 * t89;
t268 = qJD(4) - t77;
t122 = qJDD(1) * qJ(2);
t124 = qJD(1) * qJD(2);
t198 = 0.2e1 * t124;
t265 = t198 + 0.2e1 * t122;
t118 = g(3) * t130;
t168 = t276 * pkin(3) + qJ(4) * t189 + t122 + t124;
t184 = qJD(3) * pkin(8) - qJD(4);
t206 = qJ(4) * qJDD(1);
t14 = pkin(8) * t201 + (t184 * qJD(1) - t206) * t133 + t168;
t135 = -pkin(3) - pkin(8);
t86 = t136 * qJDD(1) + qJDD(2);
t183 = -t133 * t86 + qJDD(4);
t212 = qJD(3) * t130;
t72 = t89 * t212;
t169 = t183 + t72;
t20 = pkin(4) * t266 + t135 * qJDD(3) + t169;
t208 = qJD(5) * t132;
t224 = pkin(4) * t214 + t268;
t36 = t135 * qJD(3) + t224;
t200 = t129 * t20 + t132 * t14 + t36 * t208;
t113 = t133 * qJ(4);
t254 = pkin(8) * t130;
t164 = -t113 + t254;
t39 = t164 * qJD(1) + t222;
t232 = qJD(5) * t39;
t225 = t133 * t134;
t58 = t129 * t225 + t131 * t132;
t227 = t131 * t133;
t60 = -t129 * t227 + t132 * t134;
t264 = -g(1) * t60 - g(2) * t58 - (t232 + t118) * t129 + t200;
t233 = qJD(3) * t70;
t243 = t129 * t67;
t258 = t97 ^ 2;
t263 = t132 * t258 + t233 - t243;
t121 = qJDD(3) * qJ(4);
t123 = qJD(3) * qJD(4);
t210 = qJD(3) * t133;
t75 = t130 * t86;
t29 = -t89 * t210 - t121 - t123 - t75;
t230 = qJDD(3) * pkin(3);
t30 = t169 - t230;
t185 = qJD(3) * pkin(3) - qJD(4);
t53 = -t185 - t77;
t125 = qJD(3) * qJ(4);
t76 = t130 * t89;
t55 = -t76 - t125;
t143 = -t29 * t130 - t30 * t133 + (t130 * t53 - t133 * t55) * qJD(3);
t202 = qJDD(3) * t136;
t221 = t130 * pkin(3) - t113;
t79 = qJ(2) + t221;
t262 = (qJD(1) * t79 + t54) * qJD(3) + t202;
t244 = qJ(6) * t67;
t1 = qJD(6) * t97 - t39 * t209 + t200 - t244;
t12 = -t129 * t39 + t132 * t36;
t223 = qJD(6) - t12;
t6 = -pkin(5) * t97 + t223;
t13 = t129 * t36 + t132 * t39;
t7 = qJ(6) * t97 + t13;
t171 = t129 * t6 + t132 * t7;
t187 = t129 * t14 - t132 * t20 + t39 * t208 + t36 * t209;
t256 = pkin(5) * t67;
t2 = qJDD(6) + t187 + t256;
t261 = qJD(5) * t171 + t1 * t129 - t2 * t132;
t155 = t97 * t208 - t243;
t195 = t97 * t213;
t260 = t130 * (t26 + t195) - t133 * (t155 + t233) - t97 * t215;
t259 = t70 ^ 2;
t255 = t7 * t97;
t253 = g(3) * t133;
t252 = t13 * t97;
t250 = t70 * t68;
t74 = pkin(3) * t214 + qJ(4) * t216;
t47 = pkin(8) * t214 + t74;
t48 = -pkin(4) * t216 + t76;
t247 = t129 * t48 + t132 * t47;
t166 = pkin(5) * t132 + qJ(6) * t129;
t158 = -pkin(4) - t166;
t150 = t133 * t158;
t246 = -qJD(1) * t150 + t166 * qJD(5) - qJD(6) * t132 + t268;
t177 = -t221 - t254;
t63 = qJ(2) - t177;
t81 = t248 * t133;
t245 = t129 * t81 + t132 * t63;
t241 = t130 * t70;
t240 = t132 * t26;
t239 = t135 * t67;
t238 = pkin(1) * qJDD(1);
t41 = t125 + t48;
t15 = pkin(5) * t68 - qJ(6) * t70 + t41;
t235 = qJD(3) * t15;
t229 = t130 * t131;
t228 = t130 * t134;
t226 = t132 * t133;
t220 = t134 * pkin(1) + t131 * qJ(2);
t127 = t130 ^ 2;
t128 = t133 ^ 2;
t218 = t127 - t128;
t137 = qJD(3) ^ 2;
t217 = t137 + t138;
t207 = qJD(5) * t135;
t204 = qJDD(3) * t130;
t203 = qJDD(3) * t133;
t199 = g(1) * t229;
t197 = g(1) * (pkin(3) * t227 + qJ(4) * t229);
t193 = t133 * t138 * t130;
t192 = g(1) * t227 - g(2) * t225 - t118;
t190 = pkin(3) * t210 + qJ(4) * t212 + qJD(2);
t182 = (-t127 - t128) * qJDD(1);
t181 = pkin(3) * t229 + t134 * pkin(7) + t220;
t180 = qJDD(2) - t238;
t57 = t129 * t131 - t132 * t225;
t59 = t129 * t134 + t131 * t226;
t179 = g(1) * t57 - g(2) * t59;
t178 = g(1) * t58 - g(2) * t60;
t175 = g(1) * t134 + g(2) * t131;
t172 = t129 * t7 - t132 * t6;
t167 = -pkin(3) * t133 - qJ(4) * t130;
t165 = pkin(5) * t129 - qJ(6) * t132;
t114 = t134 * qJ(2);
t157 = pkin(3) * t228 - qJ(4) * t225 + t114;
t156 = -t97 * t207 - t253;
t78 = qJ(4) + t165;
t38 = t184 * t133 + t190;
t65 = t248 * t212;
t154 = -t129 * t65 + t132 * t38 + t81 * t208 - t63 * t209;
t21 = -t276 * pkin(4) - t29;
t3 = pkin(5) * t27 + qJ(6) * t26 - qJD(6) * t70 + t21;
t153 = -t156 - t3;
t152 = 0.2e1 * qJ(2) * t205 + t202;
t149 = -t136 * t137 - t175;
t148 = t54 * t214 + t183 + t192;
t147 = t97 * t15 - t239;
t144 = g(1) * t59 + g(2) * t57 - t132 * t118 - t187;
t142 = t149 + t265;
t22 = (-qJD(1) * qJD(4) - t206) * t133 + t168;
t45 = -qJD(4) * t133 + t190;
t141 = -qJD(1) * t45 - qJDD(1) * t79 - t149 - t22;
t140 = t15 * t70 + qJDD(6) - t144;
t139 = t67 * t226 + t130 * t27 + t68 * t210 + (t130 * t211 + (qJD(5) * t133 + qJD(1)) * t129) * t97;
t106 = t130 * t136;
t95 = t136 * t210;
t82 = g(2) * t129 * t228;
t80 = -pkin(4) * t130 + t106;
t66 = -pkin(4) * t210 + t95;
t62 = -t217 * t130 + t203;
t61 = t217 * t133 + t204;
t37 = t158 * t130 + t106;
t31 = pkin(5) * t70 + qJ(6) * t68;
t25 = -pkin(5) * t133 + t129 * t63 - t132 * t81;
t24 = qJ(6) * t133 + t245;
t19 = pkin(5) * t216 + t129 * t47 - t132 * t48;
t18 = -qJ(6) * t216 + t247;
t11 = t95 + (t165 * qJD(5) - qJD(6) * t129) * t130 + qJD(3) * t150;
t5 = pkin(5) * t212 + t245 * qJD(5) + t129 * t38 + t132 * t65;
t4 = -qJ(6) * t212 + qJD(6) * t133 + t154;
t8 = [qJDD(1), t267, t175, qJDD(2) - t267 - 0.2e1 * t238, -t175 + t265, -t180 * pkin(1) - g(1) * (-pkin(1) * t131 + t114) - g(2) * t220 + (t198 + t122) * qJ(2), qJDD(1) * t128 - 0.2e1 * t130 * t188, -0.2e1 * t130 * t111 + 0.2e1 * t218 * t205, -t130 * t137 + t203, -t133 * t137 - t204, 0, t130 * t142 + t133 * t152, -t130 * t152 + t133 * t142, t136 * t182 - t143 + t267, t141 * t130 - t133 * t262, t130 * t262 + t141 * t133, t22 * t79 + t54 * t45 - g(1) * (t136 * t131 + t157) - g(2) * (-t131 * t113 + t181) + t143 * t136, t208 * t241 + (-t130 * t26 + t70 * t210) * t129 (-t129 * t68 + t132 * t70) * t210 + (-t129 * t27 - t240 + (-t129 * t70 - t132 * t68) * qJD(5)) * t130 (-t26 + t195) * t133 + (t155 - t233) * t130 (t97 * t211 - t27) * t133 + (t234 + t273) * t130, -t133 * t67 - t97 * t212, -t187 * t133 + t66 * t68 + t80 * t27 + ((-qJD(5) * t81 - t38) * t97 + t63 * t67) * t129 + ((-qJD(5) * t63 - t65) * t97 - t81 * t67 - t41 * t210) * t132 + (-qJD(3) * t12 - t21 * t132 + t41 * t209) * t130 + t178, -t154 * t97 + t245 * t67 + t66 * t70 - t80 * t26 + ((qJD(3) * t41 + t232) * t129 - t200) * t133 + (qJD(3) * t13 + t21 * t129 + t41 * t208) * t130 - t179, t11 * t68 + t25 * t67 + t27 * t37 - t5 * t97 + (-t15 * t211 - t2) * t133 + (qJD(3) * t6 - t132 * t3 + t15 * t209) * t130 + t178, -t24 * t27 - t25 * t26 - t4 * t68 + t5 * t70 + t171 * t210 + (-qJD(5) * t172 + t1 * t132 + t129 * t2 - t175) * t130, -t11 * t70 - t24 * t67 + t26 * t37 + t4 * t97 + (-t15 * t213 + t1) * t133 + (-qJD(3) * t7 - t129 * t3 - t15 * t208) * t130 + t179, t1 * t24 + t7 * t4 + t3 * t37 + t15 * t11 + t2 * t25 + t6 * t5 - g(1) * (-pkin(5) * t58 + pkin(8) * t228 - qJ(6) * t57 + t157) - g(2) * (pkin(4) * t134 + pkin(5) * t60 + qJ(6) * t59 + t181) + (g(1) * t248 - g(2) * t164) * t131; 0, 0, 0, qJDD(1), -t138, t180 - t274, 0, 0, 0, 0, 0, t62, -t61, t182, -t62, t61, t143 + t275, 0, 0, 0, 0, 0, t139, -t260, t139 (-t70 * t212 + qJD(1) * t68 + (qJD(5) * t68 - t26) * t133) * t132 + (-t68 * t212 - qJD(1) * t70 + (t27 - t231) * t133) * t129, t260, -t171 * qJD(1) + (qJD(3) * t172 + t3) * t130 + (t235 - t261) * t133 - t267; 0, 0, 0, 0, 0, 0, t193, -t218 * t138, t111, -t201, qJDD(3) (t86 - t237) * t133 - t192, t274 * t130 + t253 - t75, t167 * qJDD(1) + ((-t55 - t125) * t133 + (t185 + t53) * t130) * qJD(1), t74 * t216 + t148 - 0.2e1 * t230, 0.2e1 * t121 + 0.2e1 * t123 + t75 + (qJD(1) * t74 - g(3)) * t133 + t275 * t130, -t30 * pkin(3) + g(3) * t221 - t29 * qJ(4) - t167 * t119 - t268 * t55 - t53 * t76 - t54 * t74 - t197, -t242 * t70 - t240 (-t27 - t249) * t132 + (t26 + t251) * t129 (-t133 * t242 + t241) * qJD(1) + t273 (-t130 * t68 - t97 * t226) * qJD(1) - t155, t97 * t216, t12 * t216 + qJ(4) * t27 + t82 + t224 * t68 + (-t239 + (t41 - t48) * t97) * t132 + (-t199 - t253 + t21 + (t47 - t207) * t97) * t129, -qJ(4) * t26 + t247 * t97 - t13 * t216 + t224 * t70 + (-t97 * t41 + t239) * t129 + (t156 + t21 - t269) * t132, -t6 * t216 + t19 * t97 + t27 * t78 + t82 + t246 * t68 + t147 * t132 + (-t153 - t199) * t129, t18 * t68 - t19 * t70 + (-t7 * t214 + t135 * t26 + t2 + (-t135 * t68 - t7) * qJD(5)) * t132 + (-t6 * t214 - t135 * t27 - t1 + (t135 * t70 - t6) * qJD(5)) * t129 - t192, t7 * t216 - t18 * t97 + t26 * t78 - t246 * t70 + t147 * t129 + (t153 + t269) * t132, t3 * t78 - t7 * t18 - t6 * t19 - t197 - g(3) * (t133 * t165 + t177) + t246 * t15 - (pkin(8) * t133 + t130 * t165) * t120 + t261 * t135 - (-t78 * t130 + t135 * t133) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, qJDD(3) - t193, -t128 * t138 - t137, qJD(3) * t55 + t148 - t230 + t72, 0, 0, 0, 0, 0, t270, -t263, t270, t129 * t271 + t132 * t272, t263, -t235 + (-t2 + t255) * t132 + (t6 * t97 + t1) * t129 + t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t250, -t68 ^ 2 + t259, -t272, t271, -t67, -t41 * t70 + t144 + t252, t12 * t97 + t41 * t68 - t264, -t31 * t68 - t140 + t252 - 0.2e1 * t256, pkin(5) * t26 - qJ(6) * t27 + (-t13 + t7) * t70 + (t6 - t223) * t68, -0.2e1 * t244 - t15 * t68 + t31 * t70 + (0.2e1 * qJD(6) - t12) * t97 + t264, t1 * qJ(6) - t2 * pkin(5) - t15 * t31 - t6 * t13 - g(1) * (-pkin(5) * t59 + qJ(6) * t60) - g(2) * (-pkin(5) * t57 + qJ(6) * t58) + t223 * t7 - t166 * t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 + t250, -t272, -t258 - t259, t140 - t255 + t256;];
tau_reg  = t8;
