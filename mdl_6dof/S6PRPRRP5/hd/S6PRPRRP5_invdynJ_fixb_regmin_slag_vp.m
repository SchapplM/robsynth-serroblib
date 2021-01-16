% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRRP5
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
% Datum: 2021-01-16 01:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:51:50
% EndTime: 2021-01-16 01:52:05
% DurationCPUTime: 4.13s
% Computational Cost: add. (2966->445), mult. (6159->600), div. (0->0), fcn. (4531->10), ass. (0->235)
t135 = sin(qJ(5));
t138 = cos(qJ(5));
t136 = sin(qJ(4));
t216 = t138 * qJD(4);
t191 = t136 * t216;
t139 = cos(qJ(4));
t218 = qJD(5) * t139;
t193 = t135 * t218;
t208 = t139 * qJDD(2);
t26 = (t191 + t193) * qJD(2) - qJD(5) * t216 - t135 * qJDD(4) - t138 * t208;
t226 = qJD(2) * t136;
t119 = qJD(5) + t226;
t225 = qJD(2) * t139;
t95 = t135 * t225 - t216;
t271 = t119 * t95;
t294 = -t26 - t271;
t131 = sin(pkin(6));
t140 = cos(qJ(2));
t137 = sin(qJ(2));
t241 = t136 * t137;
t62 = (-t135 * t241 + t138 * t140) * t131;
t175 = pkin(4) * t139 + pkin(9) * t136;
t92 = t175 * qJD(4) + qJD(3);
t293 = -qJD(1) * t62 + t138 * t92;
t101 = pkin(4) * t136 - pkin(9) * t139 + qJ(3);
t141 = pkin(2) + pkin(8);
t221 = qJD(4) * t141;
t195 = t139 * t221;
t219 = qJD(5) * t138;
t238 = t137 * t138;
t63 = (t135 * t140 + t136 * t238) * t131;
t292 = qJD(1) * t63 - t101 * t219 - t135 * t92 + t138 * t195;
t215 = qJD(2) * qJD(4);
t291 = -t136 * t215 + t208;
t184 = t135 * t141 + pkin(5);
t217 = qJD(6) * t138;
t220 = qJD(5) * t135;
t240 = t136 * t138;
t111 = t141 * t240;
t256 = t135 * t101 - t111;
t276 = qJ(6) * t191 - t256 * qJD(5) + (qJ(6) * t220 + t184 * qJD(4) - t217) * t139 + t293;
t192 = t138 * t218;
t290 = -qJ(6) * t192 + (-qJD(6) * t139 + (qJ(6) * qJD(4) + qJD(5) * t141) * t136) * t135 - t292;
t224 = qJD(4) * t135;
t97 = t138 * t225 + t224;
t252 = qJD(5) * t97;
t27 = -t138 * qJDD(4) + t291 * t135 + t252;
t270 = t119 * t97;
t289 = t27 + t270;
t133 = cos(pkin(6));
t230 = qJD(1) * t133;
t231 = qJD(1) * t131;
t199 = t140 * t231;
t172 = qJD(3) - t199;
t84 = -t141 * qJD(2) + t172;
t288 = -t136 * t230 + t139 * t84;
t130 = sin(pkin(10));
t132 = cos(pkin(10));
t244 = t133 * t137;
t161 = t130 * t244 - t132 * t140;
t79 = t130 * t140 + t132 * t244;
t176 = -g(1) * t161 + g(2) * t79;
t234 = qJDD(1) - g(3);
t248 = t131 * t137;
t287 = -t234 * t248 + t176;
t200 = t137 * t231;
t229 = qJD(2) * qJ(3);
t100 = t200 + t229;
t286 = (-t100 + t200 - t229) * qJD(4) + qJDD(4) * t141;
t262 = t138 * t161;
t239 = t136 * t140;
t83 = -t131 * t239 + t133 * t139;
t45 = t131 * t238 - t135 * t83;
t281 = g(3) * t45;
t247 = t131 * t139;
t243 = t133 * t140;
t80 = t130 * t243 + t132 * t137;
t40 = t130 * t247 + t136 * t80;
t78 = t130 * t137 - t132 * t243;
t42 = t132 * t247 - t78 * t136;
t69 = t79 * t138;
t285 = -g(1) * (-t135 * t40 - t262) - t281 - g(2) * (t135 * t42 + t69);
t284 = t97 ^ 2;
t283 = pkin(5) * t95;
t282 = g(1) * t80;
t46 = t135 * t248 + t138 * t83;
t280 = g(3) * t46;
t236 = t139 * t140;
t82 = t131 * t236 + t133 * t136;
t279 = g(3) * t82;
t189 = t139 * t215;
t209 = t136 * qJDD(2);
t91 = qJDD(5) + t189 + t209;
t278 = t91 * pkin(5);
t115 = t139 * t230;
t49 = t136 * t84 + t115;
t37 = qJD(4) * pkin(9) + t49;
t58 = t101 * qJD(2) + t200;
t19 = -t135 * t37 + t138 * t58;
t11 = -qJ(6) * t97 + t19;
t8 = pkin(5) * t119 + t11;
t277 = t11 - t8;
t134 = -qJ(6) - pkin(9);
t99 = t175 * qJD(2);
t274 = t135 * t99 + t138 * t288;
t182 = qJD(5) * t134;
t196 = t135 * t226;
t273 = -qJ(6) * t196 + t135 * t182 + t217 - t274;
t86 = t138 * t99;
t272 = t138 * t182 - t86 - (pkin(5) * t139 + qJ(6) * t240) * qJD(2) + (-qJD(6) + t288) * t135;
t20 = t135 * t58 + t138 * t37;
t12 = -qJ(6) * t95 + t20;
t269 = t12 * t119;
t268 = t135 * t26;
t266 = t135 * t79;
t265 = t135 * t161;
t264 = t135 * t91;
t263 = t138 * t80;
t261 = t138 * t91;
t260 = t138 * t97;
t259 = t139 * t26;
t67 = t78 * t135;
t257 = t80 * t135;
t255 = qJ(6) * t139;
t254 = qJD(4) * t95;
t253 = qJD(4) * t97;
t251 = qJDD(2) * pkin(2);
t250 = t130 * t133;
t249 = t131 * t136;
t246 = t131 * t140;
t245 = t132 * t133;
t242 = t133 * t141;
t237 = t137 * t139;
t143 = qJD(2) ^ 2;
t235 = t140 * t143;
t129 = t139 ^ 2;
t233 = t136 ^ 2 - t129;
t142 = qJD(4) ^ 2;
t232 = -t142 - t143;
t228 = qJD(2) * t100;
t227 = qJD(2) * t131;
t223 = qJD(4) * t136;
t222 = qJD(4) * t139;
t214 = qJDD(1) * t131;
t213 = qJDD(1) * t133;
t212 = qJDD(2) * qJ(3);
t211 = qJDD(4) * t136;
t206 = -qJD(4) * t115 - t136 * t213 - t84 * t223;
t205 = qJ(3) * t245;
t185 = -qJD(6) - t283;
t36 = -qJD(4) * pkin(4) - t288;
t25 = -t185 + t36;
t204 = t25 * t219;
t178 = t139 * t200;
t202 = -g(3) * t63 + t95 * t178;
t201 = -g(3) * t62 + t97 * t178;
t198 = t137 * t227;
t197 = t140 * t227;
t194 = t119 * t220;
t188 = t137 * t214;
t187 = t140 * t214;
t186 = t139 * t213;
t183 = pkin(5) * t135 + t141;
t107 = qJD(1) * t198;
t158 = qJDD(3) + t107 - t187;
t51 = -t141 * qJDD(2) + t158;
t16 = qJDD(4) * pkin(9) + qJD(4) * t288 + t136 * t51 + t186;
t24 = t188 + t101 * qJDD(2) + (t92 + t199) * qJD(2);
t181 = -t135 * t24 - t138 * t16 - t58 * t219 + t37 * t220;
t180 = -t51 + t228;
t179 = qJD(5) * t136 + qJD(2);
t39 = t130 * t249 - t80 * t139;
t41 = t132 * t249 + t139 * t78;
t177 = g(1) * t39 - g(2) * t41;
t174 = t12 * t138 - t135 * t8;
t173 = t12 * t135 + t138 * t8;
t171 = t100 * t140 + t137 * (-qJD(2) * pkin(2) + t172);
t121 = pkin(5) * t138 + pkin(4);
t170 = t136 * t121 + t139 * t134;
t169 = -g(2) * t78 + g(3) * t246 - t282;
t168 = qJDD(2) * t137 + t235;
t166 = t27 * qJ(6) + t181;
t163 = t119 * t219 + t264;
t162 = -t194 + t261;
t160 = t133 * t236 - t249;
t159 = t133 * t239 + t247;
t157 = t177 + t279;
t156 = g(1) * t40 - g(2) * t42 + g(3) * t83;
t17 = -qJDD(4) * pkin(4) - t139 * t51 - t206;
t7 = pkin(5) * t27 + qJDD(6) + t17;
t155 = t157 - t7;
t152 = g(3) * t248 + t176;
t151 = -pkin(9) * t91 + t119 * t36;
t150 = qJD(5) * pkin(9) * t119 + g(1) * (t160 * t130 + t132 * t237) - g(2) * (-t130 * t237 + t160 * t132) + t17;
t149 = -t169 + t187;
t23 = t138 * t24;
t148 = -t20 * qJD(5) - t135 * t16 + t23;
t147 = qJDD(3) - t149;
t145 = t26 * qJ(6) + t148;
t52 = t188 + t212 + (qJD(3) + t199) * qJD(2);
t144 = t172 * qJD(2) + t141 * t142 - t152 + t212 + t52;
t126 = qJDD(4) * t139;
t124 = t132 * qJ(3);
t122 = t130 * qJ(3);
t116 = qJ(3) * t248;
t113 = qJ(3) * t250;
t104 = t134 * t138;
t103 = t134 * t135;
t93 = t183 * t139;
t90 = t95 ^ 2;
t89 = t138 * t101;
t77 = t168 * t131;
t76 = (-qJDD(2) * t140 + t137 * t143) * t131;
t68 = t78 * t138;
t65 = t135 * t279;
t61 = t79 * t136;
t60 = t161 * t136;
t59 = t158 - t251;
t53 = -t136 * t221 + (-t135 * t223 + t192) * pkin(5);
t44 = t83 * qJD(4) - t139 * t198;
t43 = -t82 * qJD(4) + t136 * t198;
t38 = -t135 * t255 + t256;
t33 = -t130 * t241 + t159 * t132;
t31 = t159 * t130 + t132 * t241;
t30 = t184 * t136 - t138 * t255 + t89;
t28 = -pkin(5) * t196 + t49;
t14 = t45 * qJD(5) + t135 * t197 + t138 * t43;
t13 = -t46 * qJD(5) - t135 * t43 + t138 * t197;
t6 = -t139 * t27 + (t254 - t264) * t136 + (-t135 * t222 - t179 * t138) * t119;
t5 = t259 + (t253 - t261) * t136 + (t179 * t135 - t139 * t216) * t119;
t4 = -t119 * t14 - t26 * t82 + t44 * t97 - t46 * t91;
t3 = t119 * t13 + t27 * t82 + t44 * t95 + t45 * t91;
t2 = -t95 * qJD(6) - t166;
t1 = -t97 * qJD(6) + t145 + t278;
t9 = [t234, 0, -t76, -t77, t76, t77, qJDD(1) * t133 ^ 2 - g(3) + (t171 * qJD(2) + t137 * t52 - t140 * t59) * t131, 0, 0, 0, 0, 0, -qJD(4) * t44 - qJDD(4) * t82 + (t168 * t136 + t137 * t189) * t131, -qJD(4) * t43 - qJDD(4) * t83 + (t291 * t137 + t139 * t235) * t131, 0, 0, 0, 0, 0, t3, t4, t3, t4, -t13 * t97 - t14 * t95 + t26 * t45 - t27 * t46, t1 * t45 + t12 * t14 + t13 * t8 + t2 * t46 + t25 * t44 + t7 * t82 - g(3); 0, qJDD(2), t149, t287, t147 - 0.2e1 * t251, 0.2e1 * qJD(2) * qJD(3) + 0.2e1 * t212 - t287, t52 * qJ(3) + t100 * qJD(3) - t59 * pkin(2) - g(1) * (-(pkin(2) * t132 + t113) * t137 + (-pkin(2) * t250 + t124) * t140) - g(2) * (-(pkin(2) * t130 - t205) * t137 + (pkin(2) * t245 + t122) * t140) - g(3) * (pkin(2) * t246 + t116) - t171 * t231, qJDD(2) * t129 - 0.2e1 * t136 * t189, -0.2e1 * t136 * t208 + 0.2e1 * t233 * t215, -t136 * t142 + t126, -t139 * t142 - t211, 0, t144 * t136 - t286 * t139, t286 * t136 + t144 * t139, -t97 * t193 + (-t223 * t97 - t259) * t138, (t135 * t97 + t138 * t95) * t223 + (t268 - t138 * t27 + (t135 * t95 - t260) * qJD(5)) * t139, (-t119 * t216 - t26) * t136 + (t162 + t253) * t139, (t119 * t224 - t27) * t136 + (-t163 - t254) * t139, t119 * t222 + t136 * t91, t89 * t91 + g(1) * t138 * t60 - g(2) * (t138 * t61 - t67) + (-t219 * t37 - t221 * t95 + t23) * t136 + (qJD(5) * t111 + t293) * t119 + (t19 * qJD(4) + t141 * t27 + t219 * t36) * t139 + ((-qJD(5) * t101 + t195) * t119 + t17 * t139 + t282 + (-t36 * qJD(4) - qJD(5) * t58 + t141 * t91 - t16) * t136) * t135 + t202, -t256 * t91 - g(1) * (t135 * t60 - t263) - g(2) * (-t135 * t61 - t68) + t292 * t119 + (-t141 * t194 + (-t36 * t138 - t141 * t97) * qJD(4) + t181) * t136 + (-t20 * qJD(4) + t17 * t138 - t141 * t26 - t220 * t36) * t139 + t201, g(1) * t257 + g(2) * t67 + t93 * t27 + t30 * t91 + t53 * t95 + (-t138 * t176 - t224 * t25 + t1) * t136 + t276 * t119 + (t8 * qJD(4) + t7 * t135 + t204) * t139 + t202, g(1) * t263 + g(2) * t68 - t93 * t26 - t38 * t91 + t53 * t97 + (t135 * t176 - t216 * t25 - t2) * t136 - t290 * t119 + (-t12 * qJD(4) + t7 * t138 - t220 * t25) * t139 + t201, t26 * t30 - t27 * t38 - t276 * t97 - t290 * t95 + t173 * t223 + (-qJD(5) * t174 - t1 * t138 - t135 * t2 + t152) * t139, t2 * t38 + t1 * t30 + t7 * t93 - g(1) * (-pkin(5) * t257 - (t132 * t141 + t113) * t137 + (-t130 * t242 + t124) * t140 - t170 * t161) - g(2) * (-pkin(5) * t67 - (t130 * t141 - t205) * t137 + (t132 * t242 + t122) * t140 + t170 * t79) - g(3) * (t116 + (t137 * t170 + t140 * t183) * t131) + t276 * t8 + (t53 + t178) * t25 + t290 * t12; 0, 0, 0, 0, qJDD(2), -t143, t107 + t147 - t228 - t251, 0, 0, 0, 0, 0, t232 * t136 + t126, t232 * t139 - t211, 0, 0, 0, 0, 0, t6, t5, t6, t5, (-t95 * t222 + qJD(2) * t97 + (-t27 + t252) * t136) * t138 + (t97 * t222 + qJD(2) * t95 + (qJD(5) * t95 - t26) * t136) * t135, -t173 * qJD(2) + (qJD(4) * t174 - t7) * t139 + (qJD(4) * t25 - qJD(5) * t173 - t1 * t135 + t138 * t2) * t136 + t169; 0, 0, 0, 0, 0, 0, 0, t139 * t143 * t136, -t233 * t143, t208, -t209, qJDD(4), qJD(4) * t49 - t180 * t139 + t157 + t206, t180 * t136 + t156 - t186, t119 * t260 - t268, -t289 * t135 + t294 * t138, (t119 * t240 - t139 * t97) * qJD(2) + t163, (-t119 * t135 * t136 + t139 * t95) * qJD(2) + t162, -t119 * t225, -t19 * t225 - pkin(4) * t27 - t86 * t119 - t49 * t95 + (t119 * t288 + t151) * t135 + (-t150 + t279) * t138, pkin(4) * t26 + t274 * t119 + t150 * t135 + t151 * t138 + t20 * t225 - t49 * t97 - t65, -t8 * t225 + t103 * t91 - t121 * t27 - t28 * t95 + t272 * t119 + (t25 * t226 + (t25 + t283) * qJD(5)) * t135 + t155 * t138, t204 + t104 * t91 + t121 * t26 - t28 * t97 - t65 - t273 * t119 + (t12 * t139 + t240 * t25) * qJD(2) + (pkin(5) * t252 - t177 + t7) * t135, t103 * t26 + t104 * t27 - t272 * t97 - t273 * t95 + (-t119 * t8 + t2) * t138 + (-t1 - t269) * t135 - t156, -t2 * t104 + t1 * t103 - t7 * t121 - g(1) * (-t121 * t39 - t134 * t40) - g(2) * (t121 * t41 + t134 * t42) - g(3) * (-t121 * t82 - t134 * t83) + t272 * t8 + (pkin(5) * t220 - t28) * t25 + t273 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97 * t95, -t90 + t284, -t26 + t271, t270 - t27, t91, t20 * t119 - t36 * t97 - g(1) * (-t135 * t31 - t262) - g(2) * (t135 * t33 + t69) - t281 + t148, t19 * t119 + t36 * t95 - g(1) * (-t138 * t31 + t265) - g(2) * (t138 * t33 - t266) + t280 + t181, 0.2e1 * t278 + t269 + (t185 - t25) * t97 + t145 + t285, t11 * t119 - t284 * pkin(5) - g(1) * (-t138 * t40 + t265) - g(2) * (t138 * t42 - t266) + t280 + (qJD(6) + t25) * t95 + t166, pkin(5) * t26 + t277 * t95, -t277 * t12 + (-t25 * t97 + t1 + t285) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t289, t294, -t90 - t284, t12 * t95 + t8 * t97 - t155;];
tau_reg = t9;
