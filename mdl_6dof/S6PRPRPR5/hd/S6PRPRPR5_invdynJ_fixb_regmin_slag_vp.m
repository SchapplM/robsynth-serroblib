% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:45:06
% EndTime: 2019-03-08 19:45:15
% DurationCPUTime: 3.36s
% Computational Cost: add. (2721->381), mult. (6434->492), div. (0->0), fcn. (5375->14), ass. (0->204)
t135 = sin(pkin(11));
t138 = cos(pkin(11));
t141 = sin(qJ(4));
t262 = cos(qJ(4));
t101 = t262 * t135 + t141 * t138;
t279 = t101 * qJD(2);
t284 = qJD(6) + t279;
t140 = sin(qJ(6));
t286 = t140 * t284;
t143 = cos(qJ(6));
t220 = qJD(4) * t141;
t202 = t135 * t220;
t207 = t262 * t138;
t190 = qJD(2) * t207;
t198 = qJDD(2) * t262;
t211 = t138 * qJDD(2);
t208 = qJD(4) * t190 + t135 * t198 + t141 * t211;
t53 = qJD(2) * t202 - t208;
t47 = -qJDD(6) + t53;
t41 = t143 * t47;
t173 = -t284 * t286 - t41;
t222 = qJD(2) * t135;
t204 = t141 * t222;
t92 = -t190 + t204;
t71 = qJD(4) * t140 - t143 * t92;
t287 = t284 * t71;
t134 = pkin(11) + qJ(4);
t129 = sin(t134);
t136 = sin(pkin(10));
t142 = sin(qJ(2));
t144 = cos(qJ(2));
t239 = cos(pkin(10));
t240 = cos(pkin(6));
t185 = t240 * t239;
t88 = t136 * t142 - t144 * t185;
t197 = t136 * t240;
t90 = t239 * t142 + t144 * t197;
t189 = g(1) * t90 + g(2) * t88;
t137 = sin(pkin(6));
t232 = t137 * t144;
t281 = -g(3) * t232 + t189;
t157 = t281 * t129;
t256 = pkin(8) + qJ(3);
t108 = t256 * t135;
t109 = t256 * t138;
t181 = t207 * t232;
t201 = qJD(4) * t262;
t205 = qJD(1) * t232;
t235 = t135 * t141;
t253 = qJD(1) * t181 + (qJD(3) * t135 + qJD(4) * t109) * t141 - t205 * t235 - qJD(3) * t207 + t108 * t201;
t68 = -t141 * t108 + t262 * t109;
t285 = t253 * qJD(4) - qJDD(4) * t68 - t157;
t223 = qJD(1) * t142;
t206 = t137 * t223;
t105 = qJD(2) * qJ(3) + t206;
t194 = qJD(1) * t240;
t119 = t138 * t194;
t250 = pkin(8) * qJD(2);
t63 = t119 + (-t105 - t250) * t135;
t76 = t138 * t105 + t135 * t194;
t64 = t138 * t250 + t76;
t251 = -t141 * t64 + t262 * t63;
t278 = qJD(5) - t251;
t212 = t135 * qJDD(2);
t183 = -t138 * t198 + t141 * t212;
t283 = 0.2e1 * qJD(4) * t279 + t183;
t158 = t101 * t232;
t252 = -qJD(1) * t158 + t101 * qJD(3) + t68 * qJD(4);
t233 = t137 * t142;
t146 = qJD(2) ^ 2;
t164 = (qJDD(2) * t144 - t142 * t146) * t137;
t231 = t140 * t144;
t114 = t137 * t231;
t86 = -t135 * t233 + t240 * t138;
t87 = t240 * t135 + t138 * t233;
t39 = t141 * t87 - t262 * t86;
t282 = t143 * t39 + t114;
t280 = t284 - qJD(6);
t184 = qJD(3) - t205;
t221 = qJD(2) * t142;
t200 = qJD(1) * t221;
t277 = t137 * t200 + qJDD(3);
t23 = t141 * t63 + t262 * t64;
t16 = -qJD(4) * qJ(5) - t23;
t265 = pkin(5) * t92;
t11 = -t16 - t265;
t267 = pkin(4) + pkin(9);
t276 = t267 * t47 + (t11 - t23 + t265) * t284;
t215 = qJDD(1) * t137;
t199 = t144 * t215;
t166 = -t199 + t277;
t238 = qJDD(2) * pkin(2);
t77 = t166 - t238;
t179 = t189 - t77;
t261 = g(3) * t144;
t275 = (t200 - t261) * t137 + t179 + t238;
t40 = t141 * t86 + t262 * t87;
t20 = qJD(2) * t158 + t40 * qJD(4);
t97 = t101 * qJD(4);
t54 = qJD(2) * t97 + t183;
t274 = (-t144 * t54 + t92 * t221) * t137 - qJD(4) * t20 - qJDD(4) * t39;
t19 = -qJD(2) * t181 - t86 * t201 + (qJD(4) * t87 + t222 * t232) * t141;
t273 = (-t144 * t53 - t221 * t279) * t137 - qJD(4) * t19 + qJDD(4) * t40;
t182 = t135 * (-t105 * t135 + t119) - t138 * t76;
t272 = t182 * t144 - (-qJD(2) * pkin(2) + t184) * t142;
t130 = cos(t134);
t67 = t262 * t108 + t141 * t109;
t270 = -t252 * qJD(4) - qJDD(4) * t67 + t130 * t281;
t269 = t92 ^ 2;
t268 = t279 ^ 2;
t266 = pkin(4) * t54;
t100 = -t207 + t235;
t128 = pkin(3) * t138 + pkin(2);
t174 = -qJ(5) * t101 - t128;
t29 = t267 * t100 + t174;
t260 = t29 * t47;
t259 = t71 * t92;
t73 = qJD(4) * t143 + t140 * t92;
t258 = t73 * t92;
t257 = t279 * t92;
t255 = -pkin(5) * t97 - t253;
t96 = -t138 * t201 + t202;
t254 = -t96 * pkin(5) + t252;
t249 = qJ(5) * t92;
t248 = t11 * t100;
t247 = t140 * t47;
t245 = t140 * t97;
t218 = qJD(6) * t143;
t219 = qJD(6) * t140;
t17 = -qJD(4) * t219 + t143 * qJDD(4) + t140 * t54 + t92 * t218;
t244 = t143 * t17;
t242 = t143 * t284;
t192 = qJDD(1) * t240;
t214 = qJDD(2) * qJ(3);
t74 = t142 * t215 + t214 + (qJD(3) + t205) * qJD(2);
t45 = t135 * t192 + t138 * t74;
t237 = qJDD(4) * pkin(4);
t234 = t136 * t137;
t230 = t144 * t146;
t229 = t23 * qJD(4);
t227 = pkin(5) * t279 + t278;
t225 = qJDD(1) - g(3);
t224 = t135 ^ 2 + t138 ^ 2;
t213 = qJDD(4) * qJ(5);
t210 = g(3) * t233;
t209 = t143 * t232;
t203 = t137 * t221;
t196 = t137 * t239;
t117 = t138 * t192;
t35 = t117 + (-pkin(8) * qJDD(2) - t74) * t135;
t36 = pkin(8) * t211 + t45;
t195 = -t141 * t35 - t63 * t201 + t64 * t220 - t262 * t36;
t191 = qJDD(4) * t140 - t143 * t54;
t89 = t136 * t144 + t142 * t185;
t91 = -t142 * t197 + t239 * t144;
t188 = g(1) * t91 + g(2) * t89;
t180 = qJ(5) * t96 - qJD(5) * t101;
t28 = pkin(4) * t97 + t180;
t187 = -t28 + t206;
t10 = -t267 * qJD(4) + t227;
t83 = -t128 * qJD(2) + t184;
t152 = -qJ(5) * t279 + t83;
t14 = t267 * t92 + t152;
t4 = t10 * t143 - t14 * t140;
t5 = t10 * t140 + t14 * t143;
t169 = -t140 * t39 + t209;
t168 = t141 * t36 + t64 * t201 + t63 * t220 - t262 * t35;
t56 = -t129 * t196 + t130 * t89;
t58 = t129 * t234 + t130 * t91;
t81 = t240 * t129 + t130 * t233;
t165 = -g(1) * t58 - g(2) * t56 - g(3) * t81;
t163 = -t242 * t284 + t247;
t162 = qJDD(5) + t168;
t156 = t199 + t281;
t44 = -t135 * t74 + t117;
t155 = -t135 * t44 + t138 * t45 - t188;
t6 = -qJD(4) * qJD(5) + t195 - t213;
t154 = t165 - t195;
t3 = -pkin(5) * t54 - t6;
t37 = t101 * pkin(5) + t67;
t153 = -t3 * t100 - t11 * t97 - t37 * t47 - t188;
t65 = -t128 * qJDD(2) + t166;
t55 = t89 * t129 + t130 * t196;
t57 = t129 * t91 - t130 * t234;
t80 = t129 * t233 - t240 * t130;
t151 = g(1) * t57 + g(2) * t55 + g(3) * t80 - t168;
t150 = t3 + (t284 * t267 + t249) * t284 + t165;
t149 = qJ(5) * t53 + t65;
t26 = pkin(4) * t92 + t152;
t148 = t26 * t279 + qJDD(5) - t151;
t147 = -qJD(5) * t279 + t149;
t84 = qJD(4) * t92;
t48 = pkin(4) * t100 + t174;
t46 = pkin(4) * t279 + t249;
t38 = -t100 * pkin(5) + t68;
t21 = t267 * t97 + t180;
t18 = t73 * qJD(6) + t191;
t15 = -qJD(4) * pkin(4) + t278;
t9 = t147 + t266;
t8 = t267 * t54 + t147;
t7 = t162 - t237;
t2 = -t53 * pkin(5) - t267 * qJDD(4) + t162;
t1 = t143 * t2;
t12 = [t225, 0, t164 (-qJDD(2) * t142 - t230) * t137, t138 * t164, -t135 * t164, t224 * t137 * t230 + (-t135 * t86 + t138 * t87) * qJDD(2), t44 * t86 + t45 * t87 - g(3) + (-t272 * qJD(2) - t144 * t77) * t137, 0, 0, 0, 0, 0, t274, -t273, t19 * t92 + t20 * t279 - t39 * t53 - t40 * t54, -t274, t273, t15 * t20 + t16 * t19 + t39 * t7 - t40 * t6 - g(3) + (-t144 * t9 + t26 * t221) * t137, 0, 0, 0, 0, 0 (qJD(6) * t169 - t140 * t203 + t143 * t20) * t284 - t282 * t47 - t19 * t71 + t40 * t18 -(t282 * qJD(6) + t140 * t20 + t143 * t203) * t284 - t169 * t47 - t19 * t73 + t40 * t17; 0, qJDD(2), t156, -t225 * t233 + t188, t275 * t138, -t275 * t135, -t210 + t155 + (t184 * qJD(2) + t214) * t224, -t182 * qJD(3) + t179 * pkin(2) + t155 * qJ(3) + (-g(3) * (pkin(2) * t144 + qJ(3) * t142) + t272 * qJD(1)) * t137, -t101 * t53 - t279 * t96, t100 * t53 - t101 * t54 - t279 * t97 + t92 * t96, -qJD(4) * t96 + qJDD(4) * t101, -qJD(4) * t97 - qJDD(4) * t100, 0, t100 * t65 - t128 * t54 - t92 * t206 + t83 * t97 + t270, t101 * t65 + t128 * t53 - t206 * t279 - t83 * t96 + t285, t100 * t6 + t101 * t7 - t15 * t96 + t16 * t97 + t252 * t279 + t253 * t92 - t53 * t67 - t54 * t68 - t188 - t210, -t100 * t9 + t187 * t92 - t26 * t97 - t48 * t54 - t270, -t101 * t9 + t187 * t279 + t26 * t96 + t48 * t53 - t285, t26 * t28 + t9 * t48 - t6 * t68 + t7 * t67 - t188 * t256 + t253 * t16 + t252 * t15 + (-g(3) * t256 - qJD(1) * t26) * t233 + (-t261 * t137 + t189) * (pkin(4) * t130 + qJ(5) * t129 + t128) t73 * t245 + (t140 * t17 + t218 * t73) * t100 (-t140 * t71 + t143 * t73) * t97 + (-t140 * t18 + t244 + (-t140 * t73 - t143 * t71) * qJD(6)) * t100, t284 * t245 + t101 * t17 - t73 * t96 + (t218 * t284 - t247) * t100, t97 * t242 - t101 * t18 + t71 * t96 + (-t219 * t284 - t41) * t100, -t101 * t47 - t284 * t96, t1 * t101 + t38 * t18 - t4 * t96 + t255 * t71 + (-t101 * t8 + t129 * t189 - t21 * t284 + t260) * t140 + (t254 * t284 + t153) * t143 + (t223 * t286 - g(3) * (t129 * t231 + t142 * t143)) * t137 + ((-t140 * t37 - t143 * t29) * t284 - t5 * t101 + t140 * t248) * qJD(6), t38 * t17 + t5 * t96 + t255 * t73 + (t260 - (qJD(6) * t10 + t8) * t101 + qJD(6) * t248 + (-qJD(6) * t37 + t206 - t21) * t284 + t157) * t143 + (-(-qJD(6) * t14 + t2) * t101 + t210 + (qJD(6) * t29 - t254) * t284 - t153) * t140; 0, 0, 0, 0, -t211, t212, -t224 * t146, t182 * qJD(2) - t156 - t238 + t277, 0, 0, 0, 0, 0, t283 (-t92 - t204) * qJD(4) + t208, -t268 - t269, -t283, t53 + t84, t266 - t16 * t92 + (-qJD(5) - t15) * t279 + t149 - t281, 0, 0, 0, 0, 0, t163 + t259, -t173 + t258; 0, 0, 0, 0, 0, 0, 0, 0, t257, t268 - t269 (t92 - t204) * qJD(4) + t208, -t183, qJDD(4), -t279 * t83 + t151 + t229, qJD(4) * t251 + t83 * t92 - t154, pkin(4) * t53 - qJ(5) * t54 + (-t16 - t23) * t279 + (t15 - t278) * t92, t46 * t92 + t148 - t229 - 0.2e1 * t237, 0.2e1 * t213 - t26 * t92 + t46 * t279 + (0.2e1 * qJD(5) - t251) * qJD(4) + t154, -t6 * qJ(5) - t7 * pkin(4) - t26 * t46 - t15 * t23 - g(1) * (-pkin(4) * t57 + qJ(5) * t58) - g(2) * (-pkin(4) * t55 + qJ(5) * t56) - g(3) * (-pkin(4) * t80 + qJ(5) * t81) - t278 * t16, -t286 * t73 + t244 (-t284 * t73 - t18) * t143 + (-t17 + t287) * t140, t173 + t258, t163 - t259, t284 * t92, qJ(5) * t18 + t150 * t140 + t143 * t276 + t227 * t71 + t4 * t92, qJ(5) * t17 - t140 * t276 + t150 * t143 + t227 * t73 - t5 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53 + t84, qJDD(4) - t257, -qJD(4) ^ 2 - t268, t16 * qJD(4) + t148 - t237, 0, 0, 0, 0, 0, -qJD(4) * t71 + t173, -qJD(4) * t73 + t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t71, -t71 ^ 2 + t73 ^ 2, t17 + t287, t280 * t73 - t191, -t47, -t140 * t8 + t1 - t11 * t73 - g(1) * (-t140 * t90 + t143 * t57) - g(2) * (-t140 * t88 + t143 * t55) - g(3) * (t143 * t80 + t114) + t280 * t5, -t143 * t8 - t140 * t2 + t11 * t71 - g(1) * (-t140 * t57 - t143 * t90) - g(2) * (-t140 * t55 - t143 * t88) - g(3) * (-t140 * t80 + t209) + t280 * t4;];
tau_reg  = t12;
