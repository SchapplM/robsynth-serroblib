% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:41:20
% EndTime: 2019-03-08 19:41:30
% DurationCPUTime: 4.09s
% Computational Cost: add. (4483->430), mult. (10733->598), div. (0->0), fcn. (9090->18), ass. (0->218)
t184 = cos(pkin(11));
t301 = cos(qJ(4));
t248 = t301 * t184;
t181 = sin(pkin(11));
t187 = sin(qJ(4));
t265 = t187 * t181;
t207 = t248 - t265;
t182 = sin(pkin(6));
t190 = cos(qJ(2));
t266 = t182 * t190;
t197 = t207 * t266;
t296 = pkin(8) + qJ(3);
t149 = t296 * t181;
t151 = t296 * t184;
t208 = -t149 * t301 - t151 * t187;
t314 = qJD(1) * t197 - qJD(3) * t207 - qJD(4) * t208;
t137 = t207 * qJD(4);
t143 = t181 * t301 + t184 * t187;
t138 = t143 * qJD(4);
t188 = sin(qJ(2));
t260 = qJD(1) * t182;
t247 = t188 * t260;
t313 = pkin(4) * t138 - qJ(5) * t137 - qJD(5) * t143 - t247;
t145 = qJD(2) * qJ(3) + t247;
t280 = cos(pkin(6));
t237 = qJD(1) * t280;
t160 = t184 * t237;
t290 = pkin(8) * qJD(2);
t103 = t160 + (-t145 - t290) * t181;
t118 = t145 * t184 + t181 * t237;
t104 = t184 * t290 + t118;
t46 = t103 * t187 + t104 * t301;
t312 = t46 * qJD(4);
t164 = qJD(2) * t248;
t246 = qJD(2) * t265;
t132 = -t164 + t246;
t124 = qJD(6) + t132;
t311 = t124 - qJD(6);
t134 = t143 * qJD(2);
t180 = sin(pkin(12));
t183 = cos(pkin(12));
t112 = -qJD(4) * t183 + t134 * t180;
t114 = qJD(4) * t180 + t134 * t183;
t186 = sin(qJ(6));
t189 = cos(qJ(6));
t48 = t112 * t189 + t114 * t186;
t310 = t124 * t48;
t142 = t180 * t189 + t183 * t186;
t136 = t142 * qJD(6);
t283 = t132 * t142 + t136;
t220 = t112 * t186 - t114 * t189;
t309 = t124 * t220;
t292 = t180 * t314 + t183 * t313;
t291 = t180 * t313 - t183 * t314;
t109 = -t149 * t187 + t151 * t301;
t198 = t143 * t266;
t281 = -qJD(1) * t198 + qJD(3) * t143 + qJD(4) * t109;
t45 = t103 * t301 - t104 * t187;
t179 = pkin(11) + qJ(4);
t173 = sin(t179);
t279 = cos(pkin(10));
t229 = t280 * t279;
t278 = sin(pkin(10));
t128 = t188 * t278 - t190 * t229;
t228 = t280 * t278;
t130 = t188 * t279 + t190 * t228;
t234 = g(1) * t130 + g(2) * t128;
t307 = g(3) * t266 - t234;
t308 = t307 * t173;
t191 = qJD(2) ^ 2;
t202 = (qJDD(2) * t190 - t188 * t191) * t182;
t244 = t190 * t260;
t227 = qJD(3) - t244;
t259 = qJD(2) * t188;
t243 = qJD(1) * t259;
t306 = t182 * t243 + qJDD(3);
t305 = -qJDD(4) * pkin(4) + qJDD(5);
t140 = t180 * t186 - t183 * t189;
t284 = t124 * t140;
t240 = qJDD(2) * t301;
t254 = t181 * qJDD(2);
t226 = -t184 * t240 + t187 * t254;
t92 = qJD(2) * t138 + t226;
t87 = qJDD(6) + t92;
t304 = t124 * t284 - t142 * t87;
t256 = qJDD(1) * t182;
t242 = t190 * t256;
t211 = -t242 + t306;
t277 = qJDD(2) * pkin(2);
t119 = t211 - t277;
t217 = -t119 + t234;
t303 = t182 * (-g(3) * t190 + t243) + t217 + t277;
t219 = (-t145 * t181 + t160) * t181 - t118 * t184;
t302 = t190 * t219 - (-qJD(2) * pkin(2) + t227) * t188;
t253 = t184 * qJDD(2);
t251 = qJD(4) * t164 + t181 * t240 + t187 * t253;
t91 = -qJD(4) * t246 + t251;
t70 = -qJDD(4) * t183 + t180 * t91;
t71 = qJDD(4) * t180 + t183 * t91;
t15 = -qJD(6) * t220 + t186 * t71 + t189 * t70;
t125 = t132 ^ 2;
t300 = pkin(9) * t183;
t297 = g(3) * t182;
t295 = pkin(9) + qJ(5);
t255 = qJDD(2) * qJ(3);
t116 = t188 * t256 + t255 + (qJD(3) + t244) * qJD(2);
t236 = qJDD(1) * t280;
t158 = t184 * t236;
t72 = t158 + (-pkin(8) * qJDD(2) - t116) * t181;
t85 = t116 * t184 + t181 * t236;
t73 = pkin(8) * t253 + t85;
t214 = t187 * t72 + t301 * t73;
t11 = qJDD(4) * qJ(5) + (qJD(5) + t45) * qJD(4) + t214;
t168 = pkin(3) * t184 + pkin(2);
t106 = -qJDD(2) * t168 + t211;
t24 = t92 * pkin(4) - t91 * qJ(5) - t134 * qJD(5) + t106;
t7 = t11 * t183 + t180 * t24;
t294 = pkin(5) * t138 - t137 * t300 + t292;
t274 = t137 * t180;
t293 = pkin(9) * t274 - t291;
t42 = qJD(4) * qJ(5) + t46;
t123 = -qJD(2) * t168 + t227;
t55 = t132 * pkin(4) - t134 * qJ(5) + t123;
t19 = t180 * t55 + t183 * t42;
t86 = pkin(4) * t134 + qJ(5) * t132;
t26 = t180 * t86 + t183 * t45;
t289 = t134 * t48;
t287 = t180 * t92;
t286 = t183 * t92;
t285 = t220 * t134;
t88 = -pkin(4) * t207 - qJ(5) * t143 - t168;
t40 = t109 * t183 + t180 * t88;
t282 = pkin(5) * t274 + t281;
t275 = t132 * t180;
t273 = t143 * t180;
t272 = t143 * t183;
t178 = pkin(12) + qJ(6);
t172 = sin(t178);
t175 = cos(t179);
t270 = t172 * t175;
t174 = cos(t178);
t269 = t174 * t175;
t268 = t175 * t190;
t267 = t182 * t188;
t264 = t190 * t191;
t38 = -qJD(4) * pkin(4) + qJD(5) - t45;
t263 = -qJD(5) + t38;
t262 = qJDD(1) - g(3);
t261 = t181 ^ 2 + t184 ^ 2;
t258 = qJD(6) * t186;
t257 = qJD(6) * t189;
t252 = g(3) * t267;
t6 = -t11 * t180 + t183 * t24;
t2 = pkin(5) * t92 - pkin(9) * t71 + t6;
t5 = -pkin(9) * t70 + t7;
t250 = -t186 * t5 + t189 * t2;
t245 = t182 * t259;
t18 = -t180 * t42 + t183 * t55;
t25 = -t180 * t45 + t183 * t86;
t239 = t182 * t279;
t238 = t182 * t278;
t39 = -t109 * t180 + t183 * t88;
t235 = -t124 * t283 - t140 * t87;
t129 = t188 * t229 + t190 * t278;
t131 = -t188 * t228 + t190 * t279;
t233 = g(1) * t131 + g(2) * t129;
t232 = -t7 * t180 - t6 * t183;
t231 = t186 * t2 + t189 * t5;
t13 = -pkin(9) * t112 + t19;
t9 = pkin(5) * t132 - pkin(9) * t114 + t18;
t4 = t13 * t189 + t186 * t9;
t230 = t13 * t186 - t189 * t9;
t225 = -t18 * t180 + t183 * t19;
t29 = -pkin(5) * t207 - pkin(9) * t272 + t39;
t31 = -pkin(9) * t273 + t40;
t224 = -t186 * t31 + t189 * t29;
t223 = t186 * t29 + t189 * t31;
t126 = -t181 * t267 + t184 * t280;
t127 = t181 * t280 + t184 * t267;
t77 = t126 * t187 + t127 * t301;
t56 = -t180 * t77 - t183 * t266;
t57 = -t180 * t266 + t183 * t77;
t222 = -t186 * t57 + t189 * t56;
t221 = t186 * t56 + t189 * t57;
t213 = t187 * t73 - t301 * t72 + t312;
t209 = t126 * t301 - t127 * t187;
t148 = t295 * t180;
t206 = pkin(9) * t275 - qJD(5) * t183 + qJD(6) * t148 + t26;
t150 = t295 * t183;
t205 = pkin(5) * t134 + qJD(5) * t180 + qJD(6) * t150 + t132 * t300 + t25;
t14 = -t112 * t257 - t114 * t258 - t186 * t70 + t189 * t71;
t204 = g(1) * (t131 * t173 - t175 * t238) + g(2) * (t129 * t173 + t175 * t239) + g(3) * (t173 * t267 - t175 * t280);
t122 = t173 * t280 + t175 * t267;
t95 = t129 * t175 - t173 * t239;
t97 = t131 * t175 + t173 * t238;
t203 = g(1) * t97 + g(2) * t95 + g(3) * t122;
t12 = t213 + t305;
t201 = -t12 + t204;
t196 = t307 * t175;
t195 = -t307 + t242;
t194 = t12 * t143 + t137 * t38 - t233;
t84 = -t116 * t181 + t158;
t193 = -t84 * t181 + t85 * t184 - t233;
t192 = t204 - t213;
t169 = -pkin(5) * t183 - pkin(4);
t81 = t140 * t143;
t80 = t142 * t143;
t67 = pkin(5) * t273 - t208;
t44 = qJD(2) * t198 + qJD(4) * t77;
t43 = qJD(2) * t197 + qJD(4) * t209;
t36 = t180 * t245 + t183 * t43;
t35 = -t180 * t43 + t183 * t245;
t34 = t137 * t142 + t257 * t272 - t258 * t273;
t33 = -t136 * t143 - t137 * t140;
t32 = -pkin(5) * t275 + t46;
t30 = pkin(5) * t112 + t38;
t8 = pkin(5) * t70 + t12;
t1 = [t262, 0, t202 (-qJDD(2) * t188 - t264) * t182, t184 * t202, -t181 * t202, t261 * t182 * t264 + (-t126 * t181 + t127 * t184) * qJDD(2), t84 * t126 + t85 * t127 - g(3) + (-qJD(2) * t302 - t119 * t190) * t182, 0, 0, 0, 0, 0, -t44 * qJD(4) + t209 * qJDD(4) + (t132 * t259 - t190 * t92) * t182, -t43 * qJD(4) - t77 * qJDD(4) + (t134 * t259 - t190 * t91) * t182, t112 * t44 + t132 * t35 - t209 * t70 + t56 * t92, t114 * t44 - t132 * t36 - t209 * t71 - t57 * t92, -t112 * t36 - t114 * t35 - t56 * t71 - t57 * t70, -t12 * t209 + t18 * t35 + t19 * t36 + t38 * t44 + t56 * t6 + t57 * t7 - g(3), 0, 0, 0, 0, 0 (-qJD(6) * t221 - t186 * t36 + t189 * t35) * t124 + t222 * t87 + t44 * t48 - t209 * t15 -(qJD(6) * t222 + t186 * t35 + t189 * t36) * t124 - t221 * t87 - t44 * t220 - t209 * t14; 0, qJDD(2), t195, -t262 * t267 + t233, t303 * t184, -t303 * t181, -t252 + t193 + (qJD(2) * t227 + t255) * t261, -t219 * qJD(3) + t217 * pkin(2) + t193 * qJ(3) + (-g(3) * (pkin(2) * t190 + qJ(3) * t188) + t302 * qJD(1)) * t182, t134 * t137 + t143 * t91, -t132 * t137 - t134 * t138 - t143 * t92 + t207 * t91, qJD(4) * t137 + qJDD(4) * t143, -qJD(4) * t138 + qJDD(4) * t207, 0, -qJD(4) * t281 + qJDD(4) * t208 - t106 * t207 + t123 * t138 - t132 * t247 - t168 * t92 - t196, qJD(4) * t314 - t109 * qJDD(4) + t106 * t143 + t123 * t137 - t134 * t247 - t168 * t91 + t308, -t208 * t70 + t18 * t138 - t6 * t207 + t39 * t92 - t183 * t196 + (t194 - t252) * t180 + t292 * t132 + t281 * t112, -t208 * t71 - t19 * t138 + t7 * t207 - t40 * t92 - t234 * t180 * t175 + t194 * t183 - (-t180 * t268 + t183 * t188) * t297 - t291 * t132 + t281 * t114, -t39 * t71 - t40 * t70 + t232 * t143 + (-t18 * t183 - t180 * t19) * t137 - t292 * t114 - t291 * t112 - t308, -t12 * t208 + t292 * t18 + t291 * t19 + t281 * t38 + t6 * t39 + t7 * t40 - (t188 * t297 + t233) * t296 + (-t190 * t297 + t234) * (pkin(4) * t175 + qJ(5) * t173 + t168) -t14 * t81 - t220 * t33, -t14 * t80 + t15 * t81 + t220 * t34 - t33 * t48, t124 * t33 - t138 * t220 - t14 * t207 - t81 * t87, -t124 * t34 - t138 * t48 + t15 * t207 - t80 * t87, t124 * t138 - t207 * t87, t224 * t87 - t250 * t207 - t230 * t138 + t67 * t15 + t8 * t80 + t30 * t34 - g(1) * (-t130 * t269 + t131 * t172) - g(2) * (-t128 * t269 + t129 * t172) + t282 * t48 - (t172 * t188 + t174 * t268) * t297 + (t186 * t293 + t189 * t294) * t124 + (-t124 * t223 + t207 * t4) * qJD(6), -t223 * t87 + t231 * t207 - t4 * t138 + t67 * t14 - t8 * t81 + t30 * t33 - g(1) * (t130 * t270 + t131 * t174) - g(2) * (t128 * t270 + t129 * t174) - t282 * t220 - (-t172 * t268 + t174 * t188) * t297 + (-t186 * t294 + t189 * t293) * t124 + (-t124 * t224 - t207 * t230) * qJD(6); 0, 0, 0, 0, -t253, t254, -t261 * t191, qJD(2) * t219 - t195 - t277 + t306, 0, 0, 0, 0, 0, 0.2e1 * t134 * qJD(4) + t226 (-t132 - t246) * qJD(4) + t251, -t112 * t134 - t125 * t180 + t286, -t114 * t134 - t125 * t183 - t287, -t180 * t70 - t183 * t71 + (-t112 * t183 + t114 * t180) * t132, t132 * t225 - t38 * t134 - t232 + t307, 0, 0, 0, 0, 0, t235 - t289, t285 + t304; 0, 0, 0, 0, 0, 0, 0, 0, t134 * t132, t134 ^ 2 - t125 (t132 - t246) * qJD(4) + t251, -t226, qJDD(4), -t123 * t134 + t192 + t312, t123 * t132 + t203 - t214, -qJ(5) * t287 - pkin(4) * t70 - t46 * t112 - t18 * t134 + (t180 * t263 - t25) * t132 + t201 * t183, -qJ(5) * t286 - pkin(4) * t71 - t46 * t114 + t19 * t134 + (t183 * t263 + t26) * t132 - t201 * t180, t26 * t112 + t25 * t114 + (-qJ(5) * t70 - qJD(5) * t112 - t132 * t18 + t7) * t183 + (qJ(5) * t71 + qJD(5) * t114 - t132 * t19 - t6) * t180 - t203, -t18 * t25 - t19 * t26 - t38 * t46 + t225 * qJD(5) + t201 * pkin(4) + (-t6 * t180 + t7 * t183 - t203) * qJ(5), t14 * t142 + t220 * t284, -t14 * t140 - t142 * t15 + t220 * t283 + t284 * t48, t285 - t304, t235 + t289, -t124 * t134 (-t148 * t189 - t150 * t186) * t87 + t169 * t15 + t8 * t140 + t230 * t134 - t32 * t48 + t283 * t30 + (t186 * t206 - t189 * t205) * t124 + t204 * t174 -(-t148 * t186 + t150 * t189) * t87 + t169 * t14 + t8 * t142 + t4 * t134 + t32 * t220 - t284 * t30 + (t186 * t205 + t189 * t206) * t124 - t204 * t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114 * t132 + t70, -t112 * t132 + t71, -t112 ^ 2 - t114 ^ 2, t112 * t19 + t114 * t18 - t192 + t305, 0, 0, 0, 0, 0, t15 - t309, t14 - t310; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220 * t48, t220 ^ 2 - t48 ^ 2, t14 + t310, -t15 - t309, t87, t30 * t220 - g(1) * (t130 * t174 - t172 * t97) - g(2) * (t128 * t174 - t172 * t95) - g(3) * (-t122 * t172 - t174 * t266) + t250 + t311 * t4, t30 * t48 - g(1) * (-t130 * t172 - t174 * t97) - g(2) * (-t128 * t172 - t174 * t95) - g(3) * (-t122 * t174 + t172 * t266) - t231 - t311 * t230;];
tau_reg  = t1;
