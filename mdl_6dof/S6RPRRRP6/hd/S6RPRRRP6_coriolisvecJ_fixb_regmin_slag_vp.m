% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:16:44
% EndTime: 2019-03-09 06:16:52
% DurationCPUTime: 3.43s
% Computational Cost: add. (6834->358), mult. (17553->475), div. (0->0), fcn. (13591->8), ass. (0->187)
t185 = sin(pkin(10));
t186 = cos(pkin(10));
t189 = sin(qJ(3));
t191 = cos(qJ(3));
t153 = t191 * t185 + t189 * t186;
t188 = sin(qJ(4));
t275 = pkin(8) + pkin(9);
t229 = qJD(4) * t275;
t241 = t191 * t186;
t173 = qJD(1) * t241;
t247 = t189 * t185;
t224 = qJD(1) * t247;
t141 = t173 - t224;
t254 = t141 * t188;
t142 = t153 * qJD(1);
t104 = t142 * pkin(3) - t141 * pkin(8);
t190 = cos(qJ(4));
t272 = pkin(7) + qJ(2);
t164 = t272 * t185;
t154 = qJD(1) * t164;
t165 = t272 * t186;
t155 = qJD(1) * t165;
t286 = -t191 * t154 - t189 * t155;
t262 = t188 * t104 + t190 * t286;
t300 = -pkin(9) * t254 + t188 * t229 + t262;
t94 = t190 * t104;
t299 = t142 * pkin(4) - t188 * t286 + t94 + (-pkin(9) * t141 + t229) * t190;
t119 = t190 * qJD(3) - t188 * t142;
t120 = t188 * qJD(3) + t190 * t142;
t187 = sin(qJ(5));
t273 = cos(qJ(5));
t206 = t187 * t119 + t273 * t120;
t276 = t206 ^ 2;
t64 = -t273 * t119 + t187 * t120;
t61 = t64 ^ 2;
t298 = -t61 + t276;
t236 = qJD(4) * t188;
t297 = t236 - t254;
t296 = t206 * t64;
t295 = t64 * qJ(6);
t227 = t273 * t190;
t250 = t187 * t188;
t205 = t227 - t250;
t281 = qJD(4) + qJD(5);
t222 = t273 * qJD(5);
t283 = t273 * qJD(4) + t222;
t260 = t205 * t141 - t283 * t190 + t281 * t250;
t228 = t273 * t188;
t157 = t187 * t190 + t228;
t113 = t281 * t157;
t259 = -t157 * t141 + t113;
t294 = t153 * qJD(2);
t170 = qJD(3) * t173;
t208 = qJD(3) * t224 - t170;
t293 = qJD(3) * qJD(4) - t208;
t235 = qJD(4) * t190;
t225 = t153 * t235;
t152 = -t241 + t247;
t143 = t152 * qJD(3);
t248 = t188 * t143;
t292 = t225 - t248;
t136 = qJD(4) - t141;
t131 = qJD(5) + t136;
t230 = t142 * t235 + t188 * t293;
t234 = qJD(5) * t187;
t76 = -t142 * t236 + t190 * t293;
t28 = -t119 * t222 + t120 * t234 + t187 * t230 - t273 * t76;
t291 = t64 * t131 - t28;
t144 = t153 * qJD(3);
t132 = qJD(1) * t144;
t108 = -t189 * t154 + t191 * t155;
t103 = qJD(3) * pkin(8) + t108;
t177 = -t186 * pkin(2) - pkin(1);
t163 = t177 * qJD(1) + qJD(2);
t80 = -t141 * pkin(3) - t142 * pkin(8) + t163;
t56 = t190 * t103 + t188 * t80;
t198 = t152 * qJD(2);
t69 = -qJD(1) * t198 + qJD(3) * t286;
t88 = t132 * pkin(3) + t208 * pkin(8);
t85 = t190 * t88;
t197 = -t56 * qJD(4) - t188 * t69 + t85;
t14 = t132 * pkin(4) - t76 * pkin(9) + t197;
t202 = -t103 * t236 + t188 * t88 + t190 * t69 + t80 * t235;
t20 = -t230 * pkin(9) + t202;
t55 = -t188 * t103 + t190 * t80;
t42 = -t120 * pkin(9) + t55;
t32 = t136 * pkin(4) + t42;
t43 = t119 * pkin(9) + t56;
t217 = -t187 * t14 - t273 * t20 - t32 * t222 + t43 * t234;
t102 = -qJD(3) * pkin(3) - t286;
t60 = -t119 * pkin(4) + t102;
t290 = t60 * t64 + t217;
t35 = t64 * pkin(5) + qJD(6) + t60;
t289 = t206 * t35;
t288 = t297 * pkin(4) - t108;
t287 = qJ(6) * t206;
t166 = t275 * t188;
t167 = t275 * t190;
t240 = -t187 * t166 + t273 * t167;
t285 = -t240 * qJD(5) + t300 * t187 - t299 * t273;
t284 = t166 * t222 + t167 * t234 + t299 * t187 + t300 * t273;
t282 = -t191 * t164 - t189 * t165;
t39 = t273 * t43;
t17 = t187 * t32 + t39;
t196 = -t17 * qJD(5) + t273 * t14 - t187 * t20;
t280 = -t60 * t206 + t196;
t29 = t206 * qJD(5) + t187 * t76 + t273 * t230;
t279 = t131 * t206 - t29;
t278 = -t260 * t131 + t157 * t132;
t277 = -t205 * t28 - t206 * t259;
t37 = t187 * t43;
t16 = t273 * t32 - t37;
t6 = t16 - t287;
t5 = t131 * pkin(5) + t6;
t274 = t5 - t6;
t271 = t273 * t42 - t37;
t270 = -qJ(6) * t259 + t205 * qJD(6) - t284;
t269 = -t142 * pkin(5) + qJ(6) * t260 - t157 * qJD(6) + t285;
t118 = -t189 * t164 + t191 * t165;
t252 = t153 * t190;
t106 = t152 * pkin(3) - t153 * pkin(8) + t177;
t99 = t190 * t106;
t51 = t152 * pkin(4) - pkin(9) * t252 - t188 * t118 + t99;
t253 = t153 * t188;
t111 = t190 * t118;
t261 = t188 * t106 + t111;
t57 = -pkin(9) * t253 + t261;
t267 = t187 * t51 + t273 * t57;
t266 = t142 * t64;
t264 = t206 * t142;
t263 = t76 * t188;
t258 = t119 * t136;
t257 = t119 * t142;
t256 = t120 * t136;
t255 = t120 * t142;
t249 = t188 * t132;
t123 = t190 * t132;
t245 = t190 * t143;
t239 = t185 ^ 2 + t186 ^ 2;
t238 = qJD(3) * t189;
t237 = qJD(3) * t191;
t233 = qJD(1) * qJD(2);
t182 = -t190 * pkin(4) - pkin(3);
t226 = t153 * t236;
t221 = -t187 * t42 - t39;
t219 = -t187 * t57 + t273 * t51;
t216 = t239 * qJD(1) ^ 2;
t215 = -t273 * t166 - t187 * t167;
t214 = t136 * t190;
t70 = qJD(1) * t294 - t154 * t238 + t155 * t237;
t82 = -t164 * t238 + t165 * t237 + t294;
t213 = -t157 * t29 + t260 * t64;
t212 = -t131 * t259 + t205 * t132;
t83 = pkin(4) * t253 - t282;
t210 = 0.2e1 * t239 * t233;
t209 = -t297 * t136 + t123;
t59 = pkin(4) * t292 + t82;
t81 = t282 * qJD(3) - t198;
t105 = t144 * pkin(3) + t143 * pkin(8);
t95 = t190 * t105;
t24 = pkin(9) * t245 + t144 * pkin(4) - t188 * t81 + t95 + (-t111 + (pkin(9) * t153 - t106) * t188) * qJD(4);
t201 = t188 * t105 + t106 * t235 - t118 * t236 + t190 * t81;
t26 = -pkin(9) * t292 + t201;
t207 = t187 * t24 + t51 * t222 - t57 * t234 + t273 * t26;
t203 = -t226 - t245;
t200 = -pkin(8) * t132 + t136 * t102;
t48 = t230 * pkin(4) + t70;
t11 = t29 * pkin(5) + t48;
t195 = -t267 * qJD(5) - t187 * t26 + t273 * t24;
t181 = t273 * pkin(4) + pkin(5);
t109 = t132 * t152;
t97 = t205 * t153;
t96 = t157 * t153;
t91 = qJ(6) * t205 + t240;
t90 = -t157 * qJ(6) + t215;
t34 = -t143 * t228 - t187 * t226 - t234 * t253 + (-t143 * t187 + t283 * t153) * t190;
t33 = t113 * t153 + t143 * t227 - t187 * t248;
t21 = -t96 * qJ(6) + t267;
t18 = t152 * pkin(5) - t97 * qJ(6) + t219;
t9 = t271 - t287;
t8 = t221 + t295;
t7 = t17 - t295;
t4 = -t34 * qJ(6) - t96 * qJD(6) + t207;
t3 = t144 * pkin(5) + t33 * qJ(6) - t97 * qJD(6) + t195;
t2 = -qJ(6) * t29 - qJD(6) * t64 - t217;
t1 = t132 * pkin(5) + t28 * qJ(6) - qJD(6) * t206 + t196;
t10 = [0, 0, 0, 0, 0, t210, qJ(2) * t210, -t142 * t143 - t208 * t153, -t153 * t132 - t143 * t141 - t142 * t144 + t208 * t152, -t143 * qJD(3), -t144 * qJD(3), 0, -t82 * qJD(3) + t177 * t132 + t163 * t144, -t81 * qJD(3) - t163 * t143 - t177 * t208, t203 * t120 + t76 * t252 -(t190 * t119 - t120 * t188) * t143 + (-t190 * t230 - t263 + (-t188 * t119 - t120 * t190) * qJD(4)) * t153, t120 * t144 + t153 * t123 + t203 * t136 + t76 * t152, t119 * t144 - t136 * t292 - t230 * t152 - t153 * t249, t136 * t144 + t109 (-t118 * t235 + t95) * t136 + t99 * t132 + (-t103 * t235 + t85) * t152 + t55 * t144 - t82 * t119 - t282 * t230 + t102 * t225 + ((-qJD(4) * t106 - t81) * t136 - t118 * t132 + (-qJD(4) * t80 - t69) * t152 + t70 * t153 - t102 * t143) * t188, t203 * t102 + t82 * t120 - t261 * t132 - t201 * t136 - t56 * t144 - t202 * t152 + t70 * t252 - t282 * t76, -t206 * t33 - t28 * t97, -t206 * t34 + t28 * t96 - t97 * t29 + t33 * t64, -t33 * t131 + t97 * t132 + t144 * t206 - t28 * t152, -t34 * t131 - t96 * t132 - t64 * t144 - t29 * t152, t131 * t144 + t109, t131 * t195 + t132 * t219 + t16 * t144 + t152 * t196 + t83 * t29 + t60 * t34 + t48 * t96 + t59 * t64, -t207 * t131 - t267 * t132 - t17 * t144 + t217 * t152 + t206 * t59 - t83 * t28 - t60 * t33 + t48 * t97, -t1 * t97 + t18 * t28 - t2 * t96 - t206 * t3 - t21 * t29 + t5 * t33 - t7 * t34 - t4 * t64, t2 * t21 + t7 * t4 + t1 * t18 + t5 * t3 + t11 * (t96 * pkin(5) + t83) + t35 * (t34 * pkin(5) + t59); 0, 0, 0, 0, 0, -t216, -qJ(2) * t216, 0, 0, 0, 0, 0, 0.2e1 * t142 * qJD(3), t170 + (t141 - t224) * qJD(3), 0, 0, 0, 0, 0, t209 + t257, -t136 ^ 2 * t190 - t249 - t255, 0, 0, 0, 0, 0, t212 - t266, -t264 - t278, t213 - t277, t1 * t205 - t35 * t142 + t2 * t157 - t259 * t5 - t260 * t7; 0, 0, 0, 0, 0, 0, 0, -t142 * t141, -t141 ^ 2 + t142 ^ 2, t170 + (-t141 - t224) * qJD(3), 0, 0, t108 * qJD(3) - t163 * t142 - t70, -t163 * t141 + t152 * t233, t120 * t214 + t263 (t76 + t258) * t190 + (-t230 - t256) * t188, t136 * t214 + t249 - t255, t209 - t257, -t136 * t142, -pkin(3) * t230 - t70 * t190 - t55 * t142 + t108 * t119 + (-pkin(8) * t235 - t94) * t136 + (t136 * t286 + t200) * t188, -pkin(3) * t76 - t108 * t120 + t56 * t142 + t70 * t188 + (pkin(8) * t236 + t262) * t136 + t200 * t190, -t28 * t157 - t206 * t260, t213 + t277, -t264 + t278, t212 + t266, -t131 * t142, t285 * t131 + t132 * t215 - t16 * t142 + t182 * t29 - t205 * t48 + t259 * t60 + t288 * t64, t284 * t131 - t240 * t132 + t17 * t142 + t48 * t157 - t182 * t28 + t288 * t206 - t260 * t60, -t1 * t157 + t2 * t205 - t206 * t269 - t259 * t7 + t260 * t5 - t270 * t64 + t90 * t28 - t91 * t29, t2 * t91 + t1 * t90 + t11 * (-pkin(5) * t205 + t182) + t270 * t7 + t269 * t5 + (t259 * pkin(5) + t288) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120 * t119, -t119 ^ 2 + t120 ^ 2, t76 - t258, -t230 + t256, t132, -t102 * t120 + t56 * t136 + t197, -t102 * t119 + t55 * t136 - t202, t296, t298, t291, t279, t132, -t221 * t131 + (-t120 * t64 - t131 * t234 + t273 * t132) * pkin(4) + t280, t271 * t131 + (-t120 * t206 - t131 * t222 - t187 * t132) * pkin(4) + t290, t181 * t28 - t5 * t64 + t7 * t206 + t9 * t64 + t8 * t206 + (-t187 * t29 + (t187 * t206 - t273 * t64) * qJD(5)) * pkin(4), -pkin(5) * t289 + t1 * t181 - t5 * t8 - t7 * t9 + (-t35 * t120 + t2 * t187 + (-t187 * t5 + t273 * t7) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t296, t298, t291, t279, t132, t17 * t131 + t280, t16 * t131 + t290, pkin(5) * t28 - t274 * t64, t274 * t7 + (t1 - t289) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61 - t276, t206 * t5 + t7 * t64 + t11;];
tauc_reg  = t10;
