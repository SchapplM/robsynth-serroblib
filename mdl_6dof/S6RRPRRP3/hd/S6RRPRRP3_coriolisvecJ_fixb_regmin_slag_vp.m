% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:50:34
% EndTime: 2019-03-09 11:50:44
% DurationCPUTime: 3.53s
% Computational Cost: add. (7225->372), mult. (18186->505), div. (0->0), fcn. (13651->8), ass. (0->204)
t193 = sin(qJ(4));
t190 = sin(pkin(10));
t182 = t190 * pkin(2) + pkin(8);
t282 = pkin(9) + t182;
t224 = qJD(4) * t282;
t194 = sin(qJ(2));
t247 = qJD(1) * t194;
t191 = cos(pkin(10));
t196 = cos(qJ(2));
t258 = t191 * t196;
t154 = qJD(1) * t258 - t190 * t247;
t263 = t154 * t193;
t281 = -qJ(3) - pkin(7);
t175 = t281 * t196;
t171 = qJD(1) * t175;
t159 = t190 * t171;
t174 = t281 * t194;
t170 = qJD(1) * t174;
t111 = t191 * t170 + t159;
t195 = cos(qJ(4));
t166 = t190 * t196 + t191 * t194;
t156 = t166 * qJD(1);
t95 = pkin(2) * t247 + t156 * pkin(3) - t154 * pkin(8);
t270 = t195 * t111 + t193 * t95;
t306 = -pkin(9) * t263 + t193 * t224 + t270;
t88 = t195 * t95;
t305 = t156 * pkin(4) - t193 * t111 + t88 + (-pkin(9) * t154 + t224) * t195;
t121 = t195 * qJD(2) - t193 * t156;
t122 = t193 * qJD(2) + t195 * t156;
t192 = sin(qJ(5));
t283 = cos(qJ(5));
t211 = t192 * t121 + t283 * t122;
t285 = t211 ^ 2;
t63 = -t283 * t121 + t192 * t122;
t61 = t63 ^ 2;
t304 = -t61 + t285;
t246 = qJD(4) * t193;
t303 = t246 - t263;
t302 = t211 * t63;
t301 = t63 * qJ(6);
t236 = t283 * t195;
t257 = t192 * t193;
t210 = t236 - t257;
t290 = qJD(4) + qJD(5);
t232 = t283 * qJD(5);
t291 = t283 * qJD(4) + t232;
t269 = t210 * t154 - t291 * t195 + t290 * t257;
t237 = t283 * t193;
t169 = t192 * t195 + t237;
t115 = t290 * t169;
t268 = -t169 * t154 + t115;
t243 = qJD(1) * qJD(2);
t230 = t196 * t243;
t231 = t194 * t243;
t213 = -t190 * t231 + t191 * t230;
t300 = qJD(2) * qJD(4) + t213;
t245 = qJD(4) * t195;
t234 = t166 * t245;
t165 = t190 * t194 - t258;
t158 = t165 * qJD(2);
t255 = t193 * t158;
t299 = t234 - t255;
t148 = qJD(4) - t154;
t141 = qJD(5) + t148;
t239 = t156 * t245 + t193 * t300;
t244 = qJD(5) * t192;
t80 = -t156 * t246 + t195 * t300;
t28 = -t121 * t232 + t122 * t244 + t192 * t239 - t283 * t80;
t298 = t63 * t141 - t28;
t155 = t166 * qJD(2);
t142 = qJD(1) * t155;
t275 = qJD(2) * pkin(2);
t162 = t170 + t275;
t259 = t191 * t171;
t108 = t190 * t162 - t259;
t103 = qJD(2) * pkin(8) + t108;
t238 = -t196 * pkin(2) - pkin(1);
t216 = t238 * qJD(1);
t172 = qJD(3) + t216;
t84 = -t154 * pkin(3) - t156 * pkin(8) + t172;
t55 = t195 * t103 + t193 * t84;
t180 = pkin(2) * t231;
t83 = t142 * pkin(3) - t213 * pkin(8) + t180;
t74 = t195 * t83;
t225 = qJD(2) * t281;
t152 = t196 * qJD(3) + t194 * t225;
t128 = t152 * qJD(1);
t153 = -t194 * qJD(3) + t196 * t225;
t129 = t153 * qJD(1);
t79 = t191 * t128 + t190 * t129;
t203 = -t55 * qJD(4) - t193 * t79 + t74;
t14 = t142 * pkin(4) - t80 * pkin(9) + t203;
t207 = -t103 * t246 + t193 * t83 + t195 * t79 + t84 * t245;
t20 = -t239 * pkin(9) + t207;
t54 = -t193 * t103 + t195 * t84;
t45 = -t122 * pkin(9) + t54;
t32 = t148 * pkin(4) + t45;
t46 = t121 * pkin(9) + t55;
t223 = -t192 * t14 - t283 * t20 - t32 * t232 + t46 * t244;
t107 = t191 * t162 + t159;
t102 = -qJD(2) * pkin(3) - t107;
t60 = -t121 * pkin(4) + t102;
t297 = t60 * t63 + t223;
t296 = -0.2e1 * t243;
t37 = t63 * pkin(5) + qJD(6) + t60;
t295 = t211 * t37;
t110 = t190 * t170 - t259;
t218 = t303 * pkin(4) - t110;
t294 = qJ(6) * t211;
t163 = t282 * t193;
t164 = t282 * t195;
t249 = -t192 * t163 + t283 * t164;
t293 = -t249 * qJD(5) + t306 * t192 - t305 * t283;
t292 = t163 * t232 + t164 * t244 + t305 * t192 + t306 * t283;
t43 = t283 * t46;
t17 = t192 * t32 + t43;
t202 = -t17 * qJD(5) + t283 * t14 - t192 * t20;
t289 = -t60 * t211 + t202;
t29 = t211 * qJD(5) + t192 * t80 + t283 * t239;
t288 = t141 * t211 - t29;
t287 = -t269 * t141 + t169 * t142;
t286 = -t210 * t28 - t211 * t268;
t41 = t192 * t46;
t16 = t283 * t32 - t41;
t6 = t16 - t294;
t5 = t141 * pkin(5) + t6;
t284 = t5 - t6;
t280 = -qJ(6) * t268 + t210 * qJD(6) - t292;
t279 = -t156 * pkin(5) + qJ(6) * t269 - t169 * qJD(6) + t293;
t278 = t283 * t45 - t41;
t106 = t165 * pkin(3) - t166 * pkin(8) + t238;
t101 = t195 * t106;
t120 = t190 * t174 - t191 * t175;
t261 = t166 * t195;
t51 = t165 * pkin(4) - pkin(9) * t261 - t193 * t120 + t101;
t113 = t195 * t120;
t250 = t193 * t106 + t113;
t262 = t166 * t193;
t58 = -pkin(9) * t262 + t250;
t276 = t192 * t51 + t283 * t58;
t274 = t156 * t63;
t272 = t211 * t156;
t271 = t80 * t193;
t267 = t121 * t148;
t266 = t121 * t156;
t265 = t122 * t148;
t264 = t122 * t156;
t256 = t193 * t142;
t125 = t195 * t142;
t254 = t195 * t158;
t198 = qJD(1) ^ 2;
t253 = t196 * t198;
t197 = qJD(2) ^ 2;
t252 = t197 * t194;
t251 = t197 * t196;
t248 = t194 ^ 2 - t196 ^ 2;
t241 = t194 * t275;
t183 = -t191 * pkin(2) - pkin(3);
t235 = t166 * t246;
t228 = -t192 * t45 - t43;
t227 = -t192 * t58 + t283 * t51;
t222 = pkin(1) * t296;
t78 = t190 * t128 - t191 * t129;
t93 = t190 * t152 - t191 * t153;
t221 = -t283 * t163 - t192 * t164;
t119 = -t191 * t174 - t190 * t175;
t220 = t148 * t195;
t219 = -t169 * t29 + t269 * t63;
t217 = -t141 * t268 + t210 * t142;
t92 = pkin(4) * t262 + t119;
t215 = -t120 * t142 + t78 * t166;
t173 = -t195 * pkin(4) + t183;
t59 = pkin(4) * t299 + t93;
t214 = -t303 * t148 + t125;
t96 = t155 * pkin(3) + t158 * pkin(8) + t241;
t89 = t195 * t96;
t94 = t191 * t152 + t190 * t153;
t24 = pkin(9) * t254 + t155 * pkin(4) - t193 * t94 + t89 + (-t113 + (pkin(9) * t166 - t106) * t193) * qJD(4);
t206 = t106 * t245 - t120 * t246 + t193 * t96 + t195 * t94;
t26 = -pkin(9) * t299 + t206;
t212 = t192 * t24 + t51 * t232 - t58 * t244 + t283 * t26;
t208 = -t235 - t254;
t204 = t148 * t102 - t182 * t142;
t52 = t239 * pkin(4) + t78;
t15 = t29 * pkin(5) + t52;
t201 = -t276 * qJD(5) - t192 * t26 + t283 * t24;
t187 = t283 * pkin(4) + pkin(5);
t109 = t142 * t165;
t99 = t210 * t166;
t98 = t169 * t166;
t77 = qJ(6) * t210 + t249;
t76 = -t169 * qJ(6) + t221;
t34 = -t158 * t237 - t192 * t235 - t244 * t262 + (-t158 * t192 + t291 * t166) * t195;
t33 = t115 * t166 + t158 * t236 - t192 * t255;
t21 = -t98 * qJ(6) + t276;
t19 = t165 * pkin(5) - t99 * qJ(6) + t227;
t9 = t278 - t294;
t8 = t228 + t301;
t7 = t17 - t301;
t4 = -t34 * qJ(6) - t98 * qJD(6) + t212;
t3 = t155 * pkin(5) + t33 * qJ(6) - t99 * qJD(6) + t201;
t2 = -t29 * qJ(6) - t63 * qJD(6) - t223;
t1 = t142 * pkin(5) + t28 * qJ(6) - qJD(6) * t211 + t202;
t10 = [0, 0, 0, 0.2e1 * t194 * t230, t248 * t296, t251, -t252, 0, -pkin(7) * t251 + t194 * t222, pkin(7) * t252 + t196 * t222, t107 * t158 - t108 * t155 + t119 * t213 + t94 * t154 + t93 * t156 - t79 * t165 + t215, -t107 * t93 + t108 * t94 + t78 * t119 + t79 * t120 + (t172 + t216) * t241, t208 * t122 + t80 * t261 -(t195 * t121 - t122 * t193) * t158 + (-t195 * t239 - t271 + (-t193 * t121 - t122 * t195) * qJD(4)) * t166, t122 * t155 + t166 * t125 + t208 * t148 + t80 * t165, t121 * t155 - t148 * t299 - t239 * t165 - t166 * t256, t148 * t155 + t109 (-t120 * t245 + t89) * t148 + t101 * t142 + (-t103 * t245 + t74) * t165 + t54 * t155 - t93 * t121 + t119 * t239 + t102 * t234 + ((-qJD(4) * t106 - t94) * t148 + (-qJD(4) * t84 - t79) * t165 - t102 * t158 + t215) * t193, t208 * t102 + t119 * t80 + t93 * t122 - t250 * t142 - t206 * t148 - t55 * t155 - t207 * t165 + t78 * t261, -t211 * t33 - t28 * t99, -t211 * t34 + t28 * t98 - t99 * t29 + t33 * t63, -t33 * t141 + t99 * t142 + t155 * t211 - t28 * t165, -t34 * t141 - t98 * t142 - t63 * t155 - t29 * t165, t141 * t155 + t109, t141 * t201 + t142 * t227 + t16 * t155 + t165 * t202 + t92 * t29 + t60 * t34 + t52 * t98 + t59 * t63, -t212 * t141 - t276 * t142 - t17 * t155 + t223 * t165 + t211 * t59 - t92 * t28 - t60 * t33 + t52 * t99, -t1 * t99 + t19 * t28 - t2 * t98 - t21 * t29 - t211 * t3 + t5 * t33 - t7 * t34 - t4 * t63, t2 * t21 + t7 * t4 + t1 * t19 + t5 * t3 + t15 * (t98 * pkin(5) + t92) + t37 * (t34 * pkin(5) + t59); 0, 0, 0, -t194 * t253, t248 * t198, 0, 0, 0, t198 * pkin(1) * t194, pkin(1) * t253 (t108 - t110) * t156 + (-t111 + t107) * t154 + (-t190 * t142 - t191 * t213) * pkin(2), t107 * t110 - t108 * t111 + (-t172 * t247 + t190 * t79 - t191 * t78) * pkin(2), t122 * t220 + t271 (t80 + t267) * t195 + (-t239 - t265) * t193, t148 * t220 + t256 - t264, t214 - t266, -t148 * t156, t183 * t239 - t78 * t195 - t54 * t156 + t110 * t121 + (-t182 * t245 - t88) * t148 + (t111 * t148 + t204) * t193, -t110 * t122 + t55 * t156 + t183 * t80 + t78 * t193 + (t182 * t246 + t270) * t148 + t204 * t195, -t28 * t169 - t211 * t269, t219 + t286, -t272 + t287, t217 + t274, -t141 * t156, t293 * t141 + t142 * t221 - t16 * t156 + t173 * t29 - t210 * t52 + t218 * t63 + t268 * t60, t292 * t141 - t249 * t142 + t17 * t156 + t52 * t169 - t173 * t28 + t218 * t211 - t269 * t60, -t1 * t169 + t2 * t210 - t211 * t279 - t268 * t7 + t269 * t5 + t76 * t28 - t280 * t63 - t77 * t29, t2 * t77 + t1 * t76 + t15 * (-pkin(5) * t210 + t173) + t280 * t7 + t279 * t5 + (t268 * pkin(5) + t218) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154 ^ 2 - t156 ^ 2, t107 * t156 - t108 * t154 + t180, 0, 0, 0, 0, 0, t214 + t266, -t148 ^ 2 * t195 - t256 - t264, 0, 0, 0, 0, 0, t217 - t274, -t272 - t287, t219 - t286, t1 * t210 - t37 * t156 + t2 * t169 - t268 * t5 - t269 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t122 * t121, -t121 ^ 2 + t122 ^ 2, t80 - t267, -t239 + t265, t142, -t102 * t122 + t55 * t148 + t203, -t102 * t121 + t54 * t148 - t207, t302, t304, t298, t288, t142, -t228 * t141 + (-t122 * t63 - t141 * t244 + t283 * t142) * pkin(4) + t289, t278 * t141 + (-t122 * t211 - t141 * t232 - t192 * t142) * pkin(4) + t297, t187 * t28 - t5 * t63 + t7 * t211 + t9 * t63 + t8 * t211 + (-t192 * t29 + (t192 * t211 - t283 * t63) * qJD(5)) * pkin(4), -pkin(5) * t295 + t1 * t187 - t5 * t8 - t7 * t9 + (-t37 * t122 + t2 * t192 + (-t192 * t5 + t283 * t7) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, t304, t298, t288, t142, t17 * t141 + t289, t16 * t141 + t297, pkin(5) * t28 - t284 * t63, t284 * t7 + (t1 - t295) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61 - t285, t211 * t5 + t7 * t63 + t15;];
tauc_reg  = t10;
