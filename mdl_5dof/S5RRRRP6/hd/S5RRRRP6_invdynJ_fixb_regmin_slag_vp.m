% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:11
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:10:24
% EndTime: 2021-01-16 00:10:36
% DurationCPUTime: 3.53s
% Computational Cost: add. (5087->406), mult. (11323->512), div. (0->0), fcn. (7962->10), ass. (0->225)
t293 = cos(qJ(3));
t223 = qJD(3) * t293;
t166 = cos(qJ(2));
t168 = pkin(7) + pkin(6);
t130 = t168 * t166;
t119 = qJD(1) * t130;
t162 = sin(qJ(3));
t110 = t162 * t119;
t163 = sin(qJ(2));
t129 = t168 * t163;
t117 = qJD(1) * t129;
t77 = -t293 * t117 - t110;
t308 = -pkin(2) * t223 + t77;
t159 = qJ(2) + qJ(3);
t154 = sin(t159);
t167 = cos(qJ(1));
t255 = t154 * t167;
t164 = sin(qJ(1));
t256 = t154 * t164;
t304 = g(1) * t255 + g(2) * t256;
t224 = qJD(1) * t293;
t254 = t162 * t166;
t109 = -qJD(1) * t254 - t163 * t224;
t161 = sin(qJ(4));
t165 = cos(qJ(4));
t116 = t293 * t163 + t254;
t227 = t293 * t166;
t191 = -t162 * t163 + t227;
t239 = qJD(2) + qJD(3);
t80 = t239 * t191;
t173 = t80 * qJD(1);
t172 = t116 * qJDD(1) + t173;
t212 = t165 * t239;
t238 = qJDD(2) + qJDD(3);
t244 = qJD(4) * t161;
t190 = qJD(4) * t212 + t109 * t244 + t161 * t238 + t165 * t172;
t246 = qJD(1) * t163;
t107 = t162 * t246 - t166 * t224;
t103 = qJD(4) + t107;
t87 = -t161 * t109 - t212;
t267 = t87 * t103;
t307 = t190 - t267;
t206 = g(1) * t167 + g(2) * t164;
t243 = qJD(4) * t165;
t220 = t87 * pkin(4) + qJD(5);
t278 = qJD(2) * pkin(2);
t112 = -t117 + t278;
t73 = t293 * t112 - t110;
t62 = -t239 * pkin(3) - t73;
t42 = t220 + t62;
t245 = qJD(3) * t162;
t242 = qJD(1) * qJD(2);
t221 = t166 * t242;
t241 = t163 * qJDD(1);
t83 = qJDD(2) * pkin(2) + t168 * (-t221 - t241);
t222 = t163 * t242;
t240 = t166 * qJDD(1);
t86 = t168 * (-t222 + t240);
t216 = t112 * t245 + t119 * t223 + t162 * t86 - t293 * t83;
t26 = -t238 * pkin(3) + t216;
t187 = t165 * t109 - t161 * t239;
t31 = -qJD(4) * t187 + t161 * t172 - t165 * t238;
t8 = t31 * pkin(4) + qJDD(5) + t26;
t306 = t8 * t161 + t42 * t243;
t266 = t187 * t103;
t305 = t31 - t266;
t111 = t293 * t119;
t76 = -t162 * t117 + t111;
t210 = pkin(2) * t245 - t76;
t287 = t165 * pkin(4);
t150 = pkin(3) + t287;
t155 = cos(t159);
t160 = -qJ(5) - pkin(8);
t215 = t155 * t150 - t154 * t160;
t286 = t166 * pkin(2);
t152 = pkin(1) + t286;
t303 = -t116 * pkin(8) - t152;
t261 = t107 * t161;
t302 = -qJ(5) * t261 + t165 * qJD(5);
t70 = -t109 * pkin(3) + t107 * pkin(8);
t61 = pkin(2) * t246 + t70;
t301 = t161 * t61 + t308 * t165;
t128 = t152 * qJD(1);
t59 = t107 * pkin(3) + t109 * pkin(8) - t128;
t74 = t162 * t112 + t111;
t63 = t239 * pkin(8) + t74;
t33 = t161 * t59 + t165 * t63;
t300 = -t33 * t109 + t26 * t161 + t62 * t243;
t289 = g(3) * t161;
t250 = t167 * t165;
t253 = t164 * t161;
t96 = t155 * t253 + t250;
t251 = t167 * t161;
t252 = t164 * t165;
t98 = -t155 * t251 + t252;
t299 = -g(1) * t98 + g(2) * t96 + t154 * t289;
t19 = -t87 * qJ(5) + t33;
t260 = t107 * t165;
t298 = -t19 * t109 + t42 * t260 + t306;
t297 = t187 ^ 2;
t203 = -qJDD(1) * t227 + t162 * t241;
t81 = t239 * t116;
t49 = qJD(1) * t81 + t203;
t48 = qJDD(4) + t49;
t294 = t48 * pkin(4);
t147 = g(3) * t154;
t290 = g(3) * t155;
t32 = -t161 * t63 + t165 * t59;
t18 = qJ(5) * t187 + t32;
t14 = t103 * pkin(4) + t18;
t285 = -t18 + t14;
t284 = t161 * t70 + t165 * t73;
t149 = t162 * pkin(2) + pkin(8);
t249 = -qJ(5) - t149;
t214 = qJD(4) * t249;
t282 = t161 * t214 - t301 + t302;
t156 = t165 * qJ(5);
t204 = -t109 * pkin(4) + t107 * t156;
t58 = t165 * t61;
t281 = t165 * t214 - t204 - t58 + (-qJD(5) + t308) * t161;
t72 = -pkin(3) * t191 + t303;
t92 = -t162 * t129 + t293 * t130;
t84 = t165 * t92;
t280 = t161 * t72 + t84;
t279 = pkin(4) * qJD(4);
t277 = t161 * t80;
t276 = t165 * t48;
t275 = t165 * t80;
t274 = t165 * t187;
t272 = t190 * qJ(5);
t271 = t190 * t161;
t270 = t31 * qJ(5);
t268 = t62 * t107;
t218 = qJD(4) * t160;
t265 = t161 * t218 - t284 + t302;
t66 = t165 * t70;
t264 = t165 * t218 - t204 - t66 + (-qJD(5) + t73) * t161;
t233 = pkin(4) * t244;
t93 = pkin(4) * t261;
t263 = t233 + t93 + t210;
t262 = t103 * t109;
t259 = t109 * t107;
t258 = t116 * t161;
t248 = t304 * t165;
t157 = t163 ^ 2;
t247 = -t166 ^ 2 + t157;
t236 = t163 * t278;
t41 = t81 * pkin(3) - t80 * pkin(8) + t236;
t228 = qJD(2) * t168;
t118 = t163 * t228;
t120 = t166 * t228;
t192 = -t293 * t129 - t162 * t130;
t44 = t192 * qJD(3) - t293 * t118 - t162 * t120;
t237 = t161 * t41 + t165 * t44 + t72 * t243;
t234 = qJD(4) * pkin(8) * t103;
t35 = t42 * t244;
t54 = t62 * t244;
t230 = t206 * t155 + t147;
t229 = -t8 - t290;
t226 = t116 * t243;
t225 = -t26 - t290;
t219 = pkin(4) * t161 + t168;
t145 = pkin(2) * t222;
t17 = t49 * pkin(3) - pkin(8) * t173 + t303 * qJDD(1) + t145;
t177 = t112 * t223 - t119 * t245 + t162 * t83 + t293 * t86;
t25 = t238 * pkin(8) + t177;
t217 = t161 * t17 + t165 * t25 + t59 * t243 - t63 * t244;
t213 = t103 * t165;
t151 = -t293 * pkin(2) - pkin(3);
t209 = -g(1) * t96 - g(2) * t98;
t97 = -t155 * t252 + t251;
t99 = t155 * t250 + t253;
t208 = -g(1) * t97 - g(2) * t99;
t207 = -pkin(8) * t48 + t268;
t205 = g(1) * t164 - g(2) * t167;
t202 = -t149 * t48 + t268;
t201 = -t14 * t165 - t161 * t19;
t199 = -qJ(5) * t80 - qJD(5) * t116;
t198 = t150 * t154 + t155 * t160;
t197 = t32 * t109 + t248 + t54;
t196 = t206 * t154;
t195 = -0.2e1 * pkin(1) * t242 - pkin(6) * qJDD(2);
t194 = t226 + t277;
t193 = t152 + t215;
t133 = t155 * t289;
t184 = -t161 * t196 + t133;
t183 = g(1) * t99 - g(2) * t97 + t165 * t147 - t217;
t169 = qJD(2) ^ 2;
t182 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t169 + t205;
t170 = qJD(1) ^ 2;
t181 = pkin(1) * t170 - pkin(6) * qJDD(1) + t206;
t16 = t165 * t17;
t180 = -qJD(4) * t33 - t161 * t25 + t16;
t179 = t14 * t109 + t229 * t165 + t248 + t35;
t178 = -t128 * t109 - t216 - t290 + t304;
t45 = t92 * qJD(3) - t162 * t118 + t293 * t120;
t1 = qJD(5) * t187 + t180 - t272 + t294;
t3 = -t87 * qJD(5) + t217 - t270;
t176 = qJD(4) * t201 - t1 * t161 - t14 * t260 + t3 * t165 - t19 * t261 - t230;
t175 = t180 + t299;
t174 = -t128 * t107 - t177 + t230;
t127 = t165 * pkin(8) + t156;
t126 = t160 * t161;
t125 = t151 - t287;
t114 = t165 * t149 + t156;
t113 = t249 * t161;
t104 = -t152 * qJDD(1) + t145;
t85 = t87 ^ 2;
t68 = t165 * t72;
t60 = pkin(4) * t258 - t192;
t51 = -t107 ^ 2 + t109 ^ 2;
t47 = -t93 + t74;
t40 = -t203 + (-qJD(1) * t116 - t109) * t239;
t39 = t107 * t239 + t172;
t38 = t165 * t41;
t34 = -qJ(5) * t258 + t280;
t28 = -pkin(4) * t191 - t116 * t156 - t161 * t92 + t68;
t22 = pkin(4) * t194 + t45;
t11 = t103 * t213 - t109 * t187 + t161 * t48;
t10 = -t103 ^ 2 * t161 - t87 * t109 + t276;
t9 = -t187 * t213 + t271;
t6 = -qJ(5) * t226 + (-qJD(4) * t92 + t199) * t161 + t237;
t5 = t81 * pkin(4) - t161 * t44 + t38 + t199 * t165 + (-t84 + (qJ(5) * t116 - t72) * t161) * qJD(4);
t4 = -t305 * t161 + t307 * t165;
t2 = [qJDD(1), t205, t206, t157 * qJDD(1) + 0.2e1 * t163 * t221, 0.2e1 * t163 * t240 - 0.2e1 * t247 * t242, qJDD(2) * t163 + t169 * t166, qJDD(2) * t166 - t169 * t163, 0, t163 * t195 + t166 * t182, -t163 * t182 + t166 * t195, -t109 * t80 + t116 * t172, -t80 * t107 + t109 * t81 - t116 * t49 + t172 * t191, t116 * t238 + t239 * t80, t191 * t238 - t239 * t81, 0, -t104 * t191 + t107 * t236 - t128 * t81 - t152 * t49 + t155 * t205 + t192 * t238 - t239 * t45, -g(1) * t256 + g(2) * t255 + t104 * t116 - t109 * t236 - t128 * t80 - t152 * t172 - t238 * t92 - t239 * t44, -t80 * t274 + (t165 * t190 + t187 * t244) * t116, (t161 * t187 - t165 * t87) * t80 + (-t271 - t165 * t31 + (t161 * t87 + t274) * qJD(4)) * t116, t116 * t276 - t190 * t191 - t187 * t81 + (-t116 * t244 + t275) * t103, -t103 * t194 + t191 * t31 - t258 * t48 - t87 * t81, t103 * t81 - t191 * t48, (-t243 * t92 + t38) * t103 + t68 * t48 - (-t243 * t63 + t16) * t191 + t32 * t81 + t45 * t87 - t192 * t31 + t62 * t226 + ((-qJD(4) * t72 - t44) * t103 - t92 * t48 - (-qJD(4) * t59 - t25) * t191 + t26 * t116 + t62 * t80) * t161 + t208, -(-t244 * t92 + t237) * t103 - t280 * t48 + t217 * t191 - t33 * t81 - t45 * t187 - t192 * t190 + t62 * t275 + (t26 * t165 - t54) * t116 + t209, -t1 * t191 + t5 * t103 + t306 * t116 + t14 * t81 + t22 * t87 + t42 * t277 + t28 * t48 + t60 * t31 + t208, t42 * t275 - t6 * t103 + t3 * t191 - t19 * t81 - t22 * t187 + t60 * t190 - t34 * t48 + (t8 * t165 - t35) * t116 + t209, -t28 * t190 - t34 * t31 + t5 * t187 - t6 * t87 + t201 * t80 + t205 * t154 + (-t1 * t165 - t161 * t3 + (t14 * t161 - t165 * t19) * qJD(4)) * t116, t1 * t28 + t14 * t5 + t19 * t6 + t42 * t22 + t3 * t34 + t8 * t60 + (-g(1) * t219 - g(2) * t193) * t167 + (g(1) * t193 - g(2) * t219) * t164; 0, 0, 0, -t163 * t170 * t166, t247 * t170, t241, t240, qJDD(2), -g(3) * t166 + t163 * t181, g(3) * t163 + t166 * t181, -t259, t51, t39, t40, t238, t76 * t239 + (-t107 * t246 + t293 * t238 - t239 * t245) * pkin(2) + t178, t77 * t239 + (t109 * t246 - t162 * t238 - t223 * t239) * pkin(2) + t174, t9, t4, t11, t10, t262, t151 * t31 + t210 * t87 + t225 * t165 + t202 * t161 + (-t149 * t243 + t308 * t161 - t58) * t103 + t197, t151 * t190 - t210 * t187 + t202 * t165 + (t149 * t244 + t301) * t103 + t184 + t300, t103 * t281 + t113 * t48 + t125 * t31 + t261 * t42 + t263 * t87 + t179, -t103 * t282 - t114 * t48 + t125 * t190 - t187 * t263 + t184 + t298, -t113 * t190 - t114 * t31 + t187 * t281 - t282 * t87 + t176, t3 * t114 + t1 * t113 + t8 * t125 - g(3) * (t215 + t286) + t263 * t42 + t282 * t19 + t281 * t14 + t206 * (pkin(2) * t163 + t198); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t259, t51, t39, t40, t238, t239 * t74 + t178, t239 * t73 + t174, t9, t4, t11, t10, t262, -pkin(3) * t31 - t66 * t103 - t74 * t87 + (t73 * t103 + t207) * t161 + (t225 - t234) * t165 + t197, -pkin(3) * t190 + t284 * t103 + t74 * t187 + t133 + t207 * t165 + (-t196 + t234) * t161 + t300, t126 * t48 - t150 * t31 - t47 * t87 + (t107 * t42 + t279 * t87) * t161 + t264 * t103 + t179, -t127 * t48 - t150 * t190 + t47 * t187 + t133 - t265 * t103 + (-t187 * t279 - t196) * t161 + t298, -t126 * t190 - t127 * t31 + t187 * t264 - t265 * t87 + t176, t3 * t127 + t1 * t126 - t8 * t150 - g(3) * t215 + (-t47 + t233) * t42 + t265 * t19 + t264 * t14 + t206 * t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t187 * t87, -t85 + t297, t190 + t267, -t31 - t266, t48, t33 * t103 + t187 * t62 + t175, t32 * t103 + t62 * t87 + t183, 0.2e1 * t294 - t272 + t19 * t103 - (-t220 - t42) * t187 + t175, -t297 * pkin(4) + t270 + t18 * t103 + (qJD(5) + t42) * t87 + t183, -pkin(4) * t190 - t285 * t87, t285 * t19 + (t187 * t42 + t1 + t299) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t305, t307, -t85 - t297, -t14 * t187 + t19 * t87 - t229 - t304;];
tau_reg = t2;
