% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% tau_reg [6x35]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRRR14V3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR14V3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14V3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_invdynJ_fixb_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:10:12
% EndTime: 2019-04-12 15:10:38
% DurationCPUTime: 8.60s
% Computational Cost: add. (3544->486), mult. (8525->701), div. (0->0), fcn. (6853->10), ass. (0->222)
t113 = sin(qJ(2));
t216 = qJD(4) * t113;
t305 = -qJD(1) * t216 + qJDD(2);
t112 = sin(qJ(4));
t304 = t305 * t112;
t117 = cos(qJ(4));
t101 = t113 * qJDD(1);
t118 = cos(qJ(2));
t207 = t118 * qJD(1);
t129 = (qJD(4) + t207) * qJD(2) + t101;
t123 = t129 * t112 - t117 * t305;
t222 = qJD(2) * t112;
t223 = qJD(1) * t113;
t76 = t117 * t223 + t222;
t93 = -qJD(4) + t207;
t262 = t76 * t93;
t303 = t123 - t262;
t302 = qJ(3) * t123;
t301 = (t117 * t129 + t304) * qJ(3);
t114 = sin(qJ(1));
t119 = cos(qJ(1));
t169 = g(1) * t119 + g(2) * t114;
t276 = g(3) * t118;
t131 = -t169 * t113 + t276;
t297 = qJ(3) * t93 ^ 2 + qJDD(3) + t131;
t111 = sin(qJ(5));
t116 = cos(qJ(5));
t106 = t117 * qJD(2);
t284 = -t112 * t223 + t106;
t69 = t284 * qJ(3);
t159 = t116 * qJD(3) - t111 * t69;
t296 = qJD(5) * t159;
t295 = t284 * qJD(3);
t110 = sin(qJ(6));
t115 = cos(qJ(6));
t238 = t113 * t114;
t226 = t119 * t112;
t228 = t117 * t118;
t67 = t114 * t228 - t226;
t44 = t111 * t238 + t116 * t67;
t233 = t114 * t118;
t66 = t112 * t233 + t117 * t119;
t293 = t110 * t44 - t115 * t66;
t292 = t110 * t66 + t115 * t44;
t212 = qJD(5) * t111;
t213 = qJD(4) * t117;
t70 = t111 * t113 + t116 * t228;
t57 = t70 * qJD(1);
t291 = -t112 * t212 + t116 * t213 - t57;
t153 = t111 * t93;
t162 = -t116 * t76 + t153;
t73 = qJD(5) - t284;
t19 = -t110 * t162 - t115 * t73;
t290 = t19 * t73;
t78 = t116 * t93;
t48 = t111 * t76 + t78;
t289 = t48 * t93;
t32 = qJDD(5) + t123;
t144 = -t116 * t32 + t73 * t212;
t263 = t284 * t93;
t206 = qJD(1) * qJD(2);
t188 = t118 * t206;
t33 = qJD(4) * t106 + (t101 + t188) * t117 + t304;
t288 = t33 - t263;
t287 = qJD(3) * t76;
t278 = g(2) * t119;
t168 = g(1) * t114 - t278;
t285 = t168 * t113;
t191 = t113 * t213;
t220 = qJD(2) * t118;
t283 = -t112 * t220 - t191;
t210 = qJD(5) * t116;
t104 = t118 * qJDD(1);
t74 = t113 * t206 + qJDD(4) - t104;
t11 = t111 * t74 + t116 * t33 - t93 * t210 - t76 * t212;
t21 = t110 * t73 - t115 * t162;
t4 = qJD(6) * t21 + t11 * t110 - t115 * t32;
t277 = g(3) * t113;
t282 = t169 * t118 + t277;
t259 = qJ(3) * t74;
t280 = qJD(3) * (qJD(4) + t93) - t259;
t14 = t295 - t302;
t7 = qJDD(3) * t111 + t116 * t14 + t296;
t12 = -qJD(5) * t162 + t111 * t33 - t116 * t74;
t15 = t287 + t301;
t52 = qJD(3) * t111 + t116 * t69;
t63 = t76 * qJ(3);
t163 = t110 * t52 - t115 * t63;
t1 = -t163 * qJD(6) + t110 * t15 + t115 * t7;
t108 = t113 ^ 2;
t224 = -t118 ^ 2 + t108;
t279 = 0.2e1 * t113 * t104 - 0.2e1 * t224 * t206;
t208 = qJD(6) * t115;
t209 = qJD(6) * t110;
t3 = t115 * t11 + t110 * t32 + t162 * t209 + t73 * t208;
t275 = t110 * t3;
t40 = qJD(6) + t48;
t274 = t19 * t40;
t273 = t21 * t40;
t272 = t21 * t284;
t271 = t3 * t111;
t246 = t110 * t116;
t35 = -t115 * t76 + t246 * t284;
t270 = t35 * t40;
t232 = t115 * t116;
t36 = t110 * t76 + t232 * t284;
t269 = t36 * t40;
t268 = t4 * t111;
t267 = t48 * t73;
t266 = t48 * t76;
t265 = t162 * t73;
t264 = t162 * t76;
t214 = qJD(4) * t116;
t176 = -qJD(6) + t214;
t177 = -qJD(6) * t116 + qJD(4);
t211 = qJD(5) * t115;
t189 = t111 * t211;
t193 = t112 * t207;
t231 = t115 * t117;
t261 = t176 * t231 + (t177 * t110 - t189) * t112 - t110 * t193 - t115 * t57;
t215 = qJD(4) * t115;
t260 = t115 * t193 - t117 * t209 + (t208 * t116 - t215) * t112 + t291 * t110;
t258 = qJ(3) * t93;
t257 = t11 * t111;
t10 = qJDD(6) + t12;
t256 = t110 * t10;
t255 = t110 * t40;
t253 = t111 * t32;
t252 = t111 * t67;
t251 = t111 * t73;
t250 = t112 * t33;
t249 = t115 * t10;
t182 = t115 * t40;
t247 = t117 * t93;
t245 = t111 * t112;
t244 = t111 * t115;
t243 = t111 * t117;
t242 = t112 * t113;
t241 = t112 * t115;
t240 = t112 * t116;
t239 = t112 * t118;
t237 = t113 * t117;
t236 = t113 * t119;
t121 = qJD(1) ^ 2;
t235 = t113 * t121;
t234 = t114 * t116;
t230 = t116 * t117;
t229 = t116 * t118;
t227 = t118 * t119;
t221 = qJD(2) * t113;
t219 = qJD(3) * t113;
t218 = qJD(4) * t111;
t217 = qJD(4) * t112;
t205 = qJD(1) * qJD(3);
t204 = 0.2e1 * t205;
t203 = 0.2e1 * qJD(3) * qJD(2);
t202 = t73 * t218;
t201 = t73 * t214;
t200 = qJD(5) * t255;
t198 = t40 * t211;
t197 = t110 * t242;
t196 = t113 * t241;
t195 = t118 * t235;
t183 = -t116 * qJDD(3) + t111 * t14;
t181 = t116 * t73;
t180 = t224 * t121;
t120 = qJD(2) ^ 2;
t179 = -t108 * t121 - t120;
t175 = -g(1) * t236 - g(2) * t238 + t276;
t71 = -t114 * t117 + t118 * t226;
t171 = g(1) * t71 + g(2) * t66;
t72 = t112 * t114 + t117 * t227;
t170 = -g(1) * t72 - g(2) * t67;
t165 = -qJD(4) * t258 + qJDD(3);
t164 = -t40 * t63 - t296;
t23 = t110 * t63 + t115 * t52;
t160 = -qJD(3) * t73 + qJD(4) * t63;
t158 = t93 * t162;
t157 = t93 * t73;
t155 = t251 * t284 - t144;
t154 = -qJDD(3) - t175;
t151 = qJDD(2) * t118 - t113 * t120;
t149 = t117 * t48 + t73 * t245;
t148 = -t117 * t162 + t73 * t240;
t147 = t112 * t74 - t93 * t213;
t146 = t117 * t74 + t93 * t217;
t145 = -t73 * t210 - t253;
t143 = -t40 * t208 - t256;
t142 = -t40 * t209 + t249;
t65 = -t111 * t118 + t113 * t230;
t42 = t115 * t65 + t197;
t68 = -t110 * t117 + t112 * t232;
t62 = t110 * t240 + t231;
t46 = -t111 * t72 + t116 * t236;
t64 = t111 * t237 + t229;
t140 = -g(1) * t46 - g(2) * (t113 * t234 - t252) + g(3) * t64;
t139 = -t177 + t207;
t8 = qJD(5) * t52 + t183;
t138 = -t140 + t8;
t137 = g(3) * t242 + t171;
t136 = qJDD(1) * t108 + 0.2e1 * t113 * t188;
t56 = -t116 * t223 + t207 * t243;
t135 = t111 * t213 + t112 * t210 - t56;
t134 = t137 - t15;
t133 = t40 * (-qJD(6) - t78);
t2 = -qJD(6) * t23 - t110 * t7 + t115 * t15;
t125 = -t21 * t245 + t40 * t68;
t124 = -t19 * t245 + t40 * t62;
t122 = qJ(3) ^ 2;
t95 = g(1) * t233;
t54 = t65 * t119;
t53 = t65 * t114;
t47 = t111 * t236 + t116 * t72;
t41 = t110 * t65 - t196;
t27 = (-qJD(5) + t106) * t229 + (-t112 * t214 + (-qJD(5) * t117 + qJD(2)) * t111) * t113;
t26 = -t216 * t245 - t118 * t212 - t116 * t221 + (t111 * t220 + t113 * t210) * t117;
t18 = t110 * t71 + t115 * t47;
t17 = -t110 * t47 + t115 * t71;
t16 = t21 * t212;
t6 = (qJD(6) * t242 + t27) * t115 + (-qJD(6) * t65 - t283) * t110;
t5 = qJD(6) * t42 + t110 * t27 + t115 * t283;
t9 = [qJDD(1), t168, t169, t136, t279, qJDD(2) * t113 + t118 * t120, t151, 0, -g(2) * t227 + t95, -t285, t95 + (t113 * t204 - t278) * t118 + qJ(3) * t279, t151 * qJ(3) + qJDD(3) * t113 + t118 * t203 - t169, 0.2e1 * (qJ(3) * qJDD(1) + t205) * t108 + (0.4e1 * qJ(3) * t188 + t168) * t113, t136 * t122 + (t108 * t204 + t285) * qJ(3), t118 * t76 * t106 + (t117 * t33 - t217 * t76) * t113 (-t112 * t76 + t117 * t284) * t220 + (-t250 - t117 * t123 + (-t112 * t284 - t117 * t76) * qJD(4)) * t113 (-t106 * t93 - t33) * t118 + (qJD(2) * t76 + t146) * t113 (t222 * t93 + t123) * t118 + (qJD(2) * t284 - t147) * t113, -t118 * t74 - t221 * t93, g(1) * t67 - g(2) * t72 + (t15 + (qJ(3) * t247 + qJD(3) * t112) * qJD(2)) * t118 + (-qJD(2) * t63 + t165 * t112 + t117 * t280) * t113, -g(1) * t66 + g(2) * t71 + (t14 + (qJD(3) * t117 - t112 * t258) * qJD(2)) * t118 + (-qJD(2) * t69 - t112 * t280 + t165 * t117) * t113, t11 * t65 - t162 * t27, -t11 * t64 - t12 * t65 + t162 * t26 - t27 * t48, -t162 * t191 + t27 * t73 + t32 * t65 + (t11 * t113 - t162 * t220) * t112, -t48 * t191 - t26 * t73 - t32 * t64 + (-t113 * t12 - t220 * t48) * t112, t73 * t191 + (t113 * t32 + t220 * t73) * t112, g(1) * t44 - g(2) * t47 + t15 * t64 + t26 * t63 + (qJ(3) * t149 + t112 * t159) * t220 + (t159 * t213 - t112 * t8 + t149 * qJD(3) + ((t12 + t202) * t117 + (-qJD(4) * t48 - t145) * t112) * qJ(3)) * t113, -g(1) * t252 - g(2) * t46 + t15 * t65 + t63 * t27 + (qJ(3) * t148 - t112 * t52) * t220 + (g(1) * t234 - t52 * t213 - t7 * t112 + t148 * qJD(3) + ((t11 + t201) * t117 + (qJD(4) * t162 - t144) * t112) * qJ(3)) * t113, t21 * t6 + t3 * t42, -t19 * t6 - t21 * t5 - t3 * t41 - t4 * t42, t10 * t42 + t21 * t26 + t3 * t64 + t40 * t6, -t10 * t41 - t19 * t26 - t4 * t64 - t40 * t5, t10 * t64 + t26 * t40, t2 * t64 - t163 * t26 + t8 * t41 - t159 * t5 + g(1) * t292 - g(2) * t18 + t124 * t219 + (t124 * t220 + ((t176 * t255 - t19 * t218 + t249) * t117 + ((-t110 * t212 - t215) * t40 - t268 + (-qJD(5) * t19 - t143) * t116) * t112) * t113) * qJ(3), -t1 * t64 - t23 * t26 + t8 * t42 - t159 * t6 - g(1) * t293 - g(2) * t17 + t125 * t219 + (t125 * t220 + ((t176 * t182 - t21 * t218 - t256) * t117 + (-(-qJD(4) * t110 + t189) * t40 - t271 + (-qJD(5) * t21 + t142) * t116) * t112) * t113) * qJ(3); 0, 0, 0, -t195, t180, t101, t104, qJDD(2), -t175, t282, qJ(3) * t180 + t154, qJ(3) * t104, -t277 + 0.2e1 * qJ(3) * qJDD(2) + t203 + (-0.2e1 * qJ(3) * t235 - t169) * t118 (qJDD(2) - t195) * t122 + (-t282 + t203) * qJ(3), -t247 * t76 + t250, -t303 * t112 + t288 * t117 (-t113 * t76 + t228 * t93) * qJD(1) + t147 (-t113 * t284 - t239 * t93) * qJD(1) + t146, t93 * t223, -t259 * t112 - t297 * t117 + t63 * t223, t297 * t112 - t259 * t117 + t69 * t223, t11 * t240 - t291 * t162, t48 * t57 - t162 * t56 + (t111 * t162 - t116 * t48) * t213 + (-t257 - t116 * t12 + (t111 * t48 + t116 * t162) * qJD(5)) * t112, -t57 * t73 + (-t11 + t201) * t117 + (t158 - t144) * t112, t56 * t73 + (t12 - t202) * t117 + (t145 + t289) * t112, -t112 * t157 - t117 * t32, g(1) * t54 + g(2) * t53 - g(3) * t70 - t56 * t63 + (t8 + t160 * t111 + (t145 - t289) * qJ(3)) * t117 + (t63 * t210 + qJD(3) * t48 + t111 * t15 - t93 * t159 + (-t111 * t157 + t12) * qJ(3)) * t112, -t63 * t57 - t282 * t116 + (t7 + t160 * t116 + t131 * t111 + (t158 + t144) * qJ(3)) * t117 + (-t63 * t212 - qJD(3) * t162 + t15 * t116 + t93 * t52 + (-t73 * t78 + t11) * qJ(3)) * t112, t21 * t261 + t3 * t68, -t19 * t261 - t21 * t260 - t3 * t62 - t4 * t68, t10 * t68 + t135 * t21 + t245 * t3 + t261 * t40, -t10 * t62 - t135 * t19 - t245 * t4 - t260 * t40, t10 * t245 + t135 * t40, t2 * t245 + t8 * t62 - g(1) * (-t115 * t54 - t119 * t197) - g(2) * (-t114 * t197 - t115 * t53) - g(3) * (t110 * t239 + t115 * t70) - t260 * t159 - t135 * t163 + ((-t110 * t230 + t241) * t40 + t19 * t243) * qJD(3) + ((t110 * t133 + t153 * t19 + t249) * t112 + (-t10 * t246 + t268 + (t111 * t255 + t116 * t19) * qJD(5) - t139 * t182) * t117) * qJ(3), -t1 * t245 + t8 * t68 - g(1) * (t110 * t54 - t119 * t196) - g(2) * (t110 * t53 - t114 * t196) - g(3) * (-t110 * t70 + t115 * t239) - t261 * t159 - t135 * t23 + (-(t110 * t112 + t115 * t230) * t40 + t21 * t243) * qJD(3) + ((t115 * t133 + t153 * t21 - t256) * t112 + (-t10 * t232 + t271 + (t116 * t21 + t244 * t40) * qJD(5) + t139 * t255) * t117) * qJ(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) - t195, t101, t179, qJ(3) * t179 - t154, 0, 0, 0, 0, 0, t303, t288, 0, 0, 0, 0, 0, t155 - t266, -t116 * t73 ^ 2 - t253 + t264, 0, 0, 0, 0, 0, t270 + (-t4 - t200) * t116 + (t143 + t290) * t111, t269 + t16 + (-t3 - t198) * t116 + (-t142 - t272) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76 * t284, -t284 ^ 2 + t76 ^ 2, t33 + t263, -t123 - t262, t74, -t69 * t93 + t137 - 0.2e1 * t287 - t301, g(3) * t237 + t63 * t93 - t170 - 0.2e1 * t295 + t302, -t162 * t181 + t257 (t11 - t267) * t116 + (-t12 + t265) * t111, t181 * t73 + t253 + t264, t155 + t266, -t73 * t76, t116 * t134 - t159 * t76 - t48 * t69, -t111 * t134 + t162 * t69 + t52 * t76, t3 * t244 + (-t111 * t209 + t115 * t210 - t36) * t21, t19 * t36 + t21 * t35 + (-t110 * t21 - t115 * t19) * t210 + (-t275 - t115 * t4 + (t110 * t19 - t115 * t21) * qJD(6)) * t111, -t269 + t16 + (-t3 + t198) * t116 + (t142 - t272) * t111, t270 + (t4 - t200) * t116 + (t143 - t290) * t111, -t116 * t10 + t251 * t40, -t69 * t182 + t159 * t35 + t170 * t110 + (t110 * t164 + t115 * t171 - t2) * t116 + t68 * t277 + (t8 * t110 - t159 * t208 - t163 * t73 + t63 * t19) * t111, t69 * t255 + t159 * t36 + t170 * t115 + (-t110 * t171 + t115 * t164 + t1) * t116 - t62 * t277 + (t8 * t115 + t159 * t209 + t63 * t21 - t23 * t73) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162 * t48, t162 ^ 2 - t48 ^ 2, t11 + t267, -t12 - t265, t32, t162 * t63 + t140 - t183 + (-qJD(5) + t73) * t52, g(1) * t47 + g(2) * t44 + g(3) * t65 + t159 * t73 + t48 * t63 - t7, t182 * t21 + t275 (t3 - t274) * t115 + (-t4 - t273) * t110, t162 * t21 + t182 * t40 + t256, -t110 * t40 ^ 2 - t162 * t19 + t249, t40 * t162, -t115 * t138 - t162 * t163 - t19 * t52, t110 * t138 - t162 * t23 - t21 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t19, -t19 ^ 2 + t21 ^ 2, t3 + t274, t273 - t4, t10, -g(1) * t17 + g(2) * t293 + g(3) * t41 + t159 * t21 + t23 * t40 + t2, g(1) * t18 + g(2) * t292 + g(3) * t42 - t159 * t19 - t163 * t40 - t1;];
tau_reg  = t9;