% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:57:09
% EndTime: 2019-03-09 09:57:19
% DurationCPUTime: 3.19s
% Computational Cost: add. (5351->415), mult. (13261->523), div. (0->0), fcn. (9229->6), ass. (0->208)
t185 = cos(qJ(2));
t252 = qJD(1) * t185;
t167 = -qJD(4) + t252;
t179 = sin(pkin(9));
t183 = sin(qJ(2));
t253 = qJD(1) * t183;
t234 = t179 * t253;
t180 = cos(pkin(9));
t248 = qJD(2) * t180;
t133 = -t234 + t248;
t238 = t180 * t253;
t249 = qJD(2) * t179;
t134 = t238 + t249;
t182 = sin(qJ(4));
t184 = cos(qJ(4));
t85 = -t184 * t133 + t134 * t182;
t270 = t167 * t85;
t243 = qJD(1) * qJD(2);
t232 = t185 * t243;
t220 = t179 * t232;
t265 = t184 * t180;
t240 = t185 * t265;
t221 = qJD(2) * t240;
t244 = qJD(4) * t184;
t43 = t182 * (qJD(4) * t134 + t220) - qJD(1) * t221 - t133 * t244;
t191 = t43 - t270;
t239 = t179 * t252;
t245 = qJD(4) * t182;
t291 = -t179 * t245 + t180 * t244;
t257 = -qJD(1) * t240 + t182 * t239 + t291;
t289 = t85 ^ 2;
t171 = t183 * t243;
t298 = 0.2e1 * t171;
t283 = pkin(8) + qJ(3);
t148 = t283 * t179;
t149 = t283 * t180;
t213 = pkin(2) * t183 - qJ(3) * t185;
t139 = t213 * qJD(1);
t103 = pkin(7) * t234 + t180 * t139;
t266 = t180 * t185;
t207 = pkin(3) * t183 - pkin(8) * t266;
t74 = qJD(1) * t207 + t103;
t124 = t179 * t139;
t267 = t180 * t183;
t268 = t179 * t185;
t201 = -pkin(7) * t267 - pkin(8) * t268;
t90 = qJD(1) * t201 + t124;
t297 = qJD(3) * t265 - t148 * t244 - t184 * t90 + (-qJD(3) * t179 - qJD(4) * t149 - t74) * t182;
t137 = t179 * t184 + t180 * t182;
t123 = t137 * qJD(4);
t200 = t137 * t185;
t256 = -qJD(1) * t200 + t123;
t211 = t133 * t182 + t134 * t184;
t144 = -pkin(2) * t185 - qJ(3) * t183 - pkin(1);
t126 = t144 * qJD(1);
t173 = pkin(7) * t252;
t151 = qJD(2) * qJ(3) + t173;
t93 = t180 * t126 - t151 * t179;
t52 = -pkin(3) * t252 - pkin(8) * t134 + t93;
t94 = t179 * t126 + t180 * t151;
t60 = pkin(8) * t133 + t94;
t23 = t182 * t60 - t184 * t52;
t208 = pkin(5) * t211 + t23;
t260 = qJD(5) + t208;
t83 = t211 ^ 2;
t296 = -0.2e1 * t243;
t280 = qJ(5) * t253 - t297;
t101 = -t148 * t182 + t149 * t184;
t295 = qJD(3) * t137 + qJD(4) * t101 - t182 * t90 + t184 * t74;
t162 = t167 ^ 2;
t294 = -t83 - t162;
t272 = t211 * t167;
t153 = qJD(5) * t167;
t163 = qJ(5) * t171;
t293 = t163 - t153;
t127 = pkin(3) * t239 + t173;
t292 = t257 * qJ(5) + t137 * qJD(5) + t127;
t290 = qJD(4) * t211;
t288 = pkin(5) * t85;
t172 = pkin(7) * t253;
t273 = qJD(2) * pkin(2);
t225 = qJD(3) - t273;
t143 = t172 + t225;
t102 = -t133 * pkin(3) + t143;
t194 = -qJ(5) * t211 + t102;
t284 = pkin(4) + qJ(6);
t16 = t284 * t85 + t194;
t287 = t16 * t211;
t27 = pkin(4) * t85 + t194;
t286 = t27 * t211;
t285 = t211 * t85;
t233 = t284 * t183;
t282 = t257 * pkin(5) + qJD(1) * t233 + t295;
t281 = t256 * pkin(5) + t280;
t279 = -pkin(4) * t253 - t295;
t136 = t179 * t182 - t265;
t278 = t136 * qJD(6) + t256 * t284 - t292;
t277 = -t256 * pkin(4) + t292;
t24 = t182 * t52 + t184 * t60;
t132 = t180 * t144;
t92 = -pkin(8) * t267 + t132 + (-pkin(7) * t179 - pkin(3)) * t185;
t165 = pkin(7) * t266;
t108 = t179 * t144 + t165;
t269 = t179 * t183;
t99 = -pkin(8) * t269 + t108;
t275 = t182 * t92 + t184 * t99;
t196 = qJD(2) * t200;
t44 = qJD(1) * t196 + t290;
t274 = qJ(5) * t44;
t271 = t85 * qJ(5);
t187 = qJD(1) ^ 2;
t264 = t185 * t187;
t186 = qJD(2) ^ 2;
t263 = t186 * t183;
t262 = t186 * t185;
t261 = -qJD(5) - t23;
t15 = t24 - t288;
t259 = -qJD(6) - t15;
t120 = qJD(2) * t213 - t183 * qJD(3);
t111 = t120 * qJD(1);
t142 = (qJD(3) - t172) * qJD(2);
t73 = t179 * t111 + t180 * t142;
t247 = qJD(2) * t183;
t241 = pkin(7) * t247;
t97 = t180 * t120 + t179 * t241;
t166 = pkin(7) * t232;
t119 = pkin(3) * t220 + t166;
t246 = qJD(2) * t185;
t174 = pkin(7) * t246;
t237 = t179 * t246;
t128 = pkin(3) * t237 + t174;
t140 = pkin(3) * t269 + t183 * pkin(7);
t254 = t183 ^ 2 - t185 ^ 2;
t100 = t148 * t184 + t149 * t182;
t251 = qJD(2) * t100;
t250 = qJD(2) * t101;
t242 = pkin(7) * t268;
t170 = -t180 * pkin(3) - pkin(2);
t231 = -t171 + t285;
t229 = -t182 * t99 + t184 * t92;
t198 = t207 * qJD(2);
t72 = t180 * t111 - t179 * t142;
t47 = qJD(1) * t198 + t72;
t59 = -pkin(8) * t220 + t73;
t228 = -t182 * t47 - t184 * t59 - t52 * t244 + t60 * t245;
t227 = t182 * t59 - t184 * t47 + t60 * t244 + t52 * t245;
t226 = pkin(1) * t296;
t224 = -t133 + t248;
t223 = -t134 + t249;
t222 = pkin(4) * t171;
t19 = qJ(5) * t167 - t24;
t33 = qJ(5) * t185 - t275;
t34 = t185 * pkin(4) - t229;
t219 = qJD(2) * t233;
t218 = -t153 - t228;
t217 = -t83 - t289;
t216 = -t143 + t225;
t117 = -t182 * t269 + t183 * t265;
t214 = -t117 * qJ(5) + t140;
t64 = t198 + t97;
t112 = t179 * t120;
t75 = qJD(2) * t201 + t112;
t212 = -t182 * t75 + t184 * t64 - t99 * t244 - t92 * t245;
t206 = -t137 * qJ(5) + t170;
t205 = -pkin(5) * t44 - t228;
t204 = -t43 * pkin(5) + t227;
t203 = -t167 * t24 - t227;
t202 = t182 * t64 + t184 * t75 + t92 * t244 - t99 * t245;
t68 = t123 * t183 + t182 * t237 - t221;
t199 = t68 * qJ(5) - t117 * qJD(5) + t128;
t195 = -t16 * t85 + t205;
t5 = -t222 + t227;
t9 = pkin(4) * t44 + qJ(5) * t43 - qJD(5) * t211 + t119;
t3 = t44 * qJ(6) + t85 * qJD(6) + t9;
t192 = -qJD(1) * t219 + t204;
t10 = -qJ(5) * t247 + qJD(5) * t185 - t202;
t25 = -t43 - t270;
t190 = -t137 * t232 - t290;
t189 = t190 - t272;
t188 = t44 - t272;
t158 = 0.2e1 * t163;
t116 = t137 * t183;
t107 = t132 - t242;
t104 = -pkin(7) * t238 + t124;
t98 = -t180 * t241 + t112;
t82 = pkin(4) * t136 + t206;
t69 = t291 * t183 + t196;
t67 = -pkin(5) * t136 + t101;
t66 = pkin(5) * t137 + t100;
t57 = t284 * t136 + t206;
t53 = pkin(4) * t116 + t214;
t36 = pkin(4) * t211 + t271;
t35 = t284 * t116 + t214;
t28 = -pkin(5) * t116 - t33;
t26 = pkin(5) * t117 + qJ(6) * t185 + t34;
t22 = t211 * t284 + t271;
t18 = pkin(4) * t167 - t261;
t17 = pkin(4) * t69 + t199;
t13 = qJD(6) - t19 - t288;
t12 = t284 * t167 + t260;
t11 = -pkin(4) * t247 - t212;
t8 = t116 * qJD(6) + t284 * t69 + t199;
t7 = -pkin(5) * t69 - t10;
t6 = -t68 * pkin(5) + t185 * qJD(6) - t212 - t219;
t4 = -t163 - t218;
t2 = t205 + t293;
t1 = qJD(6) * t167 + t192;
t14 = [0, 0, 0, t185 * t298, t254 * t296, t262, -t263, 0, -pkin(7) * t262 + t183 * t226, pkin(7) * t263 + t185 * t226 (-qJD(1) * t97 - t72) * t185 + ((-pkin(7) * t133 + t143 * t179) * t185 + (t93 + (t107 + 0.2e1 * t242) * qJD(1)) * t183) * qJD(2) (qJD(1) * t98 + t73) * t185 + ((pkin(7) * t134 + t143 * t180) * t185 + (-t94 + (-t108 + 0.2e1 * t165) * qJD(1)) * t183) * qJD(2), t98 * t133 - t97 * t134 + (-t179 * t73 - t180 * t72) * t183 + (-t179 * t94 - t180 * t93 + (-t107 * t180 - t108 * t179) * qJD(1)) * t246, t72 * t107 + t73 * t108 + t93 * t97 + t94 * t98 + (t143 + t172) * t174, -t117 * t43 - t211 * t68, t116 * t43 - t117 * t44 - t211 * t69 + t68 * t85, t68 * t167 + t43 * t185 + (qJD(1) * t117 + t211) * t247, t69 * t167 + t44 * t185 + (-qJD(1) * t116 - t85) * t247 (-t167 - t252) * t247, -t212 * t167 + t227 * t185 + t128 * t85 + t140 * t44 + t119 * t116 + t102 * t69 + (qJD(1) * t229 - t23) * t247, t202 * t167 - t228 * t185 + t128 * t211 - t140 * t43 + t119 * t117 - t102 * t68 + (-t275 * qJD(1) - t24) * t247, t10 * t85 + t11 * t211 + t116 * t4 + t117 * t5 - t18 * t68 + t19 * t69 + t33 * t44 - t34 * t43, -t11 * t167 - t9 * t116 - t17 * t85 - t5 * t185 - t27 * t69 - t53 * t44 + (qJD(1) * t34 + t18) * t247, t10 * t167 - t9 * t117 - t17 * t211 + t4 * t185 + t27 * t68 + t53 * t43 + (-qJD(1) * t33 - t19) * t247, t10 * t19 + t11 * t18 + t17 * t27 + t33 * t4 + t34 * t5 + t53 * t9, t1 * t117 - t116 * t2 - t12 * t68 - t13 * t69 + t211 * t6 - t26 * t43 - t28 * t44 - t7 * t85, -t3 * t117 + t16 * t68 - t7 * t167 - t2 * t185 + t35 * t43 - t8 * t211 + (qJD(1) * t28 + t13) * t247, t1 * t185 + t3 * t116 + t16 * t69 + t6 * t167 + t35 * t44 + t8 * t85 + (-qJD(1) * t26 - t12) * t247, t1 * t26 + t12 * t6 + t13 * t7 + t16 * t8 + t2 * t28 + t3 * t35; 0, 0, 0, -t183 * t264, t254 * t187, 0, 0, 0, t187 * pkin(1) * t183, pkin(1) * t264 ((-qJ(3) * t249 - t93) * t183 + (-pkin(7) * t224 + t179 * t216 + t103) * t185) * qJD(1) ((-qJ(3) * t248 + t94) * t183 + (pkin(7) * t223 + t180 * t216 - t104) * t185) * qJD(1), t103 * t134 - t104 * t133 + (qJD(3) * t133 + t93 * t252 + t73) * t180 + (qJD(3) * t134 + t94 * t252 - t72) * t179, -t93 * t103 - t94 * t104 + (-t179 * t93 + t180 * t94) * qJD(3) + (-t72 * t179 + t73 * t180) * qJ(3) + (-t143 - t273) * t173, -t43 * t137 + t211 * t257, t43 * t136 - t137 * t44 - t211 * t256 - t257 * t85, -t257 * t167 + (qJD(2) * t137 - t211) * t253, t256 * t167 + (-qJD(2) * t136 + t85) * t253, t167 * t253, t119 * t136 - t127 * t85 + t170 * t44 + t295 * t167 + t256 * t102 + (t23 - t251) * t253, t119 * t137 - t127 * t211 - t170 * t43 + t297 * t167 + t257 * t102 + (t24 - t250) * t253, -t100 * t43 - t101 * t44 + t136 * t4 + t137 * t5 + t257 * t18 + t256 * t19 - t211 * t279 + t280 * t85, -t136 * t9 - t44 * t82 + t277 * t85 - t256 * t27 + t279 * t167 + (-t18 + t251) * t253, -t137 * t9 + t43 * t82 + t277 * t211 - t257 * t27 + t280 * t167 + (t19 + t250) * t253, t5 * t100 - t4 * t101 - t279 * t18 + t280 * t19 - t277 * t27 + t9 * t82, t1 * t137 + t257 * t12 - t256 * t13 - t136 * t2 + t211 * t282 + t281 * t85 - t43 * t66 - t44 * t67, -t137 * t3 + t43 * t57 - t278 * t211 + t281 * t167 - t257 * t16 + (qJD(2) * t67 - t13) * t253, t136 * t3 + t44 * t57 + t278 * t85 + t282 * t167 + t256 * t16 + (-qJD(2) * t66 + t12) * t253, t1 * t66 + t282 * t12 - t281 * t13 + t278 * t16 + t2 * t67 + t3 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t223 * t252, t224 * t252, -t133 ^ 2 - t134 ^ 2, -t133 * t94 + t134 * t93 + t166, 0, 0, 0, 0, 0, t188, -t191, t217, t190 + t272, t191, -t18 * t211 - t19 * t85 + t9, t217, t191, t188, -t12 * t211 + t13 * t85 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, t83 - t289, t25, t189, t171, -t102 * t211 + t203, t102 * t85 + t167 * t23 + t228, pkin(4) * t43 - t274 + (-t19 - t24) * t211 + (t18 + t261) * t85, t36 * t85 - t203 - 0.2e1 * t222 + t286, t167 * t261 + t211 * t36 - t27 * t85 + t158 + t218, -t5 * pkin(4) - t4 * qJ(5) - t18 * t24 + t19 * t261 - t27 * t36, -t274 + t284 * t43 + (t13 + t259) * t211 + (t12 - t260) * t85, -t167 * t208 + t211 * t22 - 0.2e1 * t153 + t158 + t195, -t287 - t22 * t85 + (-0.2e1 * qJD(6) - t15) * t167 + t284 * t298 - t204, t2 * qJ(5) - t1 * t284 + t12 * t259 + t13 * t260 - t16 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t231, t294, -t167 * t19 + t286 + t5, t25, t294, t231, t287 + (qJD(6) + t13) * t167 + t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, t171 + t285, -t162 - t289, -t12 * t167 + t195 + t293;];
tauc_reg  = t14;
