% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [5x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 18:44
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR14_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR14_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR14_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR14_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 18:43:37
% EndTime: 2024-09-27 18:43:40
% DurationCPUTime: 2.16s
% Computational Cost: add. (4356->265), mult. (9216->381), div. (0->0), fcn. (6791->10), ass. (0->212)
t162 = sin(qJ(4));
t163 = sin(qJ(3));
t166 = cos(qJ(4));
t167 = cos(qJ(3));
t293 = -t162 * t163 + t166 * t167;
t159 = sin(pkin(5));
t292 = t293 * t159;
t161 = sin(qJ(5));
t242 = qJD(5) * t161;
t156 = qJD(1) + qJD(2);
t160 = cos(pkin(5));
t257 = t160 * t156;
t143 = qJD(3) + t257;
t168 = cos(qJ(2));
t268 = pkin(1) * qJD(1);
t272 = t156 * pkin(2);
t138 = t168 * t268 + t272;
t255 = t160 * t167;
t124 = t138 * t255;
t164 = sin(qJ(2));
t236 = t164 * t268;
t260 = t156 * t159;
t127 = pkin(8) * t260 + t236;
t201 = pkin(9) * t260 + t127;
t81 = -t201 * t163 + t124;
t68 = t143 * pkin(3) + t81;
t256 = t160 * t163;
t231 = t138 * t256;
t82 = t201 * t167 + t231;
t78 = t166 * t82;
t195 = -t162 * t68 - t78;
t286 = t292 * t156;
t275 = t286 * pkin(10);
t28 = -t195 + t275;
t26 = t28 * t242;
t165 = cos(qJ(5));
t104 = t165 * t286;
t258 = t159 * t167;
t232 = t156 * t258;
t233 = t163 * t260;
t114 = -t162 * t232 - t166 * t233;
t66 = t161 * t114 + t104;
t259 = t156 * t167;
t110 = (-pkin(3) * t259 - t138) * t159;
t73 = -pkin(4) * t286 + t110;
t291 = -t73 * t66 + t26;
t245 = qJD(3) * t167;
t225 = t160 * t245;
t250 = t167 * t168;
t252 = t163 * t164;
t290 = (-t160 * t252 + t250) * t268 - pkin(2) * t225;
t235 = t159 * (-pkin(8) - pkin(9));
t208 = t163 * t235;
t289 = -qJD(3) * t208 + t290;
t176 = pkin(1) * (-t163 * t168 - t164 * t255);
t117 = qJD(1) * t176;
t238 = pkin(2) * t256;
t288 = (t167 * t235 - t238) * qJD(3) - t117;
t139 = qJD(4) + t143;
t136 = qJD(5) + t139;
t267 = pkin(1) * qJD(2);
t234 = qJD(1) * t267;
t207 = t164 * t234;
t193 = t160 * t207;
t206 = t168 * t234;
t249 = t138 * t225 + t167 * t206;
t48 = (-t201 * qJD(3) - t193) * t163 + t249;
t175 = qJD(2) * t176;
t174 = qJD(1) * t175;
t49 = -t82 * qJD(3) + t174;
t218 = -t162 * t48 + t166 * t49;
t172 = t195 * qJD(4) + t218;
t281 = qJD(3) + qJD(4);
t71 = t286 * t281;
t9 = -t71 * pkin(10) + t172;
t287 = (-t28 * t136 - t9) * t161 + t291;
t191 = t165 * t114 - t161 * t286;
t271 = t191 * t66;
t285 = t191 * t73;
t17 = t191 ^ 2 - t66 ^ 2;
t190 = t162 * t167 + t163 * t166;
t279 = t159 * t281;
t88 = t190 * t279;
t72 = t156 * t88;
t20 = qJD(5) * t104 + t114 * t242 - t161 * t72 + t165 * t71;
t13 = -t66 * t136 + t20;
t170 = t191 * qJD(5) - t161 * t71 - t165 * t72;
t14 = -t136 * t191 + t170;
t153 = t160 * pkin(3);
t107 = pkin(2) * t255 + t153 + t208;
t149 = pkin(9) * t258;
t179 = -pkin(8) * t258 - t238;
t116 = t149 - t179;
t243 = qJD(4) * t166;
t244 = qJD(4) * t162;
t284 = -t107 * t243 + t116 * t244 - t288 * t162 + t289 * t166;
t192 = -t162 * t107 - t166 * t116;
t283 = t192 * qJD(4) + t289 * t162 + t288 * t166;
t246 = qJD(3) * t163;
t282 = t159 * (-pkin(3) * t246 + t236);
t106 = t114 * pkin(10);
t76 = t162 * t82;
t217 = t166 * t68 - t76;
t27 = t106 + t217;
t241 = qJD(3) - t143;
t280 = t241 * t127 + t193;
t87 = t293 * t279;
t278 = t87 * pkin(10);
t277 = pkin(3) * t136;
t276 = pkin(3) * t167;
t274 = t114 * pkin(4);
t273 = t292 * pkin(4);
t269 = t166 * t81 - t76;
t266 = t165 * t28;
t265 = t110 * t114;
t264 = t114 * t286;
t155 = t159 ^ 2;
t263 = t155 * t156;
t262 = t155 * t163;
t261 = t155 * t167;
t253 = t162 * t165;
t151 = t168 * pkin(1) + pkin(2);
t248 = t151 * t225 + t250 * t267;
t226 = t159 * t246;
t145 = pkin(3) * t226;
t115 = t156 * t145 + t159 * t207;
t247 = t163 ^ 2 - t167 ^ 2;
t24 = t139 * pkin(4) + t27;
t196 = -t161 * t24 - t266;
t215 = -t162 * t49 - t166 * t48 - t68 * t243 + t82 * t244;
t8 = -t72 * pkin(10) - t215;
t227 = -t161 * t8 + t165 * t9;
t173 = t196 * qJD(5) + t227;
t122 = t190 * t159;
t85 = t165 * t122 + t161 * t292;
t30 = t85 * qJD(5) + t161 * t87 + t165 * t88;
t42 = t72 * pkin(4) + t115;
t84 = t161 * t122 - t165 * t292;
t240 = t173 * t160 + t73 * t30 + t42 * t84;
t239 = t110 * t88 - t115 * t292 + t172 * t160;
t237 = t164 * t267;
t230 = t151 * t256;
t224 = -pkin(4) * t136 - t24;
t74 = t88 * pkin(4) + t145;
t223 = -t138 - t272;
t142 = t164 * pkin(1) + t159 * pkin(8);
t222 = -pkin(9) * t159 - t142;
t221 = t160 * pkin(4) - t122 * pkin(10);
t220 = qJD(5) * t24 + t8;
t216 = -t162 * t81 - t78;
t214 = -((-qJD(3) * t127 - t193) * t163 + t249) * t160 + t207 * t262;
t213 = t143 + t257;
t212 = pkin(1) * t156 * t252;
t211 = pkin(3) * t233;
t210 = 0.2e1 * qJD(3) * t263;
t209 = t160 * t237;
t29 = -t84 * qJD(5) - t161 * t88 + t165 * t87;
t203 = -(t161 * t9 + t220 * t165 - t26) * t160 + t73 * t29 + t42 * t85;
t202 = t110 * t87 + t115 * t122 + t160 * t215;
t86 = t88 * pkin(10);
t200 = -qJD(5) * (t166 * t107 - t162 * t116 + t221) + t86 + t284;
t119 = t292 * pkin(10);
t199 = qJD(5) * (t119 - t192) + t278 - t283;
t132 = (-pkin(2) - t276) * t159;
t198 = (-qJD(2) + t156) * t268;
t197 = (-qJD(1) - t156) * t267;
t125 = (-t151 - t276) * t159;
t94 = t151 * t255 + t222 * t163 + t153;
t181 = -t167 * t142 - t230;
t95 = t149 - t181;
t194 = -t162 * t94 - t166 * t95;
t188 = qJD(3) * (-t151 * t156 - t138);
t187 = -t159 * t236 + t74;
t186 = t227 + t285;
t185 = t164 * t197;
t184 = t164 * t198;
t183 = -t110 * t286 + t215;
t69 = (t222 * qJD(3) - t209) * t163 + t248;
t70 = t175 + (t222 * t167 - t230) * qJD(3);
t180 = t162 * t70 + t166 * t69 + t94 * t243 - t95 * t244;
t171 = t194 * qJD(4) - t162 * t69 + t166 * t70;
t154 = t156 ^ 2;
t150 = t166 * pkin(3) + pkin(4);
t146 = t159 * t237;
t126 = t163 * t167 * t210;
t123 = t146 + t145;
t111 = t247 * t210;
t98 = t213 * t159 * t245;
t97 = t213 * t226;
t96 = t132 - t273;
t93 = t125 - t273;
t92 = t211 - t274;
t54 = t146 + t74;
t52 = ((-t167 * t127 - t231) * qJD(3) + t174) * t160;
t43 = t114 ^ 2 - t286 ^ 2;
t40 = -t281 * t260 * t190 - t114 * t139;
t39 = -t139 * t286 + t71;
t38 = t119 - t194;
t37 = -t162 * t95 + t166 * t94 + t221;
t36 = -t88 * t139 - t72 * t160;
t35 = t87 * t139 + t71 * t160;
t32 = t106 + t269;
t31 = t216 - t275;
t25 = -t114 * t87 + t71 * t122;
t16 = t171 - t278;
t15 = t180 - t86;
t10 = t114 * t88 - t122 * t72 + t286 * t87 + t292 * t71;
t6 = -t30 * t136 + t160 * t170;
t5 = t29 * t136 + t20 * t160;
t4 = -t191 * t29 + t20 * t85;
t3 = t170 * t85 + t191 * t30 - t20 * t84 + t29 * t66;
t1 = [0, 0, 0, 0, t185, t168 * t197, t126, -t111, t98, -t97, 0, (t181 * qJD(3) + t175) * t143 + t52 + (t163 * t188 + t167 * t185) * t155, -((-qJD(3) * t142 - t209) * t163 + t248) * t143 + (qJD(2) * t212 + t167 * t188) * t155 + t214, t25, t10, t35, t36, 0, -t123 * t286 + t125 * t72 + t171 * t139 + t239, -t123 * t114 + t125 * t71 - t180 * t139 + t202, t4, t3, t5, t6, 0, (-t161 * t15 + t165 * t16 + (-t161 * t37 - t165 * t38) * qJD(5)) * t136 - t54 * t66 - t93 * t170 + t240, -(t165 * t15 + t161 * t16 + (-t161 * t38 + t165 * t37) * qJD(5)) * t136 - t54 * t191 + t93 * t20 + t203; 0, 0, 0, 0, t184, t168 * t198, t126, -t111, t98, -t97, 0, -t117 * t143 + t52 + t184 * t261 + (t179 * t143 + t223 * t262) * qJD(3), (pkin(8) * t226 + t290) * t143 + (-qJD(1) * t212 + t223 * t245) * t155 + t214, t25, t10, t35, t36, 0, t132 * t72 + t283 * t139 + t282 * t286 + t239, t114 * t282 + t132 * t71 + t284 * t139 + t202, t4, t3, t5, t6, 0, -t96 * t170 - t187 * t66 + (t200 * t161 - t199 * t165) * t136 + t240, t96 * t20 - t187 * t191 + (t199 * t161 + t200 * t165) * t136 + t203; 0, 0, 0, 0, 0, 0, -t163 * t154 * t261, t247 * t155 * t154, t241 * t232, -t241 * t233, 0, -t280 * t167 + (-t206 + (-t241 * t160 + t263) * t138) * t163, t155 * t138 * t259 + t124 * t143 + t280 * t163 - t249, t264, t43, t39, t40, 0, -t216 * t139 + t286 * t211 + t265 + (-t78 + (-pkin(3) * t139 - t68) * t162) * qJD(4) + t218, t269 * t139 + (t114 * t233 - t139 * t243) * pkin(3) + t183, t271, t17, t13, t14, 0, -(-t161 * t32 + t165 * t31) * t136 + t92 * t66 + (-t161 * t166 - t253) * qJD(4) * t277 + ((-pkin(3) * t253 - t150 * t161) * t136 + t196) * qJD(5) + t186, t92 * t191 + (t31 * t136 - t9 - (-qJD(4) - qJD(5)) * t162 * t277) * t161 + ((-pkin(3) * t243 - qJD(5) * t150 + t32) * t136 - t220) * t165 + t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t264, t43, t39, t40, 0, -t195 * t139 + t172 + t265, t217 * t139 + t183, t271, t17, t13, t14, 0, -(-t161 * t27 - t266) * t136 - t66 * t274 + (t224 * t161 - t266) * qJD(5) + t186, -t191 * t274 + (t224 * qJD(5) + t27 * t136 - t8) * t165 + t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t271, t17, t13, t14, 0, -t196 * t136 + t173 + t285, (-t8 + (-qJD(5) + t136) * t24) * t165 + t287;];
tauc_reg = t1;
