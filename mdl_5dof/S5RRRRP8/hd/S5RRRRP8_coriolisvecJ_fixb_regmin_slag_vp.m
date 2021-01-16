% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:22
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:21:20
% EndTime: 2021-01-16 00:21:32
% DurationCPUTime: 3.24s
% Computational Cost: add. (4073->359), mult. (10134->496), div. (0->0), fcn. (6884->6), ass. (0->185)
t157 = sin(qJ(3));
t255 = pkin(7) + pkin(8);
t202 = qJD(3) * t255;
t158 = sin(qJ(2));
t159 = cos(qJ(3));
t226 = t158 * t159;
t160 = cos(qJ(2));
t227 = t157 * t160;
t178 = pkin(2) * t158 - pkin(7) * t160;
t113 = t178 * qJD(1);
t92 = t157 * t113;
t276 = -t157 * t202 - t92 - (-pkin(6) * t226 - pkin(8) * t227) * qJD(1);
t225 = t159 * t160;
t175 = pkin(3) * t158 - pkin(8) * t225;
t217 = qJD(1) * t158;
t197 = t157 * t217;
t236 = pkin(6) * t197 + t159 * t113;
t275 = t175 * qJD(1) + t159 * t202 + t236;
t208 = t159 * qJD(2);
t106 = -t197 + t208;
t209 = t157 * qJD(2);
t107 = t159 * t217 + t209;
t156 = sin(qJ(4));
t253 = cos(qJ(4));
t195 = t253 * qJD(4);
t213 = qJD(3) * t158;
t193 = qJD(1) * t213;
t206 = qJD(1) * qJD(2);
t194 = t160 * t206;
t269 = qJD(2) * qJD(3) + t194;
t203 = t269 * t157 + t159 * t193;
t211 = qJD(4) * t156;
t76 = -t157 * t193 + t269 * t159;
t21 = -t106 * t195 + t107 * t211 + t156 * t203 - t253 * t76;
t207 = t160 * qJD(1);
t141 = -qJD(3) + t207;
t130 = -qJD(4) + t141;
t59 = -t253 * t106 + t156 * t107;
t241 = t59 * t130;
t274 = -t21 - t241;
t172 = t156 * t106 + t253 * t107;
t256 = t172 ^ 2;
t57 = t59 ^ 2;
t273 = -t57 + t256;
t229 = t156 * t157;
t171 = t253 * t159 - t229;
t258 = qJD(3) + qJD(4);
t259 = t253 * qJD(3) + t195;
t245 = -t259 * t159 + t171 * t207 + t258 * t229;
t109 = t156 * t159 + t253 * t157;
t71 = t258 * t109;
t244 = -t109 * t207 + t71;
t272 = t172 * t59;
t271 = t59 * qJ(5);
t22 = qJD(4) * t172 + t156 * t76 + t253 * t203;
t240 = t172 * t130;
t270 = -t22 - t240;
t198 = t160 * t209;
t212 = qJD(3) * t159;
t199 = t158 * t212;
t268 = t198 + t199;
t144 = t158 * t206;
t181 = pkin(6) * t144;
t116 = t178 * qJD(2);
t99 = qJD(1) * t116;
t237 = -t157 * t181 - t159 * t99;
t151 = pkin(6) * t207;
t124 = qJD(2) * pkin(7) + t151;
t118 = -t160 * pkin(2) - t158 * pkin(7) - pkin(1);
t98 = t118 * qJD(1);
t243 = t157 * t98;
t67 = t159 * t124 + t243;
t169 = -t67 * qJD(3) - t237;
t19 = pkin(3) * t144 - t76 * pkin(8) + t169;
t214 = qJD(3) * t157;
t176 = -t124 * t214 + t157 * t99 + t98 * t212;
t168 = -t159 * t181 + t176;
t24 = -t203 * pkin(8) + t168;
t66 = -t157 * t124 + t159 * t98;
t43 = -t107 * pkin(8) + t66;
t37 = -t141 * pkin(3) + t43;
t44 = t106 * pkin(8) + t67;
t186 = -t156 * t19 - t37 * t195 + t44 * t211 - t253 * t24;
t123 = -qJD(2) * pkin(2) + pkin(6) * t217;
t77 = -t106 * pkin(3) + t123;
t267 = t77 * t59 + t186;
t266 = -0.2e1 * t206;
t192 = -t59 * pkin(4) - qJD(5);
t34 = -t192 + t77;
t265 = t34 * t172;
t264 = qJ(5) * t172;
t263 = -t151 + (-t157 * t207 + t214) * pkin(3);
t125 = t255 * t157;
t126 = t255 * t159;
t220 = -t156 * t125 + t253 * t126;
t139 = pkin(4) * t144;
t262 = -t172 * qJD(5) + t139;
t261 = t220 * qJD(4) + t276 * t156 + t275 * t253;
t260 = -t125 * t195 - t126 * t211 - t275 * t156 + t276 * t253;
t42 = t253 * t44;
t12 = t156 * t37 + t42;
t191 = -t156 * t24 + t253 * t19;
t167 = -qJD(4) * t12 + t191;
t257 = -t77 * t172 + t167;
t40 = t156 * t44;
t11 = t253 * t37 - t40;
t6 = t11 - t264;
t5 = -t130 * pkin(4) + t6;
t254 = t5 - t6;
t252 = pkin(3) * t130;
t153 = t158 * pkin(6);
t251 = pkin(4) * t217 - t245 * qJ(5) + t109 * qJD(5) + t261;
t250 = -t244 * qJ(5) + t171 * qJD(5) + t260;
t249 = t253 * t43 - t40;
t248 = t244 * pkin(4) + t263;
t105 = t159 * t118;
t65 = -pkin(8) * t226 + t105 + (-pkin(6) * t157 - pkin(3)) * t160;
t143 = pkin(6) * t225;
t219 = t157 * t118 + t143;
t228 = t157 * t158;
t72 = -pkin(8) * t228 + t219;
t246 = t156 * t65 + t253 * t72;
t242 = t21 * qJ(5);
t239 = t76 * t157;
t238 = t157 * t116 + t118 * t212;
t235 = t159 * t116 + t209 * t153;
t234 = t106 * t141;
t233 = t107 * t141;
t232 = t123 * t157;
t231 = t123 * t159;
t230 = t141 * t159;
t162 = qJD(1) ^ 2;
t224 = t160 * t162;
t161 = qJD(2) ^ 2;
t223 = t161 * t158;
t222 = t161 * t160;
t117 = pkin(3) * t228 + t153;
t154 = t158 ^ 2;
t218 = -t160 ^ 2 + t154;
t216 = qJD(2) * t158;
t215 = qJD(2) * t160;
t210 = t123 * qJD(3);
t78 = t268 * pkin(3) + pkin(6) * t215;
t149 = -t159 * pkin(3) - pkin(2);
t201 = t157 * t213;
t200 = t160 * t214;
t190 = -t156 * t43 - t42;
t188 = -t156 * t72 + t253 * t65;
t185 = pkin(1) * t266;
t184 = -t253 * t125 - t156 * t126;
t183 = -t106 + t208;
t182 = -t107 + t209;
t180 = t253 * t215;
t179 = t156 * t144;
t177 = qJD(1) * t154 - t141 * t160;
t56 = t203 * pkin(3) + pkin(6) * t194;
t174 = t22 * qJ(5) + t186;
t29 = t175 * qJD(2) + (-t143 + (pkin(8) * t158 - t118) * t157) * qJD(3) + t235;
t31 = -t268 * pkin(8) + (-t158 * t208 - t200) * pkin(6) + t238;
t173 = t156 * t29 + t65 * t195 - t72 * t211 + t253 * t31;
t10 = t22 * pkin(4) + t56;
t2 = -t59 * qJD(5) - t174;
t166 = -t246 * qJD(4) - t156 * t31 + t253 * t29;
t163 = t167 + t242;
t148 = t253 * pkin(3) + pkin(4);
t102 = t195 * t252;
t87 = t171 * t158;
t86 = t109 * t158;
t81 = -pkin(4) * t171 + t149;
t68 = t86 * pkin(4) + t117;
t48 = qJ(5) * t171 + t220;
t47 = -t109 * qJ(5) + t184;
t39 = t107 * pkin(3) + pkin(4) * t172;
t33 = t157 * t180 - t156 * t201 - t211 * t228 + (t156 * t215 + t259 * t158) * t159;
t32 = t156 * t198 + t71 * t158 - t159 * t180;
t26 = t33 * pkin(4) + t78;
t25 = -t86 * qJ(5) + t246;
t23 = -t160 * pkin(4) - t87 * qJ(5) + t188;
t9 = t249 - t264;
t8 = t190 + t271;
t7 = t12 - t271;
t4 = -t33 * qJ(5) - t86 * qJD(5) + t173;
t3 = pkin(4) * t216 + t32 * qJ(5) - t87 * qJD(5) + t166;
t1 = t163 + t262;
t13 = [0, 0, 0, 0.2e1 * t160 * t144, t218 * t266, t222, -t223, 0, -pkin(6) * t222 + t158 * t185, pkin(6) * t223 + t160 * t185, t76 * t226 + (t160 * t208 - t201) * t107, (t159 * t106 - t107 * t157) * t215 + (-t159 * t203 - t239 + (-t157 * t106 - t107 * t159) * qJD(3)) * t158, t141 * t201 - t76 * t160 + (t107 * t158 + t159 * t177) * qJD(2), t141 * t199 + t203 * t160 + (t106 * t158 - t157 * t177) * qJD(2), (-t141 - t207) * t216, -(-t118 * t214 + t235) * t141 + (pkin(6) * t203 + t159 * t210 + (t105 * qJD(1) + t66) * qJD(2)) * t158 + ((-pkin(6) * t106 + t232) * qJD(2) + (t243 + (pkin(6) * t141 + t124) * t159) * qJD(3) + t237) * t160, (-pkin(6) * t200 + t238) * t141 + t176 * t160 + (pkin(6) * t76 - t157 * t210) * t158 + ((pkin(6) * t107 + t231) * t160 + (-pkin(6) * t230 - qJD(1) * t219 - t67) * t158) * qJD(2), -t172 * t32 - t21 * t87, -t172 * t33 + t21 * t86 - t87 * t22 + t32 * t59, t32 * t130 + t21 * t160 + (qJD(1) * t87 + t172) * t216, t33 * t130 + t22 * t160 + (-qJD(1) * t86 - t59) * t216, (-t130 - t207) * t216, t11 * t216 + t117 * t22 - t130 * t166 + t144 * t188 - t160 * t167 + t77 * t33 + t56 * t86 + t78 * t59, t173 * t130 - t186 * t160 + t78 * t172 - t117 * t21 + t56 * t87 - t77 * t32 + (-t246 * qJD(1) - t12) * t216, -t1 * t160 + t10 * t86 - t3 * t130 + t68 * t22 + t26 * t59 + t34 * t33 + (qJD(1) * t23 + t5) * t216, t10 * t87 + t4 * t130 + t2 * t160 - t68 * t21 + t26 * t172 - t34 * t32 + (-qJD(1) * t25 - t7) * t216, -t1 * t87 - t172 * t3 - t2 * t86 + t23 * t21 - t25 * t22 + t5 * t32 - t7 * t33 - t4 * t59, t1 * t23 + t10 * t68 + t2 * t25 + t34 * t26 + t5 * t3 + t7 * t4; 0, 0, 0, -t158 * t224, t218 * t162, 0, 0, 0, t162 * pkin(1) * t158, pkin(1) * t224, -t107 * t230 + t239, (t76 - t234) * t159 + (-t203 + t233) * t157, -t141 * t212 + (t141 * t225 + t158 * t182) * qJD(1), t141 * t214 + (-t141 * t227 + t158 * t183) * qJD(1), t141 * t217, -pkin(2) * t203 + t236 * t141 + (pkin(7) * t230 + t232) * qJD(3) + ((-pkin(7) * t209 - t66) * t158 + (-pkin(6) * t183 - t232) * t160) * qJD(1), -pkin(2) * t76 - t92 * t141 + (-t157 * pkin(7) * t141 + t231) * qJD(3) + (-t123 * t225 + (-pkin(7) * t208 + t67) * t158 + (t141 * t226 + t160 * t182) * pkin(6)) * qJD(1), -t21 * t109 - t172 * t245, -t109 * t22 - t171 * t21 - t172 * t244 + t245 * t59, t245 * t130 + (qJD(2) * t109 - t172) * t217, t244 * t130 + (qJD(2) * t171 + t59) * t217, t130 * t217, -t11 * t217 + t261 * t130 + t144 * t184 + t149 * t22 - t171 * t56 + t244 * t77 + t263 * t59, t56 * t109 - t149 * t21 - t245 * t77 + t263 * t172 + t260 * t130 + (-qJD(2) * t220 + t12) * t217, -t10 * t171 + t81 * t22 + t248 * t59 + t244 * t34 + t251 * t130 + (qJD(2) * t47 - t5) * t217, t10 * t109 - t81 * t21 + t248 * t172 - t245 * t34 + t250 * t130 + (-qJD(2) * t48 + t7) * t217, -t1 * t109 + t171 * t2 + t172 * t251 + t47 * t21 - t48 * t22 - t244 * t7 + t245 * t5 - t250 * t59, t1 * t47 + t10 * t81 + t2 * t48 + t248 * t34 + t250 * t7 - t251 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107 * t106, -t106 ^ 2 + t107 ^ 2, t76 + t234, -t203 - t233, t144, -t123 * t107 - t67 * t141 + t169, -t123 * t106 - t66 * t141 - t168, t272, t273, t274, t270, t144, t190 * t130 + (-t107 * t59 + t130 * t211 + t253 * t144) * pkin(3) + t257, t102 - t249 * t130 + (-t107 * t172 - t179) * pkin(3) + t267, t148 * t144 + t242 + t8 * t130 - t265 - t39 * t59 + (-t42 + (-t37 + t252) * t156) * qJD(4) + t191 + t262, -pkin(3) * t179 - t9 * t130 - t172 * t39 + t34 * t59 + t102 - t2, t148 * t21 - t5 * t59 + t7 * t172 + t9 * t59 + t8 * t172 + (-t156 * t22 + (t156 * t172 - t253 * t59) * qJD(4)) * pkin(3), t1 * t148 - t34 * t39 - t5 * t8 - t7 * t9 + (t156 * t2 + (-t156 * t5 + t253 * t7) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t272, t273, t274, t270, t144, -t12 * t130 + t257, -t11 * t130 + t267, -t7 * t130 + 0.2e1 * t139 + (t192 - t34) * t172 + t163, -t256 * pkin(4) - t6 * t130 + (qJD(5) + t34) * t59 + t174, t21 * pkin(4) - t254 * t59, t254 * t7 + (t1 - t265) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22 - t240, -t21 + t241, -t57 - t256, t172 * t5 + t7 * t59 + t10;];
tauc_reg = t13;
