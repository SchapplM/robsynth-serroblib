% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRRP6
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
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:34
% EndTime: 2019-12-31 21:54:42
% DurationCPUTime: 2.90s
% Computational Cost: add. (5946->392), mult. (14750->502), div. (0->0), fcn. (10147->6), ass. (0->211)
t162 = cos(qJ(2));
t270 = cos(qJ(3));
t205 = t270 * qJD(3);
t206 = t270 * qJD(2);
t174 = t162 * (t206 + t205);
t159 = sin(qJ(3));
t160 = sin(qJ(2));
t221 = qJD(1) * qJD(2);
t204 = t160 * t221;
t225 = qJD(1) * t160;
t208 = t159 * t225;
t227 = qJD(3) * t208 + t159 * t204;
t285 = -qJD(1) * t174 + t227;
t284 = -0.2e1 * t221;
t158 = sin(qJ(4));
t272 = -pkin(7) - pkin(6);
t210 = qJD(2) * t272;
t192 = qJD(1) * t210;
t126 = t162 * t192;
t112 = t270 * t126;
t125 = t160 * t192;
t198 = t159 * t125 - t112;
t141 = t272 * t162;
t134 = qJD(1) * t141;
t122 = t270 * t134;
t140 = t272 * t160;
t132 = qJD(1) * t140;
t263 = qJD(2) * pkin(2);
t124 = t132 + t263;
t96 = t159 * t124 - t122;
t55 = t96 * qJD(3) + t198;
t232 = t159 * t162;
t130 = t270 * t160 + t232;
t120 = t130 * qJD(1);
t161 = cos(qJ(4));
t220 = qJD(2) + qJD(3);
t107 = t161 * t120 + t158 * t220;
t167 = t158 * t285;
t61 = t107 * qJD(4) - t167;
t21 = t61 * pkin(4) + t55;
t222 = qJD(4) * t161;
t105 = t120 * t158 - t161 * t220;
t200 = t105 * pkin(4) + qJD(5);
t121 = t159 * t134;
t95 = t270 * t124 + t121;
t82 = -t220 * pkin(3) - t95;
t59 = t200 + t82;
t283 = t21 * t158 + t59 * t222;
t282 = t55 * t158 + t82 * t222;
t110 = t159 * t140 - t270 * t141;
t103 = t161 * t110;
t209 = t270 * t162;
t233 = t159 * t160;
t129 = -t209 + t233;
t153 = -pkin(2) * t162 - pkin(1);
t94 = pkin(3) * t129 - pkin(8) * t130 + t153;
t63 = t158 * t94 + t103;
t194 = pkin(2) * t205;
t216 = pkin(2) * t225;
t118 = -qJD(1) * t209 + t208;
t92 = pkin(3) * t120 + pkin(8) * t118;
t81 = t92 + t216;
t98 = t270 * t132 + t121;
t47 = t158 * t81 + t161 * t98;
t281 = -t161 * t194 + t47;
t223 = qJD(4) * t158;
t280 = t161 * t21 - t59 * t223;
t116 = qJD(4) + t118;
t246 = t105 * t116;
t195 = qJD(4) * t220;
t60 = t120 * t223 + (t285 - t195) * t161;
t279 = -t60 - t246;
t278 = t55 * t161 - t82 * t223;
t224 = qJD(3) * t159;
t97 = t132 * t159 - t122;
t189 = pkin(2) * t224 - t97;
t240 = t118 * t158;
t277 = (t223 + t240) * pkin(4);
t276 = t270 * t140 + t159 * t141;
t275 = -qJ(5) * t240 + t161 * qJD(5);
t139 = qJD(1) * t153;
t79 = pkin(3) * t118 - pkin(8) * t120 + t139;
t83 = t220 * pkin(8) + t96;
t39 = -t158 * t83 + t161 * t79;
t40 = t158 * t79 + t161 * t83;
t274 = -t158 * t39 + t161 * t40;
t273 = t107 ^ 2;
t102 = t220 * t130;
t91 = t102 * qJD(1);
t271 = pkin(4) * t91;
t269 = pkin(2) * t159;
t268 = t161 * pkin(4);
t267 = -qJ(5) - pkin(8);
t26 = -qJ(5) * t107 + t39;
t25 = pkin(4) * t116 + t26;
t266 = t25 - t26;
t150 = pkin(8) + t269;
t228 = -qJ(5) - t150;
t197 = qJD(4) * t228;
t265 = -t158 * t197 - t275 + t281;
t51 = t158 * t92 + t161 * t95;
t264 = pkin(2) * qJD(3);
t262 = t276 * t55;
t213 = t160 * t263;
t44 = t227 * pkin(8) + t91 * pkin(3) + (-pkin(8) * t174 + t213) * qJD(1);
t54 = t124 * t205 + t270 * t125 + t159 * t126 + t134 * t224;
t201 = -t158 * t44 - t161 * t54 - t79 * t222 + t83 * t223;
t10 = t201 * t161;
t261 = t129 * t91;
t259 = t158 * t60;
t258 = t158 * t91;
t255 = t161 * t61;
t254 = t161 * t91;
t155 = t161 * qJ(5);
t187 = t120 * pkin(4) + t118 * t155;
t46 = -t158 * t98 + t161 * t81;
t252 = (-t194 - qJD(5)) * t158 + t161 * t197 - t187 - t46;
t202 = qJD(4) * t267;
t251 = t158 * t202 + t275 - t51;
t50 = -t158 * t95 + t161 * t92;
t250 = -t158 * qJD(5) + t161 * t202 - t187 - t50;
t249 = t277 + t189;
t191 = t162 * t206;
t101 = -t162 * t205 + t220 * t233 - t191;
t248 = t101 * t158;
t247 = t101 * t161;
t245 = t105 * t158;
t244 = t107 * t105;
t243 = t107 * t116;
t242 = t107 * t161;
t241 = t116 * t120;
t239 = t118 * t161;
t238 = t120 * t118;
t237 = t130 * t158;
t236 = t130 * t161;
t235 = t139 * t120;
t164 = qJD(1) ^ 2;
t231 = t162 * t164;
t163 = qJD(2) ^ 2;
t230 = t163 * t160;
t229 = t163 * t162;
t226 = t160 ^ 2 - t162 ^ 2;
t219 = -t39 * t239 - t40 * t240 - t10;
t58 = pkin(3) * t102 + pkin(8) * t101 + t213;
t133 = t160 * t210;
t65 = t276 * qJD(3) + t270 * t133 + t210 * t232;
t218 = t158 * t58 + t161 * t65 + t94 * t222;
t217 = t270 * pkin(2);
t211 = t160 * t231;
t207 = t130 * t222;
t203 = -t158 * t65 + t161 * t58;
t62 = -t110 * t158 + t161 * t94;
t199 = pkin(1) * t284;
t196 = t116 * t161;
t151 = -t217 - pkin(3);
t193 = t40 * t120 + t282;
t190 = t162 * t204;
t188 = -t96 + t277;
t186 = t118 * t82 - t150 * t91;
t27 = -qJ(5) * t105 + t40;
t185 = t158 * t27 + t161 * t25;
t184 = t158 * t40 + t161 * t39;
t182 = qJ(5) * t101 - qJD(5) * t130;
t181 = -t39 * t120 - t278;
t180 = qJ(5) * t61 + t201;
t179 = t207 - t248;
t178 = -t130 * t223 - t247;
t177 = t120 * t27 + t59 * t239 + t283;
t175 = -t120 * t25 + t59 * t240 - t280;
t12 = -qJD(4) * t40 - t158 * t54 + t161 * t44;
t173 = -t184 * qJD(4) - t12 * t158;
t170 = qJ(5) * t60 + t12;
t1 = -qJD(5) * t107 + t170 + t271;
t4 = -qJD(5) * t105 - t180;
t172 = -t185 * qJD(4) - t1 * t158 + t4 * t161 - t25 * t239 - t27 * t240;
t171 = t173 - t10;
t66 = t110 * qJD(3) + t159 * t133 - t272 * t191;
t166 = t139 * t118 - t54;
t165 = -t120 * t222 - t158 * t195 + t167;
t152 = -pkin(3) - t268;
t138 = pkin(8) * t161 + t155;
t137 = t267 * t158;
t136 = t151 - t268;
t128 = t150 * t161 + t155;
t127 = t228 * t158;
t104 = t105 ^ 2;
t80 = pkin(4) * t237 - t276;
t71 = -t118 ^ 2 + t120 ^ 2;
t67 = t118 * t220 - t285;
t45 = -t104 + t273;
t43 = -qJ(5) * t237 + t63;
t37 = t102 * t116 + t261;
t36 = pkin(4) * t129 - t130 * t155 + t62;
t32 = t179 * pkin(4) + t66;
t31 = t165 + t243;
t30 = -t60 + t246;
t23 = -t107 * t120 + t116 * t196 + t258;
t22 = -t116 ^ 2 * t158 + t105 * t120 + t254;
t19 = t116 * t245 - t255;
t18 = t107 * t196 - t259;
t16 = t179 * t105 + t61 * t237;
t15 = t178 * t107 - t60 * t236;
t14 = -t63 * qJD(4) + t203;
t13 = -t110 * t223 + t218;
t9 = -t102 * t105 - t179 * t116 - t129 * t61 - t91 * t237;
t8 = t102 * t107 + t178 * t116 - t129 * t60 + t91 * t236;
t7 = -qJ(5) * t207 + (-qJD(4) * t110 + t182) * t158 + t218;
t6 = t279 * t161 + (-t61 - t243) * t158;
t5 = pkin(4) * t102 + t182 * t161 + (-t103 + (qJ(5) * t130 - t94) * t158) * qJD(4) + t203;
t2 = (t105 * t161 + t107 * t158) * t101 + (t259 - t255 + (-t242 + t245) * qJD(4)) * t130;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t190, t226 * t284, t229, -0.2e1 * t190, -t230, 0, -pkin(6) * t229 + t160 * t199, pkin(6) * t230 + t162 * t199, 0, 0, -t120 * t101 - t130 * t285, t101 * t118 - t120 * t102 + t129 * t285 - t130 * t91, -t101 * t220, t102 * t118 + t261, -t102 * t220, 0, t153 * t91 + t139 * t102 - t66 * t220 + (qJD(1) * t129 + t118) * t213, pkin(2) * t130 * t204 - t139 * t101 + t120 * t213 - t153 * t285 - t220 * t65, t95 * t101 - t96 * t102 - t110 * t91 - t65 * t118 + t66 * t120 - t54 * t129 + t55 * t130 + t276 * t285, t110 * t54 + 0.2e1 * t139 * t213 + t65 * t96 - t66 * t95 - t262, t15, t2, t8, t16, t9, t37, t102 * t39 + t105 * t66 + t116 * t14 + t12 * t129 + t282 * t130 - t82 * t248 - t276 * t61 + t62 * t91, -t102 * t40 + t107 * t66 - t116 * t13 + t129 * t201 + t278 * t130 - t82 * t247 + t276 * t60 - t63 * t91, -t105 * t13 - t107 * t14 + t60 * t62 - t61 * t63 + t184 * t101 + (-t274 * qJD(4) - t12 * t161 + t158 * t201) * t130, t12 * t62 + t13 * t40 + t14 * t39 - t201 * t63 + t66 * t82 - t262, t15, t2, t8, t16, t9, t37, t1 * t129 + t102 * t25 + t105 * t32 + t116 * t5 + t283 * t130 - t59 * t248 + t36 * t91 + t61 * t80, -t102 * t27 + t107 * t32 - t116 * t7 - t129 * t4 + t280 * t130 - t59 * t247 - t43 * t91 - t60 * t80, -t105 * t7 - t107 * t5 + t36 * t60 - t43 * t61 + t185 * t101 + (-t1 * t161 - t158 * t4 + (t158 * t25 - t161 * t27) * qJD(4)) * t130, t1 * t36 + t21 * t80 + t25 * t5 + t27 * t7 + t32 * t59 + t4 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211, t226 * t164, 0, t211, 0, 0, t164 * pkin(1) * t160, pkin(1) * t231, 0, 0, t238, t71, t67, -t238, 0, 0, t134 * t205 + t112 - t118 * t216 - t235 + t97 * t220 + (-qJD(3) * t124 - t220 * t264 - t125) * t159, t98 * t220 + (-t120 * t225 - t205 * t220) * pkin(2) + t166, t285 * t217 - t91 * t269 + (t96 + t189) * t120 + (-t194 - t95 + t98) * t118, t95 * t97 - t96 * t98 + (-t139 * t225 - t270 * t55 + t159 * t54 + (-t159 * t95 + t270 * t96) * qJD(3)) * pkin(2), t18, t6, t23, t19, t22, -t241, t151 * t61 + t186 * t158 + t189 * t105 + (-t150 * t222 - t158 * t194 - t46) * t116 + t181, -t151 * t60 + t186 * t161 + t189 * t107 + (t150 * t223 + t281) * t116 + t193, t47 * t105 + t46 * t107 + (-t105 * t194 - t150 * t61 + (t107 * t150 - t39) * qJD(4)) * t161 + (t107 * t194 - t150 * t60 - t12 + (t105 * t150 - t40) * qJD(4)) * t158 + t219, t55 * t151 - t39 * t46 - t40 * t47 - t82 * t97 + (t159 * t82 + t274 * t270) * t264 + t171 * t150, t18, t6, t23, t19, t22, -t241, t105 * t249 + t116 * t252 + t127 * t91 + t136 * t61 + t175, t249 * t107 + t265 * t116 - t128 * t91 - t136 * t60 + t177, t265 * t105 - t252 * t107 + t127 * t60 - t128 * t61 + t172, t1 * t127 + t128 * t4 + t136 * t21 + t249 * t59 + t252 * t25 - t265 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t238, t71, t67, -t238, 0, 0, t96 * qJD(2) - t198 - t235, t220 * t95 + t166, 0, 0, t18, t6, t23, t19, t22, -t241, t82 * t240 - pkin(3) * t61 - t105 * t96 - t116 * t50 + (-t116 * t222 - t258) * pkin(8) + t181, t82 * t239 + pkin(3) * t60 - t107 * t96 + t116 * t51 + (t116 * t223 - t254) * pkin(8) + t193, t105 * t51 + t107 * t50 + (-t259 - t255 + (t242 + t245) * qJD(4)) * pkin(8) + t173 + t219, -pkin(3) * t55 + pkin(8) * t171 - t39 * t50 - t40 * t51 - t82 * t96, t18, t6, t23, t19, t22, -t241, t105 * t188 + t116 * t250 + t137 * t91 + t152 * t61 + t175, t107 * t188 - t116 * t251 - t138 * t91 - t152 * t60 + t177, -t105 * t251 - t107 * t250 + t137 * t60 - t138 * t61 + t172, t1 * t137 + t138 * t4 + t152 * t21 + t188 * t59 + t25 * t250 + t251 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, t45, t30, -t244, t31, t91, -t107 * t82 + t116 * t40 + t12, t105 * t82 + t116 * t39 + t201, 0, 0, t244, t45, t30, -t244, t31, t91, 0.2e1 * t271 + t116 * t27 + (-t200 - t59) * t107 + t170, -pkin(4) * t273 + t116 * t26 + (qJD(5) + t59) * t105 + t180, pkin(4) * t60 - t266 * t105, t266 * t27 + (-t107 * t59 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165 + t243, t279, -t104 - t273, t27 * t105 + t25 * t107 + t21;];
tauc_reg = t3;
