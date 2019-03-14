% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:21:22
% EndTime: 2019-03-09 02:21:32
% DurationCPUTime: 3.70s
% Computational Cost: add. (9423->389), mult. (22347->514), div. (0->0), fcn. (16647->10), ass. (0->200)
t175 = sin(pkin(11));
t177 = cos(pkin(11));
t181 = sin(qJ(4));
t279 = cos(qJ(4));
t150 = t279 * t175 + t181 * t177;
t142 = t150 * qJD(1);
t180 = sin(qJ(5));
t182 = cos(qJ(5));
t120 = -t182 * qJD(4) + t142 * t180;
t122 = qJD(4) * t180 + t142 * t182;
t179 = sin(qJ(6));
t278 = cos(qJ(6));
t194 = -t179 * t120 + t278 * t122;
t66 = t278 * t120 + t122 * t179;
t277 = t66 * t194;
t222 = t279 * t177;
t160 = qJD(1) * t222;
t236 = t181 * t175;
t219 = qJD(1) * t236;
t140 = -t160 + t219;
t239 = t179 * t180;
t193 = t278 * t182 - t239;
t288 = qJD(5) + qJD(6);
t215 = t278 * qJD(6);
t290 = t278 * qJD(5) + t215;
t256 = -t193 * t140 - t290 * t182 + t288 * t239;
t221 = t278 * t180;
t152 = t179 * t182 + t221;
t110 = t288 * t152;
t255 = t152 * t140 + t110;
t280 = -pkin(9) - pkin(8);
t224 = qJD(5) * t280;
t246 = t140 * t180;
t103 = pkin(4) * t142 + pkin(8) * t140;
t164 = sin(pkin(10)) * pkin(1) + qJ(3);
t156 = t164 * qJD(1);
t171 = t177 * qJD(2);
t268 = pkin(7) * qJD(1);
t118 = t171 + (-t156 - t268) * t175;
t131 = t175 * qJD(2) + t177 * t156;
t119 = t177 * t268 + t131;
t292 = t279 * t118 - t181 * t119;
t43 = t180 * t103 + t182 * t292;
t302 = -pkin(9) * t246 + t180 * t224 - t43;
t42 = t182 * t103 - t180 * t292;
t301 = pkin(5) * t142 + t42 + (pkin(9) * t140 - t224) * t182;
t159 = qJD(4) * t160;
t218 = qJD(4) * t236;
t133 = qJD(1) * t218 - t159;
t300 = -qJD(4) * qJD(5) + t133;
t232 = qJD(5) * t180;
t299 = t232 + t246;
t298 = t194 ^ 2 - t66 ^ 2;
t137 = qJD(5) + t140;
t132 = qJD(6) + t137;
t231 = qJD(5) * t182;
t225 = t142 * t231 - t300 * t180;
t230 = qJD(6) * t179;
t80 = t142 * t232 + t300 * t182;
t26 = t120 * t215 + t122 * t230 + t179 * t225 + t278 * t80;
t297 = t132 * t66 - t26;
t64 = t181 * t118 + t279 * t119;
t58 = qJD(4) * pkin(8) + t64;
t155 = -cos(pkin(10)) * pkin(1) - pkin(3) * t177 - pkin(2);
t138 = t155 * qJD(1) + qJD(3);
t79 = t140 * pkin(4) - pkin(8) * t142 + t138;
t38 = -t180 * t58 + t182 * t79;
t30 = -pkin(9) * t122 + t38;
t28 = pkin(5) * t137 + t30;
t39 = t180 * t79 + t182 * t58;
t31 = -pkin(9) * t120 + t39;
t145 = t150 * qJD(4);
t134 = qJD(1) * t145;
t195 = t222 - t236;
t188 = t195 * qJD(3);
t289 = qJD(1) * t188;
t50 = qJD(4) * t292 + t289;
t90 = pkin(4) * t134 + pkin(8) * t133;
t15 = -qJD(5) * t39 - t180 * t50 + t182 * t90;
t6 = pkin(5) * t134 + pkin(9) * t80 + t15;
t14 = t180 * t90 + t182 * t50 + t79 * t231 - t58 * t232;
t9 = -t225 * pkin(9) + t14;
t187 = -t179 * t6 - t28 * t215 + t31 * t230 - t278 * t9;
t57 = -qJD(4) * pkin(4) - t292;
t46 = t120 * pkin(5) + t57;
t296 = t46 * t66 + t187;
t200 = -t137 * t38 + t14;
t294 = t137 * t39 + t15;
t211 = t137 * t180;
t293 = t122 * t211;
t276 = pkin(7) + t164;
t146 = t276 * t175;
t147 = t276 * t177;
t291 = -t279 * t146 - t181 * t147;
t228 = t278 * t31;
t8 = t179 * t28 + t228;
t2 = -t8 * qJD(6) - t179 * t9 + t278 * t6;
t287 = -t46 * t194 + t2;
t27 = qJD(6) * t194 - t179 * t80 + t278 * t225;
t286 = t132 * t194 - t27;
t285 = -t256 * t132 + t152 * t134;
t125 = t182 * t134;
t217 = t150 * t232;
t144 = -qJD(4) * t222 + t218;
t244 = t144 * t182;
t191 = t217 + t244;
t284 = t150 * t125 - t137 * t191;
t283 = t145 * t120 - t195 * t225;
t282 = -t193 * t26 - t194 * t255;
t281 = t142 ^ 2;
t157 = t280 * t180;
t158 = t280 * t182;
t117 = t179 * t157 - t278 * t158;
t275 = t117 * qJD(6) + t302 * t179 + t301 * t278;
t116 = t278 * t157 + t179 * t158;
t274 = -t116 * qJD(6) + t301 * t179 - t302 * t278;
t35 = t110 * t150 + t193 * t144;
t99 = t193 * t150;
t273 = -t99 * t27 + t35 * t66;
t242 = t150 * t180;
t36 = -t144 * t221 - t179 * t217 - t230 * t242 + (-t144 * t179 + t290 * t150) * t182;
t98 = t152 * t150;
t272 = -t36 * t132 - t98 * t134;
t271 = t145 * t194 + t195 * t26;
t206 = t225 * t182;
t270 = t120 * t244 - t150 * t206;
t269 = t122 * t145 + t195 * t80;
t101 = -t181 * t146 + t279 * t147;
t92 = t182 * t101;
t93 = -pkin(4) * t195 - pkin(8) * t150 + t155;
t48 = t180 * t93 + t92;
t189 = t150 * qJD(3);
t51 = qJD(1) * t189 + t64 * qJD(4);
t267 = t291 * t51;
t264 = t142 * t66;
t263 = t195 * t51;
t261 = t179 * t31;
t260 = t180 * t51;
t259 = t51 * t182;
t258 = t194 * t142;
t257 = t80 * t180;
t76 = t180 * t225;
t254 = -t120 * t231 - t76;
t253 = t120 * t140;
t252 = t120 * t142;
t251 = t120 * t180;
t250 = t122 * t120;
t249 = t122 * t142;
t248 = t122 * t182;
t105 = t134 * t195;
t247 = t140 * t142;
t245 = t144 * t180;
t241 = t150 * t182;
t238 = t180 * t134;
t235 = -t134 * t150 + t140 * t144;
t234 = t175 ^ 2 + t177 ^ 2;
t233 = qJD(4) * t144;
t227 = t122 * t245;
t104 = pkin(4) * t145 + pkin(8) * t144;
t71 = t291 * qJD(4) + t188;
t213 = t182 * t104 - t180 * t71;
t47 = -t101 * t180 + t182 * t93;
t212 = qJD(1) * t234;
t210 = t137 * t182;
t209 = -t152 * t27 + t256 * t66;
t208 = t299 * pkin(5) - t64;
t207 = -t255 * t132 + t193 * t134;
t205 = t194 * t36 - t26 * t98;
t204 = -t132 * t35 + t134 * t99;
t203 = -t145 * t66 + t195 * t27;
t202 = t180 * t39 + t182 * t38;
t201 = t180 * t38 - t182 * t39;
t199 = (-t156 * t175 + t171) * t175 - t131 * t177;
t198 = t133 * t195 + t142 * t145;
t197 = -t299 * t137 + t125;
t40 = -pkin(5) * t195 - pkin(9) * t241 + t47;
t41 = -pkin(9) * t242 + t48;
t17 = -t179 * t41 + t278 * t40;
t18 = t179 * t40 + t278 * t41;
t192 = t150 * t231 - t245;
t21 = -t101 * t232 + t180 * t104 + t182 * t71 + t93 * t231;
t190 = -pkin(8) * t134 + t137 * t57;
t185 = -t137 * t192 - t150 * t238;
t184 = -t202 * qJD(5) + t14 * t182 - t15 * t180;
t72 = qJD(4) * t101 + t189;
t169 = -pkin(5) * t182 - pkin(4);
t139 = t140 ^ 2;
t136 = t145 * qJD(4);
t73 = pkin(5) * t242 - t291;
t44 = pkin(5) * t192 + t72;
t34 = t225 * pkin(5) + t51;
t22 = -t48 * qJD(5) + t213;
t19 = -pkin(9) * t192 + t21;
t16 = pkin(9) * t244 + pkin(5) * t145 + (-t92 + (pkin(9) * t150 - t93) * t180) * qJD(5) + t213;
t11 = t278 * t30 - t261;
t10 = -t179 * t30 - t228;
t7 = t278 * t28 - t261;
t4 = -t18 * qJD(6) + t278 * t16 - t179 * t19;
t3 = t17 * qJD(6) + t179 * t16 + t278 * t19;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(3) * t212 (t164 * t212 - t199) * qJD(3), -t133 * t150 - t142 * t144, -t198 + t235, -t233, t140 * t145 - t105, -t136, 0, -qJD(4) * t72 + t134 * t155 + t138 * t145, -qJD(4) * t71 - t133 * t155 - t138 * t144, -t101 * t134 + t133 * t291 - t140 * t71 + t142 * t72 + t144 * t292 - t145 * t64 + t150 * t51 + t195 * t50, t101 * t50 - t292 * t72 + t64 * t71 - t267, -t122 * t191 - t80 * t241, t227 + (t257 + (-t248 + t251) * qJD(5)) * t150 + t270, t269 + t284, t192 * t120 + t150 * t76, t185 - t283, t137 * t145 - t105, t22 * t137 + t47 * t134 - t15 * t195 + t38 * t145 + t72 * t120 - t291 * t225 - t57 * t245 + (t57 * t231 + t260) * t150, -t57 * t244 + t291 * t80 + t122 * t72 - t134 * t48 - t137 * t21 + t14 * t195 - t145 * t39 + (-t57 * t232 + t259) * t150, -t21 * t120 - t48 * t225 - t22 * t122 + t47 * t80 + t202 * t144 + (qJD(5) * t201 - t14 * t180 - t15 * t182) * t150, t14 * t48 + t15 * t47 + t21 * t39 + t22 * t38 + t57 * t72 - t267, -t194 * t35 - t26 * t99, -t205 + t273, t204 + t271, t27 * t98 + t36 * t66, t203 + t272, t132 * t145 - t105, t132 * t4 + t134 * t17 + t145 * t7 - t195 * t2 + t27 * t73 + t34 * t98 + t36 * t46 + t44 * t66, -t132 * t3 - t134 * t18 - t145 * t8 - t187 * t195 + t194 * t44 - t26 * t73 + t34 * t99 - t35 * t46, t17 * t26 - t18 * t27 + t187 * t98 - t194 * t4 - t2 * t99 - t3 * t66 + t35 * t7 - t36 * t8, t17 * t2 - t18 * t187 + t3 * t8 + t34 * t73 + t4 * t7 + t44 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, t233, t198 + t235, -t144 * t64 - t145 * t292 + t150 * t50 - t263, 0, 0, 0, 0, 0, 0, t185 + t283, t269 - t284, -t227 + (-t257 + (t248 + t251) * qJD(5)) * t150 + t270, t144 * t201 + t145 * t57 + t150 * t184 - t263, 0, 0, 0, 0, 0, 0, -t203 + t272, -t204 + t271, t205 + t273, t145 * t46 - t187 * t99 - t195 * t34 - t2 * t98 - t35 * t8 - t36 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234 * qJD(1) ^ 2, t199 * qJD(1), 0, 0, 0, 0, 0, 0, 0.2e1 * t142 * qJD(4), t159 + (-t140 - t219) * qJD(4), -t139 - t281, t140 * t64 + t142 * t292, 0, 0, 0, 0, 0, 0, t197 - t252, -t137 ^ 2 * t182 - t238 - t249 (t80 - t253) * t182 + t293 + t254, -t142 * t57 + t200 * t180 + t294 * t182, 0, 0, 0, 0, 0, 0, t207 - t264, -t258 - t285, t209 - t282, -t142 * t46 - t152 * t187 + t193 * t2 - t255 * t7 - t256 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t247, -t139 + t281, t159 + (t140 - t219) * qJD(4), -t247, 0, 0 -(qJD(3) + t138) * t142, t138 * t140 - t289, 0, 0, t122 * t210 - t257 (-t80 - t253) * t182 - t293 + t254, t137 * t210 + t238 - t249, t120 * t211 - t206, t197 + t252, -t137 * t142, -pkin(4) * t225 - t259 - t38 * t142 - t64 * t120 + (-pkin(8) * t231 - t42) * t137 + t190 * t180, pkin(4) * t80 - t122 * t64 + t142 * t39 + t260 + (pkin(8) * t232 + t43) * t137 + t190 * t182, t43 * t120 + t42 * t122 + ((qJD(5) * t122 - t225) * pkin(8) + t200) * t182 + ((qJD(5) * t120 - t80) * pkin(8) - t294) * t180, -pkin(4) * t51 + pkin(8) * t184 - t38 * t42 - t39 * t43 - t57 * t64, -t26 * t152 - t194 * t256, t209 + t282, -t258 + t285, -t193 * t27 + t255 * t66, t207 + t264, -t132 * t142, t116 * t134 - t275 * t132 - t142 * t7 + t169 * t27 - t193 * t34 + t208 * t66 + t255 * t46, -t117 * t134 + t274 * t132 + t142 * t8 + t152 * t34 - t169 * t26 + t194 * t208 - t256 * t46, t116 * t26 - t117 * t27 - t152 * t2 - t187 * t193 + t194 * t275 - t255 * t8 + t256 * t7 + t274 * t66, t116 * t2 - t117 * t187 + t169 * t34 + t208 * t46 - t274 * t8 - t275 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t250, -t120 ^ 2 + t122 ^ 2, t120 * t137 - t80, -t250, t122 * t137 - t225, t134, -t122 * t57 + t294, t120 * t57 - t200, 0, 0, t277, t298, t297, -t277, t286, t134, -t10 * t132 + (-t122 * t66 - t132 * t230 + t278 * t134) * pkin(5) + t287, t11 * t132 + (-t122 * t194 - t132 * t215 - t134 * t179) * pkin(5) + t296, t10 * t194 + t11 * t66 + t8 * t194 - t7 * t66 + (t278 * t26 - t179 * t27 + (t179 * t194 - t278 * t66) * qJD(6)) * pkin(5), -t7 * t10 - t8 * t11 + (t278 * t2 - t187 * t179 - t122 * t46 + (-t179 * t7 + t278 * t8) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t277, t298, t297, -t277, t286, t134, t8 * t132 + t287, t7 * t132 + t296, 0, 0;];
tauc_reg  = t1;