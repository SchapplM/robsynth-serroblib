% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRRP7
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
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:44
% EndTime: 2019-12-31 21:57:54
% DurationCPUTime: 3.36s
% Computational Cost: add. (5984->405), mult. (14880->505), div. (0->0), fcn. (10182->6), ass. (0->214)
t150 = cos(qJ(4));
t223 = qJD(4) * t150;
t148 = sin(qJ(3));
t149 = sin(qJ(2));
t226 = qJD(1) * t149;
t210 = t148 * t226;
t151 = cos(qJ(2));
t272 = cos(qJ(3));
t211 = t272 * t151;
t113 = -qJD(1) * t211 + t210;
t240 = t113 * t150;
t294 = t223 + t240;
t233 = t148 * t151;
t124 = t149 * t272 + t233;
t115 = t124 * qJD(1);
t147 = sin(qJ(4));
t221 = qJD(2) + qJD(3);
t100 = t115 * t147 - t150 * t221;
t102 = t150 * t115 + t147 * t221;
t243 = t102 * t150;
t246 = t100 * t147;
t206 = t272 * qJD(3);
t207 = t272 * qJD(2);
t165 = t151 * (t207 + t206);
t222 = qJD(1) * qJD(2);
t205 = t149 * t222;
t228 = qJD(3) * t210 + t148 * t205;
t292 = -qJD(1) * t165 + t228;
t157 = t147 * t292;
t62 = qJD(4) * t102 - t157;
t255 = t150 * t62;
t199 = qJD(4) * t221;
t224 = qJD(4) * t147;
t61 = t115 * t224 + (t292 - t199) * t150;
t258 = t147 * t61;
t195 = t151 * t207;
t234 = t148 * t149;
t97 = -t151 * t206 + t221 * t234 - t195;
t293 = t124 * (qJD(4) * (-t243 + t246) - t255 + t258) + (t100 * t150 + t102 * t147) * t97;
t111 = qJD(4) + t113;
t244 = t102 * t111;
t248 = t100 * t111;
t3 = (t62 + t244) * t147 + (t61 + t248) * t150;
t291 = -0.2e1 * t222;
t98 = t221 * t124;
t88 = t98 * qJD(1);
t266 = qJ(5) * t88;
t265 = qJD(2) * pkin(2);
t216 = t149 * t265;
t43 = t228 * pkin(8) + t88 * pkin(3) + (-pkin(8) * t165 + t216) * qJD(1);
t274 = -pkin(7) - pkin(6);
t131 = t274 * t149;
t126 = qJD(1) * t131;
t119 = t126 + t265;
t214 = qJD(2) * t274;
t196 = qJD(1) * t214;
t120 = t149 * t196;
t121 = t151 * t196;
t132 = t274 * t151;
t128 = qJD(1) * t132;
t225 = qJD(3) * t148;
t55 = t119 * t206 + t120 * t272 + t148 * t121 + t128 * t225;
t143 = -pkin(2) * t151 - pkin(1);
t130 = qJD(1) * t143;
t79 = pkin(3) * t113 - pkin(8) * t115 + t130;
t117 = t272 * t128;
t94 = t119 * t148 - t117;
t82 = pkin(8) * t221 + t94;
t9 = t147 * t43 + t150 * t55 + t223 * t79 - t224 * t82;
t2 = qJD(5) * t111 + t266 + t9;
t204 = t147 * t55 - t150 * t43 + t223 * t82 + t224 * t79;
t273 = pkin(4) * t88;
t6 = t204 - t273;
t290 = t6 * t147 + t2 * t150;
t254 = t150 * t88;
t289 = (t111 * t224 - t254) * pkin(8);
t108 = t272 * t121;
t202 = t148 * t120 - t108;
t56 = qJD(3) * t94 + t202;
t116 = t148 * t128;
t93 = t119 * t272 + t116;
t81 = -pkin(3) * t221 - t93;
t288 = t147 * t56 + t223 * t81;
t40 = t147 * t79 + t150 * t82;
t29 = qJ(5) * t111 + t40;
t36 = pkin(4) * t100 - qJ(5) * t102 + t81;
t286 = -t29 * t115 - t240 * t36;
t12 = t62 * pkin(4) + t61 * qJ(5) - t102 * qJD(5) + t56;
t285 = -t12 * t150 + t224 * t36;
t282 = t56 * t150 - t224 * t81;
t95 = t126 * t148 - t117;
t193 = pkin(2) * t225 - t95;
t241 = t113 * t147;
t280 = -qJD(5) * t147 - t294 * qJ(5) + (t224 + t241) * pkin(4);
t279 = t131 * t272 + t132 * t148;
t123 = -t211 + t234;
t256 = t147 * t97;
t173 = t124 * t223 - t256;
t238 = t124 * t147;
t278 = t100 * t98 + t111 * t173 + t123 * t62 + t238 * t88;
t105 = t131 * t148 - t132 * t272;
t92 = pkin(3) * t123 - pkin(8) * t124 + t143;
t268 = t105 * t150 + t147 * t92;
t60 = pkin(3) * t98 + pkin(8) * t97 + t216;
t127 = t149 * t214;
t67 = qJD(3) * t279 + t127 * t272 + t214 * t233;
t15 = -qJD(4) * t268 - t147 * t67 + t150 * t60;
t276 = t102 ^ 2;
t275 = t111 ^ 2;
t271 = pkin(2) * t148;
t270 = pkin(4) * t115;
t8 = t9 * t150;
t141 = pkin(8) + t271;
t198 = pkin(2) * t206;
t184 = t150 * t198;
t269 = -t100 * t184 - t141 * t255;
t89 = pkin(3) * t115 + pkin(8) * t113;
t53 = t147 * t89 + t150 * t93;
t218 = pkin(2) * t226;
t80 = t89 + t218;
t96 = t126 * t272 + t116;
t47 = t147 * t80 + t150 * t96;
t267 = pkin(2) * qJD(3);
t264 = t102 * t36;
t263 = t279 * t56;
t262 = t12 * t147;
t260 = t123 * t88;
t259 = t141 * t88;
t257 = t147 * t88;
t253 = t150 * t97;
t250 = t280 + t193;
t249 = -t94 + t280;
t247 = t100 * t141;
t245 = t102 * t100;
t242 = t111 * t115;
t239 = t115 * t113;
t237 = t124 * t150;
t236 = t130 * t115;
t153 = qJD(1) ^ 2;
t232 = t151 * t153;
t152 = qJD(2) ^ 2;
t231 = t152 * t149;
t230 = t152 * t151;
t39 = -t147 * t82 + t150 * t79;
t229 = qJD(5) - t39;
t227 = t149 ^ 2 - t151 ^ 2;
t220 = -t240 * t39 - t241 * t40 + t8;
t219 = t272 * pkin(2);
t215 = t149 * t232;
t213 = t147 * t272;
t212 = t150 * t272;
t209 = t141 * t223;
t208 = t100 ^ 2 - t276;
t203 = pkin(1) * t291;
t201 = t111 * t147;
t200 = t111 * t150;
t197 = t40 * t115 + t288;
t194 = t151 * t205;
t191 = pkin(4) * t150 + qJ(5) * t147;
t190 = pkin(4) * t147 - qJ(5) * t150;
t189 = t113 * t81 - t259;
t28 = -pkin(4) * t111 + t229;
t188 = t147 * t29 - t150 * t28;
t187 = t147 * t40 + t150 * t39;
t46 = -t147 * t96 + t150 * t80;
t52 = -t147 * t93 + t150 * t89;
t183 = -t100 * t115 - t254;
t63 = -t105 * t147 + t150 * t92;
t179 = t243 + t246;
t129 = -pkin(3) - t191;
t178 = t28 * t115 + t285;
t177 = -t39 * t115 - t282;
t176 = t111 * t40 - t204;
t175 = t223 * t36 + t262;
t172 = -t124 * t224 - t253;
t14 = -t105 * t224 + t147 * t60 + t150 * t67 + t223 * t92;
t171 = -t29 * t241 + t28 * t294 + t290;
t170 = (-t111 * t223 - t257) * pkin(8);
t169 = t100 * t201 - t255;
t167 = t102 * t198 - t141 * t61;
t164 = t141 * t224 - t184;
t163 = -t147 * t198 - t209;
t162 = t100 * t173 + t238 * t62;
t161 = -qJD(4) * t188 + t290;
t160 = -qJD(4) * t187 + t147 * t204 + t8;
t68 = qJD(3) * t105 + t148 * t127 - t195 * t274;
t156 = t130 * t113 - t55;
t154 = t115 * t223 + t147 * t199 - t157 - t244;
t142 = -t219 - pkin(3);
t122 = -t219 + t129;
t109 = t115 * qJ(5);
t72 = -t113 ^ 2 + t115 ^ 2;
t70 = t113 * t221 - t292;
t69 = pkin(4) * t102 + qJ(5) * t100;
t66 = t124 * t190 - t279;
t58 = pkin(8) * t255;
t45 = -pkin(4) * t123 - t63;
t44 = qJ(5) * t123 + t268;
t38 = -t52 - t270;
t37 = t109 + t53;
t35 = -t46 - t270;
t34 = t109 + t47;
t33 = t111 * t98 + t260;
t27 = -t61 + t248;
t23 = -t102 * t115 + t111 * t200 + t257;
t22 = -t147 * t275 - t183;
t21 = t111 * t201 + t183;
t18 = t102 * t200 - t258;
t17 = -t190 * t97 + (qJD(4) * t191 - qJD(5) * t150) * t124 + t68;
t16 = t102 * t172 - t237 * t61;
t13 = -pkin(4) * t98 - t15;
t11 = qJ(5) * t98 + qJD(5) * t123 + t14;
t7 = t102 * t98 + t111 * t172 - t123 * t61 + t237 * t88;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t194, t227 * t291, t230, -0.2e1 * t194, -t231, 0, -pkin(6) * t230 + t149 * t203, pkin(6) * t231 + t151 * t203, 0, 0, -t115 * t97 - t124 * t292, t97 * t113 - t115 * t98 + t123 * t292 - t124 * t88, -t97 * t221, t113 * t98 + t260, -t98 * t221, 0, t143 * t88 + t130 * t98 - t68 * t221 + (qJD(1) * t123 + t113) * t216, pkin(2) * t124 * t205 + t115 * t216 - t130 * t97 - t143 * t292 - t221 * t67, -t105 * t88 - t67 * t113 + t68 * t115 - t55 * t123 + t56 * t124 + t279 * t292 + t93 * t97 - t94 * t98, t105 * t55 + 0.2e1 * t130 * t216 + t67 * t94 - t68 * t93 - t263, t16, t293, t7, t162, -t278, t33, t100 * t68 + t111 * t15 - t123 * t204 + t124 * t288 - t256 * t81 - t279 * t62 + t39 * t98 + t63 * t88, t102 * t68 - t111 * t14 - t123 * t9 + t124 * t282 - t253 * t81 - t268 * t88 + t279 * t61 - t40 * t98, -t100 * t14 - t102 * t15 + t61 * t63 - t62 * t268 + t187 * t97 + (t204 * t150 - t147 * t9 + (t147 * t39 - t150 * t40) * qJD(4)) * t124, t14 * t40 + t15 * t39 - t204 * t63 + t268 * t9 + t68 * t81 - t263, t16, t7, -t293, t33, t278, t162, t100 * t17 - t111 * t13 - t123 * t6 + t124 * t175 - t256 * t36 - t28 * t98 - t45 * t88 + t62 * t66, -t100 * t11 + t102 * t13 - t44 * t62 - t45 * t61 + t188 * t97 + (-t147 * t2 + t150 * t6 + (-t147 * t28 - t150 * t29) * qJD(4)) * t124, -t102 * t17 + t11 * t111 + t123 * t2 + t124 * t285 + t253 * t36 + t29 * t98 + t44 * t88 + t61 * t66, t11 * t29 + t12 * t66 + t13 * t28 + t17 * t36 + t2 * t44 + t45 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t215, t227 * t153, 0, t215, 0, 0, t153 * pkin(1) * t149, pkin(1) * t232, 0, 0, t239, t72, t70, -t239, 0, 0, t128 * t206 + t108 - t113 * t218 - t236 + t95 * t221 + (-qJD(3) * t119 - t221 * t267 - t120) * t148, t96 * t221 + (-t115 * t226 - t206 * t221) * pkin(2) + t156, t292 * t219 - t88 * t271 + (t94 + t193) * t115 + (-t198 - t93 + t96) * t113, t93 * t95 - t94 * t96 + (-t130 * t226 - t272 * t56 + t148 * t55 + (-t148 * t93 + t272 * t94) * qJD(3)) * pkin(2), t18, -t3, t23, t169, t22, -t242, t142 * t62 + t189 * t147 + t193 * t100 + (-t46 + t163) * t111 + t177, -t142 * t61 + t189 * t150 + t193 * t102 + (t47 + t164) * t111 + t197, t47 * t100 + t46 * t102 + (t102 * t141 - t39) * t223 + (t204 + (-t40 + t247) * qJD(4) + t167) * t147 + t220 + t269, t56 * t142 - t39 * t46 - t40 * t47 - t81 * t95 + (t148 * t81 + t212 * t40 - t213 * t39) * t267 + t160 * t141, t18, t23, t3, -t242, t21, t169, t122 * t62 + (t113 * t36 - t259) * t147 + t250 * t100 + (t35 + t163) * t111 + t178, t34 * t100 + (-t35 + t209) * t102 + ((-t29 + t247) * qJD(4) + t167) * t147 + t171 + t269, -t262 + t122 * t61 + (-qJD(4) * t36 + t259) * t150 - t250 * t102 + (-t34 - t164) * t111 + t286, t12 * t122 - t28 * t35 - t29 * t34 + t250 * t36 + (t212 * t29 + t213 * t28) * t267 + t161 * t141; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t239, t72, t70, -t239, 0, 0, qJD(2) * t94 - t202 - t236, t221 * t93 + t156, 0, 0, t18, -t3, t23, t169, t22, -t242, -pkin(3) * t62 - t100 * t94 - t111 * t52 + t241 * t81 + t170 + t177, pkin(3) * t61 - t102 * t94 + t111 * t53 + t240 * t81 + t197 + t289, t100 * t53 + t102 * t52 - t58 + (-pkin(8) * t61 + t204) * t147 + (pkin(8) * t179 - t187) * qJD(4) + t220, -pkin(3) * t56 + pkin(8) * t160 - t39 * t52 - t40 * t53 - t81 * t94, t18, t23, t3, -t242, t21, t169, t100 * t249 + t111 * t38 + t129 * t62 + t241 * t36 + t170 + t178, -t29 * t224 + t100 * t37 - t102 * t38 - t58 + (qJD(4) * t179 - t258) * pkin(8) + t171, -t102 * t249 - t111 * t37 + t129 * t61 - t175 + t286 - t289, pkin(8) * t161 + t12 * t129 + t249 * t36 - t28 * t38 - t29 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t245, -t208, t27, -t245, -t154, t88, -t102 * t81 + t176, t100 * t81 + t111 * t39 - t9, 0, 0, t245, t27, t208, t88, t154, -t245, -t100 * t69 + t176 - t264 + 0.2e1 * t273, pkin(4) * t61 - qJ(5) * t62 + (t29 - t40) * t102 + (t28 - t229) * t100, 0.2e1 * t266 - t100 * t36 + t102 * t69 + (0.2e1 * qJD(5) - t39) * t111 + t9, -pkin(4) * t6 + qJ(5) * t2 + t229 * t29 - t28 * t40 - t36 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115 * t221 + t245, t27, -t275 - t276, -t111 * t29 + t264 + t6;];
tauc_reg = t1;
