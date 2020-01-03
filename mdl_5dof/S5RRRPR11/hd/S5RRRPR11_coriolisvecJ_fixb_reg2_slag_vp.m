% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR11_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR11_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:35:14
% EndTime: 2019-12-31 21:35:26
% DurationCPUTime: 4.32s
% Computational Cost: add. (4397->440), mult. (10483->579), div. (0->0), fcn. (6584->6), ass. (0->223)
t163 = sin(qJ(2));
t166 = cos(qJ(2));
t126 = -t166 * pkin(2) - t163 * pkin(7) - pkin(1);
t105 = t126 * qJD(1);
t232 = t166 * qJD(1);
t156 = pkin(6) * t232;
t133 = qJD(2) * pkin(7) + t156;
t162 = sin(qJ(3));
t165 = cos(qJ(3));
t62 = t165 * t105 - t162 * t133;
t246 = qJD(4) - t62;
t243 = qJD(1) * t163;
t222 = t165 * t243;
t242 = qJD(2) * t162;
t114 = t222 + t242;
t312 = -t232 + qJD(3);
t260 = t114 * t312;
t223 = t162 * t243;
t233 = t165 * qJD(2);
t112 = t223 - t233;
t263 = t112 * t312;
t231 = qJD(1) * qJD(2);
t211 = t166 * t231;
t239 = qJD(3) * t162;
t218 = t163 * t239;
t230 = qJD(2) * qJD(3);
t79 = qJD(1) * t218 + (-t211 - t230) * t165;
t240 = qJD(2) * t166;
t221 = t162 * t240;
t238 = qJD(3) * t165;
t176 = t163 * t238 + t221;
t80 = qJD(1) * t176 + t162 * t230;
t319 = (t80 + t260) * t162 + (t79 + t263) * t165;
t155 = pkin(6) * t243;
t276 = qJD(2) * pkin(2);
t132 = t155 - t276;
t183 = t114 * qJ(4) - t132;
t291 = pkin(3) + pkin(4);
t41 = -t291 * t112 + t183;
t161 = sin(qJ(5));
t164 = cos(qJ(5));
t60 = t112 * t161 + t114 * t164;
t318 = t41 * t60;
t191 = -t164 * t112 + t114 * t161;
t284 = t60 * t191;
t199 = pkin(2) * t163 - pkin(7) * t166;
t119 = t199 * qJD(1);
t103 = t162 * t119;
t153 = qJ(4) * t243;
t252 = t163 * t165;
t253 = t162 * t166;
t290 = pkin(7) - pkin(8);
t317 = t290 * t239 + t103 + t153 + (-pkin(6) * t252 + pkin(8) * t253) * qJD(1);
t135 = t290 * t165;
t289 = pkin(6) * t162;
t224 = -pkin(3) - t289;
t250 = t165 * t166;
t175 = -pkin(8) * t250 + (-pkin(4) + t224) * t163;
t251 = t165 * t119;
t316 = qJD(1) * t175 - qJD(3) * t135 - t251;
t313 = -t114 * pkin(8) + t246;
t122 = t199 * qJD(2);
t106 = qJD(1) * t122;
t151 = t163 * t231;
t202 = pkin(6) * t151;
t206 = t105 * t239 - t165 * t106 + t133 * t238 - t162 * t202;
t63 = t105 * t162 + t133 * t165;
t181 = t312 * t63 - t206;
t311 = -t191 ^ 2 + t60 ^ 2;
t140 = qJD(5) - t312;
t16 = qJD(5) * t60 - t161 * t79 - t164 * t80;
t303 = t60 * t140 - t16;
t267 = t80 * t165;
t268 = t79 * t162;
t310 = ((t112 * t162 - t114 * t165) * qJD(3) - t267 + t268) * t163 - (t112 * t165 + t114 * t162) * t240;
t205 = t112 + t233;
t309 = (t163 * t205 + t253 * t312) * qJD(1) - t312 * t239;
t159 = t163 ^ 2;
t187 = qJD(1) * t159 + t166 * t312;
t217 = t312 * t238;
t308 = (t112 * t163 + t162 * t187) * qJD(2) + t163 * t217 - t166 * t80;
t306 = -0.2e1 * t231;
t136 = t312 * qJD(4);
t142 = qJ(4) * t151;
t172 = -t105 * t238 - t162 * t106 + t133 * t239 + t165 * t202;
t20 = t136 + t142 - t172;
t10 = pkin(8) * t80 + t20;
t11 = t79 * pkin(8) - t291 * t151 + t206;
t30 = -t291 * t312 + t313;
t138 = t312 * qJ(4);
t47 = pkin(8) * t112 + t63;
t39 = t138 + t47;
t6 = t161 * t30 + t164 * t39;
t2 = -t6 * qJD(5) - t161 * t10 + t164 * t11;
t305 = t6 * t140 + t2;
t235 = qJD(5) * t164;
t236 = qJD(5) * t161;
t15 = -t112 * t235 + t114 * t236 - t161 * t80 + t164 * t79;
t304 = -t140 * t191 + t15;
t203 = pkin(3) * t151;
t27 = -t203 + t206;
t52 = t138 + t63;
t302 = -t312 * t52 + t27;
t301 = t80 - t260;
t299 = t162 * qJD(4) + t156;
t254 = t162 * t164;
t188 = t161 * t165 - t254;
t116 = t161 * t162 + t164 * t165;
t298 = t116 * qJD(5) - t161 * t239 - t164 * t238;
t1 = t164 * t10 + t161 * t11 + t30 * t235 - t39 * t236;
t297 = t191 * t41 - t1;
t256 = t162 * qJ(4);
t294 = -t291 * t165 - t256;
t292 = t114 ^ 2;
t288 = pkin(7) * t114;
t287 = pkin(7) * t312;
t5 = -t161 * t39 + t164 * t30;
t286 = t140 * t5;
t124 = t164 * qJ(4) - t161 * t291;
t283 = qJD(5) * t124 + t313 * t161 + t164 * t47;
t123 = -t161 * qJ(4) - t164 * t291;
t282 = -qJD(5) * t123 + t161 * t47 - t313 * t164;
t134 = t290 * t162;
t76 = t134 * t161 + t135 * t164;
t281 = qJD(5) * t76 - t317 * t161 + t316 * t164;
t75 = t134 * t164 - t135 * t161;
t280 = -qJD(5) * t75 + t316 * t161 + t317 * t164;
t265 = qJ(4) * t165;
t182 = -t291 * t162 + t265;
t279 = t312 * t182 + t299;
t180 = t166 * t116;
t278 = qJD(1) * t180 + t298;
t277 = t161 * t238 + t162 * t235 - t164 * t239 - t165 * t236 - t188 * t232;
t54 = pkin(3) * t112 - t183;
t275 = t114 * t54;
t173 = -pkin(6) * t211 - t79 * qJ(4) + t114 * qJD(4);
t24 = pkin(3) * t80 - t173;
t271 = t24 * t162;
t270 = t24 * t165;
t195 = pkin(3) * t162 - t265;
t266 = t312 * t195 - t299;
t264 = t112 * qJ(4);
t261 = t114 * t112;
t259 = t132 * t162;
t258 = t132 * t165;
t255 = t162 * t163;
t169 = qJD(1) ^ 2;
t249 = t166 * t169;
t168 = qJD(2) ^ 2;
t248 = t168 * t163;
t247 = t168 * t166;
t245 = t162 * t122 + t126 * t238;
t150 = pkin(6) * t250;
t87 = t162 * t126 + t150;
t244 = -t166 ^ 2 + t159;
t241 = qJD(2) * t163;
t237 = qJD(4) * t165;
t229 = t162 * t287;
t228 = t165 * t287;
t149 = pkin(6) * t253;
t227 = pkin(7) * t233;
t225 = t163 * t249;
t220 = t166 * t233;
t216 = t312 * t243;
t213 = t112 ^ 2 - t292;
t208 = pkin(1) * t306;
t86 = t126 * t165 - t149;
t207 = t140 ^ 2;
t204 = -t114 + t242;
t201 = t166 * t151;
t200 = t163 * t224;
t77 = -qJ(4) * t166 + t87;
t198 = qJD(3) * t150 - t165 * t122 + t126 * t239;
t197 = (qJD(3) * t112 - t79) * pkin(7);
t196 = t165 * pkin(3) + t256;
t158 = t166 * pkin(3);
t55 = t166 * pkin(4) + t149 + t158 + (-pkin(8) * t163 - t126) * t165;
t61 = pkin(8) * t255 + t77;
t21 = -t161 * t61 + t164 * t55;
t22 = t161 * t55 + t164 * t61;
t51 = -pkin(3) * t312 + t246;
t194 = -t162 * t52 + t165 * t51;
t193 = -t162 * t63 - t165 * t62;
t82 = -pkin(6) * t222 + t103;
t184 = pkin(6) + t195;
t178 = -pkin(6) + t182;
t174 = t162 * t263 - t267;
t44 = (-t163 * t233 - t166 * t239) * pkin(6) + t245;
t171 = t112 * t176 + t80 * t255;
t170 = t312 * t62 + t172;
t154 = qJ(4) * t241;
t125 = -pkin(2) - t196;
t109 = pkin(2) - t294;
t99 = t116 * t163;
t98 = t161 * t252 - t163 * t254;
t92 = t184 * t163;
t85 = (t312 - t232) * t241;
t81 = pkin(6) * t223 + t251;
t78 = t158 - t86;
t74 = t178 * t163;
t70 = pkin(7) * t267;
t69 = pkin(3) * t114 + t264;
t68 = qJD(1) * t200 - t251;
t67 = t153 + t82;
t49 = -t291 * t114 - t264;
t48 = -t79 + t263;
t45 = t241 * t289 - t198;
t43 = (qJD(3) * t196 - t237) * t163 + t184 * t240;
t42 = t217 + (t204 * t163 - t250 * t312) * qJD(1);
t40 = qJD(2) * t200 + t198;
t38 = -t166 * qJD(4) + t154 + t44;
t35 = qJD(2) * t180 + (qJD(3) - qJD(5)) * t163 * t188;
t34 = t161 * t220 + t298 * t163 - t164 * t221;
t31 = (t294 * qJD(3) + t237) * t163 + t178 * t240;
t29 = t165 * t260 - t268;
t26 = -t79 * t252 + (-t218 + t220) * t114;
t25 = t154 + (-pkin(6) * qJD(2) + pkin(8) * qJD(3)) * t252 + (-qJD(4) + (-pkin(6) * qJD(3) + pkin(8) * qJD(2)) * t162) * t166 + t245;
t23 = pkin(8) * t218 + qJD(2) * t175 + t198;
t19 = -t312 * t218 + t79 * t166 + (t114 * t163 + t165 * t187) * qJD(2);
t12 = -t291 * t80 + t173;
t4 = -qJD(5) * t22 - t161 * t25 + t164 * t23;
t3 = qJD(5) * t21 + t161 * t23 + t164 * t25;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t201, t244 * t306, t247, -0.2e1 * t201, -t248, 0, -pkin(6) * t247 + t163 * t208, pkin(6) * t248 + t166 * t208, 0, 0, t26, t310, t19, t171, -t308, t85, t45 * t312 + t206 * t166 + (pkin(6) * t80 + t132 * t238) * t163 + ((pkin(6) * t112 + t259) * t166 + (t62 + (t86 + t149) * qJD(1)) * t163) * qJD(2), -t44 * t312 - t172 * t166 + (-pkin(6) * t79 - t132 * t239) * t163 + ((pkin(6) * t114 + t258) * t166 + (-t63 + (-t87 + t150) * qJD(1)) * t163) * qJD(2), -t44 * t112 - t45 * t114 + t86 * t79 - t87 * t80 + t193 * t240 + (t162 * t172 + t165 * t206 + (t162 * t62 - t165 * t63) * qJD(3)) * t163, -t172 * t87 - t206 * t86 + t63 * t44 + t62 * t45 + (t132 + t155) * pkin(6) * t240, t26, t19, -t310, t85, t308, t171, t43 * t112 - t40 * t312 + t92 * t80 + (t242 * t54 + t27) * t166 + (t54 * t238 + t271 + (-qJD(1) * t78 - t51) * qJD(2)) * t163, -t38 * t112 + t40 * t114 - t77 * t80 - t78 * t79 + t194 * t240 + (-t162 * t20 + t165 * t27 + (-t162 * t51 - t165 * t52) * qJD(3)) * t163, -t43 * t114 + t38 * t312 + t92 * t79 + (-t233 * t54 - t20) * t166 + (t54 * t239 - t270 + (qJD(1) * t77 + t52) * qJD(2)) * t163, t20 * t77 + t24 * t92 + t27 * t78 + t38 * t52 + t40 * t51 + t43 * t54, -t15 * t99 + t35 * t60, t15 * t98 - t16 * t99 - t191 * t35 - t34 * t60, t35 * t140 - t15 * t166 + (-qJD(1) * t99 - t60) * t241, t16 * t98 + t191 * t34, -t34 * t140 - t16 * t166 + (qJD(1) * t98 + t191) * t241, (-t140 - t232) * t241, t12 * t98 + t4 * t140 + t74 * t16 + t2 * t166 + t31 * t191 + t41 * t34 + (-qJD(1) * t21 - t5) * t241, -t1 * t166 + t12 * t99 - t3 * t140 - t74 * t15 + t31 * t60 + t41 * t35 + (qJD(1) * t22 + t6) * t241, -t1 * t98 + t15 * t21 - t16 * t22 - t191 * t3 - t2 * t99 - t34 * t6 - t35 * t5 - t4 * t60, t1 * t22 + t12 * t74 + t2 * t21 + t3 * t6 + t31 * t41 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t225, t244 * t169, 0, t225, 0, 0, t169 * pkin(1) * t163, pkin(1) * t249, 0, 0, t29, -t319, t42, t174, t309, -t216, -pkin(2) * t80 - t81 * t312 + (-t228 + t259) * qJD(3) + ((-pkin(7) * t242 - t62) * t163 + (-pkin(6) * t205 - t259) * t166) * qJD(1), pkin(2) * t79 + t82 * t312 + (t229 + t258) * qJD(3) + ((t63 - t227) * t163 + (pkin(6) * t204 - t258) * t166) * qJD(1), t82 * t112 + t81 * t114 - t70 + (t62 * t232 - t172 + (-t62 + t288) * qJD(3)) * t165 + (t197 - t181) * t162, -t62 * t81 - t63 * t82 + (-t132 - t276) * t156 + (qJD(3) * t193 + t162 * t206 - t165 * t172) * pkin(7), t29, t42, t319, -t216, -t309, t174, t125 * t80 + t68 * t312 - t270 + t266 * t112 + (t162 * t54 - t228) * qJD(3) + (t163 * t51 + (-pkin(7) * t241 - t166 * t54) * t162) * qJD(1), t67 * t112 - t68 * t114 - t70 + (-t51 * t232 + t20 + (t51 + t288) * qJD(3)) * t165 + (t197 + t302) * t162, t125 * t79 - t67 * t312 - t271 - t266 * t114 + (-t165 * t54 - t229) * qJD(3) + (t54 * t250 + (-t52 + t227) * t163) * qJD(1), t24 * t125 - t51 * t68 - t52 * t67 + t266 * t54 + (qJD(3) * t194 + t27 * t162 + t20 * t165) * pkin(7), t15 * t188 - t278 * t60, t15 * t116 + t16 * t188 + t191 * t278 - t277 * t60, -t278 * t140 + (qJD(2) * t188 + t60) * t243, t16 * t116 + t191 * t277, -t277 * t140 + (qJD(2) * t116 - t191) * t243, t140 * t243, t109 * t16 + t12 * t116 + t279 * t191 + t277 * t41 - t281 * t140 + (-qJD(2) * t75 + t5) * t243, -t109 * t15 - t12 * t188 + t279 * t60 - t278 * t41 + t280 * t140 + (qJD(2) * t76 - t6) * t243, -t1 * t116 + t75 * t15 - t76 * t16 + t188 * t2 + t191 * t280 - t277 * t6 + t278 * t5 + t281 * t60, t1 * t76 + t12 * t109 + t2 * t75 + t279 * t41 - t280 * t6 - t281 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t261, -t213, t48, -t261, -t301, t151, -t114 * t132 + t181, t112 * t132 + t170, 0, 0, t261, t48, t213, t151, t301, -t261, -t112 * t69 + t181 + 0.2e1 * t203 - t275, pkin(3) * t79 - t80 * qJ(4) + (t52 - t63) * t114 + (t51 - t246) * t112, -t112 * t54 + t114 * t69 + 0.2e1 * t136 + 0.2e1 * t142 - t170, -t27 * pkin(3) + t20 * qJ(4) + t246 * t52 - t51 * t63 - t54 * t69, -t284, -t311, t304, t284, -t303, t151, -t123 * t151 - t283 * t140 - t191 * t49 - t2 + t318, t124 * t151 + t282 * t140 - t49 * t60 - t297, t123 * t15 - t124 * t16 + (t282 + t5) * t191 + (t283 - t6) * t60, t1 * t124 + t2 * t123 - t282 * t6 - t283 * t5 - t41 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151 + t261, t48, -t312 ^ 2 - t292, t275 + t302, 0, 0, 0, 0, 0, 0, -t114 * t191 - t151 * t164 - t161 * t207, -t114 * t60 + t151 * t161 - t164 * t207, t303 * t161 + t304 * t164, -t41 * t114 + t305 * t164 + (t1 - t286) * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t284, t311, -t304, -t284, t303, -t151, t305 - t318, t286 + t297, 0, 0;];
tauc_reg = t7;
