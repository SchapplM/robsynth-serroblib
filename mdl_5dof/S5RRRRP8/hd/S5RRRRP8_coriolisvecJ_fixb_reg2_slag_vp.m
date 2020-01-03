% Calculate inertial parameters regressor of coriolis joint torque vector for
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
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP8_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP8_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:02:09
% EndTime: 2019-12-31 22:02:21
% DurationCPUTime: 4.28s
% Computational Cost: add. (5992->430), mult. (15020->573), div. (0->0), fcn. (10129->6), ass. (0->215)
t201 = sin(qJ(2));
t203 = cos(qJ(2));
t214 = pkin(2) * t201 - pkin(7) * t203;
t154 = t214 * qJD(1);
t200 = sin(qJ(3));
t133 = t200 * t154;
t299 = -pkin(8) - pkin(7);
t241 = qJD(3) * t299;
t202 = cos(qJ(3));
t269 = t201 * t202;
t270 = t200 * t203;
t322 = t200 * t241 - t133 - (-pkin(6) * t269 - pkin(8) * t270) * qJD(1);
t259 = qJD(1) * t201;
t233 = t200 * t259;
t113 = pkin(6) * t233 + t202 * t154;
t268 = t202 * t203;
t211 = pkin(3) * t201 - pkin(8) * t268;
t321 = t211 * qJD(1) - t202 * t241 + t113;
t249 = t203 * qJD(1);
t184 = -qJD(3) + t249;
t173 = -qJD(4) + t184;
t256 = qJD(2) * t202;
t147 = -t233 + t256;
t238 = t202 * t259;
t258 = qJD(2) * t200;
t148 = t238 + t258;
t199 = sin(qJ(4));
t298 = cos(qJ(4));
t93 = -t298 * t147 + t148 * t199;
t285 = t173 * t93;
t253 = qJD(3) * t201;
t229 = qJD(1) * t253;
t248 = qJD(1) * qJD(2);
t230 = t203 * t248;
t313 = qJD(2) * qJD(3) + t230;
t112 = t200 * t229 - t313 * t202;
t231 = t298 * qJD(4);
t242 = t313 * t200 + t202 * t229;
t251 = qJD(4) * t199;
t44 = t298 * t112 - t147 * t231 + t148 * t251 + t199 * t242;
t320 = -t44 - t285;
t209 = t199 * t147 + t298 * t148;
t300 = t209 ^ 2;
t91 = t93 ^ 2;
t319 = -t91 + t300;
t318 = pkin(3) * t251;
t317 = t298 * pkin(3);
t316 = t93 * qJ(5);
t315 = t93 * t209;
t239 = t200 * t249;
t240 = t298 * t202;
t272 = t199 * t200;
t302 = qJD(3) + qJD(4);
t303 = t298 * qJD(3) + t231;
t263 = -t199 * t239 - t303 * t202 + t240 * t249 + t302 * t272;
t150 = t199 * t202 + t298 * t200;
t105 = t302 * t150;
t262 = -t150 * t249 + t105;
t283 = t209 * t173;
t45 = t209 * qJD(4) - t199 * t112 + t298 * t242;
t314 = -t45 - t283;
t255 = qJD(2) * t203;
t234 = t200 * t255;
t252 = qJD(3) * t202;
t236 = t201 * t252;
t312 = t234 + t236;
t193 = pkin(6) * t259;
t286 = qJD(2) * pkin(2);
t165 = t193 - t286;
t115 = -pkin(3) * t147 + t165;
t187 = t201 * t248;
t159 = -pkin(2) * t203 - pkin(7) * t201 - pkin(1);
t139 = t159 * qJD(1);
t194 = pkin(6) * t249;
t166 = qJD(2) * pkin(7) + t194;
t101 = t139 * t200 + t166 * t202;
t157 = t214 * qJD(2);
t140 = qJD(1) * t157;
t220 = pkin(6) * t187;
t60 = -qJD(3) * t101 + t202 * t140 + t200 * t220;
t42 = pkin(3) * t187 + pkin(8) * t112 + t60;
t254 = qJD(3) * t200;
t59 = t139 * t252 + t200 * t140 - t166 * t254 - t202 * t220;
t47 = -t242 * pkin(8) + t59;
t100 = t202 * t139 - t166 * t200;
t76 = -pkin(8) * t148 + t100;
t68 = -pkin(3) * t184 + t76;
t77 = pkin(8) * t147 + t101;
t226 = -t199 * t42 - t68 * t231 + t77 * t251 - t298 * t47;
t311 = t115 * t93 + t226;
t310 = -0.2e1 * t248;
t167 = t299 * t200;
t168 = t299 * t202;
t111 = t199 * t167 - t298 * t168;
t289 = t111 * qJD(4) + t322 * t199 + t321 * t298;
t288 = -t167 * t231 - t168 * t251 + t321 * t199 - t322 * t298;
t228 = -pkin(4) * t93 - qJD(5);
t63 = t115 - t228;
t309 = t63 * t209;
t186 = pkin(6) * t268;
t119 = t200 * t159 + t186;
t271 = t200 * t201;
t106 = -pkin(8) * t271 + t119;
t146 = t202 * t159;
t297 = pkin(6) * t200;
t99 = -pkin(8) * t269 + t146 + (-pkin(3) - t297) * t203;
t56 = t298 * t106 + t199 * t99;
t307 = t100 * t184 + t59;
t306 = t101 * t184 - t60;
t305 = qJ(5) * t209;
t216 = -t194 + (-t239 + t254) * pkin(3);
t182 = pkin(4) * t187;
t304 = -t209 * qJD(5) + t182;
t227 = -t199 * t47 + t298 * t42;
t73 = t298 * t77;
t31 = t199 * t68 + t73;
t6 = -t31 * qJD(4) + t227;
t301 = -t115 * t209 + t6;
t71 = t199 * t77;
t30 = t298 * t68 - t71;
t18 = t30 - t305;
t17 = -pkin(4) * t173 + t18;
t293 = t17 - t18;
t292 = pkin(4) * t259 - t263 * qJ(5) + t150 * qJD(5) + t289;
t149 = -t240 + t272;
t291 = -t262 * qJ(5) - qJD(5) * t149 - t288;
t221 = pkin(3) * t231;
t290 = -t199 * pkin(3) * t45 - t93 * t221;
t38 = t298 * t76 - t71;
t287 = t262 * pkin(4) + t216;
t284 = t44 * qJ(5);
t280 = t112 * t200;
t279 = t147 * t184;
t278 = t148 * t147;
t277 = t148 * t184;
t276 = t165 * t200;
t275 = t165 * t202;
t274 = t184 * t200;
t273 = t184 * t202;
t205 = qJD(1) ^ 2;
t267 = t203 * t205;
t204 = qJD(2) ^ 2;
t266 = t204 * t201;
t265 = t204 * t203;
t257 = qJD(2) * t201;
t261 = t202 * t157 + t257 * t297;
t158 = pkin(3) * t271 + t201 * pkin(6);
t197 = t201 ^ 2;
t260 = -t203 ^ 2 + t197;
t250 = t165 * qJD(3);
t246 = pkin(6) * t270;
t195 = pkin(6) * t255;
t244 = t199 * t271;
t243 = t201 * t267;
t116 = t312 * pkin(3) + t195;
t192 = -pkin(3) * t202 - pkin(2);
t237 = t200 * t253;
t235 = t173 * t259;
t37 = -t199 * t76 - t73;
t55 = -t106 * t199 + t298 * t99;
t225 = pkin(1) * t310;
t110 = t298 * t167 + t168 * t199;
t223 = -t147 + t256;
t222 = -t148 + t258;
t219 = t298 * t255;
t218 = t199 * t187;
t217 = t203 * t187;
t215 = t242 * t202;
t213 = -t100 * t202 - t101 * t200;
t212 = qJD(1) * t197 - t184 * t203;
t90 = t242 * pkin(3) + pkin(6) * t230;
t210 = qJ(5) * t45 + t226;
t54 = t211 * qJD(2) + (-t186 + (pkin(8) * t201 - t159) * t200) * qJD(3) + t261;
t74 = t200 * t157 + t159 * t252 + (-t201 * t256 - t203 * t254) * pkin(6);
t58 = -pkin(8) * t312 + t74;
t13 = -t106 * t251 + t199 * t54 + t99 * t231 + t298 * t58;
t27 = t45 * pkin(4) + t90;
t3 = -qJD(5) * t93 - t210;
t14 = -t56 * qJD(4) - t199 * t58 + t298 * t54;
t206 = t6 + t284;
t191 = pkin(4) + t317;
t143 = t173 * t221;
t128 = t201 * t240 - t244;
t127 = t150 * t201;
t122 = pkin(4) * t149 + t192;
t118 = t146 - t246;
t117 = (-t173 - t249) * t257;
t114 = -pkin(6) * t238 + t133;
t102 = pkin(4) * t127 + t158;
t81 = -qJ(5) * t149 + t111;
t80 = -qJ(5) * t150 + t110;
t75 = -t119 * qJD(3) + t261;
t70 = pkin(3) * t148 + pkin(4) * t209;
t62 = t200 * t219 - t199 * t237 - qJD(4) * t244 + (t199 * t255 + t303 * t201) * t202;
t61 = t105 * t201 + t199 * t234 - t202 * t219;
t49 = pkin(4) * t62 + t116;
t48 = -qJ(5) * t127 + t56;
t46 = -pkin(4) * t203 - qJ(5) * t128 + t55;
t29 = t262 * t173 + (-qJD(2) * t149 + t93) * t259;
t28 = t263 * t173 + (qJD(2) * t150 - t209) * t259;
t21 = t38 - t305;
t20 = t37 + t316;
t19 = t31 - t316;
t16 = t127 * t45 + t62 * t93;
t15 = -t128 * t44 - t209 * t61;
t12 = t149 * t45 + t262 * t93;
t11 = -t150 * t44 - t209 * t263;
t10 = t173 * t62 + t203 * t45 + (-qJD(1) * t127 - t93) * t257;
t9 = t173 * t61 + t203 * t44 + (qJD(1) * t128 + t209) * t257;
t8 = -qJ(5) * t62 - qJD(5) * t127 + t13;
t7 = pkin(4) * t257 + t61 * qJ(5) - t128 * qJD(5) + t14;
t4 = t127 * t44 - t128 * t45 - t209 * t62 + t61 * t93;
t2 = t206 + t304;
t1 = t149 * t44 - t150 * t45 - t209 * t262 + t263 * t93;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t217, t260 * t310, t265, -0.2e1 * t217, -t266, 0, -pkin(6) * t265 + t201 * t225, pkin(6) * t266 + t203 * t225, 0, 0, -t112 * t269 + (t202 * t255 - t237) * t148, (t147 * t202 - t148 * t200) * t255 + (-t215 + t280 + (-t147 * t200 - t148 * t202) * qJD(3)) * t201, t184 * t237 + t112 * t203 + (t148 * t201 + t202 * t212) * qJD(2), -t147 * t312 + t242 * t271, t184 * t236 + t242 * t203 + (t147 * t201 - t200 * t212) * qJD(2), (-t184 - t249) * t257, -t75 * t184 - t60 * t203 + (pkin(6) * t242 + t202 * t250) * t201 + ((-pkin(6) * t147 + t276) * t203 + (t100 + (t118 + t246) * qJD(1)) * t201) * qJD(2), t184 * t74 + t203 * t59 + (-pkin(6) * t112 - t200 * t250) * t201 + ((pkin(6) * t148 + t275) * t203 + (-t101 + (-t119 + t186) * qJD(1)) * t201) * qJD(2), -t75 * t148 + t118 * t112 + t74 * t147 - t119 * t242 + t213 * t255 + (-t59 * t200 - t60 * t202 + (t100 * t200 - t101 * t202) * qJD(3)) * t201, t100 * t75 + t101 * t74 + t118 * t60 + t119 * t59 + (t165 + t193) * t195, t15, t4, t9, t16, t10, t117, t115 * t62 + t116 * t93 + t127 * t90 - t14 * t173 + t158 * t45 - t203 * t6 + (qJD(1) * t55 + t30) * t257, -t115 * t61 + t116 * t209 + t128 * t90 + t13 * t173 - t158 * t44 - t203 * t226 + (-qJD(1) * t56 - t31) * t257, t127 * t226 - t128 * t6 - t13 * t93 - t14 * t209 + t30 * t61 - t31 * t62 + t44 * t55 - t45 * t56, t115 * t116 + t13 * t31 + t14 * t30 + t158 * t90 - t226 * t56 + t55 * t6, t15, t4, t9, t16, t10, t117, t102 * t45 + t127 * t27 - t173 * t7 - t2 * t203 + t49 * t93 + t62 * t63 + (qJD(1) * t46 + t17) * t257, -t102 * t44 + t128 * t27 + t173 * t8 + t203 * t3 + t49 * t209 - t61 * t63 + (-qJD(1) * t48 - t19) * t257, -t127 * t3 - t128 * t2 + t17 * t61 - t19 * t62 - t209 * t7 + t44 * t46 - t45 * t48 - t8 * t93, t102 * t27 + t17 * t7 + t19 * t8 + t2 * t46 + t3 * t48 + t49 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, t260 * t205, 0, t243, 0, 0, t205 * pkin(1) * t201, pkin(1) * t267, 0, 0, -t148 * t273 - t280, (-t112 - t279) * t202 + (-t242 + t277) * t200, -t184 * t252 + (t184 * t268 + t201 * t222) * qJD(1), t147 * t274 - t215, t184 * t254 + (-t184 * t270 + t201 * t223) * qJD(1), t184 * t259, -pkin(2) * t242 + t113 * t184 + (pkin(7) * t273 + t276) * qJD(3) + ((-pkin(7) * t258 - t100) * t201 + (-pkin(6) * t223 - t276) * t203) * qJD(1), pkin(2) * t112 - t114 * t184 + (-pkin(7) * t274 + t275) * qJD(3) + ((-pkin(7) * t256 + t101) * t201 + (pkin(6) * t222 - t275) * t203) * qJD(1), t113 * t148 - t114 * t147 + ((qJD(3) * t148 - t242) * pkin(7) + t307) * t202 + ((-qJD(3) * t147 - t112) * pkin(7) + t306) * t200, -t100 * t113 - t101 * t114 + (-t165 - t286) * t194 + (qJD(3) * t213 - t60 * t200 + t59 * t202) * pkin(7), t11, t1, t28, t12, t29, t235, t149 * t90 + t192 * t45 + t216 * t93 + t289 * t173 + t262 * t115 + (qJD(2) * t110 - t30) * t259, t150 * t90 - t192 * t44 + t216 * t209 - t288 * t173 - t263 * t115 + (-qJD(2) * t111 + t31) * t259, t110 * t44 - t111 * t45 + t149 * t226 - t150 * t6 + t209 * t289 - t262 * t31 + t263 * t30 + t288 * t93, t110 * t6 - t111 * t226 + t216 * t115 + t192 * t90 - t288 * t31 - t289 * t30, t11, t1, t28, t12, t29, t235, t122 * t45 + t149 * t27 + t287 * t93 + t262 * t63 + t292 * t173 + (qJD(2) * t80 - t17) * t259, -t122 * t44 + t150 * t27 + t287 * t209 - t263 * t63 + t291 * t173 + (-qJD(2) * t81 + t19) * t259, -t149 * t3 - t150 * t2 + t263 * t17 - t262 * t19 + t209 * t292 - t291 * t93 + t44 * t80 - t45 * t81, t122 * t27 - t292 * t17 + t291 * t19 + t2 * t80 + t287 * t63 + t3 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t278, -t147 ^ 2 + t148 ^ 2, -t112 + t279, t278, -t242 - t277, t187, -t148 * t165 - t306, -t147 * t165 - t307, 0, 0, t315, t319, t320, -t315, t314, t187, t37 * t173 + (-t148 * t93 + t173 * t251 + t298 * t187) * pkin(3) + t301, -t173 * t38 + t143 + (-t148 * t209 - t218) * pkin(3) + t311, t44 * t317 + t290 + (t31 + t37 + t318) * t209 + (-t30 + t38) * t93, -t30 * t37 - t31 * t38 + (t298 * t6 - t115 * t148 - t199 * t226 + (-t199 * t30 + t298 * t31) * qJD(4)) * pkin(3), t315, t319, t320, -t315, t314, t187, t191 * t187 + t284 + t20 * t173 - t309 - t70 * t93 + (-t73 + (pkin(3) * t173 - t68) * t199) * qJD(4) + t227 + t304, -pkin(3) * t218 - t173 * t21 - t209 * t70 + t63 * t93 + t143 - t3, t191 * t44 + t290 + (t19 + t20 + t318) * t209 + (-t17 + t21) * t93, -t17 * t20 - t19 * t21 + t2 * t191 - t63 * t70 + (t199 * t3 + (-t17 * t199 + t298 * t19) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t315, t319, t320, -t315, t314, t187, -t31 * t173 + t301, -t173 * t30 + t311, 0, 0, t315, t319, t320, -t315, t314, t187, -t19 * t173 + 0.2e1 * t182 + (t228 - t63) * t209 + t206, -pkin(4) * t300 - t173 * t18 + (qJD(5) + t63) * t93 + t210, pkin(4) * t44 - t293 * t93, t293 * t19 + (t2 - t309) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 - t283, -t44 + t285, -t91 - t300, t17 * t209 + t19 * t93 + t27;];
tauc_reg = t5;
