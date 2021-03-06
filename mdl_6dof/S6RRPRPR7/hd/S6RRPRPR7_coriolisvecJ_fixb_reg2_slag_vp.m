% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:47:43
% EndTime: 2019-03-09 10:47:58
% DurationCPUTime: 5.74s
% Computational Cost: add. (10623->476), mult. (24947->615), div. (0->0), fcn. (17296->8), ass. (0->252)
t192 = qJD(2) - qJD(4);
t199 = sin(qJ(6));
t202 = cos(qJ(6));
t200 = sin(qJ(4));
t201 = sin(qJ(2));
t203 = cos(qJ(4));
t204 = cos(qJ(2));
t141 = t200 * t201 + t203 * t204;
t125 = t141 * qJD(1);
t260 = qJD(1) * t204;
t261 = qJD(1) * t201;
t127 = -t200 * t260 + t203 * t261;
t197 = sin(pkin(10));
t282 = cos(pkin(10));
t219 = -t197 * t125 + t282 * t127;
t225 = t192 * t199 - t202 * t219;
t315 = t282 * t125 + t197 * t127;
t325 = qJD(6) + t315;
t339 = t199 * t325;
t345 = t225 * t339;
t304 = t225 * t219;
t217 = t141 * qJD(4);
t253 = qJD(1) * qJD(2);
t246 = t204 * t253;
t247 = t201 * t253;
t82 = qJD(1) * t217 - t200 * t247 - t203 * t246;
t256 = qJD(4) * t203;
t257 = qJD(4) * t200;
t258 = qJD(2) * t204;
t337 = t200 * t258 + t201 * t256 - t204 * t257;
t83 = t337 * qJD(1) - t203 * t247;
t48 = -t197 * t82 + t282 * t83;
t291 = t199 * t48;
t338 = t202 * t325;
t326 = -t325 * t338 - t291;
t344 = t304 - t326;
t221 = t197 * t83 + t282 * t82;
t254 = qJD(6) * t202;
t255 = qJD(6) * t199;
t27 = t192 * t254 + t202 * t221 + t219 * t255;
t28 = -qJD(6) * t225 - t199 * t221;
t63 = t202 * t192 + t199 * t219;
t300 = -t199 * t28 - t63 * t254;
t306 = t63 * t315;
t343 = (t27 + t306) * t202 - t300;
t292 = t199 * t27;
t342 = t225 * t338 + t292;
t280 = qJ(5) * t125;
t183 = pkin(7) * t261;
t148 = pkin(8) * t261 - t183;
t205 = -pkin(2) - pkin(3);
t250 = t205 * qJD(2);
t109 = qJD(3) + t250 - t148;
t184 = pkin(7) * t260;
t150 = -pkin(8) * t260 + t184;
t194 = qJD(2) * qJ(3);
t128 = t150 + t194;
t71 = t109 * t200 + t128 * t203;
t57 = t71 - t280;
t293 = t197 * t57;
t279 = qJ(5) * t127;
t70 = t203 * t109 - t128 * t200;
t56 = t70 - t279;
t53 = -pkin(4) * t192 + t56;
t29 = t282 * t53 - t293;
t25 = t192 * pkin(5) - t29;
t341 = t25 * t325;
t285 = t315 * t192;
t336 = t221 + t285;
t302 = t315 ^ 2;
t303 = t219 ^ 2;
t335 = t302 - t303;
t235 = t201 * t250;
t187 = t201 * qJD(3);
t264 = qJ(3) * t246 + qJD(1) * t187;
t97 = qJD(1) * t235 + t264;
t52 = t83 * pkin(4) + t97;
t14 = t48 * pkin(5) + pkin(9) * t221 + t52;
t259 = qJD(2) * t201;
t313 = pkin(7) - pkin(8);
t149 = t313 * t259;
t193 = qJD(2) * qJD(3);
t114 = -qJD(1) * t149 + t193;
t173 = pkin(7) * t246;
t134 = -pkin(8) * t246 + t173;
t44 = -t71 * qJD(4) - t114 * t200 + t203 * t134;
t209 = qJ(5) * t82 - qJD(5) * t127 + t44;
t43 = t109 * t256 + t203 * t114 - t128 * t257 + t200 * t134;
t23 = -qJ(5) * t83 - qJD(5) * t125 + t43;
t6 = t197 * t209 + t282 * t23;
t54 = t282 * t57;
t30 = t197 * t53 + t54;
t26 = -pkin(9) * t192 + t30;
t129 = -qJD(1) * pkin(1) - pkin(2) * t260 - qJ(3) * t261;
t108 = pkin(3) * t260 - t129;
t75 = pkin(4) * t125 + qJD(5) + t108;
t33 = pkin(5) * t315 - pkin(9) * t219 + t75;
t8 = t199 * t33 + t202 * t26;
t2 = -qJD(6) * t8 + t202 * t14 - t199 * t6;
t334 = t325 * t8 + t2;
t226 = t199 * t26 - t202 * t33;
t1 = -qJD(6) * t226 + t14 * t199 + t202 * t6;
t333 = t226 * t325 + t1;
t332 = t219 * t8;
t331 = t29 * t315;
t307 = t63 * t219;
t330 = t325 * t219;
t294 = t192 * t219;
t329 = -t48 - t294;
t317 = t125 * t192 + t82;
t328 = t219 * t226;
t327 = t315 * t219;
t324 = -t315 * t8 - t2;
t5 = t197 * t23 - t282 * t209;
t323 = t219 * t75 + t5;
t322 = -t315 * t75 + t6;
t321 = t226 * t315 + t1;
t320 = pkin(5) * t219 + pkin(9) * t315;
t318 = -t192 * t71 + t44;
t316 = t127 * t192 + t83;
t161 = t313 * t201;
t162 = t313 * t204;
t101 = t200 * t161 + t203 * t162;
t153 = -qJ(3) * t200 + t203 * t205;
t100 = t203 * t161 - t162 * t200;
t142 = -t200 * t204 + t201 * t203;
t222 = -qJ(5) * t142 + t100;
t73 = -qJ(5) * t141 + t101;
t41 = t197 * t73 - t282 * t222;
t312 = t41 * t5;
t87 = -t197 * t141 + t282 * t142;
t311 = t5 * t87;
t310 = pkin(4) * t127;
t218 = -t197 * t200 + t282 * t203;
t309 = t218 * t5;
t86 = t282 * t141 + t142 * t197;
t308 = t48 * t86;
t305 = t225 * t63;
t301 = t87 * t48;
t112 = qJD(3) * t203 + qJD(4) * t153;
t154 = qJ(3) * t203 + t200 * t205;
t113 = -qJD(3) * t200 - qJD(4) * t154;
t90 = -t148 * t200 + t150 * t203;
t215 = t90 - t280;
t91 = t203 * t148 + t200 * t150;
t66 = t91 + t279;
t299 = (-t113 + t215) * t282 + (t112 - t66) * t197;
t38 = t197 * t215 + t282 * t66;
t68 = t282 * t112 + t197 * t113;
t298 = -t68 + t38;
t297 = qJD(2) * pkin(2);
t296 = t192 * t70;
t290 = t199 * t63;
t289 = t199 * t225;
t288 = t202 * t28;
t45 = t202 * t48;
t287 = t202 * t63;
t286 = t202 * t225;
t284 = t112 - t91;
t283 = t113 - t90;
t278 = qJD(6) * t325;
t277 = t108 * t127;
t275 = t127 * t125;
t207 = qJD(1) ^ 2;
t272 = t204 * t207;
t206 = qJD(2) ^ 2;
t271 = t206 * t201;
t188 = t206 * t204;
t139 = t197 * t203 + t282 * t200;
t268 = t192 * t139;
t122 = t218 * qJD(2);
t123 = t218 * qJD(4);
t267 = -t123 + t122;
t147 = -pkin(4) + t153;
t93 = t197 * t147 + t282 * t154;
t263 = qJ(3) * t258 + t187;
t195 = t201 ^ 2;
t262 = t204 ^ 2 - t195;
t156 = -t204 * pkin(2) - t201 * qJ(3) - pkin(1);
t252 = t87 * t255;
t251 = t87 * t254;
t248 = t125 ^ 2 - t127 ^ 2;
t240 = -0.2e1 * pkin(1) * t253;
t239 = qJD(3) - t297;
t136 = t204 * pkin(3) - t156;
t237 = qJD(1) * t156 + t129;
t236 = t192 ^ 2;
t89 = -pkin(9) + t93;
t234 = t89 * t278 - t5;
t233 = t201 * t246;
t176 = pkin(4) * t197 + pkin(9);
t232 = t176 * t278 + t5;
t94 = -t203 * t259 + t337;
t95 = qJD(2) * t141 - t217;
t50 = -t197 * t94 + t282 * t95;
t231 = t25 * t50 + t311;
t230 = -t199 * t8 + t202 * t226;
t229 = -t199 * t226 - t202 * t8;
t228 = t325 * t50 + t301;
t96 = pkin(4) * t141 + t136;
t40 = pkin(5) * t86 - pkin(9) * t87 + t96;
t42 = t197 * t222 + t282 * t73;
t18 = -t199 * t42 + t202 * t40;
t19 = t199 * t40 + t202 * t42;
t224 = t45 + (-t199 * t315 - t255) * t325;
t180 = qJ(3) * t260;
t117 = t205 * t261 + t180;
t107 = pkin(2) * t247 - t264;
t118 = pkin(2) * t259 - t263;
t223 = -pkin(7) * t206 - qJD(1) * t118 - t107;
t220 = -t176 * t48 + t341;
t92 = t282 * t147 - t197 * t154;
t103 = t235 + t263;
t151 = qJD(2) * t162;
t58 = -t203 * t149 + t200 * t151 + t161 * t256 - t162 * t257;
t84 = t117 - t310;
t213 = -t325 * t68 - t48 * t89 - t341;
t61 = pkin(4) * t94 + t103;
t212 = -t108 * t125 + t43;
t211 = qJD(6) * t230 + t1 * t202 - t199 * t2;
t59 = -qJD(4) * t101 + t149 * t200 + t151 * t203;
t152 = -pkin(7) * t247 + t193;
t155 = t183 + t239;
t158 = t184 + t194;
t210 = t152 * t204 + (t155 * t204 + (-t158 + t184) * t201) * qJD(2);
t208 = -qJ(5) * t95 - qJD(5) * t142 + t59;
t177 = -t282 * pkin(4) - pkin(5);
t170 = t201 * t272;
t160 = -0.2e1 * t233;
t159 = 0.2e1 * t233;
t157 = t262 * t207;
t146 = pkin(2) * t261 - t180;
t140 = t262 * t253;
t99 = t122 * t202 + t199 * t261;
t98 = -t122 * t199 + t202 * t261;
t88 = pkin(5) - t92;
t49 = t197 * t95 + t282 * t94;
t39 = t310 + t320;
t36 = -qJ(5) * t94 - qJD(5) * t141 + t58;
t34 = -t320 + t84;
t32 = t282 * t56 - t293;
t31 = t197 * t56 + t54;
t17 = pkin(5) * t49 - pkin(9) * t50 + t61;
t16 = t197 * t208 + t282 * t36;
t15 = t197 * t36 - t282 * t208;
t13 = t199 * t34 + t202 * t38;
t12 = -t199 * t38 + t202 * t34;
t10 = t199 * t39 + t202 * t32;
t9 = -t199 * t32 + t202 * t39;
t4 = -qJD(6) * t19 - t16 * t199 + t17 * t202;
t3 = qJD(6) * t18 + t16 * t202 + t17 * t199;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, 0.2e1 * t140, t188, t160, -t271, 0, -pkin(7) * t188 + t201 * t240, pkin(7) * t271 + t204 * t240, 0, 0, t159, t188, -0.2e1 * t140, 0, t271, t160, t204 * t223 + t237 * t259, t210, t201 * t223 - t237 * t258, pkin(7) * t210 + t107 * t156 + t118 * t129, t127 * t95 - t142 * t82, -t125 * t95 - t127 * t94 + t141 * t82 - t142 * t83, -t95 * t192, t125 * t94 + t141 * t83, t94 * t192, 0, t103 * t125 + t108 * t94 + t136 * t83 + t141 * t97 - t192 * t59, t103 * t127 + t108 * t95 - t136 * t82 + t142 * t97 + t192 * t58, t100 * t82 - t101 * t83 - t125 * t58 - t127 * t59 - t141 * t43 - t142 * t44 - t70 * t95 - t71 * t94, t100 * t44 + t101 * t43 + t103 * t108 + t136 * t97 + t58 * t71 + t59 * t70, t219 * t50 - t221 * t87, -t219 * t49 + t221 * t86 - t315 * t50 - t301, -t50 * t192, t315 * t49 + t308, t49 * t192, 0, t15 * t192 + t315 * t61 + t48 * t96 + t49 * t75 + t52 * t86, t16 * t192 + t219 * t61 - t221 * t96 + t75 * t50 + t52 * t87, t15 * t219 - t16 * t315 - t221 * t41 - t29 * t50 - t30 * t49 - t42 * t48 - t6 * t86 + t311, -t15 * t29 + t16 * t30 + t42 * t6 + t52 * t96 + t61 * t75 + t312, t225 * t252 + (-t225 * t50 - t27 * t87) * t202 (-t287 + t289) * t50 + (t292 - t288 + (t286 + t290) * qJD(6)) * t87, t202 * t228 - t225 * t49 - t252 * t325 - t27 * t86, t63 * t251 + (t28 * t87 + t50 * t63) * t199, -t199 * t228 - t251 * t325 - t28 * t86 - t49 * t63, t325 * t49 + t308, t15 * t63 + t18 * t48 + t199 * t231 + t2 * t86 - t226 * t49 + t25 * t251 + t28 * t41 + t325 * t4, -t1 * t86 - t15 * t225 - t19 * t48 + t202 * t231 - t25 * t252 - t27 * t41 - t3 * t325 - t49 * t8, t18 * t27 - t19 * t28 - t3 * t63 + t4 * t225 + t230 * t50 + (qJD(6) * t229 - t1 * t199 - t2 * t202) * t87, t1 * t19 + t15 * t25 + t18 * t2 - t226 * t4 + t3 * t8 + t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170, -t157, 0, t170, 0, 0, t207 * pkin(1) * t201, pkin(1) * t272, 0, 0, -t170, 0, t157, 0, 0, t170 (-t129 * t201 + t146 * t204) * qJD(1) ((t158 - t194) * t201 + (-t155 + t239) * t204) * qJD(1), 0.2e1 * t193 + (t129 * t204 + t146 * t201) * qJD(1), qJ(3) * t152 + qJD(3) * t158 - t129 * t146 + (t158 * t201 + (-t155 - t297) * t204) * qJD(1) * pkin(7), -t275, t248, t317, t275, t316, 0, -t117 * t125 - t192 * t283 + t277 - t44, -t117 * t127 + t192 * t284 + t212, t153 * t82 - t154 * t83 + (-t71 - t283) * t127 + (t70 - t284) * t125, -t108 * t117 + t153 * t44 + t154 * t43 + t283 * t70 + t284 * t71, -t327, t335, t336, t327, -t329, 0, t192 * t299 - t315 * t84 + t323, -t192 * t298 - t219 * t84 + t322, t221 * t92 + t298 * t315 - t93 * t48 + t331 + (t299 - t30) * t219, -t29 * t299 - t298 * t30 - t5 * t92 + t6 * t93 - t75 * t84, t342, t343 - t345, -t344, -t290 * t325 + t288, t325 * t339 - t307 - t45, t330, -t12 * t325 + t199 * t213 - t202 * t234 + t28 * t88 + t299 * t63 - t328, t13 * t325 + t199 * t234 + t202 * t213 - t225 * t299 - t27 * t88 - t332, -t12 * t225 + t13 * t63 + (-t28 * t89 - t63 * t68 + (-t225 * t89 - t226) * qJD(6) - t321) * t202 + (-t27 * t89 - t225 * t68 + (t63 * t89 + t8) * qJD(6) - t324) * t199, t12 * t226 - t13 * t8 + t211 * t89 - t229 * t68 + t25 * t299 + t5 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t170, 0, -t195 * t207 - t206, -qJD(2) * t158 + t129 * t261 + t173, 0, 0, 0, 0, 0, 0, -t125 * t261 - t200 * t236, -t127 * t261 - t203 * t236, -t200 * t316 + t203 * t317, -t108 * t261 + t318 * t203 + (t43 + t296) * t200, 0, 0, 0, 0, 0, 0, -t192 * t268 - t261 * t315, -t192 * t267 - t219 * t261, -t139 * t48 + t218 * t221 - t219 * t268 + t267 * t315, t139 * t6 - t261 * t75 - t267 * t30 + t268 * t29 - t309, 0, 0, 0, 0, 0, 0, -t139 * t291 - t218 * t28 - t268 * t63 + (-t123 * t199 - t139 * t254 - t98) * t325, -t139 * t45 + t218 * t27 + t268 * t225 + (-t123 * t202 + t139 * t255 + t99) * t325, t63 * t99 - t225 * t98 + (-t287 - t289) * t123 + (-t292 - t288 + (-t286 + t290) * qJD(6)) * t139, -t123 * t229 + t139 * t211 + t226 * t98 - t25 * t268 - t8 * t99 - t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, -t248, -t317, -t275, -t316, 0, -t277 + t318, -t212 - t296, 0, 0, t327, -t335, -t336, -t327, t329, 0, -t192 * t31 - t310 * t315 - t323, -t192 * t32 - t219 * t310 - t322, -t331 + t32 * t315 + (-t197 * t48 + t221 * t282) * pkin(4) + (t30 - t31) * t219, t29 * t31 - t30 * t32 + (-t127 * t75 + t197 * t6 - t282 * t5) * pkin(4), -t342, t325 * t289 - t343, t344, t339 * t63 - t288, t224 + t307, -t330, t177 * t28 + t199 * t220 - t202 * t232 - t31 * t63 - t325 * t9 + t328, t10 * t325 - t177 * t27 + t199 * t232 + t202 * t220 + t225 * t31 + t332, t10 * t63 - t225 * t9 + (-t176 * t28 + (-t176 * t225 + t226) * qJD(6) + t321) * t202 + (-t176 * t27 + (t176 * t63 - t8) * qJD(6) + t324) * t199, -t10 * t8 + t176 * t211 + t177 * t5 + t226 * t9 - t25 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 - t294, -t221 + t285, -t302 - t303, t219 * t29 + t30 * t315 + t52, 0, 0, 0, 0, 0, 0, t224 - t307, t304 + t326 (t27 - t306) * t202 - t345 + t300, t333 * t199 + t334 * t202 - t25 * t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t305, t225 ^ 2 - t63 ^ 2, t325 * t63 - t27, t305, -t225 * t325 - t28, t48, t225 * t25 + t334, t25 * t63 - t333, 0, 0;];
tauc_reg  = t7;
