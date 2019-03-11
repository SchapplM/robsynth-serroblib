% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% tauc_reg [6x33]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:14:29
% EndTime: 2019-03-09 13:14:41
% DurationCPUTime: 4.54s
% Computational Cost: add. (10680->349), mult. (28294->474), div. (0->0), fcn. (23059->10), ass. (0->224)
t198 = cos(qJ(6));
t254 = qJD(6) * t198;
t195 = sin(qJ(5));
t199 = cos(qJ(5));
t192 = sin(pkin(11));
t193 = cos(pkin(11));
t197 = sin(qJ(2));
t201 = cos(qJ(2));
t169 = -t192 * t197 + t193 * t201;
t159 = t169 * qJD(1);
t170 = t192 * t201 + t193 * t197;
t161 = t170 * qJD(1);
t196 = sin(qJ(4));
t200 = cos(qJ(4));
t218 = t196 * t159 + t200 * t161;
t219 = t200 * t159 - t196 * t161;
t72 = t195 * t218 - t199 * t219;
t336 = t198 * t72;
t340 = t254 + t336;
t259 = qJD(4) * t196;
t287 = -qJ(3) - pkin(7);
t238 = qJD(2) * t287;
t156 = t201 * qJD(3) + t197 * t238;
t138 = t156 * qJD(1);
t157 = -t197 * qJD(3) + t201 * t238;
t139 = t157 * qJD(1);
t100 = -t192 * t138 + t193 * t139;
t253 = qJD(1) * qJD(2);
t246 = t201 * t253;
t247 = t197 * t253;
t151 = -t192 * t247 + t193 * t246;
t83 = -t151 * pkin(8) + t100;
t101 = t193 * t138 + t192 * t139;
t160 = t170 * qJD(2);
t150 = qJD(1) * t160;
t84 = -t150 * pkin(8) + t101;
t179 = t287 * t201;
t175 = qJD(1) * t179;
t164 = t192 * t175;
t178 = t287 * t197;
t174 = qJD(1) * t178;
t281 = qJD(2) * pkin(2);
t168 = t174 + t281;
t120 = t193 * t168 + t164;
t293 = t161 * pkin(8);
t94 = qJD(2) * pkin(3) + t120 - t293;
t267 = t193 * t175;
t121 = t192 * t168 - t267;
t294 = t159 * pkin(8);
t99 = t121 + t294;
t210 = -(qJD(4) * t94 + t84) * t200 - t196 * t83 + t99 * t259;
t69 = qJD(4) * t218 + t200 * t150 + t196 * t151;
t17 = -t69 * pkin(9) - t210;
t224 = -t196 * t94 - t200 * t99;
t208 = qJD(4) * t224 - t196 * t84 + t200 * t83;
t258 = qJD(4) * t200;
t68 = -t196 * t150 + t200 * t151 + t159 * t258 - t161 * t259;
t18 = -t68 * pkin(9) + t208;
t257 = qJD(5) * t195;
t313 = pkin(9) * t219;
t49 = -t224 + t313;
t242 = t195 * t18 - t49 * t257;
t189 = qJD(2) + qJD(4);
t239 = -t196 * t99 + t200 * t94;
t314 = pkin(9) * t218;
t48 = t239 - t314;
t46 = t189 * pkin(4) + t48;
t2 = (qJD(5) * t46 + t17) * t199 + t242;
t249 = -t201 * pkin(2) - pkin(1);
t226 = t249 * qJD(1);
t176 = qJD(3) + t226;
t126 = -t159 * pkin(3) + t176;
t85 = -pkin(4) * t219 + t126;
t326 = t85 * t72;
t339 = t326 - t2;
t188 = qJD(5) + t189;
t271 = t72 * t188;
t256 = qJD(5) * t199;
t34 = -t195 * t69 + t199 * t68 - t218 * t257 + t219 * t256;
t318 = t34 + t271;
t194 = sin(qJ(6));
t255 = qJD(6) * t194;
t299 = t195 * t219 + t199 * t218;
t14 = t188 * t254 + t198 * t34 - t255 * t299;
t60 = t194 * t188 + t198 * t299;
t15 = qJD(6) * t60 + t194 * t34;
t58 = -t198 * t188 + t194 * t299;
t338 = t14 * t198 - t194 * t15 - t340 * t58;
t12 = t14 * t194;
t331 = t340 * t60 + t12;
t35 = qJD(5) * t299 + t195 * t68 + t199 * t69;
t31 = t194 * t35;
t335 = -qJD(6) - t72;
t63 = t335 * t254;
t284 = t31 - t63;
t290 = t60 * t299;
t330 = -t335 * t336 + t284 - t290;
t277 = t195 * t49;
t23 = t199 * t46 - t277;
t21 = -t188 * pkin(5) - t23;
t327 = t21 * t72;
t273 = t199 * t49;
t24 = t195 * t46 + t273;
t243 = t195 * t17 - t199 * t18;
t3 = qJD(5) * t24 + t243;
t311 = t85 * t299;
t337 = -t3 - t311;
t324 = t299 * t72;
t272 = t299 * t188;
t300 = -t35 + t272;
t124 = -t192 * t174 + t267;
t102 = t124 - t294;
t125 = t193 * t174 + t164;
t103 = t125 - t293;
t183 = t193 * pkin(2) + pkin(3);
t296 = pkin(2) * t192;
t215 = t200 * t183 - t196 * t296;
t334 = -t215 * qJD(4) + t196 * t102 + t200 * t103;
t155 = t196 * t183 + t200 * t296;
t333 = -t155 * qJD(4) - t200 * t102 + t196 * t103;
t332 = qJD(6) + t335;
t319 = t299 ^ 2 - t72 ^ 2;
t45 = pkin(5) * t299 + t72 * pkin(10);
t289 = t299 * t58;
t321 = t194 * t335;
t33 = t198 * t35;
t329 = -t321 * t335 + t289 + t33;
t328 = t321 * t60 + t338;
t323 = t314 - t334;
t322 = -t313 - t333;
t312 = t335 * t299;
t269 = t218 * t189;
t320 = -t69 + t269;
t22 = t188 * pkin(10) + t24;
t38 = t72 * pkin(5) - pkin(10) * t299 + t85;
t9 = t194 * t38 + t198 * t22;
t305 = t3 * t194 + t21 * t254 + t299 * t9;
t225 = t194 * t22 - t198 * t38;
t306 = t21 * t255 + t225 * t299;
t316 = -0.2e1 * t253;
t295 = t218 * pkin(4);
t310 = -t255 * t335 - t33;
t268 = t219 * t189;
t309 = t68 - t268;
t308 = t218 * t219;
t307 = t218 ^ 2 - t219 ^ 2;
t302 = -t126 * t218 + t208;
t301 = -t126 * t219 + t210;
t128 = t193 * t178 + t192 * t179;
t113 = -t170 * pkin(8) + t128;
t129 = t192 * t178 - t193 * t179;
t114 = t169 * pkin(8) + t129;
t123 = t196 * t169 + t200 * t170;
t52 = -t123 * pkin(9) + t200 * t113 - t196 * t114;
t217 = t200 * t169 - t196 * t170;
t223 = -t196 * t113 - t200 * t114;
t53 = pkin(9) * t217 - t223;
t37 = t195 * t52 + t199 * t53;
t111 = -t192 * t156 + t193 * t157;
t163 = t169 * qJD(2);
t90 = -t163 * pkin(8) + t111;
t112 = t193 * t156 + t192 * t157;
t91 = -t160 * pkin(8) + t112;
t214 = t113 * t258 - t114 * t259 + t196 * t90 + t200 * t91;
t79 = qJD(4) * t123 + t200 * t160 + t196 * t163;
t29 = -t79 * pkin(9) + t214;
t207 = qJD(4) * t223 - t196 * t91 + t200 * t90;
t78 = qJD(4) * t217 - t196 * t160 + t200 * t163;
t30 = -t78 * pkin(9) + t207;
t36 = t195 * t53 - t199 * t52;
t4 = -qJD(5) * t36 + t195 * t30 + t199 * t29;
t80 = t195 * t123 - t199 * t217;
t40 = -qJD(5) * t80 - t195 * t79 + t199 * t78;
t81 = t199 * t123 + t195 * t217;
t141 = -t169 * pkin(3) + t249;
t95 = -pkin(4) * t217 + t141;
t43 = t80 * pkin(5) - t81 * pkin(10) + t95;
t298 = t21 * t40 + (qJD(6) * t43 + t4) * t335 - (qJD(6) * t38 + t2) * t80 + t3 * t81 - t37 * t35;
t292 = t21 * t81;
t291 = t43 * t35;
t288 = t81 * t35;
t154 = pkin(4) + t215;
t221 = t199 * t154 - t195 * t155;
t285 = -qJD(5) * t221 + t322 * t195 - t323 * t199;
t220 = t195 * t154 + t199 * t155;
t283 = qJD(5) * t220 + t323 * t195 + t322 * t199;
t279 = t194 * t60;
t203 = qJD(1) ^ 2;
t266 = t201 * t203;
t202 = qJD(2) ^ 2;
t265 = t202 * t197;
t264 = t202 * t201;
t261 = t197 ^ 2 - t201 ^ 2;
t260 = qJD(1) * t197;
t187 = t197 * t281;
t250 = t81 * t255;
t248 = -pkin(4) * t188 - t46;
t182 = pkin(2) * t247;
t127 = t150 * pkin(3) + t182;
t134 = t160 * pkin(3) + t187;
t133 = pkin(2) * t260 + t161 * pkin(3);
t232 = pkin(1) * t316;
t107 = pkin(10) + t220;
t86 = t133 + t295;
t231 = qJD(6) * t107 + t45 + t86;
t184 = t195 * pkin(4) + pkin(10);
t230 = qJD(6) * t184 + t295 + t45;
t25 = t195 * t48 + t273;
t228 = pkin(4) * t257 - t25;
t227 = -t335 * t40 + t288;
t216 = t321 * t72 - t310;
t54 = t69 * pkin(4) + t127;
t57 = t79 * pkin(4) + t134;
t212 = -t107 * t35 - t285 * t335 + t327;
t26 = t199 * t48 - t277;
t204 = -t184 * t35 + t327 - (-pkin(4) * t256 + t26) * t335;
t185 = -t199 * pkin(4) - pkin(5);
t106 = -pkin(5) - t221;
t41 = qJD(5) * t81 + t195 * t78 + t199 * t79;
t10 = t41 * pkin(5) - t40 * pkin(10) + t57;
t7 = t35 * pkin(5) - t34 * pkin(10) + t54;
t6 = t198 * t7;
t5 = qJD(5) * t37 + t195 * t29 - t199 * t30;
t1 = [0, 0, 0, 0.2e1 * t197 * t246, t261 * t316, t264, -t265, 0, -pkin(7) * t264 + t197 * t232, pkin(7) * t265 + t201 * t232, -t100 * t170 + t101 * t169 - t111 * t161 + t112 * t159 - t120 * t163 - t121 * t160 - t128 * t151 - t129 * t150, t100 * t128 + t101 * t129 + t120 * t111 + t121 * t112 + (t176 + t226) * t187, t68 * t123 + t218 * t78, -t123 * t69 + t217 * t68 - t218 * t79 + t219 * t78, t78 * t189, -t79 * t189, 0, t126 * t79 - t127 * t217 - t134 * t219 + t141 * t69 + t189 * t207, t127 * t123 + t126 * t78 + t134 * t218 + t141 * t68 - t189 * t214, t299 * t40 + t34 * t81, -t299 * t41 - t34 * t80 - t40 * t72 - t288, t40 * t188, -t41 * t188, 0, -t5 * t188 + t95 * t35 + t85 * t41 + t54 * t80 + t57 * t72, -t4 * t188 + t299 * t57 + t95 * t34 + t85 * t40 + t54 * t81, -t60 * t250 + (t14 * t81 + t40 * t60) * t198 (-t198 * t58 - t279) * t40 + (-t12 - t15 * t198 + (t194 * t58 - t198 * t60) * qJD(6)) * t81, t14 * t80 + t198 * t227 + t250 * t335 + t60 * t41, -t15 * t80 - t194 * t227 - t58 * t41 + t63 * t81, -t335 * t41 + t35 * t80, t36 * t15 - t225 * t41 + t5 * t58 + t6 * t80 + (-t10 * t335 + t291 + (-t22 * t80 + t335 * t37 + t292) * qJD(6)) * t198 + t298 * t194, t36 * t14 - t9 * t41 + t5 * t60 + ((-qJD(6) * t37 + t10) * t335 - t291 - (-qJD(6) * t22 + t7) * t80 - qJD(6) * t292) * t194 + t298 * t198; 0, 0, 0, -t197 * t266, t261 * t203, 0, 0, 0, t203 * pkin(1) * t197, pkin(1) * t266 (t121 + t124) * t161 + (t120 - t125) * t159 + (-t150 * t192 - t151 * t193) * pkin(2), -t120 * t124 - t121 * t125 + (t100 * t193 + t101 * t192 - t176 * t260) * pkin(2), -t308, t307, t309, t320, 0, t133 * t219 + t333 * t189 + t302, -t133 * t218 + t334 * t189 + t301, t324, t319, t318, t300, 0, -t188 * t283 - t86 * t72 + t337, t188 * t285 - t299 * t86 + t339, t331, t279 * t335 + t338, t330, t216 + t289, t312, t106 * t15 + t283 * t58 + (t231 * t335 - t3) * t198 + t212 * t194 + t306, t106 * t14 + t198 * t212 - t231 * t321 + t283 * t60 + t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t159 ^ 2 - t161 ^ 2, t120 * t161 - t121 * t159 + t182, 0, 0, 0, 0, 0, t69 + t269, t68 + t268, 0, 0, 0, 0, 0, t35 + t272, t34 - t271, 0, 0, 0, 0, 0, t216 - t289, -t198 * t335 ^ 2 - t290 - t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t308, t307, t309, t320, 0, -t189 * t224 + t302, t189 * t239 + t301, t324, t319, t318, t300, 0, -t72 * t295 + t25 * t188 - t311 + (t195 * t248 - t273) * qJD(5) - t243, -t299 * t295 + t26 * t188 + t326 + (qJD(5) * t248 - t17) * t199 - t242, t331, t328, t330, t329, t312, t185 * t15 + t228 * t58 + (t230 * t335 - t3) * t198 + t204 * t194 + t306, t185 * t14 + t198 * t204 + t228 * t60 - t230 * t321 + t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t324, t319, t318, t300, 0, t24 * t188 + t337, t23 * t188 + t339, t331, t328, t330, t329, t312, -pkin(5) * t15 - t3 * t198 + (-t194 * t23 + t198 * t45) * t335 - t24 * t58 + t194 * t327 - t284 * pkin(10) + t306, -pkin(5) * t14 - (t194 * t45 + t198 * t23) * t335 - t24 * t60 + t21 * t336 + t310 * pkin(10) + t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t58, -t58 ^ 2 + t60 ^ 2, -t335 * t58 + t14, -t335 * t60 - t15, t35, -t194 * t2 - t21 * t60 - t332 * t9 + t6, -t194 * t7 - t198 * t2 + t21 * t58 + t332 * t225;];
tauc_reg  = t1;
