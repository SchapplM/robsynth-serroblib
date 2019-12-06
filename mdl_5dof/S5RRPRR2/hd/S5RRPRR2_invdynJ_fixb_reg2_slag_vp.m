% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:28:32
% EndTime: 2019-12-05 18:28:41
% DurationCPUTime: 4.68s
% Computational Cost: add. (8773->410), mult. (21556->527), div. (0->0), fcn. (16457->16), ass. (0->209)
t220 = sin(pkin(9));
t221 = cos(pkin(9));
t229 = cos(qJ(2));
t279 = t221 * t229;
t263 = qJD(1) * t279;
t225 = sin(qJ(2));
t275 = qJD(1) * t225;
t156 = t220 * t275 - t263;
t168 = t220 * t229 + t221 * t225;
t159 = t168 * qJD(1);
t224 = sin(qJ(4));
t228 = cos(qJ(4));
t108 = -t228 * t156 - t159 * t224;
t223 = sin(qJ(5));
t227 = cos(qJ(5));
t158 = t168 * qJD(2);
t268 = t229 * qJDD(1);
t269 = t225 * qJDD(1);
t249 = t220 * t269 - t221 * t268;
t119 = qJD(1) * t158 + t249;
t270 = qJD(1) * qJD(2);
t261 = t225 * t270;
t243 = qJDD(1) * t168 - t220 * t261;
t260 = t229 * t270;
t120 = t221 * t260 + t243;
t247 = t156 * t224 - t228 * t159;
t239 = qJD(4) * t247 - t228 * t119 - t224 * t120;
t273 = qJD(4) * t228;
t274 = qJD(4) * t224;
t244 = -t224 * t119 + t228 * t120 - t156 * t273 - t159 * t274;
t271 = qJD(5) * t227;
t272 = qJD(5) * t223;
t10 = -t108 * t271 - t223 * t239 - t227 * t244 - t247 * t272;
t216 = qJD(2) + qJD(4);
t209 = qJD(5) + t216;
t59 = t108 * t227 + t223 * t247;
t293 = t209 * t59;
t319 = -t10 - t293;
t55 = t108 * t223 - t227 * t247;
t240 = -qJD(5) * t55 - t223 * t244 + t227 * t239;
t290 = t55 * t209;
t314 = t240 + t290;
t299 = t59 ^ 2;
t300 = t55 ^ 2;
t320 = -t299 + t300;
t298 = t59 * t55;
t217 = qJ(2) + pkin(9);
t210 = qJ(4) + t217;
t202 = qJ(5) + t210;
t194 = sin(t202);
t195 = cos(t202);
t324 = pkin(8) * t247;
t297 = qJ(3) + pkin(6);
t186 = t297 * t229;
t174 = qJD(1) * t186;
t162 = t220 * t174;
t185 = t297 * t225;
t173 = qJD(1) * t185;
t294 = qJD(2) * pkin(2);
t166 = -t173 + t294;
t117 = t221 * t166 - t162;
t307 = pkin(7) * t159;
t80 = qJD(2) * pkin(3) + t117 - t307;
t280 = t221 * t174;
t118 = t220 * t166 + t280;
t308 = pkin(7) * t156;
t86 = t118 - t308;
t44 = -t224 * t86 + t228 * t80;
t29 = t44 + t324;
t27 = pkin(4) * t216 + t29;
t325 = pkin(8) * t108;
t45 = t224 * t80 + t228 * t86;
t30 = t45 + t325;
t291 = t227 * t30;
t13 = t223 * t27 + t291;
t214 = qJDD(2) + qJDD(4);
t256 = qJD(2) * t297;
t153 = -t225 * qJD(3) - t229 * t256;
t116 = qJDD(2) * pkin(2) + t153 * qJD(1) - qJDD(1) * t185;
t152 = t229 * qJD(3) - t225 * t256;
t123 = qJD(1) * t152 + qJDD(1) * t186;
t64 = t221 * t116 - t220 * t123;
t43 = qJDD(2) * pkin(3) - pkin(7) * t120 + t64;
t65 = t220 * t116 + t221 * t123;
t46 = -pkin(7) * t119 + t65;
t9 = -qJD(4) * t45 - t224 * t46 + t228 * t43;
t6 = t214 * pkin(4) - pkin(8) * t244 + t9;
t255 = -t224 * t43 - t228 * t46 - t80 * t273 + t86 * t274;
t7 = pkin(8) * t239 - t255;
t2 = -qJD(5) * t13 - t223 * t7 + t227 * t6;
t226 = sin(qJ(1));
t230 = cos(qJ(1));
t251 = g(1) * t230 + g(2) * t226;
t212 = t229 * pkin(2);
t203 = t212 + pkin(1);
t180 = -qJD(1) * t203 + qJD(3);
t126 = t156 * pkin(3) + t180;
t70 = -pkin(4) * t108 + t126;
t317 = -g(3) * t195 + t251 * t194 - t55 * t70 + t2;
t1 = (qJD(5) * t27 + t7) * t227 + t223 * t6 - t30 * t272;
t318 = g(3) * t194 + t251 * t195 - t59 * t70 - t1;
t292 = t223 * t30;
t12 = t227 * t27 - t292;
t330 = t12 * t59 + t13 * t55;
t201 = pkin(2) * t221 + pkin(3);
t310 = pkin(2) * t220;
t150 = t228 * t201 - t224 * t310;
t124 = t173 * t220 - t280;
t87 = t124 + t308;
t125 = -t221 * t173 - t162;
t88 = t125 - t307;
t289 = t150 * qJD(4) - t224 * t87 - t228 * t88;
t151 = t201 * t224 + t228 * t310;
t288 = -t151 * qJD(4) + t224 * t88 - t228 * t87;
t285 = t247 * t216;
t329 = t239 - t285;
t283 = t108 * t216;
t328 = t244 - t283;
t284 = t108 ^ 2;
t286 = t247 ^ 2;
t327 = -t284 + t286;
t199 = sin(t210);
t200 = cos(t210);
t326 = -g(3) * t200 + t251 * t199;
t323 = -t324 + t289;
t322 = t325 + t288;
t282 = t108 * t247;
t316 = g(3) * t199 - t108 * t126 + t251 * t200 + t255;
t315 = t126 * t247 + t326 + t9;
t313 = g(1) * t226 - g(2) * t230;
t312 = t159 ^ 2;
t311 = t239 * pkin(4);
t309 = pkin(2) * t225;
t302 = g(3) * t229;
t301 = t119 * pkin(3);
t148 = pkin(4) + t150;
t94 = t148 * t227 - t151 * t223;
t296 = qJD(5) * t94 + t322 * t223 + t323 * t227;
t95 = t148 * t223 + t151 * t227;
t295 = -qJD(5) * t95 - t323 * t223 + t322 * t227;
t127 = -t221 * t185 - t186 * t220;
t101 = -pkin(7) * t168 + t127;
t128 = -t220 * t185 + t221 * t186;
t167 = t220 * t225 - t279;
t102 = -pkin(7) * t167 + t128;
t53 = t224 * t101 + t228 * t102;
t287 = pkin(6) * qJDD(1);
t281 = t159 * t156;
t100 = t221 * t152 + t220 * t153;
t208 = cos(t217);
t278 = pkin(3) * t208 + t212;
t218 = t225 ^ 2;
t219 = t229 ^ 2;
t277 = t218 - t219;
t276 = t218 + t219;
t215 = -pkin(7) - t297;
t205 = t225 * t294;
t232 = qJD(1) ^ 2;
t266 = t225 * t232 * t229;
t265 = pkin(4) * t200 + t278;
t130 = pkin(3) * t158 + t205;
t129 = pkin(2) * t275 + pkin(3) * t159;
t52 = t228 * t101 - t102 * t224;
t99 = -t152 * t220 + t221 * t153;
t132 = pkin(3) * t167 - t203;
t253 = t225 * t260;
t207 = sin(t217);
t252 = -pkin(3) * t207 - t309;
t122 = -t167 * t224 + t168 * t228;
t33 = -pkin(8) * t122 + t52;
t121 = t228 * t167 + t168 * t224;
t34 = -pkin(8) * t121 + t53;
t20 = -t223 * t34 + t227 * t33;
t21 = t223 * t33 + t227 * t34;
t69 = -t121 * t223 + t122 * t227;
t246 = -0.2e1 * pkin(1) * t270 - pkin(6) * qJDD(2);
t161 = t167 * qJD(2);
t76 = pkin(7) * t161 + t99;
t77 = -pkin(7) * t158 + t100;
t24 = t101 * t273 - t102 * t274 + t224 * t76 + t228 * t77;
t147 = pkin(2) * t261 - qJDD(1) * t203 + qJDD(3);
t231 = qJD(2) ^ 2;
t242 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t231 + t313;
t241 = pkin(1) * t232 + t251 - t287;
t82 = t147 + t301;
t25 = -t53 * qJD(4) - t224 * t77 + t228 * t76;
t238 = t147 - t313;
t236 = t238 + t301;
t211 = -pkin(8) + t215;
t206 = qJDD(5) + t214;
t172 = pkin(1) + t278;
t154 = t156 ^ 2;
t138 = pkin(1) + t265;
t81 = pkin(4) * t121 + t132;
t71 = -pkin(4) * t247 + t129;
t68 = t227 * t121 + t122 * t223;
t67 = qJD(4) * t122 + t228 * t158 - t224 * t161;
t66 = t224 * t158 + t228 * t161 + t167 * t273 + t168 * t274;
t49 = pkin(4) * t67 + t130;
t26 = t82 - t311;
t23 = qJD(5) * t69 - t223 * t66 + t227 * t67;
t22 = t121 * t271 + t122 * t272 + t223 * t67 + t227 * t66;
t19 = t66 * pkin(8) + t25;
t18 = -pkin(8) * t67 + t24;
t15 = t227 * t29 - t292;
t14 = -t223 * t29 - t291;
t4 = -qJD(5) * t21 - t223 * t18 + t227 * t19;
t3 = qJD(5) * t20 + t227 * t18 + t223 * t19;
t5 = [0, 0, 0, 0, 0, qJDD(1), t313, t251, 0, 0, qJDD(1) * t218 + 0.2e1 * t253, 0.2e1 * t225 * t268 - 0.2e1 * t277 * t270, qJDD(2) * t225 + t229 * t231, qJDD(1) * t219 - 0.2e1 * t253, qJDD(2) * t229 - t225 * t231, 0, t225 * t246 + t229 * t242, -t225 * t242 + t229 * t246, 0.2e1 * t276 * t287 - t251, -g(1) * (-pkin(1) * t226 + pkin(6) * t230) - g(2) * (pkin(1) * t230 + pkin(6) * t226) + (t276 * pkin(6) ^ 2 + pkin(1) ^ 2) * qJDD(1), t120 * t168 - t159 * t161, -t119 * t168 - t120 * t167 + t156 * t161 - t158 * t159, -qJD(2) * t161 + qJDD(2) * t168, t119 * t167 + t156 * t158, -qJD(2) * t158 - qJDD(2) * t167, 0, t127 * qJDD(2) - t203 * t119 + t147 * t167 + t180 * t158 + t313 * t208 + (t156 * t309 + t99) * qJD(2), -t128 * qJDD(2) - t203 * t120 + t147 * t168 - t180 * t161 - t313 * t207 + (t159 * t309 - t100) * qJD(2), -t100 * t156 + t117 * t161 - t118 * t158 - t119 * t128 - t120 * t127 - t159 * t99 - t167 * t65 - t168 * t64 - t251, t65 * t128 + t118 * t100 + t64 * t127 + t117 * t99 - t147 * t203 + t180 * t205 - g(1) * (-t203 * t226 + t230 * t297) - g(2) * (t203 * t230 + t226 * t297), t122 * t244 + t247 * t66, -t108 * t66 - t121 * t244 + t122 * t239 + t247 * t67, t122 * t214 - t216 * t66, -t108 * t67 - t121 * t239, -t121 * t214 - t216 * t67, 0, -t108 * t130 + t82 * t121 + t126 * t67 - t132 * t239 + t200 * t313 + t52 * t214 + t25 * t216, t82 * t122 - t126 * t66 - t130 * t247 + t132 * t244 - t199 * t313 - t53 * t214 - t24 * t216, t108 * t24 + t121 * t255 - t122 * t9 + t239 * t53 - t244 * t52 + t247 * t25 + t44 * t66 - t45 * t67 - t251, -t255 * t53 + t45 * t24 + t9 * t52 + t44 * t25 + t82 * t132 + t126 * t130 - g(1) * (-t172 * t226 - t215 * t230) - g(2) * (t172 * t230 - t215 * t226), -t10 * t69 - t22 * t55, t10 * t68 - t22 * t59 - t23 * t55 + t240 * t69, t206 * t69 - t209 * t22, -t23 * t59 - t240 * t68, -t206 * t68 - t209 * t23, 0, t195 * t313 + t20 * t206 + t4 * t209 + t70 * t23 - t240 * t81 + t26 * t68 - t49 * t59, -t81 * t10 - t194 * t313 - t21 * t206 - t3 * t209 - t70 * t22 + t26 * t69 + t49 * t55, -t1 * t68 + t10 * t20 + t12 * t22 - t13 * t23 - t2 * t69 + t21 * t240 + t3 * t59 - t4 * t55 - t251, t1 * t21 + t13 * t3 + t2 * t20 + t12 * t4 + t26 * t81 + t70 * t49 - g(1) * (-t138 * t226 - t211 * t230) - g(2) * (t138 * t230 - t211 * t226); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t266, t277 * t232, t269, t266, t268, qJDD(2), t225 * t241 - t302, g(3) * t225 + t229 * t241, 0, 0, t281, -t154 + t312, (t156 + t263) * qJD(2) + t243, -t281, -t249, qJDD(2), -g(3) * t208 - t124 * qJD(2) - t180 * t159 + t251 * t207 + (qJDD(2) * t221 - t156 * t275) * pkin(2) + t64, g(3) * t207 + t125 * qJD(2) + t180 * t156 + t251 * t208 + (-qJDD(2) * t220 - t159 * t275) * pkin(2) - t65, (t118 + t124) * t159 + (-t117 + t125) * t156 + (-t119 * t220 - t120 * t221) * pkin(2), -t117 * t124 - t118 * t125 + (-t302 + t220 * t65 + t221 * t64 + (-qJD(1) * t180 + t251) * t225) * pkin(2), t282, t327, t328, -t282, t329, t214, t108 * t129 + t150 * t214 + t216 * t288 + t315, t129 * t247 - t151 * t214 - t216 * t289 + t316, -t150 * t244 + t151 * t239 + (t288 - t45) * t247 + (t289 + t44) * t108, -g(3) * t278 - t126 * t129 + t9 * t150 - t151 * t255 - t251 * t252 + t288 * t44 + t289 * t45, -t298, t320, t319, t298, t314, t206, t206 * t94 + t295 * t209 + t59 * t71 + t317, -t206 * t95 - t296 * t209 - t55 * t71 + t318, t94 * t10 + t240 * t95 - t295 * t55 + t296 * t59 + t330, t1 * t95 + t2 * t94 - t70 * t71 - g(3) * t265 - t251 * (-pkin(4) * t199 + t252) + t296 * t13 + t295 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t159 * qJD(2) + t249, (-t156 + t263) * qJD(2) + t243, -t154 - t312, t117 * t159 + t118 * t156 + t238, 0, 0, 0, 0, 0, 0, -t239 - t285, t244 + t283, -t284 - t286, -t45 * t108 - t247 * t44 + t236, 0, 0, 0, 0, 0, 0, -t240 + t290, -t10 + t293, -t299 - t300, t12 * t55 - t13 * t59 + t236 - t311; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t282, t327, t328, -t282, t329, t214, t45 * t216 + t315, t216 * t44 + t316, 0, 0, -t298, t320, t319, t298, t314, t206, -t14 * t209 + (t206 * t227 - t209 * t272 - t247 * t59) * pkin(4) + t317, t15 * t209 + (-t206 * t223 - t209 * t271 + t247 * t55) * pkin(4) + t318, t14 * t55 - t15 * t59 + (t10 * t227 + t240 * t223 + (t223 * t55 + t227 * t59) * qJD(5)) * pkin(4) + t330, -t12 * t14 - t13 * t15 + (t1 * t223 + t247 * t70 + t2 * t227 + (-t12 * t223 + t13 * t227) * qJD(5) + t326) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t298, t320, t319, t298, t314, t206, t13 * t209 + t317, t12 * t209 + t318, 0, 0;];
tau_reg = t5;
