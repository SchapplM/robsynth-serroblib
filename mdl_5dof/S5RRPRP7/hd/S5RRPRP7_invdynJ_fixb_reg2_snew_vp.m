% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRP7
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRP7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:32
% EndTime: 2019-12-31 20:01:45
% DurationCPUTime: 5.13s
% Computational Cost: add. (10414->376), mult. (23932->481), div. (0->0), fcn. (16471->8), ass. (0->232)
t199 = sin(pkin(8));
t200 = cos(pkin(8));
t203 = sin(qJ(2));
t206 = cos(qJ(2));
t253 = qJD(1) * t206;
t254 = qJD(1) * t203;
t173 = t199 * t254 - t200 * t253;
t170 = qJD(4) + t173;
t289 = t170 ^ 2;
t175 = t199 * t253 + t200 * t254;
t202 = sin(qJ(4));
t205 = cos(qJ(4));
t156 = -t205 * qJD(2) + t175 * t202;
t290 = t156 ^ 2;
t131 = t290 - t289;
t189 = t203 * qJDD(1);
t246 = qJD(1) * qJD(2);
t240 = t206 * t246;
t180 = t189 + t240;
t191 = t206 * qJDD(1);
t241 = t203 * t246;
t181 = t191 - t241;
t232 = t180 * t199 - t200 * t181;
t147 = qJDD(4) + t232;
t158 = qJD(2) * t202 + t175 * t205;
t262 = t158 * t156;
t96 = -t262 - t147;
t279 = t202 * t96;
t71 = -t131 * t205 - t279;
t135 = t170 * t158;
t150 = t180 * t200 + t181 * t199;
t234 = t205 * qJDD(2) - t202 * t150;
t218 = qJD(4) * t158 - t234;
t80 = -t135 + t218;
t350 = t203 * (t199 * t80 + t200 * t71) + t206 * (t199 * t71 - t200 * t80);
t155 = t158 ^ 2;
t302 = -t155 - t289;
t62 = t205 * t302 + t279;
t349 = pkin(1) * t62;
t348 = pkin(2) * t62;
t347 = pkin(3) * t62;
t346 = pkin(7) * t62;
t273 = t205 * t96;
t64 = -t202 * t302 + t273;
t345 = pkin(7) * t64;
t344 = t199 * t64;
t343 = t200 * t64;
t301 = t155 - t290;
t225 = -qJDD(2) * t202 - t205 * t150;
t216 = -t156 * qJD(4) - t225;
t263 = t156 * t170;
t298 = -t263 + t216;
t280 = t202 * t298;
t303 = t135 + t218;
t48 = t205 * t303 + t280;
t340 = t203 * (-t199 * t301 + t200 * t48) + t206 * (t199 * t48 + t200 * t301);
t299 = -t262 + t147;
t278 = t202 * t299;
t296 = -t289 - t290;
t305 = t205 * t296 - t278;
t321 = t199 * t305 - t200 * t303;
t337 = qJ(3) * t321;
t336 = -t131 * t202 + t273;
t335 = pkin(2) * t321 + pkin(7) * t305;
t272 = t205 * t299;
t306 = t202 * t296 + t272;
t320 = t199 * t303 + t200 * t305;
t334 = -pkin(2) * t306 + qJ(3) * t320;
t333 = pkin(6) * (-t203 * t321 + t206 * t320) - pkin(1) * t306;
t297 = t263 + t216;
t132 = -t155 + t289;
t322 = -t132 * t202 + t272;
t332 = t203 * (t199 * t297 + t200 * t322) + t206 * (t199 * t322 - t200 * t297);
t329 = pkin(3) * t306;
t327 = pkin(7) * t306;
t326 = qJ(5) * t298;
t323 = t205 * t132 + t278;
t300 = t155 + t290;
t319 = pkin(3) * t300;
t148 = t175 * t173;
t295 = qJDD(2) - t148;
t318 = t199 * t295;
t316 = t199 * t300;
t313 = t200 * t295;
t311 = t200 * t300;
t252 = qJD(2) * t175;
t124 = t232 + t252;
t304 = -t202 * t303 + t205 * t298;
t168 = qJD(2) * t173;
t126 = t150 - t168;
t208 = qJD(1) ^ 2;
t204 = sin(qJ(1));
t288 = cos(qJ(1));
t223 = t288 * g(1) + t204 * g(2);
t270 = qJDD(1) * pkin(6);
t214 = -t208 * pkin(1) - t223 + t270;
t160 = -t203 * g(3) + t206 * t214;
t197 = t206 ^ 2;
t194 = t197 * t208;
t220 = qJD(2) * pkin(2) - qJ(3) * t254;
t115 = -pkin(2) * t194 + t181 * qJ(3) - qJD(2) * t220 + t160;
t213 = t203 * t214;
t258 = t203 * t208;
t209 = -t213 - t180 * qJ(3) + qJDD(2) * pkin(2) + (pkin(2) * t258 + qJ(3) * t246 - g(3)) * t206;
t68 = -0.2e1 * qJD(3) * t173 + t200 * t115 + t199 * t209;
t294 = pkin(4) * t218 - t326;
t116 = pkin(4) * t156 - qJ(5) * t158;
t139 = pkin(3) * t173 - pkin(7) * t175;
t207 = qJD(2) ^ 2;
t54 = -pkin(3) * t207 + qJDD(2) * pkin(7) - t139 * t173 + t68;
t237 = g(1) * t204 - t288 * g(2);
t222 = qJDD(1) * pkin(1) + t237;
t119 = pkin(2) * t181 + (qJ(3) * t197 + pkin(6)) * t208 - t220 * t254 - qJDD(3) + t222;
t61 = t124 * pkin(3) - t126 * pkin(7) - t119;
t30 = t202 * t61 + t205 * t54;
t236 = -t147 * qJ(5) + t156 * t116 - t30;
t293 = -pkin(4) * (t302 + t289) - qJ(5) * t96 - t236;
t260 = t170 * t205;
t243 = t156 * t260;
t219 = t202 * t218 + t243;
t244 = t200 * t262;
t245 = t199 * t262;
t292 = t203 * (t200 * t219 - t245) + t206 * (t199 * t219 + t244);
t261 = t170 * t202;
t129 = t158 * t261;
t227 = t129 - t243;
t291 = t203 * (t147 * t199 + t200 * t227) + t206 * (-t200 * t147 + t199 * t227);
t171 = t173 ^ 2;
t172 = t175 ^ 2;
t287 = pkin(3) * t199;
t286 = pkin(4) * t205;
t235 = t115 * t199 - t200 * t209;
t224 = -qJDD(2) * pkin(3) - t207 * pkin(7) + t235;
t53 = (0.2e1 * qJD(3) + t139) * t175 + t224;
t283 = t202 * t53;
t281 = t202 * t297;
t248 = qJD(3) * t175;
t67 = t235 + 0.2e1 * t248;
t36 = t199 * t68 - t200 * t67;
t277 = t203 * t36;
t276 = t205 * t53;
t274 = t205 * t297;
t271 = qJ(5) * t205;
t269 = t119 * t199;
t268 = t119 * t200;
t143 = qJDD(2) + t148;
t266 = t143 * t199;
t265 = t143 * t200;
t188 = t206 * t258;
t259 = t203 * (qJDD(2) + t188);
t257 = t206 * (qJDD(2) - t188);
t251 = qJD(2) * t199;
t250 = qJD(2) * t200;
t247 = qJD(5) * t170;
t242 = -pkin(3) * t200 - pkin(2);
t239 = -qJ(5) * t202 - pkin(3);
t37 = t199 * t67 + t200 * t68;
t29 = t202 * t54 - t205 * t61;
t16 = t202 * t29 + t205 * t30;
t159 = t206 * g(3) + t213;
t233 = t203 * t159 + t206 * t160;
t161 = 0.2e1 * t247;
t231 = t161 - t236;
t23 = -pkin(4) * t289 + t231;
t24 = -t147 * pkin(4) - qJ(5) * t289 + t116 * t158 + qJDD(5) + t29;
t230 = -pkin(4) * t24 + qJ(5) * t23;
t229 = -pkin(4) * t297 - qJ(5) * t80;
t228 = t156 * t261 - t205 * t218;
t15 = t202 * t30 - t205 * t29;
t125 = -t232 + t252;
t217 = (-t156 * t202 - t158 * t205) * t170;
t78 = t205 * t216 - t129;
t212 = t203 * (t200 * t78 + t245) + t206 * (t199 * t78 - t244);
t167 = -0.2e1 * t248;
t211 = 0.2e1 * qJD(5) * t158 - t139 * t175 + t167 - t224 - t294;
t210 = pkin(4) * t299 + qJ(5) * t296 - t24;
t196 = t203 ^ 2;
t192 = t196 * t208;
t182 = t191 - 0.2e1 * t241;
t179 = t189 + 0.2e1 * t240;
t177 = pkin(6) * t208 + t222;
t164 = -t172 - t207;
t163 = -t172 + t207;
t162 = t171 - t207;
t140 = -t207 - t171;
t127 = t150 + t168;
t120 = -t171 - t172;
t109 = -t164 * t199 - t265;
t108 = t164 * t200 - t266;
t99 = t140 * t200 - t318;
t98 = t140 * t199 + t313;
t89 = t125 * t200 + t127 * t199;
t88 = t125 * t199 - t127 * t200;
t87 = (qJD(4) + t170) * t156 + t225;
t82 = (-qJD(4) + t170) * t158 + t234;
t77 = t158 * t260 + t202 * t216;
t51 = t205 * t82 + t281;
t49 = -t205 * t80 + t281;
t47 = t202 * t82 - t274;
t46 = -t202 * t80 - t274;
t44 = -t199 * t87 + t343;
t42 = t200 * t87 + t344;
t40 = -t199 * t298 - t343;
t38 = t200 * t298 - t344;
t35 = t200 * t51 - t316;
t34 = t200 * t49 - t316;
t33 = t199 * t51 + t311;
t32 = t199 * t49 + t311;
t31 = t276 - t346;
t27 = t283 - t327;
t26 = (pkin(4) * t170 - 0.2e1 * qJD(5)) * t158 + t53 + t294;
t25 = -pkin(3) * t46 - t229;
t22 = t30 - t347;
t21 = t29 - t329;
t20 = (-t303 - t135) * pkin(4) + t211;
t19 = -pkin(4) * t135 + t211 + t326;
t18 = qJ(5) * t300 + t24;
t17 = (t300 - t289) * pkin(4) + t231;
t14 = -t210 - t329;
t13 = -0.2e1 * t247 - t293 + t347;
t12 = -t20 * t202 - t271 * t303 - t327;
t11 = -pkin(4) * t280 + t19 * t205 + t346;
t9 = t16 * t199 - t200 * t53;
t8 = -pkin(7) * t47 - t15;
t7 = t202 * t24 + t205 * t23;
t6 = t202 * t23 - t205 * t24;
t5 = -pkin(7) * t46 - t17 * t202 + t18 * t205;
t4 = t199 * t26 + t200 * t7;
t3 = t199 * t7 - t200 * t26;
t2 = -pkin(7) * t6 + (pkin(4) * t202 - t271) * t26;
t1 = -pkin(3) * t6 - t230;
t10 = [0, 0, 0, 0, 0, qJDD(1), t237, t223, 0, 0, (t180 + t240) * t203, t179 * t206 + t182 * t203, t259 + t206 * (-t192 + t207), (t181 - t241) * t206, t203 * (t194 - t207) + t257, 0, t206 * t177 + pkin(1) * t182 + pkin(6) * (t206 * (-t194 - t207) - t259), -t203 * t177 - pkin(1) * t179 + pkin(6) * (-t257 - t203 * (-t192 - t207)), pkin(1) * (t192 + t194) + (t196 + t197) * t270 + t233, pkin(1) * t177 + pkin(6) * t233, t203 * (t150 * t200 - t175 * t251) + t206 * (t150 * t199 + t175 * t250), t203 * (-t124 * t200 - t126 * t199) + t206 * (-t124 * t199 + t126 * t200), t203 * (-t163 * t199 + t313) + t206 * (t163 * t200 + t318), t203 * (t173 * t250 + t199 * t232) + t206 * (t173 * t251 - t200 * t232), t203 * (t162 * t200 - t266) + t206 * (t162 * t199 + t265), (t203 * (-t173 * t200 + t175 * t199) + t206 * (-t173 * t199 - t175 * t200)) * qJD(2), t203 * (-qJ(3) * t98 - t269) + t206 * (-pkin(2) * t124 + qJ(3) * t99 + t268) - pkin(1) * t124 + pkin(6) * (-t203 * t98 + t206 * t99), t203 * (-qJ(3) * t108 - t268) + t206 * (-pkin(2) * t126 + qJ(3) * t109 - t269) - pkin(1) * t126 + pkin(6) * (-t108 * t203 + t109 * t206), t203 * (-qJ(3) * t88 - t36) + t206 * (-pkin(2) * t120 + qJ(3) * t89 + t37) - pkin(1) * t120 + pkin(6) * (-t203 * t88 + t206 * t89), -qJ(3) * t277 + t206 * (pkin(2) * t119 + qJ(3) * t37) + pkin(1) * t119 + pkin(6) * (t206 * t37 - t277), t212, -t340, t332, t292, -t350, t291, t203 * (-t199 * t21 + t200 * t27 - t337) + t206 * (t199 * t27 + t200 * t21 + t334) + t333, t203 * (-qJ(3) * t42 - t199 * t22 + t200 * t31) + t206 * (qJ(3) * t44 + t199 * t31 + t200 * t22 - t348) - t349 + pkin(6) * (-t203 * t42 + t206 * t44), t203 * (-qJ(3) * t33 + t200 * t8) + t206 * (qJ(3) * t35 + t199 * t8) + pkin(6) * (-t203 * t33 + t206 * t35) + (t203 * t287 + t206 * t242 - pkin(1)) * t47, (t203 * (-pkin(7) * t200 + t287) + t206 * (-pkin(7) * t199 + t242) - pkin(1)) * t15 + (pkin(6) + qJ(3)) * ((t16 * t200 + t199 * t53) * t206 - t203 * t9), t212, t332, t340, t291, t350, t292, t203 * (t12 * t200 - t14 * t199 - t337) + t206 * (t12 * t199 + t14 * t200 + t334) + t333, t203 * (-qJ(3) * t32 - t199 * t25 + t200 * t5) + t206 * (-pkin(2) * t46 + qJ(3) * t34 + t199 * t5 + t200 * t25) - pkin(1) * t46 + pkin(6) * (-t203 * t32 + t206 * t34), t203 * (-qJ(3) * t38 + t11 * t200 - t13 * t199) + t206 * (qJ(3) * t40 + t11 * t199 + t13 * t200 + t348) + t349 + pkin(6) * (-t203 * t38 + t206 * t40), t203 * (-qJ(3) * t3 - t1 * t199 + t2 * t200) + t206 * (-pkin(2) * t6 + qJ(3) * t4 + t1 * t200 + t199 * t2) - pkin(1) * t6 + pkin(6) * (-t203 * t3 + t206 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t188, t192 - t194, t189, t188, t191, qJDD(2), -t159, -t160, 0, 0, t148, t172 - t171, t127, -t148, t125, qJDD(2), pkin(2) * t98 + t167 - t235, pkin(2) * t108 - t68, pkin(2) * t88, pkin(2) * t36, t77, t304, t323, t228, -t336, t217, -pkin(3) * t303 - t276 + t335, pkin(2) * t42 + pkin(3) * t87 + t283 + t345, pkin(2) * t33 + pkin(7) * t51 + t16 + t319, pkin(2) * t9 - pkin(3) * t53 + pkin(7) * t16, t77, t323, -t304, t217, t336, t228, t20 * t205 + t239 * t303 + t335, pkin(2) * t32 + pkin(7) * t49 + t17 * t205 + t18 * t202 + t319, pkin(2) * t38 - t345 + t19 * t202 + (pkin(3) + t286) * t298, pkin(2) * t3 + pkin(7) * t7 + (t239 - t286) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t126, t120, -t119, 0, 0, 0, 0, 0, 0, t306, t62, t47, t15, 0, 0, 0, 0, 0, 0, t306, t46, -t62, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t262, t301, t297, -t262, -t80, t147, -t29, -t30, 0, 0, t262, t297, -t301, t147, t80, -t262, t210, t229, t161 + t293, t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299, t297, t302, t24;];
tauJ_reg = t10;
