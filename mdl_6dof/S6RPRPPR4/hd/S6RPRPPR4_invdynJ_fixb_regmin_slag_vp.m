% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:48:26
% EndTime: 2019-03-09 02:48:35
% DurationCPUTime: 4.01s
% Computational Cost: add. (4837->455), mult. (11599->564), div. (0->0), fcn. (9036->12), ass. (0->226)
t187 = sin(pkin(9));
t189 = cos(pkin(9));
t192 = sin(qJ(3));
t310 = cos(qJ(3));
t146 = t187 * t310 + t189 * t192;
t186 = sin(pkin(10));
t188 = cos(pkin(10));
t191 = sin(qJ(6));
t194 = cos(qJ(6));
t324 = -t186 * t194 + t188 * t191;
t80 = t324 * t146;
t185 = pkin(9) + qJ(3);
t181 = cos(t185);
t172 = g(3) * t181;
t270 = qJD(1) * qJD(2);
t304 = pkin(7) + qJ(2);
t314 = qJDD(1) * t304 + t270;
t112 = t314 * t187;
t113 = t314 * t189;
t155 = t304 * t187;
t147 = qJD(1) * t155;
t157 = t304 * t189;
t148 = qJD(1) * t157;
t95 = -t147 * t192 + t148 * t310;
t336 = qJD(3) * t95;
t250 = -t310 * t112 - t192 * t113 - t336;
t208 = qJDD(3) * pkin(3) - qJDD(4) + t250;
t259 = t208 - t172;
t180 = sin(t185);
t193 = sin(qJ(1));
t195 = cos(qJ(1));
t243 = g(1) * t195 + g(2) * t193;
t329 = t243 * t180;
t198 = -t329 - t259;
t262 = t310 * t189;
t168 = qJD(1) * t262;
t282 = t192 * t187;
t261 = qJD(1) * t282;
t132 = -t168 + t261;
t271 = qJD(6) - t132;
t323 = t146 * qJD(1);
t103 = -qJD(3) * t188 + t186 * t323;
t224 = qJD(3) * t186 + t188 * t323;
t54 = -t103 * t194 + t191 * t224;
t339 = t271 * t54;
t56 = t103 * t191 + t194 * t224;
t338 = t271 * t56;
t273 = qJD(6) * t194;
t274 = qJD(6) * t191;
t295 = t132 * t324 + t186 * t273 - t188 * t274;
t102 = -t155 * t192 + t157 * t310;
t70 = qJD(2) * t146 + qJD(3) * t102;
t337 = t188 * t146 * qJD(5) - t70;
t335 = t103 * t132;
t292 = qJDD(1) * pkin(1);
t326 = g(1) * t193 - g(2) * t195;
t222 = -qJDD(2) + t292 + t326;
t334 = qJD(3) * t323;
t333 = t224 ^ 2;
t176 = pkin(2) * t189 + pkin(1);
t152 = -qJD(1) * t176 + qJD(2);
t68 = pkin(3) * t132 - qJ(4) * t323 + t152;
t86 = qJD(3) * qJ(4) + t95;
t40 = t186 * t68 + t188 * t86;
t26 = qJ(5) * t132 + t40;
t151 = -qJDD(1) * t176 + qJDD(2);
t254 = qJDD(1) * t310;
t268 = t189 * qJDD(1);
t265 = qJD(3) * t168 + t187 * t254 + t192 * t268;
t91 = -qJD(3) * t261 + t265;
t138 = t146 * qJD(3);
t269 = t187 * qJDD(1);
t236 = -t189 * t254 + t192 * t269;
t92 = qJD(1) * t138 + t236;
t32 = pkin(3) * t92 - qJ(4) * t91 - qJD(4) * t323 + t151;
t215 = t147 * t310 + t148 * t192;
t216 = -t192 * t112 + t113 * t310;
t41 = qJDD(3) * qJ(4) + (qJD(4) - t215) * qJD(3) + t216;
t10 = -t186 * t41 + t188 * t32;
t239 = qJDD(5) - t10;
t8 = -pkin(4) * t92 + t239;
t332 = -t132 * t26 + t8;
t328 = -t155 * t310 - t192 * t157;
t327 = pkin(3) * t181 + qJ(4) * t180;
t275 = qJD(5) * t224;
t74 = qJDD(3) * t186 + t188 * t91;
t325 = qJ(5) * t74 + t275;
t322 = -qJD(6) + t271;
t321 = qJ(2) * qJDD(1);
t73 = -qJDD(3) * t188 + t186 * t91;
t15 = qJD(6) * t56 + t191 * t74 - t194 * t73;
t143 = t186 * t191 + t188 * t194;
t296 = t271 * t143;
t89 = -qJDD(6) + t92;
t319 = t271 * t296 - t324 * t89;
t127 = t132 ^ 2;
t318 = -t103 * t323 - t127 * t186 + t188 * t92;
t297 = t186 * t92;
t317 = t127 * t188 + t224 * t323 + t297;
t305 = g(3) * t180;
t316 = t181 * t243 + t305;
t293 = qJ(5) * t188;
t312 = pkin(4) + pkin(5);
t315 = t186 * t312 - t293;
t199 = t208 + t325;
t7 = -t312 * t73 + t199;
t313 = t7 + t329;
t311 = pkin(4) * t73;
t309 = pkin(8) * t186;
t303 = -pkin(8) + qJ(4);
t11 = t186 * t32 + t188 * t41;
t214 = t262 - t282;
t137 = t214 * qJD(3);
t58 = pkin(3) * t138 - qJ(4) * t137 - qJD(4) * t146;
t69 = qJD(2) * t214 + qJD(3) * t328;
t30 = t186 * t58 + t188 * t69;
t88 = pkin(3) * t323 + qJ(4) * t132;
t49 = t186 * t88 - t188 * t215;
t90 = -pkin(3) * t214 - qJ(4) * t146 - t176;
t52 = t102 * t188 + t186 * t90;
t300 = t323 * t54;
t299 = t323 * t56;
t294 = qJ(4) * t188;
t291 = t180 * t193;
t290 = t180 * t195;
t289 = t181 * t195;
t288 = t186 * t193;
t285 = t188 * t195;
t284 = t304 * t195;
t281 = t193 * t188;
t280 = t195 * t186;
t221 = qJD(3) * pkin(3) - qJD(4) - t215;
t279 = -qJD(4) - t221;
t205 = qJ(5) * t224 + t221;
t38 = pkin(4) * t103 - t205;
t278 = qJD(4) - t38;
t277 = (g(1) * t285 + g(2) * t281) * t180;
t276 = t187 ^ 2 + t189 ^ 2;
t272 = t188 * qJD(4);
t33 = qJ(5) * t323 + t49;
t44 = -qJ(5) * t214 + t52;
t267 = t92 * t294;
t266 = pkin(3) * t289 + qJ(4) * t290 + t176 * t195;
t9 = -t199 + t311;
t264 = -t9 - t172;
t2 = -pkin(8) * t74 - t312 * t92 + t239;
t6 = qJ(5) * t92 + qJD(5) * t132 + t11;
t3 = pkin(8) * t73 + t6;
t263 = -t191 * t3 + t194 * t2;
t258 = g(3) * t327;
t257 = qJ(5) * t186 + pkin(3);
t61 = t186 * t69;
t29 = t188 * t58 - t61;
t39 = -t186 * t86 + t188 * t68;
t84 = t186 * t215;
t48 = t188 * t88 + t84;
t96 = t186 * t102;
t51 = t188 * t90 - t96;
t253 = t276 * qJD(1) ^ 2;
t252 = qJD(5) * t186 - t132 * t315 + t95;
t251 = t271 ^ 2;
t17 = qJ(5) * t138 - qJD(5) * t214 + t30;
t249 = 0.2e1 * t276;
t248 = t143 * t89 - t271 * t295;
t247 = g(1) * t291 - g(2) * t290;
t246 = qJD(5) - t39;
t116 = t181 * t288 + t285;
t118 = t181 * t280 - t281;
t245 = -g(1) * t116 + g(2) * t118;
t117 = t181 * t281 - t280;
t119 = t181 * t285 + t288;
t244 = g(1) * t117 - g(2) * t119;
t241 = -t186 * t6 + t188 * t8;
t240 = t191 * t2 + t194 * t3;
t238 = t137 * t38 + t146 * t9;
t237 = pkin(4) * t186 - t293;
t235 = -t10 * t188 - t11 * t186;
t234 = -t137 * t221 - t146 * t208;
t16 = -pkin(8) * t224 - t132 * t312 + t246;
t18 = pkin(8) * t103 + t26;
t4 = t16 * t194 - t18 * t191;
t5 = t16 * t191 + t18 * t194;
t233 = -t186 * t39 + t188 * t40;
t22 = t96 + (-pkin(8) * t146 - t90) * t188 + t312 * t214;
t31 = t146 * t309 + t44;
t232 = -t191 * t31 + t194 * t22;
t231 = t191 * t22 + t194 * t31;
t228 = qJ(4) * t74 + qJD(4) * t224;
t227 = -t103 ^ 2 - t333;
t226 = t116 * t194 - t117 * t191;
t225 = t116 * t191 + t117 * t194;
t223 = pkin(4) * t188 + t257;
t220 = t132 * t224 + t73;
t154 = t303 * t186;
t213 = qJD(6) * t154 + t132 * t309 + t272 - t33;
t156 = t303 * t188;
t212 = qJD(4) * t186 - qJD(6) * t156 + t84 - (pkin(8) * t132 - t88) * t188 + t312 * t323;
t211 = -t103 * t273 - t191 * t73 - t194 * t74 + t224 * t274;
t210 = g(3) * t324;
t81 = t143 * t146;
t209 = t74 - t335;
t206 = t222 + t292;
t204 = -t188 * t74 - t186 * t73 + (-t103 * t188 + t186 * t224) * t132;
t202 = -t103 * t272 - t294 * t73 - t316;
t201 = (-g(1) * (-t176 - t327) - g(2) * t304) * t193;
t200 = t249 * t270 - t243;
t160 = qJ(4) * t289;
t158 = t193 * t181 * qJ(4);
t139 = t188 * t312 + t257;
t76 = t118 * t191 + t119 * t194;
t75 = t118 * t194 - t119 * t191;
t53 = t146 * t237 - t328;
t50 = -t132 * t237 + t95;
t47 = -t146 * t315 + t328;
t46 = pkin(4) * t214 - t51;
t43 = qJD(6) * t81 + t137 * t324;
t42 = -qJD(6) * t80 + t137 * t143;
t35 = -pkin(4) * t323 - t48;
t28 = t137 * t237 - t337;
t25 = -pkin(4) * t132 + t246;
t23 = -t103 * t312 + t205;
t21 = -pkin(4) * t138 - t29;
t20 = -t137 * t315 + t337;
t13 = t137 * t309 + t17;
t12 = t61 + (-pkin(8) * t137 - t58) * t188 - t312 * t138;
t1 = [qJDD(1), t326, t243, t206 * t189, -t206 * t187, t249 * t321 + t200, pkin(1) * t222 + (t276 * t321 + t200) * qJ(2), t137 * t323 + t146 * t91, -t132 * t137 - t138 * t323 - t146 * t92 + t214 * t91, qJD(3) * t137 + qJDD(3) * t146, -qJD(3) * t138 + qJDD(3) * t214, 0, -qJD(3) * t70 + qJDD(3) * t328 + t138 * t152 - t151 * t214 - t176 * t92 + t181 * t326, -qJD(3) * t69 - qJDD(3) * t102 + t137 * t152 + t146 * t151 - t176 * t91 - t247, -t10 * t214 + t103 * t70 + t132 * t29 + t138 * t39 + t186 * t234 - t328 * t73 + t51 * t92 + t244, t11 * t214 - t132 * t30 - t138 * t40 + t188 * t234 + t224 * t70 - t328 * t74 - t52 * t92 + t245, -t103 * t30 - t224 * t29 - t51 * t74 - t52 * t73 + t235 * t146 + (-t186 * t40 - t188 * t39) * t137 + t247, -g(1) * t284 - g(2) * t266 + t10 * t51 + t11 * t52 + t208 * t328 - t221 * t70 + t39 * t29 + t40 * t30 + t201, t103 * t28 - t132 * t21 - t138 * t25 + t186 * t238 + t214 * t8 - t46 * t92 + t53 * t73 + t244, -t103 * t17 + t224 * t21 - t44 * t73 + t46 * t74 + t241 * t146 + (-t186 * t26 + t188 * t25) * t137 + t247, t132 * t17 + t138 * t26 - t188 * t238 - t214 * t6 - t224 * t28 + t44 * t92 - t53 * t74 - t245, t6 * t44 + t26 * t17 + t9 * t53 + t38 * t28 + t8 * t46 + t25 * t21 - g(1) * (-pkin(4) * t117 - qJ(5) * t116 + t284) - g(2) * (pkin(4) * t119 + qJ(5) * t118 + t266) + t201, -t211 * t81 + t42 * t56, -t15 * t81 + t211 * t80 - t42 * t54 - t43 * t56, -t138 * t56 - t211 * t214 + t271 * t42 - t81 * t89, t138 * t54 - t15 * t214 - t271 * t43 + t80 * t89, -t138 * t271 - t214 * t89 (t12 * t194 - t13 * t191) * t271 - t232 * t89 + t263 * t214 - t4 * t138 + t20 * t54 + t47 * t15 + t7 * t80 + t23 * t43 + g(1) * t225 - g(2) * t76 + (-t214 * t5 - t231 * t271) * qJD(6) -(t12 * t191 + t13 * t194) * t271 + t231 * t89 - t240 * t214 + t5 * t138 + t20 * t56 - t47 * t211 + t7 * t81 + t23 * t42 + g(1) * t226 - g(2) * t75 + (-t214 * t4 - t232 * t271) * qJD(6); 0, 0, 0, -t268, t269, -t253, -qJ(2) * t253 - t222, 0, 0, 0, 0, 0, t236 + 0.2e1 * t334 (-t132 - t261) * qJD(3) + t265, t318, -t317, t204, t132 * t233 + t221 * t323 - t235 - t326, t318, t204, t317, -t323 * t38 + (t186 * t25 + t188 * t26) * t132 - t241 - t326, 0, 0, 0, 0, 0, t248 + t300, t299 + t319; 0, 0, 0, 0, 0, 0, 0, t323 * t132, t323 ^ 2 - t127 (t132 - t261) * qJD(3) + t265, -t236, qJDD(3), -t152 * t323 - t172 + t250 + t329 + t336, t152 * t132 - t216 + t316, -qJ(4) * t297 - pkin(3) * t73 - t103 * t95 - t323 * t39 + t259 * t188 + (t186 * t279 - t48) * t132 + t277, -t267 - pkin(3) * t74 - t224 * t95 + t323 * t40 + (t188 * t279 + t49) * t132 + t198 * t186, t103 * t49 + t224 * t48 + (-t132 * t39 + t11) * t188 + (-t132 * t40 - t10 + t228) * t186 + t202, t208 * pkin(3) - t40 * t49 - t39 * t48 + t221 * t95 - g(1) * (-pkin(3) * t290 + t160) - g(2) * (-pkin(3) * t291 + t158) - t258 + t233 * qJD(4) + (-t10 * t186 + t11 * t188) * qJ(4), -t103 * t50 + t132 * t35 + t323 * t25 - t223 * t73 + t264 * t188 + (-qJ(4) * t92 - qJD(5) * t103 - t132 * t278) * t186 + t277, t103 * t33 - t224 * t35 + (t132 * t25 + t6) * t188 + (t228 + t332) * t186 + t202, t267 + t224 * t50 - t323 * t26 + t223 * t74 + (t188 * t278 - t33) * t132 + (t329 + t264 + t275) * t186, -t26 * t33 - t38 * t50 - t25 * t35 - g(1) * t160 - g(2) * t158 - t258 + (-pkin(4) * t172 + qJ(4) * t6 + qJD(4) * t26) * t188 + (qJ(4) * t8 - qJ(5) * t172 + qJD(4) * t25 - qJD(5) * t38) * t186 + (-t9 + t329) * t223, t211 * t324 - t296 * t56, t143 * t211 + t15 * t324 - t295 * t56 + t296 * t54, t299 - t319, t248 - t300, t271 * t323 -(t154 * t194 - t156 * t191) * t89 + t139 * t15 + t4 * t323 + t252 * t54 + t295 * t23 - (t191 * t213 - t194 * t212) * t271 + (-t172 + t313) * t143 (t154 * t191 + t156 * t194) * t89 - t139 * t211 - t5 * t323 + t252 * t56 - t296 * t23 + t181 * t210 - (t191 * t212 + t194 * t213) * t271 - t313 * t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t220, t209, t227, t103 * t40 + t224 * t39 + t198, t220, t227, -t209, t103 * t26 - t224 * t25 + t198 + t311 - t325, 0, 0, 0, 0, 0, -t15 - t338, t211 + t339; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103 * t224 - t236 - t334, t74 + t335, -t127 - t333, -g(1) * t118 - g(2) * t116 - t186 * t305 + t224 * t38 + t332, 0, 0, 0, 0, 0, -t191 * t251 - t194 * t89 - t224 * t54, t191 * t89 - t194 * t251 - t224 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t54, -t54 ^ 2 + t56 ^ 2, -t211 + t339, -t15 + t338, -t89, -g(1) * t75 - g(2) * t226 + t180 * t210 - t23 * t56 + t322 * t5 + t263, g(1) * t76 + g(2) * t225 + t143 * t305 + t23 * t54 + t322 * t4 - t240;];
tau_reg  = t1;
