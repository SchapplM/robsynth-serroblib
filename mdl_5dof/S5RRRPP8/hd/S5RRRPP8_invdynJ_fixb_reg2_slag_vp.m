% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:30
% EndTime: 2019-12-31 21:09:40
% DurationCPUTime: 4.94s
% Computational Cost: add. (3402->510), mult. (7610->589), div. (0->0), fcn. (4802->6), ass. (0->243)
t154 = sin(qJ(3));
t157 = cos(qJ(3));
t158 = cos(qJ(2));
t269 = qJD(1) * t158;
t120 = -qJD(3) + t269;
t155 = sin(qJ(2));
t270 = qJD(1) * t155;
t247 = t157 * t270;
t268 = qJD(2) * t154;
t92 = t247 + t268;
t298 = t120 * t92;
t248 = t154 * t270;
t262 = t157 * qJD(2);
t90 = t248 - t262;
t299 = t120 * t90;
t261 = qJD(1) * qJD(2);
t242 = t158 * t261;
t260 = t155 * qJDD(1);
t264 = qJD(3) * t155;
t336 = qJD(1) * t264 - qJDD(2);
t38 = -qJD(3) * t262 + (-t242 - t260) * t157 + t336 * t154;
t39 = ((qJD(3) + t269) * qJD(2) + t260) * t154 + t336 * t157;
t341 = (t38 - t299) * t157 + (t39 - t298) * t154;
t309 = pkin(3) + qJ(5);
t142 = t158 * qJDD(1);
t330 = -t155 * t261 + t142;
t88 = qJDD(3) - t330;
t249 = t309 * t88;
t144 = t155 * pkin(7);
t148 = t158 * pkin(2);
t252 = -pkin(1) - t148;
t203 = -t252 + t144;
t187 = t203 * qJD(1);
t340 = -t330 * pkin(6) - qJDD(2) * pkin(7) + qJD(3) * t187;
t23 = t39 + t298;
t159 = cos(qJ(1));
t287 = t155 * t159;
t156 = sin(qJ(1));
t289 = t155 * t156;
t337 = g(1) * t287 + g(2) * t289;
t137 = pkin(6) * t269;
t105 = qJD(2) * pkin(7) + t137;
t45 = t154 * t105 + t157 * t187;
t201 = t92 * pkin(4) + t45;
t280 = qJD(4) + t201;
t266 = qJD(2) * t158;
t295 = t39 * t157;
t296 = t38 * t154;
t163 = (qJD(3) * (t154 * t90 - t157 * t92) - t295 + t296) * t155 - (t154 * t92 + t157 * t90) * t266;
t265 = qJD(3) * t154;
t335 = pkin(6) * (t155 * t262 + t158 * t265);
t86 = t90 ^ 2;
t87 = t92 ^ 2;
t256 = -t87 + t86;
t108 = qJD(4) * t120;
t80 = t88 * qJ(4);
t334 = t108 - t80;
t116 = t120 ^ 2;
t333 = -t87 - t116;
t290 = t154 * t158;
t259 = pkin(3) * t290;
t331 = pkin(3) * t265 - qJD(1) * t259 - t154 * qJD(4) - t137;
t221 = g(1) * t159 + g(2) * t156;
t329 = -pkin(4) * t39 + qJDD(5);
t192 = t120 * t265 + t157 * t88;
t19 = qJD(1) * (t120 * t290 - t155 * t90) - t192;
t263 = qJD(3) * t157;
t193 = t120 * t263 - t154 * t88;
t285 = t157 * t158;
t20 = qJD(1) * (t120 * t285 - t155 * t92) - t193;
t6 = (qJD(2) * t92 + t192) * t155 - (t120 * t262 - t38) * t158;
t326 = (-qJD(2) * t90 + t193) * t155 + (t120 * t268 + t39) * t158;
t325 = -0.2e1 * pkin(1);
t324 = 0.2e1 * t80;
t323 = pkin(4) + pkin(7);
t321 = pkin(7) * t88;
t320 = t38 * pkin(4);
t319 = t88 * pkin(3);
t318 = t90 * pkin(4);
t317 = pkin(6) * t154;
t316 = g(1) * t156;
t313 = g(2) * t159;
t312 = g(3) * t155;
t311 = g(3) * t158;
t310 = t90 * t92;
t292 = qJ(4) * t157;
t206 = qJ(5) * t154 - t292;
t189 = t206 * t158;
t308 = -qJD(1) * t189 + qJD(3) * t206 - t157 * qJD(5) + t331;
t307 = -qJ(4) * t263 + t269 * t292 + t331;
t222 = pkin(2) * t155 - pkin(7) * t158;
t97 = t222 * qJD(2);
t274 = t148 + t144;
t99 = -pkin(1) - t274;
t306 = t154 * t97 + t99 * t263;
t239 = pkin(6) * t157 - qJ(4);
t178 = -pkin(4) * t290 - t239 * t155;
t94 = t222 * qJD(1);
t78 = t154 * t94;
t305 = -qJD(1) * t178 - t323 * t265 - t78;
t107 = t323 * t157;
t251 = -pkin(3) - t317;
t173 = pkin(4) * t285 + (-qJ(5) + t251) * t155;
t297 = t157 * t94;
t304 = -qJD(1) * t173 + qJD(3) * t107 + t297;
t303 = qJ(4) * t39;
t46 = t157 * t105 - t154 * t187;
t35 = t120 * qJ(4) - t46;
t302 = t120 * t35;
t301 = t120 * t45;
t300 = t120 * t46;
t294 = t90 * qJ(4);
t125 = pkin(6) * t285;
t60 = t154 * t99 + t125;
t293 = pkin(6) * qJDD(1);
t291 = t154 * t155;
t288 = t155 * t157;
t286 = t156 * t158;
t284 = t158 * qJ(5);
t283 = t158 * t159;
t282 = t159 * t154;
t281 = -qJD(4) - t45;
t31 = t46 - t318;
t279 = -qJD(5) - t31;
t278 = t337 * t154;
t277 = t337 * t157;
t275 = pkin(3) * t291 - qJ(4) * t288;
t273 = t159 * pkin(1) + t156 * pkin(6);
t151 = t155 ^ 2;
t152 = t158 ^ 2;
t272 = t151 - t152;
t271 = t151 + t152;
t267 = qJD(2) * t155;
t123 = pkin(6) * t290;
t161 = qJD(1) ^ 2;
t255 = t155 * t161 * t158;
t104 = -qJD(2) * pkin(2) + pkin(6) * t270;
t197 = -t92 * qJ(4) + t104;
t25 = t309 * t90 + t197;
t254 = t25 * t265;
t253 = t25 * t263;
t135 = pkin(6) * t260;
t73 = -qJDD(2) * pkin(2) + pkin(6) * t242 + t135;
t172 = t38 * qJ(4) - t92 * qJD(4) + t73;
t2 = t90 * qJD(5) + t309 * t39 + t172;
t250 = -t2 - t311;
t246 = t154 * t264;
t245 = t155 * t263;
t244 = t120 * t270;
t240 = -t154 * qJ(4) - pkin(2);
t74 = t154 * t286 + t157 * t159;
t75 = t156 * t285 - t282;
t238 = -t74 * pkin(3) + qJ(4) * t75;
t76 = -t156 * t157 + t158 * t282;
t77 = t154 * t156 + t157 * t283;
t237 = -t76 * pkin(3) + qJ(4) * t77;
t50 = qJD(1) * t97 - qJDD(1) * t203;
t236 = t105 * t265 - t154 * t50 + t340 * t157;
t235 = t105 * t263 - t340 * t154 - t157 * t50;
t59 = t157 * t99 - t123;
t232 = pkin(3) * t245 + pkin(6) * t266 + qJ(4) * t246 + qJD(2) * t259;
t231 = pkin(3) * t285 + qJ(4) * t290 + t274;
t230 = pkin(2) * t283 + pkin(7) * t287 + t273;
t229 = t155 * t242;
t228 = g(1) * t74 - g(2) * t76;
t227 = g(1) * t75 - g(2) * t77;
t149 = t159 * pkin(6);
t226 = -t75 * pkin(3) - t74 * qJ(4) + t149;
t225 = qJD(3) * t125 - t157 * t97 + t99 * t265;
t54 = qJ(4) * t158 - t60;
t224 = -t158 * qJD(4) + t306;
t223 = t251 * t155;
t218 = (qJD(3) * t90 - t38) * pkin(7);
t217 = (qJD(3) * t92 - t39) * pkin(7);
t216 = -qJDD(4) - t235;
t215 = pkin(6) * t90 + t104 * t154;
t214 = pkin(6) * t92 + t104 * t157;
t213 = qJD(3) * t104 - t321;
t18 = t309 * t120 + t280;
t24 = qJD(5) - t35 - t318;
t212 = t154 * t24 - t157 * t18;
t34 = pkin(3) * t120 - t281;
t211 = t154 * t35 + t157 * t34;
t209 = -t154 * t46 + t157 * t45;
t4 = t236 + t334;
t202 = pkin(3) * t157 - t240;
t200 = -qJ(5) * t291 - t275;
t199 = pkin(7) * qJD(3) * t120 - t311;
t194 = -pkin(6) * qJDD(2) + t261 * t325;
t5 = pkin(3) * t39 + t172;
t191 = t199 - t5;
t190 = t199 - t73;
t188 = pkin(1) * t161 + t221;
t186 = t309 * t157 - t240;
t160 = qJD(2) ^ 2;
t185 = pkin(6) * t160 + qJDD(1) * t325 + t313;
t184 = t77 * pkin(3) + qJ(4) * t76 + t230;
t182 = t203 * t316;
t37 = pkin(3) * t90 + t197;
t181 = t120 * t37 + t321;
t179 = t88 - t310;
t175 = g(1) * t76 + g(2) * t74 + g(3) * t291 - t235;
t174 = -t221 * t158 - t312;
t14 = -t157 * t298 - t296;
t15 = -t154 * t299 - t295;
t171 = -qJDD(4) + t175;
t170 = g(1) * t77 + g(2) * t75 + g(3) * t288 + t236;
t169 = t38 + t299;
t168 = t37 * t92 - t171;
t12 = t90 * t245 + (t155 * t39 + t90 * t266) * t154;
t11 = t92 * t158 * t262 + (-t38 * t157 - t92 * t265) * t155;
t167 = t25 * t92 - t171 - t320;
t166 = -t25 * t90 - t170 + t329;
t147 = t158 * pkin(3);
t145 = t155 * pkin(6);
t131 = g(1) * t289;
t128 = pkin(7) * t283;
t124 = pkin(7) * t286;
t106 = t323 * t154;
t65 = t145 + t275;
t57 = -pkin(6) * t247 + t78;
t56 = pkin(6) * t248 + t297;
t55 = t147 - t59;
t53 = t145 - t200;
t51 = -t120 * t267 - t158 * t88;
t49 = pkin(3) * t92 + t294;
t48 = qJD(1) * t223 - t297;
t47 = t239 * t270 - t78;
t43 = -pkin(4) * t291 - t54;
t40 = t284 + t123 + t147 + (pkin(4) * t155 - t99) * t157;
t32 = t309 * t92 + t294;
t29 = t267 * t317 - t225;
t28 = t306 - t335;
t27 = (-qJ(4) * t266 - qJD(4) * t155) * t157 + t232;
t26 = qJD(2) * t223 + t225;
t21 = -qJ(4) * t267 - t224 + t335;
t17 = qJD(2) * t189 + (qJD(5) * t154 + (qJ(5) * qJD(3) - qJD(4)) * t157) * t155 + t232;
t16 = (-pkin(4) * t288 - t123) * qJD(3) + t178 * qJD(2) + t224;
t13 = -pkin(4) * t246 + qJD(2) * t173 + t158 * qJD(5) + t225;
t8 = -t216 - t319;
t3 = -t4 + t329;
t1 = qJD(5) * t120 - t216 - t249 - t320;
t7 = [0, 0, 0, 0, 0, qJDD(1), -t313 + t316, t221, 0, 0, qJDD(1) * t151 + 0.2e1 * t229, 0.2e1 * t155 * t142 - 0.2e1 * t272 * t261, qJDD(2) * t155 + t158 * t160, qJDD(1) * t152 - 0.2e1 * t229, qJDD(2) * t158 - t155 * t160, 0, t194 * t155 + (-t185 + t316) * t158, t155 * t185 + t158 * t194 - t131, 0.2e1 * t271 * t293 - t221, -g(1) * (-t156 * pkin(1) + t149) - g(2) * t273 + (t271 * pkin(6) ^ 2 + pkin(1) ^ 2) * qJDD(1), t11, t163, t6, t12, t326, t51, -t29 * t120 + t59 * t88 + (qJD(2) * t215 + t235) * t158 + (pkin(6) * t39 - qJD(2) * t45 + t104 * t263 + t73 * t154) * t155 + t227, t28 * t120 - t60 * t88 + (qJD(2) * t214 - t236) * t158 + (-pkin(6) * t38 - qJD(2) * t46 - t104 * t265 + t73 * t157) * t155 - t228, -t28 * t90 - t29 * t92 + t59 * t38 - t60 * t39 + t131 + t209 * t266 + (-t313 + t235 * t157 + t154 * t236 + (-t154 * t45 - t157 * t46) * qJD(3)) * t155, -t236 * t60 + t46 * t28 - t235 * t59 - t45 * t29 - g(1) * t149 - g(2) * t230 + t182 + (t104 * t266 + t155 * t73) * pkin(6), t51, -t6, -t326, t11, t163, t12, t21 * t90 + t26 * t92 - t55 * t38 + t54 * t39 + t131 + t211 * t266 + (-t313 + t154 * t4 + t157 * t8 + (-t154 * t34 + t157 * t35) * qJD(3)) * t155, -t26 * t120 - t27 * t90 - t65 * t39 + t55 * t88 + (-t268 * t37 - t8) * t158 + (qJD(2) * t34 - t5 * t154 - t263 * t37) * t155 - t227, t21 * t120 - t27 * t92 + t65 * t38 - t54 * t88 + (-t262 * t37 + t4) * t158 + (-qJD(2) * t35 - t5 * t157 + t265 * t37) * t155 + t228, -g(1) * t226 - g(2) * t184 + t35 * t21 + t34 * t26 + t37 * t27 + t4 * t54 + t5 * t65 + t8 * t55 + t182, t51, -t326, t6, t12, -t163, t11, t13 * t92 - t16 * t90 - t40 * t38 - t43 * t39 + t131 - t212 * t266 + (-t313 + t1 * t157 - t154 * t3 + (-t154 * t18 - t157 * t24) * qJD(3)) * t155, -t16 * t120 - t17 * t92 + t53 * t38 + t43 * t88 + (-t25 * t262 - t3) * t158 + (qJD(2) * t24 - t2 * t157 + t254) * t155 + t228, t13 * t120 + t17 * t90 + t53 * t39 - t40 * t88 + (t25 * t268 + t1) * t158 + (-qJD(2) * t18 + t2 * t154 + t253) * t155 + t227, t2 * t53 + t25 * t17 + t1 * t40 + t18 * t13 + t3 * t43 + t24 * t16 - g(1) * (-t75 * qJ(5) + t226) - g(2) * (pkin(4) * t287 + qJ(5) * t77 + t184) - (-t323 * t155 + t252) * t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t255, t272 * t161, t260, t255, t142, qJDD(2), t155 * t188 - t135 - t311, t312 + (t188 - t293) * t158, 0, 0, t14, -t341, t20, t15, -t19, t244, -pkin(2) * t39 + t56 * t120 + t213 * t154 + t190 * t157 + (t155 * t45 - t158 * t215) * qJD(1) + t277, pkin(2) * t38 - t57 * t120 + t213 * t157 - t190 * t154 + (t155 * t46 - t158 * t214) * qJD(1) - t278, t56 * t92 + t57 * t90 + (t217 - t236 - t301) * t157 + (t235 + t218 + t300) * t154 + t174, -t73 * pkin(2) - t46 * t57 + t45 * t56 - t104 * t137 - g(1) * (-pkin(2) * t287 + t128) - g(2) * (-pkin(2) * t289 + t124) - g(3) * t274 + (qJD(3) * t209 + t154 * t235 - t157 * t236) * pkin(7), t244, -t20, t19, t14, -t341, t15, -t47 * t90 - t48 * t92 + (-t120 * t34 + t217 - t4) * t157 + (t218 + t8 - t302) * t154 + t174, t48 * t120 + t154 * t181 - t157 * t191 + t202 * t39 - t270 * t34 - t307 * t90 - t277, -t47 * t120 + t154 * t191 + t157 * t181 - t202 * t38 + t270 * t35 - t307 * t92 + t278, -t35 * t47 - t34 * t48 - g(1) * t128 - g(2) * t124 - g(3) * t231 + t307 * t37 + (qJD(3) * t211 + t8 * t154 - t4 * t157) * pkin(7) + (t221 * t155 - t5) * t202, t244, t19, t20, t15, t341, t14, -t312 + t1 * t154 - t106 * t38 - t107 * t39 + t157 * t3 + t304 * t92 - t305 * t90 - t212 * qJD(3) + (qJD(1) * t212 - t221) * t158, -t253 + t107 * t88 - t186 * t38 - t308 * t92 + t250 * t154 - t305 * t120 + (-t155 * t24 + t25 * t285) * qJD(1) + t278, t254 - t106 * t88 - t186 * t39 + t308 * t90 + t250 * t157 + t304 * t120 + (t155 * t18 - t25 * t290) * qJD(1) + t277, -t2 * t186 + t1 * t106 + t3 * t107 - g(1) * (pkin(4) * t283 + t128) - g(2) * (pkin(4) * t286 + t124) - g(3) * (t157 * t284 + t231) + t308 * t25 + t305 * t24 + t304 * t18 + (-g(3) * pkin(4) + t221 * t186) * t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t310, -t256, -t169, -t310, -t23, t88, -t104 * t92 + t175 - t300, t104 * t90 + t170 + t301, 0, 0, t88, t169, t23, t310, -t256, -t310, pkin(3) * t38 - t303 + (-t35 - t46) * t92 + (t34 + t281) * t90, t49 * t90 + t168 + t300 - 0.2e1 * t319, t120 * t281 - t37 * t90 + t49 * t92 - t108 - t170 + t324, -t8 * pkin(3) - g(1) * t237 - g(2) * t238 + g(3) * t275 - t4 * qJ(4) + t281 * t35 - t34 * t46 - t37 * t49, t88, t23, -t169, -t310, t256, t310, -t303 + t309 * t38 + (t24 + t279) * t92 + (t18 - t280) * t90, -t120 * t201 + t32 * t92 - 0.2e1 * t108 + t166 + t324, -t32 * t90 + (-0.2e1 * qJD(5) - t31) * t120 + 0.2e1 * t249 - t167, -t1 * t309 + t3 * qJ(4) - t25 * t32 - g(1) * (-qJ(5) * t76 + t237) - g(2) * (-qJ(5) * t74 + t238) - g(3) * t200 + t280 * t24 + t279 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169, t179, t333, t168 - t302 - t319, 0, 0, 0, 0, 0, 0, -t169, t333, -t179, -t249 + (qJD(5) + t24) * t120 + t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t88 + t310, -t86 - t116, -t120 * t18 + t166 - t334;];
tau_reg = t7;
