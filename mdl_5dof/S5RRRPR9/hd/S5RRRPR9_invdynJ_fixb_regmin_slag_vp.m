% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:23:55
% EndTime: 2021-01-15 23:24:21
% DurationCPUTime: 5.70s
% Computational Cost: add. (5155->465), mult. (11820->640), div. (0->0), fcn. (8459->14), ass. (0->225)
t227 = cos(qJ(2));
t286 = qJD(1) * t227;
t193 = -qJD(3) + t286;
t185 = -qJD(5) + t193;
t221 = sin(qJ(5));
t222 = sin(qJ(3));
t223 = sin(qJ(2));
t287 = qJD(1) * t223;
t269 = t222 * t287;
t226 = cos(qJ(3));
t276 = t226 * qJD(2);
t169 = t269 - t276;
t283 = qJD(2) * t222;
t171 = t226 * t287 + t283;
t218 = sin(pkin(9));
t219 = cos(pkin(9));
t246 = -t169 * t218 + t219 * t171;
t225 = cos(qJ(5));
t98 = t219 * t169 + t171 * t218;
t307 = t225 * t98;
t49 = t221 * t246 + t307;
t310 = t185 * t49;
t277 = qJD(5) * t221;
t274 = qJD(1) * qJD(2);
t262 = t227 * t274;
t273 = t223 * qJDD(1);
t279 = qJD(3) * t223;
t331 = -qJD(1) * t279 + qJDD(2);
t90 = qJD(3) * t276 + (t262 + t273) * t226 + t331 * t222;
t91 = t222 * (qJD(2) * (qJD(3) + t286) + t273) - t331 * t226;
t38 = t218 * t90 + t219 * t91;
t39 = -t218 * t91 + t219 * t90;
t8 = -qJD(5) * t307 - t221 * t38 + t225 * t39 - t246 * t277;
t343 = t8 - t310;
t330 = -t221 * t98 + t225 * t246;
t342 = t330 * t49;
t341 = t330 ^ 2 - t49 ^ 2;
t215 = qJ(3) + pkin(9);
t212 = qJ(5) + t215;
t197 = sin(t212);
t198 = cos(t212);
t228 = cos(qJ(1));
t224 = sin(qJ(1));
t297 = t224 * t227;
t118 = t197 * t228 - t198 * t297;
t295 = t227 * t228;
t120 = t197 * t224 + t198 * t295;
t178 = -pkin(2) * t227 - pkin(7) * t223 - pkin(1);
t158 = t178 * qJD(1);
t205 = pkin(6) * t286;
t182 = qJD(2) * pkin(7) + t205;
t108 = t158 * t222 + t182 * t226;
t72 = -qJ(4) * t169 + t108;
t309 = t219 * t72;
t107 = t226 * t158 - t182 * t222;
t71 = -qJ(4) * t171 + t107;
t63 = -pkin(3) * t193 + t71;
t28 = t218 * t63 + t309;
t333 = pkin(8) * t98;
t20 = t28 - t333;
t18 = t20 * t277;
t211 = t227 * qJDD(1);
t326 = -t223 * t274 + t211;
t161 = qJDD(3) - t326;
t251 = pkin(2) * t223 - pkin(7) * t227;
t173 = t251 * qJD(2);
t109 = qJD(1) * t173 + qJDD(1) * t178;
t103 = t226 * t109;
t141 = pkin(6) * t326 + qJDD(2) * pkin(7);
t16 = t161 * pkin(3) - t90 * qJ(4) - qJD(3) * t108 - t171 * qJD(4) - t222 * t141 + t103;
t278 = qJD(3) * t226;
t280 = qJD(3) * t222;
t238 = t222 * t109 + t226 * t141 + t158 * t278 - t182 * t280;
t21 = -qJ(4) * t91 - qJD(4) * t169 + t238;
t4 = t219 * t16 - t218 * t21;
t2 = pkin(4) * t161 - pkin(8) * t39 + t4;
t316 = g(3) * t223;
t181 = -qJD(2) * pkin(2) + pkin(6) * t287;
t113 = pkin(3) * t169 + qJD(4) + t181;
t57 = pkin(4) * t98 + t113;
t340 = g(1) * t120 - g(2) * t118 + t198 * t316 - t221 * t2 + t57 * t49 + t18;
t306 = t330 * t185;
t9 = qJD(5) * t330 + t221 * t39 + t225 * t38;
t338 = -t9 - t306;
t296 = t226 * t227;
t320 = pkin(3) * t223;
t245 = -qJ(4) * t296 + t320;
t220 = qJ(4) + pkin(7);
t258 = qJD(3) * t220;
t172 = t251 * qJD(1);
t290 = pkin(6) * t269 + t226 * t172;
t337 = -qJD(1) * t245 - t222 * qJD(4) - t226 * t258 - t290;
t153 = t222 * t172;
t275 = t226 * qJD(4);
t298 = t223 * t226;
t299 = t222 * t227;
t336 = t153 + (-pkin(6) * t298 - qJ(4) * t299) * qJD(1) + t222 * t258 - t275;
t203 = pkin(6) * t273;
t142 = -qJDD(2) * pkin(2) + pkin(6) * t262 + t203;
t250 = g(1) * t228 + g(2) * t224;
t315 = g(3) * t227;
t234 = t223 * t250 - t315;
t335 = qJD(3) * pkin(7) * t193 - t142 + t234;
t117 = t197 * t297 + t198 * t228;
t119 = -t197 * t295 + t198 * t224;
t5 = t218 * t16 + t219 * t21;
t3 = -pkin(8) * t38 + t5;
t270 = t225 * t2 - t221 * t3;
t334 = -g(1) * t119 + g(2) * t117 + t197 * t316 - t57 * t330 + t270;
t163 = t218 * t226 + t219 * t222;
t146 = t163 * qJD(3);
t294 = t163 * t286 - t146;
t162 = t218 * t222 - t219 * t226;
t332 = t193 * t162;
t329 = pkin(8) * t246;
t312 = t336 * t218 + t337 * t219;
t311 = t337 * t218 - t336 * t219;
t321 = pkin(3) * t222;
t252 = pkin(3) * t280 - t286 * t321 - t205;
t149 = t222 * t297 + t226 * t228;
t151 = -t222 * t295 + t224 * t226;
t324 = -g(1) * t151 + g(2) * t149;
t322 = pkin(3) * t218;
t319 = pkin(6) * t222;
t95 = t225 * t162 + t163 * t221;
t314 = -qJD(5) * t95 + t294 * t221 + t332 * t225;
t96 = -t162 * t221 + t163 * t225;
t313 = qJD(5) * t96 + t332 * t221 - t294 * t225;
t195 = pkin(6) * t296;
t282 = qJD(2) * t223;
t291 = t226 * t173 + t282 * t319;
t44 = -t223 * t275 + t245 * qJD(2) + (-t195 + (qJ(4) * t223 - t178) * t222) * qJD(3) + t291;
t292 = t222 * t173 + t178 * t278;
t53 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t298 + (-qJD(4) * t223 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t227) * t222 + t292;
t15 = t218 * t44 + t219 * t53;
t67 = t218 * t72;
t32 = t219 * t71 - t67;
t27 = t219 * t63 - t67;
t17 = -pkin(4) * t193 + t27 - t329;
t308 = t225 * t17;
t305 = t90 * t222;
t165 = t226 * t178;
t104 = -qJ(4) * t298 + t165 + (-pkin(3) - t319) * t227;
t289 = t222 * t178 + t195;
t300 = t222 * t223;
t110 = -qJ(4) * t300 + t289;
t55 = t218 * t104 + t219 * t110;
t304 = -t294 * pkin(4) + t252;
t303 = t169 * t193;
t302 = t171 * t193;
t301 = t171 * t226;
t179 = t220 * t222;
t180 = t220 * t226;
t112 = -t218 * t179 + t219 * t180;
t174 = pkin(3) * t300 + t223 * pkin(6);
t216 = t223 ^ 2;
t288 = -t227 ^ 2 + t216;
t285 = qJD(2) * t169;
t284 = qJD(2) * t171;
t281 = qJD(2) * t227;
t267 = t222 * t281;
t116 = pkin(3) * t267 + pkin(6) * t281 + t278 * t320;
t202 = pkin(3) * t226 + pkin(2);
t268 = t193 * t276;
t266 = t227 * t276;
t265 = t193 * t280;
t264 = t193 * t278;
t260 = qJD(5) * t17 + t3;
t14 = -t218 * t53 + t219 * t44;
t31 = -t218 * t71 - t309;
t54 = t219 * t104 - t110 * t218;
t111 = -t219 * t179 - t180 * t218;
t257 = t202 * t227 + t220 * t223;
t256 = -qJD(3) * t158 - t141;
t83 = -pkin(8) * t162 + t112;
t254 = pkin(4) * t287 + t332 * pkin(8) + qJD(5) * t83 - t312;
t82 = -pkin(8) * t163 + t111;
t253 = t294 * pkin(8) + qJD(5) * t82 + t311;
t249 = g(1) * t224 - g(2) * t228;
t248 = t182 * t278 - t103;
t7 = t221 * t17 + t225 * t20;
t247 = -pkin(7) * t161 + qJD(3) * t181;
t135 = t163 * t223;
t136 = t162 * t223;
t74 = t225 * t135 - t136 * t221;
t75 = -t135 * t221 - t136 * t225;
t199 = pkin(3) * t219 + pkin(4);
t244 = t199 * t221 + t225 * t322;
t243 = t199 * t225 - t221 * t322;
t241 = -0.2e1 * pkin(1) * t274 - pkin(6) * qJDD(2);
t240 = t222 * t161 - t264;
t239 = t226 * t161 + t265;
t230 = qJD(1) ^ 2;
t237 = pkin(1) * t230 + t250;
t62 = pkin(3) * t91 + qJDD(4) + t142;
t229 = qJD(2) ^ 2;
t232 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t229 + t249;
t208 = cos(t215);
t207 = sin(t215);
t201 = pkin(6) + t321;
t157 = qJDD(5) + t161;
t152 = t222 * t224 + t226 * t295;
t150 = t222 * t228 - t224 * t296;
t140 = pkin(1) + t257;
t132 = t207 * t224 + t208 * t295;
t131 = -t207 * t295 + t208 * t224;
t130 = t207 * t228 - t208 * t297;
t129 = t207 * t297 + t208 * t228;
t124 = pkin(4) * t162 - t202;
t105 = pkin(4) * t135 + t174;
t77 = t146 * t223 + t218 * t267 - t219 * t266;
t76 = t162 * t279 - t163 * t281;
t64 = pkin(3) * t171 + pkin(4) * t246;
t56 = -pkin(4) * t76 + t116;
t34 = -pkin(8) * t135 + t55;
t33 = -pkin(4) * t227 + pkin(8) * t136 + t54;
t26 = t32 - t329;
t25 = t31 + t333;
t24 = qJD(5) * t75 - t221 * t77 - t225 * t76;
t23 = -qJD(5) * t74 + t221 * t76 - t225 * t77;
t22 = pkin(4) * t38 + t62;
t11 = pkin(8) * t76 + t15;
t10 = pkin(4) * t282 + pkin(8) * t77 + t14;
t6 = -t20 * t221 + t308;
t1 = [qJDD(1), t249, t250, qJDD(1) * t216 + 0.2e1 * t223 * t262, 0.2e1 * t211 * t223 - 0.2e1 * t274 * t288, qJDD(2) * t223 + t227 * t229, qJDD(2) * t227 - t223 * t229, 0, t223 * t241 + t227 * t232, -t223 * t232 + t227 * t241, t90 * t298 + (-t222 * t279 + t266) * t171, (-t169 * t226 - t171 * t222) * t281 + (-t305 - t226 * t91 + (t169 * t222 - t301) * qJD(3)) * t223, (-t90 - t268) * t227 + (t239 + t284) * t223, (t193 * t283 + t91) * t227 + (-t240 - t285) * t223, -t161 * t227 - t193 * t282, -(-t178 * t280 + t291) * t193 + t165 * t161 - g(1) * t150 - g(2) * t152 + ((t264 + t285) * pkin(6) + (-pkin(6) * t161 + qJD(2) * t181 - t256) * t222 + t248) * t227 + (pkin(6) * t91 + qJD(2) * t107 + t142 * t222 + t181 * t278) * t223, t292 * t193 - t289 * t161 - g(1) * t149 - g(2) * t151 + (t181 * t276 + (-t265 + t284) * pkin(6) + t238) * t227 + (-t181 * t280 - t108 * qJD(2) + t142 * t226 + (t90 - t268) * pkin(6)) * t223, -g(1) * t130 - g(2) * t132 - t113 * t76 + t116 * t98 + t135 * t62 - t14 * t193 + t161 * t54 + t174 * t38 - t227 * t4 + t27 * t282, -g(1) * t129 - g(2) * t131 - t113 * t77 + t116 * t246 - t136 * t62 + t15 * t193 - t161 * t55 + t174 * t39 + t227 * t5 - t28 * t282, -t135 * t5 + t136 * t4 - t14 * t246 - t15 * t98 + t223 * t249 + t27 * t77 + t28 * t76 - t38 * t55 - t39 * t54, t5 * t55 + t28 * t15 + t4 * t54 + t27 * t14 + t62 * t174 + t113 * t116 - g(1) * (-t140 * t224 + t201 * t228) - g(2) * (t140 * t228 + t201 * t224), t23 * t330 + t75 * t8, -t23 * t49 - t24 * t330 - t74 * t8 - t75 * t9, t157 * t75 - t185 * t23 - t227 * t8 + t282 * t330, -t157 * t74 + t185 * t24 + t227 * t9 - t282 * t49, -t157 * t227 - t185 * t282, -(t225 * t10 - t221 * t11) * t185 + (-t221 * t34 + t225 * t33) * t157 - t270 * t227 + t6 * t282 + t56 * t49 + t105 * t9 + t22 * t74 + t57 * t24 - g(1) * t118 - g(2) * t120 + (-(-t221 * t33 - t225 * t34) * t185 + t7 * t227) * qJD(5), -t7 * t282 - g(1) * t117 - g(2) * t119 + t105 * t8 - t18 * t227 + t22 * t75 + t57 * t23 + t56 * t330 + ((-qJD(5) * t34 + t10) * t185 - t33 * t157 + t2 * t227) * t221 + ((qJD(5) * t33 + t11) * t185 - t34 * t157 + t260 * t227) * t225; 0, 0, 0, -t223 * t230 * t227, t288 * t230, t273, t211, qJDD(2), t223 * t237 - t203 - t315, t316 + (-pkin(6) * qJDD(1) + t237) * t227, -t193 * t301 + t305, (t90 + t303) * t226 + (-t91 + t302) * t222, (-t171 * t223 + t193 * t296) * qJD(1) + t240, (t169 * t223 - t193 * t299) * qJD(1) + t239, t193 * t287, -pkin(2) * t91 + t290 * t193 + t247 * t222 + (-t107 * t223 + (-pkin(6) * t169 - t181 * t222) * t227) * qJD(1) + t335 * t226, -pkin(2) * t90 - t153 * t193 + t247 * t226 + (-t181 * t296 + t108 * t223 + (-t171 * t227 + t193 * t298) * pkin(6)) * qJD(1) - t335 * t222, t111 * t161 - t113 * t294 + t162 * t62 - t193 * t312 - t202 * t38 + t208 * t234 + t252 * t98 - t27 * t287, -t112 * t161 + t113 * t332 + t163 * t62 + t193 * t311 - t202 * t39 - t207 * t234 + t246 * t252 + t28 * t287, -t111 * t39 - t112 * t38 - t162 * t5 - t163 * t4 - t227 * t250 - t246 * t312 - t27 * t332 + t28 * t294 - t311 * t98 - t316, t5 * t112 + t4 * t111 - t62 * t202 - g(3) * t257 + t311 * t28 + t312 * t27 - t250 * (-t202 * t223 + t220 * t227) + t252 * t113, t314 * t330 + t8 * t96, -t313 * t330 - t314 * t49 - t8 * t95 - t96 * t9, t96 * t157 - t314 * t185 - t287 * t330, -t95 * t157 + t313 * t185 + t49 * t287, t185 * t287, (-t221 * t83 + t225 * t82) * t157 + t124 * t9 + t22 * t95 - t6 * t287 + t313 * t57 + t304 * t49 + (t221 * t253 + t225 * t254) * t185 + t234 * t198, -(t221 * t82 + t225 * t83) * t157 + t124 * t8 + t22 * t96 + t7 * t287 + t314 * t57 + t304 * t330 + (-t221 * t254 + t225 * t253) * t185 - t234 * t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171 * t169, -t169 ^ 2 + t171 ^ 2, t90 - t303, -t302 - t91, t161, -t108 * t193 - t181 * t171 + (t256 + t316) * t222 - t248 + t324, g(1) * t152 - g(2) * t150 + g(3) * t298 - t107 * t193 + t169 * t181 - t238, t207 * t316 - g(1) * t131 + g(2) * t129 - t113 * t246 + t31 * t193 + (t161 * t219 - t171 * t98) * pkin(3) + t4, t208 * t316 + g(1) * t132 - g(2) * t130 + t113 * t98 - t32 * t193 + (-t161 * t218 - t171 * t246) * pkin(3) - t5, (-t218 * t38 - t219 * t39) * pkin(3) + (t31 + t28) * t246 + (-t27 + t32) * t98, -t27 * t31 - t28 * t32 + (g(3) * t300 - t113 * t171 + t5 * t218 + t4 * t219 + t324) * pkin(3), t342, t341, t343, t338, t157, t243 * t157 + (-t221 * t26 + t225 * t25) * t185 - t64 * t49 + (t185 * t244 - t7) * qJD(5) + t334, -t244 * t157 - t225 * t3 - (t221 * t25 + t225 * t26) * t185 - t64 * t330 + (t185 * t243 - t308) * qJD(5) + t340; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t193 * t246 + t38, t193 * t98 + t39, -t246 ^ 2 - t98 ^ 2, t246 * t27 + t28 * t98 - t234 + t62, 0, 0, 0, 0, 0, t9 - t306, t8 + t310; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t342, t341, t343, t338, t157, (-qJD(5) - t185) * t7 + t334, -t6 * t185 - t225 * t260 + t340;];
tau_reg = t1;
