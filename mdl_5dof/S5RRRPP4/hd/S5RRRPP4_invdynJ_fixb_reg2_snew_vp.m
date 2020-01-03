% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:47
% EndTime: 2019-12-31 20:55:59
% DurationCPUTime: 5.93s
% Computational Cost: add. (14049->363), mult. (32103->463), div. (0->0), fcn. (22692->8), ass. (0->230)
t229 = sin(qJ(2));
t232 = cos(qJ(2));
t228 = sin(qJ(3));
t231 = cos(qJ(3));
t198 = (t228 * t232 + t229 * t231) * qJD(1);
t216 = t229 * qJDD(1);
t264 = qJD(1) * qJD(2);
t259 = t232 * t264;
t204 = t216 + t259;
t217 = t232 * qJDD(1);
t260 = t229 * t264;
t205 = t217 - t260;
t252 = t204 * t228 - t231 * t205;
t154 = -qJD(3) * t198 - t252;
t268 = qJD(1) * t229;
t196 = -t231 * t232 * qJD(1) + t228 * t268;
t248 = t204 * t231 + t205 * t228;
t155 = -qJD(3) * t196 + t248;
t226 = sin(pkin(8));
t227 = cos(pkin(8));
t117 = t154 * t226 + t155 * t227;
t177 = t227 * t196 + t198 * t226;
t223 = qJD(2) + qJD(3);
t276 = t223 * t177;
t313 = t117 - t276;
t222 = qJDD(2) + qJDD(3);
t179 = -t196 * t226 + t198 * t227;
t280 = t179 * t177;
t132 = -t280 - t222;
t286 = t132 * t226;
t175 = t179 ^ 2;
t303 = t223 ^ 2;
t310 = -t175 - t303;
t70 = -t227 * t310 - t286;
t285 = t132 * t227;
t72 = t226 * t310 - t285;
t41 = t228 * t72 + t231 * t70;
t59 = t228 * t70 - t231 * t72;
t367 = pkin(1) * t313 - pkin(6) * (t229 * t41 + t232 * t59);
t366 = pkin(2) * t41;
t365 = pkin(7) * t41;
t364 = pkin(2) * t313 - pkin(7) * t59;
t304 = t177 ^ 2;
t162 = t304 - t303;
t104 = -t162 * t226 + t285;
t108 = -t162 * t227 - t286;
t362 = t229 * (t104 * t228 - t108 * t231) - t232 * (t104 * t231 + t108 * t228);
t361 = pkin(3) * t70;
t360 = qJ(4) * t70;
t359 = qJ(4) * t72;
t168 = t223 * t179;
t254 = -t227 * t154 + t155 * t226;
t311 = t168 + t254;
t50 = -t226 * t311 + t227 * t313;
t293 = t226 * t313;
t52 = t227 * t311 + t293;
t356 = t229 * (t228 * t50 + t231 * t52) + t232 * (t228 * t52 - t231 * t50);
t312 = t117 + t276;
t87 = -t254 + t168;
t327 = t226 * t312 + t227 * t87;
t328 = t226 * t87 - t227 * t312;
t341 = t228 * t327 + t231 * t328;
t355 = pkin(2) * t341;
t354 = pkin(7) * t341;
t111 = -t304 - t175;
t342 = -t228 * t328 + t231 * t327;
t353 = -pkin(2) * t111 + pkin(7) * t342;
t350 = pkin(6) * (-t229 * t341 + t232 * t342) - pkin(1) * t111;
t308 = -t280 + t222;
t284 = t308 * t226;
t307 = -t303 - t304;
t314 = t227 * t307 - t284;
t122 = t227 * t308;
t315 = t226 * t307 + t122;
t325 = t228 * t314 + t231 * t315;
t349 = pkin(2) * t325;
t47 = pkin(3) * t328;
t348 = pkin(7) * t325;
t345 = qJ(4) * t328;
t326 = -t228 * t315 + t231 * t314;
t344 = -pkin(2) * t311 + pkin(7) * t326;
t343 = -pkin(3) * t111 + qJ(4) * t327;
t340 = pkin(6) * (-t229 * t325 + t232 * t326) - pkin(1) * t311;
t163 = -t175 + t303;
t330 = t227 * t163 + t284;
t331 = -t163 * t226 + t122;
t339 = t229 * (-t228 * t330 + t231 * t331) + t232 * (t228 * t331 + t231 * t330);
t335 = qJ(4) * t314;
t334 = qJ(4) * t315;
t266 = qJD(4) * t179;
t329 = pkin(3) * t315 - 0.2e1 * t266;
t185 = t198 * t196;
t309 = -t185 + t222;
t318 = t228 * t309;
t317 = t231 * t309;
t316 = t313 * qJ(5);
t192 = t223 * t196;
t146 = t155 + t192;
t134 = t175 - t304;
t234 = qJD(1) ^ 2;
t270 = t229 * t234;
t230 = sin(qJ(1));
t301 = cos(qJ(1));
t246 = t301 * g(1) + t230 * g(2);
t287 = qJDD(1) * pkin(6);
t200 = -t234 * pkin(1) - t246 + t287;
t277 = t200 * t229;
t150 = qJDD(2) * pkin(2) - pkin(7) * t204 - t277 + (pkin(2) * t270 + pkin(7) * t264 - g(3)) * t232;
t189 = -t229 * g(3) + t232 * t200;
t225 = t232 ^ 2;
t219 = t225 * t234;
t244 = qJD(2) * pkin(2) - pkin(7) * t268;
t151 = -pkin(2) * t219 + t205 * pkin(7) - qJD(2) * t244 + t189;
t113 = -t231 * t150 + t151 * t228;
t65 = t309 * pkin(3) - t146 * qJ(4) - t113;
t114 = t228 * t150 + t231 * t151;
t194 = t196 ^ 2;
t249 = pkin(3) * t223 - qJ(4) * t198;
t67 = -t194 * pkin(3) + t154 * qJ(4) - t223 * t249 + t114;
t39 = -0.2e1 * qJD(4) * t177 + t226 * t65 + t227 * t67;
t240 = (-t177 * t226 - t179 * t227) * t223;
t275 = t223 * t226;
t161 = t179 * t275;
t274 = t223 * t227;
t261 = t177 * t274;
t250 = t161 - t261;
t306 = t229 * (-t228 * t240 + t231 * t250) + t232 * (t228 * t250 + t231 * t240);
t242 = t226 * t254 + t261;
t251 = t177 * t275 - t227 * t254;
t305 = t229 * (-t228 * t251 + t231 * t242) + t232 * (t228 * t242 + t231 * t251);
t195 = t198 ^ 2;
t302 = 2 * qJD(5);
t300 = pkin(4) * t254;
t299 = pkin(4) * t227;
t256 = t226 * t67 - t227 * t65;
t38 = t256 + 0.2e1 * t266;
t19 = t226 * t39 - t227 * t38;
t297 = t19 * t228;
t296 = t19 * t231;
t257 = g(1) * t230 - t301 * g(2);
t243 = qJDD(1) * pkin(1) + t257;
t157 = pkin(2) * t205 - t244 * t268 + (pkin(7) * t225 + pkin(6)) * t234 + t243;
t74 = pkin(3) * t154 + qJ(4) * t194 - t198 * t249 - qJDD(4) + t157;
t295 = t226 * t74;
t291 = t227 * t74;
t61 = -t113 * t231 + t114 * t228;
t289 = t229 * t61;
t288 = qJ(5) * t227;
t283 = t157 * t228;
t282 = t157 * t231;
t182 = t185 + t222;
t279 = t182 * t228;
t278 = t182 * t231;
t273 = t223 * t228;
t272 = t223 * t231;
t211 = t232 * t270;
t271 = t229 * (qJDD(2) + t211);
t269 = t232 * (qJDD(2) - t211);
t265 = qJD(3) + t223;
t133 = pkin(4) * t177 - qJ(5) * t179;
t245 = t222 * qJ(5) - t177 * t133 + t223 * t302 + t39;
t28 = -pkin(4) * t303 + t245;
t241 = -t222 * pkin(4) - qJ(5) * t303 + qJDD(5) + t256;
t30 = (0.2e1 * qJD(4) + t133) * t179 + t241;
t12 = t226 * t28 - t227 * t30;
t263 = pkin(3) * t12 - pkin(4) * t30 + qJ(5) * t28;
t262 = -pkin(4) * t312 + qJ(5) * t87 + t47;
t258 = -qJ(5) * t226 - pkin(3);
t20 = t226 * t38 + t227 * t39;
t255 = -t39 - t361;
t62 = t113 * t228 + t231 * t114;
t188 = g(3) * t232 + t277;
t253 = t229 * t188 + t232 * t189;
t247 = -t256 + t329;
t239 = (-qJD(3) + t223) * t198 - t252;
t238 = -pkin(4) * t310 - qJ(5) * t132 + t28 + t361;
t237 = pkin(4) * t308 + qJ(5) * t307 - t133 * t179 - t241 + t329;
t236 = -pkin(4) * t168 + t179 * t302 + t74;
t235 = t236 + t316;
t233 = qJD(2) ^ 2;
t224 = t229 ^ 2;
t218 = t224 * t234;
t206 = t217 - 0.2e1 * t260;
t203 = t216 + 0.2e1 * t259;
t199 = pkin(6) * t234 + t243;
t191 = -t195 + t303;
t190 = t194 - t303;
t187 = -t195 - t303;
t184 = t195 - t194;
t180 = -t303 - t194;
t156 = -t194 - t195;
t148 = -t187 * t228 - t278;
t147 = t187 * t231 - t279;
t145 = t155 - t192;
t144 = -t265 * t196 + t248;
t141 = t265 * t198 + t252;
t137 = t180 * t231 - t318;
t136 = t180 * t228 + t317;
t94 = t146 * t228 + t231 * t239;
t93 = -t146 * t231 + t228 * t239;
t80 = t117 * t227 - t161;
t79 = t117 * t226 + t179 * t274;
t56 = -t291 + t360;
t45 = -t295 - t334;
t40 = -pkin(3) * t313 - t295 - t359;
t37 = -pkin(3) * t311 + t291 + t335;
t35 = t235 - t300;
t26 = t229 * (-t228 * t79 + t231 * t80) + t232 * (t228 * t80 + t231 * t79);
t25 = t235 + (-t254 - t311) * pkin(4);
t24 = t236 - t300 + 0.2e1 * t316;
t23 = -qJ(5) * t111 + t30;
t22 = (-t111 - t303) * pkin(4) + t245;
t21 = -t226 * t25 - t288 * t311 - t334;
t18 = pkin(3) * t19;
t17 = -pkin(4) * t293 + t227 * t24 - t360;
t16 = t227 * t25 + t258 * t311 + t335;
t15 = t359 + t226 * t24 + (pkin(3) + t299) * t313;
t14 = pkin(3) * t74 + qJ(4) * t20;
t13 = t226 * t30 + t227 * t28;
t10 = -t19 - t345;
t9 = t20 + t343;
t8 = -t22 * t226 + t227 * t23 - t345;
t7 = t22 * t227 + t226 * t23 + t343;
t6 = t20 * t231 - t297;
t5 = t20 * t228 + t296;
t4 = -qJ(4) * t12 + (-pkin(4) * t226 + t288) * t35;
t3 = -t12 * t228 + t13 * t231;
t2 = t12 * t231 + t13 * t228;
t1 = qJ(4) * t13 + (-t258 + t299) * t35;
t11 = [0, 0, 0, 0, 0, qJDD(1), t257, t246, 0, 0, (t204 + t259) * t229, t203 * t232 + t206 * t229, t271 + t232 * (-t218 + t233), (t205 - t260) * t232, t229 * (t219 - t233) + t269, 0, t232 * t199 + pkin(1) * t206 + pkin(6) * (t232 * (-t219 - t233) - t271), -t229 * t199 - pkin(1) * t203 + pkin(6) * (-t269 - t229 * (-t218 - t233)), pkin(1) * (t218 + t219) + (t224 + t225) * t287 + t253, pkin(1) * t199 + pkin(6) * t253, t229 * (t155 * t231 - t198 * t273) + t232 * (t155 * t228 + t198 * t272), t229 * (-t141 * t231 - t145 * t228) + t232 * (-t141 * t228 + t145 * t231), t229 * (-t191 * t228 + t317) + t232 * (t191 * t231 + t318), t229 * (-t154 * t228 + t196 * t272) + t232 * (t154 * t231 + t196 * t273), t229 * (t190 * t231 - t279) + t232 * (t190 * t228 + t278), (t229 * (-t196 * t231 + t198 * t228) + t232 * (-t196 * t228 - t198 * t231)) * t223, t229 * (-pkin(7) * t136 - t283) + t232 * (-pkin(2) * t141 + pkin(7) * t137 + t282) - pkin(1) * t141 + pkin(6) * (-t136 * t229 + t137 * t232), t229 * (-pkin(7) * t147 - t282) + t232 * (-pkin(2) * t144 + pkin(7) * t148 - t283) - pkin(1) * t144 + pkin(6) * (-t147 * t229 + t148 * t232), t229 * (-pkin(7) * t93 - t61) + t232 * (-pkin(2) * t156 + pkin(7) * t94 + t62) - pkin(1) * t156 + pkin(6) * (-t229 * t93 + t232 * t94), -pkin(7) * t289 + t232 * (pkin(2) * t157 + pkin(7) * t62) + pkin(1) * t157 + pkin(6) * (t232 * t62 - t289), t26, -t356, t339, t305, t362, t306, t229 * (-t228 * t37 + t231 * t45 - t348) + t232 * (t228 * t45 + t231 * t37 + t344) + t340, t229 * (-t228 * t40 + t231 * t56 + t365) + t232 * (t228 * t56 + t231 * t40 - t364) - t367, t229 * (t10 * t231 - t228 * t9 - t354) + t232 * (t10 * t228 + t231 * t9 + t353) + t350, t229 * (-pkin(7) * t5 - qJ(4) * t296 - t14 * t228) + t232 * (pkin(2) * t74 + pkin(7) * t6 - qJ(4) * t297 + t14 * t231) + pkin(1) * t74 + pkin(6) * (-t229 * t5 + t232 * t6), t26, t339, t356, t306, -t362, t305, t229 * (-t16 * t228 + t21 * t231 - t348) + t232 * (t16 * t231 + t21 * t228 + t344) + t340, t229 * (-t228 * t7 + t231 * t8 - t354) + t232 * (t228 * t8 + t231 * t7 + t353) + t350, t229 * (-t15 * t228 + t17 * t231 - t365) + t232 * (t15 * t231 + t17 * t228 + t364) + t367, t229 * (-pkin(7) * t2 - t1 * t228 + t231 * t4) + t232 * (pkin(2) * t35 + pkin(7) * t3 + t1 * t231 + t228 * t4) + pkin(1) * t35 + pkin(6) * (-t2 * t229 + t232 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t211, -t219 + t218, t216, t211, t217, qJDD(2), -t188, -t189, 0, 0, t185, t184, t146, -t185, t239, t222, pkin(2) * t136 - t113, pkin(2) * t147 - t114, pkin(2) * t93, pkin(2) * t61, t280, t134, t312, -t280, t87, t222, t247 + t349, t255 - t366, t47 + t355, pkin(2) * t5 + t18, t280, t312, -t134, t222, -t87, -t280, t237 + t349, t262 + t355, t238 + t366, pkin(2) * t2 + t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, t184, t146, -t185, t239, t222, -t113, -t114, 0, 0, t280, t134, t312, -t280, t87, t222, t247, t255, t47, t18, t280, t312, -t134, t222, -t87, -t280, t237, t262, t238, t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t311, t313, t111, -t74, 0, 0, 0, 0, 0, 0, t311, t111, -t313, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t308, t312, t310, t30;];
tauJ_reg = t11;
