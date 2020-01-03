% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRP5
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:58
% EndTime: 2019-12-31 19:55:09
% DurationCPUTime: 5.98s
% Computational Cost: add. (12672->372), mult. (30102->480), div. (0->0), fcn. (21245->8), ass. (0->226)
t217 = sin(qJ(2));
t220 = cos(qJ(2));
t214 = sin(pkin(8));
t215 = cos(pkin(8));
t219 = cos(qJ(4));
t210 = qJDD(2) + qJDD(4);
t185 = (-t217 * t214 + t220 * t215) * qJD(1);
t187 = (t220 * t214 + t217 * t215) * qJD(1);
t216 = sin(qJ(4));
t161 = -t219 * t185 + t216 * t187;
t163 = t216 * t185 + t219 * t187;
t275 = t163 * t161;
t123 = -t275 - t210;
t266 = t216 * t123;
t160 = t163 ^ 2;
t211 = qJD(2) + qJD(4);
t293 = t211 ^ 2;
t299 = -t160 - t293;
t86 = t219 * t299 + t266;
t260 = t219 * t123;
t88 = -t216 * t299 + t260;
t41 = t214 * t86 - t215 * t88;
t52 = t214 * t88 + t215 * t86;
t342 = pkin(6) * (t217 * t52 + t220 * t41);
t341 = qJ(3) * t41;
t340 = qJ(3) * t52;
t339 = -pkin(2) * t52 - pkin(3) * t86;
t294 = t161 ^ 2;
t150 = t294 - t293;
t100 = -t219 * t150 - t266;
t96 = -t216 * t150 + t260;
t337 = t217 * (t215 * t100 - t214 * t96) + t220 * (t214 * t100 + t215 * t96);
t335 = pkin(7) * t86;
t334 = pkin(7) * t88;
t204 = t217 * qJDD(1);
t249 = qJD(1) * qJD(2);
t246 = t220 * t249;
t193 = t204 + t246;
t205 = t220 * qJDD(1);
t247 = t217 * t249;
t194 = t205 - t247;
t171 = -t214 * t193 + t215 * t194;
t172 = t215 * t193 + t214 * t194;
t234 = t216 * t171 + t219 * t172;
t106 = -t161 * qJD(4) + t234;
t273 = t211 * t161;
t301 = t106 - t273;
t156 = t211 * t163;
t241 = -t219 * t171 + t216 * t172;
t229 = t163 * qJD(4) + t241;
t77 = t156 + t229;
t46 = -t216 * t77 + t219 * t301;
t283 = t216 * t301;
t48 = t219 * t77 + t283;
t331 = t217 * (t214 * t46 + t215 * t48) + t220 * (t214 * t48 - t215 * t46);
t298 = -t275 + t210;
t265 = t216 * t298;
t297 = -t293 - t294;
t302 = t219 * t297 - t265;
t113 = t219 * t298;
t303 = t216 * t297 + t113;
t308 = t214 * t302 + t215 * t303;
t326 = qJ(3) * t308;
t309 = -t214 * t303 + t215 * t302;
t325 = qJ(3) * t309;
t324 = pkin(2) * t308 + pkin(3) * t303;
t323 = pkin(6) * (-t217 * t308 + t220 * t309);
t151 = -t160 + t293;
t310 = t219 * t151 + t265;
t311 = -t216 * t151 + t113;
t322 = t217 * (-t214 * t310 + t215 * t311) + t220 * (t214 * t311 + t215 * t310);
t102 = -t294 - t160;
t321 = pkin(1) * t102;
t320 = pkin(2) * t102;
t319 = pkin(3) * t102;
t317 = pkin(7) * t302;
t316 = pkin(7) * t303;
t274 = t187 * t185;
t304 = qJDD(2) + t274;
t315 = t214 * t304;
t314 = t215 * t304;
t305 = t301 * qJ(5);
t256 = qJD(2) * t185;
t145 = t172 - t256;
t300 = t106 + t273;
t125 = t160 - t294;
t222 = qJD(1) ^ 2;
t261 = t217 * t222;
t218 = sin(qJ(1));
t291 = cos(qJ(1));
t233 = t291 * g(1) + t218 * g(2);
t276 = qJDD(1) * pkin(6);
t190 = -t222 * pkin(1) - t233 + t276;
t263 = t217 * t190;
t137 = qJDD(2) * pkin(2) - t193 * qJ(3) - t263 + (pkin(2) * t261 + qJ(3) * t249 - g(3)) * t220;
t175 = -t217 * g(3) + t220 * t190;
t213 = t220 ^ 2;
t207 = t213 * t222;
t257 = qJD(1) * t217;
t231 = qJD(2) * pkin(2) - qJ(3) * t257;
t138 = -pkin(2) * t207 + t194 * qJ(3) - qJD(2) * t231 + t175;
t237 = -0.2e1 * qJD(3) * t187 + t215 * t137 - t214 * t138;
t92 = 0.2e1 * qJD(3) * t185 + t214 * t137 + t215 * t138;
t228 = (-t161 * t216 - t163 * t219) * t211;
t272 = t211 * t216;
t149 = t163 * t272;
t271 = t211 * t219;
t248 = t161 * t271;
t238 = t149 - t248;
t296 = t217 * (-t214 * t228 + t215 * t238) + t220 * (t214 * t238 + t215 * t228);
t230 = t216 * t229 + t248;
t239 = t161 * t272 - t219 * t229;
t295 = t217 * (-t214 * t239 + t215 * t230) + t220 * (t214 * t230 + t215 * t239);
t183 = t185 ^ 2;
t184 = t187 ^ 2;
t292 = 2 * qJD(5);
t290 = pkin(4) * t219;
t289 = t229 * pkin(4);
t124 = t161 * pkin(4) - t163 * qJ(5);
t225 = pkin(3) * t304 - t145 * pkin(7) + t237;
t235 = qJD(2) * pkin(3) - t187 * pkin(7);
t62 = -t183 * pkin(3) + t171 * pkin(7) - qJD(2) * t235 + t92;
t37 = t216 * t225 + t219 * t62;
t236 = t210 * qJ(5) - t161 * t124 + t211 * t292 + t37;
t30 = -pkin(4) * t293 + t236;
t36 = t216 * t62 - t219 * t225;
t32 = -t210 * pkin(4) - qJ(5) * t293 + t163 * t124 + qJDD(5) + t36;
t288 = -pkin(4) * t32 + qJ(5) * t30;
t78 = -t156 + t229;
t287 = -pkin(4) * t300 - qJ(5) * t78;
t18 = t216 * t37 - t219 * t36;
t286 = t214 * t18;
t285 = t215 * t18;
t282 = t216 * t300;
t244 = t218 * g(1) - t291 * g(2);
t232 = qJDD(1) * pkin(1) + t244;
t140 = t194 * pkin(2) - qJDD(3) - t231 * t257 + (qJ(3) * t213 + pkin(6)) * t222 + t232;
t90 = t171 * pkin(3) + t183 * pkin(7) - t187 * t235 + t140;
t281 = t216 * t90;
t57 = t214 * t92 + t215 * t237;
t280 = t217 * t57;
t278 = t219 * t90;
t277 = qJ(5) * t219;
t270 = t214 * t140;
t167 = qJDD(2) - t274;
t269 = t214 * t167;
t268 = t215 * t140;
t267 = t215 * t167;
t199 = t220 * t261;
t262 = t217 * (qJDD(2) + t199);
t259 = t220 * (qJDD(2) - t199);
t255 = qJD(2) * t187;
t254 = qJD(2) * t214;
t253 = qJD(2) * t215;
t250 = qJD(4) + t211;
t245 = -qJ(5) * t216 - pkin(3);
t58 = -t214 * t237 + t215 * t92;
t19 = t216 * t36 + t219 * t37;
t174 = t220 * g(3) + t263;
t240 = t217 * t174 + t220 * t175;
t143 = t171 + t255;
t227 = -pkin(4) * t299 - qJ(5) * t123 + t30;
t226 = pkin(4) * t298 + qJ(5) * t297 - t32;
t224 = -pkin(4) * t156 + t163 * t292 + t90;
t223 = t224 + t305;
t221 = qJD(2) ^ 2;
t212 = t217 ^ 2;
t206 = t212 * t222;
t195 = t205 - 0.2e1 * t247;
t192 = t204 + 0.2e1 * t246;
t189 = t222 * pkin(6) + t232;
t178 = -t184 - t221;
t177 = -t184 + t221;
t176 = t183 - t221;
t165 = -t221 - t183;
t144 = t172 + t256;
t142 = -t171 + t255;
t141 = -t183 - t184;
t133 = -t214 * t178 - t267;
t132 = t215 * t178 - t269;
t131 = t215 * t165 - t315;
t130 = t214 * t165 + t314;
t108 = t215 * t143 + t214 * t145;
t107 = t214 * t143 - t215 * t145;
t80 = -t250 * t161 + t234;
t79 = (-qJD(4) + t211) * t163 - t241;
t76 = t250 * t163 + t241;
t73 = t219 * t300;
t72 = t219 * t106 - t149;
t71 = t216 * t106 + t163 * t271;
t56 = -t278 - t335;
t51 = t219 * t79 + t282;
t49 = -t219 * t78 + t282;
t47 = t216 * t79 - t73;
t45 = -t216 * t78 - t73;
t43 = -t281 - t316;
t38 = -pkin(3) * t80 - t281 + t334;
t34 = -pkin(3) * t77 + t278 + t317;
t33 = t223 - t289;
t28 = -t214 * t47 + t215 * t51;
t27 = -t214 * t45 + t215 * t49;
t26 = t214 * t51 + t215 * t47;
t25 = t214 * t49 + t215 * t45;
t24 = t223 + (-t229 - t76) * pkin(4);
t23 = t224 - t289 + 0.2e1 * t305;
t22 = t217 * (-t214 * t71 + t215 * t72) + t220 * (t214 * t72 + t215 * t71);
t21 = -qJ(5) * t102 + t32;
t20 = (-t102 - t293) * pkin(4) + t236;
t17 = -t216 * t24 - t76 * t277 - t316;
t16 = -pkin(4) * t283 + t219 * t23 + t335;
t15 = pkin(3) * t90 + pkin(7) * t19;
t14 = t219 * t24 + t245 * t76 + t317;
t13 = -t334 + t216 * t23 + (pkin(3) + t290) * t301;
t12 = t216 * t32 + t219 * t30;
t11 = t216 * t30 - t219 * t32;
t10 = -pkin(7) * t47 - t18;
t9 = pkin(7) * t51 + t19 - t319;
t8 = -pkin(7) * t45 - t216 * t20 + t219 * t21;
t7 = t215 * t19 - t286;
t6 = t214 * t19 + t285;
t5 = pkin(7) * t49 + t219 * t20 + t216 * t21 - t319;
t4 = -pkin(7) * t11 + (-pkin(4) * t216 + t277) * t33;
t3 = -t214 * t11 + t215 * t12;
t2 = t215 * t11 + t214 * t12;
t1 = pkin(7) * t12 + (-t245 + t290) * t33;
t29 = [0, 0, 0, 0, 0, qJDD(1), t244, t233, 0, 0, (t193 + t246) * t217, t220 * t192 + t217 * t195, t262 + t220 * (-t206 + t221), (t194 - t247) * t220, t217 * (t207 - t221) + t259, 0, t220 * t189 + pkin(1) * t195 + pkin(6) * (t220 * (-t207 - t221) - t262), -t217 * t189 - pkin(1) * t192 + pkin(6) * (-t259 - t217 * (-t206 - t221)), pkin(1) * (t206 + t207) + (t212 + t213) * t276 + t240, pkin(1) * t189 + pkin(6) * t240, t217 * (t215 * t172 - t187 * t254) + t220 * (t214 * t172 + t187 * t253), t217 * (-t215 * t142 - t214 * t144) + t220 * (-t214 * t142 + t215 * t144), t217 * (-t214 * t177 + t314) + t220 * (t215 * t177 + t315), t217 * (-t214 * t171 - t185 * t253) + t220 * (t215 * t171 - t185 * t254), t217 * (t215 * t176 - t269) + t220 * (t214 * t176 + t267), (t217 * (t185 * t215 + t187 * t214) + t220 * (t185 * t214 - t187 * t215)) * qJD(2), t217 * (-qJ(3) * t130 - t270) + t220 * (-pkin(2) * t142 + qJ(3) * t131 + t268) - pkin(1) * t142 + pkin(6) * (-t217 * t130 + t220 * t131), t217 * (-qJ(3) * t132 - t268) + t220 * (-pkin(2) * t144 + qJ(3) * t133 - t270) - pkin(1) * t144 + pkin(6) * (-t217 * t132 + t220 * t133), t217 * (-qJ(3) * t107 - t57) + t220 * (-pkin(2) * t141 + qJ(3) * t108 + t58) - pkin(1) * t141 + pkin(6) * (-t217 * t107 + t220 * t108), -qJ(3) * t280 + t220 * (pkin(2) * t140 + qJ(3) * t58) + pkin(1) * t140 + pkin(6) * (t220 * t58 - t280), t22, -t331, t322, t295, -t337, t296, t217 * (-t214 * t34 + t215 * t43 - t326) + t220 * (-pkin(2) * t77 + t214 * t43 + t215 * t34 + t325) - pkin(1) * t77 + t323, t217 * (-t214 * t38 + t215 * t56 - t340) + t220 * (-pkin(2) * t80 + t214 * t56 + t215 * t38 - t341) - pkin(1) * t80 - t342, t217 * (-qJ(3) * t26 + t215 * t10 - t214 * t9) + t220 * (qJ(3) * t28 + t214 * t10 + t215 * t9 - t320) - t321 + pkin(6) * (-t217 * t26 + t220 * t28), t217 * (-pkin(7) * t285 - qJ(3) * t6 - t214 * t15) + t220 * (pkin(2) * t90 - pkin(7) * t286 + qJ(3) * t7 + t215 * t15) + pkin(1) * t90 + pkin(6) * (-t217 * t6 + t220 * t7), t22, t322, t331, t296, t337, t295, t217 * (-t214 * t14 + t215 * t17 - t326) + t220 * (-pkin(2) * t76 + t215 * t14 + t214 * t17 + t325) - pkin(1) * t76 + t323, t217 * (-qJ(3) * t25 - t214 * t5 + t215 * t8) + t220 * (qJ(3) * t27 + t214 * t8 + t215 * t5 - t320) - t321 + pkin(6) * (-t217 * t25 + t220 * t27), t217 * (-t214 * t13 + t215 * t16 + t340) + t220 * (pkin(2) * t301 + t215 * t13 + t214 * t16 + t341) + pkin(1) * t301 + t342, t217 * (-qJ(3) * t2 - t214 * t1 + t215 * t4) + t220 * (pkin(2) * t33 + qJ(3) * t3 + t215 * t1 + t214 * t4) + pkin(1) * t33 + pkin(6) * (-t217 * t2 + t220 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t199, -t207 + t206, t204, t199, t205, qJDD(2), -t174, -t175, 0, 0, -t274, t184 - t183, t145, t274, t143, qJDD(2), pkin(2) * t130 + t237, pkin(2) * t132 - t92, pkin(2) * t107, pkin(2) * t57, t275, t125, t300, -t275, -t78, t210, t324 - t36, -t339 - t37, pkin(2) * t26 + pkin(3) * t47, pkin(2) * t6 + pkin(3) * t18, t275, t300, -t125, t210, t78, -t275, t226 + t324, pkin(2) * t25 + pkin(3) * t45 + t287, t227 + t339, pkin(2) * t2 + pkin(3) * t11 + t288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, t144, t141, -t140, 0, 0, 0, 0, 0, 0, t77, t80, t102, -t90, 0, 0, 0, 0, 0, 0, t76, t102, -t301, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, t125, t300, -t275, -t78, t210, -t36, -t37, 0, 0, t275, t300, -t125, t210, t78, -t275, t226, t287, t227, t288; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t298, t300, t299, t32;];
tauJ_reg = t29;
