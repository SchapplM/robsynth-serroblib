% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:13:58
% EndTime: 2019-03-09 02:14:04
% DurationCPUTime: 2.70s
% Computational Cost: add. (2756->303), mult. (6065->406), div. (0->0), fcn. (4203->6), ass. (0->129)
t280 = sin(pkin(9));
t281 = cos(pkin(9));
t284 = sin(qJ(4));
t286 = cos(qJ(4));
t256 = t280 * t286 + t281 * t284;
t292 = qJD(1) * t256;
t356 = qJD(5) + t292;
t363 = pkin(5) * t356;
t285 = cos(qJ(5));
t283 = sin(qJ(5));
t321 = qJD(4) * t283;
t324 = qJD(1) * t284;
t308 = t280 * t324;
t323 = qJD(1) * t286;
t312 = t281 * t323;
t359 = t308 - t312;
t234 = -t285 * t359 + t321;
t362 = t234 * t356;
t317 = qJD(5) * t285;
t204 = -t359 * t317 + (qJD(5) - t292) * t321;
t293 = t256 * qJD(4);
t290 = qJD(1) * t293;
t316 = t285 * qJD(4);
t318 = qJD(5) * t283;
t203 = -qJD(5) * t316 + t285 * t290 - t318 * t359;
t263 = qJD(4) * t308;
t319 = qJD(4) * t286;
t311 = t281 * t319;
t243 = qJD(1) * t311 - t263;
t282 = -pkin(1) - qJ(3);
t351 = t282 * qJD(1);
t264 = qJD(2) + t351;
t306 = -pkin(7) * qJD(1) + t264;
t244 = t306 * t280;
t245 = t306 * t281;
t215 = t286 * t244 + t284 * t245;
t211 = qJD(4) * pkin(8) + t215;
t279 = qJD(1) * qJ(2);
t273 = qJD(3) + t279;
t275 = t280 * pkin(3);
t259 = qJD(1) * t275 + t273;
t212 = pkin(4) * t292 + pkin(8) * t359 + t259;
t193 = t211 * t285 + t212 * t283;
t294 = t256 * qJD(3);
t352 = -t284 * t244 + t245 * t286;
t197 = -qJD(1) * t294 + qJD(4) * t352;
t278 = qJD(1) * qJD(2);
t213 = t243 * pkin(4) + pkin(8) * t290 + t278;
t208 = t285 * t213;
t289 = -t193 * qJD(5) - t197 * t283 + t208;
t181 = pkin(5) * t243 + qJ(6) * t203 - qJD(6) * t234 + t289;
t232 = -t283 * t359 - t316;
t187 = -qJ(6) * t232 + t193;
t361 = t187 * t356 + t181;
t295 = t285 * t197 - t211 * t318 + t212 * t317 + t283 * t213;
t182 = -qJ(6) * t204 - qJD(6) * t232 + t295;
t192 = -t211 * t283 + t285 * t212;
t186 = -qJ(6) * t234 + t192;
t185 = t186 + t363;
t360 = -t185 * t356 + t182;
t325 = t280 ^ 2 + t281 ^ 2;
t358 = t325 * qJD(3);
t355 = MDP(6) * qJ(2) + MDP(7) * t280 + MDP(8) * t281 + MDP(5);
t349 = -pkin(7) + t282;
t257 = t349 * t280;
t258 = t349 * t281;
t228 = t284 * t257 - t258 * t286;
t350 = t234 ^ 2;
t348 = -qJ(6) - pkin(8);
t344 = t203 * t283;
t343 = t204 * t285;
t342 = t232 * t292;
t341 = t232 * t283;
t340 = t232 * t285;
t339 = t234 * t283;
t338 = t234 * t285;
t336 = t292 * t283;
t255 = t280 * t284 - t286 * t281;
t335 = t255 * t283;
t334 = t255 * t285;
t229 = t257 * t286 + t258 * t284;
t227 = t285 * t229;
t239 = t285 * t243;
t269 = qJ(2) + t275;
t332 = t185 - t186;
t331 = -t283 * t204 - t232 * t317;
t226 = -pkin(4) * t359 + pkin(8) * t292;
t330 = t283 * t226 + t285 * t352;
t329 = -t336 * t356 + t239;
t225 = pkin(4) * t256 + pkin(8) * t255 + t269;
t328 = t283 * t225 + t227;
t307 = qJD(5) * t348;
t327 = -qJ(6) * t336 + qJD(6) * t285 + t283 * t307 - t330;
t222 = t285 * t226;
t326 = pkin(5) * t359 - t222 + (-qJ(6) * t292 + t307) * t285 + (-qJD(6) + t352) * t283;
t320 = qJD(4) * t284;
t252 = -t280 * t319 - t281 * t320;
t322 = qJD(4) * t252;
t205 = -t228 * qJD(4) - t294;
t253 = -t280 * t320 + t311;
t223 = pkin(4) * t253 - pkin(8) * t252 + qJD(2);
t313 = t285 * t205 + t283 * t223 + t225 * t317;
t310 = t255 * t317;
t305 = qJD(1) * t325;
t304 = t285 * t356;
t303 = -t185 * t285 - t187 * t283;
t302 = -t185 * t283 + t187 * t285;
t300 = -qJ(6) * t252 + qJD(6) * t255;
t299 = -MDP(23) * t283 - MDP(24) * t285;
t210 = -qJD(4) * pkin(4) - t352;
t298 = -t252 * t283 + t310;
t297 = t252 * t285 + t255 * t318;
t296 = -pkin(8) * t243 + t210 * t356;
t291 = qJD(3) * t359 - t244 * t319 - t245 * t320;
t189 = pkin(5) * t204 - t291;
t206 = -t255 * qJD(3) + t229 * qJD(4);
t288 = (t338 + t341) * MDP(25) + t303 * MDP(26) + (-MDP(23) * t285 + MDP(24) * t283) * t356;
t287 = qJD(1) ^ 2;
t261 = t348 * t285;
t260 = t348 * t283;
t231 = t232 ^ 2;
t220 = t285 * t225;
t218 = t285 * t223;
t195 = pkin(5) * t232 + qJD(6) + t210;
t194 = qJ(6) * t335 + t328;
t191 = pkin(5) * t256 + qJ(6) * t334 - t229 * t283 + t220;
t184 = qJ(6) * t310 + (-qJD(5) * t229 + t300) * t283 + t313;
t183 = pkin(5) * t253 - t205 * t283 + t218 + t300 * t285 + (-t227 + (-qJ(6) * t255 - t225) * t283) * qJD(5);
t1 = [0.2e1 * qJD(3) * MDP(9) * t305 + ((t273 + t279) * qJD(2) + (-t264 - t351) * t358) * MDP(10) + (-t252 * t359 + t255 * t290) * MDP(11) + (t255 * t243 - t252 * t292 + t253 * t359 + t256 * t290) * MDP(12) + MDP(13) * t322 - t253 * qJD(4) * MDP(14) + (0.2e1 * t292 * qJD(2) - qJD(4) * t206 + t243 * t269 + t253 * t259) * MDP(16) + (-qJD(2) * t359 - t205 * qJD(4) + t259 * t252 + (-qJD(2) * t255 - t269 * t293) * qJD(1)) * MDP(17) + (t203 * t334 + t297 * t234) * MDP(18) + ((-t339 - t340) * t252 + (-t344 + t343 + (t338 - t341) * qJD(5)) * t255) * MDP(19) + (-t203 * t256 + t234 * t253 - t255 * t239 + t297 * t356) * MDP(20) + (-t204 * t256 - t232 * t253 + t243 * t335 + t298 * t356) * MDP(21) + (t243 * t256 + t253 * t356) * MDP(22) + ((-t229 * t317 + t218) * t356 + t220 * t243 + (-t211 * t317 + t208) * t256 + t192 * t253 + t206 * t232 + t228 * t204 - t210 * t310 + ((-qJD(5) * t225 - t205) * t356 - t229 * t243 + (-qJD(5) * t212 - t197) * t256 + t291 * t255 + t210 * t252) * t283) * MDP(23) + (-(-t229 * t318 + t313) * t356 - t328 * t243 - t295 * t256 - t193 * t253 + t206 * t234 - t228 * t203 + t291 * t334 + t297 * t210) * MDP(24) + (-t183 * t234 - t184 * t232 + t191 * t203 - t194 * t204 + t303 * t252 + (t302 * qJD(5) + t181 * t285 + t182 * t283) * t255) * MDP(25) + (t182 * t194 + t187 * t184 + t181 * t191 + t185 * t183 + t189 * (-pkin(5) * t335 + t228) + t195 * (-t298 * pkin(5) + t206)) * MDP(26) + 0.2e1 * t355 * t278; MDP(16) * t322 + (t204 * t255 - t232 * t252) * MDP(23) + (-t203 * t255 - t234 * t252) * MDP(24) + (t189 * t255 - t195 * t252) * MDP(26) - t355 * t287 + (-qJD(4) * MDP(17) + (t339 - t340) * MDP(25) + t302 * MDP(26) + t299 * t356) * t253 + ((-t273 - t358) * MDP(10) - t292 * MDP(16) + t359 * MDP(17) + t288) * qJD(1) + ((-t343 - t344) * MDP(25) + (-t181 * t283 + t182 * t285) * MDP(26) + t299 * t243 + t288 * qJD(5)) * t256; (t264 * t305 + t278) * MDP(10) - t263 * MDP(16) + t329 * MDP(23) + t331 * MDP(25) - t325 * MDP(9) * t287 - (-MDP(23) * t232 - MDP(24) * t234 - MDP(26) * t195) * t359 + (-MDP(23) * qJD(5) * t356 - t243 * MDP(24) + MDP(25) * t362 + MDP(26) * t360) * t283 + ((t203 - t342) * MDP(25) + t361 * MDP(26) - t356 ^ 2 * MDP(24)) * t285 + ((-t359 + t312) * MDP(16) + (-t280 * t323 - t281 * t324 - t292) * MDP(17)) * qJD(4); t359 ^ 2 * MDP(12) + (t263 + (-t359 - t312) * qJD(4)) * MDP(14) + (qJD(4) * t215 + t259 * t359 + t291) * MDP(16) + (t234 * t304 - t344) * MDP(18) + ((-t203 - t342) * t285 - t356 * t339 + t331) * MDP(19) + (t234 * t359 + t243 * t283 + t304 * t356) * MDP(20) + (-t232 * t359 - t318 * t356 + t329) * MDP(21) + t356 * t359 * MDP(22) + (-pkin(4) * t204 + t192 * t359 + t291 * t285 - t215 * t232 + (-pkin(8) * t317 - t222) * t356 + (t352 * t356 + t296) * t283) * MDP(23) + (pkin(4) * t203 - t193 * t359 - t291 * t283 - t215 * t234 + (pkin(8) * t318 + t330) * t356 + t296 * t285) * MDP(24) + (t203 * t260 + t204 * t261 - t327 * t232 - t326 * t234 - t283 * t361 + t360 * t285) * MDP(25) + (-t182 * t261 + t181 * t260 + t189 * (-pkin(5) * t285 - pkin(4)) + (t283 * t363 - t215) * t195 + t327 * t187 + t326 * t185) * MDP(26) + (-t359 * MDP(11) + (qJD(3) + t259) * MDP(17) - MDP(12) * t292) * t292; t234 * t232 * MDP(18) + (-t231 + t350) * MDP(19) + (t232 * t356 - t203) * MDP(20) + (-t204 + t362) * MDP(21) + t243 * MDP(22) + (t193 * t356 - t210 * t234 + t289) * MDP(23) + (t192 * t356 + t210 * t232 - t295) * MDP(24) + (pkin(5) * t203 - t332 * t232) * MDP(25) + (t332 * t187 + (-t195 * t234 + t181) * pkin(5)) * MDP(26); (-t231 - t350) * MDP(25) + (t185 * t234 + t187 * t232 + t189) * MDP(26);];
tauc  = t1;
