% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:55:00
% EndTime: 2019-12-31 19:55:05
% DurationCPUTime: 2.25s
% Computational Cost: add. (2503->246), mult. (6520->316), div. (0->0), fcn. (4676->6), ass. (0->132)
t307 = sin(pkin(8));
t308 = cos(pkin(8));
t311 = cos(qJ(2));
t341 = qJD(1) * qJD(2);
t336 = t311 * t341;
t310 = sin(qJ(2));
t337 = t310 * t341;
t271 = -t307 * t337 + t308 * t336;
t330 = t307 * t310 - t308 * t311;
t280 = t330 * qJD(1);
t288 = t307 * t311 + t308 * t310;
t281 = t288 * qJD(1);
t309 = sin(qJ(4));
t319 = qJD(2) * t281;
t365 = cos(qJ(4));
t338 = qJD(4) * t365;
t343 = qJD(4) * t309;
t324 = t365 * t271 - t280 * t338 - t281 * t343 - t309 * t319;
t236 = -t365 * t280 - t281 * t309;
t304 = qJD(2) + qJD(4);
t353 = t236 * t304;
t194 = t324 - t353;
t323 = qJD(2) * t288;
t317 = t365 * t323;
t327 = t309 * t280 - t365 * t281;
t315 = -qJD(1) * t317 + t327 * qJD(4) - t309 * t271;
t354 = t236 ^ 2;
t355 = t327 * t304;
t376 = t327 ^ 2;
t378 = t236 * t327;
t379 = t194 * MDP(15) + MDP(13) * t378 + (t315 - t355) * MDP(16) + (-t354 + t376) * MDP(14);
t339 = -t311 * pkin(2) - pkin(1);
t331 = t339 * qJD(1);
t295 = qJD(3) + t331;
t245 = pkin(3) * t280 + t295;
t199 = -pkin(4) * t236 + qJ(5) * t327 + t245;
t375 = t199 * t236;
t374 = t199 * t327;
t372 = t245 * t236;
t371 = t245 * t327;
t208 = -pkin(4) * t327 - qJ(5) * t236;
t370 = -0.2e1 * t341;
t369 = MDP(4) * t310;
t368 = (t310 ^ 2 - t311 ^ 2) * MDP(5);
t361 = -qJ(3) - pkin(6);
t296 = t361 * t310;
t292 = qJD(1) * t296;
t297 = t361 * t311;
t293 = qJD(1) * t297;
t352 = t308 * t293;
t243 = -t292 * t307 + t352;
t363 = pkin(7) * t280;
t227 = t243 + t363;
t283 = t307 * t293;
t244 = t308 * t292 + t283;
t362 = pkin(7) * t281;
t228 = t244 - t362;
t300 = pkin(2) * t308 + pkin(3);
t364 = pkin(2) * t307;
t340 = t309 * t364;
t367 = qJD(4) * t340 + t309 * t227 + t365 * t228 - t300 * t338;
t360 = qJD(2) * pkin(2);
t335 = qJD(2) * t361;
t277 = qJD(3) * t311 + t310 * t335;
t278 = -qJD(3) * t310 + t311 * t335;
t229 = -t307 * t277 + t308 * t278;
t322 = t330 * qJD(2);
t219 = pkin(7) * t322 + t229;
t230 = t308 * t277 + t307 * t278;
t220 = -pkin(7) * t323 + t230;
t247 = t308 * t296 + t297 * t307;
t231 = -pkin(7) * t288 + t247;
t248 = t307 * t296 - t308 * t297;
t232 = -t330 * pkin(7) + t248;
t328 = t365 * t231 - t309 * t232;
t191 = t328 * qJD(4) + t309 * t219 + t365 * t220;
t359 = t191 * t304;
t205 = t309 * t231 + t365 * t232;
t192 = t205 * qJD(4) - t365 * t219 + t309 * t220;
t358 = t192 * t304;
t287 = t292 + t360;
t239 = t308 * t287 + t283;
t223 = qJD(2) * pkin(3) + t239 - t362;
t240 = t307 * t287 - t352;
t224 = t240 - t363;
t198 = t309 * t223 + t365 * t224;
t357 = t198 * t304;
t312 = qJD(2) ^ 2;
t351 = t310 * t312;
t350 = t311 * t312;
t313 = qJD(1) ^ 2;
t349 = t311 * t313;
t321 = t309 * t300 + t365 * t364;
t348 = -t321 * qJD(4) - t365 * t227 + t309 * t228;
t347 = -qJD(5) + t367;
t258 = t277 * qJD(1);
t259 = t278 * qJD(1);
t226 = t308 * t258 + t307 * t259;
t344 = qJD(1) * t310;
t197 = t365 * t223 - t309 * t224;
t342 = qJD(5) - t197;
t302 = t310 * t360;
t252 = pkin(2) * t344 + pkin(3) * t281;
t334 = pkin(1) * t370;
t225 = -t258 * t307 + t308 * t259;
t215 = -pkin(7) * t271 + t225;
t216 = -pkin(7) * t319 + t226;
t332 = t309 * t215 + t365 * t216 + t223 * t338 - t224 * t343;
t189 = -t365 * t215 + t309 * t216 + t223 * t343 + t224 * t338;
t303 = t304 * qJD(5);
t188 = t303 + t332;
t326 = t197 * t304 - t332;
t325 = t189 - t374;
t320 = t365 * t330;
t318 = t348 * t304 - t189;
t253 = pkin(3) * t323 + t302;
t261 = t330 * pkin(3) + t339;
t299 = pkin(2) * t337;
t246 = pkin(3) * t319 + t299;
t242 = t365 * t288 - t309 * t330;
t190 = -pkin(4) * t315 - qJ(5) * t324 + qJD(5) * t327 + t246;
t275 = -t365 * t300 - pkin(4) + t340;
t274 = qJ(5) + t321;
t241 = t288 * t309 + t320;
t210 = t242 * qJD(4) - t309 * t322 + t317;
t209 = t288 * t343 + t304 * t320 + t309 * t323;
t203 = t241 * pkin(4) - t242 * qJ(5) + t261;
t202 = t208 + t252;
t196 = t304 * qJ(5) + t198;
t195 = -t304 * pkin(4) + t342;
t193 = t210 * pkin(4) + t209 * qJ(5) - t242 * qJD(5) + t253;
t1 = [0.2e1 * t336 * t369 + t368 * t370 + MDP(6) * t350 - MDP(7) * t351 + (-pkin(6) * t350 + t310 * t334) * MDP(9) + (pkin(6) * t351 + t311 * t334) * MDP(10) + (-t230 * t280 - t226 * t330 - t229 * t281 - t247 * t271 - t225 * t288 + (t239 * t330 - t240 * t288 - t248 * t281) * qJD(2)) * MDP(11) + (t225 * t247 + t226 * t248 + t239 * t229 + t240 * t230 + (t295 + t331) * t302) * MDP(12) + (t209 * t327 + t242 * t324) * MDP(13) + (-t209 * t236 + t210 * t327 - t241 * t324 + t242 * t315) * MDP(14) + (t210 * t245 - t236 * t253 + t241 * t246 - t261 * t315 - t358) * MDP(18) + (-t209 * t245 + t242 * t246 - t253 * t327 + t261 * t324 - t359) * MDP(19) + (t190 * t241 - t193 * t236 + t199 * t210 - t203 * t315 - t358) * MDP(20) + (-t188 * t241 + t189 * t242 + t191 * t236 - t192 * t327 - t195 * t209 - t196 * t210 + t205 * t315 - t324 * t328) * MDP(21) + (-t190 * t242 + t193 * t327 + t199 * t209 - t203 * t324 + t359) * MDP(22) + (t188 * t205 - t189 * t328 + t190 * t203 + t191 * t196 + t192 * t195 + t193 * t199) * MDP(23) + (-t209 * MDP(15) - t210 * MDP(16)) * t304; -t349 * t369 + t313 * t368 + ((t240 + t243) * t281 - (-t244 + t239) * t280 + (-t308 * t271 - t307 * t319) * pkin(2)) * MDP(11) + (-t239 * t243 - t240 * t244 + (t225 * t308 + t226 * t307 - t295 * t344) * pkin(2)) * MDP(12) + (t236 * t252 + t318 + t371) * MDP(18) + (t252 * t327 + t367 * t304 - t332 - t372) * MDP(19) + (t202 * t236 + t318 + t374) * MDP(20) + (t315 * t274 + t275 * t324 + (-t196 + t348) * t327 + (-t195 - t347) * t236) * MDP(21) + (-t202 * t327 - t347 * t304 + t188 + t375) * MDP(22) + (t188 * t274 + t189 * t275 - t348 * t195 - t347 * t196 - t199 * t202) * MDP(23) + (t313 * t310 * MDP(9) + MDP(10) * t349) * pkin(1) + t379; (-t280 ^ 2 - t281 ^ 2) * MDP(11) + (t239 * t281 + t240 * t280 + t299) * MDP(12) + (-t354 - t376) * MDP(21) + (t195 * t327 - t196 * t236 + t190) * MDP(23) + (MDP(19) - MDP(22)) * (t324 + t353) + (MDP(18) + MDP(20)) * (-t315 - t355); (-t189 + t357 + t371) * MDP(18) + (t326 - t372) * MDP(19) + (t208 * t236 - t325 + t357) * MDP(20) + (-pkin(4) * t324 + qJ(5) * t315 - (t196 - t198) * t327 - (t195 - t342) * t236) * MDP(21) + (-t208 * t327 + 0.2e1 * t303 - t326 + t375) * MDP(22) + (-pkin(4) * t189 + qJ(5) * t188 - t195 * t198 + t342 * t196 - t199 * t208) * MDP(23) + t379; MDP(20) * t378 + t194 * MDP(21) + (-t304 ^ 2 - t376) * MDP(22) + (-t196 * t304 + t325) * MDP(23);];
tauc = t1;
