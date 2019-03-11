% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPPRRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:29:05
% EndTime: 2019-03-09 02:29:10
% DurationCPUTime: 2.23s
% Computational Cost: add. (1843->282), mult. (3859->385), div. (0->0), fcn. (2548->6), ass. (0->126)
t267 = qJD(4) + qJD(5);
t277 = cos(qJ(6));
t275 = sin(qJ(5));
t278 = cos(qJ(5));
t279 = cos(qJ(4));
t324 = qJD(1) * t279;
t276 = sin(qJ(4));
t325 = qJD(1) * t276;
t232 = t275 * t325 - t278 * t324;
t274 = sin(qJ(6));
t334 = t232 * t274;
t222 = -t277 * t267 - t334;
t237 = t275 * t279 + t276 * t278;
t233 = t237 * qJD(1);
t354 = qJD(6) + t233;
t357 = t222 * t354;
t289 = t232 * t277 - t267 * t274;
t356 = t289 * t354;
t320 = qJD(5) * t275;
t322 = qJD(4) * t276;
t355 = t275 * t322 + t276 * t320;
t333 = t233 * t267;
t353 = t274 * t354;
t265 = qJ(2) * qJD(1) + qJD(3);
t244 = -pkin(7) * qJD(1) + t265;
t321 = qJD(4) * t279;
t323 = qJD(2) * t276;
t221 = t244 * t321 + (-pkin(8) * t321 + t323) * qJD(1);
t230 = -pkin(8) * t324 + t279 * t244;
t226 = qJD(4) * pkin(4) + t230;
t352 = (qJD(5) * t226 + t221) * t278;
t351 = MDP(11) * (t276 ^ 2 - t279 ^ 2);
t273 = pkin(1) + qJ(3);
t350 = qJD(1) * t273;
t349 = MDP(6) * qJ(2) + MDP(5) + MDP(7);
t314 = qJD(1) * qJD(2);
t258 = t279 * t314;
t220 = t258 + (pkin(8) * qJD(1) - t244) * t322;
t229 = -pkin(8) * t325 + t244 * t276;
t304 = t220 * t275 - t229 * t320;
t186 = t304 + t352;
t336 = t229 * t278;
t206 = t226 * t275 + t336;
t305 = -t278 * t220 + t221 * t275;
t187 = qJD(5) * t206 + t305;
t272 = -pkin(7) + qJ(2);
t347 = pkin(8) - t272;
t240 = t347 * t276;
t241 = t347 * t279;
t218 = -t240 * t275 + t241 * t278;
t227 = qJD(2) * t279 + t322 * t347;
t228 = -qJD(4) * t241 + t323;
t192 = -qJD(5) * t218 + t227 * t275 + t228 * t278;
t337 = t229 * t275;
t205 = t226 * t278 - t337;
t201 = -pkin(5) * t267 - t205;
t245 = -qJD(2) + t350;
t235 = pkin(4) * t325 + t245;
t204 = pkin(5) * t233 + pkin(9) * t232 + t235;
t330 = t278 * t279;
t293 = t267 * t330;
t328 = t355 * qJD(1);
t212 = qJD(1) * t293 - t328;
t236 = t275 * t276 - t330;
t255 = t276 * pkin(4) + t273;
t213 = pkin(5) * t237 + pkin(9) * t236 + t255;
t319 = qJD(5) * t278;
t216 = -t275 * t321 - t276 * t319 - t278 * t322 - t279 * t320;
t219 = -t240 * t278 - t241 * t275;
t348 = -(qJD(6) * t213 + t192) * t354 - t219 * t212 - (qJD(6) * t204 + t186) * t237 - t187 * t236 + t201 * t216;
t318 = qJD(6) * t274;
t317 = qJD(6) * t277;
t329 = t267 * t317 - t277 * t333;
t195 = t232 * t318 + t329;
t345 = t195 * t236;
t344 = t195 * t274;
t343 = t201 * t233;
t342 = t201 * t236;
t341 = t333 * t274;
t340 = t213 * t212;
t339 = t222 * t232;
t338 = t289 * t232;
t335 = t232 * t267;
t332 = t274 * t212;
t331 = t277 * t212;
t268 = qJD(3) * qJD(1);
t313 = qJD(1) * qJD(4);
t307 = t279 * t313;
t239 = pkin(4) * t307 + t268;
t316 = qJD(2) - t245;
t247 = pkin(4) * t321 + qJD(3);
t312 = pkin(4) * t324;
t310 = 0.2e1 * t268;
t306 = -pkin(4) * t267 - t226;
t303 = t277 * t354;
t214 = -pkin(5) * t232 + pkin(9) * t233;
t260 = pkin(4) * t275 + pkin(9);
t299 = qJD(6) * t260 + t214 + t312;
t298 = t278 * t267;
t297 = MDP(23) * t267;
t207 = t230 * t275 + t336;
t296 = pkin(4) * t320 - t207;
t208 = t230 * t278 - t337;
t295 = -pkin(4) * t319 + t208;
t202 = pkin(9) * t267 + t206;
t189 = t202 * t277 + t204 * t274;
t294 = t187 * t274 - t189 * t232 + t201 * t317;
t292 = qJD(2) + t245 + t350;
t291 = -t212 * t260 + t343;
t290 = t202 * t274 - t204 * t277;
t280 = qJD(4) ^ 2;
t288 = -t272 * t280 + t310;
t287 = -t187 * t277 + t201 * t318 - t232 * t290;
t286 = t232 * t235 - t305;
t285 = t233 * t235 - t304;
t284 = t216 * t277 + t236 * t318;
t196 = -qJD(6) * t289 - t341;
t282 = ((t195 - t357) * t277 + (-t196 + t356) * t274) * MDP(25) + (-t289 * t303 + t344) * MDP(24) + (-t353 * t354 + t331 - t339) * MDP(27) + (t303 * t354 + t332 - t338) * MDP(26) + (-t298 * t324 + t328 - t335) * MDP(20) + (t232 ^ 2 - t233 ^ 2) * MDP(18) + (-MDP(17) * t233 + MDP(28) * t354) * t232;
t281 = qJD(1) ^ 2;
t261 = -pkin(4) * t278 - pkin(5);
t217 = t293 - t355;
t194 = pkin(5) * t217 - pkin(9) * t216 + t247;
t193 = qJD(5) * t219 - t227 * t278 + t228 * t275;
t191 = pkin(5) * t212 + pkin(9) * t333 + t239;
t190 = t277 * t191;
t1 = [(qJD(2) * t265 + qJD(3) * t245 + (qJ(2) * qJD(2) + qJD(3) * t273) * qJD(1)) * MDP(9) + 0.2e1 * t313 * t351 + t292 * t321 * MDP(15) + (t279 * t288 - t292 * t322) * MDP(16) + (-t216 * t232 + t236 * t333) * MDP(17) + (t212 * t236 - t216 * t233 + t217 * t232 + t237 * t333) * MDP(18) + (t212 * t255 + t217 * t235 + t233 * t247 + t237 * t239) * MDP(22) + (t216 * t235 - t232 * t247 - t236 * t239 - t255 * t333) * MDP(23) + (-t277 * t345 - t284 * t289) * MDP(24) + ((-t222 * t277 + t274 * t289) * t216 + (t344 + t196 * t277 + (-t222 * t274 - t277 * t289) * qJD(6)) * t236) * MDP(25) + (t195 * t237 - t217 * t289 - t236 * t331 + t284 * t354) * MDP(26) + (t236 * t332 - t196 * t237 - t217 * t222 + (-t216 * t274 + t236 * t317) * t354) * MDP(27) + (t212 * t237 + t217 * t354) * MDP(28) + (-t290 * t217 + t190 * t237 + t193 * t222 + t218 * t196 + (t194 * t354 + t340 + (-t202 * t237 - t219 * t354 - t342) * qJD(6)) * t277 + t348 * t274) * MDP(29) + (-t189 * t217 - t193 * t289 + t218 * t195 + (-(-qJD(6) * t219 + t194) * t354 - t340 - (-qJD(6) * t202 + t191) * t237 + qJD(6) * t342) * t274 + t348 * t277) * MDP(30) + MDP(8) * t310 - t280 * t279 * MDP(13) + (-0.2e1 * MDP(10) * t307 - t280 * MDP(12) + MDP(15) * t288) * t276 + (MDP(19) * t216 - MDP(20) * t217 - MDP(22) * t193 - MDP(23) * t192) * t267 + 0.2e1 * t349 * t314; -t268 * MDP(9) + (t328 + t335) * MDP(22) + MDP(23) * t333 + (-t331 - t339) * MDP(29) + (t332 + t338) * MDP(30) + (MDP(29) * t353 + MDP(30) * t303) * t354 - t349 * t281 + (-t265 * MDP(9) + (0.2e1 * qJD(4) * MDP(16) + t278 * t297) * t276 + (-0.2e1 * qJD(4) * MDP(15) - MDP(22) * t298 + t275 * t297) * t279) * qJD(1); -t281 * MDP(8) + (t196 * t236 - t216 * t222 - t237 * t332) * MDP(29) + (t216 * t289 - t237 * t331 + t345) * MDP(30) + (MDP(22) * t216 - MDP(23) * t217) * t267 + (-t233 * MDP(22) + t232 * MDP(23) + t316 * MDP(9)) * qJD(1) + ((-qJD(1) * t277 - t217 * t274 - t237 * t317) * MDP(29) + (qJD(1) * t274 - t217 * t277 + t237 * t318) * MDP(30)) * t354 + (MDP(15) * t276 + MDP(16) * t279) * (-t280 - t281); (-t245 * t324 + t258) * MDP(15) + (-t233 * t312 + t207 * t267 + (t275 * t306 - t336) * qJD(5) + t286) * MDP(22) + (t232 * t312 + t208 * t267 + (qJD(5) * t306 - t221) * t278 + t285) * MDP(23) + (t261 * t196 + t291 * t274 + t296 * t222 + (t274 * t295 - t277 * t299) * t354 + t287) * MDP(29) + t282 + (t261 * t195 + t291 * t277 - t296 * t289 + (t274 * t299 + t277 * t295) * t354 + t294) * MDP(30) - t316 * MDP(16) * t325 + (t279 * t276 * MDP(10) - t351) * t281; (t286 + (-qJD(5) + t267) * t206) * MDP(22) + (t205 * t267 + t285 - t352) * MDP(23) + (-pkin(5) * t196 - (-t205 * t274 + t214 * t277) * t354 - t206 * t222 + t274 * t343 + (-t317 * t354 - t332) * pkin(9) + t287) * MDP(29) + (-pkin(5) * t195 + (t205 * t277 + t214 * t274) * t354 + t206 * t289 + t277 * t343 + (t318 * t354 - t331) * pkin(9) + t294) * MDP(30) + t282; -t289 * t222 * MDP(24) + (-t222 ^ 2 + t289 ^ 2) * MDP(25) + (t329 + t357) * MDP(26) + (t341 - t356) * MDP(27) + t212 * MDP(28) + (-t186 * t274 + t189 * t354 + t201 * t289 + t190) * MDP(29) + (-t186 * t277 - t191 * t274 + t201 * t222 - t290 * t354) * MDP(30) + (MDP(26) * t334 + MDP(27) * t289 - MDP(29) * t189 + MDP(30) * t290) * qJD(6);];
tauc  = t1;
