% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S6RPRPRP1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:03:07
% EndTime: 2019-03-09 03:03:14
% DurationCPUTime: 2.43s
% Computational Cost: add. (3048->305), mult. (7189->421), div. (0->0), fcn. (4902->8), ass. (0->138)
t290 = sin(pkin(9)) * pkin(1) + pkin(7);
t354 = qJ(4) + t290;
t298 = sin(pkin(10));
t300 = cos(pkin(10));
t303 = sin(qJ(3));
t305 = cos(qJ(3));
t282 = t298 * t305 + t300 * t303;
t274 = t282 * qJD(1);
t304 = cos(qJ(5));
t336 = qJD(1) * qJD(3);
t329 = t305 * t336;
t330 = t303 * t336;
t316 = -t298 * t330 + t300 * t329;
t337 = t304 * qJD(3);
t302 = sin(qJ(5));
t340 = qJD(5) * t302;
t227 = -qJD(5) * t337 + t274 * t340 - t304 * t316;
t258 = qJD(3) * t302 + t274 * t304;
t273 = t282 * qJD(3);
t268 = qJD(1) * t273;
t324 = t354 * qJD(1);
t260 = t305 * qJD(2) - t324 * t303;
t370 = qJD(3) * pkin(3);
t253 = t260 + t370;
t261 = t303 * qJD(2) + t324 * t305;
t358 = t300 * t261;
t215 = t298 * t253 + t358;
t212 = qJD(3) * pkin(8) + t215;
t292 = -cos(pkin(9)) * pkin(1) - pkin(2);
t318 = -pkin(3) * t305 + t292;
t311 = t318 * qJD(1);
t271 = qJD(4) + t311;
t341 = qJD(1) * t305;
t342 = qJD(1) * t303;
t272 = -t298 * t342 + t300 * t341;
t224 = -pkin(4) * t272 - pkin(8) * t274 + t271;
t202 = t212 * t304 + t224 * t302;
t335 = qJD(1) * qJD(4);
t247 = t260 * qJD(3) + t305 * t335;
t372 = -t261 * qJD(3) - t303 * t335;
t208 = t300 * t247 + t372 * t298;
t287 = pkin(3) * t330;
t229 = t268 * pkin(4) - t316 * pkin(8) + t287;
t226 = t304 * t229;
t308 = -t202 * qJD(5) - t208 * t302 + t226;
t191 = pkin(5) * t268 + qJ(6) * t227 - qJD(6) * t258 + t308;
t256 = t274 * t302 - t337;
t197 = -qJ(6) * t256 + t202;
t270 = qJD(5) - t272;
t378 = t270 * t197 + t191;
t228 = t258 * qJD(5) + t302 * t316;
t339 = qJD(5) * t304;
t313 = t304 * t208 - t212 * t340 + t224 * t339 + t302 * t229;
t192 = -qJ(6) * t228 - qJD(6) * t256 + t313;
t201 = -t212 * t302 + t304 * t224;
t196 = -qJ(6) * t258 + t201;
t195 = pkin(5) * t270 + t196;
t377 = -t270 * t195 + t192;
t376 = MDP(5) * t303;
t373 = (t303 ^ 2 - t305 ^ 2) * MDP(6);
t371 = t258 ^ 2;
t369 = t227 * t302;
t368 = t256 * t272;
t367 = t256 * t302;
t366 = t258 * t270;
t365 = t258 * t304;
t364 = t270 * t302;
t363 = t272 * t302;
t281 = t298 * t303 - t300 * t305;
t276 = t281 * qJD(3);
t362 = t276 * t302;
t361 = t276 * t304;
t360 = t282 * t302;
t359 = t282 * t304;
t250 = t298 * t261;
t357 = t302 * t268;
t306 = qJD(3) ^ 2;
t356 = t303 * t306;
t278 = t354 * t303;
t280 = t354 * t305;
t245 = -t278 * t298 + t280 * t300;
t241 = t304 * t245;
t266 = t304 * t268;
t355 = t305 * t306;
t289 = pkin(3) * t298 + pkin(8);
t353 = qJ(6) + t289;
t352 = t195 - t196;
t351 = -t228 * t359 + t256 * t361;
t217 = t260 * t300 - t250;
t238 = pkin(3) * t342 + pkin(4) * t274 - pkin(8) * t272;
t350 = t304 * t217 + t302 * t238;
t349 = -t302 * t228 - t256 * t339;
t348 = -t227 * t281 + t258 * t273;
t242 = pkin(4) * t281 - pkin(8) * t282 + t318;
t347 = t302 * t242 + t241;
t346 = t270 * t363 + t266;
t326 = qJD(5) * t353;
t345 = qJ(6) * t363 + qJD(6) * t304 - t302 * t326 - t350;
t232 = t304 * t238;
t344 = -pkin(5) * t274 - t232 + (qJ(6) * t272 - t326) * t304 + (-qJD(6) + t217) * t302;
t284 = qJD(1) * t292;
t334 = t303 * t370;
t333 = t258 * t362;
t327 = qJD(3) * t354;
t262 = qJD(4) * t305 - t303 * t327;
t263 = -qJD(4) * t303 - t305 * t327;
t223 = t262 * t300 + t263 * t298;
t239 = pkin(4) * t273 + pkin(8) * t276 + t334;
t332 = t304 * t223 + t302 * t239 + t242 * t339;
t291 = -pkin(3) * t300 - pkin(4);
t331 = t282 * t339;
t207 = t247 * t298 - t300 * t372;
t214 = t253 * t300 - t250;
t216 = t298 * t260 + t358;
t222 = t262 * t298 - t300 * t263;
t244 = t300 * t278 + t280 * t298;
t325 = t270 * t304;
t323 = -t195 * t304 - t197 * t302;
t322 = t207 * t282 - t245 * t268;
t321 = -t228 * t281 - t256 * t273;
t319 = qJ(6) * t276 - qJD(6) * t282;
t317 = 0.2e1 * qJD(3) * t284;
t199 = pkin(5) * t228 + t207;
t211 = -qJD(3) * pkin(4) - t214;
t315 = t331 - t362;
t314 = -t282 * t340 - t361;
t312 = t270 * t211 - t289 * t268;
t279 = t353 * t304;
t277 = t353 * t302;
t255 = t256 ^ 2;
t237 = t304 * t242;
t233 = t304 * t239;
t205 = pkin(5) * t256 + qJD(6) + t211;
t204 = -qJ(6) * t360 + t347;
t203 = pkin(5) * t281 - qJ(6) * t359 - t245 * t302 + t237;
t194 = -qJ(6) * t331 + (-qJD(5) * t245 + t319) * t302 + t332;
t193 = pkin(5) * t273 - t223 * t302 + t233 + t319 * t304 + (-t241 + (qJ(6) * t282 - t242) * t302) * qJD(5);
t1 = [0.2e1 * t329 * t376 - 0.2e1 * t336 * t373 + MDP(7) * t355 - MDP(8) * t356 + (-t290 * t355 + t303 * t317) * MDP(10) + (t290 * t356 + t305 * t317) * MDP(11) + (-t208 * t281 + t214 * t276 - t215 * t273 + t222 * t274 + t223 * t272 + t244 * t316 + t322) * MDP(12) + (t207 * t244 + t208 * t245 - t214 * t222 + t215 * t223 + (t271 + t311) * t334) * MDP(13) + (-t227 * t359 + t314 * t258) * MDP(14) + (t333 + (t369 + (-t365 + t367) * qJD(5)) * t282 + t351) * MDP(15) + (t282 * t266 + t314 * t270 + t348) * MDP(16) + (-t315 * t270 - t282 * t357 + t321) * MDP(17) + (t268 * t281 + t270 * t273) * MDP(18) + ((-t245 * t339 + t233) * t270 + t237 * t268 + (-t212 * t339 + t226) * t281 + t201 * t273 + t222 * t256 + t244 * t228 + t211 * t331 + ((-qJD(5) * t242 - t223) * t270 + (-qJD(5) * t224 - t208) * t281 - t211 * t276 + t322) * t302) * MDP(19) + (-(-t245 * t340 + t332) * t270 - t347 * t268 - t313 * t281 - t202 * t273 + t222 * t258 - t244 * t227 + t207 * t359 + t314 * t211) * MDP(20) + (-t193 * t258 - t194 * t256 + t203 * t227 - t204 * t228 - t323 * t276 + (-t191 * t304 - t192 * t302 + (t195 * t302 - t197 * t304) * qJD(5)) * t282) * MDP(21) + (t192 * t204 + t197 * t194 + t191 * t203 + t195 * t193 + t199 * (pkin(5) * t360 + t244) + t205 * (t315 * pkin(5) + t222)) * MDP(22); (-t276 * t272 + t273 * t274 + t281 * t316) * MDP(12) + (t207 * t281 - t214 * t273 - t215 * t276) * MDP(13) + (t270 * t362 - t321) * MDP(19) + (t270 * t361 + t348) * MDP(20) + (-t333 + t351) * MDP(21) + (t195 * t362 - t197 * t361 + t199 * t281 + t205 * t273) * MDP(22) + (-MDP(10) * t303 - MDP(11) * t305) * t306 + (t208 * MDP(13) - MDP(21) * t369 + (-t191 * t302 + t192 * t304) * MDP(22) + (-MDP(19) * t302 - MDP(20) * t304 - MDP(12)) * t268 + ((t365 + t367) * MDP(21) + t323 * MDP(22) + (-MDP(19) * t304 + MDP(20) * t302) * t270) * qJD(5)) * t282; ((t215 - t216) * t274 + (-t217 + t214) * t272 + (-t298 * t268 - t300 * t316) * pkin(3)) * MDP(12) + (t214 * t216 - t215 * t217 + (-t207 * t300 + t208 * t298 - t271 * t342) * pkin(3)) * MDP(13) + (t258 * t325 - t369) * MDP(14) + ((-t227 + t368) * t304 - t258 * t364 + t349) * MDP(15) + (-t258 * t274 + t270 * t325 + t357) * MDP(16) + (t256 * t274 - t270 * t340 + t346) * MDP(17) - t270 * t274 * MDP(18) + (-t201 * t274 - t207 * t304 - t216 * t256 + t291 * t228 + (-t289 * t339 - t232) * t270 + (t217 * t270 + t312) * t302) * MDP(19) + (t202 * t274 + t207 * t302 - t216 * t258 - t291 * t227 + (t289 * t340 + t350) * t270 + t312 * t304) * MDP(20) + (-t227 * t277 - t228 * t279 - t345 * t256 - t344 * t258 - t302 * t378 + t304 * t377) * MDP(21) + (t192 * t279 - t191 * t277 + t199 * (-pkin(5) * t304 + t291) + (pkin(5) * t364 - t216) * t205 + t345 * t197 + t344 * t195) * MDP(22) + (-t305 * t376 + t373) * qJD(1) ^ 2 + (-MDP(10) * t342 - MDP(11) * t341) * t284; -t272 ^ 2 * MDP(12) + (-t215 * t272 + t287) * MDP(13) + t346 * MDP(19) + t349 * MDP(21) + (-MDP(12) * t274 + MDP(13) * t214 - MDP(19) * t256 - MDP(20) * t258 - MDP(22) * t205) * t274 + (-qJD(5) * t270 * MDP(19) - t268 * MDP(20) + MDP(21) * t366 + MDP(22) * t377) * t302 + ((t227 + t368) * MDP(21) + t378 * MDP(22) - t270 ^ 2 * MDP(20)) * t304; t258 * t256 * MDP(14) + (-t255 + t371) * MDP(15) + (t256 * t270 - t227) * MDP(16) + (-t228 + t366) * MDP(17) + t268 * MDP(18) + (t202 * t270 - t211 * t258 + t308) * MDP(19) + (t201 * t270 + t211 * t256 - t313) * MDP(20) + (pkin(5) * t227 - t352 * t256) * MDP(21) + (t352 * t197 + (-t205 * t258 + t191) * pkin(5)) * MDP(22); (-t255 - t371) * MDP(21) + (t195 * t258 + t197 * t256 + t199) * MDP(22);];
tauc  = t1;
