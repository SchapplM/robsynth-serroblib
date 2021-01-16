% Calculate Coriolis joint torque vector for
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
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:19
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:18:25
% EndTime: 2021-01-15 20:18:33
% DurationCPUTime: 2.65s
% Computational Cost: add. (2577->269), mult. (6720->342), div. (0->0), fcn. (4796->6), ass. (0->134)
t319 = sin(pkin(8));
t322 = sin(qJ(2));
t347 = qJD(1) * qJD(2);
t344 = t322 * t347;
t307 = t319 * t344;
t320 = cos(pkin(8));
t323 = cos(qJ(2));
t343 = t323 * t347;
t278 = t320 * t343 - t307;
t321 = sin(qJ(4));
t377 = cos(qJ(4));
t345 = qJD(4) * t377;
t350 = qJD(4) * t321;
t297 = t319 * t323 + t320 * t322;
t353 = qJD(1) * t297;
t362 = t320 * t323;
t296 = t319 * t322 - t362;
t354 = qJD(1) * t296;
t356 = t319 * t343 + t320 * t344;
t333 = t377 * t278 - t321 * t356 - t345 * t354 - t350 * t353;
t243 = -t321 * t353 - t354 * t377;
t316 = qJD(2) + qJD(4);
t364 = t243 * t316;
t201 = t333 - t364;
t336 = t321 * t354 - t353 * t377;
t327 = qJD(4) * t336 - t321 * t278 - t377 * t356;
t365 = t243 ^ 2;
t366 = t336 * t316;
t389 = t336 ^ 2;
t391 = t243 * t336;
t392 = t201 * MDP(17) + MDP(15) * t391 + (t327 - t366) * MDP(18) + (-t365 + t389) * MDP(16);
t388 = 0.2e1 * pkin(2) * t322;
t312 = -pkin(2) * t323 - pkin(1);
t352 = qJD(1) * t312;
t304 = qJD(3) + t352;
t252 = pkin(3) * t354 + t304;
t206 = -pkin(4) * t243 + qJ(5) * t336 + t252;
t387 = t206 * t243;
t386 = t206 * t336;
t384 = t252 * t243;
t383 = t252 * t336;
t215 = -pkin(4) * t336 - qJ(5) * t243;
t382 = -0.2e1 * t347;
t381 = MDP(4) * t322;
t380 = MDP(5) * (t322 ^ 2 - t323 ^ 2);
t372 = -qJ(3) - pkin(6);
t305 = t372 * t322;
t301 = qJD(1) * t305;
t306 = t372 * t323;
t302 = qJD(1) * t306;
t363 = t320 * t302;
t250 = -t301 * t319 + t363;
t374 = pkin(7) * t354;
t234 = t250 + t374;
t291 = t319 * t302;
t251 = t320 * t301 + t291;
t373 = pkin(7) * t353;
t235 = t251 - t373;
t311 = pkin(2) * t320 + pkin(3);
t376 = pkin(2) * t319;
t346 = t321 * t376;
t379 = qJD(4) * t346 + t321 * t234 + t235 * t377 - t311 * t345;
t371 = qJD(2) * pkin(2);
t342 = qJD(2) * t372;
t284 = qJD(3) * t323 + t322 * t342;
t285 = -qJD(3) * t322 + t323 * t342;
t236 = -t319 * t284 + t320 * t285;
t331 = t296 * qJD(2);
t226 = pkin(7) * t331 + t236;
t237 = t320 * t284 + t319 * t285;
t332 = t297 * qJD(2);
t227 = -pkin(7) * t332 + t237;
t254 = t320 * t305 + t306 * t319;
t238 = -pkin(7) * t297 + t254;
t255 = t319 * t305 - t320 * t306;
t239 = -pkin(7) * t296 + t255;
t337 = t238 * t377 - t321 * t239;
t198 = qJD(4) * t337 + t321 * t226 + t227 * t377;
t370 = t198 * t316;
t212 = t321 * t238 + t239 * t377;
t199 = qJD(4) * t212 - t226 * t377 + t321 * t227;
t369 = t199 * t316;
t295 = t301 + t371;
t246 = t320 * t295 + t291;
t230 = qJD(2) * pkin(3) + t246 - t373;
t247 = t319 * t295 - t363;
t231 = t247 - t374;
t205 = t321 * t230 + t231 * t377;
t368 = t205 * t316;
t324 = qJD(2) ^ 2;
t361 = t322 * t324;
t360 = t323 * t324;
t325 = qJD(1) ^ 2;
t359 = t323 * t325;
t330 = t321 * t311 + t376 * t377;
t358 = -t330 * qJD(4) - t234 * t377 + t321 * t235;
t357 = -qJD(5) + t379;
t265 = t284 * qJD(1);
t266 = t285 * qJD(1);
t233 = t320 * t265 + t319 * t266;
t351 = qJD(1) * t322;
t204 = t230 * t377 - t321 * t231;
t348 = qJD(5) - t204;
t314 = t322 * t371;
t313 = pkin(2) * t351;
t259 = pkin(3) * t353 + t313;
t341 = pkin(1) * t382;
t232 = -t265 * t319 + t320 * t266;
t222 = -pkin(7) * t278 + t232;
t223 = -pkin(7) * t356 + t233;
t339 = t321 * t222 + t377 * t223 + t230 * t345 - t231 * t350;
t196 = -t377 * t222 + t321 * t223 + t230 * t350 + t231 * t345;
t310 = pkin(2) * t344;
t253 = pkin(3) * t356 + t310;
t315 = t316 * qJD(5);
t195 = t315 + t339;
t268 = pkin(3) * t296 + t312;
t249 = -t321 * t296 + t297 * t377;
t335 = t204 * t316 - t339;
t334 = t196 - t386;
t329 = t316 * t358 - t196;
t260 = pkin(3) * t332 + t314;
t197 = -pkin(4) * t327 - qJ(5) * t333 + qJD(5) * t336 + t253;
t282 = -t311 * t377 - pkin(4) + t346;
t281 = qJ(5) + t330;
t248 = t296 * t377 + t297 * t321;
t217 = qJD(4) * t249 - t321 * t331 + t332 * t377;
t216 = t296 * t345 + t297 * t350 + t321 * t332 + t331 * t377;
t210 = pkin(4) * t248 - qJ(5) * t249 + t268;
t209 = t215 + t259;
t203 = t316 * qJ(5) + t205;
t202 = -t316 * pkin(4) + t348;
t200 = t217 * pkin(4) + t216 * qJ(5) - t249 * qJD(5) + t260;
t1 = [0.2e1 * t343 * t381 + t380 * t382 + MDP(6) * t360 - MDP(7) * t361 + (-pkin(6) * t360 + t322 * t341) * MDP(9) + (pkin(6) * t361 + t323 * t341) * MDP(10) + (t312 * t356 + (t304 * t297 + t354 * t388 + t236) * qJD(2)) * MDP(11) + (t312 * t278 + (-t304 * t296 + t353 * t388 - t237) * qJD(2)) * MDP(12) + (-t237 * t354 - t255 * t356 - t233 * t296 - t236 * t353 - t254 * t278 - t232 * t297 + (t246 * t296 - t247 * t297) * qJD(2)) * MDP(13) + (t232 * t254 + t233 * t255 + t236 * t246 + t237 * t247 + (t304 + t352) * t314) * MDP(14) + (t216 * t336 + t249 * t333) * MDP(15) + (-t216 * t243 + t217 * t336 - t248 * t333 + t249 * t327) * MDP(16) + (t217 * t252 - t243 * t260 + t248 * t253 - t268 * t327 - t369) * MDP(20) + (-t216 * t252 + t249 * t253 - t260 * t336 + t268 * t333 - t370) * MDP(21) + (t197 * t248 - t200 * t243 + t206 * t217 - t210 * t327 - t369) * MDP(22) + (-t195 * t248 + t196 * t249 + t198 * t243 - t199 * t336 - t202 * t216 - t203 * t217 + t212 * t327 - t333 * t337) * MDP(23) + (-t197 * t249 + t200 * t336 + t206 * t216 - t210 * t333 + t370) * MDP(24) + (t195 * t212 - t196 * t337 + t197 * t210 + t198 * t203 + t199 * t202 + t200 * t206) * MDP(25) + (-t216 * MDP(17) - t217 * MDP(18)) * t316; -t359 * t381 + t325 * t380 + (-qJD(2) * t250 - t304 * t353 - t313 * t354 + t232) * MDP(11) + (qJD(2) * t251 + t304 * t354 - t313 * t353 - t233) * MDP(12) + ((t247 + t250) * t353 + (t251 - t246) * t354 + (-t320 * t278 - t319 * t356) * pkin(2)) * MDP(13) + (-t246 * t250 - t247 * t251 + (t232 * t320 + t233 * t319 - t304 * t351) * pkin(2)) * MDP(14) + (t243 * t259 + t329 + t383) * MDP(20) + (t259 * t336 + t379 * t316 - t339 - t384) * MDP(21) + (t209 * t243 + t329 + t386) * MDP(22) + (t327 * t281 + t282 * t333 + (-t203 + t358) * t336 + (-t202 - t357) * t243) * MDP(23) + (-t209 * t336 - t316 * t357 + t195 + t387) * MDP(24) + (t195 * t281 + t196 * t282 - t202 * t358 - t203 * t357 - t206 * t209) * MDP(25) + (MDP(9) * t322 * t325 + MDP(10) * t359) * pkin(1) + t392; (qJD(2) * t353 + t356) * MDP(11) + (-t307 + (qJD(1) * t362 - t354) * qJD(2)) * MDP(12) + (-t353 ^ 2 - t354 ^ 2) * MDP(13) + (t246 * t353 + t247 * t354 + t310) * MDP(14) + (-t365 - t389) * MDP(23) + (t202 * t336 - t203 * t243 + t197) * MDP(25) + (MDP(21) - MDP(24)) * (t333 + t364) + (MDP(20) + MDP(22)) * (-t327 - t366); (-t196 + t368 + t383) * MDP(20) + (t335 - t384) * MDP(21) + (t215 * t243 - t334 + t368) * MDP(22) + (-pkin(4) * t333 + qJ(5) * t327 - (t203 - t205) * t336 - (t202 - t348) * t243) * MDP(23) + (-t215 * t336 + 0.2e1 * t315 - t335 + t387) * MDP(24) + (-pkin(4) * t196 + qJ(5) * t195 - t202 * t205 + t203 * t348 - t206 * t215) * MDP(25) + t392; MDP(22) * t391 + t201 * MDP(23) + (-t316 ^ 2 - t389) * MDP(24) + (-t203 * t316 + t334) * MDP(25);];
tauc = t1;
