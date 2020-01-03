% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:13
% EndTime: 2019-12-31 19:33:22
% DurationCPUTime: 3.81s
% Computational Cost: add. (2439->312), mult. (6425->449), div. (0->0), fcn. (4696->8), ass. (0->144)
t335 = cos(qJ(2));
t382 = cos(pkin(8));
t352 = t382 * t335;
t320 = qJD(1) * t352;
t330 = sin(pkin(8));
t333 = sin(qJ(2));
t359 = t333 * qJD(1);
t294 = t330 * t359 - t320;
t290 = qJD(5) + t294;
t310 = t330 * t335 + t333 * t382;
t297 = t310 * qJD(1);
t329 = sin(pkin(9));
t331 = cos(pkin(9));
t278 = t331 * qJD(2) - t297 * t329;
t334 = cos(qJ(5));
t390 = t278 * t334;
t332 = sin(qJ(5));
t311 = t329 * t334 + t331 * t332;
t364 = t290 * t311;
t277 = qJD(2) * t329 + t297 * t331;
t236 = t277 * t334 + t278 * t332;
t358 = qJD(1) * qJD(2);
t389 = -0.2e1 * t358;
t388 = MDP(5) * (t333 ^ 2 - t335 ^ 2);
t296 = t310 * qJD(2);
t287 = qJD(1) * t296;
t369 = t334 * t331;
t371 = t329 * t332;
t309 = -t369 + t371;
t365 = t290 * t309;
t387 = -t287 * t311 + t290 * t365;
t386 = pkin(7) * t331;
t322 = pkin(2) * t330 + qJ(4);
t385 = pkin(7) + t322;
t384 = -qJ(3) - pkin(6);
t383 = qJD(2) * pkin(2);
t234 = t277 * t332 - t390;
t381 = t234 * t297;
t380 = t236 * t297;
t354 = qJD(2) * t384;
t291 = qJD(3) * t335 + t333 * t354;
t283 = t291 * qJD(1);
t292 = -qJD(3) * t333 + t335 * t354;
t284 = t292 * qJD(1);
t241 = t283 * t330 - t382 * t284;
t317 = t384 * t333;
t318 = t384 * t335;
t273 = -t382 * t317 - t318 * t330;
t379 = t241 * t273;
t355 = t333 * t358;
t288 = qJD(2) * t320 - t330 * t355;
t377 = t288 * t329;
t376 = t288 * t331;
t375 = t294 * t329;
t340 = -t330 * t333 + t352;
t299 = t340 * qJD(2);
t374 = t299 * t329;
t373 = t310 * t329;
t372 = t310 * t331;
t314 = qJD(1) * t318;
t302 = t330 * t314;
t336 = qJD(2) ^ 2;
t370 = t333 * t336;
t368 = t335 * t336;
t337 = qJD(1) ^ 2;
t367 = t335 * t337;
t321 = pkin(2) * t355;
t229 = pkin(3) * t287 - qJ(4) * t288 - qJD(4) * t297 + t321;
t242 = t382 * t283 + t330 * t284;
t239 = qJD(2) * qJD(4) + t242;
t207 = t329 * t229 + t331 * t239;
t357 = t333 * t383;
t238 = pkin(3) * t296 - qJ(4) * t299 - qJD(4) * t310 + t357;
t255 = t291 * t382 + t330 * t292;
t214 = t329 * t238 + t331 * t255;
t356 = -pkin(2) * t335 - pkin(1);
t348 = t356 * qJD(1);
t316 = qJD(3) + t348;
t245 = pkin(3) * t294 - qJ(4) * t297 + t316;
t313 = qJD(1) * t317;
t307 = t313 + t383;
t353 = t382 * t314;
t265 = t330 * t307 - t353;
t261 = qJD(2) * qJ(4) + t265;
t218 = t329 * t245 + t331 * t261;
t253 = pkin(2) * t359 + pkin(3) * t297 + qJ(4) * t294;
t269 = t313 * t382 + t302;
t223 = t329 * t253 + t331 * t269;
t263 = -pkin(3) * t340 - qJ(4) * t310 + t356;
t274 = t330 * t317 - t318 * t382;
t226 = t329 * t263 + t331 * t274;
t360 = qJD(5) * t334;
t366 = t278 * t360 + t288 * t369;
t210 = pkin(7) * t278 + t218;
t362 = qJD(5) * t210;
t361 = qJD(5) * t310;
t351 = pkin(1) * t389;
t206 = t331 * t229 - t239 * t329;
t213 = t331 * t238 - t255 * t329;
t217 = t331 * t245 - t261 * t329;
t222 = t331 * t253 - t269 * t329;
t225 = t331 * t263 - t274 * t329;
t254 = t291 * t330 - t382 * t292;
t268 = t313 * t330 - t353;
t203 = -pkin(7) * t377 + t207;
t205 = pkin(4) * t294 - pkin(7) * t277 + t217;
t350 = -qJD(5) * t205 - t203;
t325 = -pkin(2) * t382 - pkin(3);
t349 = -t309 * t287 - t364 * t290;
t264 = t307 * t382 + t302;
t199 = t205 * t334 - t210 * t332;
t200 = t205 * t332 + t210 * t334;
t215 = -pkin(4) * t340 - pkin(7) * t372 + t225;
t219 = -pkin(7) * t373 + t226;
t347 = t215 * t334 - t219 * t332;
t346 = t215 * t332 + t219 * t334;
t345 = -t217 * t329 + t218 * t331;
t344 = t241 * t310 + t273 * t288;
t343 = -qJD(5) * t277 - t377;
t306 = t385 * t331;
t342 = pkin(4) * t297 + qJD(4) * t329 + qJD(5) * t306 + t294 * t386 + t222;
t305 = t385 * t329;
t341 = pkin(7) * t375 - qJD(4) * t331 + qJD(5) * t305 + t223;
t256 = -qJD(2) * pkin(3) + qJD(4) - t264;
t339 = t256 * t299 + t344;
t338 = -t287 * t322 + t288 * t325 + (-qJD(4) + t256) * t294;
t212 = qJD(5) * t236 + t288 * t311;
t315 = -t331 * pkin(4) + t325;
t293 = t294 ^ 2;
t258 = t309 * t310;
t257 = t311 * t310;
t251 = pkin(4) * t373 + t273;
t240 = -pkin(4) * t375 + t268;
t231 = pkin(4) * t374 + t254;
t230 = -pkin(4) * t278 + t256;
t224 = pkin(4) * t377 + t241;
t221 = t299 * t311 + t360 * t372 - t361 * t371;
t220 = -t299 * t309 - t311 * t361;
t211 = t343 * t332 + t366;
t208 = -pkin(7) * t374 + t214;
t204 = pkin(4) * t296 - t299 * t386 + t213;
t202 = pkin(4) * t287 - pkin(7) * t376 + t206;
t201 = t334 * t202;
t1 = [0.2e1 * t335 * MDP(4) * t355 + t388 * t389 + MDP(6) * t368 - MDP(7) * t370 + (-pkin(6) * t368 + t333 * t351) * MDP(9) + (pkin(6) * t370 + t335 * t351) * MDP(10) + (t242 * t340 + t254 * t297 - t255 * t294 - t264 * t299 - t265 * t296 - t274 * t287 + t344) * MDP(11) + (t379 + t242 * t274 - t264 * t254 + t265 * t255 + (t316 + t348) * t357) * MDP(12) + (-t206 * t340 + t213 * t294 + t217 * t296 + t225 * t287 - t254 * t278 + t329 * t339) * MDP(13) + (t207 * t340 - t214 * t294 - t218 * t296 - t226 * t287 + t254 * t277 + t331 * t339) * MDP(14) + (-t213 * t277 + t214 * t278 + (-t206 * t310 - t217 * t299 - t225 * t288) * t331 + (-t207 * t310 - t218 * t299 - t226 * t288) * t329) * MDP(15) + (t206 * t225 + t207 * t226 + t213 * t217 + t214 * t218 + t254 * t256 + t379) * MDP(16) + (-t211 * t258 + t220 * t236) * MDP(17) + (-t211 * t257 + t212 * t258 - t220 * t234 - t221 * t236) * MDP(18) + (-t211 * t340 + t220 * t290 + t236 * t296 - t258 * t287) * MDP(19) + (t212 * t340 - t221 * t290 - t234 * t296 - t257 * t287) * MDP(20) + (-t287 * t340 + t290 * t296) * MDP(21) + ((t204 * t334 - t208 * t332) * t290 + t347 * t287 - (-t203 * t332 + t201) * t340 + t199 * t296 + t231 * t234 + t251 * t212 + t224 * t257 + t230 * t221 + (t200 * t340 - t290 * t346) * qJD(5)) * MDP(22) + (-(t204 * t332 + t208 * t334) * t290 - t346 * t287 + (t202 * t332 + t203 * t334) * t340 - t200 * t296 + t231 * t236 + t251 * t211 - t224 * t258 + t230 * t220 + (t199 * t340 - t290 * t347) * qJD(5)) * MDP(23); -t333 * MDP(4) * t367 + t337 * t388 + ((t265 - t268) * t297 + (-t264 + t269) * t294 + (-t287 * t330 - t288 * t382) * pkin(2)) * MDP(11) + (t264 * t268 - t265 * t269 + (-t241 * t382 + t242 * t330 - t316 * t359) * pkin(2)) * MDP(12) + (-t217 * t297 - t222 * t294 - t241 * t331 + t268 * t278 + t329 * t338) * MDP(13) + (t218 * t297 + t223 * t294 + t241 * t329 - t268 * t277 + t331 * t338) * MDP(14) + (t222 * t277 - t223 * t278 + (qJD(4) * t278 - t217 * t294 + t207) * t331 + (qJD(4) * t277 - t218 * t294 - t206) * t329) * MDP(15) + (-t217 * t222 - t218 * t223 + t241 * t325 - t256 * t268 + (-t206 * t329 + t207 * t331) * t322 + t345 * qJD(4)) * MDP(16) + (t211 * t311 - t236 * t365) * MDP(17) + (-t211 * t309 - t212 * t311 + t234 * t365 - t236 * t364) * MDP(18) + (-t380 - t387) * MDP(19) + (t349 + t381) * MDP(20) - t290 * t297 * MDP(21) + ((-t305 * t334 - t306 * t332) * t287 + t315 * t212 + t224 * t309 - t199 * t297 - t240 * t234 + (t332 * t341 - t334 * t342) * t290 + t364 * t230) * MDP(22) + (-(-t305 * t332 + t306 * t334) * t287 + t315 * t211 + t224 * t311 + t200 * t297 - t240 * t236 + (t332 * t342 + t334 * t341) * t290 - t365 * t230) * MDP(23) + (MDP(9) * t333 * t337 + MDP(10) * t367) * pkin(1); (-t297 ^ 2 - t293) * MDP(11) + (t264 * t297 + t321) * MDP(12) + (t278 * t297 + t287 * t331) * MDP(13) + (-t277 * t297 - t287 * t329 - t293 * t331) * MDP(14) + (t206 * t331 + t207 * t329 - t256 * t297) * MDP(16) + (t349 - t381) * MDP(22) + (-t380 + t387) * MDP(23) + (-t329 ^ 2 - t331 ^ 2) * MDP(15) * t288 + (t265 * MDP(12) + (t277 * t329 + t278 * t331) * MDP(15) + t345 * MDP(16) - MDP(13) * t375) * t294; (t277 * t294 + t377) * MDP(13) + (t278 * t294 + t376) * MDP(14) + (-t277 ^ 2 - t278 ^ 2) * MDP(15) + (t217 * t277 - t218 * t278 + t241) * MDP(16) + (t236 * t290 + t212) * MDP(22) + (t290 * t390 + (-t277 * t290 + t343) * t332 + t366) * MDP(23); -t234 ^ 2 * MDP(18) + (t234 * t290 + t366) * MDP(19) + t287 * MDP(21) + (t200 * t290 + t201) * MDP(22) + (t199 * t290 + t230 * t234) * MDP(23) + (MDP(17) * t234 + t236 * MDP(18) + t290 * MDP(20) - t230 * MDP(22)) * t236 + (MDP(20) * t343 - MDP(22) * t362 + MDP(23) * t350) * t334 + (t343 * MDP(19) + (-qJD(5) * t278 - t376) * MDP(20) + t350 * MDP(22) + (-t202 + t362) * MDP(23)) * t332;];
tauc = t1;
