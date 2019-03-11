% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPPRPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:49:25
% EndTime: 2019-03-09 01:49:32
% DurationCPUTime: 2.72s
% Computational Cost: add. (1541->315), mult. (3205->439), div. (0->0), fcn. (1985->6), ass. (0->142)
t291 = sin(qJ(4));
t352 = qJD(1) * t291;
t278 = qJD(6) + t352;
t286 = sin(pkin(9));
t293 = cos(qJ(4));
t351 = qJD(1) * t293;
t332 = t286 * t351;
t287 = cos(pkin(9));
t340 = t287 * qJD(4);
t253 = t332 - t340;
t331 = t287 * t351;
t345 = qJD(4) * t286;
t255 = t331 + t345;
t290 = sin(qJ(6));
t292 = cos(qJ(6));
t305 = t253 * t290 - t255 * t292;
t376 = t278 * t305;
t288 = -pkin(7) + qJ(2);
t343 = qJD(4) * t293;
t350 = qJD(2) * t291;
t375 = t288 * t343 + t350;
t341 = qJD(6) * t292;
t342 = qJD(6) * t290;
t374 = -t286 * t342 + t287 * t341;
t284 = t291 ^ 2;
t373 = MDP(11) * (-t293 ^ 2 + t284);
t289 = pkin(1) + qJ(3);
t372 = qJD(1) * t289;
t371 = qJD(6) * t278;
t258 = t286 * t292 + t287 * t290;
t300 = t258 * qJD(6);
t370 = qJ(2) * MDP(6) + MDP(5) + MDP(7);
t369 = pkin(8) + qJ(5);
t363 = t255 * t290;
t212 = t292 * t253 + t363;
t367 = t212 * t278;
t281 = qJD(1) * qJ(2) + qJD(3);
t274 = -pkin(7) * qJD(1) + t281;
t344 = qJD(4) * t291;
t259 = t274 * t344;
t338 = qJD(1) * qJD(2);
t236 = -t293 * t338 + t259;
t366 = t236 * t286;
t365 = t236 * t287;
t364 = t236 * t293;
t362 = t274 * t293;
t361 = t286 * t293;
t360 = t287 * t293;
t359 = t288 * t291;
t264 = t291 * t274;
t358 = t292 * t287;
t315 = pkin(4) * t293 + qJ(5) * t291;
t240 = qJD(4) * t315 - qJD(5) * t293 + qJD(3);
t225 = t240 * qJD(1);
t228 = t291 * t338 + (qJD(5) + t362) * qJD(4);
t200 = t286 * t225 + t287 * t228;
t265 = pkin(4) * t291 - qJ(5) * t293 + t289;
t235 = qJD(1) * t265 - qJD(2);
t248 = qJD(4) * qJ(5) + t264;
t205 = t286 * t235 + t287 * t248;
t333 = t286 * t352;
t357 = -t290 * t333 + t352 * t358 + t374;
t299 = t258 * qJD(1);
t356 = t291 * t299 + t300;
t260 = t315 * qJD(1);
t218 = t286 * t260 + t274 * t360;
t337 = qJD(1) * qJD(4);
t326 = t291 * t337;
t317 = t292 * t326;
t318 = t290 * t326;
t355 = -t286 * t317 - t287 * t318;
t224 = t286 * t265 + t287 * t359;
t294 = qJD(4) ^ 2;
t295 = qJD(1) ^ 2;
t353 = -t294 - t295;
t349 = qJD(2) * t293;
t257 = t286 * t290 - t358;
t348 = qJD(4) * t257;
t347 = qJD(4) * t258;
t346 = qJD(4) * t278;
t275 = -qJD(2) + t372;
t339 = qJD(2) - t275;
t336 = pkin(8) * t287 * t291;
t334 = 0.2e1 * qJD(3) * qJD(1);
t209 = t286 * t240 + t287 * t375;
t330 = t286 * t344;
t325 = MDP(25) * t351;
t324 = -t286 * t288 + pkin(5);
t323 = pkin(5) * t286 - t288;
t322 = qJD(4) * pkin(4) - qJD(5);
t199 = t287 * t225 - t228 * t286;
t298 = (pkin(5) * t293 + t336) * qJD(1);
t193 = qJD(4) * t298 + t199;
t319 = pkin(8) * t330;
t194 = qJD(1) * t319 + t200;
t321 = t292 * t193 - t194 * t290;
t204 = t287 * t235 - t248 * t286;
t320 = t339 * qJD(1);
t217 = t287 * t260 - t274 * t361;
t242 = -t322 - t362;
t316 = t288 * t351 - t242;
t313 = t193 * t290 + t194 * t292;
t195 = pkin(5) * t352 - pkin(8) * t255 + t204;
t198 = -pkin(8) * t253 + t205;
t190 = t195 * t292 - t198 * t290;
t191 = t195 * t290 + t198 * t292;
t312 = -t199 * t287 - t200 * t286;
t311 = -t199 * t286 + t200 * t287;
t310 = -t204 * t287 - t205 * t286;
t309 = t204 * t286 - t205 * t287;
t250 = t287 * t265;
t211 = -pkin(8) * t360 + t291 * t324 + t250;
t215 = -pkin(8) * t361 + t224;
t308 = t211 * t292 - t215 * t290;
t307 = t211 * t290 + t215 * t292;
t306 = t253 * t287 - t255 * t286;
t304 = -t288 * t294 + t334;
t272 = t369 * t287;
t303 = qJD(5) * t286 + qJD(6) * t272 + t217 + t298;
t271 = t369 * t286;
t302 = pkin(8) * t333 - qJD(5) * t287 + qJD(6) * t271 + t218;
t301 = qJD(1) * t257;
t297 = -t253 * t341 + t286 * t318 - t287 * t317;
t197 = -qJD(6) * t305 + t355;
t296 = -qJ(5) * t343 + (t242 + t322) * t291;
t196 = -t255 * t342 + t297;
t280 = -pkin(5) * t287 - pkin(4);
t252 = t323 * t293;
t239 = t257 * t293;
t238 = t258 * t293;
t234 = -pkin(5) * t333 + t264;
t231 = -t323 * t344 - t349;
t230 = t287 * t240;
t223 = -t286 * t359 + t250;
t219 = t259 + (-pkin(5) * t330 - t349) * qJD(1);
t216 = pkin(5) * t253 + t242;
t208 = -t286 * t375 + t230;
t207 = -t290 * t291 * t340 - t292 * t330 + t293 * t374;
t206 = t257 * t344 - t293 * t300;
t202 = t319 + t209;
t201 = -t286 * t350 + t230 + (t293 * t324 + t336) * qJD(4);
t1 = [MDP(8) * t334 + (qJD(2) * t281 + qJD(3) * t275 + (qJ(2) * qJD(2) + qJD(3) * t289) * qJD(1)) * MDP(9) + 0.2e1 * t337 * t373 + (-t208 * t255 - t209 * t253 + ((t223 * t287 + t224 * t286) * qJD(1) - t310) * t344) * MDP(19) + (-t288 * t364 + t199 * t223 + t200 * t224 + t204 * t208 + t205 * t209 + (t288 * t344 - t349) * t242) * MDP(20) + (-t196 * t239 - t206 * t305) * MDP(21) + (-t196 * t238 + t197 * t239 - t206 * t212 + t207 * t305) * MDP(22) + (t206 * t278 + (-qJD(1) * t239 - t305) * t343) * MDP(23) + (-t207 * t278 + (-qJD(1) * t238 - t212) * t343) * MDP(24) + (t278 + t352) * MDP(25) * t343 + ((t201 * t292 - t202 * t290) * t278 + t231 * t212 + t252 * t197 + t219 * t238 + t216 * t207 - t307 * t371 + (qJD(1) * t308 + t190) * t343) * MDP(26) + (-(t201 * t290 + t202 * t292) * t278 - t231 * t305 + t252 * t196 - t219 * t239 + t216 * t206 - t308 * t371 + (-qJD(1) * t307 - t191) * t343) * MDP(27) + (MDP(15) * t343 - MDP(16) * t344) * (qJD(2) + t275 + t372) + (-t294 * MDP(12) + t304 * MDP(15) + (qJD(1) * t208 + t199 + (t253 * t288 + t286 * t316) * qJD(4)) * MDP(17) + (-qJD(1) * t209 - t200 + (t255 * t288 + t287 * t316) * qJD(4)) * MDP(18) + t196 * MDP(23) - t197 * MDP(24) + (-qJD(6) * t191 + t321) * MDP(26) + (-qJD(6) * t190 - t313) * MDP(27)) * t291 + (-0.2e1 * MDP(10) * t326 - t294 * MDP(13) + t304 * MDP(16) + (-qJD(2) * t253 + t366 + (qJD(1) * t223 + t204) * qJD(4)) * MDP(17) + (-qJD(2) * t255 + t365 + (-qJD(1) * t224 - t205) * qJD(4)) * MDP(18) + t312 * MDP(19)) * t293 + 0.2e1 * t370 * t338; t312 * MDP(20) + (MDP(26) * t356 + MDP(27) * t357) * t278 + ((MDP(17) * t286 + MDP(18) * t287) * t284 - t370) * t295 + ((-qJD(3) - t281) * MDP(9) + (t306 * MDP(19) + t309 * MDP(20) + (0.2e1 * MDP(16) + (-t286 ^ 2 - t287 ^ 2) * MDP(19)) * qJD(4)) * t291 + (-0.2e1 * qJD(4) * MDP(15) + (t253 - t340) * MDP(17) + (t255 + t345) * MDP(18) + t242 * MDP(20) + (t212 + t348) * MDP(26) + (-t305 + t347) * MDP(27)) * t293) * qJD(1); -t295 * MDP(8) + MDP(9) * t320 + t353 * t293 * MDP(16) + (-t306 * t343 + (t253 * t286 + t255 * t287) * qJD(1)) * MDP(19) + (t310 * qJD(1) - t309 * t343 - t364) * MDP(20) + (t278 * t301 + (-t258 * t346 - t197) * t293) * MDP(26) + (t278 * t299 + (t257 * t346 - t196) * t293) * MDP(27) + (t353 * MDP(15) + (-t287 * t295 + (t253 - t332) * qJD(4)) * MDP(17) + (t286 * t295 + (t255 - t331) * qJD(4)) * MDP(18) + (qJD(4) * t242 + t311) * MDP(20) + (t257 * t371 + (-t258 * t351 + t212) * qJD(4)) * MDP(26) + (t278 * t300 + (t293 * t301 - t305) * qJD(4)) * MDP(27)) * t291; t293 * MDP(15) * t320 - t339 * MDP(16) * t352 + (-t253 * t264 - t365 + (-t204 * t293 - t217 * t291 + t286 * t296) * qJD(1)) * MDP(17) + (-t255 * t264 + t366 + (t205 * t293 + t218 * t291 + t287 * t296) * qJD(1)) * MDP(18) + (t217 * t255 + t218 * t253 + (-qJD(5) * t253 - t204 * t352 + t200) * t287 + (qJD(5) * t255 - t205 * t352 - t199) * t286) * MDP(19) + (-pkin(4) * t236 + qJ(5) * t311 - qJD(5) * t309 - t204 * t217 - t205 * t218 - t242 * t264) * MDP(20) + (t196 * t258 - t305 * t357) * MDP(21) + (-t196 * t257 - t197 * t258 - t212 * t357 + t305 * t356) * MDP(22) + (t357 * t278 + (t305 + t347) * t351) * MDP(23) + (-t356 * t278 + (t212 - t348) * t351) * MDP(24) - t278 * t325 + (t280 * t197 - t234 * t212 + t219 * t257 + (t290 * t302 - t292 * t303) * t278 + t356 * t216 + ((-t271 * t292 - t272 * t290) * qJD(4) - t190) * t351) * MDP(26) + (t280 * t196 + t234 * t305 + t219 * t258 + (t290 * t303 + t292 * t302) * t278 + t357 * t216 + (-(-t271 * t290 + t272 * t292) * qJD(4) + t191) * t351) * MDP(27) + (t293 * t291 * MDP(10) - t373) * t295; (-t253 ^ 2 - t255 ^ 2) * MDP(19) + (t204 * t255 + t205 * t253 + t236) * MDP(20) + (t197 - t376) * MDP(26) + (t196 - t367) * MDP(27) + ((t255 - t345) * MDP(17) + (-t253 - t340) * MDP(18)) * t352; -t305 * t212 * MDP(21) + (-t212 ^ 2 + t305 ^ 2) * MDP(22) + (t297 + t367) * MDP(23) + (-t355 - t376) * MDP(24) + qJD(4) * t325 + (t191 * t278 + t216 * t305 + t321) * MDP(26) + (t190 * t278 + t212 * t216 - t313) * MDP(27) + (-MDP(23) * t363 + MDP(24) * t305 - MDP(26) * t191 - MDP(27) * t190) * qJD(6);];
tauc  = t1;
