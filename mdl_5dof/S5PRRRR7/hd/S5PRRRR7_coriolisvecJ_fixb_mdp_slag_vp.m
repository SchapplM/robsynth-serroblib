% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:13:04
% EndTime: 2019-12-05 17:13:10
% DurationCPUTime: 2.14s
% Computational Cost: add. (1415->221), mult. (3598->318), div. (0->0), fcn. (2626->8), ass. (0->122)
t283 = cos(qJ(3));
t280 = sin(qJ(2));
t338 = qJD(1) * t280;
t263 = qJD(2) * pkin(6) + t338;
t316 = pkin(7) * qJD(2) + t263;
t242 = t316 * t283;
t282 = cos(qJ(4));
t233 = t282 * t242;
t279 = sin(qJ(3));
t241 = t316 * t279;
t234 = qJD(3) * pkin(3) - t241;
t278 = sin(qJ(4));
t303 = -t234 * t278 - t233;
t335 = qJD(2) * t283;
t320 = t282 * t335;
t336 = qJD(2) * t279;
t321 = t278 * t336;
t246 = -t320 + t321;
t349 = pkin(8) * t246;
t197 = -t303 - t349;
t281 = cos(qJ(5));
t239 = t281 * t246;
t248 = -t278 * t335 - t282 * t336;
t277 = sin(qJ(5));
t214 = t248 * t277 - t239;
t272 = -pkin(3) * t283 - pkin(2);
t284 = cos(qJ(2));
t337 = qJD(1) * t284;
t252 = qJD(2) * t272 - t337;
t223 = pkin(4) * t246 + t252;
t330 = qJD(5) * t277;
t364 = t197 * t330 - t223 * t214;
t274 = qJD(3) + qJD(4);
t326 = qJD(2) * qJD(3);
t319 = t283 * t326;
t217 = qJD(4) * t320 - t274 * t321 + t282 * t319;
t334 = qJD(3) * t279;
t226 = -t263 * t334 + (-pkin(7) * t334 + t283 * t337) * qJD(2);
t333 = qJD(3) * t283;
t227 = -t263 * t333 + (-pkin(7) * t333 - t279 * t337) * qJD(2);
t313 = -t278 * t226 + t282 * t227;
t291 = qJD(4) * t303 + t313;
t187 = -pkin(8) * t217 + t291;
t325 = -qJD(4) - qJD(5);
t273 = qJD(3) - t325;
t363 = (-t197 * t273 - t187) * t277 + t364;
t256 = t278 * t283 + t279 * t282;
t361 = t274 * t256;
t218 = t361 * qJD(2);
t362 = -0.2e1 * t326;
t332 = qJD(4) * t278;
t312 = t278 * t227 - t242 * t332;
t357 = (qJD(4) * t234 + t226) * t282;
t186 = -pkin(8) * t218 + t312 + t357;
t302 = t246 * t277 + t281 * t248;
t298 = -t277 * t186 + t281 * t187 + t223 * t302;
t188 = -qJD(5) * t239 + t281 * t217 - t277 * t218 + t248 * t330;
t289 = qJD(5) * t302 - t217 * t277 - t281 * t218;
t360 = t214 * t302 * MDP(19) + (-t214 ^ 2 + t302 ^ 2) * MDP(20) + (-t214 * t273 + t188) * MDP(21) + (-t273 * t302 + t289) * MDP(22);
t359 = MDP(5) * t279;
t358 = MDP(6) * (t279 ^ 2 - t283 ^ 2);
t255 = t278 * t279 - t282 * t283;
t244 = t255 * t280;
t352 = pkin(6) + pkin(7);
t322 = qJD(3) * t352;
t257 = t279 * t322;
t258 = t283 * t322;
t296 = t256 * t284;
t259 = t352 * t279;
t260 = t352 * t283;
t301 = t259 * t278 - t260 * t282;
t356 = qJD(1) * t296 + qJD(4) * t301 + t278 * t257 - t282 * t258;
t295 = t255 * t284;
t331 = qJD(4) * t282;
t355 = -qJD(1) * t295 + t282 * t257 + t278 * t258 + t259 * t331 + t260 * t332;
t245 = t248 * pkin(8);
t231 = t278 * t242;
t311 = t282 * t234 - t231;
t196 = t245 + t311;
t323 = pkin(3) * t334;
t299 = t323 - t338;
t354 = t279 * MDP(10) + t283 * MDP(11);
t353 = qJD(5) - t273;
t351 = pkin(3) * t273;
t350 = pkin(4) * t248;
t348 = qJD(2) * pkin(2);
t346 = t252 * t248;
t345 = t278 * t281;
t285 = qJD(3) ^ 2;
t344 = t279 * t285;
t343 = t281 * t197;
t342 = t283 * t285;
t341 = -t282 * t241 - t231;
t253 = (t323 + t338) * qJD(2);
t324 = pkin(3) * t336;
t318 = -pkin(3) * t274 - t234;
t194 = pkin(4) * t274 + t196;
t317 = -pkin(4) * t273 - t194;
t310 = t241 * t278 - t233;
t309 = pkin(4) * t361 + t299;
t306 = pkin(8) * t361 - qJD(5) * (-pkin(8) * t256 - t259 * t282 - t260 * t278) + t355;
t224 = t274 * t255;
t305 = -pkin(8) * t224 + qJD(5) * (-pkin(8) * t255 - t301) - t356;
t304 = -t277 * t194 - t343;
t221 = t255 * t281 + t256 * t277;
t222 = -t255 * t277 + t256 * t281;
t297 = t252 * t246 - t312;
t292 = -0.2e1 * qJD(3) * t348;
t288 = -t248 * t246 * MDP(12) + (t246 * t274 + t217) * MDP(14) + (-t248 * t274 - t218) * MDP(15) + (-t246 ^ 2 + t248 ^ 2) * MDP(13) + t360;
t286 = qJD(2) ^ 2;
t271 = pkin(3) * t282 + pkin(4);
t243 = t256 * t280;
t238 = pkin(4) * t255 + t272;
t229 = t324 - t350;
t204 = pkin(4) * t218 + t253;
t201 = t245 + t341;
t200 = t310 + t349;
t199 = -qJD(2) * t296 + t244 * t274;
t198 = -qJD(2) * t295 - t280 * t361;
t191 = qJD(5) * t222 - t224 * t277 + t281 * t361;
t190 = -qJD(5) * t221 - t224 * t281 - t277 * t361;
t1 = [(MDP(17) * t199 - MDP(18) * t198) * t274 + ((-t198 * t277 + t199 * t281) * MDP(24) - (t198 * t281 + t199 * t277) * MDP(25) + ((t243 * t277 + t244 * t281) * MDP(24) - (-t243 * t281 + t244 * t277) * MDP(25)) * qJD(5)) * t273 + (-t218 * MDP(17) - t217 * MDP(18) + MDP(24) * t289 - t188 * MDP(25) - t286 * MDP(4) + t354 * t362) * t284 + (-t286 * MDP(3) + (t246 * MDP(17) - t248 * MDP(18) - MDP(24) * t214 - MDP(25) * t302) * qJD(2) + (-MDP(10) * t283 + MDP(11) * t279) * (t285 + t286)) * t280; 0.2e1 * t319 * t359 + t358 * t362 + MDP(7) * t342 - MDP(8) * t344 + (-pkin(6) * t342 + t279 * t292) * MDP(10) + (pkin(6) * t344 + t283 * t292) * MDP(11) + (t217 * t256 + t224 * t248) * MDP(12) + (-t217 * t255 - t218 * t256 + t224 * t246 + t248 * t361) * MDP(13) + (t272 * t218 + t299 * t246 + t252 * t361 + t253 * t255) * MDP(17) + (t272 * t217 - t252 * t224 - t299 * t248 + t253 * t256) * MDP(18) + (t188 * t222 - t190 * t302) * MDP(19) + (-t188 * t221 + t190 * t214 + t191 * t302 + t222 * t289) * MDP(20) + (t223 * t191 + t204 * t221 - t214 * t309 - t238 * t289) * MDP(24) + (t238 * t188 + t223 * t190 + t204 * t222 - t302 * t309) * MDP(25) + (-t224 * MDP(14) - MDP(15) * t361 + MDP(17) * t356 + MDP(18) * t355) * t274 + (t190 * MDP(21) - t191 * MDP(22) + (t277 * t306 - t281 * t305) * MDP(24) + (t277 * t305 + t281 * t306) * MDP(25)) * t273; (-(t200 * t281 - t201 * t277) * t273 + t229 * t214 + (-t277 * t282 - t345) * qJD(4) * t351 + ((-pkin(3) * t345 - t271 * t277) * t273 + t304) * qJD(5) + t298) * MDP(24) + (t341 * t274 + t248 * t324 + (qJD(4) * t318 - t226) * t282 + t297) * MDP(18) + (-t310 * t274 - t246 * t324 + t346 + (t278 * t318 - t233) * qJD(4) + t313) * MDP(17) + (t229 * t302 + (-t278 * t325 * t351 + t200 * t273 - t187) * t277 + (-qJD(5) * t194 - t186 + (-pkin(3) * t331 - qJD(5) * t271 + t201) * t273) * t281 + t364) * MDP(25) + t288 + t354 * qJD(2) * t348 + (-t283 * t359 + t358) * t286; (-t274 * t303 + t291 + t346) * MDP(17) + (t274 * t311 + t297 - t357) * MDP(18) + (-(-t196 * t277 - t343) * t273 - t214 * t350 + (t277 * t317 - t343) * qJD(5) + t298) * MDP(24) + (-t302 * t350 + (qJD(5) * t317 + t196 * t273 - t186) * t281 + t363) * MDP(25) + t288; (t304 * t353 + t298) * MDP(24) + ((-t194 * t353 - t186) * t281 + t363) * MDP(25) + t360;];
tauc = t1;
