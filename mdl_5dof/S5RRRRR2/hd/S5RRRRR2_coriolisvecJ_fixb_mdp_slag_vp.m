% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [2x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S5RRRRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:35
% EndTime: 2019-03-29 15:26:39
% DurationCPUTime: 2.07s
% Computational Cost: add. (1913->228), mult. (3957->356), div. (0->0), fcn. (2870->8), ass. (0->117)
t268 = qJD(3) + qJD(4);
t269 = qJD(1) + qJD(2);
t273 = sin(qJ(4));
t274 = sin(qJ(3));
t277 = cos(qJ(3));
t358 = cos(qJ(4));
t288 = -t273 * t274 + t358 * t277;
t363 = t288 * t269;
t214 = t363 * t268;
t249 = t273 * t277 + t358 * t274;
t362 = t268 * t249;
t215 = t362 * t269;
t361 = MDP(7) * t277;
t275 = sin(qJ(2));
t360 = t275 * MDP(5);
t359 = (t274 ^ 2 - t277 ^ 2) * MDP(8);
t314 = qJD(4) * t358;
t356 = pkin(1) * qJD(1);
t321 = t277 * t356;
t278 = cos(qJ(2));
t329 = qJD(2) * t278;
t336 = t275 * t277;
t238 = (-qJD(3) * t336 - t274 * t329) * t356;
t322 = t275 * t356;
t310 = t274 * t322;
t301 = qJD(3) * t310;
t331 = -t358 * t238 - t273 * t301;
t253 = qJD(3) * pkin(2) - t310;
t341 = t253 * t273;
t201 = (t273 * t329 + t275 * t314) * t321 + qJD(4) * t341 + t331;
t357 = pkin(1) * t275;
t355 = pkin(1) * qJD(2);
t243 = t249 * t269;
t272 = sin(qJ(5));
t326 = qJD(5) * t272;
t276 = cos(qJ(5));
t325 = qJD(5) * t276;
t334 = t276 * t214 + t268 * t325;
t195 = -t243 * t326 + t334;
t354 = t195 * t272;
t353 = t214 * t272;
t352 = t215 * t272;
t351 = t215 * t276;
t350 = t215 * t277;
t345 = t243 * t272;
t221 = -t276 * t268 + t345;
t235 = qJD(5) - t363;
t349 = t221 * t235;
t294 = -t243 * t276 - t268 * t272;
t348 = t294 * t235;
t224 = t268 * t288;
t347 = t224 * t272;
t346 = t224 * t276;
t339 = t269 * t277;
t246 = -pkin(2) * t339 - t278 * t356;
t344 = t246 * t243;
t343 = t249 * t272;
t342 = t249 * t276;
t340 = t269 * t274;
t279 = qJD(3) ^ 2;
t335 = t277 * t279;
t320 = qJD(1) * t355;
t308 = t275 * t320;
t328 = qJD(3) * t274;
t316 = t269 * t328;
t244 = pkin(2) * t316 + t308;
t333 = t246 * t224 + t244 * t249;
t332 = -t244 * t288 + t246 * t362;
t327 = qJD(3) * t278;
t324 = qJD(5) * t277;
t323 = -qJD(1) - t269;
t315 = t274 * t327;
t309 = t275 * t321;
t303 = t273 * t309;
t228 = -t358 * t253 + t303;
t219 = t228 * t326;
t220 = t228 * t325;
t313 = t228 * (-t235 - t363);
t302 = t277 * t278 * t320;
t200 = -qJD(4) * t303 + t273 * t238 + t253 * t314 + (-t301 + t302) * t358;
t312 = -t200 * t272 + t276 * t244;
t311 = t235 * t276;
t255 = t358 * t309;
t229 = t255 + t341;
t210 = t229 * t276 + t246 * t272;
t297 = t229 * t272 - t246 * t276;
t299 = t200 * t276 + t244 * t272;
t307 = t201 * t342 + t228 * t346 + (-t297 * qJD(5) + t299) * t288 - t210 * t362;
t306 = (-qJD(2) + t269) * t356;
t304 = t201 * t272 + t210 * t243 + t220;
t300 = -(-t210 * qJD(5) + t312) * t288 + t201 * t343 - t297 * t362 + t228 * t347 + t249 * t220;
t233 = t249 * t322;
t298 = -t228 * t363 - t233 * t235;
t240 = t288 * t357;
t262 = -pkin(1) * t278 - pkin(2) * t277;
t296 = t240 * t276 + t262 * t272;
t295 = -t240 * t272 + t262 * t276;
t293 = qJD(5) * t273 + t340;
t292 = t278 * t306;
t291 = MDP(12) * t306;
t289 = -t201 * t276 + t243 * t297 + t219;
t287 = -t249 * t326 + t346;
t286 = t278 * t249;
t285 = t278 * t288;
t196 = -t294 * qJD(5) + t353;
t283 = ((t195 - t349) * t276 + (-t196 + t348) * t272) * MDP(22) + (-t294 * t311 + t354) * MDP(21) + (-t235 ^ 2 * t272 + t221 * t243 + t351) * MDP(24) + (t235 * t311 + t243 * t294 + t352) * MDP(23) + (t243 * t268 - t215) * MDP(17) + (t243 ^ 2 - t363 ^ 2) * MDP(15) - (MDP(14) * t363 + MDP(25) * t235) * t243;
t282 = -t279 * t274 * MDP(10) + ((-t221 * t276 + t272 * t294) * t224 + (-t354 - t196 * t276 + (t221 * t272 + t276 * t294) * qJD(5)) * t249) * MDP(22) + (-t195 * t288 + t215 * t342 + t287 * t235 - t294 * t362) * MDP(23) + (-t215 * t343 + t196 * t288 - t221 * t362 + (-t249 * t325 - t347) * t235) * MDP(24) + (t214 * t288 - t215 * t249 + t224 * t363 - t243 * t362) * MDP(15) + (t195 * t342 - t287 * t294) * MDP(21) + (-t215 * t288 + t235 * t362) * MDP(25) + (t214 * t249 + t224 * t243) * MDP(14) - 0.2e1 * t269 * qJD(3) * t359 + 0.2e1 * t316 * t361 + MDP(9) * t335 + (t224 * MDP(16) - MDP(17) * t362) * t268;
t280 = -t246 * t363 - t200;
t258 = t274 * t308;
t254 = pkin(2) * t328 + t275 * t355;
t239 = t249 * t357;
t234 = t285 * t356;
t232 = t286 * t356;
t231 = t273 * t310 - t255;
t205 = (qJD(2) * t286 + t224 * t275) * pkin(1);
t204 = (qJD(2) * t285 - t275 * t362) * pkin(1);
t1 = [t282 + (-t205 * t268 + t262 * t215 - t254 * t363 + t332) * MDP(19) + (-t204 * t268 + t214 * t262 + t243 * t254 + t333) * MDP(20) + ((-t296 * qJD(5) - t204 * t272 + t254 * t276) * t235 + t295 * t215 + t205 * t221 + t239 * t196 + t300) * MDP(26) + (-(t204 * t276 + t254 * t272) * t235 - t296 * t215 - t205 * t294 + t239 * t195 + (-t228 * t343 - t295 * t235) * qJD(5) + t307) * MDP(27) + (-t275 * t335 + t323 * t315 + (t323 * t336 - t315) * qJD(2)) * pkin(1) * MDP(12) + (t258 + ((qJD(2) * t269 + t279) * t275 * t274 - 0.2e1 * t327 * t339) * pkin(1)) * MDP(13) + (t278 * MDP(6) + t360) * t323 * t355; t282 + (-t269 * t310 + t258) * MDP(13) + t291 * t336 + (-(-t234 * t272 + t276 * t322) * t235 - t232 * t221 + (-t276 * t350 + (t272 * t324 + t276 * t328) * t235) * pkin(2) + t300) * MDP(26) + t306 * t360 + MDP(6) * t292 + (t363 * t322 + t232 * t268 + (-t328 * t363 - t350) * pkin(2) + t332) * MDP(19) + (-t249 * t219 + (t234 * t276 + t272 * t322) * t235 + t232 * t294 + (t272 * t350 + (-t272 * t328 + t276 * t324) * t235) * pkin(2) + t307) * MDP(27) + (-t243 * t322 + t234 * t268 + (-t214 * t277 + t243 * t328) * pkin(2) + t333) * MDP(20); t283 + (-t233 * t268 + (-t243 * t340 - t268 * t314) * pkin(2) + t280) * MDP(20) + t278 * t274 * t291 + t277 * MDP(13) * t292 + (-t273 * t302 + pkin(2) * t363 * t340 - t231 * t268 - t344 + (-t255 + (-pkin(2) * t268 - t253) * t273) * qJD(4) - t331) * MDP(19) + (-t231 * t294 + t298 * t276 + (-t358 * t195 + (-qJD(4) * t294 - t351) * t273 + (t293 * t272 - t276 * t314) * t235) * pkin(2) + t304) * MDP(27) + (t231 * t221 + t298 * t272 + (-t358 * t196 + (qJD(4) * t221 - t352) * t273 + (-t272 * t314 - t293 * t276) * t235) * pkin(2) + t289) * MDP(26) + (-t274 * t361 + t359) * t269 ^ 2; (t229 * t268 - t201 - t344) * MDP(19) + (-t228 * t268 + t280) * MDP(20) + (-t221 * t229 + t272 * t313 + t289) * MDP(26) + (t229 * t294 + t276 * t313 + t304) * MDP(27) + t283; -t294 * t221 * MDP(21) + (-t221 ^ 2 + t294 ^ 2) * MDP(22) + (t334 + t349) * MDP(23) + (-t348 - t353) * MDP(24) + t215 * MDP(25) + (t210 * t235 + t228 * t294 + t312) * MDP(26) + (t221 * t228 - t235 * t297 - t299) * MDP(27) + (-MDP(23) * t345 + t294 * MDP(24) - t210 * MDP(26) + t297 * MDP(27)) * qJD(5);];
tauc  = t1;
