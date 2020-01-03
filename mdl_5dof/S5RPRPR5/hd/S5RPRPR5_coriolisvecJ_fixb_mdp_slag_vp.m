% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RPRPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:43:10
% EndTime: 2020-01-03 11:43:18
% DurationCPUTime: 2.49s
% Computational Cost: add. (1679->225), mult. (4661->333), div. (0->0), fcn. (3368->8), ass. (0->131)
t301 = sin(pkin(8));
t300 = sin(pkin(9));
t302 = cos(pkin(9));
t305 = sin(qJ(3));
t307 = cos(qJ(3));
t279 = t300 * t307 + t302 * t305;
t314 = qJD(1) * t279;
t262 = t301 * t314;
t306 = cos(qJ(5));
t257 = t306 * t262;
t318 = t300 * t305 - t302 * t307;
t344 = qJD(3) * t301;
t264 = t318 * t344;
t260 = qJD(1) * t264;
t268 = t279 * t344;
t261 = qJD(1) * t268;
t348 = qJD(1) * t301;
t333 = t307 * t348;
t334 = t305 * t348;
t265 = -t300 * t334 + t302 * t333;
t304 = sin(qJ(5));
t340 = qJD(5) * t304;
t201 = -qJD(5) * t257 + t304 * t260 - t306 * t261 - t265 * t340;
t320 = -t262 * t304 + t306 * t265;
t202 = qJD(5) * t320 - t306 * t260 - t261 * t304;
t226 = t265 * t304 + t257;
t303 = cos(pkin(8));
t347 = qJD(1) * t303;
t289 = -qJD(3) + t347;
t285 = -qJD(5) + t289;
t363 = t226 * t285;
t364 = t320 * t285;
t385 = t226 * MDP(17) * t320 + (-t202 - t364) * MDP(20) + (-t226 ^ 2 + t320 ^ 2) * MDP(18) + (t201 - t363) * MDP(19);
t384 = t305 * MDP(10) + t307 * MDP(11);
t296 = t301 ^ 2;
t297 = t303 ^ 2;
t383 = (t296 + t297) * (qJ(2) * MDP(7) + MDP(6));
t280 = -pkin(2) * t303 - pkin(6) * t301 - pkin(1);
t273 = qJD(1) * t280 + qJD(2);
t270 = t307 * t273;
t359 = t303 * t305;
t365 = qJ(4) * t301;
t313 = -qJ(2) * t359 - t307 * t365;
t243 = qJD(1) * t313 + t270;
t234 = -pkin(3) * t289 + t243;
t366 = qJ(2) * t307;
t336 = t303 * t366;
t244 = -qJ(4) * t334 + qJD(1) * t336 + t273 * t305;
t360 = t302 * t244;
t211 = t300 * t234 + t360;
t370 = pkin(7) * t262;
t200 = t211 - t370;
t199 = t200 * t340;
t274 = pkin(3) * t334 + qJ(2) * t348 + qJD(4);
t242 = pkin(4) * t262 + t274;
t382 = t242 * t226 + t199;
t372 = 0.2e1 * t296;
t379 = MDP(9) * (t305 ^ 2 - t307 ^ 2);
t378 = t372 + t297;
t375 = qJD(5) + t285;
t341 = qJD(4) * t301;
t309 = qJD(3) * t313 - t305 * t341;
t338 = qJD(1) * qJD(2);
t331 = t303 * t338;
t342 = qJD(3) * t307;
t356 = t273 * t342 + t307 * t331;
t220 = qJD(1) * t309 + t356;
t346 = qJD(2) * t305;
t332 = t303 * t346;
t312 = -t307 * t341 - t332;
t343 = qJD(3) * t305;
t361 = t301 * t305;
t221 = -t273 * t343 + ((qJ(4) * t361 - t336) * qJD(3) + t312) * qJD(1);
t194 = -t220 * t300 + t302 * t221;
t192 = pkin(7) * t261 + t194;
t195 = t302 * t220 + t300 * t221;
t193 = pkin(7) * t260 + t195;
t328 = t306 * t192 - t304 * t193;
t374 = -t242 * t320 + t328;
t329 = -t280 + t365;
t373 = t305 * t329 - t336;
t371 = pkin(3) * t300;
t369 = pkin(7) * t265;
t367 = qJ(2) * t305;
t308 = qJD(1) ^ 2;
t362 = t296 * t308;
t238 = t300 * t244;
t345 = qJD(2) * t307;
t355 = t280 * t342 + t303 * t345;
t232 = t309 + t355;
t233 = qJD(3) * t373 + t312;
t206 = t302 * t232 + t300 * t233;
t214 = t302 * t243 - t238;
t248 = -t329 * t307 + (-pkin(3) - t367) * t303;
t216 = t300 * t248 - t302 * t373;
t358 = t279 * qJD(3) - t303 * t314;
t357 = t289 * t318;
t287 = t301 * pkin(3) * t342;
t354 = qJD(1) * t287 + t301 * t338;
t353 = t301 * qJD(2) + t287;
t352 = pkin(3) * t361 + t301 * qJ(2);
t339 = qJD(3) + t289;
t335 = qJ(2) * t343;
t205 = -t232 * t300 + t302 * t233;
t210 = t302 * t234 - t238;
t213 = -t243 * t300 - t360;
t215 = t302 * t248 + t300 * t373;
t326 = qJD(1) * t339;
t325 = t307 * t296 * t305 * MDP(8);
t324 = qJD(5) * t279 + t358;
t323 = -qJD(5) * t318 + t357;
t198 = -pkin(4) * t289 + t210 - t369;
t321 = -t304 * t198 - t306 * t200;
t271 = t279 * t301;
t272 = t318 * t301;
t319 = -t271 * t306 + t272 * t304;
t236 = -t271 * t304 - t272 * t306;
t293 = pkin(3) * t302 + pkin(4);
t251 = pkin(3) * t333 + pkin(4) * t265;
t249 = pkin(4) * t271 + t352;
t245 = -pkin(4) * t264 + t353;
t237 = -pkin(4) * t260 + t354;
t212 = -pkin(7) * t271 + t216;
t209 = -pkin(4) * t303 + pkin(7) * t272 + t215;
t208 = qJD(5) * t236 - t264 * t306 - t268 * t304;
t207 = qJD(5) * t319 + t264 * t304 - t268 * t306;
t204 = t214 - t369;
t203 = t213 + t370;
t197 = pkin(7) * t264 + t206;
t196 = pkin(7) * t268 + t205;
t1 = [(t289 * t332 + (-(-t280 * t305 - t336) * t289 + t273 * t359) * qJD(3) + t378 * qJD(1) * (qJ(2) * t342 + t346)) * MDP(13) + ((-t303 * t335 + t355) * t289 + t356 * t303 + (-t335 * t378 + t345 * t372) * qJD(1)) * MDP(14) + (t194 * t272 - t195 * t271 - t205 * t265 - t206 * t262 + t210 * t268 + t211 * t264 + t215 * t261 + t216 * t260) * MDP(15) + (t194 * t215 + t195 * t216 + t210 * t205 + t211 * t206 + t274 * t353 + t352 * t354) * MDP(16) + (t201 * t236 + t207 * t320) * MDP(17) + (t201 * t319 - t202 * t236 - t207 * t226 - t208 * t320) * MDP(18) + (-t201 * t303 - t207 * t285) * MDP(19) + (t202 * t303 + t208 * t285) * MDP(20) + (-(t196 * t306 - t197 * t304) * t285 - t328 * t303 + t245 * t226 + t249 * t202 - t237 * t319 + t242 * t208 + (-(-t209 * t304 - t212 * t306) * t285 - t321 * t303) * qJD(5)) * MDP(22) + (-t199 * t303 + t249 * t201 + t242 * t207 + t245 * t320 + t237 * t236 + ((-qJD(5) * t212 + t196) * t285 + t192 * t303) * t304 + ((qJD(5) * t209 + t197) * t285 + (qJD(5) * t198 + t193) * t303) * t306) * MDP(23) + (t372 * t379 - 0.2e1 * t325) * qJD(1) * qJD(3) + 0.2e1 * t338 * t383 + t384 * (t289 + t347) * t344; (t260 * t279 - t261 * t318 - t262 * t357 + t265 * t358) * MDP(15) + (-t194 * t318 + t195 * t279 - t210 * t358 + t211 * t357 - t274 * t348) * MDP(16) + (-t226 * t348 + (t304 * t323 + t306 * t324) * t285) * MDP(22) + (-t320 * t348 + (-t304 * t324 + t306 * t323) * t285) * MDP(23) - t308 * t383 + (MDP(13) * t305 + MDP(14) * t307) * (-t289 ^ 2 - t362); t308 * t325 - t362 * t379 + ((-t273 * t339 - t331) * t305 + (-t303 * t326 - t362) * t366) * MDP(13) + (-t270 * t289 + (t339 * t347 + t362) * t367 - t356) * MDP(14) + ((t211 + t213) * t265 - (t210 - t214) * t262 + (t260 * t300 + t261 * t302) * pkin(3)) * MDP(15) + (-t210 * t213 - t211 * t214 + (t194 * t302 + t195 * t300 - t274 * t333) * pkin(3)) * MDP(16) + ((t203 * t306 - t204 * t304) * t285 - t251 * t226 + (-(-t293 * t304 - t306 * t371) * t285 + t321) * qJD(5) + t374) * MDP(22) + (-t306 * t193 - t304 * t192 - (t203 * t304 + t204 * t306) * t285 - t251 * t320 + ((t293 * t306 - t304 * t371) * t285 - t306 * t198) * qJD(5) + t382) * MDP(23) - t384 * t301 * t326 + t385; (-t262 ^ 2 - t265 ^ 2) * MDP(15) + (t210 * t265 + t211 * t262 + t354) * MDP(16) + (t202 - t364) * MDP(22) + (t201 + t363) * MDP(23); (t321 * t375 + t374) * MDP(22) + ((t200 * t285 - t192) * t304 + (-t198 * t375 - t193) * t306 + t382) * MDP(23) + t385;];
tauc = t1;
