% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:58:32
% EndTime: 2019-12-31 20:58:37
% DurationCPUTime: 1.96s
% Computational Cost: add. (1797->277), mult. (4491->320), div. (0->0), fcn. (2857->4), ass. (0->123)
t297 = cos(qJ(2));
t322 = pkin(2) * t297 + pkin(1);
t273 = t322 * qJD(1);
t295 = sin(qJ(3));
t296 = sin(qJ(2));
t354 = cos(qJ(3));
t264 = t295 * t297 + t296 * t354;
t334 = qJD(1) * t264;
t363 = -qJ(4) * t334 - t273;
t292 = qJD(2) + qJD(3);
t355 = -pkin(7) - pkin(6);
t274 = t355 * t296;
t275 = t355 * t297;
t231 = t295 * t274 - t354 * t275;
t316 = t354 * qJD(3);
t361 = -t354 * qJD(2) - t316;
t360 = t297 * MDP(4) - pkin(1) * MDP(9);
t268 = qJD(1) * t274;
t350 = qJD(2) * pkin(2);
t259 = t268 + t350;
t270 = qJD(1) * t275;
t342 = t295 * t270;
t224 = t354 * t259 + t342;
t359 = qJD(4) - t224;
t358 = MDP(16) + MDP(18);
t357 = pkin(1) * t297 * MDP(10) + (t296 ^ 2 - t297 ^ 2) * MDP(5);
t320 = t354 * t297;
t312 = qJD(1) * t320;
t333 = qJD(1) * t296;
t249 = t295 * t333 - t312;
t356 = t249 ^ 2;
t248 = t334 ^ 2;
t298 = -pkin(3) - pkin(4);
t229 = t292 * t264;
t221 = t229 * qJD(1);
t352 = pkin(4) * t221;
t351 = pkin(4) * t334;
t349 = qJ(4) * t221;
t210 = pkin(3) * t249 + t363;
t347 = t210 * t334;
t285 = pkin(2) * t295 + qJ(4);
t346 = t221 * t285;
t345 = t249 * qJ(5);
t344 = t249 * t292;
t341 = t295 * t296;
t227 = t354 * t268 + t342;
t244 = t334 * qJ(5);
t209 = t244 + t227;
t279 = t316 * pkin(2) + qJD(4);
t338 = -t209 + t279;
t337 = -t227 + t279;
t222 = pkin(3) * t334 + t249 * qJ(4);
t258 = t354 * t270;
t225 = t259 * t295 - t258;
t336 = t292 * t312;
t332 = qJD(3) * t295;
t331 = t334 * MDP(11);
t204 = t244 + t224;
t330 = qJD(4) - t204;
t328 = 0.2e1 * qJD(1);
t326 = t296 * t350;
t325 = pkin(2) * t332;
t315 = t295 * t292;
t202 = -t315 * t333 + t336 + t344;
t324 = t202 * MDP(13) + (t248 - t356) * MDP(12);
t309 = t292 * t341;
t220 = qJD(1) * t309 - t336;
t323 = pkin(3) * t221 + qJ(4) * t220 - qJD(4) * t334;
t321 = qJD(2) * t355;
t319 = t292 * t332;
t313 = qJD(1) * t321;
t260 = t296 * t313;
t261 = t297 * t313;
t314 = t259 * t316 + t260 * t354 + t261 * t295 + t270 * t332;
t193 = t259 * t332 + t260 * t295 - t261 * t354 - t270 * t316;
t287 = -pkin(2) * t354 - pkin(3);
t226 = t268 * t295 - t258;
t208 = t226 + t345;
t311 = -t208 + t325;
t310 = -t226 + t325;
t205 = t225 + t345;
t230 = -t274 * t354 - t275 * t295;
t212 = pkin(2) * t333 + t222;
t308 = qJ(4) * t264 + t322;
t307 = qJ(5) * t221 + qJD(5) * t249 + t314;
t186 = qJ(5) * t220 - qJD(5) * t334 + t193;
t306 = t224 * t292 - t314;
t304 = t226 * t292 - t193;
t303 = t227 * t292 - t314;
t269 = t296 * t321;
t271 = t297 * t321;
t196 = t269 * t354 + t271 * t295 + t274 * t316 + t275 * t332;
t289 = t292 * qJD(4);
t185 = t289 + t307;
t190 = qJD(1) * t326 + t323;
t195 = t249 * t298 + qJD(5) - t363;
t302 = -t195 * t334 + t186;
t228 = t297 * t361 + t309;
t301 = -qJ(4) * t228 + qJD(4) * t264 - t326;
t197 = qJD(3) * t231 + t269 * t295 - t271 * t354;
t290 = t292 * qJ(4);
t284 = 0.2e1 * t289;
t283 = -pkin(4) + t287;
t276 = pkin(2) * t319;
t272 = t279 * t292;
t263 = -t320 + t341;
t223 = pkin(3) * t263 - t308;
t215 = t290 + t225;
t214 = qJ(5) * t263 + t231;
t213 = -qJ(5) * t264 + t230;
t211 = -pkin(3) * t292 + t359;
t207 = t263 * t298 + t308;
t203 = -t222 - t351;
t201 = t205 + t290;
t198 = -t212 - t351;
t194 = t292 * t298 - t244 + t359;
t192 = t289 + t314;
t191 = pkin(3) * t229 - t301;
t189 = qJ(5) * t228 - qJD(5) * t264 + t197;
t188 = qJ(5) * t229 + qJD(5) * t263 + t196;
t187 = t229 * t298 + t301;
t184 = -t190 - t352;
t1 = [(-t220 * t264 - t228 * t334) * MDP(11) + (t220 * t263 - t221 * t264 + t228 * t249 - t229 * t334) * MDP(12) + (-t221 * t322 - t229 * t273) * MDP(16) + (t220 * t322 + t228 * t273) * MDP(17) + (t190 * t263 + t191 * t249 + t210 * t229 + t221 * t223) * MDP(18) + (-t192 * t263 + t193 * t264 - t196 * t249 + t197 * t334 - t211 * t228 - t215 * t229 - t220 * t230 - t221 * t231) * MDP(19) + (-t190 * t264 - t191 * t334 + t210 * t228 + t220 * t223) * MDP(20) + (t190 * t223 + t191 * t210 + t192 * t231 + t193 * t230 + t196 * t215 + t197 * t211) * MDP(21) + (-t184 * t263 - t187 * t249 - t195 * t229 - t207 * t221) * MDP(22) + (t184 * t264 + t187 * t334 - t195 * t228 - t207 * t220) * MDP(23) + (t185 * t263 - t186 * t264 + t188 * t249 - t189 * t334 + t194 * t228 + t201 * t229 + t213 * t220 + t214 * t221) * MDP(24) + (t184 * t207 + t185 * t214 + t186 * t213 + t187 * t195 + t188 * t201 + t189 * t194) * MDP(25) + (-t228 * MDP(13) - t229 * MDP(14) - t189 * MDP(22) + t188 * MDP(23) - t358 * t197 + (-MDP(17) + MDP(20)) * t196) * t292 + (-t357 * t328 + (t360 * t328 + ((qJD(1) * t263 + t249) * MDP(16) + 0.2e1 * t334 * MDP(17)) * pkin(2)) * t296 + (t297 * MDP(6) - t296 * MDP(7) + (MDP(10) * t296 - MDP(9) * t297) * pkin(6)) * qJD(2)) * qJD(2); t249 * t331 + (t334 * t273 + (-t249 * t333 - t319) * pkin(2) + t304) * MDP(16) + (-t273 * t249 + (-t292 * t316 - t333 * t334) * pkin(2) + t303) * MDP(17) + (-t212 * t249 - t276 + t304 - t347) * MDP(18) + (-t220 * t287 - t346 + (t215 + t310) * t334 + (t211 - t337) * t249) * MDP(19) + (-t210 * t249 + t212 * t334 + t272 + t289 - t303) * MDP(20) + (t192 * t285 + t193 * t287 - t210 * t212 + t211 * t310 + t215 * t337) * MDP(21) + (t198 * t249 + t208 * t292 - t276 - t302) * MDP(22) + (t195 * t249 - t198 * t334 - t209 * t292 + t185 + t272) * MDP(23) + (t220 * t283 + t346 + (-t201 - t311) * t334 + (-t194 + t338) * t249) * MDP(24) + (t185 * t285 + t186 * t283 + t194 * t311 - t195 * t198 + t201 * t338) * MDP(25) + t324 + (-t296 * t360 + t357) * qJD(1) ^ 2; t306 * MDP(17) + (pkin(3) * t220 - t349) * MDP(19) + (t284 - t306) * MDP(20) + (-pkin(3) * t193 + qJ(4) * t192 - t210 * t222 - t211 * t225 + t215 * t359) * MDP(21) + (t205 * t292 - t186) * MDP(22) + (-t204 * t292 + t284 + t307) * MDP(23) + (t220 * t298 + t349) * MDP(24) + (qJ(4) * t185 + t186 * t298 - t194 * t205 - t195 * t203 + t201 * t330) * MDP(25) + (t273 * MDP(16) - t210 * MDP(18) + (t215 - t225) * MDP(19) + t222 * MDP(20) + t195 * MDP(22) - t203 * MDP(23) + (-t201 + t205) * MDP(24)) * t334 + (t331 - t273 * MDP(17) - t222 * MDP(18) + (t211 - t359) * MDP(19) - t210 * MDP(20) + t203 * MDP(22) + t195 * MDP(23) + (-t194 + t330) * MDP(24)) * t249 + t324 + t358 * (t225 * t292 - t193); t202 * MDP(19) + (-t215 * t292 + t193 + t347) * MDP(21) + t220 * MDP(24) + (-t201 * t292 + t302) * MDP(25) + (-t292 * MDP(24) + (MDP(18) + MDP(22)) * t334) * t249 + (MDP(20) + MDP(23)) * (-t292 ^ 2 - t248); -MDP(22) * t292 * t334 + (t336 - t344) * MDP(23) + (-t248 - t356) * MDP(24) + (t194 * t334 - t201 * t249 - t323 - t352) * MDP(25) + (-t297 * MDP(22) * t315 + (MDP(22) * t361 - MDP(23) * t315 - MDP(25) * t350) * t296) * qJD(1);];
tauc = t1;
