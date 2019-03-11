% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RPPRPR4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:47:12
% EndTime: 2019-03-09 01:47:18
% DurationCPUTime: 2.09s
% Computational Cost: add. (1988->261), mult. (4123->369), div. (0->0), fcn. (2705->8), ass. (0->132)
t274 = cos(qJ(6));
t323 = t274 * qJD(4);
t268 = sin(pkin(10));
t273 = sin(qJ(4));
t275 = cos(qJ(4));
t352 = cos(pkin(10));
t241 = t268 * t275 + t352 * t273;
t237 = t241 * qJD(1);
t272 = sin(qJ(6));
t344 = t237 * t272;
t215 = -t323 - t344;
t313 = t352 * t275;
t249 = qJD(1) * t313;
t329 = qJD(1) * t273;
t235 = t268 * t329 - t249;
t322 = qJD(6) - t235;
t364 = t322 * t215;
t217 = qJD(4) * t272 - t237 * t274;
t363 = t322 * t217;
t309 = t322 * t274;
t236 = t241 * qJD(4);
t230 = qJD(1) * t236;
t341 = t272 * t230;
t362 = t322 * t309 - t341;
t269 = sin(pkin(9));
t277 = qJD(4) ^ 2;
t278 = qJD(1) ^ 2;
t360 = t269 * (t277 + t278);
t276 = -pkin(1) - pkin(2);
t247 = t276 * qJD(1) + qJD(2);
t270 = cos(pkin(9));
t330 = qJD(1) * t270;
t233 = qJ(2) * t330 + t269 * t247;
t226 = -qJD(1) * pkin(7) + t233;
t359 = t275 * qJD(3) - t226 * t273;
t328 = qJD(2) * t270;
t303 = -qJD(5) + t328;
t293 = t303 * t273;
t327 = qJD(3) * t273;
t358 = (-t226 * t275 - t327) * qJD(4) + (qJD(4) * t275 * qJ(5) - t293) * qJD(1);
t326 = qJD(4) * t273;
t238 = qJD(4) * t313 - t268 * t326;
t357 = MDP(15) * t273 + MDP(16) * t275;
t356 = (t273 ^ 2 - t275 ^ 2) * MDP(11);
t355 = qJ(2) * MDP(6) + t270 * MDP(8) + MDP(5);
t292 = t303 * t275;
t196 = t359 * qJD(4) + (qJ(5) * t326 + t292) * qJD(1);
t183 = t196 * t268 - t352 * t358;
t254 = pkin(4) * t268 + pkin(8);
t354 = (-pkin(4) * t329 - pkin(5) * t237 - pkin(8) * t235 + qJD(6) * t254) * t322 + t183;
t184 = t352 * t196 + t358 * t268;
t212 = qJ(5) * t329 + t359;
t209 = qJD(4) * pkin(4) + t212;
t213 = t327 + (-qJ(5) * qJD(1) + t226) * t275;
t343 = t268 * t213;
t187 = t352 * t209 - t343;
t185 = -qJD(4) * pkin(5) - t187;
t336 = t270 * qJ(2) + t269 * t276;
t243 = -pkin(7) + t336;
t340 = qJ(5) - t243;
t310 = qJD(4) * t340;
t211 = t273 * t310 + t292;
t280 = t275 * t310 - t293;
t191 = t352 * t211 + t268 * t280;
t331 = qJD(1) * t269;
t232 = -qJ(2) * t331 + t247 * t270;
t225 = qJD(1) * pkin(3) - t232;
t265 = t275 * pkin(4);
t214 = qJD(1) * t265 + qJD(5) + t225;
t193 = -pkin(5) * t235 + pkin(8) * t237 + t214;
t285 = -t268 * t273 + t313;
t311 = -t269 * qJ(2) + t270 * t276;
t242 = pkin(3) - t311;
t294 = t242 + t265;
t203 = pkin(5) * t285 + pkin(8) * t241 + t294;
t227 = t340 * t275;
t312 = t340 * t273;
t199 = -t352 * t227 + t268 * t312;
t299 = -t183 * t241 + t199 * t230;
t353 = -t185 * t238 - (qJD(6) * t203 + t191) * t322 - (qJD(6) * t193 + t184) * t285 + t299;
t324 = qJD(6) * t272;
t320 = qJD(1) * qJD(4);
t314 = t273 * t320;
t231 = qJD(4) * t249 - t268 * t314;
t337 = qJD(6) * t323 - t274 * t231;
t200 = t237 * t324 + t337;
t350 = t200 * t272;
t349 = t203 * t230;
t348 = t215 * t237;
t347 = t217 * t237;
t345 = t231 * t272;
t219 = t274 * t230;
t207 = t352 * t213;
t188 = t268 * t209 + t207;
t339 = -t270 * t237 + t238 * t269;
t338 = -t269 * t236 - t285 * t330;
t325 = qJD(6) * t241;
t321 = qJD(1) * qJD(2);
t318 = 0.2e1 * t320;
t317 = t241 * t341;
t316 = t241 * t219;
t250 = t269 * t321;
t305 = t275 * t318;
t304 = 0.2e1 * t250;
t302 = -pkin(4) * t326 + t269 * qJD(2);
t300 = -qJD(6) * t270 + t338;
t186 = qJD(4) * pkin(8) + t188;
t182 = t186 * t274 + t193 * t272;
t298 = t186 * t272 - t193 * t274;
t297 = t200 * t285 - t217 * t236;
t201 = t217 * qJD(6) - t345;
t296 = -t201 * t285 + t215 * t236;
t295 = t232 * t269 - t233 * t270;
t229 = t285 * t269;
t290 = qJD(6) * t229 + t331;
t239 = -pkin(4) * t314 + t250;
t289 = -t219 + (t235 * t272 - t324) * t322;
t288 = t238 * t272 + t274 * t325;
t287 = -t238 * t274 + t241 * t324;
t286 = -t243 * t277 + t304;
t282 = qJD(4) * (-qJD(1) * t242 - t225 - t328);
t192 = t352 * t212 - t343;
t281 = t254 * t230 + (t185 + t192) * t322;
t255 = -t352 * pkin(4) - pkin(5);
t228 = t241 * t269;
t202 = -pkin(5) * t236 + pkin(8) * t238 + t302;
t198 = -t227 * t268 - t352 * t312;
t197 = -pkin(5) * t230 + pkin(8) * t231 + t239;
t195 = t274 * t197;
t190 = t212 * t268 + t207;
t189 = t211 * t268 - t352 * t280;
t1 = [MDP(7) * t304 + ((-t269 * t311 + t270 * t336) * qJD(1) - t295) * qJD(2) * MDP(9) - t318 * t356 + (-t184 * t285 + t187 * t238 + t188 * t236 - t189 * t237 + t191 * t235 - t198 * t231 + t299) * MDP(17) + (t183 * t198 + t184 * t199 - t187 * t189 + t188 * t191 + t214 * t302 + t239 * t294) * MDP(18) + (-t200 * t241 * t274 + t287 * t217) * MDP(19) + ((t215 * t274 + t217 * t272) * t238 + (t350 + t201 * t274 + (-t215 * t272 + t217 * t274) * qJD(6)) * t241) * MDP(20) + (t287 * t322 + t297 + t316) * MDP(21) + (t288 * t322 + t296 - t317) * MDP(22) + (-t230 * t285 - t236 * t322) * MDP(23) + (t298 * t236 + t189 * t215 + t195 * t285 + t198 * t201 + (t202 * t322 - t349 + (-t185 * t241 - t186 * t285 - t199 * t322) * qJD(6)) * t274 + t353 * t272) * MDP(24) + (t182 * t236 + t189 * t217 + t198 * t200 + (-(-qJD(6) * t199 + t202) * t322 + t349 - (-qJD(6) * t186 + t197) * t285 + t185 * t325) * t272 + t353 * t274) * MDP(25) + (-t277 * MDP(12) + t286 * MDP(15) + t282 * MDP(16)) * t275 + (MDP(10) * t305 + t277 * MDP(13) + t282 * MDP(15) - t286 * MDP(16)) * t273 + 0.2e1 * t355 * t321; t295 * qJD(1) * MDP(9) + (0.2e1 * t270 * t314 - t275 * t360) * MDP(15) + (t270 * t305 + t273 * t360) * MDP(16) + (-t228 * t231 + t229 * t230 + t338 * t235 - t339 * t237) * MDP(17) + (t183 * t228 + t184 * t229 - t339 * t187 + t338 * t188 - t214 * t331 - t239 * t270) * MDP(18) + (-(-t229 * t272 - t270 * t274) * t230 + t228 * t201 - (t300 * t272 + t290 * t274) * t322 + t339 * t215) * MDP(24) + ((t229 * t274 - t270 * t272) * t230 + t228 * t200 - (-t290 * t272 + t300 * t274) * t322 + t339 * t217) * MDP(25) + (-t269 * MDP(7) - t355) * t278; (t230 * t241 + t231 * t285 + t235 * t238 - t236 * t237) * MDP(17) + (-t183 * t285 + t184 * t241 - t187 * t236 + t188 * t238) * MDP(18) + (t296 + t317) * MDP(24) + (-t297 + t316) * MDP(25) - t357 * t277 - (t288 * MDP(24) - t287 * MDP(25)) * t322; ((-t188 + t190) * t237 + (t187 - t192) * t235 + (t230 * t268 + t352 * t231) * pkin(4)) * MDP(17) + (t187 * t190 - t188 * t192 + (-t352 * t183 + t184 * t268 + t214 * t329) * pkin(4)) * MDP(18) + (t217 * t309 + t350) * MDP(19) + ((t200 - t364) * t274 + (-t201 - t363) * t272) * MDP(20) + (t347 + t362) * MDP(21) + (t289 - t348) * MDP(22) + t322 * t237 * MDP(23) + (-t190 * t215 + t255 * t201 - t237 * t298 + t281 * t272 - t354 * t274) * MDP(24) + (-t182 * t237 - t190 * t217 + t255 * t200 + t354 * t272 + t281 * t274) * MDP(25) + t357 * (t225 - t328) * qJD(1) + (-t273 * t275 * MDP(10) + t356) * t278; (-t235 ^ 2 - t237 ^ 2) * MDP(17) + (-t187 * t237 - t188 * t235 + t239) * MDP(18) + (t289 + t348) * MDP(24) + (t347 - t362) * MDP(25); t217 * t215 * MDP(19) + (-t215 ^ 2 + t217 ^ 2) * MDP(20) + (t337 + t364) * MDP(21) + (t345 + t363) * MDP(22) - t230 * MDP(23) + (t182 * t322 - t184 * t272 - t185 * t217 + t195) * MDP(24) + (-t184 * t274 + t185 * t215 - t197 * t272 - t298 * t322) * MDP(25) + (MDP(21) * t344 - t217 * MDP(22) - t182 * MDP(24) + t298 * MDP(25)) * qJD(6);];
tauc  = t1;
