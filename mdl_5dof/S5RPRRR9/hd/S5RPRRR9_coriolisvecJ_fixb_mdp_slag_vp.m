% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:08:06
% EndTime: 2019-12-31 19:08:13
% DurationCPUTime: 3.49s
% Computational Cost: add. (2990->282), mult. (7996->382), div. (0->0), fcn. (6378->8), ass. (0->130)
t295 = cos(qJ(5));
t336 = qJD(5) * t295;
t290 = sin(pkin(9));
t294 = sin(qJ(3));
t291 = cos(pkin(9));
t297 = cos(qJ(3));
t349 = t291 * t297;
t308 = t290 * t294 - t349;
t268 = t308 * qJD(1);
t350 = t290 * t297;
t274 = t291 * t294 + t350;
t269 = t274 * qJD(1);
t293 = sin(qJ(4));
t296 = cos(qJ(4));
t249 = t296 * t268 + t269 * t293;
t381 = t249 * t295;
t386 = t336 + t381;
t312 = -t268 * t293 + t296 * t269;
t292 = sin(qJ(5));
t337 = qJD(5) * t292;
t340 = qJD(3) * t297;
t280 = t291 * qJD(1) * t340;
t341 = qJD(1) * t294;
t329 = t290 * t341;
t264 = -qJD(3) * t329 + t280;
t271 = t274 * qJD(3);
t265 = qJD(1) * t271;
t338 = qJD(4) * t296;
t339 = qJD(4) * t293;
t224 = t296 * t264 - t293 * t265 - t268 * t338 - t269 * t339;
t289 = qJD(3) + qJD(4);
t344 = t295 * t224 + t289 * t336;
t206 = -t312 * t337 + t344;
t205 = t206 * t295;
t243 = t289 * t292 + t295 * t312;
t362 = t224 * t292;
t207 = t243 * qJD(5) + t362;
t352 = t312 * t292;
t241 = -t295 * t289 + t352;
t385 = -t292 * t207 - t386 * t241 + t205;
t204 = t206 * t292;
t225 = qJD(4) * t312 + t264 * t293 + t296 * t265;
t221 = t292 * t225;
t379 = -qJD(5) - t249;
t345 = -t336 * t379 + t221;
t348 = t295 * t379;
t354 = t249 * t289;
t356 = t312 * t289;
t358 = t243 * t312;
t384 = (-t225 + t356) * MDP(18) - t249 ^ 2 * MDP(16) + (-t249 * t348 + t345 - t358) * MDP(24) + (MDP(15) * t249 + t312 * MDP(16) + t379 * MDP(26)) * t312 + (t224 + t354) * MDP(17) + (t386 * t243 + t204) * MDP(22);
t366 = pkin(6) + qJ(2);
t278 = t366 * t290;
t275 = qJD(1) * t278;
t279 = t366 * t291;
t276 = qJD(1) * t279;
t373 = -t297 * t275 - t276 * t294;
t239 = -pkin(7) * t269 + t373;
t238 = qJD(3) * pkin(3) + t239;
t310 = t275 * t294 - t276 * t297;
t240 = -pkin(7) * t268 - t310;
t361 = t240 * t293;
t212 = t238 * t296 - t361;
t210 = -pkin(4) * t289 - t212;
t383 = t210 * t249;
t284 = -pkin(2) * t291 - pkin(1);
t277 = qJD(1) * t284 + qJD(2);
t255 = pkin(3) * t268 + t277;
t382 = t249 * t255;
t322 = t379 * t292;
t380 = (qJ(2) * MDP(7) + MDP(6)) * (t290 ^ 2 + t291 ^ 2);
t227 = pkin(4) * t312 + pkin(8) * t249;
t359 = t241 * t312;
t223 = t295 * t225;
t372 = -t337 * t379 - t223;
t360 = t240 * t296;
t213 = t238 * t293 + t360;
t231 = -pkin(7) * t265 - qJD(2) * t268 + qJD(3) * t373;
t304 = t274 * qJD(2);
t302 = qJD(1) * t304;
t232 = -pkin(7) * t264 + qJD(3) * t310 - t302;
t325 = t231 * t293 - t296 * t232;
t195 = qJD(4) * t213 + t325;
t211 = pkin(8) * t289 + t213;
t214 = pkin(4) * t249 - pkin(8) * t312 + t255;
t314 = t211 * t292 - t214 * t295;
t371 = -t195 * t295 + t210 * t337 + t312 * t314;
t197 = t211 * t295 + t214 * t292;
t370 = t195 * t292 + t197 * t312 + t210 * t336;
t369 = -t312 * t255 - t325;
t324 = t232 * t293 - t240 * t339;
t194 = (qJD(4) * t238 + t231) * t296 + t324;
t301 = -t278 * t340 + qJD(2) * t349 + (-qJD(2) * t290 - qJD(3) * t279) * t294;
t233 = -pkin(7) * t271 + t301;
t270 = t308 * qJD(3);
t309 = t278 * t294 - t279 * t297;
t299 = qJD(3) * t309 - t304;
t234 = pkin(7) * t270 + t299;
t245 = -pkin(7) * t274 - t278 * t297 - t279 * t294;
t246 = -pkin(7) * t308 - t309;
t313 = t245 * t296 - t246 * t293;
t198 = qJD(4) * t313 + t233 * t296 + t234 * t293;
t219 = t245 * t293 + t246 * t296;
t254 = t274 * t296 - t293 * t308;
t259 = pkin(3) * t308 + t284;
t311 = -t274 * t293 - t296 * t308;
t220 = -pkin(4) * t311 - pkin(8) * t254 + t259;
t229 = qJD(4) * t311 - t270 * t296 - t271 * t293;
t368 = t195 * t254 + t210 * t229 - t219 * t225 + (qJD(5) * t220 + t198) * t379 + (qJD(5) * t214 + t194) * t311;
t367 = pkin(3) * t269;
t364 = t210 * t254;
t363 = t220 * t225;
t357 = t243 * t292;
t331 = qJD(1) * qJD(2);
t327 = -pkin(3) * t289 - t238;
t285 = pkin(3) * t293 + pkin(8);
t318 = qJD(5) * t285 + t227 + t367;
t215 = t239 * t293 + t360;
t317 = pkin(3) * t339 - t215;
t216 = t239 * t296 - t361;
t316 = -pkin(3) * t338 + t216;
t315 = -t225 * t285 + t383;
t306 = t249 * t322 - t372;
t305 = t229 * t295 - t254 * t337;
t286 = -pkin(3) * t296 - pkin(4);
t230 = qJD(4) * t254 - t270 * t293 + t296 * t271;
t202 = pkin(3) * t271 + pkin(4) * t230 - pkin(8) * t229;
t201 = pkin(3) * t265 + pkin(4) * t225 - pkin(8) * t224;
t200 = t295 * t201;
t199 = qJD(4) * t219 + t233 * t293 - t234 * t296;
t1 = [(t264 * t274 - t269 * t270) * MDP(8) + (-t264 * t308 - t265 * t274 + t268 * t270 - t269 * t271) * MDP(9) + (t284 * t265 + t277 * t271) * MDP(13) + (t284 * t264 - t277 * t270) * MDP(14) + (t224 * t254 + t229 * t312) * MDP(15) + (t224 * t311 - t225 * t254 - t229 * t249 - t230 * t312) * MDP(16) + (t225 * t259 + t230 * t255 + (t249 * t271 - t265 * t311) * pkin(3)) * MDP(20) + (t224 * t259 + t229 * t255 + (t254 * t265 + t271 * t312) * pkin(3)) * MDP(21) + (t205 * t254 + t243 * t305) * MDP(22) + ((-t241 * t295 - t357) * t229 + (-t204 - t207 * t295 + (t241 * t292 - t243 * t295) * qJD(5)) * t254) * MDP(23) + (-t206 * t311 + t223 * t254 + t230 * t243 - t305 * t379) * MDP(24) + (-t254 * t221 + t207 * t311 - t230 * t241 - (-t229 * t292 - t254 * t336) * t379) * MDP(25) + (-t225 * t311 - t230 * t379) * MDP(26) + (-t314 * t230 + t199 * t241 - t200 * t311 - t313 * t207 + (-t202 * t379 + t363 + (t211 * t311 + t219 * t379 + t364) * qJD(5)) * t295 + t368 * t292) * MDP(27) + (-t197 * t230 + t199 * t243 - t313 * t206 + ((-qJD(5) * t219 + t202) * t379 - t363 + (-qJD(5) * t211 + t201) * t311 - qJD(5) * t364) * t292 + t368 * t295) * MDP(28) + 0.2e1 * t331 * t380 + (MDP(17) * t229 - MDP(18) * t230 - MDP(20) * t199 - MDP(21) * t198) * t289 + (-t270 * MDP(10) - t271 * MDP(11) + t299 * MDP(13) - t301 * MDP(14)) * qJD(3); t280 * MDP(14) + (t225 + t356) * MDP(20) + (t224 - t354) * MDP(21) + (t306 - t359) * MDP(27) + (-t348 * t379 - t221 - t358) * MDP(28) + ((qJD(1) * t350 + t291 * t341 + t269) * MDP(13) + (-t268 - t329) * MDP(14)) * qJD(3) - qJD(1) ^ 2 * t380; t269 * t268 * MDP(8) + (-t268 ^ 2 + t269 ^ 2) * MDP(9) + (t280 + (t268 - t329) * qJD(3)) * MDP(10) + (-t277 * t269 - t302) * MDP(13) + (t277 * t268 + t308 * t331) * MDP(14) + (-t249 * t367 + t215 * t289 + (t293 * t327 - t360) * qJD(4) + t369) * MDP(20) + (-t312 * t367 + t216 * t289 + t382 + (qJD(4) * t327 - t231) * t296 - t324) * MDP(21) + (t357 * t379 + t385) * MDP(23) + (t306 + t359) * MDP(25) + (t286 * t207 + t315 * t292 + t317 * t241 - (t292 * t316 - t295 * t318) * t379 + t371) * MDP(27) + (t286 * t206 + t315 * t295 + t317 * t243 - (t292 * t318 + t295 * t316) * t379 + t370) * MDP(28) + t384; ((-qJD(4) + t289) * t213 + t369) * MDP(20) + (t212 * t289 - t194 + t382) * MDP(21) + (t243 * t322 + t385) * MDP(23) + (-t322 * t379 + t223 + t359) * MDP(25) + (-pkin(4) * t207 + (-t212 * t292 + t227 * t295) * t379 - t213 * t241 + t292 * t383 - t345 * pkin(8) + t371) * MDP(27) + (-pkin(4) * t206 - (t212 * t295 + t227 * t292) * t379 - t213 * t243 + t210 * t381 + t372 * pkin(8) + t370) * MDP(28) + t384; t243 * t241 * MDP(22) + (-t241 ^ 2 + t243 ^ 2) * MDP(23) + (-t241 * t379 + t344) * MDP(24) + (-t243 * t379 - t362) * MDP(25) + t225 * MDP(26) + (-t194 * t292 - t197 * t379 - t210 * t243 + t200) * MDP(27) + (-t194 * t295 - t201 * t292 + t210 * t241 + t314 * t379) * MDP(28) + (-MDP(24) * t352 - MDP(25) * t243 - MDP(27) * t197 + MDP(28) * t314) * qJD(5);];
tauc = t1;
