% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RPRRR2
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
%   see S5RPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RPRRR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:12:23
% EndTime: 2019-12-05 18:12:32
% DurationCPUTime: 3.65s
% Computational Cost: add. (2728->246), mult. (7484->331), div. (0->0), fcn. (6028->8), ass. (0->125)
t289 = cos(pkin(9));
t295 = cos(qJ(3));
t336 = qJD(3) * t295;
t280 = t289 * qJD(1) * t336;
t288 = sin(pkin(9));
t292 = sin(qJ(3));
t337 = qJD(1) * t292;
t327 = t288 * t337;
t264 = -qJD(3) * t327 + t280;
t343 = t288 * t295;
t274 = t289 * t292 + t343;
t271 = t274 * qJD(3);
t265 = qJD(1) * t271;
t342 = t289 * t295;
t308 = t288 * t292 - t342;
t268 = t308 * qJD(1);
t269 = t274 * qJD(1);
t291 = sin(qJ(4));
t294 = cos(qJ(4));
t334 = qJD(4) * t294;
t335 = qJD(4) * t291;
t208 = t264 * t294 - t265 * t291 - t268 * t334 - t269 * t335;
t312 = -t268 * t291 + t269 * t294;
t209 = qJD(4) * t312 + t264 * t291 + t265 * t294;
t249 = t268 * t294 + t269 * t291;
t287 = qJD(3) + qJD(4);
t345 = t249 * t287;
t346 = t312 * t287;
t293 = cos(qJ(5));
t244 = t293 * t249;
t290 = sin(qJ(5));
t333 = qJD(5) * t290;
t191 = -qJD(5) * t244 + t208 * t293 - t209 * t290 - t312 * t333;
t215 = -t290 * t312 - t244;
t313 = t249 * t290 - t293 * t312;
t299 = qJD(5) * t313 - t208 * t290 - t209 * t293;
t328 = -qJD(4) - qJD(5);
t284 = qJD(3) - t328;
t347 = t215 * t284;
t348 = t284 * t313;
t363 = (t191 - t347) * MDP(24) + t215 * t313 * MDP(22) + (t299 - t348) * MDP(25) + (-t215 ^ 2 + t313 ^ 2) * MDP(23);
t374 = t363 + t249 * t312 * MDP(15) + (-t209 + t346) * MDP(18) + (-t249 ^ 2 + t312 ^ 2) * MDP(16) + (t208 + t345) * MDP(17);
t350 = pkin(6) + qJ(2);
t278 = t350 * t288;
t275 = qJD(1) * t278;
t279 = t350 * t289;
t276 = qJD(1) * t279;
t310 = t275 * t292 - t276 * t295;
t240 = -pkin(7) * t268 - t310;
t237 = t294 * t240;
t356 = -t275 * t295 - t276 * t292;
t239 = -pkin(7) * t269 + t356;
t238 = qJD(3) * pkin(3) + t239;
t315 = -t238 * t291 - t237;
t369 = pkin(8) * t249;
t198 = -t315 - t369;
t282 = -pkin(2) * t289 - pkin(1);
t277 = qJD(1) * t282 + qJD(2);
t255 = pkin(3) * t268 + t277;
t225 = pkin(4) * t249 + t255;
t373 = t198 * t333 - t225 * t215;
t223 = -pkin(7) * t265 - qJD(2) * t268 + qJD(3) * t356;
t304 = t274 * qJD(2);
t303 = qJD(1) * t304;
t224 = -pkin(7) * t264 + qJD(3) * t310 - t303;
t321 = t291 * t224 - t240 * t335;
t354 = (qJD(4) * t238 + t223) * t294 + t321;
t187 = -pkin(8) * t209 + t354;
t322 = -t291 * t223 + t224 * t294;
t301 = qJD(4) * t315 + t322;
t188 = -pkin(8) * t208 + t301;
t365 = -t290 * t187 + t188 * t293 + t225 * t313;
t368 = t255 * t249;
t367 = (MDP(7) * qJ(2) + MDP(6)) * (t288 ^ 2 + t289 ^ 2);
t366 = (-t198 * t284 - t188) * t290 + t373;
t351 = pkin(4) * t312;
t361 = pkin(8) * t312;
t360 = t255 * t312;
t355 = qJD(5) - t284;
t353 = pkin(3) * t269;
t352 = pkin(3) * t284;
t235 = t291 * t240;
t341 = t291 * t293;
t340 = t293 * t198;
t339 = t239 * t294 - t235;
t329 = qJD(1) * qJD(2);
t326 = -pkin(3) * t287 - t238;
t320 = t238 * t294 - t235;
t197 = t320 - t361;
t195 = pkin(4) * t287 + t197;
t325 = -pkin(4) * t284 - t195;
t319 = -t239 * t291 - t237;
t316 = -t290 * t195 - t340;
t246 = -pkin(7) * t274 - t278 * t295 - t279 * t292;
t309 = t278 * t292 - t279 * t295;
t247 = -pkin(7) * t308 - t309;
t314 = -t246 * t291 - t247 * t294;
t254 = t274 * t294 - t291 * t308;
t311 = -t274 * t291 - t294 * t308;
t221 = t254 * t290 - t293 * t311;
t222 = t254 * t293 + t290 * t311;
t259 = pkin(3) * t308 + t282;
t302 = -t278 * t336 + qJD(2) * t342 + (-qJD(2) * t288 - qJD(3) * t279) * t292;
t229 = -pkin(7) * t271 + t302;
t270 = t308 * qJD(3);
t297 = qJD(3) * t309 - t304;
t230 = pkin(7) * t270 + t297;
t305 = t229 * t294 + t230 * t291 + t246 * t334 - t247 * t335;
t300 = qJD(4) * t314 - t229 * t291 + t230 * t294;
t283 = pkin(3) * t294 + pkin(4);
t233 = -pkin(4) * t311 + t259;
t231 = t353 + t351;
t220 = qJD(4) * t254 - t270 * t291 + t271 * t294;
t219 = qJD(4) * t311 - t270 * t294 - t271 * t291;
t204 = pkin(3) * t271 + pkin(4) * t220;
t203 = pkin(3) * t265 + pkin(4) * t209;
t202 = pkin(8) * t311 - t314;
t201 = -pkin(8) * t254 + t246 * t294 - t247 * t291;
t200 = t339 - t361;
t199 = t319 + t369;
t194 = qJD(5) * t222 + t219 * t290 + t220 * t293;
t193 = -qJD(5) * t221 + t219 * t293 - t220 * t290;
t190 = -pkin(8) * t219 + t300;
t189 = -pkin(8) * t220 + t305;
t1 = [(t264 * t274 - t269 * t270) * MDP(8) + (-t264 * t308 - t265 * t274 + t268 * t270 - t269 * t271) * MDP(9) + (t282 * t265 + t277 * t271) * MDP(13) + (t282 * t264 - t277 * t270) * MDP(14) + (t208 * t254 + t219 * t312) * MDP(15) + (t208 * t311 - t209 * t254 - t219 * t249 - t220 * t312) * MDP(16) + (t259 * t209 + t255 * t220 + (t249 * t271 - t265 * t311) * pkin(3)) * MDP(20) + (t259 * t208 + t255 * t219 + (t254 * t265 + t271 * t312) * pkin(3)) * MDP(21) + (t191 * t222 - t193 * t313) * MDP(22) + (-t191 * t221 + t193 * t215 + t194 * t313 + t222 * t299) * MDP(23) + (t225 * t194 + t203 * t221 - t204 * t215 - t233 * t299) * MDP(27) + (t233 * t191 + t225 * t193 + t203 * t222 - t204 * t313) * MDP(28) + 0.2e1 * t329 * t367 + (MDP(17) * t219 - MDP(18) * t220 + MDP(20) * t300 - MDP(21) * t305) * t287 + (t193 * MDP(24) - t194 * MDP(25) + (-t189 * t290 + t190 * t293 + (-t201 * t290 - t202 * t293) * qJD(5)) * MDP(27) + (-t189 * t293 - t190 * t290 - (t201 * t293 - t202 * t290) * qJD(5)) * MDP(28)) * t284 + (-t270 * MDP(10) - t271 * MDP(11) + MDP(13) * t297 - MDP(14) * t302) * qJD(3); t280 * MDP(14) + (t209 + t346) * MDP(20) + (t208 - t345) * MDP(21) + (-t299 - t348) * MDP(27) + (t191 + t347) * MDP(28) + ((qJD(1) * t343 + t289 * t337 + t269) * MDP(13) + (-t268 - t327) * MDP(14)) * qJD(3) - qJD(1) ^ 2 * t367; t269 * t268 * MDP(8) + (-t268 ^ 2 + t269 ^ 2) * MDP(9) + (t280 + (t268 - t327) * qJD(3)) * MDP(10) + (-t277 * t269 - t303) * MDP(13) + (t277 * t268 + t308 * t329) * MDP(14) + (-t249 * t353 - t360 - t319 * t287 + (t291 * t326 - t237) * qJD(4) + t322) * MDP(20) + (-t312 * t353 + t368 + t339 * t287 + (qJD(4) * t326 - t223) * t294 - t321) * MDP(21) + (t231 * t215 - (t199 * t293 - t200 * t290) * t284 + (-t290 * t294 - t341) * qJD(4) * t352 + ((-pkin(3) * t341 - t283 * t290) * t284 + t316) * qJD(5) + t365) * MDP(27) + (t231 * t313 + (-t291 * t328 * t352 + t199 * t284 - t188) * t290 + (-qJD(5) * t195 - t187 + (-pkin(3) * t334 - qJD(5) * t283 + t200) * t284) * t293 + t373) * MDP(28) + t374; (-t287 * t315 + t301 - t360) * MDP(20) + (t287 * t320 - t354 + t368) * MDP(21) + (t215 * t351 - (-t197 * t290 - t340) * t284 + (t290 * t325 - t340) * qJD(5) + t365) * MDP(27) + (t313 * t351 + (qJD(5) * t325 + t197 * t284 - t187) * t293 + t366) * MDP(28) + t374; (t316 * t355 + t365) * MDP(27) + ((-t195 * t355 - t187) * t293 + t366) * MDP(28) + t363;];
tauc = t1;
