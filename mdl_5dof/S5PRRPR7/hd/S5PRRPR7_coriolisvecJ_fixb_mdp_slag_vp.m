% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR7_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR7_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:42
% EndTime: 2019-12-05 16:37:50
% DurationCPUTime: 2.98s
% Computational Cost: add. (1458->315), mult. (3919->481), div. (0->0), fcn. (2861->10), ass. (0->144)
t287 = sin(qJ(3));
t371 = MDP(5) * t287;
t290 = cos(qJ(3));
t370 = MDP(6) * (t287 ^ 2 - t290 ^ 2);
t307 = pkin(3) * t287 - qJ(4) * t290;
t244 = qJD(3) * t307 - qJD(4) * t287;
t284 = cos(pkin(10));
t291 = cos(qJ(2));
t283 = sin(pkin(5));
t344 = qJD(1) * t283;
t327 = t291 * t344;
t288 = sin(qJ(2));
t328 = t288 * t344;
t282 = sin(pkin(10));
t358 = t282 * t290;
t369 = t327 * t358 + (t244 - t328) * t284;
t354 = t284 * t291;
t368 = t282 * t244 - (t282 * t288 + t290 * t354) * t344;
t262 = qJD(2) * pkin(7) + t328;
t285 = cos(pkin(5));
t343 = qJD(1) * t285;
t367 = -t287 * t262 + t290 * t343;
t366 = pkin(7) * t284;
t365 = qJD(2) * pkin(2);
t364 = qJ(4) * t284;
t272 = t287 * t343;
t342 = qJD(2) * t283;
t324 = t291 * t342;
t315 = qJD(1) * t324;
t338 = qJD(3) * t290;
t209 = qJD(3) * t272 + t262 * t338 + t287 * t315;
t363 = t209 * t282;
t362 = t209 * t284;
t341 = qJD(2) * t287;
t257 = qJD(3) * t282 + t284 * t341;
t289 = cos(qJ(5));
t360 = t257 * t289;
t267 = -pkin(3) * t290 - qJ(4) * t287 - pkin(2);
t359 = t267 * t284;
t357 = t283 * t288;
t356 = t283 * t291;
t355 = t284 * t289;
t286 = sin(qJ(5));
t353 = t286 * t287;
t292 = qJD(3) ^ 2;
t352 = t287 * t292;
t351 = t289 * t290;
t350 = t290 * t292;
t207 = t290 * t315 + (qJD(4) + t367) * qJD(3);
t215 = (t244 + t328) * qJD(2);
t190 = t207 * t284 + t282 * t215;
t329 = pkin(7) * t282 + pkin(4);
t339 = qJD(3) * t287;
t349 = -t329 * t339 - t369;
t236 = t262 * t290 + t272;
t226 = qJD(3) * qJ(4) + t236;
t237 = qJD(2) * t267 - t327;
t197 = t226 * t284 + t282 * t237;
t330 = pkin(7) * t339;
t318 = t282 * t330;
t348 = t318 + t369;
t317 = t284 * t330;
t347 = -t317 + t368;
t259 = t307 * qJD(2);
t206 = t282 * t259 + t284 * t367;
t331 = qJD(2) * qJD(3);
t321 = t290 * t331;
t314 = t284 * t321;
t322 = t287 * t331;
t346 = t286 * t322 + t289 * t314;
t239 = t282 * t267 + t290 * t366;
t340 = qJD(2) * t290;
t337 = qJD(5) * t286;
t336 = qJD(5) * t289;
t196 = -t226 * t282 + t237 * t284;
t191 = pkin(4) * t340 - t196;
t335 = t191 * qJD(5);
t278 = t284 * qJD(3);
t334 = t289 * qJD(2);
t333 = t289 * qJD(3);
t255 = t282 * t341 - t278;
t332 = -qJD(5) - t255;
t325 = t288 * t342;
t323 = t286 * t340;
t320 = -qJD(3) * pkin(3) + qJD(4);
t188 = pkin(8) * t322 + t190;
t312 = pkin(4) * t282 - pkin(8) * t284;
t298 = t312 * t340;
t195 = qJD(3) * t298 + t209;
t319 = -t188 * t286 + t195 * t289;
t316 = t287 * t327;
t313 = t282 * MDP(20) * t340;
t311 = qJD(5) * t364 + t236 + t298;
t263 = -t327 - t365;
t310 = -t263 - t327;
t301 = t284 * t351 + t353;
t243 = t301 * qJD(2);
t309 = -t282 * t337 - t243;
t302 = pkin(7) + t312;
t241 = t302 * t287;
t308 = -qJD(5) * t241 - (pkin(8) - t366) * t339 - t368;
t306 = t188 * t289 + t286 * t195;
t192 = -pkin(8) * t340 + t197;
t223 = -t367 + t320;
t194 = pkin(4) * t255 - pkin(8) * t257 + t223;
t186 = t192 * t289 + t194 * t286;
t305 = t192 * t286 - t194 * t289;
t189 = -t207 * t282 + t215 * t284;
t247 = t285 * t287 + t290 * t357;
t220 = t247 * t284 - t282 * t356;
t246 = -t285 * t290 + t287 * t357;
t304 = t220 * t289 + t246 * t286;
t303 = -t220 * t286 + t246 * t289;
t205 = t259 * t284 - t282 * t367;
t249 = t284 * t353 + t351;
t229 = t257 * t286 + t290 * t334;
t264 = -pkin(4) * t284 - pkin(8) * t282 - pkin(3);
t299 = pkin(8) * t341 - qJD(4) * t284 - qJD(5) * t264 + t206;
t297 = qJD(5) * t323 - t286 * t314 + t289 * t322;
t228 = -pkin(8) * t290 + t239;
t296 = qJD(5) * t228 - t302 * t338 + t316;
t295 = qJD(3) * (-t310 - t365);
t294 = -qJ(4) * t339 + (-t223 + t320) * t290;
t293 = qJD(2) ^ 2;
t279 = t282 ^ 2;
t250 = -t286 * t290 + t287 * t355;
t242 = t284 * t323 - t287 * t334;
t238 = -pkin(7) * t358 + t359;
t231 = -t323 + t360;
t227 = t290 * t329 - t359;
t222 = qJD(3) * t247 + t287 * t324;
t221 = -qJD(3) * t246 + t290 * t324;
t219 = t247 * t282 + t283 * t354;
t211 = -t287 * t333 - t290 * t337 + (t286 * t338 + t287 * t336) * t284;
t210 = qJD(3) * t301 - qJD(5) * t249;
t204 = t221 * t284 + t282 * t325;
t203 = t221 * t282 - t284 * t325;
t201 = t257 * t336 - t297;
t200 = -qJD(5) * t229 + t346;
t198 = -pkin(4) * t341 - t205;
t187 = -pkin(4) * t322 - t189;
t184 = -qJD(5) * t186 + t319;
t183 = -qJD(5) * t305 + t306;
t1 = [t222 * t255 * MDP(12) + t222 * t257 * MDP(13) + (t203 * t257 - t204 * t255) * MDP(14) + (-t189 * t219 + t190 * t220 - t196 * t203 + t197 * t204 + t209 * t246 + t222 * t223) * MDP(15) + (-(-qJD(5) * t304 - t204 * t286 + t222 * t289) * t332 + t203 * t229 + t219 * t201) * MDP(21) + ((qJD(5) * t303 + t204 * t289 + t222 * t286) * t332 + t203 * t231 + t219 * t200) * MDP(22) + (-MDP(10) * t222 - MDP(11) * t221) * qJD(3) + (-MDP(4) * t291 + (-MDP(10) * t290 + MDP(11) * t287 - MDP(3)) * t288) * t293 * t283 + ((MDP(12) * t203 + MDP(13) * t204) * t290 + ((-MDP(10) * t356 - MDP(12) * t219 - MDP(13) * t220) * t287 + (-MDP(11) * t356 + (MDP(13) * t246 + MDP(14) * t219) * t284 + (t246 * MDP(12) - t220 * MDP(14) + MDP(21) * t303 - MDP(22) * t304) * t282) * t290) * qJD(3)) * qJD(2); 0.2e1 * t321 * t371 - 0.2e1 * t331 * t370 + MDP(7) * t350 - MDP(8) * t352 + (-pkin(7) * t350 + t287 * t295) * MDP(10) + (pkin(7) * t352 + t290 * t295) * MDP(11) + ((-t255 * t327 + t363 + (qJD(2) * t238 + t196) * qJD(3)) * t287 + (-t189 + (pkin(7) * t255 + t223 * t282) * qJD(3) + (t318 - t348) * qJD(2)) * t290) * MDP(12) + ((-t257 * t327 + t362 + (-qJD(2) * t239 - t197) * qJD(3)) * t287 + (t190 + (pkin(7) * t257 + t223 * t284) * qJD(3) + (t317 + t347) * qJD(2)) * t290) * MDP(13) + ((-t189 * t284 - t190 * t282) * t287 - t348 * t257 - t347 * t255 + (-t196 * t284 - t197 * t282 + (-t238 * t284 - t239 * t282) * qJD(2)) * t338) * MDP(14) + (-t223 * t316 + t189 * t238 + t190 * t239 + t347 * t197 + t348 * t196 + (t209 * t287 + t223 * t338) * pkin(7)) * MDP(15) + (t200 * t250 + t210 * t231) * MDP(16) + (-t200 * t249 - t201 * t250 - t210 * t229 - t211 * t231) * MDP(17) + (-t210 * t332 + (t200 * t287 + (qJD(2) * t250 + t231) * t338) * t282) * MDP(18) + (t211 * t332 + (-t201 * t287 + (-qJD(2) * t249 - t229) * t338) * t282) * MDP(19) + (t279 * t341 - t282 * t332) * MDP(20) * t338 + (t187 * t249 + t191 * t211 + t227 * t201 - (t286 * t308 - t289 * t296) * t332 + t349 * t229 + (t184 * t287 + ((-t228 * t286 + t241 * t289) * qJD(2) - t305) * t338) * t282) * MDP(21) + (t187 * t250 + t191 * t210 + t227 * t200 - (t286 * t296 + t289 * t308) * t332 + t349 * t231 + (-t183 * t287 + (-(t228 * t289 + t241 * t286) * qJD(2) - t186) * t338) * t282) * MDP(22); (qJD(3) * t236 - t263 * t341 - t209) * MDP(10) + t310 * t340 * MDP(11) + (-t362 - t236 * t255 + (-t196 * t287 + t205 * t290 + t282 * t294) * qJD(2)) * MDP(12) + (t363 - t236 * t257 + (t197 * t287 - t206 * t290 + t284 * t294) * qJD(2)) * MDP(13) + (t205 * t257 + t206 * t255 + (-qJD(4) * t255 + t196 * t340 + t190) * t284 + (qJD(4) * t257 + t197 * t340 - t189) * t282) * MDP(14) + (-pkin(3) * t209 - t196 * t205 - t197 * t206 - t223 * t236 + (-t196 * t282 + t197 * t284) * qJD(4) + (-t189 * t282 + t190 * t284) * qJ(4)) * MDP(15) + (t200 * t282 * t289 + t231 * t309) * MDP(16) + (t229 * t243 + t231 * t242 + (-t200 * t286 - t201 * t289 + (t229 * t286 - t231 * t289) * qJD(5)) * t282) * MDP(17) + (-t200 * t284 - t309 * t332 + (-t231 * t282 + t279 * t333) * t340) * MDP(18) + (t201 * t284 - (-t282 * t336 + t242) * t332 + (-qJD(3) * t279 * t286 + t229 * t282) * t340) * MDP(19) + (t332 - t278) * t313 + (-t184 * t284 - t191 * t242 - t198 * t229 - (t286 * t299 - t289 * t311) * t332 + (t289 * t335 + qJ(4) * t201 + qJD(4) * t229 + t187 * t286 + ((t264 * t289 - t286 * t364) * qJD(3) + t305) * t340) * t282) * MDP(21) + (t183 * t284 - t191 * t243 - t198 * t231 - (t286 * t311 + t289 * t299) * t332 + (-t286 * t335 + qJ(4) * t200 + qJD(4) * t231 + t187 * t289 + (-(qJ(4) * t355 + t264 * t286) * qJD(3) + t186) * t340) * t282) * MDP(22) + (-t290 * t371 + t370) * t293; -t255 ^ 2 * MDP(14) + (t197 * t255 + t209) * MDP(15) + (-MDP(14) * t257 + MDP(15) * t196 - MDP(21) * t229 - MDP(22) * t231) * t257 + (-t257 * MDP(12) + t255 * MDP(13) + (t284 * MDP(13) + (MDP(21) * t289 - MDP(22) * t286 + MDP(12)) * t282) * qJD(3)) * t340 - (MDP(21) * t286 + MDP(22) * t289) * t332 ^ 2; t231 * t229 * MDP(16) + (-t229 ^ 2 + t231 ^ 2) * MDP(17) + (-t229 * t332 + t346) * MDP(18) + (-t231 * t332 + t297) * MDP(19) + qJD(3) * t313 + (-t186 * t332 - t191 * t231 + t319) * MDP(21) + (t191 * t229 + t305 * t332 - t306) * MDP(22) + (-MDP(18) * t229 - MDP(19) * t360 - MDP(21) * t186 + MDP(22) * t305) * qJD(5);];
tauc = t1;
