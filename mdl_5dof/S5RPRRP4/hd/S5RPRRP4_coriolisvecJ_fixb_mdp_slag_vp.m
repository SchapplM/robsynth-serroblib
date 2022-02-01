% Calculate Coriolis joint torque vector for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP4_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRRP4_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:32:55
% EndTime: 2022-01-23 09:32:59
% DurationCPUTime: 2.06s
% Computational Cost: add. (2080->248), mult. (5493->334), div. (0->0), fcn. (3782->6), ass. (0->137)
t323 = qJD(3) + qJD(4);
t277 = sin(pkin(8));
t280 = sin(qJ(3));
t282 = cos(qJ(3));
t370 = (t282 * MDP(10) + t280 * MDP(9)) * t277;
t273 = t277 ^ 2;
t362 = 0.2e1 * t273;
t281 = cos(qJ(4));
t279 = sin(qJ(4));
t347 = t279 * t282;
t252 = t280 * t281 + t347;
t242 = t252 * t277;
t368 = (t280 ^ 2 - t282 ^ 2) * MDP(8);
t335 = qJD(1) * t277;
t319 = t280 * t335;
t302 = t279 * t319;
t318 = t282 * t335;
t303 = t281 * t318;
t238 = -t302 + t303;
t232 = t238 * qJ(5);
t278 = cos(pkin(8));
t255 = -pkin(2) * t278 - pkin(6) * t277 - pkin(1);
t245 = t255 * qJD(1) + qJD(2);
t241 = t282 * t245;
t349 = t278 * t280;
t360 = pkin(7) * t277;
t289 = -qJ(2) * t349 - t282 * t360;
t224 = t289 * qJD(1) + t241;
t334 = qJD(1) * t278;
t266 = -qJD(3) + t334;
t210 = -pkin(3) * t266 + t224;
t357 = qJ(2) * t282;
t322 = t278 * t357;
t225 = -pkin(7) * t319 + qJD(1) * t322 + t245 * t280;
t217 = t279 * t225;
t312 = t281 * t210 - t217;
t188 = -t232 + t312;
t274 = t278 ^ 2;
t367 = t362 + t274;
t366 = -MDP(12) * t280 - MDP(13) * t282;
t315 = -t255 + t360;
t364 = t315 * t280 - t322;
t363 = t238 ^ 2;
t261 = -qJD(4) + t266;
t361 = pkin(3) * t261;
t359 = MDP(6) * qJ(2);
t358 = qJ(2) * t280;
t215 = t323 * t242;
t208 = qJD(1) * t215;
t356 = qJ(5) * t208;
t290 = qJD(1) * t252;
t235 = t277 * t290;
t355 = qJ(5) * t235;
t354 = t235 * t261;
t353 = t238 * t261;
t246 = pkin(3) * t319 + qJ(2) * t335;
t352 = t246 * t238;
t283 = qJD(1) ^ 2;
t351 = t273 * t283;
t350 = t277 * t280;
t348 = t279 * t280;
t219 = t281 * t225;
t187 = -pkin(4) * t261 + t188;
t346 = t187 - t188;
t345 = t281 * t224 - t217;
t251 = t281 * t282 - t348;
t344 = (t323 - t334) * t251;
t230 = t323 * t252;
t343 = t278 * t290 - t230;
t325 = qJD(1) * qJD(2);
t316 = t278 * t325;
t330 = qJD(3) * t282;
t342 = t245 * t330 + t282 * t316;
t332 = qJD(2) * t282;
t341 = t255 * t330 + t278 * t332;
t340 = t323 * t302;
t264 = t277 * pkin(3) * t330;
t244 = qJD(1) * t264 + t277 * t325;
t248 = t277 * qJD(2) + t264;
t250 = pkin(3) * t350 + t277 * qJ(2);
t339 = t273 + t274;
t333 = qJD(2) * t280;
t331 = qJD(3) * t280;
t329 = qJD(4) * t279;
t328 = qJD(4) * t281;
t327 = qJD(3) + t266;
t314 = -pkin(4) * t235 - qJD(5);
t222 = -t314 + t246;
t326 = qJD(5) + t222;
t320 = qJ(2) * t331;
t317 = t278 * t333;
t288 = t289 * qJD(3);
t205 = qJD(1) * t288 + t342;
t206 = -t245 * t331 + (-t317 + (pkin(7) * t350 - t322) * qJD(3)) * qJD(1);
t313 = -t279 * t205 + t281 * t206;
t220 = t288 + t341;
t221 = t364 * qJD(3) - t317;
t311 = -t220 * t279 + t281 * t221;
t310 = -t224 * t279 - t219;
t309 = qJD(1) * t327;
t308 = t323 * t282;
t307 = pkin(3) * t318;
t306 = MDP(22) * t323;
t305 = t282 * t273 * t280 * MDP(7);
t304 = -t281 * t205 - t279 * t206 - t210 * t328 + t225 * t329;
t297 = t281 * t277 * t308;
t209 = qJD(1) * t297 - t340;
t199 = pkin(4) * t209 + t244;
t300 = -t210 * t279 - t219;
t226 = -t315 * t282 + (-pkin(3) - t358) * t278;
t299 = -t226 * t279 + t281 * t364;
t294 = qJ(5) * t209 + t304;
t293 = t246 * t235 + t304;
t233 = t235 ^ 2;
t292 = t238 * t235 * MDP(14) + (-t230 * t335 - t354) * MDP(16) + (-t323 * t303 + t340 - t353) * MDP(17) + (-t233 + t363) * MDP(15);
t291 = t281 * t220 + t279 * t221 + t226 * t328 + t329 * t364;
t287 = t326 * t235 + t294;
t286 = t300 * qJD(4) + t313;
t285 = t286 + t356;
t284 = (-t219 + (-t210 + t361) * t279) * qJD(4) + t313;
t270 = pkin(3) * t281 + pkin(4);
t249 = t328 * t361;
t243 = t251 * t277;
t228 = pkin(4) * t238 + t307;
t227 = pkin(4) * t242 + t250;
t216 = -t323 * t277 * t348 + t297;
t200 = pkin(4) * t216 + t248;
t195 = -qJ(5) * t242 - t299;
t194 = -pkin(4) * t278 - qJ(5) * t243 + t226 * t281 + t279 * t364;
t191 = -t232 + t345;
t190 = t310 + t355;
t189 = -t300 - t355;
t186 = qJ(5) * t215 + t299 * qJD(4) - qJD(5) * t243 + t311;
t185 = -qJ(5) * t216 - qJD(5) * t242 + t291;
t184 = -qJD(5) * t238 + t285;
t183 = -qJD(5) * t235 - t294;
t1 = [(t266 * t317 + t367 * qJD(1) * (qJ(2) * t330 + t333)) * MDP(12) + ((-t278 * t320 + t341) * t266 + t342 * t278 + (-t367 * t320 + t332 * t362) * qJD(1)) * MDP(13) + (-t208 * t243 - t215 * t238) * MDP(14) + (t208 * t242 - t209 * t243 + t215 * t235 - t216 * t238) * MDP(15) + (t208 * t278 + t215 * t261) * MDP(16) + (t209 * t278 + t216 * t261) * MDP(17) + (-t311 * t261 - t313 * t278 + t248 * t235 + t250 * t209 + t244 * t242 + t246 * t216 + (-t299 * t261 - t300 * t278) * qJD(4)) * MDP(19) + (-t250 * t208 - t246 * t215 + t248 * t238 + t244 * t243 + t291 * t261 - t304 * t278) * MDP(20) + (-t184 * t278 - t186 * t261 + t199 * t242 + t200 * t235 + t209 * t227 + t216 * t222) * MDP(21) + (t183 * t278 + t185 * t261 + t199 * t243 + t200 * t238 - t208 * t227 - t215 * t222) * MDP(22) + (-t183 * t242 - t184 * t243 - t185 * t235 - t186 * t238 + t187 * t215 - t189 * t216 + t194 * t208 - t195 * t209) * MDP(23) + (t183 * t195 + t184 * t194 + t185 * t189 + t186 * t187 + t199 * t227 + t200 * t222) * MDP(24) + 0.2e1 * (MDP(5) + t359) * t339 * t325 + ((-(-t255 * t280 - t322) * t266 + t245 * t349) * MDP(12) + (t362 * t368 - 0.2e1 * t305) * qJD(1) + (t266 + t334) * t370) * qJD(3); (t208 * t251 - t209 * t252 - t344 * t235 - t343 * t238) * MDP(23) + (t183 * t252 + t184 * t251 + t343 * t187 + t344 * t189 - t222 * t335) * MDP(24) + (MDP(19) + MDP(21)) * (-t235 * t335 - t343 * t261) + (MDP(20) + MDP(22)) * (-t238 * t335 + t344 * t261) + (-t274 * MDP(5) + (-MDP(5) + t366) * t273 - t339 * t359) * t283 + t366 * t266 ^ 2; t283 * t305 - t351 * t368 + ((-t327 * t245 - t316) * t280 + (-t278 * t309 - t351) * t357) * MDP(12) + (-t241 * t266 + (t327 * t334 + t351) * t358 - t342) * MDP(13) + (-t235 * t307 + t310 * t261 + t284 - t352) * MDP(19) + (-t238 * t307 - t345 * t261 + t249 + t293) * MDP(20) + (t190 * t261 - t228 * t235 - t326 * t238 + t284 + t356) * MDP(21) + (-t191 * t261 - t228 * t238 + t249 + t287) * MDP(22) + (t208 * t270 + (t189 + t190) * t238 + (-t187 + t191) * t235 + (-t209 * t279 + (-t235 * t281 + t238 * t279) * qJD(4)) * pkin(3)) * MDP(23) + (t184 * t270 - t187 * t190 - t189 * t191 - t222 * t228 + (t183 * t279 + (-t187 * t279 + t189 * t281) * qJD(4)) * pkin(3)) * MDP(24) + t292 - t309 * t370; (t300 * t261 + t286 - t352) * MDP(19) + (-t312 * t261 + t293) * MDP(20) + (-t189 * t261 + (-t222 + t314) * t238 + t285) * MDP(21) + (-pkin(4) * t363 - t188 * t261 + t287) * MDP(22) + (pkin(4) * t208 - t346 * t235) * MDP(23) + (t346 * t189 + (-t222 * t238 + t184) * pkin(4)) * MDP(24) + t292; (-t340 - t353) * MDP(21) + MDP(22) * t354 + (-t233 - t363) * MDP(23) + (t187 * t238 + t189 * t235 + t199) * MDP(24) + (-t306 * t347 + (MDP(21) * t308 - t280 * t306) * t281) * t335;];
tauc = t1;
