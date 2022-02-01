% Calculate Coriolis joint torque vector for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRPRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:03:00
% EndTime: 2022-01-20 11:03:04
% DurationCPUTime: 1.58s
% Computational Cost: add. (1738->205), mult. (3017->275), div. (0->0), fcn. (2230->8), ass. (0->112)
t281 = cos(qJ(2));
t335 = pkin(1) * qJD(1);
t314 = t281 * t335;
t298 = qJD(3) - t314;
t273 = qJD(1) + qJD(2);
t275 = cos(pkin(9));
t280 = cos(qJ(4));
t326 = t280 * t275;
t309 = t273 * t326;
t274 = sin(pkin(9));
t277 = sin(qJ(4));
t330 = t274 * t277;
t310 = t273 * t330;
t233 = -t309 + t310;
t279 = cos(qJ(5));
t225 = t279 * t233;
t257 = qJD(4) * t309;
t229 = -qJD(4) * t310 + t257;
t254 = t274 * t280 + t275 * t277;
t246 = t254 * qJD(4);
t230 = t273 * t246;
t235 = t254 * t273;
t276 = sin(qJ(5));
t316 = qJD(5) * t276;
t184 = -qJD(5) * t225 + t279 * t229 - t276 * t230 - t235 * t316;
t204 = t235 * t276 + t225;
t293 = t233 * t276 - t279 * t235;
t284 = t293 * qJD(5) - t229 * t276 - t279 * t230;
t272 = qJD(4) + qJD(5);
t332 = t204 * t272;
t333 = t293 * t272;
t351 = (-t204 ^ 2 + t293 ^ 2) * MDP(18) - t204 * MDP(17) * t293 + (t184 + t332) * MDP(19) + (t284 - t333) * MDP(20);
t334 = pkin(1) * qJD(2);
t311 = qJD(1) * t334;
t252 = t273 * qJD(3) + t281 * t311;
t319 = t274 ^ 2 + t275 ^ 2;
t350 = t319 * t252;
t288 = t254 * t252;
t278 = sin(qJ(2));
t312 = t278 * t335;
t256 = qJ(3) * t273 + t312;
t307 = pkin(7) * t273 + t256;
t222 = t307 * t274;
t223 = t307 * t275;
t294 = t222 * t277 - t223 * t280;
t181 = -pkin(8) * t229 + t294 * qJD(4) - t288;
t193 = -pkin(8) * t233 - t294;
t267 = -t275 * pkin(3) - pkin(2);
t232 = t267 * t273 + t298;
t208 = t233 * pkin(4) + t232;
t349 = t208 * t204 + t193 * t316 + (-t193 * t272 - t181) * t276;
t345 = -t280 * t222 - t223 * t277;
t253 = -t326 + t330;
t346 = t253 * t252;
t180 = -pkin(8) * t230 + t345 * qJD(4) - t346;
t348 = -t276 * t180 + t279 * t181 + t208 * t293;
t259 = (-pkin(7) - qJ(3)) * t274;
t268 = t275 * pkin(7);
t260 = qJ(3) * t275 + t268;
t291 = -t259 * t277 - t260 * t280;
t344 = t291 * qJD(4) - t298 * t254;
t343 = qJD(5) - t272;
t317 = qJD(4) * t280;
t342 = -t253 * t314 + (qJD(3) * t274 + qJD(4) * t260) * t277 - qJD(3) * t326 - t259 * t317;
t341 = t278 * t275 * MDP(7) + t281 * MDP(6);
t340 = pkin(1) * t281;
t339 = pkin(4) * t235;
t245 = t253 * qJD(4);
t338 = pkin(8) * t245;
t337 = pkin(8) * t254;
t336 = t246 * pkin(4);
t327 = t279 * t193;
t213 = t253 * t279 + t254 * t276;
t195 = -t213 * qJD(5) - t245 * t279 - t246 * t276;
t214 = -t253 * t276 + t254 * t279;
t264 = t278 * t311;
t215 = pkin(4) * t230 + t264;
t325 = t208 * t195 + t215 * t214;
t196 = t214 * qJD(5) - t245 * t276 + t279 * t246;
t324 = t208 * t196 + t215 * t213;
t323 = -t232 * t245 + t254 * t264;
t322 = t232 * t246 + t253 * t264;
t313 = t278 * t334;
t192 = -pkin(8) * t235 + t345;
t191 = qJD(4) * pkin(4) + t192;
t308 = -pkin(4) * t272 - t191;
t262 = t281 * t334 + qJD(3);
t306 = t262 * t319;
t305 = t319 * t281;
t301 = t319 * qJD(3);
t242 = t246 * pkin(8);
t297 = -qJD(5) * (t259 * t280 - t260 * t277 - t337) + t242 + t342;
t250 = t253 * pkin(8);
t296 = qJD(5) * (-t250 - t291) - t338 - t344;
t266 = pkin(1) * t278 + qJ(3);
t247 = (-pkin(7) - t266) * t274;
t248 = t266 * t275 + t268;
t292 = -t247 * t277 - t248 * t280;
t231 = t253 * pkin(4) + t267;
t289 = -t312 + t336;
t286 = (-t184 * t213 - t195 * t204 + t196 * t293 + t214 * t284) * MDP(18) + (t184 * t214 - t195 * t293) * MDP(17) + (-t229 * t253 - t230 * t254 + t233 * t245 - t235 * t246) * MDP(11) + (t229 * t254 - t235 * t245) * MDP(10) + (t195 * MDP(19) - t196 * MDP(20)) * t272 + (-t245 * MDP(12) - t246 * MDP(13)) * qJD(4);
t285 = t247 * t317 + t262 * t326 + (-qJD(4) * t248 - t262 * t274) * t277;
t283 = t292 * qJD(4) - t254 * t262;
t258 = t267 - t340;
t255 = -t273 * pkin(2) + t298;
t224 = t313 + t336;
t220 = t231 - t340;
t202 = -t250 - t292;
t201 = t247 * t280 - t248 * t277 - t337;
t187 = t283 + t338;
t186 = -t242 + t285;
t1 = [t286 + (t273 * t306 + t350) * MDP(8) + (-t224 * t293 + t220 * t184 - (t186 * t279 + t187 * t276 + (t201 * t279 - t202 * t276) * qJD(5)) * t272 + t325) * MDP(23) + (t224 * t204 - t220 * t284 + (-t186 * t276 + t187 * t279 + (-t201 * t276 - t202 * t279) * qJD(5)) * t272 + t324) * MDP(22) + (-t273 * t313 - t264) * MDP(5) + (t283 * qJD(4) + t258 * t230 + t233 * t313 + t322) * MDP(15) + (-t285 * qJD(4) + t258 * t229 + t235 * t313 + t323) * MDP(16) + (t256 * t306 + t266 * t350 + (t255 + (-pkin(2) - t340) * qJD(1)) * t313) * MDP(9) + t341 * (-qJD(1) - t273) * t334; ((-t305 * t335 + t301) * t273 + t350) * MDP(8) + (-t231 * t284 + (t297 * t276 - t296 * t279) * t272 + t289 * t204 + t324) * MDP(22) + t286 + (t231 * t184 + (t296 * t276 + t297 * t279) * t272 - t289 * t293 + t325) * MDP(23) + (qJD(4) * t342 + t267 * t229 - t235 * t312 + t323) * MDP(16) + (t273 * t312 - t264) * MDP(5) + (t344 * qJD(4) + t267 * t230 - t233 * t312 + t322) * MDP(15) + (t256 * t301 + qJ(3) * t350 + ((-pkin(2) * qJD(2) - t255) * t278 - t256 * t305) * t335) * MDP(9) + t341 * (-qJD(2) + t273) * t335; (-t319 * t273 * t256 + t264) * MDP(9) + t257 * MDP(16) + (-t284 - t333) * MDP(22) + (t184 - t332) * MDP(23) - t319 * MDP(8) * t273 ^ 2 + (0.2e1 * t235 * MDP(15) + (-t233 - t310) * MDP(16)) * qJD(4); t235 * t233 * MDP(10) + (-t233 ^ 2 + t235 ^ 2) * MDP(11) + (t257 + (t233 - t310) * qJD(4)) * MDP(12) + (-t232 * t235 - t288) * MDP(15) + (t232 * t233 + t346) * MDP(16) + (-t204 * t339 - (-t192 * t276 - t327) * t272 + (t308 * t276 - t327) * qJD(5) + t348) * MDP(22) + (t293 * t339 + (t308 * qJD(5) + t192 * t272 - t180) * t279 + t349) * MDP(23) + t351; (t343 * (-t191 * t276 - t327) + t348) * MDP(22) + ((-t191 * t343 - t180) * t279 + t349) * MDP(23) + t351;];
tauc = t1;
