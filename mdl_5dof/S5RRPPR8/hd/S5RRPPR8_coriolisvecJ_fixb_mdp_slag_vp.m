% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:19
% EndTime: 2019-12-31 19:39:24
% DurationCPUTime: 2.40s
% Computational Cost: add. (1177->265), mult. (2906->368), div. (0->0), fcn. (1870->6), ass. (0->127)
t302 = cos(pkin(8));
t305 = sin(qJ(2));
t329 = qJD(1) * qJD(2);
t326 = t305 * t329;
t278 = t302 * t326;
t301 = sin(pkin(8));
t307 = cos(qJ(2));
t325 = t307 * t329;
t243 = t301 * t325 - t278;
t261 = t301 * t305 + t302 * t307;
t254 = t261 * qJD(2);
t244 = qJD(1) * t254;
t249 = t261 * qJD(1);
t335 = qJD(1) * t307;
t336 = qJD(1) * t305;
t251 = -t301 * t335 + t302 * t336;
t304 = sin(qJ(5));
t306 = cos(qJ(5));
t330 = qJD(5) * t306;
t331 = qJD(5) * t304;
t182 = -t304 * t243 + t306 * t244 - t249 * t330 - t251 * t331;
t183 = t306 * t243 + t304 * t244 - t249 * t331 + t251 * t330;
t209 = -t249 * t304 + t251 * t306;
t296 = qJD(2) - qJD(5);
t344 = t209 * t296;
t352 = t306 * t249 + t251 * t304;
t345 = t352 * t296;
t362 = (t182 - t345) * MDP(21) - (t183 + t344) * MDP(22) + (t209 ^ 2 - t352 ^ 2) * MDP(20) + t352 * MDP(19) * t209;
t334 = qJD(2) * t305;
t347 = pkin(6) - qJ(4);
t246 = -qJD(4) * t307 - t347 * t334;
t297 = qJD(2) * qJD(3);
t222 = t246 * qJD(1) + t297;
t282 = pkin(6) * t325;
t332 = qJD(4) * t305;
t333 = qJD(2) * t307;
t231 = t282 + (-qJ(4) * t333 - t332) * qJD(1);
t322 = t222 * t301 - t302 * t231;
t186 = -pkin(7) * t244 - t322;
t193 = t302 * t222 + t301 * t231;
t187 = -pkin(7) * t243 + t193;
t303 = qJD(1) * pkin(1);
t259 = -pkin(2) * t335 - qJ(3) * t336 - t303;
t226 = pkin(3) * t335 + qJD(4) - t259;
t204 = pkin(4) * t249 + t226;
t356 = -t306 * t186 + t304 * t187 + t204 * t209;
t355 = t304 * t186 + t306 * t187 - t204 * t352;
t354 = MDP(14) * pkin(6);
t299 = t305 ^ 2;
t353 = (-t307 ^ 2 + t299) * MDP(5);
t351 = qJD(5) + t296;
t308 = -pkin(2) - pkin(3);
t350 = pkin(1) * MDP(9);
t349 = pkin(7) * t249;
t348 = pkin(7) * t251;
t346 = qJD(2) * pkin(2);
t310 = qJD(1) ^ 2;
t341 = t307 * t310;
t277 = t347 * t307;
t247 = qJD(2) * t277 - t332;
t200 = t302 * t246 + t301 * t247;
t288 = pkin(6) * t336;
t267 = qJ(4) * t336 - t288;
t327 = t308 * qJD(2);
t239 = qJD(3) + t327 - t267;
t289 = pkin(6) * t335;
t269 = -qJ(4) * t335 + t289;
t298 = qJD(2) * qJ(3);
t258 = t269 + t298;
t198 = t301 * t239 + t302 * t258;
t214 = t302 * t267 + t301 * t269;
t276 = t347 * t305;
t219 = t301 * t276 + t302 * t277;
t292 = t305 * qJD(3);
t340 = qJ(3) * t325 + qJD(1) * t292;
t339 = qJ(3) * t333 + t292;
t274 = -t307 * pkin(2) - t305 * qJ(3) - pkin(1);
t328 = t305 * t341;
t324 = MDP(12) + t354;
t323 = qJD(3) - t346;
t270 = -qJ(3) * t301 + t302 * t308;
t197 = t302 * t239 - t258 * t301;
t199 = -t246 * t301 + t302 * t247;
t213 = -t267 * t301 + t302 * t269;
t218 = t302 * t276 - t277 * t301;
t260 = t307 * pkin(3) - t274;
t232 = pkin(2) * t326 - t340;
t248 = pkin(2) * t334 - t339;
t321 = -qJD(1) * t248 - t232;
t320 = qJD(3) * t301 + t213;
t319 = qJD(3) * t302 - t214;
t318 = t305 * t327;
t188 = -qJD(2) * pkin(4) + t197 - t348;
t189 = t198 - t349;
t315 = t306 * t188 - t304 * t189;
t314 = -t304 * t188 - t306 * t189;
t313 = t197 * t301 - t198 * t302;
t262 = -t301 * t307 + t302 * t305;
t211 = t261 * t306 + t262 * t304;
t212 = -t261 * t304 + t262 * t306;
t286 = qJ(3) * t335;
t245 = t308 * t336 + t286;
t312 = t296 * (-t301 * t306 - t302 * t304);
t311 = t296 * (-t301 * t304 + t302 * t306);
t225 = t318 + t339;
t217 = qJD(1) * t318 + t340;
t309 = qJD(2) ^ 2;
t275 = t289 + t298;
t273 = t288 + t323;
t272 = -pkin(6) * t326 + t297;
t271 = qJ(3) * t302 + t301 * t308;
t268 = pkin(2) * t336 - t286;
t266 = -pkin(4) + t270;
t253 = t301 * t333 - t302 * t334;
t215 = pkin(4) * t261 + t260;
t210 = -pkin(4) * t251 + t245;
t203 = pkin(4) * t253 + t225;
t202 = -pkin(7) * t261 + t219;
t201 = -pkin(7) * t262 + t218;
t196 = pkin(4) * t243 + t217;
t195 = t214 + t348;
t194 = t213 - t349;
t191 = -pkin(7) * t253 + t200;
t190 = -pkin(7) * t254 + t199;
t185 = t212 * qJD(5) + t306 * t253 + t254 * t304;
t184 = -t211 * qJD(5) - t253 * t304 + t254 * t306;
t1 = [(t232 * t274 + t248 * t259) * MDP(14) + (t217 * t261 + t225 * t249 + t226 * t253 + t243 * t260) * MDP(15) + (t217 * t262 + t225 * t251 + t226 * t254 + t244 * t260) * MDP(16) + (-t193 * t261 - t197 * t254 - t198 * t253 - t199 * t251 - t200 * t249 - t218 * t244 - t219 * t243 + t262 * t322) * MDP(17) + (t193 * t219 + t197 * t199 + t198 * t200 + t217 * t260 - t218 * t322 + t225 * t226) * MDP(18) + (t182 * t212 + t184 * t209) * MDP(19) + (-t182 * t211 - t183 * t212 - t184 * t352 - t185 * t209) * MDP(20) + (t215 * t183 + t204 * t185 + t196 * t211 + t203 * t352) * MDP(24) + (t215 * t182 + t204 * t184 + t196 * t212 + t203 * t209) * MDP(25) + (-t184 * MDP(21) + t185 * MDP(22) + (-t190 * t306 + t191 * t304) * MDP(24) + (t190 * t304 + t191 * t306) * MDP(25) + ((t201 * t304 + t202 * t306) * MDP(24) + (t201 * t306 - t202 * t304) * MDP(25)) * qJD(5)) * t296 + (t321 * MDP(13) + (-MDP(7) + (MDP(10) - MDP(13)) * pkin(6)) * t309) * t305 + (t309 * MDP(6) + t321 * MDP(11) + t272 * MDP(12) + (t272 * MDP(14) + (-MDP(11) - MDP(9)) * t309) * pkin(6)) * t307 + (-t199 * MDP(15) + t200 * MDP(16) - 0.2e1 * qJD(1) * t353 + (-0.2e1 * MDP(10) * t303 + (-qJD(1) * t274 - t259) * MDP(13) + t324 * t273) * t307 + (t259 * MDP(11) - t324 * t275 + (t274 * MDP(11) - 0.2e1 * t350 + (t324 * pkin(6) + (2 * MDP(4))) * t307) * qJD(1)) * t305) * qJD(2); -MDP(4) * t328 + pkin(1) * MDP(10) * t341 + 0.2e1 * t297 * MDP(13) + (qJ(3) * t272 + qJD(3) * t275 - t259 * t268) * MDP(14) + (t320 * qJD(2) + t226 * t251 - t245 * t249 + t322) * MDP(15) + (t319 * qJD(2) - t226 * t249 - t245 * t251 + t193) * MDP(16) + (-t243 * t271 - t244 * t270 + (-t198 + t320) * t251 + (t197 - t319) * t249) * MDP(17) + (-t313 * qJD(3) + t193 * t271 - t197 * t213 - t198 * t214 - t226 * t245 - t270 * t322) * MDP(18) + (-t210 * t352 + (t194 * t306 - t195 * t304) * t296 - qJD(3) * t312 + (-(-t266 * t304 - t271 * t306) * t296 - t314) * qJD(5) + t356) * MDP(24) + (-t210 * t209 - (t194 * t304 + t195 * t306) * t296 + qJD(3) * t311 + ((t266 * t306 - t271 * t304) * t296 + t315) * qJD(5) + t355) * MDP(25) + (t305 * t350 + t353) * t310 + ((-t259 * t305 + t268 * t307) * MDP(11) + ((t275 - t298) * t305 + (-t273 + t323) * t307) * MDP(12) + (t259 * t307 + t268 * t305) * MDP(13) + (t275 * t305 + (-t273 - t346) * t307) * t354) * qJD(1) - t362; -MDP(11) * t328 + (-t299 * t310 - t309) * MDP(13) + (-qJD(2) * t275 + t259 * t336 + t282) * MDP(14) + (-t249 * t336 - t301 * t309) * MDP(15) + (-t251 * t336 - t302 * t309) * MDP(16) + (-t243 * t301 - t244 * t302 + (t249 * t302 - t251 * t301) * qJD(2)) * MDP(17) + (t313 * qJD(2) + t193 * t301 - t226 * t336 - t302 * t322) * MDP(18) + (t296 * t312 - t336 * t352) * MDP(24) + (-t209 * t336 - t296 * t311) * MDP(25); -t278 * MDP(15) + (-t249 ^ 2 - t251 ^ 2) * MDP(17) + (t197 * t251 + t198 * t249 + t340) * MDP(18) + (t183 - t344) * MDP(24) + (t182 + t345) * MDP(25) + (-t251 * MDP(15) + t249 * MDP(16) + ((MDP(15) * t301 + MDP(16) * t302) * t307 + (t301 * MDP(16) + t308 * MDP(18)) * t305) * qJD(1)) * qJD(2); (t351 * t314 - t356) * MDP(24) + (-t351 * t315 - t355) * MDP(25) + t362;];
tauc = t1;
