% Calculate Coriolis joint torque vector for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 20:52
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 20:52:25
% EndTime: 2021-01-15 20:52:31
% DurationCPUTime: 1.78s
% Computational Cost: add. (1357->270), mult. (3162->339), div. (0->0), fcn. (1861->4), ass. (0->127)
t282 = sin(qJ(4));
t283 = sin(qJ(2));
t284 = cos(qJ(4));
t285 = cos(qJ(2));
t293 = t283 * t282 + t285 * t284;
t332 = t283 * t284;
t347 = (MDP(17) * t293 + MDP(18) * t332) * qJD(4);
t323 = qJD(1) * t283;
t268 = pkin(6) * t323;
t242 = pkin(7) * t323 - t268;
t286 = -pkin(2) - pkin(3);
t313 = t286 * qJD(2);
t213 = qJD(3) + t313 - t242;
t322 = qJD(1) * t285;
t269 = pkin(6) * t322;
t244 = -pkin(7) * t322 + t269;
t278 = qJD(2) * qJ(3);
t227 = t244 + t278;
t295 = -t213 * t282 - t227 * t284;
t223 = t293 * qJD(1);
t336 = qJ(5) * t223;
t184 = -t295 - t336;
t276 = qJD(2) - qJD(4);
t346 = t184 * t276;
t334 = t223 * t276;
t225 = -t282 * t322 + t284 * t323;
t345 = t225 * t276;
t279 = t283 ^ 2;
t344 = (-t285 ^ 2 + t279) * MDP(5);
t304 = t284 * t213 - t227 * t282;
t335 = qJ(5) * t225;
t183 = t304 - t335;
t182 = -pkin(4) * t276 + t183;
t343 = t182 - t183;
t315 = qJD(1) * qJD(2);
t310 = t283 * t315;
t319 = qJD(4) * t282;
t311 = t285 * t319;
t327 = qJD(1) * t311 + t284 * t310;
t342 = t327 * MDP(18);
t340 = pkin(6) - pkin(7);
t318 = qJD(4) * t284;
t320 = qJD(2) * t285;
t290 = t282 * t320 + t283 * t318;
t195 = qJD(1) * t290 - t327;
t339 = pkin(4) * t195;
t338 = pkin(6) * MDP(14);
t281 = qJD(1) * pkin(1);
t337 = qJD(2) * pkin(2);
t333 = t282 * t285;
t303 = -t242 * t282 + t284 * t244;
t186 = t303 - t336;
t247 = qJ(3) * t284 + t282 * t286;
t215 = -qJD(3) * t282 - qJD(4) * t247;
t331 = -t186 + t215;
t329 = t284 * t242 + t282 * t244;
t187 = t329 + t335;
t296 = -qJ(3) * t282 + t284 * t286;
t214 = qJD(3) * t284 + qJD(4) * t296;
t330 = -t187 + t214;
t309 = t285 * t315;
t328 = t282 * t310 + t284 * t309;
t272 = t283 * qJD(3);
t326 = qJ(3) * t309 + qJD(1) * t272;
t325 = qJ(3) * t320 + t272;
t321 = qJD(2) * t283;
t228 = -pkin(2) * t322 - qJ(3) * t323 - t281;
t210 = pkin(3) * t322 - t228;
t307 = -pkin(4) * t223 - qJD(5);
t190 = t210 - t307;
t317 = t190 * MDP(25);
t316 = t228 * MDP(11);
t249 = -t285 * pkin(2) - t283 * qJ(3) - pkin(1);
t314 = t286 * t283;
t252 = t340 * t285;
t308 = qJD(2) * t282 * MDP(18);
t306 = MDP(12) + t338;
t305 = qJD(3) - t337;
t243 = t340 * t321;
t245 = qJD(2) * t252;
t302 = t243 * t282 + t284 * t245;
t236 = t285 * pkin(3) - t249;
t209 = pkin(2) * t310 - t326;
t221 = pkin(2) * t321 - t325;
t301 = -qJD(1) * t221 - t209;
t300 = t276 ^ 2;
t277 = qJD(2) * qJD(3);
t216 = -qJD(1) * t243 + t277;
t261 = pkin(6) * t309;
t235 = -pkin(7) * t309 + t261;
t299 = -t213 * t318 - t284 * t216 + t227 * t319 - t282 * t235;
t298 = -t213 * t319 - t282 * t216 - t227 * t318 + t284 * t235;
t297 = t283 * t313;
t251 = t340 * t283;
t294 = -t251 * t282 - t252 * t284;
t264 = qJ(3) * t322;
t219 = qJD(1) * t314 + t264;
t177 = -t195 * qJ(5) - t223 * qJD(5) - t299;
t194 = qJD(4) * t223 - t328;
t292 = qJ(5) * t194 + t298;
t205 = t297 + t325;
t291 = -t284 * t243 + t282 * t245 + t251 * t318 - t252 * t319;
t200 = qJD(1) * t297 + t326;
t289 = t223 * MDP(15) - t276 * MDP(18) - t210 * MDP(20);
t288 = qJD(1) ^ 2;
t287 = qJD(2) ^ 2;
t250 = t269 + t278;
t248 = t268 + t305;
t246 = -pkin(6) * t310 + t277;
t241 = -pkin(4) + t296;
t240 = pkin(2) * t323 - t264;
t239 = t332 - t333;
t222 = t223 ^ 2;
t204 = t215 * t276;
t203 = t214 * t276;
t199 = pkin(4) * t293 + t236;
t198 = t276 * t293;
t197 = -t284 * t321 + t290 - t311;
t196 = -pkin(4) * t225 + t219;
t189 = -qJ(5) * t293 - t294;
t188 = -qJ(5) * t239 + t251 * t284 - t252 * t282;
t185 = pkin(4) * t197 + t205;
t181 = t200 + t339;
t180 = -qJ(5) * t198 + qJD(4) * t294 - qJD(5) * t239 + t302;
t179 = -qJ(5) * t197 - qJD(5) * t293 + t291;
t178 = -qJD(5) * t225 + t292;
t1 = [(t209 * t249 + t221 * t228) * MDP(14) + (-t194 * t239 + t198 * t225) * MDP(15) + (t194 * t293 - t195 * t239 - t197 * t225 - t198 * t223) * MDP(16) + (t236 * t195 + t210 * t197 + t200 * t293 + t205 * t223) * MDP(20) + (-t236 * t194 + t210 * t198 + t200 * t239 + t205 * t225) * MDP(21) + (t181 * t293 + t185 * t223 + t190 * t197 + t195 * t199) * MDP(22) + (t181 * t239 + t185 * t225 + t190 * t198 - t194 * t199) * MDP(23) + (-t177 * t293 - t178 * t239 - t179 * t223 - t180 * t225 - t182 * t198 - t184 * t197 + t188 * t194 - t189 * t195) * MDP(24) + (t177 * t189 + t178 * t188 + t179 * t184 + t180 * t182 + t181 * t199 + t185 * t190) * MDP(25) + (-t198 * MDP(17) + t197 * MDP(18) + (t251 * t319 + t252 * t318 - t302) * MDP(20) + t291 * MDP(21) - t180 * MDP(22) + t179 * MDP(23)) * t276 + (t301 * MDP(13) + (-MDP(7) + (MDP(10) - MDP(13)) * pkin(6)) * t287) * t283 + (t287 * MDP(6) + t301 * MDP(11) + t246 * MDP(12) + (t246 * MDP(14) + (-MDP(11) - MDP(9)) * t287) * pkin(6)) * t285 + (-0.2e1 * qJD(1) * t344 + (-0.2e1 * MDP(10) * t281 + (-qJD(1) * t249 - t228) * MDP(13) + t306 * t248) * t285 + (t316 - t306 * t250 + (t249 * MDP(11) - 0.2e1 * pkin(1) * MDP(9) + (pkin(6) * t306 + (2 * MDP(4))) * t285) * qJD(1)) * t283) * qJD(2); 0.2e1 * t277 * MDP(13) + (qJ(3) * t246 + qJD(3) * t250 - t228 * t240) * MDP(14) - t328 * MDP(17) - t342 + (t276 * t303 - t204 - t298) * MDP(20) + (-t276 * t329 + t203 - t299) * MDP(21) + (t186 * t276 - t204 - t292) * MDP(22) + (-t187 * t276 + t177 + t203) * MDP(23) + (t194 * t241 - t195 * t247) * MDP(24) + (t177 * t247 + t178 * t241 + t182 * t331 + t184 * t330 - t190 * t196) * MDP(25) + (t276 * MDP(17) - t219 * MDP(20) - t210 * MDP(21) - t196 * MDP(22) - t190 * MDP(23) + (t182 - t330) * MDP(24) + MDP(16) * t223) * t223 + (-t219 * MDP(21) + (qJD(5) + t190) * MDP(22) - t196 * MDP(23) + (-t184 - t331) * MDP(24) - MDP(16) * t225 - t289) * t225 + (-t283 * t285 * MDP(4) + t344 + (MDP(10) * t285 + MDP(9) * t283) * pkin(1)) * t288 + ((-t316 + (t250 - t278) * MDP(12) + t240 * MDP(13)) * t283 + (t240 * MDP(11) + (-t248 + t305) * MDP(12) + t228 * MDP(13) + t308) * t285 + t347 + (t250 * t283 + (-t248 - t337) * t285) * t338) * qJD(1); (-t279 * t288 - t287) * MDP(13) + (-qJD(2) * t250 + t261) * MDP(14) + (-t288 * t285 * MDP(11) + (MDP(14) * t228 - t317) * qJD(1)) * t283 + (MDP(21) + MDP(23)) * (-t225 * t323 - t284 * t300) + (MDP(20) + MDP(22)) * (-t223 * t323 - t282 * t300) + ((t194 + t334) * MDP(24) + (t178 - t346) * MDP(25)) * t284 + ((-t195 - t345) * MDP(24) + (t182 * t276 + t177) * MDP(25)) * t282; -t222 * MDP(16) + (t328 - t334) * MDP(17) + t342 + (t276 * t295 + t298) * MDP(20) + (t210 * t223 - t276 * t304 + t299) * MDP(21) + (t292 - t346) * MDP(22) + (-t183 * t276 + t190 * t223 - t177) * MDP(23) + (pkin(4) * t194 - t223 * t343) * MDP(24) + (pkin(4) * t178 + t184 * t343) * MDP(25) + (-t285 * t308 - t347) * qJD(1) + ((-t190 + t307) * MDP(22) - pkin(4) * t317 + (-MDP(23) * pkin(4) + MDP(16)) * t225 + t289) * t225; (-t327 - t345) * MDP(22) + (t328 + t334) * MDP(23) + (-t225 ^ 2 - t222) * MDP(24) + (t182 * t225 + t184 * t223 + t326 + t339) * MDP(25) + ((MDP(22) * t332 - MDP(23) * t293) * qJD(4) + (MDP(22) * t333 + MDP(25) * t314) * qJD(2)) * qJD(1);];
tauc = t1;
