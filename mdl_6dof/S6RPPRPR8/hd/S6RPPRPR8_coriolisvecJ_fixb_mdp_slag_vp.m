% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPPRPR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPRPR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPPRPR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:56:14
% EndTime: 2019-03-09 01:56:19
% DurationCPUTime: 2.71s
% Computational Cost: add. (1762->292), mult. (3900->370), div. (0->0), fcn. (2701->6), ass. (0->125)
t262 = sin(pkin(9));
t331 = sin(qJ(4));
t294 = qJD(1) * t331;
t249 = t262 * t294;
t263 = cos(pkin(9));
t267 = cos(qJ(4));
t309 = qJD(1) * t267;
t296 = t263 * t309;
t234 = -t249 + t296;
t230 = qJD(6) + t234;
t240 = t267 * t262 + t263 * t331;
t310 = qJD(1) * t240;
t343 = qJD(4) * t310;
t311 = t262 ^ 2 + t263 ^ 2;
t342 = t311 * qJD(3);
t264 = -pkin(1) - qJ(3);
t337 = qJD(1) * t264;
t248 = qJD(2) + t337;
t291 = -pkin(7) * qJD(1) + t248;
t227 = t291 * t262;
t228 = t291 * t263;
t312 = -t331 * t227 + t267 * t228;
t341 = qJD(5) - t312;
t299 = MDP(16) - MDP(19);
t340 = MDP(6) * qJ(2) + t262 * MDP(7) + t263 * MDP(8) + MDP(5);
t339 = t230 ^ 2;
t239 = t262 * t331 - t267 * t263;
t206 = t239 * t343;
t293 = qJD(4) * t331;
t304 = qJD(4) * t267;
t236 = -t262 * t304 - t263 * t293;
t338 = t234 * t236 + t206;
t300 = pkin(5) * t234 + t341;
t203 = t267 * t227 + t331 * t228;
t198 = -qJD(4) * qJ(5) - t203;
t329 = pkin(5) * t310;
t185 = -t198 - t329;
t332 = pkin(4) + pkin(8);
t336 = t332 * t343 + (t185 - t203 + t329) * t230;
t328 = -pkin(7) + t264;
t241 = t328 * t262;
t242 = t328 * t263;
t196 = qJD(3) * t240 + t241 * t293 - t242 * t304;
t334 = t310 ^ 2;
t333 = t234 ^ 2;
t247 = qJD(4) * t249;
t295 = t263 * t304;
t226 = qJD(1) * t295 - t247;
t330 = pkin(4) * t226;
t257 = t262 * pkin(3);
t326 = qJ(5) * t310;
t253 = qJ(2) + t257;
t284 = qJ(5) * t239 + t253;
t192 = t240 * t332 + t284;
t325 = t192 * t343;
t265 = sin(qJ(6));
t303 = qJD(6) * t265;
t266 = cos(qJ(6));
t302 = qJD(6) * t266;
t313 = t265 * t226 + t302 * t310;
t193 = -qJD(4) * t303 + t313;
t324 = t193 * t266;
t261 = qJD(1) * qJ(2);
t256 = qJD(3) + t261;
t243 = qJD(1) * t257 + t256;
t281 = -qJ(5) * t234 + t243;
t199 = pkin(4) * t310 + t281;
t323 = t199 * t234;
t305 = qJD(4) * t265;
t209 = -t266 * t310 + t305;
t322 = t209 * t230;
t321 = t209 * t310;
t211 = qJD(4) * t266 + t265 * t310;
t320 = t211 * t230;
t319 = t211 * t310;
t318 = t310 * t234;
t316 = t240 * t265;
t217 = t266 * t343;
t308 = qJD(4) * t196;
t279 = t267 * t241 + t242 * t331;
t197 = -qJD(3) * t239 + qJD(4) * t279;
t307 = qJD(4) * t197;
t260 = qJD(1) * qJD(2);
t298 = -MDP(17) + MDP(20);
t292 = qJD(3) * t309;
t290 = qJ(5) * t343 + t260;
t207 = t331 * t241 - t242 * t267;
t289 = qJD(1) * t311;
t288 = t230 * t265;
t285 = qJD(3) * t294;
t287 = t227 * t293 - t228 * t304 + t262 * t292 + t263 * t285;
t180 = -qJD(4) * t332 + t300;
t186 = t310 * t332 + t281;
t173 = t180 * t266 - t186 * t265;
t174 = t180 * t265 + t186 * t266;
t278 = -t230 * t288 - t217;
t237 = -t262 * t293 + t295;
t277 = t237 * t265 + t240 * t302;
t276 = -qJD(5) * t234 + t290;
t275 = -qJ(5) * t236 + qJD(5) * t239 + qJD(2);
t181 = -qJD(4) * qJD(5) + t287;
t175 = -pkin(5) * t226 - t181;
t274 = t175 + (t230 * t332 + t326) * t230;
t200 = -pkin(5) * t239 + t207;
t273 = t175 * t240 + t185 * t237 + t200 * t343;
t272 = t265 * t343 - t339 * t266;
t184 = t227 * t304 + t228 * t293 - t262 * t285 + t263 * t292;
t195 = -qJD(4) * pkin(4) + t341;
t271 = -t181 * t240 + t184 * t239 - t195 * t236 - t198 * t237;
t270 = t203 * qJD(4) - t184;
t269 = qJD(1) ^ 2;
t218 = t266 * t226;
t205 = pkin(4) * t234 + t326;
t204 = pkin(4) * t240 + t284;
t201 = -t240 * pkin(5) + t279;
t194 = t211 * qJD(6) - t218;
t191 = pkin(4) * t237 + t275;
t187 = t276 + t330;
t183 = t236 * pkin(5) + t197;
t182 = -t237 * pkin(5) - t196;
t179 = t237 * t332 + t275;
t178 = t226 * t332 + t276;
t177 = -pkin(5) * t343 + t184;
t176 = t266 * t177;
t1 = [0.2e1 * qJD(3) * MDP(9) * t289 + ((t256 + t261) * qJD(2) + (-t248 - t337) * t342) * MDP(10) + t338 * MDP(11) + (t226 * t239 - t234 * t237 - t236 * t310 + t240 * t343) * MDP(12) + (0.2e1 * t310 * qJD(2) + t226 * t253 + t237 * t243 - t307) * MDP(16) + (t308 - t343 * t253 + t236 * t243 + (-qJD(1) * t239 + t234) * qJD(2)) * MDP(17) + (t196 * t310 + t197 * t234 - t207 * t343 - t226 * t279 - t271) * MDP(18) + (-t187 * t240 - t191 * t310 - t199 * t237 - t204 * t226 + t307) * MDP(19) + (t187 * t239 - t191 * t234 - t199 * t236 + t204 * t343 - t308) * MDP(20) + (-t181 * t279 + t184 * t207 + t187 * t204 + t191 * t199 + t195 * t197 + t196 * t198) * MDP(21) + (t193 * t316 + t211 * t277) * MDP(22) + ((-t209 * t265 + t211 * t266) * t237 + (t324 - t194 * t265 + (-t209 * t266 - t211 * t265) * qJD(6)) * t240) * MDP(23) + (-t193 * t239 + t211 * t236 + t230 * t277 - t316 * t343) * MDP(24) + (-t240 * t217 + t194 * t239 - t209 * t236 + (t237 * t266 - t240 * t303) * t230) * MDP(25) + (t230 * t236 + t206) * MDP(26) + (t173 * t236 - t176 * t239 + t182 * t209 + t201 * t194 + (t178 * t239 - t179 * t230 + t325) * t265 + (t183 * t230 - t273) * t266 + ((-t192 * t266 - t200 * t265) * t230 + t174 * t239 + t185 * t316) * qJD(6)) * MDP(27) + (-t174 * t236 + t182 * t211 + t201 * t193 + (-(qJD(6) * t200 + t179) * t230 + t325 + (qJD(6) * t180 + t178) * t239 + t185 * qJD(6) * t240) * t266 + (-(-qJD(6) * t192 + t183) * t230 + (-qJD(6) * t186 + t177) * t239 + t273) * t265) * MDP(28) + (t236 * MDP(13) - t237 * MDP(14)) * qJD(4) + 0.2e1 * t340 * t260; (-t226 * t240 - t237 * t310 - t338) * MDP(18) + t271 * MDP(21) + (t194 * t240 + t209 * t237 - t217 * t239) * MDP(27) + (t193 * t240 + t206 * t265 + t211 * t237) * MDP(28) + (t236 * t299 + t237 * t298) * qJD(4) + ((-t256 - t342) * MDP(10) - t199 * MDP(21) + t298 * t234 - t299 * t310) * qJD(1) - t340 * t269 + ((qJD(1) * t265 - t236 * t266 - t239 * t303) * MDP(27) + (qJD(1) * t266 + t236 * t265 - t239 * t302) * MDP(28)) * t230; -t311 * MDP(9) * t269 + (t248 * t289 + t260) * MDP(10) + (-t333 - t334) * MDP(18) + (t330 - t198 * t310 + (-qJD(5) - t195) * t234 + t290) * MDP(21) + (t272 + t321) * MDP(27) + (t339 * t265 + t217 + t319) * MDP(28) + t299 * (qJD(4) * (t234 + t296) - t247) + (-0.2e1 * MDP(17) + 0.2e1 * MDP(20)) * t343; MDP(11) * t318 + (t333 - t334) * MDP(12) + (t247 + (t234 - t296) * qJD(4)) * MDP(14) + (-t243 * t234 + t270) * MDP(16) + (qJD(4) * t312 + t243 * t310 + t287) * MDP(17) + (pkin(4) * t343 - qJ(5) * t226 + (-t198 - t203) * t234 + (t195 - t341) * t310) * MDP(18) + (t205 * t310 - t270 + t323) * MDP(19) + (-t199 * t310 + t205 * t234 + (0.2e1 * qJD(5) - t312) * qJD(4) - t287) * MDP(20) + (-pkin(4) * t184 - qJ(5) * t181 - t195 * t203 - t198 * t341 - t199 * t205) * MDP(21) + (-t211 * t288 + t324) * MDP(22) + ((-t194 - t320) * t266 + (-t193 + t322) * t265) * MDP(23) + (t278 + t319) * MDP(24) + (t272 - t321) * MDP(25) + t230 * t310 * MDP(26) + (qJ(5) * t194 + t173 * t310 + t300 * t209 + t274 * t265 + t336 * t266) * MDP(27) + (qJ(5) * t193 - t174 * t310 + t300 * t211 - t336 * t265 + t274 * t266) * MDP(28); -MDP(19) * t318 + (-qJD(4) ^ 2 - t333) * MDP(20) + (t198 * qJD(4) + t184 + t323) * MDP(21) + (-qJD(4) * t209 + t278) * MDP(27) + (-qJD(4) * t211 + t272) * MDP(28); t211 * t209 * MDP(22) + (-t209 ^ 2 + t211 ^ 2) * MDP(23) + (t313 + t322) * MDP(24) + (t218 + t320) * MDP(25) - t343 * MDP(26) + (t174 * t230 - t178 * t265 - t185 * t211 + t176) * MDP(27) + (t173 * t230 - t177 * t265 - t178 * t266 + t185 * t209) * MDP(28) + (-MDP(24) * t305 - MDP(25) * t211 - MDP(27) * t174 - MDP(28) * t173) * qJD(6);];
tauc  = t1;
