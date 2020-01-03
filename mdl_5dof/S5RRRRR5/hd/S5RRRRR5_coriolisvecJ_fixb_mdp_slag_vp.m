% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:26
% EndTime: 2020-01-03 12:13:30
% DurationCPUTime: 1.44s
% Computational Cost: add. (1517->201), mult. (2445->280), div. (0->0), fcn. (1484->8), ass. (0->123)
t344 = pkin(8) + pkin(9);
t257 = qJD(1) + qJD(2);
t254 = qJD(3) + t257;
t260 = sin(qJ(5));
t261 = sin(qJ(4));
t264 = cos(qJ(5));
t265 = cos(qJ(4));
t227 = t260 * t265 + t264 * t261;
t256 = qJD(4) + qJD(5);
t342 = t256 * t227;
t187 = t254 * t342;
t263 = sin(qJ(2));
t330 = pkin(1) * qJD(2);
t301 = qJD(1) * t330;
t331 = pkin(1) * qJD(1);
t305 = t263 * t331;
t343 = qJD(3) * t305 + t263 * t301;
t338 = qJD(5) - t256;
t319 = t264 * t265;
t326 = t260 * t261;
t226 = -t319 + t326;
t275 = t226 * t256;
t262 = sin(qJ(3));
t267 = cos(qJ(2));
t266 = cos(qJ(3));
t322 = t263 * t266;
t280 = t262 * t267 + t322;
t312 = qJD(3) * t262;
t288 = pkin(2) * t312 - t280 * t331;
t278 = t288 * t254;
t302 = t267 * t331;
t233 = t257 * pkin(2) + t302;
t289 = t267 * t301;
t192 = t233 * t312 + t262 * t289 + t343 * t266;
t341 = -t278 - t192;
t340 = MDP(10) * t261;
t339 = (t261 ^ 2 - t265 ^ 2) * MDP(11);
t314 = t343 * t262;
t191 = (qJD(3) * t233 + t289) * t266 - t314;
t336 = t254 * pkin(3);
t335 = t265 * pkin(4);
t334 = t266 * pkin(2);
t250 = t267 * pkin(1) + pkin(2);
t219 = pkin(1) * t322 + t262 * t250 + pkin(8);
t333 = -pkin(9) - t219;
t248 = t262 * pkin(2) + pkin(8);
t332 = -pkin(9) - t248;
t311 = qJD(3) * t266;
t198 = t250 * t312 + (t280 * qJD(2) + t263 * t311) * pkin(1);
t329 = t198 * t254;
t213 = t262 * t233 + t266 * t305;
t328 = t213 * t254;
t268 = qJD(4) ^ 2;
t327 = t248 * t268;
t324 = t262 * t263;
t323 = t263 * MDP(5);
t295 = t344 * t254 + t213;
t194 = t295 * t265;
t321 = t264 * t194;
t310 = qJD(4) * t261;
t252 = pkin(4) * t310;
t185 = t254 * t252 + t192;
t245 = t262 * t305;
t212 = t266 * t233 - t245;
t251 = -pkin(3) - t335;
t196 = t251 * t254 - t212;
t318 = t185 * t227 - t196 * t275;
t317 = t185 * t226 + t196 * t342;
t205 = -t212 - t336;
t309 = qJD(4) * t265;
t316 = t192 * t261 + t205 * t309;
t315 = t252 + t288;
t308 = t268 * MDP(13);
t307 = -qJD(1) - t257;
t306 = t261 * t254 * pkin(4);
t304 = pkin(2) * t311;
t300 = t254 * t326;
t299 = t254 * t319;
t298 = qJD(4) * t344;
t297 = t254 * t309;
t193 = t295 * t261;
t190 = qJD(4) * pkin(4) - t193;
t296 = -pkin(4) * t256 - t190;
t294 = qJD(4) * t333;
t293 = qJD(4) * t332;
t218 = pkin(1) * t324 - t266 * t250 - pkin(3);
t287 = -t213 + t252;
t285 = pkin(8) * t268 - t328;
t284 = qJD(4) * (t212 - t336);
t283 = qJD(4) * t295;
t281 = t219 * t268 + t329;
t197 = t250 * t311 + (-t263 * t312 + (t266 * t267 - t324) * qJD(2)) * pkin(1);
t279 = qJD(4) * (t218 * t254 - t197);
t174 = t265 * t191 - t261 * t283;
t175 = -t261 * t191 - t265 * t283;
t216 = t227 * t254;
t277 = -t260 * t174 + t264 * t175 - t196 * t216;
t186 = qJD(5) * t299 - t256 * t300 + t264 * t297;
t214 = -t299 + t300;
t274 = t216 * t214 * MDP(17) + (t214 * t256 + t186) * MDP(19) + (t216 * t256 - t187) * MDP(20) + (-t214 ^ 2 + t216 ^ 2) * MDP(18);
t273 = (-t186 * t226 - t227 * t187 + t214 * t275 - t216 * t342) * MDP(18) + (t186 * t227 - t216 * t275) * MDP(17) - 0.2e1 * t254 * qJD(4) * t339 + 0.2e1 * t297 * t340 + t268 * t265 * MDP(12) + (-MDP(19) * t275 - MDP(20) * t342) * t256;
t221 = t266 * t302 - t245;
t271 = qJD(4) * ((-pkin(3) - t334) * t254 + t221 - t304);
t270 = t196 * t214 + (t338 * t194 - t175) * t260;
t269 = -t261 * t308 + t273;
t255 = t265 * pkin(9);
t244 = t265 * pkin(8) + t255;
t243 = t344 * t261;
t237 = t251 - t334;
t231 = t265 * t298;
t230 = t261 * t298;
t225 = t265 * t248 + t255;
t224 = t332 * t261;
t217 = t218 - t335;
t210 = -t261 * t304 + t265 * t293;
t209 = t261 * t293 + t265 * t304;
t208 = t265 * t219 + t255;
t207 = t333 * t261;
t201 = t205 * t310;
t195 = t252 + t198;
t184 = -t261 * t197 + t265 * t294;
t183 = t265 * t197 + t261 * t294;
t1 = [(-t192 - t329) * MDP(8) + (-t197 * t254 - t233 * t311 + t314) * MDP(9) + t201 * MDP(15) + t316 * MDP(16) + (t195 * t214 + t217 * t187 + (-t260 * t183 + t264 * t184 + (-t207 * t260 - t208 * t264) * qJD(5)) * t256 + t317) * MDP(22) + (t195 * t216 + t217 * t186 - (t264 * t183 + t260 * t184 + (t207 * t264 - t208 * t260) * qJD(5)) * t256 + t318) * MDP(23) + ((-t192 - t281) * MDP(15) + MDP(16) * t279) * t265 + (MDP(15) * t279 + t281 * MDP(16) - t308) * t261 + (t307 * t323 + (-qJD(1) * t266 * MDP(9) + t307 * MDP(6)) * t267) * t330 + t273; (t201 + t261 * t271 + (-t327 + t341) * t265) * MDP(15) + ((t327 + t278) * t261 + t265 * t271 + t316) * MDP(16) + (t221 * t254 + (-t289 + (-pkin(2) * t254 - t233) * qJD(3)) * t266 + t314) * MDP(9) + t341 * MDP(8) + (t237 * t187 + (-t260 * t209 + t264 * t210 + (-t224 * t260 - t225 * t264) * qJD(5)) * t256 + t221 * t342 + t315 * t214 + t317) * MDP(22) + (t237 * t186 - (t264 * t209 + t260 * t210 + (t224 * t264 - t225 * t260) * qJD(5)) * t256 - t221 * t275 + t315 * t216 + t318) * MDP(23) + t269 + (MDP(6) * t267 + t323) * (-qJD(2) + t257) * t331; (-t192 + t328) * MDP(8) + (t212 * t254 - t191) * MDP(9) + (t201 + t261 * t284 + (-t192 - t285) * t265) * MDP(15) + (t285 * t261 + t265 * t284 + t316) * MDP(16) + (t251 * t187 + (t260 * t230 - t264 * t231 + (t243 * t260 - t244 * t264) * qJD(5)) * t256 + t287 * t214 + t212 * t342 + t317) * MDP(22) + (t251 * t186 - (-t264 * t230 - t260 * t231 + (-t243 * t264 - t244 * t260) * qJD(5)) * t256 + t287 * t216 - t212 * t275 + t318) * MDP(23) + t269; (-t214 * t306 - (t260 * t193 - t321) * t256 + (t296 * t260 - t321) * qJD(5) + t277) * MDP(22) + (-t216 * t306 + (t296 * qJD(5) - t193 * t256 - t174) * t264 + t270) * MDP(23) + t274 + (t261 * MDP(15) + t265 * MDP(16)) * (-t205 * t254 - t191) + (-t265 * t340 + t339) * t254 ^ 2; (t277 + t338 * (-t260 * t190 - t321)) * MDP(22) + ((-t190 * t338 - t174) * t264 + t270) * MDP(23) + t274;];
tauc = t1;
