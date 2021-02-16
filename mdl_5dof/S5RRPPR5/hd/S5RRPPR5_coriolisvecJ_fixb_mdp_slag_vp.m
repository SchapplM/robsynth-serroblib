% Calculate Coriolis joint torque vector for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:36:36
% EndTime: 2021-01-15 19:36:43
% DurationCPUTime: 2.37s
% Computational Cost: add. (1315->246), mult. (3434->315), div. (0->0), fcn. (2341->6), ass. (0->117)
t290 = sin(pkin(8));
t291 = cos(pkin(8));
t293 = sin(qJ(2));
t295 = cos(qJ(2));
t267 = t290 * t295 + t291 * t293;
t256 = t267 * qJD(2);
t248 = qJD(1) * t256;
t328 = qJD(2) * t293;
t318 = qJD(1) * t328;
t276 = t290 * t318;
t337 = t291 * t295;
t319 = qJD(1) * t337;
t249 = qJD(2) * t319 - t276;
t292 = sin(qJ(5));
t294 = cos(qJ(5));
t329 = qJD(1) * t293;
t254 = t290 * t329 - t319;
t330 = qJD(1) * t267;
t363 = t294 * t254 - t330 * t292;
t186 = qJD(5) * t363 + t292 * t248 + t294 * t249;
t355 = t254 * t292 + t294 * t330;
t187 = qJD(5) * t355 - t248 * t294 + t249 * t292;
t285 = qJD(2) - qJD(5);
t342 = t363 * t285;
t343 = t355 * t285;
t369 = -(t187 + t343) * MDP(22) - t363 * MDP(19) * t355 + (t355 ^ 2 - t363 ^ 2) * MDP(20) + (t186 + t342) * MDP(21);
t366 = 0.2e1 * t330;
t320 = pkin(2) * t295 + pkin(1);
t362 = qJD(1) * t320;
t273 = qJD(3) - t362;
t364 = -qJ(4) * t330 + t273;
t346 = -qJ(3) - pkin(6);
t317 = qJD(2) * t346;
t251 = qJD(3) * t295 + t293 * t317;
t238 = t251 * qJD(1);
t301 = -qJD(3) * t293 + t295 * t317;
t239 = t301 * qJD(1);
t205 = t238 * t291 + t239 * t290;
t286 = qJD(2) * qJD(4);
t202 = t286 + t205;
t191 = pkin(7) * t248 + t202;
t204 = t238 * t290 - t239 * t291;
t192 = -pkin(7) * t249 + t204;
t351 = -pkin(3) - pkin(4);
t194 = t254 * t351 - t364;
t360 = t292 * t191 - t294 * t192 + t194 * t355;
t356 = MDP(4) * t295 - MDP(9) * pkin(1);
t274 = t346 * t293;
t271 = qJD(1) * t274;
t275 = t346 * t295;
t272 = qJD(1) * t275;
t339 = t290 * t272;
t229 = t271 * t291 + t339;
t326 = qJD(4) - t229;
t354 = qJD(5) + t285;
t325 = MDP(11) + MDP(15);
t353 = pkin(1) * t295 * MDP(10) + (t293 ^ 2 - t295 ^ 2) * MDP(5);
t352 = -t294 * t191 - t292 * t192 - t194 * t363;
t253 = t330 ^ 2;
t349 = pkin(3) * t248;
t348 = pkin(7) * t254;
t347 = pkin(7) * t330;
t230 = -t274 * t291 - t275 * t290;
t344 = t204 * t230;
t338 = t291 * t272;
t212 = t251 * t291 + t290 * t301;
t265 = qJD(2) * pkin(2) + t271;
t225 = t265 * t290 - t338;
t231 = t274 * t290 - t275 * t291;
t327 = -t347 + t326;
t324 = MDP(12) - MDP(17);
t323 = 0.2e1 * qJD(1);
t321 = pkin(2) * t329;
t221 = qJD(2) * qJ(4) + t225;
t283 = -pkin(2) * t291 - pkin(3);
t280 = pkin(2) * t318;
t316 = qJ(4) * t249 - t280;
t211 = t251 * t290 - t291 * t301;
t224 = t265 * t291 + t339;
t228 = t271 * t290 - t338;
t314 = qJD(2) * t229 - t205;
t313 = qJD(4) - t224;
t195 = qJD(2) * t351 + t313 - t347;
t201 = t221 + t348;
t310 = t294 * t195 - t292 * t201;
t309 = -t292 * t195 - t294 * t201;
t266 = t290 * t293 - t337;
t306 = t266 * t294 - t267 * t292;
t227 = t266 * t292 + t267 * t294;
t305 = qJ(4) * t267 + t320;
t304 = -qJ(4) * t254 - t321;
t303 = qJD(4) * t330 + t316;
t302 = qJD(2) * t228 - t204;
t259 = qJD(2) * t337 - t290 * t328;
t300 = -pkin(2) * t328 + qJ(4) * t259 + qJD(4) * t267;
t299 = t204 * t267 + t211 * t330 - t212 * t254 + t230 * t249 - t231 * t248;
t296 = qJD(2) ^ 2;
t281 = pkin(2) * t290 + qJ(4);
t279 = -pkin(4) + t283;
t223 = pkin(3) * t266 - t305;
t215 = -qJD(2) * pkin(3) + t313;
t214 = pkin(7) * t266 + t231;
t213 = -pkin(7) * t267 + t230;
t210 = pkin(3) * t330 - t304;
t209 = pkin(3) * t254 + t364;
t206 = t228 + t348;
t203 = t266 * t351 + t305;
t200 = pkin(3) * t256 - t300;
t198 = pkin(7) * t256 + t212;
t197 = -pkin(7) * t259 + t211;
t196 = t330 * t351 + t304;
t193 = -t303 + t349;
t190 = t256 * t351 + t300;
t189 = qJD(5) * t227 - t256 * t294 + t259 * t292;
t188 = qJD(5) * t306 + t256 * t292 + t259 * t294;
t185 = t248 * t351 + t303;
t1 = [(-t248 * t320 + t256 * t273) * MDP(11) + (-t249 * t320 + t259 * t273) * MDP(12) + (-t205 * t266 - t224 * t259 - t225 * t256 + t299) * MDP(13) + (t205 * t231 - t211 * t224 + t212 * t225 + t344) * MDP(14) + (t193 * t266 + t200 * t254 + t209 * t256 + t223 * t248) * MDP(15) + (-t202 * t266 + t215 * t259 - t221 * t256 + t299) * MDP(16) + (-t193 * t267 - t200 * t330 - t209 * t259 - t223 * t249) * MDP(17) + (t193 * t223 + t200 * t209 + t202 * t231 + t211 * t215 + t212 * t221 + t344) * MDP(18) + (t186 * t227 + t188 * t355) * MDP(19) + (t186 * t306 - t187 * t227 + t188 * t363 - t189 * t355) * MDP(20) + (-t185 * t306 + t203 * t187 + t194 * t189 - t190 * t363) * MDP(24) + (t185 * t227 + t203 * t186 + t194 * t188 + t190 * t355) * MDP(25) + (-t188 * MDP(21) + t189 * MDP(22) + (-t197 * t294 + t198 * t292) * MDP(24) + (t197 * t292 + t198 * t294) * MDP(25) + ((t213 * t292 + t214 * t294) * MDP(24) + (t213 * t294 - t214 * t292) * MDP(25)) * qJD(5)) * t285 + (MDP(6) * t295 - MDP(7) * t293 + (MDP(10) * t293 - MDP(9) * t295) * pkin(6)) * t296 + (-t324 * t212 - t325 * t211 - t353 * t323 + (t356 * t323 + ((qJD(1) * t266 + t254) * MDP(11) + MDP(12) * t366 + (t273 - t362) * MDP(14)) * pkin(2)) * t293) * qJD(2); (-t254 * t321 - t273 * t330 + t302) * MDP(11) + (t254 * t273 - t321 * t330 + t314) * MDP(12) + ((t225 - t228) * t330 + (-t224 + t229) * t254 + (-t248 * t290 - t249 * t291) * pkin(2)) * MDP(13) + (t224 * t228 - t225 * t229 + (-t204 * t291 + t205 * t290 - t273 * t329) * pkin(2)) * MDP(14) + (-t209 * t330 - t210 * t254 + t302) * MDP(15) + (-t248 * t281 + t249 * t283 + (t221 - t228) * t330 + (t215 - t326) * t254) * MDP(16) + (-t209 * t254 + t210 * t330 + 0.2e1 * t286 - t314) * MDP(17) + (t202 * t281 + t204 * t283 - t209 * t210 - t215 * t228 + t221 * t326) * MDP(18) + (t196 * t363 + (t294 * t206 + t292 * t327) * t285 + (-(-t279 * t292 - t281 * t294) * t285 - t309) * qJD(5) + t360) * MDP(24) + (-t196 * t355 + (-t292 * t206 + t294 * t327) * t285 + ((t279 * t294 - t281 * t292) * t285 + t310) * qJD(5) - t352) * MDP(25) + (-t293 * t356 + t353) * qJD(1) ^ 2 - t369; (t224 * t330 + t225 * t254 + t280) * MDP(14) + (t349 + t221 * t254 + (-qJD(4) - t215) * t330 - t316) * MDP(18) + (-t187 + t343) * MDP(24) + (-t186 + t342) * MDP(25) + t324 * (-t276 + (-t254 + t319) * qJD(2)) + (MDP(13) + MDP(16)) * (-t254 ^ 2 - t253) + t325 * qJD(2) * t366; (-t276 + (t254 + t319) * qJD(2)) * MDP(16) + (-t253 - t296) * MDP(17) + (-qJD(2) * t221 + t204) * MDP(18) + (MDP(15) * t254 + MDP(18) * t209 + MDP(24) * t363 - MDP(25) * t355) * t330 + (-MDP(24) * t292 - MDP(25) * t294) * t285 ^ 2; (t309 * t354 - t360) * MDP(24) + (-t310 * t354 + t352) * MDP(25) + t369;];
tauc = t1;
