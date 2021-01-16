% Calculate Coriolis joint torque vector for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:13:25
% EndTime: 2021-01-15 21:13:34
% DurationCPUTime: 2.66s
% Computational Cost: add. (1391->248), mult. (3466->359), div. (0->0), fcn. (2303->6), ass. (0->127)
t249 = qJD(2) + qJD(4);
t254 = sin(qJ(4));
t257 = cos(qJ(2));
t329 = cos(qJ(4));
t289 = qJD(1) * t329;
t255 = sin(qJ(2));
t307 = qJD(1) * t255;
t334 = -t254 * t307 + t257 * t289;
t335 = t334 * t249;
t224 = t254 * t257 + t255 * t329;
t200 = t249 * t224;
t194 = t200 * qJD(1);
t258 = pkin(1) + pkin(2);
t216 = qJD(5) - t334;
t333 = t216 ^ 2;
t250 = t255 ^ 2;
t251 = t257 ^ 2;
t332 = (t250 - t251) * MDP(5);
t328 = pkin(3) + qJ(3);
t230 = t328 * t255;
t225 = qJD(1) * t230;
t252 = qJD(2) * pkin(1);
t212 = qJD(2) * pkin(2) - t225 + t252;
t286 = qJD(2) * t328;
t303 = qJD(3) * t255;
t214 = -t257 * t286 - t303;
t210 = t214 * qJD(1);
t213 = qJD(3) * t257 - t255 * t286;
t211 = t213 * qJD(1);
t231 = t328 * t257;
t226 = qJD(1) * t231;
t302 = qJD(4) * t254;
t265 = t254 * t210 + t211 * t329 - t226 * t302;
t288 = qJD(4) * t329;
t172 = t212 * t288 + t265;
t293 = t329 * t226;
t192 = t254 * t212 + t293;
t284 = -t329 * t210 + t254 * t211;
t173 = qJD(4) * t192 + t284;
t268 = -t230 * t329 - t254 * t231;
t176 = qJD(4) * t268 + t213 * t329 + t254 * t214;
t294 = t329 * t212;
t314 = t254 * t226;
t191 = -t294 + t314;
t306 = qJD(1) * t257;
t219 = -t254 * t306 - t255 * t289;
t232 = t258 * t257;
t222 = -qJD(1) * t232 + qJD(3);
t198 = pkin(4) * t219 + t222;
t267 = -t254 * t255 + t257 * t329;
t199 = t249 * t267;
t206 = -t254 * t230 + t231 * t329;
t208 = -pkin(4) * t224 - t232;
t331 = t173 * t224 + t191 * t199 - t206 * t194 - (qJD(5) * t208 + t176) * t216 + (qJD(5) * t198 + t172) * t267;
t330 = 0.2e1 * t251;
t253 = sin(qJ(5));
t301 = qJD(5) * t253;
t256 = cos(qJ(5));
t300 = qJD(5) * t256;
t309 = t249 * t300 + t256 * t335;
t178 = t219 * t301 + t309;
t327 = t178 * t253;
t326 = t191 * t224;
t325 = t335 * t253;
t317 = t219 * t253;
t202 = -t256 * t249 - t317;
t324 = t202 * t216;
t323 = t202 * t219;
t273 = t219 * t256 - t249 * t253;
t322 = t273 * t216;
t321 = t273 * t219;
t320 = t208 * t194;
t318 = t219 * t249;
t316 = t222 * t334;
t315 = t253 * t194;
t259 = qJD(2) ^ 2;
t313 = t255 * t259;
t260 = qJD(1) ^ 2;
t312 = t255 * t260;
t190 = t256 * t194;
t311 = t257 * t259;
t310 = t257 * t260;
t298 = qJD(1) * qJD(2);
t287 = t255 * t298;
t242 = pkin(1) * t287;
t220 = pkin(2) * t287 + t242;
t227 = t258 * t307;
t305 = qJD(2) * t255;
t228 = t258 * t305;
t304 = qJD(2) * t257;
t297 = pkin(1) * t306;
t236 = qJD(3) - t297;
t299 = -qJD(3) + t236;
t296 = 0.2e1 * t298;
t295 = t258 * t329;
t292 = qJ(3) * t304;
t229 = -qJ(3) * t307 + t252;
t290 = t229 * t304;
t285 = t191 * (-t216 - t334);
t283 = t256 * t216;
t237 = t254 * t258 + pkin(4);
t280 = -pkin(4) * t334 + qJD(5) * t237 + t227;
t279 = (-qJD(3) - t236) * qJD(1);
t278 = t257 * t296;
t186 = t249 * pkin(4) + t192;
t175 = t186 * t256 + t198 * t253;
t277 = t173 * t253 - t175 * t219 + t191 * t300;
t196 = -t254 * t225 + t293;
t276 = t258 * t302 - t196;
t275 = t186 * t253 - t198 * t256;
t274 = -t191 * t334 - t194 * t237;
t272 = t190 + (t253 * t334 - t301) * t216;
t197 = -t225 * t329 - t314;
t271 = -t258 * t288 + t197;
t270 = -t173 * t256 + t191 * t301 - t219 * t275;
t269 = t222 * t219 - t284;
t266 = t199 * t256 - t224 * t301;
t263 = -t333 * t256 - t315;
t179 = -t273 * qJD(5) + t325;
t262 = ((t178 - t324) * t256 + (-t179 + t322) * t253) * MDP(23) + (-t273 * t283 + t327) * MDP(22) + (t272 - t323) * MDP(25) + (t216 * t283 + t315 - t321) * MDP(24) + (-t318 - t194) * MDP(18) + (t219 ^ 2 - t334 ^ 2) * MDP(16) + (MDP(15) * t334 + t216 * MDP(26)) * t219;
t261 = qJ(3) ^ 2;
t221 = (-t292 - t303) * qJD(1);
t185 = -pkin(4) * t199 + t228;
t183 = -pkin(4) * t335 + t220;
t182 = t256 * t183;
t177 = qJD(4) * t206 + t254 * t213 - t214 * t329;
t1 = [t255 * MDP(4) * t278 - t296 * t332 + MDP(6) * t311 - MDP(7) * t313 + (-qJ(3) * t311 + (-0.3e1 * t297 + t299) * t305) * MDP(11) + (qJ(3) * t313 + (t299 * t257 + (0.2e1 * t250 - t251) * qJD(1) * pkin(1)) * qJD(2)) * MDP(12) + (-t290 - t221 * t255 + (-0.2e1 * t255 * t292 + (t250 + t330) * qJD(3)) * qJD(1)) * MDP(13) + ((qJD(1) * qJD(3) * t330 - t290) * qJ(3) + (-qJ(3) * t221 - qJD(3) * t229 + (-0.2e1 * t261 * t306 + (t236 - t297) * pkin(1)) * qJD(2)) * t255) * MDP(14) + (-t199 * t219 + t224 * t335) * MDP(15) + (-t194 * t224 + t199 * t334 + t200 * t219 + t267 * t335) * MDP(16) + (-t194 * t232 + t200 * t222 - t220 * t267 - t228 * t334) * MDP(20) + (t199 * t222 - t219 * t228 + t220 * t224 - t232 * t335) * MDP(21) + (t178 * t224 * t256 - t266 * t273) * MDP(22) + ((-t202 * t256 + t253 * t273) * t199 + (-t327 - t179 * t256 + (t202 * t253 + t256 * t273) * qJD(5)) * t224) * MDP(23) + (-t178 * t267 + t190 * t224 - t200 * t273 + t216 * t266) * MDP(24) + (-t224 * t315 + t179 * t267 - t200 * t202 + (-t199 * t253 - t224 * t300) * t216) * MDP(25) + (-t194 * t267 + t200 * t216) * MDP(26) + (-t275 * t200 + t177 * t202 - t268 * t179 - t182 * t267 + (t185 * t216 + t320 + (t186 * t267 - t206 * t216 + t326) * qJD(5)) * t256 + t331 * t253) * MDP(27) + (-t175 * t200 - t177 * t273 - t268 * t178 + (-(-qJD(5) * t206 + t185) * t216 - t320 + (-qJD(5) * t186 + t183) * t267 - qJD(5) * t326) * t253 + t331 * t256) * MDP(28) + (t199 * MDP(17) - t200 * MDP(18) - t177 * MDP(20) - t176 * MDP(21)) * t249; (-pkin(1) * t250 * t260 + t257 * t279) * MDP(12) + t260 * t332 + (-t179 * t295 + t274 * t253 + t276 * t202 + (t253 * t271 - t256 * t280) * t216 + t270) * MDP(27) + (-t178 * t295 + t274 * t256 - t276 * t273 + (t253 * t280 + t256 * t271) * t216 + t277) * MDP(28) + t262 + (qJ(3) * t312 + (t229 - t252) * qJD(1)) * t257 * MDP(13) + ((qJ(3) * qJD(1) * t229 + t261 * t312) * t257 + (-t236 * t307 + t221) * pkin(1)) * MDP(14) + (t197 * t249 - t316 + t227 * t219 + (-t249 * t295 - t294) * qJD(4) - t265) * MDP(21) + (t196 * t249 + t227 * t334 + (-t293 + (-t249 * t258 - t212) * t254) * qJD(4) + t269) * MDP(20) + ((pkin(1) * t310 + t279) * MDP(11) - MDP(4) * t310) * t255; 0.2e1 * MDP(11) * t287 + MDP(12) * t278 + (-t250 - t251) * t260 * MDP(13) + (-qJ(3) * t251 * t260 + t229 * t307 + t242) * MDP(14) + (t194 - t318) * MDP(20) + 0.2e1 * t335 * MDP(21) + (t272 + t323) * MDP(27) + (t263 - t321) * MDP(28); (t269 + (-qJD(4) + t249) * t192) * MDP(20) + (-t191 * t249 - t172 - t316) * MDP(21) + (pkin(4) * t263 - t192 * t202 + t253 * t285 + t270) * MDP(27) + (t192 * t273 + t256 * t285 + (t333 * t253 - t190) * pkin(4) + t277) * MDP(28) + t262; -t273 * t202 * MDP(22) + (-t202 ^ 2 + t273 ^ 2) * MDP(23) + (t309 + t324) * MDP(24) + (-t322 - t325) * MDP(25) + t194 * MDP(26) + (-t172 * t253 + t175 * t216 + t191 * t273 + t182) * MDP(27) + (-t172 * t256 - t183 * t253 + t191 * t202 - t216 * t275) * MDP(28) + (MDP(24) * t317 + MDP(25) * t273 - MDP(27) * t175 + MDP(28) * t275) * qJD(5);];
tauc = t1;
