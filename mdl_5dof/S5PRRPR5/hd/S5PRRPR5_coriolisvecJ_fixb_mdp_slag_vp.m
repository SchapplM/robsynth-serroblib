% Calculate Coriolis joint torque vector for
% S5PRRPR5
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
%   see S5PRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:04:49
% EndTime: 2021-01-15 16:04:56
% DurationCPUTime: 2.04s
% Computational Cost: add. (1458->267), mult. (3940->394), div. (0->0), fcn. (2978->10), ass. (0->133)
t249 = sin(qJ(3));
t252 = cos(qJ(3));
t334 = (t249 ^ 2 - t252 ^ 2) * MDP(6);
t250 = sin(qJ(2));
t246 = sin(pkin(5));
t305 = qJD(1) * t246;
t293 = t250 * t305;
t229 = qJD(2) * pkin(7) + t293;
t247 = cos(pkin(5));
t304 = qJD(1) * t247;
t291 = t249 * t304;
t267 = -t229 * t252 - t291;
t253 = cos(qJ(2));
t292 = t253 * t305;
t274 = qJD(4) + t292;
t288 = qJD(3) * t252 * qJ(4);
t333 = t267 * qJD(3) + (-t274 * t249 - t288) * qJD(2);
t301 = qJD(3) * t249;
t332 = pkin(3) * t301 - t293;
t245 = sin(pkin(10));
t326 = cos(pkin(10));
t284 = t326 * t249;
t227 = t245 * t252 + t284;
t222 = t227 * qJD(2);
t235 = t252 * t304;
t183 = (-t229 * t249 + t235) * qJD(3) + (-qJ(4) * t301 + t274 * t252) * qJD(2);
t169 = t183 * t245 - t326 * t333;
t283 = t326 * t252;
t234 = qJD(2) * t283;
t303 = qJD(2) * t249;
t219 = t245 * t303 - t234;
t216 = qJD(5) + t219;
t238 = pkin(3) * t245 + pkin(8);
t296 = pkin(3) * t303;
t331 = (pkin(4) * t222 + pkin(8) * t219 + qJD(5) * t238 + t296) * t216 + t169;
t285 = t326 * t183;
t170 = t333 * t245 + t285;
t281 = qJ(4) * qJD(2) + t229;
t200 = -t281 * t249 + t235;
t198 = qJD(3) * pkin(3) + t200;
t201 = t281 * t252 + t291;
t317 = t245 * t201;
t173 = t326 * t198 - t317;
t171 = -qJD(3) * pkin(4) - t173;
t241 = -pkin(3) * t252 - pkin(2);
t215 = t241 * qJD(2) + qJD(4) - t292;
t182 = pkin(4) * t219 - pkin(8) * t222 + t215;
t316 = t245 * t249;
t264 = t283 - t316;
t194 = -pkin(4) * t264 - pkin(8) * t227 + t241;
t224 = t264 * qJD(3);
t328 = qJ(4) + pkin(7);
t231 = t328 * t252;
t205 = t326 * t231 - t328 * t316;
t221 = t227 * qJD(3);
t213 = qJD(2) * t221;
t273 = t169 * t227 - t205 * t213;
t286 = qJD(3) * t328;
t217 = qJD(4) * t252 - t249 * t286;
t263 = -qJD(4) * t249 - t252 * t286;
t308 = t326 * t217 + t245 * t263 - t264 * t292;
t330 = (qJD(5) * t182 + t170) * t264 + t171 * t224 + (-qJD(5) * t194 - t308) * t216 + t273;
t329 = pkin(3) * t249;
t327 = qJD(2) * pkin(2);
t248 = sin(qJ(5));
t299 = qJD(5) * t248;
t297 = qJD(2) * qJD(3);
t287 = t249 * t297;
t232 = t245 * t287;
t214 = qJD(3) * t234 - t232;
t251 = cos(qJ(5));
t298 = t251 * qJD(3);
t307 = qJD(5) * t298 + t251 * t214;
t184 = -t222 * t299 + t307;
t325 = t184 * t248;
t324 = t194 * t213;
t318 = t222 * t248;
t206 = -t298 + t318;
t323 = t206 * t216;
t322 = t206 * t222;
t208 = qJD(3) * t248 + t222 * t251;
t321 = t208 * t216;
t320 = t208 * t222;
t319 = t214 * t248;
t315 = t246 * t250;
t314 = t246 * t253;
t255 = qJD(2) ^ 2;
t313 = t246 * t255;
t312 = t248 * t213;
t254 = qJD(3) ^ 2;
t311 = t249 * t254;
t209 = t251 * t213;
t310 = t252 * t254;
t309 = t217 * t245 - t227 * t292 - t326 * t263;
t196 = t326 * t201;
t174 = t245 * t198 + t196;
t302 = qJD(2) * t250;
t290 = t246 * t302;
t218 = pkin(3) * t287 + qJD(1) * t290;
t300 = qJD(5) * t227;
t294 = t250 * t313;
t289 = qJD(2) * t314;
t282 = t251 * t216;
t278 = t252 * t289;
t277 = t249 * t289;
t276 = pkin(4) * t221 - pkin(8) * t224 + t332;
t172 = qJD(3) * pkin(8) + t174;
t168 = t172 * t251 + t182 * t248;
t272 = t172 * t248 - t182 * t251;
t271 = t209 + (-t219 * t248 - t299) * t216;
t225 = t247 * t249 + t252 * t315;
t268 = t247 * t252 - t249 * t315;
t191 = t326 * t225 + t245 * t268;
t270 = -t191 * t248 - t251 * t314;
t269 = -t191 * t251 + t248 * t314;
t266 = t224 * t251 - t227 * t299;
t265 = qJD(2) * t327;
t262 = t225 * qJD(3);
t260 = -0.2e1 * qJD(3) * t327;
t178 = t326 * t200 - t317;
t259 = -t238 * t213 + (t171 + t178) * t216;
t257 = -t262 - t277;
t239 = -t326 * pkin(3) - pkin(4);
t204 = t231 * t245 + t328 * t284;
t199 = t268 * qJD(3) + t278;
t190 = t225 * t245 - t326 * t268;
t185 = qJD(5) * t208 + t319;
t180 = pkin(4) * t213 - pkin(8) * t214 + t218;
t179 = t251 * t180;
t177 = t326 * t199 + t245 * t257;
t176 = t200 * t245 + t196;
t175 = t199 * t245 - t326 * t257;
t1 = [-MDP(3) * t294 - t253 * MDP(4) * t313 + (-t252 * t294 + (-t262 - 0.2e1 * t277) * qJD(3)) * MDP(10) + (t249 * t294 + (-t199 - t278) * qJD(3)) * MDP(11) + (-qJD(3) * t175 + (-t213 * t253 + t219 * t302) * t246) * MDP(12) + (-qJD(3) * t177 + (-t214 * t253 + t222 * t302) * t246) * MDP(13) + (t175 * t222 - t177 * t219 + t190 * t214 - t191 * t213) * MDP(14) + (t169 * t190 + t170 * t191 - t173 * t175 + t174 * t177 + (t215 * t302 - t218 * t253) * t246) * MDP(15) + ((t269 * qJD(5) - t177 * t248 + t251 * t290) * t216 + t270 * t213 + t175 * t206 + t190 * t185) * MDP(21) + (-(t270 * qJD(5) + t177 * t251 + t248 * t290) * t216 + t269 * t213 + t175 * t208 + t190 * t184) * MDP(22); 0.2e1 * t252 * MDP(5) * t287 - 0.2e1 * t297 * t334 + MDP(7) * t310 - MDP(8) * t311 + (-pkin(7) * t310 + t249 * t260) * MDP(10) + (pkin(7) * t311 + t252 * t260) * MDP(11) + (-t219 * t293 + t213 * t241 + t215 * t221 - t218 * t264 + (t219 * t329 - t309) * qJD(3)) * MDP(12) + (-t222 * t293 + t214 * t241 + t215 * t224 + t218 * t227 + (t222 * t329 - t308) * qJD(3)) * MDP(13) + (t170 * t264 - t173 * t224 - t174 * t221 + t204 * t214 - t308 * t219 + t309 * t222 + t273) * MDP(14) + (t169 * t204 + t170 * t205 - t309 * t173 + t308 * t174 + t332 * t215 + t218 * t241) * MDP(15) + (t184 * t227 * t251 + t266 * t208) * MDP(16) + ((-t206 * t251 - t208 * t248) * t224 + (-t325 - t185 * t251 + (t206 * t248 - t208 * t251) * qJD(5)) * t227) * MDP(17) + (-t184 * t264 + t208 * t221 + t227 * t209 + t266 * t216) * MDP(18) + (-t227 * t312 + t185 * t264 - t206 * t221 + (-t224 * t248 - t251 * t300) * t216) * MDP(19) + (-t213 * t264 + t216 * t221) * MDP(20) + (-t272 * t221 - t179 * t264 + t204 * t185 + t309 * t206 + (t324 + t276 * t216 + (t171 * t227 + t172 * t264 - t205 * t216) * qJD(5)) * t251 + t330 * t248) * MDP(21) + (-t168 * t221 + t204 * t184 + t309 * t208 + (-t324 + (-qJD(5) * t172 + t180) * t264 - t171 * t300 + (qJD(5) * t205 - t276) * t216) * t248 + t330 * t251) * MDP(22); t255 * t334 + t249 * MDP(10) * t265 + (qJD(3) * t176 - t215 * t222 - t219 * t296 - t169) * MDP(12) + (-t285 + t215 * t219 + (-t245 * t267 + t178) * qJD(3) + (t245 * t288 + (-pkin(3) * t222 + t245 * t274) * t249) * qJD(2)) * MDP(13) + ((t174 - t176) * t222 + (-t173 + t178) * t219 + (-t213 * t245 - t326 * t214) * pkin(3)) * MDP(14) + (t173 * t176 - t174 * t178 + (-t326 * t169 + t170 * t245 - t215 * t303) * pkin(3)) * MDP(15) + (t208 * t282 + t325) * MDP(16) + ((t184 - t323) * t251 + (-t185 - t321) * t248) * MDP(17) + (t216 * t282 + t312 - t320) * MDP(18) + (t271 + t322) * MDP(19) - t216 * t222 * MDP(20) + (-t176 * t206 + t239 * t185 + t222 * t272 + t259 * t248 - t331 * t251) * MDP(21) + (t168 * t222 - t176 * t208 + t239 * t184 + t331 * t248 + t259 * t251) * MDP(22) + (-t249 * t255 * MDP(5) + t265 * MDP(11)) * t252; 0.2e1 * t222 * qJD(3) * MDP(12) + (-t232 + (t234 - t219) * qJD(3)) * MDP(13) + (-t219 ^ 2 - t222 ^ 2) * MDP(14) + (t173 * t222 + t174 * t219 + t218) * MDP(15) + (t271 - t322) * MDP(21) + (-t216 ^ 2 * t251 - t312 - t320) * MDP(22); t208 * t206 * MDP(16) + (-t206 ^ 2 + t208 ^ 2) * MDP(17) + (t307 + t323) * MDP(18) + (-t319 + t321) * MDP(19) + t213 * MDP(20) + (t168 * t216 - t170 * t248 - t171 * t208 + t179) * MDP(21) + (-t170 * t251 + t171 * t206 - t180 * t248 - t216 * t272) * MDP(22) + (-MDP(18) * t318 - t208 * MDP(19) - t168 * MDP(21) + t272 * MDP(22)) * qJD(5);];
tauc = t1;
