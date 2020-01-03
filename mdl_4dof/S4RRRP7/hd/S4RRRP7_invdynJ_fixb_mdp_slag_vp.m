% Calculate vector of inverse dynamics joint torques for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_invdynJ_fixb_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRRP7_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:10
% EndTime: 2019-12-31 17:21:14
% DurationCPUTime: 3.32s
% Computational Cost: add. (1379->321), mult. (3073->414), div. (0->0), fcn. (1869->6), ass. (0->127)
t245 = cos(qJ(2));
t297 = qJD(1) * t245;
t333 = qJD(3) - t297;
t242 = sin(qJ(2));
t290 = qJD(3) * t242;
t332 = qJD(1) * t290 - qJDD(2);
t241 = sin(qJ(3));
t244 = cos(qJ(3));
t284 = qJDD(1) * t242;
t179 = t241 * (qJD(2) * (qJD(3) + t297) + t284) + t332 * t244;
t243 = sin(qJ(1));
t246 = cos(qJ(1));
t271 = g(1) * t246 + g(2) * t243;
t331 = t242 * t271;
t237 = t245 * qJDD(1);
t285 = qJD(1) * qJD(2);
t330 = -t242 * t285 + t237;
t287 = t244 * qJD(2);
t298 = qJD(1) * t242;
t210 = t241 * t298 - t287;
t294 = qJD(2) * t241;
t212 = t244 * t298 + t294;
t220 = -qJD(2) * pkin(2) + pkin(5) * t298;
t177 = pkin(3) * t210 - qJ(4) * t212 + t220;
t208 = qJDD(3) - t330;
t325 = pkin(6) * t208;
t329 = -t177 * t333 + t325;
t327 = t212 ^ 2;
t326 = pkin(3) * t208;
t322 = g(3) * t242;
t321 = g(3) * t245;
t320 = pkin(6) * qJD(3);
t319 = qJ(4) * t208;
t265 = pkin(2) * t245 + pkin(6) * t242 + pkin(1);
t205 = t265 * qJD(1);
t234 = pkin(5) * t297;
t221 = qJD(2) * pkin(6) + t234;
t183 = -t205 * t241 + t221 * t244;
t176 = qJ(4) * t333 + t183;
t318 = t176 * t333;
t278 = t245 * t285;
t178 = -qJD(3) * t287 + (-t278 - t284) * t244 + t332 * t241;
t317 = t178 * t241;
t316 = t183 * t333;
t315 = t210 * t212;
t314 = t210 * t333;
t313 = t212 * t333;
t312 = t212 * t244;
t311 = t241 * t242;
t310 = t241 * t245;
t309 = t242 * t244;
t308 = t243 * t244;
t307 = t244 * t208;
t306 = t244 * t245;
t305 = t244 * t246;
t304 = t246 * t241;
t268 = pkin(3) * t241 - qJ(4) * t244;
t303 = -qJD(4) * t241 + t333 * t268 - t234;
t274 = pkin(2) * t242 - pkin(6) * t245;
t215 = t274 * qJD(2);
t289 = qJD(3) * t244;
t302 = t241 * t215 - t265 * t289;
t301 = (g(1) * t305 + g(2) * t308) * t242;
t300 = pkin(5) * t306 - t241 * t265;
t239 = t242 ^ 2;
t299 = -t245 ^ 2 + t239;
t296 = qJD(2) * t210;
t295 = qJD(2) * t212;
t293 = qJD(2) * t242;
t292 = qJD(2) * t245;
t291 = qJD(3) * t241;
t288 = qJD(3) * t245;
t182 = -t205 * t244 - t221 * t241;
t286 = qJD(4) - t182;
t283 = pkin(5) * t241 + pkin(3);
t282 = t333 * t294;
t281 = t333 * t287;
t280 = t333 * t291;
t187 = qJD(1) * t215 - qJDD(1) * t265;
t197 = pkin(5) * t330 + qJDD(2) * pkin(6);
t275 = -t244 * t187 + t241 * t197 - t205 * t291 + t221 * t289;
t199 = t243 * t310 + t305;
t201 = t245 * t304 - t308;
t273 = -g(1) * t199 + g(2) * t201;
t200 = t243 * t306 - t304;
t202 = t241 * t243 + t245 * t305;
t272 = g(1) * t200 - g(2) * t202;
t270 = g(1) * t243 - g(2) * t246;
t232 = pkin(5) * t284;
t198 = -qJDD(2) * pkin(2) + pkin(5) * t278 + t232;
t269 = pkin(3) * t244 + qJ(4) * t241;
t267 = qJD(3) * t220 - t325;
t175 = -pkin(3) * t333 + t286;
t266 = t175 * t244 - t176 * t241;
t264 = pkin(2) + t269;
t263 = pkin(5) + t268;
t262 = t320 * t333 + t321;
t260 = -0.2e1 * pkin(1) * t285 - pkin(5) * qJDD(2);
t259 = -t208 * t241 - t289 * t333;
t258 = -t280 + t307;
t257 = t215 * t244 + t265 * t291;
t169 = pkin(3) * t179 + qJ(4) * t178 - qJD(4) * t212 + t198;
t256 = -t169 - t262;
t255 = t241 * t187 + t244 * t197 - t205 * t289 - t221 * t291;
t248 = qJD(1) ^ 2;
t254 = pkin(1) * t248 + t271;
t253 = -t245 * t271 - t322;
t247 = qJD(2) ^ 2;
t252 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t247 + t270;
t251 = g(1) * t201 + g(2) * t199 + g(3) * t311 - t275;
t250 = t177 * t212 + qJDD(4) - t251;
t249 = -g(1) * t202 - g(2) * t200 - g(3) * t309 + t255;
t214 = t274 * qJD(1);
t203 = t241 * t214;
t195 = t263 * t242;
t189 = t244 * t265 + t245 * t283;
t188 = -qJ(4) * t245 + t300;
t186 = pkin(3) * t212 + qJ(4) * t210;
t185 = -t214 * t244 - t283 * t298;
t184 = t203 + (-pkin(5) * t244 + qJ(4)) * t298;
t174 = (qJD(3) * t269 - qJD(4) * t244) * t242 + t263 * t292;
t173 = -pkin(3) * t293 + (-t241 * t293 + t244 * t288) * pkin(5) - t257;
t172 = -t178 + t314;
t171 = qJ(4) * t293 - qJD(4) * t245 + (-t241 * t288 - t242 * t287) * pkin(5) + t302;
t170 = qJDD(4) + t275 - t326;
t168 = qJD(4) * t333 + t255 + t319;
t1 = [qJDD(1) * MDP(1) + t270 * MDP(2) + t271 * MDP(3) + (qJDD(1) * t239 + 0.2e1 * t242 * t278) * MDP(4) + 0.2e1 * (t237 * t242 - t285 * t299) * MDP(5) + (qJDD(2) * t242 + t245 * t247) * MDP(6) + (qJDD(2) * t245 - t242 * t247) * MDP(7) + (t242 * t260 + t245 * t252) * MDP(9) + (-t242 * t252 + t245 * t260) * MDP(10) + (-t178 * t309 + (-t241 * t290 + t245 * t287) * t212) * MDP(11) + ((-t210 * t244 - t212 * t241) * t292 + (t317 - t179 * t244 + (t210 * t241 - t312) * qJD(3)) * t242) * MDP(12) + ((t178 + t281) * t245 + (t258 + t295) * t242) * MDP(13) + ((t179 - t282) * t245 + (t259 - t296) * t242) * MDP(14) + (-t208 * t245 + t293 * t333) * MDP(15) + (t257 * t333 - t265 * t307 + (t220 * t294 + (t259 + t296) * pkin(5) + t275) * t245 + (t220 * t289 + t182 * qJD(2) + t198 * t241 + (t179 + t282) * pkin(5)) * t242 + t272) * MDP(16) + (-t302 * t333 - t300 * t208 + (t220 * t287 + (t280 + t295) * pkin(5) + t255) * t245 + (-t220 * t291 - t183 * qJD(2) + t198 * t244 + (-t178 + t281) * pkin(5)) * t242 + t273) * MDP(17) + (-t173 * t333 + t174 * t210 + t179 * t195 - t189 * t208 + (t177 * t294 + t170) * t245 + (-qJD(2) * t175 + t169 * t241 + t177 * t289) * t242 + t272) * MDP(18) + (-t171 * t210 + t173 * t212 - t178 * t189 - t179 * t188 + t266 * t292 + (-t168 * t241 + t170 * t244 + (-t175 * t241 - t176 * t244) * qJD(3) + t270) * t242) * MDP(19) + (t171 * t333 - t174 * t212 + t178 * t195 + t188 * t208 + (-t177 * t287 - t168) * t245 + (qJD(2) * t176 - t169 * t244 + t177 * t291) * t242 - t273) * MDP(20) + (t168 * t188 + t176 * t171 + t169 * t195 + t177 * t174 + t170 * t189 + t175 * t173 - g(1) * (-pkin(3) * t200 - qJ(4) * t199) - g(2) * (pkin(3) * t202 + qJ(4) * t201) + (-g(1) * pkin(5) - g(2) * t265) * t246 + (-g(2) * pkin(5) + g(1) * t265) * t243) * MDP(21); MDP(6) * t284 + MDP(7) * t237 + qJDD(2) * MDP(8) + (t242 * t254 - t232 - t321) * MDP(9) + (t322 + (-pkin(5) * qJDD(1) + t254) * t245) * MDP(10) + (t312 * t333 - t317) * MDP(11) + ((-t178 - t314) * t244 + (-t179 - t313) * t241) * MDP(12) + ((-t212 * t242 - t306 * t333) * qJD(1) - t259) * MDP(13) + ((t210 * t242 + t310 * t333) * qJD(1) + t258) * MDP(14) - t333 * MDP(15) * t298 + (-pkin(2) * t179 + t267 * t241 + (-t321 - t198 - (t214 + t320) * t333) * t244 + (-t220 * t310 - t182 * t242 + (-t210 * t245 - t311 * t333) * pkin(5)) * qJD(1) + t301) * MDP(16) + (pkin(2) * t178 + t203 * t333 + t267 * t244 + (-t220 * t306 + t183 * t242 + (-t212 * t245 - t309 * t333) * pkin(5)) * qJD(1) + (t198 + t262 - t331) * t241) * MDP(17) + (t175 * t298 - t179 * t264 + t185 * t333 + t303 * t210 - t241 * t329 + t256 * t244 + t301) * MDP(18) + (t184 * t210 - t185 * t212 + (t168 + t333 * t175 + (qJD(3) * t212 - t179) * pkin(6)) * t244 + (t170 - t318 + (qJD(3) * t210 - t178) * pkin(6)) * t241 + t253) * MDP(19) + (-t176 * t298 - t178 * t264 - t184 * t333 - t303 * t212 + t329 * t244 + (t256 + t331) * t241) * MDP(20) + (-t175 * t185 - t176 * t184 + t303 * t177 + (qJD(3) * t266 + t168 * t244 + t170 * t241 + t253) * pkin(6) + (-t169 - t321 + t331) * t264) * MDP(21) + (-MDP(4) * t242 * t245 + MDP(5) * t299) * t248; MDP(11) * t315 + (-t210 ^ 2 + t327) * MDP(12) + t172 * MDP(13) + (-t179 + t313) * MDP(14) + t208 * MDP(15) + (-t212 * t220 + t251 + t316) * MDP(16) + (t182 * t333 + t210 * t220 - t249) * MDP(17) + (-t186 * t210 - t250 + t316 + 0.2e1 * t326) * MDP(18) + (pkin(3) * t178 - qJ(4) * t179 + (t176 - t183) * t212 + (t175 - t286) * t210) * MDP(19) + (0.2e1 * t319 - t177 * t210 + t186 * t212 - (-0.2e1 * qJD(4) + t182) * t333 + t249) * MDP(20) + (t168 * qJ(4) - t170 * pkin(3) - t177 * t186 - t175 * t183 - g(1) * (-pkin(3) * t201 + qJ(4) * t202) - g(2) * (-pkin(3) * t199 + qJ(4) * t200) + t268 * t322 + t286 * t176) * MDP(21); (-t208 + t315) * MDP(18) + t172 * MDP(19) + (-t333 ^ 2 - t327) * MDP(20) + (t250 - t318 - t326) * MDP(21);];
tau = t1;
