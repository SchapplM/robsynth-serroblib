% Calculate vector of inverse dynamics joint torques for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRRP8_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:34
% EndTime: 2019-12-31 18:47:38
% DurationCPUTime: 1.80s
% Computational Cost: add. (1486->229), mult. (2129->265), div. (0->0), fcn. (1102->6), ass. (0->114)
t232 = sin(qJ(4));
t230 = t232 ^ 2;
t234 = cos(qJ(4));
t231 = t234 ^ 2;
t289 = t230 + t231;
t335 = -t289 * MDP(18) + MDP(9);
t233 = sin(qJ(3));
t235 = cos(qJ(3));
t321 = sin(qJ(1));
t322 = cos(qJ(1));
t199 = -t233 * t321 - t235 * t322;
t200 = t233 * t322 - t235 * t321;
t261 = g(1) * t199 + g(2) * t200;
t276 = qJD(1) - qJD(3);
t334 = t276 ^ 2;
t228 = qJDD(1) - qJDD(3);
t259 = -pkin(4) * t234 - qJ(5) * t232;
t255 = pkin(3) - t259;
t333 = t255 * t228;
t332 = t255 * t276;
t236 = -pkin(1) - pkin(2);
t208 = qJD(1) * t236 + qJD(2);
t288 = qJ(2) * qJD(1);
t191 = t208 * t235 - t233 * t288;
t177 = -t191 + t332;
t331 = t276 * t177;
t330 = -qJD(3) * t288 + qJDD(1) * t236 + qJDD(2);
t329 = qJD(4) * t276;
t280 = qJD(1) * qJD(2);
t281 = qJ(2) * qJDD(1);
t327 = qJD(3) * t208 + t280 + t281;
t268 = -qJ(2) * t233 + t235 * t236;
t293 = qJ(2) * t235 + t233 * t236;
t192 = t208 * t233 + t235 * t288;
t185 = -pkin(7) * t276 + t192;
t306 = t185 * t232;
t267 = qJD(5) + t306;
t314 = qJD(4) * pkin(4);
t179 = t267 - t314;
t286 = qJD(4) * qJ(5);
t305 = t185 * t234;
t180 = t286 + t305;
t250 = t233 * t330 + t235 * t327;
t175 = -pkin(7) * t228 + t250;
t174 = t234 * t175;
t279 = qJDD(4) * qJ(5);
t169 = t279 + t174 + (qJD(5) - t306) * qJD(4);
t173 = t232 * t175;
t284 = qJD(4) * t234;
t311 = qJDD(4) * pkin(4);
t170 = t185 * t284 + qJDD(5) + t173 - t311;
t257 = t169 * t234 + t170 * t232;
t244 = (t179 * t234 - t180 * t232) * qJD(4) + t257;
t239 = t244 + t261;
t326 = t179 * t232 + t180 * t234;
t262 = g(1) * t200 - g(2) * t199;
t325 = g(3) * t234 - t232 * t261;
t324 = -qJDD(4) * t233 + 0.2e1 * t235 * t329;
t258 = pkin(4) * t232 - qJ(5) * t234;
t193 = qJD(4) * t258 - qJD(5) * t232;
t187 = t255 - t268;
t189 = qJD(2) * t235 + qJD(3) * t268;
t202 = -pkin(7) + t293;
t278 = qJDD(4) * t202;
t323 = (-t187 * t276 - t177 - t189) * qJD(4) - t278;
t320 = pkin(3) * t228;
t319 = pkin(3) * t276;
t313 = pkin(1) * qJDD(1);
t312 = pkin(7) * qJDD(4);
t304 = t189 * t276;
t190 = qJD(2) * t233 + qJD(3) * t293;
t303 = t190 * t276;
t302 = t191 * t276;
t301 = t192 * t276;
t300 = t193 * t276;
t299 = t228 * t232;
t298 = t276 * t232;
t297 = t232 * t234;
t285 = qJD(4) * t232;
t295 = t191 * t285 - t234 * t301;
t294 = t193 - t192;
t292 = pkin(1) * t322 + qJ(2) * t321;
t291 = g(1) * t321 - g(2) * t322;
t290 = t230 - t231;
t275 = 0.2e1 * t280;
t274 = t334 * t297;
t184 = -t191 + t319;
t271 = t184 + t319;
t270 = t289 * t228;
t269 = t177 + t332;
t265 = qJDD(2) - t313;
t264 = -t173 + t325;
t263 = -pkin(1) * t321 + qJ(2) * t322;
t254 = t233 * t327 - t235 * t330;
t253 = g(1) * t322 + g(2) * t321;
t237 = qJD(4) ^ 2;
t252 = pkin(7) * t237 - t262;
t251 = -t202 * t237 - t262;
t201 = pkin(3) - t268;
t249 = -t278 + (-t201 * t276 - t184 - t189) * qJD(4);
t248 = -t228 * t235 + (-t237 - t334) * t233;
t176 = t254 + t320;
t247 = -t176 - t252 - t320;
t168 = t254 - t300 + t333;
t246 = -t168 - t252 - t333;
t245 = t254 - t262;
t178 = t190 - t193;
t243 = t178 * t276 + t187 * t228 + t168 + t251;
t242 = -t201 * t228 - t176 - t251 - t303;
t241 = (-t228 * t230 - 0.2e1 * t284 * t298) * MDP(10) + 0.2e1 * (-t228 * t297 + t290 * t329) * MDP(11) + (qJDD(4) * t232 + t234 * t237) * MDP(12) + (qJDD(4) * t234 - t232 * t237) * MDP(13) - t228 * MDP(7);
t240 = t250 + t261;
t238 = qJD(1) ^ 2;
t196 = t258 * t276;
t1 = [(t232 * t242 + t234 * t249) * MDP(16) + t291 * MDP(2) + (-qJDD(2) + t291 + 0.2e1 * t313) * MDP(4) + (t168 * t187 + t177 * t178 - g(1) * (-pkin(2) * t321 + pkin(7) * t199 + t200 * t255 + t263) - g(2) * (pkin(2) * t322 + t200 * pkin(7) - t199 * t255 + t292) + t326 * t189 + (t179 * t284 - t180 * t285 + t257) * t202) * MDP(20) + (t232 * t249 - t234 * t242) * MDP(15) + (-t265 * pkin(1) - g(1) * t263 - g(2) * t292 + (t275 + t281) * qJ(2)) * MDP(6) + (t232 * t323 + t243 * t234) * MDP(17) + (t243 * t232 - t234 * t323) * MDP(19) - t241 + (t228 * t293 + t240 + t304) * MDP(9) + (-t228 * t268 + t245 + t303) * MDP(8) + (-t202 * t270 - t289 * t304 - t239) * MDP(18) + t253 * MDP(3) + qJDD(1) * MDP(1) + (-t253 + t275 + 0.2e1 * t281) * MDP(5); -qJDD(1) * MDP(4) - t238 * MDP(5) + (-qJ(2) * t238 + t265 - t291) * MDP(6) - t291 * MDP(20) + (MDP(15) + MDP(17)) * (t232 * t324 + t248 * t234) + (-MDP(16) + MDP(19)) * (t248 * t232 - t234 * t324) + ((t244 - t331) * MDP(20) - MDP(8) * t334 + t335 * t228) * t233 + (-t168 * MDP(20) - t228 * MDP(8) + (-t326 * MDP(20) - t276 * t335) * t276) * t235; (-t245 - t301) * MDP(8) + (-t240 - t302) * MDP(9) + ((qJD(4) * t271 - t312) * t232 + t247 * t234 + t295) * MDP(15) + ((-t312 + (t191 + t271) * qJD(4)) * t234 + (-t247 + t301) * t232) * MDP(16) + ((qJD(4) * t269 - t312) * t232 + (t246 + t300) * t234 + t295) * MDP(17) + (-pkin(7) * t270 + t289 * t302 + t239) * MDP(18) + ((t312 + (-t191 - t269) * qJD(4)) * t234 + (t276 * t294 + t246) * t232) * MDP(19) + (t239 * pkin(7) + t294 * t177 - t326 * t191 + (-t168 + t262) * t255) * MDP(20) + t241; -MDP(10) * t274 + t290 * MDP(11) * t334 - MDP(12) * t299 - t234 * t228 * MDP(13) + qJDD(4) * MDP(14) + (t184 * t298 + t264) * MDP(15) + (-g(3) * t232 - t174 + (t184 * t276 - t261) * t234) * MDP(16) + (0.2e1 * t311 - qJDD(5) - (-t177 * t232 - t196 * t234) * t276 + t264) * MDP(17) + (t258 * t228 - ((t180 - t286) * t232 + (qJD(5) - t179 - t314) * t234) * t276) * MDP(18) + (0.2e1 * t279 + 0.2e1 * qJD(4) * qJD(5) + t174 + (t196 * t276 + g(3)) * t232 + (t261 - t331) * t234) * MDP(19) + (-t170 * pkin(4) - g(3) * t259 + t169 * qJ(5) + t177 * t196 - t179 * t305 + t267 * t180 - t258 * t261) * MDP(20); (-qJDD(4) - t274) * MDP(17) - MDP(18) * t299 + (-t230 * t334 - t237) * MDP(19) + (-qJD(4) * t180 - t177 * t298 + t170 - t325) * MDP(20);];
tau = t1;
