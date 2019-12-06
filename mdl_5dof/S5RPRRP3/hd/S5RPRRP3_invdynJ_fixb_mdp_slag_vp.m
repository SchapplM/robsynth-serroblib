% Calculate vector of inverse dynamics joint torques for
% S5RPRRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:04:02
% EndTime: 2019-12-05 18:04:06
% DurationCPUTime: 1.94s
% Computational Cost: add. (1410->224), mult. (2957->294), div. (0->0), fcn. (1922->12), ass. (0->121)
t264 = qJ(3) + qJ(4);
t255 = sin(t264);
t256 = cos(t264);
t261 = qJ(1) + pkin(8);
t251 = sin(t261);
t252 = cos(t261);
t339 = -g(2) * t251 + g(3) * t252;
t342 = -g(1) * t256 + t255 * t339;
t268 = sin(qJ(3));
t270 = cos(qJ(3));
t265 = sin(pkin(8));
t244 = pkin(1) * t265 + pkin(6);
t330 = pkin(7) + t244;
t302 = t330 * qJD(1);
t202 = qJD(2) * t268 + t302 * t270;
t267 = sin(qJ(4));
t197 = t267 * t202;
t201 = t270 * qJD(2) - t302 * t268;
t329 = qJD(3) * pkin(3);
t200 = t201 + t329;
t335 = cos(qJ(4));
t300 = t335 * t200 - t197;
t223 = t267 * t270 + t335 * t268;
t217 = t223 * qJD(1);
t326 = t217 * qJ(5);
t341 = t326 - t300;
t260 = qJD(3) + qJD(4);
t195 = t260 * t223;
t303 = qJDD(1) * t335;
t311 = qJDD(1) * t268;
t291 = t267 * t311 - t270 * t303;
t187 = t195 * qJD(1) + t291;
t312 = qJD(1) * qJD(3);
t305 = t268 * t312;
t266 = cos(pkin(8));
t245 = -pkin(1) * t266 - pkin(2);
t257 = t270 * pkin(3);
t338 = t245 - t257;
t203 = pkin(3) * t305 + qJDD(1) * t338;
t340 = pkin(4) * t187 + qJDD(5) + t203;
t220 = t330 * t268;
t221 = t330 * t270;
t318 = -t267 * t220 + t335 * t221;
t336 = t217 ^ 2;
t308 = t335 * t270;
t296 = qJD(1) * t308;
t315 = qJD(1) * t268;
t307 = t267 * t315;
t215 = -t296 + t307;
t219 = t338 * qJD(1);
t192 = pkin(4) * t215 + qJD(5) + t219;
t328 = t192 * t217;
t327 = t215 * qJ(5);
t325 = t251 * t256;
t324 = t252 * t256;
t323 = t267 * t268;
t322 = qJDD(2) - g(1);
t174 = pkin(4) * t260 - t341;
t321 = t174 + t341;
t292 = t260 * t323;
t306 = qJD(4) * t335;
t194 = -qJD(3) * t308 - t270 * t306 + t292;
t320 = -t223 * t187 + t194 * t215;
t319 = t335 * t201 - t197;
t317 = pkin(4) * t256 + t257;
t262 = t268 ^ 2;
t316 = -t270 ^ 2 + t262;
t234 = qJD(1) * t245;
t313 = qJD(4) * t267;
t310 = qJDD(1) * t270;
t309 = t268 * t329;
t199 = t335 * t202;
t304 = qJD(3) * t330;
t231 = t244 * qJDD(1);
t301 = pkin(7) * qJDD(1) + t231;
t299 = -t201 * t267 - t199;
t298 = -t335 * t220 - t221 * t267;
t297 = -t260 * t296 - t267 * t310 - t268 * t303;
t295 = g(2) * t252 + g(3) * t251;
t269 = sin(qJ(1));
t271 = cos(qJ(1));
t293 = g(2) * t271 + g(3) * t269;
t186 = qJD(1) * t292 + t297;
t222 = -t308 + t323;
t289 = -t186 * t222 + t195 * t217;
t258 = qJDD(3) + qJDD(4);
t288 = t194 * t260 - t223 * t258;
t286 = -t267 * t200 - t199;
t212 = t268 * t304;
t213 = t270 * t304;
t285 = -t335 * t212 - t267 * t213 - t220 * t306 - t221 * t313;
t214 = t215 ^ 2;
t283 = t217 * t215 * MDP(12) + (-t297 + (t215 - t307) * t260) * MDP(14) - t291 * MDP(15) + (-t214 + t336) * MDP(13) + t258 * MDP(16);
t282 = -qJD(1) * t234 - t231 + t339;
t281 = 0.2e1 * t234 * qJD(3) - qJDD(3) * t244;
t272 = qJD(3) ^ 2;
t280 = -0.2e1 * qJDD(1) * t245 - t244 * t272 + t295;
t253 = t270 * qJDD(2);
t183 = qJDD(3) * pkin(3) - t202 * qJD(3) - t301 * t268 + t253;
t184 = t201 * qJD(3) + t268 * qJDD(2) + t301 * t270;
t279 = t286 * qJD(4) + t335 * t183 - t267 * t184;
t278 = -t318 * qJD(4) + t267 * t212 - t335 * t213;
t277 = t267 * t183 + t335 * t184 + t200 * t306 - t202 * t313;
t276 = g(1) * t255 - g(2) * t325 + g(3) * t324 + t219 * t215 - t277;
t275 = -t219 * t217 + t279 + t342;
t259 = -qJ(5) - pkin(7) - pkin(6);
t250 = t335 * pkin(3) + pkin(4);
t230 = qJDD(3) * t270 - t268 * t272;
t229 = qJDD(3) * t268 + t270 * t272;
t224 = pkin(2) + t317;
t189 = -qJ(5) * t222 + t318;
t188 = -qJ(5) * t223 + t298;
t185 = -t195 * t260 - t222 * t258;
t178 = -t326 + t319;
t177 = t299 + t327;
t176 = -t286 - t327;
t173 = t194 * qJ(5) - t223 * qJD(5) + t278;
t172 = -qJ(5) * t195 - qJD(5) * t222 + t285;
t171 = -t187 * qJ(5) - t215 * qJD(5) + t277;
t170 = t258 * pkin(4) + t186 * qJ(5) - t217 * qJD(5) + t279;
t1 = [qJDD(1) * MDP(1) + t293 * MDP(2) + (-g(2) * t269 + g(3) * t271) * MDP(3) + (t293 + (t265 ^ 2 + t266 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t262 + 0.2e1 * t270 * t305) * MDP(5) + 0.2e1 * (t268 * t310 - t316 * t312) * MDP(6) + t229 * MDP(7) + t230 * MDP(8) + (t281 * t268 + t280 * t270) * MDP(10) + (-t280 * t268 + t281 * t270) * MDP(11) + (-t186 * t223 - t194 * t217) * MDP(12) + (-t289 + t320) * MDP(13) - t288 * MDP(14) + t185 * MDP(15) + (g(2) * t324 + g(3) * t325 + t187 * t338 + t219 * t195 + t203 * t222 + t215 * t309 + t298 * t258 + t278 * t260) * MDP(17) + (-t186 * t338 - t219 * t194 + t203 * t223 + t217 * t309 - t295 * t255 - t318 * t258 - t285 * t260) * MDP(18) + (-t170 * t223 - t171 * t222 - t172 * t215 - t173 * t217 + t174 * t194 - t176 * t195 + t186 * t188 - t187 * t189 - t339) * MDP(19) + (t171 * t189 + t176 * t172 + t170 * t188 + t174 * t173 + t340 * (pkin(4) * t222 + t338) + t192 * (pkin(4) * t195 + t309) - g(2) * (-pkin(1) * t271 - t224 * t252 + t251 * t259) - g(3) * (-pkin(1) * t269 - t224 * t251 - t252 * t259)) * MDP(20); t322 * MDP(4) + t230 * MDP(10) - t229 * MDP(11) + t185 * MDP(17) + t288 * MDP(18) + (t289 + t320) * MDP(19) + (-t170 * t222 + t171 * t223 - t174 * t195 - t176 * t194 - g(1)) * MDP(20); MDP(7) * t311 + MDP(8) * t310 + qJDD(3) * MDP(9) + (-g(1) * t270 + t282 * t268 + t253) * MDP(10) + (-t322 * t268 + t282 * t270) * MDP(11) + (-t299 * t260 + (-t215 * t315 + t335 * t258 - t260 * t313) * pkin(3) + t275) * MDP(17) + (t319 * t260 + (-t217 * t315 - t267 * t258 - t260 * t306) * pkin(3) + t276) * MDP(18) + (t250 * t186 + (t176 + t177) * t217 + (-t174 + t178) * t215 + (-t187 * t267 + (-t335 * t215 + t217 * t267) * qJD(4)) * pkin(3)) * MDP(19) + (t170 * t250 - t176 * t178 - t174 * t177 - pkin(4) * t328 - g(1) * t317 - t339 * (-pkin(3) * t268 - pkin(4) * t255) + (-t192 * t315 + t171 * t267 + (-t174 * t267 + t335 * t176) * qJD(4)) * pkin(3)) * MDP(20) + t283 + (-t268 * t270 * MDP(5) + t316 * MDP(6)) * qJD(1) ^ 2; (-t286 * t260 + t275) * MDP(17) + (t300 * t260 + t276) * MDP(18) + (pkin(4) * t186 - t321 * t215) * MDP(19) + (t321 * t176 + (t170 - t328 + t342) * pkin(4)) * MDP(20) + t283; (-t214 - t336) * MDP(19) + (t174 * t217 + t176 * t215 - t295 + t340) * MDP(20);];
tau = t1;
