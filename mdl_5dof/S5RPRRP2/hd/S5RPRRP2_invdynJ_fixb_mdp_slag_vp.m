% Calculate vector of inverse dynamics joint torques for
% S5RPRRP2
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
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:36
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPRRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:36:31
% EndTime: 2021-01-15 12:36:36
% DurationCPUTime: 1.38s
% Computational Cost: add. (1348->242), mult. (2232->281), div. (0->0), fcn. (1229->12), ass. (0->130)
t331 = qJDD(2) - g(1);
t248 = cos(pkin(8));
t232 = pkin(1) * t248 + pkin(2);
t247 = sin(pkin(8));
t326 = pkin(1) * t247;
t294 = qJD(1) * t326;
t330 = -qJD(3) * t294 + t232 * qJDD(1);
t218 = t232 * qJD(1);
t329 = qJD(3) * t218 + qJDD(1) * t326;
t253 = cos(qJ(4));
t251 = sin(qJ(3));
t254 = cos(qJ(3));
t195 = t218 * t251 + t254 * t294;
t243 = qJD(1) + qJD(3);
t187 = pkin(7) * t243 + t195;
t288 = qJ(5) * t243 + t187;
t273 = t288 * t253;
t244 = qJ(1) + pkin(8);
t236 = qJ(3) + t244;
t231 = cos(t236);
t226 = g(3) * t231;
t230 = sin(t236);
t227 = g(2) * t230;
t303 = t226 - t227;
t285 = t232 * t254 - t251 * t326;
t302 = t251 * t232 + t254 * t326;
t250 = sin(qJ(4));
t297 = qJD(4) * t250;
t290 = t243 * t297;
t323 = pkin(4) * t253;
t233 = pkin(3) + t323;
t242 = qJDD(1) + qJDD(3);
t314 = t233 * t242;
t328 = -pkin(4) * t290 + t314;
t327 = -t251 * t330 - t254 * t329;
t325 = pkin(3) * t242;
t245 = t250 ^ 2;
t324 = pkin(4) * t245;
t322 = g(1) * t253;
t321 = g(2) * t231;
t249 = -qJ(5) - pkin(7);
t320 = qJD(4) * pkin(4);
t319 = qJDD(4) * pkin(4);
t318 = t195 * t243;
t198 = t302 * qJD(3);
t317 = t198 * t243;
t316 = t230 * t250;
t313 = t242 * t250;
t312 = t242 * t253;
t311 = t243 * t250;
t310 = t243 * t253;
t309 = t250 * t253;
t200 = pkin(7) + t302;
t308 = -qJ(5) - t200;
t238 = t253 * qJD(2);
t176 = -t250 * t288 + t238;
t175 = t176 + t320;
t307 = t175 - t176;
t194 = t218 * t254 - t251 * t294;
t306 = t194 * t297 + t195 * t310;
t305 = t230 * t233 + t231 * t249;
t304 = g(3) * t316 + t250 * t321;
t246 = t253 ^ 2;
t301 = -t245 - t246;
t300 = t245 - t246;
t298 = qJD(4) * t243;
t296 = qJD(4) * t253;
t295 = qJD(4) * qJD(2);
t293 = pkin(4) * t297;
t183 = t187 * t297;
t173 = pkin(7) * t242 - t327;
t264 = -qJ(5) * t242 - t173 - t295;
t260 = qJD(5) * t243 - t264;
t167 = -t183 + (-qJ(5) * t298 + qJDD(2)) * t250 + t260 * t253;
t291 = t167 * t253 + t303;
t289 = qJD(4) * t249;
t287 = -t230 * t249 + t231 * t233;
t284 = qJD(4) * t308;
t283 = 0.2e1 * t243 * t296;
t279 = -t251 * t329 + t254 * t330;
t169 = qJDD(5) - t279 - t328;
t180 = -t233 * t243 + qJD(5) - t194;
t281 = t169 * t250 + t180 * t296 + t304;
t174 = -t279 - t325;
t186 = -pkin(3) * t243 - t194;
t280 = t174 * t250 + t186 * t296 + t304;
t199 = -pkin(3) - t285;
t235 = t253 * qJDD(2);
t278 = g(2) * t316 + t235 - t322;
t256 = qJD(4) ^ 2;
t277 = pkin(7) * t256 - t325;
t276 = g(3) * t230 + t321;
t252 = sin(qJ(1));
t255 = cos(qJ(1));
t275 = -g(2) * t255 - g(3) * t252;
t209 = qJDD(4) * t250 + t253 * t256;
t210 = qJDD(4) * t253 - t250 * t256;
t274 = 0.2e1 * (t242 * t309 - t298 * t300) * MDP(9) + (t242 * t245 + t250 * t283) * MDP(8) + t209 * MDP(10) + t210 * MDP(11) + t242 * MDP(5);
t177 = qJD(2) * t250 + t273;
t272 = t175 * t250 - t177 * t253;
t192 = t198 + t293;
t193 = t199 - t323;
t271 = t192 * t243 + t193 * t242;
t270 = -t169 - t276;
t269 = -t174 - t276;
t268 = -t186 * t243 - t173 - t226;
t267 = t253 * t227 - t250 * t331 + t183;
t266 = -pkin(3) * t298 - pkin(7) * qJDD(4);
t263 = t199 * t242 + t200 * t256 + t317;
t262 = t276 - t279;
t197 = t285 * qJD(3);
t261 = -qJDD(4) * t200 + (t199 * t243 - t197) * qJD(4);
t259 = -t226 + (-qJD(5) - t180) * t243 + t264;
t258 = -t303 + t327;
t241 = t243 ^ 2;
t239 = t253 * qJ(5);
t237 = t253 * qJD(5);
t221 = pkin(7) * t253 + t239;
t220 = t249 * t250;
t202 = -qJD(5) * t250 + t253 * t289;
t201 = t250 * t289 + t237;
t191 = t200 * t253 + t239;
t190 = t308 * t250;
t189 = t194 * t296;
t181 = t186 * t297;
t178 = t180 * t297;
t171 = (-qJD(5) - t197) * t250 + t253 * t284;
t170 = t197 * t253 + t250 * t284 + t237;
t166 = -t273 * qJD(4) - t250 * t260 + t235 + t319;
t1 = [qJDD(1) * MDP(1) + t275 * MDP(2) + (g(2) * t252 - g(3) * t255) * MDP(3) + (t275 + (t247 ^ 2 + t248 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (t242 * t285 - t262 - t317) * MDP(6) + (-t197 * t243 - t242 * t302 + t258) * MDP(7) + (t181 + t261 * t250 + (-t263 + t269) * t253) * MDP(13) + (t250 * t263 + t253 * t261 + t280) * MDP(14) + (qJDD(4) * t190 + t178 + (t193 * t311 + t171) * qJD(4) + (t270 - t271) * t253) * MDP(15) + (-qJDD(4) * t191 + t271 * t250 + (t193 * t310 - t170) * qJD(4) + t281) * MDP(16) + ((t170 * t243 + t191 * t242 + (-t190 * t243 - t175) * qJD(4)) * t253 + (-t171 * t243 - t190 * t242 - t166 + (-t191 * t243 - t177) * qJD(4)) * t250 + t291) * MDP(17) + (t167 * t191 + t177 * t170 + t166 * t190 + t175 * t171 + t169 * t193 + t180 * t192 - g(2) * (pkin(2) * cos(t244) + t255 * pkin(1) + t287) - g(3) * (pkin(2) * sin(t244) + t252 * pkin(1) + t305)) * MDP(18) + t274; t331 * MDP(4) + (-qJD(4) * t272 + t166 * t253 + t167 * t250 - g(1)) * MDP(18) + (MDP(13) + MDP(15)) * t210 + (-MDP(14) - MDP(16)) * t209; (-t262 + t318) * MDP(6) + (t194 * t243 + t258) * MDP(7) + (t181 + t266 * t250 + (t269 - t277) * t253 + t306) * MDP(13) + (t189 + t266 * t253 + (t277 - t318) * t250 + t280) * MDP(14) + (qJDD(4) * t220 + t178 + (-t233 * t311 + t202) * qJD(4) + (t270 + t328) * t253 + t306) * MDP(15) + (-qJDD(4) * t221 + t189 + (-t314 - t318) * t250 + (-t201 + (-t233 * t253 + t324) * t243) * qJD(4) + t281) * MDP(16) + ((-qJD(4) * t175 + t221 * t242) * t253 + (-qJD(4) * t177 - t220 * t242 - t166) * t250 + (t201 * t253 - t202 * t250 + t301 * t194 + (-t220 * t253 - t221 * t250) * qJD(4)) * t243 + t291) * MDP(17) + (t167 * t221 + t166 * t220 - t169 * t233 - g(2) * t287 - g(3) * t305 + (-t195 + t293) * t180 + (-t194 * t253 + t201) * t177 + (t194 * t250 + t202) * t175) * MDP(18) + t274; -t241 * MDP(8) * t309 + t300 * t241 * MDP(9) + MDP(10) * t313 + MDP(11) * t312 + qJDD(4) * MDP(12) + (t250 * t268 + t278) * MDP(13) + ((-t187 * t250 + t238) * qJD(4) + (t268 - t295) * t253 + t267) * MDP(14) + (0.2e1 * t319 + (t177 - t273) * qJD(4) + (t241 * t323 + t259) * t250 + t278) * MDP(15) + (-t241 * t324 + (qJ(5) * t311 + t176) * qJD(4) + t259 * t253 + t267) * MDP(16) + (-pkin(4) * t313 + (t307 - t320) * t310) * MDP(17) + (t307 * t177 + (-t322 + t166 + (-t180 * t243 - t303) * t250) * pkin(4)) * MDP(18); (0.2e1 * t290 - t312) * MDP(15) + (t283 + t313) * MDP(16) + (t243 * t272 - t270) * MDP(18) + t301 * MDP(17) * t241;];
tau = t1;
