% Calculate vector of inverse dynamics joint torques for
% S4RRPR6
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:46
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:46:27
% EndTime: 2021-01-15 10:46:34
% DurationCPUTime: 1.94s
% Computational Cost: add. (1165->249), mult. (2752->336), div. (0->0), fcn. (1963->12), ass. (0->119)
t268 = sin(pkin(7));
t269 = cos(pkin(7));
t275 = cos(qJ(2));
t304 = qJD(1) * t275;
t296 = t269 * t304;
t272 = sin(qJ(2));
t305 = qJD(1) * t272;
t224 = t268 * t305 - t296;
t274 = cos(qJ(4));
t215 = t274 * t224;
t236 = t268 * t275 + t269 * t272;
t227 = t236 * qJD(1);
t271 = sin(qJ(4));
t310 = t227 * t271;
t189 = -t215 - t310;
t264 = qJD(2) + qJD(4);
t311 = t189 * t264;
t286 = t224 * t271 - t274 * t227;
t312 = t286 * t264;
t273 = sin(qJ(1));
t276 = cos(qJ(1));
t291 = g(1) * t276 + g(2) * t273;
t270 = -qJ(3) - pkin(5);
t293 = qJD(2) * t270;
t222 = -qJD(3) * t272 + t275 * t293;
t247 = t270 * t272;
t195 = qJDD(2) * pkin(2) + t222 * qJD(1) + qJDD(1) * t247;
t221 = qJD(3) * t275 + t272 * t293;
t248 = t270 * t275;
t202 = t221 * qJD(1) - qJDD(1) * t248;
t168 = t269 * t195 - t202 * t268;
t302 = qJD(1) * qJD(2);
t295 = t272 * t302;
t249 = t268 * t295;
t294 = t275 * t302;
t199 = t236 * qJDD(1) + t269 * t294 - t249;
t166 = qJDD(2) * pkin(3) - pkin(6) * t199 + t168;
t169 = t268 * t195 + t269 * t202;
t226 = t236 * qJD(2);
t300 = qJDD(1) * t275;
t250 = t269 * t300;
t301 = qJDD(1) * t272;
t198 = qJD(1) * t226 + t268 * t301 - t250;
t167 = -pkin(6) * t198 + t169;
t240 = qJD(1) * t247;
t314 = qJD(2) * pkin(2);
t234 = t240 + t314;
t241 = qJD(1) * t248;
t308 = t269 * t241;
t197 = t268 * t234 - t308;
t319 = pkin(6) * t224;
t177 = t197 - t319;
t259 = pkin(2) * t275 + pkin(1);
t242 = -t259 * qJD(1) + qJD(3);
t205 = pkin(3) * t224 + t242;
t265 = qJ(2) + pkin(7);
t262 = qJ(4) + t265;
t255 = sin(t262);
t256 = cos(t262);
t303 = qJD(4) * t271;
t324 = g(3) * t255 - t271 * t166 - t274 * t167 + t177 * t303 - t205 * t189 + t291 * t256;
t323 = -g(3) * t256 + t274 * t166 - t271 * t167 + t205 * t286 + t291 * t255;
t263 = qJDD(2) + qJDD(4);
t322 = t263 * MDP(19) + t189 * MDP(15) * t286 + (-t189 ^ 2 + t286 ^ 2) * MDP(16);
t292 = t274 * t198 + t199 * t271;
t165 = -t286 * qJD(4) + t292;
t321 = pkin(2) * t268;
t320 = pkin(2) * t272;
t318 = pkin(6) * t227;
t315 = g(3) * t275;
t230 = t268 * t241;
t196 = t269 * t234 + t230;
t175 = qJD(2) * pkin(3) + t196 - t318;
t313 = t175 * t274;
t309 = t268 * t272;
t181 = t269 * t221 + t268 * t222;
t204 = t269 * t240 + t230;
t207 = t268 * t247 - t269 * t248;
t266 = t272 ^ 2;
t307 = -t275 ^ 2 + t266;
t299 = pkin(2) * t295 + qJDD(3);
t298 = t272 * t314;
t297 = -qJD(4) * t215 - t271 * t198 + t274 * t199;
t180 = -t221 * t268 + t269 * t222;
t203 = -t240 * t268 + t308;
t206 = t269 * t247 + t248 * t268;
t290 = g(1) * t273 - g(2) * t276;
t289 = -t175 * t271 - t177 * t274;
t182 = -pkin(6) * t236 + t206;
t235 = -t269 * t275 + t309;
t183 = -pkin(6) * t235 + t207;
t288 = t182 * t274 - t183 * t271;
t287 = t182 * t271 + t183 * t274;
t200 = t274 * t235 + t236 * t271;
t201 = -t235 * t271 + t236 * t274;
t257 = pkin(2) * t269 + pkin(3);
t285 = t257 * t271 + t274 * t321;
t284 = t257 * t274 - t271 * t321;
t283 = -0.2e1 * pkin(1) * t302 - pkin(5) * qJDD(2);
t164 = -t227 * t303 + t297;
t220 = -t259 * qJDD(1) + t299;
t277 = qJD(2) ^ 2;
t280 = 0.2e1 * qJDD(1) * pkin(1) - pkin(5) * t277 + t290;
t278 = qJD(1) ^ 2;
t279 = pkin(1) * t278 - pkin(5) * qJDD(1) + t291;
t261 = cos(t265);
t260 = sin(t265);
t229 = t235 * qJD(2);
t211 = pkin(3) * t235 - t259;
t209 = pkin(3) * t226 + t298;
t208 = pkin(2) * t305 + pkin(3) * t227;
t179 = t204 - t318;
t178 = t203 + t319;
t176 = pkin(3) * t198 + t220;
t173 = -pkin(6) * t226 + t181;
t172 = pkin(6) * t229 + t180;
t171 = t201 * qJD(4) + t274 * t226 - t229 * t271;
t170 = -t200 * qJD(4) - t226 * t271 - t229 * t274;
t1 = [qJDD(1) * MDP(1) + t290 * MDP(2) + t291 * MDP(3) + (qJDD(1) * t266 + 0.2e1 * t272 * t294) * MDP(4) + 0.2e1 * (t272 * t300 - t307 * t302) * MDP(5) + (qJDD(2) * t272 + t275 * t277) * MDP(6) + (qJDD(2) * t275 - t272 * t277) * MDP(7) + (t283 * t272 + t280 * t275) * MDP(9) + (-t280 * t272 + t283 * t275) * MDP(10) + (qJDD(2) * t206 - t198 * t259 + t220 * t235 + t226 * t242 + t290 * t261 + (t224 * t320 + t180) * qJD(2)) * MDP(11) + (-qJDD(2) * t207 - t199 * t259 + t220 * t236 - t229 * t242 - t290 * t260 + (t227 * t320 - t181) * qJD(2)) * MDP(12) + (-t168 * t236 - t169 * t235 - t180 * t227 - t181 * t224 + t196 * t229 - t197 * t226 - t198 * t207 - t199 * t206 - t291) * MDP(13) + (t169 * t207 + t197 * t181 + t168 * t206 + t196 * t180 - t220 * t259 + t242 * t298 - g(1) * (-t259 * t273 - t270 * t276) - g(2) * (t259 * t276 - t270 * t273)) * MDP(14) + (t164 * t201 - t170 * t286) * MDP(15) + (-t164 * t200 - t165 * t201 + t170 * t189 + t171 * t286) * MDP(16) + (t170 * t264 + t201 * t263) * MDP(17) + (-t171 * t264 - t200 * t263) * MDP(18) + (-t209 * t189 + t211 * t165 + t176 * t200 + t205 * t171 + (-t287 * qJD(4) + t172 * t274 - t173 * t271) * t264 + t288 * t263 + t290 * t256) * MDP(20) + (-t209 * t286 + t211 * t164 + t176 * t201 + t205 * t170 - (t288 * qJD(4) + t172 * t271 + t173 * t274) * t264 - t287 * t263 - t290 * t255) * MDP(21); MDP(6) * t301 + MDP(7) * t300 + qJDD(2) * MDP(8) + (t279 * t272 - t315) * MDP(9) + (g(3) * t272 + t279 * t275) * MDP(10) + (-g(3) * t261 - qJD(2) * t203 - t227 * t242 + t291 * t260 + (qJDD(2) * t269 - t224 * t305) * pkin(2) + t168) * MDP(11) + (g(3) * t260 + qJD(2) * t204 + t224 * t242 + t291 * t261 + (-qJDD(2) * t268 - t227 * t305) * pkin(2) - t169) * MDP(12) + ((t197 + t203) * t227 + (-t196 + t204) * t224 + (-t198 * t268 - t199 * t269) * pkin(2)) * MDP(13) + (-t196 * t203 - t197 * t204 + (-t315 + t168 * t269 + t169 * t268 + (-qJD(1) * t242 + t291) * t272) * pkin(2)) * MDP(14) + (t164 - t311) * MDP(17) + (-t165 - t312) * MDP(18) + (t284 * t263 + t208 * t189 - (t178 * t274 - t179 * t271) * t264 + (-t285 * t264 + t289) * qJD(4) + t323) * MDP(20) + (-t285 * t263 + t208 * t286 + (t178 * t271 + t179 * t274) * t264 + (-t284 * t264 - t313) * qJD(4) + t324) * MDP(21) + (-t272 * t275 * MDP(4) + t307 * MDP(5)) * t278 + t322; -t250 * MDP(11) - t249 * MDP(12) + (-t224 ^ 2 - t227 ^ 2) * MDP(13) + (t196 * t227 + t197 * t224 - t290 + t299) * MDP(14) + (t165 - t312) * MDP(20) + (t164 + t311) * MDP(21) + (MDP(11) * t309 + t236 * MDP(12) - t259 * MDP(14)) * qJDD(1) + ((t268 * t304 + t269 * t305 + t227) * MDP(11) + (-t224 + t296) * MDP(12)) * qJD(2); (t297 - t311) * MDP(17) + (-t292 - t312) * MDP(18) + (-t289 * t264 + t323) * MDP(20) + ((-t177 * t271 + t313) * t264 + t324) * MDP(21) + (-MDP(17) * t310 + t286 * MDP(18) + t289 * MDP(20) - MDP(21) * t313) * qJD(4) + t322;];
tau = t1;
