% Calculate vector of inverse dynamics joint torques for
% S5RPPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRR12_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR12_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR12_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR12_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPPRR12_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:27
% EndTime: 2019-12-31 18:07:31
% DurationCPUTime: 2.85s
% Computational Cost: add. (1384->281), mult. (2788->361), div. (0->0), fcn. (1992->10), ass. (0->130)
t251 = sin(pkin(8));
t255 = sin(qJ(4));
t298 = qJD(1) * t255;
t285 = t251 * t298;
t252 = cos(pkin(8));
t258 = cos(qJ(4));
t297 = qJD(1) * t258;
t287 = t252 * t297;
t209 = -t285 + t287;
t254 = sin(qJ(5));
t257 = cos(qJ(5));
t194 = -t257 * qJD(4) + t209 * t254;
t215 = t251 * t258 + t252 * t255;
t207 = t215 * qJD(1);
t332 = qJD(5) + t207;
t337 = t194 * t332;
t196 = qJD(4) * t254 + t209 * t257;
t336 = t196 * t332;
t303 = t258 * t252;
t223 = qJDD(1) * t303;
t289 = qJDD(1) * t251;
t272 = -t255 * t289 + t223;
t190 = -qJD(4) * t207 + t272;
t335 = qJD(4) * qJD(5) + t190;
t253 = -pkin(1) - qJ(3);
t222 = qJD(1) * t253 + qJD(2);
t283 = -pkin(6) * qJD(1) + t222;
t199 = t283 * t251;
t200 = t283 * t252;
t182 = t199 * t258 + t200 * t255;
t323 = -qJD(1) * qJD(3) + qJDD(1) * t253;
t216 = qJDD(2) + t323;
t279 = -pkin(6) * qJDD(1) + t216;
t197 = t279 * t251;
t198 = t279 * t252;
t270 = t197 * t255 - t198 * t258;
t169 = -qJDD(4) * pkin(4) + qJD(4) * t182 + t270;
t247 = pkin(8) + qJ(4);
t234 = sin(t247);
t235 = cos(t247);
t256 = sin(qJ(1));
t259 = cos(qJ(1));
t328 = g(1) * t256 - g(2) * t259;
t263 = g(3) * t234 - t235 * t328;
t334 = t263 - t332 * (pkin(4) * t209 + t332 * pkin(7)) - t169;
t278 = t257 * t332;
t295 = qJD(4) * t258;
t286 = t252 * t295;
t221 = qJD(4) * t285;
t324 = -t215 * qJDD(1) + t221;
t191 = qJD(1) * t286 - t324;
t187 = qJDD(5) + t191;
t308 = t254 * t187;
t333 = -t278 * t332 - t308;
t299 = t251 ^ 2 + t252 ^ 2;
t330 = t222 * t299;
t248 = qJDD(1) * qJ(2);
t249 = qJD(1) * qJD(2);
t327 = t248 + t249;
t220 = qJDD(3) + t327;
t274 = g(1) * t259 + g(2) * t256;
t329 = t220 - t274;
t326 = t251 * MDP(7) + t252 * MDP(8);
t318 = -pkin(6) + t253;
t217 = t318 * t251;
t218 = t318 * t252;
t192 = t217 * t255 - t218 * t258;
t176 = -qJD(3) * t215 - qJD(4) * t192;
t269 = t199 * t255 - t200 * t258;
t178 = -qJD(4) * pkin(4) + t269;
t214 = t251 * t255 - t303;
t239 = t251 * pkin(3);
t228 = qJ(2) + t239;
t188 = pkin(4) * t215 + pkin(7) * t214 + t228;
t193 = t217 * t258 + t218 * t255;
t296 = qJD(4) * t255;
t210 = -t251 * t295 - t252 * t296;
t271 = t197 * t258 + t198 * t255;
t168 = qJDD(4) * pkin(7) - qJD(4) * t269 + t271;
t233 = qJD(1) * qJ(2) + qJD(3);
t219 = qJD(1) * t239 + t233;
t180 = pkin(4) * t207 - pkin(7) * t209 + t219;
t277 = qJD(5) * t180 + t168;
t322 = -t169 * t214 + t178 * t210 - t193 * t187 - (qJD(5) * t188 + t176) * t332 - t215 * t277;
t321 = 0.2e1 * t249;
t319 = g(3) * t235;
t317 = pkin(1) * qJDD(1);
t288 = t254 * qJDD(4) + t335 * t257;
t294 = qJD(5) * t254;
t174 = -t209 * t294 + t288;
t316 = t174 * t214;
t315 = t174 * t254;
t314 = t178 * t214;
t313 = t188 * t187;
t312 = t194 * t209;
t311 = t196 * t209;
t307 = t254 * t256;
t306 = t254 * t259;
t305 = t256 * t257;
t183 = t257 * t187;
t304 = t257 * t259;
t302 = t210 * qJD(4) - t214 * qJDD(4);
t301 = t259 * pkin(1) + t256 * qJ(2);
t293 = qJD(5) * t257;
t284 = qJDD(2) - t328;
t281 = t299 * MDP(9);
t280 = t299 * t216;
t213 = pkin(3) * t289 + t220;
t173 = pkin(4) * t191 - pkin(7) * t190 + t213;
t179 = qJD(4) * pkin(7) + t182;
t276 = qJD(5) * t179 - t173;
t211 = -t251 * t296 + t286;
t268 = -qJD(4) * t211 - qJDD(4) * t215;
t267 = t183 + (-t207 * t254 - t294) * t332;
t266 = t210 * t257 + t214 * t294;
t262 = -pkin(7) * t187 + (t178 - t269) * t332;
t260 = qJD(1) ^ 2;
t241 = t259 * qJ(2);
t237 = t257 * qJDD(4);
t206 = t234 * t304 - t307;
t205 = t234 * t306 + t305;
t204 = t234 * t305 + t306;
t203 = -t234 * t307 + t304;
t185 = pkin(4) * t211 - pkin(7) * t210 + qJD(2);
t177 = -qJD(3) * t214 + qJD(4) * t193;
t175 = t196 * qJD(5) + t190 * t254 - t237;
t172 = t257 * t173;
t171 = t179 * t257 + t180 * t254;
t170 = -t179 * t254 + t180 * t257;
t1 = [qJDD(1) * MDP(1) + t328 * MDP(2) + t274 * MDP(3) + (t284 - 0.2e1 * t317) * MDP(4) + (0.2e1 * t248 + t321 - t274) * MDP(5) + (-(qJDD(2) - t317) * pkin(1) - g(1) * (-pkin(1) * t256 + t241) - g(2) * t301 + (t248 + t321) * qJ(2)) * MDP(6) + (t328 + t299 * (-t216 - t323)) * MDP(9) + (t220 * qJ(2) + t233 * qJD(2) - g(1) * (t253 * t256 + t241) - g(2) * (qJ(3) * t259 + t301) + t253 * t280 - qJD(3) * t330) * MDP(10) + (-t190 * t214 + t209 * t210) * MDP(11) + (-t190 * t215 + t191 * t214 - t207 * t210 - t209 * t211) * MDP(12) + t302 * MDP(13) + t268 * MDP(14) + (qJD(2) * t207 - qJD(4) * t177 - qJDD(4) * t192 + t191 * t228 + t211 * t219 + t213 * t215 - t234 * t274) * MDP(16) + (qJD(2) * t209 - qJD(4) * t176 - qJDD(4) * t193 + t190 * t228 + t210 * t219 - t213 * t214 - t235 * t274) * MDP(17) + (t196 * t266 - t257 * t316) * MDP(18) + ((-t194 * t257 - t196 * t254) * t210 + (t315 + t175 * t257 + (-t194 * t254 + t196 * t257) * qJD(5)) * t214) * MDP(19) + (t174 * t215 - t183 * t214 + t196 * t211 + t266 * t332) * MDP(20) + (t214 * t308 - t175 * t215 - t194 * t211 + (-t210 * t254 + t214 * t293) * t332) * MDP(21) + (t187 * t215 + t211 * t332) * MDP(22) + (-g(1) * t206 - g(2) * t204 + t170 * t211 + t172 * t215 + t192 * t175 + t177 * t194 + (t185 * t332 + t313 + (-t179 * t215 - t193 * t332 - t314) * qJD(5)) * t257 + t322 * t254) * MDP(23) + (g(1) * t205 - g(2) * t203 - t171 * t211 + t192 * t174 + t177 * t196 + (-(-qJD(5) * t193 + t185) * t332 - t313 + t276 * t215 + qJD(5) * t314) * t254 + t322 * t257) * MDP(24) + t326 * (t327 + t329); t284 * MDP(6) + (-qJD(1) * t233 + t280 - t328) * MDP(10) + (-qJD(1) * t207 + t302) * MDP(16) + (-qJD(1) * t209 + t268) * MDP(17) + (t175 * t214 - t194 * t210 - t215 * t308) * MDP(23) + (-t215 * t183 - t196 * t210 + t316) * MDP(24) + (-MDP(6) * qJ(2) - MDP(5) - t326) * t260 + (-pkin(1) * MDP(6) + MDP(4) - t281) * qJDD(1) + ((-qJD(1) * t257 - t211 * t254 - t215 * t293) * MDP(23) + (qJD(1) * t254 - t211 * t257 + t215 * t294) * MDP(24)) * t332; (qJD(1) * t330 + t329) * MDP(10) - t221 * MDP(16) + t223 * MDP(17) + (t267 - t312) * MDP(23) + (-t311 + t333) * MDP(24) - t260 * t281 + ((MDP(16) * t255 + MDP(8)) * t252 + (MDP(16) * t258 - MDP(17) * t255 + MDP(7)) * t251) * qJDD(1) + ((t209 + t287) * MDP(16) + (-t251 * t297 - t252 * t298 - t207) * MDP(17)) * qJD(4); t209 * t207 * MDP(11) + (-t207 ^ 2 + t209 ^ 2) * MDP(12) + t272 * MDP(13) + ((t209 - t287) * qJD(4) + t324) * MDP(14) + qJDD(4) * MDP(15) + (-t209 * t219 + t263 - t270) * MDP(16) + (t207 * t219 + t234 * t328 - t271 + t319) * MDP(17) + (t196 * t278 + t315) * MDP(18) + ((t174 - t337) * t257 + (-t175 - t336) * t254) * MDP(19) + (-t311 - t333) * MDP(20) + (t267 + t312) * MDP(21) - t332 * t209 * MDP(22) + (-pkin(4) * t175 - t170 * t209 - t182 * t194 + t262 * t254 + t334 * t257) * MDP(23) + (-pkin(4) * t174 + t171 * t209 - t182 * t196 - t334 * t254 + t262 * t257) * MDP(24); t196 * t194 * MDP(18) + (-t194 ^ 2 + t196 ^ 2) * MDP(19) + (t288 + t337) * MDP(20) + (t237 + t336) * MDP(21) + t187 * MDP(22) + (-g(1) * t203 - g(2) * t205 + t171 * t332 - t178 * t196 + t172) * MDP(23) + (g(1) * t204 - g(2) * t206 + t170 * t332 + t178 * t194) * MDP(24) + ((-t168 + t319) * MDP(24) + (-MDP(21) * t209 - MDP(23) * t179 - MDP(24) * t180) * qJD(5)) * t257 + (-qJD(5) * t209 * MDP(20) - t335 * MDP(21) + (-t277 + t319) * MDP(23) + t276 * MDP(24)) * t254;];
tau = t1;
