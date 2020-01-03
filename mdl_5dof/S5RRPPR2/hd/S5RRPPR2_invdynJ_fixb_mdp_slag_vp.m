% Calculate vector of inverse dynamics joint torques for
% S5RRPPR2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RRPPR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:57:41
% EndTime: 2020-01-03 11:57:44
% DurationCPUTime: 1.51s
% Computational Cost: add. (1255->236), mult. (2000->327), div. (0->0), fcn. (1204->14), ass. (0->133)
t250 = qJDD(1) + qJDD(2);
t265 = cos(qJ(2));
t336 = pkin(1) * t265;
t240 = qJDD(1) * t336;
t262 = sin(qJ(2));
t331 = pkin(1) * qJD(2);
t303 = qJD(1) * t331;
t203 = pkin(2) * t250 - t262 * t303 + t240;
t258 = sin(pkin(8));
t260 = cos(pkin(8));
t306 = qJDD(1) * t262;
t339 = pkin(1) * t306 + t265 * t303;
t186 = t260 * t203 - t339 * t258;
t273 = qJDD(4) - t186;
t183 = -pkin(3) * t250 + t273;
t256 = qJ(1) + qJ(2);
t244 = pkin(8) + t256;
t233 = sin(t244);
t234 = cos(t244);
t288 = g(2) * t234 + g(3) * t233;
t275 = -t183 - t288;
t187 = t258 * t203 + t339 * t260;
t253 = qJD(1) + qJD(2);
t181 = qJ(4) * t250 + qJD(4) * t253 + t187;
t257 = sin(pkin(9));
t259 = cos(pkin(9));
t180 = qJDD(3) * t257 + t181 * t259;
t179 = -t259 * qJDD(3) + t181 * t257;
t329 = t179 * t257;
t338 = t180 * t259 + t329;
t261 = sin(qJ(5));
t264 = cos(qJ(5));
t277 = MDP(17) * t261 + MDP(18) * t264;
t332 = pkin(1) * qJD(1);
t305 = t262 * t332;
t223 = t258 * t305;
t304 = t265 * t332;
t208 = t260 * t304 - t223;
t308 = qJD(4) - t208;
t337 = pkin(1) * t262;
t335 = pkin(2) * t260;
t330 = MDP(8) * t259;
t212 = pkin(2) * t253 + t304;
t224 = t260 * t305;
t196 = t258 * t212 + t224;
t192 = qJ(4) * t253 + t196;
t189 = -t259 * qJD(3) + t192 * t257;
t328 = t189 * t257;
t278 = -pkin(4) * t259 - pkin(7) * t257 - pkin(3);
t230 = t258 * t337;
t239 = pkin(2) + t336;
t293 = t239 * t260 - t230;
t193 = t278 - t293;
t327 = t193 * t264;
t276 = t260 * t265 * t331 - qJD(2) * t230;
t202 = qJD(4) + t276;
t326 = t202 * t253;
t325 = t202 * t261;
t324 = t250 * t259;
t323 = t250 * t261;
t322 = t253 * t259;
t321 = t257 * t261;
t320 = t259 * t261;
t319 = t259 * t264;
t318 = t260 * t262;
t317 = t261 * t264;
t316 = pkin(1) * t318 + t258 * t239;
t251 = t257 ^ 2;
t252 = t259 ^ 2;
t315 = t251 + t252;
t255 = t264 ^ 2;
t314 = t261 ^ 2 - t255;
t311 = qJD(5) * t261;
t310 = qJD(5) * t264;
t217 = -qJDD(5) + t324;
t309 = t217 * MDP(16);
t220 = -qJD(5) + t322;
t307 = -qJD(5) - t220;
t301 = t275 * t257;
t246 = cos(t256);
t238 = pkin(2) * t246;
t300 = t234 * pkin(3) + t233 * qJ(4) + t238;
t299 = t253 * t310;
t245 = sin(t256);
t298 = g(2) * t245 - g(3) * t246;
t297 = t315 * t250;
t176 = t278 * t250 + t273;
t296 = t264 * t176 - t261 * t180;
t195 = t212 * t260 - t223;
t295 = t217 - t324;
t294 = t217 + t324;
t292 = t308 * t264;
t291 = qJD(1) * (-qJD(2) + t253);
t290 = qJD(2) * (-qJD(1) - t253);
t287 = -g(2) * t246 - g(3) * t245;
t237 = pkin(2) * t245;
t286 = t233 * pkin(3) - qJ(4) * t234 + t237;
t285 = qJD(4) - t195;
t284 = t261 * t176 + t264 * t180;
t184 = t278 * t253 + t285;
t190 = qJD(3) * t257 + t192 * t259;
t283 = t184 * t264 - t190 * t261;
t282 = -t184 * t261 - t190 * t264;
t281 = t190 * t259 + t328;
t205 = -pkin(3) - t293;
t207 = (t258 * t265 + t318) * t331;
t280 = t205 * t250 + t207 * t253;
t206 = t258 * t304 + t224;
t235 = -pkin(3) - t335;
t279 = -t206 * t253 + t235 * t250;
t274 = t240 + t287;
t209 = t278 - t335;
t232 = pkin(2) * t258 + qJ(4);
t272 = t209 * t264 - t232 * t320;
t197 = -t233 * t320 - t234 * t264;
t199 = -t233 * t264 + t234 * t320;
t271 = g(2) * t199 - g(3) * t197 + (t283 * qJD(5) + t284) * t259 + t264 * t329;
t198 = t233 * t319 - t234 * t261;
t200 = t233 * t261 + t234 * t319;
t270 = -g(2) * t200 - g(3) * t198 + t179 * t321 + t310 * t328;
t269 = -g(2) * t233 + g(3) * t234 + t338;
t268 = t232 * t310 + t308 * t261;
t211 = t257 * t311 * t322;
t267 = t259 * t309 + (t211 + (t220 * t311 - t294 * t264) * t257) * MDP(14) + (t294 * t261 + (t220 + t322) * t310) * t257 * MDP(15) + t250 * MDP(4) + (0.2e1 * (t314 * t253 * qJD(5) - t250 * t317) * MDP(13) + (t250 * t255 - 0.2e1 * t261 * t299) * MDP(12)) * t251;
t266 = cos(qJ(1));
t263 = sin(qJ(1));
t249 = t253 ^ 2;
t248 = t266 * pkin(1);
t247 = t263 * pkin(1);
t204 = qJ(4) + t316;
t191 = -pkin(3) * t253 + t285;
t171 = t282 * qJD(5) + t296;
t1 = [(t280 * t257 - t301) * MDP(9) + qJDD(1) * MDP(1) + (t275 - t280) * t330 + (t183 * t205 + t191 * t207 - g(2) * (t248 + t300) - g(3) * (t247 + t286) + t338 * t204 + t281 * t202) * MDP(11) + ((t202 * t319 + t207 * t261) * t220 + (t193 * t261 + t204 * t319) * t217 + (t204 * t250 + t326) * t264 * t251 + (t220 * t327 + (-t328 + (-t220 * t259 - t251 * t253) * t204) * t261) * qJD(5) + t271) * MDP(18) + (-g(2) * t266 - g(3) * t263) * MDP(2) + (g(2) * t263 - g(3) * t266) * MDP(3) + (t187 * t316 + t196 * t276 + t186 * t293 - t195 * t207 - g(2) * (t238 + t248) - g(3) * (t237 + t247)) * MDP(7) + t267 + (((-qJDD(1) - t250) * t262 + t265 * t290) * pkin(1) + t298) * MDP(6) + ((t250 * t265 + t262 * t290) * pkin(1) + t274) * MDP(5) + (t204 * t297 + t315 * t326 + t269) * MDP(10) + (-(-t193 * t311 + t207 * t264) * t220 - t217 * t327 + (-(-t204 * t310 - t325) * t220 + t204 * t261 * t217 - t171) * t259 + (t253 * t325 + (t299 + t323) * t204) * t251 + t270) * MDP(17); (t291 * t337 + t274) * MDP(5) + ((t265 * t291 - t306) * pkin(1) + t298) * MDP(6) + (t195 * t206 - t196 * t208 + (t186 * t260 + t187 * t258 + t287) * pkin(2)) * MDP(7) + (t275 - t279) * t330 + (t279 * t257 - t301) * MDP(9) + (t308 * t253 * t315 + t232 * t297 + t269) * MDP(10) + (t183 * t235 - t191 * t206 - g(2) * t300 - g(3) * t286 + (t180 * t232 + t308 * t190) * t259 + (t179 * t232 + t308 * t189) * t257) * MDP(11) + (-t272 * t217 - t171 * t259 + (t264 * t206 + t209 * t311 + t268 * t259) * t220 + (t232 * t323 + t268 * t253) * t251 + t270) * MDP(17) + ((t209 * t261 + t232 * t319) * t217 + (-t261 * t206 + t259 * t292) * t220 + (-t189 * t321 + t272 * t220) * qJD(5) + (t232 * t264 * t250 + (-t232 * t311 + t292) * t253) * t251 + t271) * MDP(18) + t267; (qJDD(3) - g(1)) * MDP(7) + (-t179 * t259 - g(1)) * MDP(11) + t211 * MDP(18) + (t180 * MDP(11) + (-qJD(5) * t220 * MDP(18) + t295 * MDP(17)) * t261 + (t295 * MDP(18) + (t220 - t322) * MDP(17) * qJD(5)) * t264) * t257; (-t281 * t253 + t273 + t288) * MDP(11) + (-MDP(11) * pkin(3) + MDP(9) * t257 - t330) * t250 + (-MDP(17) * t264 + MDP(18) * t261) * t217 + (-t252 * MDP(10) + (-MDP(10) - t277) * t251) * t249 - t277 * t220 ^ 2; -t309 + (-g(2) * t197 - g(3) * t199 + t282 * t220 + t296) * MDP(17) + (g(2) * t198 - g(3) * t200 - t283 * t220 - t284) * MDP(18) + (t282 * MDP(17) - t283 * MDP(18)) * qJD(5) + (MDP(12) * t317 - t314 * MDP(13)) * t251 * t249 + ((MDP(14) * t264 - MDP(15) * t261) * t250 + t277 * g(1) + ((t307 * MDP(15) - t189 * MDP(17)) * t264 + (t307 * MDP(14) + t189 * MDP(18)) * t261) * t253) * t257;];
tau = t1;
