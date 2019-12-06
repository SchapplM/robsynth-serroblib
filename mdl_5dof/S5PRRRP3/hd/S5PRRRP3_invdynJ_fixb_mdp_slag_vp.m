% Calculate vector of inverse dynamics joint torques for
% S5PRRRP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:19
% EndTime: 2019-12-05 16:44:22
% DurationCPUTime: 1.45s
% Computational Cost: add. (1194->215), mult. (2563->279), div. (0->0), fcn. (1702->8), ass. (0->113)
t252 = cos(qJ(3));
t313 = pkin(6) + pkin(7);
t220 = t313 * t252;
t251 = sin(qJ(3));
t293 = qJD(1) * t251;
t192 = qJD(2) * t220 + t293;
t250 = sin(qJ(4));
t186 = t250 * t192;
t285 = qJD(2) * t313;
t191 = t252 * qJD(1) - t251 * t285;
t307 = qJD(3) * pkin(3);
t189 = t191 + t307;
t312 = cos(qJ(4));
t278 = t312 * t189 - t186;
t204 = t250 * t252 + t312 * t251;
t197 = t204 * qJD(2);
t304 = t197 * qJ(5);
t317 = t304 - t278;
t249 = qJ(3) + qJ(4);
t240 = sin(t249);
t241 = cos(t249);
t245 = pkin(8) + qJ(2);
t236 = sin(t245);
t237 = cos(t245);
t272 = g(1) * t237 + g(2) * t236;
t316 = -g(3) * t241 + t272 * t240;
t246 = qJD(3) + qJD(4);
t242 = t252 * pkin(3);
t308 = pkin(2) + t242;
t219 = t313 * t251;
t296 = -t250 * t219 + t312 * t220;
t314 = t197 ^ 2;
t283 = t312 * t252;
t273 = qJD(2) * t283;
t292 = qJD(2) * t251;
t282 = t250 * t292;
t195 = -t273 + t282;
t218 = t308 * qJD(2);
t183 = pkin(4) * t195 + qJD(5) - t218;
t306 = t183 * t197;
t305 = t195 * qJ(5);
t303 = t236 * t241;
t302 = t237 * t241;
t301 = t250 * t251;
t300 = qJDD(1) - g(3);
t165 = pkin(4) * t246 - t317;
t299 = t165 + t317;
t182 = t246 * t204;
t279 = qJDD(2) * t312;
t289 = qJDD(2) * t251;
t269 = t250 * t289 - t252 * t279;
t173 = t182 * qJD(2) + t269;
t270 = t246 * t301;
t281 = qJD(4) * t312;
t181 = -qJD(3) * t283 - t252 * t281 + t270;
t298 = -t204 * t173 + t181 * t195;
t297 = t312 * t191 - t186;
t295 = pkin(4) * t241 + t242;
t247 = t251 ^ 2;
t294 = -t252 ^ 2 + t247;
t291 = qJD(4) * t250;
t290 = qJD(2) * qJD(3);
t288 = qJDD(2) * t252;
t287 = t251 * t307;
t284 = qJD(3) * t313;
t188 = t312 * t192;
t280 = t251 * t290;
t277 = -t191 * t250 - t188;
t276 = -t312 * t219 - t220 * t250;
t275 = -t246 * t273 - t250 * t288 - t251 * t279;
t271 = g(1) * t236 - g(2) * t237;
t172 = qJD(2) * t270 + t275;
t203 = -t283 + t301;
t268 = -t172 * t203 + t182 * t197;
t243 = qJDD(3) + qJDD(4);
t267 = t181 * t246 - t204 * t243;
t266 = -0.2e1 * pkin(2) * t290 - pkin(6) * qJDD(3);
t265 = -t250 * t189 - t188;
t193 = pkin(3) * t280 - qJDD(2) * t308;
t209 = t251 * t284;
t210 = t252 * t284;
t264 = -t312 * t209 - t250 * t210 - t219 * t281 - t220 * t291;
t194 = t195 ^ 2;
t263 = t197 * t195 * MDP(12) + (-t275 + (t195 - t282) * t246) * MDP(14) - t269 * MDP(15) + (-t194 + t314) * MDP(13) + t243 * MDP(16);
t253 = qJD(3) ^ 2;
t262 = 0.2e1 * qJDD(2) * pkin(2) - pkin(6) * t253 + t271;
t254 = qJD(2) ^ 2;
t261 = pkin(2) * t254 - pkin(6) * qJDD(2) + t272;
t260 = t173 * pkin(4) + qJDD(5) + t193;
t238 = t252 * qJDD(1);
t177 = qJDD(3) * pkin(3) + t238 - qJDD(2) * t219 + (-t252 * t285 - t293) * qJD(3);
t180 = t191 * qJD(3) + t251 * qJDD(1) + qJDD(2) * t220;
t259 = t265 * qJD(4) + t312 * t177 - t250 * t180;
t258 = -t296 * qJD(4) + t250 * t209 - t312 * t210;
t257 = t250 * t177 + t312 * t180 + t189 * t281 - t192 * t291;
t256 = g(1) * t302 + g(2) * t303 + g(3) * t240 - t218 * t195 - t257;
t255 = t218 * t197 + t259 + t316;
t244 = -qJ(5) - t313;
t234 = t312 * pkin(3) + pkin(4);
t217 = qJDD(3) * t252 - t251 * t253;
t216 = qJDD(3) * t251 + t252 * t253;
t208 = pkin(2) + t295;
t179 = -qJ(5) * t203 + t296;
t178 = -qJ(5) * t204 + t276;
t171 = -t182 * t246 - t203 * t243;
t170 = -t304 + t297;
t169 = t277 + t305;
t168 = -t265 - t305;
t162 = t181 * qJ(5) - t204 * qJD(5) + t258;
t161 = -qJ(5) * t182 - qJD(5) * t203 + t264;
t160 = -t173 * qJ(5) - t195 * qJD(5) + t257;
t159 = t243 * pkin(4) + t172 * qJ(5) - t197 * qJD(5) + t259;
t1 = [t300 * MDP(1) + t217 * MDP(10) - t216 * MDP(11) + t171 * MDP(17) + t267 * MDP(18) + (t268 + t298) * MDP(19) + (-t159 * t203 + t160 * t204 - t165 * t182 - t168 * t181 - g(3)) * MDP(20); qJDD(2) * MDP(2) + t271 * MDP(3) + t272 * MDP(4) + (qJDD(2) * t247 + 0.2e1 * t252 * t280) * MDP(5) + 0.2e1 * (t251 * t288 - t294 * t290) * MDP(6) + t216 * MDP(7) + t217 * MDP(8) + (t266 * t251 + t262 * t252) * MDP(10) + (-t262 * t251 + t266 * t252) * MDP(11) + (-t172 * t204 - t181 * t197) * MDP(12) + (-t268 + t298) * MDP(13) - t267 * MDP(14) + t171 * MDP(15) + (g(1) * t303 - g(2) * t302 - t173 * t308 - t218 * t182 + t193 * t203 + t195 * t287 + t276 * t243 + t258 * t246) * MDP(17) + (t172 * t308 + t218 * t181 + t193 * t204 + t197 * t287 - t271 * t240 - t296 * t243 - t264 * t246) * MDP(18) + (-t159 * t204 - t160 * t203 - t161 * t195 - t162 * t197 + t165 * t181 - t168 * t182 + t172 * t178 - t173 * t179 - t272) * MDP(19) + (t160 * t179 + t168 * t161 + t159 * t178 + t165 * t162 + t260 * (pkin(4) * t203 - t308) + t183 * (pkin(4) * t182 + t287) - g(1) * (-t208 * t236 - t237 * t244) - g(2) * (t208 * t237 - t236 * t244)) * MDP(20); MDP(7) * t289 + MDP(8) * t288 + qJDD(3) * MDP(9) + (-g(3) * t252 + t261 * t251 + t238) * MDP(10) + (-t300 * t251 + t261 * t252) * MDP(11) + (-t277 * t246 + (-t195 * t292 + t312 * t243 - t246 * t291) * pkin(3) + t255) * MDP(17) + (t297 * t246 + (-t197 * t292 - t250 * t243 - t246 * t281) * pkin(3) + t256) * MDP(18) + (t234 * t172 + (t168 + t169) * t197 + (-t165 + t170) * t195 + (-t173 * t250 + (-t312 * t195 + t197 * t250) * qJD(4)) * pkin(3)) * MDP(19) + (t159 * t234 - t168 * t170 - t165 * t169 - pkin(4) * t306 - g(3) * t295 - t272 * (-pkin(3) * t251 - pkin(4) * t240) + (-t183 * t292 + t160 * t250 + (-t165 * t250 + t312 * t168) * qJD(4)) * pkin(3)) * MDP(20) + t263 + (-t251 * t252 * MDP(5) + t294 * MDP(6)) * t254; (-t265 * t246 + t255) * MDP(17) + (t278 * t246 + t256) * MDP(18) + (pkin(4) * t172 - t299 * t195) * MDP(19) + (t299 * t168 + (t159 - t306 + t316) * pkin(4)) * MDP(20) + t263; (-t194 - t314) * MDP(19) + (t165 * t197 + t168 * t195 + t260 - t271) * MDP(20);];
tau = t1;
