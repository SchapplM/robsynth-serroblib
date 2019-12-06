% Calculate vector of inverse dynamics joint torques for
% S5PPRRR3
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PPRRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:17:07
% EndTime: 2019-12-05 15:17:10
% DurationCPUTime: 2.16s
% Computational Cost: add. (745->224), mult. (1669->328), div. (0->0), fcn. (1325->12), ass. (0->117)
t236 = cos(qJ(4));
t220 = -pkin(4) * t236 - pkin(3);
t234 = sin(qJ(3));
t237 = cos(qJ(3));
t228 = sin(pkin(9));
t292 = qJD(1) * t228;
t312 = qJD(2) * t237 - t234 * t292;
t182 = qJD(3) * t220 - t312;
t315 = t182 + t312;
t224 = qJD(4) + qJD(5);
t309 = pkin(6) + pkin(7);
t232 = sin(qJ(5));
t233 = sin(qJ(4));
t235 = cos(qJ(5));
t199 = t232 * t236 + t233 * t235;
t313 = t224 * t199;
t314 = qJD(3) * t313;
t276 = qJDD(3) * t236;
t277 = qJDD(3) * t233;
t165 = t232 * t277 - t235 * t276 + t314;
t198 = t232 * t233 - t235 * t236;
t249 = t224 * t198;
t311 = t224 * t236;
t229 = sin(pkin(8));
t231 = cos(pkin(8));
t230 = cos(pkin(9));
t273 = t230 * qJDD(1);
t310 = (g(1) * t231 + g(2) * t229) * t228 + t273;
t308 = qJD(3) * pkin(3);
t307 = qJDD(3) * pkin(3);
t291 = qJD(1) * t230;
t209 = t233 * t291;
t271 = t237 * t292;
t197 = qJD(2) * t234 + t271;
t265 = qJD(3) * t309 + t197;
t171 = t236 * t265 - t209;
t306 = t171 * t235;
t223 = qJDD(4) + qJDD(5);
t305 = t198 * t223;
t304 = t199 * t223;
t303 = t228 * t229;
t302 = t228 * t231;
t301 = t228 * t237;
t300 = t229 * t234;
t299 = t229 * t237;
t298 = t231 * t234;
t297 = t231 * t237;
t239 = qJD(3) ^ 2;
t296 = t237 * t239;
t295 = qJDD(1) - g(3);
t225 = t233 ^ 2;
t294 = -t236 ^ 2 + t225;
t238 = qJD(4) ^ 2;
t293 = t238 + t239;
t289 = qJD(3) * t233;
t288 = qJD(3) * t234;
t287 = qJD(3) * t236;
t286 = qJD(3) * t237;
t285 = qJD(4) * t233;
t283 = qJD(5) * t232;
t282 = qJD(5) * t235;
t280 = qJD(3) * qJD(4);
t279 = qJDD(1) * t228;
t278 = qJDD(2) * t234;
t275 = qJDD(3) * t237;
t274 = qJDD(4) * t233;
t272 = qJD(4) * t309;
t270 = t228 * t288;
t269 = t232 * t289;
t268 = t235 * t287;
t267 = t233 * t280;
t266 = t236 * t280;
t177 = qJDD(3) * pkin(6) + qJD(3) * t312 + t237 * t279 + t278;
t264 = pkin(7) * qJDD(3) + t177;
t263 = pkin(4) * t285 - t197;
t186 = t230 * t299 - t298;
t188 = t230 * t297 + t300;
t262 = g(1) * t188 + g(2) * t186;
t170 = -t233 * t265 - t236 * t291;
t169 = qJD(4) * pkin(4) + t170;
t259 = -t169 * t232 - t306;
t191 = -t230 * t236 - t233 * t301;
t192 = -t230 * t233 + t236 * t301;
t258 = t191 * t235 - t192 * t232;
t257 = t191 * t232 + t192 * t235;
t202 = t309 * t233;
t203 = t309 * t236;
t256 = -t202 * t235 - t203 * t232;
t255 = -t202 * t232 + t203 * t235;
t254 = -t234 * t239 + t275;
t253 = qJDD(3) * t234 + t296;
t252 = qJD(2) * t288 + qJD(3) * t271 - qJDD(2) * t237 + t234 * t279;
t248 = t267 - t276;
t164 = qJD(5) * t268 - t224 * t269 + t232 * t276 + (t266 + t277) * t235;
t193 = -t268 + t269;
t195 = -t232 * t287 - t235 * t289;
t247 = -t195 * t193 * MDP(13) + (t193 * t224 + t164) * MDP(15) + (-t195 * t224 - t165) * MDP(16) + (-t193 ^ 2 + t195 ^ 2) * MDP(14) + t223 * MDP(17);
t246 = g(3) * t228 * t234 - g(1) * (-t230 * t298 + t299) - g(2) * (-t230 * t300 - t297);
t189 = -t312 - t308;
t245 = -qJD(3) * t189 - t177 + t262;
t176 = t252 - t307;
t244 = -pkin(6) * qJDD(4) + (t189 + t312 - t308) * qJD(4);
t243 = qJD(3) * t197 + t246;
t160 = qJDD(4) * pkin(4) + t209 * qJD(4) - t264 * t233 + (-qJD(4) * t265 - t273) * t236;
t227 = qJ(4) + qJ(5);
t221 = sin(t227);
t222 = cos(t227);
t242 = t171 * t283 + t182 * t193 + (-t171 * t224 - t160) * t232 - g(1) * (-t188 * t222 - t221 * t302) - g(2) * (-t186 * t222 - t221 * t303) - g(3) * (t221 * t230 - t222 * t301);
t241 = -pkin(6) * t238 - t176 + t243 + t307;
t161 = qJD(4) * t170 - t233 * t273 + t236 * t264;
t240 = -g(1) * (-t188 * t221 + t222 * t302) - g(2) * (-t186 * t221 + t222 * t303) - g(3) * (-t221 * t301 - t222 * t230) + t259 * qJD(5) + t235 * t160 - t232 * t161 + t182 * t195;
t201 = t236 * t272;
t200 = t233 * t272;
t181 = qJD(4) * t191 - t236 * t270;
t180 = -qJD(4) * t192 + t233 * t270;
t167 = pkin(4) * t248 + t176;
t1 = [t295 * MDP(1) + (-g(3) + (t228 ^ 2 + t230 ^ 2) * qJDD(1)) * MDP(2) - t253 * t228 * MDP(4) - t254 * t228 * MDP(5) + (qJD(4) * t180 + qJDD(4) * t191 + (t234 * t248 - t236 * t296) * t228) * MDP(11) + (-qJD(4) * t181 - qJDD(4) * t192 + (t233 * t253 + t234 * t266) * t228) * MDP(12) + ((-qJD(5) * t257 + t180 * t235 - t181 * t232) * t224 + t258 * t223 + (t165 * t234 + t193 * t286) * t228) * MDP(18) + (-(qJD(5) * t258 + t180 * t232 + t181 * t235) * t224 - t257 * t223 + (t164 * t234 - t195 * t286) * t228) * MDP(19); (-g(1) * t229 + g(2) * t231 + qJDD(2)) * MDP(2) + t254 * MDP(4) - t253 * MDP(5) + ((-0.2e1 * t267 + t276) * t237 + (-t236 * t293 - t274) * t234) * MDP(11) + ((-qJDD(4) * t234 - 0.2e1 * t237 * t280) * t236 + (t234 * t293 - t275) * t233) * MDP(12) + ((-t165 - t314) * t237 + ((t232 * t285 + t233 * t283 - t235 * t311) * t224 - t304 + qJD(3) * t193) * t234) * MDP(18) + ((qJD(3) * t249 - t164) * t237 + (-(-t232 * t311 - t233 * t282 - t235 * t285) * t224 + t305 - qJD(3) * t195) * t234) * MDP(19); qJDD(3) * MDP(3) + (t243 - t252) * MDP(4) + (-t295 * t301 + t262 - t278) * MDP(5) + (qJDD(3) * t225 + 0.2e1 * t233 * t266) * MDP(6) + 0.2e1 * (t233 * t276 - t280 * t294) * MDP(7) + (t236 * t238 + t274) * MDP(8) + (qJDD(4) * t236 - t233 * t238) * MDP(9) + (t233 * t244 + t236 * t241) * MDP(11) + (-t233 * t241 + t236 * t244) * MDP(12) + (t164 * t199 + t195 * t249) * MDP(13) + (-t164 * t198 - t165 * t199 + t193 * t249 + t195 * t313) * MDP(14) + (-t224 * t249 + t304) * MDP(15) + (-t224 * t313 - t305) * MDP(16) + ((-qJD(5) * t255 + t200 * t232 - t201 * t235) * t224 + t256 * t223 + t220 * t165 + t167 * t198 + t263 * t193 + t246 * t222 + t315 * t313) * MDP(18) + (-(qJD(5) * t256 - t200 * t235 - t201 * t232) * t224 - t255 * t223 + t220 * t164 + t167 * t199 - t263 * t195 - t246 * t221 - t315 * t249) * MDP(19); MDP(8) * t277 + MDP(9) * t276 + qJDD(4) * MDP(10) + (-g(3) * t191 + t245 * t233 - t236 * t310) * MDP(11) + (g(3) * t192 + t233 * t310 + t245 * t236) * MDP(12) + (-(-t170 * t232 - t306) * t224 + (-t193 * t289 + t223 * t235 - t224 * t283) * pkin(4) + t240) * MDP(18) + ((-qJD(5) * t169 + t170 * t224 - t161) * t235 + (t195 * t289 - t223 * t232 - t224 * t282) * pkin(4) + t242) * MDP(19) + t247 + (-MDP(6) * t233 * t236 + MDP(7) * t294) * t239; (-t224 * t259 + t240) * MDP(18) + ((-t161 + (-qJD(5) + t224) * t169) * t235 + t242) * MDP(19) + t247;];
tau = t1;
