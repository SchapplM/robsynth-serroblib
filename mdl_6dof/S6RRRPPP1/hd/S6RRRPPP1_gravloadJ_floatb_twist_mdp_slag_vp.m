% Calculate Gravitation load on the joints for
% S6RRRPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha4,d1,d2,d3,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPP1_convert_par2_MPV_fixb.m
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPP1_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPP1_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPP1_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRRPPP1_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:18:31
% EndTime: 2019-03-09 15:18:36
% DurationCPUTime: 1.56s
% Computational Cost: add. (611->170), mult. (1670->247), div. (0->0), fcn. (1936->10), ass. (0->96)
t237 = sin(pkin(10));
t239 = cos(pkin(10));
t238 = sin(pkin(6));
t240 = cos(pkin(6));
t245 = cos(qJ(2));
t241 = sin(qJ(3));
t242 = sin(qJ(2));
t297 = t241 * t242;
t261 = t238 * t245 + t240 * t297;
t244 = cos(qJ(3));
t294 = t242 * t244;
t318 = -t237 * t261 + t239 * t294;
t323 = MDP(19) - MDP(24) - MDP(27);
t322 = MDP(23) - MDP(18) - MDP(28);
t243 = sin(qJ(1));
t246 = cos(qJ(1));
t271 = g(1) * t246 + g(2) * t243;
t321 = t242 * t271;
t317 = MDP(20) + MDP(22) + MDP(26);
t290 = t245 * t246;
t211 = t243 * t241 + t244 * t290;
t315 = g(1) * t211;
t314 = g(1) * t243;
t289 = t246 * t241;
t291 = t244 * t245;
t209 = t243 * t291 - t289;
t312 = g(2) * t209;
t234 = t245 * pkin(2);
t308 = qJ(4) * t238;
t292 = t243 * t245;
t208 = t241 * t292 + t244 * t246;
t307 = t208 * t238;
t210 = t243 * t244 - t245 * t289;
t306 = t210 * t238;
t305 = t237 * t240;
t304 = t238 * t242;
t302 = t239 * t240;
t301 = t239 * t241;
t300 = t240 * t242;
t299 = t240 * t244;
t298 = t240 * t245;
t296 = t241 * t245;
t295 = t242 * t243;
t293 = t242 * t246;
t282 = qJ(4) * t298;
t288 = pkin(9) * t292 + t243 * t282;
t287 = pkin(9) * t290 + t246 * t282;
t286 = MDP(25) + MDP(29);
t285 = g(3) * t294;
t284 = t209 * t308;
t283 = t211 * t308;
t280 = t238 * t296;
t279 = t238 * t295;
t276 = t239 * t304;
t225 = qJ(4) * t300;
t274 = -pkin(3) * t297 + t294 * t308;
t179 = -t208 * t237 + t209 * t302;
t180 = -t208 * t239 - t209 * t305;
t203 = t208 * pkin(3);
t269 = t180 * pkin(4) + qJ(5) * t179 - t203;
t181 = t210 * t237 + t211 * t302;
t182 = t210 * t239 - t211 * t305;
t205 = t210 * pkin(3);
t268 = t182 * pkin(4) + qJ(5) * t181 + t205;
t267 = pkin(3) * t291 + t242 * pkin(9) + qJ(4) * t280 + t225 + t234;
t266 = -t285 - t312;
t265 = -t209 * pkin(3) + t246 * pkin(8) - qJ(4) * t307;
t249 = t237 * t294 + t239 * t261;
t185 = t249 * t243;
t186 = t318 * t243;
t264 = -t186 * pkin(4) - qJ(5) * t185 + t288;
t187 = t249 * t246;
t188 = t318 * t246;
t263 = -t188 * pkin(4) - t187 * qJ(5) + t287;
t262 = t240 * t295 + t307;
t260 = t238 * t297 - t298;
t173 = t208 * t305 - t209 * t239 - t237 * t279;
t199 = (t237 * t241 - t239 * t299) * t242;
t200 = (t237 * t299 + t301) * t242;
t255 = -t200 * pkin(4) - qJ(5) * t199 + t274;
t190 = -t276 + (t237 * t244 + t240 * t301) * t245;
t191 = t237 * t304 + t239 * t291 - t296 * t305;
t254 = t191 * pkin(4) + qJ(5) * t190 + t267;
t172 = -t208 * t302 - t209 * t237 + t243 * t276;
t252 = t173 * pkin(4) + qJ(5) * t172 + t265;
t251 = pkin(2) * t290 + t211 * pkin(3) + t243 * pkin(8) + pkin(9) * t293 - qJ(4) * t306 + (pkin(1) + t225) * t246;
t250 = (-t234 - pkin(1) + (-qJ(4) * t240 - pkin(9)) * t242) * t314;
t174 = -t210 * t302 + t211 * t237 - t246 * t276;
t175 = t211 * t239 + (t210 * t240 + t238 * t293) * t237;
t248 = t175 * pkin(4) + qJ(5) * t174 + t251;
t247 = (pkin(3) * t244 + t241 * t308 + pkin(2)) * t321;
t207 = t280 + t300;
t202 = t260 * t246;
t201 = t260 * t243;
t193 = t240 * t293 - t306;
t1 = [t271 * MDP(3) + (g(1) * t209 - g(2) * t211) * MDP(16) + (-g(1) * t208 - g(2) * t210) * MDP(17) + (-g(1) * t265 - g(2) * t251 - t250) * MDP(21) + (-g(1) * t252 - g(2) * t248 - t250) * MDP(25) + (-g(1) * (-pkin(5) * t262 + qJ(6) * t173 + t252) - g(2) * (pkin(5) * t193 + qJ(6) * t175 + t248) - t250) * MDP(29) + (-MDP(10) * t242 + MDP(9) * t245 + MDP(2)) * (-g(2) * t246 + t314) + t317 * (g(1) * t262 - g(2) * t193) + t323 * (g(1) * t172 + g(2) * t174) + t322 * (g(1) * t173 + g(2) * t175); (g(3) * t242 + t245 * t271) * MDP(10) + (-g(1) * t287 - g(2) * t288 - g(3) * t267 + t247) * MDP(21) + (-g(1) * t263 - g(2) * t264 - g(3) * t254 + t247) * MDP(25) + (-g(1) * (-t202 * pkin(5) - t188 * qJ(6) + t263) - g(2) * (-pkin(5) * t201 - qJ(6) * t186 + t264) - g(3) * (pkin(5) * t207 + qJ(6) * t191 + t254) + t247) * MDP(29) + t317 * (g(1) * t202 + g(2) * t201 - g(3) * t207) - t322 * (g(1) * t188 + g(2) * t186 - g(3) * t191) - t323 * (g(1) * t187 + g(2) * t185 - g(3) * t190) + (t244 * MDP(16) - t241 * MDP(17) + MDP(9)) * (-g(3) * t245 + t321); (-g(1) * t210 + g(2) * t208 + g(3) * t297) * MDP(16) + (-g(1) * (t205 + t283) - g(2) * (-t203 + t284) - g(3) * t274) * MDP(21) + (-g(1) * (t268 + t283) - g(2) * (t269 + t284) - g(3) * t255) * MDP(25) + (-g(1) * (qJ(6) * t182 + t268) - g(2) * (qJ(6) * t180 + t269) - g(3) * (-qJ(6) * t200 + t255) + (-pkin(5) * t285 + (-t312 - t315) * (pkin(5) + qJ(4))) * t238) * MDP(29) + t323 * (g(1) * t181 + g(2) * t179 - g(3) * t199) + t322 * (g(1) * t182 + g(2) * t180 - g(3) * t200) + (-t238 * t317 + MDP(17)) * (-t266 + t315); (MDP(21) + t286) * (-g(1) * t193 - g(2) * t262 - g(3) * t260); t286 * (-g(1) * t174 + t266 * t237 + (-g(2) * (t208 * t240 - t279) - g(3) * t261) * t239); (-g(1) * t175 + g(2) * t173 - g(3) * t318) * MDP(29);];
taug  = t1;
