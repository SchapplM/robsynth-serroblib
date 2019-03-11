% Calculate joint inertia matrix for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR4_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:09:22
% EndTime: 2019-03-09 22:09:26
% DurationCPUTime: 1.08s
% Computational Cost: add. (1799->211), mult. (3379->296), div. (0->0), fcn. (3893->10), ass. (0->110)
t257 = sin(qJ(4));
t261 = cos(qJ(4));
t254 = sin(pkin(11));
t255 = cos(pkin(11));
t226 = -t254 * t257 + t255 * t261;
t227 = t254 * t261 + t255 * t257;
t256 = sin(qJ(6));
t260 = cos(qJ(6));
t198 = -t260 * t226 + t256 * t227;
t199 = t256 * t226 + t260 * t227;
t284 = t199 * MDP(29) - t198 * MDP(30);
t308 = t257 * MDP(20) + t261 * MDP(21) + t284;
t306 = pkin(8) + pkin(7);
t258 = sin(qJ(3));
t244 = t258 * pkin(2) + pkin(9);
t224 = (-qJ(5) - t244) * t257;
t251 = t261 * qJ(5);
t225 = t261 * t244 + t251;
t189 = t255 * t224 - t254 * t225;
t298 = t227 * pkin(10);
t173 = t189 - t298;
t190 = t254 * t224 + t255 * t225;
t220 = t226 * pkin(10);
t174 = t220 + t190;
t154 = t260 * t173 - t256 * t174;
t155 = t256 * t173 + t260 * t174;
t305 = t154 * MDP(32) - t155 * MDP(33);
t237 = (-qJ(5) - pkin(9)) * t257;
t238 = t261 * pkin(9) + t251;
t203 = t255 * t237 - t254 * t238;
t180 = t203 - t298;
t204 = t254 * t237 + t255 * t238;
t181 = t220 + t204;
t161 = t260 * t180 - t256 * t181;
t162 = t256 * t180 + t260 * t181;
t304 = t161 * MDP(32) - t162 * MDP(33);
t285 = t198 * MDP(32) + t199 * MDP(33);
t268 = t261 * MDP(23) - t257 * MDP(24);
t300 = cos(qJ(2));
t247 = -t300 * pkin(2) - pkin(1);
t303 = 0.2e1 * t247;
t302 = 2 * MDP(25);
t301 = -2 * MDP(28);
t299 = pkin(4) * t254;
t262 = cos(qJ(3));
t297 = t262 * pkin(2);
t295 = MDP(26) * pkin(4);
t259 = sin(qJ(2));
t239 = t306 * t259;
t240 = t306 * t300;
t205 = t262 * t239 + t258 * t240;
t294 = t205 * t261;
t234 = t258 * t300 + t262 * t259;
t293 = t234 * t257;
t290 = t257 * t261;
t206 = -t258 * t239 + t262 * t240;
t289 = t261 * t206;
t233 = t258 * t259 - t262 * t300;
t197 = t233 * pkin(3) - t234 * pkin(9) + t247;
t171 = t261 * t197 - t257 * t206;
t158 = t233 * pkin(4) - t234 * t251 + t171;
t168 = t289 + (-qJ(5) * t234 + t197) * t257;
t147 = t255 * t158 - t254 * t168;
t148 = t254 * t158 + t255 * t168;
t288 = -t147 * t227 + t148 * t226;
t287 = -t189 * t227 + t190 * t226;
t286 = -t203 * t227 + t204 * t226;
t185 = t227 * t234;
t186 = t226 * t234;
t167 = -t256 * t185 + t260 * t186;
t282 = MDP(27) * t167;
t166 = t260 * t185 + t256 * t186;
t164 = t166 * MDP(30);
t165 = t167 * MDP(29);
t242 = t255 * pkin(4) + pkin(5);
t214 = t260 * t242 - t256 * t299;
t281 = t214 * MDP(32);
t215 = t256 * t242 + t260 * t299;
t280 = t215 * MDP(33);
t277 = 0.2e1 * t300;
t276 = MDP(22) + MDP(31);
t275 = t233 * MDP(31) - t164 + t165;
t246 = -t261 * pkin(4) - pkin(3);
t274 = MDP(19) * t290;
t142 = t233 * pkin(5) - t186 * pkin(10) + t147;
t144 = -t185 * pkin(10) + t148;
t139 = t260 * t142 - t256 * t144;
t273 = -pkin(3) * t234 - pkin(9) * t233;
t179 = pkin(4) * t293 + t205;
t272 = (t226 * t254 - t227 * t255) * pkin(4) * MDP(25) + t308;
t252 = t257 ^ 2;
t271 = t252 * MDP(18) + MDP(15) + 0.2e1 * t274 + (MDP(27) * t199 + t198 * t301) * t199;
t140 = t256 * t142 + t260 * t144;
t245 = -pkin(3) - t297;
t270 = -t233 * t244 + t234 * t245;
t208 = -t226 * pkin(5) + t246;
t269 = t261 * MDP(20) - t257 * MDP(21);
t267 = -t257 * MDP(23) - t261 * MDP(24);
t266 = (t262 * MDP(16) - t258 * MDP(17)) * pkin(2);
t265 = 0.2e1 * t285;
t253 = t261 ^ 2;
t264 = (-t199 * t166 - t167 * t198) * MDP(28) + t199 * t282 - t205 * MDP(16) - t206 * MDP(17) + ((-t252 + t253) * MDP(19) + MDP(18) * t290 + MDP(13)) * t234 + (-MDP(14) + t308) * t233;
t236 = t246 - t297;
t207 = t208 - t297;
t200 = t205 * t257;
t172 = t257 * t197 + t289;
t169 = t185 * pkin(5) + t179;
t151 = t169 * t199;
t150 = t169 * t198;
t1 = [(t147 ^ 2 + t148 ^ 2 + t179 ^ 2) * MDP(26) + pkin(1) * MDP(9) * t277 + t234 * MDP(17) * t303 + MDP(1) + t276 * t233 ^ 2 + (t166 * t301 + t282) * t167 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t259 + MDP(5) * t277) * t259 + (t253 * MDP(18) + MDP(11) - 0.2e1 * t274) * t234 ^ 2 + (MDP(16) * t303 + 0.2e1 * t165 - 0.2e1 * t164 + 0.2e1 * (-MDP(12) + t269) * t234) * t233 + 0.2e1 * (-t172 * t233 + t234 * t294) * MDP(24) + 0.2e1 * (t171 * t233 + t205 * t293) * MDP(23) + (-t147 * t186 - t148 * t185) * t302 + 0.2e1 * (t139 * t233 + t169 * t166) * MDP(32) + 0.2e1 * (-t140 * t233 + t169 * t167) * MDP(33); (-t300 * MDP(10) - t259 * MDP(9)) * pkin(7) + (t154 * t233 + t207 * t166 + t150) * MDP(32) + (-t190 * t185 - t189 * t186 + t288) * MDP(25) + (t270 * t257 - t294) * MDP(23) + (t270 * t261 + t200) * MDP(24) + t300 * MDP(7) + (-t155 * t233 + t207 * t167 + t151) * MDP(33) + (t147 * t189 + t148 * t190 + t179 * t236) * MDP(26) + t259 * MDP(6) + t264; MDP(8) + t287 * t302 + (t189 ^ 2 + t190 ^ 2 + t236 ^ 2) * MDP(26) + t207 * t265 + t271 - 0.2e1 * t268 * t245 + 0.2e1 * t266; (t161 * t233 + t208 * t166 + t150) * MDP(32) + (-t162 * t233 + t208 * t167 + t151) * MDP(33) + (-t204 * t185 - t203 * t186 + t288) * MDP(25) + (t273 * t257 - t294) * MDP(23) + (t273 * t261 + t200) * MDP(24) + (t147 * t203 + t148 * t204 + t179 * t246) * MDP(26) + t264; (t286 + t287) * MDP(25) + (t189 * t203 + t190 * t204 + t236 * t246) * MDP(26) + t266 + t271 + t268 * (pkin(3) - t245) + t285 * (t207 + t208); t286 * t302 + (t203 ^ 2 + t204 ^ 2 + t246 ^ 2) * MDP(26) + t208 * t265 + 0.2e1 * t268 * pkin(3) + t271; t233 * MDP(22) + t171 * MDP(23) - t172 * MDP(24) + (t214 * t233 + t139) * MDP(32) + (-t215 * t233 - t140) * MDP(33) + t269 * t234 + ((-t185 * t254 - t186 * t255) * MDP(25) + (t147 * t255 + t148 * t254) * MDP(26)) * pkin(4) + t275; t267 * t244 + (t189 * t255 + t190 * t254) * t295 + t272 + t305; t267 * pkin(9) + (t203 * t255 + t204 * t254) * t295 + t272 + t304; (t254 ^ 2 + t255 ^ 2) * MDP(26) * pkin(4) ^ 2 + 0.2e1 * t281 - 0.2e1 * t280 + t276; t179 * MDP(26) + t166 * MDP(32) + t167 * MDP(33); t236 * MDP(26) + t285; t246 * MDP(26) + t285; 0; MDP(26); t139 * MDP(32) - t140 * MDP(33) + t275; t284 + t305; t284 + t304; MDP(31) - t280 + t281; 0; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
