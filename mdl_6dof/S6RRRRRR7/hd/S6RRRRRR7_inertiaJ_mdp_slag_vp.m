% Calculate joint inertia matrix for
% S6RRRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR7_inertiaJ_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:41:11
% EndTime: 2019-03-10 04:41:17
% DurationCPUTime: 2.35s
% Computational Cost: add. (3007->317), mult. (6815->436), div. (0->0), fcn. (7800->12), ass. (0->156)
t264 = sin(qJ(3));
t269 = cos(qJ(3));
t242 = -pkin(3) * t269 - pkin(10) * t264 - pkin(2);
t268 = cos(qJ(4));
t237 = t268 * t242;
t329 = pkin(11) * t264;
t263 = sin(qJ(4));
t332 = pkin(9) * t263;
t204 = -t268 * t329 + t237 + (-pkin(4) - t332) * t269;
t330 = pkin(9) * t269;
t296 = t268 * t330;
t208 = t296 + (t242 - t329) * t263;
t262 = sin(qJ(5));
t267 = cos(qJ(5));
t182 = t267 * t204 - t208 * t262;
t183 = t204 * t262 + t208 * t267;
t239 = t262 * t268 + t263 * t267;
t225 = t239 * t264;
t218 = t225 * MDP(28);
t238 = t262 * t263 - t267 * t268;
t226 = t238 * t264;
t219 = t226 * MDP(27);
t172 = -pkin(5) * t269 + pkin(12) * t226 + t182;
t177 = -pkin(12) * t225 + t183;
t261 = sin(qJ(6));
t266 = cos(qJ(6));
t159 = t266 * t172 - t177 * t261;
t160 = t172 * t261 + t177 * t266;
t194 = t266 * t225 - t226 * t261;
t192 = t194 * MDP(35);
t195 = -t225 * t261 - t226 * t266;
t193 = t195 * MDP(34);
t359 = t159 * MDP(37) - t160 * MDP(38) - t192 + t193;
t361 = t182 * MDP(30) - t183 * MDP(31) - t218 - t219 + t359;
t215 = -t263 * t330 + t237;
t216 = t242 * t263 + t296;
t360 = t215 * MDP(23) - t216 * MDP(24) + t361;
t340 = pkin(10) + pkin(11);
t243 = t340 * t263;
t244 = t340 * t268;
t211 = -t267 * t243 - t244 * t262;
t212 = -t243 * t262 + t244 * t267;
t196 = -pkin(12) * t239 + t211;
t197 = -pkin(12) * t238 + t212;
t201 = t266 * t238 + t239 * t261;
t202 = -t238 * t261 + t239 * t266;
t291 = t202 * MDP(34) - t201 * MDP(35) + (t196 * t266 - t197 * t261) * MDP(37) - (t196 * t261 + t197 * t266) * MDP(38);
t274 = t239 * MDP(27) - t238 * MDP(28) + t211 * MDP(30) - t212 * MDP(31) + t291;
t358 = t274 - (t263 * MDP(23) + t268 * MDP(24)) * pkin(10) + t263 * MDP(20) + t268 * MDP(21);
t357 = -MDP(14) + t358;
t355 = -0.2e1 * t264;
t300 = t267 * MDP(30);
t354 = pkin(4) * (-MDP(31) * t262 + t300);
t315 = MDP(38) * t261;
t275 = (MDP(37) * t266 - t315) * pkin(5);
t259 = sin(pkin(6));
t270 = cos(qJ(2));
t320 = t259 * t270;
t260 = cos(pkin(6));
t265 = sin(qJ(2));
t321 = t259 * t265;
t231 = t260 * t264 + t269 * t321;
t206 = t231 * t263 + t268 * t320;
t207 = t231 * t268 - t263 * t320;
t184 = t267 * t206 + t207 * t262;
t185 = -t206 * t262 + t207 * t267;
t167 = t266 * t184 + t185 * t261;
t168 = -t184 * t261 + t185 * t266;
t353 = t168 * MDP(34) - t167 * MDP(35);
t352 = t185 * MDP(27) - t184 * MDP(28);
t246 = pkin(8) * t321;
t338 = pkin(1) * t270;
t220 = t246 + (-pkin(2) - t338) * t260;
t230 = -t260 * t269 + t264 * t321;
t187 = t230 * pkin(3) - t231 * pkin(10) + t220;
t297 = pkin(8) * t320;
t339 = pkin(1) * t265;
t221 = t297 + (pkin(9) + t339) * t260;
t222 = (-pkin(2) * t270 - pkin(9) * t265 - pkin(1)) * t259;
t191 = t269 * t221 + t264 * t222;
t189 = -pkin(10) * t320 + t191;
t169 = t268 * t187 - t189 * t263;
t170 = t187 * t263 + t189 * t268;
t351 = -t169 * MDP(23) + t170 * MDP(24);
t348 = 0.2e1 * MDP(23);
t347 = 0.2e1 * MDP(24);
t346 = -2 * MDP(26);
t345 = 0.2e1 * MDP(30);
t344 = 0.2e1 * MDP(31);
t343 = -2 * MDP(33);
t342 = 0.2e1 * MDP(37);
t341 = 0.2e1 * MDP(38);
t337 = pkin(4) * t230;
t336 = pkin(4) * t262;
t335 = pkin(4) * t267;
t334 = pkin(5) * t230;
t333 = pkin(5) * t266;
t331 = pkin(9) * t268;
t328 = MDP(17) * pkin(2);
t327 = pkin(2) * MDP(16);
t326 = pkin(9) * MDP(17);
t162 = -pkin(11) * t207 + t169 + t337;
t164 = -pkin(11) * t206 + t170;
t324 = t164 * t267;
t158 = t162 * t262 + t324;
t156 = -pkin(12) * t184 + t158;
t325 = t156 * t266;
t190 = -t264 * t221 + t222 * t269;
t188 = pkin(3) * t320 - t190;
t323 = t188 * t263;
t322 = t188 * t268;
t319 = t260 * MDP(8);
t318 = t265 * MDP(6);
t241 = (pkin(4) * t263 + pkin(9)) * t264;
t316 = MDP(15) * t270;
t312 = t168 * MDP(32);
t307 = t185 * MDP(25);
t306 = t195 * MDP(32);
t305 = t207 * MDP(18);
t304 = t226 * MDP(25);
t251 = pkin(5) + t335;
t245 = t266 * t251;
t227 = -t261 * t336 + t245;
t303 = t227 * MDP(37);
t228 = t251 * t261 + t266 * t336;
t302 = t228 * MDP(38);
t301 = t231 * MDP(13);
t299 = t268 * MDP(18);
t298 = MDP(29) + MDP(36);
t295 = t230 * MDP(36) + t353;
t252 = -pkin(4) * t268 - pkin(3);
t294 = MDP(19) * t263 * t268;
t293 = MDP(22) + t298;
t292 = pkin(9) * MDP(16) - MDP(13);
t157 = t267 * t162 - t164 * t262;
t155 = -pkin(12) * t185 + t157 + t334;
t152 = t266 * t155 - t156 * t261;
t153 = t155 * t261 + t325;
t290 = t207 * MDP(20) - t206 * MDP(21);
t289 = t268 * MDP(20) - t263 * MDP(21);
t285 = t157 * MDP(30) - t158 * MDP(31);
t282 = -t152 * MDP(37) + t153 * MDP(38);
t279 = t230 * MDP(29) + t295 + t352;
t278 = -MDP(12) + t289;
t277 = MDP(36) - t302 + t303;
t178 = pkin(4) * t206 + t188;
t272 = t231 * MDP(12) - t290 - t352 - t353;
t257 = t268 ^ 2;
t255 = t263 ^ 2;
t254 = t259 ^ 2;
t234 = t260 * t339 + t297;
t233 = t260 * t338 - t246;
t217 = pkin(5) * t238 + t252;
t205 = pkin(5) * t225 + t241;
t163 = pkin(5) * t184 + t178;
t1 = [t254 * t265 ^ 2 * MDP(4) + t231 ^ 2 * MDP(11) + MDP(1) + (0.2e1 * t259 * t318 + t319) * t260 + (-0.2e1 * t206 * MDP(19) + t305) * t207 + (t184 * t346 + t307) * t185 + (t167 * t343 + t312) * t168 + t293 * t230 ^ 2 + (0.2e1 * MDP(5) * t265 + t316) * t254 * t270 + 0.2e1 * (t233 * t260 + t254 * t338) * MDP(9) + 0.2e1 * (t191 * t320 + t220 * t231) * MDP(17) + 0.2e1 * (-t190 * t320 + t220 * t230) * MDP(16) + 0.2e1 * (-t234 * t260 - t254 * t339) * MDP(10) + (-t170 * t230 + t188 * t207) * t347 + (-t158 * t230 + t178 * t185) * t344 + (-t153 * t230 + t163 * t168) * t341 + (t169 * t230 + t188 * t206) * t348 + (t157 * t230 + t178 * t184) * t345 + (t152 * t230 + t163 * t167) * t342 + 0.2e1 * (MDP(7) * t260 - t301) * t320 + 0.2e1 * (MDP(14) * t320 - t272) * t230; -t231 * t328 - t185 * t304 + t168 * t306 + (-t167 * t195 - t168 * t194) * MDP(33) + (t184 * t226 - t185 * t225) * MDP(26) + t233 * MDP(9) - t234 * MDP(10) + t319 + (t163 * t194 + t167 * t205) * MDP(37) + (t163 * t195 + t205 * t168) * MDP(38) + (t178 * t225 + t184 * t241) * MDP(30) + (-t178 * t226 + t185 * t241) * MDP(31) + (t270 * MDP(7) + t318) * t259 + (-t327 + t360) * t230 + (-t220 * MDP(16) + (-MDP(14) + t326) * t320 - t293 * t230 + t272 + t282 - t285 + t351) * t269 + (t207 * t299 + t231 * MDP(11) + (-t206 * t268 - t207 * t263) * MDP(19) + (pkin(9) * t206 + t323) * MDP(23) + (pkin(9) * t207 + t322) * MDP(24) + t220 * MDP(17) + t292 * t320 + t278 * t230) * t264; t194 * t205 * t342 + t225 * t241 * t345 + t328 * t355 + MDP(8) - (t225 * t346 + t241 * t344 - t304) * t226 + (t194 * t343 + t205 * t341 + t306) * t195 + (MDP(18) * t257 + t331 * t347 + t332 * t348 + MDP(11) - 0.2e1 * t294) * t264 ^ 2 + (-t159 * t342 + t160 * t341 - t182 * t345 + t183 * t344 - t215 * t348 + t216 * t347 + t269 * t293 + t278 * t355 + 0.2e1 * t192 - 0.2e1 * t193 + 0.2e1 * t218 + 0.2e1 * t219 + 0.2e1 * t327) * t269; t301 - t259 * t316 + t190 * MDP(16) - t191 * MDP(17) + t263 * t305 + (-t206 * t263 + t207 * t268) * MDP(19) + (-pkin(3) * t206 - t322) * MDP(23) + (-pkin(3) * t207 + t323) * MDP(24) + t239 * t307 + (-t184 * t239 - t185 * t238) * MDP(26) + (t178 * t238 + t184 * t252) * MDP(30) + (t178 * t239 + t185 * t252) * MDP(31) + t202 * t312 + (-t167 * t202 - t168 * t201) * MDP(33) + (t163 * t201 + t167 * t217) * MDP(37) + (t163 * t202 + t168 * t217) * MDP(38) + t357 * t230; -t239 * t304 + (-t225 * t239 + t226 * t238) * MDP(26) + (t225 * t252 + t238 * t241) * MDP(30) + (-t226 * t252 + t239 * t241) * MDP(31) + t202 * t306 + (-t194 * t202 - t195 * t201) * MDP(33) + (t194 * t217 + t201 * t205) * MDP(37) + (t195 * t217 + t202 * t205) * MDP(38) + (t263 * t299 + (-t255 + t257) * MDP(19) + (-pkin(3) * t263 - t331) * MDP(23) + (-pkin(3) * t268 + t332) * MDP(24) - t292) * t264 + (-t326 - t357) * t269; 0.2e1 * t294 + t238 * t252 * t345 + t201 * t217 * t342 + MDP(18) * t255 + MDP(15) + 0.2e1 * (MDP(23) * t268 - MDP(24) * t263) * pkin(3) + (MDP(25) * t239 + t238 * t346 + t252 * t344) * t239 + (MDP(32) * t202 + t201 * t343 + t217 * t341) * t202; t230 * MDP(22) + (t230 * t335 + t157) * MDP(30) + (-t324 + (-t162 - t337) * t262) * MDP(31) + (t227 * t230 + t152) * MDP(37) + (-t228 * t230 - t153) * MDP(38) + t279 + t290 - t351; t289 * t264 + (-MDP(22) - MDP(29) - t277 - t354) * t269 + t360; t358; t293 - 0.2e1 * t302 + 0.2e1 * t303 + 0.2e1 * t354; (t230 * t333 + t152) * MDP(37) + (-t325 + (-t155 - t334) * t261) * MDP(38) + t279 + t285; (-t298 - t275) * t269 + t361; t274; (t245 + t333) * MDP(37) + (-pkin(5) - t251) * t315 + (t300 + (-MDP(37) * t261 - MDP(38) * t266 - MDP(31)) * t262) * pkin(4) + t298; 0.2e1 * t275 + t298; -t282 + t295; -t269 * MDP(36) + t359; t291; t277; MDP(36) + t275; MDP(36);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
