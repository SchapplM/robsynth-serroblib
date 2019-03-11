% Calculate joint inertia matrix for
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP12_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP12_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP12_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:20:45
% EndTime: 2019-03-10 03:20:53
% DurationCPUTime: 2.58s
% Computational Cost: add. (5069->430), mult. (12985->611), div. (0->0), fcn. (14622->12), ass. (0->168)
t250 = sin(pkin(7));
t358 = 0.2e1 * t250;
t255 = sin(qJ(4));
t357 = -0.2e1 * t255;
t253 = cos(pkin(6));
t257 = sin(qJ(2));
t251 = sin(pkin(6));
t261 = cos(qJ(2));
t323 = t251 * t261;
t230 = t253 * t257 * pkin(1) + pkin(9) * t323;
t252 = cos(pkin(7));
t322 = t252 * t261;
t290 = t251 * t322;
t206 = (t250 * t253 + t290) * pkin(10) + t230;
t256 = sin(qJ(3));
t260 = cos(qJ(3));
t346 = pkin(1) * t261;
t241 = t253 * t346;
t324 = t251 * t257;
t209 = t253 * pkin(2) + t241 + (-pkin(10) * t252 - pkin(9)) * t324;
t216 = (-pkin(10) * t250 * t257 - pkin(2) * t261 - pkin(1)) * t251;
t279 = t209 * t252 + t216 * t250;
t189 = -t256 * t206 + t279 * t260;
t356 = 2 * MDP(16);
t355 = 2 * MDP(17);
t354 = 2 * MDP(23);
t353 = 2 * MDP(24);
t352 = -2 * MDP(26);
t351 = 2 * MDP(30);
t350 = 2 * MDP(31);
t349 = 2 * MDP(32);
t348 = 2 * MDP(33);
t347 = 2 * MDP(34);
t345 = pkin(2) * t256;
t344 = pkin(2) * t260;
t326 = t250 * t256;
t208 = t253 * t326 + (t256 * t322 + t257 * t260) * t251;
t224 = -t250 * t323 + t252 * t253;
t259 = cos(qJ(4));
t197 = t208 * t255 - t224 * t259;
t343 = pkin(5) * t197;
t225 = -t259 * t252 + t255 * t326;
t342 = pkin(5) * t225;
t254 = sin(qJ(5));
t341 = pkin(11) * t254;
t258 = cos(qJ(5));
t340 = pkin(11) * t258;
t339 = pkin(3) * MDP(23);
t338 = pkin(3) * MDP(24);
t337 = qJ(6) * t197;
t336 = qJ(6) * t225;
t193 = -t209 * t250 + t252 * t216;
t325 = t250 * t260;
t207 = -t253 * t325 + t256 * t324 - t260 * t290;
t182 = pkin(3) * t207 - pkin(11) * t208 + t193;
t190 = t260 * t206 + t279 * t256;
t185 = pkin(11) * t224 + t190;
t175 = t182 * t259 - t255 * t185;
t173 = -pkin(4) * t207 - t175;
t335 = t173 * t254;
t334 = t173 * t258;
t198 = t208 * t259 + t224 * t255;
t191 = t198 * t254 - t207 * t258;
t333 = t191 * t258;
t192 = t198 * t258 + t207 * t254;
t332 = t192 * t254;
t291 = pkin(10) * t325;
t221 = t291 + (pkin(11) + t345) * t252;
t222 = (-pkin(3) * t260 - pkin(11) * t256 - pkin(2)) * t250;
t202 = -t255 * t221 + t222 * t259;
t200 = pkin(4) * t325 - t202;
t331 = t200 * t254;
t330 = t200 * t258;
t226 = t252 * t255 + t259 * t326;
t210 = t226 * t254 + t258 * t325;
t329 = t210 * t258;
t211 = t226 * t258 - t254 * t325;
t328 = t211 * t254;
t246 = t251 ^ 2;
t327 = t246 * t257;
t321 = t253 * MDP(8);
t320 = t254 * t259;
t176 = t182 * t255 + t185 * t259;
t174 = pkin(12) * t207 + t176;
t184 = -t224 * pkin(3) - t189;
t179 = pkin(4) * t197 - pkin(12) * t198 + t184;
t170 = t258 * t174 + t254 * t179;
t239 = pkin(10) * t326;
t220 = t239 + (-pkin(3) - t344) * t252;
t199 = pkin(4) * t225 - pkin(12) * t226 + t220;
t203 = t221 * t259 + t222 * t255;
t201 = -pkin(12) * t325 + t203;
t187 = t254 * t199 + t258 * t201;
t235 = -pkin(4) * t259 - pkin(12) * t255 - pkin(3);
t218 = t254 * t235 + t259 * t340;
t247 = t254 ^ 2;
t249 = t258 ^ 2;
t319 = t247 + t249;
t318 = MDP(11) * t208;
t317 = MDP(12) * t208;
t316 = MDP(13) * t224;
t315 = MDP(13) * t256;
t314 = MDP(14) * t224;
t313 = MDP(14) * t252;
t312 = MDP(22) * t260;
t311 = MDP(25) * t211;
t310 = MDP(25) * t254;
t309 = MDP(25) * t258;
t308 = MDP(26) * t258;
t283 = -pkin(5) * t258 - qJ(6) * t254;
t234 = -pkin(4) + t283;
t307 = MDP(35) * t234;
t306 = t191 * MDP(28);
t305 = t192 * MDP(27);
t304 = t197 * MDP(29);
t303 = t198 * MDP(19);
t302 = t198 * MDP(20);
t301 = t207 * MDP(22);
t300 = t210 * MDP(28);
t299 = t211 * MDP(27);
t298 = t225 * MDP(29);
t297 = t226 * MDP(18);
t296 = t226 * MDP(19);
t295 = t226 * MDP(20);
t294 = t252 * MDP(15);
t293 = t254 * MDP(27);
t292 = MDP(31) - MDP(34);
t289 = pkin(11) * MDP(23) - MDP(20);
t288 = pkin(11) * MDP(24) - MDP(21);
t287 = -MDP(35) * pkin(5) - MDP(32);
t286 = t174 * t254 - t258 * t179;
t285 = -t258 * t199 + t201 * t254;
t284 = -0.2e1 * qJ(6) * MDP(34) - MDP(29);
t282 = -pkin(5) * t254 + qJ(6) * t258;
t167 = t170 + t337;
t168 = t286 - t343;
t281 = t167 * t258 + t168 * t254;
t180 = t187 + t336;
t181 = t285 - t342;
t280 = t180 * t258 + t181 * t254;
t213 = -qJ(6) * t259 + t218;
t232 = t258 * t235;
t214 = -t232 + (pkin(5) + t341) * t259;
t278 = t213 * t258 + t214 * t254;
t227 = t252 * t344 - t239;
t229 = t252 * t345 + t291;
t277 = t227 * MDP(16) - t229 * MDP(17);
t276 = t258 * MDP(27) - t254 * MDP(28);
t275 = t258 * MDP(28) + t293;
t274 = (-pkin(11) * t320 + t232) * MDP(30) - t218 * MDP(31);
t273 = t254 * MDP(32) - t258 * MDP(34);
t272 = -MDP(19) + t276;
t271 = (MDP(6) * t257 + MDP(7) * t261) * t251;
t270 = t202 * MDP(23) - t203 * MDP(24) + t295;
t269 = t208 * MDP(13) + t224 * MDP(15) + t189 * MDP(16) - t190 * MDP(17);
t268 = t175 * MDP(23) - t176 * MDP(24) + t301 + t302;
t267 = MDP(32) * t214 - MDP(34) * t213 - t274;
t266 = -t267 - t339;
t265 = -MDP(30) * t286 - t170 * MDP(31) + t304 + t305 - t306;
t264 = -MDP(30) * t285 - t187 * MDP(31) + t298 + t299 - t300;
t263 = -MDP(21) + (-t292 * t258 + (-MDP(30) - MDP(32)) * t254) * pkin(12) + t275;
t245 = t250 ^ 2;
t242 = pkin(12) * t320;
t228 = -pkin(9) * t324 + t241;
t223 = (pkin(11) - t282) * t255;
t188 = pkin(5) * t210 - qJ(6) * t211 + t200;
t171 = pkin(5) * t191 - qJ(6) * t192 + t173;
t1 = [t224 ^ 2 * MDP(15) + t198 ^ 2 * MDP(18) + MDP(1) + (t167 ^ 2 + t168 ^ 2 + t171 ^ 2) * MDP(35) + (MDP(4) * t257 + 0.2e1 * MDP(5) * t261) * t327 + (0.2e1 * t316 + t318) * t208 + (MDP(25) * t192 + t191 * t352) * t192 + (0.2e1 * t271 + t321) * t253 + (t301 + 0.2e1 * t302 - 0.2e1 * t314 - 0.2e1 * t317) * t207 + (-0.2e1 * MDP(21) * t207 - 0.2e1 * t303 + t304 + 0.2e1 * t305 - 0.2e1 * t306) * t197 + 0.2e1 * (-pkin(1) * t327 - t230 * t253) * MDP(10) + (t189 * t224 + t193 * t207) * t356 + (-t190 * t224 + t193 * t208) * t355 + (t175 * t207 + t184 * t197) * t354 + (-t176 * t207 + t184 * t198) * t353 + (t173 * t191 - t197 * t286) * t351 + (t167 * t197 - t171 * t192) * t347 + (-t168 * t197 + t171 * t191) * t349 + (-t170 * t197 + t173 * t192) * t350 + (-t167 * t191 + t168 * t192) * t348 + 0.2e1 * (t228 * t253 + t246 * t346) * MDP(9); t197 * t298 + t198 * t297 + t192 * t311 + (t184 * t226 + t198 * t220) * MDP(24) + (t184 * t225 + t197 * t220) * MDP(23) + t321 - t230 * MDP(10) + (-t197 * t226 - t198 * t225) * MDP(19) + t228 * MDP(9) + (t192 * t225 + t197 * t211) * MDP(27) + (t167 * t225 - t171 * t211 + t180 * t197 - t188 * t192) * MDP(34) + (t173 * t210 + t191 * t200 - t197 * t285 - t225 * t286) * MDP(30) + (-t191 * t225 - t197 * t210) * MDP(28) + (-t170 * t225 + t173 * t211 - t187 * t197 + t192 * t200) * MDP(31) + (-t168 * t225 + t171 * t210 - t181 * t197 + t188 * t191) * MDP(32) + (-t191 * t211 - t192 * t210) * MDP(26) + (-t167 * t210 + t168 * t211 - t180 * t191 + t181 * t192) * MDP(33) + (t167 * t180 + t168 * t181 + t171 * t188) * MDP(35) + t277 * t224 + t269 * t252 + t271 + (-t225 * MDP(21) + t270 - t313) * t207 + ((-MDP(16) * t207 - MDP(17) * t208) * pkin(2) + (-MDP(12) * t207 + MDP(17) * t193 + t316 + t318) * t256 + (-MDP(16) * t193 + t197 * MDP(21) - t268 + t314 + t317) * t260) * t250; t245 * t256 ^ 2 * MDP(11) + t226 ^ 2 * MDP(18) + (t180 ^ 2 + t181 ^ 2 + t188 ^ 2) * MDP(35) + MDP(8) + (t315 * t358 + t294) * t252 + (t210 * t352 + t311) * t211 + ((-t295 + t313) * t358 + (0.2e1 * MDP(12) * t256 + t312) * t245) * t260 + (0.2e1 * MDP(21) * t325 - 0.2e1 * t296 + t298 + 0.2e1 * t299 - 0.2e1 * t300) * t225 + (-t229 * t252 - t245 * t345) * t355 + (t180 * t225 - t188 * t211) * t347 + (t200 * t210 - t225 * t285) * t351 + (-t181 * t225 + t188 * t210) * t349 + (-t187 * t225 + t200 * t211) * t350 + (-t180 * t210 + t181 * t211) * t348 + (-t202 * t325 + t220 * t225) * t354 + (t203 * t325 + t220 * t226) * t353 + (t227 * t252 + t245 * t344) * t356; -t207 * MDP(14) - t198 * t338 + (-t191 * t213 + t192 * t214) * MDP(33) + (t167 * t213 + t168 * t214) * MDP(35) + (MDP(32) * t191 - MDP(34) * t192 + MDP(35) * t171) * t223 + t266 * t197 + (-t184 * MDP(23) + t168 * MDP(32) - t167 * MDP(34) - t288 * t207 - t265 + t303) * t259 + (t198 * MDP(18) + t184 * MDP(24) + t192 * t309 + (-t332 - t333) * MDP(26) + (pkin(11) * t191 + t335) * MDP(30) + (pkin(11) * t192 + t334) * MDP(31) + (-t167 * t254 + t168 * t258) * MDP(33) - t289 * t207 + t273 * t171 + t272 * t197) * t255 + t269; t294 - t226 * t338 + (-t210 * t213 + t211 * t214) * MDP(33) + (t180 * t213 + t181 * t214) * MDP(35) + (MDP(14) * t260 + t315) * t250 + (MDP(32) * t210 - MDP(34) * t211 + MDP(35) * t188) * t223 + t266 * t225 + (-t220 * MDP(23) + t181 * MDP(32) - t180 * MDP(34) + t288 * t325 - t264 + t296) * t259 + (t297 + t220 * MDP(24) + t211 * t309 + (-t328 - t329) * MDP(26) + (pkin(11) * t210 + t331) * MDP(30) + (pkin(11) * t211 + t330) * MDP(31) + (-t180 * t254 + t181 * t258) * MDP(33) + t289 * t325 + t273 * t188 + t272 * t225) * t255 + t277; MDP(15) + t338 * t357 + (t213 ^ 2 + t214 ^ 2 + t223 ^ 2) * MDP(35) + 0.2e1 * ((-t213 * t254 + t214 * t258) * MDP(33) + t273 * t223) * t255 + (t259 * MDP(29) + t272 * t357 + 0.2e1 * t267 + (2 * t339)) * t259 + (MDP(25) * t249 - 0.2e1 * t254 * t308 + MDP(18) + 0.2e1 * (t254 * MDP(30) + t258 * MDP(31)) * pkin(11)) * t255 ^ 2; t192 * t310 + (-t191 * t254 + t192 * t258) * MDP(26) + (-pkin(4) * t191 - t334) * MDP(30) + (-pkin(4) * t192 + t335) * MDP(31) + (-t171 * t258 + t191 * t234) * MDP(32) + t281 * MDP(33) + (-t171 * t254 - t192 * t234) * MDP(34) + t171 * t307 + ((t332 - t333) * MDP(33) + t281 * MDP(35)) * pkin(12) + t263 * t197 + t268; -t250 * t312 + t211 * t310 + (-t210 * t254 + t211 * t258) * MDP(26) + (-pkin(4) * t210 - t330) * MDP(30) + (-pkin(4) * t211 + t331) * MDP(31) + (-t188 * t258 + t210 * t234) * MDP(32) + t280 * MDP(33) + (-t188 * t254 - t211 * t234) * MDP(34) + t188 * t307 + ((t328 - t329) * MDP(33) + t280 * MDP(35)) * pkin(12) + t263 * t225 + t270; t242 * MDP(30) + (-t223 * t258 + t242) * MDP(32) + t278 * MDP(33) - t223 * t254 * MDP(34) + (pkin(12) * t278 + t223 * t234) * MDP(35) + (-t293 + (pkin(12) * t292 - MDP(28)) * t258 - t288) * t259 + (t254 * t309 + (-t247 + t249) * MDP(26) + (-pkin(4) * t254 - t340) * MDP(30) + (-pkin(4) * t258 + t341) * MDP(31) + t273 * t234 - t289) * t255; MDP(22) + t247 * MDP(25) + (t319 * pkin(12) ^ 2 + t234 ^ 2) * MDP(35) + t319 * pkin(12) * t348 + 0.2e1 * (MDP(30) * pkin(4) - MDP(32) * t234) * t258 + 0.2e1 * (-MDP(31) * pkin(4) - MDP(34) * t234 + t308) * t254; (-t286 + 0.2e1 * t343) * MDP(32) + (-pkin(5) * t192 - qJ(6) * t191) * MDP(33) + (t170 + 0.2e1 * t337) * MDP(34) + (-pkin(5) * t168 + qJ(6) * t167) * MDP(35) + t265; (-t285 + 0.2e1 * t342) * MDP(32) + (-pkin(5) * t211 - qJ(6) * t210) * MDP(33) + (t187 + 0.2e1 * t336) * MDP(34) + (-pkin(5) * t181 + qJ(6) * t180) * MDP(35) + t264; t232 * MDP(32) + t218 * MDP(34) + (-pkin(5) * t214 + qJ(6) * t213) * MDP(35) + ((-0.2e1 * pkin(5) - t341) * MDP(32) + t284) * t259 + (MDP(33) * t283 + t276) * t255 + t274; t282 * MDP(33) + ((MDP(35) * qJ(6) - t292) * t258 + (-MDP(30) + t287) * t254) * pkin(12) + t275; pkin(5) * t349 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(35) - t284; -MDP(32) * t197 + MDP(33) * t192 + MDP(35) * t168; -MDP(32) * t225 + MDP(33) * t211 + MDP(35) * t181; MDP(33) * t255 * t258 + MDP(32) * t259 + MDP(35) * t214; (MDP(35) * pkin(12) + MDP(33)) * t254; t287; MDP(35);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
