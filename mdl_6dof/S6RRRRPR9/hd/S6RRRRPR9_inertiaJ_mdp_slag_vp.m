% Calculate joint inertia matrix for
% S6RRRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR9_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:55:49
% EndTime: 2019-03-09 22:55:56
% DurationCPUTime: 2.04s
% Computational Cost: add. (3105->325), mult. (6819->460), div. (0->0), fcn. (7784->12), ass. (0->144)
t284 = sin(qJ(4));
t266 = pkin(3) * t284 + qJ(5);
t279 = sin(pkin(12));
t281 = cos(pkin(12));
t329 = t279 ^ 2 + t281 ^ 2;
t331 = t329 * t266;
t283 = sin(qJ(6));
t287 = cos(qJ(6));
t248 = t279 * t283 - t287 * t281;
t249 = t279 * t287 + t281 * t283;
t333 = t249 * MDP(31) - t248 * MDP(32);
t282 = cos(pkin(6));
t280 = sin(pkin(6));
t289 = cos(qJ(2));
t336 = t280 * t289;
t314 = pkin(8) * t336;
t286 = sin(qJ(2));
t347 = pkin(1) * t286;
t233 = t314 + (pkin(9) + t347) * t282;
t234 = (-pkin(2) * t289 - pkin(9) * t286 - pkin(1)) * t280;
t285 = sin(qJ(3));
t288 = cos(qJ(3));
t205 = -t233 * t285 + t288 * t234;
t337 = t280 * t286;
t236 = t285 * t282 + t288 * t337;
t315 = pkin(3) * t336;
t192 = -pkin(10) * t236 + t205 - t315;
t206 = t288 * t233 + t285 * t234;
t261 = t285 * t337;
t309 = t282 * t288 - t261;
t196 = t309 * pkin(10) + t206;
t348 = cos(qJ(4));
t312 = t348 * t196;
t181 = t284 * t192 + t312;
t177 = -qJ(5) * t336 + t181;
t210 = t236 * t284 - t348 * t309;
t211 = t348 * t236 + t284 * t309;
t262 = pkin(8) * t337;
t270 = -pkin(3) * t288 - pkin(2);
t346 = pkin(1) * t289;
t218 = t261 * pkin(3) + t262 + (t270 - t346) * t282;
t185 = t210 * pkin(4) - t211 * qJ(5) + t218;
t167 = -t177 * t279 + t281 * t185;
t168 = t281 * t177 + t279 * t185;
t361 = -t167 * t279 + t168 * t281;
t360 = t211 * MDP(20) - t210 * MDP(21);
t359 = t248 * MDP(34) + t249 * MDP(35);
t358 = MDP(25) * t281 - t279 * MDP(26);
t357 = -(MDP(16) * t285 + MDP(17) * t288) * pkin(9) + t285 * MDP(13) + t288 * MDP(14);
t349 = -pkin(10) - pkin(9);
t259 = t349 * t285;
t260 = t349 * t288;
t229 = -t348 * t259 - t260 * t284;
t230 = t284 * t259 - t348 * t260;
t250 = t284 * t285 - t348 * t288;
t251 = t284 * t288 + t348 * t285;
t356 = t251 * MDP(20) - t250 * MDP(21) - t229 * MDP(23) - t230 * MDP(24);
t355 = 2 * MDP(12);
t354 = 0.2e1 * MDP(16);
t353 = -2 * MDP(19);
t352 = 0.2e1 * MDP(23);
t351 = 2 * MDP(27);
t350 = -2 * MDP(30);
t345 = t281 * pkin(5);
t274 = t281 * pkin(11);
t343 = pkin(4) * MDP(28);
t334 = -t348 * t192 + t284 * t196;
t178 = pkin(4) * t336 + t334;
t341 = t178 * t281;
t340 = t229 * t281;
t339 = t251 * t279;
t338 = t266 * t281;
t335 = t282 * t286;
t219 = pkin(4) * t250 - qJ(5) * t251 + t270;
t194 = t279 * t219 + t281 * t230;
t330 = t329 * qJ(5);
t313 = t348 * pkin(3);
t269 = -t313 - pkin(4);
t327 = MDP(28) * t269;
t326 = t178 * MDP(28);
t203 = t211 * t279 + t281 * t336;
t204 = t211 * t281 - t279 * t336;
t186 = t287 * t203 + t204 * t283;
t325 = t186 * MDP(32);
t187 = -t203 * t283 + t204 * t287;
t324 = t187 * MDP(29);
t323 = t187 * MDP(31);
t322 = t210 * MDP(33);
t212 = t249 * t251;
t321 = t212 * MDP(32);
t213 = t248 * t251;
t320 = t213 * MDP(29);
t319 = t213 * MDP(31);
t318 = t250 * MDP(33);
t317 = t270 * MDP(24);
t316 = t285 * MDP(11);
t311 = MDP(22) + (MDP(29) * t249 + t248 * t350) * t249;
t193 = t281 * t219 - t230 * t279;
t308 = MDP(28) * t329;
t307 = -pkin(4) * t251 - qJ(5) * t250;
t305 = -t193 * t279 + t194 * t281;
t304 = -t203 * t281 + t204 * t279;
t303 = -t250 * t266 + t251 * t269;
t300 = 0.2e1 * t358;
t299 = t279 * MDP(25) + t281 * MDP(26);
t189 = pkin(5) * t250 - t251 * t274 + t193;
t190 = -pkin(11) * t339 + t194;
t173 = t189 * t287 - t190 * t283;
t174 = t189 * t283 + t190 * t287;
t298 = t173 * MDP(34) - t174 * MDP(35);
t297 = t212 * MDP(34) - t213 * MDP(35);
t296 = -t358 + t359;
t295 = 0.2e1 * t359;
t294 = (t348 * MDP(23) - t284 * MDP(24)) * pkin(3);
t293 = -t236 * MDP(13) - t309 * MDP(14);
t292 = t305 * MDP(27) + (-t212 * t249 + t213 * t248) * MDP(30) - t249 * t320 + t356 + t333 * t250;
t291 = -MDP(22) * t336 + (-t186 * t249 - t187 * t248) * MDP(30) + t249 * t324 + t360 + t333 * t210;
t276 = t280 ^ 2;
t267 = -pkin(4) - t345;
t255 = qJ(5) * t281 + t274;
t254 = (-pkin(11) - qJ(5)) * t279;
t253 = t269 - t345;
t245 = t274 + t338;
t244 = (-pkin(11) - t266) * t279;
t243 = pkin(1) * t335 + t314;
t242 = t282 * t346 - t262;
t232 = t262 + (-pkin(2) - t346) * t282;
t226 = t254 * t283 + t255 * t287;
t225 = t254 * t287 - t255 * t283;
t222 = t229 * t279;
t217 = t244 * t283 + t245 * t287;
t216 = t244 * t287 - t245 * t283;
t207 = pkin(5) * t339 + t229;
t198 = t207 * t249;
t197 = t207 * t248;
t176 = t178 * t279;
t172 = pkin(5) * t203 + t178;
t171 = t172 * t249;
t170 = t172 * t248;
t165 = -pkin(11) * t203 + t168;
t164 = pkin(5) * t210 - pkin(11) * t204 + t167;
t163 = t164 * t283 + t165 * t287;
t162 = t164 * t287 - t165 * t283;
t1 = [MDP(1) + (t167 ^ 2 + t168 ^ 2 + t178 ^ 2) * MDP(28) + t211 ^ 2 * MDP(18) + t282 ^ 2 * MDP(8) + (t236 * MDP(11) + t309 * t355) * t236 + (t186 * t350 + t324) * t187 + ((MDP(4) * t286 + 0.2e1 * MDP(5) * t289) * t286 + (MDP(15) + MDP(22)) * t289 ^ 2) * t276 + (t211 * t353 + t322 + 0.2e1 * t323 - 0.2e1 * t325) * t210 + 0.2e1 * (MDP(6) * t335 + (t282 * MDP(7) + t293 - t360) * t289) * t280 + 0.2e1 * (t242 * t282 + t276 * t346) * MDP(9) + (-t205 * t336 - t232 * t309) * t354 + 0.2e1 * (t181 * t336 + t211 * t218) * MDP(24) + 0.2e1 * (t206 * t336 + t232 * t236) * MDP(17) + (t210 * t218 + t334 * t336) * t352 + 0.2e1 * (-t243 * t282 - t276 * t347) * MDP(10) + (-t167 * t204 - t168 * t203) * t351 + 0.2e1 * (t167 * t210 + t178 * t203) * MDP(25) + 0.2e1 * (t162 * t210 + t172 * t186) * MDP(34) + 0.2e1 * (-t168 * t210 + t178 * t204) * MDP(26) + 0.2e1 * (-t163 * t210 + t172 * t187) * MDP(35); t236 * t316 + t210 * t318 - t187 * t320 + (t213 * t186 - t187 * t212) * MDP(30) + (t167 * t193 + t168 * t194 + t178 * t229) * MDP(28) + t242 * MDP(9) - t243 * MDP(10) + (t162 * t250 + t172 * t212 + t173 * t210 + t207 * t186) * MDP(34) + (-t186 * t250 - t212 * t210) * MDP(32) + (t187 * t250 - t213 * t210) * MDP(31) + (-t163 * t250 - t172 * t213 - t174 * t210 + t207 * t187) * MDP(35) + (t167 * t250 + t193 * t210 + t229 * t203) * MDP(25) + (-t168 * t250 - t194 * t210 + t229 * t204) * MDP(26) + (-t193 * t204 - t194 * t203) * MDP(27) + t282 * MDP(8) + (t236 * t288 + t285 * t309) * MDP(12) + (t270 * t210 + t218 * t250) * MDP(23) + (-pkin(2) * t236 + t232 * t285) * MDP(17) + (pkin(2) * t309 - t232 * t288) * MDP(16) + (-t250 * MDP(19) + t317) * t211 + (t211 * MDP(18) - t210 * MDP(19) + (-t167 * t281 - t168 * t279) * MDP(27) + t218 * MDP(24) + t299 * t178) * t251 + (MDP(6) * t286 + (MDP(7) - t356 - t357) * t289) * t280; MDP(8) + pkin(2) * t288 * t354 + (t193 ^ 2 + t194 ^ 2 + t229 ^ 2) * MDP(28) + (MDP(18) * t251 + 0.2e1 * t317) * t251 - (t212 * t350 - t320) * t213 + (-0.2e1 * pkin(2) * MDP(17) + t288 * t355 + t316) * t285 + (t251 * t353 + t270 * t352 + t318 - 0.2e1 * t319 - 0.2e1 * t321) * t250 + 0.2e1 * t297 * t207 + 0.2e1 * (t193 * MDP(25) - t194 * MDP(26) + t298) * t250 + 0.2e1 * ((-t193 * t281 - t194 * t279) * MDP(27) + t299 * t229) * t251; t291 + (t304 * t266 + t361) * MDP(27) + (t178 * t269 + t266 * t361) * MDP(28) - MDP(15) * t336 + t205 * MDP(16) - t206 * MDP(17) + (t186 * t253 + t210 * t216 + t170) * MDP(34) + (t187 * t253 - t210 * t217 + t171) * MDP(35) + (-t210 * t266 * t279 + t203 * t269 - t341) * MDP(25) + (t204 * t269 - t210 * t338 + t176) * MDP(26) + (-t313 * t336 - t334) * MDP(23) + (-t312 + (-t192 + t315) * t284) * MDP(24) - t293; t292 + (t303 * t279 - t340) * MDP(25) + (t303 * t281 + t222) * MDP(26) + (t229 * t269 + t305 * t266) * MDP(28) + (t212 * t253 + t216 * t250 + t197) * MDP(34) + (-t213 * t253 - t217 * t250 + t198) * MDP(35) + t357; MDP(15) + t331 * t351 + t266 ^ 2 * t308 + (-t300 + t327) * t269 + t253 * t295 + 0.2e1 * t294 + t311; -t334 * MDP(23) - t181 * MDP(24) + (-pkin(4) * t203 - t341) * MDP(25) + (-pkin(4) * t204 + t176) * MDP(26) + t361 * MDP(27) - pkin(4) * t326 + (t186 * t267 + t210 * t225 + t170) * MDP(34) + (t187 * t267 - t210 * t226 + t171) * MDP(35) + (t304 * MDP(27) + MDP(28) * t361 - t299 * t210) * qJ(5) + t291; (t307 * t279 - t340) * MDP(25) + (t307 * t281 + t222) * MDP(26) + (-pkin(4) * t229 + t305 * qJ(5)) * MDP(28) + (t212 * t267 + t225 * t250 + t197) * MDP(34) + (-t213 * t267 - t226 * t250 + t198) * MDP(35) + t292; (t330 + t331) * MDP(27) + (-pkin(4) * t269 + qJ(5) * t331) * MDP(28) + t294 + t311 + t358 * (pkin(4) - t269) + t359 * (t253 + t267); t330 * t351 + qJ(5) ^ 2 * t308 + t267 * t295 + (t300 + t343) * pkin(4) + t311; t203 * MDP(25) + t204 * MDP(26) + t186 * MDP(34) + t187 * MDP(35) + t326; t229 * MDP(28) + t299 * t251 + t297; t296 + t327; t296 - t343; MDP(28); t162 * MDP(34) - t163 * MDP(35) + t322 + t323 - t325; t298 + t318 - t319 - t321; MDP(34) * t216 - MDP(35) * t217 + t333; MDP(34) * t225 - t226 * MDP(35) + t333; 0; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
