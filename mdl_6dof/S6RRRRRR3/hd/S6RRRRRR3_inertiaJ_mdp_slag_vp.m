% Calculate joint inertia matrix for
% S6RRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR3_inertiaJ_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:42:43
% EndTime: 2019-03-10 03:42:47
% DurationCPUTime: 1.28s
% Computational Cost: add. (1855->229), mult. (3514->300), div. (0->0), fcn. (4101->10), ass. (0->129)
t260 = sin(qJ(5));
t261 = sin(qJ(4));
t265 = cos(qJ(5));
t266 = cos(qJ(4));
t233 = t260 * t261 - t265 * t266;
t235 = t260 * t266 + t265 * t261;
t259 = sin(qJ(6));
t264 = cos(qJ(6));
t197 = t264 * t233 + t259 * t235;
t198 = -t259 * t233 + t264 * t235;
t301 = t198 * MDP(34) - t197 * MDP(35);
t281 = t235 * MDP(27) - t233 * MDP(28) + t301;
t325 = t261 * MDP(20) + t266 * MDP(21);
t328 = pkin(7) + pkin(8);
t262 = sin(qJ(3));
t248 = t262 * pkin(2) + pkin(9);
t229 = (-pkin(10) - t248) * t261;
t256 = t266 * pkin(10);
t230 = t266 * t248 + t256;
t190 = t265 * t229 - t260 * t230;
t313 = t235 * pkin(11);
t175 = t190 - t313;
t191 = t260 * t229 + t265 * t230;
t228 = t233 * pkin(11);
t176 = t191 - t228;
t154 = t264 * t175 - t259 * t176;
t155 = t259 * t175 + t264 * t176;
t327 = t154 * MDP(37) - t155 * MDP(38);
t240 = (-pkin(9) - pkin(10)) * t261;
t242 = t266 * pkin(9) + t256;
t207 = t265 * t240 - t260 * t242;
t183 = t207 - t313;
t209 = t260 * t240 + t265 * t242;
t184 = t209 - t228;
t161 = t264 * t183 - t259 * t184;
t162 = t259 * t183 + t264 * t184;
t326 = t161 * MDP(37) - t162 * MDP(38);
t277 = t266 * MDP(23) - t261 * MDP(24);
t324 = t233 * MDP(30) + t235 * MDP(31);
t323 = t197 * MDP(37) + t198 * MDP(38);
t319 = cos(qJ(2));
t252 = -t319 * pkin(2) - pkin(1);
t322 = 0.2e1 * t252;
t321 = -2 * MDP(26);
t320 = -2 * MDP(33);
t318 = cos(qJ(3));
t317 = pkin(4) * t260;
t316 = t233 * pkin(5);
t263 = sin(qJ(2));
t234 = t262 * t263 - t318 * t319;
t315 = t234 * pkin(4);
t314 = t234 * pkin(5);
t312 = t264 * pkin(5);
t311 = t265 * pkin(4);
t310 = t266 * pkin(4);
t241 = t328 * t263;
t243 = t328 * t319;
t208 = t318 * t241 + t262 * t243;
t308 = t208 * t266;
t236 = t262 * t319 + t318 * t263;
t307 = t236 * t261;
t306 = t236 * t266;
t305 = t261 * t266;
t196 = t234 * pkin(3) - t236 * pkin(9) + t252;
t210 = -t262 * t241 + t318 * t243;
t170 = t266 * t196 - t261 * t210;
t158 = -pkin(10) * t306 + t170 + t315;
t302 = t266 * t210;
t165 = t302 + (-pkin(10) * t236 + t196) * t261;
t303 = t265 * t165;
t148 = t260 * t158 + t303;
t185 = t235 * t236;
t146 = -t185 * pkin(11) + t148;
t304 = t264 * t146;
t186 = t233 * t236;
t298 = MDP(25) * t186;
t167 = -t259 * t185 - t264 * t186;
t297 = MDP(32) * t167;
t166 = t264 * t185 - t259 * t186;
t163 = t166 * MDP(35);
t164 = t167 * MDP(34);
t181 = t185 * MDP(28);
t182 = t186 * MDP(27);
t249 = pkin(5) + t311;
t244 = t264 * t249;
t218 = -t259 * t317 + t244;
t294 = t218 * MDP(37);
t219 = t259 * t249 + t264 * t317;
t293 = t219 * MDP(38);
t290 = t259 * MDP(38);
t288 = t265 * MDP(30);
t286 = 0.2e1 * t319;
t285 = MDP(29) + MDP(36);
t284 = t234 * MDP(36) - t163 + t164;
t251 = -pkin(3) - t310;
t283 = MDP(19) * t305;
t282 = MDP(22) + t285;
t147 = t265 * t158 - t260 * t165;
t144 = t186 * pkin(11) + t147 + t314;
t141 = t264 * t144 - t259 * t146;
t250 = -t318 * pkin(2) - pkin(3);
t280 = -pkin(3) * t236 - pkin(9) * t234;
t142 = t259 * t144 + t304;
t279 = -t234 * t248 + t236 * t250;
t278 = t266 * MDP(20) - t261 * MDP(21);
t276 = -MDP(23) * t261 - MDP(24) * t266;
t179 = pkin(4) * t307 + t208;
t275 = t234 * MDP(29) - t181 - t182 + t284;
t239 = t250 - t310;
t274 = (t264 * MDP(37) - t290) * pkin(5);
t257 = t261 ^ 2;
t273 = t257 * MDP(18) + MDP(15) + 0.2e1 * t283 + (MDP(25) * t235 + t233 * t321) * t235 + (MDP(32) * t198 + t197 * t320) * t198;
t272 = 0.2e1 * t324;
t271 = 0.2e1 * t323;
t270 = (t318 * MDP(16) - t262 * MDP(17)) * pkin(2);
t269 = t190 * MDP(30) - t191 * MDP(31) + t281 + t327;
t268 = t207 * MDP(30) - t209 * MDP(31) + t281 + t326;
t258 = t266 ^ 2;
t267 = (-t198 * t166 - t167 * t197) * MDP(33) + t198 * t297 + (-t235 * t185 + t186 * t233) * MDP(26) - t235 * t298 - t208 * MDP(16) - t210 * MDP(17) + ((-t257 + t258) * MDP(19) + MDP(18) * t305 + MDP(13)) * t236 + (-MDP(14) + t325 + t281) * t234;
t212 = t251 + t316;
t211 = t239 + t316;
t199 = t208 * t261;
t173 = t179 * t235;
t172 = t179 * t233;
t171 = t261 * t196 + t302;
t168 = t185 * pkin(5) + t179;
t151 = t168 * t198;
t150 = t168 * t197;
t1 = [pkin(1) * MDP(9) * t286 + MDP(1) - (t185 * t321 - t298) * t186 + (t166 * t320 + t297) * t167 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t263 + MDP(5) * t286) * t263 + t282 * t234 ^ 2 + (MDP(16) * t322 - 0.2e1 * t163 + 0.2e1 * t164 - 0.2e1 * t181 - 0.2e1 * t182) * t234 + 0.2e1 * (-t171 * t234 + t208 * t306) * MDP(24) + 0.2e1 * (t170 * t234 + t208 * t307) * MDP(23) + 0.2e1 * (-t148 * t234 - t179 * t186) * MDP(31) + 0.2e1 * (-t142 * t234 + t168 * t167) * MDP(38) + 0.2e1 * (t147 * t234 + t179 * t185) * MDP(30) + 0.2e1 * (t141 * t234 + t168 * t166) * MDP(37) + (MDP(17) * t322 + 0.2e1 * (-MDP(12) + t278) * t234 + (t258 * MDP(18) + MDP(11) - 0.2e1 * t283) * t236) * t236; (-t319 * MDP(10) - t263 * MDP(9)) * pkin(7) + (t279 * t266 + t199) * MDP(24) + t267 + (t279 * t261 - t308) * MDP(23) + (t154 * t234 + t211 * t166 + t150) * MDP(37) + (-t155 * t234 + t211 * t167 + t151) * MDP(38) + (t239 * t185 + t190 * t234 + t172) * MDP(30) + (-t239 * t186 - t191 * t234 + t173) * MDP(31) + t263 * MDP(6) + t319 * MDP(7); t211 * t271 + t239 * t272 - 0.2e1 * t250 * t277 + MDP(8) + 0.2e1 * t270 + t273; (t280 * t261 - t308) * MDP(23) + t267 + (t280 * t266 + t199) * MDP(24) + (t161 * t234 + t212 * t166 + t150) * MDP(37) + (-t162 * t234 + t212 * t167 + t151) * MDP(38) + (t251 * t185 + t207 * t234 + t172) * MDP(30) + (-t251 * t186 - t209 * t234 + t173) * MDP(31); t270 + t273 + t277 * (pkin(3) - t250) + t323 * (t211 + t212) + t324 * (t239 + t251); 0.2e1 * pkin(3) * t277 + t212 * t271 + t251 * t272 + t273; t234 * MDP(22) + t170 * MDP(23) - t171 * MDP(24) + (t234 * t311 + t147) * MDP(30) + (-t303 + (-t158 - t315) * t260) * MDP(31) + (t218 * t234 + t141) * MDP(37) + (-t219 * t234 - t142) * MDP(38) + t278 * t236 + t275; t248 * t276 + t269 + t325; pkin(9) * t276 + t268 + t325; 0.2e1 * (-t260 * MDP(31) + t288) * pkin(4) + 0.2e1 * t294 - 0.2e1 * t293 + t282; t147 * MDP(30) - t148 * MDP(31) + (t234 * t312 + t141) * MDP(37) + (-t304 + (-t144 - t314) * t259) * MDP(38) + t275; t269; t268; (t244 + t312) * MDP(37) + (-pkin(5) - t249) * t290 + (t288 + (-MDP(37) * t259 - MDP(38) * t264 - MDP(31)) * t260) * pkin(4) + t285; 0.2e1 * t274 + t285; t141 * MDP(37) - t142 * MDP(38) + t284; t301 + t327; t301 + t326; MDP(36) - t293 + t294; MDP(36) + t274; MDP(36);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
