% Calculate joint inertia matrix for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR13_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR13_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR13_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:04:09
% EndTime: 2019-03-10 00:04:14
% DurationCPUTime: 1.69s
% Computational Cost: add. (1880->310), mult. (4152->439), div. (0->0), fcn. (4453->10), ass. (0->135)
t222 = sin(qJ(3));
t226 = cos(qJ(3));
t201 = -pkin(3) * t226 - pkin(10) * t222 - pkin(2);
t221 = sin(qJ(4));
t225 = cos(qJ(4));
t298 = pkin(9) * t225;
t179 = t221 * t201 + t226 * t298;
t212 = t226 * pkin(4);
t284 = t221 * t226;
t281 = pkin(9) * t284 - t225 * t201;
t175 = t212 + t281;
t283 = t222 * t225;
t166 = pkin(5) * t226 - pkin(11) * t283 + t175;
t174 = -qJ(5) * t226 + t179;
t167 = pkin(11) * t221 * t222 + t174;
t220 = sin(qJ(6));
t224 = cos(qJ(6));
t240 = -(t166 * t220 + t167 * t224) * MDP(35) + MDP(34) * (t166 * t224 - t167 * t220);
t313 = -MDP(23) * t281 - t179 * MDP(24) - t240;
t219 = cos(pkin(6));
t218 = sin(pkin(6));
t223 = sin(qJ(2));
t288 = t218 * t223;
t189 = t219 * t222 + t226 * t288;
t227 = cos(qJ(2));
t287 = t218 * t227;
t170 = t189 * t225 - t221 * t287;
t187 = -t219 * t226 + t222 * t288;
t205 = pkin(8) * t288;
t299 = pkin(1) * t227;
t180 = t205 + (-pkin(2) - t299) * t219;
t161 = pkin(3) * t187 - pkin(10) * t189 + t180;
t253 = pkin(8) * t287;
t300 = pkin(1) * t223;
t181 = t253 + (pkin(9) + t300) * t219;
t182 = (-pkin(2) * t227 - pkin(9) * t223 - pkin(1)) * t218;
t165 = t181 * t226 + t182 * t222;
t163 = -pkin(10) * t287 + t165;
t282 = -t225 * t161 + t221 * t163;
t302 = pkin(4) + pkin(5);
t145 = -pkin(11) * t170 - t187 * t302 + t282;
t151 = t221 * t161 + t225 * t163;
t186 = t187 * qJ(5);
t148 = t186 + t151;
t169 = t189 * t221 + t225 * t287;
t146 = pkin(11) * t169 + t148;
t143 = t145 * t224 - t146 * t220;
t144 = t145 * t220 + t146 * t224;
t156 = t169 * t220 + t170 * t224;
t267 = t156 * MDP(31);
t155 = -t169 * t224 + t170 * t220;
t268 = t155 * MDP(32);
t312 = t143 * MDP(34) - t144 * MDP(35) + t267 - t268;
t311 = -0.2e1 * t222;
t257 = MDP(22) + MDP(33);
t310 = t187 * t257;
t264 = t170 * MDP(20);
t266 = t169 * MDP(21);
t308 = -MDP(23) * t282 - t151 * MDP(24) + t264 - t266 - t312;
t307 = 2 * MDP(25);
t306 = 2 * MDP(26);
t305 = -2 * MDP(30);
t304 = 0.2e1 * MDP(34);
t303 = 0.2e1 * MDP(35);
t301 = pkin(10) - pkin(11);
t297 = pkin(2) * MDP(16);
t296 = pkin(2) * MDP(17);
t295 = pkin(4) * MDP(28);
t294 = qJ(5) * t221;
t293 = qJ(5) * t225;
t164 = -t222 * t181 + t226 * t182;
t162 = pkin(3) * t287 - t164;
t292 = t162 * t221;
t291 = t162 * t225;
t290 = t169 * t225;
t289 = t170 * t221;
t286 = t219 * MDP(8);
t285 = t221 * t224;
t214 = t221 ^ 2;
t216 = t225 ^ 2;
t280 = t214 + t216;
t279 = MDP(15) * t227;
t238 = qJ(5) * t170 - t162;
t152 = pkin(4) * t169 - t238;
t278 = MDP(28) * t152;
t195 = t220 * t221 + t224 * t225;
t185 = t195 * t222;
t277 = MDP(29) * t185;
t196 = -t220 * t225 + t285;
t276 = MDP(29) * t196;
t274 = MDP(34) * (qJ(5) * t220 + t224 * t302);
t272 = MDP(35) * (t224 * qJ(5) - t220 * t302);
t271 = qJ(5) * MDP(27);
t265 = t170 * MDP(18);
t184 = t220 * t283 - t222 * t285;
t263 = t184 * MDP(32);
t262 = t185 * MDP(31);
t261 = t189 * MDP(12);
t260 = t189 * MDP(13);
t259 = t225 * MDP(18);
t258 = t225 * MDP(19);
t256 = -MDP(23) - MDP(25);
t255 = MDP(24) - MDP(27);
t254 = pkin(4) * t307;
t252 = t301 * t221;
t251 = pkin(9) * MDP(16) - MDP(13);
t250 = pkin(9) * MDP(17) - MDP(14);
t249 = -pkin(4) * t225 - t294;
t248 = -pkin(4) * t221 + t293;
t149 = -pkin(4) * t187 + t282;
t247 = t148 * t225 + t149 * t221;
t246 = t174 * t225 + t175 * t221;
t245 = t225 * MDP(20) - t221 * MDP(21);
t243 = t221 * MDP(25) - t225 * MDP(27);
t241 = -t262 + t263;
t239 = MDP(34) * t224 - MDP(35) * t220;
t237 = -MDP(12) + t245;
t236 = MDP(25) + t239;
t235 = -MDP(33) - t272 - t274;
t234 = MDP(22) - t235;
t202 = t301 * t225;
t233 = t196 * MDP(31) - t195 * MDP(32) - (t202 * t220 - t224 * t252) * MDP(34) - (t224 * t202 + t220 * t252) * MDP(35);
t232 = -t221 * MDP(20) + t233;
t230 = t225 * MDP(21) - t232;
t229 = MDP(25) * t175 - MDP(27) * t174 - t313;
t213 = t218 ^ 2;
t208 = pkin(10) * t284;
t200 = -pkin(3) + t249;
t192 = t225 * t302 + pkin(3) + t294;
t191 = t219 * t300 + t253;
t190 = t219 * t299 - t205;
t183 = (pkin(9) - t248) * t222;
t171 = (-t221 * t302 - pkin(9) + t293) * t222;
t147 = -t169 * t302 + t238;
t1 = [t213 * t223 ^ 2 * MDP(4) + (t148 ^ 2 + t149 ^ 2 + t152 ^ 2) * MDP(28) + t189 ^ 2 * MDP(11) + MDP(1) + (0.2e1 * MDP(6) * t288 + t286) * t219 + (-0.2e1 * t169 * MDP(19) + t265) * t170 + (MDP(29) * t156 + t155 * t305) * t156 + (0.2e1 * (MDP(7) * t219 - t260) * t218 + (0.2e1 * MDP(5) * t223 + t279) * t213) * t227 + (0.2e1 * MDP(14) * t287 - 0.2e1 * t261 + 0.2e1 * t264 - 0.2e1 * t266 - 0.2e1 * t267 + 0.2e1 * t268 + t310) * t187 + 0.2e1 * (-t191 * t219 - t213 * t300) * MDP(10) + 0.2e1 * (t148 * t187 - t152 * t170) * MDP(27) + 0.2e1 * (t162 * t169 - t187 * t282) * MDP(23) + (t144 * t187 + t147 * t156) * t303 + 0.2e1 * (-t151 * t187 + t162 * t170) * MDP(24) + (-t149 * t187 + t152 * t169) * t307 + (-t143 * t187 + t147 * t155) * t304 + (-t148 * t169 + t149 * t170) * t306 + 0.2e1 * (t165 * t287 + t180 * t189) * MDP(17) + 0.2e1 * (-t164 * t287 + t180 * t187) * MDP(16) + 0.2e1 * (t190 * t219 + t213 * t299) * MDP(9); t156 * t277 - t189 * t296 + (t148 * t174 + t149 * t175) * MDP(28) + (-t155 * t185 - t156 * t184) * MDP(30) + t190 * MDP(9) - t191 * MDP(10) + t286 + (-t169 * t174 + t170 * t175) * MDP(26) + (t147 * t185 + t156 * t171) * MDP(35) + (t147 * t184 + t155 * t171) * MDP(34) + (MDP(6) * t223 + MDP(7) * t227) * t218 + (t169 * MDP(25) - t170 * MDP(27) + t278) * t183 + (-t229 + t241 - t297) * t187 + (-t180 * MDP(16) + t149 * MDP(25) - t148 * MDP(27) + t250 * t287 + t261 - t308 - t310) * t226 + (t189 * MDP(11) + (-t289 - t290) * MDP(19) + (-t148 * t221 + t149 * t225) * MDP(26) + t170 * t259 + (pkin(9) * t169 + t292) * MDP(23) + (pkin(9) * t170 + t291) * MDP(24) + t180 * MDP(17) + t251 * t287 + t243 * t152 + t237 * t187) * t222; MDP(8) + t296 * t311 + (t174 ^ 2 + t175 ^ 2 + t183 ^ 2) * MDP(28) + (t184 * t305 + t277) * t185 + 0.2e1 * (MDP(34) * t184 + MDP(35) * t185) * t171 + 0.2e1 * ((-t174 * t221 + t175 * t225) * MDP(26) + t243 * t183) * t222 + (MDP(18) * t216 - 0.2e1 * t221 * t258 + MDP(11) + 0.2e1 * (t221 * MDP(23) + t225 * MDP(24)) * pkin(9)) * t222 ^ 2 + (t226 * t257 + t237 * t311 + 0.2e1 * t229 + 0.2e1 * t262 - 0.2e1 * t263 + 0.2e1 * t297) * t226; t260 - t218 * t279 + t164 * MDP(16) - t165 * MDP(17) + t221 * t265 + (-t169 * t221 + t170 * t225) * MDP(19) + (-pkin(3) * t169 - t291) * MDP(23) + (-pkin(3) * t170 + t292) * MDP(24) + (-t152 * t225 + t169 * t200) * MDP(25) + t247 * MDP(26) + (-t152 * t221 - t170 * t200) * MDP(27) + t200 * t278 + t156 * t276 + (-t155 * t196 - t156 * t195) * MDP(30) + (t147 * t195 + t155 * t192) * MDP(34) + (t147 * t196 + t156 * t192) * MDP(35) + ((t289 - t290) * MDP(26) + t247 * MDP(28)) * pkin(10) + (-MDP(14) + (t221 * t256 - t225 * t255) * pkin(10) + t230) * t187; t208 * MDP(23) + (-t183 * t225 + t208) * MDP(25) + t246 * MDP(26) - t183 * t221 * MDP(27) + (pkin(10) * t246 + t183 * t200) * MDP(28) + t185 * t276 + (-t184 * t196 - t185 * t195) * MDP(30) + (t171 * t195 + t184 * t192) * MDP(34) + (t171 * t196 + t185 * t192) * MDP(35) + ((pkin(10) * t255 - MDP(21)) * t225 + t232 - t250) * t226 + (t221 * t259 + (-t214 + t216) * MDP(19) + (-pkin(3) * t221 - t298) * MDP(23) + (-pkin(3) * t225 + pkin(9) * t221) * MDP(24) + t243 * t200 - t251) * t222; MDP(15) + t214 * MDP(18) + (pkin(10) ^ 2 * t280 + t200 ^ 2) * MDP(28) + t192 * t195 * t304 + t280 * pkin(10) * t306 + (t192 * t303 + t195 * t305 + t276) * t196 + 0.2e1 * (MDP(23) * pkin(3) - MDP(25) * t200) * t225 + 0.2e1 * (-MDP(24) * pkin(3) - MDP(27) * t200 + t258) * t221; -t282 * MDP(25) + (-pkin(4) * t170 - qJ(5) * t169) * MDP(26) + (0.2e1 * t186 + t151) * MDP(27) + (-pkin(4) * t149 + qJ(5) * t148) * MDP(28) + (t234 + t254) * t187 + t308; (-0.2e1 * t212 - t281) * MDP(25) + t179 * MDP(27) + (-pkin(4) * t175 + qJ(5) * t174) * MDP(28) + (-t234 - 0.2e1 * t271) * t226 + (MDP(26) * t249 + t245) * t222 + t241 + t313; t248 * MDP(26) + ((MDP(28) * qJ(5) - t255) * t225 + (t256 - t295) * t221) * pkin(10) + t230; t254 + 0.2e1 * t271 + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(28) + 0.2e1 * t274 + 0.2e1 * t272 + t257; t170 * MDP(26) + MDP(28) * t149 - t187 * t236; MDP(26) * t283 + MDP(28) * t175 + t226 * t236; (MDP(28) * pkin(10) + MDP(26)) * t221; -t236 - t295; MDP(28); -t187 * MDP(33) + t312; t226 * MDP(33) + t240 - t241; t233; t235; t239; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
