% Calculate joint inertia matrix for
% S6RRRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP10_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:27:31
% EndTime: 2019-03-10 02:27:39
% DurationCPUTime: 2.41s
% Computational Cost: add. (2907->350), mult. (6373->476), div. (0->0), fcn. (6989->10), ass. (0->129)
t226 = sin(qJ(4));
t229 = cos(qJ(4));
t285 = -pkin(11) - pkin(10);
t205 = t285 * t229;
t225 = sin(qJ(5));
t284 = cos(qJ(5));
t249 = t284 * t226;
t177 = -t205 * t225 - t249 * t285;
t268 = t225 * t226;
t178 = -t205 * t284 + t268 * t285;
t248 = t284 * t229;
t200 = -t248 + t268;
t201 = t225 * t229 + t249;
t254 = MDP(31) - MDP(34);
t236 = t201 * MDP(27) - t200 * MDP(28) - (MDP(30) + MDP(32)) * t177 - t178 * t254;
t302 = -pkin(10) * (t226 * MDP(23) + t229 * MDP(24)) + t226 * MDP(20) + t229 * MDP(21) + t236;
t227 = sin(qJ(3));
t300 = -0.2e1 * t227;
t235 = MDP(30) * t284 - t225 * MDP(31);
t299 = t235 * pkin(4);
t223 = sin(pkin(6));
t231 = cos(qJ(2));
t270 = t223 * t231;
t224 = cos(pkin(6));
t230 = cos(qJ(3));
t228 = sin(qJ(2));
t271 = t223 * t228;
t194 = t224 * t227 + t230 * t271;
t170 = t194 * t226 + t229 * t270;
t171 = t194 * t229 - t226 * t270;
t155 = t170 * t284 + t171 * t225;
t156 = -t225 * t170 + t171 * t284;
t298 = t156 * MDP(27) - t155 * MDP(28);
t204 = -pkin(3) * t230 - pkin(10) * t227 - pkin(2);
t199 = t229 * t204;
t278 = pkin(9) * t230;
t180 = -t226 * t278 + t199;
t252 = t229 * t278;
t181 = t204 * t226 + t252;
t188 = t201 * t227;
t182 = t188 * MDP(28);
t267 = t226 * t227;
t189 = -t225 * t267 + t227 * t248;
t183 = t189 * MDP(27);
t266 = t183 - t182;
t296 = t180 * MDP(23) - t181 * MDP(24) + t266;
t207 = pkin(8) * t271;
t282 = pkin(1) * t231;
t184 = t207 + (-pkin(2) - t282) * t224;
t193 = -t224 * t230 + t227 * t271;
t158 = t193 * pkin(3) - t194 * pkin(10) + t184;
t253 = pkin(8) * t270;
t283 = pkin(1) * t228;
t185 = t253 + (pkin(9) + t283) * t224;
t186 = (-pkin(2) * t231 - pkin(9) * t228 - pkin(1)) * t223;
t163 = t230 * t185 + t227 * t186;
t160 = -pkin(10) * t270 + t163;
t146 = t229 * t158 - t160 * t226;
t147 = t158 * t226 + t160 * t229;
t295 = t146 * MDP(23) - t147 * MDP(24);
t293 = 0.2e1 * MDP(23);
t292 = 0.2e1 * MDP(24);
t291 = -2 * MDP(26);
t290 = 0.2e1 * MDP(30);
t289 = 0.2e1 * MDP(31);
t288 = 0.2e1 * MDP(32);
t287 = 2 * MDP(33);
t286 = 2 * MDP(34);
t281 = pkin(5) * t230;
t280 = pkin(9) * t226;
t279 = pkin(9) * t229;
t277 = pkin(11) * t227;
t191 = t193 * pkin(5);
t216 = t225 * pkin(4);
t276 = MDP(17) * pkin(2);
t275 = pkin(2) * MDP(16);
t274 = qJ(6) * t230;
t162 = -t227 * t185 + t186 * t230;
t159 = pkin(3) * t270 - t162;
t273 = t159 * t226;
t272 = t159 * t229;
t269 = t224 * MDP(8);
t143 = pkin(4) * t193 - pkin(11) * t171 + t146;
t145 = -pkin(11) * t170 + t147;
t139 = t225 * t143 + t284 * t145;
t168 = -t229 * t277 + t199 + (-pkin(4) - t280) * t230;
t172 = t252 + (t204 - t277) * t226;
t154 = t225 * t168 + t284 * t172;
t203 = pkin(4) * t267 + t227 * pkin(9);
t265 = MDP(15) * t231;
t264 = MDP(25) * t189;
t263 = MDP(25) * t201;
t153 = t284 * t168 - t225 * t172;
t260 = t153 * MDP(30);
t259 = t171 * MDP(18);
t258 = t194 * MDP(13);
t257 = t229 * MDP(18);
t256 = MDP(22) + MDP(29);
t251 = t284 * pkin(4);
t190 = t193 * qJ(6);
t136 = t190 + t139;
t250 = t193 * MDP(29) + t298;
t215 = -pkin(4) * t229 - pkin(3);
t247 = t226 * t229 * MDP(19);
t246 = pkin(9) * MDP(16) - MDP(13);
t245 = pkin(9) * MDP(17) - MDP(14);
t244 = pkin(5) * t288 + MDP(29);
t138 = t284 * t143 - t225 * t145;
t137 = -t138 - t191;
t243 = t171 * MDP(20) - t170 * MDP(21);
t242 = MDP(20) * t229 - MDP(21) * t226;
t238 = t138 * MDP(30) - t139 * MDP(31);
t237 = -t154 * MDP(31) + t260;
t234 = -MDP(12) + t242;
t148 = pkin(4) * t170 + t159;
t232 = t194 * MDP(12) - t243 - t298;
t221 = t229 ^ 2;
t219 = t226 ^ 2;
t218 = t223 ^ 2;
t213 = t251 + pkin(5);
t211 = t216 + qJ(6);
t196 = t224 * t283 + t253;
t195 = t224 * t282 - t207;
t166 = pkin(5) * t200 - qJ(6) * t201 + t215;
t161 = pkin(5) * t188 - qJ(6) * t189 + t203;
t150 = -t153 + t281;
t149 = t154 - t274;
t140 = pkin(5) * t155 - qJ(6) * t156 + t148;
t1 = [(t136 ^ 2 + t137 ^ 2 + t140 ^ 2) * MDP(35) + t218 * t228 ^ 2 * MDP(4) + MDP(1) + t194 ^ 2 * MDP(11) + (0.2e1 * MDP(6) * t271 + t269) * t224 + t256 * t193 ^ 2 + (-0.2e1 * t170 * MDP(19) + t259) * t171 + (MDP(25) * t156 + t155 * t291) * t156 + (0.2e1 * MDP(5) * t228 + t265) * t218 * t231 + (t146 * t193 + t159 * t170) * t293 + (t138 * t193 + t148 * t155) * t290 + (t136 * t193 - t140 * t156) * t286 + (-t147 * t193 + t159 * t171) * t292 + (-t137 * t193 + t140 * t155) * t288 + (-t139 * t193 + t148 * t156) * t289 + (-t136 * t155 + t137 * t156) * t287 + 0.2e1 * (-t196 * t224 - t218 * t283) * MDP(10) + 0.2e1 * (-t162 * t270 + t184 * t193) * MDP(16) + 0.2e1 * (t163 * t270 + t184 * t194) * MDP(17) + 0.2e1 * (t195 * t224 + t218 * t282) * MDP(9) + 0.2e1 * (MDP(7) * t224 - t258) * t270 + 0.2e1 * (MDP(14) * t270 - t232) * t193; t156 * t264 - t194 * t276 + (t136 * t149 + t137 * t150 + t140 * t161) * MDP(35) + (-t189 * t155 - t156 * t188) * MDP(26) + (-t136 * t188 + t137 * t189 - t149 * t155 + t150 * t156) * MDP(33) + t195 * MDP(9) - t196 * MDP(10) + t269 + (-t140 * t189 - t161 * t156) * MDP(34) + (t140 * t188 + t161 * t155) * MDP(32) + (t148 * t188 + t203 * t155) * MDP(30) + (t148 * t189 + t203 * t156) * MDP(31) + (MDP(6) * t228 + MDP(7) * t231) * t223 + (-MDP(32) * t150 + MDP(34) * t149 + t237 - t275 + t296) * t193 + (-t184 * MDP(16) + t137 * MDP(32) - t136 * MDP(34) - t193 * t256 + t245 * t270 + t232 - t238 - t295) * t230 + (t171 * t257 + t194 * MDP(11) + (-t170 * t229 - t171 * t226) * MDP(19) + (pkin(9) * t170 + t273) * MDP(23) + (pkin(9) * t171 + t272) * MDP(24) + t184 * MDP(17) + t246 * t270 + t234 * t193) * t227; MDP(8) + t276 * t300 + (t149 ^ 2 + t150 ^ 2 + t161 ^ 2) * MDP(35) + (MDP(18) * t221 + t279 * t292 + t280 * t293 + MDP(11) - 0.2e1 * t247) * t227 ^ 2 + (-t149 * t287 + t161 * t288 + t203 * t290) * t188 + (t150 * t287 - t161 * t286 + t188 * t291 + t203 * t289 + t264) * t189 + (-t149 * t286 + t150 * t288 - t153 * t290 + t154 * t289 - t180 * t293 + t181 * t292 + t256 * t230 + t234 * t300 + 0.2e1 * t182 - 0.2e1 * t183 + 0.2e1 * t275) * t230; t258 - t223 * t265 + t162 * MDP(16) - t163 * MDP(17) + t226 * t259 + (-t170 * t226 + t171 * t229) * MDP(19) + (-pkin(3) * t170 - t272) * MDP(23) + (-pkin(3) * t171 + t273) * MDP(24) + t156 * t263 + (-t155 * t201 - t156 * t200) * MDP(26) + (t148 * t200 + t155 * t215) * MDP(30) + (t148 * t201 + t156 * t215) * MDP(31) + (t140 * t200 + t155 * t166) * MDP(32) + (-t136 * t200 + t137 * t201 - t155 * t178 + t156 * t177) * MDP(33) + (-t140 * t201 - t156 * t166) * MDP(34) + (t136 * t178 + t137 * t177 + t140 * t166) * MDP(35) + (-MDP(14) + t302) * t193; t189 * t263 + (-t188 * t201 - t189 * t200) * MDP(26) + (t188 * t215 + t200 * t203) * MDP(30) + (t189 * t215 + t201 * t203) * MDP(31) + (t161 * t200 + t166 * t188) * MDP(32) + (-t149 * t200 + t150 * t201 + t177 * t189 - t178 * t188) * MDP(33) + (-t161 * t201 - t166 * t189) * MDP(34) + (t149 * t178 + t150 * t177 + t161 * t166) * MDP(35) + (-t245 - t302) * t230 + (t226 * t257 + (-t219 + t221) * MDP(19) + (-pkin(3) * t226 - t279) * MDP(23) + (-pkin(3) * t229 + t280) * MDP(24) - t246) * t227; MDP(15) + t219 * MDP(18) + 0.2e1 * t247 + (t166 ^ 2 + t177 ^ 2 + t178 ^ 2) * MDP(35) + (-0.2e1 * MDP(34) * t166 + t177 * t287 + t200 * t291 + t215 * t289 + t263) * t201 + 0.2e1 * (MDP(23) * t229 - MDP(24) * t226) * pkin(3) + 0.2e1 * (MDP(30) * t215 + MDP(32) * t166 - MDP(33) * t178) * t200; t193 * MDP(22) + (t193 * t251 + t138) * MDP(30) + (-t193 * t216 - t139) * MDP(31) + (t193 * t213 - t137) * MDP(32) + (-t155 * t211 - t156 * t213) * MDP(33) + (t193 * t211 + t136) * MDP(34) + (t136 * t211 - t137 * t213) * MDP(35) + t243 + t250 + t295; t260 + t153 * MDP(32) + (-t188 * t211 - t189 * t213) * MDP(33) + (t149 * t211 - t150 * t213) * MDP(35) + t242 * t227 + ((-pkin(5) - t213) * MDP(32) + (-qJ(6) - t211) * MDP(34) - t299 - t256) * t230 - t254 * t154 + t296; (-t200 * t211 - t201 * t213) * MDP(33) + (-t177 * t213 + t178 * t211) * MDP(35) + t302; (t211 ^ 2 + t213 ^ 2) * MDP(35) + 0.2e1 * t299 + t213 * t288 + t211 * t286 + t256; (-t137 + t191) * MDP(32) + (-pkin(5) * t156 - qJ(6) * t155) * MDP(33) + (0.2e1 * t190 + t139) * MDP(34) + (-pkin(5) * t137 + qJ(6) * t136) * MDP(35) + t238 + t250; -t230 * MDP(29) + (t153 - 0.2e1 * t281) * MDP(32) + (-pkin(5) * t189 - qJ(6) * t188) * MDP(33) + (t154 - 0.2e1 * t274) * MDP(34) + (-pkin(5) * t150 + qJ(6) * t149) * MDP(35) + t237 + t266; (-pkin(5) * t201 - qJ(6) * t200) * MDP(33) + (-pkin(5) * t177 + qJ(6) * t178) * MDP(35) + t236; (0.2e1 * qJ(6) + t216) * MDP(34) + (pkin(5) * t213 + qJ(6) * t211) * MDP(35) + (MDP(32) * t284 + t235) * pkin(4) + t244; qJ(6) * t286 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(35) + t244; -t193 * MDP(32) + t156 * MDP(33) + MDP(35) * t137; MDP(32) * t230 + MDP(33) * t189 + MDP(35) * t150; MDP(33) * t201 + MDP(35) * t177; -MDP(35) * t213 - MDP(32); -MDP(35) * pkin(5) - MDP(32); MDP(35);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
