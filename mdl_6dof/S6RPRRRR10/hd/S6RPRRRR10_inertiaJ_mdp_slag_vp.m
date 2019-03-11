% Calculate joint inertia matrix for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR10_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:34:14
% EndTime: 2019-03-09 07:34:19
% DurationCPUTime: 1.64s
% Computational Cost: add. (3151->255), mult. (8313->387), div. (0->0), fcn. (9777->14), ass. (0->126)
t233 = sin(qJ(6));
t237 = cos(qJ(6));
t276 = t233 * MDP(31) + t237 * MDP(32);
t250 = MDP(34) * t237 - MDP(35) * t233;
t232 = sin(pkin(6));
t288 = cos(pkin(13));
t289 = cos(pkin(7));
t253 = t289 * t288;
t231 = sin(pkin(7));
t290 = cos(pkin(6));
t255 = t290 * t231;
t304 = t232 * t253 + t255;
t230 = sin(pkin(13));
t236 = sin(qJ(3));
t239 = cos(qJ(3));
t187 = t236 * t255 + (t230 * t239 + t236 * t253) * t232;
t256 = t232 * t288;
t200 = t231 * t256 - t290 * t289;
t235 = sin(qJ(4));
t238 = cos(qJ(4));
t170 = t187 * t235 + t200 * t238;
t171 = t187 * t238 - t200 * t235;
t234 = sin(qJ(5));
t294 = cos(qJ(5));
t162 = t294 * t170 + t171 * t234;
t303 = MDP(33) * t162;
t293 = pkin(1) * t230;
t203 = qJ(2) * t256 + t290 * t293;
t182 = t304 * pkin(9) + t203;
t261 = pkin(1) * t288;
t218 = t290 * t261;
t280 = t230 * t232;
t188 = t290 * pkin(2) + t218 + (-t289 * pkin(9) - qJ(2)) * t280;
t195 = (-pkin(9) * t230 * t231 - t288 * pkin(2) - pkin(1)) * t232;
t247 = t289 * t188 + t195 * t231;
t164 = -t236 * t182 + t247 * t239;
t295 = pkin(10) + pkin(11);
t214 = t295 * t238;
t259 = t294 * t235;
t193 = t214 * t234 + t295 * t259;
t277 = t234 * t235;
t194 = t294 * t214 - t295 * t277;
t209 = -t294 * t238 + t277;
t210 = t234 * t238 + t259;
t302 = t210 * MDP(24) - t209 * MDP(25) - t193 * MDP(27) - t194 * MDP(28);
t298 = 0.2e1 * MDP(27);
t297 = 0.2e1 * MDP(34);
t296 = 0.2e1 * MDP(35);
t186 = t236 * t280 - t304 * t239;
t292 = pkin(4) * t186;
t168 = -t188 * t231 + t289 * t195;
t155 = pkin(3) * t186 - pkin(10) * t187 + t168;
t165 = t239 * t182 + t247 * t236;
t157 = -t200 * pkin(10) + t165;
t148 = t238 * t155 - t157 * t235;
t145 = -pkin(11) * t171 + t148 + t292;
t149 = t155 * t235 + t157 * t238;
t147 = -pkin(11) * t170 + t149;
t141 = t294 * t145 - t234 * t147;
t139 = -pkin(5) * t186 - t141;
t287 = t139 * t237;
t163 = -t234 * t170 + t294 * t171;
t153 = t163 * t237 + t186 * t233;
t286 = t153 * t233;
t285 = t162 * t233;
t284 = t162 * t237;
t283 = t193 * t237;
t282 = t210 * t233;
t281 = t210 * t237;
t279 = t231 * t236;
t278 = t231 * t239;
t275 = MDP(15) * t235;
t274 = MDP(20) * t238;
t273 = MDP(28) * t210;
t272 = MDP(29) * t237;
t152 = t163 * t233 - t186 * t237;
t269 = t152 * MDP(32);
t268 = t153 * MDP(31);
t160 = t162 * MDP(25);
t267 = t163 * MDP(23);
t161 = t163 * MDP(24);
t266 = t170 * MDP(18);
t265 = t171 * MDP(17);
t264 = t187 * MDP(10);
t263 = t209 * MDP(33);
t262 = t294 * pkin(4);
t223 = -pkin(4) * t238 - pkin(3);
t260 = t294 * t147;
t258 = t233 * t237 * MDP(30);
t228 = t233 ^ 2;
t257 = t228 * MDP(29) + MDP(26) + 0.2e1 * t258;
t254 = -pkin(5) * t210 - pkin(12) * t209;
t221 = pkin(4) * t234 + pkin(12);
t222 = -t262 - pkin(5);
t252 = -t209 * t221 + t210 * t222;
t251 = MDP(31) * t237 - MDP(32) * t233;
t249 = MDP(34) * t233 + MDP(35) * t237;
t142 = t234 * t145 + t260;
t246 = -MDP(23) + t251;
t204 = -t235 * t279 + t238 * t289;
t205 = t235 * t289 + t238 * t279;
t177 = -t294 * t204 + t205 * t234;
t178 = t234 * t204 + t294 * t205;
t245 = -t178 * MDP(28) + (-MDP(27) - t250) * t177;
t244 = (-t152 * t233 + t153 * t237) * MDP(30) + MDP(29) * t286 + t161 + t186 * MDP(26) - t160 + t276 * t162;
t243 = (t294 * MDP(27) - t234 * MDP(28)) * pkin(4);
t229 = t237 ^ 2;
t242 = (-t228 + t229) * t210 * MDP(30) + t272 * t282 + t302 + t276 * t209;
t156 = t200 * pkin(3) - t164;
t241 = t235 * MDP(17) + t238 * MDP(18) + (-MDP(20) * t235 - MDP(21) * t238) * pkin(10);
t140 = t186 * pkin(12) + t142;
t150 = t170 * pkin(4) + t156;
t143 = t162 * pkin(5) - t163 * pkin(12) + t150;
t136 = -t140 * t233 + t143 * t237;
t137 = t140 * t237 + t143 * t233;
t240 = t136 * MDP(34) - t137 * MDP(35) + t268 - t269 + t303;
t227 = t232 ^ 2;
t202 = -qJ(2) * t280 + t218;
t189 = t193 * t233;
t183 = pkin(5) * t209 - pkin(12) * t210 + t223;
t173 = t237 * t178 - t233 * t278;
t172 = -t233 * t178 - t237 * t278;
t167 = t183 * t233 + t194 * t237;
t166 = t183 * t237 - t194 * t233;
t138 = t139 * t233;
t1 = [0.2e1 * (t202 * t290 + t227 * t261) * MDP(4) + 0.2e1 * (-t203 * t290 - t227 * t293) * MDP(5) + t200 ^ 2 * MDP(12) - 0.2e1 * t171 * t170 * MDP(16) + 0.2e1 * (t148 * t186 + t156 * t170) * MDP(20) + 0.2e1 * (-t149 * t186 + t156 * t171) * MDP(21) + 0.2e1 * (-t142 * t186 + t150 * t163) * MDP(28) + 0.2e1 * (t165 * t200 + t168 * t187) * MDP(14) + 0.2e1 * (-t164 * t200 + t168 * t186) * MDP(13) + MDP(1) - 0.2e1 * t200 * t264 + t139 * t152 * t297 + t141 * t186 * t298 + 0.2e1 * (-t202 * t230 + t288 * t203) * MDP(6) * t232 + t163 ^ 2 * MDP(22) + t171 ^ 2 * MDP(15) + t187 ^ 2 * MDP(8) + (pkin(1) ^ 2 * t227 + t202 ^ 2 + t203 ^ 2) * MDP(7) - 0.2e1 * (t269 + t267) * t162 + (MDP(19) + MDP(26)) * t186 ^ 2 + (t136 * t297 - t137 * t296 + t150 * t298 + 0.2e1 * t268 + t303) * t162 + (MDP(29) * t153 - 0.2e1 * t152 * MDP(30) + t139 * t296) * t153 - 0.2e1 * (t187 * MDP(9) + t160 + t266) * t186 + 0.2e1 * (t200 * MDP(11) + t161 + t265) * t186; (t289 * t186 - t200 * t278) * MDP(13) + (t289 * t187 + t200 * t279) * MDP(14) + (-t170 * t278 + t204 * t186) * MDP(20) + (-t171 * t278 - t205 * t186) * MDP(21) + (-t162 * t278 - t177 * t186) * MDP(27) + (-t163 * t278 - t178 * t186) * MDP(28) + (t152 * t177 + t162 * t172) * MDP(34) + (t153 * t177 - t162 * t173) * MDP(35) + (-t288 * MDP(4) + t230 * MDP(5) - pkin(1) * MDP(7)) * t232; MDP(7); t264 - t200 * MDP(12) + t164 * MDP(13) - t165 * MDP(14) + (-t170 * t235 + t171 * t238) * MDP(16) + (-pkin(3) * t170 - t156 * t238) * MDP(20) + (-pkin(3) * t171 + t156 * t235) * MDP(21) + (t152 * t193 + t162 * t166) * MDP(34) + (t153 * t193 - t162 * t167) * MDP(35) + t171 * t275 + (t162 * MDP(27) + t163 * MDP(28)) * t223 + (t150 * MDP(27) + t240 - t267) * t209 + (t153 * t272 + t150 * MDP(28) + (-t152 * t237 - t286) * MDP(30) + t163 * MDP(22) + t249 * t139 + t246 * t162) * t210 + (-MDP(11) + t241 + t302) * t186; (t172 * t209 + t177 * t282) * MDP(34) + (-t173 * t209 + t177 * t281) * MDP(35) + (-t236 * MDP(14) + (-MDP(21) * t235 - MDP(27) * t209 + MDP(13) - t273 + t274) * t239) * t231; 0.2e1 * pkin(3) * t274 + 0.2e1 * t223 * t273 + MDP(12) + (0.2e1 * t246 * t210 + t223 * t298 + t263) * t209 + (t166 * t209 + t193 * t282) * t297 + (-t167 * t209 + t193 * t281) * t296 + (0.2e1 * MDP(16) * t238 - 0.2e1 * MDP(21) * pkin(3) + t275) * t235 + (MDP(29) * t229 + MDP(22) - 0.2e1 * t258) * t210 ^ 2; t265 - t266 + t186 * MDP(19) + t148 * MDP(20) - t149 * MDP(21) + (t186 * t262 + t141) * MDP(27) + (-t260 + (-t145 - t292) * t234) * MDP(28) + (t152 * t222 - t221 * t285 - t287) * MDP(34) + (t153 * t222 - t221 * t284 + t138) * MDP(35) + t244; MDP(20) * t204 - MDP(21) * t205 + t245; (t252 * t233 - t283) * MDP(34) + (t252 * t237 + t189) * MDP(35) + t241 + t242; -0.2e1 * t222 * t250 + MDP(19) + 0.2e1 * t243 + t257; t141 * MDP(27) - t142 * MDP(28) + (-pkin(5) * t152 - pkin(12) * t285 - t287) * MDP(34) + (-pkin(5) * t153 - pkin(12) * t284 + t138) * MDP(35) + t244; t245; (t254 * t233 - t283) * MDP(34) + (t254 * t237 + t189) * MDP(35) + t242; t243 + t257 + t250 * (pkin(5) - t222); 0.2e1 * pkin(5) * t250 + t257; t240; MDP(34) * t172 - MDP(35) * t173; MDP(34) * t166 - MDP(35) * t167 + t251 * t210 + t263; -t249 * t221 + t276; -t249 * pkin(12) + t276; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
