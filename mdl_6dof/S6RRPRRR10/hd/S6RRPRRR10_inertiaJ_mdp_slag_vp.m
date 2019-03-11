% Calculate joint inertia matrix for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR10_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:25:27
% EndTime: 2019-03-09 14:25:32
% DurationCPUTime: 1.70s
% Computational Cost: add. (2401->270), mult. (5400->381), div. (0->0), fcn. (6290->12), ass. (0->135)
t223 = sin(pkin(12));
t225 = cos(pkin(12));
t229 = sin(qJ(4));
t282 = cos(qJ(4));
t203 = t223 * t229 - t225 * t282;
t204 = t223 * t282 + t229 * t225;
t215 = -pkin(3) * t225 - pkin(2);
t178 = pkin(4) * t203 - pkin(10) * t204 + t215;
t276 = pkin(9) + qJ(3);
t207 = t276 * t223;
t208 = t276 * t225;
t180 = -t229 * t207 + t208 * t282;
t228 = sin(qJ(5));
t232 = cos(qJ(5));
t158 = t232 * t178 - t180 * t228;
t270 = t204 * t232;
t278 = pkin(5) * t203;
t155 = -pkin(11) * t270 + t158 + t278;
t272 = t180 * t232;
t156 = t272 + (-pkin(11) * t204 + t178) * t228;
t227 = sin(qJ(6));
t231 = cos(qJ(6));
t142 = t231 * t155 - t156 * t227;
t273 = t156 * t231;
t143 = t155 * t227 + t273;
t206 = t227 * t232 + t228 * t231;
t174 = t206 * t204;
t169 = t174 * MDP(32);
t205 = t227 * t228 - t231 * t232;
t175 = t205 * t204;
t296 = t175 * MDP(31);
t249 = t203 * MDP(33) - t169 - t296;
t303 = MDP(34) * t142 - MDP(35) * t143 + t249;
t226 = cos(pkin(6));
t224 = sin(pkin(6));
t230 = sin(qJ(2));
t269 = t224 * t230;
t192 = t223 * t226 + t225 * t269;
t209 = t223 * t269;
t246 = t225 * t226 - t209;
t173 = t192 * t282 + t229 * t246;
t233 = cos(qJ(2));
t268 = t224 * t233;
t163 = t173 * t228 + t232 * t268;
t164 = t173 * t232 - t228 * t268;
t152 = t231 * t163 + t164 * t227;
t150 = t152 * MDP(32);
t153 = -t163 * t227 + t164 * t231;
t151 = t153 * MDP(31);
t302 = t151 - t150;
t253 = pkin(8) * t268;
t281 = pkin(1) * t230;
t188 = t253 + (qJ(3) + t281) * t226;
t189 = (-pkin(2) * t233 - qJ(3) * t230 - pkin(1)) * t224;
t165 = -t188 * t223 + t225 * t189;
t157 = -pkin(3) * t268 - pkin(9) * t192 + t165;
t166 = t225 * t188 + t223 * t189;
t161 = pkin(9) * t246 + t166;
t147 = t229 * t157 + t161 * t282;
t145 = -pkin(10) * t268 + t147;
t172 = t192 * t229 - t246 * t282;
t212 = pkin(8) * t269;
t280 = pkin(1) * t233;
t177 = t209 * pkin(3) + t212 + (t215 - t280) * t226;
t149 = t172 * pkin(4) - t173 * pkin(10) + t177;
t139 = -t145 * t228 + t232 * t149;
t140 = t145 * t232 + t149 * t228;
t243 = t164 * MDP(24) - t163 * MDP(25);
t300 = t139 * MDP(27) - t140 * MDP(28) + t243;
t240 = -MDP(27) * t228 - MDP(28) * t232;
t283 = pkin(10) + pkin(11);
t210 = t283 * t228;
t211 = t283 * t232;
t245 = t206 * MDP(31) - t205 * MDP(32) + (-t210 * t231 - t211 * t227) * MDP(34) - (-t210 * t227 + t211 * t231) * MDP(35);
t299 = t228 * MDP(24) + t232 * MDP(25) + pkin(10) * t240 + t245;
t298 = t173 * MDP(16);
t297 = t173 * MDP(21);
t295 = t180 * MDP(21);
t279 = pkin(5) * t172;
t137 = -pkin(11) * t164 + t139 + t279;
t138 = -pkin(11) * t163 + t140;
t134 = t231 * t137 - t138 * t227;
t274 = t138 * t231;
t135 = t137 * t227 + t274;
t293 = t134 * MDP(34) - t135 * MDP(35);
t159 = t178 * t228 + t272;
t291 = t203 * MDP(26) + t158 * MDP(27) - t159 * MDP(28);
t290 = 2 * MDP(13);
t289 = 2 * MDP(20);
t288 = 0.2e1 * MDP(27);
t287 = 0.2e1 * MDP(28);
t286 = -2 * MDP(30);
t285 = 0.2e1 * MDP(34);
t284 = 0.2e1 * MDP(35);
t277 = pkin(5) * t231;
t275 = pkin(2) * MDP(14);
t271 = t204 * t228;
t267 = t226 * MDP(8);
t196 = t205 * MDP(34);
t266 = -t206 * MDP(35) - t196;
t264 = MDP(11) * t225;
t263 = MDP(12) * t223;
t262 = MDP(19) * t233;
t261 = MDP(21) * t204;
t260 = MDP(22) * t164;
t259 = MDP(22) * t232;
t258 = MDP(29) * t175;
t257 = MDP(29) * t206;
t256 = t173 * MDP(17);
t254 = MDP(26) + MDP(33);
t252 = qJ(3) * t268;
t251 = MDP(6) * t269;
t250 = t172 * MDP(33) + t302;
t248 = MDP(23) * t228 * t232;
t247 = MDP(18) * t268;
t179 = t282 * t207 + t208 * t229;
t146 = t157 * t282 - t229 * t161;
t244 = -t165 * t223 + t166 * t225;
t242 = MDP(24) * t232 - MDP(25) * t228;
t241 = MDP(27) * t232 - MDP(28) * t228;
t144 = pkin(4) * t268 - t146;
t239 = -MDP(20) - t241;
t238 = (MDP(34) * t231 - MDP(35) * t227) * pkin(5);
t237 = -t239 + t266;
t235 = -MDP(18) + t299;
t222 = t232 ^ 2;
t221 = t228 ^ 2;
t219 = t224 ^ 2;
t217 = -pkin(5) * t232 - pkin(4);
t195 = t226 * t281 + t253;
t194 = t226 * t280 - t212;
t191 = t212 + (-pkin(2) - t280) * t226;
t167 = pkin(5) * t271 + t179;
t141 = t163 * pkin(5) + t144;
t1 = [(t165 ^ 2 + t166 ^ 2 + t191 ^ 2) * MDP(14) + t173 ^ 2 * MDP(15) + t219 * t230 ^ 2 * MDP(4) + MDP(1) + (0.2e1 * t251 + t267) * t226 + (-0.2e1 * t163 * MDP(23) + t260) * t164 + (MDP(29) * t153 + t152 * t286) * t153 + (0.2e1 * (MDP(7) * t226 - t256) * t224 + (0.2e1 * MDP(5) * t230 + t262) * t219) * t233 + 0.2e1 * (t194 * t226 + t219 * t280) * MDP(9) + 0.2e1 * (-t165 * t268 - t191 * t246) * MDP(11) + 0.2e1 * (t166 * t268 + t191 * t192) * MDP(12) + 0.2e1 * (t147 * t268 + t173 * t177) * MDP(21) - t146 * t268 * t289 + (-t165 * t192 + t166 * t246) * t290 + 0.2e1 * (-t195 * t226 - t219 * t281) * MDP(10) + (t163 * t288 + t164 * t287) * t144 + (t152 * t285 + t153 * t284) * t141 + (t134 * t285 - t135 * t284 + t139 * t288 - t140 * t287 + t254 * t172 + t177 * t289 - 0.2e1 * t150 + 0.2e1 * t151 + 0.2e1 * t243 + 0.2e1 * t247 - 0.2e1 * t298) * t172; t267 + t251 + (t144 * t270 + t164 * t179) * MDP(28) + (t144 * t271 + t163 * t179) * MDP(27) + (pkin(2) * t246 - t191 * t225 + t223 * t252) * MDP(11) + (-pkin(2) * t192 + t191 * t223 + t225 * t252) * MDP(12) + ((t223 * t192 + t225 * t246) * qJ(3) + t244) * MDP(13) + (-pkin(2) * t191 + qJ(3) * t244) * MDP(14) + (-t141 * t175 + t153 * t167) * MDP(35) + t194 * MDP(9) - t195 * MDP(10) + (t141 * t174 + t152 * t167) * MDP(34) + (t152 * t175 - t153 * t174) * MDP(30) + t215 * t297 - t153 * t258 + (t179 * MDP(20) + MDP(7) + t295) * t268 + (t177 * MDP(20) + t247 + t293 - t298 + t300 + t302) * t203 + (t215 * MDP(20) + MDP(24) * t270 - MDP(25) * t271 + t291 + t303) * t172 + (t177 * MDP(21) - MDP(17) * t268 - t172 * MDP(16) + (-t163 * t232 - t164 * t228) * MDP(23) + t173 * MDP(15) + t164 * t259) * t204; 0.2e1 * t215 * t261 + MDP(8) - (t174 * t286 - t258) * t175 + (-0.2e1 * t263 + 0.2e1 * t264 + t275) * pkin(2) + (MDP(22) * t222 + MDP(15) - 0.2e1 * t248) * t204 ^ 2 + (t270 * t287 + t271 * t288) * t179 + (t174 * t285 - t175 * t284) * t167 + (t215 * t289 - 0.2e1 * t296 - 0.2e1 * t169 + 0.2e1 * (-MDP(16) + t242) * t204 + t158 * t288 - t159 * t287 + t142 * t285 - t143 * t284 + t254 * t203) * t203 + (MDP(14) * qJ(3) + t290) * (t223 ^ 2 + t225 ^ 2) * qJ(3); -MDP(11) * t246 + t192 * MDP(12) + t191 * MDP(14) + t172 * t237 + t297; t203 * t237 + t261 + t263 - t264 - t275; MDP(14); t256 - t224 * t262 + t146 * MDP(20) - t147 * MDP(21) + t228 * t260 + (-t163 * t228 + t164 * t232) * MDP(23) + (-pkin(4) * t163 - t144 * t232) * MDP(27) + (-pkin(4) * t164 + t144 * t228) * MDP(28) + t153 * t257 + (-t152 * t206 - t153 * t205) * MDP(30) + (t141 * t205 + t152 * t217) * MDP(34) + (t141 * t206 + t153 * t217) * MDP(35) + t235 * t172; -t295 - t175 * t257 + (-t174 * t206 + t175 * t205) * MDP(30) + (t167 * t205 + t174 * t217) * MDP(34) + (t167 * t206 - t175 * t217) * MDP(35) + t239 * t179 + (MDP(17) + t228 * t259 + (-t221 + t222) * MDP(23) + t240 * pkin(4)) * t204 + t235 * t203; 0; 0.2e1 * t248 + 0.2e1 * t217 * t196 + MDP(22) * t221 + MDP(19) + 0.2e1 * t241 * pkin(4) + (t205 * t286 + t217 * t284 + t257) * t206; t172 * MDP(26) + (t172 * t277 + t134) * MDP(34) + (-t274 + (-t137 - t279) * t227) * MDP(35) + t250 + t300; (t203 * t277 + t142) * MDP(34) + (-t273 + (-t155 - t278) * t227) * MDP(35) + t242 * t204 + t249 + t291; t241 + t266; t299; 0.2e1 * t238 + t254; t250 + t293; t303; t266; t245; MDP(33) + t238; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
