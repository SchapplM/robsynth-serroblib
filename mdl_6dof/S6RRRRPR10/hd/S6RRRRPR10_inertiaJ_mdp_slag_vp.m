% Calculate joint inertia matrix for
% S6RRRRPR10
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
%   see S6RRRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR10_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:08:35
% EndTime: 2019-03-09 23:08:40
% DurationCPUTime: 1.71s
% Computational Cost: add. (1830->270), mult. (3996->370), div. (0->0), fcn. (4397->10), ass. (0->125)
t218 = sin(qJ(3));
t277 = -pkin(10) - pkin(9);
t196 = t277 * t218;
t222 = cos(qJ(3));
t197 = t277 * t222;
t217 = sin(qJ(4));
t221 = cos(qJ(4));
t175 = -t221 * t196 - t197 * t217;
t176 = t196 * t217 - t197 * t221;
t193 = t217 * t218 - t221 * t222;
t194 = t217 * t222 + t218 * t221;
t248 = MDP(24) - MDP(27);
t292 = MDP(26) - MDP(23);
t296 = t194 * MDP(20) - t193 * MDP(21) + t292 * t175 - t176 * t248;
t220 = cos(qJ(6));
t209 = t220 * MDP(31);
t216 = sin(qJ(6));
t293 = -t216 * MDP(32) + t209;
t215 = cos(pkin(6));
t214 = sin(pkin(6));
t219 = sin(qJ(2));
t269 = t214 * t219;
t185 = -t215 * t222 + t218 * t269;
t186 = t215 * t218 + t222 * t269;
t166 = t221 * t185 + t186 * t217;
t223 = cos(qJ(2));
t268 = t214 * t223;
t153 = t166 * t220 + t216 * t268;
t154 = t166 * t216 - t220 * t268;
t167 = -t185 * t217 + t186 * t221;
t259 = t154 * MDP(29);
t246 = pkin(8) * t268;
t276 = pkin(1) * t219;
t181 = t246 + (pkin(9) + t276) * t215;
t182 = (-pkin(2) * t223 - pkin(9) * t219 - pkin(1)) * t214;
t156 = -t181 * t218 + t222 * t182;
t148 = -pkin(3) * t268 - pkin(10) * t186 + t156;
t157 = t181 * t222 + t182 * t218;
t151 = -pkin(10) * t185 + t157;
t264 = -t221 * t148 + t217 * t151;
t288 = t167 * MDP(20) - t166 * MDP(21);
t291 = -t264 * MDP(23) + (t153 * t220 - t154 * t216) * MDP(30) + t220 * t259 + t288 + t293 * t167;
t290 = 0.2e1 * t167;
t289 = pkin(3) * (t221 * MDP(23) - t217 * MDP(24));
t286 = -(MDP(16) * t218 + MDP(17) * t222) * pkin(9) + t218 * MDP(13) + t222 * MDP(14);
t285 = 0.2e1 * MDP(16);
t284 = 2 * MDP(25);
t283 = -2 * MDP(26);
t282 = 2 * MDP(26);
t281 = 2 * MDP(27);
t280 = 2 * MDP(34);
t279 = 2 * MDP(35);
t278 = pkin(4) + pkin(11);
t275 = pkin(1) * t223;
t274 = (MDP(28) * pkin(4));
t273 = qJ(5) * t193;
t165 = -pkin(5) * t193 + t176;
t159 = t165 * t216;
t160 = t165 * t220;
t206 = -pkin(3) * t221 - pkin(4);
t201 = -pkin(11) + t206;
t272 = t167 * t201;
t271 = t167 * t278;
t203 = pkin(3) * t217 + qJ(5);
t270 = t193 * t203;
t267 = t216 * t220;
t266 = t219 * MDP(6);
t265 = qJ(5) + t203;
t141 = t217 * t148 + t221 * t151;
t263 = MDP(28) * t206;
t262 = MDP(34) * t216;
t261 = MDP(35) * t220;
t260 = t153 * MDP(32);
t258 = t166 * MDP(19);
t207 = -pkin(3) * t222 - pkin(2);
t232 = -qJ(5) * t194 + t207;
t170 = pkin(4) * t193 + t232;
t257 = t170 * MDP(26);
t256 = t186 * MDP(11);
t255 = t194 * MDP(32);
t254 = t207 * MDP(23);
t252 = t220 * MDP(34);
t251 = MDP(15) + MDP(22);
t250 = MDP(18) + MDP(33);
t247 = 0.2e1 * t194;
t245 = qJ(5) * t268;
t200 = pkin(4) * t268;
t139 = t200 + t264;
t244 = MDP(30) * t267;
t212 = t220 ^ 2;
t243 = t212 * MDP(29) + MDP(22) - 0.2e1 * t244;
t242 = t194 * t278 + t273;
t241 = -t194 * t201 + t270;
t240 = t186 * MDP(13) - t185 * MDP(14);
t236 = t207 * MDP(24) - t170 * MDP(27);
t235 = MDP(31) * t216 + MDP(32) * t220;
t155 = t193 * t278 + t232;
t164 = pkin(5) * t194 + t175;
t144 = -t155 * t216 + t164 * t220;
t145 = t155 * t220 + t164 * t216;
t234 = t144 * MDP(34) - t145 * MDP(35);
t233 = -t216 * MDP(35) + t252;
t138 = t245 - t141;
t199 = pkin(8) * t269;
t180 = t199 + (-pkin(2) - t275) * t215;
t231 = -MDP(19) + t235;
t230 = MDP(25) + t233;
t229 = t281 + 0.2e1 * t261 + 0.2e1 * t262;
t168 = pkin(3) * t185 + t180;
t133 = pkin(5) * t167 + pkin(11) * t268 + t139;
t226 = -qJ(5) * t167 + t168;
t137 = t166 * t278 + t226;
t131 = t133 * t220 - t137 * t216;
t132 = t133 * t216 + t137 * t220;
t227 = t154 * MDP(31) + t131 * MDP(34) - t132 * MDP(35) + t260;
t211 = t216 ^ 2;
t225 = t160 * MDP(35) + t194 * t209 + ((-t211 + t212) * MDP(30) + MDP(29) * t267) * t193 + t296;
t210 = t214 ^ 2;
t188 = t215 * t276 + t246;
t187 = t215 * t275 - t199;
t142 = pkin(4) * t166 + t226;
t136 = -pkin(5) * t166 - t138;
t135 = t136 * t220;
t134 = t136 * t216;
t1 = [t215 ^ 2 * MDP(8) + (t138 ^ 2 + t139 ^ 2 + t142 ^ 2) * MDP(28) + MDP(1) + (-0.2e1 * t185 * MDP(12) + t256) * t186 + t250 * t167 ^ 2 + (0.2e1 * t153 * MDP(30) + MDP(31) * t290 + t259) * t154 + ((MDP(4) * t219 + 0.2e1 * MDP(5) * t223) * t219 + t251 * t223 ^ 2) * t210 + 0.2e1 * t215 * t266 * t214 + (-t139 * t268 - t142 * t166) * t282 + (t138 * t268 - t142 * t167) * t281 + 0.2e1 * (t166 * t168 + t264 * t268) * MDP(23) + 0.2e1 * (t141 * t268 + t167 * t168) * MDP(24) + (-t156 * t268 + t180 * t185) * t285 + 0.2e1 * (t157 * t268 + t180 * t186) * MDP(17) + 0.2e1 * (t187 * t215 + t210 * t275) * MDP(9) + 0.2e1 * (-t188 * t215 - t210 * t276) * MDP(10) + (t131 * t167 - t136 * t153) * t280 + (t138 * t166 + t139 * t167) * t284 + (-t132 * t167 + t136 * t154) * t279 + (-t258 + t260) * t290 + 0.2e1 * (MDP(7) * t215 - t240 - t288) * t268; (-t138 * t176 + t139 * t175 + t142 * t170) * MDP(28) + t187 * MDP(9) - t188 * MDP(10) + t215 * MDP(8) + (-t185 * t218 + t186 * t222) * MDP(12) + (-pkin(2) * t185 - t180 * t222) * MDP(16) + (-pkin(2) * t186 + t180 * t218) * MDP(17) + t218 * t256 + (-t153 * MDP(34) + t154 * MDP(35)) * t165 + (t175 * MDP(25) + t234 + t236) * t167 + (-t176 * MDP(25) + t254 - t257) * t166 + (t168 * MDP(24) + t139 * MDP(25) - t142 * MDP(27) + t167 * t250 + t227 - t258) * t194 + (t138 * MDP(25) + (t153 * t216 + t154 * t220) * MDP(30) - t142 * MDP(26) + t168 * MDP(23) + t216 * t259 - t233 * t136 + t231 * t167) * t193 + (t266 + (MDP(7) - t286 - t296) * t223) * t214; MDP(8) + pkin(2) * t222 * t285 + (t170 ^ 2 + t175 ^ 2 + t176 ^ 2) * MDP(28) + t236 * t247 + (MDP(11) * t218 + 0.2e1 * t222 * MDP(12) - 0.2e1 * pkin(2) * MDP(17)) * t218 + (t144 * t280 - t145 * t279 + t175 * t284 + t194 * t250) * t194 + (t159 * t279 - t160 * t280 - t176 * t284 + t231 * t247 + 0.2e1 * t254 - 0.2e1 * t257 + (t211 * MDP(29) + 0.2e1 * t244) * t193) * t193; (-t166 * t203 + t167 * t206) * MDP(25) + (-t138 * t203 + t139 * t206) * MDP(28) + (t154 * t203 - t216 * t272 + t135) * MDP(35) + (-t153 * t203 + t220 * t272 + t134) * MDP(34) + t156 * MDP(16) - t157 * MDP(17) + (-t206 * MDP(26) - MDP(27) * t265 - t251 - t289) * t268 + t139 * MDP(26) + t240 - t248 * t141 + t291; t225 + (t194 * t206 - t270) * MDP(25) + (t175 * t206 + t176 * t203) * MDP(28) + (-t220 * t241 + t159) * MDP(34) + (MDP(35) * t241 - t255) * t216 + t286; MDP(15) + (t282 + t263) * t206 + 0.2e1 * t289 + (MDP(28) * t203 + t229) * t203 + t243; -MDP(22) * t268 - t141 * MDP(24) + (-pkin(4) * t167 - qJ(5) * t166) * MDP(25) + (0.2e1 * t200 + t264) * MDP(26) + (-0.2e1 * t245 + t141) * MDP(27) + (-pkin(4) * t139 - qJ(5) * t138) * MDP(28) + (-qJ(5) * t153 - t220 * t271 + t134) * MDP(34) + (qJ(5) * t154 + t216 * t271 + t135) * MDP(35) + t291; t225 + (-pkin(4) * t175 + qJ(5) * t176) * MDP(28) + (-pkin(4) * t194 - t273) * MDP(25) + (-t220 * t242 + t159) * MDP(34) + (MDP(35) * t242 - t255) * t216; pkin(4) * t283 + qJ(5) * t281 + (-pkin(4) * t206 + qJ(5) * t203) * MDP(28) + (-t217 * t248 - t221 * t292) * pkin(3) + t243 + (t261 + t262) * t265; (t283 + t274) * pkin(4) + (MDP(28) * qJ(5) + t229) * qJ(5) + t243; -MDP(26) * t268 + t139 * MDP(28) + t167 * t230; t175 * MDP(28) + t194 * t230; MDP(26) + t263; MDP(26) - t274; MDP(28); t167 * MDP(33) + t227; t194 * MDP(33) + t193 * t235 + t234; t201 * t233 + t293; -t278 * t252 + t209 + (MDP(35) * t278 - MDP(32)) * t216; t233; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
