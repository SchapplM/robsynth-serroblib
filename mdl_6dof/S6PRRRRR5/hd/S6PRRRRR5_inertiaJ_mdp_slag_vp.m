% Calculate joint inertia matrix for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_inertiaJ_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR5_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:10:05
% EndTime: 2019-03-09 01:10:09
% DurationCPUTime: 1.45s
% Computational Cost: add. (1378->256), mult. (3364->378), div. (0->0), fcn. (3862->14), ass. (0->121)
t210 = sin(qJ(5));
t215 = cos(qJ(5));
t209 = sin(qJ(6));
t214 = cos(qJ(6));
t187 = t209 * t210 - t214 * t215;
t188 = t209 * t215 + t214 * t210;
t271 = pkin(11) + pkin(12);
t193 = t271 * t210;
t194 = t271 * t215;
t233 = t188 * MDP(28) - t187 * MDP(29) + (-t214 * t193 - t209 * t194) * MDP(31) - (-t209 * t193 + t214 * t194) * MDP(32);
t284 = -(MDP(24) * t210 + MDP(25) * t215) * pkin(11) + t210 * MDP(21) + t215 * MDP(22) + t233;
t211 = sin(qJ(4));
t216 = cos(qJ(4));
t191 = -t216 * pkin(4) - t211 * pkin(11) - pkin(3);
t186 = t215 * t191;
t256 = t211 * t215;
t268 = pkin(10) * t210;
t156 = -pkin(12) * t256 + t186 + (-pkin(5) - t268) * t216;
t266 = pkin(10) * t216;
t238 = t215 * t266;
t161 = t238 + (-pkin(12) * t211 + t191) * t210;
t142 = t214 * t156 - t209 * t161;
t143 = t209 * t156 + t214 * t161;
t176 = t188 * t211;
t169 = t176 * MDP(29);
t177 = t187 * t211;
t171 = t177 * MDP(28);
t283 = t142 * MDP(31) - t143 * MDP(32) - t169 - t171;
t282 = -MDP(15) + t284;
t167 = -t210 * t266 + t186;
t168 = t210 * t191 + t238;
t281 = t167 * MDP(24) - t168 * MDP(25) + t283;
t221 = (MDP(31) * t214 - MDP(32) * t209) * pkin(5);
t205 = sin(pkin(7));
t217 = cos(qJ(3));
t259 = t205 * t217;
t207 = cos(pkin(7));
t212 = sin(qJ(3));
t260 = t205 * t212;
t181 = t211 * t207 + t216 * t260;
t159 = t210 * t181 + t215 * t259;
t160 = t215 * t181 - t210 * t259;
t144 = t214 * t159 + t209 * t160;
t145 = -t209 * t159 + t214 * t160;
t279 = t145 * MDP(28) - t144 * MDP(29);
t195 = pkin(9) * t260;
t269 = pkin(2) * t217;
t172 = t195 + (-pkin(3) - t269) * t207;
t180 = -t216 * t207 + t211 * t260;
t149 = t180 * pkin(4) - t181 * pkin(11) + t172;
t237 = pkin(9) * t259;
t270 = pkin(2) * t212;
t173 = t237 + (pkin(10) + t270) * t207;
t174 = (-pkin(3) * t217 - pkin(10) * t212 - pkin(2)) * t205;
t153 = t216 * t173 + t211 * t174;
t151 = -pkin(11) * t259 + t153;
t135 = t215 * t149 - t210 * t151;
t136 = t210 * t149 + t215 * t151;
t277 = -t135 * MDP(24) + t136 * MDP(25);
t276 = 0.2e1 * MDP(24);
t275 = 0.2e1 * MDP(25);
t274 = -2 * MDP(27);
t273 = 0.2e1 * MDP(31);
t272 = 0.2e1 * MDP(32);
t267 = pkin(10) * t215;
t265 = t180 * pkin(5);
t264 = pkin(3) * MDP(17);
t263 = pkin(10) * MDP(18);
t152 = -t211 * t173 + t216 * t174;
t150 = pkin(4) * t259 - t152;
t262 = t150 * t210;
t261 = t150 * t215;
t258 = t207 * MDP(9);
t218 = cos(qJ(2));
t257 = t207 * t218;
t255 = t212 * MDP(7);
t134 = -t159 * pkin(12) + t136;
t254 = t214 * t134;
t206 = sin(pkin(6));
t208 = cos(pkin(6));
t213 = sin(qJ(2));
t158 = t208 * t260 + (t212 * t257 + t213 * t217) * t206;
t178 = -t206 * t218 * t205 + t208 * t207;
t148 = t158 * t216 + t178 * t211;
t157 = -t208 * t259 + (t212 * t213 - t217 * t257) * t206;
t137 = -t148 * t210 + t157 * t215;
t138 = t148 * t215 + t157 * t210;
t130 = t214 * t137 - t209 * t138;
t131 = t209 * t137 + t214 * t138;
t253 = t130 * MDP(31) - t131 * MDP(32);
t251 = MDP(16) * t217;
t246 = t160 * MDP(19);
t245 = t177 * MDP(26);
t244 = t181 * MDP(14);
t243 = t187 * MDP(31);
t242 = t188 * MDP(26);
t241 = t210 * MDP(19);
t240 = t211 * MDP(18);
t239 = MDP(23) + MDP(30);
t236 = t180 * MDP(30) + t279;
t235 = t210 * t215 * MDP(20);
t234 = pkin(10) * MDP(17) - MDP(14);
t133 = -t160 * pkin(12) + t135 + t265;
t126 = t214 * t133 - t209 * t134;
t232 = t160 * MDP(21) - t159 * MDP(22);
t231 = t215 * MDP(21) - t210 * MDP(22);
t228 = t215 * MDP(24) - t210 * MDP(25);
t127 = t209 * t133 + t254;
t226 = -t126 * MDP(31) + t127 * MDP(32);
t223 = -MDP(13) + t231;
t220 = t181 * MDP(13) - t232 - t279;
t203 = t215 ^ 2;
t201 = t210 ^ 2;
t200 = t205 ^ 2;
t199 = -t215 * pkin(5) - pkin(4);
t190 = (pkin(5) * t210 + pkin(10)) * t211;
t183 = t207 * t270 + t237;
t182 = t207 * t269 - t195;
t147 = t158 * t211 - t178 * t216;
t139 = t159 * pkin(5) + t150;
t1 = [MDP(1); (-t157 * t207 - t178 * t259) * MDP(10) + (-t158 * t207 + t178 * t260) * MDP(11) + (t147 * t259 + t157 * t180) * MDP(17) + (t148 * t259 + t157 * t181) * MDP(18) + (t137 * t180 + t147 * t159) * MDP(24) + (-t138 * t180 + t147 * t160) * MDP(25) + (t130 * t180 + t147 * t144) * MDP(31) + (-t131 * t180 + t147 * t145) * MDP(32) + (t218 * MDP(3) - t213 * MDP(4)) * t206; t200 * t212 ^ 2 * MDP(5) + t181 ^ 2 * MDP(12) + MDP(2) + (0.2e1 * t205 * t255 + t258) * t207 + t239 * t180 ^ 2 + (-0.2e1 * t159 * MDP(20) + t246) * t160 + (t145 * MDP(26) + t144 * t274) * t145 + (0.2e1 * MDP(6) * t212 + t251) * t200 * t217 + 0.2e1 * (-t152 * t259 + t172 * t180) * MDP(17) + 0.2e1 * (t153 * t259 + t172 * t181) * MDP(18) + 0.2e1 * (t182 * t207 + t200 * t269) * MDP(10) + 0.2e1 * (-t183 * t207 - t200 * t270) * MDP(11) + (t135 * t180 + t150 * t159) * t276 + (t126 * t180 + t139 * t144) * t273 + (-t136 * t180 + t150 * t160) * t275 + (-t127 * t180 + t139 * t145) * t272 + 0.2e1 * (MDP(8) * t207 - t244) * t259 + 0.2e1 * (MDP(15) * t259 - t220) * t180; -t158 * MDP(11) + (t147 * t210 * t211 - t137 * t216) * MDP(24) + (t138 * t216 + t147 * t256) * MDP(25) + (-t130 * t216 + t147 * t176) * MDP(31) + (t131 * t216 - t147 * t177) * MDP(32) + (-t216 * MDP(17) - MDP(10) + t240) * t157; (t177 * t144 - t145 * t176) * MDP(27) + t182 * MDP(10) - t183 * MDP(11) + t258 - t145 * t245 + (t139 * t176 + t190 * t144) * MDP(31) + (-t139 * t177 + t190 * t145) * MDP(32) - pkin(3) * t181 * MDP(18) + (t217 * MDP(8) + t255) * t205 + (-t264 + t281) * t180 + (-t172 * MDP(17) + (-MDP(15) + t263) * t259 - t239 * t180 + t220 + t226 + t277) * t216 + (t181 * MDP(12) + t215 * t246 + (-t159 * t215 - t160 * t210) * MDP(20) + (pkin(10) * t159 + t262) * MDP(24) + (pkin(10) * t160 + t261) * MDP(25) + t172 * MDP(18) + t234 * t259 + t223 * t180) * t211; t190 * t176 * t273 - 0.2e1 * pkin(3) * t240 + MDP(9) - (t176 * t274 + t190 * t272 - t245) * t177 + (t203 * MDP(19) + t267 * t275 + t268 * t276 + MDP(12) - 0.2e1 * t235) * t211 ^ 2 + (-t142 * t273 + t143 * t272 - t167 * t276 + t168 * t275 - 0.2e1 * t211 * t223 + t216 * t239 + 0.2e1 * t169 + 0.2e1 * t171 + 0.2e1 * t264) * t216; -t148 * MDP(18) + (t188 * MDP(32) - MDP(17) - t228 + t243) * t147; t244 - t205 * t251 + t152 * MDP(17) - t153 * MDP(18) + t160 * t241 + (-t210 * t159 + t160 * t215) * MDP(20) + (-pkin(4) * t159 - t261) * MDP(24) + (-pkin(4) * t160 + t262) * MDP(25) + t145 * t242 + (-t188 * t144 - t145 * t187) * MDP(27) + (t139 * t187 + t199 * t144) * MDP(31) + (t139 * t188 + t199 * t145) * MDP(32) + t282 * t180; -t177 * t242 + (-t188 * t176 + t177 * t187) * MDP(27) + (t199 * t176 + t190 * t187) * MDP(31) + (-t199 * t177 + t190 * t188) * MDP(32) + (-t263 - t282) * t216 + (t215 * t241 + (-t201 + t203) * MDP(20) + (-pkin(4) * t210 - t267) * MDP(24) + (-pkin(4) * t215 + t268) * MDP(25) - t234) * t211; 0.2e1 * t235 + 0.2e1 * t199 * t243 + t201 * MDP(19) + MDP(16) + 0.2e1 * t228 * pkin(4) + (t187 * t274 + t199 * t272 + t242) * t188; t137 * MDP(24) - t138 * MDP(25) + t253; t180 * MDP(23) + (t214 * t265 + t126) * MDP(31) + (-t254 + (-t133 - t265) * t209) * MDP(32) + t232 + t236 - t277; (-t239 - t221) * t216 + t231 * t211 + t281; t284; 0.2e1 * t221 + t239; t253; -t226 + t236; -t216 * MDP(30) + t283; t233; MDP(30) + t221; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
