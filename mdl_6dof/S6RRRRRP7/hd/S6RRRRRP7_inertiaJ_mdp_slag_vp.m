% Calculate joint inertia matrix for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRRP7_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:43:13
% EndTime: 2019-03-10 01:43:18
% DurationCPUTime: 1.40s
% Computational Cost: add. (2207->278), mult. (4813->391), div. (0->0), fcn. (5378->10), ass. (0->119)
t218 = sin(qJ(5));
t222 = cos(qJ(5));
t257 = t218 * MDP(27) + t222 * MDP(28);
t274 = 2 * MDP(32);
t217 = cos(pkin(6));
t220 = sin(qJ(3));
t223 = cos(qJ(3));
t216 = sin(pkin(6));
t221 = sin(qJ(2));
t262 = t216 * t221;
t181 = -t217 * t223 + t220 * t262;
t182 = t217 * t220 + t223 * t262;
t219 = sin(qJ(4));
t272 = cos(qJ(4));
t161 = t181 * t272 + t182 * t219;
t162 = -t219 * t181 + t182 * t272;
t281 = t162 * MDP(20) - t161 * MDP(21);
t273 = -pkin(10) - pkin(9);
t197 = t273 * t220;
t198 = t273 * t223;
t170 = -t272 * t197 - t198 * t219;
t171 = t219 * t197 - t198 * t272;
t280 = -t170 * MDP(23) - t171 * MDP(24);
t279 = MDP(31) * t218;
t278 = -(t220 * MDP(16) + t223 * MDP(17)) * pkin(9) + t220 * MDP(13) + t223 * MDP(14);
t277 = 0.2e1 * MDP(16);
t276 = 0.2e1 * MDP(23);
t275 = 0.2e1 * MDP(24);
t271 = pkin(1) * t221;
t224 = cos(qJ(2));
t270 = pkin(1) * t224;
t269 = pkin(11) * t222;
t268 = t222 * pkin(5);
t242 = t272 * pkin(3);
t205 = -t242 - pkin(4);
t267 = pkin(4) - t205;
t261 = t216 * t224;
t243 = pkin(8) * t261;
t175 = t243 + (pkin(9) + t271) * t217;
t176 = (-pkin(2) * t224 - pkin(9) * t221 - pkin(1)) * t216;
t154 = -t220 * t175 + t223 * t176;
t244 = pkin(3) * t261;
t146 = -t182 * pkin(10) + t154 - t244;
t155 = t175 * t223 + t176 * t220;
t150 = -pkin(10) * t181 + t155;
t258 = -t272 * t146 + t219 * t150;
t135 = pkin(4) * t261 + t258;
t266 = t135 * t222;
t265 = t161 * t218;
t264 = t171 * t222;
t204 = pkin(3) * t219 + pkin(11);
t263 = t204 * t222;
t260 = t218 * t222;
t259 = t221 * MDP(6);
t211 = t222 * qJ(6);
t256 = MDP(28) * t218;
t191 = t219 * t220 - t223 * t272;
t255 = MDP(29) * t191;
t254 = MDP(30) * t222;
t192 = t219 * t223 + t220 * t272;
t253 = MDP(32) * t192;
t152 = t218 * t162 + t222 * t261;
t252 = t152 * MDP(28);
t153 = t162 * t222 - t218 * t261;
t251 = t153 * MDP(25);
t250 = t153 * MDP(27);
t249 = t161 * MDP(29);
t248 = t162 * MDP(19);
t247 = t170 * MDP(30);
t246 = t182 * MDP(11);
t245 = t218 * MDP(32);
t207 = -pkin(3) * t223 - pkin(2);
t241 = t272 * t150;
t240 = MDP(26) * t260;
t213 = t218 ^ 2;
t239 = t213 * MDP(25) + MDP(22) + 0.2e1 * t240;
t138 = t219 * t146 + t241;
t136 = -pkin(11) * t261 + t138;
t201 = pkin(8) * t262;
t174 = t201 + (-pkin(2) - t270) * t217;
t164 = t181 * pkin(3) + t174;
t140 = t161 * pkin(4) - t162 * pkin(11) + t164;
t131 = -t136 * t218 + t222 * t140;
t128 = pkin(5) * t161 - qJ(6) * t153 + t131;
t132 = t136 * t222 + t140 * t218;
t130 = -qJ(6) * t152 + t132;
t238 = -t128 * t218 + t130 * t222;
t166 = pkin(4) * t191 - pkin(11) * t192 + t207;
t148 = t222 * t166 - t171 * t218;
t237 = -pkin(4) * t192 - pkin(11) * t191;
t236 = -t191 * t204 + t192 * t205;
t235 = t182 * MDP(13) - t181 * MDP(14);
t149 = t166 * t218 + t264;
t232 = t148 * MDP(30) - t149 * MDP(31);
t231 = t254 - t279;
t230 = t218 * MDP(30) + t222 * MDP(31);
t229 = MDP(27) * t222 - MDP(19) - t256;
t228 = (MDP(23) * t272 - t219 * MDP(24)) * pkin(3);
t227 = -MDP(22) * t261 + (-t152 * t218 + t153 * t222) * MDP(26) + t218 * t251 + t281 + t257 * t161;
t144 = t264 + (-qJ(6) * t192 + t166) * t218;
t214 = t222 ^ 2;
t226 = t144 * t222 * MDP(32) + t170 * t279 + t280 + ((-t213 + t214) * MDP(26) + MDP(25) * t260 + MDP(20)) * t192 + (-MDP(21) + t257) * t191;
t225 = t131 * MDP(30) - t132 * MDP(31) + t249 + t250 - t252;
t212 = t216 ^ 2;
t206 = -pkin(4) - t268;
t196 = t211 + t269;
t195 = (-qJ(6) - pkin(11)) * t218;
t194 = t205 - t268;
t190 = t196 * t222;
t188 = t211 + t263;
t187 = (-qJ(6) - t204) * t218;
t184 = t217 * t271 + t243;
t183 = t217 * t270 - t201;
t178 = t188 * t222;
t158 = pkin(5) * t192 * t218 + t170;
t142 = pkin(5) * t191 - t192 * t211 + t148;
t134 = t135 * t218;
t133 = pkin(5) * t152 + t135;
t1 = [t162 ^ 2 * MDP(18) + t217 ^ 2 * MDP(8) + (t128 ^ 2 + t130 ^ 2 + t133 ^ 2) * MDP(33) + MDP(1) + (-0.2e1 * t181 * MDP(12) + t246) * t182 + (-0.2e1 * t152 * MDP(26) + t251) * t153 + ((MDP(4) * t221 + 0.2e1 * MDP(5) * t224) * t221 + (MDP(15) + MDP(22)) * t224 ^ 2) * t212 + (-0.2e1 * t248 + t249 + 0.2e1 * t250 - 0.2e1 * t252) * t161 + 0.2e1 * (t217 * t259 + (MDP(7) * t217 - t235 - t281) * t224) * t216 + 0.2e1 * (t183 * t217 + t212 * t270) * MDP(9) + (-t154 * t261 + t174 * t181) * t277 + 0.2e1 * (t155 * t261 + t174 * t182) * MDP(17) + (t164 * t161 + t258 * t261) * t276 + (t138 * t261 + t164 * t162) * t275 + 0.2e1 * (-t184 * t217 - t212 * t271) * MDP(10) + 0.2e1 * (t131 * t161 + t135 * t152) * MDP(30) + 0.2e1 * (-t132 * t161 + t135 * t153) * MDP(31) + (-t128 * t153 - t130 * t152) * t274; (t128 * t142 + t130 * t144 + t133 * t158) * MDP(33) + t183 * MDP(9) - t184 * MDP(10) + t217 * MDP(8) + t220 * t246 + (-t149 * t161 + t153 * t170) * MDP(31) + (-t142 * t153 - t144 * t152) * MDP(32) + (-t181 * t220 + t182 * t223) * MDP(12) + (-pkin(2) * t181 - t174 * t223) * MDP(16) + (-pkin(2) * t182 + t174 * t220) * MDP(17) + (t148 * t161 + t152 * t170) * MDP(30) + (t161 * MDP(23) + t162 * MDP(24)) * t207 + (t259 + (MDP(7) - t278 - t280) * t224) * t216 + (MDP(21) * t261 + t164 * MDP(23) + t225 - t248) * t191 + (t222 * t251 - MDP(20) * t261 + t162 * MDP(18) + (-t152 * t222 - t153 * t218) * MDP(26) + (-t128 * t222 - t130 * t218) * MDP(32) + t164 * MDP(24) + t230 * t135 + t229 * t161) * t192; MDP(8) + pkin(2) * t223 * t277 + (t142 ^ 2 + t144 ^ 2 + t158 ^ 2) * MDP(33) + (MDP(11) * t220 + 0.2e1 * t223 * MDP(12) - 0.2e1 * pkin(2) * MDP(17)) * t220 + (t207 * t275 + (-t142 * t222 - t144 * t218) * t274 + 0.2e1 * t230 * t170 + (t214 * MDP(25) + MDP(18) - 0.2e1 * t240) * t192) * t192 + (0.2e1 * t192 * t229 + t207 * t276 + 0.2e1 * t232 + t255) * t191; t154 * MDP(16) - t155 * MDP(17) + (t128 * t187 + t130 * t188 + t133 * t194) * MDP(33) + (-t241 + (-t146 + t244) * t219) * MDP(24) + t227 - MDP(15) * t261 + (t152 * t205 - t204 * t265 - t266) * MDP(30) + (t153 * t205 - t161 * t263 + t134) * MDP(31) + (-t242 * t261 - t258) * MDP(23) + (-t152 * t188 - t153 * t187 + t238) * MDP(32) + t235; (t142 * t187 + t144 * t188 + t158 * t194) * MDP(33) + t226 + (t236 * MDP(30) + (-t188 * t192 - t142) * MDP(32)) * t218 + (MDP(31) * t236 - t187 * t253 - t247) * t222 + t278; MDP(15) + (-t187 * t218 + t178) * t274 + (t187 ^ 2 + t188 ^ 2 + t194 ^ 2) * MDP(33) + t239 - 0.2e1 * t205 * t231 + 0.2e1 * t228; -t258 * MDP(23) - t138 * MDP(24) + (-pkin(4) * t152 - pkin(11) * t265 - t266) * MDP(30) + (-pkin(4) * t153 - t161 * t269 + t134) * MDP(31) + (-t152 * t196 - t153 * t195 + t238) * MDP(32) + (t128 * t195 + t130 * t196 + t133 * t206) * MDP(33) + t227; (t142 * t195 + t144 * t196 + t158 * t206) * MDP(33) + (MDP(31) * t237 - t195 * t253 - t247) * t222 + (t237 * MDP(30) + (-t192 * t196 - t142) * MDP(32)) * t218 + t226; (t178 + t190) * MDP(32) + (t187 * t195 + t188 * t196 + t194 * t206) * MDP(33) + t267 * t254 + t228 + (-t267 * MDP(31) + (-t187 - t195) * MDP(32)) * t218 + t239; (-t195 * t218 + t190) * t274 + (t195 ^ 2 + t196 ^ 2 + t206 ^ 2) * MDP(33) + 0.2e1 * t231 * pkin(4) + t239; (-t153 * MDP(32) + t128 * MDP(33)) * pkin(5) + t225; MDP(33) * pkin(5) * t142 + t255 + (-t256 + (-MDP(32) * pkin(5) + MDP(27)) * t222) * t192 + t232; -t230 * t204 + (t187 * MDP(33) - t245) * pkin(5) + t257; -t230 * pkin(11) + (MDP(33) * t195 - t245) * pkin(5) + t257; MDP(33) * pkin(5) ^ 2 + MDP(29); t133 * MDP(33); t158 * MDP(33); t194 * MDP(33); t206 * MDP(33); 0; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
