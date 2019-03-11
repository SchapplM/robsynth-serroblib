% Calculate joint inertia matrix for
% S6RRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRRP4_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:15:44
% EndTime: 2019-03-10 01:15:48
% DurationCPUTime: 1.25s
% Computational Cost: add. (1893->253), mult. (3417->334), div. (0->0), fcn. (3778->8), ass. (0->109)
t230 = sin(qJ(5));
t231 = sin(qJ(4));
t234 = cos(qJ(5));
t235 = cos(qJ(4));
t205 = t230 * t235 + t234 * t231;
t289 = -t230 * t231 + t234 * t235;
t262 = t205 * MDP(27) + MDP(28) * t289;
t293 = t231 * MDP(20) + t235 * MDP(21);
t292 = MDP(30) + MDP(32);
t291 = -MDP(31) + MDP(34);
t290 = pkin(8) + pkin(7);
t243 = t235 * MDP(23) - t231 * MDP(24);
t226 = t230 * pkin(4);
t214 = t226 + qJ(6);
t218 = pkin(4) * t234 + pkin(5);
t288 = (-t205 * t218 + t214 * t289) * MDP(33) + t293;
t282 = cos(qJ(2));
t222 = -t282 * pkin(2) - pkin(1);
t287 = 0.2e1 * t222;
t286 = -2 * MDP(26);
t285 = 2 * MDP(32);
t284 = 0.2e1 * MDP(33);
t283 = 2 * MDP(34);
t281 = cos(qJ(3));
t232 = sin(qJ(3));
t233 = sin(qJ(2));
t204 = t232 * t233 - t281 * t282;
t280 = t204 * pkin(4);
t199 = t204 * pkin(5);
t279 = t235 * pkin(4);
t217 = pkin(2) * t232 + pkin(9);
t227 = t235 * pkin(10);
t200 = t217 * t235 + t227;
t250 = (-pkin(10) - t217) * t231;
t165 = t200 * t230 - t234 * t250;
t277 = t165 * t204;
t166 = t234 * t200 + t230 * t250;
t276 = t166 * t204;
t211 = pkin(9) * t235 + t227;
t253 = (-pkin(10) - pkin(9)) * t231;
t181 = t211 * t230 - t234 * t253;
t275 = t181 * t204;
t210 = t290 * t233;
t212 = t290 * t282;
t182 = t281 * t210 + t232 * t212;
t274 = t182 * t235;
t183 = t234 * t211 + t230 * t253;
t273 = t183 * t204;
t184 = -t232 * t210 + t281 * t212;
t272 = t184 * t235;
t206 = t232 * t282 + t281 * t233;
t271 = t206 * t231;
t270 = t206 * t235;
t268 = t231 * t235;
t169 = t204 * pkin(3) - t206 * pkin(9) + t222;
t144 = t235 * t169 - t231 * t184;
t141 = -pkin(10) * t270 + t144 + t280;
t143 = t272 + (-pkin(10) * t206 + t169) * t231;
t134 = t230 * t141 + t234 * t143;
t198 = t204 * qJ(6);
t131 = t198 + t134;
t248 = -t234 * t141 + t143 * t230;
t132 = -t199 + t248;
t266 = t131 * t289 + t132 * t205;
t265 = t165 * t205 + t166 * t289;
t264 = t181 * t205 + t183 * t289;
t221 = -pkin(3) - t279;
t168 = -pkin(5) * t289 - t205 * qJ(6) + t221;
t254 = t281 * pkin(2);
t157 = -t254 + t168;
t263 = t157 + t168;
t220 = -t254 - pkin(3);
t209 = t220 - t279;
t261 = t209 + t221;
t159 = t289 * t206;
t260 = MDP(25) * t159;
t158 = t205 * t206;
t153 = t158 * MDP(28);
t154 = t159 * MDP(27);
t259 = t230 * MDP(31);
t256 = 0.2e1 * t282;
t255 = MDP(22) + MDP(29);
t252 = t204 * MDP(29) - t153 + t154;
t251 = MDP(19) * t268;
t249 = pkin(5) * t285 + MDP(29);
t247 = -pkin(3) * t206 - pkin(9) * t204;
t228 = t231 ^ 2;
t246 = t228 * MDP(18) + MDP(15) + 0.2e1 * t251 + (MDP(25) * t205 - t286 * t289) * t205;
t245 = -t204 * t217 + t206 * t220;
t244 = MDP(20) * t235 - MDP(21) * t231;
t242 = -MDP(23) * t231 - MDP(24) * t235;
t241 = -0.2e1 * t205 * MDP(34) - t285 * t289;
t151 = pkin(4) * t271 + t182;
t240 = -t165 * t292 + t166 * t291 + t262;
t239 = -t181 * t292 + t183 * t291 + t262;
t238 = -0.2e1 * MDP(30) * t289 + 0.2e1 * t205 * MDP(31);
t237 = (t281 * MDP(16) - t232 * MDP(17)) * pkin(2);
t229 = t235 ^ 2;
t236 = (-t158 * t205 + t159 * t289) * MDP(26) + t205 * t260 - t182 * MDP(16) - t184 * MDP(17) + ((-t228 + t229) * MDP(19) + MDP(18) * t268 + MDP(13)) * t206 + (-MDP(14) + t293 + t262) * t204;
t195 = t205 * MDP(33);
t171 = t182 * t231;
t170 = (-pkin(5) * t205 + qJ(6) * t289) * MDP(33);
t147 = t151 * t205;
t146 = t151 * t289;
t145 = t231 * t169 + t272;
t137 = t158 * pkin(5) - t159 * qJ(6) + t151;
t136 = t137 * t205;
t135 = t137 * t289;
t1 = [(t131 ^ 2 + t132 ^ 2 + t137 ^ 2) * MDP(35) + MDP(1) + pkin(1) * MDP(9) * t256 + t255 * t204 ^ 2 + (t158 * t286 + t260) * t159 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t233 + MDP(5) * t256) * t233 + (MDP(16) * t287 - 0.2e1 * t153 + 0.2e1 * t154) * t204 + 0.2e1 * (-t145 * t204 + t182 * t270) * MDP(24) + 0.2e1 * (t144 * t204 + t182 * t271) * MDP(23) + 0.2e1 * (t151 * t158 - t204 * t248) * MDP(30) + (t131 * t204 - t137 * t159) * t283 + (-t132 * t204 + t137 * t158) * t285 + 0.2e1 * (-t134 * t204 + t151 * t159) * MDP(31) + (-t131 * t158 + t132 * t159) * t284 + (MDP(17) * t287 + 0.2e1 * (-MDP(12) + t244) * t204 + (t229 * MDP(18) + MDP(11) - 0.2e1 * t251) * t206) * t206; (t245 * t235 + t171) * MDP(24) + t236 + (t245 * t231 - t274) * MDP(23) + (-t158 * t166 + t159 * t165 + t266) * MDP(33) + (-t157 * t159 - t136 + t276) * MDP(34) + (t157 * t158 - t135 - t277) * MDP(32) + (t158 * t209 - t146 - t277) * MDP(30) + (t159 * t209 + t147 - t276) * MDP(31) + t233 * MDP(6) + (t131 * t166 + t132 * t165 + t137 * t157) * MDP(35) + t282 * MDP(7) + (-t282 * MDP(10) - t233 * MDP(9)) * pkin(7); MDP(8) + t265 * t284 + (t165 ^ 2 + t166 ^ 2) * MDP(35) + t209 * t238 + (t157 * MDP(35) + t241) * t157 + t246 - 0.2e1 * t243 * t220 + 0.2e1 * t237; (t131 * t183 + t132 * t181 + t137 * t168) * MDP(35) + (t247 * t235 + t171) * MDP(24) + (t247 * t231 - t274) * MDP(23) + t236 + (-t158 * t183 + t159 * t181 + t266) * MDP(33) + (-t159 * t168 - t136 + t273) * MDP(34) + (t158 * t168 - t135 - t275) * MDP(32) + (t158 * t221 - t146 - t275) * MDP(30) + (t159 * t221 + t147 - t273) * MDP(31); (t264 + t265) * MDP(33) + (t157 * t168 + t165 * t181 + t166 * t183) * MDP(35) + t237 + (t261 * MDP(31) - t263 * MDP(34)) * t205 - (t261 * MDP(30) + t263 * MDP(32)) * t289 + t246 + t243 * (pkin(3) - t220); t264 * t284 + (t181 ^ 2 + t183 ^ 2) * MDP(35) + t221 * t238 + (t168 * MDP(35) + t241) * t168 + 0.2e1 * t243 * pkin(3) + t246; t204 * MDP(22) + t144 * MDP(23) - t145 * MDP(24) + (t234 * t280 - t248) * MDP(30) + (-t204 * t226 - t134) * MDP(31) + (t204 * t218 - t132) * MDP(32) + (-t158 * t214 - t159 * t218) * MDP(33) + (t204 * t214 + t131) * MDP(34) + (t131 * t214 - t132 * t218) * MDP(35) + t244 * t206 + t252; (-t165 * t218 + t166 * t214) * MDP(35) + t242 * t217 + t240 + t288; (-t181 * t218 + t183 * t214) * MDP(35) + t242 * pkin(9) + t239 + t288; (t214 ^ 2 + t218 ^ 2) * MDP(35) + 0.2e1 * (MDP(30) * t234 - t259) * pkin(4) + t218 * t285 + t214 * t283 + t255; -t248 * MDP(30) - t134 * MDP(31) + (-t132 + t199) * MDP(32) + (-pkin(5) * t159 - qJ(6) * t158) * MDP(33) + (0.2e1 * t198 + t134) * MDP(34) + (-pkin(5) * t132 + qJ(6) * t131) * MDP(35) + t252; t170 + (-pkin(5) * t165 + qJ(6) * t166) * MDP(35) + t240; t170 + (-pkin(5) * t181 + qJ(6) * t183) * MDP(35) + t239; (0.2e1 * qJ(6) + t226) * MDP(34) + (pkin(5) * t218 + qJ(6) * t214) * MDP(35) + (t234 * t292 - t259) * pkin(4) + t249; qJ(6) * t283 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(35) + t249; -t204 * MDP(32) + t159 * MDP(33) + t132 * MDP(35); MDP(35) * t165 + t195; MDP(35) * t181 + t195; -MDP(35) * t218 - MDP(32); -MDP(35) * pkin(5) - MDP(32); MDP(35);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
