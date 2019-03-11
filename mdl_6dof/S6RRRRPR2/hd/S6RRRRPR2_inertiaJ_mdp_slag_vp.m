% Calculate joint inertia matrix for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRRRPR2_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:59:27
% EndTime: 2019-03-09 21:59:29
% DurationCPUTime: 1.03s
% Computational Cost: add. (1908->223), mult. (3477->293), div. (0->0), fcn. (4008->10), ass. (0->115)
t229 = cos(qJ(3));
t212 = pkin(2) * t229 + pkin(3);
t224 = sin(qJ(4));
t228 = cos(qJ(4));
t225 = sin(qJ(3));
t274 = pkin(2) * t225;
t182 = -t212 * t224 - t228 * t274;
t179 = qJ(5) - t182;
t221 = sin(pkin(11));
t222 = cos(pkin(11));
t263 = t221 ^ 2 + t222 ^ 2;
t268 = t263 * t179;
t209 = pkin(3) * t224 + qJ(5);
t282 = t263 * t209;
t223 = sin(qJ(6));
t227 = cos(qJ(6));
t193 = t221 * t223 - t227 * t222;
t194 = t221 * t227 + t222 * t223;
t267 = t194 * MDP(31) - t193 * MDP(32);
t281 = t193 * MDP(34) + t194 * MDP(35);
t262 = MDP(25) * t222;
t280 = -t221 * MDP(26) + t262;
t230 = cos(qJ(2));
t213 = -t230 * pkin(2) - pkin(1);
t226 = sin(qJ(2));
t241 = t225 * t226 - t229 * t230;
t177 = t241 * pkin(3) + t213;
t279 = 0.2e1 * t177;
t278 = 0.2e1 * t213;
t277 = 2 * MDP(27);
t276 = -2 * MDP(30);
t275 = pkin(7) + pkin(8);
t273 = pkin(3) * t228;
t272 = pkin(4) * t221;
t217 = t222 * pkin(10);
t271 = pkin(4) * MDP(28);
t195 = t225 * t230 + t226 * t229;
t203 = t275 * t226;
t204 = t275 * t230;
t249 = -t229 * t203 - t204 * t225;
t152 = -pkin(9) * t195 + t249;
t242 = t225 * t203 - t229 * t204;
t153 = -t241 * pkin(9) - t242;
t144 = -t152 * t228 + t153 * t224;
t270 = t144 * t222;
t161 = t228 * t195 - t224 * t241;
t269 = t161 * t221;
t160 = t195 * t224 + t228 * t241;
t143 = t160 * pkin(4) - t161 * qJ(5) + t177;
t145 = t152 * t224 + t153 * t228;
t132 = t221 * t143 + t222 * t145;
t265 = -t228 * t212 + t224 * t274;
t264 = qJ(5) * t263;
t147 = t193 * t161;
t261 = MDP(29) * t147;
t260 = MDP(33) * t160;
t146 = t194 * t161;
t259 = t146 * MDP(32);
t258 = t147 * MDP(31);
t180 = -pkin(4) + t265;
t257 = t180 * MDP(28);
t256 = t265 * MDP(23);
t255 = t182 * MDP(24);
t211 = -pkin(4) - t273;
t254 = t211 * MDP(28);
t253 = t224 * MDP(24);
t252 = t229 * MDP(16);
t210 = -pkin(5) * t222 - pkin(4);
t251 = MDP(22) + (MDP(29) * t194 + t193 * t276) * t194;
t131 = t222 * t143 - t145 * t221;
t248 = t263 * MDP(28);
t247 = MDP(15) + t251;
t246 = -pkin(4) * t161 - qJ(5) * t160;
t245 = -t131 * t221 + t132 * t222;
t244 = -t160 * t179 + t161 * t180;
t243 = -t160 * t209 + t161 * t211;
t240 = 0.2e1 * t280;
t239 = t221 * MDP(25) + t222 * MDP(26);
t128 = pkin(5) * t160 - t161 * t217 + t131;
t129 = -pkin(10) * t269 + t132;
t238 = (t128 * t227 - t129 * t223) * MDP(34) - (t128 * t223 + t129 * t227) * MDP(35);
t237 = t146 * MDP(34) - t147 * MDP(35);
t236 = -t280 + t281;
t235 = (t228 * MDP(23) - t253) * pkin(3);
t234 = 0.2e1 * t281;
t233 = t245 * MDP(27) + (-t146 * t194 + t147 * t193) * MDP(30) - t194 * t261 - t144 * MDP(23) - t145 * MDP(24) + t161 * MDP(20) + (-MDP(21) + t267) * t160;
t232 = t195 * MDP(13) - t241 * MDP(14) + t249 * MDP(16) + t242 * MDP(17) + t233;
t218 = pkin(4) * t222;
t202 = t211 * t221;
t199 = qJ(5) * t222 + t217;
t198 = (-pkin(10) - qJ(5)) * t221;
t197 = t210 - t273;
t190 = t209 * t222 + t217;
t189 = (-pkin(10) - t209) * t221;
t176 = t210 * t194;
t175 = t210 * t193;
t174 = t180 * t221;
t171 = t210 + t265;
t170 = t197 * t194;
t169 = t197 * t193;
t168 = t179 * t222 + t217;
t167 = (-pkin(10) - t179) * t221;
t164 = t198 * t223 + t199 * t227;
t163 = t198 * t227 - t199 * t223;
t157 = t171 * t194;
t156 = t171 * t193;
t155 = t189 * t223 + t190 * t227;
t154 = t189 * t227 - t190 * t223;
t149 = t167 * t223 + t168 * t227;
t148 = t167 * t227 - t168 * t223;
t140 = t144 * t221;
t135 = pkin(5) * t269 + t144;
t134 = t135 * t194;
t133 = t135 * t193;
t1 = [t226 ^ 2 * MDP(4) + (t131 ^ 2 + t132 ^ 2 + t144 ^ 2) * MDP(28) + 0.2e1 * t226 * t230 * MDP(5) + t241 * MDP(16) * t278 + MDP(1) + (MDP(18) * t161 + MDP(24) * t279) * t161 - (t146 * t276 - t261) * t147 + 0.2e1 * (-t226 * MDP(10) + t230 * MDP(9)) * pkin(1) + (MDP(11) * t195 - 0.2e1 * t241 * MDP(12) + MDP(17) * t278) * t195 + (-0.2e1 * t161 * MDP(19) + MDP(23) * t279 - 0.2e1 * t258 - 0.2e1 * t259 + t260) * t160 + 0.2e1 * t237 * t135 + 0.2e1 * (t131 * MDP(25) - t132 * MDP(26) + t238) * t160 + 0.2e1 * ((-t131 * t222 - t132 * t221) * MDP(27) + t239 * t144) * t161; (t146 * t171 + t148 * t160 + t133) * MDP(34) + (-t147 * t171 - t149 * t160 + t134) * MDP(35) + (-t230 * MDP(10) - t226 * MDP(9)) * pkin(7) + (t144 * t180 + t245 * t179) * MDP(28) + (t244 * t221 - t270) * MDP(25) + (t244 * t222 + t140) * MDP(26) + t232 + t226 * MDP(6) + t230 * MDP(7); MDP(8) + t179 ^ 2 * t248 + (-t240 + t257) * t180 + t171 * t234 + 0.2e1 * (-t225 * MDP(17) + t252) * pkin(2) - 0.2e1 * t256 + 0.2e1 * t255 + t268 * t277 + t247; (t144 * t211 + t245 * t209) * MDP(28) + (t243 * t221 - t270) * MDP(25) + (t243 * t222 + t140) * MDP(26) + t232 + (t146 * t197 + t154 * t160 + t133) * MDP(34) + (-t147 * t197 - t155 * t160 + t134) * MDP(35); (-t265 + t273) * MDP(23) + (t202 + t174) * MDP(26) + (t282 + t268) * MDP(27) + (t179 * t282 + t180 * t211) * MDP(28) + (t169 + t156) * MDP(34) + (t170 + t157) * MDP(35) + (-pkin(3) - t212) * t253 + (-t180 - t211) * t262 + (t252 + (-MDP(24) * t228 - MDP(17)) * t225) * pkin(2) + t247; t282 * t277 + t209 ^ 2 * t248 + (-t240 + t254) * t211 + t197 * t234 + 0.2e1 * t235 + t247; (t246 * t221 - t270) * MDP(25) + (t246 * t222 + t140) * MDP(26) + (-pkin(4) * t144 + t245 * qJ(5)) * MDP(28) + (t146 * t210 + t160 * t163 + t133) * MDP(34) + (-t147 * t210 - t160 * t164 + t134) * MDP(35) + t233; -t256 + t255 + (-t180 * t222 + t218) * MDP(25) + (t174 - t272) * MDP(26) + (t264 + t268) * MDP(27) + (-pkin(4) * t180 + qJ(5) * t268) * MDP(28) + (t175 + t156) * MDP(34) + (t176 + t157) * MDP(35) + t251; (-t211 * t222 + t218) * MDP(25) + (t202 - t272) * MDP(26) + (t264 + t282) * MDP(27) + (-pkin(4) * t211 + qJ(5) * t282) * MDP(28) + (t175 + t169) * MDP(34) + (t176 + t170) * MDP(35) + t235 + t251; t264 * t277 + qJ(5) ^ 2 * t248 + t210 * t234 + (t240 + t271) * pkin(4) + t251; t144 * MDP(28) + t239 * t161 + t237; t236 + t257; t236 + t254; t236 - t271; MDP(28); t238 - t258 - t259 + t260; t148 * MDP(34) - t149 * MDP(35) + t267; t154 * MDP(34) - t155 * MDP(35) + t267; t163 * MDP(34) - t164 * MDP(35) + t267; 0; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
