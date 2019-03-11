% Calculate joint inertia matrix for
% S6RRRRRR2
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
%   see S6RRRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR2_inertiaJ_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:35:51
% EndTime: 2019-03-10 03:35:54
% DurationCPUTime: 0.99s
% Computational Cost: add. (1684->215), mult. (3127->281), div. (0->0), fcn. (3696->10), ass. (0->117)
t229 = sin(qJ(5));
t234 = cos(qJ(5));
t228 = sin(qJ(6));
t233 = cos(qJ(6));
t202 = t228 * t229 - t233 * t234;
t203 = t228 * t234 + t233 * t229;
t266 = t203 * MDP(34) - t202 * MDP(35);
t254 = t229 * MDP(27) + t234 * MDP(28) + t266;
t264 = MDP(30) * t234;
t241 = -0.2e1 * MDP(31) * t229 + 0.2e1 * t264;
t236 = cos(qJ(3));
t218 = t236 * pkin(2) + pkin(3);
t230 = sin(qJ(4));
t235 = cos(qJ(4));
t231 = sin(qJ(3));
t276 = pkin(2) * t231;
t193 = -t230 * t218 - t235 * t276;
t191 = pkin(10) - t193;
t179 = (-pkin(11) - t191) * t229;
t224 = t234 * pkin(11);
t180 = t234 * t191 + t224;
t152 = t233 * t179 - t228 * t180;
t153 = t228 * t179 + t233 * t180;
t283 = t152 * MDP(37) - t153 * MDP(38);
t216 = t230 * pkin(3) + pkin(10);
t200 = (-pkin(11) - t216) * t229;
t201 = t234 * t216 + t224;
t163 = t233 * t200 - t228 * t201;
t164 = t228 * t200 + t233 * t201;
t282 = t163 * MDP(37) - t164 * MDP(38);
t207 = (-pkin(10) - pkin(11)) * t229;
t209 = t234 * pkin(10) + t224;
t181 = t233 * t207 - t228 * t209;
t182 = t228 * t207 + t233 * t209;
t281 = t181 * MDP(37) - t182 * MDP(38);
t237 = cos(qJ(2));
t220 = -t237 * pkin(2) - pkin(1);
t232 = sin(qJ(2));
t248 = t231 * t232 - t236 * t237;
t187 = t248 * pkin(3) + t220;
t280 = 0.2e1 * t187;
t279 = 0.2e1 * t220;
t278 = -2 * MDP(33);
t277 = pkin(7) + pkin(8);
t275 = pkin(4) * t229;
t204 = t231 * t237 + t236 * t232;
t172 = t230 * t204 + t235 * t248;
t274 = t172 * pkin(5);
t273 = t235 * pkin(3);
t208 = t277 * t232;
t210 = t277 * t237;
t255 = -t236 * t208 - t231 * t210;
t157 = -t204 * pkin(9) + t255;
t249 = t231 * t208 - t236 * t210;
t158 = -t248 * pkin(9) - t249;
t143 = -t235 * t157 + t230 * t158;
t272 = t143 * t234;
t173 = t235 * t204 - t230 * t248;
t271 = t173 * t229;
t270 = t173 * t234;
t269 = t229 * t234;
t142 = t172 * pkin(4) - t173 * pkin(10) + t187;
t144 = t230 * t157 + t235 * t158;
t267 = t234 * t144;
t130 = t267 + (-pkin(11) * t173 + t142) * t229;
t268 = t233 * t130;
t265 = -t235 * t218 + t230 * t276;
t148 = t202 * t173;
t263 = MDP(32) * t148;
t147 = t203 * t173;
t145 = t147 * MDP(35);
t146 = t148 * MDP(34);
t262 = t265 * MDP(23);
t261 = t193 * MDP(24);
t260 = t230 * MDP(24);
t259 = t236 * MDP(16);
t258 = MDP(29) + MDP(36);
t257 = t172 * MDP(36) - t145 - t146;
t219 = -t234 * pkin(5) - pkin(4);
t256 = MDP(26) * t269;
t132 = t234 * t142 - t229 * t144;
t129 = -pkin(11) * t270 + t132 + t274;
t126 = t233 * t129 - t228 * t130;
t253 = -pkin(4) * t173 - pkin(10) * t172;
t226 = t229 ^ 2;
t252 = t226 * MDP(25) + MDP(22) + 0.2e1 * t256 + (MDP(32) * t203 + t202 * t278) * t203;
t190 = -pkin(4) + t265;
t251 = -t172 * t191 + t173 * t190;
t217 = -pkin(4) - t273;
t250 = -t172 * t216 + t173 * t217;
t247 = t234 * MDP(27) - t229 * MDP(28);
t245 = -MDP(30) * t229 - MDP(31) * t234;
t244 = MDP(15) + t252;
t243 = (t235 * MDP(23) - t260) * pkin(3);
t242 = (MDP(37) * t233 - MDP(38) * t228) * pkin(5);
t240 = 0.2e1 * MDP(37) * t202 + 0.2e1 * MDP(38) * t203;
t227 = t234 ^ 2;
t239 = (-t203 * t147 + t148 * t202) * MDP(33) - t203 * t263 - t143 * MDP(23) - t144 * MDP(24) + ((-t226 + t227) * MDP(26) + MDP(25) * t269 + MDP(20)) * t173 + (-MDP(21) + t254) * t172;
t238 = t204 * MDP(13) - t248 * MDP(14) + t255 * MDP(16) + t249 * MDP(17) + t239;
t225 = pkin(4) * t234;
t212 = t217 * t229;
t206 = t219 - t273;
t189 = t219 * t203;
t188 = t219 * t202;
t186 = t190 * t229;
t185 = t219 + t265;
t184 = t206 * t203;
t183 = t206 * t202;
t166 = t185 * t203;
t165 = t185 * t202;
t139 = t143 * t229;
t136 = pkin(5) * t271 + t143;
t135 = t136 * t203;
t134 = t136 * t202;
t133 = t229 * t142 + t267;
t127 = t228 * t129 + t268;
t1 = [t248 * MDP(16) * t279 + 0.2e1 * t232 * t237 * MDP(5) + MDP(1) + t232 ^ 2 * MDP(4) + t258 * t172 ^ 2 - (t147 * t278 - t263) * t148 + 0.2e1 * (-t232 * MDP(10) + t237 * MDP(9)) * pkin(1) + (MDP(11) * t204 - 0.2e1 * t248 * MDP(12) + MDP(17) * t279) * t204 + (MDP(23) * t280 - 0.2e1 * t145 - 0.2e1 * t146) * t172 + 0.2e1 * (-t133 * t172 + t143 * t270) * MDP(31) + 0.2e1 * (t132 * t172 + t143 * t271) * MDP(30) + 0.2e1 * (t126 * t172 + t136 * t147) * MDP(37) + 0.2e1 * (-t127 * t172 - t136 * t148) * MDP(38) + (MDP(24) * t280 + 0.2e1 * (-MDP(19) + t247) * t172 + (t227 * MDP(25) + MDP(18) - 0.2e1 * t256) * t173) * t173; (t251 * t229 - t272) * MDP(30) + t238 + (-t237 * MDP(10) - t232 * MDP(9)) * pkin(7) + (t185 * t147 + t152 * t172 + t134) * MDP(37) + (-t185 * t148 - t153 * t172 + t135) * MDP(38) + t232 * MDP(6) + t237 * MDP(7) + (t251 * t234 + t139) * MDP(31); MDP(8) - t190 * t241 + t185 * t240 + 0.2e1 * (-t231 * MDP(17) + t259) * pkin(2) - 0.2e1 * t262 + 0.2e1 * t261 + t244; (t250 * t234 + t139) * MDP(31) + (t250 * t229 - t272) * MDP(30) + t238 + (t206 * t147 + t163 * t172 + t134) * MDP(37) + (-t206 * t148 - t164 * t172 + t135) * MDP(38); (-t265 + t273) * MDP(23) + (t212 + t186) * MDP(31) + (t183 + t165) * MDP(37) + (t184 + t166) * MDP(38) + (-t190 - t217) * t264 + (-pkin(3) - t218) * t260 + (t259 + (-MDP(24) * t235 - MDP(17)) * t231) * pkin(2) + t244; t206 * t240 - t217 * t241 + 0.2e1 * t243 + t244; (t253 * t229 - t272) * MDP(30) + t239 + (t253 * t234 + t139) * MDP(31) + (t219 * t147 + t181 * t172 + t134) * MDP(37) + (-t219 * t148 - t182 * t172 + t135) * MDP(38); -t262 + t261 + (-t190 * t234 + t225) * MDP(30) + (t186 - t275) * MDP(31) + (t188 + t165) * MDP(37) + (t189 + t166) * MDP(38) + t252; (-t217 * t234 + t225) * MDP(30) + (t212 - t275) * MDP(31) + (t188 + t183) * MDP(37) + (t189 + t184) * MDP(38) + t243 + t252; pkin(4) * t241 + t219 * t240 + t252; t172 * MDP(29) + t132 * MDP(30) - t133 * MDP(31) + (t233 * t274 + t126) * MDP(37) + (-t268 + (-t129 - t274) * t228) * MDP(38) + t247 * t173 + t257; t245 * t191 + t254 + t283; t245 * t216 + t254 + t282; t245 * pkin(10) + t254 + t281; 0.2e1 * t242 + t258; t126 * MDP(37) - t127 * MDP(38) + t257; t266 + t283; t266 + t282; t266 + t281; MDP(36) + t242; MDP(36);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
