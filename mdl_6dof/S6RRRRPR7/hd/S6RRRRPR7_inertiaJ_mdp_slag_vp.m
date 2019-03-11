% Calculate joint inertia matrix for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_inertiaJ_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR7_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:34:12
% EndTime: 2019-03-09 22:34:17
% DurationCPUTime: 1.45s
% Computational Cost: add. (2785->254), mult. (6139->380), div. (0->0), fcn. (7107->12), ass. (0->119)
t221 = sin(qJ(6));
t225 = cos(qJ(6));
t263 = t221 * MDP(29) + t225 * MDP(30);
t222 = sin(qJ(4));
t226 = cos(qJ(4));
t235 = (MDP(23) * t226 - MDP(24) * t222) * pkin(3);
t220 = cos(pkin(6));
t223 = sin(qJ(3));
t227 = cos(qJ(3));
t218 = sin(pkin(6));
t224 = sin(qJ(2));
t268 = t218 * t224;
t193 = -t220 * t227 + t223 * t268;
t194 = t220 * t223 + t227 * t268;
t171 = t193 * t226 + t194 * t222;
t172 = -t193 * t222 + t194 * t226;
t288 = t172 * MDP(20) - t171 * MDP(21);
t287 = -MDP(32) * t225 + MDP(33) * t221;
t228 = cos(qJ(2));
t267 = t218 * t228;
t249 = pkin(8) * t267;
t275 = pkin(1) * t224;
t186 = t249 + (pkin(9) + t275) * t220;
t187 = (-pkin(2) * t228 - pkin(9) * t224 - pkin(1)) * t218;
t163 = -t186 * t223 + t227 * t187;
t161 = -pkin(3) * t267 - pkin(10) * t194 + t163;
t164 = t186 * t227 + t187 * t223;
t162 = -pkin(10) * t193 + t164;
t144 = t226 * t161 - t162 * t222;
t145 = t161 * t222 + t162 * t226;
t217 = sin(pkin(12));
t219 = cos(pkin(12));
t157 = -t171 * t217 + t172 * t219;
t147 = t157 * t221 + t225 * t267;
t148 = t157 * t225 - t221 * t267;
t156 = t219 * t171 + t172 * t217;
t256 = t148 * MDP(27);
t286 = (-t147 * t221 + t148 * t225) * MDP(28) + t221 * t256 + t288 + t263 * t156 + t144 * MDP(23) - t145 * MDP(24);
t285 = -(MDP(16) * t223 + MDP(17) * t227) * pkin(9) + t223 * MDP(13) + t227 * MDP(14);
t276 = pkin(9) + pkin(10);
t247 = t276 * t227;
t248 = t276 * t223;
t181 = -t222 * t247 - t226 * t248;
t182 = -t222 * t248 + t226 * t247;
t199 = t222 * t223 - t226 * t227;
t200 = t222 * t227 + t223 * t226;
t284 = t200 * MDP(20) - t199 * MDP(21) + t181 * MDP(23) - t182 * MDP(24);
t283 = 0.2e1 * MDP(16);
t282 = -2 * MDP(19);
t281 = 0.2e1 * MDP(23);
t280 = 0.2e1 * MDP(24);
t279 = 2 * MDP(25);
t278 = 0.2e1 * MDP(32);
t277 = 0.2e1 * MDP(33);
t274 = pkin(1) * t228;
t273 = pkin(3) * t222;
t139 = -pkin(4) * t267 - qJ(5) * t172 + t144;
t143 = -qJ(5) * t171 + t145;
t134 = t139 * t219 - t143 * t217;
t132 = pkin(5) * t267 - t134;
t272 = t132 * t225;
t170 = -t199 * qJ(5) + t182;
t231 = -t200 * qJ(5) + t181;
t153 = t170 * t217 - t219 * t231;
t149 = t153 * t221;
t271 = t153 * t225;
t270 = t156 * t221;
t269 = t156 * t225;
t266 = t220 * t224;
t265 = t221 * t225;
t135 = t217 * t139 + t219 * t143;
t208 = pkin(3) * t226 + pkin(4);
t191 = t217 * t208 + t219 * t273;
t258 = t147 * MDP(28);
t257 = t147 * MDP(30);
t255 = t156 * MDP(29);
t254 = t156 * MDP(31);
t177 = t219 * t199 + t200 * t217;
t253 = t177 * MDP(31);
t252 = t200 * MDP(18);
t251 = t223 * MDP(11);
t250 = MDP(15) + MDP(22);
t209 = -pkin(3) * t227 - pkin(2);
t246 = MDP(28) * t265;
t214 = t221 ^ 2;
t245 = t214 * MDP(27) + MDP(22) + 0.2e1 * t246;
t178 = -t199 * t217 + t200 * t219;
t190 = t208 * t219 - t217 * t273;
t188 = -pkin(5) - t190;
t189 = pkin(11) + t191;
t244 = -t177 * t189 + t178 * t188;
t206 = pkin(4) * t217 + pkin(11);
t207 = -pkin(4) * t219 - pkin(5);
t243 = -t177 * t206 + t178 * t207;
t184 = pkin(4) * t199 + t209;
t242 = t194 * MDP(13) - t193 * MDP(14);
t237 = -t221 * MDP(32) - t225 * MDP(33);
t203 = pkin(8) * t268;
t185 = t203 + (-pkin(2) - t274) * t220;
t234 = (t225 * MDP(29) - t221 * MDP(30)) * t178;
t233 = -0.2e1 * t287;
t215 = t225 ^ 2;
t232 = t284 + ((-t214 + t215) * MDP(28) + MDP(27) * t265) * t178 + t263 * t177;
t175 = pkin(3) * t193 + t185;
t158 = pkin(4) * t171 + t175;
t133 = -pkin(11) * t267 + t135;
t136 = pkin(5) * t156 - pkin(11) * t157 + t158;
t129 = -t133 * t221 + t136 * t225;
t130 = t133 * t225 + t136 * t221;
t230 = t148 * MDP(29) + t129 * MDP(32) - t130 * MDP(33) + t254 - t257;
t213 = t218 ^ 2;
t196 = pkin(1) * t266 + t249;
t195 = t220 * t274 - t203;
t155 = t219 * t170 + t217 * t231;
t152 = pkin(5) * t177 - pkin(11) * t178 + t184;
t142 = t152 * t221 + t155 * t225;
t141 = t152 * t225 - t155 * t221;
t131 = t132 * t221;
t1 = [t220 ^ 2 * MDP(8) + (t134 ^ 2 + t135 ^ 2 + t158 ^ 2) * MDP(26) + MDP(1) + (t194 * MDP(11) - 0.2e1 * t193 * MDP(12)) * t194 + (t172 * MDP(18) + t171 * t282) * t172 + (t254 - 0.2e1 * t257) * t156 + (0.2e1 * t255 + t256 - 0.2e1 * t258) * t148 + ((MDP(4) * t224 + 0.2e1 * MDP(5) * t228) * t224 + t250 * t228 ^ 2) * t213 + 0.2e1 * (MDP(6) * t266 + (MDP(7) * t220 - t242 - t288) * t228) * t218 + 0.2e1 * (t164 * t267 + t185 * t194) * MDP(17) + (-t163 * t267 + t185 * t193) * t283 + (t145 * t267 + t172 * t175) * t280 + (-t144 * t267 + t171 * t175) * t281 + 0.2e1 * (t195 * t220 + t213 * t274) * MDP(9) + 0.2e1 * (-t196 * t220 - t213 * t275) * MDP(10) + (-t130 * t156 + t132 * t148) * t277 + (-t134 * t157 - t135 * t156) * t279 + (t129 * t156 + t132 * t147) * t278; t194 * t251 + t172 * t252 + (-t223 * t193 + t194 * t227) * MDP(12) + (-pkin(2) * t193 - t185 * t227) * MDP(16) + (-pkin(2) * t194 + t185 * t223) * MDP(17) + (t172 * t209 + t175 * t200) * MDP(24) + (t171 * t209 + t175 * t199) * MDP(23) + (-t142 * t156 + t148 * t153) * MDP(33) + t220 * MDP(8) + (t141 * t156 + t147 * t153) * MDP(32) + (-t171 * t200 - t172 * t199) * MDP(19) + t195 * MDP(9) - t196 * MDP(10) + (-t134 * t153 + t135 * t155 + t158 * t184) * MDP(26) + (t153 * t157 - t155 * t156) * MDP(25) + (-t135 * MDP(25) + t230) * t177 + (-t134 * MDP(25) + (-t148 * MDP(28) - t156 * MDP(30) + t132 * MDP(32)) * t221 + (t132 * MDP(33) + t255 + t256 - t258) * t225) * t178 + (MDP(6) * t224 + (MDP(7) - t284 - t285) * t228) * t218; MDP(8) + pkin(2) * t227 * t283 + t209 * t199 * t281 + (t153 ^ 2 + t155 ^ 2 + t184 ^ 2) * MDP(26) + (MDP(27) * t215 - 0.2e1 * t246) * t178 ^ 2 + (0.2e1 * t227 * MDP(12) - 0.2e1 * pkin(2) * MDP(17) + t251) * t223 + (t199 * t282 + t209 * t280 + t252) * t200 + (0.2e1 * t234 + t253) * t177 + (t153 * t178 - t155 * t177) * t279 + (t141 * t177 + t178 * t149) * t278 + (-t142 * t177 + t178 * t271) * t277; t163 * MDP(16) - t164 * MDP(17) + (-t156 * t191 - t157 * t190) * MDP(25) + (t134 * t190 + t135 * t191) * MDP(26) + (t147 * t188 - t189 * t270 - t272) * MDP(32) + (t148 * t188 - t189 * t269 + t131) * MDP(33) + (-t250 - t235) * t267 + t242 + t286; (-t177 * t191 - t178 * t190) * MDP(25) + (-t153 * t190 + t155 * t191) * MDP(26) + (t244 * t221 - t271) * MDP(32) + (t244 * t225 + t149) * MDP(33) + t232 + t285; MDP(15) + (t190 ^ 2 + t191 ^ 2) * MDP(26) - t188 * t233 + 0.2e1 * t235 + t245; -MDP(22) * t267 + (t147 * t207 - t206 * t270 - t272) * MDP(32) + (t148 * t207 - t206 * t269 + t131) * MDP(33) + ((-t156 * t217 - t157 * t219) * MDP(25) + (t134 * t219 + t135 * t217) * MDP(26)) * pkin(4) + t286; (t243 * t221 - t271) * MDP(32) + (t243 * t225 + t149) * MDP(33) + ((-t177 * t217 - t178 * t219) * MDP(25) + (-t153 * t219 + t155 * t217) * MDP(26)) * pkin(4) + t232; (t190 * t219 + t191 * t217) * MDP(26) * pkin(4) + t235 + t245 + t287 * (t188 + t207); (t217 ^ 2 + t219 ^ 2) * MDP(26) * pkin(4) ^ 2 - t207 * t233 + t245; MDP(26) * t158 - t156 * t287; MDP(26) * t184 - t177 * t287; 0; 0; MDP(26); t230; t141 * MDP(32) - t142 * MDP(33) + t234 + t253; t237 * t189 + t263; t237 * t206 + t263; -t287; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
