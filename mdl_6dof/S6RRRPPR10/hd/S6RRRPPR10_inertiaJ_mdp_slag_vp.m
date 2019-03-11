% Calculate joint inertia matrix for
% S6RRRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPPR10_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:28:02
% EndTime: 2019-03-09 16:28:06
% DurationCPUTime: 1.52s
% Computational Cost: add. (1668->260), mult. (3592->374), div. (0->0), fcn. (3826->10), ass. (0->112)
t210 = cos(qJ(3));
t207 = sin(qJ(3));
t234 = -qJ(4) * t207 - pkin(2);
t256 = pkin(3) + qJ(5);
t176 = -t210 * t256 + t234;
t265 = pkin(4) + pkin(9);
t185 = t265 * t207;
t203 = cos(pkin(11));
t181 = t203 * t185;
t201 = sin(pkin(11));
t151 = pkin(5) * t207 + t181 + (pkin(10) * t210 - t176) * t201;
t155 = t203 * t176 + t201 * t185;
t249 = t203 * t210;
t152 = -pkin(10) * t249 + t155;
t209 = cos(qJ(6));
t206 = sin(qJ(6));
t252 = t201 * t206;
t167 = -t209 * t249 + t210 * t252;
t178 = t201 * t209 + t203 * t206;
t168 = t178 * t210;
t267 = (t151 * t209 - t152 * t206) * MDP(31) - (t151 * t206 + t152 * t209) * MDP(32) - t168 * MDP(28) + t167 * MDP(29);
t266 = 2 * MDP(24);
t154 = -t176 * t201 + t181;
t184 = -pkin(3) * t210 + t234;
t264 = t210 * MDP(12) - pkin(2) * MDP(17) - t184 * MDP(20) + t154 * MDP(22) - t155 * MDP(23) + t267;
t204 = cos(pkin(6));
t202 = sin(pkin(6));
t208 = sin(qJ(2));
t251 = t202 * t208;
t173 = t204 * t207 + t210 * t251;
t263 = 0.2e1 * t173;
t211 = cos(qJ(2));
t250 = t202 * t211;
t260 = -2 * MDP(27);
t259 = 0.2e1 * MDP(32);
t258 = pkin(1) * t208;
t257 = pkin(1) * t211;
t255 = -pkin(10) - t256;
t254 = pkin(9) * MDP(18);
t253 = pkin(9) * MDP(21);
t248 = t204 * MDP(8);
t247 = t208 * MDP(6);
t236 = pkin(8) * t250;
t164 = t236 + (pkin(9) + t258) * t204;
t165 = (-pkin(2) * t211 - pkin(9) * t208 - pkin(1)) * t202;
t149 = -t207 * t164 + t165 * t210;
t190 = pkin(3) * t250;
t148 = -t149 + t190;
t139 = pkin(4) * t173 + qJ(5) * t250 + t148;
t172 = -t204 * t210 + t207 * t251;
t189 = pkin(8) * t251;
t163 = t189 + (-pkin(2) - t257) * t204;
t221 = -qJ(4) * t173 + t163;
t141 = t172 * t256 + t221;
t133 = t201 * t139 + t203 * t141;
t150 = t210 * t164 + t207 * t165;
t186 = t265 * t210;
t188 = t201 ^ 2 + t203 ^ 2;
t246 = MDP(15) * t211;
t142 = t154 * t203 + t155 * t201;
t245 = MDP(25) * t142;
t244 = MDP(31) * t178;
t235 = qJ(4) * t250;
t147 = t235 - t150;
t143 = -pkin(4) * t172 - t147;
t243 = t143 * MDP(25);
t158 = t172 * t203 + t201 * t250;
t159 = -t172 * t201 + t203 * t250;
t144 = -t209 * t158 - t159 * t206;
t242 = t144 * MDP(29);
t145 = t158 * t206 - t159 * t209;
t241 = t145 * MDP(26);
t240 = t168 * MDP(26);
t239 = t172 * MDP(12);
t238 = MDP(11) + MDP(30);
t237 = MDP(17) - MDP(20);
t233 = -pkin(3) * MDP(21) + MDP(19);
t132 = t203 * t139 - t141 * t201;
t232 = pkin(2) * MDP(16) + t184 * MDP(19);
t231 = t132 * t203 + t133 * t201;
t230 = t203 * MDP(22) - t201 * MDP(23);
t229 = MDP(22) * t201 + MDP(23) * t203;
t226 = -t167 * MDP(31) - t168 * MDP(32);
t179 = t203 * t209 - t252;
t225 = MDP(31) * t179 - t178 * MDP(32);
t224 = MDP(18) + t230;
t223 = qJ(4) * MDP(25) + t229;
t222 = -t158 * MDP(22) - t159 * MDP(23) + t243;
t130 = pkin(5) * t173 + pkin(10) * t159 + t132;
t131 = pkin(10) * t158 + t133;
t128 = t130 * t209 - t131 * t206;
t129 = t130 * t206 + t131 * t209;
t219 = t145 * MDP(28) + t128 * MDP(31) - t129 * MDP(32) - t242;
t182 = t255 * t201;
t183 = t255 * t203;
t218 = t179 * MDP(28) - t178 * MDP(29) + (-t182 * t206 + t183 * t209) * MDP(31) - (t182 * t209 + t183 * t206) * MDP(32);
t217 = t224 + t225;
t215 = (t158 * t201 + t159 * t203) * MDP(24) + t231 * MDP(25);
t214 = -pkin(3) * MDP(18) - t230 * t256 + MDP(13) + t218;
t213 = pkin(9) ^ 2;
t212 = qJ(4) ^ 2;
t200 = t210 ^ 2;
t199 = t207 ^ 2;
t196 = t202 ^ 2;
t191 = pkin(5) * t201 + qJ(4);
t177 = t188 * t256;
t175 = t204 * t258 + t236;
t174 = t204 * t257 - t189;
t171 = pkin(5) * t249 + t186;
t146 = pkin(3) * t172 + t221;
t134 = -pkin(5) * t158 + t143;
t1 = [(t146 ^ 2 + t147 ^ 2 + t148 ^ 2) * MDP(21) + (t132 ^ 2 + t133 ^ 2 + t143 ^ 2) * MDP(25) + t196 * t208 ^ 2 * MDP(4) + MDP(1) + (0.2e1 * t202 * t247 + t248) * t204 + t238 * t173 ^ 2 + (MDP(28) * t263 + t144 * t260 + t241) * t145 + (0.2e1 * MDP(5) * t208 + t246) * t196 * t211 + 0.2e1 * (t150 * t250 + t163 * t173) * MDP(17) + 0.2e1 * (-t146 * t173 + t147 * t250) * MDP(20) + 0.2e1 * (-t149 * t250 + t163 * t172) * MDP(16) + 0.2e1 * (-t146 * t172 - t148 * t250) * MDP(19) + 0.2e1 * (t174 * t204 + t196 * t257) * MDP(9) + 0.2e1 * (-t175 * t204 - t196 * t258) * MDP(10) + 0.2e1 * (t147 * t172 + t148 * t173) * MDP(18) + 0.2e1 * (t132 * t173 - t143 * t158) * MDP(22) + 0.2e1 * (t128 * t173 + t134 * t144) * MDP(31) + 0.2e1 * (-t133 * t173 - t143 * t159) * MDP(23) + (-t129 * t173 + t134 * t145) * t259 + (t132 * t159 + t133 * t158) * t266 + (-t239 - t242) * t263 + 0.2e1 * (-t173 * MDP(13) + t172 * MDP(14) + MDP(7) * t204) * t250; (t132 * t154 + t133 * t155) * MDP(25) + t174 * MDP(9) - t175 * MDP(10) + t146 * t184 * MDP(21) + (t154 * t159 + t155 * t158) * MDP(24) + (t144 * t168 + t145 * t167) * MDP(27) - t145 * t240 + t248 + (-t134 * t167 + t144 * t171) * MDP(31) + (-t134 * t168 + t145 * t171) * MDP(32) + t222 * t186 - t232 * t172 + (-t147 * t253 + (-pkin(9) * t172 - t147) * MDP(18) + (t132 * t201 - t133 * t203) * MDP(24) + t146 * MDP(19) - t163 * MDP(16) + t230 * t143) * t210 + t264 * t173 + (-t239 + t163 * MDP(17) - t146 * MDP(20) + t132 * MDP(22) - t133 * MDP(23) + (MDP(18) + t253) * t148 + (t238 + t254) * t173 + t219) * t207 + (t247 + (-t207 * MDP(13) - t210 * MDP(14) + MDP(7) + (t237 * t210 + (MDP(16) - MDP(19)) * t207) * pkin(9)) * t211) * t202; MDP(8) + (t184 ^ 2 + t200 * t213) * MDP(21) + (t154 ^ 2 + t155 ^ 2 + t186 ^ 2) * MDP(25) + (MDP(21) * t213 + t238) * t199 - (0.2e1 * t167 * MDP(27) - t240) * t168 + 0.2e1 * t226 * t171 + 0.2e1 * (t199 + t200) * t254 + 0.2e1 * t264 * t207 + 0.2e1 * ((t154 * t201 - t155 * t203) * MDP(24) + t230 * t186 + t232) * t210; -t202 * t246 + t149 * MDP(16) - t150 * MDP(17) + (-t149 + 0.2e1 * t190) * MDP(19) + (-0.2e1 * t235 + t150) * MDP(20) + (-pkin(3) * t148 - qJ(4) * t147) * MDP(21) + (-qJ(4) * t158 + t143 * t201) * MDP(22) + (-qJ(4) * t159 + t143 * t203) * MDP(23) - t231 * MDP(24) + qJ(4) * t243 + t179 * t241 + (-t144 * t179 - t145 * t178) * MDP(27) + (t134 * t178 + t144 * t191) * MDP(31) + (t134 * t179 + t145 * t191) * MDP(32) - t215 * t256 + (-MDP(18) * qJ(4) - MDP(14)) * t172 + t214 * t173; -t142 * MDP(24) - t179 * t240 + (t167 * t179 + t168 * t178) * MDP(27) + (-t167 * t191 + t171 * t178) * MDP(31) + (-t168 * t191 + t171 * t179) * MDP(32) - t256 * t245 + t223 * t186 + (t224 * qJ(4) + MDP(14)) * t210 + t214 * t207 + ((MDP(21) * qJ(4) - t237) * t210 + (-MDP(16) + t233) * t207) * pkin(9); MDP(15) - 0.2e1 * pkin(3) * MDP(19) + (pkin(3) ^ 2 + t212) * MDP(21) + t177 * t266 + (t188 * t256 ^ 2 + t212) * MDP(25) + 0.2e1 * t191 * t244 + (MDP(26) * t179 + t178 * t260 + t191 * t259) * t179 + 0.2e1 * (MDP(20) + t229) * qJ(4); -MDP(19) * t250 + t148 * MDP(21) + t173 * t217 + t215; t245 + (t217 + t253) * t207; -MDP(24) * t188 - MDP(25) * t177 + t233; MDP(25) * t188 + MDP(21); t144 * MDP(31) + t145 * MDP(32) + t222; t186 * MDP(25) + t210 * t230 + t226; MDP(32) * t179 + t223 + t244; 0; MDP(25); t173 * MDP(30) + t219; t207 * MDP(30) + t267; t218; t225; 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
