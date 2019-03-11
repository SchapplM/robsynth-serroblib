% Calculate joint inertia matrix for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP8_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:38:33
% EndTime: 2019-03-09 21:38:38
% DurationCPUTime: 1.62s
% Computational Cost: add. (1931->330), mult. (4150->444), div. (0->0), fcn. (4320->8), ass. (0->108)
t196 = sin(qJ(3));
t261 = 0.2e1 * t196;
t193 = sin(pkin(6));
t197 = sin(qJ(2));
t246 = t193 * t197;
t181 = pkin(8) * t246;
t194 = cos(pkin(6));
t200 = cos(qJ(2));
t254 = pkin(1) * t200;
t161 = t181 + (-pkin(2) - t254) * t194;
t199 = cos(qJ(3));
t168 = -t194 * t199 + t196 * t246;
t169 = t194 * t196 + t199 * t246;
t144 = pkin(3) * t168 - pkin(10) * t169 + t161;
t245 = t193 * t200;
t224 = pkin(8) * t245;
t255 = pkin(1) * t197;
t162 = t224 + (pkin(9) + t255) * t194;
t163 = (-pkin(2) * t200 - pkin(9) * t197 - pkin(1)) * t193;
t148 = t162 * t199 + t163 * t196;
t146 = -pkin(10) * t245 + t148;
t195 = sin(qJ(4));
t198 = cos(qJ(4));
t137 = t198 * t144 - t195 * t146;
t136 = -t168 * pkin(4) - t137;
t138 = t195 * t144 + t198 * t146;
t153 = t169 * t198 - t195 * t245;
t232 = t153 * MDP(20);
t152 = t169 * t195 + t198 * t245;
t234 = t152 * MDP(21);
t260 = t137 * MDP(23) - t138 * MDP(24) - t136 * MDP(25) + t232 - t234;
t201 = pkin(4) + pkin(5);
t247 = qJ(5) * t198;
t259 = t195 * t201 - t247;
t258 = -(pkin(5) + t201) * MDP(29) - MDP(22);
t257 = 2 * MDP(26);
t256 = 0.2e1 * MDP(29);
t252 = pkin(10) - qJ(6);
t251 = pkin(2) * MDP(16);
t250 = pkin(2) * MDP(17);
t249 = pkin(4) * MDP(25);
t248 = qJ(5) * t152;
t244 = t194 * MDP(8);
t243 = t195 * t199;
t242 = t196 * t198;
t241 = t199 * qJ(5);
t147 = -t196 * t162 + t199 * t163;
t176 = -pkin(3) * t199 - pkin(10) * t196 - pkin(2);
t240 = pkin(9) * t243 - t198 * t176;
t160 = t198 * t199 * pkin(9) + t195 * t176;
t190 = t195 ^ 2;
t192 = t198 ^ 2;
t239 = t190 + t192;
t238 = MDP(15) * t200;
t237 = MDP(19) * t198;
t235 = t152 * MDP(19);
t233 = t153 * MDP(18);
t218 = -pkin(4) * t195 + t247;
t164 = (pkin(9) - t218) * t196;
t231 = t164 * MDP(27);
t230 = t168 * MDP(22);
t229 = t169 * MDP(12);
t228 = t169 * MDP(13);
t227 = MDP(24) - MDP(27);
t226 = MDP(25) + MDP(29);
t225 = MDP(26) - MDP(31);
t166 = t168 * qJ(5);
t135 = t166 + t138;
t145 = pkin(3) * t245 - t147;
t223 = 0.2e1 * t166 + t138;
t188 = t199 * pkin(4);
t156 = t188 + t240;
t222 = -0.2e1 * t241 + t160;
t221 = qJ(5) * t195 + pkin(3);
t220 = pkin(9) * MDP(17) - MDP(14);
t219 = -pkin(4) * MDP(28) - MDP(25);
t155 = t160 - t241;
t217 = t135 * t198 + t136 * t195;
t216 = t155 * t198 + t156 * t195;
t215 = -t240 * MDP(23) - t160 * MDP(24);
t214 = -MDP(29) * t195 + MDP(30) * t198;
t213 = qJ(5) * t153 - t145;
t212 = qJ(6) * t153 - t136;
t177 = t252 * t195;
t178 = t252 * t198;
t211 = -t195 * MDP(20) + t177 * MDP(29) - t178 * MDP(30);
t172 = t201 * t198 + t221;
t210 = MDP(29) * t198 + MDP(30) * t195 + MDP(32) * t172;
t175 = -pkin(4) * t198 - t221;
t209 = MDP(23) * pkin(3) - MDP(25) * t175 + MDP(29) * t172 - MDP(31) * t178;
t208 = -MDP(24) * pkin(3) - MDP(27) * t175 + MDP(30) * t172 - MDP(31) * t177;
t206 = t198 * MDP(21) - t211;
t149 = pkin(5) * t199 - qJ(6) * t242 + t156;
t180 = t195 * t196 * qJ(6);
t150 = t155 + t180;
t205 = MDP(25) * t156 - MDP(27) * t155 + MDP(29) * t149 - MDP(30) * t150 - t215;
t203 = qJ(5) ^ 2;
t189 = t193 ^ 2;
t184 = pkin(10) * t243;
t171 = t194 * t255 + t224;
t170 = t194 * t254 - t181;
t154 = (-pkin(9) - t259) * t196;
t151 = t152 * qJ(6);
t139 = pkin(4) * t152 - t213;
t134 = -t201 * t152 + t213;
t133 = t151 + t135;
t132 = -pkin(5) * t168 - t212;
t1 = [(t132 ^ 2 + t133 ^ 2 + t134 ^ 2) * MDP(32) + (t135 ^ 2 + t136 ^ 2 + t139 ^ 2) * MDP(28) + t169 ^ 2 * MDP(11) + t189 * t197 ^ 2 * MDP(4) + MDP(1) + (0.2e1 * MDP(6) * t246 + t244) * t194 + (t233 - 0.2e1 * t235) * t153 + (0.2e1 * (MDP(7) * t194 - t228) * t193 + (0.2e1 * MDP(5) * t197 + t238) * t189) * t200 + (0.2e1 * MDP(14) * t245 - 0.2e1 * t229 + t230 + 0.2e1 * t232 - 0.2e1 * t234) * t168 + 0.2e1 * (-t147 * t245 + t161 * t168) * MDP(16) + 0.2e1 * (t148 * t245 + t161 * t169) * MDP(17) + 0.2e1 * (t170 * t194 + t189 * t254) * MDP(9) + 0.2e1 * (-t171 * t194 - t189 * t255) * MDP(10) + 0.2e1 * (t137 * t168 + t145 * t152) * MDP(23) + 0.2e1 * (t133 * t168 + t134 * t153) * MDP(30) + 0.2e1 * (t135 * t168 - t139 * t153) * MDP(27) + 0.2e1 * (-t136 * t168 + t139 * t152) * MDP(25) + (-t132 * t168 - t134 * t152) * t256 + 0.2e1 * (-t138 * t168 + t145 * t153) * MDP(24) + (-t135 * t152 + t136 * t153) * t257 + 0.2e1 * (-t132 * t153 + t133 * t152) * MDP(31); (t132 * t149 + t133 * t150 + t134 * t154) * MDP(32) + (t135 * t155 + t136 * t156 + t139 * t164) * MDP(28) + t170 * MDP(9) - t171 * MDP(10) + t244 - t169 * t250 + (MDP(6) * t197 + MDP(7) * t200) * t193 + (t156 * MDP(26) + t154 * MDP(30) - MDP(31) * t149 - t231) * t153 + (t164 * MDP(25) - t155 * MDP(26) - t154 * MDP(29) + MDP(31) * t150) * t152 + (-t205 - t251) * t168 + (-t161 * MDP(16) - t135 * MDP(27) + t132 * MDP(29) - t133 * MDP(30) + t220 * t245 + t229 - t230 - t260) * t199 + (-MDP(13) * t245 + t169 * MDP(11) - t168 * MDP(12) + t161 * MDP(17) + (MDP(16) * t245 + t152 * MDP(23) + t153 * MDP(24)) * pkin(9) + (-t153 * MDP(19) - t168 * MDP(21) + t145 * MDP(23) + t139 * MDP(25) - t135 * MDP(26) - t134 * MDP(29) + t133 * MDP(31)) * t195 + (t168 * MDP(20) + t145 * MDP(24) + t136 * MDP(26) - t139 * MDP(27) + t134 * MDP(30) - t132 * MDP(31) + t233 - t235) * t198) * t196; MDP(8) + (t155 ^ 2 + t156 ^ 2 + t164 ^ 2) * MDP(28) + (t149 ^ 2 + t150 ^ 2 + t154 ^ 2) * MDP(32) + ((-t155 * t195 + t156 * t198) * MDP(26) + (-t149 * t198 + t150 * t195) * MDP(31) + (t195 * MDP(25) - t198 * MDP(27)) * t164 + t214 * t154) * t261 + (0.2e1 * t251 + (-t198 * MDP(20) + t195 * MDP(21) + MDP(12)) * t261 + 0.2e1 * t205 + t199 * MDP(22)) * t199 + (-0.2e1 * t250 + (MDP(18) * t192 - 0.2e1 * t195 * t237 + MDP(11) + 0.2e1 * (t195 * MDP(23) + t198 * MDP(24)) * pkin(9)) * t196) * t196; t228 - t193 * t238 + t147 * MDP(16) - t148 * MDP(17) + t195 * t233 + (-t152 * t195 + t153 * t198) * MDP(19) + (-pkin(3) * t152 - t145 * t198) * MDP(23) + (-pkin(3) * t153 + t145 * t195) * MDP(24) + (-t139 * t198 + t152 * t175) * MDP(25) + t217 * MDP(26) + (-t139 * t195 - t153 * t175) * MDP(27) + t139 * t175 * MDP(28) + (t134 * t198 - t152 * t172) * MDP(29) + (t134 * t195 + t153 * t172) * MDP(30) + (-t132 * t195 - t133 * t198 + t152 * t178 - t153 * t177) * MDP(31) + (t132 * t177 + t133 * t178 + t134 * t172) * MDP(32) + ((-t152 * t198 + t153 * t195) * MDP(26) + t217 * MDP(28)) * pkin(10) + (-MDP(14) + (-t227 * t198 + (-MDP(23) - MDP(25)) * t195) * pkin(10) + t206) * t168; t184 * MDP(23) + (-t164 * t198 + t184) * MDP(25) + t216 * MDP(26) - t195 * t231 + (t216 * pkin(10) + t164 * t175) * MDP(28) + (-t149 * t195 - t150 * t198) * MDP(31) + (t149 * t177 + t150 * t178) * MDP(32) + t210 * t154 + ((t227 * pkin(10) - MDP(21)) * t198 + t211 - t220) * t199 + (MDP(13) - pkin(9) * MDP(16) + (-t190 + t192) * MDP(19) + (-MDP(23) * pkin(9) + t208) * t198 + (MDP(18) * t198 + MDP(24) * pkin(9) - t209) * t195) * t196; MDP(15) + t190 * MDP(18) + (t239 * pkin(10) ^ 2 + t175 ^ 2) * MDP(28) + (t172 ^ 2 + t177 ^ 2 + t178 ^ 2) * MDP(32) + t239 * pkin(10) * t257 + 0.2e1 * t209 * t198 + 0.2e1 * (t208 + t237) * t195; (-pkin(4) * t153 - t248) * MDP(26) + t223 * MDP(27) + (-pkin(4) * t136 + qJ(5) * t135) * MDP(28) + t212 * MDP(29) + (t151 + t223) * MDP(30) + (t153 * t201 + t248) * MDP(31) + (qJ(5) * t133 - t132 * t201) * MDP(32) + (t249 - t258) * t168 + t260; (-0.2e1 * t188 - t240) * MDP(25) + t222 * MDP(27) + (-pkin(4) * t156 + qJ(5) * t155) * MDP(28) - t156 * MDP(29) + (t180 + t222) * MDP(30) + (qJ(5) * t150 - t149 * t201) * MDP(32) + t258 * t199 + ((-t225 * qJ(5) - MDP(21)) * t195 + (-MDP(26) * pkin(4) + MDP(29) * qJ(6) + MDP(31) * t201 + MDP(20)) * t198) * t196 + t215; t218 * MDP(26) + t259 * MDP(31) + (qJ(5) * t178 - t177 * t201) * MDP(32) + ((MDP(28) * qJ(5) - t227) * t198 + (-MDP(23) + t219) * t195) * pkin(10) + t206; MDP(22) + 0.2e1 * t249 + (pkin(4) ^ 2 + t203) * MDP(28) + t201 * t256 + (t201 ^ 2 + t203) * MDP(32) + 0.2e1 * (MDP(27) + MDP(30)) * qJ(5); MDP(28) * t136 + MDP(32) * t132 + t225 * t153 - t226 * t168; MDP(28) * t156 + MDP(32) * t149 + t226 * t199 + t225 * t242; MDP(32) * t177 + (MDP(28) * pkin(10) + t225) * t195; -MDP(32) * t201 - MDP(29) + t219; MDP(28) + MDP(32); -t152 * MDP(29) + t153 * MDP(30) + MDP(32) * t134; MDP(32) * t154 + t214 * t196; t210; 0; 0; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
