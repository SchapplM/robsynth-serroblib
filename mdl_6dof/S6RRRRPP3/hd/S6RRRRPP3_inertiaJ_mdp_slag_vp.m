% Calculate joint inertia matrix for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP3_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:56:33
% EndTime: 2019-03-09 20:56:36
% DurationCPUTime: 0.97s
% Computational Cost: add. (1251->226), mult. (2152->280), div. (0->0), fcn. (2181->6), ass. (0->95)
t196 = sin(qJ(3));
t178 = pkin(2) * t196 + pkin(9);
t195 = sin(qJ(4));
t191 = t195 ^ 2;
t198 = cos(qJ(4));
t192 = t198 ^ 2;
t232 = t191 + t192;
t235 = t232 * t178;
t258 = t195 * MDP(20) + t198 * MDP(21);
t197 = sin(qJ(2));
t245 = cos(qJ(3));
t246 = cos(qJ(2));
t162 = t196 * t197 - t245 * t246;
t257 = 0.2e1 * t162;
t225 = 0.2e1 * t195;
t224 = 0.2e1 * t198;
t256 = pkin(8) + pkin(7);
t163 = t196 * t246 + t245 * t197;
t180 = -t246 * pkin(2) - pkin(1);
t138 = t162 * pkin(3) - t163 * pkin(9) + t180;
t170 = t256 * t197;
t172 = t256 * t246;
t145 = -t196 * t170 + t245 * t172;
t132 = t195 * t138 + t198 * t145;
t157 = t162 * qJ(5);
t128 = -t157 - t132;
t131 = t198 * t138 - t195 * t145;
t255 = t131 * MDP(23) - t132 * MDP(24) - t128 * MDP(27);
t144 = t245 * t170 + t196 * t172;
t243 = t163 * t195;
t210 = pkin(4) * t243 + t144;
t242 = t163 * t198;
t134 = -qJ(5) * t242 + t210;
t254 = t134 * MDP(26);
t253 = t144 * MDP(24);
t252 = MDP(26) - MDP(31);
t193 = pkin(4) + qJ(6);
t241 = t195 * qJ(5);
t251 = -t193 * t198 - t241;
t250 = -0.2e1 * pkin(4) * MDP(26) + 0.2e1 * t193 * MDP(31) + MDP(22);
t249 = 0.2e1 * t180;
t248 = 2 * MDP(25);
t247 = 2 * MDP(29);
t223 = t245 * pkin(2);
t179 = -t223 - pkin(3);
t244 = pkin(3) - t179;
t187 = qJ(5) * t198;
t240 = t195 * t198;
t153 = -pkin(3) + t251;
t147 = -t223 + t153;
t239 = -t147 - t153;
t158 = (pkin(5) + t178) * t195;
t190 = t198 * pkin(5);
t159 = t178 * t198 + t190;
t238 = t158 * t195 + t159 * t198;
t217 = -t198 * pkin(4) - t241;
t168 = -pkin(3) + t217;
t156 = -t223 + t168;
t237 = t156 + t168;
t169 = (pkin(5) + pkin(9)) * t195;
t171 = pkin(9) * t198 + t190;
t236 = t169 * t195 + t171 * t198;
t182 = t195 * MDP(29);
t234 = t195 * MDP(25) + t182;
t233 = t232 * pkin(9);
t231 = MDP(28) * t134;
t230 = MDP(28) * t195;
t129 = -pkin(4) * t162 - t131;
t228 = t129 * MDP(28);
t227 = 0.2e1 * t246;
t226 = -MDP(24) + MDP(27);
t222 = MDP(19) * t240;
t221 = t191 * MDP(18) + MDP(15) + 0.2e1 * t222;
t219 = -pkin(4) * MDP(28) + MDP(26);
t218 = (-t193 * t195 + t187) * MDP(29) + (-pkin(4) * t195 + t187) * MDP(25) + t258;
t216 = t198 * MDP(20) - t195 * MDP(21);
t130 = (qJ(6) * t195 - t187) * t163 + t210;
t214 = -t144 * MDP(23) - t130 * MDP(31);
t213 = -t134 * MDP(27) - t130 * MDP(30);
t212 = t159 * MDP(30) - t158 * MDP(31);
t211 = t171 * MDP(30) - t169 * MDP(31);
t209 = -t128 * MDP(28) + t226 * t162;
t208 = t228 + (-MDP(23) + MDP(26)) * t162;
t207 = (t245 * MDP(16) - t196 * MDP(17)) * pkin(2);
t206 = MDP(23) * pkin(3) + MDP(26) * t168 - MDP(31) * t153;
t205 = -MDP(24) * pkin(3) - MDP(27) * t168 - MDP(30) * t153;
t204 = -MDP(23) * t179 + MDP(26) * t156 - MDP(31) * t147;
t203 = MDP(24) * t179 - MDP(27) * t156 - MDP(30) * t147;
t202 = (MDP(28) * qJ(5) + t226) * t198 + (-MDP(23) + t219) * t195;
t125 = pkin(5) * t242 - t193 * t162 - t131;
t127 = -pkin(5) * t243 - t128;
t201 = (-t128 * t198 + t129 * t195) * MDP(25) + t198 * t254 + t195 * t253 + (t125 * t195 + t127 * t198) * MDP(29) - t144 * MDP(16) - t145 * MDP(17) + ((-t191 + t192) * MDP(19) + MDP(18) * t240 + MDP(13)) * t163 + (-MDP(14) + t258) * t162;
t199 = qJ(5) ^ 2;
t185 = t198 * MDP(29);
t1 = [pkin(1) * MDP(9) * t227 + (t128 ^ 2 + t129 ^ 2 + t134 ^ 2) * MDP(28) + (t125 ^ 2 + t127 ^ 2 + t130 ^ 2) * MDP(32) + MDP(1) + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t197 + MDP(5) * t227) * t197 + (MDP(16) * t249 + t162 * MDP(22)) * t162 + (t129 * MDP(26) + t127 * MDP(30) - t125 * MDP(31) + t255) * t257 + (MDP(17) * t249 + (-MDP(12) + t216) * t257 + (t129 * MDP(25) + t125 * MDP(29) + t213 + t253) * t224 + (t128 * MDP(25) - t127 * MDP(29) - t214 - t254) * t225 + (t192 * MDP(18) + MDP(11) - 0.2e1 * t222) * t163) * t163; (t125 * t158 + t127 * t159 + t130 * t147) * MDP(32) + (-t246 * MDP(10) - t197 * MDP(9)) * pkin(7) + t156 * t231 + t212 * t162 + (t209 * t178 + (MDP(29) * t158 + t203) * t163 + t214) * t198 + (t208 * t178 + (-MDP(29) * t159 - t204) * t163 + t213) * t195 + t197 * MDP(6) + t201 + t246 * MDP(7); MDP(8) + (t232 * t178 ^ 2 + t156 ^ 2) * MDP(28) + (t147 ^ 2 + t158 ^ 2 + t159 ^ 2) * MDP(32) + 0.2e1 * t207 + t235 * t248 + t238 * t247 + t204 * t224 + t203 * t225 + t221; (t125 * t169 + t127 * t171 + t130 * t153) * MDP(32) + t201 + t168 * t231 + t211 * t162 + (t209 * pkin(9) + (MDP(29) * t169 + t205) * t163 + t214) * t198 + (t208 * pkin(9) + (-MDP(29) * t171 - t206) * t163 + t213) * t195; (t233 + t235) * MDP(25) + (pkin(9) * t235 + t156 * t168) * MDP(28) + (t236 + t238) * MDP(29) + (t147 * t153 + t158 * t169 + t159 * t171) * MDP(32) + t207 + (t244 * MDP(23) + t237 * MDP(26) + t239 * MDP(31)) * t198 + (-t244 * MDP(24) - t237 * MDP(27) + t239 * MDP(30)) * t195 + t221; (t232 * pkin(9) ^ 2 + t168 ^ 2) * MDP(28) + (t153 ^ 2 + t169 ^ 2 + t171 ^ 2) * MDP(32) + t233 * t248 + t236 * t247 + t206 * t224 + t205 * t225 + t221; (-pkin(4) * t129 - qJ(5) * t128) * MDP(28) + (0.2e1 * t157 + t132) * MDP(30) + (qJ(5) * t127 - t125 * t193) * MDP(32) + (qJ(5) * MDP(27) + t250) * t162 + (t217 * MDP(25) + t251 * MDP(29) + (-MDP(30) * t195 - MDP(31) * t198) * pkin(5) + t216) * t163 - t252 * t131 + t255; (qJ(5) * t159 - t158 * t193) * MDP(32) + t202 * t178 + t212 + t218; (qJ(5) * t171 - t169 * t193) * MDP(32) + t202 * pkin(9) + t211 + t218; (pkin(4) ^ 2 + t199) * MDP(28) + (t193 ^ 2 + t199) * MDP(32) + 0.2e1 * (MDP(27) + MDP(30)) * qJ(5) + t250; t228 + MDP(32) * t125 + (MDP(25) + MDP(29)) * t242 + t252 * t162; MDP(32) * t158 + t178 * t230 + t234; MDP(32) * t169 + pkin(9) * t230 + t234; -MDP(32) * t193 - MDP(31) + t219; MDP(28) + MDP(32); t162 * MDP(30) + MDP(32) * t127 - t163 * t182; MDP(32) * t159 + t185; MDP(32) * t171 + t185; MDP(32) * qJ(5) + MDP(30); 0; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
