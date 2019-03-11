% Calculate joint inertia matrix for
% S6RRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRPRR6_inertiaJ_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:29:33
% EndTime: 2019-03-09 18:29:36
% DurationCPUTime: 1.01s
% Computational Cost: add. (1842->198), mult. (3729->286), div. (0->0), fcn. (4203->10), ass. (0->100)
t198 = sin(qJ(3));
t202 = cos(qJ(3));
t235 = -qJ(4) - pkin(8);
t182 = t235 * t198;
t183 = t235 * t202;
t194 = sin(pkin(11));
t195 = cos(pkin(11));
t161 = t195 * t182 + t194 * t183;
t175 = t194 * t202 + t195 * t198;
t148 = -t175 * pkin(9) + t161;
t162 = t194 * t182 - t195 * t183;
t174 = -t194 * t198 + t195 * t202;
t149 = t174 * pkin(9) + t162;
t197 = sin(qJ(5));
t201 = cos(qJ(5));
t129 = t201 * t148 - t197 * t149;
t130 = t197 * t148 + t201 * t149;
t155 = -t201 * t174 + t197 * t175;
t156 = t197 * t174 + t201 * t175;
t119 = -t156 * pkin(10) + t129;
t120 = -t155 * pkin(10) + t130;
t196 = sin(qJ(6));
t200 = cos(qJ(6));
t135 = t200 * t155 + t196 * t156;
t136 = -t196 * t155 + t200 * t156;
t213 = t136 * MDP(29) - t135 * MDP(30) + (t200 * t119 - t196 * t120) * MDP(32) - (t196 * t119 + t200 * t120) * MDP(33);
t206 = t156 * MDP(22) - t155 * MDP(23) + t129 * MDP(25) - t130 * MDP(26) + t213;
t251 = t206 - (t198 * MDP(16) + t202 * MDP(17)) * pkin(8) + t198 * MDP(13) + t202 * MDP(14);
t208 = (t200 * MDP(32) - t196 * MDP(33)) * pkin(5);
t199 = sin(qJ(2));
t168 = t175 * t199;
t232 = t199 * t202;
t234 = t198 * t199;
t169 = -t194 * t234 + t195 * t232;
t144 = t201 * t168 + t197 * t169;
t145 = -t197 * t168 + t201 * t169;
t123 = t200 * t144 + t196 * t145;
t124 = -t196 * t144 + t200 * t145;
t230 = t124 * MDP(29) - t123 * MDP(30);
t250 = t145 * MDP(22) - t144 * MDP(23);
t187 = t195 * pkin(3) + pkin(4);
t239 = pkin(3) * t194;
t171 = t201 * t187 - t197 * t239;
t170 = pkin(5) + t171;
t172 = t197 * t187 + t201 * t239;
t231 = t200 * t172;
t221 = (t196 * t170 + t231) * MDP(33);
t146 = t200 * t170 - t196 * t172;
t222 = t146 * MDP(32);
t249 = t222 - t221;
t203 = cos(qJ(2));
t181 = -t203 * pkin(2) - t199 * pkin(8) - pkin(1);
t176 = t202 * t181;
t238 = pkin(7) * t198;
t157 = -qJ(4) * t232 + t176 + (-pkin(3) - t238) * t203;
t236 = pkin(7) * t203;
t216 = t202 * t236;
t160 = t216 + (-qJ(4) * t199 + t181) * t198;
t137 = t195 * t157 - t194 * t160;
t128 = -t203 * pkin(4) - t169 * pkin(9) + t137;
t138 = t194 * t157 + t195 * t160;
t131 = -t168 * pkin(9) + t138;
t117 = t201 * t128 - t197 * t131;
t118 = t197 * t128 + t201 * t131;
t115 = -t203 * pkin(5) - t145 * pkin(10) + t117;
t116 = -t144 * pkin(10) + t118;
t108 = t200 * t115 - t196 * t116;
t109 = t196 * t115 + t200 * t116;
t247 = t108 * MDP(32) - t109 * MDP(33) + t230;
t248 = t117 * MDP(25) - t118 * MDP(26) + t247 + t250;
t244 = 2 * MDP(18);
t243 = -2 * MDP(21);
t242 = 0.2e1 * MDP(26);
t241 = -2 * MDP(28);
t240 = 0.2e1 * MDP(33);
t237 = pkin(7) * t202;
t233 = t198 * t202;
t180 = pkin(3) * t234 + t199 * pkin(7);
t225 = t124 * MDP(27);
t224 = t135 * MDP(32);
t223 = t145 * MDP(20);
t220 = t155 * MDP(25);
t219 = t171 * MDP(25);
t218 = t172 * MDP(26);
t217 = MDP(24) + MDP(31);
t188 = -t202 * pkin(3) - pkin(2);
t215 = MDP(12) * t233;
t214 = MDP(15) + t217;
t158 = t168 * pkin(4) + t180;
t167 = -t174 * pkin(4) + t188;
t212 = t202 * MDP(13) - t198 * MDP(14);
t207 = t217 - t218 + t219;
t192 = t202 ^ 2;
t191 = t199 ^ 2;
t190 = t198 ^ 2;
t166 = t198 * t181 + t216;
t165 = -t198 * t236 + t176;
t139 = t155 * pkin(5) + t167;
t132 = t144 * pkin(5) + t158;
t1 = [(t137 ^ 2 + t138 ^ 2 + t180 ^ 2) * MDP(19) - 0.2e1 * pkin(1) * t199 * MDP(10) + MDP(1) + (t144 * t243 + t223) * t145 + (t123 * t241 + t225) * t124 + t214 * t203 ^ 2 + (t192 * MDP(11) + MDP(4) - 0.2e1 * t215) * t191 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t212) * t199 - t250 - t230) * t203 + 0.2e1 * (-t165 * t203 + t191 * t238) * MDP(16) + 0.2e1 * (t166 * t203 + t191 * t237) * MDP(17) + (t118 * t203 + t158 * t145) * t242 + 0.2e1 * (-t117 * t203 + t158 * t144) * MDP(25) + 0.2e1 * (-t108 * t203 + t132 * t123) * MDP(32) + (t109 * t203 + t132 * t124) * t240 + (-t137 * t169 - t138 * t168) * t244; (-t137 * t175 + t138 * t174 - t161 * t169 - t162 * t168) * MDP(18) + (t137 * t161 + t138 * t162 + t180 * t188) * MDP(19) + (-t156 * t144 - t145 * t155) * MDP(21) + (t167 * t144 + t158 * t155) * MDP(25) + (t167 * t145 + t158 * t156) * MDP(26) + (-t136 * t123 - t124 * t135) * MDP(28) + (t139 * t123 + t132 * t135) * MDP(32) + (t139 * t124 + t132 * t136) * MDP(33) + t156 * t223 + t136 * t225 + (MDP(6) + (-t190 + t192) * MDP(12) + (-pkin(2) * t198 - t237) * MDP(16) + (-pkin(2) * t202 + t238) * MDP(17) - pkin(7) * MDP(9) + MDP(11) * t233) * t199 + (-pkin(7) * MDP(10) + MDP(7) - t251) * t203; MDP(8) + t190 * MDP(11) + 0.2e1 * t215 + (-t161 * t175 + t162 * t174) * t244 + (t161 ^ 2 + t162 ^ 2 + t188 ^ 2) * MDP(19) + 0.2e1 * t167 * t220 + 0.2e1 * t139 * t224 + 0.2e1 * (t202 * MDP(16) - t198 * MDP(17)) * pkin(2) + (MDP(20) * t156 + t155 * t243 + t167 * t242) * t156 + (MDP(27) * t136 + t135 * t241 + t139 * t240) * t136; t165 * MDP(16) - t166 * MDP(17) + t212 * t199 + (-MDP(15) - t207 - t249) * t203 + ((-t168 * t194 - t169 * t195) * MDP(18) + (t137 * t195 + t138 * t194) * MDP(19)) * pkin(3) + t248; ((t174 * t194 - t175 * t195) * MDP(18) + (t161 * t195 + t162 * t194) * MDP(19)) * pkin(3) + t251; (t194 ^ 2 + t195 ^ 2) * MDP(19) * pkin(3) ^ 2 + 0.2e1 * t219 - 0.2e1 * t218 + 0.2e1 * t222 - 0.2e1 * t221 + t214; t180 * MDP(19) + t144 * MDP(25) + t145 * MDP(26) + t123 * MDP(32) + t124 * MDP(33); t188 * MDP(19) + t156 * MDP(26) + t136 * MDP(33) + t220 + t224; 0; MDP(19); (-t217 - t208) * t203 + t248; t206; (t200 * pkin(5) + t146) * MDP(32) + (-t231 + (-pkin(5) - t170) * t196) * MDP(33) + t207; 0; 0.2e1 * t208 + t217; -t203 * MDP(31) + t247; t213; MDP(31) + t249; 0; MDP(31) + t208; MDP(31);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
