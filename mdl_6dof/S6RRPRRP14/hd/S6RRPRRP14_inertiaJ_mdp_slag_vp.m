% Calculate joint inertia matrix for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP14_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP14_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP14_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:09:54
% EndTime: 2019-03-09 13:09:58
% DurationCPUTime: 1.37s
% Computational Cost: add. (1417->262), mult. (2986->358), div. (0->0), fcn. (3018->8), ass. (0->99)
t171 = sin(pkin(6));
t178 = cos(qJ(2));
t225 = t171 * t178;
t175 = sin(qJ(2));
t226 = t171 * t175;
t160 = pkin(8) * t226;
t172 = cos(pkin(6));
t235 = pkin(1) * t178;
t205 = -pkin(2) - t235;
t134 = pkin(3) * t226 + t160 + (-pkin(9) + t205) * t172;
t179 = -pkin(2) - pkin(9);
t204 = -qJ(3) * t175 - pkin(1);
t140 = (t179 * t178 + t204) * t171;
t174 = sin(qJ(4));
t177 = cos(qJ(4));
t129 = t134 * t177 - t174 * t140;
t130 = t134 * t174 + t140 * t177;
t149 = t172 * t177 - t174 * t225;
t243 = t149 * MDP(17) + t129 * MDP(20) - t130 * MDP(21);
t242 = 0.2e1 * t172;
t156 = pkin(4) * t174 - pkin(10) * t177 + qJ(3);
t173 = sin(qJ(5));
t176 = cos(qJ(5));
t223 = t176 * t179;
t143 = t173 * t156 + t174 * t223;
t138 = qJ(6) * t174 + t143;
t153 = t176 * t156;
t224 = t173 * t179;
t141 = -t153 + (-pkin(5) + t224) * t174;
t189 = (-t174 * t224 + t153) * MDP(27) - t143 * MDP(28);
t241 = qJ(3) * MDP(20) - t141 * MDP(29) + t138 * MDP(31) + t189;
t233 = MDP(32) * pkin(10);
t240 = MDP(30) + t233;
t206 = -MDP(28) + MDP(31);
t207 = MDP(27) + MDP(29);
t186 = -t207 * t173 + t206 * t176;
t238 = 0.2e1 * t177;
t237 = 0.2e1 * MDP(29);
t236 = 0.2e1 * MDP(31);
t148 = t172 * t174 + t177 * t225;
t234 = pkin(5) * t148;
t232 = pkin(2) * MDP(14);
t231 = qJ(6) * t148;
t127 = -pkin(4) * t226 - t129;
t230 = t127 * t173;
t229 = t127 * t176;
t136 = t149 * t173 - t176 * t226;
t228 = t136 * t176;
t137 = t149 * t176 + t173 * t226;
t227 = t137 * t173;
t128 = pkin(10) * t226 + t130;
t151 = t172 * t175 * pkin(1) + pkin(8) * t225;
t164 = t172 * qJ(3);
t145 = -t164 - t151;
t139 = pkin(3) * t225 - t145;
t133 = pkin(4) * t148 - pkin(10) * t149 + t139;
t124 = t176 * t128 + t173 * t133;
t166 = t173 ^ 2;
t169 = t176 ^ 2;
t222 = t166 + t169;
t220 = MDP(16) * t149;
t125 = pkin(5) * t136 - qJ(6) * t137 + t127;
t218 = t125 * MDP(32);
t217 = t136 * MDP(25);
t216 = t137 * MDP(22);
t215 = t137 * MDP(24);
t214 = t148 * MDP(26);
t212 = t149 * MDP(21);
t197 = -pkin(5) * t176 - qJ(6) * t173;
t157 = -pkin(4) + t197;
t211 = t157 * MDP(32);
t210 = t176 * MDP(22);
t209 = t176 * MDP(23);
t208 = t179 * MDP(21);
t203 = -MDP(32) * pkin(5) - MDP(29);
t202 = t128 * t173 - t176 * t133;
t201 = t222 * MDP(30);
t199 = t179 * MDP(20) + MDP(17);
t198 = qJ(6) * t236 + MDP(26);
t196 = -pkin(5) * t173 + qJ(6) * t176;
t121 = t124 + t231;
t122 = t202 - t234;
t195 = t121 * t176 + t122 * t173;
t194 = t138 * t176 + t141 * t173;
t193 = -t151 * MDP(10) + (t172 * t235 - t160) * MDP(9);
t191 = t176 * MDP(24) - t173 * MDP(25);
t190 = t173 * MDP(24) + t176 * MDP(25);
t188 = t173 * MDP(29) - t176 * MDP(31);
t187 = -MDP(16) + t191;
t184 = -MDP(27) * t202 - t124 * MDP(28) + t214 + t215 - t217;
t183 = (MDP(32) * qJ(6) + t206) * t176 + (-MDP(27) + t203) * t173;
t182 = (t227 - t228) * MDP(30) + t195 * MDP(32);
t181 = t186 * pkin(10) - MDP(18) + t190;
t170 = t177 ^ 2;
t167 = t174 ^ 2;
t147 = t205 * t172 + t160;
t146 = (-pkin(2) * t178 + t204) * t171;
t144 = (-t179 - t196) * t177;
t1 = [t172 ^ 2 * MDP(8) + t149 ^ 2 * MDP(15) + (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) * MDP(14) + (t121 ^ 2 + t122 ^ 2 + t125 ^ 2) * MDP(32) + MDP(1) + (-0.2e1 * t136 * MDP(23) + t216) * t137 + (-0.2e1 * MDP(18) * t226 + t214 + 0.2e1 * t215 - 0.2e1 * t217 - 0.2e1 * t220) * t148 + 0.2e1 * (-t121 * t136 + t122 * t137) * MDP(30) + 0.2e1 * (t127 * t136 - t148 * t202) * MDP(27) + (t121 * t148 - t125 * t137) * t236 + (-t122 * t148 + t125 * t136) * t237 + 0.2e1 * (-t124 * t148 + t127 * t137) * MDP(28) + (t147 * MDP(12) - t145 * MDP(13) + t193) * t242 + 0.2e1 * (t148 * MDP(20) + t212) * t139 + ((t175 * MDP(6) + t178 * MDP(7)) * t242 + 0.2e1 * (-t145 * MDP(11) + t146 * MDP(12)) * t178 + ((MDP(19) + MDP(4)) * t175 ^ 2 + 0.2e1 * (-MDP(10) * t175 + MDP(9) * t178) * pkin(1)) * t171 + 0.2e1 * (t147 * MDP(11) - t146 * MDP(13) + MDP(5) * t225 + t243) * t175) * t171; (-t136 * t138 + t137 * t141) * MDP(30) + t160 * MDP(12) + qJ(3) * t212 + (0.2e1 * t164 + t151) * MDP(13) + (t121 * t138 + t122 * t141) * MDP(32) + (-pkin(2) * t147 - qJ(3) * t145) * MDP(14) + ((-0.2e1 * pkin(2) - t235) * MDP(12) + MDP(8)) * t172 + (t136 * MDP(29) - t137 * MDP(31) + t218) * t144 + t241 * t148 + (MDP(20) * t139 - MDP(29) * t122 + MDP(31) * t121 + t184 - t220) * t174 + ((qJ(3) * MDP(11) + MDP(7)) * t178 + (-pkin(2) * MDP(11) + MDP(6) + (-MDP(18) - t208) * t174) * t175) * t171 + ((-t227 - t228) * MDP(23) + (-t121 * t173 + t122 * t176) * MDP(30) + (-t137 * t179 + t229) * MDP(28) + (-t136 * t179 + t230) * MDP(27) + t139 * MDP(21) + t149 * MDP(15) + t137 * t210 + t199 * t226 + t188 * t125 + t187 * t148) * t177 + t193; MDP(8) + t167 * MDP(26) + (t138 ^ 2 + t141 ^ 2 + t144 ^ 2) * MDP(32) + (-0.2e1 * MDP(12) + t232) * pkin(2) + (MDP(14) * qJ(3) + MDP(21) * t238 + 0.2e1 * MDP(13)) * qJ(3) + ((-t138 * t173 + t141 * t176) * MDP(30) + t188 * t144) * t238 + 0.2e1 * (t187 * t177 + t241) * t174 + (MDP(22) * t169 - 0.2e1 * t173 * t209 + MDP(15) + 0.2e1 * (-t173 * MDP(27) - t176 * MDP(28)) * t179) * t170; MDP(11) * t226 + t172 * MDP(12) + t147 * MDP(14) + (MDP(20) * t226 - t207 * t136 + t206 * t137 - t218) * t177 + (-MDP(21) * t226 + t186 * t148 + t182) * t174; MDP(12) - t232 + (-t144 * t177 + t194 * t174) * MDP(32) + t186 * (t167 + t170); MDP(14) + (t222 * t167 + t170) * MDP(32); MDP(19) * t226 + t173 * t216 + (-t136 * t173 + t137 * t176) * MDP(23) + (-pkin(4) * t136 - t229) * MDP(27) + (-pkin(4) * t137 + t230) * MDP(28) + (-t125 * t176 + t136 * t157) * MDP(29) + t195 * MDP(30) + (-t125 * t173 - t137 * t157) * MDP(31) + t125 * t211 + t182 * pkin(10) + t181 * t148 + t243; (-t176 * MDP(29) - t173 * MDP(31) + t211) * t144 + (t181 - t208) * t174 + (t173 * t210 + (-t166 + t169) * MDP(23) + (-pkin(4) * t173 + t223) * MDP(27) + (-pkin(4) * t176 - t224) * MDP(28) + t188 * t157 + t199) * t177 + t240 * t194; (t222 * t233 - MDP(21) + t201) * t174 + (t206 * t173 + t207 * t176 + MDP(20) - t211) * t177; MDP(19) + t166 * MDP(22) + (t222 * pkin(10) ^ 2 + t157 ^ 2) * MDP(32) + 0.2e1 * pkin(10) * t201 + 0.2e1 * (MDP(27) * pkin(4) - MDP(29) * t157) * t176 + 0.2e1 * (-MDP(28) * pkin(4) - MDP(31) * t157 + t209) * t173; (-t202 + 0.2e1 * t234) * MDP(29) + (-pkin(5) * t137 - qJ(6) * t136) * MDP(30) + (t124 + 0.2e1 * t231) * MDP(31) + (-pkin(5) * t122 + qJ(6) * t121) * MDP(32) + t184; t153 * MDP(29) + t143 * MDP(31) + (-pkin(5) * t141 + qJ(6) * t138) * MDP(32) + ((0.2e1 * pkin(5) - t224) * MDP(29) + t198) * t174 + (t197 * MDP(30) + t191) * t177 + t189; t183 * t174; t196 * MDP(30) + t183 * pkin(10) + t190; pkin(5) * t237 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(32) + t198; -t148 * MDP(29) + t137 * MDP(30) + t122 * MDP(32); MDP(30) * t176 * t177 - MDP(29) * t174 + MDP(32) * t141; t173 * t174 * MDP(32); t240 * t173; t203; MDP(32);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
