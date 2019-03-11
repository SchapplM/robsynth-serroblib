% Calculate joint inertia matrix for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPPR7_inertiaJ_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:01:08
% EndTime: 2019-03-09 16:01:12
% DurationCPUTime: 1.29s
% Computational Cost: add. (1082->221), mult. (1975->310), div. (0->0), fcn. (1964->8), ass. (0->88)
t180 = sin(pkin(10));
t181 = cos(pkin(10));
t183 = sin(qJ(3));
t186 = cos(qJ(3));
t149 = t180 * t183 + t181 * t186;
t184 = sin(qJ(2));
t145 = t149 * t184;
t187 = cos(qJ(2));
t175 = t187 * pkin(3);
t162 = -pkin(2) * t187 - pkin(8) * t184 - pkin(1);
t226 = t183 * t187;
t224 = pkin(7) * t226 - t186 * t162;
t139 = t175 + t224;
t225 = t184 * t186;
t126 = pkin(4) * t187 - qJ(5) * t225 + t139;
t230 = pkin(7) * t186;
t142 = t183 * t162 + t187 * t230;
t138 = -qJ(4) * t187 + t142;
t227 = t183 * t184;
t130 = qJ(5) * t227 + t138;
t204 = -t181 * t126 + t130 * t180;
t113 = pkin(5) * t187 - pkin(9) * t145 - t204;
t118 = t180 * t126 + t181 * t130;
t144 = t180 * t225 - t181 * t227;
t114 = -pkin(9) * t144 + t118;
t182 = sin(qJ(6));
t185 = cos(qJ(6));
t121 = t185 * t144 + t145 * t182;
t122 = -t144 * t182 + t145 * t185;
t242 = -(t113 * t182 + t114 * t185) * MDP(32) + (t113 * t185 - t114 * t182) * MDP(31) + t122 * MDP(28) - t121 * MDP(29);
t241 = 2 * MDP(19);
t240 = -t224 * MDP(16) - t142 * MDP(17) + t204 * MDP(22) + t118 * MDP(23) - t242;
t239 = 0.2e1 * t187;
t237 = -t186 * pkin(3) - t183 * qJ(4);
t161 = -pkin(2) + t237;
t148 = t186 * pkin(4) - t161;
t133 = pkin(5) * t149 + t148;
t235 = 0.2e1 * t133;
t234 = 2 * MDP(24);
t233 = -2 * MDP(27);
t232 = -pkin(3) - pkin(4);
t231 = pkin(3) * t183;
t229 = pkin(8) - qJ(5);
t228 = pkin(3) * MDP(21);
t163 = t229 * t186;
t205 = t229 * t183;
t136 = t181 * t163 + t180 * t205;
t176 = t183 ^ 2;
t178 = t186 ^ 2;
t223 = t176 + t178;
t222 = MDP(25) * t148;
t221 = qJ(4) * MDP(20);
t217 = t122 * MDP(26);
t151 = -t180 * t186 + t181 * t183;
t128 = t185 * t149 + t151 * t182;
t215 = t128 * MDP(31);
t160 = t181 * qJ(4) + t180 * t232;
t158 = qJ(4) * t180 - t181 * t232;
t203 = -pkin(5) - t158;
t214 = (t160 * t182 - t185 * t203) * MDP(31);
t213 = (t185 * t160 + t182 * t203) * MDP(32);
t212 = t149 * MDP(22);
t211 = t151 * MDP(23);
t210 = t158 * MDP(22);
t209 = t160 * MDP(23);
t208 = t186 * MDP(12);
t207 = MDP(17) - MDP(20);
t206 = MDP(30) + MDP(15);
t134 = t163 * t180 - t181 * t205;
t202 = t138 * t186 + t139 * t183;
t201 = t186 * MDP(13) - t183 * MDP(14);
t199 = t183 * MDP(18) - t186 * MDP(20);
t198 = t144 * MDP(22) + t145 * MDP(23);
t195 = t121 * MDP(31) + t122 * MDP(32);
t194 = -(t180 * t182 - t185 * t181) * MDP(31) - (t180 * t185 + t181 * t182) * MDP(32);
t193 = -pkin(9) * t151 - t134;
t192 = -MDP(30) - t213 - t214;
t165 = qJ(4) * t225;
t137 = t165 + (t232 * t183 - pkin(7)) * t184;
t123 = -pkin(9) * t149 + t136;
t129 = -t149 * t182 + t151 * t185;
t191 = t129 * MDP(28) - t128 * MDP(29) - (t123 * t182 - t185 * t193) * MDP(31) - (t185 * t123 + t182 * t193) * MDP(32);
t190 = t181 * MDP(22) - t180 * MDP(23) + MDP(18) + t194;
t189 = -t183 * MDP(13) - t134 * MDP(22) - t136 * MDP(23) + t191;
t167 = pkin(8) * t226;
t143 = -t165 + (pkin(7) + t231) * t184;
t119 = pkin(5) * t144 + t137;
t1 = [(t138 ^ 2 + t139 ^ 2 + t143 ^ 2) * MDP(21) + (t118 ^ 2 + t137 ^ 2 + t204 ^ 2) * MDP(25) + MDP(1) + t206 * t187 ^ 2 + (t121 * t233 + t217) * t122 + (-t118 * t144 + t145 * t204) * t234 + 0.2e1 * t198 * t137 + 0.2e1 * t195 * t119 + (t139 * MDP(18) - t138 * MDP(20) + pkin(1) * MDP(9) - t240) * t239 + (-0.2e1 * pkin(1) * MDP(10) + (MDP(5) - t201) * t239 + (-t138 * t183 + t139 * t186) * t241 + 0.2e1 * t199 * t143 + (t178 * MDP(11) - 0.2e1 * t183 * t208 + MDP(4) + 0.2e1 * (t183 * MDP(16) + t186 * MDP(17)) * pkin(7)) * t184) * t184; t167 * MDP(16) + (-t143 * t186 + t167) * MDP(18) + t202 * MDP(19) - t143 * t183 * MDP(20) + (t202 * pkin(8) + t143 * t161) * MDP(21) + (t137 * t149 + t144 * t148) * MDP(22) + (t137 * t151 + t145 * t148) * MDP(23) + (-t118 * t149 + t134 * t145 - t136 * t144 + t151 * t204) * MDP(24) + (t118 * t136 + t134 * t204 + t137 * t148) * MDP(25) + (-t121 * t129 - t122 * t128) * MDP(27) + (t119 * t128 + t121 * t133) * MDP(31) + (t119 * t129 + t122 * t133) * MDP(32) + t129 * t217 + (-pkin(7) * MDP(10) + MDP(7) + (t207 * pkin(8) - MDP(14)) * t186 + t189) * t187 + (MDP(6) + (-t176 + t178) * MDP(12) + (-pkin(2) * t183 - t230) * MDP(16) + (-pkin(2) * t186 + pkin(7) * t183) * MDP(17) + t186 * t183 * MDP(11) - pkin(7) * MDP(9) + t199 * t161) * t184; MDP(8) + t176 * MDP(11) + (t223 * pkin(8) ^ 2 + t161 ^ 2) * MDP(21) + (t134 ^ 2 + t136 ^ 2) * MDP(25) + t215 * t235 + (0.2e1 * t211 + 0.2e1 * t212 + t222) * t148 + (t134 * t151 - t136 * t149) * t234 + t223 * pkin(8) * t241 + (MDP(26) * t129 + MDP(32) * t235 + t128 * t233) * t129 + 0.2e1 * (pkin(2) * MDP(16) - t161 * MDP(18)) * t186 + 0.2e1 * (-pkin(2) * MDP(17) - t161 * MDP(20) + t208) * t183; (-0.2e1 * t175 - t224) * MDP(18) + t142 * MDP(20) + (-pkin(3) * t139 + qJ(4) * t138) * MDP(21) + (-t144 * t160 + t145 * t158) * MDP(24) + (t118 * t160 + t158 * t204) * MDP(25) + (-MDP(15) + t192 - t209 - t210 - 0.2e1 * t221) * t187 + (t237 * MDP(19) + t201) * t184 + t240; t186 * MDP(14) + (qJ(4) * t186 - t231) * MDP(19) + (-t149 * t160 + t151 * t158) * MDP(24) + (t134 * t158 + t136 * t160) * MDP(25) + ((MDP(21) * qJ(4) - t207) * t186 + (-MDP(16) - MDP(18) - t228) * t183) * pkin(8) - t189; 0.2e1 * pkin(3) * MDP(18) + 0.2e1 * t221 + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(21) + (t158 ^ 2 + t160 ^ 2) * MDP(25) + 0.2e1 * t210 + 0.2e1 * t209 + 0.2e1 * t214 + 0.2e1 * t213 + t206; MDP(19) * t225 + t139 * MDP(21) + (-t144 * t180 - t145 * t181) * MDP(24) + (t118 * t180 - t181 * t204) * MDP(25) + t190 * t187; (-t149 * t180 - t151 * t181) * MDP(24) + (-t134 * t181 + t136 * t180) * MDP(25) + (MDP(21) * pkin(8) + MDP(19)) * t183; -t228 + (-t158 * t181 + t160 * t180) * MDP(25) - t190; MDP(21) + (t180 ^ 2 + t181 ^ 2) * MDP(25); t137 * MDP(25) + t195 + t198; t129 * MDP(32) + t211 + t212 + t215 + t222; 0; 0; MDP(25); t187 * MDP(30) + t242; t191; t192; t194; 0; MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
