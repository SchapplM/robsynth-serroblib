% Calculate joint inertia matrix for
% S6RRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR6_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR6_inertiaJ_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:53:27
% EndTime: 2019-03-09 13:53:30
% DurationCPUTime: 0.69s
% Computational Cost: add. (799->151), mult. (1377->196), div. (0->0), fcn. (1496->8), ass. (0->79)
t151 = sin(qJ(6));
t155 = cos(qJ(6));
t180 = t155 * MDP(34);
t166 = -t151 * MDP(35) + t180;
t153 = sin(qJ(4));
t157 = cos(qJ(4));
t152 = sin(qJ(5));
t156 = cos(qJ(5));
t124 = t152 * t153 - t156 * t157;
t127 = t152 * t157 + t156 * t153;
t163 = -t127 * MDP(28) + (-MDP(27) - t166) * t124;
t206 = t157 * MDP(20) - t153 * MDP(21) + t163;
t191 = t151 * MDP(31) + t155 * MDP(32);
t154 = sin(qJ(2));
t158 = cos(qJ(2));
t125 = t154 * t153 + t158 * t157;
t126 = -t158 * t153 + t154 * t157;
t109 = t156 * t125 + t152 * t126;
t110 = -t152 * t125 + t156 * t126;
t147 = t151 ^ 2;
t149 = t155 ^ 2;
t192 = t151 * t155;
t196 = pkin(7) - pkin(8);
t133 = t196 * t154;
t134 = t196 * t158;
t116 = -t157 * t133 + t153 * t134;
t101 = -t126 * pkin(9) - t116;
t117 = t153 * t133 + t157 * t134;
t102 = -t125 * pkin(9) + t117;
t96 = -t156 * t101 + t152 * t102;
t97 = t152 * t101 + t156 * t102;
t205 = -t96 * MDP(27) - t97 * MDP(28) + (MDP(29) * t192 + MDP(24) - (t147 - t149) * MDP(30)) * t110 + (t191 - MDP(25)) * t109;
t201 = MDP(34) * t151 + MDP(35) * t155;
t199 = t126 * MDP(17) - t125 * MDP(18) - t116 * MDP(20) - t117 * MDP(21) + t205;
t132 = -t158 * pkin(2) - t154 * qJ(3) - pkin(1);
t123 = t158 * pkin(3) - t132;
t115 = t125 * pkin(4) + t123;
t198 = 0.2e1 * t115;
t197 = 0.2e1 * t123;
t195 = t156 * pkin(4);
t93 = t96 * t151;
t193 = t96 * t155;
t148 = t154 ^ 2;
t190 = t158 ^ 2 + t148;
t187 = t109 * MDP(33);
t159 = -pkin(2) - pkin(3);
t130 = t153 * qJ(3) - t157 * t159;
t129 = -pkin(4) - t130;
t131 = t157 * qJ(3) + t153 * t159;
t113 = -t156 * t129 + t152 * t131;
t186 = t113 * MDP(27);
t114 = t152 * t129 + t156 * t131;
t185 = t114 * MDP(28);
t184 = t130 * MDP(20);
t183 = t131 * MDP(21);
t178 = t147 * MDP(29) + MDP(26);
t176 = MDP(30) * t192;
t135 = 0.2e1 * t176;
t175 = t135 + t178;
t174 = MDP(19) + t178;
t173 = -pkin(2) * MDP(14) - MDP(11);
t136 = -0.2e1 * t176;
t172 = t135 + t174;
t171 = -pkin(5) * t110 - pkin(10) * t109;
t111 = pkin(5) + t113;
t112 = -pkin(10) + t114;
t170 = -t109 * t112 + t110 * t111;
t144 = t152 * pkin(4);
t137 = t144 + pkin(10);
t138 = -pkin(5) - t195;
t168 = -t109 * t137 + t110 * t138;
t167 = MDP(31) * t155 - MDP(32) * t151;
t164 = (t156 * MDP(27) - t152 * MDP(28)) * pkin(4);
t162 = 0.2e1 * t166;
t107 = t111 * t151;
t92 = t109 * pkin(5) - t110 * pkin(10) + t115;
t91 = t151 * t92 + t155 * t97;
t90 = -t151 * t97 + t155 * t92;
t1 = [t148 * MDP(4) + (pkin(7) ^ 2 * t190 + t132 ^ 2) * MDP(14) + t125 * MDP(20) * t197 + t110 * MDP(28) * t198 + MDP(1) + (MDP(15) * t126 - 0.2e1 * t125 * MDP(16) + MDP(21) * t197) * t126 + (t149 * MDP(29) + MDP(22) + t136) * t110 ^ 2 + (MDP(27) * t198 + t187) * t109 + 0.2e1 * (t90 * t109 + t110 * t93) * MDP(34) + 0.2e1 * (-t91 * t109 + t110 * t193) * MDP(35) + 0.2e1 * t190 * MDP(12) * pkin(7) + 0.2e1 * (-t132 * MDP(11) + pkin(1) * MDP(9)) * t158 + 0.2e1 * (-pkin(1) * MDP(10) - t132 * MDP(13) + t158 * MDP(5)) * t154 + 0.2e1 * (-MDP(23) + t167) * t110 * t109; ((MDP(14) * qJ(3) - MDP(10) + MDP(13)) * t158 + (-MDP(9) + t173) * t154) * pkin(7) + (t155 * t170 - t93) * MDP(35) + (t151 * t170 + t193) * MDP(34) + t154 * MDP(6) + t158 * MDP(7) + (-t154 * pkin(2) + t158 * qJ(3)) * MDP(12) - t199; MDP(8) + 0.2e1 * pkin(2) * MDP(11) + 0.2e1 * qJ(3) * MDP(13) + (pkin(2) ^ 2 + qJ(3) ^ 2) * MDP(14) + t111 * t162 + 0.2e1 * t184 + 0.2e1 * t183 + 0.2e1 * t186 + 0.2e1 * t185 + t172; (pkin(7) * MDP(14) + MDP(12)) * t154 + t201 * (-t109 * t127 + t110 * t124); t173 - t206; MDP(14); (t151 * t168 - t193) * MDP(34) + (t155 * t168 + t93) * MDP(35) + t199; -t184 - t183 + (-t113 - t195) * MDP(27) + (-t114 + t144) * MDP(28) + t136 + (-t138 * t151 + t107) * MDP(35) + (-t111 + t138) * t180 - t174; t206; -0.2e1 * t138 * t166 + 0.2e1 * t164 + t172; (t151 * t171 - t193) * MDP(34) + (t155 * t171 + t93) * MDP(35) + t205; -t186 - t185 + t136 + (pkin(5) * t151 + t107) * MDP(35) + (-pkin(5) - t111) * t180 - t178; t163; t164 + t175 + t166 * (pkin(5) - t138); pkin(5) * t162 + t175; t90 * MDP(34) - t91 * MDP(35) + t110 * t167 + t187; -t112 * t201 - t191; -t201 * t127; -t137 * t201 + t191; -pkin(10) * t201 + t191; MDP(33);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
