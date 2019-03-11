% Calculate joint inertia matrix for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPRP5_inertiaJ_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:51:36
% EndTime: 2019-03-09 16:51:39
% DurationCPUTime: 1.01s
% Computational Cost: add. (1744->233), mult. (3368->331), div. (0->0), fcn. (3592->8), ass. (0->89)
t179 = MDP(25) + MDP(27);
t178 = MDP(26) - MDP(29);
t163 = sin(qJ(3));
t165 = cos(qJ(3));
t210 = t163 * MDP(13) + t165 * MDP(14);
t209 = t163 * MDP(16);
t208 = MDP(15) + MDP(24);
t206 = 2 * MDP(18);
t205 = -2 * MDP(21);
t204 = 2 * MDP(26);
t203 = 2 * MDP(27);
t202 = 2 * MDP(28);
t201 = 2 * MDP(29);
t200 = cos(qJ(5));
t160 = sin(pkin(10));
t199 = pkin(3) * t160;
t166 = cos(qJ(2));
t198 = pkin(5) * t166;
t197 = pkin(7) * t163;
t196 = pkin(7) * t166;
t164 = sin(qJ(2));
t155 = t164 * pkin(7);
t195 = -qJ(4) - pkin(8);
t194 = MDP(18) * pkin(3);
t193 = MDP(19) * pkin(3);
t192 = qJ(6) * t166;
t148 = t195 * t163;
t149 = t195 * t165;
t161 = cos(pkin(10));
t131 = t160 * t148 - t161 * t149;
t171 = t160 * t163 - t161 * t165;
t119 = -t171 * pkin(9) + t131;
t162 = sin(qJ(5));
t130 = t148 * t161 + t149 * t160;
t142 = t160 * t165 + t161 * t163;
t168 = -pkin(9) * t142 + t130;
t109 = t119 * t162 - t200 * t168;
t191 = t109 * t166;
t110 = t200 * t119 + t162 * t168;
t190 = t110 * t166;
t189 = t163 * t164;
t188 = t164 * t165;
t187 = t165 * t166;
t147 = -pkin(2) * t166 - pkin(8) * t164 - pkin(1);
t143 = t165 * t147;
t126 = -qJ(4) * t188 + t143 + (-pkin(3) - t197) * t166;
t177 = pkin(7) * t187;
t129 = t177 + (-qJ(4) * t164 + t147) * t163;
t112 = t161 * t126 - t129 * t160;
t137 = -t160 * t189 + t161 * t188;
t106 = -pkin(4) * t166 - pkin(9) * t137 + t112;
t113 = t160 * t126 + t161 * t129;
t169 = t142 * t164;
t111 = -pkin(9) * t169 + t113;
t98 = t162 * t106 + t200 * t111;
t146 = pkin(3) * t189 + t155;
t125 = t200 * t142 - t162 * t171;
t186 = t125 * MDP(20);
t153 = pkin(3) * t161 + pkin(4);
t139 = t200 * t153 - t162 * t199;
t185 = t139 * MDP(25);
t140 = t162 * t153 + t200 * t199;
t184 = t140 * MDP(26);
t181 = t166 * MDP(22);
t180 = t166 * MDP(23);
t154 = -t165 * pkin(3) - pkin(2);
t176 = MDP(14) * t189;
t175 = t163 * t165 * MDP(12);
t174 = MDP(13) * t188;
t173 = -t200 * t106 + t111 * t162;
t117 = t137 * t162 + t200 * t169;
t118 = t200 * t137 - t162 * t169;
t172 = t118 * MDP(22) - t117 * MDP(23) - t166 * MDP(24);
t124 = t142 * t162 + t200 * t171;
t170 = t125 * MDP(22) - t124 * MDP(23) - t179 * t109 - t178 * t110;
t135 = t171 * pkin(4) + t154;
t127 = pkin(4) * t169 + t146;
t158 = t165 ^ 2;
t157 = t164 ^ 2;
t156 = t163 ^ 2;
t138 = -pkin(5) - t139;
t136 = qJ(6) + t140;
t134 = t147 * t163 + t177;
t133 = -t163 * t196 + t143;
t107 = t124 * pkin(5) - t125 * qJ(6) + t135;
t99 = t117 * pkin(5) - t118 * qJ(6) + t127;
t96 = t173 + t198;
t95 = t98 - t192;
t1 = [-0.2e1 * pkin(1) * t164 * MDP(10) + 0.2e1 * (pkin(7) * t157 * t165 + t134 * t166) * MDP(17) + 0.2e1 * (t117 * t127 + t166 * t173) * MDP(25) + t158 * t157 * MDP(11) + (t117 * t99 + t166 * t96) * t203 + t166 * t98 * t204 + MDP(1) + (-t112 * t137 - t113 * t169) * t206 - 0.2e1 * t166 * t174 - 0.2e1 * t157 * t175 + 0.2e1 * t117 * t180 + t157 * MDP(4) + (t112 ^ 2 + t113 ^ 2 + t146 ^ 2) * MDP(19) + 0.2e1 * (-t133 * t166 + t157 * t197) * MDP(16) + (t96 ^ 2 + t99 ^ 2) * MDP(30) + (t95 * MDP(30) - t117 * t202 - t166 * t201) * t95 + t208 * t166 ^ 2 + 0.2e1 * (t164 * MDP(5) + pkin(1) * MDP(9) + t176) * t166 + (MDP(20) * t118 + t117 * t205 + t127 * t204 - t99 * t201 + t96 * t202 - 0.2e1 * t181) * t118; -MDP(9) * t155 - MDP(10) * t196 + pkin(8) * t187 * MDP(17) + (-t112 * t142 - t113 * t171 - t130 * t137 - t131 * t169) * MDP(18) + (t112 * t130 + t113 * t131 + t146 * t154) * MDP(19) + t118 * t186 + (-t117 * t125 - t118 * t124) * MDP(21) - t125 * t181 + t124 * t180 + (t117 * t135 + t124 * t127 + t191) * MDP(25) + (t118 * t135 + t125 * t127 + t190) * MDP(26) + (t107 * t117 + t124 * t99 + t191) * MDP(27) + (t109 * t118 - t110 * t117 - t124 * t95 + t125 * t96) * MDP(28) + (-t107 * t118 - t125 * t99 - t190) * MDP(29) + (t107 * t99 + t109 * t96 + t110 * t95) * MDP(30) + (t163 * MDP(11) - pkin(7) * MDP(16)) * t188 + (MDP(6) + (-t156 + t158) * MDP(12) - pkin(2) * t209 + (-pkin(2) * t165 + t197) * MDP(17)) * t164 + (pkin(8) * t209 + MDP(7) - t210) * t166; MDP(8) + t156 * MDP(11) + 0.2e1 * t175 + (-t130 * t142 - t131 * t171) * t206 + (t130 ^ 2 + t131 ^ 2 + t154 ^ 2) * MDP(19) + (t107 ^ 2 + t109 ^ 2 + t110 ^ 2) * MDP(30) + (-0.2e1 * t107 * MDP(29) + t109 * t202 + t124 * t205 + t135 * t204 + t186) * t125 + 0.2e1 * (t165 * MDP(16) - t163 * MDP(17)) * pkin(2) + 0.2e1 * (t135 * MDP(25) + t107 * MDP(27) - t110 * MDP(28)) * t124; t174 - t176 - t166 * MDP(15) + t133 * MDP(16) - t134 * MDP(17) + (-t161 * t137 - t160 * t169) * t194 + (t112 * t161 + t113 * t160) * t193 + (-t139 * t166 - t173) * MDP(25) + (t140 * t166 - t98) * MDP(26) + ((-pkin(5) + t138) * t166 - t173) * MDP(27) + (-t117 * t136 + t118 * t138) * MDP(28) + ((-qJ(6) - t136) * t166 + t98) * MDP(29) + (t136 * t95 + t138 * t96) * MDP(30) + t172; (-t161 * t142 - t160 * t171) * t194 + (t130 * t161 + t131 * t160) * t193 + (-t124 * t136 + t125 * t138) * MDP(28) + (t109 * t138 + t110 * t136) * MDP(30) + t170 + (-t165 * MDP(17) - t209) * pkin(8) + t210; (t136 ^ 2 + t138 ^ 2) * MDP(30) + (t160 ^ 2 + t161 ^ 2) * MDP(19) * pkin(3) ^ 2 + 0.2e1 * t185 - 0.2e1 * t184 - 0.2e1 * MDP(27) * t138 + t136 * t201 + t208; MDP(19) * t146 + t99 * MDP(30) + t179 * t117 + t178 * t118; MDP(19) * t154 + t107 * MDP(30) + t179 * t124 + t178 * t125; 0; MDP(19) + MDP(30); -t173 * MDP(25) - t98 * MDP(26) + (-t173 - 0.2e1 * t198) * MDP(27) + (-pkin(5) * t118 - qJ(6) * t117) * MDP(28) + (t98 - 0.2e1 * t192) * MDP(29) + (-pkin(5) * t96 + qJ(6) * t95) * MDP(30) + t172; (-pkin(5) * t125 - qJ(6) * t124) * MDP(28) + (-pkin(5) * t109 + qJ(6) * t110) * MDP(30) + t170; MDP(24) + t185 - t184 + (0.2e1 * pkin(5) + t139) * MDP(27) + (0.2e1 * qJ(6) + t140) * MDP(29) + (-pkin(5) * t138 + qJ(6) * t136) * MDP(30); 0; MDP(24) + pkin(5) * t203 + qJ(6) * t201 + (pkin(5) ^ 2 + qJ(6) ^ 2) * MDP(30); t166 * MDP(27) + t118 * MDP(28) + t96 * MDP(30); t125 * MDP(28) + t109 * MDP(30); MDP(30) * t138 - MDP(27); 0; -MDP(30) * pkin(5) - MDP(27); MDP(30);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
