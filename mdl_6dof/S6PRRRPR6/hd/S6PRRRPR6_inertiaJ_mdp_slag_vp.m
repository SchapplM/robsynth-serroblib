% Calculate joint inertia matrix for
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_inertiaJ_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRPR6_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:36:53
% EndTime: 2019-03-08 23:36:55
% DurationCPUTime: 0.80s
% Computational Cost: add. (539->174), mult. (1113->252), div. (0->0), fcn. (1113->10), ass. (0->80)
t137 = sin(qJ(6));
t139 = sin(qJ(3));
t138 = sin(qJ(4));
t141 = cos(qJ(6));
t183 = t141 * t138;
t142 = cos(qJ(4));
t184 = t139 * t142;
t110 = t137 * t184 - t139 * t183;
t117 = t137 * t138 + t141 * t142;
t111 = t117 * t139;
t143 = cos(qJ(3));
t131 = t143 * pkin(4);
t123 = -t143 * pkin(3) - t139 * pkin(9) - pkin(2);
t185 = t138 * t143;
t180 = pkin(8) * t185 - t142 * t123;
t105 = t131 + t180;
t96 = t143 * pkin(5) - pkin(10) * t184 + t105;
t193 = pkin(8) * t142;
t108 = t138 * t123 + t143 * t193;
t104 = -t143 * qJ(5) + t108;
t186 = t138 * t139;
t97 = pkin(10) * t186 + t104;
t203 = t111 * MDP(25) - t110 * MDP(26) - (t137 * t97 - t141 * t96) * MDP(28) - (t137 * t96 + t141 * t97) * MDP(29);
t202 = t180 * MDP(17) + t108 * MDP(18) + t203;
t199 = MDP(22) * pkin(9) + MDP(20);
t198 = MDP(17) + MDP(19);
t187 = t138 * qJ(5);
t195 = pkin(4) + pkin(5);
t114 = t195 * t142 + pkin(3) + t187;
t197 = 0.2e1 * t114;
t196 = -2 * MDP(24);
t194 = pkin(9) - pkin(10);
t191 = pkin(4) * MDP(22);
t190 = cos(pkin(6));
t136 = sin(pkin(6));
t140 = sin(qJ(2));
t189 = t136 * t140;
t144 = cos(qJ(2));
t188 = t136 * t144;
t182 = t142 * qJ(5);
t132 = t138 ^ 2;
t134 = t142 ^ 2;
t179 = t132 + t134;
t178 = MDP(11) * t139;
t177 = MDP(20) * t139;
t176 = qJ(5) * MDP(21);
t174 = t111 * MDP(23);
t172 = t117 * MDP(28);
t171 = (t137 * qJ(5) + t141 * t195) * MDP(28);
t170 = (t141 * qJ(5) - t137 * t195) * MDP(29);
t169 = t142 * MDP(13);
t168 = MDP(16) + MDP(27);
t166 = MDP(18) - MDP(21);
t165 = t194 * t138;
t163 = -t198 - t191;
t162 = -t142 * pkin(4) - t187;
t161 = -pkin(4) * t138 + t182;
t158 = MDP(22) * qJ(5) - t166;
t113 = t190 * t139 + t143 * t189;
t100 = t113 * t142 - t138 * t188;
t99 = t113 * t138 + t142 * t188;
t93 = t100 * t137 - t99 * t141;
t94 = t100 * t141 + t99 * t137;
t156 = -t93 * MDP(28) - t94 * MDP(29);
t155 = t104 * t142 + t105 * t138;
t154 = t142 * MDP(14) - t138 * MDP(15);
t152 = t138 * MDP(19) - t142 * MDP(21);
t150 = t141 * MDP(28) - t137 * MDP(29);
t149 = MDP(19) + t150;
t148 = -MDP(27) - t170 - t171;
t118 = -t137 * t142 + t183;
t124 = t194 * t142;
t147 = t118 * MDP(25) - t117 * MDP(26) - (t137 * t124 - t141 * t165) * MDP(28) - (t141 * t124 + t137 * t165) * MDP(29);
t146 = -t138 * MDP(14) + t147;
t127 = pkin(9) * t185;
t122 = -pkin(3) + t162;
t112 = t139 * t189 - t190 * t143;
t109 = (pkin(8) - t161) * t139;
t101 = (-t195 * t138 - pkin(8) + t182) * t139;
t1 = [MDP(1) + (t100 ^ 2 + t112 ^ 2 + t99 ^ 2) * MDP(22); (t100 * t104 + t99 * t105 + t112 * t109) * MDP(22) + (-t112 * t110 - t93 * t143) * MDP(28) + (-t112 * t111 - t94 * t143) * MDP(29) + t166 * (t100 * t143 + t112 * t184) + (-t100 * t138 + t142 * t99) * t177 + (-t140 * MDP(4) + (MDP(10) * t143 + MDP(3) - t178) * t144) * t136 + t198 * (t112 * t186 + t99 * t143); MDP(2) - 0.2e1 * pkin(2) * t178 + (t104 ^ 2 + t105 ^ 2 + t109 ^ 2) * MDP(22) + t168 * t143 ^ 2 + (t110 * t196 + t174) * t111 + 0.2e1 * (t110 * MDP(28) + t111 * MDP(29)) * t101 + 0.2e1 * ((-t104 * t138 + t105 * t142) * MDP(20) + t152 * t109) * t139 + (t134 * MDP(12) - 0.2e1 * t138 * t169 + MDP(5) + 0.2e1 * (t138 * MDP(17) + t142 * MDP(18)) * pkin(8)) * t139 ^ 2 + 0.2e1 * (pkin(2) * MDP(10) + (MDP(6) - t154) * t139 + t105 * MDP(19) - t104 * MDP(21) + t202) * t143; -t113 * MDP(11) + (t122 * MDP(22) - t118 * MDP(29) + t166 * t138 - t142 * t198 - MDP(10) - t172) * t112 + t199 * (t100 * t142 + t99 * t138); t127 * MDP(17) + (-t109 * t142 + t127) * MDP(19) + t155 * MDP(20) - t109 * t138 * MDP(21) + (t155 * pkin(9) + t109 * t122) * MDP(22) + t118 * t174 + (-t118 * t110 - t111 * t117) * MDP(24) + (t101 * t117 + t114 * t110) * MDP(28) + (t101 * t118 + t114 * t111) * MDP(29) + (-pkin(8) * MDP(11) + MDP(8) + (t166 * pkin(9) - MDP(15)) * t142 + t146) * t143 + (MDP(7) - pkin(8) * MDP(10) + t142 * t138 * MDP(12) + (-t132 + t134) * MDP(13) + (-pkin(3) * t138 - t193) * MDP(17) + (-pkin(3) * t142 + pkin(8) * t138) * MDP(18) + t152 * t122) * t139; MDP(9) + t132 * MDP(12) + (t179 * pkin(9) ^ 2 + t122 ^ 2) * MDP(22) + t172 * t197 + 0.2e1 * t179 * MDP(20) * pkin(9) + (MDP(23) * t118 + MDP(29) * t197 + t117 * t196) * t118 + 0.2e1 * (pkin(3) * MDP(17) - t122 * MDP(19)) * t142 + 0.2e1 * (-pkin(3) * MDP(18) - t122 * MDP(21) + t169) * t138; t158 * t100 + t163 * t99 - t156; (-0.2e1 * t131 - t180) * MDP(19) + t108 * MDP(21) + (-t105 * pkin(4) + t104 * qJ(5)) * MDP(22) + (-MDP(16) + t148 - 0.2e1 * t176) * t143 + (t162 * MDP(20) + t154) * t139 - t202; t142 * MDP(15) + t161 * MDP(20) + (t163 * t138 + t158 * t142) * pkin(9) - t146; 0.2e1 * pkin(4) * MDP(19) + 0.2e1 * t176 + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(22) + 0.2e1 * t171 + 0.2e1 * t170 + t168; t99 * MDP(22); t105 * MDP(22) + t142 * t177 + t149 * t143; t199 * t138; -t149 - t191; MDP(22); t156; t143 * MDP(27) + t203; t147; t148; t150; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
