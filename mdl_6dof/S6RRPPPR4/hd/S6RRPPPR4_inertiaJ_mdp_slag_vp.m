% Calculate joint inertia matrix for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_inertiaJ_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RRPPPR4_inertiaJ_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:19:46
% EndTime: 2019-03-09 08:19:48
% DurationCPUTime: 0.56s
% Computational Cost: add. (511->150), mult. (855->198), div. (0->0), fcn. (743->6), ass. (0->65)
t176 = pkin(3) + pkin(7);
t131 = sin(pkin(9));
t132 = cos(pkin(9));
t117 = t131 ^ 2 + t132 ^ 2;
t169 = sin(qJ(6));
t170 = cos(qJ(6));
t134 = sin(qJ(2));
t135 = cos(qJ(2));
t164 = t131 * t135;
t171 = -pkin(4) - pkin(5);
t133 = -pkin(2) - qJ(4);
t150 = -qJ(3) * t134 - pkin(1);
t103 = t133 * t135 + t150;
t113 = t176 * t134;
t91 = -t131 * t103 + t113 * t132;
t87 = pkin(8) * t164 + t171 * t134 - t91;
t163 = t132 * t135;
t92 = t132 * t103 + t131 * t113;
t89 = t134 * qJ(5) + t92;
t88 = pkin(8) * t163 + t89;
t141 = t169 * t131 + t170 * t132;
t97 = t141 * t135;
t106 = -t170 * t131 + t169 * t132;
t98 = t106 * t135;
t175 = (-t169 * t88 + t170 * t87) * MDP(28) - (t169 * t87 + t170 * t88) * MDP(29) + t98 * MDP(25) + t97 * MDP(26);
t90 = -pkin(4) * t134 - t91;
t85 = t131 * t89 - t132 * t90;
t86 = t131 * t92 + t132 * t91;
t172 = MDP(18) * t86 + MDP(22) * t85;
t168 = pkin(8) + t133;
t165 = qJ(5) * t131;
t162 = t117 * t133 ^ 2;
t114 = t176 * t135;
t161 = MDP(16) * t132;
t148 = qJ(5) * t132 - qJ(3);
t111 = pkin(4) * t131 - t148;
t160 = MDP(22) * t111;
t159 = MDP(23) * t141;
t158 = MDP(28) * t106;
t157 = t114 * MDP(18);
t156 = MDP(15) + MDP(19);
t155 = MDP(16) - MDP(21);
t154 = -MDP(17) - MDP(20);
t153 = MDP(18) + MDP(22);
t149 = -pkin(2) * MDP(14) + MDP(12);
t147 = t155 * t131;
t144 = -t97 * MDP(28) + t98 * MDP(29);
t143 = t132 * MDP(15) - t131 * MDP(16);
t142 = -MDP(28) * t141 + MDP(29) * t106;
t140 = t170 * MDP(28) - t169 * MDP(29);
t139 = t156 * t132 - t147;
t109 = t168 * t131;
t110 = t168 * t132;
t138 = -t141 * MDP(25) + t106 * MDP(26) - (-t169 * t109 - t170 * t110) * MDP(28) + (t170 * t109 - t169 * t110) * MDP(29);
t137 = pkin(7) ^ 2;
t136 = qJ(3) ^ 2;
t130 = t135 ^ 2;
t129 = t134 ^ 2;
t115 = t132 * t133 * t134;
t112 = -pkin(2) * t135 + t150;
t104 = t117 * t133;
t102 = t171 * t131 + t148;
t96 = (pkin(4) * t132 + t165) * t135 + t114;
t93 = (t171 * t132 - t165) * t135 - t114;
t1 = [(t112 ^ 2 + t130 * t137) * MDP(14) + (t114 ^ 2 + t91 ^ 2 + t92 ^ 2) * MDP(18) + (t89 ^ 2 + t90 ^ 2 + t96 ^ 2) * MDP(22) + MDP(1) - (-MDP(23) * t98 - 0.2e1 * t97 * MDP(24)) * t98 + (t137 * MDP(14) + MDP(27) + MDP(4)) * t129 + 0.2e1 * t144 * t93 + 0.2e1 * (t129 + t130) * MDP(11) * pkin(7) + 0.2e1 * (-pkin(1) * MDP(10) - t112 * MDP(13) + t91 * MDP(15) - t92 * MDP(16) - t90 * MDP(19) + t89 * MDP(21) + t135 * MDP(5) - t175) * t134 + 0.2e1 * ((t131 * t91 - t132 * t92) * MDP(17) + (-t131 * t90 - t132 * t89) * MDP(20) + (t132 * MDP(19) + t131 * MDP(21)) * t96 + t143 * t114 + t112 * MDP(12) + pkin(1) * MDP(9)) * t135; t135 * MDP(7) + (t114 * t131 + t115) * MDP(15) + t114 * t161 - t86 * MDP(17) + (t111 * t163 + t131 * t96 + t115) * MDP(19) - t85 * MDP(20) + (t111 * t164 - t132 * t96) * MDP(21) + t96 * t160 + t98 * t159 + (-t106 * t98 + t141 * t97) * MDP(24) + (-t102 * t97 + t106 * t93) * MDP(28) + (t102 * t98 + t141 * t93) * MDP(29) + (t157 + (MDP(11) + t143) * t135) * qJ(3) + t172 * t133 + (-pkin(2) * MDP(11) - t133 * t147 + MDP(6) + t138) * t134 + ((MDP(14) * qJ(3) - MDP(10) + MDP(13)) * t135 + (-MDP(9) + t149) * t134) * pkin(7); MDP(8) - 0.2e1 * pkin(2) * MDP(12) + (pkin(2) ^ 2 + t136) * MDP(14) + (t136 + t162) * MDP(18) + t162 * MDP(22) + 0.2e1 * t102 * t158 + (0.2e1 * MDP(19) * t131 - 0.2e1 * MDP(21) * t132 + t160) * t111 + 0.2e1 * t154 * t104 - (0.2e1 * MDP(24) * t106 - 0.2e1 * MDP(29) * t102 - t159) * t141 + 0.2e1 * (MDP(15) * t131 + MDP(13) + t161) * qJ(3); (MDP(14) * pkin(7) + MDP(11) + t139 - t142) * t134 + t172; t153 * t104 + t154 * t117 + t149; t153 * t117 + MDP(14); MDP(22) * t96 + t139 * t135 - t144 + t157; qJ(3) * MDP(18) - MDP(29) * t141 + t156 * t131 + t155 * t132 - t158 + t160; 0; t153; -MDP(20) * t164 + t90 * MDP(22) + (-MDP(19) - t140) * t134; (-MDP(22) * t133 + MDP(20)) * t132; -t132 * MDP(22); 0; MDP(22); -MDP(27) * t134 + t175; -t138; t142; 0; t140; MDP(27);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
