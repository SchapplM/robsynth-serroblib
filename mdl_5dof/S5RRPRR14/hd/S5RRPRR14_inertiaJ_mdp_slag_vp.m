% Calculate joint inertia matrix for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR14_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR14_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR14_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:38:40
% EndTime: 2019-12-31 20:38:42
% DurationCPUTime: 0.61s
% Computational Cost: add. (921->174), mult. (2120->259), div. (0->0), fcn. (2340->10), ass. (0->88)
t131 = sin(pkin(5));
t182 = 0.2e1 * t131;
t181 = 2 * MDP(13);
t180 = 2 * MDP(20);
t179 = 2 * MDP(27);
t178 = 2 * MDP(28);
t177 = cos(qJ(4));
t136 = sin(qJ(2));
t176 = pkin(1) * t136;
t138 = cos(qJ(2));
t175 = pkin(1) * t138;
t174 = pkin(8) + qJ(3);
t173 = pkin(2) * MDP(14);
t130 = sin(pkin(10));
t132 = cos(pkin(10));
t133 = cos(pkin(5));
t169 = t131 * t136;
t112 = t130 * t169 - t132 * t133;
t113 = t130 * t133 + t132 * t169;
t135 = sin(qJ(4));
t100 = -t135 * t112 + t177 * t113;
t134 = sin(qJ(5));
t137 = cos(qJ(5));
t168 = t131 * t138;
t96 = t100 * t137 - t134 * t168;
t172 = MDP(22) * t96;
t99 = t177 * t112 + t113 * t135;
t171 = MDP(26) * t99;
t120 = t174 * t132;
t148 = t174 * t130;
t103 = t120 * t135 + t177 * t148;
t119 = t177 * t130 + t135 * t132;
t170 = t103 * t119;
t167 = t133 * MDP(8);
t166 = t136 * MDP(6);
t95 = t100 * t134 + t137 * t168;
t165 = t95 * MDP(23);
t164 = t95 * MDP(25);
t163 = t96 * MDP(24);
t150 = pkin(7) * t168;
t108 = t150 + (qJ(3) + t176) * t133;
t109 = (-pkin(2) * t138 - qJ(3) * t136 - pkin(1)) * t131;
t98 = t132 * t108 + t130 * t109;
t161 = MDP(11) * t132;
t160 = MDP(12) * t130;
t122 = pkin(7) * t169;
t111 = t122 + (-pkin(2) - t175) * t133;
t159 = MDP(14) * t111;
t158 = MDP(19) * t138;
t157 = MDP(21) * t119;
t156 = MDP(22) * t134;
t155 = t100 * MDP(16);
t154 = t100 * MDP(17);
t153 = t100 * MDP(21);
t104 = t177 * t120 - t135 * t148;
t152 = t104 * MDP(21);
t118 = t130 * t135 - t177 * t132;
t151 = t118 * MDP(26);
t124 = -pkin(3) * t132 - pkin(2);
t149 = MDP(23) * t134 * t137;
t97 = -t108 * t130 + t132 * t109;
t147 = -t130 * t97 + t132 * t98;
t146 = MDP(24) * t137 - MDP(25) * t134;
t145 = MDP(27) * t137 - MDP(28) * t134;
t144 = -MDP(27) * t134 - MDP(28) * t137;
t91 = -pkin(3) * t168 - pkin(8) * t113 + t97;
t94 = -pkin(8) * t112 + t98;
t88 = -t135 * t94 + t177 * t91;
t89 = t135 * t91 + t177 * t94;
t143 = MDP(20) + t145;
t101 = pkin(3) * t112 + t111;
t142 = t134 * MDP(24) + t137 * MDP(25) + t144 * pkin(9);
t87 = -pkin(9) * t168 + t89;
t90 = pkin(4) * t99 - pkin(9) * t100 + t101;
t84 = -t134 * t87 + t137 * t90;
t85 = t134 * t90 + t137 * t87;
t141 = t84 * MDP(27) - t85 * MDP(28) + t163 - t164 + t171;
t140 = -MDP(18) + t142;
t129 = t137 ^ 2;
t128 = t134 ^ 2;
t126 = t131 ^ 2;
t115 = t133 * t176 + t150;
t114 = t133 * t175 - t122;
t102 = pkin(4) * t118 - pkin(9) * t119 + t124;
t93 = t102 * t134 + t104 * t137;
t92 = t102 * t137 - t104 * t134;
t86 = pkin(4) * t168 - t88;
t1 = [(t111 ^ 2 + t97 ^ 2 + t98 ^ 2) * MDP(14) + t100 ^ 2 * MDP(15) + t126 * t136 ^ 2 * MDP(4) + MDP(1) + (-0.2e1 * t165 + t172) * t96 + (t166 * t182 + t167) * t133 + (-0.2e1 * t155 + 0.2e1 * t163 - 0.2e1 * t164 + t171) * t99 + ((0.2e1 * MDP(5) * t136 + t158) * t126 + (t99 * MDP(18) + MDP(7) * t133 - t154) * t182) * t138 + (t84 * t99 + t86 * t95) * t179 + (-t85 * t99 + t86 * t96) * t178 + (-t112 * t98 - t113 * t97) * t181 + 0.2e1 * (-t115 * t133 - t126 * t176) * MDP(10) + (t101 * t99 - t88 * t168) * t180 + 0.2e1 * (t100 * t101 + t89 * t168) * MDP(21) + 0.2e1 * (t111 * t112 - t97 * t168) * MDP(11) + 0.2e1 * (t111 * t113 + t98 * t168) * MDP(12) + 0.2e1 * (t114 * t133 + t126 * t175) * MDP(9); t167 + t114 * MDP(9) - t115 * MDP(10) + (-pkin(2) * t112 - t111 * t132) * MDP(11) + (-pkin(2) * t113 + t111 * t130) * MDP(12) + t147 * MDP(13) - pkin(2) * t159 + (t103 * t95 + t92 * t99) * MDP(27) + (t103 * t96 - t93 * t99) * MDP(28) + (t99 * MDP(20) + t153) * t124 + (t166 + (t103 * MDP(20) + MDP(7) + t152) * t138) * t131 + ((-t112 * t132 + t113 * t130) * MDP(13) + t147 * MDP(14) + (MDP(11) * t130 + MDP(12) * t132) * t168) * qJ(3) + (MDP(18) * t168 + t101 * MDP(20) + t141 - t155) * t118 + (-MDP(17) * t168 + t100 * MDP(15) - t99 * MDP(16) + t101 * MDP(21) + (-t96 * MDP(23) - t99 * MDP(25) + MDP(27) * t86) * t134 + (t99 * MDP(24) + MDP(28) * t86 - t165 + t172) * t137) * t119; 0.2e1 * t124 * t157 + MDP(8) + (-0.2e1 * t160 + 0.2e1 * t161 + t173) * pkin(2) + (MDP(22) * t129 + MDP(15) - 0.2e1 * t149) * t119 ^ 2 + (t124 * t180 + t151 + 0.2e1 * (-MDP(16) + t146) * t119) * t118 + (t118 * t92 + t134 * t170) * t179 + (-t118 * t93 + t137 * t170) * t178 + (MDP(14) * qJ(3) + t181) * (t130 ^ 2 + t132 ^ 2) * qJ(3); t112 * MDP(11) + t113 * MDP(12) + t143 * t99 + t153 + t159; t143 * t118 + t157 + t160 - t161 - t173; MDP(14); t154 - t131 * t158 + t88 * MDP(20) - t89 * MDP(21) + t96 * t156 + (-t134 * t95 + t137 * t96) * MDP(23) + (-pkin(4) * t95 - t137 * t86) * MDP(27) + (-pkin(4) * t96 + t134 * t86) * MDP(28) + t140 * t99; -t152 - t143 * t103 + t140 * t118 + (MDP(17) + t137 * t156 + (-t128 + t129) * MDP(23) + t144 * pkin(4)) * t119; 0; MDP(22) * t128 + 0.2e1 * pkin(4) * t145 + MDP(19) + 0.2e1 * t149; t141; MDP(27) * t92 - MDP(28) * t93 + t146 * t119 + t151; t145; t142; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
