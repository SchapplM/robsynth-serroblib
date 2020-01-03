% Calculate joint inertia matrix for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPR10_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR10_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S5RRRPR10_inertiaJ_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:29:54
% EndTime: 2019-12-31 21:29:56
% DurationCPUTime: 0.65s
% Computational Cost: add. (960->167), mult. (2182->264), div. (0->0), fcn. (2380->10), ass. (0->80)
t133 = sin(pkin(5));
t176 = 0.2e1 * t133;
t137 = sin(qJ(3));
t140 = cos(qJ(3));
t175 = -(MDP(16) * t137 + MDP(17) * t140) * pkin(8) + t137 * MDP(13) + t140 * MDP(14);
t174 = 0.2e1 * MDP(16);
t173 = 2 * MDP(18);
t172 = 2 * MDP(25);
t171 = 2 * MDP(26);
t138 = sin(qJ(2));
t170 = pkin(1) * t138;
t141 = cos(qJ(2));
t169 = pkin(1) * t141;
t168 = -qJ(4) - pkin(8);
t132 = sin(pkin(10));
t134 = cos(pkin(10));
t135 = cos(pkin(5));
t165 = t133 * t141;
t153 = pkin(7) * t165;
t111 = t153 + (pkin(8) + t170) * t135;
t112 = (-pkin(2) * t141 - pkin(8) * t138 - pkin(1)) * t133;
t100 = -t137 * t111 + t140 * t112;
t166 = t133 * t138;
t115 = t135 * t137 + t140 * t166;
t93 = -pkin(3) * t165 - t115 * qJ(4) + t100;
t101 = t140 * t111 + t137 * t112;
t114 = -t135 * t140 + t137 * t166;
t97 = -t114 * qJ(4) + t101;
t90 = t132 * t93 + t134 * t97;
t122 = t168 * t140;
t151 = t168 * t137;
t106 = -t132 * t122 - t134 * t151;
t120 = t132 * t140 + t134 * t137;
t167 = t106 * t120;
t164 = t135 * MDP(8);
t163 = t138 * MDP(6);
t103 = -t132 * t114 + t134 * t115;
t136 = sin(qJ(5));
t139 = cos(qJ(5));
t98 = t136 * t103 + t139 * t165;
t162 = t98 * MDP(21);
t161 = t98 * MDP(23);
t99 = t139 * t103 - t136 * t165;
t160 = t99 * MDP(20);
t159 = t99 * MDP(22);
t158 = MDP(15) * t141;
t157 = MDP(20) * t136;
t102 = t134 * t114 + t132 * t115;
t156 = t102 * MDP(24);
t155 = t115 * MDP(11);
t119 = t132 * t137 - t134 * t140;
t154 = t119 * MDP(24);
t128 = -t140 * pkin(3) - pkin(2);
t152 = t136 * t139 * MDP(21);
t89 = -t132 * t97 + t134 * t93;
t150 = t115 * MDP(13) - t114 * MDP(14);
t147 = MDP(25) * t139 - MDP(26) * t136;
t146 = MDP(25) * t136 + MDP(26) * t139;
t123 = pkin(7) * t166;
t110 = t123 + (-pkin(2) - t169) * t135;
t145 = (MDP(22) * t139 - MDP(23) * t136) * t120;
t104 = t114 * pkin(3) + t110;
t144 = t136 * MDP(22) + t139 * MDP(23) - t146 * (t132 * pkin(3) + pkin(9));
t88 = -pkin(9) * t165 + t90;
t91 = t102 * pkin(4) - t103 * pkin(9) + t104;
t85 = -t136 * t88 + t139 * t91;
t86 = t136 * t91 + t139 * t88;
t143 = MDP(25) * t85 - MDP(26) * t86 + t156 + t159 - t161;
t131 = t139 ^ 2;
t130 = t136 ^ 2;
t129 = t133 ^ 2;
t127 = -t134 * pkin(3) - pkin(4);
t117 = t135 * t170 + t153;
t116 = t135 * t169 - t123;
t108 = -t134 * t122 + t132 * t151;
t105 = t119 * pkin(4) - t120 * pkin(9) + t128;
t95 = t136 * t105 + t139 * t108;
t94 = t139 * t105 - t136 * t108;
t87 = pkin(4) * t165 - t89;
t1 = [(t104 ^ 2 + t89 ^ 2 + t90 ^ 2) * MDP(19) + t129 * t138 ^ 2 * MDP(4) + MDP(1) + (t160 - 0.2e1 * t162) * t99 + (t163 * t176 + t164) * t135 + (-0.2e1 * t114 * MDP(12) + t155) * t115 + (t156 + 0.2e1 * t159 - 0.2e1 * t161) * t102 + ((0.2e1 * MDP(5) * t138 + t158) * t129 + (MDP(7) * t135 - t150) * t176) * t141 + 0.2e1 * (t116 * t135 + t129 * t169) * MDP(9) + 0.2e1 * (-t117 * t135 - t129 * t170) * MDP(10) + (-t100 * t165 + t110 * t114) * t174 + 0.2e1 * (t101 * t165 + t110 * t115) * MDP(17) + (-t90 * t102 - t89 * t103) * t173 + (t85 * t102 + t87 * t98) * t172 + (-t86 * t102 + t87 * t99) * t171; t164 + t116 * MDP(9) - t117 * MDP(10) + t137 * t155 + (-t137 * t114 + t115 * t140) * MDP(12) + (-pkin(2) * t114 - t110 * t140) * MDP(16) + (-pkin(2) * t115 + t110 * t137) * MDP(17) + (-t102 * t108 + t103 * t106) * MDP(18) + (t104 * t128 - t106 * t89 + t108 * t90) * MDP(19) + (t102 * t94 + t106 * t98) * MDP(25) + (-t102 * t95 + t106 * t99) * MDP(26) + (-MDP(18) * t90 + t143) * t119 + (-t89 * MDP(18) + (-MDP(21) * t99 - MDP(23) * t102 + MDP(25) * t87) * t136 + (MDP(22) * t102 + MDP(26) * t87 + t160 - t162) * t139) * t120 + (t163 + (MDP(7) - t175) * t141) * t133; MDP(8) + pkin(2) * t140 * t174 + (t106 ^ 2 + t108 ^ 2 + t128 ^ 2) * MDP(19) + (MDP(20) * t131 - 0.2e1 * t152) * t120 ^ 2 + (MDP(11) * t137 + 0.2e1 * MDP(12) * t140 - 0.2e1 * MDP(17) * pkin(2)) * t137 + (0.2e1 * t145 + t154) * t119 + (-t108 * t119 + t167) * t173 + (t119 * t94 + t136 * t167) * t172 + (-t119 * t95 + t139 * t167) * t171; -t133 * t158 + t100 * MDP(16) - t101 * MDP(17) + t99 * t157 + (-t136 * t98 + t139 * t99) * MDP(21) + (t127 * t98 - t139 * t87) * MDP(25) + (t127 * t99 + t136 * t87) * MDP(26) + t144 * t102 + ((-t102 * t132 - t103 * t134) * MDP(18) + (t132 * t90 + t134 * t89) * MDP(19)) * pkin(3) + t150; -t147 * t106 + t144 * t119 + (t139 * t157 + (-t130 + t131) * MDP(21) + t146 * t127) * t120 + ((-t119 * t132 - t120 * t134) * MDP(18) + (-t106 * t134 + t108 * t132) * MDP(19)) * pkin(3) + t175; 0.2e1 * t152 + t130 * MDP(20) + MDP(15) + (t132 ^ 2 + t134 ^ 2) * MDP(19) * pkin(3) ^ 2 - 0.2e1 * t147 * t127; t104 * MDP(19) + t147 * t102; t128 * MDP(19) + t147 * t119; 0; MDP(19); t143; t94 * MDP(25) - t95 * MDP(26) + t145 + t154; t144; t147; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
