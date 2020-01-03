% Calculate joint inertia matrix for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP9_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP9_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRRRP9_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:06:06
% EndTime: 2019-12-31 22:06:08
% DurationCPUTime: 0.76s
% Computational Cost: add. (658->166), mult. (1243->232), div. (0->0), fcn. (1173->6), ass. (0->64)
t115 = sin(qJ(3));
t118 = cos(qJ(3));
t129 = MDP(23) + MDP(25);
t151 = MDP(24) - MDP(27);
t144 = -pkin(8) - pkin(7);
t100 = t144 * t118;
t114 = sin(qJ(4));
t117 = cos(qJ(4));
t127 = t144 * t115;
t84 = -t114 * t100 - t117 * t127;
t85 = -t117 * t100 + t114 * t127;
t135 = t117 * t118;
t95 = t114 * t115 - t135;
t96 = t114 * t118 + t117 * t115;
t124 = t96 * MDP(20) - t95 * MDP(21) - t129 * t84 - t151 * t85;
t157 = t124 - (t115 * MDP(16) + t118 * MDP(17)) * pkin(7) + t115 * MDP(13) + t118 * MDP(14);
t131 = t114 * MDP(24);
t154 = (t117 * MDP(23) - t131) * pkin(3);
t116 = sin(qJ(2));
t90 = t96 * t116;
t137 = t115 * t116;
t91 = -t114 * t137 + t116 * t135;
t138 = t91 * MDP(20) - t90 * MDP(21);
t119 = cos(qJ(2));
t140 = pkin(8) * t116;
t143 = pkin(6) * t115;
t99 = -t119 * pkin(2) - t116 * pkin(7) - pkin(1);
t94 = t118 * t99;
t77 = -t118 * t140 + t94 + (-pkin(3) - t143) * t119;
t141 = pkin(6) * t119;
t128 = t118 * t141;
t79 = t128 + (t99 - t140) * t115;
t71 = -t114 * t79 + t117 * t77;
t152 = t71 * MDP(23) + t138;
t149 = -2 * MDP(19);
t148 = 0.2e1 * MDP(24);
t147 = 2 * MDP(25);
t146 = 2 * MDP(26);
t145 = 2 * MDP(27);
t142 = pkin(6) * t118;
t139 = t119 * pkin(4);
t72 = t114 * t77 + t117 * t79;
t136 = t115 * t118;
t134 = t119 * qJ(5);
t132 = t91 * MDP(18);
t98 = pkin(3) * t137 + t116 * pkin(6);
t130 = MDP(15) + MDP(22);
t107 = -t118 * pkin(3) - pkin(2);
t126 = MDP(12) * t136;
t125 = pkin(4) * t147 + MDP(22);
t123 = t118 * MDP(13) - t115 * MDP(14);
t112 = t118 ^ 2;
t111 = t116 ^ 2;
t110 = t115 ^ 2;
t108 = t114 * pkin(3);
t105 = t117 * pkin(3) + pkin(4);
t103 = t108 + qJ(5);
t87 = t115 * t99 + t128;
t86 = -t115 * t141 + t94;
t76 = t95 * pkin(4) - t96 * qJ(5) + t107;
t73 = t90 * pkin(4) - t91 * qJ(5) + t98;
t70 = -t71 + t139;
t69 = -t134 + t72;
t1 = [MDP(1) - 0.2e1 * pkin(1) * t116 * MDP(10) + (t69 ^ 2 + t70 ^ 2 + t73 ^ 2) * MDP(28) + (t90 * t149 + t132) * t91 + t130 * t119 ^ 2 + (t112 * MDP(11) + MDP(4) - 0.2e1 * t126) * t111 + 0.2e1 * (pkin(1) * MDP(9) + (MDP(5) - t123) * t116 - t138) * t119 + 0.2e1 * (t111 * t143 - t86 * t119) * MDP(16) + 0.2e1 * (t111 * t142 + t87 * t119) * MDP(17) + 0.2e1 * (-t71 * t119 + t98 * t90) * MDP(23) + (t72 * t119 + t98 * t91) * t148 + (t70 * t119 + t73 * t90) * t147 + (-t69 * t90 + t70 * t91) * t146 + (-t69 * t119 - t73 * t91) * t145; t96 * t132 + (-t96 * t90 - t91 * t95) * MDP(19) + (t107 * t90 + t98 * t95) * MDP(23) + (t107 * t91 + t98 * t96) * MDP(24) + (t73 * t95 + t76 * t90) * MDP(25) + (-t69 * t95 + t70 * t96 + t84 * t91 - t85 * t90) * MDP(26) + (-t73 * t96 - t76 * t91) * MDP(27) + (t69 * t85 + t70 * t84 + t73 * t76) * MDP(28) + (-pkin(6) * MDP(10) + MDP(7) - t157) * t119 + (MDP(6) - pkin(6) * MDP(9) + MDP(11) * t136 + (-t110 + t112) * MDP(12) + (-pkin(2) * t115 - t142) * MDP(16) + (-pkin(2) * t118 + t143) * MDP(17)) * t116; MDP(8) + t110 * MDP(11) + 0.2e1 * t126 + (t76 ^ 2 + t84 ^ 2 + t85 ^ 2) * MDP(28) + (MDP(18) * t96 - 0.2e1 * t76 * MDP(27) + t107 * t148 + t84 * t146 + t95 * t149) * t96 + 0.2e1 * (t118 * MDP(16) - t115 * MDP(17)) * pkin(2) + 0.2e1 * (t107 * MDP(23) + t76 * MDP(25) - t85 * MDP(26)) * t95; t86 * MDP(16) - t87 * MDP(17) + t71 * MDP(25) + (-t103 * t90 - t105 * t91) * MDP(26) + (t69 * t103 - t70 * t105) * MDP(28) + t123 * t116 + ((-pkin(4) - t105) * MDP(25) + (-qJ(5) - t103) * MDP(27) - t154 - t130) * t119 - t151 * t72 + t152; (-t103 * t95 - t105 * t96) * MDP(26) + (t85 * t103 - t84 * t105) * MDP(28) + t157; (t103 ^ 2 + t105 ^ 2) * MDP(28) + 0.2e1 * t154 + t105 * t147 + t103 * t145 + t130; -t119 * MDP(22) - t72 * MDP(24) + (t71 - 0.2e1 * t139) * MDP(25) + (-pkin(4) * t91 - t90 * qJ(5)) * MDP(26) + (-0.2e1 * t134 + t72) * MDP(27) + (-t70 * pkin(4) + t69 * qJ(5)) * MDP(28) + t152; (-pkin(4) * t96 - t95 * qJ(5)) * MDP(26) + (-t84 * pkin(4) + t85 * qJ(5)) * MDP(28) + t124; (0.2e1 * qJ(5) + t108) * MDP(27) + (t105 * pkin(4) + t103 * qJ(5)) * MDP(28) + (t129 * t117 - t131) * pkin(3) + t125; qJ(5) * t145 + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(28) + t125; t119 * MDP(25) + t91 * MDP(26) + t70 * MDP(28); t96 * MDP(26) + t84 * MDP(28); -t105 * MDP(28) - MDP(25); -MDP(28) * pkin(4) - MDP(25); MDP(28);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
