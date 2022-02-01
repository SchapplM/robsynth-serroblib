% Calculate joint inertia matrix for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRRRP2_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:49:29
% EndTime: 2022-01-20 11:49:29
% DurationCPUTime: 0.25s
% Computational Cost: add. (454->103), mult. (765->126), div. (0->0), fcn. (733->6), ass. (0->58)
t108 = sin(qJ(4));
t109 = sin(qJ(3));
t111 = cos(qJ(4));
t112 = cos(qJ(3));
t91 = t108 * t109 - t111 * t112;
t92 = t108 * t112 + t111 * t109;
t137 = t91 * MDP(21) + t92 * MDP(22);
t117 = t112 * MDP(12) - t109 * MDP(13);
t106 = t111 * pkin(3);
t100 = t106 + pkin(4);
t143 = pkin(3) * t108;
t147 = t112 * MDP(10) + t109 * MDP(9) + (-t100 * t92 - t91 * t143) * MDP(23);
t146 = 0.2e1 * MDP(21);
t145 = 0.2e1 * MDP(23);
t144 = t91 * pkin(4);
t113 = cos(qJ(2));
t142 = t113 * pkin(1);
t110 = sin(qJ(2));
t99 = t110 * pkin(1) + pkin(7);
t89 = (-pkin(8) - t99) * t109;
t107 = t112 * pkin(8);
t90 = t112 * t99 + t107;
t125 = -t108 * t90 + t111 * t89;
t133 = t92 * qJ(5);
t68 = t125 - t133;
t122 = -t108 * t89 - t111 * t90;
t134 = t91 * qJ(5);
t69 = -t122 - t134;
t140 = -t68 * t92 - t69 * t91;
t95 = (-pkin(7) - pkin(8)) * t109;
t96 = t112 * pkin(7) + t107;
t124 = -t108 * t96 + t111 * t95;
t72 = t124 - t133;
t121 = -t108 * t95 - t111 * t96;
t73 = -t121 - t134;
t139 = -t72 * t92 - t73 * t91;
t102 = -t112 * pkin(3) - pkin(2);
t94 = t102 - t142;
t80 = t94 + t144;
t81 = t102 + t144;
t138 = t80 + t81;
t136 = t92 * MDP(16) - t91 * MDP(17);
t135 = MDP(24) * pkin(4);
t132 = t102 + t94;
t131 = t80 * MDP(24);
t130 = t81 * MDP(24);
t129 = t100 * MDP(24);
t127 = t111 * MDP(19);
t123 = MDP(4) + (MDP(14) * t92 - 0.2e1 * MDP(15) * t91) * t92 + (MDP(7) * t109 + 0.2e1 * MDP(8) * t112) * t109;
t120 = t125 * MDP(19) + t122 * MDP(20) + t68 * MDP(21) - t69 * MDP(22) + t136;
t119 = t124 * MDP(19) + t121 * MDP(20) + t72 * MDP(21) - t73 * MDP(22) + t136;
t118 = 0.2e1 * t137;
t116 = -t109 * MDP(12) - t112 * MDP(13);
t115 = (t113 * MDP(5) - t110 * MDP(6)) * pkin(1);
t114 = 0.2e1 * t91 * MDP(19) + 0.2e1 * t92 * MDP(20);
t101 = -pkin(2) - t142;
t83 = t92 * pkin(4) * MDP(23);
t1 = [MDP(1) + t140 * t145 + (t68 ^ 2 + t69 ^ 2) * MDP(24) + t94 * t114 + (t118 + t131) * t80 + t123 - 0.2e1 * t117 * t101 + 0.2e1 * t115; (t139 + t140) * MDP(23) + (t68 * t72 + t69 * t73 + t80 * t81) * MDP(24) + t115 + (t132 * MDP(20) + t138 * MDP(22)) * t92 + (t132 * MDP(19) + t138 * MDP(21)) * t91 + t123 + t117 * (pkin(2) - t101); t139 * t145 + (t72 ^ 2 + t73 ^ 2) * MDP(24) + (t118 + t130) * t81 + t102 * t114 + 0.2e1 * t117 * pkin(2) + t123; (t68 * t100 + t69 * t143) * MDP(24) + t116 * t99 + t120 + t147; (t72 * t100 + t73 * t143) * MDP(24) + t116 * pkin(7) + t119 + t147; MDP(11) + MDP(18) + (t146 + t129) * t100 + (0.2e1 * t127 + (MDP(24) * t143 - 0.2e1 * MDP(20) - 0.2e1 * MDP(22)) * t108) * pkin(3); t68 * t135 + t120 - t83; t72 * t135 + t119 - t83; MDP(18) + (0.2e1 * pkin(4) + t106) * MDP(21) + pkin(4) * t129 + (t127 + (-MDP(20) - MDP(22)) * t108) * pkin(3); MDP(18) + (t146 + t135) * pkin(4); t131 + t137; t130 + t137; 0; 0; MDP(24);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
