% Calculate joint inertia matrix for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP4_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:24:51
% EndTime: 2021-01-15 22:24:53
% DurationCPUTime: 0.33s
% Computational Cost: add. (589->109), mult. (1051->155), div. (0->0), fcn. (1080->6), ass. (0->47)
t107 = MDP(18) + MDP(22);
t120 = MDP(19) - MDP(24);
t115 = 2 * MDP(18);
t119 = -2 * MDP(19);
t118 = MDP(21) + MDP(25);
t96 = cos(qJ(2));
t89 = -t96 * pkin(2) - pkin(1);
t117 = 0.2e1 * t89;
t116 = 0.2e1 * t96;
t114 = 2 * MDP(22);
t113 = 2 * MDP(24);
t112 = pkin(6) + pkin(7);
t93 = sin(qJ(3));
t111 = pkin(2) * t93;
t95 = cos(qJ(3));
t88 = t95 * pkin(2) + pkin(3);
t91 = sin(pkin(8));
t92 = cos(pkin(8));
t75 = t92 * t111 + t91 * t88;
t94 = sin(qJ(2));
t79 = t93 * t94 - t95 * t96;
t80 = t93 * t96 + t95 * t94;
t67 = t92 * t79 + t91 * t80;
t68 = -t91 * t79 + t92 * t80;
t71 = t79 * pkin(3) + t89;
t60 = t67 * pkin(4) - t68 * qJ(5) + t71;
t110 = t60 * MDP(25);
t109 = t71 * MDP(21);
t108 = t95 * MDP(16);
t90 = t91 * pkin(3);
t105 = -t90 - t75;
t104 = t112 * t94;
t103 = t112 * t96;
t82 = t92 * t88;
t74 = -t91 * t111 + t82;
t100 = -t93 * t103 - t95 * t104;
t99 = -t95 * t103 + t93 * t104;
t66 = -t79 * qJ(4) - t99;
t98 = -t80 * qJ(4) + t100;
t62 = t91 * t66 - t92 * t98;
t64 = t92 * t66 + t91 * t98;
t101 = t80 * MDP(13) - t79 * MDP(14) + t100 * MDP(16) + t99 * MDP(17) - t107 * t62 - t120 * t64;
t85 = t92 * pkin(3) + pkin(4);
t84 = t90 + qJ(5);
t73 = -pkin(4) - t74;
t72 = qJ(5) + t75;
t1 = [MDP(1) + pkin(1) * MDP(9) * t116 + t79 * MDP(16) * t117 + (0.2e1 * t68 * MDP(19) + t67 * t115 + t109) * t71 + (-0.2e1 * t68 * MDP(24) + t67 * t114 + t110) * t60 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t94 + MDP(5) * t116) * t94 + (MDP(11) * t80 - 0.2e1 * t79 * MDP(12) + MDP(17) * t117) * t80 + t118 * (t62 ^ 2 + t64 ^ 2) + 0.2e1 * (MDP(20) + MDP(23)) * (t62 * t68 - t64 * t67); t94 * MDP(6) + t96 * MDP(7) + (-t75 * t67 - t74 * t68) * MDP(20) + (-t62 * t74 + t64 * t75) * MDP(21) + (-t72 * t67 + t73 * t68) * MDP(23) + (t62 * t73 + t64 * t72) * MDP(25) + (-t96 * MDP(10) - t94 * MDP(9)) * pkin(6) + t101; MDP(8) + MDP(15) + (t74 ^ 2 + t75 ^ 2) * MDP(21) + (t72 ^ 2 + t73 ^ 2) * MDP(25) + 0.2e1 * (-t93 * MDP(17) + t108) * pkin(2) + t74 * t115 + t75 * t119 - 0.2e1 * t73 * MDP(22) + t72 * t113; (-t84 * t67 - t85 * t68) * MDP(23) + (-t62 * t85 + t64 * t84) * MDP(25) + ((-t67 * t91 - t68 * t92) * MDP(20) + (-t62 * t92 + t64 * t91) * MDP(21)) * pkin(3) + t101; MDP(15) + t82 * MDP(18) + t105 * MDP(19) + (0.2e1 * pkin(4) + t82) * MDP(22) + (0.2e1 * qJ(5) - t105) * MDP(24) + (t72 * t84 - t73 * t85) * MDP(25) + (t75 * t91 * MDP(21) + (t74 * MDP(21) + t107) * t92) * pkin(3) + (t108 + (-t107 * t91 - MDP(17)) * t93) * pkin(2); MDP(15) + (t84 ^ 2 + t85 ^ 2) * MDP(25) + t85 * t114 + t84 * t113 + (t92 * t115 + t91 * t119 + (t91 ^ 2 + t92 ^ 2) * MDP(21) * pkin(3)) * pkin(3); t107 * t67 + t120 * t68 + t109 + t110; 0; 0; t118; t68 * MDP(23) + t62 * MDP(25); t73 * MDP(25) - MDP(22); -t85 * MDP(25) - MDP(22); 0; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
