% Calculate joint inertia matrix for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR15_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR15_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR15_inertiaJ_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:25
% EndTime: 2019-12-31 18:37:26
% DurationCPUTime: 0.35s
% Computational Cost: add. (274->103), mult. (519->155), div. (0->0), fcn. (472->6), ass. (0->53)
t116 = qJ(4) * MDP(17);
t82 = sin(pkin(8));
t83 = cos(pkin(8));
t106 = t82 ^ 2 + t83 ^ 2;
t115 = t106 * t116;
t87 = cos(qJ(3));
t114 = 0.2e1 * t87;
t113 = -2 * MDP(19);
t112 = 2 * MDP(24);
t111 = pkin(7) * t87;
t110 = (pkin(1) * MDP(6));
t88 = -pkin(1) - pkin(6);
t109 = t82 * t88;
t85 = sin(qJ(3));
t108 = t85 * t88;
t107 = pkin(7) + qJ(4);
t72 = t85 * pkin(3) - t87 * qJ(4) + qJ(2);
t62 = t83 * t108 + t82 * t72;
t80 = t85 ^ 2;
t81 = t87 ^ 2;
t105 = -t80 - t81;
t104 = pkin(3) * MDP(17);
t103 = MDP(17) * t85;
t84 = sin(qJ(5));
t86 = cos(qJ(5));
t70 = t84 * t82 - t86 * t83;
t66 = t70 * t87;
t102 = t66 * MDP(18);
t101 = t70 * MDP(23);
t100 = t82 * MDP(15);
t99 = t83 * MDP(14);
t98 = t88 * MDP(17);
t97 = t106 * MDP(16);
t71 = t86 * t82 + t84 * t83;
t94 = -t82 * MDP(14) - t83 * MDP(15);
t64 = t71 * t87;
t93 = -t66 * MDP(20) - t64 * MDP(21);
t92 = -t100 + t99 + t104;
t73 = t107 * t82;
t74 = t107 * t83;
t91 = t71 * MDP(20) - t70 * MDP(21) + (-t86 * t73 - t84 * t74) * MDP(23) - (-t84 * t73 + t86 * t74) * MDP(24);
t90 = t71 * MDP(24) + t101 - t92;
t77 = -t83 * pkin(4) - pkin(3);
t69 = (pkin(4) * t82 - t88) * t87;
t68 = t83 * t72;
t65 = t70 * t85;
t63 = t71 * t85;
t61 = -t82 * t108 + t68;
t58 = -t82 * t111 + t62;
t57 = -t83 * t111 + t68 + (pkin(4) - t109) * t85;
t56 = t84 * t57 + t86 * t58;
t55 = t86 * t57 - t84 * t58;
t1 = [MDP(1) + t81 * MDP(7) + (t81 * t88 ^ 2 + t61 ^ 2 + t62 ^ 2) * MDP(17) + t80 * MDP(22) - (t64 * t113 - t102) * t66 + ((-2 * MDP(4) + t110) * pkin(1)) + (MDP(13) * t114 + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2) + 0.2e1 * (qJ(2) * MDP(12) - t87 * MDP(8) + t93) * t85 + 0.2e1 * (-t81 * t109 + t61 * t85) * MDP(14) + 0.2e1 * (-t81 * t88 * t83 - t62 * t85) * MDP(15) + 0.2e1 * (t55 * t85 + t69 * t64) * MDP(23) + (-t56 * t85 - t69 * t66) * t112 + (-t61 * t83 - t62 * t82) * MDP(16) * t114; MDP(4) - t110 + t81 * t98 + (-t63 * t85 - t87 * t64) * MDP(23) + (t65 * t85 + t87 * t66) * MDP(24) + (t105 * MDP(15) + t62 * t103) * t83 + (t105 * MDP(14) - t61 * t103) * t82; MDP(6) + (t106 * t80 + t81) * MDP(17); -t71 * t102 + (-t71 * t64 + t66 * t70) * MDP(19) + (t77 * t64 + t69 * t70) * MDP(23) + (-t77 * t66 + t69 * t71) * MDP(24) + (MDP(9) + t94 * pkin(3) + (MDP(12) + t92) * t88) * t87 + (-t88 * MDP(13) + t94 * qJ(4) - MDP(10) + t91) * t85 + (MDP(16) + t116) * (-t61 * t82 + t62 * t83); (-MDP(13) + t97 + t115) * t85 + (MDP(12) - t90) * t87; 0.2e1 * t77 * t101 + MDP(11) + (-0.2e1 * t100 + 0.2e1 * t99 + t104) * pkin(3) + (MDP(18) * t71 + t77 * t112 + t70 * t113) * t71 + (0.2e1 * t97 + t115) * qJ(4); t64 * MDP(23) - t66 * MDP(24) + (-t94 - t98) * t87; -t87 * MDP(17); t90; MDP(17); t85 * MDP(22) + t55 * MDP(23) - t56 * MDP(24) + t93; -t63 * MDP(23) + t65 * MDP(24); t91; 0; MDP(22);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
