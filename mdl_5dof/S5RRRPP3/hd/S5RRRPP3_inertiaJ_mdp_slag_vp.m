% Calculate joint inertia matrix for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRRPP3_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:42
% EndTime: 2019-12-31 20:53:43
% DurationCPUTime: 0.31s
% Computational Cost: add. (260->104), mult. (398->127), div. (0->0), fcn. (252->4), ass. (0->42)
t92 = sin(qJ(3));
t88 = t92 ^ 2;
t94 = cos(qJ(3));
t107 = t94 ^ 2 + t88;
t93 = sin(qJ(2));
t76 = t93 * pkin(1) + pkin(7);
t110 = t107 * t76;
t104 = 0.2e1 * t94;
t118 = 2 * MDP(14);
t117 = 2 * MDP(18);
t95 = cos(qJ(2));
t116 = t95 * pkin(1);
t77 = -pkin(2) - t116;
t115 = pkin(2) - t77;
t90 = pkin(3) + qJ(5);
t102 = -t92 * qJ(4) - pkin(2);
t61 = -t90 * t94 + t102;
t58 = t61 - t116;
t114 = -t58 - t61;
t63 = (pkin(4) + t76) * t92;
t87 = t94 * pkin(4);
t64 = t94 * t76 + t87;
t113 = t63 * t92 + t64 * t94;
t69 = -t94 * pkin(3) + t102;
t62 = t69 - t116;
t112 = t62 + t69;
t70 = (pkin(4) + pkin(7)) * t92;
t71 = t94 * pkin(7) + t87;
t111 = t70 * t92 + t71 * t94;
t109 = (MDP(18) + MDP(14)) * t92;
t108 = t107 * pkin(7);
t106 = MDP(17) * t92;
t105 = 0.2e1 * t92;
t103 = t92 * MDP(8) * t104 + t88 * MDP(7) + MDP(4);
t84 = t94 * qJ(4);
t101 = (-t90 * t92 + t84) * MDP(18) + (-t92 * pkin(3) + t84) * MDP(14) + t94 * MDP(10) + t92 * MDP(9);
t100 = -pkin(3) * MDP(17) + MDP(15);
t99 = (t95 * MDP(5) - t93 * MDP(6)) * pkin(1);
t98 = (qJ(4) * MDP(17) - MDP(13) + MDP(16)) * t94 + (-MDP(12) + t100) * t92;
t96 = qJ(4) ^ 2;
t80 = t94 * MDP(18);
t1 = [MDP(1) + (t107 * t76 ^ 2 + t62 ^ 2) * MDP(17) + (t58 ^ 2 + t63 ^ 2 + t64 ^ 2) * MDP(21) + 0.2e1 * t99 + t110 * t118 + t113 * t117 + (-t77 * MDP(12) + t62 * MDP(15) - t58 * MDP(20)) * t104 + (t77 * MDP(13) - t62 * MDP(16) - t58 * MDP(19)) * t105 + t103; (t108 + t110) * MDP(14) + (pkin(7) * t110 + t62 * t69) * MDP(17) + (t111 + t113) * MDP(18) + (t58 * t61 + t63 * t70 + t64 * t71) * MDP(21) + t99 + (t115 * MDP(12) + t112 * MDP(15) + t114 * MDP(20)) * t94 + (-t115 * MDP(13) - t112 * MDP(16) + t114 * MDP(19)) * t92 + t103; (t107 * pkin(7) ^ 2 + t69 ^ 2) * MDP(17) + (t61 ^ 2 + t70 ^ 2 + t71 ^ 2) * MDP(21) + t108 * t118 + t111 * t117 + (pkin(2) * MDP(12) + t69 * MDP(15) - t61 * MDP(20)) * t104 + (-pkin(2) * MDP(13) - t69 * MDP(16) - t61 * MDP(19)) * t105 + t103; t64 * MDP(19) - t63 * MDP(20) + (t64 * qJ(4) - t63 * t90) * MDP(21) + t98 * t76 + t101; t71 * MDP(19) - t70 * MDP(20) + (t71 * qJ(4) - t70 * t90) * MDP(21) + t98 * pkin(7) + t101; MDP(11) - 0.2e1 * pkin(3) * MDP(15) + (pkin(3) ^ 2 + t96) * MDP(17) + 0.2e1 * t90 * MDP(20) + (t90 ^ 2 + t96) * MDP(21) + 0.2e1 * (MDP(16) + MDP(19)) * qJ(4); t63 * MDP(21) + t76 * t106 + t109; t70 * MDP(21) + pkin(7) * t106 + t109; -t90 * MDP(21) - MDP(20) + t100; MDP(17) + MDP(21); t64 * MDP(21) + t80; t71 * MDP(21) + t80; MDP(21) * qJ(4) + MDP(19); 0; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
