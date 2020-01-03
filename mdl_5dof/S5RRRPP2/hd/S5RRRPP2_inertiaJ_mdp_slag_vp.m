% Calculate joint inertia matrix for
% S5RRRPP2
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
%   see S5RRRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RRRPP2_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:59
% EndTime: 2019-12-31 20:52:01
% DurationCPUTime: 0.32s
% Computational Cost: add. (265->107), mult. (401->129), div. (0->0), fcn. (256->4), ass. (0->39)
t83 = sin(qJ(3));
t81 = t83 ^ 2;
t85 = cos(qJ(3));
t100 = t85 ^ 2 + t81;
t84 = sin(qJ(2));
t67 = pkin(1) * t84 + pkin(7);
t103 = t100 * t67;
t97 = 0.2e1 * t85;
t109 = t85 * pkin(3) + t83 * qJ(4);
t108 = 2 * MDP(15);
t107 = 2 * MDP(20);
t86 = cos(qJ(2));
t68 = -pkin(1) * t86 - pkin(2);
t106 = pkin(2) - t68;
t54 = t68 - t109;
t79 = t85 * pkin(4);
t53 = -t54 + t79;
t96 = pkin(2) + t109;
t55 = t79 + t96;
t105 = t53 + t55;
t104 = -t54 + t96;
t102 = t85 * MDP(18) + t83 * MDP(19);
t101 = t100 * pkin(7);
t99 = qJ(4) * t85;
t98 = 0.2e1 * t83;
t95 = t83 * MDP(8) * t97 + t81 * MDP(7) + MDP(4);
t87 = pkin(3) + pkin(4);
t94 = (t83 * t87 - t99) * MDP(20) + (-pkin(3) * t83 + t99) * MDP(15) + t85 * MDP(10) + t83 * MDP(9);
t93 = -pkin(3) * MDP(17) - MDP(14);
t92 = (t86 * MDP(5) - t84 * MDP(6)) * pkin(1);
t91 = (MDP(17) * qJ(4) - MDP(13) + MDP(16)) * t85 + (-MDP(12) + t93) * t83;
t89 = qJ(4) ^ 2;
t74 = t83 * qJ(5);
t70 = t83 * MDP(15);
t62 = (pkin(7) - qJ(5)) * t85;
t61 = pkin(7) * t83 - t74;
t57 = (-qJ(5) + t67) * t85;
t56 = t67 * t83 - t74;
t1 = [MDP(1) + (t100 * t67 ^ 2 + t54 ^ 2) * MDP(17) + (t53 ^ 2 + t56 ^ 2 + t57 ^ 2) * MDP(21) + 0.2e1 * t92 + (-t68 * MDP(12) - t54 * MDP(14) + t53 * MDP(18)) * t97 + (t68 * MDP(13) - t54 * MDP(16) + t53 * MDP(19)) * t98 + t103 * t108 + (-t56 * t83 - t57 * t85) * t107 + t95; (t101 + t103) * MDP(15) + (pkin(7) * t103 - t54 * t96) * MDP(17) + (t53 * t55 + t56 * t61 + t57 * t62) * MDP(21) + t92 + (t106 * MDP(12) + t104 * MDP(14) + t105 * MDP(18) + (-t57 - t62) * MDP(20)) * t85 + (-t106 * MDP(13) + t104 * MDP(16) + t105 * MDP(19) + (-t56 - t61) * MDP(20)) * t83 + t95; (t100 * pkin(7) ^ 2 + t96 ^ 2) * MDP(17) + (t55 ^ 2 + t61 ^ 2 + t62 ^ 2) * MDP(21) + (MDP(12) * pkin(2) + MDP(14) * t96 + MDP(18) * t55) * t97 + (-MDP(13) * pkin(2) + MDP(16) * t96 + MDP(19) * t55) * t98 + t101 * t108 + (-t61 * t83 - t62 * t85) * t107 + t95; -t56 * MDP(18) + t57 * MDP(19) + (qJ(4) * t57 - t56 * t87) * MDP(21) + t91 * t67 + t94; -t61 * MDP(18) + t62 * MDP(19) + (qJ(4) * t62 - t61 * t87) * MDP(21) + t91 * pkin(7) + t94; MDP(11) + 0.2e1 * pkin(3) * MDP(14) + (pkin(3) ^ 2 + t89) * MDP(17) + 0.2e1 * t87 * MDP(18) + (t87 ^ 2 + t89) * MDP(21) + 0.2e1 * (MDP(16) + MDP(19)) * qJ(4); MDP(21) * t56 + t70 + (MDP(17) * t67 - MDP(20)) * t83; MDP(21) * t61 + t70 + (MDP(17) * pkin(7) - MDP(20)) * t83; -MDP(21) * t87 - MDP(18) + t93; MDP(17) + MDP(21); t53 * MDP(21) + t102; t55 * MDP(21) + t102; 0; 0; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
