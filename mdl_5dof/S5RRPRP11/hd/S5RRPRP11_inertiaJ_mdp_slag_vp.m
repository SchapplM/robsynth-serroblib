% Calculate joint inertia matrix for
% S5RRPRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRP11_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP11_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRPRP11_inertiaJ_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:14:05
% EndTime: 2019-12-31 20:14:06
% DurationCPUTime: 0.37s
% Computational Cost: add. (303->120), mult. (500->155), div. (0->0), fcn. (370->4), ass. (0->45)
t107 = pkin(3) + pkin(6);
t79 = sin(qJ(4));
t74 = t79 ^ 2;
t81 = cos(qJ(4));
t76 = t81 ^ 2;
t71 = t74 + t76;
t106 = MDP(25) * t71;
t104 = 2 * MDP(22);
t83 = -pkin(2) - pkin(7);
t80 = sin(qJ(2));
t103 = pkin(4) * t80;
t102 = t80 * t83;
t82 = cos(qJ(2));
t93 = -qJ(3) * t80 - pkin(1);
t61 = t83 * t82 + t93;
t68 = t107 * t80;
t58 = t81 * t61 + t79 * t68;
t69 = t107 * t82;
t101 = pkin(2) * MDP(14);
t100 = qJ(5) * t80;
t88 = pkin(4) * t79 - qJ(5) * t81;
t65 = qJ(3) + t88;
t99 = MDP(25) * t65;
t98 = MDP(25) * t83;
t97 = pkin(6) ^ 2 * MDP(14);
t96 = MDP(14) * qJ(3);
t95 = -MDP(21) + MDP(24);
t94 = t81 * t79 * MDP(16);
t92 = t61 * t79 - t81 * t68;
t91 = MDP(12) - t101;
t90 = -MDP(25) * pkin(4) - MDP(22);
t89 = (MDP(20) + MDP(22)) * t81;
t67 = pkin(4) * t81 + qJ(5) * t79;
t87 = -t79 * MDP(17) - t81 * MDP(18);
t86 = -MDP(20) * t92 - t58 * MDP(21);
t85 = t95 * t79 + t89;
t77 = t82 ^ 2;
t75 = t80 ^ 2;
t66 = -pkin(2) * t82 + t93;
t64 = t71 * t83;
t59 = t67 * t82 + t69;
t56 = t92 - t103;
t55 = t100 + t58;
t54 = t55 * t79 - t56 * t81;
t1 = [MDP(1) + (t55 ^ 2 + t56 ^ 2 + t59 ^ 2) * MDP(25) + t66 ^ 2 * MDP(14) + (t74 * MDP(15) + 0.2e1 * t94 + t97) * t77 + (MDP(19) + MDP(4) + t97) * t75 + 0.2e1 * (t75 + t77) * MDP(11) * pkin(6) + 0.2e1 * (pkin(1) * MDP(9) + t66 * MDP(12) + (-t55 * t81 - t56 * t79) * MDP(23) + (t81 * MDP(20) - t79 * MDP(21)) * t69 + (t81 * MDP(22) + t79 * MDP(24)) * t59) * t82 + 0.2e1 * (-pkin(1) * MDP(10) - t66 * MDP(13) + (MDP(5) + t87) * t82 - t56 * MDP(22) + t55 * MDP(24) + t86) * t80; t59 * t99 - t54 * MDP(23) + (-MDP(11) * pkin(2) + MDP(6)) * t80 + t89 * t102 + (MDP(7) + qJ(3) * MDP(11) + (t74 - t76) * MDP(16)) * t82 + (-t56 * t98 + t80 * MDP(17) + t69 * MDP(21) - t59 * MDP(24) + (MDP(20) * qJ(3) + MDP(22) * t65) * t82) * t81 + (-t82 * t81 * MDP(15) - t80 * MDP(18) + t69 * MDP(20) + (-qJ(3) * t82 - t102) * MDP(21) + t59 * MDP(22) + (t65 * t82 + t102) * MDP(24) + t55 * t98) * t79 + ((-MDP(10) + MDP(13) + t96) * t82 + (-MDP(9) + t91) * t80) * pkin(6); -0.2e1 * t94 - 0.2e1 * t64 * MDP(23) + t76 * MDP(15) + MDP(8) + t83 ^ 2 * t106 + (-0.2e1 * MDP(24) * t81 + t79 * t104 + t99) * t65 + (-0.2e1 * MDP(12) + t101) * pkin(2) + (0.2e1 * MDP(20) * t79 + 0.2e1 * MDP(21) * t81 + 0.2e1 * MDP(13) + t96) * qJ(3); MDP(25) * t54 + (MDP(14) * pkin(6) + MDP(11) + t85) * t80; -MDP(23) * t71 + MDP(25) * t64 + t91; MDP(14) + t106; t80 * MDP(19) + (-t92 + 0.2e1 * t103) * MDP(22) + (0.2e1 * t100 + t58) * MDP(24) + (-pkin(4) * t56 + qJ(5) * t55) * MDP(25) + (t88 * MDP(23) + t87) * t82 + t86; t81 * MDP(17) - t79 * MDP(18) - t67 * MDP(23) + ((MDP(20) - t90) * t81 + (MDP(25) * qJ(5) + t95) * t79) * t83; MDP(25) * t67 + t85; MDP(19) + pkin(4) * t104 + 0.2e1 * qJ(5) * MDP(24) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(25); -t79 * t82 * MDP(23) - t80 * MDP(22) + MDP(25) * t56; (MDP(23) - t98) * t81; -t81 * MDP(25); t90; MDP(25);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
