% Calculate joint inertia matrix for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:13
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(4,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_inertiaJ_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S5RRPRR1_inertiaJ_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:13:20
% EndTime: 2021-01-15 21:13:22
% DurationCPUTime: 0.22s
% Computational Cost: add. (209->75), mult. (417->111), div. (0->0), fcn. (382->6), ass. (0->45)
t79 = sin(qJ(5));
t82 = cos(qJ(5));
t101 = t79 * MDP(24) + t82 * MDP(25);
t108 = 2 * pkin(1);
t80 = sin(qJ(4));
t83 = cos(qJ(4));
t85 = pkin(1) + pkin(2);
t93 = t82 * MDP(27) - t79 * MDP(28);
t90 = MDP(20) + t93;
t107 = (-MDP(21) * t80 + t90 * t83) * t85;
t84 = cos(qJ(2));
t68 = t85 * t84;
t106 = -0.2e1 * t68;
t81 = sin(qJ(2));
t64 = t80 * t81 - t83 * t84;
t105 = pkin(4) * t64;
t102 = pkin(3) + qJ(3);
t66 = t102 * t81;
t67 = t102 * t84;
t55 = t83 * t66 + t80 * t67;
t52 = t55 * t79;
t104 = t55 * t82;
t103 = t79 * t82;
t99 = t64 * MDP(26);
t65 = t80 * t84 + t83 * t81;
t98 = t65 * MDP(21);
t97 = t81 * MDP(12);
t96 = MDP(23) * t103;
t75 = t79 ^ 2;
t95 = t75 * MDP(22) + MDP(19) + 0.2e1 * t96;
t94 = MDP(24) * t82 - MDP(25) * t79;
t92 = -MDP(27) * t79 - MDP(28) * t82;
t70 = t80 * t85 + pkin(4);
t91 = -t65 * t83 * t85 - t64 * t70;
t56 = -t80 * t66 + t83 * t67;
t77 = t82 ^ 2;
t89 = -t55 * MDP(20) - t56 * MDP(21) + ((-t75 + t77) * MDP(23) + MDP(22) * t103 + MDP(17)) * t65 + (-MDP(18) + t101) * t64;
t87 = pkin(1) ^ 2;
t86 = qJ(3) ^ 2;
t78 = t84 ^ 2;
t76 = t81 ^ 2;
t57 = -t65 * pkin(4) - t68;
t50 = t82 * t56 + t79 * t57;
t49 = -t79 * t56 + t82 * t57;
t1 = [MDP(1) + t76 * MDP(4) + 0.2e1 * t81 * t84 * MDP(5) + (t76 * t86 + (t86 + t87) * t78) * MDP(14) + t98 * t106 + (t77 * MDP(22) + MDP(15) - 0.2e1 * t96) * t65 ^ 2 + (MDP(20) * t106 + t99) * t64 + 0.2e1 * (t49 * t64 + t65 * t52) * MDP(27) + 0.2e1 * (t65 * t104 - t50 * t64) * MDP(28) + 0.2e1 * (t76 + t78) * MDP(13) * qJ(3) + (t78 * MDP(11) - t84 * t97) * t108 + 0.2e1 * (-MDP(16) + t94) * t65 * t64; (t91 * t79 - t104) * MDP(27) + (t91 * t82 + t52) * MDP(28) + (-qJ(3) * MDP(12) + MDP(7)) * t84 + (-qJ(3) * MDP(11) + MDP(6) + (-MDP(14) * qJ(3) - MDP(13)) * pkin(1)) * t81 + t89; MDP(11) * t108 + t87 * MDP(14) + MDP(8) + 0.2e1 * t107 + t95; t97 + t98 + (-MDP(14) * pkin(1) - MDP(11)) * t84 + t90 * t64; 0; MDP(14); (-t79 * t105 - t104) * MDP(27) + (-t82 * t105 + t52) * MDP(28) + t89; t107 + t95; 0; t95; t49 * MDP(27) - t50 * MDP(28) + t94 * t65 + t99; t92 * t70 + t101; t93; t92 * pkin(4) + t101; MDP(26);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
