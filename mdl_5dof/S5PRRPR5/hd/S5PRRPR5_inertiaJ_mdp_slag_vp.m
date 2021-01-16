% Calculate joint inertia matrix for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR5_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:04:45
% EndTime: 2021-01-15 16:04:46
% DurationCPUTime: 0.23s
% Computational Cost: add. (278->92), mult. (601->154), div. (0->0), fcn. (648->10), ass. (0->48)
t80 = cos(qJ(3));
t71 = -t80 * pkin(3) - pkin(2);
t104 = 0.2e1 * t71;
t103 = sin(qJ(2));
t74 = sin(pkin(10));
t76 = cos(pkin(10));
t78 = sin(qJ(3));
t65 = t74 * t80 + t76 * t78;
t77 = sin(qJ(5));
t102 = t65 * t77;
t79 = cos(qJ(5));
t101 = t65 * t79;
t75 = sin(pkin(5));
t81 = cos(qJ(2));
t100 = t75 * t81;
t99 = t77 * t79;
t98 = qJ(4) + pkin(7);
t97 = MDP(15) * pkin(3);
t96 = cos(pkin(5));
t95 = MDP(10) * t80;
t64 = t74 * t78 - t76 * t80;
t94 = t64 * MDP(20);
t93 = t65 * MDP(13);
t92 = MDP(17) * t99;
t91 = t75 * t103;
t90 = t98 * t78;
t89 = t71 * MDP(15) + t93;
t88 = t79 * MDP(21) - t77 * MDP(22);
t87 = MDP(21) * t77 + MDP(22) * t79;
t86 = MDP(12) + t88;
t85 = (MDP(18) * t79 - MDP(19) * t77) * t65;
t84 = -t78 * t91 + t96 * t80;
t83 = t77 * MDP(18) + t79 * MDP(19) - t87 * (t74 * pkin(3) + pkin(8));
t73 = t79 ^ 2;
t72 = t77 ^ 2;
t70 = -t76 * pkin(3) - pkin(4);
t67 = t98 * t80;
t62 = t96 * t78 + t80 * t91;
t60 = t76 * t67 - t74 * t90;
t58 = t74 * t67 + t76 * t90;
t57 = t64 * pkin(4) - t65 * pkin(8) + t71;
t56 = t76 * t62 + t74 * t84;
t54 = t74 * t62 - t76 * t84;
t53 = -t77 * t100 + t79 * t56;
t52 = -t79 * t100 - t77 * t56;
t51 = t77 * t57 + t79 * t60;
t50 = t79 * t57 - t77 * t60;
t1 = [MDP(1) + (t75 ^ 2 * t81 ^ 2 + t54 ^ 2 + t56 ^ 2) * MDP(15); (t54 * t65 - t56 * t64) * MDP(14) + (t54 * t58 + t56 * t60) * MDP(15) + (t54 * t102 + t52 * t64) * MDP(21) + (t54 * t101 - t53 * t64) * MDP(22) + (-t103 * MDP(4) + (-MDP(11) * t78 - t64 * MDP(12) + MDP(3) - t89 + t95) * t81) * t75; MDP(2) + 0.2e1 * pkin(2) * t95 + t93 * t104 + (t58 ^ 2 + t60 ^ 2 + t71 ^ 2) * MDP(15) + (t73 * MDP(16) - 0.2e1 * t92) * t65 ^ 2 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t78 + 0.2e1 * t80 * MDP(6)) * t78 + (MDP(12) * t104 + 0.2e1 * t85 + t94) * t64 + 0.2e1 * (t58 * t65 - t60 * t64) * MDP(14) + 0.2e1 * (t58 * t102 + t50 * t64) * MDP(21) + 0.2e1 * (t58 * t101 - t51 * t64) * MDP(22); t84 * MDP(10) - t62 * MDP(11) + (t74 * t97 - MDP(13)) * t56 + (-t76 * t97 - t86) * t54; -t60 * MDP(13) + t78 * MDP(7) + t80 * MDP(8) - t86 * t58 + t83 * t64 + (-t78 * MDP(10) - t80 * MDP(11)) * pkin(7) + (MDP(16) * t99 + (-t72 + t73) * MDP(17) + t87 * t70) * t65 + ((-t64 * t74 - t65 * t76) * MDP(14) + (-t58 * t76 + t60 * t74) * MDP(15)) * pkin(3); 0.2e1 * t92 + t72 * MDP(16) + MDP(9) + (t74 ^ 2 + t76 ^ 2) * MDP(15) * pkin(3) ^ 2 - 0.2e1 * t88 * t70 + 0.2e1 * (t76 * MDP(12) - t74 * MDP(13)) * pkin(3); -MDP(15) * t100; t86 * t64 + t89; 0; MDP(15); t52 * MDP(21) - t53 * MDP(22); t50 * MDP(21) - t51 * MDP(22) + t85 + t94; t83; t88; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
