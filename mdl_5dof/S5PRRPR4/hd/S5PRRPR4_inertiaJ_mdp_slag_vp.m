% Calculate joint inertia matrix for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRRPR4_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:52:43
% EndTime: 2021-01-15 15:52:44
% DurationCPUTime: 0.20s
% Computational Cost: add. (234->77), mult. (481->119), div. (0->0), fcn. (509->8), ass. (0->38)
t89 = MDP(15) * pkin(3);
t66 = sin(pkin(9));
t67 = cos(pkin(9));
t69 = sin(qJ(3));
t72 = cos(qJ(3));
t58 = t66 * t69 - t67 * t72;
t65 = -t72 * pkin(3) - pkin(2);
t88 = 0.2e1 * t58 * pkin(4) + 0.2e1 * t65;
t87 = pkin(3) * t66;
t86 = qJ(4) + pkin(6);
t59 = t66 * t72 + t67 * t69;
t70 = sin(qJ(2));
t53 = t59 * t70;
t54 = t58 * t70;
t68 = sin(qJ(5));
t71 = cos(qJ(5));
t85 = (-t71 * t53 + t68 * t54) * MDP(21) + (t68 * t53 + t71 * t54) * MDP(22);
t48 = t71 * t58 + t68 * t59;
t84 = t48 * MDP(21);
t64 = t67 * pkin(3) + pkin(4);
t83 = (t71 * t64 - t68 * t87) * MDP(21);
t82 = (-t68 * t64 - t71 * t87) * MDP(22);
t81 = t58 * MDP(12);
t80 = t59 * MDP(13);
t79 = t65 * MDP(15);
t78 = t72 * MDP(10);
t61 = t86 * t69;
t62 = t86 * t72;
t50 = -t67 * t61 - t66 * t62;
t44 = -t59 * pkin(7) + t50;
t51 = -t66 * t61 + t67 * t62;
t45 = -t58 * pkin(7) + t51;
t49 = -t68 * t58 + t71 * t59;
t77 = t49 * MDP(18) - t48 * MDP(19) + (t71 * t44 - t68 * t45) * MDP(21) + (-t68 * t44 - t71 * t45) * MDP(22);
t76 = -t69 * MDP(10) - t72 * MDP(11);
t75 = t49 * MDP(22) + t79 + t80 + t81 + t84;
t73 = cos(qJ(2));
t1 = [MDP(1) + (t53 ^ 2 + t54 ^ 2 + t73 ^ 2) * MDP(15); -t70 * MDP(4) + (t53 * t59 + t54 * t58) * MDP(14) + (-t53 * t50 - t54 * t51) * MDP(15) + (-t69 * MDP(11) + MDP(3) - t75 + t78) * t73; MDP(2) + 0.2e1 * pkin(2) * t78 + 0.2e1 * (-t50 * t59 - t51 * t58) * MDP(14) + (t50 ^ 2 + t51 ^ 2) * MDP(15) + t84 * t88 + (t79 + 0.2e1 * t80 + 0.2e1 * t81) * t65 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t69 + 0.2e1 * t72 * MDP(6)) * t69 + (MDP(16) * t49 - 0.2e1 * t48 * MDP(17) + MDP(22) * t88) * t49; -t53 * MDP(12) + t54 * MDP(13) + t76 * t70 + (-t53 * t67 - t54 * t66) * t89 + t85; t50 * MDP(12) - t51 * MDP(13) + t69 * MDP(7) + t72 * MDP(8) + t76 * pkin(6) + ((-t58 * t66 - t59 * t67) * MDP(14) + (t50 * t67 + t51 * t66) * MDP(15)) * pkin(3) + t77; MDP(9) + MDP(20) + 0.2e1 * t83 + 0.2e1 * t82 + (0.2e1 * t67 * MDP(12) - 0.2e1 * t66 * MDP(13) + (t66 ^ 2 + t67 ^ 2) * t89) * pkin(3); -t73 * MDP(15); t75; 0; MDP(15); t85; t77; MDP(20) + t82 + t83; 0; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
