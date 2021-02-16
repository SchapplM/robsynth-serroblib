% Calculate joint inertia matrix for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S4RRPR7_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:56:39
% EndTime: 2021-01-15 10:56:40
% DurationCPUTime: 0.19s
% Computational Cost: add. (177->63), mult. (355->104), div. (0->0), fcn. (335->6), ass. (0->33)
t59 = cos(qJ(2));
t51 = -t59 * pkin(2) - pkin(1);
t74 = 0.2e1 * t51;
t73 = 0.2e1 * t59;
t70 = -qJ(3) - pkin(5);
t47 = t70 * t59;
t54 = sin(pkin(7));
t55 = cos(pkin(7));
t57 = sin(qJ(2));
t66 = t70 * t57;
t40 = -t54 * t47 - t55 * t66;
t45 = t54 * t59 + t55 * t57;
t72 = t40 * t45;
t56 = sin(qJ(4));
t58 = cos(qJ(4));
t71 = t56 * t58;
t44 = t54 * t57 - t55 * t59;
t69 = t44 * MDP(19);
t68 = t45 * MDP(12);
t67 = MDP(16) * t71;
t65 = t58 * MDP(20) - t56 * MDP(21);
t64 = MDP(20) * t56 + MDP(21) * t58;
t63 = MDP(11) + t65;
t62 = (MDP(17) * t58 - MDP(18) * t56) * t45;
t61 = t56 * MDP(17) + t58 * MDP(18) - (t54 * pkin(2) + pkin(6)) * t64;
t53 = t58 ^ 2;
t52 = t56 ^ 2;
t50 = -t55 * pkin(2) - pkin(3);
t42 = -t55 * t47 + t54 * t66;
t39 = t44 * pkin(3) - t45 * pkin(6) + t51;
t38 = t56 * t39 + t58 * t42;
t37 = t58 * t39 - t56 * t42;
t1 = [MDP(1) + pkin(1) * MDP(9) * t73 + t68 * t74 + (t40 ^ 2 + t42 ^ 2 + t51 ^ 2) * MDP(14) + (t53 * MDP(15) - 0.2e1 * t67) * t45 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t57 + MDP(5) * t73) * t57 + (MDP(11) * t74 + 0.2e1 * t62 + t69) * t44 + 0.2e1 * (-t42 * t44 + t72) * MDP(13) + 0.2e1 * (t37 * t44 + t56 * t72) * MDP(20) + 0.2e1 * (-t38 * t44 + t58 * t72) * MDP(21); -t42 * MDP(12) + t57 * MDP(6) + t59 * MDP(7) - t63 * t40 + t61 * t44 + (-t59 * MDP(10) - t57 * MDP(9)) * pkin(5) + (MDP(15) * t71 + (-t52 + t53) * MDP(16) + t64 * t50) * t45 + ((-t44 * t54 - t45 * t55) * MDP(13) + (-t40 * t55 + t42 * t54) * MDP(14)) * pkin(2); 0.2e1 * t67 + t52 * MDP(15) + MDP(8) + (t54 ^ 2 + t55 ^ 2) * MDP(14) * pkin(2) ^ 2 - 0.2e1 * t65 * t50 + 0.2e1 * (t55 * MDP(11) - t54 * MDP(12)) * pkin(2); t51 * MDP(14) + t44 * t63 + t68; 0; MDP(14); t37 * MDP(20) - t38 * MDP(21) + t62 + t69; t61; t65; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
