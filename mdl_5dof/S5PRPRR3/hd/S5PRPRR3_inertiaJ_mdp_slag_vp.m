% Calculate joint inertia matrix for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRR3_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:47:34
% EndTime: 2019-12-05 15:47:35
% DurationCPUTime: 0.12s
% Computational Cost: add. (107->42), mult. (207->69), div. (0->0), fcn. (213->8), ass. (0->28)
t54 = sin(qJ(4));
t57 = cos(qJ(4));
t64 = t57 * MDP(11);
t53 = sin(qJ(5));
t56 = cos(qJ(5));
t44 = t53 * t54 - t56 * t57;
t35 = t44 * MDP(18);
t45 = t53 * t57 + t54 * t56;
t65 = -t45 * MDP(19) - t35;
t70 = -MDP(12) * t54 + t64 + t65;
t52 = cos(pkin(9));
t49 = -pkin(2) * t52 - pkin(3);
t69 = -0.2e1 * pkin(4) * t57 + 0.2e1 * t49;
t51 = sin(pkin(9));
t48 = pkin(2) * t51 + pkin(6);
t68 = pkin(7) + t48;
t67 = MDP(5) * pkin(2);
t55 = sin(qJ(2));
t58 = cos(qJ(2));
t43 = t51 * t58 + t52 * t55;
t66 = (-MDP(18) * t45 + MDP(19) * t44) * t43;
t39 = t68 * t54;
t40 = t68 * t57;
t63 = t45 * MDP(15) - t44 * MDP(16) + (-t39 * t56 - t40 * t53) * MDP(18) + (t39 * t53 - t40 * t56) * MDP(19);
t61 = -MDP(11) * t54 - MDP(12) * t57;
t60 = (MDP(18) * t56 - MDP(19) * t53) * pkin(4);
t41 = t51 * t55 - t52 * t58;
t1 = [MDP(1) + (t41 ^ 2 + t43 ^ 2) * MDP(5); t43 * t51 * t67 + t58 * MDP(3) - t55 * MDP(4) + (-t52 * t67 - t70) * t41; -0.2e1 * t49 * t64 + t35 * t69 + MDP(2) + (t51 ^ 2 + t52 ^ 2) * MDP(5) * pkin(2) ^ 2 + (0.2e1 * MDP(12) * t49 + MDP(6) * t54 + 0.2e1 * MDP(7) * t57) * t54 + (MDP(13) * t45 - 0.2e1 * MDP(14) * t44 + MDP(19) * t69) * t45; 0; 0; MDP(5); t43 * t61 + t66; t54 * MDP(8) + t57 * MDP(9) + t48 * t61 + t63; t70; MDP(10) + MDP(17) + 0.2e1 * t60; t66; t63; t65; MDP(17) + t60; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
