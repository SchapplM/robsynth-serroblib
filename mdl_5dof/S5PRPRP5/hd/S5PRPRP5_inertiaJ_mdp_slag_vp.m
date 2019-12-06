% Calculate joint inertia matrix for
% S5PRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRP5_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:38:42
% EndTime: 2019-12-05 15:38:43
% DurationCPUTime: 0.16s
% Computational Cost: add. (173->61), mult. (348->90), div. (0->0), fcn. (335->6), ass. (0->30)
t47 = sin(pkin(8));
t48 = cos(pkin(8));
t49 = sin(qJ(4));
t68 = cos(qJ(4));
t39 = t68 * t47 + t49 * t48;
t42 = -t48 * pkin(3) - pkin(2);
t54 = -t49 * t47 + t68 * t48;
t32 = -pkin(4) * t54 - t39 * qJ(5) + t42;
t60 = MDP(15) - MDP(18);
t63 = t48 * MDP(5);
t64 = t47 * MDP(6);
t67 = pkin(2) * MDP(8);
t69 = t32 * MDP(19) - (MDP(14) + MDP(16)) * t54 + t60 * t39 - t63 + t64 - t67;
t66 = pkin(6) + qJ(3);
t65 = t47 ^ 2 + t48 ^ 2;
t62 = MDP(8) + MDP(19);
t59 = t66 * t47;
t58 = t65 * qJ(3);
t57 = -pkin(4) * MDP(19) - MDP(16);
t56 = -MDP(14) + t57;
t55 = MDP(19) * qJ(5) - t60;
t51 = cos(qJ(2));
t50 = sin(qJ(2));
t46 = t51 ^ 2;
t40 = t66 * t48;
t36 = t54 * t50;
t35 = t39 * t50;
t34 = t68 * t40 - t49 * t59;
t33 = t49 * t40 + t68 * t59;
t1 = [MDP(1) + (t65 * t50 ^ 2 + t46) * MDP(8) + (t35 ^ 2 + t36 ^ 2 + t46) * MDP(19); (t35 * t39 + t36 * t54) * MDP(17) + (t35 * t33 + t36 * t34) * MDP(19) + (t65 * MDP(7) + MDP(8) * t58 - MDP(4)) * t50 + (MDP(3) - t69) * t51; MDP(2) + (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) * MDP(19) + t65 * MDP(8) * qJ(3) ^ 2 + 0.2e1 * MDP(7) * t58 + (0.2e1 * t63 - 0.2e1 * t64 + t67) * pkin(2) - 0.2e1 * (t42 * MDP(14) + t32 * MDP(16) - MDP(17) * t34) * t54 + (0.2e1 * MDP(10) * t54 + 0.2e1 * t42 * MDP(15) + 0.2e1 * MDP(17) * t33 - 0.2e1 * t32 * MDP(18) + MDP(9) * t39) * t39; -t62 * t51; t69; t62; t56 * t35 + t55 * t36; t39 * MDP(11) + t54 * MDP(12) + (-pkin(4) * t39 + qJ(5) * t54) * MDP(17) + t55 * t34 + t56 * t33; 0; MDP(13) + 0.2e1 * pkin(4) * MDP(16) + 0.2e1 * qJ(5) * MDP(18) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(19); t35 * MDP(19); t39 * MDP(17) + t33 * MDP(19); 0; t57; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
