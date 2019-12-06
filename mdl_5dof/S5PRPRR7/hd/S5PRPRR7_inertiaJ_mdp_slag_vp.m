% Calculate joint inertia matrix for
% S5PRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR7_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR7_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRPRR7_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:00:44
% EndTime: 2019-12-05 16:00:45
% DurationCPUTime: 0.12s
% Computational Cost: add. (91->47), mult. (159->59), div. (0->0), fcn. (144->6), ass. (0->25)
t36 = sin(qJ(5));
t37 = sin(qJ(4));
t39 = cos(qJ(5));
t40 = cos(qJ(4));
t31 = t36 * t40 + t37 * t39;
t45 = t36 * t37 - t39 * t40;
t57 = MDP(20) * t45 + t31 * MDP(21);
t56 = 0.2e1 * pkin(4) * t37 + (2 * qJ(3));
t42 = -pkin(2) - pkin(6);
t55 = -pkin(7) + t42;
t54 = (pkin(2) * MDP(7));
t41 = cos(qJ(2));
t53 = t57 * t41;
t51 = qJ(3) * MDP(7);
t50 = t31 * MDP(20);
t49 = t37 * MDP(13);
t48 = t40 * MDP(14);
t47 = MDP(5) - t54;
t33 = t55 * t37;
t34 = t55 * t40;
t46 = -t45 * MDP(17) - t31 * MDP(18) + (-t33 * t36 + t34 * t39) * MDP(20) + (-t33 * t39 - t34 * t36) * MDP(21);
t44 = -MDP(13) * t40 + MDP(14) * t37;
t43 = (MDP(20) * t39 - MDP(21) * t36) * pkin(4);
t38 = sin(qJ(2));
t1 = [MDP(1) + (t38 ^ 2 + t41 ^ 2) * MDP(7); (MDP(3) - t47) * t41 + (-MDP(21) * t45 - MDP(4) + MDP(6) + t48 + t49 + t50 + t51) * t38; t50 * t56 + MDP(2) + (MDP(8) * t40 - 0.2e1 * MDP(9) * t37) * t40 + ((-2 * MDP(5) + t54) * pkin(2)) - (-MDP(15) * t45 - 0.2e1 * MDP(16) * t31 + MDP(21) * t56) * t45 + (0.2e1 * MDP(6) + 0.2e1 * t48 + 0.2e1 * t49 + t51) * qJ(3); -t41 * MDP(7); t47; MDP(7); t41 * t44 + t53; (MDP(13) * t42 + MDP(10)) * t40 + (-MDP(14) * t42 - MDP(11)) * t37 + t46; -t44 - t57; MDP(12) + MDP(19) + 0.2e1 * t43; t53; t46; -t57; MDP(19) + t43; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
