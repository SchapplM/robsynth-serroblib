% Calculate joint inertia matrix for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPPPR5_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:31
% EndTime: 2019-12-31 17:46:31
% DurationCPUTime: 0.14s
% Computational Cost: add. (136->49), mult. (197->70), div. (0->0), fcn. (161->6), ass. (0->30)
t53 = sin(pkin(7));
t55 = cos(pkin(7));
t58 = -pkin(1) - pkin(2);
t43 = t55 * qJ(2) + t53 * t58;
t39 = -qJ(4) + t43;
t52 = sin(pkin(8));
t54 = cos(pkin(8));
t65 = t52 ^ 2 + t54 ^ 2;
t60 = t65 * MDP(13);
t69 = t39 * t60;
t62 = t54 * MDP(10);
t63 = t52 * MDP(11);
t41 = t53 * qJ(2) - t55 * t58;
t40 = pkin(3) + t41;
t64 = t40 * MDP(13);
t56 = sin(qJ(5));
t57 = cos(qJ(5));
t37 = t56 * t52 - t57 * t54;
t34 = t37 * MDP(19);
t38 = -t57 * t52 - t56 * t54;
t66 = t38 * MDP(20) - t34;
t68 = t62 - t63 + t64 + t66;
t67 = pkin(6) - t39;
t61 = t65 * MDP(12);
t51 = t55 ^ 2;
t49 = t53 ^ 2;
t33 = t54 * pkin(4) + t40;
t32 = t67 * t54;
t31 = t67 * t52;
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + (t41 ^ 2 + t43 ^ 2) * MDP(9) - 0.2e1 * t33 * t34 + (0.2e1 * t62 - 0.2e1 * t63 + t64) * t40 + (MDP(14) * t38 + 0.2e1 * t37 * MDP(15) + 0.2e1 * t33 * MDP(20)) * t38 + 0.2e1 * t41 * MDP(7) + 0.2e1 * t43 * MDP(8) + (-0.2e1 * t61 + t69) * t39; -pkin(1) * MDP(6) - MDP(4) + (t43 * MDP(9) + MDP(8) - t61 + t69) * t53 + (-t41 * MDP(9) - MDP(7) - t68) * t55; MDP(6) + (t49 + t51) * MDP(9) + (t65 * t49 + t51) * MDP(13); 0; 0; MDP(9) + t60; t68; -t55 * MDP(13); 0; MDP(13); t38 * MDP(16) + t37 * MDP(17) + (t57 * t31 + t56 * t32) * MDP(19) + (-t56 * t31 + t57 * t32) * MDP(20); (t38 * MDP(19) + t37 * MDP(20)) * t53; t66; 0; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
