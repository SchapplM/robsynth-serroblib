% Calculate joint inertia matrix for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRPRP4_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:36:07
% EndTime: 2019-12-05 15:36:07
% DurationCPUTime: 0.12s
% Computational Cost: add. (108->48), mult. (206->74), div. (0->0), fcn. (174->6), ass. (0->24)
t43 = sin(pkin(8));
t44 = cos(pkin(8));
t46 = sin(qJ(2));
t48 = cos(qJ(2));
t34 = t43 * t46 - t44 * t48;
t59 = t34 ^ 2;
t58 = MDP(5) * pkin(2);
t45 = sin(qJ(4));
t41 = t45 ^ 2;
t47 = cos(qJ(4));
t57 = t47 ^ 2 + t41;
t56 = MDP(12) - MDP(15);
t40 = -t44 * pkin(2) - pkin(3);
t55 = t57 * MDP(14);
t54 = t57 * MDP(16);
t53 = -pkin(4) * MDP(16) - MDP(13);
t52 = MDP(11) - t53;
t51 = MDP(16) * qJ(5) - t56;
t50 = -t52 * t45 + t51 * t47;
t39 = t43 * pkin(2) + pkin(6);
t36 = t43 * t48 + t44 * t46;
t33 = t36 ^ 2;
t32 = -t47 * pkin(4) - t45 * qJ(5) + t40;
t1 = [MDP(1) + (t33 + t59) * MDP(5) + (t57 * t33 + t59) * MDP(16); t48 * MDP(3) - t46 * MDP(4) + (t39 * t54 + t43 * t58 + t55) * t36 + (-t44 * t58 + t32 * MDP(16) + (-MDP(11) - MDP(13)) * t47 + t56 * t45) * t34; MDP(2) + t41 * MDP(6) + (t57 * t39 ^ 2 + t32 ^ 2) * MDP(16) + (t43 ^ 2 + t44 ^ 2) * MDP(5) * pkin(2) ^ 2 + 0.2e1 * t39 * t55 + 0.2e1 * (-t40 * MDP(11) - t32 * MDP(13)) * t47 + 0.2e1 * (t40 * MDP(12) - t32 * MDP(15) + t47 * MDP(7)) * t45; 0; 0; MDP(5) + t54; t50 * t36; t45 * MDP(8) + t47 * MDP(9) + (-t45 * pkin(4) + t47 * qJ(5)) * MDP(14) + t50 * t39; t51 * t45 + t52 * t47; MDP(10) + 0.2e1 * pkin(4) * MDP(13) + 0.2e1 * qJ(5) * MDP(15) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(16); t45 * t36 * MDP(16); (t39 * MDP(16) + MDP(14)) * t45; -t47 * MDP(16); t53; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
