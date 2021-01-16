% Calculate joint inertia matrix for
% S5PRPRP3
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
%   see S5PRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRPRP3_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:14:02
% EndTime: 2021-01-15 15:14:03
% DurationCPUTime: 0.12s
% Computational Cost: add. (110->50), mult. (204->72), div. (0->0), fcn. (181->6), ass. (0->29)
t45 = sin(qJ(4));
t54 = MDP(12) + MDP(14);
t63 = t54 * t45;
t43 = sin(pkin(8));
t44 = cos(pkin(8));
t46 = sin(qJ(2));
t48 = cos(qJ(2));
t34 = t43 * t46 - t44 * t48;
t62 = t34 ^ 2;
t61 = MDP(5) * pkin(2);
t41 = t45 ^ 2;
t47 = cos(qJ(4));
t60 = t47 ^ 2 + t41;
t59 = MDP(16) * pkin(4);
t39 = t43 * pkin(2) + pkin(6);
t58 = qJ(5) + t39;
t53 = -t44 * pkin(2) - pkin(3);
t37 = -t47 * pkin(4) + t53;
t57 = t37 * MDP(16);
t56 = t45 * MDP(14);
t55 = t47 * MDP(13);
t52 = MDP(13) + t59;
t31 = t58 * t45;
t32 = t58 * t47;
t51 = t31 * t45 + t32 * t47;
t50 = MDP(11) + t52;
t36 = t43 * t48 + t44 * t46;
t33 = t36 ^ 2;
t1 = [MDP(1) + (t33 + t62) * MDP(5) + (t60 * t33 + t62) * MDP(16); t48 * MDP(3) - t46 * MDP(4) + (t60 * MDP(15) + t51 * MDP(16) + t43 * t61) * t36 + (-t44 * t61 + t57 + (-MDP(11) - MDP(13)) * t47 + t63) * t34; MDP(2) + t41 * MDP(6) + 0.2e1 * t45 * t47 * MDP(7) + 0.2e1 * t51 * MDP(15) + (t31 ^ 2 + t32 ^ 2) * MDP(16) + (t43 ^ 2 + t44 ^ 2) * MDP(5) * pkin(2) ^ 2 + 0.2e1 * (-t47 * MDP(11) + t45 * MDP(12)) * t53 + (-0.2e1 * t55 + 0.2e1 * t56 + t57) * t37; 0; (-t31 * t47 + t32 * t45) * MDP(16); t60 * MDP(16) + MDP(5); (-t50 * t45 - t54 * t47) * t36; -t32 * MDP(14) + (-MDP(12) * t39 + MDP(9)) * t47 - t52 * t31 + (-MDP(11) * t39 - MDP(15) * pkin(4) + MDP(8)) * t45; t50 * t47 - t63; MDP(10) + (0.2e1 * MDP(13) + t59) * pkin(4); t34 * MDP(16); -t55 + t56 + t57; 0; 0; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
