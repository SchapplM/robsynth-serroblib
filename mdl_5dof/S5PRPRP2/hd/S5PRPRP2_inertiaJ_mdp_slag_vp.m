% Calculate joint inertia matrix for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:05
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRPRP2_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:04:52
% EndTime: 2021-01-15 15:04:53
% DurationCPUTime: 0.16s
% Computational Cost: add. (139->65), mult. (287->91), div. (0->0), fcn. (222->4), ass. (0->27)
t44 = sin(qJ(4));
t52 = MDP(14) + MDP(16);
t59 = t52 * t44;
t45 = cos(qJ(4));
t41 = t45 ^ 2;
t58 = t44 ^ 2 + t41;
t57 = MDP(18) * pkin(4);
t56 = qJ(3) * t44;
t42 = sin(pkin(8));
t55 = qJ(5) * t42;
t54 = t44 * MDP(11);
t53 = -MDP(13) - MDP(15);
t43 = cos(pkin(8));
t51 = t45 * t43 * qJ(3);
t37 = -t43 * pkin(3) - t42 * pkin(6) - pkin(2);
t35 = t45 * t37;
t31 = -t45 * t55 + t35 + (-pkin(4) - t56) * t43;
t32 = t51 + (t37 - t55) * t44;
t50 = t31 * t45 + t32 * t44;
t49 = -t53 + t57;
t48 = t44 * MDP(15) + t45 * MDP(16);
t47 = (-t43 * t56 + t35) * MDP(13) - (t44 * t37 + t51) * MDP(14) - t32 * MDP(16);
t46 = qJ(3) ^ 2;
t39 = t43 ^ 2;
t38 = t42 ^ 2;
t36 = (pkin(4) * t44 + qJ(3)) * t42;
t1 = [MDP(1) + (t38 + t39) * MDP(7) + (t58 * t38 + t39) * MDP(18); (-t43 * t36 + (-t31 * t44 + t32 * t45) * t42) * MDP(18); MDP(2) + (pkin(2) ^ 2 + t39 * t46) * MDP(7) + t39 * MDP(12) + (t31 ^ 2 + t32 ^ 2 + t36 ^ 2) * MDP(18) + (-0.2e1 * t45 * t44 * MDP(9) + t46 * MDP(7) + t41 * MDP(8)) * t38 + 0.2e1 * (-t50 * MDP(17) + t48 * t36) * t42 + 0.2e1 * (t39 * MDP(6) + (t44 * MDP(13) + t45 * MDP(14) + MDP(6)) * t38) * qJ(3) + 0.2e1 * (pkin(2) * MDP(5) + (-t45 * MDP(10) + t54) * t42 - t31 * MDP(15) - t47) * t43; 0; -pkin(2) * MDP(7) + t50 * MDP(18) - t58 * MDP(17) * t42 + (t53 * t45 - MDP(5) + t59) * t43; t58 * MDP(18) + MDP(7); (-t49 * t44 - t52 * t45) * t42; t31 * t57 + t35 * MDP(15) + (-MDP(12) + (-0.2e1 * pkin(4) - t56) * MDP(15)) * t43 + (-t54 + (-MDP(15) * qJ(5) - MDP(17) * pkin(4) + MDP(10)) * t45) * t42 + t47; t49 * t45 - t59; MDP(12) + (0.2e1 * MDP(15) + t57) * pkin(4); -t43 * MDP(18); t36 * MDP(18) + t48 * t42; 0; 0; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
