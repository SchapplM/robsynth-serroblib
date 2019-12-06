% Calculate joint inertia matrix for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRRRP1_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:40:13
% EndTime: 2019-12-05 16:40:14
% DurationCPUTime: 0.11s
% Computational Cost: add. (101->52), mult. (177->72), div. (0->0), fcn. (118->4), ass. (0->28)
t56 = 2 * MDP(15);
t45 = cos(qJ(3));
t55 = t45 * pkin(2);
t35 = -pkin(3) - t55;
t54 = pkin(3) - t35;
t42 = sin(qJ(4));
t44 = cos(qJ(4));
t53 = t42 * MDP(10) + t44 * MDP(11);
t52 = t42 * MDP(14);
t51 = t42 * MDP(15);
t50 = t44 * MDP(13);
t41 = t42 ^ 2;
t49 = 0.2e1 * t42 * t44 * MDP(9) + t41 * MDP(8) + MDP(5);
t36 = -t44 * pkin(4) - pkin(3);
t48 = t50 - t52;
t47 = -MDP(13) * t42 - MDP(14) * t44;
t43 = sin(qJ(3));
t46 = (t45 * MDP(6) - t43 * MDP(7)) * pkin(2);
t40 = t44 * qJ(5);
t34 = t43 * pkin(2) + pkin(7);
t32 = t44 * pkin(7) + t40;
t31 = (-qJ(5) - pkin(7)) * t42;
t30 = t36 - t55;
t29 = t32 * t44;
t28 = t44 * t34 + t40;
t27 = (-qJ(5) - t34) * t42;
t26 = t28 * t44;
t1 = [MDP(1) + (t44 ^ 2 + t41) * MDP(16); (t44 * t27 + t42 * t28) * MDP(16); MDP(2) + (-t27 * t42 + t26) * t56 + (t27 ^ 2 + t28 ^ 2 + t30 ^ 2) * MDP(16) + t49 - 0.2e1 * t48 * t35 + 0.2e1 * t46; (t44 * t31 + t42 * t32) * MDP(16); (t26 + t29) * MDP(15) + (t27 * t31 + t28 * t32 + t30 * t36) * MDP(16) + t54 * t50 + t46 + (-t54 * MDP(14) + (-t27 - t31) * MDP(15)) * t42 + t49; (-t31 * t42 + t29) * t56 + (t31 ^ 2 + t32 ^ 2 + t36 ^ 2) * MDP(16) + 0.2e1 * t48 * pkin(3) + t49; -t52 + (MDP(16) * pkin(4) + MDP(13)) * t44; t47 * t34 + (t27 * MDP(16) - t51) * pkin(4) + t53; t47 * pkin(7) + (MDP(16) * t31 - t51) * pkin(4) + t53; pkin(4) ^ 2 * MDP(16) + MDP(12); 0; t30 * MDP(16); t36 * MDP(16); 0; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
