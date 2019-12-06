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
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRPRP2_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:31:05
% EndTime: 2019-12-05 15:31:05
% DurationCPUTime: 0.13s
% Computational Cost: add. (106->55), mult. (224->81), div. (0->0), fcn. (168->4), ass. (0->25)
t44 = cos(qJ(4));
t40 = t44 ^ 2;
t43 = sin(qJ(4));
t56 = t43 ^ 2 + t40;
t55 = MDP(17) * pkin(4);
t54 = qJ(3) * t43;
t41 = sin(pkin(8));
t53 = qJ(5) * t41;
t52 = MDP(12) * t43;
t51 = t43 * MDP(15);
t50 = t44 * MDP(15);
t42 = cos(pkin(8));
t49 = t44 * t42 * qJ(3);
t48 = MDP(14) + t55;
t36 = -t42 * pkin(3) - t41 * pkin(6) - pkin(2);
t34 = t44 * t36;
t30 = -t44 * t53 + t34 + (-pkin(4) - t54) * t42;
t31 = t49 + (t36 - t53) * t43;
t47 = t30 * t44 + t31 * t43;
t46 = -(-t42 * t54 + t34) * MDP(14) + (t43 * t36 + t49) * MDP(15);
t45 = qJ(3) ^ 2;
t38 = t42 ^ 2;
t37 = t41 ^ 2;
t35 = (pkin(4) * t43 + qJ(3)) * t41;
t1 = [MDP(1) + (t37 + t38) * MDP(8) + (t56 * t37 + t38) * MDP(17); (-t42 * t35 + (-t30 * t43 + t31 * t44) * t41) * MDP(17); MDP(2) + (pkin(2) ^ 2 + t38 * t45) * MDP(8) + t38 * MDP(13) + (t30 ^ 2 + t31 ^ 2 + t35 ^ 2) * MDP(17) + (-0.2e1 * t44 * t43 * MDP(10) + t45 * MDP(8) + t40 * MDP(9)) * t37 + 0.2e1 * (t38 * MDP(7) + (t43 * MDP(14) + MDP(7) + t50) * t37) * qJ(3) + 0.2e1 * (pkin(2) * MDP(5) + t46) * t42 + 0.2e1 * (-pkin(2) * MDP(6) + (-MDP(11) * t44 + t52) * t42 - t47 * MDP(16)) * t41; 0; -pkin(2) * MDP(8) + t47 * MDP(17) + (-t44 * MDP(14) - MDP(5) + t51) * t42 + (-t56 * MDP(16) + MDP(6)) * t41; t56 * MDP(17) + MDP(8); (-t48 * t43 - t50) * t41; t30 * t55 - t42 * MDP(13) + (-t52 + (-MDP(16) * pkin(4) + MDP(11)) * t44) * t41 - t46; t48 * t44 - t51; pkin(4) ^ 2 * MDP(17) + MDP(13); -t42 * MDP(17); t35 * MDP(17); 0; 0; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
