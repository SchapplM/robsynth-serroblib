% Calculate joint inertia matrix for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5RPPRP1_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:36:25
% EndTime: 2019-12-05 17:36:26
% DurationCPUTime: 0.15s
% Computational Cost: add. (147->59), mult. (270->88), div. (0->0), fcn. (209->6), ass. (0->29)
t51 = sin(pkin(7));
t44 = t51 * pkin(1) + qJ(3);
t54 = sin(qJ(4));
t67 = t44 * t54;
t55 = cos(qJ(4));
t49 = t55 ^ 2;
t66 = t54 ^ 2 + t49;
t65 = MDP(17) * pkin(4);
t50 = sin(pkin(8));
t64 = qJ(5) * t50;
t63 = MDP(12) * t54;
t62 = t54 * MDP(15);
t61 = t55 * MDP(15);
t52 = cos(pkin(8));
t60 = t55 * t52 * t44;
t53 = cos(pkin(7));
t45 = -t53 * pkin(1) - pkin(2);
t59 = MDP(14) + t65;
t42 = -t52 * pkin(3) - t50 * pkin(6) + t45;
t40 = t55 * t42;
t36 = -t55 * t64 + t40 + (-pkin(4) - t67) * t52;
t37 = t60 + (t42 - t64) * t54;
t58 = t36 * t55 + t37 * t54;
t57 = -(-t52 * t67 + t40) * MDP(14) + (t54 * t42 + t60) * MDP(15);
t47 = t52 ^ 2;
t46 = t50 ^ 2;
t43 = t44 ^ 2;
t41 = (pkin(4) * t54 + t44) * t50;
t1 = [MDP(1) + (t47 * t43 + t45 ^ 2) * MDP(8) + t47 * MDP(13) + (t36 ^ 2 + t37 ^ 2 + t41 ^ 2) * MDP(17) + (t51 ^ 2 + t53 ^ 2) * MDP(4) * pkin(1) ^ 2 + (-0.2e1 * t55 * t54 * MDP(10) + t43 * MDP(8) + t49 * MDP(9)) * t46 + 0.2e1 * (t47 * MDP(7) + (t54 * MDP(14) + MDP(7) + t61) * t46) * t44 + 0.2e1 * (-t45 * MDP(5) + t57) * t52 + 0.2e1 * (t45 * MDP(6) + (-MDP(11) * t55 + t63) * t52 - t58 * MDP(16)) * t50; (-t41 * t52 + (-t36 * t54 + t37 * t55) * t50) * MDP(17); MDP(4) + (t46 + t47) * MDP(8) + (t66 * t46 + t47) * MDP(17); t45 * MDP(8) + t58 * MDP(17) + (-t55 * MDP(14) - MDP(5) + t62) * t52 + (-t66 * MDP(16) + MDP(6)) * t50; 0; t66 * MDP(17) + MDP(8); t36 * t65 - t52 * MDP(13) + (-t63 + (-MDP(16) * pkin(4) + MDP(11)) * t55) * t50 - t57; (-t59 * t54 - t61) * t50; t59 * t55 - t62; pkin(4) ^ 2 * MDP(17) + MDP(13); t41 * MDP(17); -t52 * MDP(17); 0; 0; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
