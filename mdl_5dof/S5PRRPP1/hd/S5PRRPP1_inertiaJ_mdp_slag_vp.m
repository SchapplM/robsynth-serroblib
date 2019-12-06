% Calculate joint inertia matrix for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRPP1_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:53
% EndTime: 2019-12-05 16:06:55
% DurationCPUTime: 0.19s
% Computational Cost: add. (185->58), mult. (340->90), div. (0->0), fcn. (326->4), ass. (0->23)
t69 = MDP(13) + MDP(17);
t68 = -qJ(4) - pkin(6);
t52 = sin(pkin(8));
t53 = cos(pkin(8));
t54 = sin(qJ(3));
t55 = cos(qJ(3));
t42 = t52 * t54 - t53 * t55;
t44 = t52 * t55 + t53 * t54;
t51 = -t55 * pkin(3) - pkin(2);
t35 = t42 * pkin(4) - t44 * qJ(5) + t51;
t67 = t35 * MDP(17);
t66 = t42 * MDP(14);
t65 = t44 * MDP(16);
t64 = t55 * MDP(10);
t46 = t68 * t55;
t61 = t68 * t54;
t38 = -t52 * t46 - t53 * t61;
t40 = -t53 * t46 + t52 * t61;
t63 = t38 ^ 2 + t40 ^ 2;
t57 = t65 - t66;
t49 = t53 * pkin(3) + pkin(4);
t47 = t52 * pkin(3) + qJ(5);
t1 = [MDP(1) + t69 * (t42 ^ 2 + t44 ^ 2); t69 * (t42 * t38 + t44 * t40); MDP(2) + 0.2e1 * pkin(2) * t64 + (t51 ^ 2 + t63) * MDP(13) + t63 * MDP(17) + (-0.2e1 * t65 + 0.2e1 * t66 + t67) * t35 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t54 + 0.2e1 * t55 * MDP(6)) * t54 + 0.2e1 * (MDP(12) + MDP(15)) * (t38 * t44 - t40 * t42); t64 - t54 * MDP(11) + (-t42 * t49 + t44 * t47) * MDP(17) + (-t42 * t53 + t44 * t52) * MDP(13) * pkin(3) + t57; t54 * MDP(7) + t55 * MDP(8) - t38 * MDP(14) + (-t47 * t42 - t49 * t44) * MDP(15) + t40 * MDP(16) + (-t38 * t49 + t40 * t47) * MDP(17) + (-t54 * MDP(10) - t55 * MDP(11)) * pkin(6) + ((-t42 * t52 - t44 * t53) * MDP(12) + (-t38 * t53 + t40 * t52) * MDP(13)) * pkin(3); MDP(9) + (t47 ^ 2 + t49 ^ 2) * MDP(17) + (t52 ^ 2 + t53 ^ 2) * MDP(13) * pkin(3) ^ 2 + 0.2e1 * t49 * MDP(14) + 0.2e1 * t47 * MDP(16); 0; t51 * MDP(13) - t57 + t67; 0; t69; t42 * MDP(17); t44 * MDP(15) + t38 * MDP(17); -t49 * MDP(17) - MDP(14); 0; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
