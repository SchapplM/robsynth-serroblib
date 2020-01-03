% Calculate joint inertia matrix for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5RPRRP2_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:29
% EndTime: 2020-01-03 11:45:30
% DurationCPUTime: 0.12s
% Computational Cost: add. (151->58), mult. (258->81), div. (0->0), fcn. (192->6), ass. (0->35)
t71 = 2 * MDP(15);
t52 = sin(pkin(8));
t70 = pkin(1) * t52;
t56 = cos(qJ(4));
t69 = t56 * pkin(4);
t53 = cos(pkin(8));
t45 = t53 * pkin(1) + pkin(2);
t55 = sin(qJ(3));
t57 = cos(qJ(3));
t38 = t57 * t45 - t55 * t70;
t36 = -pkin(3) - t38;
t68 = pkin(3) - t36;
t54 = sin(qJ(4));
t67 = t54 * MDP(10) + t56 * MDP(11);
t66 = t38 * MDP(6);
t39 = -t55 * t45 - t57 * t70;
t65 = t39 * MDP(7);
t64 = t54 * MDP(14);
t63 = t54 * MDP(15);
t62 = t56 * MDP(13);
t51 = t54 ^ 2;
t61 = 0.2e1 * t54 * t56 * MDP(9) + t51 * MDP(8) + MDP(5);
t60 = t62 - t64;
t59 = -MDP(13) * t54 - MDP(14) * t56;
t50 = t56 * qJ(5);
t46 = -pkin(3) - t69;
t42 = t56 * pkin(7) + t50;
t41 = (-qJ(5) - pkin(7)) * t54;
t40 = t42 * t56;
t37 = pkin(7) - t39;
t35 = t36 - t69;
t34 = t56 * t37 + t50;
t33 = (-qJ(5) - t37) * t54;
t32 = t34 * t56;
t1 = [MDP(1) + (t33 ^ 2 + t34 ^ 2 + t35 ^ 2) * MDP(16) + (t52 ^ 2 + t53 ^ 2) * MDP(4) * pkin(1) ^ 2 - 0.2e1 * t60 * t36 + 0.2e1 * t66 + 0.2e1 * t65 + (-t33 * t54 + t32) * t71 + t61; (t33 * t56 + t34 * t54) * MDP(16); MDP(4) + (t56 ^ 2 + t51) * MDP(16); t66 + t65 + (t32 + t40) * MDP(15) + (t33 * t41 + t34 * t42 + t35 * t46) * MDP(16) + t68 * t62 + (-t68 * MDP(14) + (-t33 - t41) * MDP(15)) * t54 + t61; (t41 * t56 + t42 * t54) * MDP(16); (-t41 * t54 + t40) * t71 + (t41 ^ 2 + t42 ^ 2 + t46 ^ 2) * MDP(16) + 0.2e1 * t60 * pkin(3) + t61; t59 * t37 + (MDP(16) * t33 - t63) * pkin(4) + t67; -t64 + (MDP(16) * pkin(4) + MDP(13)) * t56; t59 * pkin(7) + (MDP(16) * t41 - t63) * pkin(4) + t67; MDP(16) * pkin(4) ^ 2 + MDP(12); t35 * MDP(16); 0; t46 * MDP(16); 0; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
