% Calculate joint inertia matrix for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRRP3_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:14
% EndTime: 2019-12-05 16:44:14
% DurationCPUTime: 0.13s
% Computational Cost: add. (149->53), mult. (276->85), div. (0->0), fcn. (266->4), ass. (0->24)
t58 = pkin(6) + pkin(7);
t45 = sin(qJ(4));
t57 = pkin(3) * t45;
t56 = cos(qJ(3));
t55 = cos(qJ(4));
t46 = sin(qJ(3));
t36 = t45 * t46 - t55 * t56;
t38 = t45 * t56 + t55 * t46;
t54 = -t36 * MDP(17) - t38 * MDP(18);
t53 = MDP(20) * pkin(4);
t40 = t58 * t46;
t41 = t58 * t56;
t52 = -t55 * t40 - t45 * t41;
t51 = MDP(17) * t55;
t50 = t56 * MDP(10);
t47 = t45 * t40 - t55 * t41;
t49 = t38 * MDP(14) - t36 * MDP(15) + t52 * MDP(17) + t47 * MDP(18);
t48 = -t56 * pkin(3) - pkin(2);
t43 = t55 * pkin(3) + pkin(4);
t35 = t38 ^ 2;
t30 = t36 * pkin(4) + t48;
t27 = -t36 * qJ(5) - t47;
t26 = -t38 * qJ(5) + t52;
t1 = [MDP(1) + (t36 ^ 2 + t35) * MDP(20); (-t36 * t26 + t38 * t27) * MDP(20); MDP(2) + 0.2e1 * pkin(2) * t50 + t35 * MDP(12) - 0.2e1 * t38 * t36 * MDP(13) + 0.2e1 * (-t26 * t38 - t27 * t36) * MDP(19) + (t26 ^ 2 + t27 ^ 2 + t30 ^ 2) * MDP(20) - 0.2e1 * t54 * t48 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t46 + 0.2e1 * t56 * MDP(6)) * t46; t50 - t46 * MDP(11) + (-t36 * t43 + t38 * t57) * MDP(20) + t54; t46 * MDP(7) + t56 * MDP(8) + (-t36 * t57 - t43 * t38) * MDP(19) + (t26 * t43 + t27 * t57) * MDP(20) + (-t46 * MDP(10) - t56 * MDP(11)) * pkin(6) + t49; t43 ^ 2 * MDP(20) + MDP(16) + MDP(9) + (0.2e1 * t51 + (MDP(20) * t57 - 0.2e1 * MDP(18)) * t45) * pkin(3); -t36 * t53 + t54; (-t38 * MDP(19) + t26 * MDP(20)) * pkin(4) + t49; t43 * t53 + MDP(16) + (-MDP(18) * t45 + t51) * pkin(3); pkin(4) ^ 2 * MDP(20) + MDP(16); 0; t30 * MDP(20); 0; 0; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
