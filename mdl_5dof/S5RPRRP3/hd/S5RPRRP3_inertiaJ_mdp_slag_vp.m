% Calculate joint inertia matrix for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% MDP [20x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5RPRRP3_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:54
% EndTime: 2019-12-05 18:03:55
% DurationCPUTime: 0.15s
% Computational Cost: add. (185->57), mult. (317->92), div. (0->0), fcn. (302->6), ass. (0->28)
t52 = sin(pkin(8));
t48 = t52 * pkin(1) + pkin(6);
t68 = pkin(7) + t48;
t54 = sin(qJ(4));
t67 = pkin(3) * t54;
t66 = cos(qJ(3));
t65 = cos(qJ(4));
t55 = sin(qJ(3));
t44 = t54 * t55 - t65 * t66;
t46 = t54 * t66 + t65 * t55;
t64 = -t44 * MDP(17) - t46 * MDP(18);
t63 = MDP(20) * pkin(4);
t53 = cos(pkin(8));
t49 = -t53 * pkin(1) - pkin(2);
t41 = t68 * t55;
t42 = t68 * t66;
t62 = -t65 * t41 - t54 * t42;
t61 = MDP(17) * t65;
t60 = t66 * MDP(10);
t58 = t54 * t41 - t65 * t42;
t59 = t46 * MDP(14) - t44 * MDP(15) + t62 * MDP(17) + t58 * MDP(18);
t57 = -t66 * pkin(3) + t49;
t51 = t65 * pkin(3) + pkin(4);
t43 = t46 ^ 2;
t35 = t44 * pkin(4) + t57;
t32 = -t44 * qJ(5) - t58;
t31 = -t46 * qJ(5) + t62;
t1 = [MDP(1) - 0.2e1 * t49 * t60 + t43 * MDP(12) - 0.2e1 * t46 * t44 * MDP(13) + 0.2e1 * (-t31 * t46 - t32 * t44) * MDP(19) + (t31 ^ 2 + t32 ^ 2 + t35 ^ 2) * MDP(20) + (t52 ^ 2 + t53 ^ 2) * MDP(4) * pkin(1) ^ 2 - 0.2e1 * t64 * t57 + (0.2e1 * t49 * MDP(11) + MDP(5) * t55 + 0.2e1 * t66 * MDP(6)) * t55; (-t31 * t44 + t32 * t46) * MDP(20); MDP(4) + (t44 ^ 2 + t43) * MDP(20); t55 * MDP(7) + t66 * MDP(8) + (-t44 * t67 - t51 * t46) * MDP(19) + (t31 * t51 + t32 * t67) * MDP(20) + (-t55 * MDP(10) - t66 * MDP(11)) * t48 + t59; t60 - t55 * MDP(11) + (-t44 * t51 + t46 * t67) * MDP(20) + t64; t51 ^ 2 * MDP(20) + MDP(16) + MDP(9) + (0.2e1 * t61 + (MDP(20) * t67 - 0.2e1 * MDP(18)) * t54) * pkin(3); (-t46 * MDP(19) + t31 * MDP(20)) * pkin(4) + t59; -t44 * t63 + t64; t51 * t63 + MDP(16) + (-MDP(18) * t54 + t61) * pkin(3); pkin(4) ^ 2 * MDP(20) + MDP(16); t35 * MDP(20); 0; 0; 0; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
