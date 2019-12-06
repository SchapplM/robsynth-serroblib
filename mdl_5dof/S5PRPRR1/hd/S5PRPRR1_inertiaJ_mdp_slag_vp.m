% Calculate joint inertia matrix for
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% MDP [22x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(22,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [22 1]), ...
  'S5PRPRR1_inertiaJ_mdp_slag_vp: MDP has to be [22x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:05
% EndTime: 2019-12-05 15:43:06
% DurationCPUTime: 0.12s
% Computational Cost: add. (156->45), mult. (295->66), div. (0->0), fcn. (320->6), ass. (0->32)
t54 = sin(pkin(9));
t55 = cos(pkin(9));
t57 = sin(qJ(4));
t73 = cos(qJ(4));
t47 = t73 * t54 + t57 * t55;
t46 = t57 * t54 - t73 * t55;
t66 = t46 * MDP(14);
t56 = sin(qJ(5));
t58 = cos(qJ(5));
t41 = t58 * t46 + t56 * t47;
t37 = t41 * MDP(21);
t42 = -t56 * t46 + t58 * t47;
t70 = -t42 * MDP(22) - t37;
t76 = -t47 * MDP(15) - t66 + t70;
t51 = -t55 * pkin(3) - pkin(2);
t75 = 0.2e1 * t46 * pkin(4) + 0.2e1 * t51;
t74 = 0.2e1 * t51;
t72 = pkin(2) * MDP(8);
t71 = pkin(6) + qJ(3);
t69 = t54 ^ 2 + t55 ^ 2;
t68 = t54 * MDP(6);
t67 = t55 * MDP(5);
t48 = t71 * t54;
t49 = t71 * t55;
t65 = -t73 * t48 - t57 * t49;
t64 = t69 * MDP(8);
t35 = -t47 * pkin(7) + t65;
t61 = t57 * t48 - t73 * t49;
t36 = -t46 * pkin(7) - t61;
t63 = t42 * MDP(18) - t41 * MDP(19) + (t58 * t35 - t56 * t36) * MDP(21) + (-t56 * t35 - t58 * t36) * MDP(22);
t60 = (MDP(21) * t58 - MDP(22) * t56) * pkin(4);
t1 = [MDP(1) + t64; 0; t66 * t74 + t37 * t75 + MDP(2) + (0.2e1 * t67 - 0.2e1 * t68 + t72) * pkin(2) + (-0.2e1 * t46 * MDP(10) + MDP(15) * t74 + MDP(9) * t47) * t47 + (MDP(16) * t42 - 0.2e1 * t41 * MDP(17) + MDP(22) * t75) * t42 + (0.2e1 * t69 * MDP(7) + t64 * qJ(3)) * qJ(3); 0; -t67 + t68 - t72 - t76; MDP(8); t76; t47 * MDP(11) - t46 * MDP(12) + t65 * MDP(14) + t61 * MDP(15) + t63; 0; MDP(13) + MDP(20) + 0.2e1 * t60; t70; t63; 0; MDP(20) + t60; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
