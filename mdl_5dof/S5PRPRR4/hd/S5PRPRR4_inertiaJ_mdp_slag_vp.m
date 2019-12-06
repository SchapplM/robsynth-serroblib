% Calculate joint inertia matrix for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_inertiaJ_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PRPRR4_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:36
% EndTime: 2019-12-05 15:51:36
% DurationCPUTime: 0.16s
% Computational Cost: add. (122->65), mult. (295->108), div. (0->0), fcn. (288->10), ass. (0->38)
t60 = sin(qJ(5));
t63 = cos(qJ(5));
t68 = MDP(18) * t60 + MDP(19) * t63;
t80 = t60 * MDP(15) + t63 * MDP(16) - t68 * pkin(8);
t56 = sin(pkin(10));
t57 = sin(pkin(5));
t58 = cos(pkin(10));
t62 = sin(qJ(2));
t65 = cos(qJ(2));
t48 = (t56 * t65 + t58 * t62) * t57;
t59 = cos(pkin(5));
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t44 = t48 * t61 - t59 * t64;
t79 = t44 * t61;
t51 = t56 * pkin(2) + pkin(7);
t78 = t51 * t60;
t77 = t51 * t63;
t76 = t51 * t64;
t75 = t60 * t63;
t74 = t61 * MDP(12);
t73 = t64 * MDP(17);
t72 = MDP(14) * t75;
t52 = -t58 * pkin(2) - pkin(3);
t71 = MDP(15) * t63 - MDP(16) * t60;
t69 = t63 * MDP(18) - t60 * MDP(19);
t67 = MDP(11) + t69;
t55 = t63 ^ 2;
t54 = t61 ^ 2;
t53 = t60 ^ 2;
t49 = -t64 * pkin(4) - t61 * pkin(8) + t52;
t46 = (t56 * t62 - t58 * t65) * t57;
t45 = t48 * t64 + t59 * t61;
t43 = t60 * t49 + t63 * t76;
t42 = t63 * t49 - t60 * t76;
t41 = t45 * t63 + t46 * t60;
t40 = -t45 * t60 + t46 * t63;
t1 = [MDP(1) + (t46 ^ 2 + t48 ^ 2 + t59 ^ 2) * MDP(5); (-t40 * t64 + t60 * t79) * MDP(18) + (t41 * t64 + t63 * t79) * MDP(19) + (t65 * MDP(3) - t62 * MDP(4)) * t57 + (-t64 * MDP(11) + t74) * t46 + (-t46 * t58 + t48 * t56) * MDP(5) * pkin(2); MDP(2) + (t56 ^ 2 + t58 ^ 2) * MDP(5) * pkin(2) ^ 2 + (-0.2e1 * t52 * MDP(11) + t73) * t64 + (t55 * MDP(13) + MDP(6) - 0.2e1 * t72) * t54 + 0.2e1 * (-t42 * t64 + t54 * t78) * MDP(18) + 0.2e1 * (t43 * t64 + t54 * t77) * MDP(19) + 0.2e1 * (t52 * MDP(12) + (MDP(7) - t71) * t64) * t61; t59 * MDP(5); 0; MDP(5); -t45 * MDP(12) - t67 * t44; (-t51 * MDP(12) + MDP(9) - t80) * t64 + (MDP(8) - t51 * MDP(11) + MDP(13) * t75 + (-t53 + t55) * MDP(14) + (-pkin(4) * t60 - t77) * MDP(18) + (-pkin(4) * t63 + t78) * MDP(19)) * t61; t67 * t64 - t74; t53 * MDP(13) + 0.2e1 * pkin(4) * t69 + MDP(10) + 0.2e1 * t72; t40 * MDP(18) - t41 * MDP(19); t42 * MDP(18) - t43 * MDP(19) + t71 * t61 - t73; -t68 * t61; t80; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
