% Calculate joint inertia matrix for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR8_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5PRPRR8_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:27
% EndTime: 2019-12-05 16:04:29
% DurationCPUTime: 0.17s
% Computational Cost: add. (106->69), mult. (235->101), div. (0->0), fcn. (201->8), ass. (0->37)
t52 = sin(qJ(5));
t55 = cos(qJ(5));
t78 = t52 * MDP(20) + t55 * MDP(21);
t77 = (pkin(2) * MDP(7));
t51 = cos(pkin(5));
t53 = sin(qJ(4));
t56 = cos(qJ(4));
t50 = sin(pkin(5));
t57 = cos(qJ(2));
t74 = t50 * t57;
t42 = t51 * t53 + t56 * t74;
t76 = t42 * t56;
t54 = sin(qJ(2));
t75 = t50 * t54;
t73 = t52 * t55;
t58 = -pkin(2) - pkin(7);
t72 = t52 * t58;
t71 = t55 * t58;
t69 = qJ(3) * MDP(7);
t67 = t53 * MDP(13);
t65 = MDP(16) * t73;
t64 = MDP(5) - t77;
t63 = MDP(17) * t55 - MDP(18) * t52;
t62 = t55 * MDP(20) - t52 * MDP(21);
t60 = MDP(13) + t62;
t59 = t52 * MDP(17) + t55 * MDP(18) - pkin(8) * t78;
t49 = t56 ^ 2;
t48 = t55 ^ 2;
t47 = t53 ^ 2;
t46 = t52 ^ 2;
t44 = t53 * pkin(4) - t56 * pkin(8) + qJ(3);
t43 = t51 * t56 - t53 * t74;
t41 = t52 * t44 + t53 * t71;
t40 = t55 * t44 - t53 * t72;
t39 = t43 * t55 + t52 * t75;
t38 = -t43 * t52 + t55 * t75;
t1 = [MDP(1) + (t51 ^ 2 + (t54 ^ 2 + t57 ^ 2) * t50 ^ 2) * MDP(7); (t38 * t53 + t52 * t76) * MDP(20) + (-t39 * t53 + t55 * t76) * MDP(21) + ((MDP(3) - t64) * t57 + (t56 * MDP(14) - MDP(4) + MDP(6) + t67 + t69) * t54) * t50; t47 * MDP(19) + MDP(2) + ((-2 * MDP(5) + t77) * pkin(2)) + (0.2e1 * MDP(6) + 0.2e1 * t67 + t69) * qJ(3) + (t48 * MDP(15) + MDP(8) - 0.2e1 * t65) * t49 + 0.2e1 * (t40 * t53 - t49 * t72) * MDP(20) + 0.2e1 * (-t41 * t53 - t49 * t71) * MDP(21) + 0.2e1 * (qJ(3) * MDP(14) + (-MDP(9) + t63) * t53) * t56; -MDP(7) * t74; t64 + t78 * (-t47 - t49); MDP(7); -t43 * MDP(14) - t60 * t42; (-t58 * MDP(14) - MDP(11) + t59) * t53 + (MDP(10) + t58 * MDP(13) + MDP(15) * t73 + (-t46 + t48) * MDP(16) + (-pkin(4) * t52 + t71) * MDP(20) + (-pkin(4) * t55 - t72) * MDP(21)) * t56; -t53 * MDP(14) + t60 * t56; t46 * MDP(15) + 0.2e1 * pkin(4) * t62 + MDP(12) + 0.2e1 * t65; t38 * MDP(20) - t39 * MDP(21); t53 * MDP(19) + t40 * MDP(20) - t41 * MDP(21) + t63 * t56; -t78 * t53; t59; MDP(19);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
