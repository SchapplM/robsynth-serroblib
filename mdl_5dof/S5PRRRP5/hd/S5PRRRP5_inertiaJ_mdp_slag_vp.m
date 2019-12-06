% Calculate joint inertia matrix for
% S5PRRRP5
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
%   see S5PRRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(20,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [20 1]), ...
  'S5PRRRP5_inertiaJ_mdp_slag_vp: MDP has to be [20x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:49:18
% EndTime: 2019-12-05 16:49:19
% DurationCPUTime: 0.17s
% Computational Cost: add. (163->63), mult. (325->98), div. (0->0), fcn. (318->6), ass. (0->30)
t67 = pkin(6) + pkin(7);
t49 = sin(qJ(4));
t66 = pkin(3) * t49;
t50 = sin(qJ(3));
t52 = cos(qJ(4));
t53 = cos(qJ(3));
t42 = t49 * t53 + t52 * t50;
t51 = sin(qJ(2));
t37 = t42 * t51;
t41 = t49 * t50 - t52 * t53;
t38 = t41 * t51;
t65 = -t37 * MDP(17) + t38 * MDP(18);
t64 = MDP(20) * pkin(4);
t60 = -t53 * pkin(3) - pkin(2);
t34 = t41 * pkin(4) + t60;
t63 = t34 * MDP(20);
t62 = t52 * MDP(17);
t61 = t53 * MDP(10);
t44 = t67 * t50;
t45 = t67 * t53;
t59 = -t52 * t44 - t49 * t45;
t57 = t49 * t44 - t52 * t45;
t58 = t42 * MDP(14) - t41 * MDP(15) + t59 * MDP(17) + t57 * MDP(18);
t56 = -t50 * MDP(10) - t53 * MDP(11);
t55 = t41 * MDP(17) + t42 * MDP(18);
t54 = cos(qJ(2));
t47 = t52 * pkin(3) + pkin(4);
t31 = -t41 * qJ(5) - t57;
t30 = -t42 * qJ(5) + t59;
t1 = [MDP(1) + (t37 ^ 2 + t38 ^ 2 + t54 ^ 2) * MDP(20); -t51 * MDP(4) + (t37 * t42 + t38 * t41) * MDP(19) + (-t37 * t30 - t38 * t31) * MDP(20) + (-t50 * MDP(11) + MDP(3) - t55 + t61 - t63) * t54; MDP(2) + 0.2e1 * pkin(2) * t61 + t42 ^ 2 * MDP(12) - 0.2e1 * t42 * t41 * MDP(13) + 0.2e1 * (-t30 * t42 - t31 * t41) * MDP(19) + (t30 ^ 2 + t31 ^ 2 + t34 ^ 2) * MDP(20) + 0.2e1 * t55 * t60 + (-0.2e1 * pkin(2) * MDP(11) + MDP(5) * t50 + 0.2e1 * t53 * MDP(6)) * t50; (-t37 * t47 - t38 * t66) * MDP(20) + t56 * t51 + t65; t50 * MDP(7) + t53 * MDP(8) + (-t41 * t66 - t47 * t42) * MDP(19) + (t30 * t47 + t31 * t66) * MDP(20) + t56 * pkin(6) + t58; t47 ^ 2 * MDP(20) + MDP(16) + MDP(9) + (0.2e1 * t62 + (MDP(20) * t66 - 0.2e1 * MDP(18)) * t49) * pkin(3); -t37 * t64 + t65; (-t42 * MDP(19) + t30 * MDP(20)) * pkin(4) + t58; t47 * t64 + MDP(16) + (-t49 * MDP(18) + t62) * pkin(3); pkin(4) ^ 2 * MDP(20) + MDP(16); -t54 * MDP(20); t63; 0; 0; MDP(20);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
