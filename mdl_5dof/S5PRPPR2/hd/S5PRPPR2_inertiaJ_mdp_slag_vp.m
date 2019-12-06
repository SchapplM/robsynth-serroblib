% Calculate joint inertia matrix for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PRPPR2_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:24:42
% EndTime: 2019-12-05 15:24:43
% DurationCPUTime: 0.12s
% Computational Cost: add. (102->39), mult. (202->65), div. (0->0), fcn. (199->8), ass. (0->32)
t48 = sin(pkin(8));
t41 = t48 * pkin(2) + qJ(4);
t47 = sin(pkin(9));
t49 = cos(pkin(9));
t64 = t47 ^ 2 + t49 ^ 2;
t58 = t64 * MDP(9);
t69 = t41 * t58;
t50 = cos(pkin(8));
t52 = sin(qJ(2));
t54 = cos(qJ(2));
t34 = t48 * t52 - t50 * t54;
t68 = t34 ^ 2;
t44 = -t50 * pkin(2) - pkin(3);
t67 = -0.2e1 * t49 * pkin(4) + 0.2e1 * t44;
t66 = pkin(6) + t41;
t65 = MDP(5) * pkin(2);
t63 = t44 * MDP(9);
t62 = t47 * MDP(7);
t61 = t49 * MDP(6);
t51 = sin(qJ(5));
t53 = cos(qJ(5));
t36 = t51 * t47 - t53 * t49;
t60 = t36 * MDP(15);
t59 = t64 * MDP(8);
t38 = t53 * t47 + t51 * t49;
t57 = -t38 * MDP(16) - t60;
t56 = -t57 - t61 + t62 + t63;
t37 = t48 * t54 + t50 * t52;
t33 = t37 ^ 2;
t32 = t66 * t49;
t31 = t66 * t47;
t1 = [MDP(1) + (t33 + t68) * MDP(5) + (t64 * t33 + t68) * MDP(9); t54 * MDP(3) - t52 * MDP(4) + (t48 * t65 + t59 + t69) * t37 + (-t50 * t65 + t56) * t34; t60 * t67 + MDP(2) + (t48 ^ 2 + t50 ^ 2) * MDP(5) * pkin(2) ^ 2 + (-0.2e1 * t61 + 0.2e1 * t62 + t63) * t44 + (MDP(10) * t38 - 0.2e1 * t36 * MDP(11) + MDP(16) * t67) * t38 + (0.2e1 * t59 + t69) * t41; 0; 0; MDP(5) + t58; t34 * MDP(9); t56; 0; MDP(9); (-t38 * MDP(15) + t36 * MDP(16)) * t37; t38 * MDP(12) - t36 * MDP(13) + (-t53 * t31 - t51 * t32) * MDP(15) + (t51 * t31 - t53 * t32) * MDP(16); t57; 0; MDP(14);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
