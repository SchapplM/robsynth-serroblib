% Calculate joint inertia matrix for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S4RRPR7_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:45
% EndTime: 2019-12-31 17:06:45
% DurationCPUTime: 0.14s
% Computational Cost: add. (157->56), mult. (317->95), div. (0->0), fcn. (303->6), ass. (0->30)
t56 = cos(qJ(2));
t68 = 0.2e1 * t56;
t65 = -qJ(3) - pkin(5);
t44 = t65 * t56;
t51 = sin(pkin(7));
t52 = cos(pkin(7));
t54 = sin(qJ(2));
t62 = t65 * t54;
t37 = -t51 * t44 - t52 * t62;
t42 = t51 * t56 + t52 * t54;
t67 = t37 * t42;
t53 = sin(qJ(4));
t55 = cos(qJ(4));
t66 = t53 * t55;
t41 = t51 * t54 - t52 * t56;
t64 = t41 * MDP(17);
t63 = MDP(14) * t66;
t48 = -t56 * pkin(2) - pkin(1);
t61 = t55 * MDP(18) - t53 * MDP(19);
t60 = MDP(18) * t53 + MDP(19) * t55;
t59 = (MDP(15) * t55 - MDP(16) * t53) * t42;
t58 = t53 * MDP(15) + t55 * MDP(16) - t60 * (t51 * pkin(2) + pkin(6));
t50 = t55 ^ 2;
t49 = t53 ^ 2;
t47 = -t52 * pkin(2) - pkin(3);
t39 = -t52 * t44 + t51 * t62;
t36 = t41 * pkin(3) - t42 * pkin(6) + t48;
t35 = t53 * t36 + t55 * t39;
t34 = t55 * t36 - t53 * t39;
t1 = [MDP(1) + pkin(1) * MDP(9) * t68 + (t37 ^ 2 + t39 ^ 2 + t48 ^ 2) * MDP(12) + (t50 * MDP(13) - 0.2e1 * t63) * t42 ^ 2 + (-0.2e1 * pkin(1) * MDP(10) + MDP(4) * t54 + MDP(5) * t68) * t54 + (0.2e1 * t59 + t64) * t41 + 0.2e1 * (-t39 * t41 + t67) * MDP(11) + 0.2e1 * (t34 * t41 + t53 * t67) * MDP(18) + 0.2e1 * (-t35 * t41 + t55 * t67) * MDP(19); t54 * MDP(6) + t56 * MDP(7) - t61 * t37 + t58 * t41 + (-t56 * MDP(10) - t54 * MDP(9)) * pkin(5) + (MDP(13) * t66 + (-t49 + t50) * MDP(14) + t60 * t47) * t42 + ((-t41 * t51 - t42 * t52) * MDP(11) + (-t37 * t52 + t39 * t51) * MDP(12)) * pkin(2); 0.2e1 * t63 + t49 * MDP(13) + MDP(8) + (t51 ^ 2 + t52 ^ 2) * MDP(12) * pkin(2) ^ 2 - 0.2e1 * t61 * t47; t48 * MDP(12) + t61 * t41; 0; MDP(12); t34 * MDP(18) - t35 * MDP(19) + t59 + t64; t58; t61; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
