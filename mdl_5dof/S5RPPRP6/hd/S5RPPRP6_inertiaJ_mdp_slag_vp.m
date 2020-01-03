% Calculate joint inertia matrix for
% S5RPPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
% MDP [21x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(21,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [21 1]), ...
  'S5RPPRP6_inertiaJ_mdp_slag_vp: MDP has to be [21x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:17
% EndTime: 2019-12-31 17:55:17
% DurationCPUTime: 0.14s
% Computational Cost: add. (189->60), mult. (287->82), div. (0->0), fcn. (262->4), ass. (0->27)
t48 = sin(pkin(7));
t49 = cos(pkin(7));
t51 = sin(qJ(4));
t60 = t51 * t49;
t61 = cos(qJ(4));
t38 = t61 * t48 + t60;
t64 = 0.2e1 * t38;
t63 = 2 * MDP(20);
t50 = -pkin(1) - qJ(3);
t62 = -pkin(6) + t50;
t41 = t48 ^ 2 + t49 ^ 2;
t42 = t48 * pkin(3) + qJ(2);
t59 = -MDP(17) + MDP(20);
t57 = t61 * t49;
t37 = t51 * t48 - t57;
t58 = t37 ^ 2 + t38 ^ 2;
t56 = -pkin(4) * MDP(21) - MDP(18);
t55 = -MDP(16) + t56;
t54 = t48 * MDP(7) + t49 * MDP(8);
t53 = MDP(21) * qJ(5) + t59;
t52 = qJ(2) ^ 2;
t40 = t62 * t48;
t36 = t41 * t50;
t34 = t61 * t40 + t62 * t60;
t33 = t51 * t40 - t62 * t57;
t32 = t38 * pkin(4) + t37 * qJ(5) + t42;
t1 = [MDP(1) - (2 * pkin(1) * MDP(4)) + ((pkin(1) ^ 2) + t52) * MDP(6) - 0.2e1 * t36 * MDP(9) + (t41 * t50 ^ 2 + t52) * MDP(10) + (t32 ^ 2 + t33 ^ 2 + t34 ^ 2) * MDP(21) + (MDP(11) * t37 + MDP(12) * t64 - 0.2e1 * t42 * MDP(17) - 0.2e1 * MDP(19) * t33 + t32 * t63) * t37 + (t42 * MDP(16) + t32 * MDP(18) - MDP(19) * t34) * t64 + 0.2e1 * (MDP(5) + t54) * qJ(2); MDP(4) - pkin(1) * MDP(6) - t41 * MDP(9) + t36 * MDP(10) - t58 * MDP(19) + (t33 * t37 + t34 * t38) * MDP(21); t41 * MDP(10) + t58 * MDP(21) + MDP(6); qJ(2) * MDP(10) + t32 * MDP(21) + (MDP(16) + MDP(18)) * t38 + t59 * t37 + t54; 0; MDP(10) + MDP(21); -t37 * MDP(13) - t38 * MDP(14) + (t37 * pkin(4) - t38 * qJ(5)) * MDP(19) + t53 * t34 + t55 * t33; t55 * t37 + t53 * t38; 0; MDP(15) + 0.2e1 * pkin(4) * MDP(18) + (qJ(5) * t63) + (pkin(4) ^ 2 + (qJ(5) ^ 2)) * MDP(21); -t37 * MDP(19) + t33 * MDP(21); t37 * MDP(21); 0; t56; MDP(21);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
