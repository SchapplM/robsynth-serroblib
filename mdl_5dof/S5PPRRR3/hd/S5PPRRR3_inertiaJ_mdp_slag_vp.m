% Calculate joint inertia matrix for
% S5PPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% MDP [19x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PPRRR3_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:17:00
% EndTime: 2019-12-05 15:17:01
% DurationCPUTime: 0.12s
% Computational Cost: add. (88->42), mult. (195->69), div. (0->0), fcn. (201->8), ass. (0->26)
t54 = cos(qJ(4));
t66 = -0.2e1 * pkin(4) * t54 - (2 * pkin(3));
t65 = pkin(6) + pkin(7);
t48 = sin(pkin(9));
t55 = cos(qJ(3));
t64 = t48 * t55;
t49 = cos(pkin(9));
t51 = sin(qJ(4));
t38 = -t49 * t54 - t51 * t64;
t39 = -t49 * t51 + t54 * t64;
t50 = sin(qJ(5));
t53 = cos(qJ(5));
t63 = (t38 * t53 - t39 * t50) * MDP(18) + (-t38 * t50 - t39 * t53) * MDP(19);
t42 = t50 * t51 - t53 * t54;
t43 = t50 * t54 + t51 * t53;
t52 = sin(qJ(3));
t62 = (-MDP(18) * t43 + MDP(19) * t42) * t52;
t61 = MDP(11) * t54;
t60 = MDP(18) * t42;
t44 = t65 * t51;
t45 = t65 * t54;
t59 = t43 * MDP(15) - t42 * MDP(16) + (-t44 * t53 - t45 * t50) * MDP(18) + (t44 * t50 - t45 * t53) * MDP(19);
t58 = -MDP(11) * t51 - MDP(12) * t54;
t57 = (MDP(18) * t53 - MDP(19) * t50) * pkin(4);
t56 = -MDP(12) * t51 - MDP(19) * t43 + MDP(4) - t60 + t61;
t1 = [MDP(1) + (t48 ^ 2 + t49 ^ 2) * MDP(2); 0; MDP(2); (-MDP(5) * t55 - t52 * t56) * t48; -t52 * MDP(5) + t55 * t56; 0.2e1 * pkin(3) * t61 + t60 * t66 + MDP(3) + (-0.2e1 * MDP(12) * pkin(3) + MDP(6) * t51 + 0.2e1 * MDP(7) * t54) * t51 + (MDP(13) * t43 - 0.2e1 * MDP(14) * t42 + MDP(19) * t66) * t43; MDP(11) * t38 - MDP(12) * t39 + t63; t52 * t58 + t62; t51 * MDP(8) + t54 * MDP(9) + pkin(6) * t58 + t59; MDP(10) + MDP(17) + 0.2e1 * t57; t63; t62; t59; MDP(17) + t57; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
