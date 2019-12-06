% Calculate joint inertia matrix for
% S5PPRRR2
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
%   see S5PPRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(19,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [19 1]), ...
  'S5PPRRR2_inertiaJ_mdp_slag_vp: MDP has to be [19x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:49
% EndTime: 2019-12-05 15:14:50
% DurationCPUTime: 0.10s
% Computational Cost: add. (87->36), mult. (177->57), div. (0->0), fcn. (187->8), ass. (0->24)
t45 = sin(qJ(4));
t48 = cos(qJ(4));
t54 = t48 * MDP(11);
t44 = sin(qJ(5));
t47 = cos(qJ(5));
t36 = t44 * t45 - t47 * t48;
t30 = t36 * MDP(18);
t37 = t44 * t48 + t47 * t45;
t55 = -t37 * MDP(19) - t30;
t59 = -t45 * MDP(12) + t54 + t55;
t58 = -0.2e1 * t48 * pkin(4) - (2 * pkin(3));
t57 = pkin(6) + pkin(7);
t42 = sin(pkin(9));
t43 = cos(pkin(9));
t46 = sin(qJ(3));
t49 = cos(qJ(3));
t35 = t49 * t42 + t46 * t43;
t56 = (-MDP(18) * t37 + MDP(19) * t36) * t35;
t38 = t57 * t45;
t39 = t57 * t48;
t53 = t37 * MDP(15) - t36 * MDP(16) + (-t47 * t38 - t44 * t39) * MDP(18) + (t44 * t38 - t47 * t39) * MDP(19);
t51 = -MDP(11) * t45 - MDP(12) * t48;
t50 = (MDP(18) * t47 - MDP(19) * t44) * pkin(4);
t1 = [MDP(1) + (t42 ^ 2 + t43 ^ 2) * MDP(2); 0; MDP(2); -t35 * MDP(5) + (-MDP(4) - t59) * (t46 * t42 - t49 * t43); 0; 0.2e1 * pkin(3) * t54 + t30 * t58 + MDP(3) + (-0.2e1 * MDP(12) * pkin(3) + MDP(6) * t45 + 0.2e1 * MDP(7) * t48) * t45 + (MDP(13) * t37 - 0.2e1 * MDP(14) * t36 + MDP(19) * t58) * t37; t51 * t35 + t56; t59; t45 * MDP(8) + t48 * MDP(9) + t51 * pkin(6) + t53; MDP(10) + MDP(17) + 0.2e1 * t50; t56; t55; t53; MDP(17) + t50; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
