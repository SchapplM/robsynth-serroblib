% Calculate joint inertia matrix for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PRPRR2_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:13
% EndTime: 2019-12-05 15:45:14
% DurationCPUTime: 0.08s
% Computational Cost: add. (101->35), mult. (196->52), div. (0->0), fcn. (195->8), ass. (0->24)
t45 = sin(qJ(5));
t48 = cos(qJ(5));
t54 = t48 * MDP(14) - t45 * MDP(15);
t43 = sin(pkin(9));
t62 = pkin(2) * t43;
t60 = t45 * MDP(11) + t48 * MDP(12);
t44 = cos(pkin(9));
t39 = t44 * pkin(2) + pkin(3);
t46 = sin(qJ(4));
t49 = cos(qJ(4));
t34 = t49 * t39 - t46 * t62;
t59 = t34 * MDP(7);
t35 = -t46 * t39 - t49 * t62;
t58 = t35 * MDP(8);
t55 = MDP(6) + (0.2e1 * MDP(10) * t48 + MDP(9) * t45) * t45;
t53 = -MDP(14) * t45 - MDP(15) * t48;
t47 = sin(qJ(2));
t50 = cos(qJ(2));
t36 = -t43 * t47 + t44 * t50;
t37 = t43 * t50 + t44 * t47;
t31 = t46 * t36 + t49 * t37;
t52 = -t31 * MDP(8) + (-MDP(7) - t54) * (-t49 * t36 + t46 * t37);
t32 = -pkin(4) - t34;
t1 = [MDP(1) + (t36 ^ 2 + t37 ^ 2) * MDP(5); t50 * MDP(3) - t47 * MDP(4) + (t36 * t44 + t37 * t43) * MDP(5) * pkin(2) + t52; MDP(2) + (t43 ^ 2 + t44 ^ 2) * MDP(5) * pkin(2) ^ 2 - 0.2e1 * t54 * t32 + 0.2e1 * t59 + 0.2e1 * t58 + t55; 0; 0; MDP(5); t52; t55 + t58 + t59 + t54 * (pkin(4) - t32); 0; 0.2e1 * pkin(4) * t54 + t55; t53 * t31; t53 * (pkin(7) - t35) + t60; t54; t53 * pkin(7) + t60; MDP(13);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
