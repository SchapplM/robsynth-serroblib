% Calculate joint inertia matrix for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S5PPRPR1_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:01:20
% EndTime: 2019-12-05 15:01:21
% DurationCPUTime: 0.08s
% Computational Cost: add. (79->33), mult. (169->58), div. (0->0), fcn. (170->8), ass. (0->26)
t38 = cos(pkin(9));
t55 = -0.2e1 * t38 * pkin(4) - (2 * pkin(3));
t54 = pkin(3) * MDP(9);
t53 = pkin(6) + qJ(4);
t36 = sin(pkin(9));
t52 = t36 ^ 2 + t38 ^ 2;
t51 = t36 * MDP(7);
t50 = t38 * MDP(6);
t40 = sin(qJ(5));
t42 = cos(qJ(5));
t26 = t40 * t36 - t42 * t38;
t49 = t26 * MDP(15);
t48 = t52 * MDP(9);
t47 = t52 * qJ(4);
t28 = t42 * t36 + t40 * t38;
t46 = -t28 * MDP(16) - t49;
t45 = -t46 - t50 + t51 - t54;
t43 = cos(qJ(3));
t41 = sin(qJ(3));
t39 = cos(pkin(8));
t37 = sin(pkin(8));
t31 = t53 * t38;
t30 = t53 * t36;
t29 = t43 * t37 + t41 * t39;
t27 = t41 * t37 - t43 * t39;
t1 = [MDP(1) + (t37 ^ 2 + t39 ^ 2) * MDP(2) + (t52 * t29 ^ 2 + t27 ^ 2) * MDP(9); 0; MDP(2) + t48; (t52 * MDP(8) + MDP(9) * t47 - MDP(5)) * t29 + (-MDP(4) + t45) * t27; 0; t49 * t55 + MDP(3) + qJ(4) ^ 2 * t48 + 0.2e1 * MDP(8) * t47 + (0.2e1 * t50 - 0.2e1 * t51 + t54) * pkin(3) + (MDP(10) * t28 - 0.2e1 * t26 * MDP(11) + MDP(16) * t55) * t28; t27 * MDP(9); 0; t45; MDP(9); (-t28 * MDP(15) + t26 * MDP(16)) * t29; t46; t28 * MDP(12) - t26 * MDP(13) + (-t42 * t30 - t40 * t31) * MDP(15) + (t40 * t30 - t42 * t31) * MDP(16); 0; MDP(14);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
