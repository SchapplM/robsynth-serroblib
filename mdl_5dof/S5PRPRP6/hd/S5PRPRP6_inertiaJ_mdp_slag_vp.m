% Calculate joint inertia matrix for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5PRPRP6_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:41:19
% EndTime: 2019-12-05 15:41:20
% DurationCPUTime: 0.12s
% Computational Cost: add. (92->54), mult. (156->69), div. (0->0), fcn. (93->4), ass. (0->27)
t36 = cos(qJ(4));
t31 = t36 ^ 2;
t34 = sin(qJ(4));
t28 = t34 ^ 2 + t31;
t53 = t28 * MDP(18);
t41 = -pkin(4) * MDP(18) - MDP(15);
t43 = MDP(14) - MDP(17);
t52 = t34 * (MDP(18) * qJ(5) - t43) + t36 * (MDP(13) - t41);
t51 = -0.2e1 * t36;
t50 = 0.2e1 * MDP(15);
t49 = (pkin(2) * MDP(7));
t38 = -pkin(2) - pkin(6);
t48 = MDP(18) * t38;
t47 = qJ(3) * MDP(7);
t26 = t34 * pkin(4) - t36 * qJ(5) + qJ(3);
t46 = t26 * MDP(18);
t45 = t36 * MDP(18);
t44 = MDP(13) + MDP(15);
t42 = -MDP(5) + t49;
t37 = cos(qJ(2));
t35 = sin(qJ(2));
t32 = t37 ^ 2;
t30 = t35 ^ 2;
t27 = t36 * pkin(4) + t34 * qJ(5);
t25 = t28 * t38;
t24 = t28 * t37;
t1 = [MDP(1) + (t30 + t32) * MDP(7) + (t28 * t32 + t30) * MDP(18); t24 * MDP(16) + (-t28 * t48 + MDP(3) + t42) * t37 + (t44 * t34 + t43 * t36 - MDP(4) + MDP(6) + t46 + t47) * t35; -0.2e1 * t25 * MDP(16) + t34 * MDP(9) * t51 + t31 * MDP(8) + MDP(2) + t38 ^ 2 * t53 + (MDP(17) * t51 + t34 * t50 + t46) * t26 + ((-2 * MDP(5) + t49) * pkin(2)) + (0.2e1 * t34 * MDP(13) + 0.2e1 * t36 * MDP(14) + 0.2e1 * MDP(6) + t47) * qJ(3); -t24 * MDP(18) - t37 * MDP(7); -t28 * MDP(16) + t25 * MDP(18) - t42; MDP(7) + t53; -t52 * t37; t36 * MDP(10) - t34 * MDP(11) - t27 * MDP(16) + t52 * t38; t27 * MDP(18) - t43 * t34 + t44 * t36; MDP(12) + pkin(4) * t50 + 0.2e1 * qJ(5) * MDP(17) + (pkin(4) ^ 2 + qJ(5) ^ 2) * MDP(18); t37 * t45; (MDP(16) - t48) * t36; -t45; t41; MDP(18);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
