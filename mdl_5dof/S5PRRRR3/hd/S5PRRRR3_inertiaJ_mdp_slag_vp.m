% Calculate joint inertia matrix for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S5PRRRR3_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:19
% EndTime: 2019-12-05 17:06:19
% DurationCPUTime: 0.11s
% Computational Cost: add. (104->49), mult. (182->55), div. (0->0), fcn. (124->6), ass. (0->31)
t44 = sin(qJ(5));
t47 = cos(qJ(5));
t56 = t47 * MDP(16);
t53 = -t44 * MDP(17) + t56;
t50 = 0.2e1 * t53;
t46 = sin(qJ(3));
t64 = pkin(2) * t46;
t63 = pkin(4) * t44;
t48 = cos(qJ(4));
t62 = t48 * pkin(3);
t61 = t44 * MDP(13) + t47 * MDP(14);
t49 = cos(qJ(3));
t39 = t49 * pkin(2) + pkin(3);
t34 = t48 * t39;
t45 = sin(qJ(4));
t32 = -t45 * t64 + t34;
t60 = t32 * MDP(9);
t59 = t49 * MDP(6);
t33 = -t45 * t39 - t48 * t64;
t58 = t33 * MDP(10);
t57 = t45 * MDP(10);
t55 = MDP(8) + (MDP(11) * t44 + 0.2e1 * MDP(12) * t47) * t44;
t54 = MDP(5) + t55;
t52 = -MDP(16) * t44 - MDP(17) * t47;
t51 = (t48 * MDP(9) - t57) * pkin(3);
t43 = pkin(4) * t47;
t38 = -pkin(4) - t62;
t35 = t38 * t44;
t30 = -pkin(4) - t32;
t29 = t30 * t44;
t1 = [MDP(1); 0; MDP(2) - t30 * t50 + 0.2e1 * (-t46 * MDP(7) + t59) * pkin(2) + 0.2e1 * t58 + 0.2e1 * t60 + t54; 0; (t34 + t62) * MDP(9) + (t35 + t29) * MDP(17) + (-t30 - t38) * t56 + (-pkin(3) - t39) * t57 + (t59 + (-MDP(10) * t48 - MDP(9) * t45 - MDP(7)) * t46) * pkin(2) + t54; -t38 * t50 + 0.2e1 * t51 + t54; 0; t60 + t58 + (-t30 * t47 + t43) * MDP(16) + (t29 - t63) * MDP(17) + t55; (-t38 * t47 + t43) * MDP(16) + (t35 - t63) * MDP(17) + t51 + t55; pkin(4) * t50 + t55; t53; t52 * (pkin(8) - t33) + t61; t52 * (t45 * pkin(3) + pkin(8)) + t61; t52 * pkin(8) + t61; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
