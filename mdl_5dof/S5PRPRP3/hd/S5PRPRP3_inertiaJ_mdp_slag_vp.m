% Calculate joint inertia matrix for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S5PRPRP3_inertiaJ_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:41
% EndTime: 2019-12-05 15:33:42
% DurationCPUTime: 0.08s
% Computational Cost: add. (85->40), mult. (164->65), div. (0->0), fcn. (145->6), ass. (0->26)
t39 = sin(pkin(8));
t40 = cos(pkin(8));
t42 = sin(qJ(2));
t44 = cos(qJ(2));
t30 = t39 * t42 - t40 * t44;
t56 = t30 ^ 2;
t55 = MDP(5) * pkin(2);
t41 = sin(qJ(4));
t37 = t41 ^ 2;
t43 = cos(qJ(4));
t54 = t43 ^ 2 + t37;
t53 = MDP(14) * pkin(4);
t35 = t39 * pkin(2) + pkin(6);
t52 = qJ(5) + t35;
t49 = -t40 * pkin(2) - pkin(3);
t33 = -t43 * pkin(4) + t49;
t51 = t33 * MDP(14);
t50 = t41 * MDP(12);
t48 = MDP(11) + t53;
t27 = t52 * t41;
t28 = t52 * t43;
t47 = t27 * t41 + t28 * t43;
t46 = -t43 * MDP(11) + t50;
t32 = t39 * t44 + t40 * t42;
t29 = t32 ^ 2;
t1 = [MDP(1) + (t29 + t56) * MDP(5) + (t54 * t29 + t56) * MDP(14); t44 * MDP(3) - t42 * MDP(4) + (-t40 * t55 + t46 + t51) * t30 + (t54 * MDP(13) + t47 * MDP(14) + t39 * t55) * t32; MDP(2) + t37 * MDP(6) + 0.2e1 * t41 * t43 * MDP(7) + 0.2e1 * t47 * MDP(13) + (t27 ^ 2 + t28 ^ 2 + t33 ^ 2) * MDP(14) + (t39 ^ 2 + t40 ^ 2) * MDP(5) * pkin(2) ^ 2 + 0.2e1 * t46 * t49; 0; (-t27 * t43 + t28 * t41) * MDP(14); t54 * MDP(14) + MDP(5); (-MDP(12) * t43 - t48 * t41) * t32; -t27 * t53 + (-MDP(12) * t35 + MDP(9)) * t43 + (-MDP(11) * t35 - MDP(13) * pkin(4) + MDP(8)) * t41; t48 * t43 - t50; pkin(4) ^ 2 * MDP(14) + MDP(10); t30 * MDP(14); t51; 0; 0; MDP(14);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
