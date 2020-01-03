% Calculate joint inertia matrix for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% MDP [18x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(18,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'S5RPPPR4_inertiaJ_mdp_slag_vp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:16
% EndTime: 2019-12-31 17:45:16
% DurationCPUTime: 0.08s
% Computational Cost: add. (88->38), mult. (137->55), div. (0->0), fcn. (105->6), ass. (0->25)
t42 = sin(pkin(8));
t44 = cos(pkin(8));
t54 = t42 ^ 2 + t44 ^ 2;
t58 = MDP(11) * t54;
t45 = cos(pkin(7));
t39 = -t45 * pkin(1) - pkin(2);
t57 = t39 * MDP(7);
t43 = sin(pkin(7));
t37 = t43 * pkin(1) + qJ(3);
t56 = 0.2e1 * t42 * pkin(4) + 0.2e1 * t37;
t36 = -qJ(4) + t39;
t55 = -pkin(6) + t36;
t53 = t42 * MDP(8);
t52 = t44 * MDP(9);
t51 = MDP(7) + t58;
t46 = sin(qJ(5));
t47 = cos(qJ(5));
t30 = t47 * t42 + t46 * t44;
t50 = t30 * MDP(17);
t31 = -t46 * t42 + t47 * t44;
t49 = -t31 * MDP(18) - t50;
t29 = t55 * t44;
t28 = t55 * t42;
t27 = t54 * t36;
t1 = [t50 * t56 + MDP(1) + (t43 ^ 2 + t45 ^ 2) * MDP(4) * pkin(1) ^ 2 + t36 ^ 2 * t58 - 0.2e1 * t27 * MDP(10) + (MDP(12) * t31 - 0.2e1 * t30 * MDP(13) + MDP(18) * t56) * t31 + ((2 * MDP(6)) + 0.2e1 * t53 + 0.2e1 * t52 + (MDP(7) + MDP(11)) * t37) * t37 + ((2 * MDP(5)) + t57) * t39; 0; MDP(4) + t51; -MDP(10) * t54 + t27 * MDP(11) + MDP(5) + t57; 0; t51; t37 * MDP(11) - t49 + t52 + t53; 0; 0; MDP(11); t31 * MDP(14) - t30 * MDP(15) + (-t46 * t28 + t47 * t29) * MDP(17) + (-t47 * t28 - t46 * t29) * MDP(18); t49; t31 * MDP(17) - t30 * MDP(18); 0; MDP(16);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
