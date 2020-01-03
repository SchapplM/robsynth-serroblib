% Calculate joint inertia matrix for
% S4PRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4PRPR5_inertiaJ_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:18
% EndTime: 2019-12-31 16:23:18
% DurationCPUTime: 0.04s
% Computational Cost: add. (30->20), mult. (64->37), div. (0->0), fcn. (57->6), ass. (0->13)
t25 = cos(qJ(4));
t30 = t25 * MDP(11);
t23 = sin(qJ(4));
t29 = -t23 * MDP(12) + t30;
t28 = -MDP(11) * t23 - MDP(12) * t25;
t26 = cos(qJ(2));
t24 = sin(qJ(2));
t22 = cos(pkin(7));
t21 = sin(pkin(7));
t20 = -t22 * pkin(2) - pkin(3);
t17 = t21 * t26 + t22 * t24;
t15 = t21 * t24 - t22 * t26;
t1 = [MDP(1) + (t15 ^ 2 + t17 ^ 2) * MDP(5); t26 * MDP(3) - t24 * MDP(4) - t29 * t15 + (-t15 * t22 + t17 * t21) * MDP(5) * pkin(2); -0.2e1 * t20 * t30 + MDP(2) + (t21 ^ 2 + t22 ^ 2) * MDP(5) * pkin(2) ^ 2 + (0.2e1 * t20 * MDP(12) + MDP(6) * t23 + 0.2e1 * t25 * MDP(7)) * t23; 0; 0; MDP(5); t28 * t17; t23 * MDP(8) + t25 * MDP(9) + t28 * (t21 * pkin(2) + pkin(5)); t29; MDP(10);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
