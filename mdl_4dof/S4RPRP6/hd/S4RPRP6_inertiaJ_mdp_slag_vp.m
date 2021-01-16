% Calculate joint inertia matrix for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% MDP [17x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:27
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RPRP6_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:27:41
% EndTime: 2021-01-15 10:27:41
% DurationCPUTime: 0.05s
% Computational Cost: add. (61->37), mult. (90->47), div. (0->0), fcn. (49->2), ass. (0->17)
t29 = (pkin(1) * MDP(6));
t28 = MDP(17) * pkin(3);
t22 = -pkin(1) - pkin(5);
t27 = -qJ(4) + t22;
t20 = sin(qJ(3));
t18 = t20 * pkin(3) + qJ(2);
t26 = t18 * MDP(17);
t25 = t20 * MDP(14);
t21 = cos(qJ(3));
t24 = t21 * MDP(15);
t23 = MDP(14) + t28;
t19 = t21 ^ 2;
t17 = t20 ^ 2 + t19;
t16 = t27 * t21;
t15 = t27 * t20;
t14 = t15 * t20 + t16 * t21;
t1 = [MDP(1) + t19 * MDP(7) - 0.2e1 * t21 * t20 * MDP(8) - 0.2e1 * t14 * MDP(16) + (t15 ^ 2 + t16 ^ 2) * MDP(17) + (0.2e1 * t24 + 0.2e1 * t25 + t26) * t18 + ((-2 * MDP(4) + t29) * pkin(1)) + (0.2e1 * t20 * MDP(12) + 0.2e1 * t21 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2); -t17 * MDP(16) + t14 * MDP(17) + MDP(4) - t29; t17 * MDP(17) + MDP(6); -t15 * MDP(15) + (-MDP(13) * t22 - MDP(10)) * t20 + t23 * t16 + (MDP(12) * t22 - MDP(16) * pkin(3) + MDP(9)) * t21; (-MDP(13) - MDP(15)) * t20 + (MDP(12) + t23) * t21; MDP(11) + (0.2e1 * MDP(14) + t28) * pkin(3); t24 + t25 + t26; 0; 0; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
