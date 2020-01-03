% Calculate joint inertia matrix for
% S4RPRP7
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
%   see S4RPRP7_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP7_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(17,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'S4RPRP7_inertiaJ_mdp_slag_vp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:17
% EndTime: 2019-12-31 16:47:17
% DurationCPUTime: 0.06s
% Computational Cost: add. (62->39), mult. (94->53), div. (0->0), fcn. (42->2), ass. (0->15)
t23 = cos(qJ(3));
t20 = t23 ^ 2;
t22 = sin(qJ(3));
t18 = t22 ^ 2 + t20;
t30 = t18 * MDP(17);
t29 = -0.2e1 * t23;
t28 = 2 * MDP(14);
t27 = (pkin(1) * MDP(6));
t26 = -MDP(13) + MDP(16);
t25 = -pkin(3) * MDP(17) - MDP(14);
t24 = -pkin(1) - pkin(5);
t17 = t23 * pkin(3) + t22 * qJ(4);
t16 = t22 * pkin(3) - t23 * qJ(4) + qJ(2);
t15 = t18 * t24;
t1 = [-0.2e1 * t15 * MDP(15) + t22 * MDP(8) * t29 + t20 * MDP(7) + MDP(1) + t24 ^ 2 * t30 + (MDP(16) * t29 + MDP(17) * t16 + t22 * t28) * t16 + ((-2 * MDP(4) + t27) * pkin(1)) + (0.2e1 * t22 * MDP(12) + 0.2e1 * t23 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2); -t18 * MDP(15) + t15 * MDP(17) + MDP(4) - t27; MDP(6) + t30; -t22 * MDP(10) - t17 * MDP(15) + t23 * MDP(9) + ((MDP(12) - t25) * t23 + (MDP(17) * qJ(4) + t26) * t22) * t24; t17 * MDP(17) + (MDP(12) + MDP(14)) * t23 + t26 * t22; MDP(11) + pkin(3) * t28 + 0.2e1 * qJ(4) * MDP(16) + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(17); (-MDP(17) * t24 + MDP(15)) * t23; -t23 * MDP(17); t25; MDP(17);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
