% Calculate joint inertia matrix for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% MDP [16x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(16,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [16 1]), ...
  'S4RPPR5_inertiaJ_mdp_slag_vp: MDP has to be [16x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:49
% EndTime: 2019-12-31 16:39:49
% DurationCPUTime: 0.04s
% Computational Cost: add. (50->29), mult. (72->42), div. (0->0), fcn. (45->4), ass. (0->12)
t23 = sin(pkin(6));
t24 = cos(pkin(6));
t27 = -pkin(1) - pkin(2);
t19 = t24 * qJ(2) + t23 * t27;
t26 = cos(qJ(4));
t30 = t26 * MDP(15);
t17 = t23 * qJ(2) - t24 * t27;
t25 = sin(qJ(4));
t29 = -t25 * MDP(16) + t30;
t28 = -MDP(15) * t25 - MDP(16) * t26;
t15 = pkin(3) + t17;
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + (t17 ^ 2 + t19 ^ 2) * MDP(9) + 0.2e1 * t15 * t30 + 0.2e1 * t17 * MDP(7) + 0.2e1 * t19 * MDP(8) + (MDP(10) * t25 + 0.2e1 * t26 * MDP(11) - 0.2e1 * t15 * MDP(16)) * t25; -pkin(1) * MDP(6) - MDP(4) + (t19 * MDP(9) + MDP(8)) * t23 + (-t17 * MDP(9) - MDP(7) - t29) * t24; MDP(6) + (t23 ^ 2 + t24 ^ 2) * MDP(9); 0; 0; MDP(9); -t25 * MDP(12) - t26 * MDP(13) + t28 * (-pkin(5) + t19); t28 * t23; t29; MDP(14);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
