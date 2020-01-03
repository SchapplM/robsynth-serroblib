% Calculate joint inertia matrix for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S4PRRP5_inertiaJ_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:29:22
% EndTime: 2019-12-31 16:29:22
% DurationCPUTime: 0.08s
% Computational Cost: add. (39->28), mult. (83->43), div. (0->0), fcn. (56->4), ass. (0->15)
t31 = -qJ(4) - pkin(5);
t22 = sin(qJ(3));
t19 = t22 ^ 2;
t24 = cos(qJ(3));
t30 = t24 ^ 2 + t19;
t29 = MDP(13) * pkin(3);
t18 = -t24 * pkin(3) - pkin(2);
t28 = t18 * MDP(13);
t16 = t31 * t22;
t17 = t31 * t24;
t27 = -t16 * t22 - t17 * t24;
t26 = t24 * MDP(10) - t22 * MDP(11);
t25 = cos(qJ(2));
t23 = sin(qJ(2));
t1 = [MDP(1) + (t30 * t23 ^ 2 + t25 ^ 2) * MDP(13); (MDP(3) + t26 - t28) * t25 + (t30 * MDP(12) + t27 * MDP(13) - MDP(4)) * t23; MDP(2) + t19 * MDP(5) + 0.2e1 * t22 * t24 * MDP(6) + 0.2e1 * t27 * MDP(12) + (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) * MDP(13) + 0.2e1 * t26 * pkin(2); (-MDP(11) * t24 + (-MDP(10) - t29) * t22) * t23; t16 * t29 + (-MDP(11) * pkin(5) + MDP(8)) * t24 + (-MDP(10) * pkin(5) - MDP(12) * pkin(3) + MDP(7)) * t22; MDP(13) * pkin(3) ^ 2 + MDP(9); -t25 * MDP(13); t28; 0; MDP(13);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
