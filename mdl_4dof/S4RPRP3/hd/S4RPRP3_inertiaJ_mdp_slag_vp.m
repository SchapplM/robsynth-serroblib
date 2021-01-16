% Calculate joint inertia matrix for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 10:20
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPRP3_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 10:20:37
% EndTime: 2021-01-15 10:20:37
% DurationCPUTime: 0.08s
% Computational Cost: add. (61->34), mult. (103->49), div. (0->0), fcn. (66->4), ass. (0->17)
t33 = MDP(15) * pkin(3);
t22 = sin(pkin(6));
t19 = t22 * pkin(1) + pkin(5);
t32 = qJ(4) + t19;
t25 = cos(qJ(3));
t23 = cos(pkin(6));
t28 = -t23 * pkin(1) - pkin(2);
t18 = -t25 * pkin(3) + t28;
t31 = t18 * MDP(15);
t24 = sin(qJ(3));
t30 = t24 * MDP(13);
t29 = t25 * MDP(12);
t27 = MDP(12) + t33;
t21 = t24 ^ 2;
t17 = t32 * t25;
t16 = t32 * t24;
t1 = [MDP(1) + t21 * MDP(5) + 0.2e1 * t24 * t25 * MDP(6) + 0.2e1 * (t16 * t24 + t17 * t25) * MDP(14) + (t16 ^ 2 + t17 ^ 2) * MDP(15) + (t22 ^ 2 + t23 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * (-t25 * MDP(10) + t24 * MDP(11)) * t28 + (-0.2e1 * t29 + 0.2e1 * t30 + t31) * t18; (-t16 * t25 + t17 * t24) * MDP(15); MDP(4) + (t25 ^ 2 + t21) * MDP(15); -t17 * MDP(13) + (-MDP(11) * t19 + MDP(8)) * t25 - t27 * t16 + (-MDP(10) * t19 - MDP(14) * pkin(3) + MDP(7)) * t24; (-MDP(11) - MDP(13)) * t24 + (MDP(10) + t27) * t25; MDP(9) + (0.2e1 * MDP(12) + t33) * pkin(3); -t29 + t30 + t31; 0; 0; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
