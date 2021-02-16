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
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP5_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:36
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP5_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP5_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:36:08
% EndTime: 2021-01-14 22:36:09
% DurationCPUTime: 0.06s
% Computational Cost: add. (54->38), mult. (109->50), div. (0->0), fcn. (74->4), ass. (0->18)
t37 = qJ(4) + pkin(5);
t25 = sin(qJ(3));
t22 = t25 ^ 2;
t27 = cos(qJ(3));
t36 = t27 ^ 2 + t22;
t35 = MDP(15) * pkin(3);
t21 = -t27 * pkin(3) - pkin(2);
t34 = t21 * MDP(15);
t33 = t25 * MDP(13);
t32 = t27 * MDP(12);
t31 = -MDP(11) - MDP(13);
t30 = MDP(12) + t35;
t19 = t37 * t25;
t20 = t37 * t27;
t29 = t19 * t25 + t20 * t27;
t28 = cos(qJ(2));
t26 = sin(qJ(2));
t1 = [MDP(1) + (t36 * t26 ^ 2 + t28 ^ 2) * MDP(15); (t36 * MDP(14) + t29 * MDP(15) - MDP(4)) * t26 + (-t34 + MDP(3) + (MDP(10) + MDP(12)) * t27 + t31 * t25) * t28; MDP(2) + t22 * MDP(5) + 0.2e1 * t25 * t27 * MDP(6) + 0.2e1 * t29 * MDP(14) + (t19 ^ 2 + t20 ^ 2) * MDP(15) + (-0.2e1 * t32 + 0.2e1 * t33 + t34) * t21 + 0.2e1 * (t27 * MDP(10) - t25 * MDP(11)) * pkin(2); (t31 * t27 + (-MDP(10) - t30) * t25) * t26; -t20 * MDP(13) + (-MDP(11) * pkin(5) + MDP(8)) * t27 - t30 * t19 + (-MDP(10) * pkin(5) - MDP(14) * pkin(3) + MDP(7)) * t25; MDP(9) + (0.2e1 * MDP(12) + t35) * pkin(3); -t28 * MDP(15); -t32 + t33 + t34; 0; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
