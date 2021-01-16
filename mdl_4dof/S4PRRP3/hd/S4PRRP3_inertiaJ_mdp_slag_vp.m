% Calculate joint inertia matrix for
% S4PRRP3
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
%   see S4PRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 22:27
% Revision: beb2ba9bd8c5bd556f42a244985f3dab86917626 (2021-01-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP3_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 22:27:06
% EndTime: 2021-01-14 22:27:06
% DurationCPUTime: 0.05s
% Computational Cost: add. (42->30), mult. (79->42), div. (0->0), fcn. (47->2), ass. (0->13)
t21 = -qJ(4) - pkin(5);
t20 = MDP(15) * pkin(3);
t15 = cos(qJ(3));
t12 = -t15 * pkin(3) - pkin(2);
t19 = t12 * MDP(15);
t14 = sin(qJ(3));
t18 = t14 * MDP(13);
t17 = t15 * MDP(12);
t16 = MDP(12) + t20;
t13 = t14 ^ 2;
t11 = t21 * t15;
t10 = t21 * t14;
t1 = [MDP(1) + (t15 ^ 2 + t13) * MDP(15); (t15 * t10 - t14 * t11) * MDP(15); MDP(2) + t13 * MDP(5) + 0.2e1 * t14 * t15 * MDP(6) + 0.2e1 * (-t10 * t14 - t11 * t15) * MDP(14) + (t10 ^ 2 + t11 ^ 2) * MDP(15) + (-0.2e1 * t17 + 0.2e1 * t18 + t19) * t12 + 0.2e1 * (t15 * MDP(10) - t14 * MDP(11)) * pkin(2); (-MDP(11) - MDP(13)) * t14 + (MDP(10) + t16) * t15; t11 * MDP(13) + (-MDP(11) * pkin(5) + MDP(8)) * t15 + t16 * t10 + (-MDP(10) * pkin(5) - MDP(14) * pkin(3) + MDP(7)) * t14; MDP(9) + (0.2e1 * MDP(12) + t20) * pkin(3); 0; -t17 + t18 + t19; 0; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
