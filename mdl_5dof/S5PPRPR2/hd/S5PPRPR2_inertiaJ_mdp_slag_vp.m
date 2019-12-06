% Calculate joint inertia matrix for
% S5PPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S5PPRPR2_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:21
% EndTime: 2019-12-05 15:03:22
% DurationCPUTime: 0.08s
% Computational Cost: add. (42->27), mult. (77->34), div. (0->0), fcn. (71->6), ass. (0->17)
t28 = (MDP(8) * pkin(3));
t27 = qJ(4) * MDP(8);
t17 = sin(qJ(5));
t26 = t17 * MDP(14);
t19 = cos(qJ(5));
t25 = t19 * MDP(15);
t24 = MDP(6) - t28;
t23 = t19 * MDP(14) - t17 * MDP(15);
t22 = -t25 - t26;
t21 = -pkin(3) - pkin(6);
t20 = cos(qJ(3));
t18 = sin(qJ(3));
t16 = cos(pkin(8));
t15 = sin(pkin(8));
t13 = t20 * t15 + t18 * t16;
t12 = t18 * t15 - t20 * t16;
t1 = [MDP(1) + (t15 ^ 2 + t16 ^ 2) * MDP(2) + (t12 ^ 2 + t13 ^ 2) * MDP(8); 0; MDP(2) + MDP(8); (-MDP(4) + t24) * t12 + (-MDP(5) + MDP(7) - t22 + t27) * t13; 0; MDP(3) + (-0.2e1 * t17 * MDP(10) + MDP(9) * t19) * t19 + ((-2 * MDP(6) + t28) * pkin(3)) + (0.2e1 * MDP(7) + 0.2e1 * t25 + 0.2e1 * t26 + t27) * qJ(4); t12 * MDP(8); 0; t24; MDP(8); t23 * t12; t22; (MDP(14) * t21 + MDP(11)) * t19 + (-MDP(15) * t21 - MDP(12)) * t17; t23; MDP(13);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
