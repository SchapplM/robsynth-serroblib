% Calculate joint inertia matrix for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRRR1_inertiaJ_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:25
% EndTime: 2019-03-08 18:25:25
% DurationCPUTime: 0.03s
% Computational Cost: add. (32->22), mult. (56->25), div. (0->0), fcn. (34->4), ass. (0->14)
t19 = sin(qJ(3));
t28 = pkin(2) * t19;
t21 = cos(qJ(3));
t17 = t21 * pkin(2) + pkin(3);
t20 = cos(qJ(4));
t16 = t20 * t17;
t18 = sin(qJ(4));
t27 = (-t18 * t28 + t16) * MDP(9);
t26 = t21 * MDP(6);
t25 = (-t18 * t17 - t20 * t28) * MDP(10);
t24 = t18 * MDP(10);
t23 = MDP(5) + MDP(8);
t22 = (t20 * MDP(9) - t24) * pkin(3);
t1 = [MDP(1); 0; MDP(2) + 0.2e1 * (-t19 * MDP(7) + t26) * pkin(2) + 0.2e1 * t25 + 0.2e1 * t27 + t23; 0; (t20 * pkin(3) + t16) * MDP(9) + (-pkin(3) - t17) * t24 + (t26 + (-MDP(10) * t20 - MDP(9) * t18 - MDP(7)) * t19) * pkin(2) + t23; 0.2e1 * t22 + t23; 0; MDP(8) + t25 + t27; MDP(8) + t22; MDP(8);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
