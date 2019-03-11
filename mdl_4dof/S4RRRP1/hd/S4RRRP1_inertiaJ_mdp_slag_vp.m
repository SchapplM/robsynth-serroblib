% Calculate joint inertia matrix for
% S4RRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RRRP1_inertiaJ_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:04
% EndTime: 2019-03-08 18:36:04
% DurationCPUTime: 0.05s
% Computational Cost: add. (57->34), mult. (97->44), div. (0->0), fcn. (59->4), ass. (0->20)
t24 = sin(qJ(2));
t34 = pkin(1) * t24;
t33 = MDP(10) * pkin(2);
t26 = cos(qJ(2));
t21 = t26 * pkin(1) + pkin(2);
t25 = cos(qJ(3));
t19 = t25 * t21;
t23 = sin(qJ(3));
t17 = -t23 * t34 + t19;
t32 = t17 * MDP(8);
t18 = t23 * t21 + t25 * t34;
t31 = t18 * MDP(9);
t30 = t25 * MDP(8);
t29 = t26 * MDP(5);
t22 = t25 * pkin(2);
t20 = t22 + pkin(3);
t28 = MDP(10) * t20;
t27 = MDP(4) + MDP(7);
t16 = pkin(3) + t17;
t1 = [MDP(1) + (t16 ^ 2 + t18 ^ 2) * MDP(10) + 0.2e1 * (-t24 * MDP(6) + t29) * pkin(1) + 0.2e1 * t32 - 0.2e1 * t31 + t27; (t19 + t22) * MDP(8) + t16 * t28 + ((-pkin(2) - t21) * MDP(9) + t18 * t33) * t23 + (t29 + (-MDP(8) * t23 - MDP(9) * t25 - MDP(6)) * t24) * pkin(1) + t27; t20 ^ 2 * MDP(10) + (0.2e1 * t30 + (t23 * t33 - 0.2e1 * MDP(9)) * t23) * pkin(2) + t27; t16 * pkin(3) * MDP(10) + MDP(7) - t31 + t32; pkin(3) * t28 + MDP(7) + (-t23 * MDP(9) + t30) * pkin(2); pkin(3) ^ 2 * MDP(10) + MDP(7); 0; 0; 0; MDP(10);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
