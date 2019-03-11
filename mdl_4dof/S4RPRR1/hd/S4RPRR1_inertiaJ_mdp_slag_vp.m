% Calculate joint inertia matrix for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPRR1_inertiaJ_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:32:00
% EndTime: 2019-03-08 18:32:00
% DurationCPUTime: 0.04s
% Computational Cost: add. (66->26), mult. (113->31), div. (0->0), fcn. (84->6), ass. (0->20)
t28 = sin(pkin(7));
t42 = pkin(1) * t28;
t29 = cos(pkin(7));
t27 = t29 * pkin(1) + pkin(2);
t31 = sin(qJ(3));
t33 = cos(qJ(3));
t25 = t31 * t27 + t33 * t42;
t32 = cos(qJ(4));
t41 = t32 * t25;
t24 = t33 * t27 - t31 * t42;
t23 = pkin(3) + t24;
t30 = sin(qJ(4));
t20 = t32 * t23 - t30 * t25;
t40 = t20 * MDP(9);
t39 = t24 * MDP(6);
t38 = t25 * MDP(7);
t37 = (-t30 * t23 - t41) * MDP(10);
t36 = MDP(5) + MDP(8);
t35 = (-t30 * MDP(10) + t32 * MDP(9)) * pkin(3);
t1 = [MDP(1) + (t28 ^ 2 + t29 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * t37 + 0.2e1 * t39 - 0.2e1 * t38 + 0.2e1 * t40 + t36; 0; MDP(4); t39 - t38 + (t32 * pkin(3) + t20) * MDP(9) + (-t41 + (-pkin(3) - t23) * t30) * MDP(10) + t36; 0; 0.2e1 * t35 + t36; MDP(8) + t37 + t40; 0; MDP(8) + t35; MDP(8);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
