% Calculate joint inertia matrix for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% MDP [12x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RRPP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [12 1]), ...
  'S4RRPP2_inertiaJ_mdp_slag_vp: MDP has to be [12x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:34:05
% EndTime: 2019-03-08 18:34:06
% DurationCPUTime: 0.05s
% Computational Cost: add. (63->40), mult. (81->42), div. (0->0), fcn. (22->2), ass. (0->18)
t40 = 2 * qJ(3);
t39 = 2 * MDP(7);
t38 = 2 * MDP(10);
t29 = sin(qJ(2));
t37 = t29 * MDP(6);
t36 = MDP(7) + MDP(10);
t35 = MDP(8) + MDP(11);
t30 = cos(qJ(2));
t26 = t30 * pkin(1) + pkin(2);
t34 = 0.2e1 * pkin(2);
t33 = qJ(3) ^ 2;
t31 = pkin(2) + pkin(3);
t28 = t29 * pkin(1);
t24 = t28 + qJ(3);
t22 = pkin(3) + t26;
t21 = t24 ^ 2;
t20 = t24 * qJ(3);
t1 = [MDP(1) + MDP(4) + (t26 ^ 2 + t21) * MDP(9) + (t22 ^ 2 + t21) * MDP(12) + 0.2e1 * (t30 * MDP(5) - t37) * pkin(1) + t22 * t38 + t26 * t39 + 0.2e1 * t35 * t24; MDP(4) + t34 * MDP(7) + (t26 * pkin(2) + t20) * MDP(9) + (0.2e1 * pkin(3) + t34) * MDP(10) + (t22 * t31 + t20) * MDP(12) + t35 * (t40 + t28) + (-t37 + (MDP(5) + t36) * t30) * pkin(1); MDP(4) + pkin(2) * t39 + (pkin(2) ^ 2 + t33) * MDP(9) + t31 * t38 + (t31 ^ 2 + t33) * MDP(12) + t35 * t40; -t22 * MDP(12) - t26 * MDP(9) - t36; -t31 * MDP(12) - pkin(2) * MDP(9) - t36; MDP(9) + MDP(12); 0; 0; 0; MDP(12);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
