% Calculate joint inertia matrix for
% S4PRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
% MDP [8x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [8 1]), ...
  'S4PRRP2_inertiaJ_mdp_slag_vp: MDP has to be [8x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:24:10
% EndTime: 2019-03-08 18:24:10
% DurationCPUTime: 0.03s
% Computational Cost: add. (30->20), mult. (55->31), div. (0->0), fcn. (47->4), ass. (0->12)
t15 = sin(qJ(3));
t22 = pkin(2) * t15;
t21 = MDP(8) * pkin(3);
t16 = sin(qJ(2));
t17 = cos(qJ(3));
t18 = cos(qJ(2));
t12 = -t15 * t16 + t17 * t18;
t13 = t15 * t18 + t17 * t16;
t20 = t12 * MDP(6) - t13 * MDP(7);
t19 = t17 * MDP(6);
t14 = t17 * pkin(2) + pkin(3);
t1 = [MDP(1) + (t12 ^ 2 + t13 ^ 2) * MDP(8); t18 * MDP(3) - t16 * MDP(4) + (t12 * t14 + t13 * t22) * MDP(8) + t20; t14 ^ 2 * MDP(8) + MDP(2) + MDP(5) + (0.2e1 * t19 + (MDP(8) * t22 - 0.2e1 * MDP(7)) * t15) * pkin(2); t12 * t21 + t20; t14 * t21 + MDP(5) + (-t15 * MDP(7) + t19) * pkin(2); pkin(3) ^ 2 * MDP(8) + MDP(5); 0; 0; 0; MDP(8);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
