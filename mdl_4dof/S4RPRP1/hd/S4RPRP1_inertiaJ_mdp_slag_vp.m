% Calculate joint inertia matrix for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPRP1_inertiaJ_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:29:48
% EndTime: 2019-03-08 18:29:48
% DurationCPUTime: 0.03s
% Computational Cost: add. (57->29), mult. (94->36), div. (0->0), fcn. (54->4), ass. (0->15)
t36 = 2 * pkin(3);
t35 = 2 * qJ(4);
t26 = sin(pkin(6));
t34 = pkin(1) * t26;
t27 = cos(pkin(6));
t25 = t27 * pkin(1) + pkin(2);
t28 = sin(qJ(3));
t29 = cos(qJ(3));
t33 = t28 * t25 + t29 * t34;
t20 = t29 * t25 - t28 * t34;
t32 = t20 * MDP(6);
t31 = t33 * MDP(7);
t19 = -pkin(3) - t20;
t18 = qJ(4) + t33;
t1 = [MDP(1) + MDP(5) + (t18 ^ 2 + t19 ^ 2) * MDP(10) + (t26 ^ 2 + t27 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * t32 - 0.2e1 * t31 - 0.2e1 * t19 * MDP(8) + 0.2e1 * t18 * MDP(9); 0; MDP(4) + MDP(10); MDP(5) + t32 - t31 + (t36 + t20) * MDP(8) + (t35 + t33) * MDP(9) + (-t19 * pkin(3) + t18 * qJ(4)) * MDP(10); 0; MDP(5) + MDP(8) * t36 + MDP(9) * t35 + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(10); t19 * MDP(10) - MDP(8); 0; -pkin(3) * MDP(10) - MDP(8); MDP(10);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
