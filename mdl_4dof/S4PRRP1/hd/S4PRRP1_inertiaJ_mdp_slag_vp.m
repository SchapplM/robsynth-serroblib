% Calculate joint inertia matrix for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4PRRP1_inertiaJ_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:23:04
% EndTime: 2019-03-08 18:23:04
% DurationCPUTime: 0.03s
% Computational Cost: add. (29->23), mult. (44->27), div. (0->0), fcn. (13->2), ass. (0->10)
t19 = 2 * MDP(8);
t18 = 2 * qJ(4);
t14 = sin(qJ(3));
t17 = t14 * MDP(7);
t16 = pkin(3) * t19 + MDP(5);
t15 = cos(qJ(3));
t13 = t14 * pkin(2);
t11 = t15 * pkin(2) + pkin(3);
t10 = t13 + qJ(4);
t1 = [MDP(1) + MDP(10); 0; MDP(2) + MDP(5) + (t10 ^ 2 + t11 ^ 2) * MDP(10) + 0.2e1 * (t15 * MDP(6) - t17) * pkin(2) + t11 * t19 + 0.2e1 * t10 * MDP(9); 0; (t18 + t13) * MDP(9) + (t11 * pkin(3) + t10 * qJ(4)) * MDP(10) + (-t17 + (MDP(6) + MDP(8)) * t15) * pkin(2) + t16; MDP(9) * t18 + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(10) + t16; 0; -t11 * MDP(10) - MDP(8); -pkin(3) * MDP(10) - MDP(8); MDP(10);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
