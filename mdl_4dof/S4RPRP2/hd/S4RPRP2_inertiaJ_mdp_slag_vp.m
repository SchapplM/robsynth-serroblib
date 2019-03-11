% Calculate joint inertia matrix for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% MDP [10x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP2_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [10 1]), ...
  'S4RPRP2_inertiaJ_mdp_slag_vp: MDP has to be [10x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:53
% EndTime: 2019-03-08 18:30:53
% DurationCPUTime: 0.03s
% Computational Cost: add. (43->26), mult. (53->32), div. (0->0), fcn. (27->2), ass. (0->11)
t12 = sin(qJ(3));
t13 = cos(qJ(3));
t14 = -pkin(1) - pkin(2);
t9 = t12 * qJ(2) - t13 * t14;
t18 = t9 * MDP(8);
t17 = MDP(10) * pkin(3);
t10 = t13 * qJ(2) + t12 * t14;
t16 = t10 * MDP(9);
t15 = t12 * MDP(9);
t8 = -pkin(3) - t9;
t1 = [MDP(1) + (2 * pkin(1) * MDP(4)) + 0.2e1 * qJ(2) * MDP(5) + ((pkin(1) ^ 2) + qJ(2) ^ 2) * MDP(6) + MDP(7) + (t10 ^ 2 + t8 ^ 2) * MDP(10) + 0.2e1 * t18 + 0.2e1 * t16; -MDP(4) - pkin(1) * MDP(6) - t13 * MDP(8) + t15 + (t10 * t12 + t8 * t13) * MDP(10); MDP(6) + (t12 ^ 2 + t13 ^ 2) * MDP(10); t8 * t17 - MDP(7) - t16 - t18; -t15 + (MDP(8) + t17) * t13; pkin(3) ^ 2 * MDP(10) + MDP(7); 0; 0; 0; MDP(10);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1) t1(2) t1(4) t1(7); t1(2) t1(3) t1(5) t1(8); t1(4) t1(5) t1(6) t1(9); t1(7) t1(8) t1(9) t1(10);];
Mq  = res;
