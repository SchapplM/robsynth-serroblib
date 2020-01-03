% Calculate joint inertia matrix for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP6_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP6_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(5,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_inertiaJ_mdp_slag_vp: pkin has to be [5x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4RPRP6_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:14
% EndTime: 2019-12-31 16:46:14
% DurationCPUTime: 0.04s
% Computational Cost: add. (46->29), mult. (70->42), div. (0->0), fcn. (37->2), ass. (0->13)
t24 = (pkin(1) * MDP(6));
t23 = MDP(15) * pkin(3);
t21 = -pkin(1) - pkin(5);
t22 = -qJ(4) + t21;
t20 = cos(qJ(3));
t19 = sin(qJ(3));
t18 = t20 ^ 2;
t17 = t19 * pkin(3) + qJ(2);
t16 = t19 ^ 2 + t18;
t15 = t22 * t20;
t14 = t22 * t19;
t13 = t14 * t19 + t15 * t20;
t1 = [MDP(1) + t18 * MDP(7) - 0.2e1 * t20 * t19 * MDP(8) - 0.2e1 * t13 * MDP(14) + (t14 ^ 2 + t15 ^ 2 + t17 ^ 2) * MDP(15) + ((-2 * MDP(4) + t24) * pkin(1)) + (0.2e1 * t19 * MDP(12) + 0.2e1 * t20 * MDP(13) + MDP(6) * qJ(2) + (2 * MDP(5))) * qJ(2); -t16 * MDP(14) + t13 * MDP(15) + MDP(4) - t24; t16 * MDP(15) + MDP(6); t15 * t23 + (-MDP(13) * t21 - MDP(10)) * t19 + (MDP(12) * t21 - MDP(14) * pkin(3) + MDP(9)) * t20; -t19 * MDP(13) + (MDP(12) + t23) * t20; MDP(15) * pkin(3) ^ 2 + MDP(11); t17 * MDP(15); 0; 0; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
