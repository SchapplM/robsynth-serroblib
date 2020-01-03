% Calculate joint inertia matrix for
% S4PRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [15x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(15,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP4_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP4_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'S4PRRP4_inertiaJ_mdp_slag_vp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:27:58
% EndTime: 2019-12-31 16:27:58
% DurationCPUTime: 0.05s
% Computational Cost: add. (40->29), mult. (78->44), div. (0->0), fcn. (38->2), ass. (0->9)
t17 = sin(qJ(3));
t15 = t17 ^ 2;
t18 = cos(qJ(3));
t23 = t18 ^ 2 + t15;
t22 = -pkin(3) * MDP(15) - MDP(12);
t21 = MDP(10) - t22;
t20 = MDP(15) * qJ(4) - MDP(11) + MDP(14);
t14 = -t18 * pkin(3) - t17 * qJ(4) - pkin(2);
t1 = [t23 * MDP(15) + MDP(1); 0; MDP(2) + t15 * MDP(5) + (t23 * pkin(5) ^ 2 + t14 ^ 2) * MDP(15) + 0.2e1 * t23 * MDP(13) * pkin(5) + 0.2e1 * (pkin(2) * MDP(10) - t14 * MDP(12)) * t18 + 0.2e1 * (-pkin(2) * MDP(11) - t14 * MDP(14) + t18 * MDP(6)) * t17; t20 * t17 + t21 * t18; t17 * MDP(7) + t18 * MDP(8) + (-t17 * pkin(3) + t18 * qJ(4)) * MDP(13) + (-t21 * t17 + t20 * t18) * pkin(5); MDP(9) + 0.2e1 * pkin(3) * MDP(12) + 0.2e1 * qJ(4) * MDP(14) + (pkin(3) ^ 2 + qJ(4) ^ 2) * MDP(15); -t18 * MDP(15); (MDP(15) * pkin(5) + MDP(13)) * t17; t22; MDP(15);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
