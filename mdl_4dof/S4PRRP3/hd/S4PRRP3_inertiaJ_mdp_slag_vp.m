% Calculate joint inertia matrix for
% S4PRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4PRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S4PRRP3_inertiaJ_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:56
% EndTime: 2019-12-31 16:26:56
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->22), mult. (59->36), div. (0->0), fcn. (35->2), ass. (0->10)
t18 = -qJ(4) - pkin(5);
t17 = MDP(13) * pkin(3);
t14 = sin(qJ(3));
t16 = t14 * MDP(11);
t15 = cos(qJ(3));
t13 = t14 ^ 2;
t12 = -t15 * pkin(3) - pkin(2);
t11 = t18 * t15;
t10 = t18 * t14;
t1 = [MDP(1) + (t15 ^ 2 + t13) * MDP(13); (t15 * t10 - t14 * t11) * MDP(13); MDP(2) + t13 * MDP(5) + 0.2e1 * t14 * t15 * MDP(6) + 0.2e1 * (-t10 * t14 - t11 * t15) * MDP(12) + (t10 ^ 2 + t11 ^ 2 + t12 ^ 2) * MDP(13) + 0.2e1 * (t15 * MDP(10) - t16) * pkin(2); -t16 + (MDP(10) + t17) * t15; t10 * t17 + (-MDP(11) * pkin(5) + MDP(8)) * t15 + (-MDP(10) * pkin(5) - MDP(12) * pkin(3) + MDP(7)) * t14; MDP(13) * pkin(3) ^ 2 + MDP(9); 0; t12 * MDP(13); 0; MDP(13);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
