% Calculate joint inertia matrix for
% S4RPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPPR4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S4RPPR4_inertiaJ_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:54
% EndTime: 2019-12-31 16:38:54
% DurationCPUTime: 0.04s
% Computational Cost: add. (28->19), mult. (47->26), div. (0->0), fcn. (24->4), ass. (0->11)
t16 = cos(pkin(6));
t14 = -t16 * pkin(1) - pkin(2);
t23 = t14 * MDP(7);
t17 = sin(qJ(4));
t22 = t17 * MDP(13);
t18 = cos(qJ(4));
t21 = t18 * MDP(14);
t20 = t18 * MDP(13) - t17 * MDP(14);
t15 = sin(pkin(6));
t12 = t15 * pkin(1) + qJ(3);
t1 = [MDP(1) + (t15 ^ 2 + t16 ^ 2) * MDP(4) * pkin(1) ^ 2 + (MDP(8) * t18 - 0.2e1 * t17 * MDP(9)) * t18 + ((2 * MDP(5)) + t23) * t14 + (MDP(7) * t12 + (2 * MDP(6)) + 0.2e1 * t21 + 0.2e1 * t22) * t12; 0; MDP(4) + MDP(7); MDP(5) + t23; 0; MDP(7); t18 * MDP(10) - t17 * MDP(11) + t20 * (-pkin(5) + t14); -t21 - t22; t20; MDP(12);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
