% Calculate joint inertia matrix for
% S4RPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S4RPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [4x4]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRP3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP3_inertiaJ_mdp_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP3_inertiaJ_mdp_slag_vp: pkin has to be [6x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S4RPRP3_inertiaJ_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:42:48
% EndTime: 2019-12-31 16:42:48
% DurationCPUTime: 0.04s
% Computational Cost: add. (44->26), mult. (79->43), div. (0->0), fcn. (50->4), ass. (0->14)
t29 = MDP(13) * pkin(3);
t21 = sin(pkin(6));
t18 = t21 * pkin(1) + pkin(5);
t28 = qJ(4) + t18;
t23 = sin(qJ(3));
t27 = t23 * MDP(11);
t22 = cos(pkin(6));
t26 = -t22 * pkin(1) - pkin(2);
t24 = cos(qJ(3));
t20 = t23 ^ 2;
t17 = -t24 * pkin(3) + t26;
t16 = t28 * t24;
t15 = t28 * t23;
t1 = [MDP(1) + t20 * MDP(5) + 0.2e1 * t23 * t24 * MDP(6) + 0.2e1 * (t15 * t23 + t16 * t24) * MDP(12) + (t15 ^ 2 + t16 ^ 2 + t17 ^ 2) * MDP(13) + (t21 ^ 2 + t22 ^ 2) * MDP(4) * pkin(1) ^ 2 + 0.2e1 * (-t24 * MDP(10) + t27) * t26; (-t15 * t24 + t16 * t23) * MDP(13); MDP(4) + (t24 ^ 2 + t20) * MDP(13); -t15 * t29 + (-MDP(11) * t18 + MDP(8)) * t24 + (-MDP(10) * t18 - MDP(12) * pkin(3) + MDP(7)) * t23; -t27 + (MDP(10) + t29) * t24; MDP(13) * pkin(3) ^ 2 + MDP(9); t17 * MDP(13); 0; 0; MDP(13);];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
