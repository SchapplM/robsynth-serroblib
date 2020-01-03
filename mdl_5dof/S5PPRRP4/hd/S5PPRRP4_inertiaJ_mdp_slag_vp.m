% Calculate joint inertia matrix for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP4_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_inertiaJ_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S5PPRRP4_inertiaJ_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:41
% EndTime: 2019-12-31 17:34:42
% DurationCPUTime: 0.05s
% Computational Cost: add. (47->32), mult. (97->48), div. (0->0), fcn. (67->4), ass. (0->17)
t34 = -qJ(5) - pkin(6);
t23 = sin(qJ(4));
t20 = t23 ^ 2;
t25 = cos(qJ(4));
t33 = t25 ^ 2 + t20;
t32 = MDP(14) * pkin(4);
t19 = -t25 * pkin(4) - pkin(3);
t31 = t19 * MDP(14);
t30 = t23 * MDP(12);
t29 = -MDP(11) - t32;
t17 = t34 * t23;
t18 = t34 * t25;
t28 = -t17 * t23 - t18 * t25;
t27 = t25 * MDP(11) - t30;
t26 = cos(qJ(3));
t24 = sin(qJ(3));
t1 = [t33 * MDP(14) + MDP(1) + MDP(2); 0; MDP(2) + (t33 * t24 ^ 2 + t26 ^ 2) * MDP(14); (-t25 * t17 + t23 * t18) * MDP(14); (MDP(4) + t27 - t31) * t26 + (t33 * MDP(13) + t28 * MDP(14) - MDP(5)) * t24; MDP(3) + t20 * MDP(6) + 0.2e1 * t23 * t25 * MDP(7) + 0.2e1 * t28 * MDP(13) + (t17 ^ 2 + t18 ^ 2 + t19 ^ 2) * MDP(14) + 0.2e1 * t27 * pkin(3); t29 * t25 + t30; (-MDP(12) * t25 + t29 * t23) * t24; t17 * t32 + (-MDP(12) * pkin(6) + MDP(9)) * t25 + (-MDP(11) * pkin(6) - MDP(13) * pkin(4) + MDP(8)) * t23; MDP(14) * pkin(4) ^ 2 + MDP(10); 0; -t26 * MDP(14); t31; 0; MDP(14);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
