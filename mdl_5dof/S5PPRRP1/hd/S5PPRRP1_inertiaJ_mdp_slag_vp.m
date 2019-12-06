% Calculate joint inertia matrix for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP1_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_inertiaJ_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S5PPRRP1_inertiaJ_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:14
% EndTime: 2019-12-05 15:07:15
% DurationCPUTime: 0.06s
% Computational Cost: add. (63->34), mult. (132->55), div. (0->0), fcn. (117->6), ass. (0->21)
t42 = -qJ(5) - pkin(6);
t31 = sin(qJ(4));
t27 = t31 ^ 2;
t33 = cos(qJ(4));
t41 = t33 ^ 2 + t27;
t40 = MDP(14) * pkin(4);
t26 = -t33 * pkin(4) - pkin(3);
t39 = t26 * MDP(14);
t38 = t31 * MDP(12);
t37 = MDP(11) + t40;
t24 = t42 * t31;
t25 = t42 * t33;
t36 = -t24 * t31 - t25 * t33;
t35 = t33 * MDP(11) - t38;
t34 = cos(qJ(3));
t32 = sin(qJ(3));
t30 = cos(pkin(8));
t29 = sin(pkin(8));
t23 = t34 * t29 + t32 * t30;
t22 = t32 * t29 - t34 * t30;
t1 = [MDP(1) + (t29 ^ 2 + t30 ^ 2) * MDP(2) + (t41 * t23 ^ 2 + t22 ^ 2) * MDP(14); 0; t41 * MDP(14) + MDP(2); (-MDP(4) - t35 + t39) * t22 + (t41 * MDP(13) + t36 * MDP(14) - MDP(5)) * t23; (t33 * t24 - t31 * t25) * MDP(14); MDP(3) + t27 * MDP(6) + 0.2e1 * t31 * t33 * MDP(7) + 0.2e1 * t36 * MDP(13) + (t24 ^ 2 + t25 ^ 2 + t26 ^ 2) * MDP(14) + 0.2e1 * t35 * pkin(3); (-MDP(12) * t33 - t37 * t31) * t23; t37 * t33 - t38; t24 * t40 + (-MDP(12) * pkin(6) + MDP(9)) * t33 + (-MDP(11) * pkin(6) - MDP(13) * pkin(4) + MDP(8)) * t31; MDP(14) * pkin(4) ^ 2 + MDP(10); t22 * MDP(14); 0; t39; 0; MDP(14);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
