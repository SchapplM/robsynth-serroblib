% Calculate joint inertia matrix for
% S5PPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% MDP [13x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PPRPR3_convert_par2_MPV_fixb.m
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR3_inertiaJ_mdp_slag_vp(qJ, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_inertiaJ_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_inertiaJ_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [13 1]), ...
  'S5PPRPR3_inertiaJ_mdp_slag_vp: MDP has to be [13x1] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:05:22
% EndTime: 2019-12-05 15:05:23
% DurationCPUTime: 0.06s
% Computational Cost: add. (57->34), mult. (138->62), div. (0->0), fcn. (143->8), ass. (0->19)
t43 = MDP(6) * pkin(3);
t37 = cos(qJ(5));
t42 = t37 * MDP(12);
t31 = sin(pkin(9));
t33 = cos(pkin(9));
t36 = sin(qJ(3));
t38 = cos(qJ(3));
t26 = t31 * t38 + t33 * t36;
t24 = t31 * t36 - t33 * t38;
t35 = sin(qJ(5));
t41 = -t35 * MDP(13) + t42;
t40 = -t35 * MDP(12) - t37 * MDP(13);
t34 = cos(pkin(8));
t32 = sin(pkin(8));
t30 = t34 ^ 2;
t29 = -t33 * pkin(3) - pkin(4);
t23 = t24 * t32;
t21 = t26 * t32;
t1 = [MDP(1) + (t32 ^ 2 + t30) * MDP(2) + (t21 ^ 2 + t23 ^ 2 + t30) * MDP(6); (t21 * t24 - t23 * t26) * MDP(6); MDP(2) + (t24 ^ 2 + t26 ^ 2) * MDP(6); (-t36 * MDP(4) - t38 * MDP(5)) * t32 - t41 * t21 + (-t21 * t33 - t23 * t31) * t43; t38 * MDP(4) - t36 * MDP(5) - t41 * t24 + (-t24 * t33 + t26 * t31) * t43; -0.2e1 * t29 * t42 + MDP(3) + (t31 ^ 2 + t33 ^ 2) * MDP(6) * pkin(3) ^ 2 + (0.2e1 * t29 * MDP(13) + MDP(7) * t35 + 0.2e1 * t37 * MDP(8)) * t35; -t34 * MDP(6); 0; 0; MDP(6); (t35 * t23 - t34 * t37) * MDP(12) + (t37 * t23 + t34 * t35) * MDP(13); t40 * t26; t37 * MDP(10) + t35 * MDP(9) + t40 * (t31 * pkin(3) + pkin(6)); t41; MDP(11);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
