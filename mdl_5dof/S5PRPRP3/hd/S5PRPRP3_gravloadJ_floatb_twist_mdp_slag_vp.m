% Calculate Gravitation load on the joints for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% MDP [14x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRPRP3_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [14 1]), ...
  'S5PRPRP3_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [14x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:41
% EndTime: 2019-12-05 15:33:42
% DurationCPUTime: 0.17s
% Computational Cost: add. (80->33), mult. (120->49), div. (0->0), fcn. (99->8), ass. (0->20)
t32 = sin(pkin(7));
t33 = cos(pkin(7));
t43 = g(1) * t33 + g(2) * t32;
t31 = qJ(2) + pkin(8);
t29 = sin(t31);
t30 = cos(t31);
t53 = -g(3) * t30 + t29 * t43;
t50 = g(3) * t29;
t35 = sin(qJ(4));
t48 = t32 * t35;
t37 = cos(qJ(4));
t47 = t32 * t37;
t46 = t33 * t35;
t45 = t33 * t37;
t44 = MDP(14) + MDP(5);
t38 = cos(qJ(2));
t36 = sin(qJ(2));
t34 = -qJ(5) - pkin(6);
t28 = pkin(4) * t37 + pkin(3);
t1 = [(-MDP(1) - t44) * g(3); (g(3) * t36 + t38 * t43) * MDP(4) + (-t30 * t43 - t50) * MDP(13) + (-g(3) * (pkin(2) * t38 + t28 * t30 - t29 * t34) + t43 * (pkin(2) * t36 + t28 * t29 + t30 * t34)) * MDP(14) + (MDP(5) * pkin(2) + MDP(3)) * (-g(3) * t38 + t36 * t43) + (MDP(11) * t37 - MDP(12) * t35) * t53; t44 * (-g(1) * t32 + g(2) * t33); (-g(1) * (-t30 * t45 - t48) - g(2) * (-t30 * t47 + t46) + t37 * t50) * MDP(12) + (MDP(14) * pkin(4) + MDP(11)) * (-g(1) * (-t30 * t46 + t47) - g(2) * (-t30 * t48 - t45) + t35 * t50); -t53 * MDP(14);];
taug = t1;
