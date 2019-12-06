% Calculate Gravitation load on the joints for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5PRRRR7_convert_par2_MPV_fixb.m
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRR7_gravloadJ_floatb_twist_mdp_slag_vp(qJ, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_gravloadJ_floatb_twist_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_gravloadJ_floatb_twist_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5PRRRR7_gravloadJ_floatb_twist_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From gravload_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:59
% EndTime: 2019-12-05 17:13:00
% DurationCPUTime: 0.16s
% Computational Cost: add. (179->42), mult. (199->68), div. (0->0), fcn. (202->10), ass. (0->21)
t35 = qJ(3) + qJ(4);
t34 = qJ(5) + t35;
t30 = sin(t34);
t31 = cos(t34);
t32 = sin(t35);
t33 = cos(t35);
t38 = sin(qJ(3));
t40 = cos(qJ(3));
t56 = t40 * MDP(10) - t38 * MDP(11) + t33 * MDP(17) - t32 * MDP(18) + t31 * MDP(24) - t30 * MDP(25) + MDP(3);
t39 = sin(qJ(2));
t55 = g(3) * t39;
t36 = sin(pkin(9));
t41 = cos(qJ(2));
t54 = t36 * t41;
t37 = cos(pkin(9));
t53 = t37 * t41;
t52 = t38 * t41;
t51 = t40 * t41;
t50 = (-g(1) * (-t30 * t53 + t36 * t31) - g(2) * (-t30 * t54 - t37 * t31) + t30 * t55) * MDP(24) + (-g(1) * (-t36 * t30 - t31 * t53) - g(2) * (t37 * t30 - t31 * t54) + t31 * t55) * MDP(25);
t43 = (-g(1) * (-t32 * t53 + t36 * t33) - g(2) * (-t32 * t54 - t37 * t33) + t32 * t55) * MDP(17) + (-g(1) * (-t36 * t32 - t33 * t53) - g(2) * (t37 * t32 - t33 * t54) + t33 * t55) * MDP(18) + t50;
t1 = [-g(3) * MDP(1); (t39 * MDP(4) - t56 * t41) * g(3) + (MDP(4) * t41 + t56 * t39) * (g(1) * t37 + g(2) * t36); (-g(1) * (t36 * t40 - t37 * t52) - g(2) * (-t36 * t52 - t37 * t40) + t38 * t55) * MDP(10) + (-g(1) * (-t36 * t38 - t37 * t51) - g(2) * (-t36 * t51 + t37 * t38) + t40 * t55) * MDP(11) + t43; t43; t50;];
taug = t1;
